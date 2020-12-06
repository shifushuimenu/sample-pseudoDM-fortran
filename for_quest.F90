module types 
implicit none 

integer, parameter :: dp = 8
integer, parameter :: ZERO = 0

end module types 


module square_lattice_FT
    use types 
    implicit none 
    private 
    public zexpiqr
    public lq, nk, init
    public listk 

    public calc_FT_coeff 
    public FT_GreensFunction

    real(dp), parameter :: pi = dacos(-1.d0)
    integer :: l
    integer :: a1_p(2), a2_p(2)
    real(dp) :: b1_p(2), b2_p(2)
    integer, allocatable :: listk(:,:), list(:,:)
    integer, allocatable :: invlistk(:,:)

    integer :: lq, nk
    complex(dp), allocatable :: zexpiqr(:,:)
    ! Are the FT coeff's initialized ? 
    logical, save :: init = .false.           

    contains 
    subroutine calc_FT_coeff(Nsites)
        implicit none 
        ! Purpose:
        ! --------
        ! Calculate the coefficients of the Fourier transform 
        !     C(i,j) = exp(i vec{q}_i * \vec{r}_j)
        !
        ! Arguments:
        ! ----------
        integer, intent(in) :: Nsites  ! Number of sites of square lattice 

        ! ... local variables ...
        integer :: i, j, ir, ix, iy, iq
        real(dp) :: qvec(2), ri(2)

        lq = Nsites 
        l = nint( sqrt(dble(lq)) )
        nk = l*l ! (l+1)*(l+1)
        a1_p(1) = 1 ; a1_p(2) = 0
        a2_p(1) = 0 ; a2_p(2) = 1
        b1_p(1) = 2.0_dp*pi/dble(l) ; b1_p(2) = 0.0_dp
        b2_p(1) = 0.0_dp            ; b2_p(2) = 2.0_dp*pi/dble(l)

        allocate( zexpiqr(lq,nk) )
        allocate( listk(nk,2), list(lq,2) )
        allocate( invlistk( 1-l/2 : l/2, 1-l/2 : l/2 ) )

        ir = 0
        do ix = 1,l
            do iy = 1,l 
                ir = ir + 1
                list(ir,1) = ix
                list(ir,2) = iy
            enddo
        enddo

        nk = 0
        do j = 1,l ! 0, l
            do i = 1,l ! 0, l 
                nk = nk + 1
                listk(nk,1) = i - l/2 
                listk(nk,2) = j - l/2
                invlistk(listk(nk,1), listk(nk,2)) = nk
            enddo 
        enddo 

        do iq = 1, nk
            qvec = dble(listk(iq,1))*b1_p + dble(listk(iq,2))*b2_p
            do i = 1, lq
                ri = dble(list(i,1))*a1_p + dble(list(i,2))*a2_p
                zexpiqr(i,iq) = exp( dcmplx( 0.0_dp, qvec(1)*ri(1) + qvec(2)*ri(2) )  ) / sqrt(dble(lq))
            enddo 
        enddo 

        init = .true.

    end subroutine calc_FT_coeff

    subroutine FT_GreensFunction(G, Gk)
        implicit none 
        ! Purpose:
        ! --------
        ! Fourier transform the Green's function from real to momentum space. 
        !
        ! Arguments:
        ! ----------
        real(dp), intent(in) :: G(:,:)       ! Green's function in real space 
        complex(dp), intent(out) :: Gk(:,:)  ! Green's function in momentum space 

        ! ... Local Variables ...
        integer :: i, j, k1, k2 

        ! ... executable ...
        if( size(G,1) /= size(G,2) ) then 
            print*, "Error: Green's function should be a square matrix."
            print*, "Exiting ..."
            stop 
        endif 

        print*, init        
        if (.not. init) call calc_FT_coeff(size(G,1))
        print*, init, lq, nk

        Gk = ZERO
        do k2 = 1, nk
            do k1 = 1, nk
                do j = 1, lq 
                    do i = 1, lq
                        Gk(k1, k2) = Gk(k1, k2) + zexpiqr(i,k1) * conjg(zexpiqr(j,k2)) * G(i,j)
                    enddo                     
                enddo 
            enddo 
        enddo

    end subroutine FT_GreensFunction

end module square_lattice_FT


module direct_sampling_pseudo_DM
! --------------------------------------------------------------------------------        
! Componentwise direct sampling of occupation number pseudo-snapshots 
! from a pseudo-density matrix as it appears naturally 
! in finite-temperature determinantal QMC simulations. 
! 
! Reference:
! ==========
!    Stephan Humeniuk and Yuan Wan: "Numerically exact quantum gas microscopy 
!        of interacting lattice fermions", arXiv:2009.07377 (2020)
!
! 
! Usage:
! ------
! The subroutine `run_sample_pseudo_DM` should be called in the 
! measurement part of the determinantal QMC (DQMC). 
! `Nsamples_per_HS` pseudo-snapshots per DQMC Monte Carlo step are output
! into separate files for spin up and spin down. These files are understood to
! be "synchronized", i.e. the snapshot for spin-up in the m-th line
! of the file 
!         `outfile_basename`//"ncpu"//`MPI_rank`_up.dat
! is associated with the snapshot for spin-dn in the m-th line of 
!         `outfile_basename`//"ncpu"//`MPI_rank`_dn.dat
! 
! In order to deal with the case that the simulated Hamiltonian exhibits itself 
! a sign problem, the sign of the Monte Carlo weight for spin up and spin down 
! should be provided by the calling DQMC outer loop so that it can be multiplied 
! to the signed reweighting factors of the corresponding pseudo-snapshots.
! 
!--------------------------------------------------------------------------------
    use types
    use square_lattice_FT
    implicit none 
    private 
    public run_sample_pseudo_DM

    contains

subroutine run_sample_pseudo_DM(G_up, G_dn, BSS_sign_up, BSS_sign_dn, sp_basis, &
    Nsamples_per_HS, outfile_basename, MPI_rank)
! Purpose:
! --------
!   This (driver) routine calls `sample_FF_GreensFunction` for both spin components 
!   and makes sure that the snapshots, signs and reweighting factors for both spin 
!   components are output into two files (one for `up` and one for `dn`) in a synchronized way.
!   Care should be taken that the snapshots for the different spin components are 
!   only combined if they come from the same Hubbard-Stratonovich (HS) sample since only in
!   that case they are statistically independent. In other words, in one HS sample 
!   the samples for different spin components can be combined arbitrarily (but across 
!   different HS samples they must not be combined).

! Arguments:
! ----------
real(dp), intent(in) :: G_up(:,:)      ! Green's function for spin up
real(dp), intent(in) :: G_dn(:,:)      ! Green's function for spin down 
real(dp), intent(in) :: BSS_sign_up    ! sign of the weight of the HS sample for spin up
real(dp), intent(in) :: BSS_sign_dn    ! sign of the weight of the HS sample for spin dn
character(len=*), intent(in) :: sp_basis  ! single-particle basis \in ["real_space", "momentum_space"]
integer, intent(in) :: Nsamples_per_HS    ! number of occupation number states generated per HS sample (10-100 is a good choice)
character(len=*), intent(in) :: outfile_basename  ! base name of the output file for snapshots
integer, intent(in) :: MPI_rank           ! MPI rank for labelling output files 

! ... Local variables ...
integer :: Nsites
integer :: Nx
integer :: Ny

integer :: i, sss, spin_idx
complex(dp) :: Green_kspace(size(G_up, dim=1), size(G_up, dim=1))
real(dp) :: Green_xspace(size(G_up, dim=1), size(G_up, dim=1))
integer  :: occ_vector(size(G_up, dim=1))
complex(dp) :: weight_phase
real(dp) :: weight_factor
real(dp) :: abs_corr(size(G_up, dim=1))
integer :: Ksites(size(G_up, dim=1))

integer, parameter :: Nspin_species = 2
real(dp) :: BSS_sign(1:Nspin_species)
character(len=3)  :: chr_spin(1:Nspin_species)
character(len=5)  :: chr_rank

logical, save :: init_outfiles = .false.

integer :: occ_vector_tmp(size(G_up, dim=1), Nsamples_per_HS)
complex(dp) :: weight_phase_tmp(Nsamples_per_HS)
real(dp) :: weight_factor_tmp(Nsamples_per_HS)

BSS_sign = (/ BSS_sign_up, BSS_sign_dn /)

write(chr_rank, "(i5.5)") MPI_rank
chr_spin(1) = "_up"
chr_spin(2) = "_dn"

! Check that MPI is working correctly 
print*, "I am rank", MPI_rank

! assuming a square system
Nsites = size(G_up, dim=1)
Nx = int(sqrt(float(Nsites)))
Ny = Nx

do spin_idx = 1, 2
    occ_vector_tmp(:,:) = ZERO
    weight_factor_tmp(:) = ZERO
    weight_phase_tmp(:) = cmplx(ZERO, ZERO, kind=dp)

    if (spin_idx == 1) Green_xspace(:,:) = G_up(:,:)
    if (spin_idx == 2) Green_xspace(:,:) = G_dn(:,:)

    if (trim(sp_basis) == "momentum_space") then 
        call FT_GreensFunction(Green_xspace, Green_kspace)
    elseif (trim(sp_basis) == "real_space") then 
        Green_kspace(:,:) = cmplx(x=Green_xspace(:,:), kind=dp)
    else 
        print*, "Error run_sample_pseudo_DM: unknown single-particle basis ", trim(sp_basis) 
        stop
    endif 

    do sss = 1, Nsamples_per_HS
        call sample_FF_GreensFunction(Green_kspace, occ_vector, abs_corr, Ksites, weight_phase, weight_factor)
        occ_vector_tmp(:, sss) = occ_vector(:)
        weight_factor_tmp(sss) = weight_factor
        weight_phase_tmp(sss)  = weight_phase
    enddo 


    open(100, file=trim(outfile_basename)//"_ncpu"//chr_rank//chr_spin(spin_idx)//".dat", status="unknown", position="append")
    do sss = 1, Nsamples_per_HS
        if (.not.init_outfiles) then 
        write(100, *) "Pseudo-snapshots for spin"//chr_spin(spin_idx)//":"
        write(100, *) "BSS_sign    real(exp(i*sampling_phase))       aimag(exp(i*sampling_phase))     reweighting_factor    occ_vector(1:Nsites)"
        write(100, *) "========================================================================================================================="
        init_outfiles = .true.
        endif 
        write(100, '(4(f24.8, 6x), *(i3))')  BSS_sign(spin_idx), real(weight_phase_tmp(sss)), aimag(weight_phase_tmp(sss)), &
                        weight_factor_tmp(sss), ( occ_vector_tmp(i, sss), i = 1, Nsites )
    enddo 
    close(100)

enddo 

end subroutine            


subroutine sample_FF_GreensFunction(G, occ_vector, abs_corr, Ksites, &
        reweighting_phase, reweighting_factor)
    ! Purpose:
    ! --------
    ! Sample pseudo-snapshots from a free-fermion pseudo-density matrix using 
    ! the single-particle Green's function for a given spin species as input. 
    ! 
    ! In general, the single-particle Green's function may be complex as is the case 
    ! in the momentum basis. In that case there is a phase problem and 
    ! every snapshot is associated a complex phase `exp(i phi) = a + 1j b` 
    ! and a reweighting factor. 
    ! 
    ! Irrespective of the factor ordering during the componentwise sampling, 
    ! the sites in the output vector of occupation numbers are ordered in the same 
    ! manner as the sites in the input Green's function matrix. 

    implicit none 
    ! Arguments:
    ! ----------
    complex(dp), intent(in)  :: G(:,:)              ! equal-time Green's function matrix. 
    integer, intent(out)     :: occ_vector(:)       ! vector of occupation numbers for a given spin species 
    real(dp), intent(out)    :: abs_corr(:)
    complex(dp), intent(out) :: reweighting_phase   ! "complex sign" of the snapshot 
    real(dp), intent(out)    :: reweighting_factor        
    integer, intent(inout)   :: Ksites(:)           ! Factor ordering of the components. 

    ! ... Local Variables ...
    integer :: D  ! dimension of the single-particle Hilbert space 
    complex(dp), allocatable :: corr(:)  ! correction term taking care of correlations between different sites 
    complex(dp), allocatable :: cond_prob(:)   ! conditional probability 
    real(dp), allocatable :: q_prob(:)         ! absolute value of the conditional probability 
    real(dp), allocatable :: phase_angle(:)    ! phase angle of the conditional probability
    complex(dp), allocatable :: phase(:)       ! exp(i phi(:))
    real(dp) :: norm
    real(dp) :: eta 
    integer :: k, occ
       
    complex(dp), allocatable :: Xinv(:,:)
    complex(dp), allocatable :: Xinv_new(:,:)
    complex(dp) :: gg
    complex(dp), allocatable :: uu(:), vv(:)
    integer :: ii, jj

    ! check inputs 
    D = size(G, 1)
    if( size(G,1) /= size(G,2) ) then 
        print*, "ERROR: Green's function must be square."
    endif 
    if( size(G,1) /= size(occ_vector,1) ) then 
        print*, "ERROR: size(G,1) /= size(occ_vector,1)"
    endif 

    allocate(corr(D))
    allocate(cond_prob(0:1))
    allocate(q_prob(0:1))
    allocate(phase_angle(0:1))
    allocate(phase(0:1))

    ! lexicographic factor ordering 
    do ii = 1, D 
        Ksites(ii) = ii 
    enddo 

    ! helper variables 
    allocate(Xinv(1:D,1:D))
    allocate(Xinv_new(1:D,1:D))
    allocate(uu(1:D))
    allocate(vv(1:D))

    corr = cmplx(ZERO, ZERO)
    cond_prob = cmplx(ZERO, ZERO)
    reweighting_factor = 1.0_dp
    reweighting_phase = cmplx(1.0_dp, 0.0_dp)

    occ_vector(:) = -1 ! initialize to invalid value 
    Xinv = cmplx(ZERO, ZERO)
    Xinv_new = cmplx(ZERO, ZERO)

    ! Componentwise direct sampling
    do k = 1, D
        ! "Correction term" due to correlations between components, i.e. sites
        if ( k == 1 ) then 
            corr(1) = cmplx(ZERO, ZERO)
        elseif ( k == 2) then
            ! matmul() does not work for scalars 
            corr(k) = G(Ksites(k), Ksites(1)) * Xinv(1,1) * G(Ksites(1), Ksites(k))
        else
            corr(k) = dot_product( G(Ksites(k), Ksites(1:k-1)), matmul( Xinv(1:k-1, 1:k-1), G(Ksites(1:k-1), Ksites(k)) ) )
        endif 
    
        cond_prob(1) = 1.0_dp - G(Ksites(k), Ksites(k)) + corr(k)
        cond_prob(0) = G(Ksites(k), Ksites(k)) - corr(k)


        ! Take care of quasi-probability distribution 
        do ii = 0,1
            phase_angle(ii) = atan2(real(cond_prob(ii)), aimag(cond_prob(ii)))
            phase(ii) = cmplx(cos(phase_angle(ii)), sin(phase_angle(ii)))
        enddo 
        
        ! Reweighting 
        q_prob(0:1) = abs(cond_prob(0:1))
        norm = sum(q_prob(0:1))
        q_prob(0:1) = q_prob(0:1) / norm

        call random_number(eta)
        if( eta < q_prob(1) ) then 
            occ = 1 
            reweighting_phase = reweighting_phase * phase(1)
        else
            occ = 0
            reweighting_phase = reweighting_phase * phase(0)
        endif 
        reweighting_factor = reweighting_factor * norm 

        occ_vector(Ksites(k)) = occ 

        if( k==1 ) then 
            Xinv(1,1) = 1.0_dp / (G(Ksites(1), Ksites(1)) - occ)
        else
            ! Low-rank update 
            gg = 1.0_dp / (G(Ksites(k), Ksites(k)) - occ - corr(k))
            if ( k==2 ) then 
                ! matmul() does not work for scalars 
                uu(1) = Xinv(1,1) * G(Ksites(1), Ksites(2))
                vv(1) = G(Ksites(2), Ksites(1)) * Xinv(1,1)
            else
                uu(1:k-1) = matmul( Xinv(1:k-1,1:k-1), G(Ksites(1:k-1), Ksites(k)) )
                vv(1:k-1) = matmul( G(Ksites(k), Ksites(1:k-1)), Xinv(1:k-1,1:k-1) )
            endif 
            do jj = 1, k-1
                do ii = 1, k-1
                    Xinv_new(ii, jj) = Xinv(ii, jj) &
                        + gg * uu(ii) * vv(jj)
                enddo 
            enddo 
            do ii = 1, k-1
                Xinv_new(k, ii) = -gg*vv(ii)
                Xinv_new(ii, k) = -gg*uu(ii)
            enddo
            Xinv_new(k,k) = gg
            Xinv(1:k, 1:k) = Xinv_new(1:k, 1:k)

        endif 
    enddo 

    abs_corr(1:D) = abs(corr(1:D))

end subroutine 

end module  direct_sampling_pseudo_DM
