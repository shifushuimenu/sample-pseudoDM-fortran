
module types
    implicit none 

    integer, parameter :: dp = 8
    integer, parameter :: ZERO = 0

end module types 

module parallelization
    implicit none 

    integer :: MPI_rank, MPI_size
    integer :: ierr
    character(len=5)  :: chr_rank
    integer, parameter :: root_rank = 0

end module parallelization 

module util
    use types
    implicit none 

    contains 

    subroutine rotate(A, l)
    ! Purpose:
    ! --------
    !    Cyclically left-shift the elements of integer array 'A'
    !    by 'l' positions. 
    ! Arguments:
    ! ----------
        integer, intent(inout) :: A(:)  
        integer, intent(in)  :: l
    ! ... Local variables ...
        integer :: ii, n, s, temp(size(A,1))
    ! ... Executable ...
        n = size(A,1)
        s = mod(l, n)
        temp(1:s) = A(1:s)
        do ii = 1, n-s
            A(ii) = A(ii+s)
        enddo 
        A(n-s+1:n) = temp(1:s)
    end subroutine rotate

    subroutine random_permutation(A)
    ! Purpose:
    ! --------
    !   Randomly permute the elements in the integer array A(:)
    !   (in place).
    !
    ! Arguments:
    ! ----------
        integer, intent(inout) :: A(:)
    ! ... Local variables ...
        integer :: n, k, l
        real(dp) :: eta
    ! ... Executable ...
        n = size(A,1)
        do k = n, 2, -1
            call random_number(eta)
            l = ceiling(eta * (k-1))
            call rotate(A(1:k), l)
        enddo 
    end subroutine 

    subroutine init_RNG(MPI_rank)
    ! Purpose:
    ! --------
    !    Initialize the standard pseudo-random number generator 
    !    with a seed obtained from the system time (at the millisecond level)
    !    and the MPI rank.
    !    Call the random number generator:
    !
    !        integer :: eta
    !        call random_number(eta)
    !
    ! Arguments:
    ! ----------
    integer, intent(in) :: MPI_rank

    ! ... Local variables ...
        integer :: n, values(1:8)
        integer, allocatable :: seed(:)

        call date_and_time(values=values)
        call random_seed(size=n)
        allocate(seed(n))
        seed(1:n) = values(8) + MPI_rank
        call random_seed(put=seed)

    end subroutine init_RNG

end module util


module square_lattice_FT
    use types 
    implicit none 
    private 
    public zexpiqr
    public lq, nk, init
    public listk 

    public calc_FT_coeff 
    public FT_GreensFunction
    public order_Ksites

    real(dp), parameter :: pi = dacos(-1.d0)
    integer :: l
    integer :: a1_p(2), a2_p(2)
    real(dp) :: b1_p(2), b2_p(2)
    integer, allocatable :: listk(:,:), list(:,:)
    integer, allocatable :: invlistk(:,:)

    integer :: lq, nk
    complex(dp), allocatable :: zexpiqr(:,:)
    ! Are the FT coeff's initialized ? 
    logical :: init = .false.           

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
        print*, "nk=", nk
        a1_p(1) = 1 ; a1_p(2) = 0
        a2_p(1) = 0 ; a2_p(2) = 1
        b1_p(1) = 2.d0*pi/dble(l) ; b1_p(2) = 0.d0
        b2_p(1) = 0.d0            ; b2_p(2) = 2.d0*pi/dble(l)

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
                zexpiqr(i,iq) = exp( dcmplx( 0.d0, qvec(1)*ri(1) + qvec(2)*ri(2) )  ) / sqrt(dble(lq))
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

    subroutine order_Ksites(Ksites, ordering)
    ! Purpose:
    ! --------
    ! Arrange momentum points in a certain order. 
    !
    ! odering == 'pair':    
    !   Order the momentum points in such a ways that momenta related 
    !   by momentum inversion come in pairs, one after the other. 
    ! ordering == 'quad':
    !   Order the momentum points such that the upper right quadrant 
    !   (kx >=0, ky>=0) comes first. 
    ! ordering == 'sqFS': 
    !   Order the momentum points such that the points on the diamond-shaped
    !   Fermi surface of the square lattice at half filling come first. 
    !
    ! Precondition:
    ! -------------
    ! The system is square, i.e. number of sites is n = l*l
    !
    ! Arguments:
    ! ----------
        use types
        integer, intent(out) :: Ksites(:)
        character(len=4) :: ordering 
    ! ...Local variables... 
        integer :: i, j, l, n
        integer :: k, k_inv
        integer :: counter 
        logical :: taken(size(Ksites,1))
        integer :: Ksites_reordered(size(Ksites,1))
        n = size(Ksites, 1)
        l = nint( sqrt(dble(n)) )
        
        taken(:) = .false.
        Ksites_reordered(:) = -1
        
        select case(ordering)

            case('pair')
                counter = 1
                do k=1,n
                    if (.not.taken(k)) then 
                        Ksites_reordered(counter) = k
                        taken(k) = .true.
                        counter = counter + 1 
                    endif 
                    i = listk(k,1)
                    j = listk(k,2)
                    if (i < l/2) then 
                        i = -i
                    endif 
                    if (j < l/2) then 
                        j = -j
                    endif 
                    k_inv = invlistk(i,j)
                    if( .not.(taken(k_inv)) .and. (k_inv /= k) ) then 
                        Ksites_reordered(counter) = k_inv
                        taken(k_inv) = .true.
                        counter = counter + 1
                    endif 
                enddo 
        
                if (counter /= (n+1)) then 
                    print*, "order_Ksites(): ERROR: counter /= n+1"
                    stop
                endif         

            case('quad')
                counter = 1
                do j=0,l/2
                    do i=0,l/2
                        k = invlistk(i,j)
                        taken(k) = .true.
                        Ksites_reordered(counter) = k
                        counter = counter + 1
                    enddo 
                enddo
                ! the remaining BZ ...
                do k=1,n
                    if( .not.taken(k) ) then 
                        taken(k) = .true.                        
                        Ksites_reordered(counter) = k
                        counter = counter + 1
                    endif 
                enddo 

                if (counter /= (n+1)) then 
                    print*, "order_Ksites(): ERROR: counter /= n+1"
                    stop
                endif     

            case('sqFS')
                counter = 1
                ! lower right edge
                do i = 1, l/2, +1
                    j = -l/2 + i
                    k = invlistk(i, j)
                    taken(k) = .true.
                    Ksites_reordered(counter) = k
                    counter = counter + 1 
                enddo 
                ! upper right edge 
                do i = 0, l/2-1, +1
                    j =  l/2 - i
                    k = invlistk(i, j)
                    taken(k) = .true.
                    Ksites_reordered(counter) = k
                    counter = counter + 1 
                enddo 
                ! lower left edge 
                do i = -1, -l/2+1, -1
                    j = -l/2 - i
                    k = invlistk(i, j)
                    taken(k) = .true.
                    Ksites_reordered(counter) = k
                    counter = counter + 1                     
                enddo 
                ! upper left edge 
                do i = -1, -l/2+1, -1
                    j = l/2 + i
                    k = invlistk(i, j)
                    taken(k) = .true.
                    Ksites_reordered(counter) = k
                    counter = counter + 1 
                enddo 
                ! the remaining BZ ...
                do k=1,n
                    if( .not.taken(k) ) then 
                        taken(k) = .true.                        
                        Ksites_reordered(counter) = k
                        counter = counter + 1
                    endif 
                enddo 

            case('diag')
                counter = 1
                do i = 1-l/2, l/2, +1
                    j = i
                    k = invlistk(i, j)
                    taken(k) = .true.
                    Ksites_reordered(counter) = k 
                    counter = counter + 1
                enddo
                ! the remaining BZ ...
                do k=1,n
                    if( .not.taken(k) ) then 
                        taken(k) = .true.                        
                        Ksites_reordered(counter) = k
                        counter = counter + 1
                    endif 
                enddo 

            case default
                do k=1,n 
                    Ksites_reordered(k) = k
                enddo 

        end select 

        Ksites(:) = Ksites_reordered(:)

    end subroutine order_Ksites

end module square_lattice_FT


module sample_pseudo_DM
    implicit none 
    private 
    public sample_FF_GreensFunction 

    contains 
    subroutine sample_FF_GreensFunction(G, occ_vector, abs_corr, Ksites, reweighting_phase, reweighting_factor)
        ! Purpose:
        ! --------
        ! Sample pseudo-snapshots from a free-fermion pseudo-density matrix using 
        ! the single-particle Green's function for a given spin species as input. 
        ! In general, the single-particle Green's function may be complex as is the case 
        ! in the momentum basis. In that case there is a phase problem and 
        ! every snapshot is associated a complex phase `exp(i phi) = a + 1j b` 
        ! and a reweighting factor. 
        use types 
        use util
        use square_lattice_FT, only: order_Ksites
        implicit none 
        ! Arguments:
        ! ----------
        complex(dp), intent(in) :: G(:,:)                  ! Green's function matrix.
        integer, intent(out) :: occ_vector(:)              ! vector of occupation numbers for given spin species.
        real(dp), intent(out) :: abs_corr(:)               
        complex(dp), intent(out) :: reweighting_phase 
        real(dp), intent(out) :: reweighting_factor        
        ! REMOVE
        integer, intent(out) :: Ksites(:)                  ! Factor ordering of the components. 
        ! REMOVE

        ! ... Local Variables ...
        integer :: D  ! dimension of the single-particle Hilbert space 
        complex(dp), allocatable :: corr(:)  ! correction term taking care of correlations between different sites 
        complex(dp), allocatable :: cond_prob(:)   ! conditional probability 
        real(dp), allocatable :: q_prob(:)      ! absolute value of the conditional probability 
        real(dp), allocatable :: phase_angle(:) ! phase angle of the conditional probability
        complex(dp), allocatable :: phase(:)    ! exp(i phi(:))
        real(dp) :: norm
        real(dp) :: eta 
        integer :: k, occ
        ! CHANGE BACK
        ! integer, allocatable :: Ksites(:)
        ! CHANGE BACK        
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
        ! CHANGE BACK
        ! allocate(Ksites(D))
        ! CHANGE BACK
        do ii = 1, D 
            Ksites(ii) = ii 
        enddo 
        !! call random_permutation(Ksites)
        call order_Ksites(Ksites, 'sqFS')
        !! REMOVE
        D = 2*int(sqrt(float(D))) - 2
        ! D = (2*int(sqrt(float(D))) - 2) / 2 + 1
        !!! REMOVE

        ! helper variables 
        allocate(Xinv(1:D,1:D))
        allocate(Xinv_new(1:D,1:D))
        allocate(uu(1:D))
        allocate(vv(1:D))

        corr = complex(ZERO, ZERO)
        cond_prob = complex(ZERO, ZERO)
        reweighting_factor = 1.0_dp
        reweighting_phase = complex(1.0_dp, 0.0_dp)

        occ_vector(:) = -1 ! initialize to invalid value 
        Xinv = complex(ZERO, ZERO)
        Xinv_new = complex(ZERO, ZERO)

        ! Componentwise direct sampling
        do k = 1, D
            ! "Correction" due to correlations between sites
            if ( k == 1 ) then 
                corr(1) = complex(ZERO, ZERO)
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
                phase(ii) = complex(cos(phase_angle(ii)), sin(phase_angle(ii)))
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

            print*, k, "reweighting_factor=", reweighting_factor 

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

end module 


program sample_kspace
    use types 
    use parallelization
    use util
    use square_lattice_FT
    use sample_pseudo_DM
    implicit none 
    
    integer :: Nsites
    integer :: Nx
    integer :: Ny

    integer :: i,j, ss, sss, spin_idx
    real(dp), allocatable :: Green_xspace(:,:)
    complex(dp), allocatable :: Green_kspace(:,:)
    integer, allocatable :: occ_vector(:)
    real(dp), allocatable :: occ_matrix(:,:)
    complex(dp) :: weight_phase
    real(dp) :: weight_factor
    real(dp), allocatable :: abs_corr(:)
    integer, allocatable :: Ksites(:), invKsites(:)
    
    character(len=3)  :: chr_spin(1:2)
    character(len=30) :: filename 
    integer :: max_HS_samples
    integer :: Nsamples_per_HS
    integer :: skip
    integer, parameter :: Nspin_species = 2

    integer, allocatable :: occ_vector_tmp(:,:)
    complex(dp), allocatable :: weight_phase_tmp(:)
    real(dp), allocatable :: weight_factor_tmp(:)

    namelist /simpara/ filename, Nsites, max_HS_samples, Nsamples_per_HS, skip 

#if defined (USE_MPI)
    include "mpif.h"
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, MPI_rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, MPI_size, ierr)
#else
    MPI_rank = 0
    MPI_size = 1
#endif 
    write(chr_rank, "(i5.5)") MPI_rank
    chr_spin(1) = "_up"
    chr_spin(2) = "_dn"

    ! Check that MPI is working correctly 
    print*, "I am rank", MPI_rank, "of ", MPI_size

    print*, "Reading parameter file..."
    open(100, file="simparams.in", status="old", action="read")
    read(100, nml=simpara)
    close(100)

    Nx = int(sqrt(float(Nsites)))
    Ny = Nx

    allocate( Green_xspace(Nsites, Nsites) )
    allocate( Green_kspace(Nsites, Nsites) )
    allocate( occ_vector(Nsites) )
    allocate( abs_corr(Nsites) )
    allocate( Ksites(Nsites) )

    allocate( occ_vector_tmp(Nsites, Nsamples_per_HS) )
    allocate( weight_factor_tmp(Nsamples_per_HS) )
    allocate( weight_phase_tmp(Nsamples_per_HS) )

    ! allocate(invKsites(Nsites))
    ! allocate(occ_matrix(1-Nx/2:Nx/2, 1-Ny/2:Ny/2))

    call init_RNG(MPI_rank)

    do spin_idx = 1, 2
        open(50+spin_idx+MPI_rank, file="Green_ncpu"//chr_rank//chr_spin(spin_idx)//".dat", status='old', action='read')
    enddo 
    do ss = 1, max_HS_samples
        do spin_idx = 1, 2
            occ_vector_tmp(:,:) = ZERO
            weight_factor_tmp(:) = ZERO
            weight_phase_tmp(:) = cmplx(ZERO, ZERO)
            if (ss > skip) then 

                print*, "MPI_rank=", MPI_rank, "of MPI_size", MPI_size, "spin_idx=", spin_idx, "HS_sample=", ss
                Green_xspace = ZERO
                read(50+spin_idx+MPI_rank, *) Green_xspace
                read(50+spin_idx+MPI_rank, *)
                read(50+spin_idx+MPI_rank, *)
                print*, Green_xspace(1,1)

                ! uncomment, if you want to sample in momentum space
                ! call FT_GreensFunction(Green_xspace, Green_kspace)
                Green_kspace = Green_xspace

                do sss = 1, Nsamples_per_HS
                    call sample_FF_GreensFunction(Green_kspace, occ_vector, abs_corr, Ksites, weight_phase, weight_factor)
                    occ_vector_tmp(:, sss) = occ_vector(:)
                    weight_factor_tmp(sss) = weight_factor
                    weight_phase_tmp(sss)  = weight_phase
                enddo 

                open(100, file="Fock_samples_ncpu"//chr_rank//chr_spin(spin_idx)//".dat", status="unknown", position="append")
                do sss = 1, Nsamples_per_HS
                    write(100, *) 1.0, 1.0, real(weight_phase_tmp(sss)), aimag(weight_phase_tmp(sss)), &
                                  weight_factor_tmp(sss), occ_vector_tmp(:, sss)
                enddo 
                close(100)

            endif 
        enddo 
    enddo

    do spin_idx = 1, 2
        close(50+spin_idx+MPI_rank)
    enddo

#if defined(USE_MPI)
    call MPI_FINALIZE(ierr)
#endif 

end program sample_kspace    
