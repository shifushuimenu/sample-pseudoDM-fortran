
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

    subroutine init_RNG
    ! Purpose:
    ! --------
    !    Initialize the standard pseudo-random number generator 
    !    with a seed obtained from the system time (at the millisecond level).
    !    Call the random number generator:
    !
    !        integer :: eta
    !        call random_number(eta)
    !
    ! ... Local variables ...
        integer :: n, values(1:8)
        integer, allocatable :: seed(:)

        call date_and_time(values=values)
        call random_seed(size=n)
        allocate(seed(n))
        seed(1:n) = values(8) 
        call random_seed(put=seed)

    end subroutine init_RNG

end module util


module sample_pseudo_DM
    implicit none 
    private 
    public sample_FF_GreensFunction 

    contains 
    subroutine sample_FF_GreensFunction(G, occ_vector, abs_corr, Ksites, reweighting_phase, reweighting_factor)
        use types 
        use util
        implicit none 
        ! Arguments:
        ! ----------
        complex(dp), intent(in) :: G(:,:)
        integer, intent(out) :: occ_vector(:)
        real(dp), intent(out) :: abs_corr(:)
        complex(dp), intent(out) :: reweighting_phase 
        real(dp), intent(out) :: reweighting_factor        
        ! REMOVE
        integer, intent(out) :: Ksites(:)
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
        ! call random_permutation(Ksites)

        ! helper variables 
        allocate(Xinv(1:D,1:D))
        allocate(Xinv_new(1:D,1:D))
        allocate(uu(1:D))
        allocate(vv(1:D))

        corr = complex(ZERO, ZERO)
        cond_prob = complex(ZERO, ZERO)
        reweighting_factor = 1.0_dp
        reweighting_phase = complex(1.0_dp, 0.0_dp)

        occ_vector = -1 ! initialize to invalid value 
        Xinv = complex(ZERO, ZERO)
        Xinv_new = complex(ZERO, ZERO)

        ! Componentwise direct sampling
        do k = 1, D
            ! "Correction" due to correlations between sites
            if ( k == 1 ) then 
                corr(1) = complex(ZERO, ZERO)
            elseif ( k == 2) then
                ! matmul() does not work for scalars 
                corr(k) = G(k, Ksites(1)) * Xinv(1,1) * G(Ksites(1), k)
            else
                corr(k) = dot_product( G(k, Ksites(1:k-1)), matmul( Xinv(1:k-1, 1:k-1), G(Ksites(1:k-1), k) ) )
            endif 
        
            cond_prob(1) = 1.0_dp - G(k,k) + corr(k)
            cond_prob(0) = G(k,k) - corr(k)


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

            occ_vector(k) = occ 

            if( k==1 ) then 
                Xinv(1,1) = 1.0_dp / (G(1,1) - occ)
            else
                ! Low-rank update 
                gg = 1.0_dp / (G(k,k) - occ - corr(k))
                if ( k==2 ) then 
                    ! matmul() does not work for scalars 
                    uu(1) = Xinv(1,1) * G(Ksites(1), 2)
                    vv(1) = G(2, Ksites(1)) * Xinv(1,1)
                else
                    uu(1:k-1) = matmul( Xinv(1:k-1,1:k-1), G(Ksites(1:k-1),k) )
                    vv(1:k-1) = matmul( G(k, Ksites(1:k-1)), Xinv(1:k-1,1:k-1) )
                endif 
                do jj = 1, k-1
                    do ii = 1, k-1
                        Xinv_new(Ksites(ii), Ksites(jj)) = Xinv(Ksites(ii), Ksites(jj)) &
                            + gg * uu(ii) * vv(jj)
                    enddo 
                enddo 
                do ii = 1, k-1
                    Xinv_new(k, Ksites(ii)) = -gg*vv(Ksites(ii))
                    Xinv_new(Ksites(ii), k) = -gg*uu(Ksites(ii))
                enddo
                Xinv_new(k,k) = gg
                Xinv(1:k, 1:k) = Xinv_new(1:k, 1:k)

            endif 
        enddo 

        abs_corr(1:D) = abs(corr(1:D))

    end subroutine 

end module 

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

    subroutine order_Ksites(Ksites)
    ! Purpose:
    ! --------
    ! Order the momentum points in such a ways that momenta related 
    ! by momentum inversion come in pairs, one after the other. 
    !
    ! Precondition:
    ! The system is square, i.e. number oif sites is n = l*l
    !
    ! Arguments:
    ! ----------
        use types
        integer, intent(out) :: Ksites(:)
        integer :: i, j, l, n
        integer :: k, k_inv
        integer :: counter 
        logical :: taken(size(Ksites,1))
        integer :: Ksites_reordered(size(Ksites,1))
        n = size(Ksites, 1)
        l = nint( sqrt(dble(n)) )
        
        taken(:) = .false.
        Ksites_reordered(:) = -1
        
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

        Ksites(:) = Ksites_reordered(:)

    end subroutine order_Ksites


end module square_lattice_FT


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
    open(100, file="simpara.in", status="old", action="read")
    read(100, nml=simpara)
    close(100)

    Nx = int(sqrt(float(Nsites)))
    Ny = Nx

    allocate( Green_xspace(Nsites, Nsites) )
    allocate( Green_kspace(Nsites, Nsites) )
    allocate(occ_vector(Nsites))
    allocate(abs_corr(Nsites))
    allocate(Ksites(Nsites))
    ! allocate(invKsites(Nsites))
    ! allocate(occ_matrix(1-Nx/2:Nx/2, 1-Ny/2:Ny/2))

    call init_RNG

    do spin_idx = 1, 2
        open(50+spin_idx+MPI_rank, file="Green_ncpu"//chr_rank//chr_spin(spin_idx)//".dat", status='old', action='read')
    enddo 
    do ss = 1, max_HS_samples
        do spin_idx = 1, 2
            if (ss > skip) then 

                print*, "MPI_rank=", MPI_rank, "of MPI_size", MPI_size, "spin_idx=", spin_idx, "HS_sample=", ss
                Green_xspace = ZERO
                read(50+spin_idx+MPI_rank, *) Green_xspace
                read(50+spin_idx+MPI_rank, *)
                read(50+spin_idx+MPI_rank, *)
                print*, Green_xspace(1,1)

                call FT_GreensFunction(Green_xspace, Green_kspace)

                ! Green_kspace = Green_xspace
                do sss = 1, Nsamples_per_HS
                    call sample_FF_GreensFunction(Green_kspace, occ_vector, abs_corr, Ksites, weight_phase, weight_factor)
                enddo 

                open(100, file="KFock_samples_ncpu"//chr_rank//chr_spin(spin_idx)//".dat", status="unknown", position="append")
                write(100, *) 1.0, 1.0, real(weight_phase), aimag(weight_phase), weight_factor, occ_vector(:)
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