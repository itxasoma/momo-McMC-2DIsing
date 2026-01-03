MODULE ising_energy
    use GLOBAL
    use lattice
    IMPLICIT NONE
CONTAINS

! Computes total energy E for the Ising model
SUBROUTINE E_ising(s, nbr, E)
    IMPLICIT NONE
    integer, intent(in):: s(:)       ! spins: +1 or -1
    integer, intent(in):: nbr(:,:)   ! nbr(i,k): k-th neighbour of site i
    integer, intent(out):: E
    integer:: i, k
    double precision:: sums

    sums = 0.d0
    do i = 1, N
        do k = 1, z
            sums = sums + s(i) * s(nbr(i,k))
        enddo
    enddo

    E = -0.5d0 * sums

END SUBROUTINE

! OLD VERSION (DO NOT REMOVE)
! SUBROUTINE h_local(z,s, nbr, hh)
!     IMPLICIT NONE
!     integer, intent(in):: s(:),z       ! spins: +1 or -1
!     integer, intent(in):: nbr(:,:)   ! nbr(i,k): k-th neighbour of site i
!     integer, intent(out):: hh(:)
!     integer:: i, k
!     print*, "z:", z
!     hh(:) = 0
!     do i = 1, N
!         do k = 1, z
!             hh(i) = hh(i) + s(nbr(i,k))
!         enddo
!     enddo

! END SUBROUTINE

SUBROUTINE h_local(z, s, nbr, hh, ispin)
    IMPLICIT NONE
    integer, intent(in):: s(:), z       ! spins: +1 or -1
    integer, intent(in):: nbr(:,:)   ! nbr(i,k): k-th neighbour of site i
    integer, intent(out):: hh
    integer:: ispin, k
    !print*, "z:", z
    hh = 0
    do k = 1, z
        hh = hh + s(nbr(ispin,k))
    enddo

END SUBROUTINE

SUBROUTINE table_exp(z, table)
! TODO : REALMENTE SOLO NECESITAS TABULAR VALORES > 0 DE DELTA E
!  CREO QUE TAMBIEN PUEDES DEFINIR TABLE(-8:8:2) (SI COGES NEGATIVOS) O TABLE(0:8:2)
! Creates table of exponentials to be read in the main program
    IMPLICIT NONE
    integer, intent(in) :: z
    integer::n
    double precision, intent(out) :: table(2*z+1)
    do n=0,size(table)-1
        table(n+1)=exp(-beta*(-2*z+2*n))  
        ! print*, n+1, -2*z+2*n, exp(-beta*(-2*z+2*n))   
    enddo  
END SUBROUTINE


! Computes total magnetization M for the Ising model
SUBROUTINE M_ising(s, M)
    IMPLICIT NONE
    integer, intent(in) :: s(:)
    integer, intent(out):: M
    integer:: i
    double precision:: sums

    sums = 0.d0

    do i = 1, N
        sums = sums + dble( s(i) )
    enddo

    M = sums

END SUBROUTINE



END MODULE
