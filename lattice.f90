MODULE lattice
    USE GLOBAL
    IMPLICIT NONE
    integer, allocatable:: s(:)  ! spins: +1 o -1
    integer, allocatable:: nbr(:,:)   ! nbr(i,k): k-essim neighb of site i
    integer, allocatable:: ina(:,:)   ! pbc

CONTAINS

SUBROUTINE alloc_lattice(N)
    integer, intent(in):: N
    allocate(s(N))
    allocate(nbr(N, z))
    allocate(ina(0:1, L))
END SUBROUTINE

SUBROUTINE dealloc_lattice()
    deallocate(s, nbr, ina)
END SUBROUTINE

! To account for periodic (toroidal) b.c. conditions we can define shift array
SUBROUTINE pbc()
    IMPLICIT NONE
    integer:: i

    do i = 1, L
        ina(0,i) = i-1
        ina(1,i) = i+1
    enddo
    ina(0,1) = L
    ina(1,L) = 1

END SUBROUTINE


! Build square lattice
SUBROUTINE square_lattice()
    IMPLICIT NONE
    integer:: i
    double precision:: x,y
    i = 0
    do y = 1, L
        do x = 1, L
            i= i+1
            nbr(i,1) = ina(1,x) + L*(y-1)
            nbr(i,2) = ina(0,x) + L*(y-1)
            nbr(i,3) = x + L*(ina(1,y)-1)
            nbr(i,4) = x + L*(ina(0,y)-1)
        enddo
    enddo
END SUBROUTINE


END MODULE
