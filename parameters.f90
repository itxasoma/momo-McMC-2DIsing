MODULE GLOBAL
  IMPLICIT NONE
  double precision :: T, beta
  integer :: L, N, z, nMCS, nmeas

CONTAINS

  SUBROUTINE init_globals()
    IMPLICIT NONE
    beta = 1.0d0 /T
    N = L * L
  END SUBROUTINE init_globals

END MODULE GLOBAL


 