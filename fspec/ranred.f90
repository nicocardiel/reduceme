!------------------------------------------------------------------------------
!Comment
!
! REAL FUNCTION RANRED(NSEED)
!
! Input: NSEED
! Output: RANRED (function)
!
! Return a random number in the range [0,1).
!
! INTEGER NSEED -> NSEED=0: RANRED returns the next random number in the
!                           sequence.
!                  NSEED<0: the generator is initialised with a seed taken
!                           from the system clock, and RANRED returns a
!                           random number in the new sequence.
!                           In this case NSEED returns 0.
!                  NSEED>0: the generator is initialised with the seed
!                           NSEED, and RANRED returns a random number in the
!                           new sequence. In this case NSEED returns 0.
!
! The generator reproduces exactly the GNU Fortran extensions RAND() and
! SRAND() used in previous versions, written in standard Fortran:
! - Park & Miller "minimal standard" generator:
!   X(n+1) = 16807 * X(n) mod (2**31-1)
!   (seed 0 is replaced by 123459876, as in SRAND);
! - returned value: REAL((X(n+1)-1)/2**9) / 2**22, i.e. the integer is
!   truncated to 22 bits, which is exactly representable in single
!   precision and always lower than 1.0.
! The only difference is the seed used when NSEED<0: SYSTEM_CLOCK instead of
! the GNU extension TIME() (both give non-reproducible sequences).
!
!Comment
        REAL FUNCTION RANRED(NSEED)
        USE ISO_FORTRAN_ENV, ONLY: INT64
        IMPLICIT NONE
        INTEGER NSEED
!
        INTEGER(KIND=INT64) IA,IM
        PARAMETER (IA=16807_INT64,IM=2147483647_INT64)
! state of the generator; RAND() starts with 1 when SRAND() is not called
        INTEGER(KIND=INT64) ISTATE
        SAVE ISTATE
        DATA ISTATE/1_INT64/
        INTEGER ICOUNT
!------------------------------------------------------------------------------
        IF(NSEED.LT.0)THEN
          CALL SYSTEM_CLOCK(ICOUNT)
          ISTATE=INT(ICOUNT,KIND=INT64)
          IF(ISTATE.LE.0_INT64) ISTATE=0_INT64     !no clock: seed as for 0
        ELSEIF(NSEED.GT.0)THEN
          ISTATE=INT(NSEED,KIND=INT64)
        END IF
        IF(NSEED.NE.0)THEN
          IF(ISTATE.EQ.0_INT64) ISTATE=123459876_INT64
          NSEED=0
        END IF
! 16807 * ISTATE < 2**46: no overflow in 64-bit integers
        ISTATE=MOD(IA*ISTATE,IM)
        RANRED=REAL((ISTATE-1_INT64)/512_INT64)/4194304.0
!
        END
