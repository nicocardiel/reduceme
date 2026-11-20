C------------------------------------------------------------------------------
Comment
C
C REAL FUNCTION RANRED(NSEED)
C
C Input: NSEED
C Output: RANRED (function)
C
C Return a random number in the range [0,1).
C
C INTEGER NSEED -> NSEED=0: RANRED returns the next random number in the
C                           sequence.
C                  NSEED<0: the generator is initialised with a seed taken
C                           from the system clock, and RANRED returns a
C                           random number in the new sequence.
C                           In this case NSEED returns 0.
C                  NSEED>0: the generator is initialised with the seed
C                           NSEED, and RANRED returns a random number in the
C                           new sequence. In this case NSEED returns 0.
C
C The generator reproduces exactly the GNU Fortran extensions RAND() and
C SRAND() used in previous versions, written in standard Fortran:
C - Park & Miller "minimal standard" generator:
C   X(n+1) = 16807 * X(n) mod (2**31-1)
C   (seed 0 is replaced by 123459876, as in SRAND);
C - returned value: REAL((X(n+1)-1)/2**9) / 2**22, i.e. the integer is
C   truncated to 22 bits, which is exactly representable in single
C   precision and always lower than 1.0.
C The only difference is the seed used when NSEED<0: SYSTEM_CLOCK instead of
C the GNU extension TIME() (both give non-reproducible sequences).
C
Comment
        REAL FUNCTION RANRED(NSEED)
        USE ISO_FORTRAN_ENV, ONLY: INT64
        IMPLICIT NONE
        INTEGER NSEED
C
        INTEGER(KIND=INT64) IA,IM
        PARAMETER (IA=16807_INT64,IM=2147483647_INT64)
C state of the generator; RAND() starts with 1 when SRAND() is not called
        INTEGER(KIND=INT64) ISTATE
        SAVE ISTATE
        DATA ISTATE/1_INT64/
        INTEGER ICOUNT
C------------------------------------------------------------------------------
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
C 16807 * ISTATE < 2**46: no overflow in 64-bit integers
        ISTATE=MOD(IA*ISTATE,IM)
        RANRED=REAL((ISTATE-1_INT64)/512_INT64)/4194304.0
C
        END
