!Comment
!
! INTEGER FUNCTION SYSTEMFUNCTION(COMANDO)
!
! Input: COMANDO
! Output: SYSTEMFUNCTION (function)
!
! Execute a shell command and return its exit status (0 when the command
! succeeds), or -1 if the command could not be executed.
!
! CHARACTER*(*) COMANDO -> command to be executed
!
! Previous versions used the GNU Fortran extension SYSTEM, which returned
! the status in a system-dependent encoding. The value is not used by the
! callers.
!
!Comment
!------------------------------------------------------------------------------
        INTEGER FUNCTION SYSTEMFUNCTION(COMANDO)
        IMPLICIT NONE
        CHARACTER*(*) COMANDO
!
        INTEGER IEXIT,ICMD
!------------------------------------------------------------------------------
        IEXIT=0
! CMDSTAT is present so that a command that cannot be executed does not
! stop the program (error termination otherwise)
        CALL EXECUTE_COMMAND_LINE(COMANDO,EXITSTAT=IEXIT,CMDSTAT=ICMD)
        IF(ICMD.NE.0)THEN
          SYSTEMFUNCTION=-1
        ELSE
          SYSTEMFUNCTION=IEXIT
        END IF
        END
