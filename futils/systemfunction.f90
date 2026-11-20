Comment
C
C INTEGER FUNCTION SYSTEMFUNCTION(COMANDO)
C
C Input: COMANDO
C Output: SYSTEMFUNCTION (function)
C
C Execute a shell command and return its exit status (0 when the command
C succeeds), or -1 if the command could not be executed.
C
C CHARACTER*(*) COMANDO -> command to be executed
C
C Previous versions used the GNU Fortran extension SYSTEM, which returned
C the status in a system-dependent encoding. The value is not used by the
C callers.
C
Comment
C------------------------------------------------------------------------------
        INTEGER FUNCTION SYSTEMFUNCTION(COMANDO)
        IMPLICIT NONE
        CHARACTER*(*) COMANDO
C
        INTEGER IEXIT,ICMD
C------------------------------------------------------------------------------
        IEXIT=0
C CMDSTAT is present so that a command that cannot be executed does not
C stop the program (error termination otherwise)
        CALL EXECUTE_COMMAND_LINE(COMANDO,EXITSTAT=IEXIT,CMDSTAT=ICMD)
        IF(ICMD.NE.0)THEN
          SYSTEMFUNCTION=-1
        ELSE
          SYSTEMFUNCTION=IEXIT
        END IF
        END
