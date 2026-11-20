C------------------------------------------------------------------------------
C Version 07-April-2003                                        file: giveerrf.f
C @ ncl & fjg
C------------------------------------------------------------------------------
Comment
C
C Program: giveerrf.f
C Classification: miscellany
C Description: Given a character string corresponding to a file name, this
C program returns the expected error frame associated to that file.
C
Comment
C
        PROGRAM GIVEERRF
        IMPLICIT NONE
C
!       INCLUDE 'redlib.inc'
C
        CHARACTER*75 INFILE,ERRFILE
C------------------------------------------------------------------------------
        IF(COMMAND_ARGUMENT_COUNT().LT.1)THEN
          WRITE(*,101) 'ERROR: you must provide a file name'//
     +     ' in the command line'
        ELSE
          CALL GET_COMMAND_ARGUMENT(1,INFILE)
          CALL GUESSEF(INFILE,ERRFILE)
          WRITE(*,101) ERRFILE
        END IF
C------------------------------------------------------------------------------
        STOP
101     FORMAT(A) 
        END
