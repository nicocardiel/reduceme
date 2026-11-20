!------------------------------------------------------------------------------
! Version 20-November-2026                                   file: giveerrf.f90
! @ ncl & fjg
!------------------------------------------------------------------------------
!Comment
!
! Program: giveerrf.f
! Classification: miscellany
! Description: Given a character string corresponding to a file name, this
! program returns the expected error frame associated to that file.
!
!Comment
!
        PROGRAM GIVEERRF
        IMPLICIT NONE
!
!       INCLUDE 'redlib.inc'
!
        CHARACTER*75 INFILE,ERRFILE
!------------------------------------------------------------------------------
        IF(COMMAND_ARGUMENT_COUNT().LT.1)THEN
          WRITE(*,101) 'ERROR: you must provide a file name'//&
     &     ' in the command line'
        ELSE
          CALL GET_COMMAND_ARGUMENT(1,INFILE)
          CALL GUESSEF(INFILE,ERRFILE)
          WRITE(*,101) ERRFILE
        END IF
!------------------------------------------------------------------------------
        STOP
101     FORMAT(A) 
        END
