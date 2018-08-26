C
C=======================================================================
C
C     String handling routines:
C     Strings are handled differently in C++ and in FORTRAN.  In order
C     to simplify the framework strings are passed from FORTRAN to C++
C     as arrays of integer characters, terminated by the EOS symbol
C     which we set to -1
C     blstr2int converts a FORTRAN string to an integer array,
C     blint2str converts an integer array to a FORTRAN string.
C      
C-----------------------------------------------------------------------
      SUBROUTINE blstr2int(iarr, n, str)
      CHARACTER*(*) str
      INTEGER n, i, j
      INTEGER iarr(n)
      INTEGER EOS
      PARAMETER (EOS=-1)
C
      IF ( n .LE. len(str) ) THEN
          call bl_abort("blstr2int: str to large for iarr")
      END IF
C     Make sure that IARR is empty
      DO J = 1, N
          iarr(J) = ichar(' ')
      END DO
      j = 1
      DO i = 1, len(str)
          iarr(j) = ichar(str(i:i))
          j = j + 1
      END DO
C     EOS
      iarr(j) = EOS
C
      END
C-----------------------------------------------------------------------
      SUBROUTINE blint2str(str, iarr, n)
      CHARACTER*(*) str
      INTEGER n
      INTEGER iarr(n)
      INTEGER EOS
      PARAMETER (EOS=-1)
      INTEGER i
C
      DO i = 1, LEN(str)
          str(i:i) = ' '
      END DO
      DO i = 1, n
          IF ( i .GT. LEN(str) ) then
             call bl_abort("blint2str: iarr to large for str")
          end if
          IF ( iarr(i) .EQ. EOS ) GO TO 100
          str(i:i) = char(iarr(i))
      END DO
 100  CONTINUE
C
      END
