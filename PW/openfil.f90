!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE openfil()
  !----------------------------------------------------------------------------
  !
  ! ... This routine opens all files needed to the self consistent run,
  ! ... sets various file names, units, record lengths
  !
  USE parameters,     ONLY :  DP
  USE io_global,      ONLY :  stdout
  USE basis,          ONLY :  natomwfc, startingwfc
  USE wvfct,          ONLY :  nbnd, npwx
  USE varie,          ONLY :  order, lneb
  USE ldaU,           ONLY :  lda_plus_U
  USE io_files,       ONLY :  prefix, &
                              iunat, iunocc, iunwfc, iunoldwfc, iunoldwfc2, &
                              iunigk, nwordwfc, nwordatwfc, iunneb
  USE restart_module, ONLY :  readfile_new
!#ifdef __PARA
  USE para,           ONLY : 
!#endif
  !
  IMPLICIT NONE
  !
  LOGICAL       :: exst
  INTEGER       :: ndr, kunittmp, ierr
  REAL(KIND=DP) :: edum(1,1), wdum(1,1)
  !
  !
  ! ... iunwfc contains the wavefunctions
  !
  iunwfc = 10
  !
  ! ... iunoldwfc contains the old wavefunctions, used in molecular dynamics
  !
  iunoldwfc = 11
  iunoldwfc2= 12
  !
  ! ... nwordwfc is the record length for the direct-access file
  ! ... containing wavefunctions
  !
  nwordwfc = 2 * nbnd * npwx
  !
  CALL diropn( iunwfc, TRIM( prefix )//'.wfc', nwordwfc, exst )
  !
  IF ( startingwfc == 'file' .AND. .NOT. exst ) THEN
     !
#if defined __NEW_PUNCH
     ndr      = 4
     kunittmp = 1
#  ifdef __PARA
     kunittmp = kunit
#  endif
     !
     CALL readfile_new( 'wave', ndr, edum, wdum, kunittmp, nwordwfc, &
                        iunwfc, ierr )
     IF ( ierr > 0 ) THEN
        !
#else
        WRITE( stdout, '(5X,"Cannot read wfc file: not found")' )
        startingwfc = 'atomic'
#endif
#if defined __NEW_PUNCH
     END IF
#endif
     !
  END IF
  !
  ! ... Needed for LDA+U
  !
  ! ... iunat contains the orthogonalized wfcs
  !
  iunat = 13
  nwordatwfc = 2 * npwx * natomwfc
  !
  IF ( lda_plus_u ) &
     CALL diropn( iunat, TRIM( prefix )//'.atwfc', nwordatwfc, exst )
  !
  ! ... iunocc contains the atomic occupations computed in new_ns
  ! ... it is opened and closed for each reading-writing operation
  !
  iunocc = 14
  !
  ! ... if extrapolation of wfc's is requested (order=2)
  ! ... another file is needed to store the "old" wfc's
  !
  IF ( order > 1 ) &
     CALL diropn( iunoldwfc, TRIM( prefix )//'.oldwfc', nwordwfc, exst )
  !
  IF ( order > 2 ) &
     CALL diropn( iunoldwfc2, TRIM( prefix )//'.oldwfc2', nwordwfc, exst )
  !
  ! ... iunigk contains the number of PW and the indices igk
  ! ... Note that unit 15 is reserved for error messages 
  !
  iunigk = 16
  CALL seqopn( iunigk, TRIM( prefix )//'.igk', 'UNFORMATTED', exst )
  !
  RETURN
  !
END SUBROUTINE openfil

