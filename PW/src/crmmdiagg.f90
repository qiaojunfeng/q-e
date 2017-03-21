!
! Copyright (C) 2017 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO ( 0._DP, 0._DP )
#define ONE  ( 1._DP, 0._DP )
!
!----------------------------------------------------------------------------
SUBROUTINE crmmdiagg( npwx, npw, nbnd, npol, psi, spsi, e, &
                      g2kin, btype, ethr, ndiis, uspp, notconv, rmm_iter )
  !----------------------------------------------------------------------------
  !
  ! ... Iterative diagonalization of a complex hermitian matrix
  ! ... through preconditioned RMM-DIIS algorithm.
  !
  USE constants, ONLY : eps14, eps16
  USE kinds,     ONLY : DP
  USE funct,     ONLY : exx_is_active
  USE mp,        ONLY : mp_sum, mp_bcast
  USE mp_bands,  ONLY : inter_bgrp_comm, intra_bgrp_comm, me_bgrp, root_bgrp, &
                        root_bgrp_id, use_bgrp_in_hpsi, set_bgrp_indices
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER,     INTENT(IN)    :: npwx, npw, nbnd, npol
  COMPLEX(DP), INTENT(INOUT) :: psi(npwx*npol,nbnd)
  COMPLEX(DP), INTENT(OUT)   :: spsi(npwx*npol,nbnd)
  REAL(DP),    INTENT(OUT)   :: e(nbnd)
  REAL(DP),    INTENT(IN)    :: g2kin(npwx)
  INTEGER,     INTENT(IN)    :: btype(nbnd)
  REAL(DP),    INTENT(IN)    :: ethr
  INTEGER,     INTENT(IN)    :: ndiis
  LOGICAL,     INTENT(IN)    :: uspp
  INTEGER,     INTENT(OUT)   :: notconv
  INTEGER,     INTENT(OUT)   :: rmm_iter
  !
  ! ... local variables
  !
  INTEGER                  :: ierr
  INTEGER                  :: idiis
  INTEGER                  :: kdim, kdmx
  INTEGER                  :: ibnd_start, ibnd_end, ibnd_size
  REAL(DP)                 :: empty_ethr
  COMPLEX(DP), ALLOCATABLE :: phi(:,:,:), hphi(:,:,:), sphi(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: hpsi(:,:), kpsi(:,:), hkpsi(:,:), skpsi(:,:)
  COMPLEX(DP), ALLOCATABLE :: hc(:,:,:), sc(:,:,:)
  REAL(DP),    ALLOCATABLE :: php(:,:), psp(:,:)
  REAL(DP),    ALLOCATABLE :: ew(:), hw(:), sw(:)
  LOGICAL,     ALLOCATABLE :: conv(:)
  !
  REAL(DP),    PARAMETER   :: SREF = 0.5_DP
  REAL(DP),    PARAMETER   :: SMIN = 0.1_DP
  REAL(DP),    PARAMETER   :: SMAX = 2.0_DP
  !
  CALL start_clock( 'crmmdiagg' )
  !
  empty_ethr = MAX( ( ethr * 5._DP ), 1.E-5_DP )
  !
  IF ( npol == 1 ) THEN
     !
     kdim = npw
     kdmx = npwx
     !
  ELSE
     !
     kdim = npwx * npol
     kdmx = npwx * npol
     !
  END IF
  !
  CALL set_bgrp_indices( nbnd, ibnd_start, ibnd_end )
  !
  ibnd_size = MAX( ibnd_end - ibnd_start + 1, 0 )
  !
  ALLOCATE( phi( kdmx, ibnd_start:ibnd_end, ndiis ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' crmmdiagg ', ' cannot allocate phi ', ABS(ierr) )
  !
  ALLOCATE( hphi( kdmx, ibnd_start:ibnd_end, ndiis ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' crmmdiagg ', ' cannot allocate hphi ', ABS(ierr) )
  !
  IF ( uspp ) THEN
     !
     ALLOCATE( sphi( kdmx, ibnd_start:ibnd_end, ndiis ), STAT=ierr )
     IF( ierr /= 0 ) CALL errore( ' crmmdiagg ', ' cannot allocate sphi ', ABS(ierr) )
     !
  END IF
  !
  ALLOCATE( hpsi( kdmx, nbnd ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' crmmdiagg ', ' cannot allocate hpsi ', ABS(ierr) )
  !
  ALLOCATE( kpsi( kdmx, nbnd ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' crmmdiagg ', ' cannot allocate kpsi ', ABS(ierr) )
  !
  ALLOCATE( hkpsi( kdmx, nbnd ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' crmmdiagg ', ' cannot allocate hkpsi ', ABS(ierr) )
  !
  IF ( uspp ) THEN
     !
     ALLOCATE( skpsi( kdmx, nbnd ), STAT=ierr )
     IF( ierr /= 0 ) CALL errore( ' crmmdiagg ', ' cannot allocate skpsi ', ABS(ierr) )
     !
  END IF
  !
  ALLOCATE( hc( ndiis, ndiis, ibnd_start:ibnd_end ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' crmmdiagg ', ' cannot allocate hc ', ABS(ierr) )
  !
  ALLOCATE( sc( ndiis, ndiis, ibnd_start:ibnd_end ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' crmmdiagg ', ' cannot allocate sc ', ABS(ierr) )
  !
  ALLOCATE( php( ibnd_start:ibnd_end, ndiis ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' crmmdiagg ', ' cannot allocate php ', ABS(ierr) )
  !
  ALLOCATE( psp( ibnd_start:ibnd_end, ndiis ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' crmmdiagg ', ' cannot allocate psp ', ABS(ierr) )
  !
  ALLOCATE( ew( nbnd ) )
  ALLOCATE( hw( nbnd ) )
  ALLOCATE( sw( nbnd ) )
  ALLOCATE( conv( nbnd ) )
  !
  phi  = ZERO
  hphi = ZERO
  IF ( uspp ) sphi = ZERO
  !
  hpsi  = ZERO
  spsi  = ZERO
  kpsi  = ZERO
  hkpsi = ZERO
  IF ( uspp ) skpsi = ZERO
  !
  hc = ZERO
  sc = ZERO
  !
  php = 0._DP
  psp = 0._DP
  !
  e    = 0._DP
  ew   = 0._DP
  hw   = 0._DP
  sw   = 0._DP
  conv = .FALSE.
  !
  ! ... Initial eigenvalues
  !
  CALL eigenvalues( .TRUE. )
  !
  ! ... RMM-DIIS's loop
  !
  rmm_iter = 0
  notconv  = 0
  !
  DO idiis = 1, ndiis
     !
     rmm_iter = rmm_iter + 1
     !
     ! ... Perform DIIS
     !
     CALL do_diis( idiis )
     !
     ! ... Line searching
     !
     CALL line_search( )
     !
     ! ... Calculate eigenvalues and check convergence
     !
     CALL eigenvalues( .FALSE. )
     !
     IF ( notconv == 0 ) EXIT
     !
  END DO
  !
  DEALLOCATE( phi )
  DEALLOCATE( hphi )
  IF ( uspp ) DEALLOCATE( sphi )
  DEALLOCATE( hpsi )
  DEALLOCATE( kpsi )
  DEALLOCATE( hkpsi )
  IF ( uspp ) DEALLOCATE( skpsi )
  DEALLOCATE( hc )
  DEALLOCATE( sc )
  DEALLOCATE( php )
  DEALLOCATE( psp )
  DEALLOCATE( ew )
  DEALLOCATE( hw )
  DEALLOCATE( sw )
  DEALLOCATE( conv )
  !
  CALL stop_clock( 'crmmdiagg' )
  !
  RETURN
  !
  !
CONTAINS
  !
  !
  SUBROUTINE do_diis( idiis )
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: idiis
    !
    INTEGER                  :: ibnd
    INTEGER                  :: kdiis
    REAL(DP)                 :: norm
    COMPLEX(DP)              :: ec
    COMPLEX(DP), ALLOCATABLE :: vec1(:)
    COMPLEX(DP), ALLOCATABLE :: vec2(:,:)
    COMPLEX(DP), ALLOCATABLE :: vc(:)
    COMPLEX(DP), ALLOCATABLE :: tc(:,:)
    !
    ALLOCATE( vec1( kdmx ) )
    ALLOCATE( vec2( kdmx, idiis ) )
    IF ( idiis > 1 ) ALLOCATE( vc( idiis ) )
    ALLOCATE( tc( idiis, ibnd_start:ibnd_end ) )
    !
    ! ... Save current wave functions and eigenvalues
    !
    CALL ZCOPY( kdmx * ibnd_size, psi (1,ibnd_start), 1, phi (1,ibnd_start,idiis), 1 )
    CALL ZCOPY( kdmx * ibnd_size, hpsi(1,ibnd_start), 1, hphi(1,ibnd_start,idiis), 1 )
    IF ( uspp ) &
    CALL ZCOPY( kdmx * ibnd_size, spsi(1,ibnd_start), 1, sphi(1,ibnd_start,idiis), 1 )
    !
    CALL DCOPY( ibnd_size, hw(ibnd_start), 1, php(ibnd_start,idiis), 1 )
    CALL DCOPY( ibnd_size, sw(ibnd_start), 1, psp(ibnd_start,idiis), 1 )
    !
    ! ... <R_i|R_j>
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       ! ... Residual vectors : |R> = (H - e S) |psi>
       !
       DO kdiis = 1, idiis
          !
          ec = CMPLX( php(ibnd,kdiis), 0._DP, kind=DP )
          CALL ZCOPY( kdim, hphi(1,ibnd,kdiis), 1, vec2(1,kdiis), 1 )
          CALL ZAXPY( kdim, -ec, sphi(1,ibnd,kdiis), 1, vec2(1,kdiis), 1 )
          !
       END DO
       !
       ec = CMPLX( php(ibnd,idiis), 0._DP, kind=DP )
       CALL ZCOPY( kdim, hphi(1,ibnd,idiis), 1, vec1(1), 1 )
       CALL ZAXPY( kdim, -ec, sphi(1,ibnd,idiis), 1, vec1(1), 1 )
       !
       CALL ZGEMV( 'C', kdim, idiis, ONE, vec2(1,1), kdmx, &
                   vec1(1), 1, ZERO, hc(1,idiis,ibnd), 1 )
       !
    END DO
    !
    tc(:,:) = hc(:,idiis,:)
    CALL mp_sum( tc, intra_bgrp_comm )
    hc(:,idiis,:) = tc(:,:)
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       hc(idiis,idiis,ibnd) = CMPLX( DBLE( hc(idiis,idiis,ibnd) ), 0._DP, kind=DP )
       hc(idiis,1:idiis,ibnd) = CONJG( hc(1:idiis,idiis,ibnd) )
       !
    END DO
    !
    ! ... <phi_i| S |phi_j>
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       DO kdiis = 1, idiis
          !
          CALL ZCOPY( kdim, phi(1,ibnd,kdiis), 1, vec2(1,kdiis), 1 )
          !
       END DO
       !
       IF ( uspp ) THEN
          !
          CALL ZCOPY( kdim, sphi(1,ibnd,idiis), 1, vec1(1), 1 )
          !
       ELSE
          !
          CALL ZCOPY( kdim, phi(1,ibnd,idiis), 1, vec1(1), 1 )
          !
       END IF
       !
       CALL ZGEMV( 'C', kdim, idiis, ONE, vec2(1,1), kdmx, &
                   vec1(1), 1, ZERO, sc(1,idiis,ibnd), 1 )
       !
    END DO
    !
    tc(:,:) = sc(:,idiis,:)
    CALL mp_sum( tc, intra_bgrp_comm )
    sc(:,idiis,:) = tc(:,:)
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       sc(idiis,idiis,ibnd) = CMPLX( DBLE( sc(idiis,idiis,ibnd) ), 0._DP, kind=DP )
       sc(idiis,1:idiis,ibnd) = CONJG( sc(1:idiis,idiis,ibnd) )
       !
    END DO
    !
    ! ... Update current wave functions and residual vectors
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( idiis > 1 ) THEN
          !
          ! ... solve Rv = eSv
          !
          IF ( me_bgrp == root_bgrp ) CALL diag_diis( ibnd, idiis, vc(:) )
          CALL mp_bcast( vc, root_bgrp, intra_bgrp_comm )
          !
          psi (:,ibnd) = ZERO
          hpsi(:,ibnd) = ZERO
          IF ( uspp ) &
          spsi(:,ibnd) = ZERO
          kpsi(:,ibnd) = ZERO
          !
          DO kdiis = 1, idiis
             !
             ! ... Wave functions
             !
             CALL ZAXPY( kdim, vc(kdiis), phi (1,ibnd,kdiis), 1, psi (1,ibnd), 1 )
             CALL ZAXPY( kdim, vc(kdiis), hphi(1,ibnd,kdiis), 1, hpsi(1,ibnd), 1 )
             IF ( uspp ) &
             CALL ZAXPY( kdim, vc(kdiis), sphi(1,ibnd,kdiis), 1, spsi(1,ibnd), 1 )
             !
             ! ... Residual vectors
             !
             ec = CMPLX( php(ibnd,kdiis), 0._DP, kind=DP )
             CALL ZCOPY( kdim, hphi(1,ibnd,kdiis), 1, vec1(1), 1 )
             CALL ZAXPY( kdim, -ec, sphi(1,ibnd,kdiis), 1, vec1(1), 1 )
             CALL ZAXPY( kdim, vc(kdiis), vec1(1), 1, kpsi(1,ibnd), 1 )
             !
          END DO
          !
       ELSE
          !
          ! ... Wave functions
          !
          norm = SQRT( sw(ibnd) )
          CALL ZDSCAL( kdim, 1._DP / norm, psi (1,ibnd), 1 )
          CALL ZDSCAL( kdim, 1._DP / norm, hpsi(1,ibnd), 1 )
          IF ( uspp ) &
          CALL ZDSCAL( kdim, 1._DP / norm, spsi(1,ibnd), 1 )
          !
          ! ... Residual vectors
          !
          ec = CMPLX( hw(ibnd), 0._DP, kind=DP )
          CALL ZCOPY( kdim, hpsi(1,ibnd), 1, kpsi(1,ibnd), 1 )
          CALL ZAXPY( kdim, -ec, spsi(1,ibnd), 1, kpsi(1,ibnd), 1 )
          !
       END IF
       !
    END DO
    !
    DEALLOCATE( vec1 )
    DEALLOCATE( vec2 )
    IF ( idiis > 1 ) DEALLOCATE( vc )
    DEALLOCATE( tc )
    !
    RETURN
    !
  END SUBROUTINE do_diis
  !
  !
  SUBROUTINE diag_diis( ibnd, idiis, vc )
    !
    IMPLICIT NONE
    !
    INTEGER,     INTENT(IN)  :: ibnd
    INTEGER,     INTENT(IN)  :: idiis
    COMPLEX(DP), INTENT(OUT) :: vc(idiis)
    !
    INTEGER                  :: info
    INTEGER                  :: ndim, ndep
    INTEGER                  :: i, imin
    REAL(DP)                 :: emin
    COMPLEX(DP), ALLOCATABLE :: h1(:,:)
    COMPLEX(DP), ALLOCATABLE :: h2(:,:)
    COMPLEX(DP), ALLOCATABLE :: h3(:,:)
    COMPLEX(DP), ALLOCATABLE :: s1(:,:)
    COMPLEX(DP), ALLOCATABLE :: s2(:,:)
    COMPLEX(DP), ALLOCATABLE :: x1(:,:)
    REAL(DP),    ALLOCATABLE :: e1(:)
    INTEGER                  :: nwork
    COMPLEX(DP), ALLOCATABLE :: work(:)
    !
    ndim  = idiis
    nwork = 3 * ndim
    !
    ALLOCATE( h1( ndim, ndim ) )
    ALLOCATE( h2( ndim, ndim ) )
    ALLOCATE( h3( ndim, ndim ) )
    ALLOCATE( s1( ndim, ndim ) )
    ALLOCATE( s2( ndim, ndim ) )
    ALLOCATE( x1( ndim, ndim ) )
    ALLOCATE( e1( ndim ) )
    ALLOCATE( work( nwork ) )
    !
    h1(1:ndim,1:ndim) = hc(1:ndim,1:ndim,ibnd)
    s1(1:ndim,1:ndim) = sc(1:ndim,1:ndim,ibnd)
    !
    CALL ZHEEV( 'V', 'U', ndim, s1, ndim, e1, work, nwork, info )
    !
    IF( info /= 0 ) CALL errore( ' crmmdiagg ', ' cannot solve diis ', ABS(info) )
    !
    ndep = 0
    !
    DO i = 1, ndim
       !
       IF ( e1(i) > eps14 ) THEN
          !
          s2(:,i) = s1(:,i) / SQRT(e1(i))
          !
       ELSE
          !
          ndep = ndep + 1
          !
          s2(:,i) = ZERO
          !
       END IF
       !
    END DO
    !
    IF ( (ndim - ndep) <= 1 ) THEN
       !
       vc        = ZERO
       vc(idiis) = ONE
       !
       GOTO 10
       !
    END IF
    !
    CALL ZGEMM( 'N', 'C', ndim, ndim, ndim, ONE, s2, ndim, s1, ndim, ZERO, x1, ndim )
    !
    CALL ZGEMM( 'N', 'N', ndim, ndim, ndim, ONE, h1, ndim, x1, ndim, ZERO, h2, ndim )
    !
    CALL ZGEMM( 'N', 'N', ndim, ndim, ndim, ONE, x1, ndim, h2, ndim, ZERO, h3, ndim )
    !
    CALL ZHEEV( 'V', 'U', ndim, h3, ndim, e1, work, nwork, info )
    !
    IF( info /= 0 ) CALL errore( ' crmmdiagg ', ' cannot solve diis ', ABS(info) )
    !
    imin = 1
    emin = e1(1)
    !
    DO i = 2, ndim
       !
       IF ( e1(i) < emin ) imin = i
       !
    END DO
    !
    CALL ZGEMV( 'N', ndim, ndim, ONE, x1, ndim, h3(:,imin), 1, ZERO, vc, 1 )
    !
10  DEALLOCATE( h1 )
    DEALLOCATE( h2 )
    DEALLOCATE( h3 )
    DEALLOCATE( s1 )
    DEALLOCATE( s2 )
    DEALLOCATE( x1 )
    DEALLOCATE( e1 )
    DEALLOCATE( work )
    !
    RETURN
    !
  END SUBROUTINE diag_diis
  !
  !
  SUBROUTINE line_search( )
    !
    IMPLICIT NONE
    !
    INTEGER               :: ibnd, ig, ipol
    REAL(DP)              :: psir, psii, psi2
    REAL(DP)              :: kdiag, k1, k2
    REAL(DP)              :: x, x2, x3, x4
    REAL(DP), ALLOCATABLE :: ekin(:)
    REAL(DP)              :: a, b
    REAL(DP)              :: ene0, ene1
    REAL(DP)              :: step, norm
    REAL(DP)              :: php, khp, khk
    REAL(DP)              :: psp, ksp, ksk
    REAL(DP), ALLOCATABLE :: hmat(:,:), smat(:,:)
    REAL(DP)              :: c1, c2
    COMPLEX(DP)           :: z1, z2
    REAL(DP), ALLOCATABLE :: coef(:,:)
    !
    COMPLEX(DP), EXTERNAL :: ZDOTC
    !
    ALLOCATE( ekin( ibnd_start:ibnd_end ) )
    ALLOCATE( hmat( 3, ibnd_start:ibnd_end ) )
    ALLOCATE( smat( 3, ibnd_start:ibnd_end ) )
    ALLOCATE( coef( 2, ibnd_start:ibnd_end ) )
    !
    ! ... Kinetic energy
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       ekin(ibnd) = 0._DP
       !
       DO ipol = 1, npol
          !
          DO ig = 1, npw
             !
             psir = DBLE ( psi(ig+(ipol-1)*npwx,ibnd) )
             psii = AIMAG( psi(ig+(ipol-1)*npwx,ibnd) )
             psi2 = psir * psir + psii * psii
             ekin(ibnd) = ekin(ibnd) + g2kin(ig) * psi2
             !
          END DO
          !
       END DO
       !
    END DO
    !
    CALL mp_sum( ekin, intra_bgrp_comm )
    !
    ! ... Preconditioning vectors : K (H - e S) |psi>
    !
    ! ... G.Kresse and J.Furthmuller, PRB 54, 11169 (1996)
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       DO ipol = 1, npol
          !
          DO ig = 1, npw
             !
             x  = g2kin(ig) / ( 1.5_DP * ekin(ibnd) )
             x2 = x * x
             x3 = x * x2
             x4 = x * x3
             !
             k1 = 27._DP + 18._DP * x + 12._DP * x2 + 8._DP * x3
             k2 = k1 + 16._DP * x4
             kdiag = ( -4._DP / 3._DP / ekin(ibnd) ) * k1 / k2
             !
             kpsi(ig+(ipol-1)*npwx,ibnd) = kdiag * kpsi(ig+(ipol-1)*npwx,ibnd)
             !
          END DO
          !
       END DO
       !
    END DO
    !
    ! ... Share kpsi for all band-groups
    !
    IF ( ( .NOT. use_bgrp_in_hpsi ) .OR. exx_is_active() ) THEN
       !
       DO ibnd = 1, ( ibnd_start - 1)
          !
          kpsi(:,ibnd) = ZERO
          !
       END DO
       !
       DO ibnd = ( ibnd_end + 1 ), nbnd
          !
          kpsi(:,ibnd) = ZERO
          !
       END DO
       !
       CALL mp_sum( kpsi, inter_bgrp_comm )
       !
    END IF
    !
    ! ... Operate the Hamiltonian : H K (H - eS) |psi>
    !
    CALL h_psi( npwx, npw, nbnd, kpsi, hkpsi )
    !
    ! ... Operate the Overlap : S K (H - eS) |psi>
    !
    IF ( uspp ) CALL s_psi( npwx, npw, nbnd, kpsi, skpsi )
    !
    ! ... Create 2 x 2 matrix
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       php = DBLE( ZDOTC( kdim, psi (1,ibnd), 1, hpsi (1,ibnd), 1 ) )
       khp = DBLE( ZDOTC( kdim, kpsi(1,ibnd), 1, hpsi (1,ibnd), 1 ) )
       khk = DBLE( ZDOTC( kdim, kpsi(1,ibnd), 1, hkpsi(1,ibnd), 1 ) )
       !
       IF ( uspp ) THEN
          !
          psp = DBLE( ZDOTC( kdim, psi (1,ibnd), 1, spsi (1,ibnd), 1 ) )
          ksp = DBLE( ZDOTC( kdim, kpsi(1,ibnd), 1, spsi (1,ibnd), 1 ) )
          ksk = DBLE( ZDOTC( kdim, kpsi(1,ibnd), 1, skpsi(1,ibnd), 1 ) )
          !
       ELSE
          !
          psp = DBLE( ZDOTC( kdim, psi (1,ibnd), 1, psi (1,ibnd), 1 ) )
          ksp = DBLE( ZDOTC( kdim, kpsi(1,ibnd), 1, psi (1,ibnd), 1 ) )
          ksk = DBLE( ZDOTC( kdim, kpsi(1,ibnd), 1, kpsi(1,ibnd), 1 ) )
          !
       END IF
       !
       hmat(1,ibnd) = php
       hmat(2,ibnd) = khp
       hmat(3,ibnd) = khk
       !
       smat(1,ibnd) = psp
       smat(2,ibnd) = ksp
       smat(3,ibnd) = ksk
       !
    END DO
    !
    CALL mp_sum( hmat, intra_bgrp_comm )
    CALL mp_sum( smat, intra_bgrp_comm )
    !
    ! ... Line searching for each band
    !
    IF ( me_bgrp == root_bgrp ) THEN
       !
       DO ibnd = ibnd_start, ibnd_end
          !
          php = hmat(1,ibnd)
          khp = hmat(2,ibnd)
          khk = hmat(3,ibnd)
          !
          psp = smat(1,ibnd)
          ksp = smat(2,ibnd)
          ksk = smat(3,ibnd)
          IF( psp <= eps16 ) CALL errore( ' crmmdiagg ', ' psp <= 0 ', 1 )
          !
          norm = psp + 2._DP * ksp * SREF + ksk * SREF * SREF
          IF( norm <= eps16 ) CALL errore( ' crmmdiagg ', ' norm <= 0 ', 1 )
          !
          ene0 = php / psp
          ene1 = ( php + 2._DP * khp * SREF + khk * SREF * SREF ) / norm
          !
          a = 2._DP * ( khp * psp - php * ksp ) / psp / psp
          b = ( ene1 - ene0 - a * SREF ) / SREF / SREF
          IF( ABS( b ) < eps16 ) CALL errore( ' crmmdiagg ', ' b == 0 ', 1 )
          !
          step  = -0.5_DP * a / b
          step  = MAX( SMIN, step )
          step  = MIN( SMAX, step )
          norm  = psp + 2._DP * ksp * step + ksk * step * step
          IF( norm <= eps16 ) CALL errore( ' crmmdiagg ', ' norm <= 0 ', 1 )
          norm  = SQRT( norm )
          !
          coef(1,ibnd) = 1._DP / norm
          coef(2,ibnd) = step  / norm
          !
          ! ... Update current matrix elements
          !
          c1 = coef(1,ibnd)
          c2 = coef(2,ibnd)
          !
          hw(ibnd) = php * c1 * c1 + 2._DP * khp * c1 * c2 + khk * c2 * c2
          sw(ibnd) = psp * c1 * c1 + 2._DP * ksp * c1 * c2 + ksk * c2 * c2
          !
       END DO
       !
    END IF
    !
    CALL mp_bcast( coef, root_bgrp, intra_bgrp_comm )
    CALL mp_bcast( hw(ibnd_start:ibnd_end), root_bgrp, intra_bgrp_comm )
    CALL mp_bcast( sw(ibnd_start:ibnd_end), root_bgrp, intra_bgrp_comm )
    !
    ! ... Update current wave functions
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       z1 = CMPLX( coef(1,ibnd), 0._DP, kind=DP )
       z2 = CMPLX( coef(2,ibnd), 0._DP, kind=DP )
       !
       CALL ZSCAL( kdim, z1, psi (1,ibnd), 1 )
       CALL ZAXPY( kdim, z2, kpsi(1,ibnd), 1, psi(1,ibnd), 1 )
       !
       CALL ZSCAL( kdim, z1, hpsi (1,ibnd), 1 )
       CALL ZAXPY( kdim, z2, hkpsi(1,ibnd), 1, hpsi(1,ibnd), 1 )
       !
       IF ( uspp ) THEN
          !
          CALL ZSCAL( kdim, z1, spsi (1,ibnd), 1 )
          CALL ZAXPY( kdim, z2, skpsi(1,ibnd), 1, hpsi(1,ibnd), 1 )
          !
       END IF
       !
    END DO
    !
    DEALLOCATE( ekin )
    DEALLOCATE( hmat )
    DEALLOCATE( smat )
    DEALLOCATE( coef )
    !
    RETURN
    !
  END SUBROUTINE line_search
  !
  !
  SUBROUTINE eigenvalues( first )
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: first
    !
    INTEGER :: ibnd
    !
    COMPLEX(DP), EXTERNAL :: ZDOTC
    !
    IF ( first ) THEN
       !
       ! ... Operate the Hamiltonian : H |psi>
       !
       CALL h_psi( npwx, npw, nbnd, psi, hpsi )
       !
       ! ... Operate the Overlap : S |psi>
       !
       IF ( uspp ) CALL s_psi( npwx, npw, nbnd, psi, spsi )
       !
       ! ... Matrix elements
       !
       DO ibnd = ibnd_start, ibnd_end
          !
          hw(ibnd) = DBLE( ZDOTC( kdim, psi(1,ibnd), 1, hpsi(1,ibnd), 1 ) )
          !
       END DO
       !
       CALL mp_sum( hw(ibnd_start:ibnd_end), intra_bgrp_comm )
       !
       IF ( uspp ) THEN
          !
          DO ibnd = ibnd_start, ibnd_end
             !
             sw(ibnd) = DBLE( ZDOTC( kdim, psi(1,ibnd), 1, spsi(1,ibnd), 1 ) )
             !
          END DO
          !
       ELSE
          !
          DO ibnd = ibnd_start, ibnd_end
             !
             sw(ibnd) = DBLE( ZDOTC( kdim, psi(1,ibnd), 1, psi(1,ibnd), 1 ) )
             !
          END DO
          !
       END IF
       !
       CALL mp_sum( sw(ibnd_start:ibnd_end), intra_bgrp_comm )
       !
    END IF
    !
    ! ... Energy eigenvalues
    !
    IF( ANY( sw(ibnd_start:ibnd_end) <= eps16 ) ) &
    CALL errore( ' crmmdiagg ', ' sw <= 0 ', 1 )
    !
    ew(1:nbnd) = 0._DP
    ew(ibnd_start:ibnd_end) = hw(ibnd_start:ibnd_end) / sw(ibnd_start:ibnd_end)
    !
    CALL mp_sum( ew, inter_bgrp_comm )
    !
    ! ... Check convergence
    !
    WHERE( btype(1:nbnd) == 1 )
       !
       conv(1:nbnd) = ( ( ABS( ew(1:nbnd) - e(1:nbnd) ) < ethr ) )
       !
    ELSEWHERE
       !
       conv(1:nbnd) = ( ( ABS( ew(1:nbnd) - e(1:nbnd) ) < empty_ethr ) )
       !
    END WHERE
    !
    CALL mp_bcast( conv, root_bgrp_id, inter_bgrp_comm )
    !
    notconv = COUNT( .NOT. conv(1:nbnd) )
    !
    e(1:nbnd) = ew(1:nbnd)
    !
    RETURN
    !
  END SUBROUTINE eigenvalues
  !
  !
END SUBROUTINE crmmdiagg
