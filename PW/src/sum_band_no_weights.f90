!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE sum_band_no_weights()
  !----------------------------------------------------------------------------
  !
  ! ... Calculates the symmetrized charge density and related quantities
  ! ... Also computes the occupations and the sum of occupied eigenvalues.
  !
  USE kinds,                ONLY : DP
  USE ener,                 ONLY : eband
  USE control_flags,        ONLY : diago_full_acc, gamma_only, lxdm
  USE cell_base,            ONLY : at, bg, omega, tpiba
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE fft_base,             ONLY : dfftp, dffts, dtgs
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE gvect,                ONLY : ngm, g, nl, nlm
  USE gvecs,                ONLY : nls, nlsm, doublegrid
  USE klist,                ONLY : nks, nkstot, wk, xk, ngk, igk_k
  USE fixed_occ,            ONLY : one_atom_occupations
  USE ldaU,                 ONLY : lda_plus_U
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE scf,                  ONLY : rho
  USE symme,                ONLY : sym_rho
  USE io_files,             ONLY : iunwfc, nwordwfc
  USE buffers,              ONLY : get_buffer
  USE uspp,                 ONLY : nkb, vkb, becsum, nhtol, nhtoj, indv, okvan
  USE uspp_param,           ONLY : upf, nh, nhm
  USE wavefunctions_module, ONLY : evc, psic, psic_nc
  USE noncollin_module,     ONLY : noncolin, npol, nspin_mag
  USE spin_orb,             ONLY : lspinorb, domag, fcoef
  USE wvfct,                ONLY : nbnd, npwx, wg, et, btype
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp_bands,             ONLY : inter_bgrp_comm, intra_bgrp_comm, set_bgrp_indices, nbgrp
  USE mp,                   ONLY : mp_sum
  USE funct,                ONLY : dft_is_meta
  USE paw_symmetry,         ONLY : PAW_symmetrize
  USE paw_variables,        ONLY : okpaw
  USE becmod,               ONLY : allocate_bec_type, deallocate_bec_type, &
                                   becp
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  INTEGER :: ir,   &! counter on 3D r points
             is,   &! counter on spin polarizations
             ig,   &! counter on g vectors
             ibnd, & ! counter on bands
             ik,   &! counter on k points
             ibnd_start, ibnd_end, this_bgrp_nbnd ! first, last and number of band in this bgrp
  REAL (DP), ALLOCATABLE :: kplusg (:)
  !
  !
  CALL start_clock( 'sum_band' )
  !
  becsum(:,:,:) = 0.D0
  rho%of_r(:,:)      = 0.D0
  rho%of_g(:,:)      = (0.D0, 0.D0)
  if ( dft_is_meta() .OR. lxdm ) then
     rho%kin_r(:,:)      = 0.D0
     rho%kin_g(:,:)      = (0.D0, 0.D0)
  end if
  eband         = 0.D0
  !
  ! ... calculates weights of Kohn-Sham orbitals used in calculation of rho
  ! already been called before, no need.
  !CALL weights ( )
  !
  IF (one_atom_occupations) CALL new_evc()
  !
  IF ( diago_full_acc ) THEN
     !
     ! ... for diagonalization purposes all the bands are considered occupied
     !
     btype(:,:) = 1
     !
  ELSE
     !
     ! ... for diagonalization purposes a band is considered empty when its
     ! ... occupation is less than 1.0 %
     !
     btype(:,:) = 1
     FORALL( ik = 1:nks, wk(ik) > 0.D0 )
        WHERE( wg(:,ik) / wk(ik) < 0.01D0 ) btype(:,ik) = 0
     END FORALL
     !
  END IF
  !
  ! ... Needed for LDA+U: compute occupations of Hubbard states
  !
  IF (lda_plus_u) THEN
     IF(noncolin) THEN
        CALL new_ns_nc(rho%ns_nc)
     ELSE
        CALL new_ns(rho%ns)
     ENDIF
  ENDIF
  !
  call set_bgrp_indices ( nbnd, ibnd_start, ibnd_end )
  this_bgrp_nbnd = ibnd_end - ibnd_start + 1
  !
  ! ... Allocate (and later deallocate) arrays needed in specific cases
  !
  IF ( okvan ) CALL allocate_bec_type (nkb,nbnd, becp,intra_bgrp_comm)
  IF (dft_is_meta() .OR. lxdm) ALLOCATE (kplusg(npwx))
  !
  ! ... specialized routines are called to sum at Gamma or for each k point 
  ! ... the contribution of the wavefunctions to the charge
  !
  IF ( gamma_only ) THEN
     !
     CALL sum_band_gamma()
     !
  ELSE
     !
     CALL sum_band_k()
     !
  END IF
  !
  IF (dft_is_meta() .OR. lxdm) DEALLOCATE (kplusg)
  !
  IF( okpaw )  THEN
     rho%bec(:,:,:) = becsum(:,:,:) ! becsum is filled in sum_band_{k|gamma}
     ! rho%bec has to be recollected and symmetrized, becsum must not, otherwise
     ! it will break stress routines.
     CALL mp_sum(rho%bec, inter_pool_comm )
     call mp_sum(rho%bec, inter_bgrp_comm )
     CALL PAW_symmetrize(rho%bec)
  ENDIF
  !
  IF ( okvan ) CALL deallocate_bec_type ( becp )
  !
  ! ... If a double grid is used, interpolate onto the fine grid
  !
  IF ( doublegrid ) THEN
     !
     DO is = 1, nspin
        !
        CALL interpolate( rho%of_r(1,is), rho%of_r(1,is), 1 )
        if (dft_is_meta() .OR. lxdm) CALL interpolate(rho%kin_r(1,is),rho%kin_r(1,is),1)
        !
     END DO
     !
  END IF
  !
  ! ... Here we add the Ultrasoft contribution to the charge
  !
  CALL addusdens(rho%of_r(:,:)) ! okvan is checked inside the routine
  !
  IF( okvan )  THEN
     ! bgrp_parallelization is done here but not in subsequent routines
     ! (in particular stress routines uses becsum). collect it across bgrp
     call mp_sum(becsum, inter_bgrp_comm )
  ENDIF
  IF ( noncolin .AND. .NOT. domag ) rho%of_r(:,2:4)=0.D0
  !
  CALL mp_sum( eband, inter_pool_comm )
  CALL mp_sum( eband, inter_bgrp_comm )
  !
  ! ... reduce charge density across pools
  !
  CALL mp_sum( rho%of_r, inter_pool_comm )
  CALL mp_sum( rho%of_r, inter_bgrp_comm )
  if (dft_is_meta() .OR. lxdm) CALL mp_sum( rho%kin_r, inter_pool_comm )
  if (dft_is_meta() .OR. lxdm) CALL mp_sum( rho%kin_r, inter_bgrp_comm )
  !
  ! ... bring the (unsymmetrized) rho(r) to G-space (use psic as work array)
  !
  DO is = 1, nspin
     psic(:) = rho%of_r(:,is)
     CALL fwfft ('Dense', psic, dfftp)
     rho%of_g(:,is) = psic(nl(:))
  END DO
  !
  ! ... symmetrize rho(G) 
  !
  CALL sym_rho ( nspin_mag, rho%of_g )
  !
  ! ... same for rho_kin(G)
  !
  IF ( dft_is_meta() .OR. lxdm) THEN
     DO is = 1, nspin
        psic(:) = rho%kin_r(:,is)
        CALL fwfft ('Dense', psic, dfftp)
        rho%kin_g(:,is) = psic(nl(:))
     END DO
     IF (.NOT. gamma_only) CALL sym_rho( nspin, rho%kin_g )
  END IF
  !
  ! ... synchronize rho%of_r to the calculated rho%of_g (use psic as work array)
  !
  DO is = 1, nspin_mag
     !
     psic(:) = ( 0.D0, 0.D0 )
     psic(nl(:)) = rho%of_g(:,is)
     IF ( gamma_only ) psic(nlm(:)) = CONJG( rho%of_g(:,is) )
     CALL invfft ('Dense', psic, dfftp)
     rho%of_r(:,is) = psic(:)
     !
  END DO
  !
  ! ... the same for rho%kin_r and rho%kin_g
  !
  IF ( dft_is_meta() .OR. lxdm) THEN
     DO is = 1, nspin
        !
        psic(:) = ( 0.D0, 0.D0 )
        psic(nl(:)) = rho%kin_g(:,is)
        IF ( gamma_only ) psic(nlm(:)) = CONJG( rho%kin_g(:,is) )
        CALL invfft ('Dense', psic, dfftp)
        rho%kin_r(:,is) = psic(:)
        !
     END DO
  END IF
  !
  CALL stop_clock( 'sum_band' )
  !
  RETURN
  !
  CONTAINS
     !
     ! ... internal procedures
     !
     !-----------------------------------------------------------------------
     SUBROUTINE sum_band_gamma()
       !-----------------------------------------------------------------------
       !
       ! ... gamma version
       !
       USE becmod,        ONLY : becp
       USE mp_bands,      ONLY : me_bgrp
       USE mp,            ONLY : mp_sum, mp_get_comm_null
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       REAL(DP) :: w1, w2
         ! weights
       INTEGER  :: npw, idx, ioff, incr, v_siz, j
       COMPLEX(DP), ALLOCATABLE :: tg_psi(:)
       REAL(DP),    ALLOCATABLE :: tg_rho(:)
       !
       !
       ! ... here we sum for each k point the contribution
       ! ... of the wavefunctions to the charge
       !
       incr = 2
       !
       IF( dtgs%have_task_groups ) THEN
          !
          IF( dft_is_meta() .OR. lxdm) &
             CALL errore( ' sum_band ', ' task groups with meta dft, not yet implemented ', 1 )
          !
          v_siz = dtgs%tg_nnr * dtgs%nogrp
          !
          ALLOCATE( tg_psi( v_siz ) )
          ALLOCATE( tg_rho( v_siz ) )
          !
          incr  = 2 * dtgs%nogrp
          !
       END IF
       !
       k_loop: DO ik = 1, nks
          !
          IF( dtgs%have_task_groups ) tg_rho = 0.0_DP
          IF ( lsda ) current_spin = isk(ik)
          !
          npw = ngk(ik)
          !
          IF ( nks > 1 ) &
             CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
          !
          IF ( nkb > 0 ) &
             CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb )
          !
          ! ... here we compute the band energy: the sum of the eigenvalues
          !
          DO ibnd = ibnd_start, ibnd_end
             !
             ! ... the sum of eband and demet is the integral for
             ! ... e < ef of e n(e) which reduces for degauss=0 to the sum of
             ! ... the eigenvalues.
             !
             eband = eband + et(ibnd,ik) * wg(ibnd,ik)
             !
          END DO
          !
          DO ibnd = ibnd_start, ibnd_end, incr
             !
             IF( dtgs%have_task_groups ) THEN
                !
                tg_psi(:) = ( 0.D0, 0.D0 )
                ioff   = 0
                !
                DO idx = 1, 2*dtgs%nogrp, 2
                   !
                   ! ... 2*dtgs%nogrp ffts at the same time
                   !
                   IF( idx + ibnd - 1 < ibnd_end ) THEN
                      DO j = 1, npw
                         tg_psi(nls (j)+ioff)=     evc(j,idx+ibnd-1)+&
                              (0.0d0,1.d0) * evc(j,idx+ibnd)
                         tg_psi(nlsm(j)+ioff)=CONJG(evc(j,idx+ibnd-1) -&
                              (0.0d0,1.d0) * evc(j,idx+ibnd) )
                      END DO
                   ELSE IF( idx + ibnd - 1 == ibnd_end ) THEN
                      DO j = 1, npw
                         tg_psi(nls (j)+ioff)=       evc(j,idx+ibnd-1)
                         tg_psi(nlsm(j)+ioff)=CONJG( evc(j,idx+ibnd-1) )
                      END DO
                   END IF

                   ioff = ioff + dtgs%tg_nnr

                END DO
                !
                CALL invfft ('Wave', tg_psi, dffts, dtgs)
                !
                ! Now the first proc of the group holds the first two bands
                ! of the 2*dffts%nogrp bands that we are processing at the same time,
                ! the second proc. holds the third and fourth band
                ! and so on
                !
                ! Compute the proper factor for each band
                !
                DO idx = 1, dtgs%nogrp
                   IF( dtgs%nolist( idx ) == me_bgrp ) EXIT
                END DO
                !
                ! Remember two bands are packed in a single array :
                ! proc 0 has bands ibnd   and ibnd+1
                ! proc 1 has bands ibnd+2 and ibnd+3
                ! ....
                !
                idx = 2 * idx - 1
                !
                IF( idx + ibnd - 1 < ibnd_end ) THEN
                   w1 = wg( idx + ibnd - 1, ik) / omega
                   w2 = wg( idx + ibnd    , ik) / omega
                ELSE IF( idx + ibnd - 1 == ibnd_end ) THEN
                   w1 = wg( idx + ibnd - 1, ik) / omega
                   w2 = w1
                ELSE
                   w1 = 0.0d0
                   w2 = w1
                END IF
                !
                CALL get_rho_gamma(tg_rho, dtgs%tg_npp( me_bgrp + 1 ) * &
                                   dffts%nr1x * dffts%nr2x, w1, w2, tg_psi)
                !
             ELSE
                !
                psic(:) = ( 0.D0, 0.D0 )
                !
                IF ( ibnd < ibnd_end ) THEN
                   !
                   ! ... two ffts at the same time
                   !
                   psic(nls(1:npw))  = evc(1:npw,ibnd) + &
                                           ( 0.D0, 1.D0 ) * evc(1:npw,ibnd+1)
                   psic(nlsm(1:npw)) = CONJG( evc(1:npw,ibnd) - &
                                           ( 0.D0, 1.D0 ) * evc(1:npw,ibnd+1) )
                   !
                ELSE
                   !
                   psic(nls (1:npw))  = evc(1:npw,ibnd)
                   psic(nlsm(1:npw)) = CONJG( evc(1:npw,ibnd) )
                   !
                END IF
                !
                CALL invfft ('Wave', psic, dffts)
                !
                w1 = wg(ibnd,ik) / omega
                !
                ! ... increment the charge density ...
                !
                IF ( ibnd < ibnd_end ) THEN
                   !
                   ! ... two ffts at the same time
                   !
                   w2 = wg(ibnd+1,ik) / omega
                   !
                ELSE
                   !
                   w2 = w1
                   !
                END IF
                !
                CALL get_rho_gamma(rho%of_r(:,current_spin), dffts%nnr, w1, w2, psic)
                !
             END IF
             !
             IF (dft_is_meta() .OR. lxdm) THEN
                DO j=1,3
                   psic(:) = ( 0.D0, 0.D0 )
                   !
                   kplusg (1:npw) = (xk(j,ik)+g(j,1:npw)) * tpiba

                   IF ( ibnd < ibnd_end ) THEN
                      ! ... two ffts at the same time
                      psic(nls (1:npw))=CMPLX(0d0, kplusg(1:npw),kind=DP) * &
                                            ( evc(1:npw,ibnd) + &
                                            ( 0.D0, 1.D0 ) * evc(1:npw,ibnd+1) )
                      psic(nlsm(1:npw)) = CMPLX(0d0, -kplusg(1:npw),kind=DP) * &
                                       CONJG( evc(1:npw,ibnd) - &
                                            ( 0.D0, 1.D0 ) * evc(1:npw,ibnd+1) )
                   ELSE
                      psic(nls(1:npw)) = CMPLX(0d0, kplusg(1:npw),kind=DP) * &
                                              evc(1:npw,ibnd)
                      psic(nlsm(1:npw)) = CMPLX(0d0, -kplusg(1:npw),kind=DP) * &
                                       CONJG( evc(1:npw,ibnd) )
                   END IF
                   !
                   CALL invfft ('Wave', psic, dffts)
                   !
                   ! ... increment the kinetic energy density ...
                   !
                   DO ir = 1, dffts%nnr
                      rho%kin_r(ir,current_spin) = &
                                           rho%kin_r(ir,current_spin) + &
                                           w1 *  DBLE( psic(ir) )**2 + &
                                           w2 * AIMAG( psic(ir) )**2
                   END DO
                   !
                END DO
             END IF
             !
             !
          END DO
          !
          IF( dtgs%have_task_groups ) THEN
             !
             ! reduce the group charge
             !
             CALL mp_sum( tg_rho, gid = dtgs%ogrp_comm )
             !
             ioff = 0
             DO idx = 1, dtgs%nogrp
                IF( me_bgrp == dtgs%nolist( idx ) ) EXIT
                ioff = ioff + dffts%nr1x * dffts%nr2x * dffts%npp( dtgs%nolist( idx ) + 1 )
             END DO
             !
             ! copy the charge back to the processor location
             !
             DO ir = 1, dffts%nnr
                rho%of_r(ir,current_spin) = rho%of_r(ir,current_spin) + tg_rho(ir+ioff)
             END DO

          END IF
          !
          ! ... If we have a US pseudopotential we compute here the becsum term
          !
          IF ( okvan ) CALL sum_bec ( ik, current_spin, ibnd_start,ibnd_end,this_bgrp_nbnd ) 
          !
       END DO k_loop
       !
       ! ... with distributed <beta|psi>, sum over bands
       !
       IF( okvan .AND. becp%comm /= mp_get_comm_null() ) CALL mp_sum( becsum, becp%comm )
       !
       IF( dtgs%have_task_groups ) THEN
          DEALLOCATE( tg_psi )
          DEALLOCATE( tg_rho )
       END IF
       !
       RETURN
       !
     END SUBROUTINE sum_band_gamma
     !
     !
     !-----------------------------------------------------------------------
     SUBROUTINE sum_band_k()
       !-----------------------------------------------------------------------
       !
       ! ... k-points version
       !
       USE mp_bands,     ONLY : me_bgrp
       USE mp,           ONLY : mp_sum, mp_get_comm_null
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       REAL(DP) :: w1
       ! weights
       INTEGER :: npw, ipol, na, np
       !
       INTEGER  :: idx, ioff, incr, v_siz, j
       COMPLEX(DP), ALLOCATABLE :: tg_psi(:), tg_psi_nc(:,:)
       REAL(DP),    ALLOCATABLE :: tg_rho(:), tg_rho_nc(:,:)
       LOGICAL  :: use_tg
       !
       ! ... here we sum for each k point the contribution
       ! ... of the wavefunctions to the charge
       !
       use_tg = ( dtgs%have_task_groups ) .AND. &
                                ( .NOT. (dft_is_meta() .OR. lxdm) )
       !
       incr = 1
       !
       IF( use_tg ) THEN
          !
          v_siz = dtgs%tg_nnr * dtgs%nogrp
          !
          IF (noncolin) THEN
             ALLOCATE( tg_psi_nc( v_siz, npol ) )
             ALLOCATE( tg_rho_nc( v_siz, nspin_mag ) )
          ELSE
             ALLOCATE( tg_psi( v_siz ) )
             ALLOCATE( tg_rho( v_siz ) )
          ENDIF
          !
          incr  = dtgs%nogrp
          !
       END IF
       !
       k_loop: DO ik = 1, nks
          !
          IF( use_tg ) THEN
            IF (noncolin) THEN
               tg_rho_nc = 0.0_DP
            ELSE
               tg_rho = 0.0_DP
            ENDIF
          ENDIF

          IF ( lsda ) current_spin = isk(ik)
          npw = ngk (ik)
          !
          IF ( nks > 1 ) &
             CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
          !
          IF ( nkb > 0 ) &
             CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb )
          !
          ! ... here we compute the band energy: the sum of the eigenvalues
          !
          DO ibnd = ibnd_start, ibnd_end, incr
             !
             IF( use_tg ) THEN
                DO idx = 1, dtgs%nogrp
                   IF( idx + ibnd - 1 <= ibnd_end ) eband = eband + et( idx + ibnd - 1, ik ) * wg( idx + ibnd - 1, ik )
                END DO
             ELSE
                eband = eband + et( ibnd, ik ) * wg( ibnd, ik )
             END IF
             !
             ! ... the sum of eband and demet is the integral for e < ef of
             ! ... e n(e) which reduces for degauss=0 to the sum of the
             ! ... eigenvalues
             w1 = wg(ibnd,ik) / omega
             !
             IF (noncolin) THEN
                IF( use_tg ) THEN
                   !
                   tg_psi_nc = ( 0.D0, 0.D0 )
                   !
                   ioff   = 0
                   !
                   DO idx = 1, dtgs%nogrp
                      !
                      ! ... dffts%nogrp ffts at the same time
                      !
                      IF( idx + ibnd - 1 <= ibnd_end ) THEN
                         DO j = 1, npw
                            tg_psi_nc( nls(igk_k(j,ik) ) + ioff, 1 ) = &
                                                       evc( j, idx+ibnd-1 )
                            tg_psi_nc( nls(igk_k(j,ik) ) + ioff, 2 ) = &
                                                       evc( j+npwx, idx+ibnd-1 )
                         END DO
                      END IF

                      ioff = ioff + dtgs%tg_nnr

                   END DO
                   !
                   CALL invfft ('Wave', tg_psi_nc(:,1), dffts, dtgs)
                   CALL invfft ('Wave', tg_psi_nc(:,2), dffts, dtgs)
                   !
                   ! Now the first proc of the group holds the first band
                   ! of the dffts%nogrp bands that we are processing at the same time,
                   ! the second proc. holds the second and so on
                   !
                   ! Compute the proper factor for each band
                   !
                   DO idx = 1, dtgs%nogrp
                      IF( dtgs%nolist( idx ) == me_bgrp ) EXIT
                   END DO
                   !
                   ! Remember
                   ! proc 0 has bands ibnd
                   ! proc 1 has bands ibnd+1
                   ! ....
                   !
                   IF( idx + ibnd - 1 <= ibnd_end ) THEN
                      w1 = wg( idx + ibnd - 1, ik) / omega
                   ELSE
                      w1 = 0.0d0
                   END IF
                   !
                   DO ipol=1,npol
                      CALL get_rho(tg_rho_nc(:,1), dtgs%tg_npp( me_bgrp + 1 ) &
                          * dffts%nr1x * dffts%nr2x, w1, tg_psi_nc(:,ipol))
                   ENDDO
                   !
                   IF (domag) CALL get_rho_domag(tg_rho_nc(:,:), &
                          dtgs%tg_npp( me_bgrp + 1 )*dffts%nr1x*dffts%nr2x, &
                          w1, tg_psi_nc(:,:))
                   !
                ELSE
!
!     Noncollinear case without task groups
!
                   psic_nc = (0.D0,0.D0)
                   DO ig = 1, npw
                      psic_nc(nls(igk_k(ig,ik)),1)=evc(ig     ,ibnd)
                      psic_nc(nls(igk_k(ig,ik)),2)=evc(ig+npwx,ibnd)
                   END DO
                   CALL invfft ('Wave', psic_nc(:,1), dffts)
                   CALL invfft ('Wave', psic_nc(:,2), dffts)
                   !
                   ! increment the charge density ...
                   !
                   DO ipol=1,npol
                      CALL get_rho(rho%of_r(:,1), dffts%nnr, w1, psic_nc(:,ipol))
                   END DO
                   !
                   ! In this case, calculate also the three
                   ! components of the magnetization (stored in rho%of_r(ir,2-4))
                   !
                   IF (domag) THEN
                      CALL get_rho_domag(rho%of_r(:,:), dffts%nnr, w1, psic_nc(:,:))
                   ELSE
                      rho%of_r(:,2:4)=0.0_DP
                   END IF
                   !
                END IF
                !
             ELSE
                !
                IF( use_tg ) THEN
                   !
!$omp parallel default(shared), private(j,ioff,idx)
!$omp do
                   DO j = 1, SIZE( tg_psi )
                      tg_psi(j) = ( 0.D0, 0.D0 )
                   END DO
!$omp end do
                   !
                   ioff   = 0
                   !
                   DO idx = 1, dtgs%nogrp
                      !
                      ! ... dffts%nogrp ffts at the same time
                      !
                      IF( idx + ibnd - 1 <= ibnd_end ) THEN
!$omp do
                         DO j = 1, npw
                            tg_psi( nls(igk_k(j,ik))+ioff ) = evc(j,idx+ibnd-1)
                         END DO
!$omp end do
                      END IF

                      ioff = ioff + dtgs%tg_nnr

                   END DO
!$omp end parallel
                   !
                   CALL invfft ('Wave', tg_psi, dffts, dtgs)
                   !
                   ! Now the first proc of the group holds the first band
                   ! of the dffts%nogrp bands that we are processing at the same time,
                   ! the second proc. holds the second and so on
                   !
                   ! Compute the proper factor for each band
                   !
                   DO idx = 1, dtgs%nogrp
                      IF( dtgs%nolist( idx ) == me_bgrp ) EXIT
                   END DO
                   !
                   ! Remember
                   ! proc 0 has bands ibnd
                   ! proc 1 has bands ibnd+1
                   ! ....
                   !
                   IF( idx + ibnd - 1 <= ibnd_end ) THEN
                      w1 = wg( idx + ibnd - 1, ik) / omega
                   ELSE
                      w1 = 0.0d0
                   END IF
                   !
                   CALL get_rho(tg_rho, dtgs%tg_npp( me_bgrp + 1 ) * dffts%nr1x * dffts%nr2x, w1, tg_psi)
                   !
                ELSE
                   !
                   psic(:) = ( 0.D0, 0.D0 )
                   !
                   psic(nls(igk_k(1:npw,ik))) = evc(1:npw,ibnd)
                   !
                   CALL invfft ('Wave', psic, dffts)
                   !
                   ! ... increment the charge density ...
                   !
                   CALL get_rho(rho%of_r(:,current_spin), dffts%nnr, w1, psic)

                END IF
                !
                IF (dft_is_meta() .OR. lxdm) THEN
                   DO j=1,3
                      psic(:) = ( 0.D0, 0.D0 )
                      !
                      kplusg (1:npw) = (xk(j,ik)+g(j,igk_k(1:npw,ik))) * tpiba
                      psic(nls(igk_k(1:npw,ik)))=CMPLX(0d0,kplusg(1:npw),kind=DP) * &
                                              evc(1:npw,ibnd)
                      !
                      CALL invfft ('Wave', psic, dffts)
                      !
                      ! ... increment the kinetic energy density ...
                      !
                      CALL get_rho(rho%kin_r(:,current_spin), dffts%nnr, w1, psic)
                   END DO
                END IF
                !
             END IF
             !
          END DO
          !
          IF( use_tg ) THEN
             !
             ! reduce the group charge
             !
             IF (noncolin) THEN
                CALL mp_sum( tg_rho_nc, gid = dtgs%ogrp_comm )
             ELSE
                CALL mp_sum( tg_rho, gid = dtgs%ogrp_comm )
             ENDIF
             !
             ioff = 0
             DO idx = 1, dtgs%nogrp
                IF( me_bgrp == dtgs%nolist( idx ) ) EXIT
                ioff = ioff + dffts%nr1x * dffts%nr2x * dffts%npp( dtgs%nolist( idx ) + 1 )
             END DO
             !
             ! copy the charge back to the proper processor location
             !
             IF (noncolin) THEN
!$omp parallel do
                DO ir = 1, dffts%nnr
                   rho%of_r(ir,1) = rho%of_r(ir,1) + &
                                               tg_rho_nc(ir+ioff,1)
                END DO
!$omp end parallel do
                IF (domag) THEN
!$omp parallel do
                   DO ipol=2,4 
                      DO ir = 1, dffts%nnr
                         rho%of_r(ir,ipol) = rho%of_r(ir,ipol) + &
                                               tg_rho_nc(ir+ioff,ipol)
                      END DO
                   END DO
!$omp end parallel do
                ENDIF 
             ELSE
!$omp parallel do
                DO ir = 1, dffts%nnr
                   rho%of_r(ir,current_spin) = rho%of_r(ir,current_spin) + tg_rho(ir+ioff)
                END DO
!$omp end parallel do
             END IF
             !
          END IF
          !
          ! ... If we have a US pseudopotential we compute here the becsum term
          !
          IF ( okvan ) CALL sum_bec ( ik, current_spin, ibnd_start,ibnd_end,this_bgrp_nbnd ) 
          !
       END DO k_loop
       !
       ! ... with distributed <beta|psi>, sum over bands
       !
       IF( okvan .AND. becp%comm /= mp_get_comm_null() ) CALL mp_sum( becsum, becp%comm )
       !

       IF( use_tg ) THEN
          IF (noncolin) THEN
             DEALLOCATE( tg_psi_nc )
             DEALLOCATE( tg_rho_nc )
          ELSE
             DEALLOCATE( tg_psi )
             DEALLOCATE( tg_rho )
          END IF
       END IF
       !
       RETURN
       !
     END SUBROUTINE sum_band_k
     !
     !
     SUBROUTINE get_rho(rho_loc, nrxxs_loc, w1_loc, psic_loc)

        IMPLICIT NONE

        INTEGER :: nrxxs_loc
        REAL(DP) :: rho_loc(nrxxs_loc)
        REAL(DP) :: w1_loc
        COMPLEX(DP) :: psic_loc(nrxxs_loc)

        INTEGER :: ir

!$omp parallel do
        DO ir = 1, nrxxs_loc
           !
           rho_loc(ir) = rho_loc(ir) + &
                         w1_loc * ( DBLE( psic_loc(ir) )**2 + &
                                   AIMAG( psic_loc(ir) )**2 )
           !
        END DO
!$omp end parallel do

     END SUBROUTINE get_rho

     SUBROUTINE get_rho_gamma(rho_loc, nrxxs_loc, w1_loc, w2_loc, psic_loc)

        IMPLICIT NONE

        INTEGER :: nrxxs_loc
        REAL(DP) :: rho_loc(nrxxs_loc)
        REAL(DP) :: w1_loc, w2_loc
        COMPLEX(DP) :: psic_loc(nrxxs_loc)

        INTEGER :: ir

!$omp parallel do
        DO ir = 1, nrxxs_loc
           !
           rho_loc(ir) = rho_loc(ir) + &
                         w1_loc * DBLE( psic_loc(ir) )**2 + &
                         w2_loc * AIMAG( psic_loc(ir) )**2
           !
        END DO
!$omp end parallel do

     END SUBROUTINE get_rho_gamma


     SUBROUTINE get_rho_domag(rho_loc, nrxxs_loc, w1_loc, psic_loc)

        IMPLICIT NONE

        INTEGER :: nrxxs_loc
        REAL(DP) :: rho_loc(:, :)
        REAL(DP) :: w1_loc
        COMPLEX(DP) :: psic_loc(:, :)

        INTEGER :: ir

!$omp parallel do
        DO ir = 1, nrxxs_loc
           !
           rho_loc(ir,2) = rho_loc(ir,2) + w1_loc*2.D0* &
                          (DBLE(psic_loc(ir,1))* DBLE(psic_loc(ir,2)) + &
                          AIMAG(psic_loc(ir,1))*AIMAG(psic_loc(ir,2)))
 
           rho_loc(ir,3) = rho_loc(ir,3) + w1_loc*2.D0* &
                          (DBLE(psic_loc(ir,1))*AIMAG(psic_loc(ir,2)) - &
                           DBLE(psic_loc(ir,2))*AIMAG(psic_loc(ir,1)))

           rho_loc(ir,4) = rho_loc(ir,4) + w1_loc* &
                          (DBLE(psic_loc(ir,1))**2+AIMAG(psic_loc(ir,1))**2 &
                          -DBLE(psic_loc(ir,2))**2-AIMAG(psic_loc(ir,2))**2)
           !
        END DO
!$omp end parallel do

     END SUBROUTINE get_rho_domag

END SUBROUTINE sum_band_no_weights

