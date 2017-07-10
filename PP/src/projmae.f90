!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
PROGRAM do_projmae
  !-----------------------------------------------------------------------
  !
  ! projects wavefunctions onto orthogonalized atomic wavefunctions,
  ! calculates Lowdin charges, spilling parameter, projected DOS
  ! or computes the LDOS in a volume given in input as function of energy
  !
  ! See files INPUT_PROJWFC.* in Doc/ directory for usage
  ! IMPORTANT: since v.5 namelist name is &projwfc and no longer &inputpp
  !
  USE parameters, ONLY : npk
  USE io_global,  ONLY : stdout, ionode, ionode_id
  USE constants,  ONLY : rytoev
  USE kinds,      ONLY : DP
  USE klist,      ONLY : nks, nkstot, xk, degauss, ngauss, lgauss, ltetra
  USE io_files,   ONLY : nd_nmbr, prefix, tmp_dir
  USE noncollin_module, ONLY : noncolin
  USE mp,         ONLY : mp_bcast
  USE spin_orb,   ONLY: lforcet
  USE mp_world,   ONLY : world_comm
  USE mp_global,  ONLY : mp_startup, nproc_ortho, nproc_pool, nproc_pool_file
  USE environment,ONLY : environment_start, environment_end
  USE wvfct,      ONLY : et, nbnd
  USE basis,      ONLY : natomwfc
  USE control_flags, ONLY: twfcollect
  USE paw_variables, ONLY : okpaw
  ! following modules needed for generation of tetrahedra
  USE ktetra,     ONLY : tetra, tetra_type, opt_tetra_init
  USE symm_base,  ONLY : nsym, s, time_reversal, t_rev
  USE cell_base,  ONLY : at, bg
  USE start_k,    ONLY : k1, k2, k3, nk1, nk2, nk3
  USE lsda_mod,   ONLY : lsda
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  CHARACTER (len=256) :: filpdos, filproj, outdir
  REAL (DP), allocatable :: xk_collect(:,:)
  REAL (DP) :: Emin, Emax, DeltaE, degauss1, ef_0
  INTEGER :: nks2, ngauss1, ios
  LOGICAL :: lwrite_overlaps, lbinary_data
  LOGICAL :: lsym, kresolveddos, tdosinboxes, plotboxes, pawproj
  INTEGER, PARAMETER :: N_MAX_BOXES = 999
  INTEGER :: n_proj_boxes, irmin(3,N_MAX_BOXES), irmax(3,N_MAX_BOXES)
  LOGICAL :: lgww  !if .true. use GW QP energies from file bands.dat
  !
  NAMELIST / projwfc / outdir, prefix, ngauss, degauss, lsym, &
             Emin, Emax, DeltaE, filpdos, filproj, lgww, &
             kresolveddos, tdosinboxes, n_proj_boxes, irmin, irmax, plotboxes, &
             lwrite_overlaps, lbinary_data, pawproj, lforcet, ef_0
  !
  ! initialise environment
  !
#if defined(__MPI)
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'PROJMAE' )
  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( trim( outdir ) == ' ' ) outdir = './'
  filproj= ' '
  filpdos= ' '
  Emin   =-1000000.d0
  Emax   =+1000000.d0
  DeltaE = 0.01d0
  ngauss = 0
  lsym   = .true.
  degauss= 0.d0
  lgww   = .false.
  pawproj= .false.
  lwrite_overlaps   = .false.
  lbinary_data = .false.
  kresolveddos = .false.
  tdosinboxes = .false.
  plotboxes   = .false.
  n_proj_boxes= 1
  irmin(:,:)  = 1
  irmax(:,:)  = 0
  !
  ios = 0
  !

  ef_0 = 0.d0
  lforcet = .false.


  IF ( ionode )  THEN
     !
     CALL input_from_file ( )
     !
     READ (5, projwfc, iostat = ios)
     !
     tmp_dir = trimcheck (outdir)
     ! save the value of degauss and ngauss: they are read from file
     degauss1=degauss
     ngauss1 = ngauss
     !
  ENDIF
  !
  CALL mp_bcast (ios, ionode_id, world_comm )
  IF (ios /= 0) CALL errore ('do_projwfc', 'reading projwfc namelist', abs (ios) )
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir,   ionode_id, world_comm )
  CALL mp_bcast( prefix,    ionode_id, world_comm )
  CALL mp_bcast( filproj,   ionode_id, world_comm )
  CALL mp_bcast( ngauss1,   ionode_id, world_comm )
  CALL mp_bcast( degauss1,  ionode_id, world_comm )
  CALL mp_bcast( DeltaE,    ionode_id, world_comm )
  CALL mp_bcast( lsym,      ionode_id, world_comm )
  CALL mp_bcast( Emin,      ionode_id, world_comm )
  CALL mp_bcast( Emax,      ionode_id, world_comm )
  CALL mp_bcast( lwrite_overlaps, ionode_id, world_comm )
  CALL mp_bcast( lbinary_data,    ionode_id, world_comm )
  CALL mp_bcast( lgww,      ionode_id, world_comm )
  CALL mp_bcast( pawproj,   ionode_id, world_comm )
  CALL mp_bcast( tdosinboxes,     ionode_id, world_comm )
  CALL mp_bcast( n_proj_boxes,    ionode_id, world_comm )
  CALL mp_bcast( irmin,     ionode_id, world_comm )
  CALL mp_bcast( irmax,     ionode_id, world_comm )
  CALL mp_bcast( ef_0, ionode_id, world_comm )
  CALL mp_bcast( lforcet, ionode_id, world_comm )


  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  CALL read_file ( )
  !
  IF (pawproj) THEN
    IF ( .NOT. okpaw ) CALL errore ('projwfc','option pawproj only for PAW',1)
    IF ( noncolin )  CALL errore ('projwfc','option pawproj and noncolinear spin not implemented',2)
  END IF
  !
  IF (nproc_pool /= nproc_pool_file .and. .not. twfcollect)  &
     CALL errore('projwfc',&
     'pw.x run with a different number of procs/pools. Use wf_collect=.true.',1)
  !
  CALL openfil_pp ( )
  !
  !   Tetrahedron method
  !
  IF ( ltetra .AND. tetra_type == 1 .OR. tetra_type == 2 ) THEN
     !
     ! info on tetrahedra is no longer saved to file and must be rebuilt
     !
     ! workaround for old xml file, to be removed
     IF(ALLOCATED(tetra)) DEALLOCATE(tetra)
     !
     ! in the lsda case, only the first half of the k points
     ! are needed in the input of "tetrahedra"
     !
     IF ( lsda ) THEN
        nks2 = nkstot / 2
     ELSE
        nks2 = nkstot
     END IF
     IF(tetra_type == 1) THEN
        WRITE( stdout,'(/5x,"Linear tetrahedron method (read from file) ")')
     ELSE
        WRITE( stdout,'(/5x,"Optimized tetrahedron method (read from file) ")')
     END IF
     !
     ! not sure this is needed
     !
     ALLOCATE(xk_collect(3,nkstot))
     CALL poolcollect(3, nks, xk, nkstot, xk_collect)
     !
     CALL opt_tetra_init(nsym, s, time_reversal, t_rev, at, bg, npk, k1,k2,k3, &
          &              nk1, nk2, nk3, nks2, xk_collect, 1)
     !
     DEALLOCATE(xk_collect)
     !
  ELSE IF (degauss1/=0.d0) THEN
     degauss=degauss1
     ngauss =ngauss1
     WRITE( stdout,'(/5x,"Gaussian broadening (read from input): ",&
          &        "ngauss,degauss=",i4,f12.6/)') ngauss,degauss
     lgauss=.true.
  ELSE IF (lgauss) THEN
     WRITE( stdout,'(/5x,"Gaussian broadening (read from file): ",&
          &        "ngauss,degauss=",i4,f12.6/)') ngauss,degauss
  ELSE
     degauss=DeltaE/rytoev
     ngauss =0
     WRITE( stdout,'(/5x,"Gaussian broadening (default values): ",&
          &        "ngauss,degauss=",i4,f12.6/)') ngauss,degauss
     lgauss=.true.
  ENDIF
  !
  IF ( filpdos == ' ') filpdos = prefix
  !
  CALL projwave_nc(filproj,lsym,lwrite_overlaps,lbinary_data,ef_0)
  !
  CALL environment_end ( 'PROJMAE' )
  !
  CALL stop_pp
  !
END PROGRAM do_projmae

!-----------------------------------------------------------------------
SUBROUTINE projwave_nc(filproj, lsym, lwrite_ovp, lbinary, ef_0 )
  !-----------------------------------------------------------------------
  !
  USE io_global,  ONLY : stdout, ionode
  USE ions_base, ONLY : zv, tau, nat, ntyp => nsp, ityp, atm
  USE run_info, ONLY: title
  USE cell_base
  USE constants, ONLY: rytoev, eps4
  USE gvect
  USE gvecs,   ONLY: dual
  USE gvecw,   ONLY: ecutwfc
  USE fft_base, ONLY : dfftp
  USE klist, ONLY: xk, nks, nkstot, nelec, ngk, igk_k
  USE lsda_mod, ONLY: nspin
  USE noncollin_module, ONLY: noncolin, npol, angle1, angle2
  USE symm_base, ONLY: nsym, irt, t_rev
  USE wvfct, ONLY: npwx, nbnd, et, wg
  USE control_flags, ONLY: gamma_only
  USE uspp, ONLY: nkb, vkb
  USE uspp_param, ONLY: upf
  USE becmod,   ONLY: bec_type, becp, calbec, allocate_bec_type, deallocate_bec_type
  USE io_files, ONLY: nd_nmbr, prefix, tmp_dir, nwordwfc, iunwfc
  USE wavefunctions_module, ONLY: evc
  USE mp_global, ONLY : intra_pool_comm
  USE mp,        ONLY : mp_sum
  USE mp_pools,             ONLY : inter_pool_comm
  !
  USE spin_orb,   ONLY: lspinorb, domag, lforcet
  USE projections
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*) :: filproj
  CHARACTER(256) :: file_eband
  LOGICAL :: lwrite_ovp, lbinary
  LOGICAL :: lsym
  LOGICAL :: freeswfcatom
  !
  INTEGER :: ik, ibnd, i, j, k, na, nb, nt, isym, ind, n, m, m1, n1, &
             n2, l, nwfc, nwfc1, lmax_wfc, is, nspin0, iunproj, npw, &
             ind0
  REAL(DP) :: jj, ef_0
  REAL(DP), ALLOCATABLE :: e (:)
  ! Some workspace for k-point calculation ...
  REAL(DP) :: psum
  INTEGER, ALLOCATABLE :: idx(:)
  !
  REAL(DP), EXTERNAL :: get_clock
  INTEGER :: glob_k
  REAL(DP), ALLOCATABLE ::  eband_k(:)
  INTEGER, EXTERNAL :: global_kpoint_index
  !
  !
  !
     IF ( lsym ) call errore('projwave_nc','Force Theorem   &
                     & implemented only with lsym=.false.',1) 
      CALL weights()
!   write(6,*) 'ef_0 = ', ef_0
!   write(6,*) wg
      ef_0 = ef_0 / rytoev
     allocate(eband_k(nkstot))
      eband_k = 0.d0


  !
  !    loop on k points
  !
  DO ik = 1, nks

     CALL davcio (evc, 2*nwordwfc, iunwfc, ik, - 1)

     DO i = 1, nbnd
         psum = wg(i,ik) * (et(i,ik)-ef_0)
         glob_k = global_kpoint_index(nkstot, ik)
         eband_k(glob_k) = eband_k(glob_k) + psum
     ENDDO 
!-- 
     ! on k-points
  ENDDO
  !
!-- Output for the Force Theorem (AlexS)
!
 CALL mp_sum( eband_k,  inter_pool_comm )

IF ( ionode ) THEN

        open(22, file=trim(filproj)//".eband_k.dat",form='formatted', status='unknown')

       eband_k = eband_k*rytoev

    write(22, '("# mae decomposed in kpt, unit is meV, already multiplied by kpt weight")')
    write(22, '("sum(eband_k) = ", f18.10, "  eV")') sum(eband_k(:))
    write(22, '("ef_0 = ", f18.10, "  ry")') ef_0
    do ik = 1, nkstot
      write(22, '(f18.10)') eband_k(ik) * 1000    ! unit is meV
    end do

       CLOSE(22)

 ENDIF

 RETURN

END SUBROUTINE projwave_nc
