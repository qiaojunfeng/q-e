!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
MODULE projmae_module
CONTAINS
!-----------------------------------------------------------------------
SUBROUTINE extract_psi (plot_files,plot_num)
  !-----------------------------------------------------------------------
  !
  !    Reads data produced by pw.x, computes the desired quantity (rho, V, ...)
  !    and writes it to a file (or multiple files) for further processing or
  !    plotting
  !
  !    On return, plot_files contains a list of all written files.
  !
  !    DESCRIPTION of the INPUT: see file Doc/INPUT_PP
  !
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : bg, omega
  USE ener,      ONLY : ef
  USE ions_base, ONLY : nat, ntyp=>nsp, ityp, tau
  USE gvect
  USE fft_base,  ONLY : dfftp
  USE klist,     ONLY : two_fermi_energies, degauss
  USE vlocal,    ONLY : strf
  USE io_files,  ONLY : tmp_dir, prefix
  USE io_global, ONLY : ionode, ionode_id
  USE mp_global,     ONLY : nproc_pool, nproc_file, nproc_pool_file
  USE control_flags, ONLY : twfcollect
  USE noncollin_module, ONLY : i_cons
  USE paw_variables, ONLY : okpaw
  USE mp,        ONLY : mp_bcast
  USE mp_world,  ONLY : world_comm, mpime
  USE constants, ONLY : rytoev
  USE parameters, ONLY : npk
  USE io_global, ONLY : stdout
  USE klist, ONLY : nkstot
  USE wvfct, ONLY : nbnd, et, wg

  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  CHARACTER(len=256), DIMENSION(:), ALLOCATABLE, INTENT(out) :: plot_files
  INTEGER, INTENT(out) :: plot_num

  CHARACTER (len=2), DIMENSION(0:3) :: spin_desc = &
       (/ '  ', '_X', '_Y', '_Z' /)

  INTEGER :: kpoint(2), kband(2), spin_component(3), ios

  REAL(DP) :: emin, emax, sample_bias, z, dz
  
  REAL(DP) :: degauss_ldos, delta_e
  CHARACTER(len=256) :: filplot
  INTEGER :: plot_nkpt, plot_nbnd, plot_nspin, nplots
  INTEGER ::  ikpt, ibnd, ispin
  LOGICAL :: lsign

  ! directory for temporary files
  CHARACTER(len=256) :: outdir

  NAMELIST / inputprojmae / outdir, prefix, plot_num, sample_bias, &
      spin_component, z, dz, emin, emax, delta_e, degauss_ldos, kpoint, kband, &
      filplot, lsign, ef_0
  !
  REAL(DP), ALLOCATABLE :: psi(:), eband_r(:)
  REAL(DP) :: ef_0, eband_tot
  !   set default values for variables in namelist
  REAL(DP), EXTERNAL :: get_clock
  !
  CALL start_clock('extract_psi')
  !
  prefix = 'pwscf'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( trim( outdir ) == ' ' ) outdir = './'
  filplot = 'tmp.pp'
  plot_num = -1
  kpoint(2) = 0
  kband(2) = 0
  spin_component = 0
  sample_bias = 0.01d0
  z = 1.d0
  dz = 0.05d0
  lsign=.false.
  emin = -999.0d0
  emax = +999.0d0
  delta_e=0.1d0
  degauss_ldos=-999.0d0
  ef_0 = 0.0_DP
  !
  ios = 0
  !
  IF ( ionode )  THEN
     !
     !     reading the namelist inputpp
     !
     READ (5, inputprojmae, iostat = ios)
     !
     tmp_dir = trimcheck ( outdir )
     !
  ENDIF
  !
  CALL mp_bcast (ios, ionode_id, world_comm)
  !
  IF ( ios /= 0) CALL errore ('postproc', 'reading inputpp namelist', abs(ios))
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir, ionode_id, world_comm )
  CALL mp_bcast( prefix, ionode_id, world_comm )
  CALL mp_bcast( plot_num, ionode_id, world_comm )
  CALL mp_bcast( sample_bias, ionode_id, world_comm )
  CALL mp_bcast( spin_component, ionode_id, world_comm )
  CALL mp_bcast( z, ionode_id, world_comm )
  CALL mp_bcast( dz, ionode_id, world_comm )
  CALL mp_bcast( emin, ionode_id, world_comm )
  CALL mp_bcast( emax, ionode_id, world_comm )
  CALL mp_bcast( degauss_ldos, ionode_id, world_comm )
  CALL mp_bcast( delta_e, ionode_id, world_comm )
  CALL mp_bcast( kband, ionode_id, world_comm )
  CALL mp_bcast( kpoint, ionode_id, world_comm )
  CALL mp_bcast( filplot, ionode_id, world_comm )
  CALL mp_bcast( lsign, ionode_id, world_comm )
  CALL mp_bcast( ef_0, ionode_id, world_comm )
  !
  !
  !
  IF (plot_num /= 7) THEN
     call errore('postproc', 'wrong plot_num', plot_num)
  elseif (spin_component(1) < 0 .or. spin_component(1) > 3) then
     CALL errore('postproc', 'wrong spin_component', 1)
  ENDIF
  !
  !   Now allocate space for pwscf variables, read and check them.
  !
     CALL read_file ( )
     IF (nproc_pool /= nproc_pool_file .and. .not. twfcollect)  &
        CALL errore('postproc', &
        'pw.x run with a different number of procs/pools. Use wf_collect=.true.',1)
     CALL openfil_pp ( )
  !
  !
  ! Set default values for emin, emax, degauss_ldos
  ! Done here because ef, degauss must be read from file
  IF (emin > emax) CALL errore('postproc','emin > emax',0)
  ! transforming all back to Ry units
  emin = emin / rytoev
  emax = emax / rytoev
  delta_e = delta_e / rytoev
  degauss_ldos = degauss_ldos / rytoev
  ef_0 = ef_0 / rytoev

  ! Number of output files depends on input
  nplots = 1

  ALLOCATE( plot_files(nplots) )
  plot_files(1) = filplot

  ! 
  kpoint(1) = 1
  kpoint(2) = nkstot
  kband(1) = 1
  kband(2) = nbnd
#if defined(__MPI)
  write(*, '("nr1x = ", i4, ", nr2x = ", i4, ", nr3x = ", i4)') &
        dfftp%nr1x, dfftp%nr2x, dfftp%nr3x
  allocate(psi(dfftp%nr1x *  dfftp%nr2x *  dfftp%nr3x))
  allocate(eband_r(dfftp%nr1x *  dfftp%nr2x *  dfftp%nr3x))
#else
  write(*, '("nnr = ", i4)') dfftp%nnr
  allocate(psi(dfftp%nnr))
  allocate(eband_r(dfftp%nnr))
#endif
  write(*, '("kpoint(1) = ", i4, "  kpoint(2) = ", i4)') kpoint(1), kpoint(2)
  write(*, '("kband(1) = ", i4, "   kband(2) = ", i4)') kband(1), kband(2)
  eband_r = 0.0_DP
  eband_tot = 0.0_DP
  ! Plot multiple KS orbitals in one go
    DO ikpt=kpoint(1), kpoint(2)
      DO ibnd=kband(1), kband(2)
        DO ispin=spin_component(1), spin_component(2)
          !WRITE( *, '("before punch_plot_psi time is ", F10.1)') get_clock('extract_psi')
          CALL punch_plot_psi (psi, plot_num, sample_bias, z, dz, &
            emin, emax, ikpt, ibnd, ispin, lsign)
          !WRITE( *, '("after punch_plot_psi time is ", F10.1)') get_clock('extract_psi')
#if defined(__MPI)
          psi = psi * omega / (dfftp%nr1x *  dfftp%nr2x *  dfftp%nr3x)
#else
          psi = psi * omega / (dfftp%nnr)
#endif
          write(*, '("     processor index ", i3, "   sum(|psi(r)|^2) = ", f18.6)') &
                mpime, sum(psi(:))  ! should be 1
          eband_r = eband_r + wg(ibnd, ikpt) * (et(ibnd, ikpt)-ef_0) * psi
          eband_tot = eband_tot + wg(ibnd, ikpt) * (et(ibnd, ikpt)-ef_0)
          WRITE( stdout, 9000 ) get_clock( 'extract_psi' )
        ENDDO
      ENDDO
    ENDDO

   eband_r = eband_r * rytoev
   write(*, '("eband_tot = ", f18.10, "   eV")') eband_tot*rytoev
   write(*, '("sum(eband_r) = ", f18.10, "   eV")') sum(eband_r(:))
   write(*, '("save real space resolved mae to ", a, ", unit is meV")') &
         TRIM(plot_files(1))
   eband_r = eband_r * 1000  ! unit is meV
   call punch_plot_save(TRIM(plot_files(1)), plot_num, eband_r)

  deallocate(psi)
  deallocate(eband_r)

9000 FORMAT( '     total cpu time spent up to now is ',F10.1,' secs' )

  CALL stop_clock('extract_psi')
  !
END SUBROUTINE extract_psi

SUBROUTINE punch_plot_psi (psi, plot_num, sample_bias, z, dz, &
     emin, emax, kpoint, kband, spin_component, lsign)
  !-----------------------------------------------------------------------
  !
  !     This subroutine writes on output several quantities
  !     in a real space 3D mesh for subsequent processing or plotting
  !     The integer variable plot_num is used to choose the output quantity
  !     See file Doc/INPUT_PP.* for a description of plotted quantities
  !
  !     The output quantity is written (formatted) on file filplot.
  !
  USE kinds,            ONLY : DP
  USE constants,        ONLY : rytoev
  USE cell_base,        ONLY : at, bg, omega, alat, celldm, ibrav
  USE ions_base,        ONLY : nat, ntyp => nsp, ityp, tau, zv, atm
  USE run_info,         ONLY : title
  USE extfield,         ONLY : tefield, dipfield
  USE fft_base,         ONLY : dfftp
  USE scatter_mod,      ONLY : gather_grid
  USE fft_interfaces,   ONLY : fwfft, invfft
  USE gvect,            ONLY : gcutm
  USE gvecs,            ONLY : dual
  USE klist,            ONLY : nks, nkstot, xk
  USE lsda_mod,         ONLY : nspin, current_spin
  USE ener,             ONLY : ehart
  USE io_global,        ONLY : stdout, ionode
  USE scf,              ONLY : rho, vltot, v
  USE wvfct,            ONLY : nbnd, wg
  USE gvecw,            ONLY : ecutwfc
  USE noncollin_module, ONLY : noncolin
  USE paw_postproc,     ONLY : PAW_make_ae_charge

  IMPLICIT NONE
  REAL(DP), INTENT(OUT) :: psi(:)
  INTEGER, INTENT(IN) :: plot_num, kpoint, kband, spin_component
  LOGICAL, INTENT(IN) :: lsign
  REAL(DP), INTENT(IN) :: sample_bias, z, dz, &
      emin, emax
  REAL(DP) :: dummy, charge
  INTEGER :: is, ipol, istates
#if defined(__MPI)
  ! auxiliary vector (parallel case)
  REAL(DP), ALLOCATABLE :: raux1 (:)
#endif
  ! auxiliary vector
  REAL(DP), ALLOCATABLE :: raux (:), raux2(:,:)
  REAL(DP), EXTERNAL :: get_clock


#if defined(__MPI)
  ALLOCATE (raux1(  dfftp%nr1x *  dfftp%nr2x *  dfftp%nr3x))
#endif

  !WRITE( stdout, '(/5x,"Calling punch_plot, plot_num = ",i3)') plot_num
      WRITE( stdout, '(/5x,"Plotting k_point = ",i3,"  band =", i3  )') &
                                                   kpoint, kband
  IF (noncolin .and. spin_component /= 0 ) &
     WRITE( stdout, '(/5x,"Plotting spin magnetization ipol = ",i3)') &
                                                          spin_component
  !
  ALLOCATE (raux(dfftp%nnr))
  !
  !     Here we decide which quantity to plot
  !

     WRITE (title, '("k_point ",i4,", band ",i4)') kpoint ,kband

     IF (noncolin) THEN
        IF (spin_component==0) THEN
  WRITE( *, '("before local_dos time is ", F10.1)') get_clock('extract_psi')
           CALL local_dos_psi ( kpoint, kband, spin_component, raux)
  WRITE( *, '("after local_dos time is ", F10.1)') get_clock('extract_psi')
        ELSE
           CALL local_dos_mag (spin_component, kpoint, kband, raux)
        ENDIF
     ELSE
        CALL errore('postproc','noncolin = .false. not implemented',0)
     ENDIF


#if defined(__MPI)
  IF (.not. (plot_num == 5 ) ) CALL gather_grid (dfftp, raux, raux1)
  IF ( ionode ) &
     psi = raux1
  DEALLOCATE (raux1)
#else
  psi = raux
#endif

  DEALLOCATE (raux)
  RETURN
END SUBROUTINE punch_plot_psi

SUBROUTINE punch_plot_save (filplot, plot_num, raux)
  !-----------------------------------------------------------------------
  !
  !     This subroutine writes on output several quantities
  !     in a real space 3D mesh for subsequent processing or plotting
  !     The integer variable plot_num is used to choose the output quantity
  !     See file Doc/INPUT_PP.* for a description of plotted quantities
  !
  !     The output quantity is written (formatted) on file filplot.
  !
  USE kinds,            ONLY : DP
  USE constants,        ONLY : rytoev
  USE cell_base,        ONLY : at, bg, omega, alat, celldm, ibrav
  USE ions_base,        ONLY : nat, ntyp => nsp, ityp, tau, zv, atm
  USE run_info,         ONLY : title
  USE fft_base,         ONLY : dfftp
  USE fft_interfaces,   ONLY : fwfft, invfft
  USE gvect,            ONLY : gcutm
  USE gvecs,            ONLY : dual
  USE klist,            ONLY : nks, nkstot, xk
  USE lsda_mod,         ONLY : nspin, current_spin
  USE io_global,        ONLY : stdout, ionode
  USE gvecw,            ONLY : ecutwfc

  IMPLICIT NONE
  CHARACTER(len=*), INTENT(IN) :: filplot
  ! auxiliary vector
  REAL(DP), INTENT(IN) :: raux (:)
  INTEGER, INTENT(IN) :: plot_num


  IF (filplot == ' ') RETURN
  !
  !
  !     Here we decide which quantity to plot
  !

#if defined(__MPI)
  IF ( ionode ) &
#endif
  CALL plot_io (filplot, title,  dfftp%nr1x,  dfftp%nr2x,  dfftp%nr3x,  &
        dfftp%nr1,  dfftp%nr2,  dfftp%nr3, nat, ntyp, ibrav, celldm, at,&
        gcutm, dual, ecutwfc, plot_num, atm, ityp, zv, tau, raux, + 1)

  RETURN
END SUBROUTINE punch_plot_save

SUBROUTINE local_dos_psi ( kpoint, kband, spin_component, dos)
  !--------------------------------------------------------------------
  !
  !     iflag=0: calculates |psi|^2 for band "kband" at point "kpoint"
  !     spin_component: for iflag=3 and LSDA calculations only
  !                     0 for up+down dos,  1 for up dos, 2 for down dos
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : omega
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE ener,                 ONLY : ef
  USE fft_base,             ONLY : dffts, dfftp
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE gvect,                ONLY : nl, ngm, g
  USE gvecs,                ONLY : nls, nlsm, doublegrid
  USE klist,                ONLY : lgauss, degauss, ngauss, nks, wk, xk, &
                                   nkstot, ngk, igk_k
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE scf,                  ONLY : rho
  USE symme,                ONLY : sym_rho, sym_rho_init, sym_rho_deallocate
  USE uspp,                 ONLY : nkb, vkb, becsum, nhtol, nhtoj, indv
  USE uspp_param,           ONLY : upf, nh, nhm
  USE wavefunctions_module, ONLY : evc, psic, psic_nc
  USE wvfct,                ONLY : nbnd, npwx, et
  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : noncolin, npol
  USE spin_orb,             ONLY : lspinorb, fcoef
  USE io_files,             ONLY : iunwfc, nwordwfc
  USE mp_global,            ONLY : me_pool, nproc_pool, my_pool_id, npool
  USE mp,                   ONLY : mp_bcast, mp_sum
  USE mp_global,            ONLY : inter_pool_comm, intra_pool_comm
  USE becmod,               ONLY : calbec
  USE control_flags,        ONLY : tqr
  USE realus,               ONLY : addusdens_r
  IMPLICIT NONE
  !
  ! input variables
  !
  INTEGER, INTENT(in) ::  kpoint, kband, spin_component
  !
  REAL(DP), INTENT(out) :: dos (dfftp%nnr)
  !
  !    local variables
  !
  ! counters for US PPs
  INTEGER :: npw, ikb, jkb, ijkb0, ih, jh, kh, na, ijh, np
  ! counters
  INTEGER :: ir, is, ig, ibnd, ik, irm, isup, isdw, ipol, kkb, is1, is2

  REAL(DP) :: w, w1, modulus
  REAL(DP), ALLOCATABLE :: maxmod(:)
  COMPLEX(DP), ALLOCATABLE ::  becp_nc(:,:,:), be1(:,:), be2(:,:)
  INTEGER :: who_calculate, iproc
  COMPLEX(DP) :: phase
  LOGICAL :: i_am_the_pool
  INTEGER :: which_pool, kpoint_pool
  REAL(DP), EXTERNAL :: get_clock
  REAL(DP) :: wg_this
  !
  ! input checks
  !
  IF ( .not. noncolin ) CALL errore('local_dos', 'not implemented', 1)
  IF (gamma_only) CALL errore('local_dos','not available',1)
  !
  IF ( ( kband < 1 .or. kband > nbnd ) ) &
       CALL errore ('local_dos', 'wrong band specified', 1)
  IF ( ( kpoint < 1 .or. kpoint > nkstot ) ) &
       CALL errore ('local_dos', 'wrong kpoint specified', 1)
  !
        ALLOCATE (becp_nc(nkb,npol,nbnd))
        IF (lspinorb) THEN
          ALLOCATE(be1(nhm,2))
          ALLOCATE(be2(nhm,2))
        ENDIF

  rho%of_r(:,:) = 0.d0
  dos(:) = 0.d0
  becsum(:,:,:) = 0.d0
  !
  !   calculate the correct weights
  !
  wg_this  = 0.d0

  IF ( npool > 1 ) THEN
     CALL xk_pool( kpoint, nkstot, kpoint_pool,  which_pool )
     IF ( kpoint_pool < 1 .or. kpoint_pool > nks ) &
        CALL errore('local_dos','problems with xk_pool',1)
     i_am_the_pool=(my_pool_id==which_pool)
  ELSE
     i_am_the_pool=.true.
     kpoint_pool=kpoint
  ENDIF

  IF (i_am_the_pool) wg_this = 1.d0
  !
  !     here we sum for each k point the contribution
  !     of the wavefunctions to the density of states
  !
  ik = kpoint_pool
     IF (i_am_the_pool) THEN
        IF (lsda) current_spin = isk (ik)
        CALL davcio (evc, 2*nwordwfc, iunwfc, ik, - 1)
        npw = ngk(ik)
        CALL init_us_2 (npw, igk_k(1,ik), xk (1, ik), vkb)

           CALL calbec ( npw, vkb, evc, becp_nc )

WRITE( *, '("before ik ibnd loop time is ", F10.1)') get_clock('extract_psi')

     !
     !     here we compute the density of states
     !
        ibnd = kband
         ! Neglect summands with relative weights below machine epsilon
         IF ( wg_this > epsilon(0.0_DP) ) THEN

                 psic_nc = (0.d0,0.d0)
                 DO ig = 1, npw
                    psic_nc(nls(igk_k(ig,ik)),1)=evc(ig     ,ibnd)
                    psic_nc(nls(igk_k(ig,ik)),2)=evc(ig+npwx,ibnd)
                 ENDDO
                 DO ipol=1,npol
                    CALL invfft ('Wave', psic_nc(:,ipol), dffts)
                 ENDDO

              w1 = wg_this / omega

                 DO ipol=1,npol
                    DO ir=1,dffts%nnr
                       rho%of_r(ir,current_spin)=rho%of_r(ir,current_spin)+&
                          w1*(dble(psic_nc(ir,ipol))**2+ &
                             aimag(psic_nc(ir,ipol))**2)
                    ENDDO
                 ENDDO

        !
        !    If we have a US pseudopotential we compute here the becsum term
        !
              w1 = wg_this
              ijkb0 = 0
              DO np = 1, ntyp
                IF (upf(np)%tvanp  ) THEN
                  DO na = 1, nat
                    IF (ityp (na) == np) THEN

                        IF (upf(np)%has_so) THEN
                          be1=(0.d0,0.d0)
                          be2=(0.d0,0.d0)
                          DO ih = 1, nh(np)
                            ikb = ijkb0 + ih
                            DO kh = 1, nh(np)
                              IF ((nhtol(kh,np)==nhtol(ih,np)).and. &
                                  (nhtoj(kh,np)==nhtoj(ih,np)).and. &
                                  (indv(kh,np)==indv(ih,np))) THEN
                                 kkb=ijkb0 + kh
                                 DO is1=1,2
                                   DO is2=1,2
                                     be1(ih,is1)=be1(ih,is1)+ &
                                           fcoef(ih,kh,is1,is2,np)* &
                                           becp_nc(kkb,is2,ibnd)
                                     be2(ih,is1)=be2(ih,is1)+ &
                                           fcoef(kh,ih,is2,is1,np)* &
                                        conjg(becp_nc(kkb,is2,ibnd))
                                   ENDDO
                                 ENDDO
                              ENDIF
                            ENDDO
                          ENDDO
                        ENDIF
                        ijh = 1
                        DO ih = 1, nh (np)
                          ikb = ijkb0 + ih
                          IF (upf(np)%has_so) THEN
                            becsum(ijh,na,1)=becsum(ijh,na,1)+ w1*    &
                               (be1(ih,1)*be2(ih,1)+be1(ih,2)*be2(ih,2))
                          ELSE
                            becsum(ijh,na,1) = becsum(ijh,na,1)+  &
                             w1*(conjg(becp_nc(ikb,1,ibnd))*      &
                                       becp_nc(ikb,1,ibnd)+       &
                                 conjg(becp_nc(ikb,2,ibnd))*      &
                                       becp_nc(ikb,2,ibnd))
                          ENDIF
                          ijh = ijh + 1
                          DO jh = ih + 1, nh (np)
                            jkb = ijkb0 + jh
                            IF (upf(np)%has_so) THEN
                              becsum(ijh,na,1)=becsum(ijh,na,1) &
                                 + w1*((be1(jh,1)*be2(ih,1)+   &
                                        be1(jh,2)*be2(ih,2))+  &
                                       (be1(ih,1)*be2(jh,1)+   &
                                        be1(ih,2)*be2(jh,2)) )
                            ELSE
                              becsum(ijh,na,1)= becsum(ijh,na,1)+ &
                                   w1*2.d0*dble(conjg(becp_nc(ikb,1,ibnd)) &
                                     *becp_nc(jkb,1,ibnd) + &
                                conjg(becp_nc(ikb,2,ibnd)) &
                                     *becp_nc(jkb,2,ibnd) )
                            ENDIF
                            ijh = ijh + 1
                          ENDDO
                        ENDDO

                      ijkb0 = ijkb0 + nh (np)
                    ENDIF
                  ENDDO
                ELSE
                  DO na = 1, nat
                    IF (ityp (na) == np) ijkb0 = ijkb0 + nh (np)
                  ENDDO
                ENDIF
              ENDDO
           ENDIF
        ! loop over bands
    ENDIF
  ! loop over k-points



        IF (lspinorb) THEN
           DEALLOCATE(be1)
           DEALLOCATE(be2)
        ENDIF
        DEALLOCATE(becp_nc)

  IF (doublegrid) THEN

       CALL interpolate(rho%of_r, rho%of_r, 1)

  ENDIF
  !
  !    Here we add the US contribution to the charge
  !
  IF ( tqr ) THEN
   CALL addusdens_r(rho%of_r(:,:),.false.)
  ELSE
  !
  CALL addusdens(rho%of_r(:,:))
  !
  ENDIF
  !
  IF (nspin == 1 .or. nspin==4) THEN
     is = 1
     dos(:) = rho%of_r (:, is)
  ELSE

        isup = 1
        isdw = 2
        dos(:) = rho%of_r (:, isup) + rho%of_r (:, isdw)

  ENDIF

#if defined(__MPI)
  !WRITE( *, '("before local_dos:mp_sum time is ", F10.1)') get_clock('extract_psi')
  CALL mp_sum( dos, inter_pool_comm )
  !WRITE( *, '("after local_dos:mp_sum time is ", F10.1)') get_clock('extract_psi')
#endif

  RETURN

END SUBROUTINE local_dos_psi

END MODULE projmae_module
!
!-----------------------------------------------------------------------
PROGRAM projmae_real
  !-----------------------------------------------------------------------
  !
  !    Program for data analysis and plotting. The two basic steps are:
  !    1) read the output file produced by pw.x, extract and calculate
  !       the desired quantity (rho, V, ...)
  !    2) write the desired quantity to file in a suitable format for
  !       various types of plotting and various plotting programs
  !    The two steps can be performed independently. Intermediate data
  !    can be saved to file in step 1 and read from file in step 2.
  !
  !    DESCRIPTION of the INPUT : see file Doc/INPUT_PP.*
  !
  USE io_global,  ONLY : ionode
  USE mp_global,  ONLY : mp_startup
  USE environment,ONLY : environment_start, environment_end
  USE chdens_module, ONLY : chdens
  USE projmae_module, ONLY : extract_psi

  !
  IMPLICIT NONE
  !
  !CHARACTER(len=256) :: filplot
  CHARACTER(len=256), DIMENSION(:), ALLOCATABLE :: plot_files
  INTEGER :: plot_num
  !
  ! initialise environment
  !
#if defined(__MPI)
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'POST-PROC' )
  !
  IF ( ionode )  CALL input_from_file ( )
  !
  CALL extract_psi (plot_files, plot_num)
  !
  !CALL chdens (plot_files, plot_num)
  !
  CALL environment_end ( 'POST-PROC' )
  !
  CALL stop_pp()
  !
END PROGRAM projmae_real
