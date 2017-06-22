!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!  
!
!-----------------------------------------------------------------------
SUBROUTINE non_scf ( )
    !-----------------------------------------------------------------------
    !
    ! ... diagonalization of the KS hamiltonian in the non-scf case
    !
    USE kinds,                ONLY : DP
    USE bp,                   ONLY : lelfield, lberry, lorbm
    USE check_stop,           ONLY : stopped_by_user
    USE control_flags,        ONLY : io_level, conv_elec, lbands
    USE ener,                 ONLY : ef
    USE io_global,            ONLY : stdout, ionode
    USE io_files,             ONLY : iunwfc, nwordwfc, iunefield
    USE buffers,              ONLY : save_buffer
    USE klist,                ONLY : xk, wk, nks, nkstot
    USE lsda_mod,             ONLY : lsda, nspin
    USE wvfct,                ONLY : nbnd, et, npwx
    USE wavefunctions_module, ONLY : evc
    !
    !!!!!!
    USE cell_base,            ONLY : at, bg, alat, omega, tpiba2
    USE ions_base,            ONLY : zv, nat, nsp, ityp, tau, compute_eextfor
    USE vlocal,               ONLY : strf
    USE gvect,                ONLY : ngm, gstart, nl, nlm, g, gg, gcutm
    USE gvecs,                ONLY : doublegrid, ngms
    USE control_flags,        ONLY : mixing_beta, tr2, ethr, niter, nmix, &
        iprint, conv_elec, &
        restart, io_level, do_makov_payne,  &
        gamma_only, iverbosity, textfor,     &
        llondon, scf_must_converge, lxdm, ts_vdw
    USE ener,                 ONLY : etot, hwf_energy, eband, deband, ehart, &
        vtxc, etxc, etxcc, ewld, demet, epaw, &
        elondon, ef_up, ef_dw, exdm
    USE scf,                  ONLY : rho, rho_core, rhog_core, v, vltot, vrs, &
        kedtau, vnew
    USE scf,                  ONLY : scf_type, scf_type_COPY, bcast_scf_type,&
        create_scf_type, destroy_scf_type, &
        open_mix_file, close_mix_file, &
        rho, rho_core, rhog_core, v, vltot, vrs, &
        kedtau, vnew
    USE fft_base,             ONLY : dfftp
    USE noncollin_module,     ONLY : noncolin, magtot_nc, i_cons,  bfield, &
        lambda, report
    USE paw_variables,        ONLY : okpaw, ddd_paw, total_core_energy, only_paw
    USE spin_orb,             ONLY : domag
    USE ldaU,                 ONLY : eth, Hubbard_U, Hubbard_lmax, &
                                     niter_with_fixed_ns, lda_plus_u
    USE extfield,             ONLY : tefield, etotefield, monopole, etotmonofield !TB
    USE xdm_module,           ONLY : energy_xdm
    USE tsvdw_module,         ONLY : EtsvdW
    USE mp,                   ONLY : mp_sum, mp_bcast
    USE mp_bands,             ONLY : intra_bgrp_comm
    USE mp_pools,             ONLY : root_pool, my_pool_id, inter_pool_comm
    USE io_files,             ONLY : iunmix, output_drho, &
                                     iunres, iunefield, seqopn
    USE fcp_variables,        ONLY : lfcpopt, lfcpdyn
    USE plugin_variables,     ONLY : plugin_etot
    USE klist,                ONLY : xk, wk, nelec, ngk, nks, nkstot, lgauss, &
                                     two_fermi_energies, tot_charge
    USE lsda_mod,             ONLY : lsda, nspin, magtot, absmag, isk
    USE paw_onecenter,        ONLY : PAW_potential
    USE paw_symmetry,         ONLY : PAW_symmetrize_ddd
    USE dfunct,               ONLY : newd
    !!!!!!
    IMPLICIT NONE
    !
    ! ... local variables
    !
    INTEGER :: iter, i
    REAL(DP), EXTERNAL :: get_clock
    !
    !!!!!
    REAL(DP) :: &
        dr2,          &! the norm of the diffence between potential
        charge,       &! the total charge
        deband_hwf,   &! deband for the Harris-Weinert-Foulkes functional
        mag           ! local magnetization
    !
    REAL(DP) :: &
      tr2_min,     &! estimated error on energy coming from diagonalization
      descf,       &! correction for variational energy
      en_el=0.0_DP,&! electric field contribution to the total energy
      eext=0.0_DP   ! external forces contribution to the total energy
    !
    ! ... auxiliary variables for calculating and storing temporary copies of
    ! ... the charge density and of the HXC-potential
    !
    type (scf_type) :: rhoin ! used to store rho_in of current/next iteration
    !
    REAL(DP) :: exxen = 0.0d0  ! current estimate of the exchange energy
    !
    ! ... external functions
    !
    REAL(DP), EXTERNAL :: ewald
    REAL(DP) :: etot_cmp_paw(nat,2,2)
    !
    INTEGER :: printout = 2
    !!!!!
    CALL start_clock( 'electrons' )
    iter = 1
    !
    !
    ! ... calculates the ewald contribution to total energy
    !
    ewld = ewald( alat, nat, nsp, ityp, zv, at, bg, tau, &
        omega, g, gg, ngm, gcutm, gstart, gamma_only, strf )
    !
    elondon = 0.d0
    !
    call create_scf_type ( rhoin )
    !
    ! ... deband = - \sum_v <\psi_v | V_h + V_xc |\psi_v> is calculated a
    ! ... first time here using the input density and potential ( to be
    ! ... used to calculate the Harris-Weinert-Foulkes energy )
    !
    deband_hwf = delta_e()
    !
    ! ... save input density to rhoin
    !
    call scf_type_COPY( rho, rhoin )
    !
    WRITE( stdout, 9002 )
    FLUSH( stdout )
    !
    IF ( lelfield) THEN
        !
        CALL c_bands_efield ( iter )
       !
    ELSE
        !
        CALL c_bands_nscf ( )
       !
    END IF
    !
    ! ... check if calculation was stopped in c_bands
    !
    IF ( stopped_by_user ) THEN
        conv_elec=.FALSE.
        RETURN
    END IF
    !
    ! ... xk, wk, isk, et, wg are distributed across pools;
    ! ... the first node has a complete copy of xk, wk, isk,
    ! ... while eigenvalues et and weights wg must be
    ! ... explicitly collected to the first node
    ! ... this is done here for et, in weights () for wg
    !
    CALL poolrecover( et, nbnd, nkstot, nks )
    !
    ! ... calculate weights of Kohn-Sham orbitals (only weights, not Ef,
    ! ... for a "bands" calculation where Ef is read from data file)
    ! ... may be needed in further calculations such as phonon
    !
    IF ( lbands ) THEN
        CALL weights_only  ( )
    ELSE
        CALL weights  ( )
    END IF
    !
    ! ... Note that if you want to use more k-points for the phonon
    ! ... calculation then those needed for self-consistency, you can,
    ! ... by performing a scf with less k-points, followed by a non-scf
    ! ... one with additional k-points, whose weight on input is set to zero
    !
    !
    ! ... the new density is computed here. For PAW:
    ! ... sum_band computes new becsum (stored in uspp modules)
    ! ... and a subtly different copy in rho%bec (scf module)
    !
    CALL sum_band_no_weights()
    !
    ! ... the Harris-Weinert-Foulkes energy is computed here using only
    ! ... quantities obtained from the input density
    !
    hwf_energy = eband + deband_hwf + (etxc - etxcc) + ewld + ehart + demet
    If ( okpaw ) hwf_energy = hwf_energy + epaw
    ! does not consider lda_plus_u
    !
    ! ... calculate total and absolute magnetization
    !
    IF ( lsda .OR. noncolin ) CALL compute_magnetization()
    !
    ! ... eband  = \sum_v \epsilon_v    is calculated by sum_band
    ! ... deband = - \sum_v <\psi_v | V_h + V_xc |\psi_v>
    ! ... eband + deband = \sum_v <\psi_v | T + Vion |\psi_v>
    !
    deband = delta_e()
    !
    ! ... mix_rho mixes several quantities: rho in g-space, tauk (for
    ! ... meta-gga), ns and ns_nc (for lda+u) and becsum (for paw)
    ! ... Results are broadcast from pool 0 to others to prevent trouble
    ! ... on machines unable to yield the same results from the same
    ! ... calculation on same data, performed on different procs
    ! ... The mixing should be done on pool 0 only as well, but inside
    ! ... mix_rho there is a call to rho_ddot that in the PAW case
    ! ... contains a hidden parallelization level on the entire image
    !
    ! IF ( my_pool_id == root_pool )
    !CALL mix_rho ( rho, rhoin, mixing_beta, dr2, tr2_min, iter, nmix, &
    !    iunmix, conv_elec )
    ! we do not need mix_rho in nscf, just copy it
    rhoin = rho
    conv_elec = .TRUE.
    CALL bcast_scf_type ( rhoin, root_pool, inter_pool_comm )
    CALL mp_bcast ( dr2, root_pool, inter_pool_comm )
    CALL mp_bcast ( conv_elec, root_pool, inter_pool_comm )
    ! ... convergence reached:
    ! ... 1) the output HXC-potential is saved in v
    ! ... 2) vnew contains V(out)-V(in) ( used to correct the forces ).
    !
    vnew%of_r(:,:) = v%of_r(:,:)
    CALL v_of_rho( rho,rho_core,rhog_core, &
        ehart, etxc, vtxc, eth, etotefield, charge, v)
    vnew%of_r(:,:) = v%of_r(:,:) - vnew%of_r(:,:)
    !
    IF (okpaw) THEN
        CALL PAW_potential(rho%bec, ddd_paw, epaw,etot_cmp_paw)
        CALL PAW_symmetrize_ddd(ddd_paw)
    ENDIF
    !
    ! ... note that rho is here the output, not mixed, charge density
    ! ... so correction for variational energy is no longer needed
    !
    descf = 0._dp
           !
    plugin_etot = 0.0_dp
    !
    CALL plugin_scf_energy(plugin_etot,rhoin)
    !
    CALL plugin_scf_potential(rhoin,conv_elec,dr2)
    !
    ! ... define the total local potential (external + scf)
    !
    CALL sum_vrs( dfftp%nnr, nspin, vltot, v%of_r, vrs )
    !
    ! ... interpolate the total local potential
    !
    CALL interpolate_vrs( dfftp%nnr, nspin, doublegrid, kedtau, v%kin_r, vrs )
    !
    ! ... in the US case we have to recompute the self-consistent
    ! ... term in the nonlocal potential
    ! ... PAW: newd contains PAW updates of NL coefficients
    !
    CALL newd()
    !
    IF ( lelfield ) en_el =  calc_pol ( )
    !
    IF ( ( MOD(iter,report) == 0 ) .OR. ( report /= 0 .AND. conv_elec ) ) THEN
        !
        IF ( (noncolin .AND. domag) .OR. i_cons==1 .OR. nspin==2) CALL report_mag()
       !
    END IF

    WRITE( stdout, 9000 ) get_clock( 'PWSCF' )
    !
    WRITE( stdout, 9102 )
    !
    ! ... write band eigenvalues (conv_elec is used in print_ks_energies)
    !
    conv_elec = .true.
    CALL print_ks_energies ( )
    !
    etot = eband + ( etxc - etxcc ) + ewld + ehart + deband + demet + descf
    ! for hybrid calculations, add the current estimate of exchange energy
    ! (it will subtracted later if exx_is_active to be replaced with a better estimate)
    etot = etot - exxen
    hwf_energy = hwf_energy - exxen ! [LP]
    !
    IF (okpaw) etot = etot + epaw
    IF ( lda_plus_u ) etot = etot + eth
    !
    IF ( lelfield ) etot = etot + en_el
    ! not sure about the HWF functional in the above case
    IF( textfor ) THEN
        eext = alat*compute_eextfor()
        etot = etot + eext
        hwf_energy = hwf_energy + eext
    END IF
    IF (llondon) THEN
        etot = etot + elondon
        hwf_energy = hwf_energy + elondon
    END IF
    ! calculate the xdm energy contribution with converged density
    if (lxdm .and. conv_elec) then
        exdm = energy_xdm()
        etot = etot + exdm
        hwf_energy = hwf_energy + exdm
    end if
    IF (ts_vdw) THEN
        ! factor 2 converts from Ha to Ry units
        etot = etot + 2.0d0*EtsvdW
        hwf_energy = hwf_energy + 2.0d0*EtsvdW
    END IF
    !
    IF ( tefield ) THEN
        etot = etot + etotefield
        hwf_energy = hwf_energy + etotefield
    END IF
    ! TB monopole energy
    IF ( monopole) THEN
        etot = etot + etotmonofield
        hwf_energy = hwf_energy + etotmonofield
    END IF
    !
    IF ( lfcpopt .or. lfcpdyn ) THEN
        etot = etot + ef * tot_charge
        hwf_energy = hwf_energy + ef * tot_charge
    ENDIF
    !
    ! ... adds possible external contribution from plugins to the energy
    !
    etot = etot + plugin_etot
    !
    CALL print_energies ( printout )
    call destroy_scf_type ( rhoin )
    !
    ! ... save converged wfc if they have not been written previously
    ! ... FIXME: it shouldn't be necessary to do this here
    !
    IF ( nks == 1 .AND. (io_level < 2) .AND. (io_level > -1) ) &
        CALL save_buffer ( evc, nwordwfc, iunwfc, nks )
    !
    ! ... do a Berry phase polarization calculation if required
    !
    IF ( lberry ) CALL c_phase()
    !
    ! ... do an orbital magnetization (Kubo terms) calculation
    !
    IF ( lorbm ) CALL orbm_kubo()
    !
    CALL stop_clock( 'electrons' )
    !
9000 FORMAT(/'     total cpu time spent up to now is ',F10.1,' secs' )
9002 FORMAT(/'     Band Structure Calculation' )
9102 FORMAT(/'     End of band structure calculation' )
!
  !
CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE compute_magnetization()
        !-----------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        INTEGER :: ir
        !
        !
        IF ( lsda ) THEN
            !
            magtot = 0.D0
            absmag = 0.D0
            !
            DO ir = 1, dfftp%nnr
                !
                mag = rho%of_r(ir,1) - rho%of_r(ir,2)
                !
                magtot = magtot + mag
                absmag = absmag + ABS( mag )
               !
            END DO
            !
            magtot = magtot * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
            absmag = absmag * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
            !
            CALL mp_sum( magtot, intra_bgrp_comm )
            CALL mp_sum( absmag, intra_bgrp_comm )
            !
            IF (two_fermi_energies.and.lgauss) bfield(3)=0.5D0*(ef_up-ef_dw)
           !
        ELSE IF ( noncolin ) THEN
            !
            magtot_nc = 0.D0
            absmag    = 0.D0
            !
            DO ir = 1,dfftp%nnr
                !
                mag = SQRT( rho%of_r(ir,2)**2 + &
                    rho%of_r(ir,3)**2 + &
                    rho%of_r(ir,4)**2 )
                !
                DO i = 1, 3
                    !
                    magtot_nc(i) = magtot_nc(i) + rho%of_r(ir,i+1)
                   !
                END DO
                !
                absmag = absmag + ABS( mag )
               !
            END DO
            !
            CALL mp_sum( magtot_nc, intra_bgrp_comm )
            CALL mp_sum( absmag, intra_bgrp_comm )
            !
            DO i = 1, 3
                !
                magtot_nc(i) = magtot_nc(i) * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
               !
            END DO
            !
            absmag = absmag * omega / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
           !
        END IF
        !
        RETURN
      !
    END SUBROUTINE compute_magnetization
    !
    !-----------------------------------------------------------------------
    FUNCTION delta_e()
        !-----------------------------------------------------------------------
        ! ... delta_e = - \int rho%of_r(r)  v%of_r(r)
        !               - \int rho%kin_r(r) v%kin_r(r) [for Meta-GGA]
        !               - \sum rho%ns       v%ns       [for LDA+U]
        !               - \sum becsum       D1_Hxc     [for PAW]
        USE funct,  ONLY : dft_is_meta
        IMPLICIT NONE
        REAL(DP) :: delta_e, delta_e_hub
        !
        delta_e = - SUM( rho%of_r(:,:)*v%of_r(:,:) )
        !
        IF ( dft_is_meta() ) &
            delta_e = delta_e - SUM( rho%kin_r(:,:)*v%kin_r(:,:) )
        !
        delta_e = omega * delta_e / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
        !
        CALL mp_sum( delta_e, intra_bgrp_comm )
        !
        if (lda_plus_u) then
            if (noncolin) then
                delta_e_hub = - SUM (rho%ns_nc(:,:,:,:)*v%ns_nc(:,:,:,:))
                delta_e = delta_e + delta_e_hub
            else
                delta_e_hub = - SUM (rho%ns(:,:,:,:)*v%ns(:,:,:,:))
                if (nspin==1) delta_e_hub = 2.d0 * delta_e_hub
                delta_e = delta_e + delta_e_hub
            endif
        end if
        !
        IF (okpaw) delta_e = delta_e - SUM(ddd_paw(:,:,:)*rho%bec(:,:,:))
        !
        RETURN
      !
    END FUNCTION delta_e
    !
    !-----------------------------------------------------------------------
    FUNCTION delta_escf()
        !-----------------------------------------------------------------------
        !
        ! ... delta_escf = - \int \delta rho%of_r(r)  v%of_r(r)
        !                  - \int \delta rho%kin_r(r) v%kin_r(r) [for Meta-GGA]
        !                  - \sum \delta rho%ns       v%ns       [for LDA+U]
        !                  - \sum \delta becsum       D1         [for PAW]
        ! ... calculates the difference between the Hartree and XC energy
        ! ... at first order in the charge density difference \delta rho(r)
        !
        USE funct,  ONLY : dft_is_meta
        IMPLICIT NONE
        REAL(DP) :: delta_escf, delta_escf_hub
        !
        delta_escf = - SUM( ( rhoin%of_r(:,:)-rho%of_r(:,:) )*v%of_r(:,:) )
        !
        IF ( dft_is_meta() ) &
            delta_escf = delta_escf - &
            SUM( (rhoin%kin_r(:,:)-rho%kin_r(:,:) )*v%kin_r(:,:))
        !
        delta_escf = omega * delta_escf / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
        !
        CALL mp_sum( delta_escf, intra_bgrp_comm )
        !
        if (lda_plus_u) then
            if (noncolin) then
                delta_escf_hub = - SUM((rhoin%ns_nc(:,:,:,:)-rho%ns_nc(:,:,:,:))*v%ns_nc(:,:,:,:))
                delta_escf = delta_escf + delta_escf_hub
            else
                delta_escf_hub = - SUM((rhoin%ns(:,:,:,:)-rho%ns(:,:,:,:))*v%ns(:,:,:,:))
                if (nspin==1) delta_escf_hub = 2.d0 * delta_escf_hub
                delta_escf = delta_escf + delta_escf_hub
            endif
        end if

        IF (okpaw) delta_escf = delta_escf - &
            SUM(ddd_paw(:,:,:)*(rhoin%bec(:,:,:)-rho%bec(:,:,:)))

        RETURN
      !
    END FUNCTION delta_escf
    !
    !-----------------------------------------------------------------------
    FUNCTION calc_pol ( ) RESULT ( en_el )
        !-----------------------------------------------------------------------
        !
        USE kinds,     ONLY : DP
        USE constants, ONLY : pi
        USE bp,        ONLY : lelfield, ion_pol, el_pol, fc_pol, l_el_pol_old, &
            el_pol_acc, el_pol_old, efield, l3dstring, gdir, &
            transform_el, efield_cart
        !
        IMPLICIT NONE
        REAL (DP) :: en_el
        !
        INTEGER :: i, j
        REAL(DP):: sca, el_pol_cart(3),  el_pol_acc_cart(3)
        !
        IF (.not.l3dstring) THEN
            CALL c_phase_field(el_pol(gdir),ion_pol(gdir),fc_pol(gdir),gdir)
            if (.not.l_el_pol_old) then
                l_el_pol_old=.true.
                el_pol_old(gdir)=el_pol(gdir)
                en_el=-efield*(el_pol(gdir)+ion_pol(gdir))
                el_pol_acc(gdir)=0.d0
            else
                sca=(el_pol(gdir)-el_pol_old(gdir))/fc_pol(gdir)
                if(sca < - pi) then
                    el_pol_acc(gdir)=el_pol_acc(gdir)+2.d0*pi*fc_pol(gdir)
                else if(sca > pi) then
                    el_pol_acc(gdir)=el_pol_acc(gdir)-2.d0*pi*fc_pol(gdir)
                endif
                en_el=-efield*(el_pol(gdir)+ion_pol(gdir)+el_pol_acc(gdir))
                el_pol_old=el_pol
            endif
        ELSE
            do i=1,3
                CALL c_phase_field(el_pol(i),ion_pol(i),fc_pol(i),i)
            enddo
            el_pol_cart(:)=0.d0
            do i=1,3
                do j=1,3
                    !el_pol_cart(i)=el_pol_cart(i)+transform_el(j,i)*el_pol(j)
                    el_pol_cart(i)=el_pol_cart(i)+at(i,j)*el_pol(j) / &
                        (sqrt(at(1,j)**2.d0+at(2,j)**2.d0+at(3,j)**2.d0))
                enddo
            enddo

            write(stdout,'( "Electronic Dipole on Cartesian axes" )')
            do i=1,3
                write(stdout,*) i, el_pol_cart(i)
            enddo

            write(stdout,'( "Ionic Dipole on Cartesian axes" )')
            do i=1,3
                write(stdout,*) i, ion_pol(i)
            enddo

            if(.not.l_el_pol_old) then
                l_el_pol_old=.true.
                el_pol_old(:)=el_pol(:)
                en_el=0.d0
                do i=1,3
                    en_el=en_el-efield_cart(i)*(el_pol_cart(i)+ion_pol(i))
                enddo
                el_pol_acc(:)=0.d0
            else
                do i=1,3
                    sca=(el_pol(i)-el_pol_old(i))/fc_pol(i)
                    if(sca < - pi) then
                        el_pol_acc(i)=el_pol_acc(i)+2.d0*pi*fc_pol(i)
                    else if(sca > pi) then
                        el_pol_acc(i)=el_pol_acc(i)-2.d0*pi*fc_pol(i)
                    endif
                enddo
                el_pol_acc_cart(:)=0.d0
                do i=1,3
                    do j=1,3
                        el_pol_acc_cart(i)=el_pol_acc_cart(i)+transform_el(j,i)*el_pol_acc(j)
                    enddo
                enddo
                en_el=0.d0
                do i=1,3
                    en_el=en_el-efield_cart(i)*(el_pol_cart(i)+ion_pol(i)+el_pol_acc_cart(i))
                enddo
                el_pol_old(:)=el_pol(:)
            endif
        ENDIF
      !
    END FUNCTION calc_pol
    !
    !-----------------------------------------------------------------------
    SUBROUTINE print_energies ( printout )
        !-----------------------------------------------------------------------
        !
        USE constants, ONLY : eps8
        INTEGER, INTENT (IN) :: printout
        !

        IF ( printout == 0 ) RETURN
        IF ( ( conv_elec .OR. MOD(iter,iprint) == 0 ) .AND. printout > 1 ) THEN
            !
            IF ( dr2 > eps8 ) THEN
                WRITE( stdout, 9081 ) etot, hwf_energy, dr2
            ELSE
                WRITE( stdout, 9083 ) etot, hwf_energy, dr2
            END IF
            IF ( only_paw ) WRITE( stdout, 9085 ) etot+total_core_energy
            !
            WRITE( stdout, 9100) eband, deband

9100    FORMAT(/'     one-electron contribution = eband + deband:',/,&
            /'     eband                     =',F17.8,' Ry' &
            /'     deband                    =',F17.8,' Ry' )

            WRITE( stdout, 9060 ) &
                ( eband + deband ), ehart, ( etxc - etxcc ), ewld
            !
            IF ( llondon ) WRITE ( stdout , 9074 ) elondon
            IF ( lxdm )    WRITE ( stdout , 9075 ) exdm
            IF ( ts_vdw )  WRITE ( stdout , 9076 ) 2.0d0*EtsvdW
            IF ( textfor)  WRITE ( stdout , 9077 ) eext
            IF ( tefield )            WRITE( stdout, 9061 ) etotefield
            IF ( monopole )           WRITE( stdout, 9062 ) etotmonofield ! TB
            IF ( lda_plus_u )         WRITE( stdout, 9065 ) eth
            IF ( ABS (descf) > eps8 ) WRITE( stdout, 9069 ) descf
            IF ( okpaw ) THEN
                WRITE( stdout, 9067 ) epaw
                ! Detailed printout of PAW energy components, if verbosity is high
                IF(iverbosity>0)THEN
                    WRITE( stdout, 9068) SUM(etot_cmp_paw(:,1,1)), &
                        SUM(etot_cmp_paw(:,1,2)), &
                        SUM(etot_cmp_paw(:,2,1)), &
                        SUM(etot_cmp_paw(:,2,2)), &
                        SUM(etot_cmp_paw(:,1,1))+SUM(etot_cmp_paw(:,1,2))+ehart, &
                        SUM(etot_cmp_paw(:,2,1))+SUM(etot_cmp_paw(:,2,2))+etxc-etxcc
                ENDIF
            ENDIF
            !
            ! ... With Fermi-Dirac population factor, etot is the electronic
            ! ... free energy F = E - TS , demet is the -TS contribution
            !
            IF ( lgauss ) WRITE( stdout, 9070 ) demet
            !
            ! ... With Fictitious charge particle (FCP), etot is the grand
            ! ... potential energy Omega = E - muN, -muN is the potentiostat
            ! ... contribution.
            !
            IF ( lfcpopt .or. lfcpdyn ) WRITE( stdout, 9072 ) ef*tot_charge
           !
        ELSE IF ( conv_elec ) THEN
            !
            IF ( dr2 > eps8 ) THEN
                WRITE( stdout, 9081 ) etot, hwf_energy, dr2
            ELSE
                WRITE( stdout, 9083 ) etot, hwf_energy, dr2
            END IF
           !
        ELSE
            !
            IF ( dr2 > eps8 ) THEN
                WRITE( stdout, 9080 ) etot, hwf_energy, dr2
            ELSE
                WRITE( stdout, 9082 ) etot, hwf_energy, dr2
            END IF
        END IF
        !
        CALL plugin_print_energies()
        !
        IF ( lsda ) WRITE( stdout, 9017 ) magtot, absmag
        !
        IF ( noncolin .AND. domag ) &
            WRITE( stdout, 9018 ) magtot_nc(1:3), absmag
        !
        IF ( i_cons == 3 .OR. i_cons == 4 )  &
            WRITE( stdout, 9071 ) bfield(1), bfield(2), bfield(3)
        IF ( i_cons /= 0 .AND. i_cons < 4 ) &
            WRITE( stdout, 9073 ) lambda
        !
        FLUSH( stdout )
        !
        RETURN
          !
9017    FORMAT(/'     total magnetization       =', F9.2,' Bohr mag/cell', &
            /'     absolute magnetization    =', F9.2,' Bohr mag/cell' )
9018    FORMAT(/'     total magnetization       =',3F9.2,' Bohr mag/cell' &
            &   ,/'     absolute magnetization    =', F9.2,' Bohr mag/cell' )
9060    FORMAT(/'     The total energy is the sum of the following terms:',/,&
            /'     one-electron contribution =',F17.8,' Ry' &
            /'     hartree contribution      =',F17.8,' Ry' &
            /'     xc contribution           =',F17.8,' Ry' &
            /'     ewald contribution        =',F17.8,' Ry' )
9061    FORMAT( '     electric field correction =',F17.8,' Ry' )
9062    FORMAT( '     monopole field correction =',F17.8,' Ry' ) ! TB
9065    FORMAT( '     Hubbard energy            =',F17.8,' Ry' )
9067    FORMAT( '     one-center paw contrib.   =',F17.8,' Ry' )
9068    FORMAT( '      -> PAW hartree energy AE =',F17.8,' Ry' &
            /'      -> PAW hartree energy PS =',F17.8,' Ry' &
            /'      -> PAW xc energy AE      =',F17.8,' Ry' &
            /'      -> PAW xc energy PS      =',F17.8,' Ry' &
            /'      -> total E_H with PAW    =',F17.8,' Ry'&
            /'      -> total E_XC with PAW   =',F17.8,' Ry' )
9069    FORMAT( '     scf correction            =',F17.8,' Ry' )
9070    FORMAT( '     smearing contrib. (-TS)   =',F17.8,' Ry' )
9071    FORMAT( '     Magnetic field            =',3F12.7,' Ry' )
9072    FORMAT( '     pot.stat. contrib. (-muN) =',F17.8,' Ry' )
9073    FORMAT( '     lambda                    =',F11.2,' Ry' )
9074    FORMAT( '     Dispersion Correction     =',F17.8,' Ry' )
9075    FORMAT( '     Dispersion XDM Correction =',F17.8,' Ry' )
9076    FORMAT( '     Dispersion T-S Correction =',F17.8,' Ry' )
9077    FORMAT( '     External forces energy    =',F17.8,' Ry' )
9080    FORMAT(/'     total energy              =',0PF17.8,' Ry' &
            /'     Harris-Foulkes estimate   =',0PF17.8,' Ry' &
            /'     estimated scf accuracy    <',0PF17.8,' Ry' )
9081    FORMAT(/'!    total energy              =',0PF17.8,' Ry' &
            /'     Harris-Foulkes estimate   =',0PF17.8,' Ry' &
            /'     estimated scf accuracy    <',0PF17.8,' Ry' )
9082    FORMAT(/'     total energy              =',0PF17.8,' Ry' &
            /'     Harris-Foulkes estimate   =',0PF17.8,' Ry' &
            /'     estimated scf accuracy    <',1PE17.1,' Ry' )
9083    FORMAT(/'!    total energy              =',0PF17.8,' Ry' &
            /'     Harris-Foulkes estimate   =',0PF17.8,' Ry' &
            /'     estimated scf accuracy    <',1PE17.1,' Ry' )
9085    FORMAT(/'     total all-electron energy =',0PF17.6,' Ry' )

    END SUBROUTINE print_energies
  !
END SUBROUTINE non_scf

