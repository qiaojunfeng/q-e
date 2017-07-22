!
! Copyright (C) 2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
PROGRAM plotproj_k
!
!  This small program is used to select the band eigenvalues whose
!  wavefunctions projected on atomic wavefunctions have projections larger
!  than a given threshold. It requires two input files. The first is a
!  file with the band eigenvalues, written in the output of pw.x.
!  The input file with the bands has the following format:
!  nbnd, nks     ! number of bands, number of k points
!  --- blank line
!  kvector coordinates
!  --- blank line
!  bands eigenvalues
!  ...
!  --- blank line
!  kvector coordinates
!  --- blank line
!  bands eigenvalues
!  ...
!
!  The second file is written by the projwfc.x program with the option
!  lsym=.false.
!
!  The input of this program is:
!  filename     ! name of the file with the band eigenvalues
!  filename1    ! name of the file with the projections
!  fileout      ! name of the output file where the bands are written
!  threshold    ! see below
!  ncri         ! number of criterions for selecting the bands
!  for each criterion
!  first_atomic_wfc, last_atomic_wfc   ! the band is selected if the
!                                        sum of the projections on
!                                        the atomic wavefunctions between
!                                        first_atomic_wfc and
!                                        last_atomic_wfc is larger than
!                                        threshold. The sum is done on
!                                        all criterions.
!
  USE constants, ONLY: eps4
  USE spin_orb,   ONLY: lspinorb
  USE projections
  USE ions_base, ONLY : ityp, atm

  IMPLICIT NONE
  !INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  REAL(DP), ALLOCATABLE :: e_par(:,:), e_per(:,:), k(:,:)
  INTEGER :: nks = 0, nbnd = 0, ios, n, i, ibnd, na, idum, nat, &
       natomwfc, nwfc, ntyp
  !LOGICAL, ALLOCATABLE :: toplot(:,:)
  CHARACTER(len=256) :: filename, filename1
  REAL(DP) :: threshold
  REAL(DP), ALLOCATABLE :: proj_par(:,:,:), proj_per(:,:,:)
  !INTEGER, ALLOCATABLE :: first_atomic_wfc(:), last_atomic_wfc(:)
  INTEGER, ALLOCATABLE :: plot_kpts(:)
  INTEGER :: n_plot_kpts, ipltk, tmp, strind, j
  character(len=128) :: tmpstr
  CHARACTER (len=3)  :: str_label(1:2)=(/'par','per'/)
  REAL(DP) :: e_thr, e_diff
  logical, allocatable :: wfctoplot(:)
  integer :: kplot_window(2)
  INTEGER :: stdin = 5, stdout = 6, stderr = 6
  real(DP) :: emin, emax
  ! min energy for having difference between par & per
  real(dp) :: emin_hasdiff = 1, tmpcmp, psum
  real(dp), allocatable :: proj1(:)
  integer, allocatable :: idx(:)
  integer :: lmax_wfc
  !real(dp) :: compute_mj

  ! only plot nearest 5 kpts
  !ALLOCATE (toplot(nbnd,5))

  do tmp = 1,2
    WRITE(stdout,'(a3, " bands.out file : ")', advance="NO") trim(str_label(tmp))
    CALL get_file ( filename )
    ! READ (stdin,'(a)') filename

    OPEN(UNIT=1,FILE=filename,FORM='formatted',status='old',iostat=ios)
    IF (ios/=0) STOP 'Error opening band file '

    read(1, '(a)', err=20, iostat=ios) tmpstr
    READ(tmpstr(13:16), *, err=20) nbnd
    READ(tmpstr(23:28), *, err=20) nks

    IF (nks <= 0 .or. nbnd <= 0 ) THEN
       STOP 'Error reading file header'
    ELSE
       PRINT '("Reading ",i4," bands at ",i6," k-points")', nbnd, nks
    ENDIF

    if (tmp == 1) then
      ALLOCATE (e_par(nbnd, nks))
      ALLOCATE (k(3,nks))
    else
      allocate (e_per(nbnd, nks))
    endif

    DO n=1,nks
       READ(1, '(10x,3f10.6)', ERR=20, IOSTAT=ios) (k(i,n), i=1,3)
       if (tmp == 1) then
         READ(1, '(10e17.9)', END=20, ERR=20) (e_par(i,n),i=1,nbnd)
       else
         READ(1, '(10e17.9)', END=20, ERR=20) (e_per(i,n),i=1,nbnd)
       endif
    ENDDO

20 IF (ios/=0) STOP "problem reading files"
    CLOSE(UNIT=1)

    WRITE(stdout,'(a3, " filproj : ")', advance="NO") trim(str_label(tmp))
    CALL get_file ( filename1 )
    ! READ (stdin,'(a)') filename1
    OPEN(UNIT=1, FILE=filename1, FORM='formatted', STATUS='old', IOSTAT=ios)
    IF (ios/=0) STOP 'Error opening projection file '
    READ(1, *, ERR=20, IOSTAT=ios)
    READ (1, '(8i8)', ERR=20, IOSTAT=ios) idum, idum, idum, idum, idum, &
         idum, nat, ntyp
    DO i=1,2+nat+ntyp
       READ(1, *, ERR=20, IOSTAT=ios)
    ENDDO
    READ (1, '(3i8)',ERR=20, IOSTAT=ios) natomwfc, nks, nbnd
    READ (1, *, ERR=20, IOSTAT=ios)

    if (tmp == 1) then
      ALLOCATE( proj_par(natomwfc,nbnd,nks) )
    else
      ALLOCATE( proj_per(natomwfc,nbnd,nks) )
    endif

    DO nwfc = 1, natomwfc
       READ(1, *, ERR=20, IOSTAT=ios)
       DO n=1,nks
          DO ibnd=1,nbnd
             if ( tmp == 1) then
               READ(1, '(2i8,f20.10)', ERR=20, IOSTAT=ios) idum,idum,proj_par(nwfc,ibnd,n)
             else
               READ(1, '(2i8,f20.10)', ERR=20, IOSTAT=ios) idum,idum,proj_per(nwfc,ibnd,n)
             endif
          ENDDO
       ENDDO
    ENDDO
    CLOSE(1)

  enddo
! end reading par per files

  WRITE(*,'("output file > ")', advance="NO")
  READ(5,'(a)', END=25, ERR=25)  filename

  IF (filename == ' ' ) THEN
     PRINT '("skipping ...")'
     GOTO 25
  ENDIF

30 IF (ios/=0) STOP "problem reading files"
  WRITE(*,'("threshold for selecting wfc > ")', advance="NO")
  READ(5, *, ERR=30, IOSTAT=ios) threshold

  WRITE(*,'("threshold for band energy diff > ")', advance="NO")
  READ(5, *, ERR=30, IOSTAT=ios) e_thr

  WRITE(*,'("number of kpts to be ploted > ")', advance="NO")
  READ(5, *, ERR=30, IOSTAT=ios) n_plot_kpts
  ALLOCATE (plot_kpts(n_plot_kpts))

  do i = 1, n_plot_kpts
    WRITE(*,'(i6, " kpt number  > ")', advance="NO") i
    READ(5, *, ERR=30, IOSTAT=ios) plot_kpts(i)
  enddo

  write(stdout, '("ploting emin, emax > ")', advance='NO')
  read(stdin, *, ERR=30) emin, emax
  write(stdout, '()')

  ! meaning of each states, please refer to kpdos output
!  !
!  ! fill structure nlmchi
!  !
!  CALL fill_nlmchi ( natomwfc, nwfc, lmax_wfc )
!  !
!  !
!  ! write on the standard output file
!  !
!  WRITE( stdout,'(/5x,"Atomic states used for projection")')
!  WRITE( stdout,'( 5x,"(read from pseudopotential files):"/)')
!  IF (lspinorb) THEN
!     DO nwfc = 1, natomwfc
!        WRITE(stdout,1000) &
!          nwfc, nlmchi(nwfc)%na, atm(ityp(nlmchi(nwfc)%na)), &
!          nlmchi(nwfc)%n, nlmchi(nwfc)%jj, nlmchi(nwfc)%l,   &
!          compute_mj(nlmchi(nwfc)%jj,nlmchi(nwfc)%l,nlmchi(nwfc)%m)
!     ENDDO
!1000    FORMAT (5x,"state #",i4,": atom ",i3," (",a3,"), wfc ",i2, &
!                " (j=",f3.1," l=",i1," m_j=",f4.1,")")
!  ELSE
!     DO nwfc = 1, natomwfc
!        WRITE(stdout,1500) &
!          nwfc, nlmchi(nwfc)%na, atm(ityp(nlmchi(nwfc)%na)), &
!          nlmchi(nwfc)%n, nlmchi(nwfc)%l, nlmchi(nwfc)%m, &
!          0.5d0-int(nlmchi(nwfc)%ind/(2*nlmchi(nwfc)%l+2))
!     ENDDO
!1500    FORMAT (5x,"state #",i4,": atom ",i3," (",a3,"), wfc ",i2, &
!                " (l=",i1," m=",i2," s_z=",f4.1,")")
!  ENDIF
!  !

  allocate( wfctoplot(natomwfc) )

do ipltk = 1, n_plot_kpts

  ! get kplot_window
  if (plot_kpts(ipltk) < 3) then
    kplot_window(1) = 1
    kplot_window(2) = 5
  elseif ( plot_kpts(ipltk) > (nks-2) ) then
    kplot_window(1) = nks-4
    kplot_window(2) = nks
  else
    kplot_window(1) = plot_kpts(ipltk)-2
    kplot_window(2) = plot_kpts(ipltk)+2
  endif

  ! find which wfc to be ploted
  do ibnd = 1, nbnd
    if (e_par(ibnd, plot_kpts(ipltk)) < emin) cycle
    e_diff = e_par(ibnd, plot_kpts(ipltk)) - e_per(ibnd, plot_kpts(ipltk))
    if (abs(e_diff) < e_thr ) cycle
    write(stdout, '()')
    write(*, '("ibnd ", i4, ", e_par-e_per = ", f10.6, " > e_thr, ", "e_par = ", f10.6)') &
       ibnd, e_diff, e_par(ibnd,plot_kpts(ipltk))
    !
    ! sort projections by magnitude, in decreasing order
    allocate( proj1(natomwfc), idx(natomwfc) )
    do tmp = 1, 2
      !
      DO nwfc = 1, natomwfc
         idx (nwfc) = 0
         if ( tmp == 1) then
           proj1 (nwfc) = - proj_par (nwfc, ibnd, plot_kpts(ipltk))
         else
           proj1 (nwfc) = - proj_per (nwfc, ibnd, plot_kpts(ipltk))
         endif
      ENDDO
      CALL hpsort_eps (natomwfc, proj1, idx, eps4)
      !
      !  only projections that are larger than 0.001 are written
      !
      DO nwfc = 1, natomwfc
        proj1 (nwfc) = - proj1(nwfc)
        !IF ( abs (proj1(nwfc)) < 0.001d0 ) GOTO 40
        IF ( abs (proj1(nwfc)) < threshold ) GOTO 40
      ENDDO
      nwfc = natomwfc + 1
  40  nwfc = nwfc -1
      !
      ! fancy (?!?) formatting
      !
      if ( tmp == 1) then
        write(stdout, '("  par:")')
      else
        write(stdout, '("  per:")')
      endif
      WRITE( stdout, '(5x,"psi = ",5(f5.3,"*[#",i4,"]+"))') &
        (proj1 (i), idx(i), i = 1, min(5,nwfc))
      DO j = 1, (nwfc-1)/5
        WRITE( stdout, '(10x,"+",5(f5.3,"*[#",i4,"]+"))') &
             (proj1 (i), idx(i), i = 5*j+1, min(5*(j+1),nwfc))
      ENDDO
      if ( tmp == 1) then
        psum = SUM ( proj_par(1:natomwfc, ibnd, plot_kpts(ipltk)) )
      else
        psum = SUM ( proj_par(1:natomwfc, ibnd, plot_kpts(ipltk)) )
      endif
      WRITE( stdout, '(4x,"|psi|^2 = ",f5.3)') psum
      !
    enddo ! tmp
    DEALLOCATE (idx, proj1)
    !
    if (e_diff > 0) then
      tmpcmp = e_per(ibnd, plot_kpts(ipltk))
    else
      tmpcmp = e_par(ibnd, plot_kpts(ipltk))
    endif
    if (emin < tmpcmp) &
      emin_hasdiff = min(emin_hasdiff, tmpcmp)

    do tmp = 1, 2
      wfctoplot = .false.
      DO nwfc = 1, natomwfc
        if ( tmp == 1 ) then
          if (proj_par(nwfc,ibnd,plot_kpts(ipltk)) > threshold) wfctoplot(nwfc) = .true.
        else
          if (proj_per(nwfc,ibnd,plot_kpts(ipltk)) > threshold) wfctoplot(nwfc) = .true.
        endif
      enddo

      DO nwfc = 1, natomwfc
        if (.not. wfctoplot(nwfc)) cycle

        tmpstr = trim(filename)//"_ikpt"
        strind = len_trim(tmpstr) + 1
        if (plot_kpts(ipltk) < 10) then
          write(tmpstr(strind:strind), '(i1)') plot_kpts(ipltk)
          strind = strind + 1
        elseif (plot_kpts(ipltk) < 100) then
          write(tmpstr(strind:strind+1), '(i2)') plot_kpts(ipltk)
          strind = strind + 2
        else
          write(tmpstr(strind:strind+2), '(i3)') plot_kpts(ipltk)
          strind = strind + 3
        endif

        tmpstr = trim(tmpstr)//"_"//trim(str_label(tmp))//"_ibnd"
        strind = len_trim(tmpstr) + 1
        if (ibnd < 10) then
          write(tmpstr(strind:strind), '(i1)') ibnd
          strind = strind + 1
        elseif (ibnd < 100) then
          write(tmpstr(strind:strind+1), '(i2)') ibnd
          strind = strind + 2
        else
          write(tmpstr(strind:strind+2), '(i3)') ibnd
          strind = strind + 3
        endif

        tmpstr = trim(tmpstr)//'_nwfc'
        strind = strind + 5
        if (nwfc < 10) then
          write(tmpstr(strind:strind), '(i1)') nwfc
          strind = strind + 1
        elseif (ibnd < 100) then
          write(tmpstr(strind:strind+1), '(i2)') nwfc
          strind = strind + 2
        elseif (ibnd < 1000) then
          write(tmpstr(strind:strind+2), '(i3)') nwfc
          strind = strind + 3
        else
          write(tmpstr(strind:strind+2), '(i4)') nwfc
          strind = strind + 4
        endif

        OPEN (UNIT=2,FILE=trim(tmpstr),FORM='formatted',STATUS='unknown',IOSTAT=ios)
        IF (ios/=0) STOP "Error opening output file "

        DO i=1,nbnd
          DO n=kplot_window(1),kplot_window(2)
            if (tmp == 1) then
              if ((e_par(i,n) > emax) .or. (e_par(i,n) < emin)) cycle
              WRITE (2,'(i6,"  ", e17.9,"  ", f6.5)') n, e_par(i,n), proj_par(nwfc,i,n)
            else
              if ((e_per(i,n) > emax) .or. (e_per(i,n) < emin)) cycle
              WRITE (2,'(i6,"  ", e17.9,"  ", f6.5)') n, e_per(i,n), proj_par(nwfc,i,n)
            endif
          ENDDO
          write(2,'()')
        ENDDO

        CLOSE (UNIT = 2)
      ENDDO ! nwfc
    enddo ! tmp
  enddo ! ibnd

  write(stdout, '("ikpt ", i6, "  emin_hasdiff = ", f17.9)') &
         plot_kpts(ipltk), emin_hasdiff
enddo ! ipltk

25 CONTINUE

END PROGRAM plotproj_k

! copied from projwfc.f90
!FUNCTION compute_mj(j,l,m)
!   !-----------------------------------------------------------------------
!   USE kinds, ONLY: DP
!   IMPLICIT NONE
!   !
!   REAL(DP) :: compute_mj, j
!   INTEGER  :: l, m
!
!   IF (abs(j-l-0.5d0)<1.d-4) THEN
!       compute_mj=m+0.5d0
!   ELSEIF (abs(j-l+0.5d0)<1.d-4) THEN
!      compute_mj=m-0.5d0
!   ELSE
!      CALL errore('compute_mj','l and j not compatible',1)
!   ENDIF
!
!   RETURN
!END FUNCTION compute_mj
