!----------------------------------------------------------------------
!
! This file contains the common blocks and the data declaration needed
! for the routines.
!
! Input variables that have been previously declared in common_c.h have to be
! re-declared here, in a consistant block. 
!
! Zdenek Johan, Winter 1991.  (Fortran 90)
!----------------------------------------------------------------------
module phcommonvars
    use iso_c_binding
    IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
    !
    !.... parameters  IF YOU CHANGE THESE YOU HAVE TO CHANGE THEM IN
    !                  common_c.h ALSO
    !
    parameter     ( MAXBLK = 50000, MAXTS = 100)
    parameter     ( MAXSH = 32, NSD = 3 )
    parameter     ( maxtask = 256 ) ! this used to be in auxmpi.h
    !
    !  The five types of region topology are  1= Tet, 2=Hex, 3= Wedge (tri-start),
    !                                         4= Wedge (quad-first) 5=pyramid
    !
    !  The two types of face topology are  1= tri, 2=quad
    !
    parameter     ( MAXTOP = 6, MAXSURF=199, MAXREGIONS=255 )
  
    ! the common block nomodule holds all the things which have been removed
    ! from different modules

    ! the next two used to be in auxmpi.h
    !----------------------------------------------------------
    integer INEWCOMM
    common /newcom/ INEWCOMM
    bind(C, name="newcom") :: /newcom/
    !----------------------------------------------------------

    !----------------------------------------------------------
    integer sevsegtype(maxtask,15)
    common /newtyp/ sevsegtype
    !----------------------------------------------------------

    !----------------------------------------------------------
    real*8            bcttimescale,ValueListResist(0:MAXSURF), &
        rhovw,thicknessvw, evw, rnuvw, rshearconstantvw, betai, &
        ValueListWallE(0:MAXREGIONS), &
        ValueListWallh(0:MAXREGIONS), &
        tissSuppStiffCoeff, tissSuppDampCoeff, &
        tissSuppRingStiffCoeff, &
        tissSuppRingDampCoeff, &
        stateFilterCoeff, &
        rescontrol, ResCriteria, heartparam(0:15), stabflux_coeff
    integer           icardio, itvn, ipvsq, &
        incp, numINCPSrfs, nsrflistINCP(0:MAXSURF),incpfile, &
        numResistSrfs, nsrflistResist(0:MAXSURF), &
        numImpSrfs, nsrflistImp(0:MAXSURF),impfile, &
        numRCRSrfs, nsrflistRCR(0:MAXSURF),ircrfile, &
        numTRCRSrfs,nsrflistTRCR(0:MAXSURF),itrcrfile, &
        numCORSrfs, nsrflistCOR(0:MAXSURF),icorfile, &
        numVisFluxSrfs, nsrflistVisFlux(0:MAXSURF), &
        numCalcSrfs, nsrflistCalc(0:MAXSURF), &
        numDirCalcSrfs,nsrflistDirCalc(0:MAXSURF), &
        Lagrange,numLagrangeSrfs,nsrflistLagrange(0:MAXSURF),iLagfile, &
        MinNumIter, &
        ideformwall, iwallmassfactor, iwallstiffactor, nProps, &
        iUseSWB, &
        iUseSWBthickonly, &
        numWallRegions, nsrflistWallRegions(0:MAXREGIONS), &
        nWallETagID, nWallhTagID, &
        iwalldamp, iwallsupp, &
        indsurf, &
        iringdamp, iringsupp, &
        imeasdist, idistancenudge, &
        iinitialprestress, iupdateprestress, &
        iuseBET, numBETFields, iestimator, iheart, &
        numControlledCoronarySrfs, indicesOfCoronarySurfaces(0:MAXSURF), &
        numNetlistLPNSrfs, indicesOfNetlistSurfaces(0:MAXSURF), &
        numLoopClosingCircuits, &
        inputHRandSP, geombcHasObservationFields, &
        geombcHasNodeTags, pureZeroDSimulation, &
        num3DConnectedComponents, surfacesOfEachConnectedComponent(0:MAXSURF), &
        hasMasterPythonControlScript
    common /nomodule/ bcttimescale,ValueListResist, &
        rhovw,thicknessvw, evw, rnuvw, rshearconstantvw, betai, &
        ValueListWallE, &
        ValueListWallh, &
        tissSuppStiffCoeff, tissSuppDampCoeff, &
        tissSuppRingStiffCoeff, &
        tissSuppRingDampCoeff, &
        stateFilterCoeff, &
        rescontrol,ResCriteria, heartparam, &
        stabflux_coeff, &
        icardio, itvn, ipvsq, &
        incp, numINCPSrfs, nsrflistINCP,incpfile, &
        numResistSrfs, nsrflistResist, &
        numImpSrfs, nsrflistImp,impfile, &
        numRCRSrfs, nsrflistRCR,ircrfile, &
        numTRCRSrfs,nsrflistTRCR,itrcrfile, &
        numCORSrfs, nsrflistCOR,icorfile, &
        numVisFluxSrfs, nsrflistVisFlux, &
        numCalcSrfs, nsrflistCalc, &
        numDirCalcSrfs,nsrflistDirCalc, &
        Lagrange,numLagrangeSrfs,nsrflistLagrange,iLagfile, &
        MinNumIter, &
        ideformwall, iwallmassfactor, iwallstiffactor, nProps, &
        iUseSWB, &
        iUseSWBthickonly, &
        numWallRegions, nsrflistWallRegions, &
        nWallETagID, nWallhTagID, &
        iwalldamp, iwallsupp, &
        indsurf, &
        iringdamp, iringsupp, &
        imeasdist, idistancenudge, &
        iinitialprestress, iupdateprestress, &
        iuseBET, numBETFields, iestimator, iheart, &
        numControlledCoronarySrfs, indicesOfCoronarySurfaces, &
        numNetlistLPNSrfs, indicesOfNetlistSurfaces, &
        numLoopClosingCircuits, &
        inputHRandSP, geombcHasObservationFields, &
        geombcHasNodeTags, pureZeroDSimulation, &
        num3DConnectedComponents, surfacesOfEachConnectedComponent, &
        hasMasterPythonControlScript
    bind(C, name="nomodule") :: /nomodule/
    !----------------------------------------------------------

    !----------------------------------------------------------
    integer              numGRCRSrfs,nsrflistGRCR(0:MAXSURF),igrcrfile
    common /grcrbccom/   numGRCRSrfs,nsrflistGRCR,igrcrfile
    bind(C, name="grcrbccom") :: /grcrbccom/
    !----------------------------------------------------------

    !----------------------------------------------------------
    integer seqsize, stepseq(100)
    common /sequence/ seqsize, stepseq
    bind(C, name="sequence") :: /sequence/
    !----------------------------------------------------------

    !----------------------------------------------------------
    integer :: master != 0
    integer numpe, myrank
    common /workfc/ master, numpe, myrank
    bind(C, name="workfc") :: /workfc/
    !----------------------------------------------------------

    !----------------------------------------------------------
    integer, target ::  maxfront
    integer, target ::  nlwork
    common /fronts/ maxfront, nlwork
    bind(C, name="fronts") :: /fronts/
    !----------------------------------------------------------

    !----------------------------------------------------------
    integer numper, nshgt, nshg0
    common /newdim/ numper, nshgt, nshg0
    bind(C, name="newdim") :: /newdim/
    !----------------------------------------------------------
 
    !----------------------------------------------------------
    real*8 birth, death, comtim
    common /timer4/ birth, death, comtim
    bind(C, name="timer4") :: /timer4/
    !----------------------------------------------------------

    !----------------------------------------------------------
    real*8 ttim(100)
    common /extrat/ ttim
    bind(C, name="extrat") :: /extrat/
    !----------------------------------------------------------

    !----------------------------------------------------------
!    real*8             zoutSponge, radSponge, zinSponge, &
!        grthOSponge,grthISponge,betamax
!    integer            spongecontinuity, spongemomentum1, spongemomentum2, &
!        spongeenergy, spongemomentum3
!    common /spongevar/ zoutSponge, radSponge, zinSponge, &
!        grthOSponge,grthISponge,betamax, &
!        spongecontinuity, spongemomentum1, spongemomentum2, &
!        spongeenergy, spongemomentum3
!    bind(C, name="spongevar") :: /spongevar/
    !----------------------------------------------------------

    !----------------------------------------------------------
!    real*8           eles,ylimit(3,9), rmutarget, pzero,  wtavei, &
!        dtavei, dke,  fwr1, flump
!    integer          ierrcalc, ihessian, itwmod, ngaussf,idim, &
!        nlist, nintf(MAXTOP), sonfathvar
!    common /turbvar/ eles,ylimit, rmutarget, pzero,  wtavei, &
!        dtavei, dke,  fwr1, flump, &
!        ierrcalc, ihessian, itwmod, ngaussf,idim, &
!        nlist, nintf
!    bind(C, name="turbvar") :: /turbvar/
    integer           ierrcalc, ihessian
    common /turbvar/ ierrcalc, ihessian
    bind(C, name="turbvar") :: /turbvar/

    !----------------------------------------------------------

    !----------------------------------------------------------
    integer           ierrsmooth !iRANS, iLES, isubmod, ifproj, i2filt, &
        !modlstats, idis, nohomog, ierrsmooth
    common /turbvari/ ierrsmooth !iRANS, iLES, isubmod, ifproj, i2filt, &
        !modlstats, idis, nohomog, ierrsmooth
    bind(C, name="turbvari") :: /turbvari/
    !----------------------------------------------------------

    !----------------------------------------------------------
!    integer          irscale, intpres
!    real*8           plandist, &
!        thetag, ds, tolerence, radcyl, rbltin, rvscal
!    common /spebcvr/ irscale, intpres, &
!        plandist, &
!        thetag, ds, tolerence, radcyl, rbltin, rvscal
!    bind(C, name="spebcvr") :: /spebcvr/
    !----------------------------------------------------------

    !----------------------------------------------------------
    real*8         scdiff(5), tdecay
    integer        nsclr, isclr, nsolt, nosource
    integer        consrv_sclr_conv_vel
    common /sclrs/ scdiff, tdecay, &
        nsclr, isclr, nsolt, nosource, &
        consrv_sclr_conv_vel
    bind(C, name="sclrs") :: /sclrs/
    !----------------------------------------------------------
    !
    !.... common blocks
    !
    parameter (MAXQPT = 125)
    !
    !.... common blocks for hierarchic basis functions
    !
    !----------------------------------------------------------
    real*8          Qpt (MAXTOP,4,MAXQPT),  Qwt (MAXTOP,MAXQPT), &
        Qptb(MAXTOP,4,MAXQPT),  Qwtb(MAXTOP,MAXQPT)
    integer         nint(MAXTOP),           nintb(MAXTOP), &
        ngauss,                 ngaussb,   intp, &
        maxnint
    common /intpt/  Qpt                  ,  Qwt                , &
        Qptb                 ,  Qwtb               , &
        nint,                   nintb, &
        ngauss,                 ngaussb,   intp, &
        maxnint
    !----------------------------------------------------------
  
    ! nsrflist is a binary switch that tells us if a given srfID should be
    ! included in the consistent flux calculation.  It starts from zero
    ! since we need to be able to handle/ignore surfaces with no srfID attached
    !
    ! flxID(numfluxes,nIDs+1)
    ! numfluxes = area, mass, fx, fy, fz, heat, scalar_flux_{1,2,3,4}
    ! nIDs currently set to MAXSURF, each surface has its own
    !

    !----------------------------------------------------------
    real*8          flxID(10,0:MAXSURF), Force(3),HFlux, &
        flxIDsclr(4,MAXSURF)
    integer         nsrflist(0:MAXSURF), isrfIM
    common /aerfrc/ flxID              , Force   ,HFlux, &
        flxIDsclr          , &
        nsrflist           , isrfIM
    bind(C, name="aerfrc") :: /aerfrc/
    !----------------------------------------------------------

    !----------------------------------------------------------
!    real*8 a(100000)
!    common /astore/ a
!    bind(C, name="astore") :: /astore/
    !----------------------------------------------------------

    !----------------------------------------------------------
    integer         lcblk(10,MAXBLK+1),      lcblkb(10,MAXBLK+1)
    common /blkdat/ lcblk             ,      lcblkb
    !----------------------------------------------------------

    !----------------------------------------------------------
    integer :: mnodeb(9,8,3)
    common /mbndnod/ mnodeb
    bind(C, name="mbndnod") :: /mbndnod/
    !----------------------------------------------------------

    !----------------------------------------------------------
    integer, target ::  numnp
    integer, target ::  numel
    integer, target ::  numelb
    integer, target ::  numpbc
    integer, target ::  nen
    integer, target ::  nwnp
    integer, target ::  nfaces
    integer, target ::  numflx
    integer, target ::  ndof
    integer, target ::  iALE
    integer, target ::  icoord
    integer, target ::  navier
    integer, target ::  irs
    integer, target ::  iexec
    integer, target ::  necho
    integer, target ::  ichem
    integer, target ::  iRK
    integer, target ::  nedof
    integer, target ::  nshg
    integer, target ::  nnz
    integer, target ::  istop
    integer, target ::  nflow
    integer, target ::  nnz_tot
    integer, target ::  idtn
    integer, target ::  nshguniq
    integer, target ::  icurrentblk
    common /conpar/ numnp,  numel,  numelb, numpbc, nen,     nwnp, &
        nfaces, numflx, ndof,   iALE,    icoord,  navier, &
        irs,    iexec,  necho,  ichem,  iRK,     nedof, &
        nshg,   nnz,    istop,  nflow,  nnz_tot, idtn, &
        nshguniq, icurrentblk
    bind(C, name="conpar") :: /conpar/
    !----------------------------------------------------------

    !----------------------------------------------------------
    real*8           epsilon_ls, epsilon_lsd, dtlset
    integer          iLSet, &
        ivconstraint, iExpLSSclr1, iExpLSSclr2
    common /levlset/ epsilon_ls, epsilon_lsd, dtlset, &
        iLSet, &
        ivconstraint, iExpLSSclr1, iExpLSSclr2
    bind (C, name="levlset") :: /levlset/
    !----------------------------------------------------------

    !----------------------------------------------------------
    integer, target :: nshape
    integer, target :: nshapeb
    integer, target :: maxshb
    integer, target :: nshl
    integer, target :: nshlb
    integer, target :: nfath
    integer, target :: ntopsh
    integer, target :: nsonmax
    common /shpdat/ nshape, nshapeb, maxshb, &
        nshl, nshlb,nfath,  ntopsh,  nsonmax
    bind (C, name="shpdat") :: /shpdat/
    !----------------------------------------------------------

    !----------------------------------------------------------
    integer         mshp,   mshgl,  mwght,  mshpb,  mshglb, mwghtb, &
        mmut,   mrhot,  mxst
    common /datpnt/ mshp,   mshgl,  mwght,  mshpb,  mshglb, mwghtb, &
        mmut,   mrhot,  mxst
    !----------------------------------------------------------

    !----------------------------------------------------------
    integer :: mcsyst
    integer :: melCat
    integer :: nenCat(8,3)
    integer :: nfaCat(8,3)
    common /melmcat/ mcsyst, melCat, nenCat, nfaCat
    bind (C, name="melmcat") :: /melmcat/
    !----------------------------------------------------------
  
    !----------------------------------------------------------
    integer, target :: lelCat
    integer, target :: lcsyst
    integer, target :: iorder
    integer, target :: nenb
    integer, target :: nelblk
    integer, target :: nelblb
    integer, target :: ndofl
    integer, target :: nsymdl
    integer, target :: nenl
    integer, target :: nfacel
    integer, target :: nenbl
    integer, target :: intind
    integer, target :: mattyp
    common /elmpar/ lelCat, lcsyst, iorder, nenb, &
        nelblk, nelblb, ndofl,  nsymdl, nenl,   nfacel, &
        nenbl,  intind, mattyp
    bind (C, name="elmpar") :: /elmpar/
    !----------------------------------------------------------

    !----------------------------------------------------------
    common /errpar/ numerr
    integer numerr
    !----------------------------------------------------------

    !----------------------------------------------------------
    real*8          E3nsd
    integer         I3nsd,  nsymdf, ndofBC, ndiBCB, ndBCB, &
        Jactyp, jump,   ires,   iprec,  iprev,  ibound, &
        idiff,  lhs,    itau,   ipord,  ipred,  lstres, &
        iepstm
    real*8          dtsfct, taucfct
    integer         ibksiz, iabc, isurf, &
        idflx
    real*8          Bo
    integer         EntropyPressure
    common /genpar/ E3nsd,  &
        I3nsd,  nsymdf, ndofBC, ndiBCB, ndBCB, &
        Jactyp, jump,   ires,   iprec,  iprev,  ibound, &
        idiff,  lhs,    itau,   ipord,  ipred,  lstres, &
        iepstm, &
        dtsfct, taucfct, &
        ibksiz, iabc, isurf, &
        idflx,  &
        Bo,     &
        EntropyPressure
    bind (C, name="genpar") :: /genpar/
    !----------------------------------------------------------

    !!----------------------------------------------------------
    character*128   mesh_filename
    real*8          epstol(8), Delt(MAXTS), CFLfl(MAXTS), CFLsl(MAXTS)
    integer         nstep(MAXTS),   niter(MAXTS), &
        impl(MAXTS)
    real*8          rhoinf(MAXTS)
    integer         LHSupd(6),  loctim(MAXTS)
    real*8          deltol(MAXTS,2)
    integer         memLSFlag
    common /inpdat/ epstol   , Delt, CFLfl, CFLsl, &
        nstep,   niter, &
        impl, &
        rhoinf, &
        LHSupd,  loctim, &
        deltol, &
        memLSFlag
    bind (C, name="inpdat") :: /inpdat/
    !----------------------------------------------------------

    !----------------------------------------------------------
    integer         intg(2,MAXTS),  intpt(3),       intptb(3)
    common /intdat/ intg,           intpt,          intptb
    bind (C, name="intdat") :: /intdat/
    !----------------------------------------------------------

    !----------------------------------------------------------
    integer :: indQpt(3,3,4),  numQpt(3,3,4), intmax != 3
    common /mintpar/ indQpt,   numQpt,        intmax
    bind (C, name="mintpar") :: /mintpar/
    !----------------------------------------------------------

    parameter (ivol= 27)
  
    !----------------------------------------------------------
    integer         iin,    igeom,  ipar,   ibndc,  imat,   iecho, &
        iout,   ichmou, irstin, irstou, ihist,  iflux, &
        ierror, itable, iforce, igraph, itime
    common /mio   / iin,    igeom,  ipar,   ibndc,  imat,   iecho, &
        iout,   ichmou, irstin, irstou, ihist,  iflux, &
        ierror, itable, iforce, igraph, itime
  
    bind (C, name="mio") :: /mio/
    !----------------------------------------------------------

    ! /*         common /andres/ fwr1,ngaussf,idim,nlist */

    !----------------------------------------------------------
    !   character :: fin = 'input.dat'
    !   character :: fgeom = 'geombc.dat'
    !   character :: fpar = 'partition.dat'
    !   character :: fbndc = 'bc.dat'
    !   character :: fmat = 'material.dat'
    !   character :: fecho = 'echo.dat'
    !   character :: frstin = 'restart'
    !   character :: frstou = 'restart'
    !   character :: fhist = 'histor.dat'
    !   character :: ferror = 'error.dat'
    !   character :: ftable = 'table.dat'
    !   character :: fforce = 'forces.dat'
    !   character :: fgraph = 'graph.dat'
    !   character :: ftime = 'time.out'
    character*80     fin,    fgeom,  fpar,  fbndc,  fmat,   fecho, &
        frstin, frstou, fhist, ferror, ftable, fforce, &
        fgraph, ftime
    common /mioname/ fin,    fgeom,  fpar,  fbndc,  fmat,   fecho, &
        frstin, frstou, fhist, ferror, ftable, fforce, &
        fgraph, ftime
    bind (C, name="mioname") :: /mioname/
    !----------------------------------------------------------
  
    !----------------------------------------------------------
    real*8          eGMRES
    integer         lGMRES, iKs,    ntotGM
    common /itrpar/ eGMRES, lGMRES, iKs,    ntotGM
    !----------------------------------------------------------

    !----------------------------------------------------------
    integer         mHBrg,  meBrg,  myBrg,  mRcos,  mRsin
    common /itrpnt/ mHBrg,  meBrg,  myBrg,  mRcos,  mRsin
    !----------------------------------------------------------

    !----------------------------------------------------------
    real*8           pr,     Planck, Stefan, Nh,     Rh,     Rgas, &
        gamma,  gamma1, s0,     const,  xN2,    xO2, &
        yN2,    yO2,    Msh(5), cpsh(5), s0sh(5), h0sh(5), &
        Rs(5),  cps(5), cvs(5), h0s(5), Trot(5), sigs(5), &
        Tvib(5), g0s(5), dofs(5), ithm
    common /mmatpar/ pr,     Planck, Stefan, Nh,     Rh,     Rgas, &
        gamma,  gamma1, s0,     const,  xN2,    xO2, &
        yN2,    yO2,    Msh, cpsh, s0sh, h0sh, &
        Rs,  cps, cvs, h0s, Trot, sigs, &
        Tvib, g0s, dofs, ithm
    bind(C, name="mmatpar") :: /mmatpar/
    !----------------------------------------------------------

    !----------------------------------------------------------
    real*8          datmat(3,7,MAXTS)
    integer         matflg(6,MAXTS)
    integer         nummat,  mexist
    common /matdat/ datmat,      &
        matflg, &
        nummat,  mexist
    bind(C, name="matdat") :: /matdat/
    !----------------------------------------------------------

    !----------------------------------------------------------
    real*8          ro,     vel,    temper, press,  entrop
    integer         ntout, &
        ioform, iowflux, iofieldv
    character*80    iotype
    integer         ioybar
    common /outpar/ ro,     vel,    temper, press,  entrop, &
        ntout, &
        ioform, iowflux, iofieldv, &
        iotype, &
        ioybar
    bind(C, name="outpar") :: /outpar/
    !----------------------------------------------------------

    !----------------------------------------------------------
    integer         mbeg,   mend,   mprec
    common /point / mbeg,   mend,   mprec
    bind(C, name="point") :: /point/
    !----------------------------------------------------------

    !----------------------------------------------------------
    real*8          epsM
    integer         iabres
    common /precis/ epsM,   iabres
    !----------------------------------------------------------

    !----------------------------------------------------------
    integer         npro
    common /propar/ npro
    !----------------------------------------------------------

    !----------------------------------------------------------
    real*8 :: resfrt != 0.00000000000000d+0
    common /resdat/ resfrt
    bind(C, name="resdat") :: /resdat/
    !----------------------------------------------------------

    !----------------------------------------------------------
    integer         imap,   ivart,  iDC,    iPcond, Kspace, nGMRES, &
        iconvflow, iconvsclr, idcsclr(2)
    common /solpar/ imap,   ivart,  iDC,    iPcond, Kspace, nGMRES, &
        iconvflow, iconvsclr, idcsclr
    bind(C, name="solpar") :: /solpar/
    !----------------------------------------------------------

    !----------------------------------------------------------
    integer          indsym(5,5)
    common /msympar/ indsym
    bind(C, name="msympar") :: /msympar/
    !----------------------------------------------------------

    !----------------------------------------------------------
    real*8          time,    CFLfld, CFLsld, Dtgl,   Dtmax,  alpha, &
        etol
    integer         lstep,  ifunc,  itseq,  istep,  iter, &
        nitr
    real*8          almi,   alfi,   gami,   flmpl,  flmpr, &
        dtol(2)
    integer         iCFLworst
    common /timdat/ time,    CFLfld, CFLsld, Dtgl,   Dtmax,  alpha, &
        etol,    &
        lstep,  ifunc,  itseq,  istep,  iter, &
        nitr,    &
        almi,   alfi,   gami,   flmpl,  flmpr, &
        dtol, &
        iCFLworst
    bind(C, name="timdat") :: /timdat/
    !----------------------------------------------------------

    !----------------------------------------------------------
    integer         LCtime, ntseq
    common /timpar/ LCtime, ntseq
    bind(C, name="timpar") :: /timpar/
    !----------------------------------------------------------

    !----------------------------------------------------------
    integer         numeqns(100), minIters, maxIters, &
        iprjFlag,     nPrjs,    ipresPrjFlag, nPresPrjs
    real*8          prestol,      statsflow(6), statssclr(6)
    integer         iverbose
    common /incomp/ numeqns, minIters, maxIters, &
        iprjFlag,     nPrjs,    ipresPrjFlag, nPresPrjs, &
        prestol,      statsflow, statssclr, &
        iverbose
    bind(C, name="incomp") :: /incomp/
    !----------------------------------------------------------

    !----------------------------------------------------------
    ! not sure about the alignment on this one... should it be:
    ! character (c_char) :: ccode(13,8)
    ! common /mtimer1/ ccode(13,8)
    ! common_c.h has ccode defined as array of 13 doubles....
    ! ?????? ONKAR
    ! no C code is using this array
    ! the strings are not NULL terminated and thus cannot be read properly from C
    !  character(8) :: ccode(13) = (/ 'Input   ', 'PrProces', 'Rezoning', 'Elm_Form', &
    !                                 'Solver  ', 'Bnd_Flux', 'Output  ', 'Mapping ', &
    !                                 'Gather  ', 'Scatter ', 'Begin   ', 'End     ', &
    !                                 'Back    ' /)
    character(8) :: ccode(13)
    common /mtimer1/ ccode
    bind(C, name="mtimer1") :: /mtimer1/
    !----------------------------------------------------------

    !----------------------------------------------------------
    ! there was a mismatch in types for flops,  gbytes, sbytes with common_c.h ? ONKAR
    real*8 :: flops,  gbytes, sbytes  ! modified from integer to real
    integer :: iclock, icd, &
        icode, icode2, icode3
    common /mtimer2/ flops,  gbytes, sbytes, iclock, icd,    icode, &
        icode2, icode3
    bind(C, name="mtimer2") :: /mtimer2/
    !----------------------------------------------------------

    !----------------------------------------------------------
    real*8          cpu(11),        cpu0(11)
    integer         nacess(11)
    common /timer3/ cpu,        cpu0,       nacess
    !----------------------------------------------------------

    !----------------------------------------------------------
    character*80    title,  ititle
    common /title / title,  ititle
    !----------------------------------------------------------


    character*8     machin
    parameter     ( machin = 'RS/6000 ' )
    parameter     ( machfl = 4 )
  
    parameter ( zero   = 0.0000000000000000000000000000000d0, &
        pt125  = 0.1250000000000000000000000000000d0, &
        pt25   = 0.2500000000000000000000000000000d0, &
        pt33   = 0.3333333333333333333333333333333d0, &
        pt39   = 0.3968502629920498686879264098181d0, &
        pt5    = 0.5000000000000000000000000000000d0, &
        pt57   = 0.5773502691896257645091487805020d0, &
        pt66   = 0.6666666666666666666666666666667d0, &
        pt75   = 0.7500000000000000000000000000000d0, &
        one    = 1.0000000000000000000000000000000d0, &
        sqt2   = 1.4142135623730950488016887242097d0, &
        onept5 = 1.5000000000000000000000000000000d0, &
        two    = 2.0000000000000000000000000000000d0, &
        three  = 3.0000000000000000000000000000000d0, &
        four   = 4.0000000000000000000000000000000d0, &
        five   = 5.0000000000000000000000000000000d0, &
        pi     = 3.1415926535897932384626433832795d0)

 
    interface
  
        subroutine getintpnts(pts, numpts) bind(C, name="getIntPnts")
            use iso_c_binding
            real (c_double), intent(inout), dimension(*) :: pts    ! array written out by reference
            integer (c_int), value, intent(in) :: numpts           ! input integer
        end subroutine getintpnts

        subroutine shp6w(p, par, N, dN) bind(C, name="shp6w")
            use iso_c_binding
            integer (c_int), value, intent(in) :: p
            real (c_double), intent(in) :: par(3)  ! array of 3 doubles
            real (c_double), intent(inout), dimension(*) :: N  ! not used in WedgeShapeAndDrv
            real (c_double), intent(inout), dimension(*) :: dN ! not used in WedgeShapeAndDrv
        end subroutine shp6w

        subroutine shphex(p, par, N, dN) bind(C, name="shphex")
            use iso_c_binding
            integer (c_int), value, intent(in) :: p
            real (c_double), intent(in) :: par(3)  ! array of 3 doubles
            real (c_double), intent(inout), dimension(*) :: N  ! array written out from HexShapeAndDrv
            real (c_double), intent(inout), dimension(*) :: dN ! array written out from HexShapeAndDrv
        end subroutine shphex

        subroutine shppyr(p, par, N, dN) bind(C, name="shppyr")
            use iso_c_binding
            integer (c_int), value, intent(in) :: p
            real (c_double), intent(in) :: par(3)  ! array of 3 doubles
            real (c_double), intent(inout), dimension(*) :: N  ! array written out
            real (c_double), intent(inout), dimension(*) :: dN ! array written out
        end subroutine shppyr

        subroutine shptet(p, par, N, dN) bind(C, name="shptet")
            use iso_c_binding
            integer (c_int), value, intent(in) :: p
            real (c_double), intent(in) :: par(3)  ! input array 2 doubles
            real (c_double), intent(inout), dimension(*) :: N  ! ouput array
            real (c_double), intent(inout), dimension(*) :: dN ! output array
        end subroutine shptet

        subroutine shptri(p, par, N, dN) bind(C, name="shptri")
            use iso_c_binding
            integer (c_int), value, intent(in) :: p
            real (c_double), intent(in) :: par(2)  ! input array 3 doubles
            real (c_double), intent(inout), dimension(*) :: N  ! array written out
            real (c_double), intent(inout), dimension(*) :: dN ! array written out
        end subroutine shptri

        subroutine symhex(p, pt, wt, err) bind(C, name="symhex")
            use iso_c_binding
            integer (c_int), value, intent(in) :: p
            real (c_double), intent(inout), dimension(*) :: pt  ! array written out
            real (c_double), intent(inout), dimension(*) :: wt  ! array written out
            integer (c_int), intent(inout) :: err ! integer written out
        end subroutine symhex
  
        subroutine sympyr(p, pt, wt, err) bind(C, name="sympyr")
            use iso_c_binding
            integer (c_int), value, intent(in) :: p
            real (c_double), intent(inout), dimension(*) :: pt  ! array written out
            real (c_double), intent(inout), dimension(*) :: wt  ! array written out
            integer (c_int), intent(inout) :: err ! integer written out
        end subroutine sympyr
  
        subroutine symquad(p, pt, wt, err) bind(C, name="symquad")
            use iso_c_binding
            integer (c_int), value, intent(in) :: p
            real (c_double), intent(inout), dimension(*) :: pt  ! array written out
            real (c_double), intent(inout), dimension(*) :: wt  ! array written out
            integer (c_int), intent(inout) :: err ! integer written out
        end subroutine symquad
  
        subroutine symquadw(p, pt, wt, err) bind(C, name="symquadw")
            use iso_c_binding
            integer (c_int), value, intent(in) :: p
            real (c_double), intent(inout), dimension(*) :: pt  ! array written out
            real (c_double), intent(inout), dimension(*) :: wt  ! array written out
            integer (c_int), intent(inout) :: err ! integer written out
        end subroutine symquadw
  
        subroutine symtet(p, pt, wt, err) bind(C, name="symtet")
            use iso_c_binding
            integer (c_int), value, intent(in) :: p
            real (c_double), intent(inout), dimension(*) :: pt  ! array written out
            real (c_double), intent(inout), dimension(*) :: wt  ! array written out
            integer (c_int), intent(inout) :: err ! integer written out
        end subroutine symtet
  
        subroutine symtri(p, pt, wt, err) bind(C, name="symtri")
            use iso_c_binding
            integer (c_int), value, intent(in) :: p
            real (c_double), intent(inout), dimension(*) :: pt  ! array written out
            real (c_double), intent(inout), dimension(*) :: wt  ! array written out
            integer (c_int), intent(inout) :: err ! integer written out
        end subroutine symtri
  
        subroutine symtripyr(p, pt, wt, err) bind(C, name="symtripyr")
            use iso_c_binding
            integer (c_int), value, intent(in) :: p
            real (c_double), intent(inout), dimension(*) :: pt  ! array written out
            real (c_double), intent(inout), dimension(*) :: wt  ! array written out
            integer (c_int), intent(inout) :: err ! integer written out
        end subroutine symtripyr
  
        subroutine symwdg(p, pt, wt, err) bind(C, name="symwdg")
            use iso_c_binding
            integer (c_int), value, intent(in) :: p
            real (c_double), intent(inout), dimension(*) :: pt  ! array written out
            real (c_double), intent(inout), dimension(*) :: wt  ! array written out
            integer (c_int), intent(inout) :: err ! integer written out
        end subroutine symwdg

        real(c_double) function tmrc() bind(C, name="tmrc")
            use iso_c_binding, only: c_double
        end function

        ! this one is needed because it is called from both C++ and Fortran
        subroutine elmdist(u, x, xdist, xdnv, df_fem) bind(C, name="elmdist")
            use iso_c_binding
            real (c_double), intent(in), dimension(*) :: u
            real (c_double), intent(in), dimension(*) :: x
            real (c_double), intent(inout), dimension(*) :: xdist
            real (c_double), intent(inout), dimension(*) :: xdnv
            real (c_double), intent(inout), dimension(*) :: df_fem
        end subroutine elmdist

        ! the following functions are implented in PhGlobalArrayTransfer.cxx

!        subroutine phglobalarrayassignpointer (uniqptr, yoldptr, acoldptr, uoldptr, &
!            coordptr, taptr, distptr, df_femptr, &
!            oyptr, oaptr, ouptr, odptr) &
!            bind(C, name='PhGlobalArrayAssignPointer')
!            use iso_c_binding
!            type(c_ptr), value :: uniqptr
!            type(c_ptr), value :: yoldptr
!            type(c_ptr), value :: acoldptr
!            type(c_ptr), value :: uoldptr
!            type(c_ptr), value :: coordptr
!            type(c_ptr), value :: taptr
!            type(c_ptr), value :: distptr
!            type(c_ptr), value :: df_femptr
!            type(c_ptr), value :: oyptr
!            type(c_ptr), value :: oaptr
!            type(c_ptr), value :: ouptr
!            type(c_ptr), value :: odptr
!        end subroutine phglobalarrayassignpointer

        subroutine phglobalblockedarrayassignpointer (npro_in, nshl_in, ien_in) &
            bind(C, name='PhGlobalBlockedArrayAssignPointer')
            use iso_c_binding
            integer(c_int), value :: npro_in
            integer(c_int), value :: nshl_in
            type(c_ptr), value :: ien_in
        end subroutine phglobalblockedarrayassignpointer

!        subroutine phgloballumpedparameterarrayassignpointer (p_ptr, q_ptr, param_ptr, pout_ptr) &
!            bind(C, name='PhGlobalLumpedParameterArrayAssignPointer')
!            use iso_c_binding
!            type(c_ptr), value :: p_ptr
!            type(c_ptr), value :: q_ptr
!            type(c_ptr), value :: param_ptr
!            type(c_ptr), value :: pout_ptr
!        end subroutine phgloballumpedparameterarrayassignpointer

        subroutine PhAssignPointerInt (ptrInt, fieldName) &
            bind(C, name='PhAssignPointerInt')
            use iso_c_binding
            type(c_ptr), value :: ptrInt
            character(kind=c_char), intent(in) :: fieldName(*)
        end subroutine PhAssignPointerInt

        subroutine PhAssignPointerDP (ptrDP, fieldName) &
            bind(C, name='PhAssignPointerDP')
            use iso_c_binding
            type(c_ptr), value :: ptrDP
            character(kind=c_char), intent(in) :: fieldName(*)
        end subroutine PhAssignPointerDP


    end interface

end module
  
subroutine initPhCommonVars
  
    use phcommonvars
    IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

    ! set default value of iestimator to 0, i.e. flowsolver only
    iestimator = int(0);

    master = 0

    mnodeb(:,:,:) = reshape((/   1,  0,  0,   0,  0,  0,   0,  0,  0,  &
        1,  0,  0,   0,  0,  0,   0,  0,  0,  &
        1,  0,  0,   0,  0,  0,   0,  0,  0,  &
        1,  0,  0,   0,  0,  0,   0,  0,  0,  &
        1,  0,  0,   0,  0,  0,   0,  0,  0,  &
        1,  0,  0,   0,  0,  0,   0,  0,  0,  &
        1,  0,  0,   0,  0,  0,   0,  0,  0,  &
        1,  0,  0,   0,  0,  0,   0,  0,  0,  &   ! 1D
        1,  2,  0,   0,  0,  0,   0,  0,  0,  &
        1,  2,  0,   0,  0,  0,   0,  0,  0,  &
        1,  2,  0,   0,  0,  0,   0,  0,  0,  &
        1,  2,  0,   0,  0,  0,   0,  0,  0,  &
        1,  2,  5,   0,  0,  0,   0,  0,  0,  &
        1,  2,  4,   0,  0,  0,   0,  0,  0,  &
        1,  2,  4,   0,  0,  0,   0,  0,  0,  &
        1,  2,  4,   0,  0,  0,   0,  0,  0,  &  ! 2D
        1,  2,  3,   4,  0,  0,   0,  0,  0,  &
        1,  2,  3,   0,  0,  0,   0,  0,  0,  &
        1,  2,  3,   0,  0,  0,   0,  0,  0,  &
        1,  2,  5,   4,  0,  0,   0,  0,  0,  &
        1,  2,  3,   4,  9, 10,  11, 12, 21,  &
        1,  2,  3,   5,  6,  9,   0,  0,  0,  &
        1,  2,  3,   7,  9,  8,   0,  0,  0,  &
        1,  2,  5,   4,  7, 10,  13, 14, 16   /), &
        shape(mnodeb))  ! 3D
                              
    !nsd = 3
                              
    mcsyst = 4
    melCat = 8
  
    nenCat(:,:) = reshape((/  2,  2,  2,  2,    3,  3,  3,  3, &        ! 1D
        4,  3,  3,  4,    9,  6,  6,  9, &        ! 2D
        8,  4,  6,  6,   27, 10, 18, 18    /), &
        shape(nenCat))     ! 3D
    nfaCat(:,:) = reshape((/  2,  2,  2,  2,    2,  2,  2,  2, &       ! 1D
        4,  3,  3,  4,    4,  3,  3,  4, &       ! 2D
        6,  4,  5,  5,    6,  4,  5,  5    /), &     ! 3D
        shape(nfaCat))
                               
    intmax = 3
  
    iin = 10
    igeom = 11
    ipar = 12
    ibndc = 13
    imat = 14
    iecho = 15
    iout = 16
    ichmou = 17
    irstin = 18
    irstou = 19
    ihist = 20
    iflux = 21
    ierror = 22
    itable = 23
    iforce = 24
    igraph = 25
    itime = 26
  
    fin = 'input.dat'
    fgeom = 'geombc.dat'
    fpar = 'partition.dat'
    fbndc = 'bc.dat'
    fmat = 'material.dat'
    fecho = 'echo.dat'
    frstin = 'restart'
    frstou = 'restart'
    fhist = 'histor.dat'
    ferror = 'error.dat'
    ftable = 'table.dat'
    fforce = 'forces.dat'
    fgraph = 'graph.dat'
    ftime = 'time.out'
  
    Planck = 6.62617600000000d-34
    Stefan = 5.66970000000000d-08
    Nh = 6.02204500000000d+23
    Rh = 8.31441000000000d+0
    gamma = 1.40000000000000d+0
    gamma1 = 0.40000000000000d+0
    xN2 = 0.79000000000000d+0
    xO2 = 0.21000000000000d+0
    Msh(:) = reshape((/ 2.80000000000000d-2,  3.20000000000000d-2, &
        3.00000000000000d-2,  1.40000000000000d-2, &
        1.60000000000000d-2 /), &
        shape(Msh))
    h0sh(:) = reshape((/ 0.00000000000000d+0,  0.00000000000000d+0, &
        8.97750000000000d+4,  4.70820000000000d+5, &
        2.46790000000000d+5 /), &
        shape(h0sh))
    Trot(:) = reshape((/ 2.87000000000000d+0,  2.08000000000000d+0, &
        2.45000000000000d+0,  0.00000000000000d+0, &
        0.00000000000000d+0 /), &
        shape(Trot))
    sigs(:) = reshape((/ 2.00000000000000d+0,  2.00000000000000d+0, &
        1.00000000000000d+0,  0.00000000000000d+0, &
        0.00000000000000d+0 /), &
        shape(sigs))
    Tvib(:) = reshape((/ 3.39350000000000d+3,  2.27356000000000d+3, &
        2.73887000000000d+3,  0.00000000000000d+0, &
        0.00000000000000d+0 /), &
        shape(Tvib))
    g0s(:) = reshape((/ 1.00000000000000d+0,  3.00000000000000d+0, &
        4.00000000000000d+0,  4.00000000000000d+0, &
        9.00000000000000d+0 /), &
        shape(g0s))
    dofs(:) = reshape((/ 5.00000000000000d+0,  5.00000000000000d+0, &
        5.00000000000000d+0,  3.00000000000000d+0, &
        3.00000000000000d+0 /), &
        shape(dofs))
  

    mbeg = 1
    mend = 100000
    mprec = 2
  
    resfrt = 0.00000000000000d+0
  
    indsym(:,:) = reshape((/ 1,  2,  4,  7, 11, &
        2,  3,  5,  8, 12, &
        4,  5,  6,  9, 13, &
        7,  8,  9, 10, 14, &
        11, 12, 13, 14, 15   /), &
        shape(indsym))
                         
    ccode(:) = (/ 'Input   ', 'PrProces', 'Rezoning', 'Elm_Form', &
        'Solver  ', 'Bnd_Flux', 'Output  ', 'Mapping ', &
        'Gather  ', 'Scatter ', 'Begin   ', 'End     ', &
        'Back    ' /)
                         
    icd = 11
  
    numerr = 10
              
end subroutine

!
!----------------------------------------------------------------------
!
!.... element pointers
!
! mmat   (MAXBLK)  : pointer to interior element material number
! mmatb  (MAXBLK)  : pointer to boundary element material number
! mien   (MAXBLK)  : pointer to ien array
! mienb  (MAXBLK)  : pointer to ienb array
! miBCB  (MAXBLK)  : pointer to iBCB array
! mDt    (MAXBLK)  : pointer to Dt array
! mDC    (MAXBLK)  : pointer to DC array
! mBCB   (MAXBLK)  : pointer to BCB array
! mstiff (MAXBLK)  : pointer to stiff array
!
!----------------------------------------------------------------------
!
!.... common /aerfrc/   : aerodynamic forces
!
! Force(3)      : components of the aerodynamic forces
! HFlux         : total heat flux
!
!----------------------------------------------------------------------
!
!.... common /astore/   : the dynamic memory allocation area
!
! a(...)        : the blank array used for front-end data storage
!
!----------------------------------------------------------------------
!
!.... common /blkdat/   : blocking data
!
! lcblk  (10,MAXBLK+1) : blocking data for the interior elements
! lcblkb (10,MAXBLK+1) : blocking data for the boundary elements
!
!----------------------------------------------------------------------
!
!.... common /bndnod/   : boundary nodes of boundary elements
!
! mnodeb (9,8,3) : boundary nodes of each element category and dimension
!
!----------------------------------------------------------------------
!
!.... common /conpar/   : input constants
!
! numnp         : number of nodal points
! numel         : number of elements
! numelb        : number of boundary elements
! numpbc        : number of nodes having a boundary condition
! nen           : maximum number of element nodes
! nfaces        : maximum number of element faces
! nsd           : number of space dimensions
! numflx        : number of flux boundary nodes
! ndof          : number of degrees of freedom per node
! iALE          : ALE formulation flag
! icoord        : coordinate system flag
! navier        : Navier-Stokes calculation flag
! irs           : restart option 
! iexec         : execute flag
! necho         : input echo parameter
! ichem         : equilibrium chemistry flag (for outchem.step dump)
! iRK           : Runge-Kutta flag
! nshg          : global number of shape functions (degrees of freedom,
!                 or equations). Computed from the specified p-order,
!                 the number of edges, and the number of faces (in the
!                 entire mesh)
!
!----------------------------------------------------------------------
!
!.... common /datpnt/   : front-end data pointers
!
! mshp          : pointer to shape-functions 
! mshgl         : pointer to local-grad-shape-functions
! mwght         : pointer to quadrature weights
! mshpb         : pointer to shape-functions of boundary elements
! mshglb        : pointer to local-grad-shape-functions of bound. elem.
! mwghtb        : pointer to quadrature weights of bound. elements
! mmut          : pointer to table mu  = mu  (p,T)
! mrhot         : pointer to table rho = rho (p,T)
! mxst          : pointer to table xs  = xs  (p,T)
!
!----------------------------------------------------------------------
!
!.... common /elmcat/   : element category information
!
! mcsyst        : maximum number of element coordinate system
! melCat        : maximum number of element categories
! nenCat (8,3)  : number of nodes for each category and dimension
! nfaCat (8,3)  : number of faces for each category and dimension
!
!----------------------------------------------------------------------
!
!.... common /elmpar/   : element parameters
!
! lelCat        : element category (P1, Q1, P2, Q2, etc.)
! lcsyst        : element coordinate system
! iorder        : element order (=k for Pk and Qk)
! nenb          : number of element nodes per boundary sides
! maxsh         : total number integration points
! maxshb        : total number integration points of boundary elements
! nelblk        : number of element blocks
! nelblb        : number of boundary element blocks
! ndofl         : number of degrees of freedom (for current block)
! nsymdl        : number of d.o.f for symm. storage (for current block)
! nenl          : number of element nodes (for current block)
! nfacel        : number of element faces (for current block)
! nenbl         : number of boundary element nodes
! intind        : integration data index
! nintg         : number of integration points
! mattyp        : material type ( = 0 for fluid; = 1 for solid )
!
!----------------------------------------------------------------------
!
!.... common /errpar/   : element parameters
!
! numerr        : number of error statistics
!
!----------------------------------------------------------------------
!
!.... common /genpar/   : control parameters
!
! E3nsd         : NSD .eq. 3 flag; 0. for 2D, 1. for 3D
! I3nsd         : NSD .eq. 3 flag; 0  for 2D, 1  for 3D
! nsymdf        : number of d.o.f.'s in symm. storage (= ndof*(ndof+1)/2)
! ndofBC        : dimension size of the boundary condition array BC
! ndiBCB        : dimension size of the boundary condition array iBCB
! ndBCB         : dimension size of the boundary condition array BCB
! Jactyp        : Jacobian type flag
! jump          : jump term computation flag
! ires          : residual type computation flag
! iprec         : block-diagonal preconditioner flag
! iprev         : ypl array allocation flag
! ibound        : boundary element flag
! idiff         : diffusive flux vector flag
!                 ( = 0 not used; = 1 global reconstruction )
! itau          : type of tau to be used
!
!----------------------------------------------------------------------
!
!.... common /inpdat/   : time sequence input data
!
! epstol (MAXTS)  : tolerance for GMRES solvers
! Delt   (MAXTS)  : global time step
! CFLfl  (MAXTS)  : CFL number for fluid flow
! CFLsl  (MAXTS)  : CFL number for structural heating
! nstep  (MAXTS)  : number of time steps
! niter  (MAXTS)  : number of iterations per time step
! impl   (MAXTS)  : solver flag
! iturb  (MAXTS)  : turbulence model flag
! rhoinf (MAXTS)  : time integration spectral radius paramter
!                             (0=Gears       1= trapezoidal rule)
! LHSupd (MAXTS)  : LHS/preconditioner update
! loctim (MAXTS)  : local time stepping flag
!
!----------------------------------------------------------------------
!
!.... common /intdat/   : integration data
!
! intg  (2,MAXTS) : integration parameters
! intpt (3)       : integration pointers
! intptb(3)       : integration pointers of boundary elements
!
!----------------------------------------------------------------------
!
!.... common /shpdat/   : hierarchic shape function quadrature data
!
! Qpt  (3,MAXQPT)  : interior element quadrature points (xi,eta,zeta)
! Qwt  (MAXQPT)    : interior element quad. weights
! Qptb (2,MAXQPT)  : boundary element quad. pnts.
! Qwtb (MAXQPT)    : boundary element quad. weights
! nshape           : number of interior element shape functions
! nshapeb          :   "    "  boundary  "        "       "
! ngauss           : number of interior element integration points
! ngaussb          :   "    "  boundary  "        "           "
!----------------------------------------------------------------------
!
!.... common /intpar/   : integration parameters
!
! Qpt   (4,*)   : xi, eta, zeta, weight of quadrature points
! indQpt(3,3,4) : index to quadrature points for a given rule
! numQpt(3,3,4) : number of quadrature points for a given rule
! intmax        : number of allowable spatial integ. points per nsd
!
!----------------------------------------------------------------------
!
!.... common /io    /   : io channels
!
! iin           : input  (main parameters)          [INPUT.DAT]
! igeom         : input  (problem geometry)         [GEOM.DAT]
! ipar          : in/out (spectral mapping)         [PARTITION.DAT]
! ibndc         : input  (problem boundary cond.)   [BC.DAT]
! imat          : input  (element material types)   [MATERIAL.DAT]
! iecho         : output (echo of input)            [ECHO.DAT]
! iout          : output (result output)            [OUTPUT.lstep]
! ichmou        : output (chemistry output)         [OUTCHM.lstep]
! irstin        : input  (input restart)            [RESTAR.INP]
! irstou        : output (output restart)           [RESTAR.OUT]
! ihist         : output (history output)           [HISTOR.DAT]
! iflux         : output (boundary flux)            [FLUX.lstep]
! ierror        : output (error messages)           [ERROR.DAT]
! itable        : input  (equilibrium chemistry)    [TABLE.DAT]
! iforce        : output (aerodynamic forces)       [FORCES.DAT]
! ivol		: output (element volumes)          [VOL.DAT]   
! istat		: output (run status file)          [STATUS.DAT]
!
!----------------------------------------------------------------------
!
!.... common /ioname/   : io file names
!
! fin           : input.dat
! fgeom         : geom.dat
! fpar          : partition.dat
! fbndc         : bc.dat
! fmat          : material.dat
! fecho         : echo.dat
! frstin        : restar.inp
! frstou        : restar.out
! fhist         : histor.dat
! ferror        : error.dat
! ftable        : table.dat
! fforce        : forces.dat
! fvol		: vol.dat
! fstat		: status.dat
!
!----------------------------------------------------------------------
!
!.... common /itrpar/   : Preconditioned GMRES parameters
!
! eGMRES        : finite difference interval
! lGMRES        : number of GMRES cycles
! iKs           : current Krylov vector
! ntotGM        : total number of GMRES iterations
!
!----------------------------------------------------------------------
!
!.... common /itrpnt/   : Preconditioned GMRES array pointers
!
! mHBrg         : pointer to Hessenberg matrix
! meBrg         : pointer to Hessenberg's RHS matrix
! myBrg         : pointer to minimize solution matrix
! mRcos         : pointer to Rotation Cosine of QR algorithm
! mRsin         : pointer to Rotation Sine   of QR algorithm
!
!----------------------------------------------------------------------
!
!.... common /matpar/   : material constants
!
! pr            : Prandtl number
! Planck        : Planck's constant
! Stefan        : Stefan's constant (for radiation)
! Nh            : Avogadro's number
! Rh            : universal gas constant
! Rgas          : specific gas constant
! gamma         : specific heat ratio
! gamma1        : gamma - 1
! s0            : reference specific entropy
! const         : special constant
! xN2           : mole fraction of diatomic nitrogen
! xO2           : mole fraction of diatomic oxygen
! yN2           : mole fraction of diatomic nitrogen
! yO2           : mole fraction of diatomic oxygen
! Msh  (5)      : molar mass of species
! cpsh (5)      : molar heat at constant pressure of species
! s0sh (5)      : molar reference entropy of species
! h0sh (5)      : molar heat of formation of species
! Rs   (5)      : specific gas constant of species
! cps  (5)      : specific heat at constant pressure of species
! cvs  (5)      : specific heat at constant volume of species
! h0s  (5)      : specific heat of formation of species
! Trot (5)      : characteristic rotational temperature of species
! sigs (5)      : symmetry factor of species
! Tvib (5)      : characteristic vibrational temperature of species
! g0s  (5)      : ground degeneracy of electronic energy
! dofs (5)      : degrees of freedom of species
! ithm          : thermodynamic property flag
!
!----------------------------------------------------------------------
!
!.... common /matdat/   : material data
!
! datmat (3,5,2) : material data
! matflg (5,100)   : material type flag
! nummat           : number of materials
! mexist           : flag indicating the presence of MATERIAL.DAT
!
!----------------------------------------------------------------------
!
!.... common /outpar/   : output parameters
!
! ro            : density     rescaling factor for output
! vel           : velocity    rescaling factor for output
! temper        : temperature rescaling factor for output
! press         : pressure    rescaling factor for output
! entrop        : entropy     rescaling factor for output
! ntout         : number of steps between consecutive printouts
! ioform        : output I/O format
!
!----------------------------------------------------------------------
!
!.... common /point /   : dynamic storage pointer management data
!
! mbeg          : pointer to the beginning of the free storage
! mend          : pointer to the end of the storage
! mprec         : precision of the floating point data
!
!----------------------------------------------------------------------
!
!.... common /precis/   : finite difference interval data
!
! epsM          : square root of machine precision
! iabres        : absolute value residual flag
!
!----------------------------------------------------------------------
!
!....common /propar/    : processor related information
!
! npro          : number of virtual processors for the current block
!
!----------------------------------------------------------------------
!
!....common /resdat/    : residual statistics data
!
! resfrt        : first residual of convergence
!
!----------------------------------------------------------------------
!
!.... common /solpar/   : solution parameters
!
! imap          : permutation mapping flag
! ivart         : variational formulation type
! iDC           : DC type
! iPcond        : type of preconditioner
! Kspace        : dimension of Krylov space
! nGMRES        : maximum number of GMRES iterations
!
!----------------------------------------------------------------------
!
!.... common /sympar/   : symmetric storage parameters
!
! indsym (5,5)  : mapping from 2D storage to symmetric one
!
!----------------------------------------------------------------------
!
!.... common /timdat/   : time data
!
! time          : current run time
! CFLfld        : CFL number for fluid flow
! CFLsld        : CFL number for structural heating
! Dtgl          : inverse of global time step
! Dtmax         : maximum delta-time
! alpha         : trapezoidal rule parameter
! etol          : epsilon tolerance for GMRES
! lstep         : current time step
! ifunc         : func. eval. counter (=niter*(lstep-lstep0) + iter)
! itseq         : sequence number
! istep         : step number (reseted at the beginning of the run)
! iter          : iteration number
! nitr          : number of multi-corrector iterations for this sequence
!
!----------------------------------------------------------------------
!
!.... common /timpar/   : time integration parameters
!
! LCtime        : local time stepping flag
! ntseq         : number of time sequences
!
!----------------------------------------------------------------------
!
!.... common /timer1/   : timer parameters
!.... common /timer2/   : timer parameters
!.... common /timer3/   : timer parameters
!
! ccode(13)     : timing entities codes
! flops         : flop counter
! gbytes        : byte counter for gather operation
! sbytes        : byte counter for scatter operation
! iclock        : wall-clock time (in milliseconds)
! icd           : number of timing entities
! icode         : current timer code
! icode2        : last timer code
! icode3        : next-to-last timer code
! cpu(11)       : cpu time of each entity
! cpu0(11)      : initial cpu time of each entity
! nacess(11)    : number of times each entity is accessed
!
!----------------------------------------------------------------------
!
!.... common /title /   : problem title
!
! title         : problem title
! ititle        : problem title (with form feed)
!
!----------------------------------------------------------------------
!
!.... common /avging / : nfath
! 
! nfath         : total number of global fathers over which certain
!                 quantities will be averaged
! 
!----------------------------------------------------------------------
!
!.... parameters        : machine data
!
! machin        : machine type
!                  (set parameter)
! machfl        : single precision floating point lenght in bytes
!                  (set parameter)
!
!----------------------------------------------------------------------
!
!.... parameters        : useful constants
!
! zero          : 0.0
! pt125         : 0.125
! pt25          : 0.25
! pt33          : 0.33 (1/3)
! pt39          : 2^(-4/3)
! pt5           : 0.5
! pt57          : 1/sqrt(3)
! pt66          : 0.66 (2/3)
! pt75          : 0.75
! one           : 1.0
! sqt2          : sqrt(2)
! onept5        : 1.5
! two           : 2.0
! three         : 3.0
! four          : 4.0
! five          : 5.0
! pi            : the magical number :-)
! 
!---------------------------------------------------------------------- 
! 
! Zdenek Johan, Winter 1991.
! 
!----------------------------------------------------------------------
