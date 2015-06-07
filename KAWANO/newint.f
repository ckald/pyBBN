      SUBROUTINE check

C-------REMARKS.
C     This is an interface subroutine,
C     a flexible module which allows user to manipulate physical quantities
C     of interest at certain key points during the computer run.
C     Included within this subroutine is a roster of all global variables
C     and their respective COMMON areas.
C     Current instructions accumulate abundances of deuterium, helium-3,
C     helium-4, and lithium-7 for eventual plotting, taking into account
C     the contribution of beryllium-7 to lithium-7 and tritium to helium-3.

      USE commons

C-------COMMON AREAS.
      COMMON /recpr0/ reacpr                       !Reaction parameter values.
      COMMON /recpr/  iform,ii,jj,kk,ll,rev,q9     !Reaction parameter names.
      COMMON /rates/  f,r                            !Reaction rates.
      COMMON /evolp1/ t9,hv,phie,y                   !Evolution parameters.
      COMMON /evolp2/ dt9,dhv,dphie,dydt             !Evolution parameters.
      COMMON /evolp3/ t90,hv0,phie0,y0               !Evolution parameters.
      COMMON /compr0/ cy0,ct0,t9i0,t9f0,ytmin0,inc0  !Default comp parameters.
      COMMON /compr/  cy,ct,t9i,t9f,ytmin,inc        !Computation parameters.
      COMMON /modpr0/ c0,cosmo0,xi0                !Default model parameters.
      COMMON /modpr/  g,tau,xnu,c,cosmo,xi           !Model parameters.
      COMMON /varpr0/ dt0,eta0                     !Default variationl params.
      COMMON /varpr/  dt1,eta1                       !Variational parameters.
      COMMON /time/   t,dt,dlt9dt                    !Time variables.
      COMMON /xtherm/  thm,hubcst                     !Dynamic variables.
      COMMON /endens/ rhone0,rhob0,rhob,rnb          !Energy densities.
      COMMON /lncoef/ a,b,yx                         !Linear eqn coefficients.
      COMMON /nucdat/ am,zm,dm                       !Nuclide data.
      COMMON /xbessel/ bl1,bl2,bl3,bl4,bl5,           !Eval function bl(z).
     |                bm1,bm2,bm3,bm4,bm5,           !Eval function bm(z).
     |                bn1,bn2,bn3,bn4,bn5            !Eval function bn(z).
      COMMON /kays/   bk0,bk1,bk2,bk3,bk4            !Coefficients K.
      COMMON /flags/  ltime,is,ip,it,mbad            !Flags,counters.
      COMMON /check1/  itime                          !Computation location.
      COMMON /outdat/ xout,thmout,t9out,tout,dtout,  !Output data.
     |                etaout,hubout
      COMMON /nupar/  t9mev,tnmev,tnu,cnorm,nu,rhonu !Neutrino parameters.
      COMMON /runopt/ irun,isize,jsize               !Run options.
      COMMON /outopt/ nout,outfile                   !Output option.


C=================DECLARATION DIVISION=====================

C-------REACTION PARAMETER VALUES.
      DOUBLE PRECISION    reacpr(nrec,8)       !Reaction parameters.

C-------REACTION PARAMETER NAMES.
      INTEGER iform(nrec)          !Reaction type code (1-11).
      INTEGER ii(nrec)             !Incoming nuclide type (1-26).
      INTEGER jj(nrec)             !Incoming light nuclide type (1-6).
      INTEGER kk(nrec)             !Outgoing light nuclide type (1-6).
      INTEGER ll(nrec)             !Outgoing nuclide type (1-26).
      DOUBLE PRECISION    rev(nrec)            !Reverse reaction coefficient.
      DOUBLE PRECISION    q9(nrec)             !Energy released in reaction (in 10**9 K).

C-------REACTION RATES.
      DOUBLE PRECISION    f(nrec)              !Forward reaction rate coefficients.
      DOUBLE PRECISION    r(nrec)              !Reverse reaction rate coefficients.

C-------EVOLUTION PARAMETERS.
      DOUBLE PRECISION    t9                   !Temperature of photons (units of 10**9 K).
      DOUBLE PRECISION    hv                   !Defined by hv = M(atomic)n(baryon)/t9**3.
      DOUBLE PRECISION    phie                 !Chemical potential of electron.
      DOUBLE PRECISION    y(nnuc)              !Relative number abundances.

C-------EVOLUTION PARAMETERS (DERIVATIVES).
      DOUBLE PRECISION    dt9                  !Change in temperature.
      DOUBLE PRECISION    dhv                  !Change in hv.
      DOUBLE PRECISION    dphie                !Change in chemical potential.
      DOUBLE PRECISION    dydt(nnuc)           !Change in relative number abundances.

C-------EVOLUTION PARAMETERS (ORIGINAL VALUES).
      DOUBLE PRECISION    y0(nnuc)             !Rel # abundances at end of 1st R-K loop.

C-------DEFAULT COMPUTATION PARAMETERS.
      DOUBLE PRECISION    cy0                  !Default cy.
      DOUBLE PRECISION    ct0                  !Default ct.
      DOUBLE PRECISION    t9i0                 !Default t9i.
      DOUBLE PRECISION    t9f0                 !Default t9f.
      DOUBLE PRECISION    ytmin0               !Default ytmin.
      INTEGER inc0                 !Default accumulation increment.

C-------COMPUTATION PARAMETERS.
      DOUBLE PRECISION    cy                   !Time step limiting constant on abundances.
      DOUBLE PRECISION    ct                  !Time step limiting constant on temperature.
      DOUBLE PRECISION    t9i                  !Initial temperature (in 10**9 K).
      DOUBLE PRECISION    t9f                  !Final temperature (in 10**9 K).
      DOUBLE PRECISION    ytmin                !Smallest abundances allowed.
      INTEGER inc                  !Accumulation increment.

C-------DEFAULT MODEL PARAMETERS.
      DOUBLE PRECISION    c0(3)                !Default c.
      DOUBLE PRECISION    cosmo0               !Default cosmological constant.
      DOUBLE PRECISION    xi0(3)               !Default neutrino degeneracy parameters.

C-------EARLY UNIVERSE MODEL PARAMETERS.
      DOUBLE PRECISION    g                    !Gravitational constant.
      DOUBLE PRECISION    tau                  !Neutron lifetime (sec).
      DOUBLE PRECISION    xnu                  !Number of neutrino species.
      DOUBLE PRECISION    c(3)               !c(1) is variation of gravitational constant.
     |                             !c(2) is neutron half-life (min).
     |                             !c(3) is number of neutrino species.
      DOUBLE PRECISION    cosmo                !Cosmological constant.
      DOUBLE PRECISION    xi(3)                !Neutrino degeneracy parameters.
     |                             !xi(1) is e neutrino degeneracy parameter.
     |                             !xi(2) is m neutrino degeneracy parameter.
     |                             !xi(3) is t neutrino degeneracy parameter.

C-------DEFAULT VARIATIONAL PARAMETERS.
      DOUBLE PRECISION    dt0                  !Default initial time step.
      DOUBLE PRECISION    eta0                 !Default baryon-to-photon ratio.

C-------VARIATIONAL PARAMETERS.
      DOUBLE PRECISION    dt1                  !Initial time step.
      DOUBLE PRECISION    eta1                 !Baryon-to-photon ratio.

C-------TIME VARIABLES.
      DOUBLE PRECISION    t                    !Time.
      DOUBLE PRECISION    dt                   !Time step.
      DOUBLE PRECISION    dlt9dt               !(1/t9)*d(t9)/d(t).

C-------DYNAMIC VARIABLES.
      DOUBLE PRECISION    thm(14)              !Thermodynamic variables (energy densities).
      DOUBLE PRECISION    hubcst               !Expansion rate of the universe.

C-------ENERGY DENSITIES.
      DOUBLE PRECISION    rhone0               !Initial electron neutrino energy density.
      DOUBLE PRECISION    rhob0                !Initial baryon energy density.
      DOUBLE PRECISION    rhob                 !Baryon energy density.
      DOUBLE PRECISION    rnb                !Baryon energy density (ratio to init value).

C-------MATRIX COEFFICIENTS FOR LINEAR EQUATION.
      DOUBLE PRECISION a(nnuc,nnuc)!Relates y(t+dt) to y(t).
      DOUBLE PRECISION    b(nnuc)              !Contains y0 in inverse order.
      DOUBLE PRECISION    yx(nnuc)             !yy in reverse order.

C-------NUCLIDE DATA.
      DOUBLE PRECISION    am(nnuc)             !Atomic number of nuclide.
      DOUBLE PRECISION    zm(nnuc)             !Charge of nuclide.
      DOUBLE PRECISION    dm(nnuc)             !Mass excess of nuclide.

C-------EVALUATION OF FUNCTIONS bl,bm,bn.
      DOUBLE PRECISION    bl1,bl2,bl3,bl4,bl5  !Evaluation of function bl(z).
      DOUBLE PRECISION    bm1,bm2,bm3,bm4,bm5  !Evaluation of function bm(z).
      DOUBLE PRECISION    bn1,bn2,bn3,bn4,bn5  !Evaluation of function bn(z).

C-------EVALUATION OF MODIFIED BESSEL FUNCTIONS.
      DOUBLE PRECISION    bk0,bk1,bk2,bk3,bk4  !Values k0(r),k1(r),k2(r),k3(r),k4(r).

C-------FLAGS AND COUNTERS.
      INTEGER ltime                !Indicates if output buffer printed.
      INTEGER is                   !# total iterations for particular model.
      INTEGER ip                   !# iterations after outputing a line.
      INTEGER it                   !# times accumulated in output buffer.
      INTEGER mbad                 !Indicates if gaussian elimination failed.

C-------COMPUTATION LOCATION.
      INTEGER itime                !Time check.

C-------OUTPUT ARRAYS.
      DOUBLE PRECISION    xout(itmax,nnuc)     !Nuclide mass fractions.
      DOUBLE PRECISION    thmout(itmax,6)      !Thermodynamic variables.
      DOUBLE PRECISION    t9out(itmax)         !Temperature (in units of 10**9 K).
      DOUBLE PRECISION    tout(itmax)          !Time.
      DOUBLE PRECISION    dtout(itmax)         !Time step.
      DOUBLE PRECISION    etaout(itmax)        !Baryon to photon ratio.
      DOUBLE PRECISION    hubout(itmax)        !Expansion rate.

C-------NEUTRINO PARAMETERS.
      DOUBLE PRECISION    t9mev                !Temperature (in units of MeV).
      DOUBLE PRECISION    tnmev                !Neutrino temperature (in units of MeV).
      DOUBLE PRECISION    tnu                  !Neutrino temperature.
      DOUBLE PRECISION    cnorm                !Normalizing constant.
      DOUBLE PRECISION    rhonu                !Neutrino energy density.
      INTEGER nu                   !Type of neutrino.

C-------RUN OPTION.
      INTEGER irun                 !Run network size.
      INTEGER isize                !Number of nuclides in computation.
      INTEGER jsize                !Number of reactions in computation.

C-------OUTPUT FILE STATUS.
      INTEGER nout                 !Number of output requests.
      LOGICAL outfile              !Indicates if output file used.
      INTEGER outunit
      DOUBLE PRECISION    shout(itmax,nnuc)     !Nuclide mass fractions.


C==================PROCEDURE DIVISION======================

C10-----OPEN FILE---------------------------------------

      IF (itime.eq.1) THEN         !Beginning of program.
        OPEN (unit=3, file='nucint.dat',  status='unknown')
      END IF

C20-----PRINTINTO FILE------------------------------------

      IF (itime.eq.8 .or. itime.eq.10) THEN         !Right after a run.

        IF (outfile) THEN
          outunit = 2
        ELSE
          outunit = 3
        END IF

        shout = xout
        shout(it,8) = shout(it,8) + shout(it,9)  !Add beryllium to lithium.
        shout(it,5) = shout(it,5) + shout(it,4)  !Add tritium to helium-3.
        shout(it,6) = shout(it,6)-0.0003
     |            !Radiative, coulomb, finite-temperature corrections (Ref 1).
        write(outunit,199) "", "eta", "D", "He3", "He4", "Li7"
 199    FORMAT(6(A13, ' '))
        write(outunit,200) "Observables:", etaout(it),shout(it,3),
     |                shout(it,5),shout(it,6),shout(it,8)
     |                             !Output eta, H2, He3, He4, and Li7.
 200    FORMAT (A13, ' ', 5(e13.5,' '))
      END IF

C30-----close FILE--------------------------------------

      IF (itime.eq.10) THEN        !End of program.
        IF (outfile) THEN
          CLOSE (unit=3, status='delete')
        ELSE
          CLOSE (unit=3, status='keep')
        END IF
      END IF
      RETURN

C-------REFERENCES--------------------------------------
C     1) D.A. Dicus, E.W. Kolb, A.M. Gleeson, E.C.G. Sudarshan, V.L. Teplitz,
C        M.S. Turner, Phys. Rev. D., 26,2694 (1982).

      END

