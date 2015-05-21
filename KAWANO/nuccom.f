C===============IDENTIFICATION DIVISION====================

      SUBROUTINE driver

C-------LINKAGES.
C     CALLED BY - [subroutine] run
C     CALLS     - [subroutine] start, derivs, accum

C-------REMARKS.
C     Runge-Kutta computational routine

      USE commons
      USE computation_parameters
      USE flags
      USE evolution_parameters

C-------COMMON AREAS.
      COMMON /time/   t,dt,dlt9dt                    !Time variables.
      COMMON /check1/  itime                          !Computation location.
      COMMON /runopt/ irun,isize,jsize               !Run options.


C=================DECLARATION DIVISION=====================

C-------TIME AND TIME STEP VARIABLES.
      DOUBLE PRECISION    t                    !Time.
      DOUBLE PRECISION    dt                   !Time step.
      DOUBLE PRECISION    dlt9dt               !(1/t9)*d(t9)/d(t).

C-------COUNTERS AND FLAGS.
      INTEGER loop                 !Counts which Runge-Kutta loop.
      INTEGER i

C-------COMPUTATION LOCATION.
      INTEGER itime                !Time check.

C-------RUN OPTION.
      INTEGER isize                !Number of nuclides in computation.

C-------TIME AND TIME STEP VARIABLES.
      DOUBLE PRECISION    dtmin                !Mininum time step.
      DOUBLE PRECISION    dtl                 !Time step from limitation on abund changes.

C-------LABELS FOR VARIABLES TO BE TIME EVOLVED.
      INTEGER mvar                 !Total number of variables to be evolved.
      DOUBLE PRECISION    v(nvar)              !Variables to be time evolved.
      DOUBLE PRECISION    dvdt(nvar)           !Time derivatives.
      DOUBLE PRECISION    v0(nvar)             !Value of variables at original point.
      DOUBLE PRECISION    dvdt0(nvar)          !Value of derivatives at original point.

C-------EQUIVALENCE STATEMENTS.
      EQUIVALENCE (v(4),y(1)),(dvdt(4),dydt(1)),(v0(4),y0(1))


C==================PROCEDURE DIVISION======================

C10-----INPUT INITIALIZATION INFORMATION, RELABEL-------------------

      ltime = 0                    !Set termination indicator to zero.
      CALL start                   !Input initialization information.
      mvar  = isize + 3            !Total number of variables to be evolved.

C20-----LOOP ONE----------------------------------------

 200  continue                     !Begin Runge-Kutta looping.
      loop = 1                     !Loop indicator.
C..........COMPUTE DERIVATIVES OF VARIABLES TO BE EVOLVED.
      CALL derivs(loop)
      itime = 4                    !Time = 1st R-K loop.
      CALL check                   !Check interface subroutine.
C..........ACCUMULATE.
      IF ((t9.le.t9f).or.                         !Low temp.
     |    (dt.lt.abs(min_time_step/dlt9dt)).or.   !Small dt.
     |    (ip.eq.inc)) CALL accum                 !Enough iterations.
C..........POSSIBLY TERMINATE COMPUTATION.
      IF (ltime.eq.1) THEN         !Return to run selection.
        RETURN
      END IF
C..........RESET COUNTERS.
      IF (ip.eq.inc) THEN          !Reset iteration counters.
        ip = 0
      END IF
      ip = ip + 1
      is = is + 1
C..........ADJUST TIME STEP.
      IF (is.gt.3) THEN            !Adjust time step after 3 iterations.
        dtmin = abs(1./dlt9dt)*ct  !Trial value for minimum time step (Ref 1).
        DO i = 1,isize             !Go through all abundance changes.
          IF ((dydt(i).ne.0.).and.(y(i).gt.ytmin)) THEN
            dtl = abs(y(i)/dydt(i))*cy
     |            *(1.+(alog10(y(i))/alog10(ytmin))**2)  !(Ref 2).
            IF (dtl.lt.dtmin) dtmin = dtl         !Find smallest time step.
          END IF
        END DO
        IF (dtmin.gt.1.5*dt) dtmin = 1.5*dt       !Limit change in time step.
        dt = dtmin                 !Set new time step.
      END IF
      t = t + dt                   !Increment time.
C..........STORE AND INCREMENT VALUES (Ref 3).
      DO i = 1,mvar
        v0(i)    = v(i)
        dvdt0(i) = dvdt(i)
        v(i)     = v0(i) + dvdt0(i)*dt
        IF ((i.ge.4).and.(v(i).lt.ytmin)) v(i) = ytmin  !Set at minimum value.
      END DO

C30-----LOOP TWO----------------------------------------

      loop = 2                     !Step up loop counter.
C..........COMPUTE DERIVATIVES OF VARIABLES TO BE EVOLVED.
      CALL derivs(loop)
      itime = 7                    !Time = 2nd R-K loop.
      CALL check                   !Check interface subroutine.
C..........INCREMENT VALUES.
      DO i = 1,mvar
        v(i) = v0(i) + .5*(dvdt(i)+dvdt0(i))*dt
        IF ((i.ge.4).and.(v(i).lt.ytmin)) v(i) = ytmin  !Set at minimum value.
      END DO

      GO TO 200

C-------REFERENCES--------------------------------------
C     1)  Constraint on dt from the requirement that
C                (d(t9)/dt)*(dt/t9) < ct
C         Wagoner, R.V. 1969, Ap J. Suppl. No. 162, 18, page 293, equation C6.
C     2)  Constraint on dt from
C                dtl < y/(dy/dt)*cy*(1+(log(y)/log(ytmin))**2)
C         Wagoner, R.V. 1969, page 293, equation C7 but with log term squared.
C     3)  Wagoner, R.V. 1969, page 292, equations C1, C2.

      END



C===============IDENTIFICATION DIVISION====================

      SUBROUTINE start

C-------LINKAGES.
C     CALLED BY - [subroutine] driver
C     CALLS     - [subroutine] rate1, bessel, rate0
C               - [function] ex

C-------REMARKS.
C     Sets initial conditions.

      USE commons
      USE model_parameters
      USE computation_parameters
      USE variational_parameters
      USE flags
      USE evolution_parameters
      USE sterile

C-------COMMON AREAS.
      COMMON /rates/  f,r                            !Reaction rates.
      COMMON /time/   t,dt,dlt9dt                    !Time variables.
      COMMON /endens/ rhone0,rhob0,rhob,rnb          !Energy densities.
      COMMON /xbessel/ bl1,bl2,bl3,bl4,bl5,          !Eval of function bl(z).
     |                bm1,bm2,bm3,bm4,bm5,           !Eval of function bm(z).
     |                bn1,bn2,bn3,bn4,bn5            !Eval of function bn(z).
      COMMON /nupar/  t9mev,tnmev,tnu,cnorm,nu,rhonu !Neutrino parameters.
      COMMON /runopt/ irun,isize,jsize               !Run options.

C=================DECLARATION DIVISION=====================

C-------REACTION RATES.
      DOUBLE PRECISION    f(nrec)              !Forward reaction rate coefficients.
      DOUBLE PRECISION r(nrec)

C-------TIME VARIABLES.
      DOUBLE PRECISION    t                    !Time.
      DOUBLE PRECISION    dt                   !Time step.

C-------ENERGY DENSITIES.
      DOUBLE PRECISION    rhone0               !Initial electron neutrino mass density.
      DOUBLE PRECISION    rhob0                !Initial baryon mass density.

C-------EVALUATION OF FUNCTIONS bl,bm,bn.
      DOUBLE PRECISION    bl1,bl2,bl3,bl4,bl5  !Evaluation of function bl(z).

C-------NEUTRINO PARAMETERS.
      DOUBLE PRECISION    tnu                  !Neutrino temperature.
      DOUBLE PRECISION    cnorm                !Normalizing constant.

C-------RUN OPTION.
      INTEGER isize                !Number of nuclides in computation.

C-------LOCAL VARIABLES.
      DOUBLE PRECISION    z                 !Defined by z = m(electron)*c**2/k*t9.

c Julien modified, 28-02-08
      DOUBLE PRECISION tr,xx,HH,t9r,dt9r,rhor,
     |     r1r,r2r,r3r,r4r,r5r,r6r !Variables used to read the data file
      DOUBLE PRECISION t_interp             !This is used to obtain the correct time values
      CHARACTER trash           !Used to remove the first line from the data
c Julien end mod 28-02-08

      INTEGER :: i, e
      CHARACTER(255) :: input_file

C==================PROCEDURE DIVISION======================

C10-----INITIALIZE FLAGS AND COUNTERS-------------------------

      ltime = 0                    !No output yet.
      is    = 1                    !First iteration coming up.
      ip    = inc                  !Set to maximum allowed # of iteration.
      it    = 0                    !No accumulation yet.
      mbad  = 0                    !No computational errors.

C20-----SETTINGS----------------------------------------

C..........COMPUTATIONAL SETTINGS.
      t9  = t9i                    !Initial temperature.
      tnu = t9                     !Initial neutrino temperature.
      t   = 1/(const1*t9)**2       !Initial time (Ref 1).
      dt  = dt1                    !Initial time step.

C..........MODEL SETTINGS.
      g   = const2*c(1)            !Modify gravitational constant.
      tau = c(2)                 !Convert n half-life (min) to lifetime (sec).
      tau = tau/0.98               !Coulomb correction (Ref 2).
      xnu = c(3)                   !Number of neutrino species.

C30-----COMPUTE INITIAL ABUNDANCES FOR NEUTRON AND PROTON--------------

      IF ((15.011/t9+xi(1)).gt.58.) THEN      !Overabundance of antineutrinos.
        y(1) = 1.e-25              !Very little of neutrons.
        y(2) = 1.                  !Essentially all protons.
      ELSE
        IF ((15.011/t9+xi(1)).lt.-58.) THEN   !Overabundance of neutrinos.
          y(1) = 1.                !Essentially all neutrons.
          y(2) = 1.e-25            !Very little of protons.
        ELSE
          y(1) = 1./(ex(15.011/t9+xi(1))+1.)  !Initial n abundance (Ref 3).
          y(2) = 1./(ex(-15.011/t9-xi(1))+1.) !Initial p abundance (Ref 3).
        END IF
      END IF
      IF (xi(1).ne.0.) THEN        !Electron neutrino degeneracy.
        cnorm = 1.
        tnu   = .00001             !Low temperature.
        CALL rate1(0.00001)        !Find normalization constant at low temp.
        cnorm = 1/tau/f(1)
        print *, cnorm
      END IF
      y0(1) = y(1)
      y0(2) = y(2)

C40-----FIND RATIO OF BARYON DENSITY TO TEMPERATURE CUBED--------------

      z      = 5.930/t9            !Inverse of temperature.
      CALL bessel(z)
      hv     = 3.3683e+4*eta1*2.75 !(Ref 4 but with final eta).
      phie   = hv*(1.784e-5*y(2))  !Chemical potential of electron (Ref 5).
     |            /(.5*z**3*(bl1-2.*bl2+3.*bl3-4.*bl4+5.*bl5))
      rhob0  = hv*t9**3            !Baryon density.
      IF ((xi(1).eq.0.).and.(xi(2).eq.0.).and.(xi(3).eq.0)) THEN  !Nondegen.
        rhone0 = 7.366*t9**4       !Electron neutrino density (Ref 6).
      END IF

C50-----SET ABUNDANCES FOR REST OF NUCLIDES----------------------

      y(3)  = y(1)*y(2)*rhob0*ex(25.82/t9)/(.471e+10*t9**1.5)  !(Ref 7).
      y0(3) = y(3)
      DO i = 4,isize
        y(i)  = ytmin              !Set rest to minimum abundance.
        y0(i) = y(i)               !Init abundances at beginning of iteration.
      END DO
      CALL rate0                   !Compute weak decay rates.

c Julien modified, 28-02-08 - 05-03-08
C60----READ OUTPUT FILE FROM OTHER PROGRAM AND SAVE DATA -------

      CALL GETENV("INPUT", input_file)
      input_file = TRIM(input_file)
      IF (input_file == '') THEN
        input_file = 's4.dat'
      END IF
      PRINT *, 'Reading custom data from: ', input_file
      OPEN(1,file=input_file, status='old') !Here give the name of the outputted file
     |                        ! from the other program
      READ(1,*, IOSTAT=e) trash !Get the number of written lines
      IF (e == 0) THEN
        READ(1,*) trash           !Remove the second line that are labels
      END IF

      nlines=0
      DO
        READ(1,*,end=10) tr,xx,t9r,dt9r,rhor,HH,r1r,r2r,r3r,r4r,r5r,r6r
        nlines=nlines + 1
        ts(nlines) = tr
        t9s(nlines) = t9r
        dt9s(nlines) = dt9r
        rho_tot(nlines) = rhor
        ratef(nlines) = (r1r + r3r + r5r) !Will divide them by tau later on
        rater(nlines) = (r2r + r4r + r6r)
      END DO
c Saves the different used datas from the file into the corresponding
c  arrays
10    CLOSE(1)

      rater = rater / ratef(nlines)
      ratef = ratef / ratef(nlines)
c These divisions are for normalisation so that for T->0 (taken to be our
c  last outputted values), the neutron decay rate (ratef) is equal to 1
c  (and lated on we'll divide by 1/tau so that this rate is equal to
c  the real neutron decay rate). Such a normalisation is already done in
c  the other program so that ratef(nlines) was already nearly equal to 1,
c  but the previous normalization is not taking into account the change that
c  could come from adding sterile neutrinos.

      call log_interp_values(t_interp, t9, t9s, ts, .true., nlines)
      ts = ts - (t_interp - t)
c In these 2 lines we do the following: the value from our time ts is
c  not the same as the one used by Kawano (not same initial time condition)
c  so we determine how to put our time in agreement with his. First, the
c  initial time is determined by the initial temperature, so we determine
c  in our data to which time this initial temperature corresponds (first
c  line). And then we subtract to our datas the value of the difference
c  between our time and the one from Kawano, so that at this temperature
c  we have exactly the same time value.
c But this time is in fact not really used in our modifications.

c Julien end mod 28-02-08 - 05-03-08


      RETURN

C-------REFERENCES--------------------------------------
C     1) Wagoner, R.V., Fowler, W.A., and Hoyle, F. 1967, Ap. J. 148,
C          page 44, equation A15.
C     2) Coulomb correction obtained by dividing by correction factor Fp(t9)
C               Fp(t9) = 1 - 0.5(pi/(137<v>/c))
C          Wagoner, R.V. 1973, Ap. J. 179, page 358.
C     3) For the nondegenerate case:
C          Wagoner, R.V., Fowler, W.A., and Hoyle, F. 1967, Ap. J. 148,
C          page 4, equation 3.
C        For the case with neutrino degeneracy:
C          Beaudet,G. and Goret,P., 1976, Astron. & Astrophys., 49,
C          page 417, equation 9.
C     4) Wagoner, R.V. 1969, Ap J. Suppl. No. 162, 18, page 250, equation 4.
C          3.3683e+4 = Mu(ng/t9**3) with Mu the atomic mass, ng the
C          photon density.  2.75 is for the 11/4 factor difference
C          between the initial and final values of eta.
C     5) Kawano, L., 1992, Fermilab preprint FERMILAB-PUB-92/04-A,
C          Kellogg Radiation Lab preprint OAP-714,
C          equation D.2.
C     6) Wagoner, R.V., Fowler, W.A., and Hoyle, F. 1967, Ap. J. 148,
C          page 43, equation A4.
C          7.366 is used instead of 14.73 as the latter is the sum total
C          for 2 neutrino species.
C     7) Initial deuterium abundance from nuclear statistical equilibrium
C          Wagoner, R.V., Fowler, W.A., and Hoyle, F. 1967, Ap. J. 148,
C          page 19, equation 17.

      END



C===============IDENTIFICATION DIVISION====================

      SUBROUTINE derivs(loop)

C-------LINKAGES.
C     CALLED BY - [subroutine] driver
C     CALLS     - [subroutine] therm, rate1, rate4, rate3, rate2, sol

C-------REMARKS.
C     Computes derivatives of
C       - Temperature
C       - hv
C       - Chemical potential
C       - abundances

      USE commons
      USE model_parameters
      USE evolution_parameters
      USE sterile

C-------COMMON AREAS.
      COMMON /time/   t,dt,dlt9dt                    !Time variables.
      COMMON /xtherm/  thm,hubcst                     !Dynamic variables.
      COMMON /endens/ rhone0,rhob0,rhob,rnb          !Energy densities.
      COMMON /nucdat/ am(nnuc),zm,dm                 !Nuclide data.
      COMMON /runopt/ irun,isize,jsize               !Run options.


C=================DECLARATION DIVISION=====================

C-------TIME VARIABLES.
      DOUBLE PRECISION    dlt9dt               !(1/t9)*d(t9)/d(t).

C-------DYNAMIC VARIABLES.
      DOUBLE PRECISION    thm(14)              !Thermodynamic variables.
      DOUBLE PRECISION    hubcst               !Expansion rate.

C-------ENERGY DENSITIES.
      DOUBLE PRECISION    rhob0                !Initial baryon mass density.
      DOUBLE PRECISION    rhob                 !Baryon mass density.
      DOUBLE PRECISION    rnb                  !Baryon mass density (ratio to init value).

C-------NUCLIDE DATA.
      DOUBLE PRECISION    zm(nnuc)             !Charge of nuclide.
      DOUBLE PRECISION    dm(nnuc)             !Mass excess of nuclide.

C-------RUN OPTION.
      INTEGER irun                 !Run network size.
      INTEGER isize                !Number of nuclides in computation.

C-------SUMS.
      DOUBLE PRECISION    sumy                 !Sum of abundances.
      DOUBLE PRECISION    sumzy                !Sum of charge*abundances.
      DOUBLE PRECISION    sumdy                !Sum of abundance flows.
      DOUBLE PRECISION    summdy               !Sum of (mass excess)*(abundance flows).
      DOUBLE PRECISION    sumzdy               !Sum of (charge)*(abundance flows).

C-------DERIVATIVES.
      DOUBLE PRECISION    dphdt9               !d(phi e)/d(t9).
      DOUBLE PRECISION    dphdln               !d(phi e)/d(h).
      DOUBLE PRECISION    dphdzy               !d(phi e)/d(sumzy).
      DOUBLE PRECISION    dlndt9               !(1/h)*d(h)/d(t9).
      DOUBLE PRECISION    bar                  !Baryon density and pressure terms.

C-------LOCAL VARIABLES.
      INTEGER loop              !Counts which Runge-Kutta loop.
c Julien modified, 05-03-08
      DOUBLE PRECISION dt9_interp           !The interpolated dt9 value from dt9s
c Julien end mod 05-03-08


C==================PROCEDURE DIVISION======================

C10-----COMPUTE DERIVATIVES FOR ABUNDANCES-----------------------

      rnb    = hv*t9*t9*t9/rhob0   !Baryon mass density (ratio to init value).
C..........VARIOUS THERMODYNAMIC QUANTITIES.

      CALL therm
      hubcst = sqrt((8./3.)*pi*g*(thm(10))+(cosmo/3.))  !Expansion rate.
      rhob   = thm(9)             !Baryon mass density.
C..........COMPUTE REACTION RATE COEFFICIENTS.

      CALL rate1(t9)
      GO TO (100,110,120), irun    !Run network selection.
 100  CONTINUE
        CALL rate4                 !Forward rate for all of reactions.
 110  CONTINUE
        CALL rate3                 !Forward rate for reactions with A < 19.
 120  CONTINUE
        CALL rate2                 !Forward rate for reactions with A < 10.
C..........SOLVE COUPLED DIFFERENTIAL EQUATIONS.
      CALL sol(loop)
      IF (mbad.gt.0) RETURN        !Abort in case matrix not invertible.

C20-----COMPUTE DERIVATIVES FOR TEMPERATURE, hv, AND CHEMICAL POTENTIAL------

C..........INITIALIZE SUMS TO ZERO.
      sumy   = 0.
      sumzy  = 0.
      sumdy  = 0.
      summdy = 0.
      sumzdy = 0.
C..........ACCUMULATE TO GET SUM.
      DO i = 1,isize
        sumy   = sumy   + y(i)           !Sum of abundance.
        sumzy  = sumzy  + zm(i)*y(i)     !Sum of charge*abundance.
        sumdy  = sumdy  + dydt(i)        !Sum of abundance flow.
        summdy = summdy + dm(i)*dydt(i) !Sum of (mass excess)*(abundance flow).
        sumzdy = sumzdy + zm(i)*dydt(i)  !Sum of (charge)*(abundance flow).
      END DO
C..........CHANGES IN TEMPERATURE, hv, AND CHEMICAL POTENTIAL.
      dphdt9 = thm(12)*(-1.070e-4*hv*sumzy/t9 - thm(11))
      dphdln = -thm(12)*3.568e-5*hv*sumzy
      dphdzy = thm(12)*3.568e-5*hv
      bar    = 9.25e-5*t9*sumy + 1.388e-4*t9*sumdy/(3.*hubcst)
     |         + summdy/(3.*hubcst)
      dlndt9 = -(thm(2) + thm(5) + thm(6)*dphdt9 + thm(9)*1.388e-4*
     |         sumy)/(thm(1) + thm(3) + thm(4) + thm(7) + thm(9)*bar
     |         + thm(6)*(dphdln + dphdzy*sumzdy/(3.*hubcst)))   !(Ref 1).

c Julien modified, 29-02-08 - 05-03-08
c  was:
c      dt9    = (3.*hubcst)/dlndt9
c  we use:
      call log_interp_values(dt9_interp, t9, t9s, dt9s, .true., nlines)
      dt9 = dt9_interp
c Julien end mod 05-03-08

      dlt9dt = dt9/t9
      dhv    = -hv*((3.*hubcst) + 3.*dlt9dt)                    !(Ref 2).
      dphie  = dphdt9*dt9 + dphdln*(3.*hubcst) + dphdzy*sumzdy  !(Ref 3).

      RETURN

C-------REFERENCES--------------------------------------
C     1)  Kawano, L., 1992, Fermilab preprint FERMILAB-PUB-92/04-A,
C          Kellogg Radiation Lab preprint OAP-714,
C          equation D.35.
C     2)  Kawano, L. 1992, preprint, equation D.19.
C     3)  Kawano, L. 1992, preprint, equation D.20.

      END



C===============IDENTIFICATION DIVISION====================

      SUBROUTINE accum

C-------LINKAGES.
C     CALLED BY - [subroutine] driver
C     CALLS     - none

C-------REMARKS.
C     Output accumulator.

      USE commons
      USE computation_parameters
      USE flags
      USE evolution_parameters

C-------COMMON AREAS.
      COMMON /time/   t,dt,dlt9dt                    !Time variables.
      COMMON /xtherm/  thm,hubcst                     !Dynamic variables.
      COMMON /nucdat/ am,zm(nnuc),dm(nnuc)           !Nuclide data.
      COMMON /outdat/ xout,thmout,t9out,tout,dtout,  !Output data.
     |                etaout,hubout
      COMMON /runopt/ irun,isize,jsize               !Run options.


C=================DECLARATION DIVISION=====================

C-------TIME PARAMETERS.
      DOUBLE PRECISION    t                    !Time.
      DOUBLE PRECISION    dt                   !Time step.

C-------DYNAMIC VARIABLES.
      DOUBLE PRECISION    thm(14)              !Thermodynamic variables.
      DOUBLE PRECISION    hubcst               !Expansion rate.

C-------NUCLIDE DATA.
      DOUBLE PRECISION    am(nnuc)             !Atomic number of nuclide.

C-------OUTPUT ARRAYS.
      DOUBLE PRECISION    xout(itmax,nnuc)     !Nuclide mass fractions.
      DOUBLE PRECISION    thmout(itmax,6)      !Thermodynamic variables.
      DOUBLE PRECISION    t9out(itmax)         !Temperature (in units of 10**9 K).
      DOUBLE PRECISION    tout(itmax)          !Time.
      DOUBLE PRECISION    dtout(itmax)         !Time step.
      DOUBLE PRECISION    etaout(itmax)        !Baryon-to-photon ratio.
      DOUBLE PRECISION    hubout(itmax)        !Expansion rate.

C-------RUN OPTION.
      INTEGER isize                !Number of nuclides in computation.


C==================PROCEDURE DIVISION======================

      it = it + 1                  !Set up accumulation counter.

C10-----SET UP OUTPUT VARIABLES-------------------------------

C..........DIVIDE NUMBER FRACTION BY THAT OF PROTON.
      DO i = 1,isize
        xout(it,i) = y(i)/y(2)
      END DO
      xout(it,2) = y(2)*am(2)      !Exception for proton.
      xout(it,6) = y(6)*am(6)      !Exception for helium.
C..........SUM UP ABUNDANCES OF HEAVY NUCLIDES.
      xout(it,10) =  xout(it,10)+xout(it,11)+xout(it,12)+xout(it,13)
     |              +xout(it,14)+xout(it,15)+xout(it,16)+xout(it,17)
     |              +xout(it,18)+xout(it,19)+xout(it,20)+xout(it,21)
     |              +xout(it,22)+xout(it,23)+xout(it,24)+xout(it,25)
     |              +xout(it,26)   !Li8 to O16.
C..........RELABEL TEMPERATURE, TIME, THERMODYNAMIC VARIABLES, ETC.
      t9out(it)    = t9            !Temperature.
      tout(it)     = t             !Time.
      thmout(it,1) = thm(1)        !rho photon.
      thmout(it,2) = thm(4)        !rho electron.
      thmout(it,3) = thm(8)        !rho neutrino.
      thmout(it,4) = thm(9)        !rho baryon.
      thmout(it,5) = phie          !Chemical potential.
      thmout(it,6) = thm(10)       !rho total.
      dtout(it)    = dt            !Time step.
      etaout(it)   = hv/(3.3683e+4)!Baryon to photon ratio.
      hubout(it)   = hubcst        !Expansion rate.

C20-----INDICATE TERMINATION OF ACCUMULATION IF APPROPRIATE------------

      IF ((it.eq.itmax).or.(ip.lt.inc)) ltime = 1
      RETURN

      END



C===============IDENTIFICATION DIVISION====================

      SUBROUTINE therm

C-------LINKAGES.
C     CALLED BY - [subroutine] derivs
C     CALLS     - [subroutine] bessel, nudens
C               - [function] ex

C-------REMARKS.
C     Computes various temperature dependent thermodynamic quantities.

      USE commons
      USE model_parameters
      USE computation_parameters
      USE evolution_parameters
      USE sterile

C-------COMMON AREAS.
      COMMON /xtherm/  thm,hubcst                     !Dynamic variables.
      COMMON /endens/ rhone0,rhob0,rhob,rnb          !Energy densities.
      COMMON /xbessel/ bl1,bl2,bl3,bl4,bl5,           !Eval of function bl(z).
     |                bm1,bm2,bm3,bm4,bm5,           !Eval of function bm(z).
     |                bn1,bn2,bn3,bn4,bn5            !Eval of function bn(z).
      COMMON /nupar/  t9mev,tnmev,tnu,cnorm,nu,rhonu !Integration parameters.


C=================DECLARATION DIVISION=====================

C-------DYNAMIC VARIABLES.
      DOUBLE PRECISION    thm(14)              !Thermodynamic variables.

C-------ENERGY DENSITIES.
      DOUBLE PRECISION    rhone0               !Initial electron neutrino mass density.
      DOUBLE PRECISION    rhob0                !Initial baryon mass density.
      DOUBLE PRECISION    rnb                  !Baryon mass density (ratio to init value).

C-------EVALUATION OF FUNCTIONS bl,bm,bn.
      DOUBLE PRECISION    bl1,bl2,bl3,bl4,bl5  !Evaluation of function bl(z).
      DOUBLE PRECISION    bm1,bm2,bm3,bm4,bm5  !Evaluation of function bm(z).
      DOUBLE PRECISION    bn1,bn2,bn3,bn4,bn5  !Evaluation of function bn(z).

C-------NEUTRINO PARAMETERS.
      DOUBLE PRECISION    tnu                  !Neutrino temperature.
      DOUBLE PRECISION    rhonu                !Neutrino energy density.
      INTEGER nu                   !Type of neutrino.

C-------LOCAL VARIABLE.
      DOUBLE PRECISION    z                    !Defined by z = m(electron)*c**2/k*t9.
c Julien modified, 28-02-08
      DOUBLE PRECISION rho,rate             !Total energy density at given time, and
     |                          ! rate
c Julien end mod 28-02-08


C==================PROCEDURE DIVISION======================

C10-----COMPUTE FACTORS------------------------------------

      z = 5.930/t9                 !z = m(electron)c**2/k(t9).
      tnu = ((rnb)**(1./3.))*t9i   !Neutrino temperature.
C..........FACTORS OF z.
      z1 = z
      z2 = z*z
      z3 = z2*z
      z4 = z3*z
      z5 = z4*z
C..........TRIGNOMETRIC FUNCTION VALUES.
      IF (phie.le.17.) THEN        !No chance of overflow.
        cosh1 = cosh(phie)
        cosh2 = cosh(2.*phie)
        cosh3 = cosh(3.*phie)
        cosh4 = cosh(4.*phie)
        cosh5 = cosh(5.*phie)
        sinh1 = sinh(phie)
        sinh2 = sinh(2.*phie)
        sinh3 = sinh(3.*phie)
        sinh4 = sinh(4.*phie)
        sinh5 = sinh(5.*phie)
      ELSE
        cosh1 = 0.
        cosh2 = 0.
        cosh3 = 0.
        cosh4 = 0.
        cosh5 = 0.
        sinh1 = 0.
        sinh2 = 0.
        sinh3 = 0.
        sinh4 = 0.
        sinh5 = 0.
      END IF
      CALL bessel(z)

C20-----COMPUTE THERMODYNAMIC VARIABLES--------------------------
      thm(1)  = 8.418*t9*t9*t9*t9                               !(Ref 1).
      thm(2)  = 4.*thm(1)/t9                                    !(Ref 2).
      thm(3)  = thm(1)/3.                                       !(Ref 3).
      thm(4)  = 3206.*(bm1*cosh1 - bm2*cosh2 + bm3*cosh3        !(Ref 4).
     |          - bm4*cosh4 + bm5*cosh5)
      thm(5)  = 3206.*(z/t9)*(bn1*cosh1 - 2.*bn2*cosh2          !(Ref 5).
     |          + 3.*bn3*cosh3 - 4.*bn4*cosh4 + 5.*bn5*cosh5)
      thm(6)  = 3206.*(bm1*sinh1 - 2.*bm2*sinh2 + 3.*bm3*sinh3  !(Ref 6).
     |          - 4.*bm4*sinh4 + 5.*bm5*sinh5)
      thm(7)  = 3206.*(bl1*cosh1/z - bl2*cosh2/(2.*z)           !(Ref 7).
     |          + bl3*cosh3/(3.*z) - bl4*cosh4/(4.*z)
     |          + bl5*cosh5/(5.*z))
      IF ((xi(1).eq.0.).and.(xi(2).eq.0.).and.(xi(3).eq.0)) THEN  !Nondegen.
        thm(8) = xnu*rhone0*(rnb**(4./3.))                      !(Ref 8).
      ELSE                         !Include effects of neutrino degeneracy.
        thm(8) = 0.
        DO nu = 1, INT(xnu)        !For every neutrino family.
          CALL nudens              !Compute neutrino energy density.
          thm(8) = thm(8) + 12.79264*rhonu  !Have 12.79264 from units change.
        END DO
      END IF
c Julien: 29-02-08: Note that I didn't change the neutrino density,
c  only the total density, so this neutrino density is not totally
c  the correct one in our case but in the code it is needed only
c  for the total density, so it doesn't matter

      thm(9)  = rhob0*rnb                                       !(Ref 9).

c Julien modified, 28-02-08
c  was:
c      thm(10) = thm(1) + thm(4) + thm(8) + thm(9)               !(Ref 10).
c  we use:
      call log_interp_values(rho, t9, t9s, rho_tot, .true., nlines)
      thm(10) = rho + thm(9)                                    !(Ref 10).
c Julien end mod 28-02-08

      thm(11) = -(z**3/t9)*(sinh1*(3.*bl1-z*bm1)-sinh2*(3.*bl2  !(Ref 11).
     |          -2.*z*bm2) + sinh3*(3.*bl3-3.*z*bm3) - sinh4
     |          *(3.*bl4-4.*z*bm4) + sinh5*(3.*bl5-5.*z*bm5))
      thm(12) = z**3*(cosh1*bl1- 2.*cosh2*bl2                   !(Ref 12).
     |          + 3.*cosh3*bl3 - 4.*cosh4*bl4 + 5.*cosh5*bl5)
      IF (thm(12).ne.0.) thm(12) = 1./thm(12)

c Julien modified, 28-02-08
c  was:
c      thm(13) = 1.000 + 0.565/z1 - 6.382/z2 + 11.108/z3 !(Ref 13).
c     |          + 36.492/z4 + 27.512/z5
c      thm(14) = (5.252/z1 - 16.229/z2 + 18.059/z3 + 34.181/z4   !(Ref 14).
c     |          + 27.617/z5)*ex(-q*z)
c  we use:
      call log_interp_values(rate, t9, t9s, ratef, .true., nlines)
      thm(13) = rate                                            !(Ref 13).
      call log_interp_values(rate, t9, t9s, rater, .true., nlines)
      thm(14) = rate                                            !(Ref 14).
c Julien end mod 28-02-08

      RETURN

C-------REFERENCES AND NOTES-------------------------------
C     1)  thm(1)  = rho photon
C         (Wagoner, R.V., Fowler, W.A., and Hoyle, F. 1967, Ap. J. 148,
C          page 43, equation A2.)
C     2)  thm(2)  = d(rho photon)/d(t9)
C     3)  thm(3)  = (p photon)/c**2
C         (Wagoner, R.V., Fowler, W.A., and Hoyle, F. 1967,
C          page 43, equation A3.)
C     4)  thm(4)  = rho electron+positron
C         (Fowler, W.A. and Hoyle, F., 1964, Ap. J. Suppl. No. 91, 9,
C          page 281, equation B44.)
C     5)  thm(5)  = d(rho electron+positron)/d(t9)
C     6)  thm(6)  = d(rho electron+positron)/d(phi e)
C     7)  thm(7)  = (p electron+positron)/c**2
C         (Fowler, W.A. and Hoyle, F., 1964, Ap. J. Suppl. No. 91, 9,
C          page 279, equation B27.)
C     8)  thm(8)  = rho neutrino
C                 = # neutrino species x rho electron neutrino (nondegenerate)
C                 = rho nu(e) + rho nu(m) + rho nu(t)          (degenerate)
C     9)  thm(9)  = rho baryon
C     10) thm(10) = rho total
C                 = rho photon + rho electron+positron + rho neutrino
C                              + rho baryon
C     11) thm(11) = d     /pi**2(hbar*c)**3(ne- - ne+)*z**3\
C                   d(t9) \  2  (mc**2)**3                 /
C     12) thm(12) = d        /pi**2(hbar*c)**3(ne- - ne+)*z**3\
C                   d(phi e) \  2  (mc**2)**3                 /
C     13) thm(13) = rate for n->p
C     14) thm(14) = rate for p->n

      END



C===============IDENTIFICATION DIVISION====================

      SUBROUTINE bessel(z)

C-------LINKAGES.
C     CALLED BY - [subroutine] start, therm
C     CALLS     - [subroutine] knux

C-------REMARKS.
C     Evaluates functions bl(z), bm(z), and bn(z) using solutions to
C     modified Bessel functions.

C-------COMMON AREAS.
      COMMON /xbessel/ bl1,bl2,bl3,bl4,bl5,           !Eval function bl(z).
     |                bm1,bm2,bm3,bm4,bm5,           !Eval function bm(z).
     |                bn1,bn2,bn3,bn4,bn5            !Eval function bn(z).
      COMMON /kays/   bk0,bk1,bk2,bk3,bk4            !Coefficients K.


C=================DECLARATION DIVISION=====================

C-------EVALUATION OF FUNCTIONS bl,bm,bn.
      DOUBLE PRECISION   bl1,bl2,bl3,bl4,bl5  !Single variables equivalenced to array blz.
      DOUBLE PRECISION   bm1,bm2,bm3,bm4,bm5  !Single variables equivalenced to array bmz.
      DOUBLE PRECISION   bn1,bn2,bn3,bn4,bn5  !Single variables equivalenced to array bnz.

C-------EVALUATIION OF MODIFIED BESSEL FUNCTIONS.
      DOUBLE PRECISION    bk0,bk1,bk2,bk3,bk4  !Values k0(r),k1(r),k2(r),k3(r),k4(r).

C-------LOCAL VARIABLES.
      DOUBLE PRECISION    blz(5)               !Array containing values from function bl.
      DOUBLE PRECISION    bmz(5)               !Array containing values from function bm.
      DOUBLE PRECISION    bnz(5)               !Array containing values from function bn.
      DOUBLE PRECISION    z                    !Defined by z = m(electron)*c**2/k*t9.
      DOUBLE PRECISION    r                    !Multiples of z.

C-------EQUIVALENCE STATEMENTS.
      EQUIVALENCE (blz(1),bl1),(blz(2),bl2),(blz(3),bl3),(blz(4),bl4),
     |            (blz(5),bl5)
      EQUIVALENCE (bmz(1),bm1),(bmz(2),bm2),(bmz(3),bm3),(bmz(4),bm4),
     |            (bmz(5),bm5)
      EQUIVALENCE (bnz(1),bn1),(bnz(2),bn2),(bnz(3),bn3),(bnz(4),bn4),
     |            (bnz(5),bn5)


C==================PROCEDURE DIVISION======================

C10-----LOCALLY DEFINED FUNCTIONS-----------------------------

      bl(z) = bk2/z                !Function bl.
      bm(z) = 0.25*(3.*bk3+bk1)/z  !Function bm.
      bn(z) = 0.5*(bk4+bk2)/z      !Function bn.

C20-----CALCULATE FOR 1 THRU 5 Z------------------------------

      DO i=1,5
        r=i*z                      !Multiples of z.
        CALL knux(r)               !Get k0(r),k1(r),k2(r),k3(r),k4(r),k(5).
        blz(i) = bl(r)             !Put value from function bl into array.
        bmz(i) = bm(r)             !Put value from function bm into array.
        bnz(i) = bn(r)             !Put value from function bn into array.
      END DO
      RETURN

      END



C===============IDENTIFICATION DIVISION====================

      SUBROUTINE knux(z)

C-------LINKAGES.
C     CALLED BY - [subroutine] bessel
C     CALLS     - [function] exp

C-------REMARKS.
C     A subroutine for modified bessel functions of the second kind
C     k-nu(z).

C-------COMMON AREAS.
      COMMON /kays/   bk0,bk1,bk2,bk3,bk4            !Coefficients K.


C==================DECLARATION DIVISION====================

C--------MODIFIED BESSEL FUNCTION VALUES.
      DOUBLE PRECISION    bk0,bk1              !Values k0(z),k1(z)
      DOUBLE PRECISION    bi0,bi1              !Values i0(z),i1(z).
      DOUBLE PRECISION    bk2,bk3,bk4          !Values k2(z),k3(z),k4(z).

C--------EXPANSION COEFFICIENTS.
      DOUBLE PRECISION    ci0(7)               !Expansion coefficients for i0 (z.le.2).
      DOUBLE PRECISION    ci1(7)               !Expansion coefficients for i1 (z.le.2).
      DOUBLE PRECISION    ck0(7)               !Expansion coefficients for k0 (z.le.2).
      DOUBLE PRECISION    ck1(7)               !Expansion coefficients for k1 (z.le.2).
      DOUBLE PRECISION    c0(7)                !Expansion coefficients for k0 (z.gt.2).
      DOUBLE PRECISION    c1(7)                !Expansion coefficients for k1 (z.gt.2).

C--------VARIABLES TO BE EVALUATED.
      DOUBLE PRECISION    z                    !Input variable.
      DOUBLE PRECISION    y                    !Expansion variable = z/2.
      DOUBLE PRECISION    t                    !Expansion variable = z/3.75.
      DOUBLE PRECISION    coeff                !Logrithmic or exponential coefficient.


C=====================DATA DIVISION========================

C-------EXPANSION COEFFICIENTS.
      DATA ci0 / 1.,
     |           3.5156229,      3.0899424,      1.2067492,
     |           0.2659732,      0.0360768,      0.0045813/
      DATA ci1 / 0.5,
     |           0.87890594,     0.51498869,     0.15084934,
     |           0.02658733,     0.00301532,     0.00032411/
      DATA ck0 /-0.57721566,
     |           0.42278420,     0.23069756,     0.03488590,
     |           0.00262698,     0.00010750,     0.00000740/
      DATA ck1 / 1.,
     |           0.15443144,    -0.67278579,    -0.18156897,
     |          -0.01919402,    -0.00110404,    -0.00004686/
      DATA c0  / 1.25331414,
     |          -0.07832358,     0.02189568,    -0.01062446,
     |           0.00587872,    -0.00251540,     0.00053208/
      DATA c1  / 1.25331414,
     |           0.23498619,    -0.03655620,     0.01504268,
     |          -0.00780353,     0.00325614,    -0.00068245/


C==================PROCEDURE DIVISION======================

C10-----COMPUTE K0 AND K1----------------------------------

      IF (z.le.2.) THEN            !(Ref. 1).
C..........COMPUTE FACTORS.
        t = (z/3.75)
        y = (z/2)
        coeff = alog(y)
C..........VALUES FOR i0(z) and i1(z).
        bi0 = ci0(1)
        bi1 = ci1(1)
        bk0 = ck0(1)
        bk1 = ck1(1)
        DO i = 2,7
          bi0 = bi0 + ci0(i)*t**(2*(i-1))
          bi1 = bi1 + ci1(i)*t**(2*(i-1))
          bk0 = bk0 + ck0(i)*y**(2*(i-1))
          bk1 = bk1 + ck1(i)*y**(2*(i-1))
        END DO
C..........VALUES FOR k0(z) and k1(z).
        bk0 = -coeff*bi0 + bk0
        bk1 = coeff*bi1*z + bk1/z
      ELSE !(z.le.2.)               !(Ref. 2).
C..........COMPUTE FACTORS.
        y = (2.0/z)
        coeff = (ex(-z)/sqrt(z))
C..........VALUES FOR k0(z) and k1(z).
        bk0 = c0(1)
        bk1 = c1(1)
        DO i = 2,7
          bk0 = bk0 + c0(i)*y**(i-1)
          bk1 = bk1 + c1(i)*y**(i-1)
        END DO
        bk0 = coeff*bk0
        bk1 = coeff*bk1
      END IF !(z.le.2.)

C20-----FIND K2, K3, AND K4 BY ITERATION (Ref. 3)-------------------

      bk2 = 2.*(bk1/z) + bk0       !k2(z).
      bk3 = 4.*(bk2/z) + bk1       !k3(z).
      bk4 = 6.*(bk3/z) + bk2       !k4(z).

      RETURN

C-------REFERENCES--------------------------------------
C     Handbook of Mathematical Functions (Abramowitz and Stegun),
C       Dover Publications, Inc., New York
C       1) Polynomial approximations for z.le.2
C         page 378, equations 9.8.1 and 9.8.3.
C         page 379, equations 9.8.5 and 9.8.7.
C       2) Polynomial approximations for z > 2
C         page 379, equations 9.8.6 and 9.8.8.
C       3) Recursion relation from 1st line of 9.6.26, page 376.

      END



C===============IDENTIFICATION DIVISION====================

      SUBROUTINE nudens

C-------LINKAGES.
C     CALLED BY - [subroutine] therm
C     CALLS     - [function] xintd, eval

C-------REMARKS.
C     Computes energy density contribution from neutrinos.

C-------COMMON AREAS.
      COMMON /nupar/  t9mev,tnmev,tnu,cnorm,nu,rhonu !Integration parameters.

C-------EXTERNAL FUNCTIONS.
      EXTERNAL func5               !Integral for neutrinos.
      EXTERNAL func6               !Integral for antineutrinos.


C=================DECLARATION DIVISION=====================

C-------MODEL PARAMETERS.
      DOUBLE PRECISION    xi(3)                !Neutrino degeneracy parameters.
      DOUBLE PRECISION func5
      DOUBLE PRECISION func6

C-------NEUTRINO PARAMETERS.
      DOUBLE PRECISION    tnu                  !Neutrino temperature (units of 10**9 K).
      DOUBLE PRECISION    rhonu                !Neutrino energy density.
      INTEGER nu                   !Which neutrino type.

C-------LOCAL VARIABLES.
      DOUBLE PRECISION    uplim1               !Upper limit for neutrino energy integral.
      DOUBLE PRECISION    uplim2               !Upper limit for antineu energy integral.


C==================PROCEDURE DIVISION======================

C10-----COMPUTE NEUTRINO ENERGY DENSITIES------------------------

      IF (abs(xi(nu)).le.0.03) THEN
C..........SMALL xi APPROXIMATION.
        rhonu = 2.*(3.14159**2/30.)*(tnu)**4
     |          *(7./8.+(15./(4*3.14159**2))*xi(nu)**2
     |          +(15./(8.*3.14159**4))*xi(nu)**4)
      ELSE
        IF (abs(xi(nu)).ge.30.) THEN
C..........LARGE xi APPROXIMATION.
          rhonu = ((tnu)**4)/(8.*3.14159**2)*xi(nu)**4
     |            *(1+12.*1.645 /xi(nu)**2)
        ELSE
C..........DO INTEGRATION
          uplim1 = (88.029+xi(nu))*tnu
          uplim2 = (88.029-xi(nu))*tnu
          IF (uplim2.le.0.) THEN
            rhonu = xintd(0.,uplim1,func5,iter)
          ELSE
            rhonu= xintd(0.,uplim1,func5,iter)
     |             + xintd(0.,uplim2,func6,iter)
          END IF
        END IF !(abs(xi(nu)).ge.30.)
      END IF !(abs(xi(nu)).le.0.03)


      RETURN

C-------REFERENCES--------------------------------------
C     Forms of the integrals involved can be found in
C       Beaudet,G. and Goret,P., 1976, Astron. & Astrophys., 49, 415.
C       Freese, K., Kolb, E.W., Turner, M.S., 1983, Phys. Rev. D, 27, 1689.

      END



C===============IDENTIFICATION DIVISION====================

      FUNCTION eval()

C-------LINKAGES.
C     CALLED BY - [subroutine] rate1, nudens
C     CALLS     - [function] ex

C-------REMARKS.
C     Contains integrands to be integrated.

C-------COMMON AREAS.
      COMMON /nupar/  t9mev,tnmev,tnu,cnorm,nu,rhonu !Integration parameters.

C=================DECLARATION DIVISION=====================

C-------MODEL PARAMETERS.
      DOUBLE PRECISION    xi(3)                !Neutrino degeneracy parameters.

C-------NEUTRINO PARAMETERS.
      DOUBLE PRECISION    t9mev                !Temperature (in units of MeV).
      DOUBLE PRECISION    tnmev                !Neutrino temperature (in units of MeV).
      DOUBLE PRECISION    tnu                  !Neutrino temperature (units of 10**9 K).
      DOUBLE PRECISION    cnorm                !Normalizing constant.

C-------FUNCTIONS TO BE EVALUATED.
      DOUBLE PRECISION   func1                 !1st part n->p rate.
      DOUBLE PRECISION   func2                 !2nd part n->p rate.
      DOUBLE PRECISION   func3                 !1st part p->n rate.
      DOUBLE PRECISION   func4                 !2nd part p->n rate.
      DOUBLE PRECISION   func5                 !Energy density for neutrinos.
      DOUBLE PRECISION   func6                 !Energy density for antineutrinos.

C-------LOCAL VARIABLES.
      DOUBLE PRECISION    x                    !Value at which function is evaluated.
      DOUBLE PRECISION    part1                !Exponential expression with photon temp.
      DOUBLE PRECISION    part2                !Exponential expression with neutrino temp.


C==================PROCEDURE DIVISION======================

C10-----1ST PART OF INTEGRAL FOR n->p RATE-----------------------

      ENTRY func1(x)
      IF (x.le.0.) THEN
        func1 = 0.
      ELSE
        part1 = 1./(1.+ex(-.511*x/(t9mev)))
        part2 = 1./(1.+ex(+(x-2.531)*(.511/(tnmev))-xi(1)))
        func1 = cnorm*x*(x-2.531)**2*(x**2-1)**.5*part1*part2
      END IF
      RETURN

C20-----2ND PART OF INTEGRAL FOR n->p RATE-----------------------

      ENTRY func2(x)
      IF (x.le.1.) THEN
        func2 = 0.
      ELSE
        part1 = 1./(1.+ex(+.511*x/(t9mev)))
        part2 = 1./(1.+ex(-(x+2.531)*(.511/(tnmev))-xi(1)))
        func2 = cnorm*x*(x+2.531)**2*(x**2-1)**.5*part1*part2
      END IF
      RETURN

C30-----1ST PART OF INTEGRAL FOR p->n RATE-----------------------

      ENTRY func3(x)
      IF (x.le.1.) THEN
        func3 = 0.
      ELSE
        part1 = 1./(1.+ex(-.511*x/(t9mev)))
        part2 = 1./(1.+ex(+(x+2.531)*(.511/(tnmev))+xi(1)))
        func3 = cnorm*x*(x+2.531)**2*(x**2-1)**.5*part1*part2
      END IF
      RETURN

C40-----2ND PART OF INTEGRAL FOR p->n RATE-----------------------

      ENTRY func4(x)
      IF (x.le.1.) THEN
        func4 = 0.
      ELSE
        part1 = 1./(1.+ex(+.511*x/(t9mev)))
        part2 = 1./(1.+ex(-(x-2.531)*(.511/(tnmev))+xi(1)))
        func4 = cnorm*x*(x-2.531)**2*(x**2-1)**.5*part1*part2
      END IF
      RETURN

C50-----INTEGRAL FOR ENERGY DENSITY OF NEUTRINO---------------------

      ENTRY func5(x)
      func5 = 1./(2*3.14159**2)*x**3/(1.+exp(x/tnu-xi(nu)))
      RETURN

C60-----INTEGRAL FOR ENERGY DENSITY OF ANTINEUTRINO-----------------

      ENTRY func6(x)
      func6 = 1./(2*3.14159**2)*x**3/(1.+exp(x/tnu+xi(nu)))
      RETURN

C-------REFERENCES--------------------------------------
C     Forms of the integrals involved can be found in
C       Scherrer,R.J., 1983, Mon.Not.R.astr.Soc., 205, 683.
C       Beaudet,G. and Goret,P., 1976, Astron. & Astrophys., 49, 415.

      END



C===============IDENTIFICATION DIVISION====================

      FUNCTION xintd (xlow,xhi,func,nq)

C-------LINKAGES.
C     CALLED BY - [subroutine] rate1, nudens
C     CALLS     - none

C-------REMARKS.
C     Computes the integral of the function "func".


C=================DECLARATION DIVISION=====================

C-------INPUT VARIABLES.
      external func
      DOUBLE PRECISION    xlow                 !Array of low limits.
      DOUBLE PRECISION    xhi                  !Array of high limits.
      INTEGER nq                   !Number of six point gaussian quads.

C-------COMPUTATION VARIABLES.
      DOUBLE PRECISION    dist                 !Size of quad interval.
      DOUBLE PRECISION    cent                 !Center of quad interval.
      DOUBLE PRECISION    x                    !Variables of integration.
      DOUBLE PRECISION    sum                  !Summation of terms.

C-------COUNTERS.
      INTEGER nint                 !Interval number.
      INTEGER npnt                 !Point number.
      INTEGER np                   !Total number of points in interval.

C-------ABSCISSAS AND WEIGHT FACTORS.
      DOUBLE PRECISION    u(6)                 !Abscissas.
      DOUBLE PRECISION    w(6)                 !Weight factor.


C=====================DATA DIVISION========================

C-------ABSCISSAS AND WEIGHT FACTORS.
      DATA u/-.93246951420315,-.66120938646627,-.23861918608320,
     |        .23861918608320, .66120938646627, .93246951420315/
      DATA w/.17132449237917,.36076157304814,.46791393457269,
     |       .46791393457269,.36076157304814,.17132449237917/
      DATA np/6/              !6 point Gaussian integration.


C==================PROCEDURE DIVISION======================

C10-----DO INTEGRATION-------------------------------------

      sum   = 0.
      dist  = (xhi-xlow)/float(nq) !Size of quad interval.
      DO nint = 1,nq
        cent = xlow+(float(nint)-0.5)*dist  !Center of interval.
        DO npnt = 1,np
          x   = cent+0.5*dist*u(npnt) !Integration point.
          f   = func(x)            !Evaluate function x(1).
          sum = sum+f*w(npnt)      !Add up sum.
        END DO
      END DO

C20-----GET INTEGRAL VALUE---------------------------------

      xintd = sum*dist*0.5         !Do integral.
      RETURN

      END



C===============IDENTIFICATION DIVISION====================

      FUNCTION ex(x)

C-------LINKAGES.
C     CALLED BY - [subroutine] start, rate2, rate3, rate4, sol
C               - [function] eval
C     CALLS     - none

C-------REMARKS.
C     Exponential function with underflow precaution.


C==================PROCEDURE DIVISION======================

      IF (x.gt.88.029) THEN        !In danger of overflow.
        ex = exp(88.029)
      ELSE
        IF (x.lt.-88.722) THEN     !In danger of underflow.
          ex = 0.
        ELSE                       !Value of x in allowable range.
          ex = exp(x)
        END IF
      END IF
      RETURN

C-------NOTE-----------------------------------------
C     The overflow limit for the VAX/VMS system is exp(88.029).
C     The underflow limit for the VAX/VMS system is exp(-88.722).

      END



C===============IDENTIFICATION DIVISION====================

      SUBROUTINE sol(loop)

C-------LINKAGES.
C     CALLED BY - [subroutine] derivs
C     CALLS     - [subroutine] eqslin
C               - [function] ex

C-------REMARKS.
C     Computes reverse strong and electromagnetic reaction rates.
C     Fills and solves matrix equation for dydt(i).

      USE commons
      USE evolution_parameters

C--------COMMON AREAS.
      COMMON /recpr/  iform,ii,jj,kk,ll,rev,q9     !Reaction parameters names.
      COMMON /rates/  f,r                            !Reaction rates.
      COMMON /time/   t,dt,dlt9dt                    !Time varying parameters.
      COMMON /xtherm/  thm(14),hubcst                 !Dynamic variables.
      COMMON /endens/ rhone0,rhob0,rhob,rnb          !Energy densities.
      COMMON /lncoef/ a,b,yx                         !Linear eqn coefficients.
      COMMON /flags/  ltime,is,ip,it,mbad            !Flags,counters.
      COMMON /runopt/ irun,isize,jsize               !Run option.


C=================DECLARATION DIVISION=====================

C-------REACTION PARAMETERS.
      INTEGER iform(nrec)          !Reaction code number (1-88).
      INTEGER ii(nrec)             !Incoming nuclide type (1-26).
      INTEGER jj(nrec)             !Incoming light nuclide type (1-6).
      INTEGER kk(nrec)             !Outgoing light nuclide type (1-6).
      INTEGER ll(nrec)             !Outgoing nuclide type (1-26).
      DOUBLE PRECISION    rev(nrec)            !Reverse reaction coefficient.
      DOUBLE PRECISION    q9(nrec)             !Energy released in reaction.

C-------REACTION RATES.
      DOUBLE PRECISION    f(nrec)              !Forward reaction rate coefficients.
      DOUBLE PRECISION    r(nrec)              !Reverse reaction rate coefficients.

C-------TIME VARIABLES.
      DOUBLE PRECISION    dt                   !Time step.

C-------DYNAMIC VARIABLES.
      DOUBLE PRECISION    hubcst               !Expansion rate.

C-------ENERGY DENSITIES.
      DOUBLE PRECISION    rhob                 !Baryon mass density.

C-------COMPONENTS OF MATRIX EQUATION.
      DOUBLE PRECISION a(nnuc,nnuc) !Relates y(t-dt) to y(t).
      DOUBLE PRECISION    b(nnuc)              !Contains y0 in inverse order.
      DOUBLE PRECISION    yx(nnuc)             !yy in reverse order.

C-------COUNTERS AND FLAGS.
      INTEGER loop                 !Counts which Runge-Kutta loop.
      INTEGER ip                   !# time steps after outputting a line.
      INTEGER mbad                 !Indicates if gaussian elimination fails.

C-------RUN OPTIONS.
      INTEGER isize                !Number of nuclides in computation.
      INTEGER isize1               !Equals isize + 1.
      INTEGER jsize                !Number of reactions in computation.

C-------EVOLUTION EQUATION COEFFICIENTS.
      INTEGER i,j,k,l              !Equate to ii,jj,kk,ll.
      DOUBLE PRECISION    ri,rj,rk,rl          !Equate to si,sj,sk,sl.
      DOUBLE PRECISION    ci,cj,ck,cl          !Coefficients of rate equation.

C-------LOCAL VARIABLES.
      DOUBLE PRECISION    yy(nnuc)             !Abundances at end of iteration.
      DOUBLE PRECISION    si(11),sj(11),sk(11),sl(11)  !# of nuclide i,j,k,l
      DOUBLE PRECISION    bdln                 !(10**(-5))*volume expansion rate.
      INTEGER ind                  !Equate to iform.
      INTEGER ierror               !Element which does not converge.


C=====================DATA DIVISION========================

C-------NUMBER OF NUCLIDES IN REACTION TYPES 1-11.
      DATA si /1.,1.,1.,1.,1.,2.,3.,2.,1.,1.,2./
      DATA sj /0.,1.,1.,0.,1.,0.,0.,1.,1.,1.,0./
      DATA sk /0.,0.,1.,0.,0.,1.,0.,0.,1.,0.,2./
      DATA sl /1.,1.,1.,2.,2.,1.,1.,1.,2.,3.,1./


C==================PROCEDURE DIVISION======================

C10-----TEMPERATURE FACTORS AND INITIAL VALUES--------------

C..........TEMPERATURE FACTORS.
      t932  = t9**1.5              !t9**(3/2).
      t9m32 = 1./t932              !t9**(-3/2).
C..........MATRIX SIZE.
      isize1 = isize + 1
C..........INITIALIZE A-MATRIX.
      DO i = 1,isize
        DO j = 1,isize
          a(j,i) = 0.d0            !Set a-matrix to zero.
        END DO
      END DO

C20-----COMPUTE FACTORS FOR THE A-MATRIX-------------------------

      DO n = 1,jsize
C..........EQUATE VARIABLES TO ARRAYS.
        ind = iform(n)             !Type of reaction.
        i = ii(n)                  !ID # of incoming nuclide i.
        j = jj(n)                  !ID # of incoming nuclide j.
        k = kk(n)                  !ID # of outgoing nuclide k.
        l = ll(n)                  !ID # of outgoing nuclide l.
        IF ((ind.ne.0).and.(i.le.isize).and.(l.le.isize)) THEN !Reaction okay.
          ri = si(ind)             !# of incoming nuclide i.
          rj = sj(ind)             !# of incoming nuclide j.
          rk = sk(ind)             !# of outgoing nuclide k.
          rl = sl(ind)             !# of outgoing nuclide l.
C..........COMPUTE DIFFERENT REACTION RATES.
          GO TO (201,202,203,204,205,206,207,208,209,210,211),ind
 201      CONTINUE                 !1-0-0-1 configuration.
            ci = f(n)              !(Ref 1).
            cj = 0.
            ck = 0.
            cl = r(n)
            GO TO 212
 202      CONTINUE                 !1-1-0-1 configuration.
            r(n) = rev(n)*1.e+10*t932*ex(-q9(n)/t9)*f(n)  !(Ref 2).
            f(n) = rhob*f(n)
            ci = y(j)*f(n)/2.
            cj = y(i)*f(n)/2.
            ck = 0.
            cl = r(n)
            GO TO 212
 203      CONTINUE                 !1-1-1-1 configuration.
            f(n) = rhob*f(n)
            r(n) = rev(n)*ex(-q9(n)/t9)*f(n)  !(Ref 3).
            ci = y(j)*f(n)/2.
            cj = y(i)*f(n)/2.
            ck = y(l)*r(n)/2.
            cl = y(k)*r(n)/2.
            GO TO 212
 204      CONTINUE                 !1-0-0-2 configuration.
            ci = f(n)
            cj = 0.
            ck = 0.
            cl = y(l)*r(n)/2.
            GO TO 212
 205      CONTINUE                 !1-1-0-2 configuration.
            f(n) = rhob*f(n)
            r(n) = rev(n)*ex(-q9(n)/t9)*f(n)  !(Ref 3).
            ci = y(j)*f(n)/2.
            cj = y(i)*f(n)/2.
            ck = 0.
            cl = y(l)*r(n)/2.
            GO TO 212
 206      CONTINUE                 !2-0-1-1 configuration.
            f(n) = rhob*f(n)
            r(n) = rev(n)*ex(-q9(n)/t9)*f(n)  !(Ref 3).
            ci = y(i)*f(n)/2.
            cj = 0.
            ck = y(l)*r(n)/2.
            cl = y(k)*r(n)/2.
            GO TO 212
 207      CONTINUE                 !3-0-0-1 configuration.
            r(n) = rev(n)*1.e+20*t932*t932*ex(-q9(n)/t9)*f(n)  !(Ref 4).
            f(n) = rhob*rhob*f(n)
            ci = y(i)*y(i)*f(n)/6.
            cj = 0.
            ck = 0.
            cl = r(n)
            GO TO 212
 208      CONTINUE                 !2-1-0-1 configuration.
            r(n) = rev(n)*1.e+20*t932*t932*ex(-q9(n)/t9)*f(n)  !(Ref 4).
            f(n) = rhob*rhob*f(n)
            ci = y(j)*y(i)*f(n)/3.
            cj = y(i)*y(i)*f(n)/6.
            ck = 0.
            cl = r(n)
            GO TO 212
 209      CONTINUE                 !1-1-1-2 configuration.
            f(n) = rhob*f(n)
            r(n) = rev(n)*1.e-10*t9m32*rhob*ex(-q9(n)/t9)*f(n)  !(Ref 5).
            ci = y(j)*f(n)/2.
            cj = y(i)*f(n)/2.
            ck = y(l)*y(l)*r(n)/6.
            cl = y(k)*y(l)*r(n)/3.
            GO TO 212
 210      CONTINUE                 !1-1-0-3 configuration.
            f(n) = rhob*f(n)
            r(n) = rev(n)*1.e-10*t9m32*rhob*ex(-q9(n)/t9)*f(n)  !(Ref 5).
            ci = y(j)*f(n)/2.
            cj = y(i)*f(n)/2.
            ck = 0.
            cl = y(l)*y(l)*r(n)/6.
            GO TO 212
 211      CONTINUE                 !2-0-2-1 configuration.
            f(n) = rhob*f(n)
            r(n) = rev(n)*1.e-10*t9m32*rhob*ex(-q9(n)/t9)*f(n)  !(Ref 5).
            ci = y(i)*f(n)/2.
            cj = 0.
            ck = y(l)*y(k)*r(n)/3.
            cl = y(k)*y(k)*r(n)/6.
 212      CONTINUE

C30-----CONSTRUCT THE A-MATRIX--------------------------------

          i = isize1 - i           !Invert i index.
          j = isize1 - j           !Invert j index.
          k = isize1 - k           !Invert k index.
          l = isize1 - l           !Invert l index.
C..........FILL I NUCLIDE COLUMN.
          IF (j.le.isize) a(j,i) = a(j,i) +  rj*ci
          IF (k.le.isize) a(k,i) = a(k,i) -  rk*ci
          a(i,i) = a(i,i) +  ri*ci
          a(l,i) = a(l,i) -  rl*ci
C..........FILL J NUCLIDE COLUMN.
          IF (j.le.isize) THEN
            a(j,j) = a(j,j) +  rj*cj
            IF (k.le.isize) a(k,j) = a(k,j) -  rk*cj
            a(i,j) = a(i,j) +  ri*cj
            a(l,j) = a(l,j) -  rl*cj
          END IF
C..........FILL K NUCLIDE COLUMN.
          IF (k.le.isize) THEN
            IF (j.le.isize) a(j,k) = a(j,k) -  rj*ck
            a(k,k) = a(k,k) +  rk*ck
            a(i,k) = a(i,k) -  ri*ck
            a(l,k) = a(l,k) +  rl*ck
          END IF
C..........FILL L NUCLIDE COLUMN.
          IF (j.le.isize) a(j,l) = a(j,l) -  rj*cl
          IF (k.le.isize) a(k,l) = a(k,l) +  rk*cl
          a(i,l) = a(i,l) -  ri*cl
          a(l,l) = a(l,l) +  rl*cl
        END IF !((ind.ne.0).and.(i.le.isize).and.(l.le.isize))
      END DO !n = 1,jsize

C40-----PUT A-MATRIX AND B-VECTOR IN FINAL FORM OF MATRIX EQUATION--------

      bdln   = 1.e-5*(3.*hubcst)   !(10**(-5))*(Expansion rate).
      DO i = 1,isize
        i1 = isize1 - i            !Invert the rows.
        DO j = 1,isize
          j1 = isize1 - j          !Invert the columns.
          IF (ABS(a(j,i)).lt.bdln*y0(j1)/y0(i1)) THEN
            a(j,i) = 0.d0          !Set 0 if tiny.
          ELSE
            a(j,i) = a(j,i)*dt     !Bring dt over to other side.
          END IF
        END DO
        a(i,i) = 1.d0 + a(i,i)     !Add identity matrix to a-matrix.
        b(i1)  = y0(i)             !Initial abundances.
      END DO

C50-----SOLVE EQUATIONS TO GET DERIVATIVE------------------------

C..........SET MONITOR FLAG AND SOLVE BY GAUSSIAN ELIMINATION.
      IF (loop.eq.1) THEN
        CALL eqslin(ip,ierror)
      ELSE
        CALL eqslin(0,ierror)
      END IF
C..........OBTAIN DERIVATIVE.
      DO i = 1,isize
        yy(i)   = yx(isize1-i)     !Abundance at t+dt.
        dydt(i) = (yy(i) - y0(i))/dt         !Take derivative.
      END DO

C60-----POSSIBLE ERROR MESSAGES AND EXIT-------------------------

      IF (mbad.ne.0) THEN          !Problem in gaussian elimination.
        IF (mbad.eq.-1) print 6000, ierror !Error message.
        IF (mbad.ge. 1) print 6002, mbad   !Error message.
 6000   FORMAT (' ','** y(', i2, ') fails to converge **')
 6002   FORMAT (' ','** ', i2, ' th diagonal term equals zero **')
      END IF
      RETURN

C-------REFERENCES--------------------------------------
C     1) The coefficients are given in general as:
C             ci = ri*(y(j)**rj)*(y(i)**(ri-1)*f(n)/
C                  ((ri+rj)*fac(ri)*fac(rj))
C             cj = rj*(y(i)**ri)*(y(j)**(rj-1)*f(n)/
C                  ((ri+rj)*fac(ri)*fac(rj))
C             ck = rk*(y(l)**rl)*(y(k)**(rk-1)*f(n)/
C                  ((rk+rl)*fac(rk)*fac(rl))
C             cl = rl*(y(k)**rk)*(y(l)**(rl-1)*f(n)/
C                  ((rk+rl)*fac(rk)*fac(rl))
C        in which fac(x) is the factorial of x.
C     2) Form of reverse rate given in
C        Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247,
C          tables 1B, 4B, 7B.
C     3) Form of reverse rate given in
C        Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247,
C          tables 2B, 3B, 5B, 6B, 8B, 9B, 10B.
C     4) Form of reverse rate given in
C        Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247,
C          table 11B.
C     5) Form of reverse rate given in
C        Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247,
C          tables 12B, 13B.


      END



C===============IDENTIFICATION DIVISION====================

      SUBROUTINE eqslin(icnvm,ierror)

C-------LINKAGES.
C     CALLED BY - [subroutine] sol
C     CALLS     - none

C-------REMARKS.
C     Solves for new abundances using gaussian elimination
C     with back substitution, no pivoting.

      USE commons

C-------COMMON AREAS.
      COMMON /compr/  cy,ct,t9i,t9f,ytmin,inc        !Computation parameters.
      COMMON /lncoef/ a,b,y                          !Lin eqn coefficients.
      COMMON /flags/  ltime,is,ip,it,mbad            !Flags, counters.
      COMMON /runopt/ irun,isize,jsize               !Run options.


C=================DECLARATION DIVISION=====================

C-------COMPUTATION PARAMETER.
      INTEGER inc                  !Accumulation increment.

C-------MATRIX COEFFICIENTS FOR LINEAR EQUATION.
      DOUBLE PRECISION a(nnuc,nnuc)!Coefficient array.
      DOUBLE PRECISION    b(nnuc)              !Right-hand vector w/o manipulation.
      DOUBLE PRECISION    y(nnuc)              !Solution vector.

C-------COUNTERS AND FLAGS.
      INTEGER mbad                 !Indicates type of error.

C-------RUN OPTION.
      INTEGER isize                !Number of nuclides in computation.

C-------LOCAL MATRICES AND VECTORS.
      DOUBLE PRECISION a0(nnuc,nnuc)!Coefficient array w/o manipulation.
      DOUBLE PRECISION x(nnuc)     !Right-hand vector.

C-------LOCAL COMPUTATION VARIABLES.
      DOUBLE PRECISION cx          !Scaling factor in triangularization.
      DOUBLE PRECISION sum         !Sum for backsubstitution.
      DOUBLE PRECISION   xdy                   !Relative error.

C-------LOCAL COUNTERS.
      INTEGER nord                 !Order of correction.
      INTEGER icnvm                !Convergence monitor.
      INTEGER ierror               !ith nuclide fails to converge.


C==================PROCEDURE DIVISION======================

C10-----INITIALIZE VECTOR----------------------------------

C..........SET COUNTERS TO ZERO.
      nord = 0                     !No corrections yet.
      mbad = 0                     !No errors yet.
C..........SET RIGHT-HAND AND SOLUTION VECTORS TO INITIAL VALUES.
      DO i = 1,isize
        x(i) = b(i)                !Right-hand vector.
        y(i) = 0.                  !Solution vector.
      END DO
C..........SAVE MATRIX.
      IF (icnvm.eq.inc) THEN       !Monitor convergence.
        DO i = 1,isize
          DO j = 1,isize
            a0(j,i) = a(j,i)       !Initial value of coefficient array.
          END DO
        END DO
      END IF

C20-----TRIANGULARIZE MATRIX AND SAVE OPERATOR----------------------

C..........CHECK TO SEE THAT THERE ARE NO ZEROES AT PIVOT POINTS.
      DO i = 1,isize-1
        IF (a(i,i).eq.0.d0) THEN   !Don't want to divide by zero.
          mbad = i                 !Position of zero coefficient.
          RETURN                   !Terminate matrix evaluation.
        END IF
C..........TRIANGULARIZE MATRIX.
        DO j = i+1,isize
          IF (a(j,i).ne.0.d0) THEN !Progress diagonally down the column.
            cx = a(j,i)/a(i,i)     !Scaling factor down the column.
            DO k = i+1,isize       !Progress diagonally along row.
              a(j,k) = a(j,k) - cx*a(i,k)  !Subtract scaled coeff along row.
            END DO
            a(j,i) = cx            !Scaled coefficient.
C..........OPERATE ON RIGHT-HAND VECTOR.
            x(j) = x(j) - cx*x(i)  !Subtract off scaled coefficient.
          END IF
        END DO
      END DO

C30-----DO BACK SUBSTITUTION-------------------------------

 300  CONTINUE
      x(isize) = x(isize)/a(isize,isize)   !Solution for ultimate position.
      y(isize) = y(isize) + x(isize)
      DO i = isize-1,1,-1          !From i = penultimate to i = 1.
        sum = 0.d0
        DO j = i+1,isize
          sum = sum + a(i,j)*x(j)  !Sum up all previous terms.
        END DO
        x(i) = (x(i) - sum)/a(i,i)
        y(i) = y(i) + x(i)         !Add difference to initial value.
      END DO

C40-----TESTS AND EXITS------------------------------------

      IF (icnvm.eq.inc) THEN
        DO i = 1,isize
          IF (y(i).ne.0.) THEN
            xdy = ABS(x(i)/y(i))  !Relative error.
            IF (xdy.gt.eps) THEN
              IF (nord.lt.mord) THEN !Continue to higher orders.
                nord = nord + 1
C..........FIND ERROR IN RIGHT-HAND VECTOR.
                DO j = 1,isize
                  r = 0.d0         !Initialize r.
                  DO k = 1,isize
                    r = r + a0(j,k)*y(k) !Left side with approximate solution.
                  END DO
                  x(j) = b(j) - r  !Subtract difference from right side.
                END DO
C..........OPERATE ON RIGHT-HAND VECTOR.
                DO j = 1,isize-1
                  DO k = j+1,isize
                   x(k) = x(k) - a(k,j)*x(j) !Subtract off scaled coefficient.
                  END DO
                END DO
                GO TO 300       !Go for another iteratiion.
              ELSE
C..........NOT ENOUGH CONVERGENCE.
                mbad = -1          !Signal error problem.
                ierror = i         !ith nuclide for which x/y checked.
                RETURN
              END IF !(nord.lt.mord)
            END IF !(xdy.gt.eps)
          END IF !(y(i).ne.0)
        END DO !i = 1,isize
      END IF !(icnvm.eq.inc)
      RETURN                       !No more iterations & relative error small.

      END
