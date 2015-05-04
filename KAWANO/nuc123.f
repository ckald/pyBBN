 
      PROGRAM nuc123

C----------LINKAGES.
C     CALLED BY - none
C     CALLS     - [subroutine] help, setcom, setmod, run, output

C----------REMARKS.
C     Control program -
C       Offers user the main menu and channels through to various options.
C     Implementation -
C       To run this program, NUC123.FOR must be linked with NUCCOM 
C       (containing the computation subroutines), NUCRAT (with the
C       reaction rates), and NUCINT (with an interface subroutine).
C       This program has been written to be compatible with
C       ANSI FORTRAN-77 with the exception of the END DO statement
C       used to limit the number of statement labels.
C       The code was developed on the VAX/VMS system.
C     Notes -
C       The program utilizes Wagoner's code as the core of the computational
C       routines.
C     Documentation -
C       Kawano, L., 1992, Fermilab preprint FERMILAB-PUB-92/04-A,
C       Kellogg Radiation Lab preprint OAP-714.
C     Copy -
C       Version 4.1 (December 1991)

C----------PARAMETERS.
      PARAMETER (ir=1)             !Input unit number.
      PARAMETER (iw=1)             !Output unit number.
      PARAMETER (nrec=88)          !Number of nuclear reactions.
      PARAMETER (nnuc=26)          !Number of nuclides in calculation.

C----------COMMON AREAS.
      COMMON /recpr0/ reacpr                        !Reaction parameter values.
      COMMON /recpr/  iform,ii,jj,kk,ll,rev,q9       !Reaction parameter names.
      COMMON /rates/  f,r                            !Reaction rates.
      COMMON /compr0/ cy0,ct0,t9i0,t9f0,ytmin0,inc0  !Default comp parameters.
      COMMON /compr/  cy,ct,t9i,t9f,ytmin,inc        !Computation parameters.
      COMMON /modpr0/ c0,cosmo0,xi0                  !Default model parameters.
      COMMON /modpr/  g,tau,xnu,c,cosmo,xi           !Model parameters.
      COMMON /varpr0/ dt0,eta0                      !Default variationl params.
      COMMON /varpr/  dt1,eta1                       !Variational parameters.
      COMMON /check1/  itime                          !Computation location.
      COMMON /runopt/ irun,isize,jsize               !Run options.
      COMMON /outopt/ nout,outfile                   !Output option.


C==========================DECLARATION DIVISION============================

C----------REACTION PARAMETERS FROM BLOCK DATA.
      REAL    reacpr(nrec,8)       !Reaction parameters.

C----------REACTION PARAMETERS.
      INTEGER iform(nrec)          !Reaction type code (1-11).
      INTEGER ii(nrec)             !Incoming nuclide type (1-26).
      INTEGER jj(nrec)             !Incoming light nuclide type (1-6).
      INTEGER kk(nrec)             !Outgoing light nuclide type (1-6).
      INTEGER ll(nrec)             !Outgoing nuclide type (1-26).
      REAL    rev(nrec)            !Reverse reaction coefficient.
      REAL    q9(nrec)             !Energy released in reaction.

C----------REACTION RATES.
      REAL    f(nrec)              !Forward reaction rate coefficients.
      REAL    r(nrec)              !Reverse reaction rate coefficients.

C----------DEFAULT COMPUTATION PARAMETERS.
      REAL    cy0                  !Default cy.
      REAL    ct0                  !Default ct.
      REAL    t9i0                 !Default t9i.
      REAL    t9f0                 !Default t9f.
      REAL    ytmin0               !Default ytmin.
      INTEGER inc0                 !Default accumulation increment.

C----------COMPUTATIONAL PARAMETERS.
      REAL    cy                   !Time step limiting constant on abundances.
      REAL    ct                   !Time step limiting constant on temperature.
      REAL    t9i                  !Initial temperature (in 10**9 K).
      REAL    t9f                  !Final temperature (in 10**9 k).
      REAL    ytmin                !Smallest abundances allowed.
      INTEGER inc                  !Accumulation increment.

C----------DEFAULT MODEL PARAMETERS.
      REAL    c0(3)                !Default c.
      REAL    cosmo0               !Default cosmological constant.
      REAL    xi0(3)               !Default neutrino degeneracy parameters.

C----------EARLY UNIVERSE MODEL PARAMETERS.
      REAL    c(3)                !c(1) is variation of gravitational constant.
     |                             !c(2) is neutron lifetime (sec).
     |                             !c(3) is number of neutrino species.
      REAL    cosmo                !Cosmological constant.
      REAL    xi(3)                !Neutrino degeneracy parameters.

C----------DEFAULT VARIATIONAL PARAMETERS.
      REAL    dt0                  !Default initial time step.
      REAL    eta0                 !Default baryon-to-photon ratio.

C----------VARIATIONAL PARAMETERS.
      REAL    dt1                  !Initial time step.
      REAL    eta1                 !Baryon-to-photon ratio.

C----------COMPUTATION LOCATION.
      INTEGER itime                !Time check.

C----------RUN OPTION.
      INTEGER irun                 !Run network size.
      INTEGER isize                !Number of nuclides in computation.
      INTEGER jsize                !Number of reactions in computation.

C----------OUTPUT FILE STATUS.
      INTEGER nout                 !Number of output requests.
      LOGICAL outfile              !Indicates if output file used.

C----------USER RESPONSE VARIABLES.
      INTEGER inum                 !Selection number.


C===========================PROCEDURE DIVISION=====================

C--------OPEN FILES AND PRINT GREETING-----------------------------

      OPEN (unit=2, file='nuc123.dat', status='unknown')  !Output file.
      itime = 1                    !Time = beginning of program.
      CALL check                   !Check interface subroutine.
      PRINT 1000
 1000 FORMAT (6(/),
     |        2(' ',4x,'NN',6x,'NN  UU',6x,'UU',4x,8('C'),6x,'11',8x,
     |        6('2'),6x,6('3'),/),
     |        2(' ',4x,'NN',6x,'NN  UU',6x,'UU  CC',12x,'1111',6x,
     |        '22',6x,'22  33',6x,'33',/),
     |        2(' ',4x,'NNNN    NN  UU',6x,'UU  CC',14x,'11',14x,
     |        '22',10x,'33',/),
     |        2(' ',4x,'NN  NN  NN  UU',6x,'UU  CC',14x,'11',12x,
     |        '22',10x,'33',/),
     |        2(' ',4x,'NN    NNNN  UU',6x,'UU  CC',14x,'11',10x,
     |        '22',14x,'33',/),
     |        2(' ',4x,'NN',6x,'NN  UU',6x,'UU  CC',14x,'11',8x,
     |        '22',8x,'33',6x,'33',/),
     |        2(' ',4x,'NN',6x,'NN  ',10('U'),4x,8('C'),4x,6('1'),4x,
     |        10('2'),4x,6('3'),/),/,
     |        ' ',26x,'WRITTEN BY LAWRENCE KAWANO',///,
     |        ' ','(Press <RETURN> to continue): ',$)

C20--------INPUT INITIALIZATION INFORMATION AND PAUSE------------------------

      DO i  = 1,nrec
C..........READ IN REACTION PARAMETERS.
        iform(i) = int(reacpr(i,2))!Reaction type.
        ii(i)    = int(reacpr(i,3))!Incoming nuclide type.
        jj(i)    = int(reacpr(i,4))!Incoming nuclide type.
        kk(i)    = int(reacpr(i,5))!Outgoing nuclide type.
        ll(i)    = int(reacpr(i,6))!Outgoing nuclide type.
        rev(i)   = reacpr(i,7)     !Reverse reaction coefficient.
        q9(i)    = reacpr(i,8)     !Energy released.
C..........INITIALIZE REACTION RATES.
        f(i)  = 0.                 !Forward rate coeff.
        r(i)  = 0.                 !Reverse rate coeff.
C..........SET RUN OPTIONS TO DEFAULT.
      END DO
      irun       = 1               !Do full run.
      isize      = nnuc            !Use all 26 nuclides.
      jsize      = nrec            !Use all 88 reactions.
C..........SET OUTPUT OPTION TO DEFAULT.
      nout    = 0                  !No output requests.
      outfile = .false.            !Output file not used.
C..........SET VALUES TO DEFAULT.
      cy    = cy0                  !Time step limiting constant on abundances.
      ct    = ct0                  !Time step limiting constant on temperature.
      t9i   = t9i0                 !Initial temperature.
      t9f   = t9f0                 !Final temperature.
      ytmin = ytmin0               !Smallest abundances allowed.
      inc   = inc0                 !Accumulation increment.
      c(1)  = c0(1)                !Variation of gravitational constant.
      c(2)  = c0(2)                !Neutron lifetime.
      c(3)  = c0(3)                !Number of neutrino species.
      cosmo = cosmo0               !Cosmological constant.
      xi(1) = xi0(1)               !Electron degeneracy parameter.
      xi(2) = xi0(2)               !Muon degeneray parameter.
      xi(3) = xi0(3)               !Tauon degeneracy parameter.
c     To force the integration,  SHH
c      xi(1) = 1.E-20
c      xi(2) = 1.E-20
c      xi(3) = 1.E-20
      dt1   = dt0                  !Initial time step.
      eta1  = eta0                 !Baryon-to-photon ratio.
C..........ACCEPT RETURN TO CONTINUE.
      READ (5,*)                  !Pause.

C30--------PRINT MENU AND AWAIT RESPONSE---------------------------------

C..........RETURN FROM LOOPING.
 300  CONTINUE
C..........DISPLAY MENU.
      PRINT 3000
 3000 FORMAT (8(/),
     |        ' ',32x,'MENU SELECTION',/,
     |        ' ',32x,'---- ---------',//,
     |        ' ',24x,'1. HELP',/,
     |        ' ',24x,'2. SET COMPUTATION PARAMETERS',/,
     |        ' ',24x,'3. SET MODEL PARAMETERS',/,
     |        ' ',24x,'4. RUN',/,
     |        ' ',24x,'5. OUTPUT',/,
     |        ' ',24x,'6. EXIT',8(/),
     |        ' ',24x,'Enter selection (1-6): ',$)
C..........READ IN SELECTION NUMBER.
      READ (5,3001) inum
 3001 FORMAT(i1)

C40--------BRANCH TO APPROPRIATE SECTION--------------------------------

      GO TO (410,420,430,440,450,460),inum
      GO TO 460                    !Improper input or <RETURN>.
 410  CONTINUE                     !Help section.
        CALL help
        GO TO 500
 420  CONTINUE                     !Set computation parameters section.
        CALL setcom
        GO TO 500
 430  CONTINUE                     !Set model parameters section.
        CALL setmod
        GO TO 500
 440  CONTINUE                     !Run section.
        itime = 2                  !Time = beginning of run section.
        CALL check                 !Check interface subroutine.
        CALL run
        itime = 9                  !Time = end of run section.
        CALL check                 !Check interface subroutine.
        GO TO 500
 450  CONTINUE                     !Output section.
        CALL output
        GO TO 500
 460  CONTINUE                     !Exit section.
        IF (outfile) THEN 
          close (unit=2,status='keep')   !Close output file.
        ELSE
          CLOSE (unit=2,status='delete') !File not used - dispose.
        END IF
ccccccccccc  CLOSE (unit=1)             !End terminal session.
        itime = 10                 !Time = end of program.
        CALL check                 !Check interface subroutine.
        STOP

C50---------GO BACK TO MENU-----------------------------------------

 500  CONTINUE
      GO TO 300

      END



C========================IDENTIFICATION DIVISION============================

      SUBROUTINE help

C----------LINKAGES.
C     CALLED BY - [program] nuc123
C     CALLS     - none

C----------REMARKS.
C     Displays description and workings of the program.

C----------PARAMETERS.
      PARAMETER (ir=1)             !Input unit number.
      PARAMETER (iw=1)             !Output unit number.


C==========================DECLARATION DIVISION===========================

C----------USER RESPONSE VARIABLES.
      INTEGER inum                 !Selection number.


C===========================PROCEDURE DIVISION=============================

C10--------PRINT HELP SELECTION-------------------------------------------

C..........RETURN FROM LOOPING.
 100  CONTINUE
C..........DISPLAY MENU.
      PRINT 1000
 1000 FORMAT (8(/),
     |        ' ',32x,'HELP SELECTION',/,
     |        ' ',32x,'---- ---------',//,
     |        ' ',24x,'1. INTRODUCTION',/,
     |        ' ',24x,'2. SETTING UP A RUN',/,
     |        ' ',24x,'3. RUNNING THE PROGRAM',/,
     |        ' ',24x,'4. OUTPUT OPTIONS',/,
     |        ' ',24x,'5. GENERAL METHOD OF COMPUTATION',/,
     |        ' ',24x,'6. USING THE INTERFACE SUBROUTINE',/,
     |        ' ',24x,'7. EXIT',7(/),
     |        ' ',24x,'Enter selection (1-7): ',$)
C..........READ IN SELECTION NUMBER.
      READ (5,1001) inum
 1001 FORMAT (i1)

C20--------BRANCH TO APPROPRIATE SECTION---------------------------------

      GO TO (210,220,230,240,250,260,270),inum
      GO TO 270                    !Improper input or <RETURN>.

C21--------INTRODUCTION SECTION------------------------------------------

 210  CONTINUE                     !Setting up a run section.
        PRINT 2100
 2100   FORMAT (/,
     |          ' ',31x,'INTRODUCTION',/,
     |          ' ',31x,'------------',2(/),
     |          ' ','Welcome to the wonderful world of primor',
     |              'dial nucleosynthesis.  NUC123 is a      ',/,
     |          ' ','FORTRAN program designed to provide the ',
     |              'early universe researcher with the tools',/,
     |          ' ','necessary for the investigation of primo',
     |              'rdial nucleosynthesis.  Its menu-driven ',/,
     |          ' ','interface allows the user to first set c',
     |              'omputation parameters (such as the time ',/,
     |          ' ','step) and model parameters (such as the ',
     |              'neutron lifetime and number of neutri-  ',/,
     |          ' ','nos) before doing single runs or multipl',
     |              'e runs (in which desired model parame-  ',/,
     |          ' ','ters are varied over a desired range.)  ',
     |              'After the run, the user can utilize the ',/,
     |          ' ','menu to either produce an output file or',
     |              ' to view the most recent run on the     ',/,
     |          ' ','screen.  The program comes with an empty',
     |              ' subroutine CHECK into which the user   ',/,
     |          ' ','may wish to put additional code to add t',
     |              'o the computation in an original manner.',10(/),
     |          ' ','(Enter <RETURN> to go back to help menu): ',$)
        READ (5,*) 
        GO TO 300

C22--------SET UP RUN SECTION------------------------------------------

 220  CONTINUE                     !Setting up a run section.
        PRINT 2200
 2200   FORMAT (/,
     |          ' ',29x,'SETTING UP A RUN',/,
     |          ' ',29x,'------- -- - ---',2(/),
     |          ' ','I. Setting computation parameters.      ',/,
     |          ' ','   The accuracy of the computation and t',
     |              'he relevant temperature region can be   ',/,
     |          ' ','   set by the following parameters:     ',/,
     |          ' ','    A. Time step limiting constant 1  (d',
     |              'efault value of 0.3)                    ',/,
     |          ' ','    B. Time step limiting constant 2  (d',
     |              'efault value of 0.03)                   ',/,
     |          ' ','    C. Initial time step  (default value',
     |              ' of 10**-4)                             ',/,
     |          ' ','    D. Initial temperature  (default val',
     |              'ue of 10**2)                            ',/,
     |          ' ','       This is the temperature at the be',
     |              'ginning of the run in units of 10**9 K  ',/,
     |          ' ','    E. Final temperature  (default value',
     |              ' of 10**-2)                             ',/,
     |          ' ','       This is the termination temperatu',
     |              're of the run in units of 10**9 K       ',/,
     |          ' ','    F. Smallest abundances allowed  (def',
     |              'ault value of 10**-25)                  ',/,
     |          ' ','       Elemental abundances are not allo',
     |              'wed to drop below this value            ',/,
     |          ' ','    G. # of iterations for each accumula',
     |              'tion  (default value of 30)             ',/,
     |          ' ','       This is the number of iterations ',
     |              'before values are put in an output array',6(/),
     |          ' ','(Enter 1 to continue, <RETURN> to end): ',$)
        READ (5,1001) inum
        IF (inum.eq.1) THEN
          PRINT 2202
 2202     FORMAT (/,
     |            ' ','II. Setting model parameters.           ',/,
     |            ' ','   Default values here give what is know',
     |                'n as the standard model with best guess ',/,
     |            ' ','   figure on the neutron lifetime of 889',
     |                '.541 seconds.  Nonstandard scenarios can',/,
     |            ' ','   be investigated by varying the follow',
     |                'ing parameters:                         ',/,
     |            ' ','    A. The gravitational constant       ',/,
     |            ' ','       (The default value of one here gi',
     |                'ves the usual 6.6720e-8 dyne*cm**2/g**2)',/,
     |            ' ','    B. Neutron life-time  (default value',
     |                ' of 889. seconds)                    ',/,
     |            ' ','    C. Number of neutrino species  (defa',
     |                'ult value of 3 light neutrinos)         ',/,
     |            ' ','    D. Final baryon-to-photon ratio  (se',
     |                't to log(eta) = -9.5)                   ',/,
     |            ' ','    E. Cosmological constant  (default v',
     |                'alue of 0)                              ',/,
     |            ' ','    F. Neutrino degeneracy parameters  (',
     |                'default values all 0)                   ',/,
     |            ' ','       There are 3 separate parameters f',
     |                'or the electron, muon, and tau neutrinos',11(/),
     |            ' ','(Enter <RETURN> to go back to help menu): ',$)
          READ (5,*)
          GO TO 300
        ELSE
          GO TO 300
        END IF !(inum.eq.1)

C23--------RUN PROGRAM SECTION------------------------------------------

 230  CONTINUE                     !Running the program section.
        PRINT 2300
 2300   FORMAT (/,
     |          ' ',28x,'RUNNING THE PROGRAM',/,
     |          ' ',28x,'------- --- -------',2(/),
     |          ' ','I. Setting run speed.                   ',/,
     |          ' ','   The code can be run at 3 different se',
     |              'ttings of speed.  The running of the    ',/,
     |          ' ','   code can be speeded up by reducing th',
     |              'e number of nuclides and reactions.  The',/,
     |          ' ','   complete computation takes into accou',
     |              'nt the following nuclides: n, p, d, t,  ',/,
     |          ' ','   He3, He4, Li6, Li7, Be7, Li8, B8, Be9',    
     |              ',B10, B11, C11, B12, C12, N12, C13, N13,',/,
     |          ' ','   C14, N14, O14, N15, O15, and O16.    ',/,
     |          ' ','   The given CPU percentages and abundan',
     |              'ce variations are with regard to a      ',/,
     |          ' ','   single run with all default parameter',
     |              ' values.                                ',/,
     |          ' ','    A. 26 nuclides, 88 reactions (defaul',
     |              't)                                      ',/,
     |          ' ','       nuclides from n to O16           ',/,
     |          ' ','    B. 18 nuclides, 60 reactions        ',/,
     |          ' ','       nuclides from n to N12           ',/,
     |          ' ','       (63% CPU time, variation = .1%)  ',/,
     |          ' ','    C.  9 nuclides, 25 reactions        ',/,
     |          ' ','       nuclides from n to Be7           ',/,
     |          ' ','       (20% CPU time, variation = .5%)  ',4(/),
     |          ' ','(Enter 1 to continue, <RETURN> to end): ',$)
        READ (5,1001) inum
        IF (inum.eq.1) THEN
          PRINT 2302
 2302     FORMAT (/,
     |            ' ','II. Do single run.                      ',/,
     |            ' ','    A. Interactive.                     ',/,
     |            ' ','       In an interactive session, the us',
     |                'er can readily input the computational  ',/,
     |            ' ','       and model parameters and begin th',
     |                'e computation process.  The run itself  ',/,
     |            ' ','       is commenced when option 2, "GO",',
     |                ' in the "RUN" section is requested.     ',//,
     |            ' ','    B. Batch.                           ',/,
     |            ' ','       To run the program in a batch mod',
     |                'e, it must be altered slightly so that  ',/,
     |            ' ','       the I/O takes place with files in',
     |                'stead of a terminal.  This is done by   ',/,
     |            ' ','       setting different values for the ',
     |                'input and output unit number parameters ',/,
     |            ' ','       "ir" and "iw" and assigning them ',
     |                'to different files in NUC123.  In the   ',/,
     |            ' ','       file assigned the "ir" unit numbe',
     |                'r, one must place the responses to the  ',/,
     |            ' ','       queries of the program.          ',10(/),
     |            ' ','(Enter 1 to continue, <RETURN> to end): ',$)
          READ (5,1001) inum
          IF (inum.eq.1) THEN
            PRINT 2304
 2304       FORMAT (/,
     |              ' ','III. Do multiple runs.                 ',/,
     |              ' ','   A wide range of early universe model',
     |                  's can be covered by doing many runs    ',/,
     |              ' ','   while one or more parameters are var',
     |                  'ied over a range of interest.  The     ',/,
     |              ' ','   parameters that can be varied are th',
     |                  'e following:                           ',/,
     |              ' ','    A. Eta                             ',
     |                  '       - Logrithmic variation          ',/,
     |              ' ','    B. Gravitational constant          ',
     |                  '       - Linear variation              ',/,
     |              ' ','    C. Neutron lifetime                ',
     |                  '       - Linear variation              ',/,
     |              ' ','    D. Number of neutrino species      ',
     |                  '       - Linear variation              ',/,
     |              ' ','    E. Cosmological constant           ',
     |                  '       - Linear variation              ',/,
     |              ' ','    F. Neutrino degeneracy parameters  ',
     |                  '       - Linear variation              ',/,
     |              ' ','        1. Electron neutrino           ',/,
     |              ' ','        2. Muon neutrino               ',/,
     |              ' ','        3. Tauon neutrino              ',/,
     |              ' ','   At most 3 parameters can be varied. ',
     |                  ' The first parameter inputted will be  ',/,
     |              ' ','   will be varied in the outermost loop',
     |                  ' and the third parameter inputted will ',/,
     |              ' ','   be varied in the innermost loop.    ',7(/),
     |              ' ','(Enter <RETURN> to go back to help menu): ',$)
            READ (5,*)
            GO TO 300
          ELSE
            GO TO 300
          END IF !(inum.eq.1)
        ELSE
          GO TO 300
        END IF !(inum.eq.1)

C24--------OUTPUT OPTIONS SECTION----------------------------------

 240  CONTINUE                     !Output options section.
        PRINT 2400
 2400   FORMAT (/,
     |          ' ',30x,'OUTPUT OPTIONS',/,
     |          ' ',30x,'------ -------',2(/),
     |          ' ','I.  Request output file.                ',/,
     |          ' ','   After a run, the user can request the',
     |              ' program to put the resulting numbers   ',/,
     |          ' ','   into an output file.  This can be don',
     |              'e as many times as desired and all the  ',/,
     |          ' ','   information will be put in one new fi',
     |              'le under the name of "NUC123.DAT."  If  ',/,
     |          ' ','   there is no request during the entire',
     |              ' running of the program, this file is   ',/,
     |          ' ','   not created.  If an output file is re',
     |              'quested after a multiple run, only the  ',/,
     |          ' ','   information from the very last run wi',
     |              'll be given.  The output file will give ',/,
     |          ' ','   the computational and model parameter',
     |              's for each run and will contain the     ',/,
     |          ' ','   following information:               ',/,
     |          ' ','    A. Temperatures in decreasing order ',/,
     |          ' ','    B. Abundances for n, p, d, t, He3, H',
     |              'e4, Li6, Li7, Be7, and Li8 & up         ',/,
     |          ' ','       (p and He4 are in mass fraction, ',
     |              'the rest in ratios to the p abundance)  ',/,
     |          ' ','    C. Time, time interval, chemical pot',
     |              'ential of the electron                  ',/,
     |          ' ','    D. Energy densities for photons, ele',
     |              'ctrons, electron neutrinos, and baryons ',/,
     |          ' ','    E. Baryon-to-photon ratio, expansion',
     |              ' rate of the universe                   ',5(/),
     |          ' ','(Enter 1 to continue, <RETURN> to end): ',$)
        READ (5,1001) inum
        IF (inum.eq.1) THEN
          PRINT 2402
 2402     FORMAT (/,
     |            ' ','II.  Request output on screen.         ',/,
     |            ' ','   Instead of waiting to print out an o',
     |                'utput file, the user can immediately   ',/,
     |            ' ','   access the results of the latest run',
     |                ' by requesting the output on the       ',/,
     |            ' ','   screen.  There are four screens on e',
     |                'ach of which are displayed the         ',/,
     |            ' ','   computational and model parameters a',
     |                'nd the temperature:                    ',/,
     |            ' ','    A. Abundances for d, t, He3, He4, a',
     |                'nd Li7                                 ',/,
     |            ' ','       (He4 in mass fraction, rest as a',
     |                ' ratio with the p abundance)           ',/,
     |            ' ','    B. Abundances for n, p, Li6, Be7, a',
     |                'nd Li8 & up                            ',/,
     |            ' ','       (p in mass fraction, rest as a r',
     |                'atio with the p abundance)             ',/,
     |            ' ','    C. Energy densities for photons, el',
     |                'ectrons, electron neutrinos, & baryons ',/,
     |            ' ','    D. Time, time interval, chemical po',
     |                'tential of the electron,               ',/,
     |            ' ','       baryon-to-photon ratio, and expa',
     |                'nsion rate of the universe             ',11(/),
     |            ' ','(Enter <RETURN> to go back to help menu): ',$)
          READ (5,*)
          GO TO 300
        ELSE
          GO TO 300
        END IF !(inum.eq.1)

C25--------METHOD OF COMPUTATION SECTION--------------------------------

 250  CONTINUE                     !General method of computation section.
        PRINT 2500
 2500   FORMAT (/,
     |          ' ',22x,'GENERAL METHOD OF COMPUTATION',/,
     |          ' ',22x,'------- ------ -- -----------',2(/),
     |          ' ','I. Time evolution algorithm.            ',/,
     |          ' ','   The program utilizes a 2-point Runge-',
     |              'Kutta scheme (located in subroutine     ',/,
     |          ' ','   DRIVER) to time-evolve the temperatur',
     |              'e, the quantity hv (the ratio of the    ',/,
     |          ' ','   baryon density to T**3), the chemical',
     |              ' potential of the electron, and the     ',/,
     |          ' ','   nuclide abundances.  In the 2-point R',
     |              'unge-Kutta routine, a variable v at time',/,
     |          ' ','   t0 (= v0) is evolved to a time t1 by ',
     |              'adding to v0 the average of the         ',/,
     |          ' ','   derivatives evaluated at t0 and at t1',
     |              ' multiplied by dt:                      ',/,
     |          ' ','       v1 = v0 + 0.5(dvdt(t0)+dvdt(t1)) ',/,
     |          ' ','   where dvdt(t1) is gotten by first fin',
     |              'ding v1'' = v0 + dvdt(t0).  The         ',/,
     |          ' ','   derivatives of the nuclide abundances',
     |              ' are first computed and these are used  ',/,
     |          ' ','   to find the derivatives of t9, hv, an',
     |              'd phie (this is done in subroutine      ',/,
     |          ' ','   DERIVS).  To compute the time derivat',
     |              'ives of the nuclide abundances, a matrix',/,
     |          ' ','   equation is set up (in subroutine SOL',
     |              ') and is solved (in subroutine EQSLIN)  ',/,
     |          ' ','   by gaussian elimination utilizing imp',
     |              'licit differentiation.                  ',6(/),
     |          ' ','(Enter 1 to continue, <RETURN> to end): ',$)
        READ (5,1001) inum
        IF (inum.eq.1) THEN
          PRINT 2502
 2502     FORMAT (/
     |            ' ','II. Hierarchy of Subroutines.   ',/,
     |            ' ','    NUC123                       ',
     |                '     Main program (main menu)    ',/,
     |            ' ','        HELP                     ',
     |                '     Help option                 ',/,
     |            ' ','        SETCOM                   ',
     |                '     Set computational parameters',/,
     |            ' ','        SETMOD                   ',
     |                '     Set model parameters        ',/,
     |            ' ','        RUN                      ',
     |                '     Run computation code        ',/,
     |            ' ','            DRIVER               ',
     |                '     Main routine (Runge-Kutta loop)    ',/,
     |            ' ','                START            ',
     |                '     Initialization routine      ',/,
     |            ' ','                    RATE0        ',
     |                '     Computes weak decay rates   ',/,
     |            ' ','                DERIVS           ',
     |                '     Computes time derivatives   ',/,
     |            ' ','                    THERM        ',
     |                '     Computes energy densities   ',/,
     |            ' ','                        BESSEL   ',
     |                '     Gives functions of Kn       ',/,
     |            ' ','                            KNUX ',
     |                '     Computes modified Bessel fcn Kn    ',/,
     |            ' ','                        NUDENS   ',
     |                '     Computes neutrino energy density   ',/,
     |            ' ','                    RATE1-4      ',
     |                '     Computes rates for reactions',/,
     |            ' ','                    SOL          ',
     |                '     Builds A matrix for eqn dy/dt = Ay ',/,
     |            ' ','                        EQSLIN   ',
     |                '     Solves dy/dt=Ay by gaussian elim   ',/,
     |            ' ','                ACCUM            ',
     |                '     Output accumulator          ',/,
     |            ' ','        OUTPUT                   ',
     |                '     Allows user to output result',4(/),
     |            ' ','(Enter <RETURN> to go back to help menu): ',$)
          READ (5,*)
          GO TO 300
        ELSE
          GO TO 300
        END IF !(inum.eq.1)

C26--------USING INTERFACE SUBROUTINE SECTION.

 260  CONTINUE                     !Using the interface subroutine section.
        PRINT 2600
 2600   FORMAT (/,
     |          ' ',22x,'USING THE INTERFACE SUBROUTINE',/,
     |          ' ',22x,'----- --- --------- ----------',2(/),
     |          ' ','I. Purpose.                             ',/,
     |          ' ','   The interface subroutine CHECK is des',
     |              'igned to be an outlet of the program    ',/,
     |          ' ','   into which alterations can be easily ',
     |              'plugged.  Programs are normally modified',/,
     |          ' ','   by searching through the program, ide',
     |              'ntifying the appropriate areas for      ',/,
     |          ' ','   alterations, and interspersing new co',
     |              'mmands while deleting some old ones.    ',/,
     |          ' ','   This process can get tricky unless on',
     |              'e actively documents the alterations:   ',/,
     |          ' ','   one might lose track of all of the mo',
     |              'difications and deletions.  Thus, it is ',/,
     |          ' ','   worthwhile to put most if not all of ',
     |              'the necessary changes into one          ',/,
     |          ' ','   subroutine which is to be called from',
     |              ' strategic locations in the main        ',/,
     |          ' ','   program.  Furthermore, by putting cha',
     |              'nges into one small subroutine, one need',/,
     |          ' ','   only to compile the subroutine CHECK ',
     |              'each time instead of the entire nucleo- ',/,
     |          ' ','   synthesis code.                      ',8(/),
     |          ' ','(Enter 1 to continue, <RETURN> to end): ',$)
        READ (5,1001) inum
        IF (inum.eq.1) THEN
          PRINT 2602
 2602     FORMAT (/,
     |            ' ','II. Description.                        ',/,
     |            ' ','   Subroutine CHECK is an empty subrouti',
     |                'ne with a large COMMON area, giving the ',/,
     |            ' ','   user ready access to all of the impor',
     |                'tant variables in the computations.  The',/,
     |            ' ','   routine is called from various locati',
     |                'ons in the main program and the location',/,
     |            ' ','   spot in the program is labeled by the',
     |                ' flag "itime".  The set call locations  ',/,
     |            ' ','   are given below:                     ',/,
     |            ' ','    A. itime = 1 (NUC123, very beginning',
     |                ' of program run)                        ',/,
     |            ' ','       (appropriate for opening files, i',
     |                'nitializing variables)                  ',/,
     |            ' ','    B. itime = 2 (NUC123, right before g',
     |                'oing into the RUN section)              ',/,
     |            ' ','    C. itime = 3 (RUN, right before goin',
     |                'g into DRIVER to do the computations)   ',/,
     |            ' ','    D. itime = 4 (DRIVER, in 1st R-K loo',
     |                'p after computing derivatives in DERIVS)',/,
     |            ' ','    E. itime = 7 (DRIVER, in 2nd R-K loo',
     |                'p after computing derivatives in DERIVS)',/,
     |            ' ','    F. itime = 8 (RUN, right after comin',
     |                'g back from DRIVER)                     ',/,
     |            ' ','    G. itime = 9 (NUC123, right after co',
     |                'ming back from the RUN section)         ',/,
     |            ' ','    H. itime =10 (NUC123, very end of pr',
     |                'ogram run)                              ',/,
     |            ' ','       (appropriate for closing files)  ',/,
     |            ' ','   The difference between the (2,9) pair',
     |                'ing and the (3,8) pairing is that for a ',/,
     |            ' ','   multiple run, the (3,8) pairing would',
     |                ' be called before and after every run   ',/,
     |            ' ','   but the (2,9) pairing would be called',
     |                ' before and after the entire sequence.  ',4(/),
     |            ' ','(Enter 1 to continue, <RETURN> to end): ',$)
          READ (5,1001) inum
          IF (inum.eq.1) THEN
            PRINT 2604
 2604       FORMAT (/,
     |              ' ','III. Implementation.                   ',/,
     |              ' ','   The additional program statements ar',
     |                  'e needed in the subroutine CHECK.  If a',/,
     |              ' ','   particular command is to be executed',
     |                  ' when the computer is at a certain     ',/,
     |              ' ','   location in the program -- say label',
     |                  'ed by itime = 8 -- then in CHECK, one  ',/,
     |              ' ','   must place the command under the sta',
     |                  'tement, IF (itime.eq.8)....  The user  ',/,
     |              ' ','   is at leisure to place his own locat',
     |                  'ion indicators (5,6) and CALL CHECK    ',/,
     |              ' ','   statements anywhere in the program a',
     |                  's long as there is a COMMON /check/    ',/,
     |              ' ','   statement in the particular subrouti',
     |                  'ne to carry the value of itime along.  ',15(/),
     |              ' ','(Enter <RETURN> to go back to help menu): ',$)
            READ (5,*)
            GO TO 300
          ELSE
            GO TO 300
          END IF !(inum.eq.1)
        ELSE
          GO TO 300
        END IF !(inum.eq.1)

C27--------EXIT SECTION----------------------------------------------

 270  CONTINUE                     !Exit section.
        RETURN

C30--------GO BACK TO MAIN MENU-------------------------------------------

 300  CONTINUE
      GO TO 100

      END



C========================IDENTIFICATION DIVISION==========================

      SUBROUTINE setcom

C----------LINKAGES.
C     CALLED BY - [program] nuc123
C     CALLS     - none

C----------REMARKS.
C     Allows resetting of computation parameters.

C----------PARAMETERS.
      PARAMETER (ir=1)             !Input unit number.
      PARAMETER (iw=1)             !Output unit number.

C----------COMMON AREAS.
      COMMON /compr0/ cy0,ct0,t9i0,t9f0,ytmin0,inc0  !Default comp parameters.
      COMMON /compr/  cy,ct,t9i,t9f,ytmin,inc        !Computation parameters.
      COMMON /varpr0/ dt0,eta0                      !Default variationl params.
      COMMON /varpr/  dt1,eta1                       !Variational parameters.


C==========================DECLARATION DIVISION=============================

C----------DEFAULT COMPUTATION PARAMETERS.
      REAL    cy0                  !Default cy.
      REAL    ct0                  !Default ct.
      REAL    t9i0                 !Default t9i.
      REAL    t9f0                 !Default t9f.
      REAL    ytmin0               !Default ytmin.
      INTEGER inc0                 !Default accumulation increment.

C----------COMPUTATION PARAMETERS.
      REAL    cy                   !Time step limiting constant on abundances.
      REAL    ct                   !Time step limiting constant on temperature.
      REAL    t9i                  !Initial temperature (in 10**9 K).
      REAL    t9f                  !Final temperature (in 10**9 K).
      REAL    ytmin                !Smallest abundances allowed.
      INTEGER inc                  !Accumulation increment.

C----------DEFAULT VARIATIONAL  PARAMETERS.
      REAL    dt0                  !Default initial dt.

C----------VARIATIONAL  PARAMETERS.
      REAL    dt1                  !Initial time step.

C----------LOCAL VARIABLES.
      INTEGER inum                 !Selection number.


C===========================PROCEDURE DIVISION===========================

C10--------PRINT RESET SELECTION AND AWAIT RESPONSE------------------------

C..........RETURN FROM LOOPING.
 100  CONTINUE
C..........DISPLAY RESET SELECTIONS.
      PRINT 1000, cy,ct,dt1,t9i,t9f,ytmin,float(inc)
 1000 FORMAT (8(/),
     |        ' ',21x,'SET COMPUTATION PARAMETERS SELECTION',/,
     |        ' ',21x,'--- ----------- ---------- ---------',//,
     |        ' ',10x,' 1. CHANGE TIME-STEP LIMITING CONSTANT 1  FROM ',
     |            f5.3,/,
     |        ' ',10x,' 2. CHANGE TIME-STEP LIMITING CONSTANT 2  FROM ',
     |            f5.3,/,
     |        ' ',10x,' 3. CHANGE INITIAL TIME-STEP              FROM ',
     |            1pe8.2,' SECONDS',/,
     |        ' ',10x,' 4. CHANGE INITIAL TEMPERATURE            FROM ',
     |            1pe8.2,' (10**9 K)',/,
     |        ' ',10x,' 5. CHANGE FINAL TEMPERATURE              FROM ',
     |            1pe8.2,' (10**9 K)',/,
     |        ' ',10x,' 6. CHANGE SMALLEST ABUNDANCES ALLOWED    FROM ',
     |            1pe8.2,/,
     |        ' ',10x,' 7. CHANGE ACCUMULATION INCREMENT         FROM ',
     |            1pe8.2,' ITERATIONS',/,
     |        ' ',10x,' 8. RESET ALL TO DEFAULT VALUES',/,
     |        ' ',10x,' 9. EXIT',5(/),
     |        ' ',10x,'Enter selection (1-9): ',$)
C..........READ IN SELECTION NUMBER.
      READ (5,1001) inum
1001  FORMAT (i1)

C20--------BRANCH TO APPROPRIATE SECTION--------------------------------

      GO TO (210,220,230,240,250,260,270,280,300),inum
      GO TO 300                    !Improper input or <RETURN>.
 210  CONTINUE                     !Change time step limiting const 1 section.
        PRINT 2100
 2100   FORMAT (' ','Enter value for time step limiting constant 1: ',$)
        READ (5,*) cy
 2101   FORMAT (f5.3)
        GO TO 400
 220  CONTINUE                     !Change time step limiting const 2 section.
        PRINT 2200
 2200   FORMAT (' ','Enter value for time step limiting constant 2: ',$)
        READ (5,*) ct
        GO TO 400
 230  CONTINUE                     !Change initial time step section.
        PRINT 2300
 2300   FORMAT (' ','Enter value for initial time step: ',$)
        READ (5,*) dt1
        GO TO 400
 240  CONTINUE                     !Change initial temperature section.
        PRINT 2400
 2400   FORMAT (' ','Enter value for initial temperature: ',$)
        READ (5,*) t9i
        GO TO 400
 250  CONTINUE                     !Change final temperature section.
        PRINT 2500
 2500   FORMAT (' ','Enter value for final temperature: ',$)
        READ (5,*) t9f
        GO TO 400
 260  CONTINUE                     !Change smallest abundances allowed section.
        PRINT 2600
 2600   FORMAT (' ','Enter value for smallest abundances allowed: ',$)
        READ (5,*) ytmin
        GO TO 400
 270  CONTINUE                     !Change accumulation increment section.
        PRINT 2700
 2700   FORMAT (' ','Enter value for accumulation increment: ',$)
        READ (5,*) inc
        GO TO 400
 280  CONTINUE                     !Reset all to default values section.
        cy    = cy0                !Time step limiting constant on abundances.
        ct    = ct0                !Time step limiting constant on temperature.
        dt1   = dt0                !Time step.
        t9i   = t9i0               !Initial temperature.
        t9f   = t9f0               !Final temperature.
        ytmin = ytmin0             !Smallest abundances allowed.
        inc   = inc0               !Accumulation increment.
        PRINT 2800
 2800   FORMAT (' ','All values reset to default - Press <RETURN> ',
     |              'to continue: ',$)
        READ (5,*)
        GO TO 400
 300  CONTINUE                     !Exit section.
        RETURN

C40--------GO BACK TO MENU----------------------------------------

 400  CONTINUE
      GO TO 100

      END



C========================IDENTIFICATION DIVISION========================

      SUBROUTINE setmod

C----------LINKAGES.
C     CALLED BY - [program] nuc123
C     CALLS     - none

C----------REMARKS.
C     Allows resetting of model parameters.

C----------PARAMETERS.
      PARAMETER (ir=1)             !Input unit number.
      PARAMETER (iw=1)             !Output unit number.

C----------COMMON AREAS.
      COMMON /modpr0/ c0,cosmo0,xi0                  !Default model parameters.
      COMMON /modpr/  g,tau,xnu,c,cosmo,xi           !Model parameters.
      COMMON /varpr0/ dt0,eta0                      !Default variationl params.
      COMMON /varpr/  dt1,eta1                       !Variational parameters.


C==========================DECLARATION DIVISION==========================

C----------DEFAULT MODEL PARAMETERS.
      REAL    c0(3)                !Default c.
      REAL    cosmo0               !Default cosmological constant.
      REAL    xi0(3)               !Default neutrino degeneracy parameters.

C----------EARLY UNIVERSE MODEL PARAMETERS.
      REAL    c(3)                 !c(1) is variation of gravitational constant.
     |                             !c(2) is neutron lifetime (sec).
     |                             !c(3) is number of neutrino species.
      REAL    cosmo                !Cosmological constant.
      REAL    xi(3)                !Neutrino degeneracy parameters.

C----------DEFAULT VARIATIONAL PARAMETERS.
      REAL    eta0                 !Default eta.

C----------VARIATIONAL PARAMETERS.
      REAL    eta1                 !Intial baryon-to-photon ratio.

C----------USER RESPONSE VARIABLES.
      INTEGER inum                 !Selection number.


C===========================PROCEDURE DIVISION============================

C10--------PRINT RESET SELECTION AND AWAIT RESPONSE------------------------

C..........RETURN FROM LOOPING.
 100  CONTINUE
C..........DISPLAY RESET SELECTIONS.
      PRINT 1000, c(1),c(2),c(3),eta1,cosmo,xi(1),xi(2),xi(3)
 1000 FORMAT (8(/),
     |        ' ',24x,'SET MODEL PARAMETERS SELECTION',/,
     |        ' ',24x,'--- ----- ---------- ---------',//,
     |        ' ',10x,' 1. CHANGE GRAVITATIONAL CONSTANT         FROM ',
     |            1pe10.3,/,
     |        ' ',10x,' 2. CHANGE NEUTRON LIFETIME               FROM ',
     |            1pe10.3,' SECONDS',/,
     |        ' ',10x,' 3. CHANGE NUMBER OF NEUTRINO SPECIES     FROM ',
     |            1pe10.3,/,
     |        ' ',10x,' 4. CHANGE FINAL BARYON-TO-PHOTON RATIO   FROM ',
     |            1pe10.3,/,
     |        ' ',10x,' 5. CHANGE COSMOLOGICAL CONSTANT          FROM ',
     |            1pe10.3,/,
     |        ' ',10x,' 6. CHANGE XI-ELECTRON                    FROM ',
     |            1pe10.3,/,
     |        ' ',10x,' 7. CHANGE XI-MUON                        FROM ',
     |            1pe10.3,/,
     |        ' ',10x,' 8. CHANGE XI-TAUON                       FROM ',
     |            1pe10.3,/,
     |        ' ',10x,' 9. RESET ALL TO DEFAULT VALUES',/,
     |        ' ',10x,'10. EXIT',4(/),
     |        ' ',10x,' Enter selection (1-10): ',$)
C..........READ IN SELECTION NUMBER.
      READ (5,1001) inum
 1001 FORMAT (i2)

C20--------BRANCH TO APPROPRIATE SECTION---------------------------------

      GO TO (210,220,230,240,250,260,270,280,290,300),inum
      GO TO 300                    !Improper input or <RETURN>.
 210  CONTINUE                     !Change gravitational constant section.
        PRINT 2100
 2100   FORMAT (' ','Enter value for variation of gravitational ',
     |              'constant: ',$)
        READ (5,*) c(1)
        GO TO 400
 220  CONTINUE                     !Change neutron lifetime section.
        print 2200
 2200   FORMAT (' ','Enter value for neutron lifetime (sec): ',$)
        READ (5,*) c(2)
        GO TO 400
 230  CONTINUE                     !Change number of neutrino species section.
        print 2300
 2300   FORMAT (' ','Enter value for number of neutrino species: ',$)
        READ (5,*) c(3)
        GO TO 400
 240  CONTINUE                     !Change baryon-to-photon ratio section.
        print 2400
 2400   FORMAT (' ','Enter value for baryon-to-photon ratio: ',$)
        READ (5,*) eta1
        GO TO 400
 250  CONTINUE                     !Change cosmological constant section.
        print 2500
 2500   FORMAT (' ','Enter value for cosmological constant: ',$)
        READ (5,*) cosmo
        GO TO 400
 260  CONTINUE                     !Change neutrino degeneracy section.
        print 2600
 2600   FORMAT (' ','Enter value for xi electron: ',$)
        READ (5,*) xi(1)
        GO TO 400
 270  CONTINUE                     !Change neutrino degeneracy section.
        print 2700
 2700   FORMAT (' ','Enter value for xi muon: ',$)
        READ (5,*) xi(2)
        GO TO 400
 280  CONTINUE                     !Change neutrino degeneracy section.
        print 2800
 2800   FORMAT (' ','Enter value for xi tauon: ',$)
        READ (5,*) xi(3)
        IF ((xi(3).ne.0.).and.(c(3).lt.3.)) THEN
          c(3) = 3.
          print 2802
 2802     FORMAT (' ','Number of neutrinos set to 3')
          print 2804
 2804     FORMAT (' ','Press <RETURN> to continue: ',$)
          READ (5,*)
        END IF
        GO TO 400
 290  CONTINUE                     !Reset all to default values section.
        c(1)   = c0(1)
        c(2)   = c0(2)
        c(3)   = c0(3)
        cosmo  = cosmo0
        xi(1)  = xi0(1)
        xi(2)  = xi0(2)
        xi(3)  = xi0(3)
        eta1   = eta0
        print 2900
 2900   FORMAT (' ','All values reset to default - Press <RETURN> ',
     |              'to continue: ',$)
        READ (5,*)
        GO TO 400
 300  CONTINUE                     !Exit section.
        RETURN

C40--------GO BACK TO MENU---------------------------------------------

 400  CONTINUE
      GO TO 100

      END



C========================IDENTIFICATION DIVISION=========================

      SUBROUTINE run

C----------LINKAGES.
C     CALLED BY - [program] nuc123
C     CALLS     - [subroutine] driver

C----------REMARKS.
C     Activates computation routine.

C----------PARAMETERS.
      PARAMETER (ir=1)             !Input unit number.
      PARAMETER (iw=1)             !Output unit number.
      PARAMETER (nrec=88)          !Number of nuclear reactions.
      PARAMETER (lrec=64)          !Total # of nuclear reactions for irun = 2.
      PARAMETER (krec=34)          !Total # of nuclear reactions for irun = 3.
      PARAMETER (nnuc=26)          !Number of nuclides in calculation.
      PARAMETER (lnuc=18)          !Total # of nuclides for irun = 2.
      PARAMETER (knuc=9)           !Total # of nuclides for irun = 3.

C----------COMMON AREAS.
      COMMON /modpr/  g,tau,xnu,c,cosmo,xi           !Model parameters.
      COMMON /varpr/  dt1,eta1                       !Variational parameters.
      COMMON /check1/  itime                          !Computation location.
      COMMON /runopt/ irun,isize,jsize               !Run options.


C==========================DECLARATION DIVISION===========================

C----------MODEL PARAMETERS.
      REAL    eta1                 !Baryon-to-photon ratio.
      REAL    c(3)                !c(1) is variation of gravitational constant.
     |                             !c(2) is neutron lifetime (sec).
     |                             !c(3) is number of neutrino species.
      REAL    cosmo                !Cosmological constant.
      REAL    xi(3)                !Neutrino degeneracy parameters.

C----------RUN OPTION.
      INTEGER irun                 !Run network size.
      INTEGER isize                !Number of nuclides in computation.
      INTEGER jsize                !Number of reactions in computation.

C----------USER INTERACTION VARIABLES.
      REAL    rnumb1               !Run parameter for outer loop.
      REAL    rnumb2               !Run parameter for middle loop.
      REAL    rnumb3               !Run parameter for inner loop.
      REAL    rnum1(3)             !Run parameter starting value.
      REAL    rnum2(3)             !Run parameter end value.
      REAL    rnum3(3)             !Run parameter increment.
      INTEGER inumb                !Selection number.
      INTEGER inum(3)              !Selection number.
      INTEGER jnum                 !Number of loopings to be done.
      INTEGER knum                 !Number of loopings rejected.
      INTEGER lnumb1               !Run parameter for outer loop.
      INTEGER lnumb2               !Run parameter for middle loop.
      INTEGER lnumb3               !Run parameter for inner loop.
      INTEGER lnum(3)              !Run parameter end value.
      INTEGER lchose               !User response (alphanumeric).

C----------FLAG AND LABELS.
      INTEGER itime                !Computation location.
      CHARACTER*22 vtype(8)        !Label for quantities being varied.

C----------EQUIVALENCE VARIABLE.
      REAL    qvary(7)             !Array set equal to c, cosmo, and xi.

C----------EQUIVALENCE STATEMENTS.
      EQUIVALENCE (qvary(1),c(1)), (qvary(4),cosmo), (qvary(5),xi(1))


C==============================DATA DIVISION=========================

C----------LABELS FOR QUANTITIES BEING VARIED.
      DATA vtype /'baryon/photon ratio   ',
     |            'gravitational constant',
     |            'neutron lifetime      ',
     |            '# of neutrino species ',
     |            'cosmological constant ',
     |            'xi-electron           ',
     |            'xi-muon               ',
     |            'xi-tauon              '/


C===========================PROCEDURE DIVISION============================

C10--------PRINT RUN SELECTION AND AWAIT RESPONSE------------------------

C..........RETURN FROM LOOPING.
 100  CONTINUE
C..........DISPLAY RUN SELECTIONS.
      print 1000
 1000 FORMAT (8(/),
     |        ' ',32x,'RUN SELECTION',/,
     |        ' ',32x,'--- ---------',//,
     |        ' ',27x,' 1. SET RUN NETWORK',/,
     |        ' ',27x,' 2. GO',/,
     |        ' ',27x,' 3. DO MULTIPLE RUNS',/,
     |        ' ',27x,' 4. EXIT',10(/),
     |        ' ',27x,' Enter selection (1-4): ',$)
C..........READ IN SELECTION NUMBER.
      READ (5,1001) inumb
 1001 FORMAT (i1)

C20--------BRANCH TO APPROPRIATE SECTION------------------------------------

      GO TO (210,220,230,240),inumb
      GO TO 240                    !Improper input or <RETURN>.

C21--------SET RUN NETWORK SECTION--------------------------------------

 210  CONTINUE
        print 2100
 2100   FORMAT (' ','Enter network size (1-26 nuclides (default); ',
     |              '2-18; 3-9): ',$)
        READ (5,*) inumb          !Read in selection number.
        IF ((inumb.ne.1).and.(inumb.ne.2).and.(inumb.ne.3)) inumb = 1 !Default.
        IF (inumb.ne.irun) THEN    !Run network changed from previously.
          irun = inumb             !Run network size selection.
        END IF
        IF (irun.eq.1) THEN        !Maximal network size.
          isize = nnuc
          jsize = nrec
        ELSE
          IF (irun.eq.2) THEN      !Abridged network size.
            isize = lnuc
            jsize = lrec
          ELSE
            IF (irun.eq.3) THEN    !Minimal network size.
              isize = knuc
              jsize = krec
            END IF
          END IF
        END IF !(irun.eq.1) 
        print 2104, irun
 2104   FORMAT (' ','Run network set to ',i1,' - Press <RETURN> ',
     |              'to continue: ',$)
        READ (5,*)
        GO TO 300

C22--------GO SECTION--------------------------------------

 220  CONTINUE
        print 2200
 2200   FORMAT (' ','Begin computation run....')
        itime = 3
        CALL check                 !Call interface subr before computation.
        CALL driver                !Do nucleosynthesis computation.
        itime = 8
        CALL check                 !Call interface subr after computation.
        print 2202
 2202   FORMAT (' ','Computation completed - Press <RETURN> to ',
     |              'continue: ',$)
        READ (5,*)
        GO TO 300

C23--------DO MULTIPLE RUNS SECTION----------------------------------------

C..........GET NUMBER OF LOOPINGS.
 230  CONTINUE
        print 2300
 2300   FORMAT (' ','Enter the number of loopings to be done (1 ',
     |              '(default); 2; 3): ',$)
        READ (5,*) jnum           !Read in number of loopings to be done.
        IF ((jnum.ne.1).and.(jnum.ne.2).and.(jnum.ne.3)) THEN
          jnum = 1                 !Default number of loopings.
        END IF
        knum = 0.                  !No loopings rejected for now.
        DO i = 1,3
          IF (i.gt.jnum) THEN
            rnum1(i) = 0.          !Initialize initial parameter.
            rnum2(i) = 0.          !Initialize terminal parameter.
            rnum3(i) = 1.          !Initialize incremental parameter.
            inum(i)  = 0           !Initialize selection number.
          ELSE
C..........OBTAIN QUANTITY TO VARY.
            print 2302
 2302       FORMAT (8(/),
     |              ' ',30x,'QUANTITY TO VARY',/,
     |              ' ',30x,'-------- -- ----',//,
     |              ' ',25x,' 1.  ETA     (LOGRITHMIC VARIATION)',/,
     |              ' ',25x,' 2.  G           (LINEAR VARIATION)',/,
     |              ' ',25x,' 3.  TAU         (LINEAR VARIATION)',/,
     |              ' ',25x,' 4.  # NEUTRINOS (LINEAR VARIATION)',/,
     |              ' ',25x,' 5.  LAMBDA      (LINEAR VARIATION)',/,
     |              ' ',25x,' 6.  XI-ELECTRON (LINEAR VARIATION)',/,
     |              ' ',25x,' 7.  XI-MUON     (LINEAR VARIATION)',/,
     |              ' ',25x,' 8.  XI-TAUON    (LINEAR VARIATION)',/,
     |              ' ',25x,' 9.  NO SELECTION',5(/),
     |              ' ',25x,' Enter selection (1-9): ',$)
            READ (5,1001) inum(i)
            IF ((inum(i).lt.1).or.(inum(i).gt.8)) THEN  !No selection made.
              print 2304
 2304         FORMAT (' ','No selection made - Reduce number of ',
     |                    'loopings by one',/,
     |                ' ','Press <RETURN> to continue: ',$)
              READ (5,*)
              knum = knum + 1   !Step up number of loopings rejected.
              rnum1(i) = 0.     !Initialize initial parameter.
              rnum2(i) = 0.     !Initialize terminal parameter.
              rnum3(i) = 1.     !Initialize incremental parameter.
              inum(i)  = 0      !Initialize selection number.
            ELSE !((inum(i).ge.1).and.(inum(i).le.8)) 
C..........INPUT RUN SPECIFICATIONS.
 231          CONTINUE
              print 2306
 2306         FORMAT (' ','Enter minimum value: ',$)
              READ (5,*) rnum1(i)  !Read in starting value.
              print 2308
 2308         FORMAT (' ','Enter maximum value: ',$)
              READ (5,*) rnum2(i)  !Read in terminating value.
 232          CONTINUE
              print 2310
 2310         FORMAT (' ','Enter increment: ',$)
              READ (5,*) rnum3(i)  !Read in incremental value.
              IF (rnum3(i).eq.0.) THEN !Trouble with 0 division later on.
                print 2312
 2312           FORMAT (' ','Zero increment not allowed: trouble with ',
     |                      'dividing by zero')
                GO TO 232
              END IF 
              print 2314, rnum1(i), rnum2(i), rnum3(i) !Display input info.
 2314         FORMAT (' ','Run from ',1pe12.5,' to ',1pe12.5,
     |                    ' in increments of ',1pe12.5)
              print 2316
 2316         FORMAT (' ','Confirm these values (1=Y or 0=N): ',$)
              READ (5,*) lchose                       !Get confirmation.
 2301         FORMAT (a1)
              IF (lchose.eq.0) GO TO 231
            END IF !((inum(i).lt.1).or.(inum(i).gt.8)) 
          END IF !(i.gt.jnum) 
        END DO !i = 1,3
        jnum = jnum-knum           !Number of valid loopings.
        IF (jnum.ne.0) THEN        !Run requested.
C..........PRINTOUT QUANTITY TO VARY, RUN SPECIFICATIONS.
          DO l = 1,jnum+knum       !Check all loopings.
            IF (inum(l).ne.0) THEN !Proper selection was made.
              print 2318, vtype(inum(l)),rnum1(l),     !Display run params.
     |                        rnum2(l), rnum3(l)
 2318         FORMAT (' ','Run ',a22,/,
     |                    '    from ',1pe12.5,' to ',1pe12.5,
     |                    ' in increments of ',1pe12.5)
C..........GET LOGS OF eta VALUES FOR LOGRITHMIC INCREMENTATION.
              IF (inum(l).eq.1) THEN  !Work with exponents for eta increments.
                rnum1(l) = log10(rnum1(l))
                rnum2(l) = log10(rnum2(l))
              END IF
            END IF
          END DO
C..........COMPUTE NUMBER OF RUNS FOR EACH LOOPING.
          DO l = 1,3
            lnum(l) = nint((rnum2(l)-rnum1(l)+rnum3(l))/rnum3(l))
          END DO
C..........DO MULTIPLE RUNS.
          print 2200          !Inform user of beginning of computation.
          DO lnumb1 = 0,lnum(1)-1  !Outer loop.
            rnumb1 = rnum1(1)+float(lnumb1)*rnum3(1)  !Value of param for run.
            IF ((inum(1).ge.1).and.(inum(1).le.8)) THEN
              IF (inum(1).eq.1) THEN
                eta1 = 10**rnumb1    !Vary baryon-to-photon ratio.
              ELSE 
                qvary(inum(1)-1) = rnumb1  !Vary other quantities.
              END IF
            END IF
            DO lnumb2 = 0,lnum(2)-1  !Middle loop.
              rnumb2 = rnum1(2)+float(lnumb2)*rnum3(2) !Value of param for run.
              IF ((inum(2).ge.1).and.(inum(2).le.8)) THEN
                IF (inum(2).eq.1) THEN
                  eta1 = 10**rnumb2  !Vary baryon-to-photon ratio.
                ELSE 
                  qvary(inum(2)-1) = rnumb2  !Vary other quantities.
                END IF
              END IF
              DO lnumb3 = 0,lnum(3)-1  !Inner loop.
                rnumb3 = rnum1(3)+float(lnumb3)*rnum3(3)  !Value of parameter.
                IF ((inum(3).ge.1).and.(inum(3).le.8)) THEN
                  IF (inum(3).eq.1) THEN
                    eta1 = 10**rnumb3  !Vary baryon-to-photon ratio.
                  ELSE 
                    qvary(inum(3)-1) = rnumb3  !Vary other quantities.
                  END IF
                END IF
                itime = 3
                CALL check         !Check interface subr before computation.
                CALL driver        !Do nucleosynthesis computation.
                itime = 8
                CALL check       !Check interface subroutine after computation.
              END DO !lnumb3 = 0,lnum(3)-1  
            END DO !lnumb2 = 0,lnum(2)-1  
          END DO !lnumb1 = 0,lnum(1)-1  
          print 2202          !Inform user of completion of computation.
        ELSE !(jnum.eq.0)
          print 2320
 2320     FORMAT (' ','No selection made - ',
     |                'Press <RETURN> to continue: ',$)   
        END IF !(jnum.ne.0) 
        READ (5,*)
        GO TO 300

C24--------EXIT SECTION--------------------------------------------------

 240  CONTINUE
        RETURN

C30--------GO BACK TO MENU-----------------------------------------------

 300  CONTINUE
      GO TO 100

       END



C========================IDENTIFICATION DIVISION==========================

      SUBROUTINE output

C----------LINKAGES.
C     CALLED BY - [program] nuc123
C     CALLS     - none

C----------REMARKS.
C     Outputs computational results either into an output file or onto 
C     the screen

C----------PARAMETERS.
      PARAMETER (ir=1)             !Input unit number.
      PARAMETER (iw=1)             !Output unit number.
      PARAMETER (nnuc=26)          !Number of nuclides in calculation.
      PARAMETER (itmax=40)         !Maximum # of line to be printed.

C----------COMMON AREAS.
      COMMON /compr/  cy,ct,t9i,t9f,ytmin,inc        !Computation parameters.
      COMMON /modpr/  g,tau,xnu,c,cosmo,xi           !Model parameters.
      COMMON /flags/  ltime,is,ip,it,mbad            !Flags, counters.
      COMMON /outdat/ xout,thmout,t9out,tout,dtout,  !Output data.
     |                etaout,hubout
      COMMON /outopt/ nout,outfile                   !Output option.


C==========================DECLARATION DIVISION=============================

C----------COMPUTATION SETTINGS.
      REAL    cy                   !Time step limiting constant on abundances.
      REAL    ct                   !Time step limiting constant on temperature.
      REAL    t9i                  !Initial temperature (in 10**9 K).
      REAL    t9f                  !Final temperature (in 10**9 K).
      REAL    ytmin                !Smallest abundances allowed.

C----------EARLY UNIVERSE MODEL PARAMETERS.
      REAL    c(3)                !c(1) is variation of gravitational constant.
     |                             !c(2) is neutron lifetime (sec).
     |                             !c(3) is number of neutrino species.
      REAL    cosmo                !Cosmological constant.
      REAL    xi(3)                !Neutrino degeneracy parameters.

C----------COUNTER.
      INTEGER it                   !# times accumulated in output buffer.

C----------OUTPUT ARRAYS.
      REAL    xout(itmax,nnuc)     !Nuclide mass fractions.
      REAL    thmout(itmax,6)      !Thermodynamic variables.
      REAL    t9out(itmax)         !Temperature (in units of 10**9 K).
      REAL    tout(itmax)          !Time.
      REAL    dtout(itmax)         !Time step.
      REAL    etaout(itmax)        !Baryon-to-photon ratio.
      REAL    hubout(itmax)        !Expansion rate.

C----------OUTPUT FILE STATUS.
      INTEGER nout                 !Number of output requests.
      LOGICAL outfile              !Indicates if output file used.

C----------USER INTERACTION VARIABLES.
      INTEGER inum                 !Selection number.


C===========================PROCEDURE DIVISION=============================

C10--------PRINT OUTPUT SELECTION AND AWAIT RESPONSE----------------------

C..........RETURN FROM LOOPING.
 100  CONTINUE
C..........DISPLAY OUTPUT SELECTIONS.
      print 1000
 1000 FORMAT (8(/),
     |        ' ',30x,'OUTPUT SELECTION',/,
     |        ' ',30x,'------ ---------',//,
     |        ' ',25x,' 1. REQUEST OUTPUT FILE',/,
     |        ' ',25x,' 2. REQUEST OUTPUT ON SCREEN',/,
     |        ' ',25x,' 3. EXIT',11(/),
     |        ' ',25x,' Enter selection (1-3): ',$)
C..........READ IN SELECTION NUMBER.
      READ (5,1001) inum
 1001 FORMAT (i1)
C..........BRANCH TO APPROPRIATE SECTION.
      GO TO (200,300,400),inum
      GO TO 400                    !Improper input or <RETURN>.

C20--------REQUEST OUTPUT SECTION--------------------------------------

 200  CONTINUE
c      DO j = 1,it                  !Temperature in MeV.
c        t9out(j) = t9out(j)*.08617
c      END DO
c      DO j = 1,it                  !Energy density as fraction of total.
c        thmout(j,1) = thmout(j,1)/thmout(j,6)  !Rhog.
c        thmout(j,2) = thmout(j,2)/thmout(j,6)  !Rhoe.
c        thmout(j,3) = thmout(j,3)/thmout(j,6)  !Rhone.
c        thmout(j,4) = thmout(j,4)/thmout(j,6)  !Rhob.
c      END DO
C..........PRINT CAPTION.         
        nout = nout + 1            !Keep track of number of output requests.
        IF (nout.eq.1) THEN
          write(2,2000)
 2000     FORMAT (54x,'NUCLIDE ABUNDANCE YIELDS',/,
     |            54x,'------- --------- ------',//)
        END IF
        write(2,2002) cy,ct,t9i,t9f,ytmin
 2002   FORMAT (' Computational parameters:',/,
     |          '   cy = ',f5.3,'/  ct = ',f5.3,
     |          '/  initial temp = ',1pe8.2,
     |          '/  final temp = ',1pe8.2,
     |          '/  smallest abundances allowed = ',1pe8.2)
        write(2,2004) c(1),c(2),c(3),cosmo,xi(1),xi(2),xi(3)
 2004   FORMAT (' Model parameters:',/,
     |          '   g = ',f5.2,'/  tau = ',f6.2,
     |          '/  # nu = ',f5.2,'/  lambda = ',1pe10.3,
     |          '/  xi-e = ',e10.3,'/  xi-m = ',e10.3,
     |          '/  xi-t = ',e10.3,/)
C..........PRINT HEADINGS, ABUNDANCES FOR NEUTRON TO LI8.         
        write(2,2006)
 2006   FORMAT (4x,'Temp',8x,'N/H',10x,'P',10x,'D/H',9x,'T/H',8x,
     |          'He3/H',8x,'He4',8x,'Li6/H',7x,'Li7/H',7x,
     |          'Be7/H',6x,'Li8/H&up',/,132('-'))
        DO j = 1,it
          write(2,2008) t9out(j),(xout(j,i),i=1,10)
 2008     FORMAT (1pe10.3,1p10e12.3)
        END DO
C..........PRINT THERMODYNAMIC QUANTITIES.         
        write(2,2010)
 2010   FORMAT (' ',/,4x,'Temp',9x,'T',10x,'rhog',8x,'rhoe',7x,
     |              'rhone',8x,'rhob',8x,'phie',9x,'dt',9x,
     |              'eta',10x,'H',/,132('-')) 
        DO j = 1,it
          write(2,2012) t9out(j),tout(j),(thmout(j,i),i=1,5),dtout(j),
     |                   etaout(j),hubout(j)
 2012     FORMAT (1pe10.3,9e12.3)
        END DO
        write(2,2014)
 2014   FORMAT (///)
        outfile = .true.           !Output file requested.      
        print 2016
 2016   FORMAT (' ','Output file requested - Press <RETURN> to ',
     |              'continue: ',$)
        READ (5,*)
        GO TO 500

C30--------REQUEST OUTPUT ON SCREEN SECTION------------------------------

C..........RETURN FROM LOOPING.
 300  CONTINUE
c      DO j = 1,it                  !Temperature in MeV.
c        t9out(j) = t9out(j)*.08617 
c      END DO
c      DO j = 1,it                  !Energy density as fraction of total.
c        thmout(j,1) = thmout(j,1)/thmout(j,6)  !Rhog.
c        thmout(j,2) = thmout(j,2)/thmout(j,6)  !Rhoe.
c        thmout(j,3) = thmout(j,3)/thmout(j,6)  !Rhone.
c        thmout(j,4) = thmout(j,4)/thmout(j,6)  !Rhob.
c      END DO
C..........DISPLAY SCREEN OUTPUT SELECTIONS.
        print 3000
 3000   FORMAT (8(/),
     |          ' ',26x,'SCREEN OUTPUT SELECTION',/,
     |          ' ',26x,'------ ------ ---------',//,
     |          ' ',25x,' 1. DISPLAY D,T,HE3,HE4,LI7',/,
     |          ' ',25x,' 2. DISPLAY N,P,LI6,BE7,LI8&UP',/,
     |          ' ',25x,' 3. DISPLAY RHOG,RHOE,RHONE,RHOB',/,
     |          ' ',25x,' 4. DISPLAY T,DT,PHIE,ETA,H',/,
     |          ' ',25x,' 5. EXIT',9(/),
     |          ' ',25x,' Enter selection (1-5): ',$)
C..........READ IN SELECTION NUMBER.
        READ (5,1001) inum
        GO TO (310,320,330,340,350),inum
        GO TO 350                  !Improper input or <RETURN>.
 310    CONTINUE                   !Display d,t,he3,he4,li7.
C..........PRINT CAPTION.
          print 2014
          print 3100, cy,ct,t9i,t9f,ytmin
 3100     FORMAT (' ','Computational parameters:',/,
     |            ' ','   cy = ',f5.3,'/ ct = ',f5.3,
     |                '/ initial temp = ',1pe8.2,
     |                '/ final temp = ',1pe8.2,/,
     |            ' ','   smallest abundances allowed = ',1pe8.2)
          print 3102, c(1),c(2),c(3),cosmo,xi(1),xi(2),xi(3)
 3102     FORMAT (' ','Model parameters:',/,
     |            ' ','   g = ',f5.2,'/ tau = ',f6.2,
     |                '/ # nu = ',f5.2,'/ lambda = ',1pe10.3,/,
     |            ' ','   xi-e = ',e10.3,'/ xi-m = ',e10.3,
     |                '/ xi-t = ',e10.3,/)
C..........PRINT HEADINGS, ABUNDANCES FOR D,T,HE3,HE4,LI7.         
          print 3104
 3104     FORMAT (4x,'Temp',8x,'D/H',9x,'T/H',8x,'He3/H',8x,
     |            'He4',8x,'Li7/H',/,' ',80('-'))
          DO j = 1,it
            print 3106, t9out(j),(xout(j,i),i=3,6),xout(j,8)
 3106       FORMAT (1pe10.3,1p5e12.3)
          END DO
          print 2014
          print 3108
 3108     FORMAT (' ','Press <RETURN> to continue: ',$)
          READ (5,*)
          GO TO 360
 320    CONTINUE                   !Display n,p,li6,be7,li8&up.
C..........PRINT CAPTION.
          print 2014
          print 3100, cy,ct,t9i,t9f,ytmin
          print 3102, c(1),c(2),c(3),cosmo,xi(1),xi(2),xi(3)
C..........PRINT HEADINGS, ABUNDANCES FOR N,P,LI6,BE7,LI8&UP.
          print 3204
 3204     FORMAT (4x,'Temp',8x,'N/H',10x,'P',9x,
     |            'Li6/H',7x,'Be7/H',6x,'Li8/H&up',/,' ',80('-'))
          DO j = 1,it
            print 3106, t9out(j),(xout(j,i),i=1,2),xout(j,7),
     |                      (xout(j,i),i=9,10)
          END DO
          print 2014
          print 3108
          READ (5,*)
          GO TO 360
 330    CONTINUE                   !Display rhog,rhoe,rhone,rhob.
C..........PRINT CAPTION.
          print 2014
          print 3100, cy,ct,t9i,t9f,ytmin
          print 3102, c(1),c(2),c(3),cosmo,xi(1),xi(2),xi(3)
C..........PRINT ENERGY DENSITIES.
          print 3304
 3304     FORMAT (4x,'Temp',8x,'rhog',8x,'rhoe',7x,'rhone',8x,'rhob',
     |            /,' ',80('-'))
          DO j = 1,it
            print 3306, t9out(j),(thmout(j,i),i=1,4)
 3306       FORMAT (1pe10.3,4e12.3)
          END DO
          print 2014
          print 3108
          READ (5,*)
          GO TO 360
 340    CONTINUE                   !Display t,dt,phie,eta,hubcst.
C..........PRINT CAPTION.
          print 2014
          print 3100, cy,ct,t9i,t9f,ytmin
          print 3102, c(1),c(2),c(3),cosmo,xi(1),xi(2),xi(3)
C..........PRINT THERMODYNAMIC QUANTITIES.         
          print 3404
 3404     FORMAT (4x,'Temp',8x,'time',8x,'phie',9x,'dt',9x,'eta',10x,
     |            'H',/,' ',80('-'))
          DO j = 1,it
            print 3406, t9out(j),tout(j),thmout(j,5),dtout(j),
     |                      etaout(j),hubout(j)
 3406       FORMAT (1pe10.3,5e12.3)
          END DO
          print 2014
          print 3108
          READ (5,*)
          GO TO 360
 350    CONTINUE                   !Exit.
          GO TO 500
 360    CONTINUE
        GO TO 300

C40--------EXIT SECTION----------------------------------------------

 400  CONTINUE
        RETURN

C50--------GO BACK TO MENU-------------------------------------------------

 500  CONTINUE
      GO TO 100

      END

C===============IDENTIFICATION DIVISION====================

      SUBROUTINE driver

C-------LINKAGES.
C     CALLED BY - [subroutine] run
C     CALLS     - [subroutine] start, derivs, accum

C-------REMARKS.
C     Runge-Kutta computational routine

C-------PARAMETERS.
      PARAMETER (nvar=29)          !Number of variables to be evolved.
      PARAMETER (nnuc=26)          !Number of nuclides in calculation.
      PARAMETER (cl=1.e-16)        !Lower limit on size of time step.

C-------COMMON AREAS.
      COMMON /evolp1/ t9,hv,phie,y                   !Evolution parameters.
      COMMON /evolp2/ dt9,dhv,dphie,dydt             !Evolution parameters.
      COMMON /evolp3/ t90,hv0,phie0,y0               !Evolution parameters.
      COMMON /compr/  cy,ct,t9i,t9f,ytmin,inc        !Computation parameters.
      COMMON /time/   t,dt,dlt9dt                    !Time variables.
      COMMON /flags/  ltime,is,ip,it,mbad            !Flags,counters.
      COMMON /check1/  itime                          !Computation location.
      COMMON /runopt/ irun,isize,jsize               !Run options.


C=================DECLARATION DIVISION=====================

C-------EVOLUTION PARAMETERS.
      REAL    t9                   !Temperature (in units of 10**9 K).
      REAL    hv                   !Defined by hv = M(atomic)n(baryon)/t9**3.
      REAL    phie                 !Chemical potential for electron.
      REAL    y(nnuc)              !Relative number abundances.

C-------EVOLUTION PARAMETERS (DERIVATIVES).
      REAL    dydt(nnuc)           !Change in rel number abundances.

C-------EVOLUTION PARAMETERS (ORIGINAL VALUES).
      REAL    y0(nnuc)             !Rel # abund at beginning of iteration.

C-------COMPUTATION PARAMETERS.
      REAL    cy                   !Time step limiting constant on abundances.
      REAL    ct                  !Time step limiting constant on temperature.
      REAL    t9f                  !Final temperature (in 10**9 K).
      REAL    ytmin                !Smallest abundances allowed.
      INTEGER inc                  !Accumulation increment.

C-------TIME AND TIME STEP VARIABLES.
      REAL    t                    !Time.
      REAL    dt                   !Time step.
      REAL    dlt9dt               !(1/t9)*d(t9)/d(t).

C-------COUNTERS AND FLAGS.
      INTEGER loop                 !Counts which Runge-Kutta loop.
      INTEGER ltime                !Indicates termination status.
      INTEGER is                   !# total time steps for particular run.
      INTEGER ip                   !# time steps after outputting a line.

C-------COMPUTATION LOCATION.
      INTEGER itime                !Time check.

C-------RUN OPTION.
      INTEGER isize                !Number of nuclides in computation.

C-------TIME AND TIME STEP VARIABLES.
      REAL    dtmin                !Mininum time step.
      REAL    dtl                 !Time step from limitation on abund changes.

C-------LABELS FOR VARIABLES TO BE TIME EVOLVED.
      INTEGER mvar                 !Total number of variables to be evolved.
      REAL    v(nvar)              !Variables to be time evolved.
      REAL    dvdt(nvar)           !Time derivatives.
      REAL    v0(nvar)             !Value of variables at original point.
      REAL    dvdt0(nvar)          !Value of derivatives at original point.

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
     |    (dt.lt.abs(cl/dlt9dt)).or.              !Small dt.
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

C-------PARAMETERS.
      PARAMETER (nrec=88)          !Number of nuclear reactions.
      PARAMETER (nnuc=26)          !Number of nuclides in calculation.
      PARAMETER (const1=0.09615)   !Relation between time and temperature.
      PARAMETER (const2=6.6700e-8) !Gravitational constant.

C-------COMMON AREAS.
      COMMON /rates/  f,r                            !Reaction rates.
      COMMON /evolp1/ t9,hv,phie,y                   !Evolution parameters.
      COMMON /evolp2/ dt9,dhv,dphie,dydt(nnuc)       !Evolution parameters.
      COMMON /evolp3/ t90,hv0,phie0,y0               !Evolution parameters.
      COMMON /compr/  cy,ct,t9i,t9f,ytmin,inc        !Computation parameters.
      COMMON /modpr/  g,tau,xnu,c,cosmo,xi           !Model parameters.
      COMMON /varpr/  dt1,eta1                       !Variational parameters.
      COMMON /time/   t,dt,dlt9dt                    !Time variables.
      COMMON /endens/ rhone0,rhob0,rhob,rnb          !Energy densities.
      COMMON /xbessel/ bl1,bl2,bl3,bl4,bl5,          !Eval of function bl(z).
     |                bm1,bm2,bm3,bm4,bm5,           !Eval of function bm(z).
     |                bn1,bn2,bn3,bn4,bn5            !Eval of function bn(z).
      COMMON /flags/  ltime,is,ip,it,mbad            !Flags,counters.
      COMMON /nupar/  t9mev,tnmev,tnu,cnorm,nu,rhonu !Neutrino parameters.
      COMMON /runopt/ irun,isize,jsize               !Run options.


C=================DECLARATION DIVISION=====================

C-------REACTION RATES.
      REAL    f(nrec)              !Forward reaction rate coefficients.
	real r(nrec)
C-------EVOLUTION PARAMETERS.
      REAL    t9                   !Temperature (in units of 10**9 K).
      REAL    hv                   !Defined by hv = M(atomic)n(baryon)/t9**3.
      REAL    phie                 !Chemical potential of electron.
      REAL    y(nnuc)              !Relative number abundances.

C-------EVOLUTION PARAMETERS (ORIGINAL VALUES).
      REAL    y0(nnuc)             !Rel # abund at start of iteration.

C-------COMPUTATION SETTINGS.
      REAL    t9i                  !Initial temperature (in 10**9 K).
      REAL    ytmin                !Smallest abundances allowed.
      INTEGER inc                  !Accumulation increment.

C-------EARLY UNIVERSE MODEL PARAMETERS.
      REAL    g                    !Gravitational constant.
      REAL    tau                  !Neutron lifetime.
      REAL    xnu                  !Number of neutrino species.
      REAL    c(3)               !c(1) is variation of gravitational constant.
     |                             !c(2) is neutron lifetime (sec).
     |                             !c(3) is number of neutrino species.
      REAL    xi(3)                !Neutrino degeneracy parameters.

C-------VARIATIONAL PARAMETERS.
      REAL    dt1                  !Initial time step.
      REAL    eta1                 !Baryon-to-photon ratio.

C-------TIME VARIABLES.
      REAL    t                    !Time.
      REAL    dt                   !Time step.

C-------ENERGY DENSITIES.
      REAL    rhone0               !Initial electron neutrino mass density.
      REAL    rhob0                !Initial baryon mass density.

C-------EVALUATION OF FUNCTIONS bl,bm,bn.
      REAL    bl1,bl2,bl3,bl4,bl5  !Evaluation of function bl(z).

C-------COUNTERS AND FLAGS.
      INTEGER ltime                !Indicates if output buffer printed.
      INTEGER is                   !# total time steps for particular run.
      INTEGER ip                   !# time steps after outputting a line.
      INTEGER it                   !# times accumulated in output buffer.
      INTEGER mbad                 !Indicates if gaussian elimination fails.

C-------NEUTRINO PARAMETERS.
      REAL    tnu                  !Neutrino temperature.
      REAL    cnorm                !Normalizing constant.

C-------RUN OPTION.
      INTEGER isize                !Number of nuclides in computation.

C-------LOCAL VARIABLES.
      REAL    z                    !Defined by z = m(electron)*c**2/k*t9.


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

C-------PARAMETERS.
      PARAMETER (nvar=29)          !Number of variables to be evolved.
      PARAMETER (nnuc=26)          !Number of nuclides in calculation.
      PARAMETER (pi=3.141593)

C-------COMMON AREAS.
      COMMON /evolp1/ t9,hv,phie,y                   !Evolution parameters.
      COMMON /evolp2/ dt9,dhv,dphie,dydt             !Evolution parameters.
      COMMON /evolp3/ t90,hv0,phie0,y0               !Evolution parameters.
      COMMON /modpr/  g,tau,xnu,c(3),cosmo,xi(3)     !Model parameters.
      COMMON /time/   t,dt,dlt9dt                    !Time variables.
      COMMON /xtherm/  thm,hubcst                     !Dynamic variables.
      COMMON /endens/ rhone0,rhob0,rhob,rnb          !Energy densities.
      COMMON /nucdat/ am(nnuc),zm,dm                 !Nuclide data.
      COMMON /flags/  ltime,is,ip,it,mbad            !Flags,counters.
      COMMON /runopt/ irun,isize,jsize               !Run options.


C=================DECLARATION DIVISION=====================

C-------EVOLUTION PARAMETERS.
      REAL    t9                   !Temperature (in units of 10**9 K).
      REAL    hv                   !Defined by hv = M(atomic)n(baryon)/t9**3.
      REAL    phie                 !Chemical potential for electron.
      REAL    y(nnuc)              !Relative number abundances.

C-------EVOLUTION PARAMETERS (DERIVATIVES).
      REAL    dt9                  !Change in temperature.
      REAL    dhv                  !Change in hv.
      REAL    dphie                !Change in chemical potential.
      REAL    dydt(nnuc)           !Change in rel number abundances.

C-------EVOLUTION PARAMETERS (ORIGINAL VALUES).
      REAL    y0(nnuc)             !Rel # abund at beginning of iteration.

C-------MODEL PARAMETERS.
      REAL    g                    !Gravitational constant.
      REAL    cosmo                !Cosmological constant.

C-------TIME VARIABLES.
      REAL    dlt9dt               !(1/t9)*d(t9)/d(t).

C-------DYNAMIC VARIABLES.
      REAL    thm(14)              !Thermodynamic variables.
      REAL    hubcst               !Expansion rate.

C-------ENERGY DENSITIES.
      REAL    rhob0                !Initial baryon mass density.
      REAL    rhob                 !Baryon mass density.
      REAL    rnb                  !Baryon mass density (ratio to init value).

C-------NUCLIDE DATA.
      REAL    zm(nnuc)             !Charge of nuclide.
      REAL    dm(nnuc)             !Mass excess of nuclide.

C-------COUNTERS AND FLAGS.
      INTEGER mbad                 !Indicates if gaussian elimination fails.

C-------RUN OPTION.
      INTEGER irun                 !Run network size.
      INTEGER isize                !Number of nuclides in computation.

C-------SUMS.
      REAL    sumy                 !Sum of abundances.
      REAL    sumzy                !Sum of charge*abundances.
      REAL    sumdy                !Sum of abundance flows.
      REAL    summdy               !Sum of (mass excess)*(abundance flows).
      REAL    sumzdy               !Sum of (charge)*(abundance flows).

C-------DERIVATIVES.
      REAL    dphdt9               !d(phi e)/d(t9).
      REAL    dphdln               !d(phi e)/d(h).
      REAL    dphdzy               !d(phi e)/d(sumzy).
      REAL    dlndt9               !(1/h)*d(h)/d(t9).
      REAL    bar                  !Baryon density and pressure terms.

C-------LOCAL VARIABLES.
      INTEGER loop                 !Counts which Runge-Kutta loop.


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
      summdy = summdy + dm(i)*dydt(i)  !Sum of (mass excess)*(abundance flow).
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
      dt9    = (3.*hubcst)/dlndt9
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

C-------PARAMETERS.
      PARAMETER (nvar=29)          !Number of variables to be evolved.
      PARAMETER (nnuc=26)          !Number of nuclides in calculation.
      PARAMETER (itmax=40)         !Maximum # of lines to be printed.

C-------COMMON AREAS.         
      COMMON /evolp1/ t9,hv,phie,y                   !Evolution parameters.
      COMMON /compr/  cy,ct,t9i,t9f,ytmin,inc        !Computation parameters.
      COMMON /time/   t,dt,dlt9dt                    !Time variables.
      COMMON /xtherm/  thm,hubcst                     !Dynamic variables.
      COMMON /nucdat/ am,zm(nnuc),dm(nnuc)           !Nuclide data.
      COMMON /flags/  ltime,is,ip,it,mbad            !Flags,counters.
      COMMON /outdat/ xout,thmout,t9out,tout,dtout,  !Output data.
     |                etaout,hubout
      COMMON /runopt/ irun,isize,jsize               !Run options.


C=================DECLARATION DIVISION=====================

C-------EVOLUTION PARAMETERS.
      REAL    t9                   !Temperature (in units of 10**9 K).
      REAL    hv                   !Defined by hv = M(atomic)n(baryon)/t9**3.
      REAL    phie                 !Chemical potential for electron.
      REAL    y(nnuc)              !Relative number abundances.

C-------COMPUTATION PARAMETERS.
      INTEGER inc                  !Accumulation increment.

C-------TIME PARAMETERS.
      REAL    t                    !Time.
      REAL    dt                   !Time step.

C-------DYNAMIC VARIABLES.
      REAL    thm(14)              !Thermodynamic variables.
      REAL    hubcst               !Expansion rate.

C-------NUCLIDE DATA.
      REAL    am(nnuc)             !Atomic number of nuclide.

C-------COUNTERS AND FLAGS.
      INTEGER ltime                !Indicates if output buffer printed.
      INTEGER it                   !# times accumulated in output buffer.
      INTEGER ip                   !# time steps after outputting a line.

C-------OUTPUT ARRAYS.
      REAL    xout(itmax,nnuc)     !Nuclide mass fractions.
      REAL    thmout(itmax,6)      !Thermodynamic variables.
      REAL    t9out(itmax)         !Temperature (in units of 10**9 K).
      REAL    tout(itmax)          !Time.
      REAL    dtout(itmax)         !Time step.
      REAL    etaout(itmax)        !Baryon-to-photon ratio.
      REAL    hubout(itmax)        !Expansion rate.

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

C-------PARAMETER.
      PARAMETER (nnuc=26)          !Number of nuclides in calculation.
      PARAMETER (q=2.531)          !(mass(neutron)-mass(proton))/m(electron)

C-------COMMON AREAS.         
      COMMON /evolp1/ t9,hv,phie,y(nnuc)             !Evolution parameters.
      COMMON /compr/  cy,ct,t9i,t9f,ytmin,inc        !Computation parameters.
      COMMON /modpr/  g,tau,xnu,c(3),cosmo,xi        !Model parameters.
      COMMON /xtherm/  thm,hubcst                     !Dynamic variables.
      COMMON /endens/ rhone0,rhob0,rhob,rnb          !Energy densities.
      COMMON /xbessel/ bl1,bl2,bl3,bl4,bl5,           !Eval of function bl(z).
     |                bm1,bm2,bm3,bm4,bm5,           !Eval of function bm(z).
     |                bn1,bn2,bn3,bn4,bn5            !Eval of function bn(z).
      COMMON /nupar/  t9mev,tnmev,tnu,cnorm,nu,rhonu !Integration parameters.


C=================DECLARATION DIVISION=====================

C-------EVOLUTION PARAMETERS.
      REAL    t9                   !Temperature (in units of 10**9 K).
      REAL    phie                 !Chemical potential for electron.

C-------COMPUTATION PARAMETERS.
      REAL    t9i                  !Initial temperature (in 10**9 K).

C-------EARLY UNIVERSE MODEL PARAMETERS.
      REAL    xnu                  !Number of neutrino species.
      REAL    xi(3)                !Neutrino degeneracy parameters.

C-------DYNAMIC VARIABLES.
      REAL    thm(14)              !Thermodynamic variables.

C-------ENERGY DENSITIES.
      REAL    rhone0               !Initial electron neutrino mass density.
      REAL    rhob0                !Initial baryon mass density.
      REAL    rnb                  !Baryon mass density (ratio to init value).

C-------EVALUATION OF FUNCTIONS bl,bm,bn.
      REAL    bl1,bl2,bl3,bl4,bl5  !Evaluation of function bl(z).
      REAL    bm1,bm2,bm3,bm4,bm5  !Evaluation of function bm(z).
      REAL    bn1,bn2,bn3,bn4,bn5  !Evaluation of function bn(z).

C-------NEUTRINO PARAMETERS.
      REAL    tnu                  !Neutrino temperature.
      REAL    rhonu                !Neutrino energy density.
      INTEGER nu                   !Type of neutrino.

C-------LOCAL VARIABLE.
      REAL    z                    !Defined by z = m(electron)*c**2/k*t9.


C==================PROCEDURE DIVISION======================

C10-----COMPUTE FACTORS------------------------------------

      z = 5.930/t9                 !z = m(electron)c**2/k(t9).
      tnu = ((rnb)**(1./3.))*t9i   !Neutrino temperature.
C..........FACTORS OF z.
      z1 = z
      z2 = z*z
      z3 = z*z*z
      z4 = z*z*z*z
      z5 = z*z*z*z*z
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
        DO nu = 1,xnu              !For every neutrino family.
          CALL nudens              !Compute neutrino energy density.
          thm(8) = thm(8) + 12.79264*rhonu  !Have 12.79264 from units change.
        END DO
      END IF
      thm(9)  = rhob0*rnb                                       !(Ref 9).
      thm(10) = thm(1) + thm(4) + thm(8) + thm(9)               !(Ref 10).
      thm(11) = -(z**3/t9)*(sinh1*(3.*bl1-z*bm1)-sinh2*(3.*bl2  !(Ref 11).
     |          -2.*z*bm2) + sinh3*(3.*bl3-3.*z*bm3) - sinh4
     |          *(3.*bl4-4.*z*bm4) + sinh5*(3.*bl5-5.*z*bm5))
      thm(12) = z**3*(cosh1*bl1- 2.*cosh2*bl2                   !(Ref 12).
     |          + 3.*cosh3*bl3 - 4.*cosh4*bl4 + 5.*cosh5*bl5)
      IF (thm(12).ne.0.) thm(12) = 1./thm(12)
      thm(13) = 1.000 + 0.565/z1 - 6.382/z2 + 11.108/z3         !(Ref 13).
     |          + 36.492/z4 + 27.512/z5
      thm(14) = (5.252/z1 - 16.229/z2 + 18.059/z3 + 34.181/z4   !(Ref 14).
     |          + 27.617/z5)*ex(-q*z)

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
      REAL   bl1,bl2,bl3,bl4,bl5  !Single variables equivalenced to array blz.
      REAL   bm1,bm2,bm3,bm4,bm5  !Single variables equivalenced to array bmz.
      REAL   bn1,bn2,bn3,bn4,bn5  !Single variables equivalenced to array bnz.

C-------EVALUATIION OF MODIFIED BESSEL FUNCTIONS.
      REAL    bk0,bk1,bk2,bk3,bk4  !Values k0(r),k1(r),k2(r),k3(r),k4(r).

C-------LOCAL VARIABLES.
      REAL    blz(5)               !Array containing values from function bl.
      REAL    bmz(5)               !Array containing values from function bm.
      REAL    bnz(5)               !Array containing values from function bn.
      REAL    z                    !Defined by z = m(electron)*c**2/k*t9.
      REAL    r                    !Multiples of z.

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
      REAL    bk0,bk1              !Values k0(z),k1(z)
      REAL    bi0,bi1              !Values i0(z),i1(z).
      REAL    bk2,bk3,bk4          !Values k2(z),k3(z),k4(z).

C--------EXPANSION COEFFICIENTS.
      REAL    ci0(7)               !Expansion coefficients for i0 (z.le.2).
      REAL    ci1(7)               !Expansion coefficients for i1 (z.le.2).
      REAL    ck0(7)               !Expansion coefficients for k0 (z.le.2).
      REAL    ck1(7)               !Expansion coefficients for k1 (z.le.2).
      REAL    c0(7)                !Expansion coefficients for k0 (z.gt.2).
      REAL    c1(7)                !Expansion coefficients for k1 (z.gt.2).

C--------VARIABLES TO BE EVALUATED.
      REAL    z                    !Input variable.
      REAL    y                    !Expansion variable = z/2.
      REAL    t                    !Expansion variable = z/3.75.
      REAL    coeff                !Logrithmic or exponential coefficient.


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

C-------PARAMTER.
      PARAMETER (iter=50)          !Number of gaussian quads.

C-------COMMON AREAS.
      COMMON /modpr/  g,tau,xnu,c(3),cosmo,xi        !Model parameters.
      COMMON /nupar/  t9mev,tnmev,tnu,cnorm,nu,rhonu !Integration parameters.

C-------EXTERNAL FUNCTIONS.
      EXTERNAL func5               !Integral for neutrinos.
      EXTERNAL func6               !Integral for antineutrinos.


C=================DECLARATION DIVISION=====================

C-------MODEL PARAMETERS.
      REAL    xi(3)                !Neutrino degeneracy parameters.
	real func5
	real func6

C-------NEUTRINO PARAMETERS.
      REAL    tnu                  !Neutrino temperature (units of 10**9 K).
      REAL    rhonu                !Neutrino energy density.
      INTEGER nu                   !Which neutrino type.

C-------LOCAL VARIABLES.
      REAL    uplim1               !Upper limit for neutrino energy integral.
      REAL    uplim2               !Upper limit for antineu energy integral.


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
      COMMON /modpr/  g,tau,xnu,c(3),cosmo,xi        !Model parameters.
      COMMON /nupar/  t9mev,tnmev,tnu,cnorm,nu,rhonu !Integration parameters.

C=================DECLARATION DIVISION=====================

C-------MODEL PARAMETERS.
      REAL    xi(3)                !Neutrino degeneracy parameters.

C-------NEUTRINO PARAMETERS.
      REAL    t9mev                !Temperature (in units of MeV).
      REAL    tnmev                !Neutrino temperature (in units of MeV).
      REAL    tnu                  !Neutrino temperature (units of 10**9 K).
      REAL    cnorm                !Normalizing constant.

C-------FUNCTIONS TO BE EVALUATED.
      REAL   func1                 !1st part n->p rate.
      REAL   func2                 !2nd part n->p rate.
      REAL   func3                 !1st part p->n rate.
      REAL   func4                 !2nd part p->n rate.
      REAL   func5                 !Energy density for neutrinos.
      REAL   func6                 !Energy density for antineutrinos.

C-------LOCAL VARIABLES.
      REAL    x                    !Value at which function is evaluated.
      REAL    part1                !Exponential expression with photon temp.
      REAL    part2                !Exponential expression with neutrino temp.


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
      REAL    xlow                 !Array of low limits.
      REAL    xhi                  !Array of high limits.
      INTEGER nq                   !Number of six point gaussian quads.

C-------COMPUTATION VARIABLES.
      REAL    dist                 !Size of quad interval.
      REAL    cent                 !Center of quad interval.
      REAL    x                    !Variables of integration.
      REAL    sum                  !Summation of terms.

C-------COUNTERS.
      INTEGER nint                 !Interval number.
      INTEGER npnt                 !Point number.
      INTEGER np                   !Total number of points in interval.

C-------ABSCISSAS AND WEIGHT FACTORS.
      REAL    u(6)                 !Abscissas.
      REAL    w(6)                 !Weight factor.


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

C-------PARAMETERS.
      PARAMETER (ir=1)             !Input unit number.
      PARAMETER (iw=1)             !Output unit number.
      PARAMETER (nrec=88)          !Number of nuclear reactions.
      PARAMETER (nnuc=26)          !Number of nuclides in calculation.

C--------COMMON AREAS.
      COMMON /recpr/  iform,ii,jj,kk,ll,rev,q9     !Reaction parameters names.
      COMMON /rates/  f,r                            !Reaction rates.
      COMMON /evolp1/ t9,hv,phie,y                   !Evolution parameters.
      COMMON /evolp2/ dt9,dhv,dphie,dydt             !Evolution parameters.
      COMMON /evolp3/ t90,hv0,phie0,y0               !Evolution parameters.
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
      REAL    rev(nrec)            !Reverse reaction coefficient.
      REAL    q9(nrec)             !Energy released in reaction.

C-------REACTION RATES.
      REAL    f(nrec)              !Forward reaction rate coefficients.
      REAL    r(nrec)              !Reverse reaction rate coefficients.

C-------EVOLUTION PARAMETERS.
      REAL    t9                   !Temperature (in units of 10**9 K).
      REAL    y(nnuc)              !Relative number abundances.

C-------EVOLUTION PARAMETERS (DERIVATIVES).
      REAL    dydt(nnuc)           !Change in rel number abundances.

C-------EVOLUTION PARAMETERS (ORIGINAL VALUES).
      REAL    y0(nnuc)             !Rel # abund at start of iteration.

C-------TIME VARIABLES.
      REAL    dt                   !Time step.

C-------DYNAMIC VARIABLES.
      REAL    hubcst               !Expansion rate.

C-------ENERGY DENSITIES.
      REAL    rhob                 !Baryon mass density.

C-------COMPONENTS OF MATRIX EQUATION.
      DOUBLE PRECISION a(nnuc,nnuc)!Relates y(t-dt) to y(t).
      REAL    b(nnuc)              !Contains y0 in inverse order.
      REAL    yx(nnuc)             !yy in reverse order.

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
      REAL    ri,rj,rk,rl          !Equate to si,sj,sk,sl.
      REAL    ci,cj,ck,cl          !Coefficients of rate equation.

C-------LOCAL VARIABLES.
      REAL    yy(nnuc)             !Abundances at end of iteration.
      REAL    si(11),sj(11),sk(11),sl(11)  !# of nuclide i,j,k,l
      REAL    bdln                 !(10**(-5))*volume expansion rate.
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
          IF (dabs(a(j,i)).lt.bdln*y0(j1)/y0(i1)) THEN
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

C-------PARAMETERS.
      PARAMETER (nnuc=26)          !Rank of matrix.
      PARAMETER (mord=1)           !Higher order in correction.
      PARAMETER (eps=2.e-4)        !Tolerance for convergence (.ge. 1.e-7).

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
      REAL    b(nnuc)              !Right-hand vector w/o manipulation.
      REAL    y(nnuc)              !Solution vector.

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
      REAL   xdy                   !Relative error.
 
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
            xdy = dabs(x(i)/y(i))  !Relative error.
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

C===============IDENTIFICATION DIVISION====================

      SUBROUTINE rate0

C-------LINKAGES.
C     CALLED BY - [subroutine] start
C     CALLS     - none

C-------REMARKS.
C     Generates weak decay rates.

C-------PARAMETER.
      PARAMETER (nrec=88)          !Number of nuclear reactions.

C-------COMMON AREA.
      COMMON /rates/  f,r(nrec)    !Reaction rates.


C=================DECLARATION DIVISION=====================

C-------REACTION RATE COEFFICIENTS.
      REAL    f(nrec)              !Forward reaction rate coefficients.


C==================PROCEDURE DIVISION======================

C10-----SET DECAY RATE COEFFICIENTS---------------------------

C.......H3 -> e- + v + He3.........(Tilly-Weller-Hasan 1987)
      f(2)  = 1.79e-9

C.......Li8 -> e- + v + 2He4.......(Ajzenberg-Selove 1988)
      f(3)  = 8.27e-1

C.......B12 -> e- + B + C12........(Ajzenberg-Selove 1990)
      f(4)  = 3.43e+1

C.......C14 -> e- + v + N14........(Ajzenberg-Selove 1986)
      f(5)  = 3.834e-12

C.......B8 -> e+ + v + 2He4........(Ajzenberg-Selove 1988)
      f(6)  = 9.00e-1

C.......C11 -> e+ + v + B11........(Ajzenberg-Selove 1990)
      f(7)  = 5.668e-4

C.......N12 -> e+ + v + C12........(Ajzenberg-Selove 1990)
      f(8)  = 6.301e+1

C.......N13 -> e+ + v + C13........(Ajzenberg-Selove 1986)
      f(9)  = 1.159e-3

C.......O14 -> e+ + v + N14........(Ajzenberg-Selove 1986)
      f(10) = 9.8171e-3

C.......O15 -> e+ + v + N15........(Ajzenberg-Selove 1986)
      f(11) = 5.6704e-3

      RETURN

C-------REFERENCES--------------------------------------
C     Ajzenberg-Selove, F., 1990, Nucl. Phys. A506, 1.
C     Ajzenberg-Selove, F., 1988, Nucl. Phys. A490, 1.
C     Ajzenberg-Selove, F., 1986, Nucl. Phys. A449, 1.
C     Tilley, D.R., Weller, H.R., and Hasan, H.H., 1987, Nucl. Phys. A474, 1.

      END



C===============IDENTIFICATION DIVISION====================

      SUBROUTINE rate1(tph)

C-------LINKAGES.
C     CALLED BY - [subroutine] start, derivs
C     CALLS     - [function] xintd, eval

C-------REMARKS.
C     Generates rate coefficients for weak n->p and p->n reactions.

C-------PARAMETERS.
      PARAMETER (nrec=88)          !Number of nuclear reactions.
      PARAMETER (iter=50)          !Number of gaussian quads.

C-------COMMON AREAS.
      COMMON /rates/  f,r                            !Reaction rates.
      COMMON /modpr/  g,tau,xnu,c(3),cosmo,xi        !Model parameters.
      COMMON /xtherm/  thm,hubcst                     !Dynamic variables.
      COMMON /nupar/  t9mev,tnmev,tnu,cnorm,nu,rhonu !Integration parameters.

C-------EXTERNAL FUNCTIONS.
      EXTERNAL func1               !Part 1 of n->p rate.
      EXTERNAL func2               !Part 2 of n->p rate.
      EXTERNAL func3               !Part 1 of p->n rate.
      EXTERNAL func4               !Part 2 of p->n rate.


C=================DECLARATION DIVISION=====================

C-------REACTION RATES.
      REAL    f(nrec)              !Forward reaction rate coefficients.
      REAL    r(nrec)              !Reverse reaction rate coefficients.

C-------EARLY UNIVERSE MODEL PARAMETERS.
      REAL    tau                  !Neutron lifetime.
      REAL    xi(3)                !Neutrino degeneracy parameters.
     |                             !xi(1) is e neutrino degeneracy parameter.
     |                             !xi(2) is m neutrino degeneracy parameter.
     |                             !xi(3) is t neutrino degeneracy parameter.

C-------DYNAMIC VARIABLES.
      REAL    thm(14)              !Thermodynamic variables (energy densities).

C-------NEUTRINO PARAMETERS.
      REAL    t9mev                !Temperature (in units of MeV).
      REAL    tnmev                !Neutrino temperature (in units of MeV).

C-------LOCAL VARIABLES.
      REAL    tph                  !Photon temperature.
      REAL    w(2),x(2),          !Upper limits for exponentials, forward rate.
     |        y(2),z(2)           !Upper limits for exponentials, reverse rate.
      REAL    uplim1,uplim2,      !Upper limits for integrals for forward rate.
     |        uplim3,uplim4       !Upper limits for integrals for reverse rate.
      REAL    part1,part2,         !Parts of integrals for forward rate.
     |        part3,part4          !Parts of integrals for reverse rate.


C==================PROCEDURE DIVISION======================

C10-----COMPUTE WEAK REACTION RATES (NONDEGENERATE)-----------------

      IF (xi(1).eq.0.) THEN
        f(1)  = thm(13)/tau        !Forward rate for weak np reaction.
        r(1)  = thm(14)/tau        !Reverse rate for weak np reaction.
      ELSE

C20-----COMPUTE WEAK REACTION RATES (DEGENERATE)--------------------

        t9mev = tph*.086171        !Convert photon temp to units of MeV.
        tnmev = tnu*.086171        !Convert neutrino temp to units of MeV.
C..........COMPUTE OVERFLOW LIMITS FOR LIMITS OF INTEGRATION (Ref 1 & 2).
        w(1) = (-(t9mev/.511)*(-88.722))
        w(2) = ((tnmev/.511)*(88.029+xi(1))+2.531)
        x(1) = ((t9mev/.511)*(88.029))
        x(2) = (-(tnmev/.511)*(-88.722+xi(1))-2.531)
        y(1) = (-(t9mev/.511)*(-88.722))
        y(2) = ((tnmev/.511)*(88.029-xi(1))-2.531)
        z(1) = ((t9mev/.511)*(88.029))
        z(2) = (-(tnmev/.511)*(-88.722-xi(1))+2.531)
C..........COMPARE LIMITS AND TAKE LARGER OF THE TWO.
        uplim1 = abs(w(1))
        uplim2 = abs(x(1))
        uplim3 = abs(y(1))
        uplim4 = abs(z(1))
        IF (uplim1.lt.abs(w(2))) uplim1 = w(2)
        IF (uplim2.lt.abs(x(2))) uplim2 = x(2)
        IF (uplim3.lt.abs(y(2))) uplim3 = y(2)
        IF (uplim4.lt.abs(z(2))) uplim4 = z(2)
C..........EVALUATE THE INTEGRALS NUMERICALLY.
        part1 = xintd(1.,uplim1,func1,iter)
        part2 = xintd(1.,uplim2,func2,iter)
        part3 = xintd(1.,uplim3,func3,iter)
        part4 = xintd(1.,uplim4,func4,iter)
        f(1) = part1 + part2       !Add 2 integrals to get forward rate.
        r(1) = part3 + part4       !Add 2 integrals to get reverse rate.
      END IF !(xi(1).eq.0.)
      RETURN

C-------REFERENCES--------------------------------------
C     1) Forms of the integrals involved can be found in
C          Scherrer,R.J., 1983, Mon.Not.R.astr.Soc., 205, 683.
C          Beaudet,G. and Goret,P., 1976, Astron. & Astrophys., 49, 415.
C
C     2) The overflow limit for the VAX/VMS system is exp(88.029).
C        The underflow limit for the VAX/VMS system is exp(-88.722).

      END



C===============IDENTIFICATION DIVISION====================

      SUBROUTINE rate2

C-------LINKAGES.
C     CALLED BY - [subroutine] derivs
C     CALLS     - [function] ex

C-------REMARKS.
C     Generates rate coefficients for reactions involving nuclides
C     up to A = 9.

C-------PARAMETER.
      PARAMETER (nrec=88)          !Number of nuclear reactions.
      PARAMETER (nnuc=26)          !Number of nuclides in calculation.

C-------COMMON AREAS.
      COMMON /rates/  f,r(nrec)           !Reaction rates.
      COMMON /evolp1/ t9,hv,phie,y(nnuc)  !Evolution parameters.


C=================DECLARATION DIVISION=====================

C-------REACTION RATES.
      REAL    f(nrec)              !Forward reaction rate coefficients.

C-------EVOLUTION PARAMETER.
      REAL    t9                   !Temperature of photons (units of 10**9 K).


C==================PROCEDURE DIVISION======================

C10-----TEMPERATURE FACTORS--------------------------------

      t913  = t9**(.33333333)      !t9**(1/3)
      t923  = t913*t913            !t9**(2/3)
      t943  = t923*t923            !t9**(4/3)
      t953  = t9*t923              !t9**(5/3)
      t912  = sqrt(t9)             !t9**(1/2)
      t932  = t9*t912              !t9**(3/2)
      t9m1  = 1/t9                 !t9**(-1)
      t9m23 = 1.0/t923             !t9**(-2/3)
      t9m32 = 1.0/t932             !t9**(-3/2)
      t9a   = t9/(1.0+13.076*t9)   !For reaction 17.
      t9a32 = t9a**(1.5)           !t9a**(3/2)
      t9b   = t9/(1.+49.18*t9)     !For reaction 18.
      t9b32 = t9b**(1.5)           !t9b**(3/2)
      IF (t9.gt.10.) THEN          !For reaction 22.
        t9c = 1.
      ELSE
        t9c = t9/(1.-9.69e-2*t9+2.84e-2*t953/(1.-9.69e-2*t9)**(2./3.))
      END IF
      t9c13 = t9c**(.3333333)      !t9c**(1/3)
      t9c56 = t9c**(.8333333)      !t9c**(5/6)
      t9d   = t9/(1.+0.759*t9)     !For reaction 24.
      t9d13 = t9d**(.3333333)      !t9d**(1/3)
      t9d56 = t9d**(.8333333)      !t9d**(5/6)
      t9e   = t9/(1.+0.1378*t9)    !For reaction 26.
      t9e13 = t9e**(.3333333)      !t9e**(1/3)
      t9e56 = t9e**(.8333333)      !t9e**(5/6)
      t9f   = t9/(1.+0.1071*t9)    !For reaction 27.
      t9f13 = t9f**(.3333333)      !t9f**(1/3)
      t9f56 = t9f**(.8333333)      !t9f**(5/6)


C20-----NEUTRON, PHOTON REACTIONS-----------------------------

C.......H(n,g)H2...................(Smith-Kawano-Malaney 1992)
      f(12)  = 4.742e+4*(1.-.8504*t912+.4895*t9-.09623*t932
     |                     +8.471e-3*t9*t9-2.80e-4*t9*t932)

C.......H2(n,g)H3..................(Wagoner 1969)
      f(13)  = 6.62e+1*(1.+18.9*t9)

C.......He3(n,g)He4................(Wagoner 1969)
      f(14)  = 6.62e+0*(1.+905.*t9)

C.......Li6(n,g)Li7................(Malaney-Fowler 1989)
      f(15)  = 5.10e+3

C30-----NEUTRON, PROTON REACTIONS-----------------------------

C.......He3(n,p)H3.................(Smith-Kawano-Malaney 1992)
      f(16)  = 7.21e+8*(1.-.508*t912+.228*t9)

C.......Be7(n,p)Li7................(Smith-Kawano-Malaney 1992)
      f(17)  = 2.675e+9*(1.-.560*t912+.179*t9-.0283*t932
     |        + 2.214e-3*t9*t9-6.851e-5*t9*t932)
     |        + 9.391e+8*t9a32*t9m32
     |        + 4.467e+7*t9m32*ex(-0.07486/t9)

C40-----NEUTRON, ALPHA REACTIONS------------------------------

C.......Li6(n,a)H3.................(Caughlan-Fowler 1988)
      f(18)  = 2.54e+9*t9m32*ex(-2.39/t9)
     |         + 1.68e+8*(1.-.261*t9b32/t932)

C.......Be7(n,a)He4................(Wagoner 1969)
      f(19)  = 2.05e+4*(1.+3760.*t9)

C50-----PROTON, PHOTON REACTIONS------------------------------

C.......H2(p,g)He3.................(Smith-Kawano-Malaney 1992)
      f(20)  = 2.65e+3*t9m23*ex(-3.720/t913)
     |         *(1.+.112*t913+1.99*t923+1.56*t9+.162*t943+.324*t953)

C.......H3(p,g)He4.................(Caughlan-Fowler 1988)
      f(21)  = 2.20e+4*t9m23*ex(-3.869/t913)
     |         *(1.+.108*t913+1.68*t923+1.26*t9+.551*t943+1.06*t953)

C.......Li6(p,g)Be7................(Caughlan-Fowler 1988)
      f(22)  = 6.69e+5*t9c56*t9m32*ex(-8.413/t9c13)

C60-----PROTON, ALPHA REACTIONS-------------------------------

C.......Li6(p,a)He3................(Caughlan-Fowler 1988)
      f(23)  = 3.73e+10*t9m23*ex(-8.413/t913-(t9/5.50)**2)
     |         *(1.+.050*t913-.061*t923-.021*t9+.006*t943+.005*t953)
     |         + 1.33e+10*t9m32*ex(-17.763/t9)
     |         + 1.29e+09*t9m1*ex(-21.820/t9)

C.......Li7(p,a)He4................(Smith-Kawano-Malaney 1992)
      f(24)  = 1.096e+9*t9m23*ex(-8.472/t913)
     |         - 4.830e+8*t9d56*t9m32*ex(-8.472/t9d13)
     |         + 1.06e+10*t9m32*ex(-30.442/t9)
     |         + 1.56e+5*t9m23*ex((-8.472/t913)-(t9/1.696)**2)
     |           *(1.+.049*t913-2.498*t923+.860*t9+3.518*t943+3.08*t953)
     |         + 1.55e+6*t9m32*ex(-4.478/t9)

C70-----ALPHA, PHOTON REACTIONS-------------------------------

C.......H2(a,g)Li6.................(Caughlan-Fowler 1988)
      f(25)  = 3.01e+01*t9m23*ex(-7.423/t913)
     |         *(1.+.056*t913-4.85*t923+8.85*t9-.585*t943-.584*t953)
     |         + 8.55e+1*t9m32*ex(-8.228/t9)

C.......H3(a,g)Li7.................(Smith-Kawano-Malaney 1992)
      f(26)  = 3.032e+5*t9m23*ex(-8.090/t913)
     |         *(1.+.0516*t913+.0229*t923+8.28e-3*t9
     |             -3.28e-4*t943-3.01e-4*t953)
     |         + 5.109e+5*t9e56*t9m32*ex(-8.068/t9e13)

C.......He3(a,g)Be7................(Smith-Kawano-Malaney 1992)
      f(27)  = 4.817e+6*t9m23*ex(-14.964/t913)
     |         *(1.+.0325*t913-1.04e-3*t923-2.37e-4*t9
     |             -8.11e-5*t943-4.69e-5*t953)
     |         + 5.938e+6*t9f56*t9m32*ex(-12.859/t9f13)

C80-----DEUTERIUM, NEUTRON AND DEUTERIUM, PROTON REACTIONS-------------

C.......H2(d,n)He3.................(Smith-Kawano-Malaney 1992)
      f(28)  = 3.95e+8*t9m23*ex(-4.259/t913)
     |         *(1.+.098*t913+.765*t923+.525*t9+9.61e-3*t943+.0167*t953)

C.......H2(d,p)H3..................(Smith-Kawano-Malaney 1992)
      f(29)  = 4.17e+8*t9m23*ex(-4.258/t913)
     |         *(1.+.098*t913+.518*t923+.355*t9-.010*t943-.018*t953)

C.......H3(d,n)He4.................(Smith-Kawano-Malaney 1992)
      f(30)  = 1.063e+11*t9m23*ex(-4.559/t913-(t9/.0754)**2)
     |         *(1.+.092*t913-.375*t923-.242*t9+33.82*t943+55.42*t953)
     |         + 8.047e+8*t9m23*ex(-0.4857/t9)

C.......He3(d,p)He4................(Smith-Kawano-Malaney 1992)
      f(31)  = 5.021e+10*t9m23*ex(-7.144/t913-(t9/.270)**2)
     |         *(1.+.058*t913+.603*t923+.245*t9+6.97*t943+7.19*t953)
     |         + 5.212e+8/t912*ex(-1.762/t9)

C90-----THREE PARTICLE REACTIONS------------------------------

C.......He3(He3,2p)He4.............(Caughlan-Fowler 1988)
      f(32)  = 6.04e+10*t9m23*ex(-12.276/t913)
     |         *(1.+.034*t913-.522*t923-.124*t9+.353*t943+.213*t953)

C.......Li7(d,na)He4...............(Caughlan-Fowler 1988)
      f(33)  = 2.92e+11*t9m23*ex(-10.259/t913)

C.......Be7(d,pa)He4...............(Caughlan-Fowler 1988)
      f(34)  = 1.07e+12*t9m23*ex(-12.428/t913)

      RETURN

C-------REFERENCES--------------------------------------
C     Smith, M., Kawano, L.H., and Malaney, R.A., 1992, submitted to Ap. J.
C     Malaney, R.A., and Fowler, W.A., 1989, Astrophys. J., 345, L5.
C     Caughlan, G.R., and Fowler, W.A., 1988, Atomic Data and Nuclear Data 
C       Tables, 40, 283.
C     Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247.

      END



C===============IDENTIFICATION DIVISION====================

      SUBROUTINE rate3

C-------LINKAGES.
C     CALLED BY - [subroutine] derivs
C     CALLS     - [function] ex

C-------REMARKS.
C     Generates rate coefficients for reactions involving nuclides 
C     up to A = 18.

C-------PARAMETER.
      PARAMETER (nrec=88)          !Number of nuclear reactions.
      PARAMETER (nnuc=26)          !Number of nuclides in calculation.

C-------COMMON AREAS.
      COMMON /rates/  f,r(nrec)           !Reaction rates.
      COMMON /evolp1/ t9,hv,phie,y(nnuc)  !Evolution parameters.


C=================DECLARATION DIVISION=====================

C-------REACTION RATES.
      REAL    f(nrec)              !Forward reaction rate coefficients.

C-------EVOLUTION PARAMETER.
      REAL    t9                   !Temperature of photons (units of 10**9 K).


C==================PROCEDURE DIVISION======================

C10-----TEMPERATURE FACTORS--------------------------------

      t913  = t9**(.33333333)      !t9**(1/3)
      t923  = t913*t913            !t9**(2/3)
      t943  = t923*t923            !t9**(4/3)
      t953  = t9*t923              !t9**(5/3)
      t912  = sqrt(t9)             !t9**(1/2)
      t932  = t9*t912              !t9**(3/2)
      t915  = t9**(.2)             !t9**(1/5)
      t954  = t9**(1.25)           !t9**(5/4)
      t9m1  = 1.0/t9               !t9**(-1)
      t9m23 = 1.0/t923             !t9**(-2/3)
      t9m32 = 1.0/t932             !t9**(-3/2)
      t9m34 = sqrt(t9m32)          !t9**(-3/4)
      t9m15 = 1.0/t915             !t9**(-1/5)
      t9m54 = 1.0/t954             !t9**(-5/4)
      t9a   = t9/(1.+t9/15.1)      !For reaction 53.
      t9a13 = t9a**(.3333333)      !t9a**(1/3)
      t9a56 = t9a**(.8333333)      !t9a**(5/6)

C20-----NEUTRON, PHOTON REACTIONS-----------------------------

C.......Li7(n,g)Li8................(Wagoner 1969)
      f(35)  = 4.90e+3 + 9.96e+3*t9m32*ex(-2.62/t9)

C.......B10(n,g)B11................(Wagoner 1969)
      f(36)  = 6.62e+4

C.......B11(n,g)B12................(Malaney-Fowler 1989)
      f(37)  = 7.29e+2 + 2.40e+3*t9m32*ex(-0.223/t9)

C30-----NEUTRON, PROTON REACTIONS-----------------------------

C.......C11(n,p)B11................(Caughlan-Fowler 1988)
      f(38)  = 1.69e+8*(1.-.048*t912+.010*t9)

C40-----NEUTRON, ALPHA REACTIONS------------------------------

C.......B10(n,a)Li7................(Caughlan-Fowler 1988)
      f(39)  = 5.07e+8

C50-----PROTON, PHOTON REACTIONS------------------------------

C.......Be7(p,g)B8.................(Caughlan-Fowler 1988)
      f(40)  = 3.11e+5*t9m23*ex(-10.262/t913)
     |         + 2.53e+3*t9m32*ex(-7.306/t9)

C.......Be9(p,g)B10................(Caughlan-Fowler 1988)
      f(41)  = 1.33e+7*t9m23*ex(-10.359/t913-(t9/.846)**2)
     |         *(1.+.040*t913+1.52*t923+.428*t9+2.15*t943+1.54*t953)
     |         + 9.64e+4*t9m32*ex(-3.445/t9)
     |         + 2.72e+6*t9m32*ex(-10.620/t9)

C.......B10(p,g)C11................(Caughlan-Fowler 1988)
      f(42)  = 4.61e+5*t9m23*ex(-12.062/t913-(t9/4.402)**2)
     |         *(1.+.035*t913+.426*t923+.103*t9+.281*t943+.173*t953)
     |         + 1.93e+5*t9m32*ex(-12.041/t9)
     |         + 1.14e+4*t9m32*ex(-16.164/t9)

C.......B11(p,g)C12................(Caughlan-Fowler 1988)
      f(43)  = 4.62e+7*t9m23*ex(-12.095/t913-(t9/.239)**2)
     |         *(1.+.035*t913+3.00*t923+.723*t9+9.91*t943+6.07*t953)
     |         + 7.89e+3*t9m32*ex(-1.733/t9)
     |         + 9.68e+4*t9m15*ex(-5.617/t9)

C.......C11(p,g)N12................(Caughlan-Fowler 1988)
      f(44)  = 4.24e+4*t9m23*ex(-13.658/t913-(t9/1.627)**2)
     |         *(1.+.031*t913+3.11*t923+.665*t9+4.61*t943+2.50*t953)
     |         + 8.84e+3*t9m32*ex(-7.021/t9)

C60-----PROTON, NEUTRON REACTIONS-----------------------------

C.......B12(p,n)C12................(Wagoner 1969)
      f(45)  = 4.02e+11*t9m23*ex(-12.12/t913)

C70-----PROTON, ALPHA REACTIONS-------------------------------

C.......Be9(p,a)Li6................(Caughlan-Fowler 1988)
      f(46)  = 2.11e+11*t9m23*ex(-10.359/t913-(t9/.520)**2)
     |         *(1.+.040*t913+1.09*t923+.307*t9+3.21*t943+2.30*t953)
     |         + 4.51e+8*t9m1*ex(-3.046/t9)
     |         + 6.70e+8*t9m34*ex(-5.160/t9)

C.......B10(p,a)Be7................(Caughlan-Fowler 1988)
      f(47)  = 1.26e+11*t9m23*ex(-12.062/t913-(t9/4.402)**2)
     |         *(1.+.035*t913-.498*t923-.121*t9+.300*t943+.184*t953)
     |         + 2.59e+9*t9m1*ex(-12.260/t9)

C.......B12(p,a)Be9................(Wagoner 1969)
      f(48)  = 2.01e+11*t9m23*ex(-12.12/t913)

C80-----ALPHA, PHOTON REACTIONS-------------------------------

C.......Li6(a,g)B10................(Caughlan-Fowler 1988)
      f(49)  = 4.06e+6*t9m23*ex(-18.790/t913-(t9/1.326)**2)
     |         *(1.+.022*t913+1.54*t923+.239*t9+2.20*t943+.869*t953)
     |         + 1.91e+3*t9m32*ex(-3.484/t9)
     |         + 1.01e+4*t9m1*ex(-7.269/t9)

C.......Li7(a,g)B11................(Caughlan-Fowler 1988)
      f(50)  = 3.55e+7*t9m23*ex(-19.161/t913-(t9/4.195)**2)
     |         *(1.+.022*t913+.775*t923+.118*t9+.884*t943+.342*t953)
     |         + 3.33e+2*t9m32*ex(-2.977/t9)
     |         + 4.10e+4*t9m1*ex(-6.227/t9)

C.......Be7(a,g)C11................(Caughlan-Fowler 1988)
      f(51)  = 8.45e+7*t9m23*ex(-23.212/t913-(t9/4.769)**2)
     |         *(1.+.018*t913+.488*t923+.061*t9+.296*t943+.095*t953)
     |         + 1.25e+4*t9m32*ex(-6.510/t9)
     |         + 1.29e+5*t9m54*ex(-10.039/t9)

C90-----ALPHA, PROTON REACTIONS-------------------------------

C.......B8(a,p)C11.................(Wagoner 1969)
      f(52)  = 1.08e+15*t9m23*ex(-27.36/t913)

C100-------ALPHA, NEUTRON REACTIONS------------------------------

C.......Li8(a,n)B11................(Malaney-Fowler 1989)
      f(53)  = 8.62e+13*t9a56*t9m32*ex(-19.461/t9a13)

C.......Be9(a,n)C12................(Caughlan-Fowler 1988)
      f(54)  = 4.62e+13*t9m23*ex(-23.870/t913-(t9/.049)**2)
     |         *(1.+.017*t913+8.57*t923+1.05*t9+74.51*t943+23.15*t953)
     |         + 7.34e-5*t9m32*ex(-1.184/t9)
     |         + 2.27e-1*t9m32*ex(-1.834/t9)
     |         + 1.26e+5*t9m32*ex(-4.179/t9)
     |         + 2.40e+8*ex(-12.732/t9)

C110-------DEUTERIUM, NEUTRON AND DEUTERIUM, PROTON REACTIONS-------------

C.......Be9(d,n)B10................(original Wagoner code)
      f(55)  = 7.16e+8*t9m23*ex(6.44-12.6/t913)

C.......B10(d,p)B11................(original Wagoner code)
      f(56)  = 9.53e+8*t9m23*ex(7.30-14.8/t913)

C.......B11(d,n)C12................(original Wagoner code)
      f(57)  = 1.41e+9*t9m23*ex(7.40-14.8/t913)

C120-------THREE PARTICLE REACTIONS------------------------------

C.......He4(an,g)Be9...............(Caughlan-Fowler 1988)
      f(58)  = (2.59e-6/((1.+.344*t9)*t9**2))*ex(-1.062/t9)

C.......He4(2a,g)C12...............(Caughlan-Fowler 1988)
      f(59)  = 2.79e-8*t9m32*t9m32*ex(-4.4027/t9)
     |         + 1.35e-8*t9m32*ex(-24.811/t9)

C.......Li8(p,na)He4...............(original Wagoner code)
      f(60)  = 8.65e+9*t9m23*ex(-8.52/t913-(t9/2.53)**2)
     |         + 2.31e+9*t9m32*ex(-4.64/t9)

C.......B8(n,pa)He4................(original Wagoner code)
      f(61)  = 4.02e+8

C.......Be9(p,da)He4...............(Caughlan-Fowler 1988)
      f(62)  = 2.11e+11*t9m23*ex(-10.359/t913-(t9/.520)**2)
     |         *(1.+.040*t913+1.09*t923+.307*t9+3.21*t943+2.30*t953)
     |         + 5.79e+8*t9m1*ex(-3.046/t9)
     |         + 8.50e+8*t9m34*ex(-5.800/t9)

C.......B11(p,2a)He4...............(Caughlan-Fowler 1988)
      f(63)  = 2.20e+12*t9m23*ex(-12.095/t913-(t9/1.644)**2)
     |         *(1.+.034*t913+.140*t923+.034*t9+.190*t943+.116*t953)
     |         + 4.03e+6*t9m32*ex(-1.734/t9)
     |         + 6.73e+9*t9m32*ex(-6.262/t9)
     |         + 3.88e+9*t9m1*ex(-14.154/t9)

C.......C11(n,2a)He4...............(Wagoner 1969)
      f(64)  = 1.58e+8

      RETURN

C-------REFERENCES--------------------------------------
C     Malaney, R.A., and Fowler, W.A., 1989, Astrophys. J., 345, L5.
C     Caughlan, G.R., and Fowler, W.A., 1988, Atomic Data and Nuclear Data 
C       Tables, 40, 283.
C     Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247.

      END



C===============IDENTIFICATION DIVISION====================

      SUBROUTINE rate4

C-------LINKAGES.
C     CALLED BY - [subroutine] derivs
C     CALLS     - [function] ex

C-------REMARKS.
C     Generates rate coefficients for rest of reactions.

C-------PARAMETER.
      PARAMETER (nrec=88)          !Number of nuclear reactions.
      PARAMETER (nnuc=26)          !Number of nuclides in calculation.

C-------COMMON AREAS.
      COMMON /rates/  f,r(nrec)           !Reaction rates.
      COMMON /evolp1/ t9,hv,phie,y(nnuc)  !Evolution parameters.


C=================DECLARATION DIVISION=====================

C-------REACTION RATES.
      REAL    f(nrec)              !Forward reaction rate coefficients.

C-------EVOLUTION PARAMETER.
      REAL    t9                   !Temperature of photons (units of 10**9 K).


C==================PROCEDURE DIVISION======================

C10-----TEMPERATURE FACTORS--------------------------------

      t913  = t9**(.33333333)      !t9**(1/3)
      t923  = t913*t913            !t9**(2/3)
      t943  = t923*t923            !t9**(4/3)
      t953  = t9*t923              !t9**(5/3)
      t912  = sqrt(t9)             !t9**(1/2)
      t932  = t9*t912              !t9**(3/2)
      t935  = t9**(.6)             !t9**(3/5)
      t965  = t9**(1.2)            !t9**(6/5)
      t938  = t9**(.375)           !t9**(3/8)
      t9m13 = 1.0/t913             !t9**(1/3)
      t9m23 = 1.0/t923             !t9**(-2/3)
      t9m32 = 1.0/t932             !t9**(-3/2)
      t9m65 = 1.0/t965             !t9**(-6/5)
      t9a   = t9                   !For reaction 82.
     |          /(1.+4.78e-2*t9+7.56e-3*t953/(1.+4.78e-2*t9)**(2./3.))  
      t9a13 = t9a**(.33333333)     !t9a**(1/3)
      t9a56 = t9a**(.83333333)     !t9a**(5/6)
      t9b   = t9                   !For reaction 84.
     |          /(1.+7.76e-2*t9+2.64e-2*t953/(1.+7.76e-2*t9)**(2./3.))
      t9b13 = t9b**(.33333333)     !t9b**(1/3)
      t9b56 = t9b**(.83333333)     !t9b**(5/6)

C20-----NEUTRON, PHOTON REACTIONS-----------------------------

C.......C12(n,g)C13................(Wagoner 1969)
      f(65)  = 4.50e+2

C.......C13(n,g)C14................(Wagoner 1969)
      f(66)  = 1.19e+2 + 2.38e+5*t9m32*ex(-1.67/t9)

C.......N14(n,g)N15................(Wagoner 1969)
      f(67)  = 9.94e+3

C30-----NEUTRON, PROTON REACTIONS-----------------------------

C.......N13(n,p)C13................(Caughlan-Fowler 1988)
      f(68)  = 1.88e+8*(1.-.167*t912+.037*t9)

C.......N14(n,p)C14................(Caughlan-Fowler 1988)
      f(69)  = 2.39e+5*(1.+.361*t912+.502*t9)
     |         + 1.112e+8/t912*ex(-4.983/t9)

C.......O15(n,p)N15................(Caughlan-Fowler 1988)
      f(70)  = 3.50e+8*(1.+.452*t912-.191*t9)

C40-----NEUTRON, ALPHA REACTIONS------------------------------

C.......O15(n,a)C12................(Caughlan-Fowler 1988)
      f(71)  = 3.50e+7*(1.+.188*t912+.015*t9)

C50-----PROTON, PHOTON REACTIONS------------------------------

C.......C12(p,g)N13................(Caughlan-Fowler 1988)
      f(72)  = 2.04e+7*t9m23*ex(-13.690/t913-(t9/1.500)**2)
     |         *(1.+.030*t913+1.19*t923+.254*t9+2.06*t943+1.12*t953)
     |         + 1.08e+5*t9m32*ex(-4.925/t9)
     |         + 2.15e+5*t9m32*ex(-18.179/t9)

C.......C13(p,g)N14................(Caughlan-Fowler 1988)
      f(73)  = 8.01e+7*t9m23*ex(-13.717/t913-(t9/2.000)**2)
     |         *(1.+.030*t913+.958*t923+.204*t9+1.39*t943+.753*t953)
     |         + 1.21e+6*t9m65*ex(-5.701/t9)

C.......C14(p,g)N15................(Caughlan-Fowler 1988)
      f(74)  = 6.80e+6*t9m23*ex(-13.741/t913-(t9/5.721)**2)
     |         *(1.+.030*t913+.503*t923+.107*t9+.213*t943+.115*t953)
     |         + 5.36e+3*t9m32*ex(-3.811/t9)
     |         + 9.82e+4*t9m13*ex(-4.739/t9)

C.......N13(p,g)O14................(Caughlan-Fowler 1988)
      f(75)  = 4.04e+7*t9m23*ex(-15.202/t913-(t9/1.191)**2)
     |         *(1.+.027*t913-.803*t923-.154*t9+5.00*t943+2.44*t953)
     |         + 2.43e+5*t9m32*ex(-6.348/t9)

C.......N14(p,g)O15................(Caughlan-Fowler 1988)
      f(76)  = 4.90e+7*t9m23*ex(-15.228/t913-(t9/3.294)**2)
     |         *(1.+.027*t913-.778*t923-.149*t9+.261*t943+.127*t953)
     |         + 2.37e+3*t9m32*ex(-3.011/t9)
     |         + 2.19e+4*ex(-12.530/t9)

C.......N15(p,g)O16................(Caughlan-Fowler 1988)
      f(77)  = 9.78e+8*t9m23*ex(-15.251/t913-(t9/.450)**2)
     |         *(1.+.027*t913+.219*t923+.042*t9+6.83*t943+3.32*t953)
     |         + 1.11e+4*t9m32*ex(-3.328/t9)
     |         + 1.49e+4*t9m32*ex(-4.665/t9)
     |         + 3.80e+6*t9m32*ex(-11.048/t9)

C60-----PROTON, ALPHA REACTIONS-------------------------------

C.......N15(p,a)C12................(Caughlan-Fowler 1988)
      f(78)  = 1.08e+12*t9m23*ex(-15.251/t913-(t9/.522)**2)
     |         *(1.+.027*t913+2.62*t923+.501*t9+5.36*t943+2.60*t953)
     |         + 1.19e+8*t9m32*ex(-3.676/t9)
     |         + 5.41e+8/t912*ex(-8.926/t9)
     |         + 4.72e+7*t9m32*ex(-7.721/t9)
     |         + 2.20e+8*t9m32*ex(-11.418/t9)

C70-----ALPHA, PHOTON REACTIONS-------------------------------

C.......C12(a,g)O16................(Caughlan-Fowler 1988)
      f(79)  = 1.04e+8/t9**2*ex(-32.120/t913-(t9/3.496)**2)
     |         /(1.+.0489*t9m23)**2
     |         + 1.76e+8/(t9)**2/(1.+.2654*t9m23)**2*ex(-32.120/t913)
     |         + 1.25e+3*t9m32*ex(-27.499/t9)
     |         + 1.43e-2*(t9)**5*ex(-15.541/t9)

C80-----ALPHA, PROTON REACTIONS-------------------------------

C.......B10(a,p)C13................(Wagoner 1969)
      f(80)  = 9.60e+14*t9m23*ex(-27.99/t913)

C.......B11(a,p)C14................(Caughlan-Fowler 1988)
      f(81)  = 5.37e+11*t9m23*ex(-28.234/t913-(t9/0.347)**2)
     |         *(1.+.015*t913+5.575*t923+.576*t9+15.888*t943+4.174*t953)
     |         + 5.44e-3*t9m32*ex(-2.827/t9)
     |         + 3.36e+2*t9m32*ex(-5.178/t9)
     |         + 5.32e+6/t938*ex(-11.617/t9)

C.......C11(a,p)N14................(Caughlan-Fowler 1988)
      f(82)  = 7.15e+15*t9a56*t9m32*ex(-31.883/t9a13)

C.......N12(a,p)O15................(Caughlan-Fowler 1988)
      f(83)  = 5.59e+16*t9m23*ex(-35.60/t913)

C.......N13(a,p)O16................(Caughlan-Fowler 1988)
      f(84)  = 3.23e+17*t9b56*t9m32*ex(-35.829/t9b13)

C90-----ALPHA, NEUTRON REACTIONS------------------------------

C.......B10(a,n)N13................(Caughlan-Fowler 1988)
      f(85)  = 1.20e+13*t9m23*ex(-27.989/t913-(t9/9.589)**2)

C.......B11(a,n)N14................(Caughlan-Fowler 1988)
      f(86)  = 6.97e+12*t9m23*ex(-28.234/t913-(t9/0.140)**2)
     |         *(1.+.015*t913+8.115*t923+.838*t9+39.804*t943
     |             +10.456*t953)
     |         + 1.79e+0*t9m32*ex(-2.827/t9)
     |         + 1.71e+3*t9m32*ex(-5.178/t9)
     |         + 4.49e+6*t935*ex(-8.596/t9)

C.......B12(a,n)N15................(Wagoner 1969)
      f(87)  = 3.04e+15*t9m23*ex(-28.45/t913)

C.......C13(a,n)O16................(Caughlan-Fowler 1988)
      f(88)  = 6.77e+15*t9m23*ex(-32.329/t913-(t9/1.284)**2)
     |         *(1.+.013*t913+2.04*t923+.184*t9)
     |         + 3.82e+5*t9m32*ex(-9.373/t9)
     |         + 1.41e+6*t9m32*ex(-11.873/t9)
     |         + 2.00e+9*t9m32*ex(-20.409/t9)
     |         + 2.92e+9*t9m32*ex(-29.283/t9)

      RETURN

C-------REFERENCES--------------------------------------
C     Caughlan, G.R., and Fowler, W.A., 1988, Atomic Data and Nuclear Data
C       Tables, 40, 283.
C     Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247.

      END



C===============IDENTIFICATION DIVISION====================

      BLOCK DATA   

C-------PARAMETERS.
      PARAMETER (nrec=88)          !Number of nuclear reactions.
      PARAMETER (nnuc=26)          !Number of nuclides in calculation.

C-------COMMON AREAS.
      COMMON /recpr0/ reacpr                         !Reaction parameter values.
      COMMON /compr0/ cy0,ct0,t9i0,t9f0,ytmin0,inc0  !Default comp parameters.
      COMMON /modpr0/ c0,cosmo0,xi0                  !Default model parameters.
      COMMON /varpr0/ dt0,eta0                       !Default variationl params.
      COMMON /nucdat/ am,zm,dm                       !Nuclide data.


C=================DECLARATION DIVISION=====================

C-------REACTION PARAMETERS VALUES.
      REAL    reacpr(nrec,8)       !Reaction parameters.

C-------DEFAULT COMPUTATION PARAMETERS.
      REAL    cy0                  !Default time step limiting constant.
      REAL    ct0                  !Default time step limiting constant.
      REAL    t9i0                 !Default initial temperature (in 10**9 K).
      REAL    t9f0                 !Default final temperature (in 10**9 K).
      REAL    ytmin0               !Default smallest abundances allowed.
      INTEGER inc0                 !Default accumulation increment.

C-------DEFAULT MODEL PARAMETERS.
      REAL    c0(3)               !c0(1) is default variation of grav constant.
     |                             !c0(2) is default neutron half-life.
     |                             !c0(3) is default number of neutrinos.
      REAL    cosmo0               !Default cosmological constant.
      REAL    xi0(3)               !Default neutrino degeneracy parameters.

C-------DEFAULT VARIATIONAL PARAMETERS.
      REAL    dt0                  !Default initial time step.
      REAL    eta0                 !Default baryon-to-photon ratio.      

C-------NUCLIDE DATA.
      REAL    am(nnuc)             !Atomic number of nuclide.
      REAL    zm(nnuc)             !Charge of nuclide.
      REAL    dm(nnuc)             !Mass excess of nuclide.


C=====================DATA DIVISION========================

C    Nuclide and corresponding number
C    --------------------
C    1) N         7) Li6      13) B10      19) C13      25) O15
C    2) P         8) Li7      14) B11      20) N13      26) O16
C    3) H2        9) Be7      15) C11      21) C14
C    4) H3       10) Li8      16) B12      22) N14
C    5) He3      11) B8       17) C12      23) O14
C    6) He4      12) Be9      18) N12      24) N15

C--------NUCLIDE DATA.
      DATA am /1.,1.,2.,3.,3.,4.,6.,7.,7.,8.,8.,9.,10.,11.,11.,12.,
     |         12.,12.,13.,13.,14.,14.,14.,15.,15.,16./
      DATA zm /0.,1.,1.,1.,2.,2.,3.,3.,4.,3.,5.,4.,5.,5.,6.,5.,
     |         6.,7.,6.,7.,6.,7.,8.,7.,8.,8./
      DATA dm /.008665,.007825,.014102,.016050,.016030,.002603,.015125,
     |         .016004,.016929,.022487,.024609,.012186,.012939,.009305,
     |         .011432,.014354,.000000,.018641,.003354,.005738,.003242,
     |         .003074,.008597,.000108,.003070,-.005085/

C-------REACTION RATE COEFFICIENTS (Ref 1).
      DATA ((reacpr(i,j),j=1,8),i=1,11) /
C              reac# type n1 n2 n3 n4 rev-coeff q-value
C              ----  ---- -- -- -- -- ------ -------
     |             1.,1., 1.,0.,0., 2., 0.0  ,   0.0 ,     !N->P         
     |             2.,1., 4.,0.,0., 5., 0.0  ,   0.0 ,     !H3->He3
     |             3.,4.,10.,0.,0., 6., 0.0  ,   0.0 ,     !Li8->2He4
     |             4.,1.,16.,0.,0.,17., 0.0  ,   0.0 ,     !B12->C12
     |             5.,1.,21.,0.,0.,22., 0.0  ,   0.0 ,     !C14->N14
     |             6.,4.,11.,0.,0., 6., 0.0  ,   0.0 ,     !B8->2He4
     |             7.,1.,15.,0.,0.,14., 0.0  ,   0.0 ,     !C11->B11
     |             8.,1.,18.,0.,0.,17., 0.0  ,   0.0 ,     !N12->C12
     |             9.,1.,20.,0.,0.,19., 0.0  ,   0.0 ,     !N13->C13
     |            10.,1.,23.,0.,0.,22., 0.0  ,   0.0 ,     !O14->N14
     |            11.,1.,25.,0.,0.,24., 0.0  ,   0.0 /     !O15->N15
      DATA ((reacpr(i,j),j=1,8),i=12,22) /
C              reac# type n1 n2 n3 n4 rev-coeff q-value
C              ----  ---- -- -- -- -- ------ -------
     |            12.,2., 2.,1.,0., 3., 0.471,  25.82,     !H(n,g)H2
     |            13.,2., 3.,1.,0., 4., 1.63 ,  72.62,     !H2(n,g)H3
     |            14.,2., 5.,1.,0., 6., 2.61 , 238.81,     !He3(n,g)He4
     |            15.,2., 7.,1.,0., 8., 1.19 ,  84.17,     !Li6(n,g)Li7
     |            16.,3., 5.,1.,2., 4., 1.002,   8.863,    !He3(n,p)H3
     |            17.,3., 9.,1.,2., 8., 0.998,  19.081,    !Be7(n,p)Li7
     |            18.,3., 7.,1.,4., 6., 1.070,  55.494,    !Li6(n,a)H3
     |            19.,5., 9.,1.,0., 6., 4.70 , 220.39,     !Be7(n,a)He4
     |            20.,2., 3.,2.,0., 5., 1.63 ,  63.750,    !H2(p,g)He3
     |            21.,2., 4.,2.,0., 6., 2.61 , 229.932,    !H3(p,g)He4
     |            22.,2., 7.,2.,0., 9., 1.19 ,  65.054/    !Li6(p,g)Be7
      DATA ((reacpr(i,j),j=1,8),i=23,33) /
C              reac# type n1 n2 n3 n4 rev-coeff q-value
C              ----  ---- -- -- -- -- ------ -------
     |            23.,3., 7.,2.,5., 6., 1.07 ,  46.631,    !Li6(p,a)He3
     |            24.,5., 8.,2.,0., 6., 4.69 , 201.291,    !Li7(p,a)He4
     |            25.,2., 6.,3.,0., 7., 1.53 ,  17.118,    !H2(a,p)Li6
     |            26.,2., 6.,4.,0., 8., 1.11 ,  28.640,    !H3(a,p)Li7
     |            27.,2., 6.,5.,0., 9., 1.11 ,  18.423,    !He3(a,p)Be7
     |            28.,6., 3.,0.,1., 5., 1.73 ,  37.935,    !H2(d,p)He3
     |            29.,6., 3.,0.,2., 4., 1.73 ,  46.798,    !H2(d,n)H3
     |            30.,3., 4.,3.,1., 6., 5.54 , 204.117,    !H3(d,n)He4
     |            31.,3., 5.,3.,2., 6., 5.55 , 212.980,    !He3(d,p)He4
     |            32.,11.,5.,0.,2., 6., 3.39 , 149.230,    !He3(He3,2p)He4
     |            33.,9., 8.,3.,1., 6., 9.95 , 175.476/    !Li7(d,na)He4
      DATA ((reacpr(i,j),j=1,8),i=34,44) /
C              reac# type n1 n2 n3 n4 rev-coeff q-value
C              ----  ---- -- -- -- -- ------ -------
     |            34.,9., 9.,3.,2., 6., 9.97 , 194.557,    !Be7(d,pa)He4
     |            35.,2., 8.,1.,0.,10., 1.31 ,  23.59,     !Li7(n,g)Li8
     |            36.,2.,13.,1.,0.,14., 3.04 , 132.95,     !B10(n,g)B11
     |            37.,2.,14.,1.,0.,16., 2.34 ,  39.10,     !B11(n,g)B12
     |            38.,3.,15.,1.,2.,14., 1.002,  32.080,    !C11(n,p)B11
     |            39.,3.,13.,1.,6., 8., 0.758,  32.382,    !B10(n,a)Li7
     |            40.,2., 9.,2.,0.,11., 1.30 ,   1.595,    !Be7(p,g)B8
     |            41.,2.,12.,2.,0.,13., 0.973,  76.427,    !Be9(p,g)B10
     |            42.,2.,13.,2.,0.,15., 3.03 , 100.840,    !B10(p,g)C11
     |            43.,2.,14.,2.,0.,17., 7.01 , 185.173,    !B11(p,g)C12
     |            44.,2.,15.,2.,0.,18., 2.33 ,   6.975/    !C11(p,g)N12
      DATA ((reacpr(i,j),j=1,8),i=45,55) /
C              reac# type n1 n2 n3 n4 rev-coeff q-value
C              ----  ---- -- -- -- -- ------ -------
     |            45.,3.,16.,2.,1.,17., 3.00 , 146.08,     !B12(p,n)C12
     |            46.,3.,12.,2.,6., 7., 0.618,  24.674,    !Be9(p,a)Li6
     |            47.,3.,13.,2.,6., 9., 0.754,  13.301,    !B10(p,a)Be7
     |            48.,3.,16.,2.,6.,12., 0.292,  79.89,     !B12(p,a)Be9
     |            49.,2., 7.,6.,0.,13., 1.58 ,  51.753,    !Li6(a,g)B10
     |            50.,2., 8.,6.,0.,14., 4.02 , 100.538,    !Li7(a,g)B11
     |            51.,2., 9.,6.,0.,15., 4.02 ,  87.539,    !Be7(a,g)C11
     |            52.,3.,11.,6.,2.,15., 3.08 ,  86.00,     !B8(a,p)C11
     |            53.,3.,10.,6.,1.,14., 3.07 ,  76.96,     !Li8(a,n)B11   
     |            54.,3.,12.,6.,1.,17.,10.3  ,  66.160,    !Be9(a,n)C12
     |            55.,3.,12.,3.,1.,13., 2.07 ,  50.63/     !Be9(d,n)B10
      DATA ((reacpr(i,j),j=1,8),i=56,66) /
C              reac# type n1 n2 n3 n4 rev-coeff q-value
C              ----  ---- -- -- -- -- ------ -------
     |            56.,3.,13.,3.,2.,14., 6.44 , 107.13,     !B10(d,p)B11
     |            57.,3.,14.,3.,1.,17.,14.9  , 159.36,     !B11(d,n)C12
     |            58.,8., 6.,1.,0.,12., 0.584,  18.260,    !He4(an,g)Be9
     |            59.,7., 6.,0.,0.,17., 2.00 ,  84.420,    !He4(2a,g)C12
     |            60.,9.,10.,2.,1., 6., 3.58 , 177.73,     !Li8(p,na)He4
     |            61.,9.,11.,1.,2., 6., 3.58 , 218.82,     !B8(n,pa)He4
     |            62.,9.,12.,2.,3., 6., 0.807,   7.555,    !Be9(p,da)He4
     |            63.,10.,14.,2.,0.,6., 3.50 , 100.753,    !B11(p,2a)Be4
     |            64.,10.,15.,1.,0.,6., 3.49 , 132.83,     !C11(n,2a)He4
     |            65.,2.,17.,1.,0.,19., 0.886,  57.41,     !C12(n,g)C13
     |            66.,2.,19.,1.,0.,21., 3.58 ,  94.88/     !C13(n,g)C14
      DATA ((reacpr(i,j),j=1,8),i=67,77) /
C              reac# type n1 n2 n3 n4 rev-coeff q-value
C              ----  ---- -- -- -- -- ------ -------
     |            67.,2.,22.,1.,0.,24., 2.71 , 125.74,     !N14(n,g)N15
     |            68.,3.,20.,1.,2.,19., 1.002,  34.846,    !N13(n,p)C13
     |            69.,3.,22.,1.,2.,21., 3.003,   7.263,    !N14(n,p)C14
     |            70.,3.,25.,1.,2.,24., 1.002,  41.037,    !O15(n,p)N15
     |            71.,3.,25.,1.,6.,17., 0.709,  98.661,    !O15(n,a)C12
     |            72.,2.,17.,2.,0.,20., 0.884,  22.553,    !C12(p,g)N13
     |            73.,2.,19.,2.,0.,22., 1.19 ,  87.621,    !C13(p,g)N14
     |            74.,2.,21.,2.,0.,24., 0.900, 118.452,    !C14(p,g)N15
     |            75.,2.,20.,2.,0.,23., 3.57 ,  53.706,    !N13(p,g)O14
     |            76.,2.,22.,2.,0.,25., 2.70 ,  84.678,    !N14(p,g)O15
     |            77.,2.,24.,2.,0.,26., 3.62 , 140.734/    !N15(p,g)O16
      DATA ((reacpr(i,j),j=1,8),i=78,88) /
C              reac# type n1 n2 n3 n4 rev-coeff q-value
C              ----  ---- -- -- -- -- ------ -------
     |            78.,3.,24.,2.,6.,17., 0.706,  57.623,    !N15(p,a)C12
     |            79.,2.,17.,6.,0.,26., 5.13 ,  83.111,    !C12(a,g)O16
     |            80.,3.,13.,6.,2.,19., 9.36 ,  47.16,     !B10(a,p)C13  
     |            81.,3.,14.,6.,2.,21.,11.0  ,   9.098,    !B11(a,p)C14  
     |            82.,3.,15.,6.,2.,22., 3.68 ,  33.915,    !C11(a,p)N14  
     |            83.,3.,18.,6.,2.,25., 4.26 , 111.87,     !N12(a,p)O15  
     |            84.,3.,20.,6.,2.,26., 5.81 ,  60.557,    !N13(a,p)O16  
     |            85.,3.,13.,6.,1.,20., 9.34 ,  12.287,    !B10(a,n)N13  
     |            86.,3.,14.,6.,1.,22., 3.67 ,   1.835,    !B11(a,n)N14  
     |            87.,3.,16.,6.,1.,24., 4.25 ,  88.47,     !B12(a,n)N15  
     |            88.,3.,19.,6.,1.,26., 5.79 ,  25.711/    !C13(a,n)O16  

C-------DEFAULT COMPUTATION PARAMETERS.
      DATA cy0    /.300/           !Default time step limiting constant.
      DATA ct0    /.030/           !Default time step limiting constant.
      DATA t9i0   /1.00e+02/       !Default initial temperature.
      DATA t9f0   /1.00e-02/       !Default final temperature.
      DATA ytmin0 /1.00e-25/       !Default smallest abundances allowed.
      DATA inc0   /30/             !Default accumulation increment.

C--------DEFAULT MODEL PARAMETERS.
      DATA c0     /1.00,887.,3.0/!Default variation of 3 parameters.
      DATA cosmo0 /0.00/           !Default cosmological constant.
      DATA xi0    /0.00,0.00,0.00/ !Default neutrino degeneracy parameter.

C--------DEFAULT VARIATIONAL PARAMETERS.
      DATA dt0    /1.00e-04/       !Default initial time step.
      DATA eta0   /3.000e-10/      !Default baryon-to-photon ratio.

      END

C===============IDENTIFICATION DIVISION====================

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

C-------PARAMETERS.
      PARAMETER (nrec=88)          !Number of nuclear reactions.
      PARAMETER (nvar=29)          !Number of variables to be evolved.
      PARAMETER (nnuc=26)          !Number of nuclides in calculation.
      PARAMETER (itmax=40)         !Maximum # of lines to be printed.

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
      REAL    reacpr(nrec,8)       !Reaction parameters.

C-------REACTION PARAMETER NAMES.
      INTEGER iform(nrec)          !Reaction type code (1-11).
      INTEGER ii(nrec)             !Incoming nuclide type (1-26).
      INTEGER jj(nrec)             !Incoming light nuclide type (1-6).
      INTEGER kk(nrec)             !Outgoing light nuclide type (1-6).
      INTEGER ll(nrec)             !Outgoing nuclide type (1-26).
      REAL    rev(nrec)            !Reverse reaction coefficient.
      REAL    q9(nrec)             !Energy released in reaction (in 10**9 K).

C-------REACTION RATES.
      REAL    f(nrec)              !Forward reaction rate coefficients.
      REAL    r(nrec)              !Reverse reaction rate coefficients.

C-------EVOLUTION PARAMETERS.
      REAL    t9                   !Temperature of photons (units of 10**9 K).
      REAL    hv                   !Defined by hv = M(atomic)n(baryon)/t9**3.
      REAL    phie                 !Chemical potential of electron.
      REAL    y(nnuc)              !Relative number abundances.

C-------EVOLUTION PARAMETERS (DERIVATIVES).
      REAL    dt9                  !Change in temperature.
      REAL    dhv                  !Change in hv.
      REAL    dphie                !Change in chemical potential.
      REAL    dydt(nnuc)           !Change in relative number abundances.

C-------EVOLUTION PARAMETERS (ORIGINAL VALUES).
      REAL    y0(nnuc)             !Rel # abundances at end of 1st R-K loop.

C-------DEFAULT COMPUTATION PARAMETERS.
      REAL    cy0                  !Default cy.
      REAL    ct0                  !Default ct.
      REAL    t9i0                 !Default t9i.
      REAL    t9f0                 !Default t9f.
      REAL    ytmin0               !Default ytmin.
      INTEGER inc0                 !Default accumulation increment.

C-------COMPUTATION PARAMETERS.
      REAL    cy                   !Time step limiting constant on abundances.
      REAL    ct                  !Time step limiting constant on temperature.
      REAL    t9i                  !Initial temperature (in 10**9 K).
      REAL    t9f                  !Final temperature (in 10**9 K).
      REAL    ytmin                !Smallest abundances allowed.
      INTEGER inc                  !Accumulation increment.

C-------DEFAULT MODEL PARAMETERS.
      REAL    c0(3)                !Default c.
      REAL    cosmo0               !Default cosmological constant.
      REAL    xi0(3)               !Default neutrino degeneracy parameters.

C-------EARLY UNIVERSE MODEL PARAMETERS.
      REAL    g                    !Gravitational constant.
      REAL    tau                  !Neutron lifetime (sec).
      REAL    xnu                  !Number of neutrino species.
      REAL    c(3)               !c(1) is variation of gravitational constant.
     |                             !c(2) is neutron half-life (min).
     |                             !c(3) is number of neutrino species.
      REAL    cosmo                !Cosmological constant.
      REAL    xi(3)                !Neutrino degeneracy parameters.
     |                             !xi(1) is e neutrino degeneracy parameter.
     |                             !xi(2) is m neutrino degeneracy parameter.
     |                             !xi(3) is t neutrino degeneracy parameter.

C-------DEFAULT VARIATIONAL PARAMETERS.
      REAL    dt0                  !Default initial time step.
      REAL    eta0                 !Default baryon-to-photon ratio.

C-------VARIATIONAL PARAMETERS.
      REAL    dt1                  !Initial time step.
      REAL    eta1                 !Baryon-to-photon ratio.

C-------TIME VARIABLES.
      REAL    t                    !Time.
      REAL    dt                   !Time step.
      REAL    dlt9dt               !(1/t9)*d(t9)/d(t).

C-------DYNAMIC VARIABLES.
      REAL    thm(14)              !Thermodynamic variables (energy densities).
      REAL    hubcst               !Expansion rate of the universe.

C-------ENERGY DENSITIES.
      REAL    rhone0               !Initial electron neutrino energy density.
      REAL    rhob0                !Initial baryon energy density.
      REAL    rhob                 !Baryon energy density.
      REAL    rnb                !Baryon energy density (ratio to init value).

C-------MATRIX COEFFICIENTS FOR LINEAR EQUATION.
      DOUBLE PRECISION a(nnuc,nnuc)!Relates y(t+dt) to y(t).
      REAL    b(nnuc)              !Contains y0 in inverse order.
      REAL    yx(nnuc)             !yy in reverse order.

C-------NUCLIDE DATA.
      REAL    am(nnuc)             !Atomic number of nuclide.
      REAL    zm(nnuc)             !Charge of nuclide.
      REAL    dm(nnuc)             !Mass excess of nuclide.

C-------EVALUATION OF FUNCTIONS bl,bm,bn.
      REAL    bl1,bl2,bl3,bl4,bl5  !Evaluation of function bl(z).
      REAL    bm1,bm2,bm3,bm4,bm5  !Evaluation of function bm(z).
      REAL    bn1,bn2,bn3,bn4,bn5  !Evaluation of function bn(z).

C-------EVALUATION OF MODIFIED BESSEL FUNCTIONS.
      REAL    bk0,bk1,bk2,bk3,bk4  !Values k0(r),k1(r),k2(r),k3(r),k4(r).

C-------FLAGS AND COUNTERS.
      INTEGER ltime                !Indicates if output buffer printed.
      INTEGER is                   !# total iterations for particular model.
      INTEGER ip                   !# iterations after outputing a line.
      INTEGER it                   !# times accumulated in output buffer.
      INTEGER mbad                 !Indicates if gaussian elimination failed.

C-------COMPUTATION LOCATION.
      INTEGER itime                !Time check.

C-------OUTPUT ARRAYS.
      REAL    xout(itmax,nnuc)     !Nuclide mass fractions.
      REAL    thmout(itmax,6)      !Thermodynamic variables.
      REAL    t9out(itmax)         !Temperature (in units of 10**9 K).
      REAL    tout(itmax)          !Time.
      REAL    dtout(itmax)         !Time step.
      REAL    etaout(itmax)        !Baryon to photon ratio.
      REAL    hubout(itmax)        !Expansion rate.

C-------NEUTRINO PARAMETERS.
      REAL    t9mev                !Temperature (in units of MeV).
      REAL    tnmev                !Neutrino temperature (in units of MeV).
      REAL    tnu                  !Neutrino temperature.
      REAL    cnorm                !Normalizing constant.
      REAL    rhonu                !Neutrino energy density.
      INTEGER nu                   !Type of neutrino.

C-------RUN OPTION.
      INTEGER irun                 !Run network size.
      INTEGER isize                !Number of nuclides in computation.
      INTEGER jsize                !Number of reactions in computation.

C-------OUTPUT FILE STATUS.
      INTEGER nout                 !Number of output requests.
      LOGICAL outfile              !Indicates if output file used.


C==================PROCEDURE DIVISION======================

C10-----OPEN FILE---------------------------------------

      IF (itime.eq.1) THEN         !Beginning of program.
        OPEN (unit=3, file='nucint.dat',  status='unknown')
      END IF

C20-----PRINTINTO FILE------------------------------------

      IF (itime.eq.8) THEN         !Right after a run.
        xout(it,8) = xout(it,8) + xout(it,9)  !Add beryllium to lithium.
        xout(it,5) = xout(it,5) + xout(it,4)  !Add tritium to helium-3.
        xout(it,6) = xout(it,6)-0.0025  
     |            !Radiative, coulomb, finite-temperature corrections (Ref 1).
        write(3,200) etaout(it),xout(it,3),
     |                xout(it,5),xout(it,6),xout(it,8)  
     |                             !Output eta, H2, He3, He4, and Li7.
 200    FORMAT (5(e13.5,' '))
      END IF

C30-----close FILE--------------------------------------

      IF (itime.eq.10) THEN        !End of program.
        CLOSE (unit=3)
      END IF
      RETURN

C-------REFERENCES--------------------------------------
C     1) D.A. Dicus, E.W. Kolb, A.M. Gleeson, E.C.G. Sudarshan, V.L. Teplitz,
C        M.S. Turner, Phys. Rev. D., 26,2694 (1982).

      END

