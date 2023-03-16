
//  *********************************************************************
//                       SUBROUTINE ESIa
//  *********************************************************************
void PenPhys::ESIa(double E_, double DELTA, double &DE, double &EP, double &CDT, double &ES, double &CDTS, int M, int &IZZ, int &ISH)
{
  //  Random sampling of inner-shell ionisation by electron impact.
  //
  //  Sternheimer-Liljequist GOS model.
  //
  //  Input arguments:
  //    E_ ....... electron energy (eV).
  //    M ....... material where electrons propagate.
  //    DELTA ... Fermi's density effect correction.
  //  Output arguments:
  //    DE ...... energy loss (eV).
  //    EP ...... energy of the scattered electron (eV).
  //    CDT ..... cosine of the polar scattering angle.
  //    ES ...... energy of the emitted secondary electron (eV).
  //    CDTS .... polar cosine of direction of the secondary electron.
  //    IZZ ..... atomic number of the element where ionisation has
  //             occurred.
  //    ISH ..... atomic electron shell that has been ionized.

  using namespace PENELOPE_mod;

  using namespace CECUTR;
////////////  using namespace CEGRID;
  using namespace CEIN;
  using namespace CESIAC;

  const double RREV = 1.0/REV;
  const double TREV = 2.0*REV;
  const double RTREV = 1.0/TREV;

  //  ****  Energy grid point.
  int JE;
  double PK = (CEGRID_.XEL-CEGRID_.DLEMP[CEGRID_.KE-1])*CEGRID_.DLFC;
  if(RAND(1.0) < PK)
  {
    JE = CEGRID_.KE+1;
  }
  else
  {
    JE = CEGRID_.KE;
  }

  //  ************  Selection of the active oscillator.

  double TST = RAND(2.0);
  int IO, JO, IT, IOSC;
  //  ****  Binary search.
  IO = 1;
  JO = NESI[M-1]+1;
  bool brkIt = false;
  while(!brkIt)
  {
    brkIt = true;
    IT = (IO+JO)/2;
    if(TST > ESIAC[M-1][JE-1][IT-1])
    {
      IO = IT;
    }
    else
    {
      JO = IT;
    }
    if(JO-IO > 1){ brkIt = false; continue;}
  }
  IOSC = IESI[M-1][IO-1];

  IZZ = KZ[M-1][IOSC-1];
  ISH = KS[M-1][IOSC-1];

  double UK = UI[M-1][IOSC-1];
  double WK = WRI[M-1][IOSC-1];
  double WTHR;
  if(UK > 1.0E-3)
  {
    WTHR = UK;
  }
  else
  {
    WTHR = WK;
  }

  if(E_ < WTHR+1.0E-6)
  {
    DE = UK;
    EP = E_-DE;
    CDT = 1.0;
    ES = 0.0;
    CDTS = 0.0;
    return;
  }

  //  ****  Trick: The resonance energy and the cutoff recoil energy of
  //        inner shells are varied to yield a smooth threshold.

  double WM = 3.0*WK-2.0*UK;
  double WKP, QKP;
  if(E_ > WM)
  {
    WKP = WK;
    QKP = UK;
  }
  else
  {
    WKP = (E_+2.0*UK)/3.0;
    QKP = UK*(E_/WM);
    WM = E_;
  }
  double EE = E_+UK;
  double WCMAX = 0.5*EE;
  double WDMAX;
  if(WM < WCMAX){ WDMAX = WM;}else{ WDMAX = WCMAX;}
      

  //  ****  Constants.

  const double RB =E_+TREV;
  const double GAM = 1.0+E_*RREV;
  const double GAM2 = GAM*GAM;
  const double BETA2 = (GAM2-1.0)/GAM2;
  const double AMOL = ((GAM-1.0)/GAM)*((GAM-1.0)/GAM);
  const double CPS = E_*RB;
  const double CP = sqrt(CPS);

  //  ************  Partial cross sections of the active oscillator.

  //  ****  Distant excitations.
      
  double CPPS =(E_-WKP)*(E_-WKP+TREV);
  double CPP = sqrt(CPPS);
  double QM;
  double RWKP, XHDL, XHDT, F0;
  if(WKP > 1.0E-6*E_)
  {
    QM = sqrt((CP-CPP)*(CP-CPP)+REV*REV)-REV;
  }
  else
  {
    QM = (WKP*WKP)/(BETA2*TREV);
    QM = QM*(1.0-QM*RTREV);
  }
  if(QM < QKP)
  {
    RWKP = 1.0/WKP;
    XHDL = log(QKP*(QM+TREV)/(QM*(QKP+TREV)))*RWKP;

    XHDT = log(GAM2)-BETA2-DELTA;
    if(XHDT < 0.0){ XHDT = 0.0;}
    
    XHDT = XHDT*RWKP;
    F0 = (WDMAX-WTHR)*(WM+WM-WDMAX-WTHR)/((WM-UK)*(WM-UK));
    XHDL = F0*XHDL;
    XHDT = F0*XHDT;
  }
  else
  {
    XHDL = 0.0;
    XHDT = 0.0;
  }
      //  ****  Close collisions.
  double RCL = WTHR/EE;
  double RL1 = 1.0-RCL;
  double RRL1 = 1.0/RL1;
  double XHC = (AMOL*(0.5-RCL)+1.0/RCL-RRL1+(1.0-AMOL)*log(RCL*RRL1))/EE;

  double XHTOT = XHC+XHDL+XHDT;
  if(XHTOT < 1.0E-35)
  {
    DE = UK;
    EP = E_-DE;
    CDT = 1.0;
    ES = 0.0;
    CDTS = 0.0;
    return;
  }

      //  ************  Sampling of final-state variables.

  TST=RAND(3.0)*XHTOT;

      //  ****  Hard close collision.

  double TS1 = XHC;
  double A, ARCL, FB, RK, RK2, RKF, PHI;
  if(TST < TS1)
  {
    A = 5.0*AMOL;
    ARCL = A*0.5*RCL;
    brkIt = false;
    while(!brkIt)
    {
      brkIt = true;
      FB = (1.0+ARCL)*RAND(4.0);
      if(FB < 1.0)
      {
        RK = RCL/(1.0-FB*(1.0-(RCL+RCL)));
      }
      else
      {
        RK = RCL+(FB-1.0)*(0.5-RCL)/ARCL;
      }
      RK2 = RK*RK;
      RKF = RK/(1.0-RK);
      PHI = 1.0+(RKF*RKF)-RKF+AMOL*(RK2+RKF);
      if(RAND(5.0)*(1.0+A*RK2) > PHI){ brkIt = false; continue;}
    }
    //  ****  Energy and scattering angle (primary electron).
    DE = RK*EE;
    EP = E_-DE;
    CDT = sqrt(EP*RB/(E_*(RB-DE)));
    //  ****  Energy and emission angle of the delta ray.
    ES = DE-UK;
    CDTS = sqrt(DE*RB/(E_*(DE+TREV)));
    return;
  }

  //  ****  Hard distant longitudinal interaction.

  double QS, Q, QTREV;
  TS1 = TS1+XHDL;
  DE = WM-sqrt((WM-WTHR)*(WM-WTHR)-RAND(7.0)*(WDMAX-WTHR)*(WM+WM-WDMAX-WTHR));
  EP = E_-DE;
  if(TST < TS1)
  {
    QS = QM/(1.0+QM*RTREV);
    Q = QS/(pow((QS/QKP)*(1.0+QKP*RTREV),RAND(6.0))-(QS*RTREV));
    QTREV = Q*(Q+TREV);
    CDT = (CPPS+CPS-QTREV)/(2.0*CP*CPP);
    if(CDT > 1.0){ CDT = 1.0;}
    //  ****  Energy and emission angle of the delta ray.
    ES = DE-UK;  // Inner shells only.
    CDTS = 0.5*(WKP*(E_+RB-WKP)+QTREV)/sqrt(CPS*QTREV);
    if(CDTS > 1.0){ CDTS = 1.0;}
    return;
  }

  //  ****  Hard distant transverse interaction.

  CDT = 1.0;
  //  ****  Energy and emission angle of the delta ray.
  ES = DE-UK;  // Inner shells only.
  CDTS = 1.0;
      
}

