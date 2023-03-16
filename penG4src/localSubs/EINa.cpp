
//  *********************************************************************
//                       SUBROUTINE EINa
//  *********************************************************************
void PenPhys::EINa(double E_, double DELTA, double &DE, double &EP, double &CDT, double &ES, double &CDTS, int M, int &IOSC)
{
  //  Random sampling of hard inelastic collisions of electrons.
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
  //    IOSC .... index of the oscillator that has been 'ionized'.

  using namespace PENELOPE_mod;

  using namespace CECUTR;
////////////////  using namespace CEGRID;
  using namespace CEIN;
  using namespace CEINAC;
  
  bool LDIST;
  const double RREV = 1.0/REV;
  const double TREV = 2.0*REV;
  const double RTREV = 1.0/TREV;
  
  double WCCM = WCC[M-1];
  if(WCCM > E_)
  {
    DE = 0.0;
    EP = E_;
    CDT = 1.0;
    ES = 0.0;
    CDTS = 0.0;
    IOSC = NO;
    return;
  }
  //  ****  Energy grid point.
  double PK = (CEGRID_.XEL-CEGRID_.DLEMP[CEGRID_.KE-1])*CEGRID_.DLFC;
  int JE;
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
  //  ****  Binary search.
  int IT;
  int IO = 1;
  int JO = NEIN[M-1]+1;
  bool brkIt = false;
  while(!brkIt)
  {
    brkIt = true;
    IT = (IO+JO)/2;
    if(TST > EINAC[M-1][JE-1][IT-1])
    {
      IO = IT;
    }
    else
    {
      JO = IT;
    }
    if(JO-IO > 1){ brkIt = false; continue;}
  }
  IOSC = IEIN[M-1][IO-1];
  double UK = UI[M-1][IOSC-1];
  double WK=WRI[M-1][IOSC-1];
  double WTHR;
  if(UK > 1.0E-3)
  {
    WTHR = WCCM;
    if(WTHR < UK){ WTHR = UK;}
  }
  else
  {
    WTHR = WCCM;
    if(WTHR < WK){ WTHR = WK;}
  }

  if(E_ < WTHR+1.0E-6)
  {
    DE = 0.0;
    EP = E_;
    CDT = 1.0;
    ES = 0.0;
    CDTS = 0.0;
    IOSC = NO;
    return;
  }

  //  ****  Trick: The resonance energy and the cutoff recoil energy of
  //        inner shells are varied to yield a smooth threshold.

  LDIST = true;
  double WM, WKP, QKP;
  double EE, WCMAX, WDMAX;
  if(UK > 1.0E-3)
  {
    WM = 3.0*WK-2.0*UK;
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
    if(WCCM > WM){ LDIST = false;}
    EE = E_+UK;
    WCMAX = 0.5*EE;
    WDMAX = WM;
    if(WDMAX > WCMAX){ WDMAX = WCMAX;}
    if(WTHR > WDMAX){ LDIST = false;}
  }
  else
  {
    if(WCCM > WK){ LDIST = false;}
    WKP = WK;
    QKP = WK;
    WM = E_;
    EE = E_;
    WCMAX = 0.5*EE;
    WDMAX = WKP+1.0;
  }

  //  ****  Constants.

  double RB = E_+TREV;
  double GAM = 1.0+E_*RREV;
  double GAM2 = GAM*GAM;
  double BETA2 = (GAM2-1.0)/GAM2;
  double AMOL = ((GAM-1.0E0)/GAM)*((GAM-1.0E0)/GAM);
  double CPS = E_*RB;
  double CP = sqrt(CPS);

  //  ************  Partial cross sections of the active oscillator.

  //  ****  Distant excitations.
  double CPP, CPPS;
  double QM, XHDL, XHDT;
  if(LDIST)
  {
      
    CPPS = (E_-WKP)*(E_-WKP+TREV);
    CPP = sqrt(CPPS);
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
      double RWKP = 1.0/WKP;
      XHDL = log(QKP*(QM+TREV)/(QM*(QKP+TREV)))*RWKP;
      XHDT = log(GAM2)-BETA2-DELTA;
      if(XHDT < 0.0){ XHDT = 0.0;}
      XHDT = XHDT*RWKP;
      if(UK > 1.0E-3)
      {
      double F0 = (WDMAX-WTHR)*(WM+WM-WDMAX-WTHR)/((WM-UK)*(WM-UK));
      XHDL = F0*XHDL;
      XHDT = F0*XHDT;
      }
    }
    else
    {
      XHDL = 0.0;
      XHDT = 0.0;
    }
  }
  else
  {
    QM = 0.0;    // Defined to avoid compilation warnings.
    CPP = 0.0;   // same
    CPPS = 0.0;  // same
    XHDL = 0.0;
    XHDT = 0.0;
  }
  //  ****  Close collisions.
  double RCL = WTHR/EE;
  double RL1, RRL1, XHC;
  if(RCL < 0.5)
  {
    RL1 = 1.0-RCL;
    RRL1 = 1.0/RL1;
    XHC = (AMOL*(0.5-RCL)+1.0/RCL-RRL1+(1.0-AMOL)*log(RCL*RRL1))/EE;
  }
  else
  {
    XHC = 0.0;
  }

  double XHTOT = XHC+XHDL+XHDT;
  if(XHTOT < 1.0E-35)
  {
    DE = 0.0;
    EP = E_;
    CDT = 1.0;
    ES = 0.0;
    CDTS = 0.0;
    IOSC = NO;
    return;
  }

  //  ************  Sampling of final-state variables.

  TST = RAND(3.0)*XHTOT;
  
  //  ****  Hard close collision.

  double TS1 = XHC;
  double A, ARCL, RK, RK2, RKF, PHI;
  if(TST < TS1)
  {
    A = 5.0*AMOL;
    ARCL = A*0.5*RCL;
    brkIt = false;
    while(!brkIt)
    {
      brkIt = true;
      double FB = (1.0+ARCL)*RAND(4.0);
      if(FB < 1.0)
      {
        RK = RCL/(1.0E0-FB*(1.0E0-(RCL+RCL)));
      }
      else
      {
        RK = RCL+(FB-1.0E0)*(0.5E0-RCL)/ARCL;
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
    if(KS[M-1][IOSC-1] < 17)
    {
      if(UK > ECUTR[M-1])
      {
        ES = DE-UK;  // Inner shells only.
      }
      else
      {
        ES = DE;
      }
    }
    else
    {
      ES = DE;
    }
    CDTS = sqrt(DE*RB/(E_*(DE+TREV)));
    return;
  }

  //  ****  Hard distant longitudinal interaction.

  TS1 = TS1+XHDL;
  if(UK > 1.0E-3)
  {
    DE = WM-sqrt((WM-WTHR)*(WM-WTHR)-RAND(7.0)*(WDMAX-WTHR)*(WM+WM-WDMAX-WTHR));
  }
  else
  {
    DE = WKP;
  }
  EP = E_-DE;
  double QS, Q, QTREV;
  if(TST < TS1)
  {
    QS = QM/(1.0+QM*RTREV);
    Q = QS/(pow((QS/QKP)*(1.0+QKP*RTREV),RAND(6.0))-(QS*RTREV));
    QTREV = Q*(Q+TREV);
    CDT = (CPPS+CPS-QTREV)/(2.0*CP*CPP);
    if(CDT > 1.0){ CDT = 1.0;}
    //  ****  Energy and emission angle of the delta ray.
    if(KS[M-1][IOSC-1] < 17)
    {
      ES = DE-UK;  // Inner shells only.
    }
    else
    {
      ES = DE;
    }
    CDTS = 0.5*(WKP*(E_+RB-WKP)+QTREV)/sqrt(CPS*QTREV);
    if(CDTS > 1.0){ CDTS = 1.0;}
    return;
  }

  //  ****  Hard distant transverse interaction.

  CDT = 1.0;
  //  ****  Energy and emission angle of the delta ray.
  if(KS[M-1][IOSC-1] < 17)
  {
    if(UK > ECUTR[M-1])
    {
      ES = DE-UK;  // Inner shells only.
    }
    else
    {
      ES = DE;
    }
  }
  else
  {
    ES = DE;
  }
  CDTS = 1.0;

}

