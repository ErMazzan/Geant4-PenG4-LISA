
//  *********************************************************************
//                       SUBROUTINE EINaT1
//  *********************************************************************
void PenInterface::EINaT1(double &E, double &UK, double& WK, double DELTA, double &WCCM, double &H0, double &H1, double &H2, double &S0, double &S1, double &S2, double &R0, double &R1, double &R2)
{
  //  Integrated cross sections for inelastic collisions of electrons with
  //  a single-shell oscillator, restricted to energy losses larger than,
  //  and smaller than, the cutoff energy loss WCCM.
  //
  //  Sternheimer-Liljequist oscillator model.
  //
  //  Input arguments:
  //    E ..... kinetic energy (eV).
  //    UK .... ionisation energy (eV).
  //    WK .... resonance energy (eV).
  //    DELTA ... Fermi's density effect correction.
  //    WCCM ... cutoff energy loss (eV).
  //
  //  Output arguments:
  //    H0 .... total cross section for hard colls. (cm**2).
  //    H1 .... stopping cross section for hard colls. (eV*cm**2).
  //    H2 .... straggling cross section for hard colls. (eV**2*cm**2).
  //    S0 .... total cross section for soft colls. (cm**2).
  //    S1 .... stopping cross section for soft colls. (eV*cm**2).
  //    S2 .... straggling cross section for soft colls. (eV**2*cm**2).
  //    R0 .... total cross section for soft colls. (cm**2).
  //    R1 .... 1st transport cross section for soft colls. (cm**2).
  //    R2 .... 2nd transport cross section for soft colls. (cm**2).

  using namespace CEIN01;
  
  const double TREV = 2.0*REV;
  const double RTREV = 1.0/TREV;
  const double PIELR2 = PI*ELRAD*ELRAD;

  H0 = 0.0;
  H1 = 0.0;
  H2 = 0.0;
  S0 = 0.0;
  S1 = 0.0;
  S2 = 0.0;
  R0 = 0.0;
  R1 = 0.0;
  R2 = 0.0;

  double WTHR;
  if(UK > 1.0E-3)
  {
    WTHR = UK;
  }
  else
  {
    WTHR = WK;
  }
  if(E < WTHR+1.0E-6){return;}

  //  ****  Constants.

  EI = E;
  const double GAM = 1.0+E/REV;
  const double GAM2 = GAM*GAM;
  const double BETA2 = (GAM2-1.0)/GAM2;
  const double CONST = PIELR2*TREV/BETA2;

  CPS = E*(E+TREV);
  double CP = sqrt(CPS);
  AMOL = (E/(E+REV))*(E/(E+REV));

  //  ****  Trick: The resonance energy and the cutoff recoil energy of
  //        inner shells are varied to yield a smooth threshold.

  double WM, WKP, QKP, WCMAX, WDMAX;
  if(UK > 1.0E-3)
  {
    WM = 3.0*WK-2.0*UK;
    if(E > WM)
    {
      WKP = WK;
      QKP = UK;
    }
    else
    {
      WKP = (E+2.0*UK)/3.0;
      QKP = UK*(E/WM);
      WM = E;
    }
    EE = E+UK;
    WCMAX = 0.5*EE;
    if(WCMAX < WM){ WDMAX = WCMAX;}
    else{ WDMAX = WM;}
  }

  else
  {
    WM = E;
    WKP = WK;
    QKP = WK;
    EE = E;
    WCMAX = 0.5*EE;
    WDMAX = WKP+1.0;
  }

  //  ****  Distant interactions.

  double SDL1 = 0.0;
  double SDT1 = 0.0;
  double CPPS, CPP;
  double A, B, BA;
  double QM;
  double RMU1;
  if(WDMAX > WTHR+1.0E-6)
  {
    CPPS = (E-WKP)*(E-WKP+TREV);
    CPP = sqrt(CPPS);
    A = 4.0*CP*CPP;
    B = ((CP-CPP)*(CP-CPP));

    if(WKP > 1.0E-6*E)
    {
      QM = sqrt(((CP-CPP)*(CP-CPP))+(REV*REV))-REV;
    }
    else
    {
      QM = WKP*WKP/(BETA2*TREV);
      QM = QM*(1.0-QM*RTREV);
    }
    if(QM < QKP)
    {
      SDL1 = log(QKP*(QM+TREV)/(QM*(QKP+TREV)));
      SDT1 = log(GAM2)-BETA2-DELTA;
      if(SDT1 < 0.0){ SDT1 = 0.0;}
    
    //  ****  Soft distant transport moments of orders 0-2.
      if(WCCM > WTHR)
      {
        BA = B/A;
        RMU1 = (QKP*(QKP+TREV)-B)/A;
        R0 = log((RMU1+BA)/BA);
        R1 = RMU1-BA*R0;
        R2 = (BA*BA)*R0+0.5*RMU1*(RMU1-2.0*BA);
        R0 = R0/WKP;
        R1 = R1/WKP;
        R2 = R2/WKP;
        R0 = R0+SDT1/WKP;
      }
    }
  }

  double F0, F1, WL, WU;
  double WU2, WL2;  //MACG to speed up pow(Wx,n)
  double SD1 = SDL1+SDT1;
  if(SD1 > 0.0)
  {
    if(UK > 1.0E-3)
    {
    //  ****  Inner-shell excitations (triangle distribution).
      F0 = 1.0/((WM-UK)*(WM-UK));
      F1 = 2.0*F0*SD1/WKP;
      if(WCCM < UK)
      {
        WL = UK;
        WU = WDMAX;
	WL2 = WL*WL;
	WU2 = WU*WU;
        H0 = F1*(WM*(WU-WL)-(WU2-WL2)/2.0);
        H1 = F1*(WM*(WU2-WL2)/2.0-((WU2*WU)-(WL2*WL))/3.0);
        H2 = F1*(WM*((WU2*WU)-(WL2*WL))/3.0-((WU2*WU2)-(WL2*WL2))/4.0);
      }
      else
      {
        if(WCCM > WDMAX)
        {
          WL = UK;
          WU = WDMAX;
	  WL2 = WL*WL;
	  WU2 = WU*WU;
          S0 = F1*(WM*(WU-WL)-(WU2-WL2)/2.0);
          S1 = F1*(WM*(WU2-WL2)/2.0-((WU2*WU)-(WL2*WL))/3.0);
          S2 = F1*(WM*((WU2*WU)-(WL2*WL))/3.0-((WU2*WU2)-(WL2*WL2))/4.0);
        }
        else
        {
          WL = WCCM;
          WU = WDMAX;
	  WL2 = WL*WL;
	  WU2 = WU*WU;
          H0 = F1*(WM*(WU-WL)-(WU2-WL2)/2.0);
          H1 = F1*(WM*(WU2-WL2)/2.0-((WU2*WU)-(WL2*WL))/3.0);
          H2 = F1*(WM*((WU2*WU)-(WL2*WL))/3.0-((WU2*WU2)-(WL2*WL2))/4.0);
          WL = UK;
          WU = WCCM;
	  WL2 = WL*WL;
	  WU2 = WU*WU;
          S0 = F1*(WM*(WU-WL)-(WU2-WL2)/2.0);
          S1 = F1*(WM*(WU2-WL2)/2.0-((WU2*WU)-(WL2*WL))/3.0);
          S2 = F1*(WM*((WU2*WU)-(WL2*WL))/3.0-((WU2*WU2)-(WL2*WL2))/4.0);
        }
        double F2 = F0*(2.0*WM*(WU-WL)-(WU2-WL2));
        R0 = F2*R0;
        R1 = F2*R1;
        R2 = F2*R2;
      }
    }
    else
    {
    //  ****  Outer-shell excitations (delta oscillator).
      if(WCCM < WKP)
      {
        H1 = SD1;
        H0 = SD1/WKP;
        H2 = SD1*WKP;
      }
      else
      {
        S1 = SD1;
        S0 = SD1/WKP;
        S2 = SD1*WKP;
      }
    }
  }

  //  ****  Close collisions (Moller's cross section).

  if(WCMAX < WTHR+1.0E-6){}
  else
  {
    double EE2;   //MACG to speed up pow(EE,2)
    if(WCCM < WTHR)   // No soft interactions.
    {
      WL = WTHR;
      WU = WCMAX;
      WL2 = WL*WL;
      WU2 = WU*WU;
      EE2 = EE*EE;
      H0 = H0+(1.0/(EE-WU))-(1.0/(EE-WL))-(1.0/WU)+(1.0/WL)+(1.0-AMOL)*log(((EE-WU)*WL)/((EE-WL)*WU))/EE+AMOL*(WU-WL)/EE2;
    
      H1 = H1+log(WU/WL)+(EE/(EE-WU))-(EE/(EE-WL))+(2.0-AMOL)*log((EE-WU)/(EE-WL))+AMOL*(WU2-WL2)/(2.0*EE2);
    
      H2 = H2+(2.0-AMOL)*(WU-WL)+(WU*(2.0*EE-WU)/(EE-WU))-(WL*(2.0*EE-WL)/(EE-WL))+(3.0-AMOL)*EE*log((EE-WU)/(EE-WL))+AMOL*((WU2*WU)-(WL2*WL))/(3.0*EE2);
    }
    else
    {
      if(WCCM > WCMAX)
      {
        WL = WTHR;
        WU = WCMAX;
	WL2 = WL*WL;
	WU2 = WU*WU;
	EE2 = EE*EE;
        S0 = S0+(1.0/(EE-WU))-(1.0/(EE-WL))-(1.0/WU)+(1.0/WL)+(1.0-AMOL)*log(((EE-WU)*WL)/((EE-WL)*WU))/EE+AMOL*(WU-WL)/EE2;
        
        S1 = S1+log(WU/WL)+(EE/(EE-WU))-(EE/(EE-WL))+(2.0-AMOL)*log((EE-WU)/(EE-WL))+AMOL*(WU2-WL2)/(2.0*EE2);
        
        S2 = S2+(2.0-AMOL)*(WU-WL)+(WU*(2.0*EE-WU)/(EE-WU))-(WL*(2.0*EE-WL)/(EE-WL))+(3.0-AMOL)*EE*log((EE-WU)/(EE-WL))+AMOL*((WU2*WU)-(WL2*WL))/(3.0*EE2);
      }
      else
      {
        WL = WCCM;
        WU = WCMAX;
	WL2 = WL*WL;
	WU2 = WU*WU;
	EE2 = EE*EE;
        H0 = H0+(1.0/(EE-WU))-(1.0/(EE-WL))-(1.0/WU)+(1.0/WL)+(1.0-AMOL)*log(((EE-WU)*WL)/((EE-WL)*WU))/EE+AMOL*(WU-WL)/EE2;
        
        H1 = H1+log(WU/WL)+(EE/(EE-WU))-(EE/(EE-WL))+(2.0-AMOL)*log((EE-WU)/(EE-WL))+AMOL*(WU2-WL2)/(2.0*EE2);
        
        H2 = H2+(2.0-AMOL)*(WU-WL)+(WU*(2.0*EE-WU)/(EE-WU))-(WL*(2.0*EE-WL)/(EE-WL))+(3.0-AMOL)*EE*log((EE-WU)/(EE-WL))+AMOL*((WU2*WU)-(WL2*WL))/(3.0*EE2);
        WL = WTHR;
        WU = WCCM;
	WL2 = WL*WL;
	WU2 = WU*WU;
        S0 = S0+(1.0/(EE-WU))-(1.0/(EE-WL))-(1.0/WU)+(1.0/WL)+(1.0-AMOL)*log(((EE-WU)*WL)/((EE-WL)*WU))/EE+AMOL*(WU-WL)/EE2;
        
        S1 = S1+log(WU/WL)+(EE/(EE-WU))-(EE/(EE-WL))+(2.0-AMOL)*log((EE-WU)/(EE-WL))+AMOL*(WU2-WL2)/(2.0*EE2);
        
        S2 = S2+(2.0-AMOL)*(WU-WL)+(WU*(2.0*EE-WU)/(EE-WU))-(WL*(2.0*EE-WL)/(EE-WL))+(3.0-AMOL)*EE*log((EE-WU)/(EE-WL))+AMOL*((WU2*WU)-(WL2*WL))/(3.0*EE2);
      }
      //  ****  Soft close transport moments of orders 0-2.
      double CP2S = (E-WL)*(E-WL+TREV);
      double CP2 = sqrt(CP2S);
      double RMU2 = (WL*(WL+TREV)-((CP-CP2)*(CP-CP2)))/(4.0*CP*CP2);
      double CP3S = (E-WU)*(E-WU+TREV);
      double CP3 = sqrt(CP3S);
      double RMU3 = (WU*(WU+TREV)-((CP-CP3)*(CP-CP3)))/(4.0*CP*CP3);
      MOM = 0;
      R0 = R0+SUMGA(EINaDS,RMU2,RMU3,1.0E-7);
      MOM = 1;
      R1 = R1+SUMGA(EINaDS,RMU2,RMU3,1.0E-7);
      MOM = 2;
      R2 = R2+SUMGA(EINaDS,RMU2,RMU3,1.0E-7);
    }
  }
      
  H0 = CONST*H0;
  H1 = CONST*H1;
  H2 = CONST*H2;
  S0 = CONST*S0;
  S1 = CONST*S1;
  S2 = CONST*S2;
  R0 = CONST*R0;
  R1 = CONST*R1;
  R2 = CONST*R2;      
}

