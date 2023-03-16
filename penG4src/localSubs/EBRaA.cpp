
//  *********************************************************************
//                       SUBROUTINE EBRaA
//  *********************************************************************
void PenPhys::EBRaA(double &E_, double &DE, double &CDT, int &M)
{
  //  Random sampling of the initial direction of bremss photons, relative
  //  to the direction of the projectile.
  //  Numerical fit/interpolation of partial-wave shape functions generated
  //  by the program BREMS of A. Poskus, Comp. Phys. Commun. (2018).

  //  Input parameters:
  //    M ..... material where the projectile moves.
  //    E_ ..... kinetic energy of the projectile.
  //    DE .... energy of the emitted photon.
  //  Output parameter:
  //    CDT ... cosine of the polar emission angle.

  using namespace PENELOPE_mod;
  using namespace CBRANG;
  
  const double TREV = 2.0*REV;
  
  //  ****  Distribution parameters.
  
  double BETA = sqrt(E_*(E_+TREV))/(E_+REV);

  //  A pure dipole distribution is used for E>1 MeV.
  if(E_ > 1.0E6)
  {
    CDT = 2.0*RAND(1.0)-1.0;
    if(RAND(2.0) > 0.75)
    {
      if(CDT < 0.0)
      {
        CDT = -pow(-CDT,0.333333333333333);
      }
      else
      {
      CDT = pow(CDT,0.333333333333333);
      }
    }
    CDT = (CDT+BETA)/(1.0+BETA*CDT);
    return;
  }

  int IE, IET, IE1;
  if(BETA > BET[6])
  {
    BETA = BET[6];
    IE = 6;
  }
  else
  {
    if(BETA < BET[0])
    {
      IE = 1;
    }
    else
    {
      IE = 1;
      IE1 = 7;
      bool brkIt = false;
      while(!brkIt)
      {
        brkIt = true;
        IET = (IE+IE1)/2;
        if(BETA > BET[IET-1])
        {
          IE = IET;
        }
        else
        {
          IE1 = IET;
        }
        if(IE1-IE > 1){ brkIt = false; continue;}
      }
    }
  }

  double RK = 1.0+40.0*DE/E_;
  int IK = int(RK);
  if(IK > 40){ IK = 40;}

  double P10 = BP1[M-1][IE-1][IK-1][0]+BETA*(BP1[M-1][IE-1][IK-1][1]+BETA*(BP1[M-1][IE-1][IK-1][2]+BETA*BP1[M-1][IE-1][IK-1][3]));
  double P11 = BP1[M-1][IE-1][IK+1-1][0]+BETA*(BP1[M-1][IE-1][IK+1-1][1]+BETA*(BP1[M-1][IE-1][IK+1-1][2]+BETA*BP1[M-1][IE-1][IK+1-1][3]));
  double P1 = std::exp(P10+(RK-IK)*(P11-P10));

  double P20 = BP2[M-1][IE-1][IK-1][0]+BETA*(BP2[M-1][IE-1][IK-1][1]+BETA*(BP2[M-1][IE-1][IK-1][2]+BETA*BP2[M-1][IE-1][IK-1][3]));
  double P21 = BP2[M-1][IE-1][IK+1-1][0]+BETA*(BP2[M-1][IE-1][IK+1-1][1]+BETA*(BP2[M-1][IE-1][IK+1-1][2]+BETA*BP2[M-1][IE-1][IK+1-1][3]));
  double P2 = P20+(RK-IK)*(P21-P20);

  //  ****  Sampling from the Lorentz-transformed dipole distributions.

  // P1 = exp(P1)/BETA;
  // if(P1 > 1.0){ P1 = 1.0;}

  // double BETAP = BETA*(1.0+P2/BETA);
  // if(BETAP < 0.0){ BETAP = 0.0;}
  // if(BETAP > 0.999999999){ BETAP = 0.999999999;}
  
  if(RAND(3.0) < P1)
  {
    bool brkIt = false;
    while(!brkIt)
    {
      brkIt = true;
      CDT = 2.0*RAND(4.0)-1.0;
      if(2.0*RAND(5.0) > 1.0+CDT*CDT){ brkIt = false; continue;}
    }
  }
  else
  {
    bool brkIt = false;
    while(!brkIt)
    {
      brkIt = true;
      CDT = 2.0*RAND(4.0)-1.0;
      if(RAND(5.0) > 1.0-CDT*CDT){ brkIt = false; continue;}
    }
  }
  CDT = (CDT+P2)/(1.0+P2*CDT);

}

