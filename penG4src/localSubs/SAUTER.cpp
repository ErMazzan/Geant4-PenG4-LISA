
//  *********************************************************************
//                       SUBROUTINE SAUTER
//  *********************************************************************
void PenPhys::SAUTER(double &ES, double &CDTS)
{
  //  Random sampling of the initial direction of photoelectrons from the
  //  Sauter distribution.

  if(ES > 1.0E9)
  {
    CDTS = 1.0;
    return;
  }
  double GAM = 1.0+ES/REV;
  double GAM2 = GAM*GAM;
  double BETA = sqrt((GAM2-1.0)/GAM2);
  double AC = 1.0/BETA-1.0;
  double A1 = 0.5*BETA*GAM*(GAM-1.0)*(GAM-2.0);
  double A2 = AC+2.0;
  double GTMAX = 2.0*(A1+1.0/AC);
  double TSAM;
  bool brkIt = false;
  while(!brkIt)
  {
    brkIt = true;
    double RU = RAND(1.0);
    TSAM = 2.0*AC*(2.0*RU+A2*sqrt(RU))/(A2*A2-4.0*RU);
    double GTR = (2.0-TSAM)*(A1+1.0/(AC+TSAM));
    if(RAND(2.0)*GTMAX > GTR){ brkIt = false; continue;}
  }
  CDTS = 1.0-TSAM;
}

