
//  *********************************************************************
//                       SUBROUTINE PANaR 
//  *********************************************************************
void PenPhys::PANaR(double &ECUT)
{
  //  Simulation of annihilation of positrons at rest. Annihilation quanta
  //  are stored in the secondary stack only when ECUT is less than REV.

  using namespace PENELOPE_mod;
  using namespace PENERROR_mod;
  //////////  using namespace TRACK_mod;

  const double TWOPI = PI+PI;
  int ILBA[5];

  if(REV < ECUT){ return;}

  double US = U;
  double VS = V;
  double WS = W;
  
  double CDT1 = -1.0+2.0*RAND(1.0);
  double DF = TWOPI*RAND(2.0);
  DIRECT(CDT1,DF,US,VS,WS);
  ILBA[0] = ILB[0]+1;
  ILBA[1] = 3;
  ILBA[2] = 6;
  ILBA[3] = 0;
  ILBA[4] = ILB[4];
  STORES(REV,X,Y,Z,US,VS,WS,WGHT,2,ILBA,0);
  if(IRETRN != 0){ return;}
  STORES(REV,X,Y,Z,-US,-VS,-WS,WGHT,2,ILBA,0);
  if(IRETRN != 0){ return;}

}

