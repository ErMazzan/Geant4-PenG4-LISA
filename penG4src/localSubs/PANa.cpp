
//  *********************************************************************
//                       SUBROUTINE PANa
//  *********************************************************************
void PenPhys::PANa(double &E_, double &E1, double &CDT1, double &E2, double &CDT2, int &M)
{
  //  Simulation of positron annihilation (either at rest or in flight) in
  //  material M. Ei and CDTi are the energies and polar direction cosines
  //  of the two annihilation photons.

  using namespace PENELOPE_mod;

  const double TREV = 2.0*REV;

  //  ****  Slow positrons (assumed at rest).

  if(E_ < EABS[2][M-1])
  {
    E1 = 0.5*(E_+TREV);
    E2 = E1;
    CDT1 = -1.0+2.0*RAND(1.0);
    CDT2 = -CDT1;
  }
  else
  {
    //  ****  Annihilation in flight (two photons with energy and directions
    //        determined from the dcs and energy-momentum conservation).
    double GAM;
    if(E_ < 1.0){ GAM = 1.0+1.0/REV;}
    else{ GAM = 1.0+E_/REV;}
      
    double GAM21 = sqrt(GAM*GAM-1.0);
    double ANI = 1.0+GAM;
    double CHIMIN = 1.0/(ANI+GAM21);
    double RCHI = (1.0-CHIMIN)/CHIMIN;
    double GT0 = ANI*ANI-2.0;
    bool brkIt = false;
    double CHI, GREJ;
    while(!brkIt)
    {
      brkIt = true;
      CHI = CHIMIN*pow(RCHI,RAND(2.0));
      GREJ = ANI*ANI*(1.0-CHI)+GAM+GAM-1.0/CHI;
      if(RAND(3.0)*GT0 > GREJ){ brkIt = false; continue;}
    }
      
    double DET = E_+TREV;
    E1 = CHI*DET;
    CDT1 = (ANI-1.0/CHI)/GAM21;
    double CHIP = 1.0-CHI;
    E2 = DET-E1;
    CDT2 = (ANI-1.0/CHIP)/GAM21;
  }

}

