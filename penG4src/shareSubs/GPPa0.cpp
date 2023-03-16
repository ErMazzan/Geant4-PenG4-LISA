
//  *********************************************************************
//                       SUBROUTINE GPPa0
//  *********************************************************************
void PenInterface::GPPa0(int M)
{
  //  Initialisation of the sampling algorithm for electron-positron pair
  //  production by photons in material M. Bethe-Heitler differential cross
  //  section.

  using namespace PENELOPE_mod;

  using namespace CADATA;
  using namespace COMPOS;
  using namespace CGPP00;
  
  //  ***  Effective atomic number.
  
  double FACT = 0.0;
  int IZZ;
  for(int I = 0; I < NELEM[M-1]; I++)
  {
    IZZ = IZ[M-1][I];
    FACT = FACT+IZZ*ATW[IZZ-1]*STF[M-1][I];
  }
  ZEQPP[M-1] = FACT/AT[M-1];
  IZZ = (int)(ZEQPP[M-1]+0.25);
  if(IZZ <= 0){ IZZ = 1;}
  if(IZZ > 99){ IZZ = 99;}
  //  ****  DBM Coulomb correction.
  double ALZ = ZEQPP[M-1]/SL;
  double A = ALZ*ALZ;
  double FC = A*(0.202059-A*(0.03693-A*(0.00835-A*(0.00201-A*(0.00049-A*(0.00012-A*0.00003)))))+1.0/(A+1.0));
  //  ****  Screening functions and low-energy correction.
  BCB[M-1] = 2.0/RSCR[IZZ-1];
  F0[M-1][0] = 4.0*log(RSCR[IZZ-1]);
  F0[M-1][1] = F0[M-1][0]-4.0*FC;
}

