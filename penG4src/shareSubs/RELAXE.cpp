
//  *********************************************************************
//                       SUBROUTINE RELAXE
//  *********************************************************************
void PenInterface::RELAXE(int &IZ, int &JS0, int &JS1, int &JS2, double &ENERGY, double &TRPROB)
{
  //  This subroutine gives the energy (ENERGY, in eV) and the probability
  //  (PROB) of the transition JS0-JS1-JS2 of an atom of the element IZ
  //  with an initial vacancy in shell JS0.

  //  The output values ENERGY=1.0D35 and TRPROB=1.0D35 indicate that the
  //  transition data are not loaded, or that the transition does not
  //  exist.

  using namespace CRELAX;

  ENERGY = 1.0E35;
  TRPROB = 1.0E35;
  if(IZ < 3 || IZ > 99 || JS0 < 1 || JS0 > 16){ return;}
  
  int N1 = IFIRST[IZ-1][JS0-1];
  if(N1 == 0){ return;}
  int NL = ILAST[IZ-1][JS0-1];
  if(NL < N1){ return;}

  for(int N = N1-1; N < NL; N++)
  {
    if(JS1 == IS1[N] && JS2 == IS2[N])
    {
      ENERGY = ET[N];
      TRPROB = P[N];
      return;
    }
  } 
}

