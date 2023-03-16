
//  *********************************************************************
//                       SUBROUTINE GCOaT
//  *********************************************************************
void PenInterface::GCOaT(double &E, double &CS, int &M)
{
  //  Total cross section for incoherent (Compton) scattering. Relativistic
  //  Impulse approximation with analytical Compton profiles.

  //  Input arguments:
  //    E ........ photon energy (eV).
  //    M ........ material where photons propagate.
  //  Output argument:
  //    CS ....... incoherent total cross section (cm**2/molecule).

  using namespace PENELOPE_mod;

  using namespace CGCO;
  using namespace CGCO00;

////////  double GCOaD(double);

  const double RREV = 1.0/REV;
  const double PIELR2 = PI*ELRAD*ELRAD;

  double SXCO[NOCO];

  EE = E;
  MM = M;

  double EK, EKS, EK2, EK1, T0, CSL, TAU, CSKN, CSU;
  CS = 0.0;
  if(E < 5.0E6)
  {
    for(int IO = 0; IO < NOSCCO[M]; IO++)
    {
      IOSC = IO+1;
      SXCO[IO] = FCO[M][IO]*PIELR2*SUMGA(GCOaD,-1.0,1.0,1.0E-6);
      CS = CS+SXCO[IO];
    }
  }
  else
  {
    //  ****  Klein-Nishina total cross section.
    EK = E*RREV;
    EKS = EK*EK;
    EK2 = 1.0+EK+EK;
    EK1 = EKS-EK2-1.0;
    T0 = 1.0/(1.0+EK+EK);
    CSL = 0.5*EKS*T0*T0+EK2*T0+EK1*log(T0)-1.0/T0;
    for(int IO = 0; IO < NOSCCO[M]; IO++)
    {
      TAU = (E-UICO[M][IO])/E;
      if(TAU < T0)
      {
        CSKN = 0.0;
      }
      else
      {
        CSU = 0.5*EKS*TAU*TAU+EK2*TAU+EK1*log(TAU)-1.0/TAU;
        CSKN = PIELR2*(CSU-CSL)/(EK*EKS);
      }
      SXCO[IO] = FCO[M][IO]*CSKN;
      CS = CS+SXCO[IO];        
    }
  }
}

