
//  *********************************************************************
//                       SUBROUTINE GRAbT
//  *********************************************************************
void PenInterface::GRAbT(double &E, double &CS, int &M)
{
  //  Total cross section for Rayleigh (coherent) photon scattering. Form-
  //  factor approximation. Not used in the simulation!

  //  Input arguments:
  //    E ........ photon energy (eV).
  //    M ........ material where photons propagate.
  //  Output argument:
  //    CS ....... coherent total cross section (cm**2/molecule).

  using namespace CGRA00;
  
/////////////  double GRAaD(double);

  const double RREV = 1.0/REV;
  const double PIELR2 = PI*ELRAD*ELRAD;

  MM = M;
  MOM = 0;
  double SUM, EC;
  if(E < 5.0E7)
  {
    FACTE = 2.0*(E*RREV)*(E*RREV);
    SUM = SUMGA(GRAaD,-1.0,0.90,1.0E-6)+SUMGA(GRAaD,0.90,0.99,1.0E-6)+SUMGA(GRAaD,0.99,0.995,1.0E-6)+SUMGA(GRAaD,0.995,1.00,1.0E-6);
    CS = PIELR2*SUM;
  }
  else
  {
    EC = 5.0E7;
    FACTE = 2.0*(EC*RREV)*(EC*RREV);
    Q2MAX = 2.0*FACTE;
    SUM = SUMGA(GRAaD,-1.0,0.90,1.0E-6)+SUMGA(GRAaD,0.90,0.99,1.0E-6)+SUMGA(GRAaD,0.99,0.995,1.0E-6)+SUMGA(GRAaD,0.995,1.00,1.0E-6);
    CS = PIELR2*SUM;
    CS = (EC/E)*(EC/E)*CS;
    FACTE = 2.0*(E*RREV)*(E*RREV);
  }     
}

