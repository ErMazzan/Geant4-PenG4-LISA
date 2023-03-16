
//  *********************************************************************
//                       SUBROUTINE GRAaTI
//  *********************************************************************
void PenInterface::GRAaTI(double &E, double &ECS, int M)
{
  //  Total cross section for Rayleigh (coherent) photon scattering, in
  //  cm**2. Interpolated from input data.

  using namespace PENELOPE_mod;

////////////////  using namespace CEGRID;
  using namespace COMPOS;
  using namespace CGRA01;
  
  double XELN = log(E);
  double XEN = 1.0+(XELN-CEGRID_.DLEMP1)*CEGRID_.DLFC;
  int    KEN = (int)XEN;
  if(KEN == 0){ KEN = 1;}

  //  ****  Binary search.

  int II = IED[KEN-1];
  int IU = IEU[KEN-1];
  bool brkIt = false;
  while(!brkIt)
  {
    brkIt = true;
    int IT=(II+IU)/2;
    if(XELN > ERA[IT-1])
    {
      II = IT;
    }
    else
    {
      IU = IT;
    }
    if(IU-II > 1){ brkIt = false; continue;}
  }
  ECS = exp(XSRA[M-1][II-1]+(XSRA[M-1][II+1-1]-XSRA[M-1][II-1])*(XELN-ERA[II-1])/(ERA[II+1-1]-ERA[II-1]))/VMOL[M-1];
}

