
//  *********************************************************************
//                       SUBROUTINE PINaT
//  *********************************************************************
void PenInterface::PINaT(double &E, double &WCCM, double &XH0, double &XH1, double &XH2, double &XS0, double &XS1, double &XS2, double &XT1, double &XT2, double &DELTA, int M)
{
  //  Integrated cross sections for inelastic collisions of positrons of
  //  energy E in material M, restricted to energy losses larger than and
  //  less than the cutoff energy WCCM.
  //
  //  Sternheimer-Liljequist GOS model.
  //
  //  Output arguments:
  //    XH0 ... total cross section for hard colls. (cm**2).
  //    XH1 ... stopping cross section for hard colls. (eV*cm**2).
  //    XH2 ... straggling cross section for hard colls. (eV**2*cm**2).
  //    XS0 ... total cross section for soft colls. (cm**2).
  //    XS1 ... stopping cross section for soft colls. (eV*cm**2)
  //    XS2 ... straggling cross section for soft colls. (eV**2*cm**2).
  //    XT1 ... 1st transport cross section for soft colls. (cm**2).
  //    XT2 ... 2nd transport cross section for soft colls. (cm**2).
  //    DELTA ... Fermi's density effect correction.

  using namespace PENELOPE_mod;

  using namespace COMPOS;
  using namespace CEIN;
  using namespace CPIN00;

  //  ****  Constants.

  const double GAM = 1.0+E/REV;
  const double GAM2 = GAM*GAM;

  //  ************  Density effect.

  //  ****  Sternheimer's resonance energy (WL2=L**2).
  bool brkIt = false;
  double TST, WL2, FDEL, WL2L, WL2U;
  TST = ZT[M-1]/(GAM2*OP2[M-1]);
  WL2 = 0.0;
  FDEL = 0.0;
  for(int I = 0; I < NOSC[M-1]; I++)
  {
    FDEL = FDEL+F[M-1][I]/((WRI[M-1][I])*(WRI[M-1][I])+WL2);
  }
  if(FDEL < TST)
  {
    DELTA = 0.0;
  }
  else
  {
    WL2 = WRI[M-1][NOSC[M-1]-1]*WRI[M-1][NOSC[M-1]-1];
    brkIt = false;
    while(!brkIt)
    {
      brkIt = true;
      WL2 = WL2+WL2;
      FDEL = 0.0;
      for(int I = 0; I < NOSC[M-1]; I++)
      {
        FDEL = FDEL+F[M-1][I]/((WRI[M-1][I])*(WRI[M-1][I])+WL2);
      }
      if(FDEL > TST){ brkIt = false; continue;}
    }
    WL2L = 0.0;
    WL2U = WL2;
    brkIt = false;
    while(!brkIt)
    {
      brkIt = true;
      WL2 = 0.5*(WL2L+WL2U);
      FDEL = 0.0;
      for(int I = 0; I < NOSC[M-1]; I++)
      {
        FDEL = FDEL+F[M-1][I]/((WRI[M-1][I])*(WRI[M-1][I])+WL2);
      }
      if(FDEL > TST)
      {
        WL2L = WL2;
      }
      else
      {
        WL2U = WL2;
      }
      if(WL2U-WL2L > 1.0E-12*WL2){ brkIt = false; continue;}
    }
      //  ****  Density effect correction (delta).
    DELTA = 0.0;
    for(int I = 0; I < NOSC[M-1]; I++)
    {
      DELTA = DELTA+F[M-1][I]*log(1.0+WL2/((WRI[M-1][I])*(WRI[M-1][I])));
    }
    DELTA = (DELTA/ZT[M-1])-WL2/(GAM2*OP2[M-1]);
  }

  //  ****  Shell-oscillator cross sections.

  for(int I = 0; I < NOSC[M-1]; I++)
  {
    SPH0[I] = 0.0;
    SPH1[I] = 0.0;
    SPH2[I] = 0.0;
    SPS0[I] = 0.0;
    SPS1[I] = 0.0;
    SPS2[I] = 0.0;
    SPT0[I] = 0.0;
    SPT1[I] = 0.0;
    SPT2[I] = 0.0;
  }
  XH0 = 0.0;
  XH1 = 0.0;
  XH2 = 0.0;
  XS0 = 0.0;
  XS1 = 0.0;
  XS2 = 0.0;
  double XT0 = 0.0;
  XT1 = 0.0;
  XT2 = 0.0;

  double UK, WK, H0, H1, H2, S0, S1, S2, R0, R1, R2;
  for(int K = 0; K < NOSC[M-1]; K++)
  {
    UK = UI[M-1][K];
    WK = WRI[M-1][K];
    PINaT1(E,UK,WK,DELTA,WCCM,H0,H1,H2,S0,S1,S2,R0,R1,R2);
    SPH0[K] = F[M-1][K]*H0;
    SPH1[K] = F[M-1][K]*H1;
    SPH2[K] = F[M-1][K]*H2;
    SPS0[K] = F[M-1][K]*S0;
    SPS1[K] = F[M-1][K]*S1;
    SPS2[K] = F[M-1][K]*S2;
    SPT0[K] = F[M-1][K]*R0;
    SPT1[K] = F[M-1][K]*2.0*R1;
    SPT2[K] = F[M-1][K]*6.0*(R1-R2);
    XH0 = XH0+SPH0[K];
    XH1 = XH1+SPH1[K];
    XH2 = XH2+SPH2[K];
    XS0 = XS0+SPS0[K];
    XS1 = XS1+SPS1[K];
    XS2 = XS2+SPS2[K];
    XT0 = XT0+SPT0[K];
    XT1 = XT1+SPT1[K];
    XT2 = XT2+SPT2[K];
  }

}

