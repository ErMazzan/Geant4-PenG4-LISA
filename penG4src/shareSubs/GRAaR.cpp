
//  *********************************************************************
//                       SUBROUTINE GRAaR
//  *********************************************************************
void PenInterface::GRAaR(int M, FILE* IRD, FILE* IWR, int &INFO)
{
  //  This subroutine reads the squared molecular form factor and the DCS
  //  for Rayleigh scattering of photons in material M. These two functions
  //  are tabulated using the same grids for all materials.
  //     The random sampling of the scattering angle is performed using the
  //  RITA algorithm.

  using namespace PENELOPE_mod;
  using namespace PENERROR_mod;

////////////////////  using namespace CEGRID;
  using namespace COMPOS;
  using namespace CADATA;
  using namespace CGIMFP;
  using namespace CGRA00;
  using namespace CGRA01;
  using namespace CGRA02;
  using namespace CGRA03;
  using namespace CRITA;

  const double RREV = 1.0/REV;
  

  double Q[NQ], F[NQ], FFI[NQ], ER[NEX], AA[NQ], BB[NQ], CC[NQ], DD[NQ];
  const int NIP =51;
  double QI[51], FUN[51], SUM[51];

  int NQQ;
  
  fscanf(IRD, "%*32c%d%*8c%d%*[^\n]", &NQQ, &NE);
  getc(IRD);
  if(INFO >= 2){ fprintf(IWR, "\n *** Rayleigh scattering.  NQ = %3d,  NE = %4d\n", NQ, NE);}
  for(int I = 0; I < NQ; I++)
  {
    fscanf(IRD, "%lf %lf%*[^\n]", &Q[I], &FFI[I]);
    getc(IRD);
    if(I == 0){ FF0[M-1] = FFI[I]*FFI[I];}
    F[I] = log(FFI[I]*FFI[I]);
    FF[M-1][I] = FFI[I];
  }
  for(int I = 0; I < NE; I++)
  {
    fscanf(IRD, "%lf %lf%*[^\n]", &ER[I], &XSRA[M-1][I]);
    getc(IRD);
    XSRA[M-1][I] = log(XSRA[M-1][I]*VMOL[M-1]);
  }
  double XS1;
  if(ER[NE-1] < CEGRID_.EU)
  {
    XS1 = XSRA[M-1][NE-1]+(XSRA[M-1][NE-1]-XSRA[M-1][NE-1-1])*(log(CEGRID_.EU)-log(ER[NE-1-1]))/(log(ER[NE-1])-log(ER[NE-1-1]));
    XSRA[M-1][NE-1] = XS1;
    ER[NE-1] = CEGRID_.EU;
  }

  if(M == 1)
  {
    QQM = Q[NQ-1]*Q[NQ-1];
    for(int I = 0; I < NQ; I++)
    {
      QQ[I] = log(Q[I]*Q[I]);
    }
    for(int I = 0; I < NE; I++)
    {
      ERA[I] = log(ER[I]);
    }
    for(int I = 0; I < NEGP; I++)
    {
      int J;
      FINDI(ERA,CEGRID_.DLEMP[I],NE,J);
      IED[I] = J;
    }
    for(int I = 0; I < NEGP-1; I++)
    {
      IEU[I] = IED[I+1]+1;
      if(IEU[I] > NE){ IEU[I] = NE;}  
    }
    IEU[NEGP-1] = IED[NEGP-1]+1;
    if(IEU[NEGP-1] > NE){ IEU[NEGP-1] = NE;}
      
  }

  SPLINE(QQ,F,AA,BB,CC,DD,0.0,0.0,NQ);
  if(IRETRN != 0){ return;}
  for(int I = 0; I < NQ; I++)
  {
    AR[M-1][I] = AA[I];
    BR[M-1][I] = BB[I];
    CR[M-1][I] = CC[I];
    DR[M-1][I] = DD[I];
  }

  //  ****  Total cross section at the simulation grid points (slightly
  //        increased to simplify interpolation).

  double EE = CEGRID_.DLEMP[0];
  int J;
  FINDI(ERA,EE,NE,J);
  XS1 = XSRA[M-1][J-1]+(XSRA[M-1][J+1-1]-XSRA[M-1][J-1])*(EE-ERA[J-1])/(ERA[J+1-1]-ERA[J-1]);
  double XSMAX;
  for(int IE = 0; IE < NEGP-1; IE++)
  {
    XSMAX = XS1;
    int J1 = J+1;
    EE = CEGRID_.DLEMP[IE+1];
    FINDI(ERA,EE,NE,J);
    XS1 = XSRA[M-1][J-1]+(XSRA[M-1][J+1-1]-XSRA[M-1][J-1])*(EE-ERA[J-1])/(ERA[J+1-1]-ERA[J-1]);
    if(XSMAX < XS1){ XSMAX = XS1;}

    if(J1 < J)
    {
      for(int I = J1-1; I < J; I++)
      {
        if(XSMAX < XSRA[M-1][I]){ XSMAX = XSRA[M-1][I];}            
      }
    }
    SGRA[M-1][IE] = exp(XSMAX);
  }
  SGRA[M-1][NEGP-1] = SGRA[M-1][NEGP-2];

  if(INFO >= 2)
  {
    fprintf(IWR, "\n   Q/me*c     Form factor\n -------------------------\n");
    for(int I = 0; I < NQ; I++)
    {
      fprintf(IWR, "%12.5E%13.5E\n", Q[I], FFI[I]);
    }
    fprintf(IWR, "\n   Energy       CS-Rayl\n    (eV)        (cm**2)\n -------------------------\n");
    for(int I = 0; I < NE; I++)
    {
      fprintf(IWR, "%12.5E%13.5E\n", ER[I], exp(XSRA[M-1][I])/VMOL[M-1]);
    }
  }
  //
  //  ****  Initialisation of the RITA algorithm for random sampling of the
  //  squared momentum transfer from the squared molecular form factor.

  MM = M;
  double Q2MIN = 0.0;
  Q2MAX = 0.0;
  int NPT = NP;
  int NU = NPT/4;
  for(int I = 1; I < NQ; I++)
  {
    if(GRAaF2(Q[I]*Q[I]) > 1.0E-35){ Q2MAX = Q[I-1]*Q[I-1]; }
  }
  double ERRM;
  int IER;
  RITAI0(GRAaF2,Q2MIN,Q2MAX,NPT,NU,ERRM,0,IER);
  if(IER != 0)
  {
    printf("IER=%d\n",IER);
    printf("GRAaR. RITA initialisation error.\n");
    ErrorFunction(1327); return;
  }
  int NPI = (*NPM1I)+1;
  if(NPI != NP)
  {
    printf("GRAaR. RITA initialisation error.\n");
    printf("The number of fixed grid points is %d\n", NPI);
    printf("The required number of grid points was %d\n", NP);
    ErrorFunction(1327); return;
  }
  if(ERRM > 1.0E-5)
  {
    printf("GRAaR. RITA interpolation error is too large.\n");
    printf("The interpolation error is %E\n", ERRM);
    ErrorFunction(1328); return;
  }

  //  ****  Upper limit of the X2 interval for the PENELOPE grid energies.

  for(int IE = 0; IE < NEGP; IE++)
  {
    double QM = 2.0*CEGRID_.ET[IE]*RREV;
    double Q2M = QM*QM;
    int I;
    if(Q2M > QTI[0])
    {
      if(Q2M < QTI[NP-1])
      {
        I = 1;
        J = NPI;
        bool brkIt = false;
        while(!brkIt)
        {
          brkIt = true;
          int IT = (I+J)/2;
          if(Q2M > QTI[IT-1])
          {
            I = IT;
          }
          else
          {
            J = IT;
          }
          if(J-I > 1){ brkIt = false; continue;}
        }

        double Q1 = QTI[I-1];
        double Q2 = Q2M;
        double DQ = (Q2-Q1)/double(NIP-1);
        double TAU, CON1, CI, CON2;
        for(int K = 0; K < NIP; K++)
        {
          QI[K] = Q1+double(K+1-1)*DQ;
          TAU = (QI[K]-QTI[I-1])/(QTI[I+1-1]-QTI[I-1]);
          CON1 = 2.0*BI[I-1]*TAU;
          CI = 1.0+AI[I-1]+BI[I-1];
          CON2 = CI-AI[I-1]*TAU;
          double ETAP;
          if(fabs(CON1) > 1.0E-16*fabs(CON2))
          {
            ETAP = CON2*(1.0-sqrt(1.0-2.0*TAU*CON1/(CON2*CON2)))/CON1;
          }
          else
          {
            ETAP = TAU/CON2;
          }
          FUN[K] = DPACI[I-1]*pow(1.0+(AI[I-1]+BI[I-1]*ETAP)*ETAP,2)/((1.0-BI[I-1]*ETAP*ETAP)*CI*(QTI[I+1-1]-QTI[I-1]));
        }
        SLAG6(DQ,FUN,SUM,NIP);
        if(IRETRN != 0){ return;}
        PMAX[IE][M-1] = PACI[I-1]+SUM[NIP-1];
      }
      else
      {
        PMAX[IE][M-1] = 1.0;
      }
    }
    else
    {
      PMAX[IE][M-1] = PACI[0];
    }
  }

  for(int I = 0; I < NP; I++)
  {
    QRA[I][M-1] = QTI[I];
    PRA[I][M-1] = PACI[I];
    DPRA[I][M-1]= DPACI[I];
    ARA[I][M-1] = AI[I];
    BRA[I][M-1] = BI[I];
    ITLRA[I][M-1] = ITLI[I];
    ITURA[I][M-1] = ITUI[I];
  }
} 

