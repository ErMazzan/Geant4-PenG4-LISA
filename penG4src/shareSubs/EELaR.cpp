
//  *********************************************************************
//                       SUBROUTINE EELaR
//  *********************************************************************
void PenInterface::EELaR(int M, FILE* IRD, FILE* IWR, int INFO)
{
  //  This subroutine reads elastic cross sections for electrons and posi-
  //  trons in material M from the material data file. It also initializes
  //  the algorithm for simulation of elastic scattering of electrons and
  //  positrons.


  using namespace PENELOPE_mod;
  using namespace PENERROR_mod;
  using namespace COMPOS;
////////////////  using namespace CEGRID;
  using namespace CEIMFP;
  using namespace CPIMFP;
  using namespace CEEL00;
  using namespace CEINTF;

  int NDATA;

  //  ****  Reading the input cross section table.
  
  fscanf(IRD, "%*59c%d%*[^\n]", &NDATA);
  getc(IRD);
  
  if(INFO >= 2)
  {
    fprintf(IWR, "\n *** Electron and positron elastic cross sections,  NDATA =%4d\n", NDATA);
    fprintf(IWR, "\n  Energy       CS0,e-      CS1,e-      CS2,e-      CS0,e+      CS1,e+      CS2,e+\n   (eV)        (cm**2)     (cm**2)     (cm**2)     (cm**2)     (cm**2)     (cm**2)\n ------------------------------------------------------------------------------------\n");
  }
  for(int I = 0; I < NDATA; I++)
  {
    fscanf(IRD, "%lf %lf %lf %lf %lf %lf %lf%*[^\n]", &EJT[I], &XE0[I], &XE1[I], &XE2[I], &XP0[I], &XP1[I], &XP2[I]);
    getc(IRD);
    
    if(INFO >= 2){ fprintf(IWR,"%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E\n", EJT[I], XE0[I], XE1[I], XE2[I], XP0[I], XP1[I], XP2[I]);}
    EJTL[I] = log(EJT[I]);
  }

  //  ****  Elastic scattering of electrons.

  double EC;
  int J;
  for(int I = 0; I < NDATA; I++)
  {
    FJL[I] = log(XE0[I]);
  }
  SPLINE(EJTL,FJL,A,B,C,D,0.0,0.0,NDATA);
  if(IRETRN != 0){ return;}
  for(int I = 0; I < NEGP; I++)
  {
    EC = CEGRID_.DLEMP[I];
    FINDI(EJTL,EC,NDATA,J);
    XE0[I] = exp(A[J-1]+EC*(B[J-1]+EC*(C[J-1]+EC*D[J-1])));
  }
  for(int I = 0; I < NDATA; I++)
  {
    FJL[I] = log(XE1[I]);
  }
  SPLINE(EJTL,FJL,A,B,C,D,0.0,0.0,NDATA);
  if(IRETRN != 0){ return;}
  for(int I = 0; I < NEGP; I++)
  {
    EC = CEGRID_.DLEMP[I];
    FINDI(EJTL,EC,NDATA,J);
    XE1[I] = exp(A[J-1]+EC*(B[J-1]+EC*(C[J-1]+EC*D[J-1])));
  }

  for(int I = 0; I < NDATA; I++)
  {
    FJL[I] = log(XE2[I]);
  }
  SPLINE(EJTL,FJL,A,B,C,D,0.0,0.0,NDATA);
  if(IRETRN != 0){ return;}
  for(int I = 0; I < NEGP; I++)
  {
    EC = CEGRID_.DLEMP[I];
    FINDI(EJTL,EC,NDATA,J);
    XE2[I] = exp(A[J-1]+EC*(B[J-1]+EC*(C[J-1]+EC*D[J-1])));
  }

  double XS0, XS1, XS2, XS0H, XS1S, XS2S;
  double AAA, BBB;
  for(int I = 0; I < NEGP; I++)
  {
    XS0 = XE0[I];
    XS1 = XE1[I];
    XS2 = XE2[I];
    double FPEL = 1.0/(XS0*VMOL[M-1]);
    double FPT1 = 1.0/(XS1*VMOL[M-1]);
    double FPST = CEGRID_.ET[I]/(CSTPE[M-1][I]+RSTPE[M-1][I]);
      
    XS0H = C1[M-1]*FPT1;
    if(XS0H > C2[M-1]*FPST){ XS0H = C2[M-1]*FPST;}
    if(XS0H < FPEL){ XS0H = FPEL;}
    XS0H = 1.0/(VMOL[M-1]*XS0H);

    double RNDC;
    EELa0(XS0,XS1,XS2,XS0H,AAA,BBB,RNDC,XS1S,XS2S);
    SEHEL[M-1][I] = XS0H*VMOL[M-1];
    RNDCE[M-1][I] = RNDC;
    AE[M-1][I] = AAA;
    BE[M-1][I] = BBB;
    T1E0[I] = XS1S;
    T1E[M-1][I] = T1EI[I]+XS1S*VMOL[M-1];
    T2E0[I] = XS2S;
    T2E[M-1][I] = T2EI[I]+XS2S*VMOL[M-1];
  }

  //  ****  Print electron elastic scattering tables.

  if(INFO >= 3)
  {
    fprintf(IWR, "\n PENELOPE >>>  Elastic scattering of electrons\n");
    fprintf(IWR, "\n   E (eV)      MFP (mtu)   TMFP1 (mtu)  MFPh (mtu)        A           B           RNDC\n ------------------------------------------------------------------------------------------\n");
  }
  double FP0, FP1;
  double HMFP;
  for(int I = 0; I < NEGP; I++)
  {
    FP0 = RHO[M-1]/(XE0[I]*VMOL[M-1]);
    FP1 = RHO[M-1]/(XE1[I]*VMOL[M-1]);
    HMFP = RHO[M-1]/SEHEL[M-1][I];
    if(INFO >= 3){ fprintf(IWR, "%12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E\n", CEGRID_.ET[I], FP0, FP1, HMFP, AE[M-1][I], BE[M-1][I], RNDCE[M-1][I]);}
    SEHEL[M-1][I] = log(SEHEL[M-1][I]);
    AE[M-1][I] = log(AE[M-1][I]);
      //  ****  Soft scattering events are switched off when T1E is too small.
    if(T1E[M-1][I] > 1.0E-6*XE1[I]*VMOL[M-1])
    {
      if(T1E[M-1][I] < 1.0E-35){ T1E[M-1][I] = 1.0E-35;}
      T1E[M-1][I] = log(T1E[M-1][I]);
      if(T2E[M-1][I] < 1.0E-35){ T2E[M-1][I] = 1.0E-35;}
      T2E[M-1][I] = log(T2E[M-1][I]);
    }
    else
    {
      T1E[M-1][I] = -100.0;
      T2E[M-1][I] = -100.0;
    }
  }

  //  ****  Elastic scattering of positrons.

  for(int I = 0; I < NDATA; I++)
  {
    FJL[I] = log(XP0[I]);
  }
  SPLINE(EJTL,FJL,A,B,C,D,0.0,0.0,NDATA);
  if(IRETRN != 0){ return;}
  for(int I = 0; I < NEGP; I++)
  {
    EC = CEGRID_.DLEMP[I];
    FINDI(EJTL,EC,NDATA,J);
    XP0[I] = exp(A[J-1]+EC*(B[J-1]+EC*(C[J-1]+EC*D[J-1])));
  }

  for(int I = 0; I < NDATA; I++)
  {
    FJL[I] = log(XP1[I]);
  }
  SPLINE(EJTL,FJL,A,B,C,D,0.0,0.0,NDATA);
  if(IRETRN != 0){ return;}
  for(int I = 0; I < NEGP; I++)
  {
    EC = CEGRID_.DLEMP[I];
    FINDI(EJTL,EC,NDATA,J);
    XP1[I] = exp(A[J-1]+EC*(B[J-1]+EC*(C[J-1]+EC*D[J-1])));
  }

  for(int I = 0; I < NDATA; I++)
  {
    FJL[I] = log(XP2[I]);
  }
  SPLINE(EJTL,FJL,A,B,C,D,0.0,0.0,NDATA);
  if(IRETRN != 0){ return;}
  for(int I = 0; I < NEGP; I++)
  {
    EC = CEGRID_.DLEMP[I];
    FINDI(EJTL,EC,NDATA,J);
    XP2[I] = exp(A[J-1]+EC*(B[J-1]+EC*(C[J-1]+EC*D[J-1])));
  }

  for(int I = 0; I < NEGP; I++)
  {
    XS0 = XP0[I];
    XS1 = XP1[I];
    XS2 = XP2[I];
    double FPEL = 1.0/(XS0*VMOL[M-1]);
    double FPT1 = 1.0/(XS1*VMOL[M-1]);
    double FPST = CEGRID_.ET[I]/(CSTPP[M-1][I]+RSTPP[M-1][I]);

    XS0H = C1[M-1]*FPT1;
    if(XS0H > C2[M-1]*FPST){ XS0H = C2[M-1]*FPST;}
    if(XS0H < FPEL){ XS0H = 1.0/(VMOL[M-1]*FPEL);}
    else{ XS0H = 1.0/(VMOL[M-1]*XS0H);}
      
    double RNDC;
    EELa0(XS0,XS1,XS2,XS0H,AAA,BBB,RNDC,XS1S,XS2S);
    SPHEL[M-1][I] = XS0H*VMOL[M-1];
    RNDCP[M-1][I] = RNDC;
    AP[M-1][I] = AAA;
    BP[M-1][I] = BBB;
    T1P0[I] = XS1S;
    T1P[M-1][I] = T1PI[I]+XS1S*VMOL[M-1];
    T2P0[I] = XS2S;
    T2P[M-1][I] = T2PI[I]+XS2S*VMOL[M-1];
  }

  //  ****  Print positron elastic scattering tables.

  if(INFO >= 3)
  {
    fprintf(IWR, "\n PENELOPE >>>  Elastic scattering of positrons\n");
    fprintf(IWR, "\n   E (eV)      MFP (mtu)   TMFP1 (mtu)  MFPh (mtu)        A           B           RNDC\n ------------------------------------------------------------------------------------------\n");
  }
  for(int I = 0; I < NEGP; I++)
  {
    FP0 = RHO[M-1]/(XP0[I]*VMOL[M-1]);
    FP1 = RHO[M-1]/(XP1[I]*VMOL[M-1]);
    HMFP = RHO[M-1]/SPHEL[M-1][I];
    if(INFO >= 3){ fprintf(IWR, "%12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E\n", CEGRID_.ET[I], FP0, FP1, HMFP, AP[M-1][I], BP[M-1][I], RNDCP[M-1][I]);}
    SPHEL[M-1][I] = log(SPHEL[M-1][I]);
    AP[M-1][I] = log(AP[M-1][I]);
      
      //  ****  Soft scattering events are switched off when T1P is too small.
    if(T1P[M-1][I] > 1.0E-6*XP1[I]*VMOL[M-1])
    {
      if(T1P[M-1][I] < 1.0E-35){ T1P[M-1][I] = log(1.0E-35);}
      else{ T1P[M-1][I] = log(T1P[M-1][I]);}
    
      if(T2P[M-1][I] < 1.0E-35){ T2P[M-1][I] = log(1.0E-35);}
      else{ T2P[M-1][I] = log(T2P[M-1][I]);}
    }
    else
    {
      T1P[M-1][I] = -100.0;
      T2P[M-1][I] = -100.0;
    }
  }
}

