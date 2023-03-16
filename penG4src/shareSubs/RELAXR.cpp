
//  *********************************************************************
//                       SUBROUTINE RELAXR
//  *********************************************************************
void PenInterface::RELAXR(FILE* IRD, FILE* IWR, int INFO)
{
  //  This subroutine reads atomic relaxation data, for a single element,
  //  from the material definition file (unit IRD) and initializes the
  //  algorithm for random sampling of de-excitation cascades of this
  //  element. (See heading comments in subroutine RELAX).

  using namespace PENERROR_mod;
  using namespace CADATA;
  using namespace CRELAX;
  
  char CH5[6];

  const int NTRAN = 2500;
  
  double PR[2500], ER[2500], WW[2500], FF[2500];
  int JS0[2500], JS1[2500], JS2[2500], KK[2500], IORD[2500], ISR[2500], IQQ[30];
  double EE[30],ALWR[30],CP0P[30];

  char LSHELL[31][3] = {"  ","K ","L1","L2","L3","M1","M2","M3","M4","M5","N1","N2","N3","N4","N5","N6","N7","O1","O2","O3","O4","O5","O6","O7","P1","P2","P3","P4","P5","Q1","X "};

  MODER = 1;  // RELAX normal operation mode.

  //  ****  Input transition data.

  int IZ, NSHR, NT;
  fscanf(IRD, "%*16c%3d%*18c%3d%*23c%5d%*[^\n]", &IZ, &NSHR, &NT);
  getc(IRD);
  if(INFO >= 2){ fprintf(IWR, "\n *** RELAX:  Z =%3d,  no. of shells =%3d,  no. of transitions =%5d", IZ, NSHR, NT);}

  if(NT > NTRAN){ ErrorFunction(1340); return;}
  if(NCUR+NT > NRX)
  {
    fprintf(IWR, "Insufficient memory storage in RELAXR.\n");
    fprintf(IWR, "Increase the value of the parameter NRX to %d\n", NCUR+NT);
    ErrorFunction(1341); return;
  }
  
  if(INFO >= 2){ fprintf(IWR, "\n\n  i   Shell    f    Ui (eV)    Gamma(1/eV)  lifetime (s)    Ji(0)\n --------------------------------------------------------------------\n");}

  for(int IS = 0; IS < NSHR; IS++)
  {
    fscanf(IRD, "%d %5s %d %lf %lf %lf%*[^\n]", &ISR[IS], CH5, &IQQ[IS], &EE[IS], &ALWR[IS], &CP0P[IS]);
    getc(IRD);
    double ALTIME;
    if(ALWR[IS] > 0.0)
    {
      ALTIME = HBAR/ALWR[IS];
    }
    else
    {
      ALTIME = 0.0;
    }
    if(INFO >= 2){ fprintf(IWR, " %3d %5s %2s  %1d %12.5E %12.5E %12.5E %12.5E\n", ISR[IS], CH5, LSHELL[ISR[IS]], IQQ[IS], EE[IS], ALWR[IS], ALTIME, CP0P[IS]);}
  }

  if(NT > 0)
  {
    if(INFO >= 2){ fprintf(IWR, "\n  S0 S1 S2   Probability     Energy (eV)\n ----------------------------------------\n");}
    for(int I = 0; I < NT; I++)
    {
      fscanf(IRD, "%d %d %d %lf %lf%*[^\n]", &JS0[I], &JS1[I], &JS2[I], &PR[I], &ER[I]);
      getc(IRD);
      if(INFO >= 2){ fprintf(IWR, "  %2s %2s %2s %15.8E %15.8E\n", LSHELL[JS0[I]], LSHELL[JS1[I]], LSHELL[JS2[I]], PR[I], ER[I]);}
      if(PR[I] < 1.0E-35)
      {
        if(INFO < 2){ fprintf(IWR, "  %2s %2s %2s %15.8E %15.8E\n", LSHELL[JS0[I]],LSHELL[JS1[I]], LSHELL[JS2[I]], PR[I], ER[I]);}
        ErrorFunction(1342); return;
      }
    }
  }

  //  ****  Check if this element's data have already been loaded.

  if(IFIRST[IZ-1][0] != 0){ return;}

  NSHT[IZ-1] = NSHR;
  for(int IS = 0; IS < NSHR; IS++)
  {
    IKS[IZ-1][IS] = ISR[IS];
    IFI[IZ-1][ISR[IS]-1] = IQQ[IS];
    EB[IZ-1][ISR[IS]-1] = EE[IS];
    if(ALWR[IS] > 0.0)
    {
      ALW[IZ-1][ISR[IS]-1] = HBAR/ALWR[IS];
    }
    else
    {
      ALW[IZ-1][ISR[IS]-1] = 0.0;
    }
    CP0[IZ][ISR[IS]-1] = CP0P[IS];
  }
  if(NT == 0)
  {
    for(int IS = 0; IS < 16; IS++)
    {
      IFIRST[IZ-1][IS] = NCUR+1;
      ILAST[IZ-1][IS] = NCUR+1;
    //  ****  The array IS0 contains the alias values.
      IS0[NCUR+1-1] = NCUR+1;
      P[NCUR+1-1] = 1.0;
      F[NCUR+1-1] = 1.0;
      ET[NCUR+1-1] = 0.0;
      IS1[NCUR+1-1] = 1;
      IS2[NCUR+1-1] = 1;
      NCUR = NCUR+1;
    }
    return;
  }

  //  ****  Walker's aliasing.

  for(int IS = 0; IS < 16; IS++)
  {
    int N = 0;
    for(int J = 0; J < NT; J++)
    {
      if(JS0[J] == IS+1)
      {
        N = N+1;
        IORD[N-1] = J+1;
        WW[N-1] = PR[J];
      }
    }
    if(N > 1)
    {
      IRND0(WW,FF,KK,N);
      IFIRST[IZ-1][IS] = NCUR+1;
      ILAST[IZ-1][IS] = NCUR+N;
      for(int LLL = 0; LLL < N; LLL++)
      {
        P[NCUR+LLL] = WW[LLL];
        F[NCUR+LLL] = FF[LLL];
        ET[NCUR+LLL] = ER[IORD[LLL]-1];
        //  ****  The array IS0 contains the alias values.
        IS0[NCUR+LLL] = IFIRST[IZ-1][IS]+KK[LLL]-1;
        IS1[NCUR+LLL] = JS1[IORD[LLL]-1];
        IS2[NCUR+LLL] = JS2[IORD[LLL]-1];
      }
      NCUR = NCUR+N;
    }
    else
    {
      NCUR = NCUR+1;
      IFIRST[IZ-1][IS] = NCUR;
      ILAST[IZ-1][IS] = NCUR;
      IS0[NCUR+1-1] = NCUR;
      P[NCUR-1] = 1.0;
      F[NCUR-1] = 1.0;
      ET[NCUR-1] = ER[0];
      IS1[NCUR-1] = JS1[0];
      IS2[NCUR-1] = JS2[0];
    }
  }

  //  ****  Verify that transition probabilities are correctly reproduced.

  double TST = 0.0;
  double PT, PPI;
  int IO, IN;
  for(int IS = 0; IS < 16; IS++)
  {
    IO = IFIRST[IZ-1][IS];
    IN = ILAST[IZ-1][IS];
    PT = 0.0;
    for(int I = IO-1; I < IN; I++)
    {
      PT = PT+P[I];
    }
    for(int I = IO-1; I < IN; I++)
    {
      PPI = 0.0;
      for(int J = IO-1; J < IN; J++)
      {
        if(IS0[J] == I+1){ PPI = PPI+(1.0-F[J]);}
      }
      PPI = (PPI+F[I])/double(IN-IO+1);
      if(TST < fabs(1.0-PPI*PT/P[I])){ TST = fabs(1.0-PPI*PT/P[I]);}
    }
  }
  if(TST > 1.0E-12)
  {
    printf("\nTST =%13.6E",TST);
    ErrorFunction(1343); return;
  }
}

