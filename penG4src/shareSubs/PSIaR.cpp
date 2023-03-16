
//  *********************************************************************
//                       SUBROUTINE PSIaR
//  *********************************************************************
void PenInterface::PSIaR(int M, FILE *IRD, FILE* IWR, int INFO)
{
  //  This subroutine reads cross sections for inner-shell ionisation by
  //  positron impact of the elements in material M and prepares simulation
  //  tables.

  using namespace PENELOPE_mod;
  using namespace PENERROR_mod;

////////////////  using namespace CEGRID;
  using namespace COMPOS;
  using namespace CECUTR;
  using namespace CADATA;
  using namespace CPSI0;

  char CS5[16][6];
  const int NDIN=800;
  double E[NDIN], XPSIR[NDIN][16], X[NDIN], Y[NDIN];

  strcpy(CS5[0], "CS-K ");
  strcpy(CS5[1], "CS-L1");
  strcpy(CS5[2], "CS-L2");
  strcpy(CS5[3], "CS-L3");
  strcpy(CS5[4], "CS-M1");
  strcpy(CS5[5], "CS-M2");
  strcpy(CS5[6], "CS-M3");
  strcpy(CS5[7], "CS-M4");
  strcpy(CS5[8], "CS-M5");
  strcpy(CS5[9], "CS-N1");
  strcpy(CS5[10], "CS-N2");
  strcpy(CS5[11], "CS-N3");
  strcpy(CS5[12], "CS-N4");
  strcpy(CS5[13], "CS-N5");
  strcpy(CS5[14], "CS-N6");
  strcpy(CS5[15], "CS-N7");

  //  ************  Read element x-section tables

  int NSHR, IZZ, NDATA;
  for(int IEL = 0; IEL < NELEM[M-1]; IEL++)
  {
    fscanf(IRD, "%*46c%3d%*11c%3d%*10c%4d%*[^\n]", &IZZ, &NSHR, &NDATA);
    getc(IRD);

    if(INFO >= 2){ fprintf(IWR, "\n *** Positron impact ionisation cross sections,  IZ =%3d,  NSHELL =%3d,  NDATA =%4d\n", IZZ, NSHR, NDATA);}

    if(IZZ != IZ[M-1][IEL]){ ErrorFunction(1310); return;}
    if(NDATA > NDIN){ ErrorFunction(1311); return;}
    if(NSHR > 16){ ErrorFunction(1312); return;}
    for(int IE = 0; IE < NDATA; IE++)
    {
      fscanf(IRD, "%lf", &E[IE]);
      for(int IS = 0; IS < NSHR; IS++)
      {
        if(IS<NSHR-1){fscanf(IRD, "%lf", &XPSIR[IE][IS]);}else{
        fscanf(IRD, "%lf%*[^\n]",&XPSIR[IE][IS]);
        getc(IRD);
      }
      }
    }

      //  ****  Remove shells with ionisation energies less than 50 eV.

    int NSHA;
    if(NSHR > 1)
    {
      NSHA = NSHR;
      bool brkIt = false;
      for(int IS = NSHA-1; IS >= 0; IS--)
      {
        if(EB[IZZ-1][IS] < 50.0)
        {
          NSHR = NSHR-1;
        }
        else
        {
          brkIt = true;
          break;
        }
      }
      if(NSHR < 1 && !brkIt){ NSHR = 1;}
    }

    double TCS;
    if(INFO >= 2)
    {
      fprintf(IWR, "\n   Energy");
      for(int IS = 0; IS < NSHR; IS++)
      {
        fprintf(IWR, "       %s", CS5[IS]);
      }
      fprintf(IWR,"\n");
      for(int IE = 0; IE < NDATA; IE++)
      { 
        TCS = 0.0;
        for(int IS = 0; IS < NSHR; IS++)
        {
          TCS = TCS+XPSIR[IE][IS];
        }
        fprintf(IWR, " %11.5E", E[IE]);
        for(int IS = 0; IS < NSHR; IS++)
        {
          fprintf(IWR, " %11.5E", XPSIR[IE][IS]);
        }
        fprintf(IWR, " %11.5E\n", TCS);
      }
    }

    double XC, DX;
    int N, IC;
    NSPSI[IZZ-1] = NSHR;
    if(IPSIF[IZZ-1] == 0)
    {
      IPSIF[IZZ-1] = NCURP+1;
      if(NCURP+NEGP > NRP)
      {
        fprintf(IWR, "\nInsufficient memory storage in PSIaR.");
        fprintf(IWR,"\nIncrease the value of the parameter NRP to %d\n", NCURP+NEGP);
        ErrorFunction(1313); return;
      }
      for(int IS = 0; IS < NSHR; IS++)
      {
        N = 0;
        for(int I = 0; I < NDATA; I++)
        {
          if(XPSIR[I][IS] > 1.0E-35)
          {
            N = N+1;
            X[N-1] = log(E[I]);
            if(N > 1)
            {
              if(X[N-1] < X[N-1-1]+1.0E-6){ X[N-1] = X[N-1-1]+1.0E-6;}
            }
            Y[N-1] = log(XPSIR[I][IS]);
          }
        }
        if(N > 4)
        {
          for(int I = 0; I < NEGP; I++)
          {
            IC = NCURP+I+1;
            XC = CEGRID_.DLEMP[I];
            if(XC > X[0])
            {
              int J;
              FINDI(X,XC,N,J);
              if(J == N){ J = N-1;}
              DX = X[J+1-1]-X[J-1];
              if(DX > 1.0E-6)
              {
                XPSI[IC-1][IS] = Y[J-1]+(XC-X[J-1])*(Y[J+1-1]-Y[J-1])/DX;
              }
              else
              {
                XPSI[IC-1][IS] = (Y[J+1-1]+Y[J-1])/2.0;
              }
            }
            else
            {
              XPSI[IC-1][IS] = -80.6;
            }
          }
        }
        else
        {
          for(int I = 0; I < NEGP; I++)
          {
            IC = NCURP+I+1;
            XPSI[IC-1][IS] = -80.6;
          }
        }
      }
      NCURP = NCURP+NEGP;
      IPSIL[IZZ-1] = NCURP;
    }
  }
}

