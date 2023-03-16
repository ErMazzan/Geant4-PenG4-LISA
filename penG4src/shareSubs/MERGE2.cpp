
//  *********************************************************************
//                       SUBROUTINE MERGE2
//  *********************************************************************
void PenInterface::MERGE2(double* X1, double* Y1, double* X2, double* Y2, double* XM, double* YM, int &N1, int &N2, int &N)
{
  //  This subroutine merges two tables (X1,Y1), (X2,Y2) of two functions
  //  to produce a table (XM,YM) of the sum of these functions, with abs-
  //  cissas in increasing order. The abscissas and function values are
  //  assumed to be positive. N1, N2 and N are the numbers of grid points
  //  in the input and merged tables. A discontinuity in the function is
  //  described by giving twice the abscissa. Log-log linear interpolation
  //  is used to interpolate the input tables.

  using namespace PENERROR_mod;

  const double EPS = 1.0E-10;
  const int NP_S = 12000;
  const int NP2_S = NP_S + NP_S;
  double X[NP2_S], Y12[NP2_S];

  if(N1 > NP_S || N2 > NP_S)
  {
    printf("NP =%7d\n", N1 > N2 ? N1 : N2);
    ErrorFunction(1901); return;
  }

  SORT2(X1,Y1,N1);
  if(IRETRN != 0){ return;}
  SORT2(X2,Y2,N2);
  if(IRETRN != 0){ return;}

  for(int I1 = 0; I1 < N1; I1++)
  {
    X[I1] = X1[I1];
  }
  double XC, TST1, TST2, TST3, TST4, TST12, TST34, TST;
  N = N1;
  for(int I2 = 0; I2 < N2; I2++)
  {
    int I1; 
    XC = X2[I2];
    FINDI(X1,XC,N1,I1);
    if(I1 == N1){ I1 = N1-1;}
    TST1 = fabs(XC-X1[I1-1]);
    TST2 = fabs(XC-X1[I1+1-1]);
    TST12 = (TST1 < TST2 ? TST1 : TST2);
    if(I2+1 > 1)
    {
      TST3 = fabs(XC-X2[I2-1]);
    }
    else
    {
      TST3 = 1.0;
    }
    if(I2+1 < N2)
    {
      TST4 = fabs(XC-X2[I2+1]);
    }
    else
    {
      TST4 = 1.0;
    }
    TST34 = (TST3 < TST4 ? TST3 : TST4);
    TST = EPS*XC;
    if(TST34 > TST)
    {
      if(TST12 > TST)
      {
        N = N+1;
        X[N-1] = XC;
      }
    }
    else
    {
      N = N+1;
      X[N-1] = XC;
    }
  }

  //  ****  Sort and clean the merged grid.

  bool brkIt = false;
  while(!brkIt)
  {
    brkIt = true;
    for(int I = 0; I < N-1; I++)
    {
      int IMIN = I+1;
      double XMIN = X[I];
      for(int J = I+1; J < N; J++)
      {
        if(X[J] < XMIN)
        {
          IMIN = J+1;
          XMIN = X[J];
        }
      }
      double SAVE = X[I];
      X[I] = X[IMIN-1];
      X[IMIN-1] = SAVE;
    }

    for(int I = 0; I < N-2; I++)
    {
      if(X[I] > X[I+2]*(1.0E0-EPS))
      {
        X[I+1] = X[N-1];
        N = N-1;
        brkIt = false;
        break;
      }
    }
  }

  for(int I = 0; I < N; I++)
  {
    XC = X[I];
    if(I+1 < N)
    {
      if(X[I] > X[I+1]*(1.0-EPS)){ XC = X[I]*(1.0-EPS);}
    }
    if(I+1 > 1)
    {
      if(X[I] < X[I-1]*(1.0+EPS)){ XC = X[I]*(1.0+EPS);}
    }
    int J;
    double YI1, YI2;
    FINDI(X1,XC,N1,J);
    if(J == N1){ J = N1-1;}
    if(X1[J+1-1] > X1[J-1]+EPS)
    {
      YI1 = exp(log(Y1[J-1])+log(XC/X1[J-1])*log(Y1[J+1-1]/Y1[J-1])/log(X1[J+1-1]/X1[J-1]));
    }
    else
    {
      YI1 = Y1[J-1];
    }
    FINDI(X2,XC,N2,J);
    if(J == N2){ J = N2-1;}
    if(X2[J+1-1] > X2[J-1]+EPS)
    {
      YI2 = exp(log(Y2[J-1])+log(XC/X2[J-1])*log(Y2[J+1-1]/Y2[J-1])/log(X2[J+1-1]/X2[J-1]));
    }
    else
    {
      YI2 = Y2[J-1];
    }
    Y12[I] = YI1+YI2;
    if(Y12[I] < 1.0E-75){ Y12[I] = 1.0E-75;}
  }

  if(N > NP_S)
  {
    printf("NP = %7d\n", N);
    ErrorFunction(1902); return;
  }
  for(int I = 0; I < N; I++)
  {
    XM[I] = X[I];
    YM[I] = Y12[I];
  }     
}

