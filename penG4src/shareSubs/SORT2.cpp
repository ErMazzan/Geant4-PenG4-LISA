
//  *********************************************************************
//                       SUBROUTINE SORT2
//  *********************************************************************
void PenInterface::SORT2(double* X, double* Y, int &N)
{
  //  This subroutine sorts a table (X,Y) of a function with n data points.
  //  A discontinuity of the function is described by giving twice the abs-
  //  cissa. It is assumed that the function is strictly positive (negative
  //  values of Y are set to zero).

  using namespace PENERROR_mod;

  const int NP_S=12000;
  int IORDER[NP_S];
  
  if(N > NP_S)
  {
    printf("NP =%7d\n", N);
    ErrorFunction(1903); return;
  }

  if(N == 1){ return;}
  for(int I = 0; I < N; I++)
  {
    IORDER[I] = I+1;
    if(Y[I] < 1.0E-75){ Y[I] = 1.0E-75;}
  }

  int IMIN;
  double XMIN;
  for(int I = 0; I < N-1; I++)
  {
    IMIN = I+1;
    XMIN = X[I];
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
    SAVE = Y[I];
    Y[I] = Y[IMIN-1];
    Y[IMIN-1] = SAVE;
    int ISAVE = IORDER[I];
    IORDER[I] = IORDER[IMIN-1];
    IORDER[IMIN-1] = ISAVE;
    if(I+1 == 1){ continue;}
    if(IORDER[I] < IORDER[I-1] && fabs(X[I]-X[I-1]) < 1.0E-15)
    {
      SAVE = X[I-1];
      X[I-1] = X[I];
      X[I] = SAVE;
      SAVE = Y[I-1];
      Y[I-1] = Y[I];
      Y[I] = SAVE;
      ISAVE = IORDER[I-1];
      IORDER[I-1] = IORDER[I];
      IORDER[I] = ISAVE;
    }
  }
  int I = N;
  if(IORDER[I-1] < IORDER[I-1-1] && fabs(X[I-1]-X[I-1-1]) < 1.0E-15)
  {
    double SAVE = X[I-1-1];
    X[I-1-1] = X[I-1];
    X[I-1] = SAVE;
    SAVE = Y[I-1-1];
    Y[I-1-1] = Y[I-1];
    Y[I-1] = SAVE;
  }
}

