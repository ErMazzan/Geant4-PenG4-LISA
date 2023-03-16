
//  *********************************************************************
//                       SUBROUTINE FINDI
//  *********************************************************************
void FINDI(double* X, double XC, int N, int &I)
{
  //     This subroutine finds the interval (X(I),X(I+1)) that contains the
  //  value XC using binary search.

  //  Input:
  //     X(I) (I=1:N) ... grid points (the X values must be in increasing
  //                      order).
  //     XC ............. point to be located.
  //     N  ............. number of grid points.
  //  Output:
  //     I .............. interval index.

  if(XC > X[N-1])
  {
    I = N;
    return;
  }
  if(XC < X[0])
  {
    I = 1;
    return;
  }
  I = 1;
  int I1 = N;
  bool brkIt = false;
  while(!brkIt)
  {
    brkIt = true;
    int IT = (I+I1)/2;
    if(XC > X[IT-1])
    {
      I = IT;
    }
    else
    {
      I1 = IT;
    }
    if(I1-I > 1){ brkIt = false; continue;}
  }
}

