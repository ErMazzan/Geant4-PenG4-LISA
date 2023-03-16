
//  *********************************************************************
//                       SUBROUTINE SINTEG
//  *********************************************************************
void PenInterface::SINTEG(double* X, double* A, double* B, double* C, double* D, double &XL, double &XU, double &SUM, int N)
{

  //  ****  Set integration limits in increasing order.
  double XLL, XUU, SIGN;
  if(XU > XL)
  {
    XLL = XL;
    XUU = XU;
    SIGN = 1.0;
  }
  else
  {
    XLL = XU;
    XUU = XL;
    SIGN = -1.0;
  }
  //  ****  Check integral limits.
  if(XLL < X[0] || XUU > X[N-1])
  {
    printf("\n     Integral limits out of range. Stop.\n");
    ErrorFunction(1906); return;
  }
  //  ****  Find involved intervals.
  SUM = 0.0;
  int IL, IU;
  FINDI(X,XLL,N,IL);
  FINDI(X,XUU,N,IU);

  double X1, X2;
  if(IL == IU)
  {
    //  ****  Only a single interval involved.
    X1 = XLL;
    X2 = XUU;
    SUM = X2*(A[IL-1]+X2*((B[IL-1]/2)+X2*((C[IL-1]/3)+X2*D[IL-1]/4)))-X1*(A[IL-1]+X1*((B[IL-1]/2)+X1*((C[IL-1]/3)+X1*D[IL-1]/4)));
  }
  else
  {
    //  ****  Contributions from several intervals.
    X1 = XLL;
    X2 = X[IL+1-1];
    SUM = X2*(A[IL-1]+X2*((B[IL-1]/2)+X2*((C[IL-1]/3)+X2*D[IL-1]/4)))-X1*(A[IL-1]+X1*((B[IL-1]/2)+X1*((C[IL-1]/3)+X1*D[IL-1]/4)));
    IL = IL+1;
    for(int I = IL-1; I < IU;  I++)
    {
      X1 = X[I];
      X2 = X[I+1];
      if((I+1) == IU){ X2 = XUU;}
      double SUMP = X2*(A[I]+X2*((B[I]/2)+X2*((C[I]/3)+X2*D[I]/4)))-X1*(A[I]+X1*((B[I]/2)+X1*((C[I]/3)+X1*D[I]/4)));
      SUM = SUM+SUMP;
    }
  }
  SUM = SIGN*SUM;     
}

