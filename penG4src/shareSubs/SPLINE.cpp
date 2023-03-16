
//  *********************************************************************
//                       SUBROUTINE SPLINE
//  *********************************************************************
void PenInterface::SPLINE(double* X, double* Y, double* A, double* B, double* C, double* D, double S1, double SN, int N)
{
  //     Cubic spline interpolation of tabulated data.
  //
  //  Input:
  //     X(I) (I=1:N) ... grid points (the X values must be in increasing
  //                      order).
  //     Y(I) (I=1:N) ... corresponding function values.
  //     S1,SN .......... second derivatives at X(1) and X(N). The natural
  //                      spline corresponds to taking S1=SN=0.
  //     N .............. number of grid points.
  //  Output:
  //     A(I),B(I),C(I),D(I) (I=1:N) ... spline coefficients.
  //
  //  The interpolating cubic polynomial in the I-th interval, from X(I) to
  //  X(I+1), is
  //               P(x) = A(I)+x*(B(I)+x*(C(I)+x*D(I)))

  //  Reference: M.J. Maron, 'Numerical Analysis: a Practical Approach',
  //             MacMillan Publ. Co., New York, 1982.

  using namespace PENERROR_mod;

  if(N < 4)
  {
    printf(" *** Error in SPLINE: interpolation cannot be performed with %4d points.\n", N);
    ErrorFunction(1904); return;
  }
  int N1 = N-1;
  int N2 = N-2;
  int K;
      //  ****  Auxiliary arrays H(=A) and DELTA(=D).
  for(int I = 0; I < N1; I++)
  {
    A[I] = X[I+1]-X[I];
    if(A[I] < 1.0E-12*( (fabs(X[I]) > fabs(X[I+1])) ? fabs(X[I]) : fabs(X[I+1]) ))
    {
       printf(" *** Error in SPLINE: X values not in increasing order.\n");
       ErrorFunction(1905); return;
    }
    D[I] = (Y[I+1]-Y[I])/A[I];
  }
      //  ****  Symmetric coefficient matrix (augmented).
  for(int I = 0; I < N2; I++)
  {
    B[I] = 2.0*(A[I]+A[I+1]);
    K = N1-(I+1)+1;
    D[K-1] = 6.0*(D[K-1]-D[K-2]);
  }
  D[1] = D[1]-A[0]*S1;
  D[N1-1] = D[N1-1]-A[N1-1]*SN;
      //  ****  Gauss solution of the tridiagonal system.
  for(int I = 1; I < N2; I++)
  {
    double R = A[I]/B[I-1];
    B[I] = B[I]-R*A[I];
    D[I+1] = D[I+1]-R*D[I];
  }
      //  ****  The SIGMA coefficients are stored in array D.
  D[N1-1] = D[N1-1]/B[N2-1];
  for(int I = 1; I < N2; I++)
  {
    K = N1-(I+1)+1;
    D[K-1] = (D[K-1]-A[K-1]*D[K+1-1])/B[K-2];
  }
  D[N-1] = SN;
      //  ****  Spline coefficients.
  double SI1 = S1;
  for(int I = 0; I < N1; I++)
  {
    double SI = SI1;
    SI1 = D[I+1];
    double H = A[I];
    double HI = 1.0/H;
    A[I] = (HI/6.0)*(SI*(X[I+1]*X[I+1]*X[I+1])-SI1*(X[I]*X[I]*X[I]))+HI*(Y[I]*X[I+1]-Y[I+1]*X[I])+(H/6.0)*(SI1*X[I]-SI*X[I+1]);
    B[I] = (HI/2.0)*(SI1*(X[I]*X[I])-SI*(X[I+1]*X[I+1]))+HI*(Y[I+1]-Y[I])+(H/6.0)*(SI-SI1);
    C[I] = (HI/2.0)*(SI*X[I+1]-SI1*X[I]);
    D[I] = (HI/6.0)*(SI1-SI);
  }
      //  ****  Natural cubic spline for X.GT.X(N).
  double FNP = B[N1-1]+X[N-1]*(2.0*C[N1-1]+X[N-1]*3.0*D[N1-1]);
  A[N-1] = Y[N-1]-X[N-1]*FNP;
  B[N-1] = FNP;
  C[N-1] = 0.0;
  D[N-1] = 0.0;

}

