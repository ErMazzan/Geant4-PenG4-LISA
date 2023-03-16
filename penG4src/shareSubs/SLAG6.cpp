
//  *********************************************************************
//                       SUBROUTINE SLAG6
//  *********************************************************************
void PenInterface::SLAG6(double &H, double* Y, double* S, int N)
{
  //     Piecewise six-point Lagrange integration of a uniformly tabulated
  //  function.

  //  Input arguments:
  //     H ............ grid spacing.
  //     Y(I) (1:N) ... array of function values (ordered abscissas).
  //     N ............ number of data points.

  //  Output argument:
  //     S(I) (1:N) ... array of integral values defined as
  //                S(I)=INTEGRAL(Y) from X(1) to X(I)=X(1)+(I-1)*H.

  using namespace PENERROR_mod;

  if(N < 6){ ErrorFunction(1907); return;}
  double HR = H/1440.0;
  double Y1 = 0.0;
  double Y2 = Y[0];
  double Y3 = Y[1];
  double Y4 = Y[2];
  S[0] = 0.0;
  S[1] = HR*(475*Y2+1427*Y3-798*Y4+482*Y[3]-173*Y[4]+27*Y[5]);
  S[2] = S[1]+HR*(-27*Y2+637*Y3+1022*Y4-258*Y[3]+77*Y[4]-11*Y[5]);
  for(int I = 3; I < N-2; I++)
  {
    Y1 = Y2;
    Y2 = Y3;
    Y3 = Y4;
    Y4 = Y[I];
    S[I] = S[I-1]+HR*(11*(Y1+Y[I+2])-93*(Y2+Y[I+1])+802*(Y3+Y4));
  }
  double Y5 = Y[N-1-1];
  double Y6 = Y[N-1];
  S[N-1-1] = S[N-2-1]+HR*(-27*Y6+637*Y5+1022*Y4-258*Y3+77*Y2-11*Y1);
  S[N-1] = S[N-1-1]+HR*(475*Y6+1427*Y5-798*Y4+482*Y3-173*Y2+27*Y1);
}

