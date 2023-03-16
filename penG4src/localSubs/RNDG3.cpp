
//  *********************************************************************
//                       FUNCTION RNDG3
//  *********************************************************************
double PenPhys::RNDG3()
{
  //  This function delivers random values in the interval (-3.0,3.0)
  //  sampled from a truncated Gaussian distribution that has zero mean and
  //  unit variance. The sampling is performed by using the RITA method.

  using namespace CRNDG3;

  double RNDG3_RETURN;
  //  ****  Selection of the interval (Walker's aliasing).
  double RN = RAND(1.0)*NPM1+1.0;
  int K = (int)RN;
  double TST = RN-(double)K;
  double RR, D;
  int I;
  if(TST < F[K-1])
  {
    I = K;
    RR = TST;
    D = F[K-1];
  }
  else
  {
    I = KA[K-1];
    RR = TST-F[K-1];
    D = 1.0-F[K-1];
  }
  //  ****  Sampling from the rational inverse cumulative distribution.
  if(RR > 1.0E-12)
  {
    RNDG3_RETURN = X_[I-1]+((1.0+A[I-1]+B[I-1])*D*RR/(D*D+(A[I-1]*D+B[I-1]*RR)*RR))*(X_[I+1-1]-X_[I-1]);
  }
  else
  {
    RNDG3_RETURN = X_[I-1]+RAND(2.0)*(X_[I+1-1]-X_[I-1]);
  }
  return RNDG3_RETURN;
}

