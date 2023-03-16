
//  *********************************************************************
//                       SUBROUTINE IRND0
//  *********************************************************************
void PenInterface::IRND0(double* W, double* F, int* K, int &N)
{
  //  Initialisation of Walker's aliasing algorithm for random sampling
  //  from discrete probability distributions.
  //
  //  Input arguments:
  //    N ........ number of different values of the random variable.
  //    W(1:N) ... corresponding point probabilities (not necessarily
  //               normalised to unity).
  //  Output arguments:
  //    F(1:N) ... cutoff values.
  //    K(1:N) ... alias values.

  //  ****  Renormalisation.
  using namespace CRITAN;

  CNORM = 0.0;
  for(int I = 0; I < N; I++)
  {
    CNORM = CNORM+(W[I]>0.0 ? W[I] : 0.0);
  }
  CNORM = 1.0/CNORM;
  double FACT = double(N)*CNORM;
  for(int I = 0; I < N; I++)
  {
    K[I] = I+1;
    F[I] = (W[I]>0.0 ? W[I] : 0.0)*FACT;
  }
  if(N == 1){return;}
  //  ****  Cutoff and alias values.
  for(int I = 0; I < N-1; I++)
  {
    double HLOW = 1.0;
    double HIGH = 1.0;
    int ILOW = 0;
    int IHIGH = 0;
    for(int J = 0; J < N; J++)
    {
      if(K[J] == J+1)
      {
        if(F[J] < HLOW)
        {
          HLOW = F[J];
          ILOW = J+1;
        }
        else if(F[J] > HIGH)
        {
          HIGH = F[J];
          IHIGH = J+1;
        }
      }
    }
    if(ILOW == 0 || IHIGH == 0){ return;}
    K[ILOW-1] = IHIGH;
    F[IHIGH-1] = HIGH+HLOW-1.0;
  }
}

