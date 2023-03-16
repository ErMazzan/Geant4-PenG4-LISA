
//  *********************************************************************
//                      FUNCTION RNDG3F
//  *********************************************************************
double RNDG3F(double X)
{
  //  Truncated Gaussian distribution, restricted to the interval (-3.0,
  //  3.0), with zero mean and unit variance.

  double RNDG3F_RETURN;
  if(fabs(X) > 3.00001)
  {
    RNDG3F_RETURN = 0.0;
  }
  else
  {
    RNDG3F_RETURN = exp(-X*X*0.5/(1.01538698*1.01538698));
  }
  return RNDG3F_RETURN;
}

