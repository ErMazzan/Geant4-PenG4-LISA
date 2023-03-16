
//  *********************************************************************
//                       SUBROUTINE RITAV
//  *********************************************************************
void PenInterface::RITAV(double &X, double &PDF, double &CDF)
{
  //  This subroutine gives the values of the (normalised) pdf, PDF, and of
  //  its cumulative distribution function, CDF, at the point X. These
  //  values are calculated from the RITA approximation.

  using namespace CRITA;

  if(X > XT[NPM1+1-1])
  {
    PDF = 0.0;
    CDF = 1.0;
  }
  else if(X < XT[0])
  {
    PDF = 0.0;
    CDF = 0.0;
  }
  else
  {
    int I = 1;
    int I1 = NPM1+1;
    bool brkIt = false;
    while(!brkIt)
	  {
	    brkIt = true;
	    int IT = (I+I1)/2;
	    if(X > XT[IT-1])
	    {
	      I = IT;
	    }
	    else
	    {
	      I1 = IT;
	    }
	    if(I1-I > 1){ brkIt = false; continue;}	  
	  }
    double TAU = (X-XT[I-1])/(XT[I+1-1]-XT[I-1]);
    double CON1 = 2.0*B[I-1]*TAU;
    double CON2 = 1.0+B[I-1]+A[I-1]*(1.0-TAU);
    double ETA;
    if(fabs(CON1) > 1.0E-10*fabs(CON2))
      {
	ETA = CON2*(1.0-sqrt(1.0-2.0*TAU*CON1/(CON2*CON2)))/CON1;
      }
    else
      {
	ETA = TAU/CON2;
      }
    PDF = DPAC[I-1]*((1.0+(A[I-1]+B[I-1]*ETA)*ETA)*(1.0+(A[I-1]+B[I-1]*ETA)*ETA))/((1.0-B[I-1]*ETA*ETA)*(1.0+A[I-1]+B[I-1])*(XT[I+1-1]-XT[I-1]));
    CDF = PAC[I-1]+ETA*DPAC[I-1];
  }
}

