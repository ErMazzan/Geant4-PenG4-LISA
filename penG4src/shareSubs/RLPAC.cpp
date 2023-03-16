
//  *********************************************************************
//                       SUBROUTINE RLPAC
//  *********************************************************************
void PenInterface::RLPAC(double* X, double* PDF, double* PAC, int NPpar)
{
  //  Cumulative distribution function of PDF(X)/X, obtained from linear
  //  interpolation on a table of PDF values.
  //  The independent variable X is assumed to take only positive values.

  //    X ..... array of values of the variable (in increasing order).
  //    PDF ... corresponding PDF values.
  //    PAC ... cumulative probability function.
  //    NP .... number of points in the table.


  const double EPS = 1.0E-35;
  PAC[0] = 0.0;
  double X1, Y1, X2, Y2, DX, DY, B, A, DS;
  for(int I = 1; I < NPpar; I++)
  {
    if(X[I-1] < EPS){ X1 = EPS;}
    else{ X1 = X[I-1];}
    Y1 = PDF[I-1];

    if(X[I] < EPS){ X2 = EPS;}
    else{ X2 = X[I];}
      
    Y2 = PDF[I];
    DX = X2-X1;
    DY = Y2-Y1;
    B = DY/DX;
    A = Y1-B*X1;
    DS = A*log(X2/X1)+B*(X2-X1);
    PAC[I] = PAC[I-1]+DS;
  }
     
}

