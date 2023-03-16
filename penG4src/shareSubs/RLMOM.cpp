
//  *********************************************************************
//                       FUNCTION RLMOM
//  *********************************************************************
double PenInterface::RLMOM(double* X, double* FCT, double XC, int NPpar, int MOM)
{
  //  Calculation of the integral of (X**MOM)*FCT(X) over the interval from
  //  X(1) to XC, obtained by linear interpolation on a table of FCT.
  //  The independent variable X is assumed to take only positive values.

  //    X ....... array of values of the variable (in increasing order).
  //    FCT ..... corresponding FCT values.
  //    NPpar ...... number of points in the table.
  //    XC ...... upper limit of the integral, X(1).LE.XC.LE.X(NPpar).
  //    MOM ..... moment order (GE.-1).
  //    RLMOM ... integral of (X**MOM)*FCT(X) over the interval from X(1)
  //              to XC.

  using namespace PENERROR_mod;

  const double EPS = 1.0E-35;

  double RLMOM_RETURN  = 0.0;
  if(MOM < -1){ ErrorFunction(1319); return RLMOM_RETURN;}
  if(NPpar < 2){ ErrorFunction(1320); return RLMOM_RETURN;}
  if(X[0] < 0.0){ ErrorFunction(1321); return RLMOM_RETURN;}
  for(int I = 1; I < NPpar; I++)
  {
    if(X[I] < 0.0){ ErrorFunction(1322); return RLMOM_RETURN;}
    if(X[I] < X[I-1]){ ErrorFunction(1323); return RLMOM_RETURN;}
  }

  RLMOM_RETURN  = 0.0;
  if(XC < X[0]){ return RLMOM_RETURN;}
  int IEND = 0;

  double XT;
  if(XC < X[NPpar-1]){ XT = XC;}
  else{ XT = X[NPpar-1];}
  
  double X1, Y1, X2, Y2, XTC, DX, DY, A, B, DS;
  for(int I = 0; I < NPpar-1; I++)
  {
    if(X[I] < EPS){ X1 = EPS;}
    else{ X1 = X[I];}
      
    Y1 = FCT[I];

    if(X[I+1] < EPS){ X2 = EPS;}
    else{ X2 = X[I+1];}

    Y2 = FCT[I+1];
    if(XT < X2)
    {
      XTC = XT;
      IEND = 1;
    }
    else
    {
      XTC = X2;
    }
    DX = X2-X1;
    DY = Y2-Y1;
    if(fabs(DX) > 1.0E-14*fabs(DY))
    {
      B = DY/DX;
      A = Y1-B*X1;
      if(MOM == -1)
      {
        DS = A*log(XTC/X1)+B*(XTC-X1);
      }
      else
      {
        DS = A*(pow(XTC,(MOM+1))-pow(X1,(MOM+1)))/double(MOM+1)+B*(pow(XTC,(MOM+2))-pow(X1,(MOM+2)))/double(MOM+2);
      }
    }
    else
    {
      DS = 0.5*(Y1+Y2)*(XTC-X1)*pow(XTC,MOM);
    }
    RLMOM_RETURN = RLMOM_RETURN+DS;
    if(IEND != 0){ return RLMOM_RETURN;}
  }
  return RLMOM_RETURN;
}

