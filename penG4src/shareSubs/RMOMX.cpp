
//  *********************************************************************
//                       FUNCTION RMOMX
//  *********************************************************************
double PenInterface::RMOMX(double* X, double* PDF, double XD, double XU, int NPpar, int MOM)
{
  //     Calculation of momenta of a pdf, PDF(X), obtained from linear
  //  log-log interpolation of the input table. The independent variable X
  //  is assumed to take only positive values.

  //     X ........ array of variable values (in increasing order).
  //     PDF ...... corresponding PDF values (must be non-negative).
  //     NP ....... number of points in the table.
  //     XD, XU ... limits of the integration interval.
  //     MOM ...... moment order.
  //     RMOM = INTEGRAL (X**N)*PDF(X) dX over the interval (XD,XU).

  using namespace PENERROR_mod;
  
  const double EPS = 1.0E-12;
  const double ZERO = 1.0E-35;

  double RMOMX_RETURN=0.0; //Valor que torna la funcio
  
  if(NPpar < 2){ ErrorFunction(1908); return RMOMX_RETURN;}
  if(X[0] < 0.0 || PDF[0] < 0.0)
  {
    printf("X(1),PDF(1) = %E %E\n",X[0],PDF[0]);
    ErrorFunction(1909); return RMOMX_RETURN;
  }
  for(int I = 1; I < NPpar; I++)
  {
    if(X[I] < 0.0 || PDF[I] < 0.0)
    {
      printf("I,X(I),PDF(I) = %d %E %E\n",I+1,X[I],PDF[I]);
      ErrorFunction(1910); return RMOMX_RETURN;
    }
    if(X[I] < X[I-1]){ ErrorFunction(1911); return RMOMX_RETURN;}
  }

  double XLOW = (X[0] > XD ? X[0] : XD);
  if(XLOW < ZERO){ XLOW = ZERO;}
  double XUP = (X[NPpar-1] < XU ? X[NPpar-1] : XU);

  if(XLOW >= XUP)
  {
    printf("\n WARNING: XLOW is greater than XUP in RMOMX.");
    printf("\n XLOW =%E,   XUP =%E", XLOW, XUP);
    RMOMX_RETURN = 0.0;
    return RMOMX_RETURN;
  }

  int IL = 1;
  int IU = NPpar-1;
  for(int I = 0; I < NPpar-1; I++)
  {
    if(X[I] < XLOW){ IL = I+1;}
    if(X[I] < XUP){ IU = I+1;}
  }

  //  ****  A single interval.

  double XIL, XFL, YIL, YFL, X1, X2, DENOM, Y1, Y2, DXL, DYL, DSUM, AP1;
  if(IU == IL)
  {
    XIL = log((X[IL-1] > ZERO ? X[IL-1] : ZERO));
    XFL = log(X[IL+1-1]);
    YIL = log((PDF[IL-1] > ZERO ? PDF[IL-1] : ZERO));
    YFL = log((PDF[IL+1-1] > ZERO ? PDF[IL+1-1] : ZERO));
    X1 = XLOW;
    X2 = XUP;
    DENOM = XFL-XIL;
    if(fabs(DENOM) > ZERO)
    {
      Y1 = exp(YIL+(YFL-YIL)*(log(X1)-XIL)/DENOM)*pow(X1,MOM);
      Y2 = exp(YIL+(YFL-YIL)*(log(X2)-XIL)/DENOM)*pow(X2,MOM);
    }
    else
    {
      Y1 = exp(YIL)*pow(X1,MOM);
      Y2 = exp(YIL)*pow(X2,MOM);
    }
    DXL = log(X2)-log(X1);
    DYL = log((Y2 > ZERO ? Y2 : ZERO))-log((Y1 > ZERO ? Y1 : ZERO));
    if(fabs(DXL) > EPS*fabs(DYL))
    {
      AP1 = 1.0+(DYL/DXL);
      if(fabs(AP1) > EPS)
      {
        DSUM = (Y2*X2-Y1*X1)/AP1;
      }
      else
      {
        DSUM = Y1*X1*DXL;
      }
    }
    else
    {
      DSUM = 0.5*(Y1+Y2)*(X2-X1);
    }
    RMOMX_RETURN = DSUM;
    return RMOMX_RETURN;
  }

      //  ****  Multiple intervals.

  XIL = log((X[IL-1] > ZERO ? X[IL-1] : ZERO));
  XFL = log(X[IL+1-1]);
  YIL = log((PDF[IL-1] > ZERO ? PDF[IL-1] : ZERO));
  YFL = log((PDF[IL+1-1] > ZERO ? PDF[IL+1-1] : ZERO));
  X1 = XLOW;
  DENOM = XFL-XIL;
  if(fabs(DENOM) > ZERO)
  {
    Y1 = exp(YIL+(YFL-YIL)*(log(X1)-XIL)/DENOM)*pow(X1,MOM);
  }
  else
  {
    Y1 = exp(YIL)*pow(X1,MOM);
  }
  X2 = X[IL+1-1];
  Y2 = (PDF[IL+1-1] > ZERO ? PDF[IL+1-1] : ZERO)*pow(X2,MOM);
  DXL = log(X2)-log(X1);
  DYL = log((Y2 > ZERO ? Y2 : ZERO))-log((Y1 > ZERO ? Y1 : ZERO));
  if(fabs(DXL) > EPS*fabs(DYL))
  {
    AP1 = 1.0+(DYL/DXL);
    if(fabs(AP1) > EPS)
    {
      DSUM = (Y2*X2-Y1*X1)/AP1;
    }
    else
    {
      DSUM = Y1*X1*DXL;
    }
  }
  else
  {
    DSUM = 0.5*(Y1+Y2)*(X2-X1);
  }
  RMOMX_RETURN = DSUM;
      
  if(IU > IL+1)
  {
    for(int I = IL; I < IU-1; I++)
    {
      X1 = X[I];
      Y1 = (PDF[I] > ZERO ? PDF[I] : ZERO)*pow(X1,MOM);
      X2 = X[I+1];
      Y2 = (PDF[I+1] > ZERO ? PDF[I+1] : ZERO)*pow(X2,MOM);
      DXL = log(X2)-log(X1);
      DYL = log((Y2 > ZERO ? Y2 : ZERO))-log((Y1 > ZERO ? Y1 : ZERO));
      if(fabs(DXL) > EPS*fabs(DYL))
      {
        AP1 = 1.0+(DYL/DXL);
        if(fabs(AP1) > EPS)
        {
          DSUM = (Y2*X2-Y1*X1)/AP1;
        }
        else
        {
          DSUM = Y1*X1*DXL;
        }
      }
      else
      {
        DSUM = 0.5*(Y1+Y2)*(X2-X1);
      }
      RMOMX_RETURN = RMOMX_RETURN+DSUM;
    }
  }

  X1 = X[IU-1];
  Y1 = (PDF[IU-1] > ZERO ? PDF[IU-1] : ZERO)*pow(X1,MOM);
  XIL = log(X[IU-1]);
  XFL = log(X[IU+1-1]);
  YIL = log((PDF[IU-1] > ZERO ? PDF[IU-1] : ZERO));
  YFL = log((PDF[IU+1-1] > ZERO ? PDF[IU+1-1] : ZERO));
  X2 = XUP;
  DENOM = XFL-XIL;
  if(fabs(DENOM) > ZERO)
  {
    Y2 = exp(YIL+(YFL-YIL)*(log(X2)-XIL)/DENOM)*pow(X2,MOM);
  }
  else
  {
    Y2 = exp(YIL)*pow(X2,MOM);
  }
  DXL = log(X2)-log(X1);
  DYL = log((Y2 > ZERO ? Y2 : ZERO))-log((Y1 > ZERO ? Y1 : ZERO));
  if(fabs(DXL) > EPS*fabs(DYL))
  {
    AP1 = 1.0+(DYL/DXL);
    if(fabs(AP1) > EPS)
    {
      DSUM = (Y2*X2-Y1*X1)/AP1;
    }
    else
    {
      DSUM = Y1*X1*DXL;
    }
  }
  else
  {
    DSUM = 0.5*(Y1+Y2)*(X2-X1);
  }
  RMOMX_RETURN = RMOMX_RETURN+DSUM;

  return RMOMX_RETURN;
}

