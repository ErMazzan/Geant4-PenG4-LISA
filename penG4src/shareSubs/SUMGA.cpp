
//  *********************************************************************
//                       FUNCTION SUMGA
//  *********************************************************************
double PenInterface::SUMGA(FCT_SUMGA FCT, double XL, double XU, double TOL)
{
  //  This function calculates the value SUMGA of the integral of the
  //  (external) function FCT over the interval (XL,XU) using the 20-point
  //  Gauss quadrature method with an adaptive-bisection scheme.

  //  TOL is the tolerance, i.e. maximum allowed relative error; it should
  //  not be less than 1.0D-13. A warning message is written in unit 6 when
  //  the required accuracy is not attained. The common block CSUMGA can be
  //  used to transfer the error flag IERGA and the number of calculated
  //  function values to the calling program.

  //                                        Francesc Salvat. 2 April, 2012.

  using namespace CSUMGA;
  using namespace CSGAWR;
  
  const int NP_S = 10;
  const int NP2_S = 2*NP_S;
  const int NP4 = 4*NP_S;
  const int NOIT = 128;
  const int NOIT5 = NOIT/5;
  const int NCALLT = 100000;
  
  double XM[NP_S], XP[NP_S];
  double S[NOIT], SN[NOIT], XR[NOIT], XRN[NOIT];
  //  Output error codes:
  //     IERGA = 0, no problem, the calculation has converged.
  //           = 1, too many open subintervals.
  //           = 2, too many function calls.
  //           = 3, subintervals are too narrow.

  //  ****  Gauss 20-point integration formula.
  //  Abscissas.
  double X[NP_S] = {7.6526521133497334E-02,2.2778585114164508E-01,3.7370608871541956E-01,5.1086700195082710E-01,6.3605368072651503E-01,7.4633190646015079E-01,8.3911697182221882E-01,9.1223442825132591E-01,9.6397192727791379E-01,9.9312859918509492E-01};
  //  Weights.
  double W[NP_S] = {1.5275338713072585E-01,1.4917298647260375E-01,1.4209610931838205E-01,1.3168863844917663E-01,1.1819453196151842E-01,1.0193011981724044E-01,8.3276741576704749E-02,6.2672048334109064E-02,4.0601429800386941E-02,1.7614007139152118E-02};

  for(int I = 0; I < NP_S; I++)
  {
    XM[I] = 1.0-X[I];
    XP[I] = 1.0+X[I];
  }
  //  ****  Global and partial tolerances.

  double TOL1; // Global tolerance.

  if(TOL < 1.0E-13){ TOL1 = 1.0E-13;}
  else{ TOL1 = TOL;}

  if( TOL1 > 1.0E-5){ TOL1 = 1.0E-5;}
  
  double TOL2 = TOL1;  // Effective tolerance.
  double TOL3 = 1.0E-13;  // Round-off protection.
  double SUMGA_RETURN = 0.0; //Valor que torna la funcio
  IERGA = 0;
  //  ****  Straight integration from XL to XU.
  double H = XU-XL;
  double HH = 0.5*H;
  double X1 = XL;
  double SP = W[0]*(FCT(X1+XM[0]*HH)+FCT(X1+XP[0]*HH));
  for(int J = 1; J < NP_S; J++)
  {
    SP = SP+W[J]*(FCT(X1+XM[J]*HH)+FCT(X1+XP[J]*HH));
  }
  S[0] = SP*HH;
  XR[0] = X1;
  NCALL = NP2_S;
  int NOI = 1;
  int IDONE = 1;  // To prevent a compilation warning.

    //  ****  Adaptive-bisection scheme.

  bool brkIt = false;
  int NOIP;
  double SUMR;
  while(!brkIt)
  {
    brkIt = true;
    H = HH;  // Subinterval length.
    HH = 0.5*H;
    double AHH = fabs(HH);
    if(TOL2 > 0.01*TOL1){ TOL2 = TOL2*0.5;}
    SUMR = 0.0;
    NOIP = NOI;
    NOI = 0;
    bool brkIt2 = false;
    for(int I = 0; I < NOIP; I++)
    {
      double SI = S[I];  // Bisect the I-th open interval.
      
      X1 = XR[I];
      if(AHH < fabs(X1)*TOL3){ IERGA = 3;}  // The interval is too narrow.
      SP = W[0]*(FCT(X1+XM[0]*HH)+FCT(X1+XP[0]*HH));
      for(int J = 1; J < NP_S; J++)
      {
        SP = SP+W[J]*(FCT(X1+XM[J]*HH)+FCT(X1+XP[J]*HH));
      }
      double S1 = SP*HH;
      
      double X2 = X1+H;
      if(AHH < fabs(X2)*TOL3){ IERGA = 3;}  // The interval is too narrow.
      SP = W[0]*(FCT(X2+XM[0]*HH)+FCT(X2+XP[0]*HH));
      for(int J = 1; J < NP_S; J++)
      {
        SP = SP+W[J]*(FCT(X2+XM[J]*HH)+FCT(X2+XP[J]*HH));
      }
      double S2 = SP*HH;
      
      IDONE = I+1;
      NCALL = NCALL+NP4;
      double S12 = S1+S2;  // Sum of integrals on the two subintervals.
      if(fabs(S12-SI) < ((TOL2*fabs(S12) > 1.0E-35) ? TOL2*fabs(S12) : 1.0E-35))
      {
      //  ****  The integral over the parent interval has converged.
        SUMGA_RETURN = SUMGA_RETURN+S12;
      }
      else
      {
        SUMR = SUMR+S12;
        NOI = NOI+2;
        if(NOI < NOIT)
        {
         //  ****  Store open intervals.
          SN[NOI-2] = S1;
          XRN[NOI-2] = X1;
          SN[NOI-1] = S2;
          XRN[NOI-1] = X2;
        }
        else
        {
        //  ****  Too many open intervals.
          IERGA = 1;
          brkIt2 = true;
          break;
        }
      }
      if(NCALL > NCALLT)
      {
       //  ****  Too many calls to FCT.
        IERGA = 2;
        brkIt2 = true;
        break;
      }
    }
    if(brkIt2){ break;}

  //  ****  Analysis of partial results and error control.

    if(IERGA == 3)  // Intervals are too narrow.
    {
      if(NOI < NOIT5)
      {
        IERGA = 0;  // The result is probably correct.
        SUMGA_RETURN = SUMGA_RETURN+SUMR;
        return SUMGA_RETURN;
      }
      break;
    }

    if(IERGA == 0)
    {
      double Aux_Double = TOL1*fabs(SUMGA_RETURN+SUMR);
      if(Aux_Double < 1.0E-35){ Aux_Double = 1.0E-35;}
      if(fabs(SUMR) < Aux_Double || NOI == 0)
      {
        return SUMGA_RETURN;
      }
      else
      {
        for(int I = 0; I < NOI; I++)
        {
          S[I] = SN[I];
          XR[I] = XRN[I];
        }
        brkIt = false;
        continue;
      }
    }
  }

    //  ****  Warning (low accuracy) message.

  if(IDONE < NOIP)
  {
    for(int I = IDONE; I < NOIP; I++)
    {
      SUMR = SUMR+S[I];
    }
    NOI = NOI+(NOIP-IDONE);
  }
  SUMGA_RETURN = SUMGA_RETURN+SUMR;
  if(ISGAW == 0){ return SUMGA_RETURN;}
  printf("  >>> SUMGA. Gauss adaptive-bisection quadrature.\n");
  printf("  XL =%15.8E, XU =%15.8E, TOL =%8.1E\n", XL, XU, TOL);
  if(fabs(SUMGA_RETURN) > 1.0E-35)
  {
    double RERR = fabs(SUMR)/fabs(SUMGA_RETURN);
    printf("  SUMGA =%22.15E, relative error =%8.1E\n", SUMGA_RETURN, RERR);
  }
  else
  {
    double AERR = fabs(SUMR);
    printf("  SUMGA =%22.15E, absolute error =%8.1E\n", SUMGA_RETURN, AERR);
  }
  printf("  NCALL =%6d, open subintervals =%4d, H =%10.3E\n", NCALL, NOI, HH);
  if(IERGA == 1)
  {
    printf("  IERGA = 1, too many open subintervals.\n");
  }
  else if(IERGA == 2)
  {
    printf("  IERGA = 2, too many function calls.\n");
  }
  else if(IERGA == 3)
  {
    printf("  IERGA = 3, subintervals are too narrow.\n");
  }
  printf("  WARNING: the required accuracy has not been attained.\n");

  return SUMGA_RETURN; //MIRAR
  
}

