
//  *********************************************************************
//                       SUBROUTINE RITAI0
//  *********************************************************************
void PenInterface::RITAI0(PDF_RITA PDF, double XLOW, double XHIGH, int N, int NU, double &ERRM, int fileidx, int &IER)
{
  //  Initialisation of the RITA algorithm for random sampling of a
  //  continuous random variable X from a probability distribution function
  //  PDF(X) defined in the interval (XLOW,XHIGH). The external function
  //  PDF(X) --not necessarily normalised-- must be provided by the user.
  //  N is the number of points in the sampling grid. These points are
  //  determined by means of an adaptive strategy that minimizes local
  //  interpolation errors. The first NU grid points are uniformly spaced
  //  in (XLOW,XHIGH); when NU is negative, the initial grid consists of
  //  -NU points logarithmically spaced (in this case, XLOW must be
  //  nonnegative).
  //
  //  ERRM is a measure of the interpolation error (the largest value of
  //  the absolute error of the rational interpolation integrated over each
  //  grid interval).
  //
  //  IER is an error code. An output value different from zero indicates
  //  that the initialization procedure has not succeeded; the reason is
  //  written in the standard output unit (6).
  //
  //
  //  ****  Interpolation coefficients and PDF tables are printed on
  //        separate files (UNIT=IWR) if IWR is greater than zero.
  //
  //  Other subprograms needed: EXTERNAL function PDF.
  //

  
  const double EPS = 1.0E-10;
  const double ZERO = 1.0E-75;
  const double ZEROT = 0.1*ZERO;

  //
  //  The information used by the sampling function RITAI is exported
  //  through the following common block,
  using namespace CRITA;

  //  where
  //    XT(I) ..... grid points, in increasing order.
  //    PAC(I) .... value of the cumulative pdf at XT(I).
  //    DPAC(I) ... probability of the I-th interval.
  //    A(I), B(I) ... rational inverse cumulative distribution parameters.
  //    IL(I) .... largest J for which PAC(J) < (I-1)/(NP-1).
  //    IU(I) .... smallest J for which PAC(J) > I/(NP-1).
  //    NPM1 ..... numner of grid points minus one, NP-1 (8.LE.NP.LE.NM).

  using namespace CRITAN;
  //    CNORM .... normalising constant of the external PDF, output.
  //
  double ERR[NM], C[NM];
  const int NIP = 51;
  double XS[NIP], YS[NIP], SUMI[NIP];
  int NUNIF, NP_RITA;
  char Dirname[]="RITAI0";

  static int naccess=0;
  if(fileidx>0){naccess++;}
  if(fileidx>0 && naccess==1)
  {
    char Command[100];
    strcpy(Command,"mkdir -p ");
    strcat(Command,Dirname);
    system(Command);   
  }

  if(N < 9)
  {
    printf(" Error in RITAI0: N must be larger than 8.\n N=%11d\n", N);
    IER=1;
    ErrorFunction(1402); exit(1402);
  }
  if(N > NM)
  {
    printf(" Error in RITAI0: N must be less than NM=512.\n N=%11d\n", N);
    IER=2;
    ErrorFunction(1403); exit(1403);
  }
  if(XLOW > XHIGH-EPS)
  {
    printf(" Error in RITAI0: XLOW must be larger than XHIGH.\n XLOW=%13.6E, XHIGH =%13.6E\n", XLOW, XHIGH);
    IER=3;
    ErrorFunction(1404); exit(1404);
  }
  IER=0;

  //  ****  We start with a grid of NUNIF points uniformly spaced in the
  //        interval (XLOW,XHIGH).

  if(NU >= 0)
  {
    NUNIF = (NU > 8 ? NU : 8);
    if(NUNIF > N/2){ NUNIF = N/2;}
      
    NP_RITA = NUNIF;
    double DX = (XHIGH-XLOW)/double(NP_RITA-1);
    XT[0] = XLOW;
    for(int I = 0; I < NP_RITA-1; I++)
    {
      XT[I+1] = XLOW+(I+1)*DX;     
    }
    XT[NP_RITA-1] = XHIGH;
  }
  else
  {
      //  ****  If NU.LT.0,the NUNIF points are logarithmically spaced.
      //        XLOW must be greater than or equal to zero.
    NUNIF = (-NU > 8 ? -NU : 8);
    if(NUNIF > N/2){ NUNIF = N/2;}
      
    NP_RITA = NUNIF;
    if(XLOW < 0.0)
    {
      printf(" Error in RITAI0: XLOW and NU are negative.\n XLOW=%14.7E,  NU=%11d\n", XLOW, NU);
      IER=4;
      ErrorFunction(1405); exit(1405);
    }
    XT[0] = XLOW;
    int I1;
    if(XLOW < 1.0E-16)
    {
      XT[1] = XLOW+1.0E-6*(XHIGH-XLOW);
      I1 = 2;
    }
    else
    {
      I1 = 1;
    }
    double FX = exp(log(XHIGH/XT[I1-1])/double(NP_RITA-I1));
    for(int I = I1-1; I < NP_RITA-1; I++)
    {
      XT[I+1] = XT[I]*FX;
    }
    XT[NP_RITA-1] = XHIGH;
  }

  for(int I = 0; I < NP_RITA-1; I++)
  {
    double DX = XT[I+1]-XT[I];
    double DXI = DX/double(NIP-1);
    double PDFMAX = 0.0;
    for(int K = 0; K < NIP; K++)
    {
      XS[K] = XT[I]+double(K+1-1)*DXI;

      YS[K] = PDF(XS[K]);
      if(YS[K] < ZEROT){ YS[K] = ZEROT;}

      if(PDFMAX < YS[K]){ PDFMAX = YS[K];}
    
    }
      //  ****  Simpson's integration.
    double CONS = DXI*3.3333333333333333E-1*0.5;
    SUMI[0] = 0.0;
    for(int K = 1; K < NIP; K++)
    {
      double XIH = XS[K]-0.5*DXI;

      double YSH = PDF(XIH);
      if(YSH < ZEROT){ YSH = ZEROT;}

      if(PDFMAX < YSH){ PDFMAX = YSH;}

      SUMI[K] = SUMI[K-1]+CONS*(YS[K-1]+4.0*YSH+YS[K]);
    }

    DPAC[I] = SUMI[NIP-1];
    double FACT = 1.0/DPAC[I];
    for(int K = 0; K < NIP; K++)
    {
      SUMI[K] = FACT*SUMI[K];
    }
      //  ****  When the PDF vanishes at one of the interval end points, its
      //        value is modified.
    if(YS[0] < ZERO){ YS[0] = 1.0E-5*PDFMAX;}
    if(YS[NIP-1] < ZERO){ YS[NIP-1] = 1.0E-5*PDFMAX;}

    double PLI = YS[0]*FACT;
    double PUI = YS[NIP-1]*FACT;
    B[I] = 1.0-1.0/(PLI*PUI*DX*DX);
    A[I] = (1.0/(PLI*DX))-1.0-B[I];
    C[I] = 1.0+A[I]+B[I];
    if(C[I] < ZERO)
    {
      A[I] = 0.0;
      B[I] = 0.0;
      C[I] = 1.0;
    }

      //  ****  ERR(I) is the integral of the absolute difference between the
      //  rational interpolation and the true PDF, extended over the interval
      //  (XT(I),XT(I+1)). Calculated using the trapezoidal rule.

    int ICASE = 1;
    bool brkIt = false;
    while(!brkIt)
    {
      brkIt = true;
      ERR[I] = 0.0;
      for(int K = 0; K < NIP; K++)
      {
        double RR = SUMI[K];
        double PAP = DPAC[I]*((1.0+(A[I]+B[I]*RR)*RR)*(1.0+(A[I]+B[I]*RR)*RR))/((1.0-B[I]*RR*RR)*C[I]*(XT[I+1]-XT[I]));
        if(K+1 == 1 || K+1 == NIP)
        {
          ERR[I] = ERR[I]+0.5*fabs(PAP-YS[K]);
        }
        else
        {
          ERR[I] = ERR[I]+fabs(PAP-YS[K]);
        }
      }
      ERR[I] = ERR[I]*DXI;
    //  ****  If ERR(I) is too large, the PDF is approximated by a uniform
    //        distribution.
      if(ERR[I] > 0.10*DPAC[I] && ICASE == 1)
      {
        B[I] = 0.0;
        A[I] = 0.0;
        C[I] = 1.0;
        ICASE = 2;
        brkIt = false;
        continue;
      }
    }
  }
  XT[NP_RITA-1] = XHIGH;
  A[NP_RITA-1] = 0.0;
  B[NP_RITA-1] = 0.0;
  C[NP_RITA-1] = 0.0;
  ERR[NP_RITA-1] = 0.0;
  DPAC[NP_RITA-1] = 0.0;

  //  ****  New grid points are added by halving the subinterval with the
  //        largest absolute error.

  bool brkIt = false;
  while(!brkIt)
  {
    brkIt = true;
    ERRM = 0.0;
    int LMAX = 1;
    for(int I = 0; I < NP_RITA-1; I++)
    {
    //  ****  ERRM is the largest of the interval errors ERR(I).
      if(ERR[I] > ERRM)
      {
        ERRM = ERR[I];
        LMAX = I+1;
      }
    }
      
    NP_RITA = NP_RITA+1;
    for(int I = NP_RITA-1; I >= LMAX+1-1; I--)
    {
      XT[I] = XT[I-1];
      A[I] = A[I-1];
      B[I] = B[I-1];
      C[I] = C[I-1];
      ERR[I] = ERR[I-1];
      DPAC[I] = DPAC[I-1];
    }
    XT[LMAX+1-1] = 0.5*(XT[LMAX-1]+XT[LMAX+2-1]);
    for(int I = LMAX-1; I < LMAX+1; I++)
    {
      double DX = XT[I+1]-XT[I];
      double DXI = (XT[I+1]-XT[I])/double(NIP-1);
      double PDFMAX = 0.0;
      for(int K = 0; K < NIP; K++)
      {
        XS[K] = XT[I]+double(K+1-1)*DXI;
        YS[K] = PDF(XS[K]);
        if(YS[K] < ZEROT){ YS[K] = ZEROT;}
        if(PDFMAX < YS[K]){ PDFMAX = YS[K];}
      }
    //  ****  Simpson's integration.
      double CONS = DXI*3.3333333333333333E-1*0.5;
      SUMI[0] = 0.0;
      for(int K = 1; K < NIP; K++)
      {
        double XIH = XS[K]-0.5*DXI;      
        double YSH = PDF(XIH);
        if(YSH < ZEROT){ YSH = ZEROT;}  
        SUMI[K] = SUMI[K-1]+CONS*(YS[K-1]+4.0*YSH+YS[K]);
      }

      DPAC[I] = SUMI[NIP-1];
      double FACT = 1.0/DPAC[I];
      for(int K = 0; K < NIP; K++)
      {
        SUMI[K] = FACT*SUMI[K];
      }

      if(YS[0] < ZERO){ YS[0] = 1.0E-5*PDFMAX;}
      if(YS[NIP-1] < ZERO){ YS[NIP-1] = 1.0E-5*PDFMAX;}
      double PLI = YS[0]*FACT;
      double PUI = YS[NIP-1]*FACT;
      B[I] = 1.0-1.0/(PLI*PUI*DX*DX);
      A[I] = (1.0/(PLI*DX))-1.0-B[I];
      C[I] = 1.0+A[I]+B[I];
      if(C[I] < ZERO)
      {
        A[I] = 0.0;
        B[I] = 0.0;
        C[I] = 1.0;
      }

      int ICASE = 1;
      bool brkIt2 = false;
      while(!brkIt2)
      {
        brkIt2 = true;
        ERR[I] = 0.0;
        for(int K = 0; K < NIP; K++)
        {
          double RR = SUMI[K];
          double PAP = DPAC[I]*((1.0+(A[I]+B[I]*RR)*RR)*(1.0+(A[I]+B[I]*RR)*RR))/((1.0-B[I]*RR*RR)*C[I]*(XT[I+1]-XT[I]));
          if(K+1 == 1 || K+1 == NIP)
          {
            ERR[I] = ERR[I]+0.5*fabs(PAP-YS[K]);
          }
          else
          {
            ERR[I] = ERR[I]+fabs(PAP-YS[K]);
          }
        }
        ERR[I] = ERR[I]*DXI;

        if(ERR[I] > 0.10*DPAC[I] && ICASE == 1)
        {
          B[I] = 0.0;
          A[I] = 0.0;
          C[I] = 1.0;
          ICASE = 2;
          brkIt2 = false;
          continue;
        }
      }
    }

    if(NP_RITA < N){ brkIt = false; continue;}
  }
  NPM1 = NP_RITA-1;

  //  ****  Renormalisation.

  CNORM = 0.0;
  for(int I = 0; I < NPM1; I++)
  {
    CNORM = CNORM+DPAC[I];
  }
  CNORM = 1.0/CNORM;
  ERRM = 0.0;
  for(int I = 0; I < NPM1; I++)
  {
    DPAC[I] = DPAC[I]*CNORM;
    ERR[I] = ERR[I]*CNORM;

    if(ERRM < ERR[I]){ ERRM = ERR[I];}
      
  }

  PAC[0] = 0.0;
  for(int I = 0; I < NPM1; I++)
  {
      PAC[I+1] = PAC[I]+DPAC[I];
  }
  PAC[NP_RITA-1] = 1.0;

  //  ****  Pre-calculated limits for the initial binary search in
  //        subroutine RITAI.

  double BIN = 1.0/double(NPM1);
  IL[0] = 1;
  for(int I = 1; I < NPM1; I++)
  {
    double PTST = (I+1-1)*BIN;
    for(int J = IL[I-1]-1; J < NP_RITA; J++)
    {
      if(PAC[J] > PTST)
      {
        IL[I] = J+1-1;
        IU[I-1] = J+1;
        break;
      }
    }
  }
  IU[NPM1-1] = NP_RITA;
  IL[NP_RITA-1] = NP_RITA-1;
  IU[NP_RITA-1] = NP_RITA;

  //  ****  Print interpolation tables (only when IWR.GT.0).

  if(fileidx > 0)
  {
    FILE* IWR;
    char Buffer[150];
    sprintf(Buffer,"%s%s%s%d%s",Dirname,"/","param",naccess,".dat");
    IWR = fopen(Buffer, "w");
    fprintf(IWR, " #  Interpolation error =%15.7E\n", ERRM);
    fprintf(IWR, " # Normalising constant =%15.7E\n", CNORM);
    fprintf(IWR, " #     X           PDF(X)          A             B             C           error\n");
    for(int I = 0; I < NPM1; I++)
    {
      double PDFE = PDF(XT[I]);
      if(PDFE < ZEROT){ PDFE = ZEROT;}

      PDFE = PDFE*CNORM;

      fprintf(IWR, "%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E\n", XT[I], PDFE, A[I], B[I], C[I], ERR[I]);
    }
    fclose(IWR);

    sprintf(Buffer,"%s%s%s%d%s",Dirname,"/","table",naccess,".dat");
    IWR = fopen(Buffer, "w");
    fprintf(IWR, " #  Interpolation error =%15.7E\n", ERRM);
    fprintf(IWR, " # Normalising constant =%15.7E\n", CNORM);
    fprintf(IWR, " #      X             PDF_ex          PDF_ap           err\n");

    for(int I = 0; I < NPM1; I++)
    {
      double DX = (XT[I+1]-XT[I])/double(NIP-1);
      for(int K = 0; K < NIP; K += 5)
      {
        double XTAU = XT[I]+double(K+1-1)*DX;

        double P1 = PDF(XTAU)*CNORM;
        if(P1 < ZEROT){ P1 = ZEROT;}
        
        //  ****  Rational interpolation.
        double TAU = (XTAU-XT[I])/(XT[I+1]-XT[I]);
        double CON1 = 2.0*B[I]*TAU;
        double CON2 = C[I]-A[I]*TAU;
        double ETA;
        if(fabs(CON1) > 1.0E-10*fabs(CON2))
        {
          ETA = CON2*(1.0-sqrt(1.0-2.0*TAU*CON1/(CON2*CON2)))/CON1;
        }
        else
        {
          ETA = TAU/CON2;
        }
        double P2 = DPAC[I]*((1.0+(A[I]+B[I]*ETA)*ETA)*(1.0+(A[I]+B[I]*ETA)*ETA))/((1.0-B[I]*ETA*ETA)*C[I]*(XT[I+1]-XT[I]));
        fprintf(IWR, " %15.8E %15.8E %15.8E %15.8E\n", XTAU, P1, P2, (P1-P2)/P1);
      }
    }
    fclose(IWR);

    sprintf(Buffer,"%s%s%s%d%s",Dirname,"/","limits",naccess,".dat");
    IWR = fopen(Buffer, "w");    
    fprintf(IWR, " #  I      PAC(ITL)         (I-1)/NPM1           I/NPM1            PAC(ITU)            PAC(I)\n");
    for(int I = 0; I < NPM1; I++)
    {
      fprintf(IWR, "%5d %18.11E %18.11E %18.11E %18.11E %18.11E\n", I+1, PAC[IL[I]-1], (I+1-1)*BIN, (I+1)*BIN, PAC[IU[I]-1], PAC[I]);
      if(PAC[IL[I]-1] > (I+1-1)*BIN+EPS || PAC[IU[I]-1] < (I+1)*BIN-EPS)
      {
        fprintf(IWR, " #  WARNING: The first four values should be in increasing order.\n");
      }
    }
    fclose(IWR);
  }
}

