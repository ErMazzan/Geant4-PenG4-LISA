
//  *********************************************************************
//                       SUBROUTINE BRaAR
//  *********************************************************************
void PenInterface::BRaAR(int M, FILE* IRD, FILE* IWR, int INFO)
{
  //  This subroutine reads bremsstrahlung angular distribution parameters
  //  of material M from the material data file. It also initializes the
  //  algorithm for generation of the initial direction of bremsstrahlung
  //  photons.

  using namespace PENELOPE_mod;
  using namespace PENERROR_mod;
  // using namespace CEBR;
  using namespace CBRANG;

  const double TREV = 2.0*REV;

  const int NE = 7;
  const int NK = 10;
  double E[NE], XK[NK], Q1R[NE][NK], Q2R[NE][NK], Q1[NE][NKT], Q2[NE][NKT];
  double X[NK], A[NK], B[NK], C[NK], D[NK];
  double ZEQ;
  int NDATA;

  E[0] = 1.00E3;
  E[1] = 5.00E3;
  E[2] = 1.00E4;
  E[3] = 5.00E4;
  E[4] = 1.00E5;
  E[5] = 5.00E5;
  E[6] = 1.00E6;

  XK[0] = 0.00;
  XK[1] = 0.10;
  XK[2] = 0.20;
  XK[3] = 0.30;
  XK[4] = 0.40;
  XK[5] = 0.50;
  XK[6] = 0.60;
  XK[7] = 0.70;
  XK[8] = 0.80;
  XK[9] = 0.95;

  for(int IE = 0; IE < NET; IE++)
  {
    BET[IE] = sqrt(E[IE]*(E[IE]+TREV))/(E[IE]+REV);
  }

  //  ****  Grid of reduced photon energies.

  double DDK = 1.0/(double(NKT-1));
  for(int IK = 0; IK < NKT; IK++)
  {
    BK[IK] = (IK+1-1)*DDK;
  }

  //  ****  Read angular distribution parameters from file (unit IRD).

  fscanf(IRD, "%*40c%lf%*10c%4d%*[^\n]", &ZEQ, &NDATA);
  getc(IRD);
  
  if(INFO >= 2){ fprintf(IWR,"\n *** Bremss angular distribution,  ZEQ =%12.5E,  NDATA =%4d\n", ZEQ, NDATA);}
  if(NDATA != 70){ ErrorFunction(1324); return;}

  if(ZEQ > 1.0){ ZBEQ[M-1] = ZEQ;}
  else{ ZBEQ[M-1] = 1.0;}
  if(ZBEQ[M-1] > 99.0){ ZBEQ[M-1] = 99.0;}

  for(int IE1 = 0; IE1 < NE; IE1++)
  {
    for(int IK1 = 0; IK1 < NK; IK1++)
    {
      int IE, IK;
      double ER,RKR,Q1RR,Q2RR;
      fscanf(IRD, "%d %d %lf %lf %lf %lf%*[^\n]", &IE, &IK, &ER, &RKR, &Q1RR, &Q2RR);
      getc(IRD);
      if((fabs(ER-E[IE-1]) < 1.0E-6) && (fabs(RKR-XK[IK-1]) < 1.0E-6))
      {
        Q1R[IE-1][IK-1] = Q1RR;
        Q2R[IE-1][IK-1] = Q2RR;
      }
      else
      {
        printf("\nCorrupt data file (pdbrang.p08).\n");
        ErrorFunction(1325); return;
      }
    }
  }

  if(INFO >= 2)
  {
    for(int IE = 0; IE < NE; IE++)
    {
      for(int IK = 0; IK < NK; IK++)
      {
        fprintf(IWR, "%10.3E %10.3E %14.7E %14.7E\n", E[IE], XK[IK], Q1R[IE][IK], Q2R[IE][IK]);
      }
    }
  }

  //  ****  Expanded K-grid of distribution parameters.

  for(int IE = 0; IE < NE; IE++)
  {
    for(int IK = 0; IK < NK; IK++)
    {
      X[IK] = log(Q1R[IE][IK]);
    }
    SPLINE(XK,X,A,B,C,D,0.0,0.0,NK);
    if(IRETRN != 0){ return;}
    for(int IK = 0; IK < NKT; IK++)
    {
      int J;
      FINDI(XK,BK[IK],NK,J);
      Q1[IE][IK] = A[J-1]+BK[IK]*(B[J-1]+BK[IK]*(C[J-1]+BK[IK]*D[J-1]));
    }
    for(int IK = 0; IK < NK; IK++)
    {
      X[IK] = Q2R[IE][IK];
    }
    SPLINE(XK,X,A,B,C,D,0.0,0.0,NK);
    if(IRETRN != 0){ return;}
    for(int IK = 0; IK < NKT; IK++)
    {
      int J;
      FINDI(XK,BK[IK],NK,J);
      Q2[IE][IK] = A[J-1]+BK[IK]*(B[J-1]+BK[IK]*(C[J-1]+BK[IK]*D[J-1]));
    }
  }

  //  ****  ... and natural cubic spline interpolations.

  for(int IK = 0; IK < NKT; IK++)
  {
    for(int IE = 0; IE < NE; IE++)
    {
      X[IE] = Q1[IE][IK];
    }
    SPLINE(BET,X,A,B,C,D,0.0,0.0,NE);
    if(IRETRN != 0){ return;}
    for(int IE = 0; IE < NE; IE++)
    {
      BP1[M-1][IE][IK][0] = A[IE];
      BP1[M-1][IE][IK][1] = B[IE];
      BP1[M-1][IE][IK][2] = C[IE];
      BP1[M-1][IE][IK][3] = D[IE];
    }
    for(int IE = 0; IE < NE; IE++)
    {
      X[IE] = Q2[IE][IK];
    }
    SPLINE(BET,X,A,B,C,D,0.0,0.0,NE);
    if(IRETRN != 0){ return;}
    for(int IE = 0; IE < NE; IE++)
    {
      BP2[M-1][IE][IK][0] = A[IE];
      BP2[M-1][IE][IK][1] = B[IE];
      BP2[M-1][IE][IK][2] = C[IE];
      BP2[M-1][IE][IK][3] = D[IE];
    }
  }
}

