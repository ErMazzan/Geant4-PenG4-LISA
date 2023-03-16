#include <cmath>

//  *********************************************************************
//                       SUBROUTINE BRaAW
//  *********************************************************************
void PenInterface::BRaAW(double &ZEQ, FILE* IWR)
{
  //  This subroutine generates the parameters of the angular distribution
  //  of bremsstrahlung photons for the element of atomic number ZEQ. In
  //  the case of compounds (and mixtures) ZEQ is the average atomic number
  //  of the elements in the molecule. The evaluated parameters are written
  //  on the material definition file. Data are read from the database file
  //  'pdbrang.p18'.

  using namespace PENELOPE_mod;
  using namespace PENERROR_mod;

  const int NZ = 13;
  const int NE =  7;
  const int NK = 10;
  
  double Z[NZ], E[NE], XK[NK], P1[NZ][NE][NK], P2[NZ][NE][NK], Q1[NE][NK], Q2[NE][NK];
  double Y1[NZ], Y2[NZ], A[NZ], B[NZ], C[NZ], D[NZ], A2[NZ], B2[NZ];
  
  Z[0] = 1.00;
  Z[1] = 2.00;
  Z[2] = 5.00;
  Z[3] = 8.00;
  Z[4] = 11.0;
  Z[5] = 13.0;
  Z[6] = 26.0;
  Z[7] = 37.0;
  Z[8] = 47.0;
  Z[9] = 64.0;
  Z[10] = 79.0;
  Z[11] = 86.0;
  Z[12] = 92.0;
  
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

  //  ****  Read database file.
  int IZ, IE, IK, I;
  double ZR, ER, RKR, P1R, P2R;
  FILE* pdfiles = fopen("./pdfiles/pdbrang.p18", "r");
  char BUFFER[256];
  fgets(BUFFER,256,pdfiles);
  for(int IZ1 = 0; IZ1 < NZ; IZ1++)
  {
    for(int IE1 = 0; IE1 < NE; IE1++)
    {
      for(int IK1 = 0; IK1 < NK; IK1++)
      {
        fscanf(pdfiles, "%i %i %i %lf %lf %lf %lf %lf%*[^\n]", &IZ, &IE, &IK, &ZR, &ER, &RKR, &P1R, &P2R);
        getc(pdfiles);
        if((fabs(ZR-Z[IZ-1]) < 1.0E-6) && (fabs(ER-E[IE-1]) < 1.0E-6) && (fabs(RKR-XK[IK-1]) < 1.0E-6))
        {
          P1[IZ-1][IE-1][IK-1] = P1R;
          P2[IZ-1][IE-1][IK-1] = P2R;
        }
        else
        {
          printf("Corrupt data file (pdbrang.p18).\n");
          ErrorFunction(1326); return;
        }
      }
    }
  }
  fclose(pdfiles);
  
  //  ****  Interpolation in Z.

  for(IE = 0; IE < NE; IE++)
  {
    for(IK = 0; IK < NK; IK++)
    {
      for(IZ = 0; IZ < NZ; IZ++)
      {
        Y1[IZ] = log(P1[IZ][IE][IK]*Z[IZ]);
        Y2[IZ] = P2[IZ][IE][IK];
      }
      for(IZ = 0; IZ < NZ-1; IZ++)
      {
	A2[IZ] = Y2[IZ]-(Y2[IZ+1]-Y2[IZ])*Z[IZ]/(Z[IZ+1]-Z[IZ]);
	B2[IZ] = (Y2[IZ+1]-Y2[IZ])/(Z[IZ+1]-Z[IZ]);
      }
      A2[NZ-1]=A2[NZ-1-1];
      B2[NZ-1]=B2[NZ-1-1];
      SPLINE(Z,Y1,A,B,C,D,0.0,0.0,NZ);
      if(IRETRN != 0){ return;}
      FINDI(Z,ZEQ,NZ,I);
      Q1[IE][IK] = exp(A[I-1]+ZEQ*(B[I-1]+ZEQ*(C[I-1]+ZEQ*D[I-1])))/ZEQ;
      Q1[IE][IK] = std::min( Q1[IE][IK], 1.0);  // Corrects wrong values.
      Q2[IE][IK] = A2[I-1]+ZEQ*B2[I-1];
    }
  }

  //  ****  Write final table of parameters.

  int NDATA = 70;
  fprintf(IWR," *** Bremss angular distribution,  ZEQ =%12.5E,  NDATA =%4d\n", ZEQ, NDATA);
  for(IE = 0; IE < NE; IE++)
  {
    for(IK = 0; IK < NK; IK++)
    {
      fprintf(IWR, " %2d %2d %10.3E %10.3E %14.7E %14.7E\n", IE+1, IK+1, E[IE], XK[IK], Q1[IE][IK], Q2[IE][IK]);
    }
  }
}
