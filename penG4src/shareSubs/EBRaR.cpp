
//  *********************************************************************
//                       SUBROUTINE EBRaR
//  *********************************************************************
void PenInterface::EBRaR(double &WCRM, int M, FILE* IRD, FILE* IWR, int INFO)
{
  //  This subroutine reads the bremss scaled cross section for electrons
  //  in material M from the material data file. It computes restricted
  //  integrated cross sections and initializes the algorithm for simula-
  //  tion of bremss emission by electrons and positrons.

  using namespace PENELOPE_mod;
  using namespace PENERROR_mod;

///////////////////  using namespace CEGRID;
  using namespace CEBR;
  using namespace CEBR01;
  using namespace CEBR02;

  double A[NEGP], B[NEGP], C[NEGP], D[NEGP], PAC[NEGP], PDF[NEGP], ZBR;
  int    NBER;

  double WB0[NBW] = {1.0E-12, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.925, 0.95, 0.97, 0.99, 0.995, 0.999, 0.9995, 0.9999, 0.99995, 0.99999, 1.0};

  //  ****  Reading the scaled cross section table.

  fscanf(IRD, "%*45c%lf%*10c%4d%*[^\n]", &ZBR, &NBER);
  getc(IRD);
  
  if(INFO >= 2){ fprintf(IWR, "\n *** Electron scaled bremss x-section,  ZEQ =%12.5E,  NDATA =%4d\n", ZBR, NBER);fflush(IWR);}
  if(NBER != NBE){ ErrorFunction(1316); return;}
  ZBR2[M-1] = ZBR*ZBR;

  for(int IE = 0; IE < NBE; IE++)
  {
    fscanf(IRD, "%lf", &EBT[IE]);
    for(int IW = 0; IW < NBW; IW++)
    { 
      if((IW+1)%5==0){
      fscanf(IRD, "%lf%*[^\n]", &XS[IE][IW]);
      getc(IRD);
      }
      else{fscanf(IRD, "%lf", &XS[IE][IW]);}
    }
    fscanf(IRD, "%lf%*[^\n]", &TXS[IE]);
    getc(IRD);
    if(INFO >= 2)
    {
      fprintf(IWR, "%9.2E", EBT[IE]);
      int contador_salt = 0;
      for(int IW = 0; IW < NBW; IW++)
      {
        fprintf(IWR, " %11.5E", XS[IE][IW]);
        contador_salt++;
        if(contador_salt == 5)
        {
          fprintf(IWR, "\n         ");
          contador_salt = 0;
        }
      }
      fprintf(IWR, "                                    %10.3E\n", TXS[IE]);
    }
    X[IE] = log(EBT[IE]);
  }

  //  ************  Initialisation of the calculation routines.

  for(int I = 0; I < NBW; I++)
  {
    WB[I] = WB0[I];
  }

  //  ****  Compute the scaled energy loss distribution and sampling
  //        parameters for the energies in the simulation grid.
  //
  //  ****  Interpolation in E.

  double ELL;
  for(int IW = 0; IW < NBW; IW++)
  {
    for(int IE = 0; IE < NBE; IE++)
    {
      Y[IE] = log(XS[IE][IW]);
    }
    SPLINE(X,Y,A,B,C,D,0.0,0.0,NBE);
    if(IRETRN != 0){ return;}
    for(int I = 0; I < NEGP; I++)
    {
      int J;
      ELL = CEGRID_.DLEMP[I];
      FINDI(X,ELL,NBE,J);
      if(ELL > X[0])
      {
        P0[M-1][I][IW] = exp(A[J-1]+ELL*(B[J-1]+ELL*(C[J-1]+ELL*D[J-1])));
      }
      else
      {
        double F1 = A[0]+X[0]*(B[0]+X[0]*(C[0]+X[0]*D[0]));
        double FP1 = B[0]+X[0]*(2.0*C[0]+X[0]*3.0*D[0]);
        P0[M-1][I][IW] = exp(F1+FP1*(ELL-X[0]));
      }
    }
  }
  for(int IE = 0; IE < NEGP; IE++)
  {
    for(int IW = 0; IW < NBW; IW++)
    {
      PDF[IW] = P0[M-1][IE][IW];
    }

    RLPAC(WB,PDF,PAC,NBW);
    for(int IW = 0; IW < NBW; IW++)
    {
      PDFB[M-1][IE][IW] = PDF[IW];
      PACB[M-1][IE][IW] = PAC[IW];
    }
    for(int IW = 0; IW < NBW-1; IW++)
    {
      DPDFB[M-1][IE][IW] = PDFB[M-1][IE][IW+1]-PDFB[M-1][IE][IW];
    }
    DPDFB[M-1][IE][NBW-1] = 0.0;
    //  ****  The cutoff scaled energy loss is slightly modified to ensure
    // that the sampling routine EBR covers the allowed energy loss interval.
    double XC;
    if(IE+1 < NEGP)
    {
      XC = WCRM/CEGRID_.ET[IE+1];
    }
    else
    {
      XC = WCRM/CEGRID_.ET[NEGP-1];
    }
    if(XC < 1.0)
    {
      PBCUT[M-1][IE] = RLMOM(WB,PDF,XC,NBW,-1);
      if(IRETRN != 0){ return;}
      WBCUT[M-1][IE] = XC;
    }
    else
    {
      PBCUT[M-1][IE] = RLMOM(WB,PDF,1.0,NBW,-1);
      if(IRETRN != 0){ return;}
      WBCUT[M-1][IE] = 1.0;
    }
  }

}

