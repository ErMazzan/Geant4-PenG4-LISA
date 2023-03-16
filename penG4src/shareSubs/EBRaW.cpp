
//  *********************************************************************
//                       SUBROUTINE EBRaW
//  *********************************************************************
void PenInterface::EBRaW(int &M, FILE* IWR)
{
  //  This subroutine generates a table of the scaled energy-loss cross
  //  section for bremsstrahlung emission by electrons in material M. Data
  //  are read from the files 'pdebrZZ.p08'.

  using namespace PENELOPE_mod;
  using namespace PENERROR_mod;

  using namespace COMPOS;
/////////////////  using namespace CEGRID;
  using namespace CEBR;
  using namespace CEBR01;
  using namespace CEBR02;
  
  const double TREV = 2.0*REV;

  char FILEN[12];
  char LDIG[10] = {'0','1','2','3','4','5','6','7','8','9'};
  char LDIG1[2], LDIG2[2];


  double A[NEGP], B[NEGP], C[NEGP], D[NEGP];
  double WB0[NBW] = {1.0E-12, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.925, 0.95, 0.97, 0.99, 0.995, 0.999, 0.9995, 0.9999, 0.99995, 0.99999, 1.0};
  double PDF[NBE];
  

  //  ****  'Equivalent' atomic number.

  double SUMZ2 = 0.0;
  double SUMS = 0.0;
  for(int IEL = 0; IEL < NELEM[M]; IEL++)
  {
    SUMZ2 = SUMZ2+STF[M][IEL]*(IZ[M][IEL]*IZ[M][IEL]);
    SUMS = SUMS+STF[M][IEL];
  }
  ZBR2[M] = SUMZ2/SUMS;

  //  ****  Building the scaled cross section table.

  for(int IE = 0; IE < NBE; IE++)
  {
    TXS[IE] = 0.0;
    for(int IW = 0; IW < NBW; IW++)
    {
      XS[IE][IW] = 0.0;
    }
  }
  int IZZ, NLD, NLD1, NLD2, IZZZ;
  double WGHT, TXSP;
  for(int IEL = 0; IEL < NELEM[M]; IEL++)
  {
    char straux[30];
    strcpy(straux,"./pdfiles/");
    IZZ = IZ[M][IEL];
    WGHT = STF[M][IEL]*IZZ*IZZ/ZBR2[M];
    NLD = IZZ;
    NLD1 = NLD-10*(NLD/10);
    NLD2 = (NLD-NLD1)/10;
    LDIG1[0] = LDIG[NLD1+1-1];LDIG1[1]='\0';
    LDIG2[0] = LDIG[NLD2+1-1];LDIG2[1]='\0';
    strcpy(FILEN, "pdebr");
    strcat(FILEN, LDIG2);
    strcat(FILEN, LDIG1);
    strcat(FILEN, ".p08");
    strcat(straux, FILEN);
    FILE* pdebr = fopen(straux, "r");
    fscanf(pdebr, "%d%*[^\n]", &IZZZ);
    getc(pdebr);
    if(IZZZ != IZZ){ ErrorFunction(1317); return;}
    for(int IE = 0; IE < NBE; IE++)
    {
      fscanf(pdebr, "%lf", &EBT[IE]);
      for(int IW = 0; IW < NBW; IW++)
      {
        if((IW+1)%5==0){
        fscanf(pdebr, "%lf%*[^\n]", &PDF[IW]);
        getc(pdebr);
      }
      else{fscanf(pdebr, "%lf", &PDF[IW]);}
      } 
      fscanf(pdebr, "%*36c %lf%*[^\n]", &TXSP);
      getc(pdebr);
      TXS[IE] = TXS[IE]+WGHT*TXSP;
      for(int IW = 0; IW < NBW; IW++)
      {
        XS[IE][IW] = XS[IE][IW]+WGHT*PDF[IW];
      }
    }
    fclose(pdebr);
  }
  
  //  ****  The energy loss spectrum is re-normalized to reproduce the
  //        total scaled cross section of Berger and Seltzer.

  double RSUM, FACT, FNORM, TST;
  for(int IE = 0; IE < NBE; IE++)
  {
    for(int IW = 0; IW < NBW; IW++)
    {
      X[IW] = WB0[IW];
      Y[IW] = XS[IE][IW];
    }
    RSUM = RLMOM(X,Y,1.0,NBW,0);
    FACT = (EBT[IE]+REV)*1.0E-27*137.03604/((ELRAD*ELRAD)*(EBT[IE]+TREV));
    FNORM = TXS[IE]/(RSUM*FACT);
    TST = 100.0*fabs(FNORM-1.0);
    if(TST > 1.0){ ErrorFunction(1318); return;}
    for(int IW = 0; IW < NBW; IW++)
    {
      XS[IE][IW] = XS[IE][IW]*FNORM;
    }
  }

  //  ****  Write output scaled x-section table.
  
  fprintf(IWR, " *** Electron scaled bremss x-section,  ZEQ =%12.5E,  NDATA =%4d\n", sqrt(ZBR2[M]), NBE);
  for(int IE = 0; IE < NBE; IE++)
  {
    fprintf(IWR,"%9.2E", EBT[IE]);
    for(int IW = 0; IW < NBW; IW++)
    {
      if(IW != 0 && IW%5==0){fprintf(IWR,"%*c",9,' ');}
      if((IW+1)%5==0){fprintf(IWR,"%12.5E\n", XS[IE][IW]);}else{fprintf(IWR,"%12.5E", XS[IE][IW]);}

    }
    fprintf(IWR,"%*c%10.3E\n",36,' ',TXS[IE]);
  }


  //  ************  Initialisation of the calculation routines.

  for(int I = 0; I < NBW; I++)
  {
    WB[I] = WB0[I];
  }

  //  ****  Compute the scaled energy loss distribution and sampling
  //        parameters for the energies in the simulation grid.

  //  ****  Interpolation in E.

  double F1, FP1, ELL;
  int J;
  for(int IE = 0; IE < NBE; IE++)
  {
    X[IE] = log(EBT[IE]);
  }
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
      ELL = CEGRID_.DLEMP[I];
      if(ELL > X[0])
      {
        FINDI(X,ELL,NBE,J);
        P0[M][I][IW] = exp(A[J-1]+ELL*(B[J-1]+ELL*(C[J-1]+ELL*D[J-1])));
      }
      else
      {
        F1 = A[0]+X[0]*(B[0]+X[0]*(C[0]+X[0]*D[0]));
        FP1 = B[0]+X[0]*(2.0*C[0]+X[0]*3.0*D[0]);
        P0[M][I][IW] = exp(F1+FP1*(ELL-X[0]));
      }
    }
  }
}

