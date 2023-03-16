
//  *********************************************************************
//                       SUBROUTINE GPPaW
//  *********************************************************************
void PenInterface::GPPaW(double* EIT, double* XGP0, double* XGT0, int &NPTAB, int &M)
{
  //  This subroutine generates a table of electron-positron pair produc-
  //  tion cross sections for photons in material M. Data are read from the
  //  files 'pdgppZZ.p11'.

  using namespace PENELOPE_mod;
  using namespace PENERROR_mod;

  using namespace COMPOS;
  
  char FILEN[12];
  char LDIG[10] = {'0','1','2','3','4','5','6','7','8','9'};
  char LDIG1[2], LDIG2[2];

  const int NEGPP = 10000;
  //double EIT[NEGPP], XGP0[NEGPP], XGT0[NEGPP];

  //  ****  Building the cross section table.

  for(int I = 0; I < NEGPP; I++)
  {
    XGP0[I] = 0.0;
    XGT0[I] = 0.0;
  }

  int IZZ, NLD, NLD1, NLD2, IZZZ;
  double WGHT;
  for(int IEL = 0; IEL < NELEM[M]; IEL++)
  {
    char straux[30];
    strcpy(straux,"./pdfiles/");
    IZZ = IZ[M][IEL];
    WGHT = STF[M][IEL]*1.0E-24;
    NLD = IZZ;
    NLD1 = NLD-10*(NLD/10);
    NLD2 = (NLD-NLD1)/10;
    LDIG1[0] = LDIG[NLD1+1-1];LDIG1[1]='\0';
    LDIG2[0] = LDIG[NLD2+1-1];LDIG2[1]='\0';

    strcpy(FILEN, "pdgpp");
    strcat(FILEN, LDIG2);
    strcat(FILEN, LDIG1);
    strcat(FILEN, ".p11");
    strcat(straux, FILEN);
    FILE* pdgpp = fopen(straux, "r");

    fscanf(pdgpp, "%d%*[^\n]", &IZZZ);
    getc(pdgpp);
    
    if(IZZZ != IZZ){ ErrorFunction(1339); return;}
    for(int I = 0; I < NEGPP; I++)
    {
      double XG0P, XG0T;
      fscanf(pdgpp, "%lf %lf %lf%*[^\n]", &EIT[I], &XG0P, &XG0T);
      getc(pdgpp);
      XGP0[I] = XGP0[I]+WGHT*XG0P;
      XGT0[I] = XGT0[I]+WGHT*XG0T;
      NPTAB = I+1;
      if(EIT[I] > 0.999E9){ break;}
    }
    fclose(pdgpp);
  }
}

