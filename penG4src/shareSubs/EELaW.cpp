
//  *********************************************************************
//                       SUBROUTINE EELaW
//  *********************************************************************
void PenInterface::EELaW(int &M, FILE* IWR)
{
  //  This subroutine generates a table of integrated cross sections for
  //  elastic scattering of electrons and positrons in material M, and
  //  writes it on the material definition file. Data are read from the
  //  files 'pdeelZZ.p08'.

  using namespace PENELOPE_mod;
  using namespace PENERROR_mod;
  
  using namespace COMPOS;
  using namespace CEEL00;
    
  char FILEN[12];
  char LDIG[] = {'0','1','2','3','4','5','6','7','8','9'};
  char LDIG1[2], LDIG2[2];

  //  ****  Building the cross section table.

  for(int I = 0; I < NEGP; I++)
  {
    XE0[I] = 0.0;
    XE1[I] = 0.0;
    XE2[I] = 0.0;
    XP0[I] = 0.0;
    XP1[I] = 0.0;
    XP2[I] = 0.0;
  }

  int IZZ, IZZZ, NLD, NLD1, NLD2, NPTAB = 0;
  double WGHT;
  for(int IEL = 0; IEL < NELEM[M]; IEL++)
  {
    char straux[30];
    strcpy(straux,"./pdfiles/");
    IZZ = IZ[M][IEL];
    WGHT = STF[M][IEL];
    NLD = IZZ;
    NLD1 = NLD-10*(NLD/10);
    NLD2 = (NLD-NLD1)/10;
    LDIG1[0] = LDIG[NLD1+1-1];LDIG1[1]='\0';
    LDIG2[0] = LDIG[NLD2+1-1];LDIG2[1]='\0';
    strcpy(FILEN, "pdeel");
    strcat(FILEN, LDIG2);
    strcat(FILEN, LDIG1);
    strcat(FILEN, ".p08");
    strcat(straux, FILEN);
    FILE* pdeel_file = fopen(straux, "r");
    fscanf(pdeel_file, "%d%*[^\n]", &IZZZ);
    getc(pdeel_file);
    
    if(IZZZ != IZZ){ ErrorFunction(1303); return;}
    for(int I = 0; I < NEGP; I++)
    {
      double XE0P, XE1P, XE2P, XP0P, XP1P, XP2P;
      fscanf(pdeel_file, "%lf %lf %lf %lf %lf %lf %lf%*[^\n]", &EJT[I], &XE0P, &XE1P, &XE2P, &XP0P, &XP1P, &XP2P);
      getc(pdeel_file);
      if(feof(pdeel_file)){break;}
      
      XE0[I] = XE0[I]+WGHT*XE0P;
      XE1[I] = XE1[I]+WGHT*XE1P;
      XE2[I] = XE2[I]+WGHT*XE2P;
      XP0[I] = XP0[I]+WGHT*XP0P;
      XP1[I] = XP1[I]+WGHT*XP1P;
      XP2[I] = XP2[I]+WGHT*XP2P;
      NPTAB = I+1;
    }
    fclose(pdeel_file);
  }

  //  ****  Write final x-section table.

  fprintf(IWR, " *** Electron and positron elastic cross sections,  NDATA =%4d\n", NPTAB);
  for(int I = 0; I < NPTAB; I++)
  {
    fprintf(IWR, "%10.3E %11.5E %11.5E %11.5E %11.5E %11.5E %11.5E\n", EJT[I], XE0[I], XE1[I], XE2[I], XP0[I], XP1[I], XP2[I]);
  }     
}

