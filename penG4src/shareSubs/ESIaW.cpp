
//  *********************************************************************
//                       SUBROUTINE ESIaW
//  *********************************************************************
void PenInterface::ESIaW(int &M, FILE* IWR)
{
  //  This subroutine generates tables of cross sections for inner-shell
  //  ionisation by electron impact for the elements in material M and
  //  writes them on the material data file.

  //  Data are read from the files 'pdesiZZ.p14'.

  using namespace PENELOPE_mod;
  using namespace PENERROR_mod;

  using namespace COMPOS;
 
  char FILEN[12];
  char LDIG[] = {'0','1','2','3','4','5','6','7','8','9'};
  char LDIG1[2], LDIG2[2];
  
  const int NES = 850;
  double E[NES], XESIR[NES][16];

  int IZZ, IZZZ, NLD, NLD1, NLD2, NSHR, NPTAB = 0;
  for(int IEL = 0; IEL < NELEM[M]; IEL++)
  {
    char straux[30];
    strcpy(straux,"./pdfiles/");
    IZZ = IZ[M][IEL];
    NLD = IZZ;
    NLD1 = NLD-10*(NLD/10);
    NLD2 = (NLD-NLD1)/10;
    LDIG1[0] = LDIG[NLD1+1-1];LDIG1[1]='\0';
    LDIG2[0] = LDIG[NLD2+1-1];LDIG2[1]='\0';
    strcpy(FILEN, "pdesi");
    strcat(FILEN, LDIG2);
    strcat(FILEN, LDIG1);
    strcat(FILEN, ".p14");
    strcat(straux, FILEN);
    FILE* pdfiles = fopen(straux, "r");
    fscanf(pdfiles, "%*16c%2d%*6c%2d%*[^\n]", &IZZZ, &NSHR);
    getc(pdfiles);
    if(IZZZ != IZZ){ ErrorFunction(1308); return;}
    if(NSHR > 16){ ErrorFunction(1309); return;}
    fscanf(pdfiles,"%*[^\n]");
    getc(pdfiles);
    fscanf(pdfiles,"%*[^\n]");
    getc(pdfiles);

    for(int IE = 0; IE < NES; IE++)
    {
      if(feof(pdfiles)){break;}
      fscanf(pdfiles, "%lf", &E[IE]);
      for(int IS = 0; IS < NSHR; IS++)
      {
    
      if(IS<NSHR-1){fscanf(pdfiles, "%lf", &XESIR[IE][IS]);}
      else{
      fscanf(pdfiles, "%lf%*[^\n]", &XESIR[IE][IS]);
      getc(pdfiles);
      }
      }

      NPTAB = IE+1;
      if(E[IE] > 0.999E9){ break;}
    }
    fclose(pdfiles);
    fprintf(IWR, " *** Electron ionisation cross sections,  IZ =%3d,  NSHELL =%3d,  NDATA =%4d\n", IZZ, NSHR, NPTAB);
    for(int IE = 0; IE < NPTAB; IE++)
    {
      fprintf(IWR, "%12.5E", E[IE]);
      for(int IS = 0; IS < NSHR; IS++)
      {
        if(IS<NSHR-1){fprintf(IWR, "%12.5E", XESIR[IE][IS]);}
        else{fprintf(IWR, "%12.5E\n", XESIR[IE][IS]);}
      }
    }
  }
}

