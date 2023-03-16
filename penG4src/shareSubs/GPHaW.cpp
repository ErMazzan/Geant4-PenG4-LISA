
//  *********************************************************************
//                       SUBROUTINE GPHaW
//  *********************************************************************
void PenInterface::GPHaW(int &M, FILE* IWR)
{
  //  This subroutine generates the table of photoelectric cross sections
  //  for photons in material M and writes it on the material data file.
  //  Data are read from the files 'pdgphZZ.p18'.

  using namespace PENELOPE_mod;
  using namespace PENERROR_mod;

  using namespace COMPOS;
  
  char FILEN[12];
  char LDIG[10] = {'0','1','2','3','4','5','6','7','8','9'};
  char LDIG1[2], LDIG2[2];

  const int NPHM = 2000;
  double XS[17], E0[NPHM], XS0[NPHM][17];
  int ISH[17];

  int IZZ, NLD, NLD1, NLD2, IZZZ, NGP, NSHR;
  for(int IEL = 0; IEL < NELEM[M]; IEL++)
  {
    IZZ = IZ[M][IEL];
    NLD = IZZ;
    NLD1 = NLD-10*(NLD/10);
    NLD2 = (NLD-NLD1)/10;
    LDIG1[0] = LDIG[NLD1+1-1];LDIG1[1] = '\0';
    LDIG2[0] = LDIG[NLD2+1-1];LDIG2[1] = '\0';
    char straux[30];
    strcpy(straux,"./pdfiles/");
    strcpy(FILEN, "pdgph");
    strcat(FILEN, LDIG2);
    strcat(FILEN, LDIG1);
    strcat(FILEN, ".p18");
    strcat(straux,FILEN);
    FILE* pdgph = fopen(straux, "r");
    fscanf(pdgph, "%*15c%2i %3i %4i%*[^\n]", &IZZZ, &NSHR, &NGP);
    getc(pdgph);
    
    if(NGP > NPHM){ ErrorFunction(1336); return;}
    if(IZZZ != IZZ){ ErrorFunction(1337); return;}
    if(NSHR > 16){ ErrorFunction(1338); return;}
    ISH[0] = 0;
    if(NSHR>0)
    {
      fscanf(pdgph, "%*25c");
      for(int IS = 1; IS < NSHR+1-1; IS++)
      {
        fscanf(pdgph, "%2d%*10c", &ISH[IS]);
      }
      fscanf(pdgph,"%2d%*[^\n]", &ISH[NSHR+1-1]);
      getc(pdgph);
    }
    else
    {
      fscanf(pdgph,"%*[^\n]");
      getc(pdgph);
    }       
    fscanf(pdgph,"%*[^\n]");
    getc(pdgph);

    int NPTAB=0;
    for(int IE = 0; IE < NGP; IE++)
    {
      double ER;
      if(feof(pdgph)){break;}
      fscanf(pdgph,"%lf",&ER);
      for(int IS = 0; IS < NSHR+1; IS++)
      {
        fscanf(pdgph, "%lf", &XS[IS]);
      }
      fscanf(pdgph,"%*[^\n]");
      getc(pdgph);
 
      if(ER > 49.9 && ER < 1.01E9)
      {
        NPTAB = NPTAB+1;
        E0[NPTAB-1] = ER;
        for(int IS = 0; IS < NSHR+1; IS++)
        {
          XS0[NPTAB-1][IS] = XS[IS]*1.0E-24;
        }
      }
    }
    fclose(pdgph);

    fprintf(IWR, " *** Photoelectric cross sections,  IZ =%3i,  NSHELL =%3i,  NDATA = %4i\n", IZZ, NSHR, NPTAB);

    fprintf(IWR,"%13c%2d",' ',0);
    for(int IS = 1; IS < NSHR+1; IS++)
    {
      fprintf(IWR,"%10c%2d",' ',ISH[IS]);
    }
    fprintf(IWR,"\n");
    for(int IE = 0; IE < NPTAB; IE++)
    {
      fprintf(IWR,"%12.5E", E0[IE]);
      for(int IS = 0; IS < NSHR+1; IS++)
      {
        fprintf(IWR, "%12.5E", XS0[IE][IS]);
      }
      fprintf(IWR,"\n");
    }
  }
}

