
//  *********************************************************************
//                       SUBROUTINE GRAaW
//  *********************************************************************
void PenInterface::GRAaW(int &M, FILE* IWR)
{
  //  This subroutine generates tables of molecular form factors and cross
  //  sections for Rayleigh scattering of photons in material M and writes
  //  them on the material data file. Data are read from the files
  //  'pdgraZZ.p08'.

  using namespace PENELOPE_mod;
  using namespace PENERROR_mod;

  using namespace COMPOS;

  char FILEN[13];
  char LDIG[10] = {'0','1','2','3','4','5','6','7','8','9'};
  char LDIG1[2], LDIG2[2];

  double Q[NQ], FF[NQ], FF2[NQ], ER[NEM], XSR[NEM];
  double EI[NEM], XS[NEM];

  //  ****  Momentum-transfer grid points and atomic form factor.
  for(int I = 0; I < NQ; I++)
  {
    Q[I] = 0.0;
    FF2[I] = 0.0;
  }
  //  ****  Energy grid points and Rayleigh x sections.
  double FACT1 = pow(10.0,(1.0/250.0));
  double FACT2 = pow(10.0,(1.0/25.0));
  double E = 10.0/FACT1;
  int NE=0;
  for(int I = 0; I < 2*NEM; I++)
  {
    if(E < 1.5848E5)
    {
      E = E*FACT1;
    }
    else
    {
      E = E*FACT2;
    }
    if(E > 49.0)
    {
      NE = NE+1;
      ER[NE-1] = E;
      XSR[NE-1] = 0.0;
      if(E > 1.2E9){ break;}
    }
  }

  int IZZ, NLD, NLD1, NLD2, IZZZ, NQI, NEI;
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
    strcpy(FILEN, "pdaff");
    strcat(FILEN, LDIG2);
    strcat(FILEN, LDIG1);
    strcat(FILEN, ".p08");
    strcat(straux, FILEN);
    FILE* pdaff = fopen(straux, "r");
    fscanf(pdaff, "%d %d%*[^\n]", &IZZZ, &NQI);
    getc(pdaff);
    
    if(IZZZ != IZZ){ ErrorFunction(1329); return;}
    if(NQI != NQ){ ErrorFunction(1330); return;}
    for(int I = 0; I < NQ; I++)
    {  
      fscanf(pdaff, "%lf %lf%*[^\n]", &Q[I], &FF[I]);
      getc(pdaff);
    }
    fclose(pdaff);

    strcpy(straux,"./pdfiles/");    
    strcpy(FILEN, "pdgra");
    strcat(FILEN, LDIG2);
    strcat(FILEN, LDIG1);
    strcat(FILEN, ".p08");
    strcat(straux, FILEN);
    FILE* pdgra = fopen(straux, "r");

    fscanf(pdgra, "%d %d%*[^\n]", &IZZZ, &NEI);
    getc(pdgra);
    if(IZZZ != IZZ){ ErrorFunction(1331); return;}
    for(int I = 0; I < NEI; I++)
    {
      double FA1, FA2;
      fscanf(pdgra, "%lf %lf %lf %lf%*[^\n]", &EI[I], &FA1, &FA2, &XS[I]);
      getc(pdgra);
      EI[I] = log(EI[I]);
      XS[I] = log(XS[I]);
    }
    fclose(pdgra);

    for(int I = 0; I < NQ; I++)
    {
      FF2[I] = FF2[I]+STF[M][IEL]*(FF[I]*FF[I]);
    }
    for(int I = 0; I < NE; I++)
    {
      int J;
      double EE = log(ER[I]);
      FINDI(EI,EE,NEI,J);
      double XSE = exp(XS[J-1]+(XS[J+1-1]-XS[J-1])*(EE-EI[J-1])/(EI[J+1-1]-EI[J-1]));
      XSR[I] = XSR[I]+STF[M][IEL]*XSE;
    }
  }

  fprintf(IWR, " *** Rayleigh scattering.  NQ = %3d,  NE = %4d\n", NQ, NE);
  for(int I = 0; I < NQ; I++)
  {
    fprintf(IWR, "%9.2E%12.5E\n", Q[I], sqrt(FF2[I]));
  }
  for(int I = 0; I < NE; I++)
  {
    fprintf(IWR, "%12.5E%12.5E\n", ER[I], XSR[I]);
  }     
}

