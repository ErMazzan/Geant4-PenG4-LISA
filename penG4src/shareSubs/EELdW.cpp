
//  *********************************************************************
//                       SUBROUTINE EELdW
//  *********************************************************************
void PenInterface::EELdW(int &M, FILE* IWR)
{
  //  This subroutine generates a table of differential cross sections for
  //  elastic scattering of electrons and positrons in material M, and
  //  writes it on the material definition file. Data are read from the
  //  ELSEPA elastic database files.

  using namespace PENELOPE_mod;
  using namespace PENERROR_mod;

  using namespace COMPOS;
  using namespace CDCSEP;
  
  const int NE = 96;

  int IZM[30];
  double STFM[30];

  for(int I = 0; I < NELEM[M]; I++)
  {
    IZM[I] = IZ[M][I];
    STFM[I] = STF[M][I];
  }
  ELINIT(IZM,STFM,NELEM[M]);
  if(IRETRN != 0){ return;}

  //  ************  Write final DCS tables.

  //  ****  Electrons.

  int IELEC = -1;
  fprintf(IWR, " *** Electron elastic differential cross sections\n");
  double ECS0, ECS1, ECS2, TCS1, TCS2, TS0, TS1, TS2, TSTE; 
  for(int IE = 0; IE < NE; IE++)
  {
    for(int K = 0; K < NA; K++)
    {
      DCSI[K] = EDCS[IE][K];
    }
    ECS0 = 4.0*PI*RMOMX(XMU,DCSI,0.0,1.0,NA,0);
    if(IRETRN != 0){ return;}
    ECS1 = 4.0*PI*RMOMX(XMU,DCSI,0.0,1.0,NA,1);
    ECS2 = 4.0*PI*RMOMX(XMU,DCSI,0.0,1.0,NA,2);
    TCS1 = 2.0*ECS1;
    TCS2 = 6.0*(ECS1-ECS2);
    fprintf(IWR, "%3d%10.3E%12.5E%12.5E%12.5E\n", IELEC, ETS[IE], ECS0, TCS1, TCS2);
    for(int K = 0; K < NA; K++)
    {
      fprintf(IWR, " %11.5E", EDCS[IE][K]);
      if((K+1)%10==0){fprintf(IWR, "\n");}
    }
    fprintf(IWR, "\n");
     //  ****  Consistency test.
    TS0 = (ECS0-ECS[IE])/ECS[IE];
    TS1 = (TCS1-ETCS1[IE])/ETCS1[IE];
    TS2 = (TCS2-ETCS2[IE])/ETCS2[IE];

    TSTE = fabs(TS0);
    if(TSTE < fabs(TS1)){ TSTE = fabs(TS1);}
    if(TSTE < fabs(TS2)){ TSTE = fabs(TS2);}
      
    if(TSTE > 1.0E-2)
    {
      fprintf(IWR, "( E=%12.5E)", ETS[IE]);
      fprintf(IWR, "(    %11.5E %11.5E %11.5E)", ECS0, TCS1, TCS2);
      fprintf(IWR, "(    %11.5E %11.5E %11.5E)", ECS[IE], ETCS1[IE], ETCS2[IE]);
      fprintf(IWR, " Electron cross section data are corrupt.\n");
      ErrorFunction(1346); return;
    }
  }

  //  ****  Positrons.
  
  IELEC = +1;
  fprintf(IWR, " *** Positron elastic differential cross sections\n");

  for(int IE = 0; IE < NE; IE++)
  {
    for(int K = 0; K < NA; K++)
    {
      DCSI[K] = PDCS[IE][K];
    }
    ECS0 = 4.0*PI*RMOMX(XMU,DCSI,0.0,1.0,NA,0);
    if(IRETRN != 0){ return;}
    ECS1 = 4.0*PI*RMOMX(XMU,DCSI,0.0,1.0,NA,1);
    ECS2 = 4.0*PI*RMOMX(XMU,DCSI,0.0,1.0,NA,2);
    TCS1 = 2.0*ECS1;
    TCS2 = 6.0*(ECS1-ECS2);
    fprintf(IWR, "%3d%10.3E%12.5E%12.5E%12.5E\n", IELEC, ETS[IE], ECS0, TCS1, TCS2);
    for(int K = 0; K < NA; K++)
    {
      fprintf(IWR, " %11.5E", PDCS[IE][K]);
      if((K+1)%10==0){fprintf(IWR, "\n");}
    }
    fprintf(IWR, "\n");
      //  ****  Consistency test.
    TS0 = (ECS0-PCS[IE])/PCS[IE];
    TS1 = (TCS1-PTCS1[IE])/PTCS1[IE];
    TS2 = (TCS2-PTCS2[IE])/PTCS2[IE];

    TSTE = fabs(TS0);
    if(TSTE < fabs(TS1)){ TSTE = fabs(TS1);}
    if(TSTE < fabs(TS2)){ TSTE = fabs(TS2);}
      
    if(TSTE > 1.0E-2)
    {
      fprintf(IWR, "( E= %11.5E)\n", ETS[IE]);
      fprintf(IWR, "    %11.5E %11.5E %11.5E\n", ECS0, TCS1, TCS2);
      fprintf(IWR, "    %11.5E %11.5E %11.5E\n", PCS[IE], PTCS1[IE], PTCS2[IE]);
      fprintf(IWR, " Positron cross section data are corrupt.\n");
      ErrorFunction(1347); return;
    }
  }
}

