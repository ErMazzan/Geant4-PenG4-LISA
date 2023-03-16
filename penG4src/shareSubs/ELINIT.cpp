          
//  *********************************************************************
//                       SUBROUTINE ELINIT
//  *********************************************************************
void PenInterface::ELINIT(int* IZ, double* STF, int &NELEM)
{
  //  This subroutine reads atomic elastic cross sections for electrons and
  //  positrons from the database files and determines the molecular cross
  //  section as the incoherent sum of atomic cross sections.
  
  //  Input arguments:
  //    IZ (1:NELEM) .... atomic numbers of the elements in the compound.
  //    STF (1:NELEM) ... stoichiometric indexes.
  //    NELEM ........... number of different elements.
  //
  
  using namespace PENERROR_mod;
  using namespace CDCSEP;
  
  const int NE = 96;
  
  double EGRD[16] = {1.0,1.25,1.50,1.75,2.00,2.50,3.00,3.50,4.00,4.50,5.00,6.00,7.00,8.00,9.00,1.00E1};
  
  char LIT10[10] = {'0','1','2','3','4','5','6','7','8','9'};
  char LIT1[2], LIT2[2], LIT3[2];
  
  char FILE1[13];
  
  //  ****  Energy mesh points (in eV).
  
  int IE = 0;
  int IGRID = 10;
  double FGRID = 10.0;
  double EV;
  bool brkIt = false;
  while(!brkIt)
  {
    brkIt = true;
    IGRID = IGRID+1;
    EV = EGRD[IGRID-1]*FGRID;
    if(IGRID == 16)
    {
      IGRID = 1;
      FGRID = 10.0*FGRID;
    }
    IE = IE+1;
    ETS[IE-1] = EV;
    ETL[IE-1] = log(ETS[IE-1]);
    if(IE < NE){ brkIt = false; continue;}
  }
  
  //  ****  Angular grid (TH in deg, XMU=(1.0D0-COS(TH))/2).
  
  int I = 1;
  TH[I-1] = 0.0;
  THR[I-1] = TH[I-1]*PI/180.0;
  XMU[I-1] = (1.0-cos(THR[I-1]))/2.0;
  XMUL[I-1] = log(1.0E-35);
  I = 2;
  TH[I-1] = 1.0E-4;
  THR[I-1] = TH[I-1]*PI/180.0;
  XMU[I-1] = (1.0-cos(THR[I-1]))/2.0;
  
  if(XMU[I-1] < 1.0E-35){ XMUL[I-1] = log(1.0E-35);}
  else{ XMUL[I-1] = log(XMU[I-1]);}
  
  brkIt = false;
  while(!brkIt)
  {
    brkIt = true;
    I = I+1;
    if(TH[I-1-1] < 0.9999E-3)
    {
      TH[I-1] = TH[I-1-1]+2.5E-5;
    }
    else if(TH[I-1-1] < 0.9999E-2)
    {
      TH[I-1] = TH[I-1-1]+2.5E-4;
    }
    else if(TH[I-1-1] < 0.9999E-1)
    {
      TH[I-1] = TH[I-1-1]+2.5E-3;
    }
    else if(TH[I-1-1] < 0.9999)
    {
      TH[I-1] = TH[I-1-1]+2.5E-2;
    }
    else if(TH[I-1-1] < 0.9999E+1)
    {
      TH[I-1] = TH[I-1-1]+1.0E-1;
    }
    else if(TH[I-1-1] < 2.4999E+1)
    {
      TH[I-1]=TH[I-1-1]+2.5E-1;
    }
    else
    {
      TH[I-1] = TH[I-1-1]+5.0E-1;
    }
    THR[I-1] = TH[I-1]*PI/180.0;
      
    XMU[I-1] = (1.0-cos(THR[I-1]))/2.0;
    if(XMU[I-1] < 1.0E-35){ XMU[I-1] = 1.0E-35;}
      
    if(XMU[I-1] < 1.0E-35){ XMUL[I-1] = log(1.0E-35);}
    else{ XMUL[I-1] = log(XMU[I-1]);}
      
    if(I < NA){ brkIt = false; continue;}
  }
  
  for(IE = 0; IE < NE; IE++)
  {
    ECS[IE] = 0.0;
    ETCS1[IE] = 0.0;
    ETCS2[IE] = 0.0;
    PCS[IE] = 0.0;
    PTCS1[IE] = 0.0;
    PTCS2[IE] = 0.0;
    for(int IA = 0; IA < NA; IA++)
    {
      EDCS[IE][IA] = 0.0;
      PDCS[IE][IA] = 0.0;
    }
  }
  
  //  ****  Read atomic DCS tables and compute the molecular DCS as the
  //        incoherent sum of atomic DCSs.
  
  for(int IEL = 0; IEL < NELEM; IEL++)
  {
    int IZZ = IZ[IEL];
    double STFF = STF[IEL];
    int NS = IZ[IEL];
    if(NS > 999){ NS = 999;}
    int NS1 = NS-10*(NS/10);
    NS = (NS-NS1)/10;
    int NS2 = NS-10*(NS/10);
    NS = (NS-NS2)/10;
    int NS3 = NS-10*(NS/10);
    LIT1[0] = LIT10[NS1+1-1];LIT1[1] = '\0';
    LIT2[0] = LIT10[NS2+1-1];LIT2[1] = '\0';
    LIT3[0] = LIT10[NS3+1-1];LIT3[1] = '\0';
    
    char straux[30];
    strcpy(straux,"./pdfiles/");
    strcpy(FILE1, "eeldx");
    strcat(FILE1, LIT3);
    strcat(FILE1, LIT2);
    strcat(FILE1, LIT1);
    strcat(FILE1, ".p08");
    strcat(straux,FILE1);
    FILE* eeldx = fopen(straux, "r");
      
    int IELEC, IZR;
    double ENR, CSE, TCS1E, TCS2E, CSP, TCS1P, TCS2P;
    for(IE = 0; IE < NE; IE++)
    {
      fscanf(eeldx, "%d %d %lf %lf %lf %lf%*[^\n]", &IELEC, &IZR, &ENR, &CSE, &TCS1E, &TCS2E);
      getc(eeldx);
      if(IELEC != -1 || IZR != IZZ || fabs(ENR-ETS[IE]) > 1.0E-3){ ErrorFunction(1352); return;}
      for(int IA = 0; IA < NA; IA++)
      {
        fscanf(eeldx, "%lf", &DCSI[IA]);
      }
      fscanf(eeldx,"%*[^\n]");
      getc(eeldx);
      ECS[IE] = ECS[IE]+STFF*CSE;
      ETCS1[IE] = ETCS1[IE]+STFF*TCS1E;
      ETCS2[IE] = ETCS2[IE]+STFF*TCS2E;
      for(int IA = 0; IA < NA; IA++)
      {
        EDCS[IE][IA] = EDCS[IE][IA]+STFF*DCSI[IA];
      }
    }
    fclose(eeldx);
      
    strcpy(straux,"./pdfiles/");  
    strcpy(FILE1, "peldx");
    strcat(FILE1, LIT3);
    strcat(FILE1, LIT2);
    strcat(FILE1, LIT1);
    strcat(FILE1, ".p08");
    strcat(straux,FILE1);    
    FILE* peldx = fopen(straux, "r");
    for(IE = 0; IE < NE; IE++)
    {
      fscanf(peldx, "%d %d %lf %lf %lf %lf%*[^\n]", &IELEC, &IZR, &ENR, &CSP, &TCS1P, &TCS2P);
      getc(peldx);
      if(IELEC != +1 || IZR != IZZ || fabs(ENR-ETS[IE]) > 1.0E-3){ErrorFunction(1353); return;}
      for(int IA = 0; IA < NA; IA++)
      {
        fscanf(peldx, "%lf", &DCSI[IA]);
      }
      fscanf(peldx,"%*[^\n]");
      getc(peldx);
      PCS[IE] = PCS[IE]+STFF*CSP;
      PTCS1[IE] = PTCS1[IE]+STFF*TCS1P;
      PTCS2[IE] = PTCS2[IE]+STFF*TCS2P;
      for(int IA = 0; IA < NA; IA++)
      {
        PDCS[IE][IA] = PDCS[IE][IA]+STFF*DCSI[IA];
      }
    }
    fclose(peldx);
  }
}

