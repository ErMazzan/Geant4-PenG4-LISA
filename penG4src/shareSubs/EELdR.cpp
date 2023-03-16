
//  *********************************************************************
//                       SUBROUTINE EELdR
//  *********************************************************************
void PenInterface::EELdR(int M, FILE* IRD, FILE* IWR, int &INFO)
{
  //     This subroutine reads elastic cross sections for electrons and
  //  positrons in material M from the elastic scattering database. It also
  //  initializes the algorithm for simulation of elastic collisions of
  //  electrons and positrons.

  using namespace PENELOPE_mod;
  using namespace PENERROR_mod;

  using namespace COMPOS;
///////////  using namespace CEGRID;
  using namespace CEINTF;
  using namespace CEIMFP;
  using namespace CPIMFP;
  using namespace CEEL01;
  using namespace CDCSEP;
  using namespace CRITA;
  using namespace CEELDB;
  using namespace CPELDB;
  using namespace CELSEP;
  
  const int NE=96;
  const int NP_P=128;

  char CTEXT[51];

  double EGRD[16] = {1.0, 1.25, 1.50, 1.75, 2.00, 2.50, 3.00, 3.50, 4.00, 4.50, 5.00, 6.00, 7.00, 8.00, 9.00, 1.00E1};

  //  ****  Energy mesh points (in eV).
  
  int IE = 0;
  int IGRID = 10;
  double FGRID = 10.0;
  bool brkIt = false;
  double EV;
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

  int I = 0;
  TH[I] = 0.0;
  THR[I] = TH[I]*PI/180.0;
  XMU[I] = (1.0-cos(THR[I]))/2.0;
  XMUL[I] = log(1.0E-35);
  I = 1;
  TH[I] = 1.0E-4;
  THR[I] = TH[I]*PI/180.0;
  XMU[I] = (1.0-cos(THR[I]))/2.0;

  XMUL[I] = XMU[I];
  if(XMUL[I] < 1.0E-35){ XMUL[I] = 1.0E-35;}
  XMUL[I]=log(XMUL[I]);

  brkIt = false;
  while(!brkIt)
  {
    brkIt = true;
    I = I+1;
    if(TH[I-1] < 0.9999E-3)
    {
      TH[I] = TH[I-1]+2.5E-5;
    }
    else if(TH[I-1] < 0.9999E-2)
    {
      TH[I] = TH[I-1]+2.5E-4;
    }
    else if(TH[I-1] < 0.9999E-1)
    {
      TH[I] = TH[I-1]+2.5E-3;
    }
    else if(TH[I-1] < 0.9999)
    {
      TH[I] = TH[I-1]+2.5E-2;
    }
    else if(TH[I-1] < 0.9999E+1)
    {
      TH[I] = TH[I-1]+1.0E-1;
    }
    else if(TH[I-1] < 2.4999E+1)
    {
      TH[I] = TH[I-1]+2.5E-1;
    }
    else
    {
      TH[I] = TH[I-1]+5.0E-1;
    }
    THR[I] = TH[I]*PI/180.0;
    XMU[I] = (1.0-cos(THR[I]))/2.0;
    if(XMU[I] < 1.0E-35){ XMU[I] = 1.0E-35;}
    if( XMU[I] > 1.0E-35){ XMUL[I] = log(XMU[I]);}
    else{ XMUL[I] = log(1.0E-35);}
      
    if((I+1) < NA){ brkIt = false; continue;}
  }

  //  ****  Read elastic DCS tables.

  if(INFO > 3){ fprintf(IWR, "\n *** Electron elastic differential cross sections\n");}
  int IELEC = -1;
  fscanf(IRD, "%50[^\n]%*[^\n]", CTEXT);
  getc(IRD);
  for(IE = 0; IE < NE; IE++)
  {
    double ETSIE;
    fscanf(IRD, "%d %lf %lf %lf %lf%*[^\n]", &IELEC, &ETSIE, &ECS[IE], &ETCS1[IE], &ETCS2[IE]);
    getc(IRD);
    if(INFO > 3){ fprintf(IWR, "%3d%10.3E%12.5E%12.5E%12.5E\n", IELEC, ETS[IE], ECS[IE], ETCS1[IE], ETCS2[IE]);}
    fflush(IWR);
    if(IELEC != -1 || fabs(ETSIE-ETS[IE]) > 0.1)
    {
      fprintf(IWR, "\n Error reading electron elastic DCS data.\n");
      ErrorFunction(1348); return;
    }
    for(int K = 0; K < NA; K++)
    {
      fscanf(IRD, "%lf", &EDCS[IE][K]);
    }
    fscanf(IRD,"%*[^\n]");
    getc(IRD);

    if(INFO > 3)
    {
      for(int K = 0; K < NA; K++)
      {
        if((K+1)%10==0 || K==NA-1){fprintf(IWR, " %11.5E\n", EDCS[IE][K]);}
        else{fprintf(IWR, " %11.5E", EDCS[IE][K]);}
      }
    }
  //  ****  Consistency test.
    for(int K = 0; K < NA; K++)
    {
      DCSI[K] = EDCS[IE][K];
    }
    double ECS0 = 4.0*PI*RMOMX(XMU,DCSI,0.0,1.0,NA,0);
    if(IRETRN != 0){ return;}
    double ECS1 = 4.0*PI*RMOMX(XMU,DCSI,0.0,1.0,NA,1);
    double ECS2 = 4.0*PI*RMOMX(XMU,DCSI,0.0,1.0,NA,2);
    double TCS1 = 2.0*ECS1;
    double TCS2 = 6.0*(ECS1-ECS2);
    double TS0 = (ECS0-ECS[IE])/ECS[IE];
    double TS1 = (TCS1-ETCS1[IE])/ETCS1[IE];
    double TS2 = (TCS2-ETCS2[IE])/ETCS2[IE];
    double TSTE = fabs(TS0);
    if(TSTE < fabs(TS1)){ TSTE = fabs(TS1);}
    if(TSTE < fabs(TS2)){ TSTE = fabs(TS2);}
    if(TSTE > 1.0E-4)
    {
      fprintf(IWR, "\n E=%12.5E\n", ETS[IE]);
      fprintf(IWR, " Electron cross section data are corrupt.\n");
      ErrorFunction(1349); return;
    }
  }

  if(INFO > 3){ fprintf(IWR, "\n *** Positron elastic differential cross sections\n");}
  IELEC = +1;
  fscanf(IRD, "%50[^\n]%*[^\n]", CTEXT);
  getc(IRD);
  for(IE = 0; IE < NE; IE++)
  {
    double ETSIE;
    fscanf(IRD, "%d %lf %lf %lf %lf%*[^\n]", &IELEC, &ETSIE, &PCS[IE], &PTCS1[IE], &PTCS2[IE]);
    getc(IRD);
    if(INFO > 3){ fprintf(IWR, "%3d%10.3E%12.5E%12.5E%12.5E\n",IELEC,ETS[IE],PCS[IE],PTCS1[IE],PTCS2[IE]);}
    if(IELEC != +1 || fabs(ETSIE-ETS[IE]) > 0.1)
    {
      fprintf(IWR, "\n Error reading positron elastic DCS data.\n");
      ErrorFunction(1350); return;
    }
    for(int K = 0; K < NA; K++)
    {
      fscanf(IRD, "%lf", &PDCS[IE][K]);
    }
    fscanf(IRD,"%*[^\n]");
    getc(IRD);
      
    if(INFO > 3)
    {
      for(int K = 0; K < NA; K++)
      {
        if((K+1)%10==0 || K==NA-1){fprintf(IWR, " %11.5E\n", PDCS[IE][K]);}
        else{fprintf(IWR, " %11.5E", PDCS[IE][K]);}
      }
    }
      //  ****  Consistency test.
    for(int K = 0; K < NA; K++)
    {
      DCSI[K] = PDCS[IE][K];
    }
    double ECS0 = 4.0*PI*RMOMX(XMU,DCSI,0.0,1.0,NA,0);
    if(IRETRN != 0){ return;}
    double ECS1 = 4.0*PI*RMOMX(XMU,DCSI,0.0,1.0,NA,1);
    double ECS2 = 4.0*PI*RMOMX(XMU,DCSI,0.0,1.0,NA,2);
    double TCS1 = 2.0*ECS1;
    double TCS2 = 6.0*(ECS1-ECS2);
    double TS0 = (ECS0-PCS[IE])/PCS[IE];
    double TS1 = (TCS1-PTCS1[IE])/PTCS1[IE];
    double TS2 = (TCS2-PTCS2[IE])/PTCS2[IE];

    double TSTE = fabs(TS0);
    if(TSTE < fabs(TS1)){ TSTE = fabs(TS1);}
    if(TSTE < fabs(TS2)){ TSTE = fabs(TS2);}

    if(TSTE > 1.0E-4)
    {
      fprintf(IWR, "\n E=%11.5E\n", ETS[IE]);
      fprintf(IWR, " Positron cross section data are corrupt.\n");
      ErrorFunction(1351); return;
    }
  }
  int NPP = NP_P;
  int NU = NPP/4;  //ttttt

  //  ************  Electrons.

  int IEME = 0;
  for(int _KE = 0; _KE < NEGP; _KE++)
  {
    double ERRM;
    if(CEGRID_.ET[_KE] > 0.999999E8){ break;}
    DCSEL0(CEGRID_.ET[_KE],-1);
    if(IRETRN != 0){ return;}
    int IER;
    RITAI0(DCSEL,0.0,1.0,NPP,NU,ERRM,0,IER);
    if(IER != 0)
    {
      printf("IER=%d\n",IER);
      printf("EEldR. RITA initialisation error.\n");
      ErrorFunction(1327); return;
    }
    for(I = 0; I < NP_P; I++)
    {
      XSE[I][_KE][M-1] = XTI[I];
      PSE[I][_KE][M-1] = PACI[I];
      ASE[I][_KE][M-1] = AI[I];
      BSE[I][_KE][M-1] = BI[I];
      ITLE[I][_KE][M-1] = ITLI[I];
      ITUE[I][_KE][M-1] = ITUI[I];
    }
    double XM0A, XM1, XM2;
    RITAM(0.0,1.0,XM0A,XM1,XM2);
    double ECS0 = CSI;
    double ECS1 = CSI*XM1/XM0A;
    double ECS2 = CSI*XM2/XM0A;
    XE0[_KE] = ECS0;
    XE1[_KE] = 2.0*ECS1;
    XE2[_KE] = 6.0*(ECS1-ECS2);

    double FPEL = 1.0/(XE0[_KE]*VMOL[M-1]);
    double FPT1 = 1.0/(XE1[_KE]*VMOL[M-1]);
    double FPST = CEGRID_.ET[_KE]/(CSTPE[M-1][_KE]+RSTPE[M-1][_KE]);

    double XS0H = C1[M-1]*FPT1;
    if(XS0H > C2[M-1]*FPST){ XS0H = C2[M-1]*FPST;}
    if(XS0H < FPEL){ XS0H = FPEL;}

    XS0H = 1.0/(VMOL[M-1]*XS0H);

    double RNDC = 1.0-XS0H/XE0[_KE];
    if(RNDC < 1.0E-10){ RNDC = 1.0E-10;}
      
    if(RNDC < 1.0E-6){ RNDC = 0.0;}
    RNDCEd[M-1][_KE] = RNDC;
      
    double RU = RNDC;
    I = 1;
    int J = NP_P;
    bool brkIt2 = false;
    while(!brkIt2)
    {
      brkIt2 = true;
      int K = (I+J)/2;   
      if(RU > PSE[K-1][_KE][M-1])
      {
        I = K;
      }
      else
      {
        J = K;
      }
      if(J-I > 1){ brkIt2 = false; continue;}
    }

    double RR = RU-PSE[I-1][_KE][M-1];
    double DPRO = PSE[I+1-1][_KE][M-1]-PSE[I-1][_KE][M-1];
    double RMUC;
    if(DPRO < 1.0E-10)
    {
      RMUC = XSE[I-1][_KE][M-1];
    }
    else
    {
      double CI = (1.0+ASE[I-1][_KE][M-1]+BSE[I-1][_KE][M-1])*DPRO;
      RMUC = XSE[I-1][_KE][M-1]+(CI*RR/((DPRO*DPRO)+(DPRO*ASE[I-1][_KE][M-1]+BSE[I-1][_KE][M-1]*RR)*RR))*(XSE[I+1-1][_KE][M-1]-XSE[I-1][_KE][M-1]);
    }

      //  ****  Moments of the PDF on the restricted interval (0,RMUC).
      //        Total and transport cross sections for soft interactions.
      
    double XM0; 
    RITAM(0.0,RMUC,XM0,XM1,XM2);
    ECS1 = CSI*XM1/XM0A;
    ECS2 = CSI*XM2/XM0A;
    double TCS1 = 2.0*ECS1;
    double TCS2 = 6.0*(ECS1-ECS2);
    SEHEL[M-1][_KE] = XS0H*VMOL[M-1];
    T1E0[_KE] = TCS1;
    T1E[M-1][_KE] = T1EI[_KE]+TCS1*VMOL[M-1];
    T2E0[_KE] = TCS2;
    T2E[M-1][_KE] = T2EI[_KE]+TCS2*VMOL[M-1];
    IEME = _KE+1;
  }

  EELMAX[M-1] = CEGRID_.ET[IEME-1]-1.0;
  if(EELMAX[M-1] > 0.999999E8){ EELMAX[M-1] = 0.999999E8;}
  
  //  ****  Print electron elastic scattering tables.

  if(INFO >= 3){ fprintf(IWR, "\n PENELOPE >>>  Elastic scattering of electrons (ELSEPA database)\n");}
  if(INFO >= 3){ fprintf(IWR, "\n   E (eV)      MFP (mtu)   TMFP1 (mtu)  MFPh (mtu)\n --------------------------------------------------\n");}
  for(I = 0; I < IEME; I++)
  {
    double FP0 = RHO[M-1]/(XE0[I]*VMOL[M-1]);
    double FP1 = RHO[M-1]/(XE1[I]*VMOL[M-1]);
    double HMFP = RHO[M-1]/SEHEL[M-1][I];
    if(INFO >= 3){ fprintf(IWR, "%12.5E %12.5E %12.5E %12.5E\n", CEGRID_.ET[I], FP0, FP1, HMFP);}
    SEHEL[M-1][I] = log(SEHEL[M-1][I]);
      //  ****  Soft scattering events are switched off when T1E is too small.
    if(T1E[M-1][I] > 1.0E-6*XE1[I]*VMOL[M-1])
    {
      if(T1E[M-1][I] < 1.0E-35){ T1E[M-1][I] = log(1.0E-35);}
      else{ T1E[M-1][I] = log(T1E[M-1][I]);}

      if(T2E[M-1][I] < 1.0E-35){ T2E[M-1][I] = log(1.0E-35);}
      else{ T2E[M-1][I] = log(T2E[M-1][I]);}
    }
    else
    {
      T1E[M-1][I] = -100.0;
      T2E[M-1][I] = -100.0;
    }
  }

  //  ************  Positrons.

  int IEMP = 0;
  for(int _KE = 0; _KE < NEGP; _KE++)
  {
    double ERRM;
    if(CEGRID_.ET[_KE] > 0.999999E8){break;}
    DCSEL0(CEGRID_.ET[_KE],+1);
    if(IRETRN != 0){ return;}
    int IER;
    RITAI0(DCSEL,0.0,1.0,NPP,NU,ERRM,0,IER);
    if(IER != 0)
    {
      printf("IER=%d\n",IER);
      printf("EEldR. RITA initialisation error.\n");
      ErrorFunction(1327); return;
    }
    for(I = 0; I < NP_P; I++)
    {
      XSP[I][_KE][M-1] = XTI[I];
      PSP[I][_KE][M-1] = PACI[I];
      ASP[I][_KE][M-1] = AI[I];
      BSP[I][_KE][M-1] = BI[I];
      ITLP[I][_KE][M-1] = ITLI[I];
      ITUP[I][_KE][M-1] = ITUI[I];
    }
    double XM0A,XM1,XM2;
    RITAM(0.0,1.0,XM0A,XM1,XM2);
    double ECS0 = CSI;
    double ECS1 = CSI*XM1/XM0A;
    double ECS2 = CSI*XM2/XM0A;
    XP0[_KE] = ECS0;
    XP1[_KE] = 2.0*ECS1;
    XP2[_KE] = 6.0*(ECS1-ECS2);

    double FPEL = 1.0/(XP0[_KE]*VMOL[M-1]);
    double FPT1 = 1.0/(XP1[_KE]*VMOL[M-1]);
    double FPST = CEGRID_.ET[_KE]/(CSTPP[M-1][_KE]+RSTPP[M-1][_KE]);

    double XS0H = C1[M-1]*FPT1;
    if(XS0H > C2[M-1]*FPST){ XS0H = C2[M-1]*FPST;}
    if(XS0H < FPEL){ XS0H = FPEL;}
    XS0H=1.0/(VMOL[M-1]*XS0H);

    double RNDC = 1.0-XS0H/XP0[_KE];
    if(RNDC < 1.0E-10){ RNDC = 1.0E-10;}

    if(RNDC < 1.0E-6){ RNDC = 0.0;}
    RNDCPd[M-1][_KE] = RNDC;

    double RU = RNDC;
    I = 1;
    int J = NP_P;
    bool brkIt2 = false;
    while(!brkIt2)
    {
      brkIt2 = true;
      int K = (I+J)/2;
      if(RU > PSP[K-1][_KE][M-1])
      {
        I = K;
      }
      else
      {
        J = K;
      }
      if(J-I > 1){ brkIt2 = false; continue;}
    }

    double RMUC, CI;
    double RR = RU-PSP[I-1][_KE][M-1];
    double DPRO = PSP[I+1-1][_KE][M-1]-PSP[I-1][_KE][M-1];
    if(DPRO < 1.0E-10)
    {
      RMUC = XSP[I-1][_KE][M-1];
    }
    else
    {
      CI = (1.0+ASP[I-1][_KE][M-1]+BSP[I-1][_KE][M-1])*DPRO;
      RMUC = XSP[I-1][_KE][M-1]+(CI*RR/((DPRO*DPRO)+(DPRO*ASP[I-1][_KE][M-1]+BSP[I-1][_KE][M-1]*RR)*RR))*(XSP[I+1-1][_KE][M-1]-XSP[I-1][_KE][M-1]);
    }

      //  ****  Moments of the PDF on the restricted interval (0,RMUC).
      //        Total and transport cross sections for soft interactions.

    double XM0;
    RITAM(0.0,RMUC,XM0,XM1,XM2);
    ECS1 = CSI*XM1/XM0A;
    ECS2 = CSI*XM2/XM0A;
    double TCS1 = 2.0*ECS1;
    double TCS2 = 6.0*(ECS1-ECS2);
    SPHEL[M-1][_KE] = XS0H*VMOL[M-1];
    T1P0[_KE] = TCS1;
    T1P[M-1][_KE] = T1PI[_KE]+TCS1*VMOL[M-1];
    T2P0[_KE] = TCS2;
    T2P[M-1][_KE] = T2PI[_KE]+TCS2*VMOL[M-1];
    IEMP = _KE+1;
  }

  PELMAX[M-1] = CEGRID_.ET[IEMP-1]-1.0;
  if(PELMAX[M-1] > 0.999999E8){ PELMAX[M-1] = 0.999999E8;}

  //  ****  Print positron elastic scattering tables.

  if(INFO >= 3){ fprintf(IWR, "\n PENELOPE >>>  Elastic scattering of positrons (ELSEPA database)\n");}
  if(INFO >= 3){ fprintf(IWR, "\n   E (eV)      MFP (mtu)   TMFP1 (mtu)  MFPh (mtu)\n --------------------------------------------------\n");}
  for(I = 0; I < IEMP; I++)
  {
    double FP0 = RHO[M-1]/(XP0[I]*VMOL[M-1]);
    double FP1 = RHO[M-1]/(XP1[I]*VMOL[M-1]);
    double HMFP = RHO[M-1]/SPHEL[M-1][I];
    if(INFO >= 3){ fprintf(IWR, "%12.5E %12.5E %12.5E %12.5E\n", CEGRID_.ET[I], FP0, FP1, HMFP);}
    SPHEL[M-1][I] = log(SPHEL[M-1][I]);
      //  ****  Soft scattering events are switched off when T1P is too small.
    if(T1P[M-1][I] > 1.0E-6*XP1[I]*VMOL[M-1])
    {
      if(T1P[M-1][I] < 1.0E-35){ T1P[M-1][I] = log(1.0E-35);}
      else{ T1P[M-1][I] = log(T1P[M-1][I]);}

      if(T2P[M-1][I] < 1.0E-35){ T2P[M-1][I] = log(1.0E-35);}
      else{ T2P[M-1][I] = log(T2P[M-1][I]);}
    }
    else
    {
      T1P[M-1][I] = -100.0;
      T2P[M-1][I] = -100.0;
    }
  } 
}

