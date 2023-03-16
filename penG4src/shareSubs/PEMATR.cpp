
//  *********************************************************************
//                       SUBROUTINE PEMATR
//  *********************************************************************
void PenInterface::PEMATR(unsigned int M, FILE *IRD, FILE *IWR, int INFO)
{

  //  This subroutine reads the definition file of material M (unit IRD)
  //  and initialises the simulation routines for this material. Informa-
  //  tion is written on unit IWR, the amount of written information is
  //  determined by the value of INFO.

  using namespace PENELOPE_mod;
  using namespace PENERROR_mod;

  using namespace COMPOS;
  using namespace CADATA;
  using namespace CECUTR;
//////////////  using namespace CEGRID;
  using namespace CRANGE;
  using namespace CEIN;
  using namespace CEINTF;
  using namespace CEIN00;
  using namespace CPIN00;
  using namespace CESI0;
  using namespace CPSI0;
  using namespace CEINAC;
  using namespace CESIAC;
  using namespace CESIN;
  using namespace CPINAC;
  using namespace CPSIAC;
  using namespace CPSIN;
  using namespace CGCO;
  using namespace CEIMFP;
  using namespace CLAS1E;
  using namespace CPIMFP;
  using namespace CLAS1P;
  using namespace CEEL00;
  using namespace CBRYLD;
  using namespace CGIMFP;
  using namespace CGPH01;
  using namespace CGPP01;
    
  char NAME[63], LNAME[63];
  const double FOURPI=4.0*PI;
  
  //  ****  Positron inelastic coll. and inner-shell ionisation tables.
  double STFI[NO], STFO[NO];
  int INOUT[NO];
  
  //  ****  Auxiliary arrays.
  double EIT[NEGP], EITL[NEGP], FL[NEGP], F1[NEGP], F2[NEGP], F3[NEGP], F4[NEGP], RADY[NEGP], RADN[NEGP];

  //  ****  Rescaling of calculated MFPs.
  //  When ISCALE is set equal to 1, the program rescales the calculated
  //  total cross sections and MFPs to reproduce the cross sections and
  //  stopping powers read from the input material data file. To avoid this
  //  rescaling, set ISCALE=0.

  int ISCALE=1;

  //  ************  Material characteristics.
  
  if(M > MAXMAT){ErrorFunction(1008);return;}

  strcpy(LNAME, " PENELOPE (v. 2018)  Material data file ...............");
  fscanf(IRD,"%55c%*[^\n]",NAME);
  //Append end of string chars
  NAME[55] = '\0';
  getc(IRD);
  if(strcmp(NAME, LNAME) != 0)
  {
    fprintf(IWR, "\n I/O error. Corrupt material data file (material number %3d).\n      The first line is: %55s\n      ... and should be: %55s\n",M ,NAME, LNAME);
    ErrorFunction(1009);return;
  }
  fprintf(IWR, "%55s\n", NAME);

  fscanf(IRD,"%*11c%62[^\n]%*[^\n]",LNAME);
  getc(IRD);
  fprintf(IWR, " Material: %62s\n", LNAME);
  
  fscanf(IRD, "%*15c%lf%*[^\n]", &RHO[M-1]);
  getc(IRD);
  
  fprintf(IWR, " Mass density =%15.8E g/cm**3\n", RHO[M-1]);
  
  fscanf(IRD, "%*37c%3d%*[^\n]", &NELEM[M-1]);
  getc(IRD);
  fprintf(IWR, " Number of elements in the molecule = %2d\n", NELEM[M-1]);fflush(IWR);
  if(NELEM[M-1] > 30){ fprintf(IWR, "To many elements in material number %3d\n", M); ErrorFunction(1010);return;}
  
  ZT[M-1] = 0.0;
  AT[M-1] = 0.0;

  for(int I = 0; I < NELEM[M-1]; I++)
  {
    fscanf(IRD, "%*18c%d,%*19c%lf%*[^\n]", &IZ[M-1][I], &STF[M-1][I]);
    getc(IRD);
    fprintf(IWR, "    Element: %s (Z=%2d), atoms/molecule =%15.8E\n", LASYMB[IZ[M-1][I]-1], IZ[M-1][I], STF[M-1][I]);fflush(IWR);
      
    ZT[M-1] = ZT[M-1]+STF[M-1][I]*IZ[M-1][I];
    int IZZ = IZ[M-1][I];
    AT[M-1] = AT[M-1]+ATW[IZZ-1]*STF[M-1][I];
  }
  VMOL[M-1] = AVOG*RHO[M-1]/AT[M-1];

  if(INFO >= 2){ fprintf(IWR, "\n Molecular density = %15.8E 1/cm**3\n", VMOL[M-1]);}
  
  OP2[M-1] = FOURPI*ZT[M-1]*VMOL[M-1]*(A0B*A0B*A0B)*(HREV*HREV);
  double OMEGA = sqrt(OP2[M-1]);

  DEN[M-1] = RHO[M-1];
  RDEN[M-1]=1.0/DEN[M-1];

  if(INFO >= 2)
  {
    fprintf(IWR, "\n *** Electron/positron inelastic scattering.\n");
    fprintf(IWR, " Plasma energy = %15.8E eV\n", OMEGA);
  }
  
  fscanf(IRD, "%*25c%lf%*[^\n]", &EXPOT[M-1]);
  getc(IRD);
  fprintf(IWR, " Mean excitation energy =%15.8E eV\n", EXPOT[M-1]);

  //  ****  E/P inelastic collisions.

  fscanf(IRD, "%*24c%3d%*[^\n]", &NOSC[M-1]);
  getc(IRD);
  if(INFO >= 2 || NOSC[M-1] > NO){ fprintf(IWR, " Number of oscillators =%4d\n", NOSC[M-1]);}
  if(NOSC[M-1] > NO){ ErrorFunction(1011);return;}
  if(INFO >= 2){ fprintf(IWR, "\n           Fi            Ui (eV)         Wi (eV)      KZ  KS\n ------------------------------------------------------------\n");}
  double EXPT = 0.0;
  for(int I = 0; I < NOSC[M-1]; I++)
  {
    fscanf(IRD , "%*4c%lf %lf %lf %d %d%*[^\n]", &F[M-1][I], &UI[M-1][I], &WRI[M-1][I], &KZ[M-1][I], &KS[M-1][I]);
    getc(IRD);
      
    if(UI[M-1][I] < 1.0E-3)
    {
      UI[M-1][I] = 0.0;
    }
    if(INFO >= 2){ fprintf(IWR, "%4d%16.8E%16.8E%16.8E%4d%4d\n", I+1, F[M-1][I], UI[M-1][I], WRI[M-1][I], KZ[M-1][I], KS[M-1][I]);}
    EXPT = EXPT+F[M-1][I]*log(WRI[M-1][I]);
  }
  EXPT = exp(EXPT/ZT[M-1]);

  if(fabs(EXPT-EXPOT[M-1]) > 1.0E-6*EXPOT[M-1])
  {
    fprintf(IWR, "EXPT      =%.5E\nEXPOT (M) =%.5E\nInconsistent oscillator data.\n", EXPT, EXPOT[M-1]);
    ErrorFunction(1012);return;
  }

  //  ****  Compton scattering.

  fscanf(IRD,"%*19c%3d%*[^\n]", &NOSCCO[M-1]);
  getc(IRD);
  if(INFO >= 2 || NOSCCO[M-1] > NOCO)
  {
    fprintf(IWR, "\n *** Compton scattering (Impulse Approximation).\n");
    fprintf(IWR, " Number of shells =%4d\n", NOSCCO[M-1]);
    if(NOSCCO[M-1] > NOCO){ ErrorFunction(1013);return;}
  }
  if(INFO >= 2){ fprintf(IWR, "\n           Fi            Ui (eV)         Ji(0)        KZ  KS\n ------------------------------------------------------------\n");}

  for(int I = 0; I < NOSCCO[M-1]; I++)
  {
    fscanf(IRD, "%*4c%lf %lf %lf %d %d%*[^\n]", &FCO[M-1][I], &UICO[M-1][I], &FJ0[M-1][I], &KZCO[M-1][I], &KSCO[M-1][I]);
    getc(IRD);
    if(INFO >= 2){ fprintf(IWR, "%4d %15.8E %15.8E %15.8E%4d%4d\n", I+1, FCO[M-1][I], UICO[M-1][I], FJ0[M-1][I], KZCO[M-1][I], KSCO[M-1][I]);}
    FJ0[M-1][I] = FJ0[M-1][I]*SL;
  }

  //  ************  Atomic relaxation data.

  for(int I = 0; I < NELEM[M-1]; I++)
  {
      RELAXR(IRD, IWR, INFO);
      if(IRETRN != 0){return;}
  }

  //  ****  Inner-shell ionisation by electron and positron impact.

  ESIaR(M,IRD,IWR,INFO);
  if(IRETRN != 0){return;}
  PSIaR(M,IRD,IWR,INFO);
  if(IRETRN != 0){return;}

  //  ****  Electron and positron interaction properties.

  //  ****  Scaled bremsstrahlung x-section.
  int IBREMS;
  double WCRM;
  if(WCR[M-1] >= 0.0)
  {
    if(WCR[M-1] < 10.0){WCR[M-1] = 10.0;}     
    WCRM = WCR[M-1];
    IBREMS = 1;
  }
  else
  {
    WCRM = 10.0;
    WCR[M-1] = 10.0;
    IBREMS = 0;
  }
  EBRaR(WCRM, M, IRD, IWR, INFO);
  if(IRETRN != 0){return;}
  //  ****  Bremsstrahlung angular distribution.
  BRaAR(M,IRD,IWR,INFO);
  if(IRETRN != 0){return;}

  //  ****  Stopping powers.
  int NDATA;
  fscanf(IRD, "%*58c%4d%*[^\n]", &NDATA);
  getc(IRD);
  
  if(INFO >= 2){ fprintf(IWR, "\n *** Stopping powers for electrons and positrons,  NDATA =%4d\n", NDATA);}
  if(NDATA > NEGP){ ErrorFunction(1014); return;}
  if(INFO >= 2){ fprintf(IWR, "\n  Energy     Scol,e-     Srad,e-     Scol,e+     Srad,e+\n   (eV)     (MeV/mtu)   (MeV/mtu)   (MeV/mtu)   (MeV/mtu)\n ----------------------------------------------------------\n");}

  for(int I = 0; I < NDATA; I++)
  {
    fscanf(IRD, "%lf %lf %lf %lf %lf%*[^\n]", &EIT[I], &F1[I], &F2[I], &F3[I], &F4[I]);
    getc(IRD);
    
    if(INFO >= 2){ fprintf(IWR, "%10.3E%12.5E%12.5E%12.5E%12.5E\n", EIT[I], F1[I], F2[I], F3[I], F4[I]);}
    EITL[I] = log(EIT[I]);
  }

  //  ************  Inelastic x-section for electrons.

  double WCCM = WCC[M-1];
  double XH0, XH1, XH2, XS0, XS1, XS2, XT1, XT2;
  double DELTA;
  double DIFMAX=0.0;
  double STPI;
  double STPC;
  
  for(int I = 0; I < NDATA; I++)
  {
    if(EIT[I] >= CEGRID_.EL && EIT[I] <= CEGRID_.EU)
    {
      STPI = F1[I]*RHO[M-1]*1.0E6;  // Collision stopping power.
      EINaT( EIT[I], WCCM, XH0, XH1, XH2, XS0, XS1, XS2, XT1, XT2, DELTA, M);
      STPC = (XS1+XH1)*VMOL[M-1];
      double AuxDouble = 100.0*fabs(STPC-STPI)/STPI;
      if(DIFMAX < AuxDouble){DIFMAX = AuxDouble;}
    }
  }

  int ICOR;
  if(DIFMAX > 1.0E-3 && ISCALE == 1)
  {
    ICOR=1;
    for(int I = 0; I < NDATA; I++)
    {
      FL[I] = log(F1[I]*RHO[M-1]*1.0E6);
    }
    SPLINE(EITL, FL, A, B, C, D, 0.0, 0.0, NDATA);
    if(IRETRN != 0){return;}
  }
  else
  {
    ICOR=0;
  }

  //  **** Classify inner and outer shells.

  int NI = 0;
  int NS = 0;
  double EBIN;
  
  for(int IO = 0; IO < NO; IO++)
  {
    STFI[IO]=0.0;
    STFO[IO]=0.0;
  }

  for(int KO = 0; KO < NOSC[M-1]; KO++)
  {
    int IZZ = KZ[M-1][KO];
    int ISH = KS[M-1][KO];
    if(IZZ < 3 || ISH > 16 || UI[M-1][KO] < 50.0)
    {
      NI = NI+1;
      IEIN[M-1][NI-1] = KO+1;
      ISIE[KO] = -NI;
      INOUT[NI-1] = 0;
      for(int IEL = 0; IEL < NELEM[M-1]; IEL++)
      {
        if(IZZ == IZ[M-1][IEL]){ STFO[NI-1] = STF[M-1][IEL];}
      }
    }
    else if(ISH <= NSESI[IZZ-1])
    {
      EBIN = EB[IZZ-1][ISH-1];
      if(EBIN > ECUTR[M-1])
      {
        NS = NS+1;
        IESI[M-1][NS-1] = KO+1;
        ISIE[KO] = NS;
        for(int IEL = 0; IEL < NELEM[M-1]; IEL++)
        {
          if(IZZ == IZ[M-1][IEL]){ STFI[NS-1] = STF[M-1][IEL];}
      
        }
      }   
      else
      {
        NI = NI+1;
        IEIN[M-1][NI-1] = KO+1;
        ISIE[KO] = -NI;
        INOUT[NI-1] = 1;
        for(int IEL = 0; IEL < NELEM[M-1]; IEL++)
        {
          if(IZZ == IZ[M-1][IEL]){ STFO[NI-1] = STF[M-1][IEL];}
        }
      }
    }
    else
    {
      NI = NI+1;
      IEIN[M-1][NI-1] = KO+1;
      ISIE[KO] = -NI;
      INOUT[NI-1] = 0;
      for(int IEL = 0; IEL < NELEM[M-1]; IEL++)
      {
        if(IZZ == IZ[M-1][IEL]){ STFO[NI-1] = STF[M-1][IEL];}
      }
    }
  }
  NEIN[M-1] = NI;
  NESI[M-1] = NS;

  //  ****  Simulation arrays.
  double EE;
  double FACT;
  double EC;
  
  double PCSI;
  double STOT;
  double DFERMI;
  
  for(int I = 0; I < NEGP; I++)
  {
    EE = CEGRID_.ET[I];       
    EINaT(EE, WCCM, XH0, XH1, XH2, XS0, XS1, XS2, XT1, XT2, DELTA, M);
    STPC = (XS1+XH1)*VMOL[M-1];
    if(ICOR == 1)
    {
      int J;    
      EC = CEGRID_.DLEMP[I];
      FINDI(EITL, EC, NDATA, J);
      STPI = exp(A[J-1]+EC*(B[J-1]+EC*(C[J-1]+EC*D[J-1])));
      FACT = STPI/STPC;
    }
    else
    {
      FACT = 1.0;
    }
    DEL[M-1][I] = DELTA;
    CSTPE[M-1][I] = STPC*FACT;

    double XS1SI = 0.0;
    double XS2SI = 0.0;
    double XT1SI = 0.0;
    double XT2SI = 0.0;
    double STPSI = 0.0;
    double STRSI = 0.0;
    double XS1IN = 0.0;
    double XS2IN = 0.0;
    double XT1IN = 0.0;
    double XT2IN = 0.0;
    double STPIN = 0.0;
    double STRIN = 0.0;

    for(int IO = 0; IO < NEIN[M-1]; IO++)
    {
      XSEIN[I][IO] = 0.0;
    }
    for(int IO = 0; IO < NESI[M-1]; IO++)
    {
      XSESI[I][IO] = 0.0;
    }

    for(int KO = 0; KO < NOSC[M-1]; KO++)
    {
      int IO = ISIE[KO];
    //  ****  Inner shells
      if(IO > 0)
      {
        int IZZ = KZ[M-1][KO];
        int ISH = KS[M-1][KO];
        EBIN = EB[IZZ-1][ISH-1];
        if(EBIN < EE)
        {
          int INDC = IESIF[IZZ-1]-1;
          PCSI = exp(XESI[INDC+I][ISH-1])*STFI[IO-1];
          if(PCSI > 1.0E-35)
          {
            double H0, H1, H2, S0, S1, S2, R0, R1, R2;
            EINaT1(EE, UI[M-1][KO], WRI[M-1][KO], 0.0, WCCM, H0, H1, H2, S0, S1, S2, R0, R1, R2);
            STOT = F[M-1][KO]*(S0+H0);
            if(STOT > 1.0E-35)
            {
              DFERMI = (SES0[KO]+SEH0[KO])/STOT;
              PCSI = PCSI*DFERMI;
              XSESI[I][IO-1] = PCSI;
              double RNFO = PCSI/(SES0[KO]+SEH0[KO]);
              STPSI = STPSI+(SES1[KO]+SEH1[KO])*RNFO;
              STRSI = STRSI+(SES2[KO]+SEH2[KO])*RNFO;
            }
            else
            {
              XSESI[I][IO-1]=PCSI;
              STPSI = STPSI+PCSI*EBIN;
              STRSI = STRSI+PCSI*(EBIN*EBIN);
            }
          }
        }
      }
      else
      {       
        //  ****  Outer shells
        IO = -IO;
        if(INOUT[IO-1] == 1)
        {
          double H0, H1, H2, S0, S1, S2, R0, R1, R2;
          int IZZ = KZ[M-1][KO];
          int ISH = KS[M-1][KO];
          EBIN = EB[IZZ-1][ISH-1];
          int INDC = IESIF[IZZ-1]-1;
          PCSI = exp(XESI[INDC+I][ISH-1])*STFO[IO-1];
          EINaT1( EE, UI[M-1][KO], WRI[M-1][KO], 0.0, WCCM, H0, H1, H2, S0, S1, S2, R0, R1, R2);
          STOT = F[M-1][KO]*(S0+H0);
          if(STOT > 1.0E-35)
          {
            DFERMI = (SES0[KO]+SEH0[KO])/STOT;
            PCSI = PCSI*DFERMI;
            double RNFO = PCSI/(SES0[KO]+SEH0[KO]);
            XS1SI = XS1SI+SES1[KO]*RNFO;
            XS2SI = XS2SI+SES2[KO]*RNFO;
            XT1SI = XT1SI+SET1[KO]*RNFO;
            XT2SI = XT2SI+SET2[KO]*RNFO;
            XSEIN[I][IO-1] = SEH0[KO]*RNFO;
            STPSI = STPSI+(SES1[KO]+SEH1[KO])*RNFO;
            STRSI = STRSI+(SES2[KO]+SEH2[KO])*RNFO;
          }
          else
          {
            XS1SI = XS1SI+SES1[KO];
            XS2SI = XS2SI+SES2[KO];
            XT1SI = XT1SI+SET1[KO];
            XT2SI = XT2SI+SET2[KO];
            XSEIN[I][IO-1] = SEH0[KO];
            STPSI = STPSI+(SES1[KO]+SEH1[KO]);
            STRSI = STRSI+(SES2[KO]+SEH2[KO]);
          }
        }
        else
        {
          XS1IN = XS1IN+SES1[KO];
          XS2IN = XS2IN+SES2[KO];
          XT1IN = XT1IN+SET1[KO];
          XT2IN = XT2IN+SET2[KO];
          XSEIN[I][IO-1] = SEH0[KO];
          STPIN = STPIN+(SES1[KO]+SEH1[KO]);
          STRIN = STRIN+(SES2[KO]+SEH2[KO]);
        }
      }
    }

    double SSIT = 0.0;
    if(NESI[M-1] > 0)
    {
      for(int IO = 0; IO < NESI[M-1]; IO++)
      {
        ESIAC[M-1][I][IO] = SSIT;
        SSIT = SSIT+XSESI[I][IO];
      }
      double Aux_Double = SSIT*VMOL[M-1];
      if(Aux_Double < 1.0E-35){ Aux_Double = 1.0E-35;}
      
      SEISI[M-1][I] = log(Aux_Double);
    }
    else
    {
      SEISI[M-1][I] = log(1.0E-35);
    }
    if(SSIT > 1.0E-35)
    {
      for(int IO = 0; IO < NESI[M-1]; IO++)
      {
        ESIAC[M-1][I][IO] = ESIAC[M-1][I][IO]/SSIT;
      }
    }
    else
    {
      for(int IO = 0; IO < NESI[M-1]; IO++)
      {
        ESIAC[M-1][I][IO] = 1.0;
      }
    }

    STPSI = STPSI*VMOL[M-1];
    STPIN = STPIN*VMOL[M-1];
    double FNORM = (CSTPE[M-1][I]-STPSI)/STPIN;
 
    double SINT = 0.0;
    if(NEIN[M-1] == 0){ ErrorFunction(1015); return;}
    for(int IO = 0; IO < NEIN[M-1]; IO++)
    {
      if(INOUT[IO] == 0){ XSEIN[I][IO] = FNORM*XSEIN[I][IO];}
      EINAC[M-1][I][IO] = SINT;
      SINT = SINT+XSEIN[I][IO];
    }
    if(SINT > 1.0E-35)
    {
      for(int IO = 0; IO < NEIN[M-1]; IO++)
      {
        EINAC[M-1][I][IO] = EINAC[M-1][I][IO]/SINT;
      }
    }
    else
    {
      for(int IO = 0; IO < NEIN[M-1]; IO++)
      {
        EINAC[M-1][I][IO] = 1.0;
      }
    }
        
    SEHIN[M-1][I] = SINT*VMOL[M-1];
    if(SEHIN[M-1][I] < 1.0E-35){ SEHIN[M-1][I] = 1.0E-35;}
    SEHIN[M-1][I] = log(SEHIN[M-1][I]); 
    TSTRE[M-1][I] = (STRSI+FNORM*STRIN)*VMOL[M-1];
    W1E[M-1][I] = (XS1SI+FNORM*XS1IN)*VMOL[M-1];
    W2E[M-1][I] = (XS2SI+FNORM*XS2IN)*VMOL[M-1];
    T1EI[I] = (XT1SI+FNORM*XT1IN)*VMOL[M-1];
    T2EI[I] = (XT2SI+FNORM*XT2IN)*VMOL[M-1];
  }

  //  ************  Inelastic x-section for positrons.

  DIFMAX = 0.0;
    
  for(int I = 0; I < NDATA; I++)
  {
    if(EIT[I] >= CEGRID_.EL && EIT[I] <= CEGRID_.EU)
    {
      STPI = F3[I]*RHO[M-1]*1.0E6;  // Collision stopping power.
      PINaT(EIT[I], WCCM, XH0, XH1, XH2, XS0, XS1, XS2, XT1, XT2, DELTA, M);
      STPC = (XS1+XH1)*VMOL[M-1];
      if(DIFMAX < 100.0*fabs(STPC-STPI)/STPI){ DIFMAX = 100.0*fabs(STPC-STPI)/STPI;}
    }
  }

  if(DIFMAX > 1.0E-3 && ISCALE == 1)
  {
    ICOR = 1;
    for(int I = 0; I < NDATA; I++)
    {
      FL[I] = log(F3[I]*RHO[M-1]*1.0E6);
    }
    SPLINE( EITL, FL, A, B, C, D, 0.0, 0.0, NDATA);
    if(IRETRN != 0){return;}
  }
  else
  {
    ICOR = 0;
  }

  //  **** Classify inner and outer shells.

  NI = 0;
  NS = 0;
  for(int IO = 0; IO < NO; IO++)
  {
    STFI[IO]=0.0;
    STFO[IO]=0.0;
  }
    
  for(int KO = 0; KO < NOSC[M-1]; KO++)
  {
    int IZZ = KZ[M-1][KO];
    int ISH = KS[M-1][KO];
    if(IZZ < 3 || ISH > 16 || UI[M-1][KO] < 50.0)
    {
      NI = NI+1;
      IPIN[M-1][NI-1] = KO+1;
      ISIP[KO] = -NI;
      INOUT[NI-1] = 0;
      for(int IEL = 0; IEL < NELEM[M-1]; IEL++)
      {
        if(IZZ == IZ[M-1][IEL]){ STFO[NI-1] = STF[M-1][IEL];}
      }
    }
    else if(ISH <= NSPSI[IZZ-1])
    { 
      EBIN = EB[IZZ-1][ISH-1];
      if(EBIN > ECUTR[M-1])
      {
        NS = NS+1;
        IPSI[M-1][NS-1] = KO+1;
        ISIP[KO] = NS;
        for(int IEL = 0; IEL < NELEM[M-1]; IEL++)
        {
          if(IZZ == IZ[M-1][IEL]){ STFI[NS-1] = STF[M-1][IEL];}
        }
      }
      else
      {
        NI = NI+1;
        IPIN[M-1][NI-1] = KO+1;
        ISIP[KO] = -NI;
        INOUT[NI-1] = 1;
        for(int IEL = 0; IEL < NELEM[M-1]; IEL++)
        {
          if(IZZ == IZ[M-1][IEL]){ STFO[NI-1] = STF[M-1][IEL];}
        }
      }
    }
    else
    {
      NI = NI+1;
      IPIN[M-1][NI-1] = KO+1;
      ISIP[KO] = -NI;
      INOUT[NI-1] = 0;
      for(int IEL = 0; IEL < NELEM[M-1]; IEL++)
      {
        if(IZZ == IZ[M-1][IEL]){ STFO[NI-1] = STF[M-1][IEL];}
      }
    }
  }
  NPIN[M-1] = NI;
  NPSI[M-1] = NS;

  //  ****  Simulation arrays.

  for(int I = 0; I < NEGP; I++)
  {
    EE = CEGRID_.ET[I];
    PINaT( EE, WCCM, XH0, XH1, XH2, XS0, XS1, XS2, XT1, XT2, DELTA, M);
    STPC = (XS1+XH1)*VMOL[M-1];
    if(ICOR == 1)
    {
      int J;
      EC = CEGRID_.DLEMP[I];
      FINDI( EITL, EC, NDATA, J);
      STPI = exp(A[J-1]+EC*(B[J-1]+EC*(C[J-1]+EC*D[J-1])));
      FACT = STPI/STPC;
    }
    else
    {
      FACT=1.0;
    }
    CSTPP[M-1][I] = STPC*FACT;

    double XS1SI=0.0;
    double XS2SI=0.0;
    double XT1SI=0.0;
    double XT2SI=0.0;
    double STPSI=0.0;
    double STRSI=0.0;
    double XS1IN=0.0;
    double XS2IN=0.0;
    double XT1IN=0.0;
    double XT2IN=0.0;
    double STPIN=0.0;
    double STRIN=0.0;

    for(int IO = 0; IO < NPIN[M-1]; IO++)
    {
      XSPIN[I][IO] = 0.0;
    }
    for(int IO = 0; IO < NPSI[M-1]; IO++)
    {
      XSPSI[I][IO] = 0.0;
    }

    for(int KO = 0; KO < NOSC[M-1]; KO++)
    {
      int IO = ISIP[KO];
    //  ****  Inner shells
      if(IO > 0)
      {       
        int IZZ = KZ[M-1][KO];
        int ISH = KS[M-1][KO];
        EBIN = EB[IZZ-1][ISH-1];
        if(EBIN < EE)
        {
          int INDC = IPSIF[IZZ-1]-1;
          PCSI = exp(XPSI[INDC+I][ISH-1])*STFI[IO-1];
          if(PCSI > 1.0E-35)
          {
            double H0, H1, H2, S0, S1, S2, R0, R1, R2;
            PINaT1( EE, UI[M-1][KO], WRI[M-1][KO], 0.0, WCCM, H0, H1, H2, S0, S1, S2, R0, R1, R2);
            STOT = F[M-1][KO]*(S0+H0);
            if(STOT > 1.0E-35)
            {
              DFERMI = (SPS0[KO]+SPH0[KO])/STOT;
              PCSI = PCSI*DFERMI;
              XSPSI[I][IO-1] = PCSI;
              double RNFO = PCSI/(SPS0[KO]+SPH0[KO]);
              STPSI = STPSI+(SPS1[KO]+SPH1[KO])*RNFO;
              STRSI = STRSI+(SPS2[KO]+SPH2[KO])*RNFO;
            }
            else
            {
              XSPSI[I][IO-1] = PCSI;
              STPSI = STPSI+PCSI*EBIN;
              STRSI = STRSI+PCSI*(EBIN*EBIN);
            }
          }
        }
      }
      else
      {
          //  ****  Outer shells
        IO = -IO;
        if(INOUT[IO-1] == 1)
        {
          int IZZ = KZ[M-1][KO];
          int ISH = KS[M-1][KO];
          double H0, H1, H2, S0, S1, S2, R0, R1, R2;
          EBIN = EB[IZZ-1][ISH-1];
          int INDC = IPSIF[IZZ-1]-1;
          PCSI = exp(XPSI[INDC+I][ISH-1])*STFO[IO-1];
          PINaT1( EE, UI[M-1][KO], WRI[M-1][KO], 0.0, WCCM, H0, H1, H2, S0, S1, S2, R0, R1, R2);
          STOT = F[M-1][KO]*(S0+H0);
          if(STOT > 1.0E-35)
          {
            DFERMI = (SPS0[KO]+SPH0[KO])/STOT;
            PCSI = PCSI*DFERMI;
            double RNFO = PCSI/(SPS0[KO]+SPH0[KO]);
            XS1SI = XS1SI+SPS1[KO]*RNFO;
            XS2SI = XS2SI+SPS2[KO]*RNFO;
            XT1SI = XT1SI+SPT1[KO]*RNFO;
            XT2SI = XT2SI+SPT2[KO]*RNFO;
            XSPIN[I][IO-1] = SPH0[KO]*RNFO;
            STPSI = STPSI+(SPS1[KO]+SPH1[KO])*RNFO;
            STRSI = STRSI+(SPS2[KO]+SPH2[KO])*RNFO;
          }
          else
          {
            XS1SI = XS1SI+SPS1[KO];
            XS2SI = XS2SI+SPS2[KO];
            XT1SI = XT1SI+SPT1[KO];
            XT2SI = XT2SI+SPT2[KO];
            XSPIN[I][IO-1] = SPH0[KO];
            STPSI = STPSI+(SPS1[KO]+SPH1[KO]);
            STRSI = STRSI+(SPS2[KO]+SPH2[KO]);
          }
        }
        else
        {
          XS1IN = XS1IN+SPS1[KO];
          XS2IN = XS2IN+SPS2[KO];
          XT1IN = XT1IN+SPT1[KO];
          XT2IN = XT2IN+SPT2[KO];
          XSPIN[I][IO-1] = SPH0[KO];
          STPIN = STPIN+(SPS1[KO]+SPH1[KO]);
          STRIN = STRIN+(SPS2[KO]+SPH2[KO]);
        }
      }
    }

    double SSIT = 0.0;
    if(NPSI[M-1] > 0)
    {
      for(int IO = 0; IO < NPSI[M-1]; IO++)
      {
        PSIAC[M-1][I][IO] = SSIT;
        SSIT = SSIT+XSPSI[I][IO];
      }
      SPISI[M-1][I] = SSIT*VMOL[M-1];
      if(SPISI[M-1][I] < 1.0E-35){ SPISI[M-1][I] = 1.0E-35;}
      SPISI[M-1][I] = log(SPISI[M-1][I]);
    }
    else
    {
      SPISI[M-1][I] = log(1.0E-35);
    }
    if(SSIT > 1.0E-35)
    {
      for(int IO = 0; IO < NPSI[M-1]; IO++)
      {
        PSIAC[M-1][I][IO] = PSIAC[M-1][I][IO]/SSIT;
      }
    }
    else
    {
      for(int IO = 0; IO < NPSI[M-1]; IO++)
      {
        PSIAC[M-1][I][IO] = 1.0;
      }
    }

    STPSI = STPSI*VMOL[M-1];
    STPIN = STPIN*VMOL[M-1];
    double FNORM = (CSTPP[M-1][I]-STPSI)/STPIN;

    double SINT = 0.0;
    if(NPIN[M-1] == 0){ ErrorFunction(1016); return;}
    
    for(int IO = 0; IO < NPIN[M-1]; IO++)
    {
      if(INOUT[IO] == 0){ XSPIN[I][IO] = FNORM*XSPIN[I][IO];}
      PINAC[M-1][I][IO] = SINT;
      SINT = SINT+XSPIN[I][IO];
    }
    if(SINT > 1.0E-35)
    {
      for(int IO = 0; IO < NPIN[M-1]; IO++)
      {
        PINAC[M-1][I][IO] = PINAC[M-1][I][IO]/SINT;
      }
    }
    else
    {
      for(int IO = 0;  IO < NPIN[M-1]; IO++)
      {
        PINAC[M-1][I][IO] = 1.0;
      }
    }
        
    SPHIN[M-1][I] = SINT*VMOL[M-1];
    if(SPHIN[M-1][I] < 1.0E-35){ SPHIN[M-1][I] = 1.0E-35;}
    SPHIN[M-1][I] = log(SPHIN[M-1][I]);
    TSTRP[M-1][I] = (STRSI+FNORM*STRIN)*VMOL[M-1];
    W1P[M-1][I] = (XS1SI+FNORM*XS1IN)*VMOL[M-1];
    W2P[M-1][I] = (XS2SI+FNORM*XS2IN)*VMOL[M-1];
    T1PI[I] = (XT1SI+FNORM*XT1IN)*VMOL[M-1];
    T2PI[I] = (XT2SI+FNORM*XT2IN)*VMOL[M-1];
  }

  //  ****  Bremsstrahlung x-section for electrons.

  double STPR;
  DIFMAX=0.0;
  for(int I = 0; I < NDATA; I++)
  {
    if(EIT[I] >= CEGRID_.EL && EIT[I] < CEGRID_.EU)
    {
      STPI = F2[I]*RHO[M-1]*1.0E6;  // Radiative stopping power.
      EBRaT(EIT[I], WCRM, XH0, XH1, XH2, XS1, XS2, M);
      if(IRETRN != 0){return;}
      STPR=(XS1+XH1)*VMOL[M-1];
      if(DIFMAX < fabs(STPR-STPI)/(0.01*STPI)){ DIFMAX = fabs(STPR-STPI)/(0.01*STPI);}
    }
  }

  if(DIFMAX > 1.0E-3 && ISCALE == 1)
  {
    ICOR = 1;
    for(int I = 0; I < NDATA; I++)
    {
      FL[I] = log(F2[I]*RHO[M-1]*1.0E6);
    }
    SPLINE( EITL, FL, A, B, C, D, 0.0, 0.0, NDATA);
    if(IRETRN != 0){return;}
  }
  else
  {
    ICOR = 0;
  }

  for(int I = 0; I < NEGP; I++)
  {
    EBRaT( CEGRID_.ET[I], WCRM, XH0, XH1, XH2, XS1, XS2, M);
    if(IRETRN != 0){return;}
    STPR = (XS1+XH1)*VMOL[M-1];
    if(ICOR == 1)
    {
      int J;
      EC = CEGRID_.DLEMP[I];
      FINDI( EITL, EC, NDATA, J);
      STPI = exp(A[J-1]+EC*(B[J-1]+EC*(C[J-1]+EC*D[J-1])));
      FACT = STPI/STPR;
    }
    else
    {
      FACT=1.0;
    }
    RSTPE[M-1][I] = STPR*FACT;
    TSTPE[M-1][I] = log(CSTPE[M-1][I]+RSTPE[M-1][I]);
    TSTRE[M-1][I] = log(TSTRE[M-1][I]+(XS2+XH2)*VMOL[M-1]*FACT);

    SEHBR[M-1][I] = XH0*VMOL[M-1]*FACT;
    if(SEHBR[M-1][I] < 1.0E-35){SEHBR[M-1][I] = 1.0E-35;}
    SEHBR[M-1][I] = log(SEHBR[M-1][I]);
    if(IBREMS == 1)
    {
      W1E[M-1][I] = W1E[M-1][I]+XS1*VMOL[M-1]*FACT;
      W2E[M-1][I] = W2E[M-1][I]+XS2*VMOL[M-1]*FACT;
    }
  }

  //  ****  Electron range as a function of energy.

  for(int I = 0; I < NEGP; I++)
  {
    F1[I]=1.0/(CSTPE[M-1][I]+RSTPE[M-1][I]);
  }
  SPLINE( CEGRID_.ET, F1, A, B, C, D, 0.0, 0.0, NEGP);
  if(IRETRN != 0){return;}
  RANGE[0][M-1][0] = 1.0E-8;
  RANGEL[0][M-1][0] = log(RANGE[0][M-1][0]);

  double XL;
  double XU;
  double DR;
  for(int I = 1; I < NEGP; I++)
  {
    XL = CEGRID_.ET[I-1];
    XU = CEGRID_.ET[I];
    SINTEG( CEGRID_.ET, A, B, C, D, XL, XU, DR, NEGP);
    if(IRETRN != 0){return;}
    RANGE[0][M-1][I] = RANGE[0][M-1][I-1]+DR;
    RANGEL[0][M-1][I] = log(RANGE[0][M-1][I]);
  }

  //  ****  Electron radiative yield as a function of energy.

  for(int I = 0; I < NEGP; I++)
  {
    F1[I] = RSTPE[M-1][I]/(CSTPE[M-1][I]+RSTPE[M-1][I]);
  }
  RADY[0] = 1.0E-35;
  EBRY[M-1][0] = log(RADY[0]/CEGRID_.ET[0]);
  for(int I = 1; I < NEGP; I++)
  {
    XL = CEGRID_.ET[I-1];
    XU = CEGRID_.ET[I];
    RADY[I] = RADY[I-1]+RMOMX( CEGRID_.ET, F1, XL, XU, NEGP, 0);
    if(IRETRN != 0){return;}
    EBRY[M-1][I] = log(RADY[I]/CEGRID_.ET[I]);
  }

  //  ****  Electron bremss. photon number yield as a function of energy.

  for(int I = 0; I < NEGP; I++)
  {
    F1[I] = exp(SEHBR[M-1][I])/(CSTPE[M-1][I]+RSTPE[M-1][I]);
  }
  RADN[0] = 0.0;
  for(int I = 1; I < NEGP; I++)
  {
    XL = CEGRID_.ET[I-1];
    XU = CEGRID_.ET[I];
    RADN[I] = RADN[I-1]+RMOMX( CEGRID_.ET, F1, XL, XU, NEGP, 0);
    if(IRETRN != 0){return;}
  }
    
  //  ****  Print electron stopping power tables.

  if(INFO >= 3)
  {
    fprintf(IWR, "\n PENELOPE >>>  Stopping powers for electrons\n");
    fprintf(IWR, "\n   Energy        Scol         Srad         range     Rad. Yield   PhotonYield    delta\n    (eV)       (eV/mtu)     (eV/mtu)       (mtu)                    (W>WCR)\n ------------------------------------------------------------------------------------------\n");
    for(int I = 0; I < NEGP; I++)
    {
      fprintf(IWR, "%12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E\n", CEGRID_.ET[I], CSTPE[M-1][I]/RHO[M-1], RSTPE[M-1][I]/RHO[M-1], RANGE[0][M-1][I]*RHO[M-1], RADY[I]/CEGRID_.ET[I], RADN[I], DEL[M-1][I]);
    }
  }
    
  //  ****  Bremss x-section for positrons.

  DIFMAX = 0.0;
  for(int I = 0; I < NDATA; I++)
  {
    if(EIT[I] >= CEGRID_.EL && EIT[I]  < CEGRID_.EU)
    {
      STPI = F4[I]*RHO[M-1]*1.0E6;  // Radiative stopping power.
      PBRaT(EIT[I], WCRM, XH0, XH1, XH2, XS1, XS2, M);
      if(IRETRN != 0){return;}
      STPR = (XS1+XH1)*VMOL[M-1];
      if(DIFMAX < fabs(STPR-STPI)/(0.01*STPI)) { DIFMAX = fabs(STPR-STPI)/(0.01*STPI);}
    }
  }

  if(DIFMAX > 1.0E-3 && ISCALE == 1)
  {
    ICOR = 1;
    for(int I = 0; I < NDATA; I++)
    {
      FL[I] = log(F4[I]*RHO[M-1]*1.0E6);
    }
    SPLINE( EITL, FL, A, B, C, D, 0.0, 0.0, NDATA);
    if(IRETRN != 0){return;}
  }
  else
  {
    ICOR = 0;
  }

  for(int I = 0; I < NEGP; I++)
  {
    PBRaT( CEGRID_.ET[I], WCRM, XH0, XH1, XH2, XS1, XS2, M);
    if(IRETRN != 0){return;}
    STPR = (XS1+XH1)*VMOL[M-1];
    if(ICOR == 1)
    {
      int J;
      EC=CEGRID_.DLEMP[I];
      FINDI(EITL, EC, NDATA, J);
      STPI = exp(A[J-1]+EC*(B[J-1]+EC*(C[J-1]+EC*D[J-1])));
      FACT = STPI/STPR;
    }
    else
    {
      FACT = 1.0;
    }
    RSTPP[M-1][I] = STPR*FACT;
    TSTPP[M-1][I] = log(CSTPP[M-1][I]+RSTPP[M-1][I]);
    TSTRP[M-1][I] = log(TSTRP[M-1][I]+(XS2+XH2)*VMOL[M-1]*FACT);
    SPHBR[M-1][I] = XH0*VMOL[M-1]*FACT;
    if(SPHBR[M-1][I] < 1.0E-35){ SPHBR[M-1][I] = 1.0E-35;}
    SPHBR[M-1][I] = log(SPHBR[M-1][I]);
    if(IBREMS == 1)
    {
      W1P[M-1][I] = W1P[M-1][I]+XS1*VMOL[M-1]*FACT;
      W2P[M-1][I] = W2P[M-1][I]+XS2*VMOL[M-1]*FACT;
    }
  }

  //  ****  Positron annihilation.

  for(int I = 0; I < NEGP; I++)
  {
    double TXAN;
    PANaT( CEGRID_.ET[I], TXAN);
    SPAN[M-1][I] = log(TXAN*ZT[M-1]*VMOL[M-1]);
  }

  //  ****  Positron range as a function of energy.

  for(int I = 0; I < NEGP; I++)
  {
    F3[I] = 1.0/(CSTPP[M-1][I]+RSTPP[M-1][I]);
  }
  SPLINE( CEGRID_.ET, F3, A, B, C, D, 0.0, 0.0, NEGP);
  if(IRETRN != 0){return;}
  RANGE[2][M-1][0] = 1.0E-8;
  RANGEL[2][M-1][0] = log(RANGE[2][M-1][0]);
  for(int I = 1; I < NEGP; I++)
  {
    XL = CEGRID_.ET[I-1];
    XU = CEGRID_.ET[I];
    SINTEG( CEGRID_.ET, A, B, C, D, XL, XU, DR, NEGP);
    if(IRETRN != 0){return;}
    RANGE[2][M-1][I] = RANGE[2][M-1][I-1]+DR;
    RANGEL[2][M-1][I] = log(RANGE[2][M-1][I]);
  }

  //  ****  Positron radiative yield as a function of energy.

  for(int I = 0; I < NEGP; I++)
  {
    F3[I] = RSTPP[M-1][I]/(CSTPP[M-1][I]+RSTPP[M-1][I]);
  }
  RADY[0] = 1.0E-35;
  PBRY[M-1][0] = log(RADY[0]/CEGRID_.ET[0]);
  for(int I = 1; I < NEGP; I++)
  {
    XL = CEGRID_.ET[I-1];
    XU = CEGRID_.ET[I];
    RADY[I] = RADY[I-1]+RMOMX( CEGRID_.ET, F3, XL, XU, NEGP, 0);
    if(IRETRN != 0){return;}
    PBRY[M-1][I] = log(RADY[I]/CEGRID_.ET[I]);
  }

  //  ****  Positron bremss. photon number yield as a function of energy.

  for(int I = 0; I < NEGP; I++)
  {
    F3[I] = exp(SPHBR[M-1][I])/(CSTPP[M-1][I]+RSTPP[M-1][I]);
  }
  RADN[0]=0.0;
  for(int I = 1; I < NEGP; I++)
  {
    XL = CEGRID_.ET[I-1];
    XU = CEGRID_.ET[I];
    RADN[I] = RADN[I-1]+RMOMX( CEGRID_.ET, F3, XL, XU, NEGP, 0);
    if(IRETRN != 0){return;}
  }

  //  ****  Print positron stopping power tables.

  if(INFO >= 3)
  {
    fprintf(IWR, "\n PENELOPE >>>  Stopping powers for positrons\n");
    fprintf(IWR, "\n   Energy        Scol         Srad         range     Rad. Yield   PhotonYield  annih. mfp\n    (eV)       (eV/mtu)     (eV/mtu)       (mtu)                    (W>WCR)      (mtu)\n ------------------------------------------------------------------------------------------\n");
    for(int I = 0; I < NEGP; I++)
    {
      fprintf(IWR, "%12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E\n", CEGRID_.ET[I], CSTPP[M-1][I]/RHO[M-1], RSTPP[M-1][I]/RHO[M-1], RANGE[2][M-1][I]*RHO[M-1], RADY[I]/CEGRID_.ET[I], RADN[I], RHO[M-1]/exp(SPAN[M-1][I]));
    }
  }

  if(INFO >= 3)
  {
    fprintf(IWR, "\n PENELOPE >>>  Soft stopping power and energy straggling\n");
    fprintf(IWR, "\n   Energy       SStp,e-      SStr,e-     STP,e-     SStp,e+      SStr,e+     Stp,e+\n    (eV)       (eV/mtu)    (eV**2/mtu)  soft/tot   (eV/mtu)    (eV**2/mtu)  soft/tot\n ------------------------------------------------------------------------------------\n");
  }
  double FMTU = 1.0/RHO[M-1];
  double FSOFTE;
  double FSOFTP;
  double W1EW;
  double W2EW;
  double W1PW;
  double W2PW;
  for(int I = 0; I < NEGP; I++)
  {
      //  ****  Soft energy-loss events are switched off when W1 is too small.
    FSOFTE = W1E[M-1][I]/(CSTPE[M-1][I]+RSTPE[M-1][I]);
    FSOFTP = W1P[M-1][I]/(CSTPP[M-1][I]+RSTPP[M-1][I]);
    if(W1E[M-1][I] > 1.0E-5*(CSTPE[M-1][I]+RSTPE[M-1][I]))
    {
      W1EW = W1E[M-1][I];
      W2EW = W2E[M-1][I];
      if(W1E[M-1][I] < 1.0E-35){ W1E[M-1][I] = 1.0E-35;}
      W1E[M-1][I] = log(W1E[M-1][I]);
      if(W2E[M-1][I] < 1.0E-35){ W2E[M-1][I] = 1.0E-35;}
      W2E[M-1][I] = log(W2E[M-1][I]);
    }
    else
    {
      W1EW = 0.0;
      W2EW = 0.0;
      W1E[M-1][I] = -100.0;
      W2E[M-1][I] = -100.0;
    }
    if(W1P[M-1][I] > 1.0E-5*(CSTPP[M-1][I]+RSTPP[M-1][I]))
    {
      W1PW = W1P[M-1][I];
      W2PW = W2P[M-1][I];
      if(W1P[M-1][I] < 1.0E-35){ W1P[M-1][I] = 1.0E-35;}
      W1P[M-1][I] = log(W1P[M-1][I]);
      if(W2P[M-1][I] < 1.0E-35){ W2P[M-1][I] = 1.0E-35;}
      W2P[M-1][I] = log(W2P[M-1][I]);
    }
    else
    {
      W1PW = 0.0;
      W2PW = 0.0;
      W1P[M-1][I] = -100.0;
      W2P[M-1][I] = -100.0;
    }
    if(INFO >= 3){ fprintf(IWR, "%12.5E%13.5E%13.5E%10.2E%13.5E%13.5E%10.2E\n", CEGRID_.ET[I], W1EW*FMTU, W2EW*FMTU, FSOFTE, W1PW*FMTU, W2PW*FMTU, FSOFTP);}
  }

  //  ****  Elastic scattering of electrons and positrons.

  EELaR( M, IRD, IWR, INFO);
  if(IRETRN != 0){return;}
  for(int I = 0; I < NEGP; I++)
  {
    TRL1E[M-1][I] = -log(XE1[I]*VMOL[M-1]);
    TRL2E[M-1][I] = -log(XE2[I]*VMOL[M-1]);
    TRL1P[M-1][I] = -log(XP1[I]*VMOL[M-1]);
    TRL2P[M-1][I] = -log(XP2[I]*VMOL[M-1]);
  }
  EELdR( M, IRD, IWR, INFO);  // Uses the ELSEPA database.
  if(IRETRN != 0){return;}
      
  if(INFO >= 3)
  {
    fprintf(IWR, "\n PENELOPE >>>  Soft angular transport coefficients\n");
    fprintf(IWR, "\n   Energy      SITMFP1,e-   SITMFP2,e-   SITMFP1,e+   SITMFP2,e+\n    (eV)        (1/mtu)      (1/mtu)      (1/mtu)      (1/mtu)\n ----------------------------------------------------------------\n");
    for(int I = 0; I < NEGP; I++)
    {
      fprintf(IWR, "%12.5E%13.5E%13.5E%13.5E%13.5E\n", CEGRID_.ET[I], exp(T1E[M-1][I])*FMTU, exp(T2E[M-1][I])*FMTU, exp(T1P[M-1][I])*FMTU, exp(T2P[M-1][I])*FMTU);
    }
  }
    
  //  ****  Print electron mean free path tables.

  if(INFO >= 3)
  {
    fprintf(IWR, "\n PENELOPE >>>  Electron mean free paths (hard events)\n *** NOTE: The MFP for inner-shell ionisation (isi) is listed only for\n           completeness. The MFP for inelastic collisions (in) has been\n           calculated by considering all inelastic events, including isi.\n");
    fprintf(IWR, "\n   Energy        MFPel        MFPin        MFPbr       MFPtot       MFPisi\n    (eV)         (mtu)        (mtu)        (mtu)        (mtu)        (mtu)\n -----------------------------------------------------------------------------\n");
  }

  double FPEL;
  double FPSI;
  double FPIN;
  double FPBR;
  double FPTOT;
  for(int I = 0; I < NEGP; I++)
  {
    FPEL = exp(SEHEL[M-1][I]);
    FPSI = exp(SEISI[M-1][I]);
    FPIN = exp(SEHIN[M-1][I])+FPSI;
    FPBR = exp(SEHBR[M-1][I]);
    FPTOT = FPEL+FPIN+FPBR;
    if(INFO >= 3){ fprintf(IWR, "%12.5E %12.5E %12.5E %12.5E %12.5E %12.5E\n", CEGRID_.ET[I], RHO[M-1]/FPEL, RHO[M-1]/FPIN, RHO[M-1]/FPBR, RHO[M-1]/FPTOT, RHO[M-1]/FPSI);}
    SEAUX[M-1][I] = log(1.0E-35);
    SETOT[M-1][I] = log(FPTOT);
  }

  for(int I = 1; I < NEGP-1; I++)
  {
    if(exp(SETOT[M-1][I]) > 1.005*exp(SETOT[M-1][I-1]) && exp(SETOT[M-1][I]) > 1.005*exp(SETOT[M-1][I+1]) && CEGRID_.ET[I] > EABS[0][M-1] && CEGRID_.ET[I] < 1.0E6)
    {
      fprintf(IWR, "\n WARNING: The electron hard IMFP has a maximum at E = %13.6E eV\n", CEGRID_.ET[I]);
    }
  }

  //  ****  Print positron mean free path tables.

  if(INFO >= 3)
  {
    fprintf(IWR, "\n PENELOPE >>>  Positron mean free paths (hard events)\n *** NOTE: The MFP for inner-shell ionisation (isi) is listed only for\n           completeness. The MFP for inelastic collisions (in) has been\n           calculated by considering all inelastic events, including isi.\n");
    fprintf(IWR, "\n   Energy        MFPel        MFPin        MFPbr        MFPan       MFPtot       MFPisi\n    (eV)         (mtu)        (mtu)        (mtu)        (mtu)        (mtu)        (mtu)\n ------------------------------------------------------------------------------------------\n");
  }
  for(int I = 0; I < NEGP; I++)
  {
    FPEL = exp(SPHEL[M-1][I]);
    FPSI = exp(SPISI[M-1][I]);
    FPIN = exp(SPHIN[M-1][I])+FPSI;
    FPBR = exp(SPHBR[M-1][I]);
    double FPAN = exp(SPAN[M-1][I]);
    FPTOT = FPEL+FPIN+FPBR+FPAN;
    if(INFO >= 3){ fprintf(IWR, "%12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E\n", CEGRID_.ET[I], RHO[M-1]/FPEL, RHO[M-1]/FPIN, RHO[M-1]/FPBR, RHO[M-1]/FPAN, RHO[M-1]/FPTOT, RHO[M-1]/FPSI);}
    SPAUX[M-1][I] = log(1.0E-35);
    SPTOT[M-1][I] = log(FPTOT);
  }

  for(int I = 1; I < NEGP-1; I++)
  {
    if(exp(SPTOT[M-1][I]) > 1.005*exp(SPTOT[M-1][I-1]) && exp(SPTOT[M-1][I]) > 1.005*exp(SPTOT[M-1][I+1]) && CEGRID_.ET[I] > EABS[2][M-1] && CEGRID_.ET[I] < 1.0E6)
    {
      fprintf(IWR, "\n WARNING: The positron hard IMFP has a maximum at E = %13.6E eV\n", CEGRID_.ET[I]);
    }
  }

  //  ****  Correction for the energy dependence of the stopping power and
  //  the energy straggling parameter of soft interactions.

  for(int I = 0; I < NEGP-1; I++)
  {
    DW1EL[M-1][I] = (W1E[M-1][I+1]-W1E[M-1][I])/(CEGRID_.ET[I+1]-CEGRID_.ET[I]);
    DW2EL[M-1][I] = 0.5*(W2E[M-1][I+1]-W2E[M-1][I])/(CEGRID_.ET[I+1]-CEGRID_.ET[I])+DW1EL[M-1][I];
    DW1EL[M-1][I] = 0.5*DW1EL[M-1][I];
    DW1PL[M-1][I] = (W1P[M-1][I+1]-W1P[M-1][I])/(CEGRID_.ET[I+1]-CEGRID_.ET[I]);
    DW2PL[M-1][I] = 0.5*(W2P[M-1][I+1]-W2P[M-1][I])/(CEGRID_.ET[I+1]-CEGRID_.ET[I])+DW1PL[M-1][I];
    DW1PL[M-1][I] = 0.5*DW1PL[M-1][I];
  }
  DW1EL[M-1][NEGP-1] = 0.0;
  DW2EL[M-1][NEGP-1] = 0.0;
  DW1PL[M-1][NEGP-1] = 0.0;
  DW2PL[M-1][NEGP-1] = 0.0;

  //  ************  Photon interaction properties.

  //  ****  Rayleigh scattering.
  GRAaR( M, IRD, IWR, INFO);
  if(IRETRN != 0){return;}
  //  ****  Compton scattering and pair production.
  fscanf(IRD, "%*57c%4d%*[^\n]", &NDATA);
  getc(IRD);
  
  if(INFO >= 2){ fprintf(IWR, "\n *** Compton and pair-production cross sections,  NDATA =%4d\n", NDATA);}
  if(NDATA > NEGP){ ErrorFunction(1017); return;}
  if(INFO >= 2){ fprintf(IWR, "\n  Energy     CS-Comp     CS-pair   CS-triplet\n   (eV)      (cm**2)     (cm**2)     (cm**2)\n ----------------------------------------------\n");}
  for(int I = 0; I < NDATA; I++)
  {
    fscanf(IRD, "%lf %lf %lf %lf%*[^\n]", &EIT[I], &F2[I], &F3[I], &F4[I]);
    getc(IRD);
    
    if(INFO >= 2){ fprintf(IWR, "%10.3E%12.5E%12.5E%12.5E\n", EIT[I], F2[I], F3[I], F4[I]);}
    EITL[I] = log(EIT[I]);
  }

  //  ****  Compton scattering.

  double VMOLL = log(VMOL[M-1]);
  for(int I = 0; I < NDATA; I++)
  {
    FL[I] = F2[I];
    if(FL[I] < 1.0E-35){ FL[I] = 1.0E-35;}
    FL[I] = log(FL[I]);
  }
  SPLINE( EITL, FL, A, B, C, D, 0.0, 0.0, NDATA);
  if(IRETRN != 0){return;}
  for(int I = 0; I < NEGP; I++)
  { 
    int J;     
    EC = CEGRID_.DLEMP[I];
    FINDI( EITL, EC, NDATA, J);
    SGCO[M-1][I] = A[J-1]+EC*(B[J-1]+EC*(C[J-1]+EC*D[J-1]))+VMOLL;
  }

  //  ****  Electron-positron pair and triplet production.

  //  ****  Pairs.
  int NP_S = 0;
  for(int I = 0; I < NDATA; I++)
  {
    if(EIT[I] < 1.023E6) { continue;}
    NP_S = NP_S+1;
    F1[NP_S-1] = EITL[I];
    FACT = pow(EIT[I]/(EIT[I]-1.022E6),3);
        
    FL[NP_S-1] = F3[I];
    if(FL[NP_S-1] < 1.0E-35){ FL[NP_S-1] = 1.0E-35;}
    FL[NP_S-1] = log(FACT*FL[NP_S-1]);
  }
  SPLINE( F1, FL, A, B, C, D, 0.0, 0.0, NP_S);
  if(IRETRN != 0){return;}
  for(int I = 0; I < NEGP; I++)
  {
    if(CEGRID_.ET[I] < 1.023E6)
    {
      SGPP[M-1][I] = -80.6;
    }
    else
    {
      int J;  
      FACT = pow(CEGRID_.ET[I]/(CEGRID_.ET[I]-1.022E6),3);
      EC = CEGRID_.DLEMP[I];
      FINDI( F1, EC, NP_S, J);
      SGPP[M-1][I] = A[J-1]+EC*(B[J-1]+EC*(C[J-1]+EC*D[J-1]))+VMOLL-log(FACT);
    }
  }
  //  ****  Triplets.
  NP_S = 0;
  for(int I = 0; I < NDATA; I++)
  {
    if(EIT[I] < 2.045E6){ continue;}
    NP_S = NP_S+1;
    FACT = pow(EIT[I]/(EIT[I]-2.044E6),3);
    F1[NP_S-1] = EITL[I];

    FL[NP_S-1] = F4[I];
    if(FL[NP_S-1] < 1.0E-35){ FL[NP_S-1] = 1.0E-35;}
    FL[NP_S-1] = log(FACT*FL[NP_S-1]);
  }
  SPLINE( F1, FL, A, B, C, D, 0.0, 0.0, NP_S);
  if(IRETRN != 0){return;}
  for(int I = 0; I < NEGP; I++)
  {
    if(CEGRID_.ET[I] < 2.045E6)
    {
      TRIP[M-1][I] = -80.6;
    }
    else
    {
      int J;
      FACT = pow(CEGRID_.ET[I]/(CEGRID_.ET[I]-2.044E6),3);
      EC = CEGRID_.DLEMP[I];
      FINDI( F1, EC, NP_S, J);
      double TRIPL = exp(A[J-1]+EC*(B[J-1]+EC*(C[J-1]+EC*D[J-1]))+VMOLL-log(FACT));
      double PAIR = exp(SGPP[M-1][I]);
      SGPP[M-1][I] = log(PAIR+TRIPL);
      TRIP[M-1][I] = TRIPL/(PAIR+TRIPL);
    }
  }

  if(NOSCCO[M-1] > 1)
  {
    PTRSH[M-1][0] = FCO[M-1][0];
    for(int I = 1; I < NOSCCO[M-1]; I++)
    {
      PTRSH[M-1][I] = PTRSH[M-1][I-1]+FCO[M-1][I];
    }
    for(int I = 0; I < NOSCCO[M-1]; I++)
    {
      PTRSH[M-1][I] = PTRSH[M-1][I]/PTRSH[M-1][NOSCCO[M-1]-1];
    }
  }
  else
  {
    PTRSH[M-1][0] = 1.0;
  }

  //  ****  Photoelectric absorption.

  GPHaR( M, IRD, IWR, INFO);
  if(IRETRN != 0){return;}

  //  ****  Photon 'range' (= mean free path, slightly underestimated).

  double PRA;
  double PCO;
  double PPP;
  double PPH;
  double PT;
  for(int _KE = 0; _KE < NEGP; _KE++)
  { 
    double ECS;
    GRAaTI( CEGRID_.ET[_KE], ECS, M);
    PRA = ECS*VMOL[M-1];
    PCO = exp(SGCO[M-1][_KE]);
    if(CEGRID_.ET[_KE] < 1.023E6)
    {
      PPP = 1.0E-35;
    }
    else
    {
      PPP = exp(SGPP[M-1][_KE]);
    }
    PPH = SGPH[M-1][_KE];
    PT = (PRA+PCO+PPP+PPH);
    RANGE[1][M-1][_KE] = 1.0/PT;
    RANGEL[1][M-1][_KE] = log(RANGE[1][M-1][_KE]);
  }

  if(INFO >= 3)
  {
    fprintf(IWR, "\n PENELOPE >>>  Photon mass attenuation coefficients\n");
    fprintf(IWR, "\n   Energy      Rayleigh      Compton    Photoelect.     Pair         Total\n    (eV)        (1/mtu)      (1/mtu)      (1/mtu)      (1/mtu)      (1/mtu)\n -----------------------------------------------------------------------------\n");
  }
  for(int I = 0; I < NPHD; I++)
  {
    if(ER[I] >= CEGRID_.EL && ER[I] <= CEGRID_.EU)
    {
      double ECS;
      CEGRID_.XEL = log(ER[I]);
      CEGRID_.XE = 1.0+(CEGRID_.XEL-CEGRID_.DLEMP1)*CEGRID_.DLFC;
      CEGRID_.KE = (int)CEGRID_.XE;
      CEGRID_.XEK = CEGRID_.XE-(double)CEGRID_.KE;
      GRAaTI( ER[I], ECS, M);
      PRA = ECS*VMOL[M-1]/RHO[M-1];
      PCO = exp(SGCO[M-1][CEGRID_.KE-1]+(SGCO[M-1][CEGRID_.KE+1-1]-SGCO[M-1][CEGRID_.KE-1])*CEGRID_.XEK)/RHO[M-1];
      if(ER[I] < 1.022E6)
      {
        PPP = 1.0E-35;
      }
      else
      {
        PPP = exp(SGPP[M-1][CEGRID_.KE-1]+(SGPP[M-1][CEGRID_.KE+1-1]-SGPP[M-1][CEGRID_.KE-1])*CEGRID_.XEK)/RHO[M-1];
      }
      PPH = XSR[I]*VMOL[M-1]/RHO[M-1];
      PT = PRA+PCO+PPP+PPH;
      if(INFO >= 3)
      {
        fprintf(IWR, "%12.5E%13.5E%13.5E%13.5E%13.5E%13.5E\n", ER[I], PRA, PCO, PPH, PPP, PT);
      }
    }
  }

  //  ****  Pair production. Initialisation routine.

  GPPa0(M);

  for(int IE = 0; IE < NEGP; IE++)
  {
    SGAUX[M-1][IE] = log(1.0E-35);
  }

  strcpy(LNAME," PENELOPE (v. 2018)  End of material data file ........");
  fscanf(IRD,"%55[^\n]%*[^\n]",NAME);
  getc(IRD);
  if(strcmp(NAME,LNAME) != 0)
  {
    fprintf(IWR, "\n I/O error. Corrupt material data file.\n");
    fprintf(IWR, "       The last line is: [%s]\n", NAME);
    fprintf(IWR, "      ... and should be: [%s]\n", LNAME);
    ErrorFunction(1018); return;
  }
  fprintf(IWR, "%55s\n", NAME);
  fflush(IWR); // Final flush is needed
}

