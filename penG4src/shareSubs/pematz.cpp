//  *********************************************************************
//                       SUBROUTINE PEMATW
//  *********************************************************************
//            void  PEMATW(int IFILE, char MFNAME[20])
void PenInterface::PEMATZ(int& midx, FILE* IWR)
{
  //  This subroutine generates the material definition file, a part of the
  //  input data file of PENELOPE.
  //
  //  midx : material ID pre-defined in PENELOPE (Table 7.1) (0<midx<281)
  //
  //  working directory has to be moved to $PENDBASE_DIR before calling
  //  this method.
  //
  using namespace PENELOPE_mod;
  using namespace PENERROR_mod;

  using namespace COMPOS;
  using namespace CADATA;
////////////////  using namespace CEGRID;
  using namespace CEIN;
  using namespace CGCO;
  using namespace CEBR;
    
  char /*PFILE[81],*/ NAME[63];
  const double TREV = 2.0*REV;
  const double FOURPI = 4.0*PI;

/////////  double FBW[30];
  double FF[NOM], UUI[NOM], FFJ0[NOM], WWRI[NOM];
  int KKZ[NOM], KKS[NOM];
  double FFT[NOM], UIT[NOM], WRIT[NOM];
  int KZT[NOM], KST[NOM];
  double FC[NOM], UIC[NOM], FJ0C[NOM];
  int KZC[NOM], KSC[NOM];

  double EGRT[17] = {1.0, 1.25, 1.50, 1.75, 2.00, 2.50, 3.00, 3.50, 4.00, 4.50, 5.00, 5.50, 6.00, 7.00, 8.00, 9.00, 1.00E1};

  //  ****  'Standard' energy grid.
  double EIT[NEGP1], ES[NEGP1], ESS[NEGP1];
  double XGP0[NEGP1], XGT0[NEGP1], PPE[NEGP1], PPT[NEGP1];
  

  //  ****  Lowest negative valences of the elements.

  double VAL[99] = { 0,  0, -1, -2, -3, -2,  0,  0,  0,
  /*1*/        0, -1, -2, -3, -2, -3, -4,  0,  0, -1,
  /*2*/       -2, -3, -2, -2, -2, -2, -2, -2, -2, -1,
  /*3*/       -2, -3, -2, -3, -4,  0,  0, -1, -2, -3,
  /*4*/       -4, -3, -6, -4, -3, -3, -2, -1, -2, -3,
  /*5*/       -2, -3, -4,  0,  0, -1, -2, -3, -3, -3,
  /*6*/       -3, -3, -2, -2, -3, -3, -3, -3, -3, -3,
  /*7*/       -2, -3, -4, -5, -6, -4, -3, -3, -2, -1,
  /*8*/       -1, -1, -2, -3, -2,  0,  0, -1, -2, -3,
  /*9*/       -4, -4, -3, -3, -3, -3, -3, -3, -3, -3};
  //           0   1   2   3   4   5   6   7   8   9
  

  int M=0;
///////  int IREAD;
///////  int IREAD2;
  int IZZ;
  int IORD;
  int NS;
  int IIF;
  bool brkIt;
  double RGROUP;
  int IDNUM = midx;

//MA  printf("      \n Select one option (1 or 2):\n   1: Enter composition data from the keyboard\n   2: Read them from the file pdcompos.pen\n");
//MA  scanf("%d", &IREAD);
//MA  if(IREAD == 1)
//MA  {
//MA ++++++++++++++++ Lines in if-block are deleted
//MA  }
//MA  else
//MA  {
    //
    //  ************  Material data read from file 'pdcompos.pen'.
    //
//MA    printf(" Enter material identification number ...\n");
//MA    scanf("%d", &IDNUM);
    if(IDNUM < 1 || IDNUM > 300){ ErrorFunction(1104); return;}

    FILE* pdcompos = fopen("./pdfiles/pdcompos.pen","r");
    if(pdcompos == NULL){ ErrorFunction(1107); return;}
    for(int I = 0; I < 15; I++) //Pasem 15 linies
    {
      char BUFFER[256];
      fgets(BUFFER,256,pdcompos);
    }
    for(int K1 = 0; K1 < 300; K1++)
    {
      fscanf(pdcompos,"%3d %62c%*[^\n]",&IORD,NAME);
      getc(pdcompos);
      //Append end of string chars
      NAME[62] = '\0';
      double HOLLOW;
      fscanf(pdcompos,"%d %lf %lf %lf%*[^\n]",&NELEM[M],&HOLLOW,&EXPOT[M],&RHO[M]);
      getc(pdcompos);
      if(NELEM[M] > 30){ ErrorFunction(1105); return;}
      if(NELEM[M] < 1){ ErrorFunction(1106); return;}
      for(int I = 0; I < NELEM[M]; I++)
      {
        fscanf(pdcompos,"%d %lf %lf%*[^\n]",&IZ[M][I],&HOLLOW,&STF[M][I]);
        getc(pdcompos);
      }
      if(IORD == IDNUM){ break;}
    }
    fclose(pdcompos);
    printf("  %3d %62s\n", IORD, NAME);

    ZT[M] = 0.0;
    AT[M] = 0.0;
    for(int I = 0; I < NELEM[M]; I++)
    {
      IZZ = IZ[M][I];
      if(IZZ < 1 || IZZ > 99)
      {
        printf("  Element:    (Z=%2d), atoms/molecule = %12.5E\n", IZZ, STF[M][I]);
        ErrorFunction(1108); return;
      }
      printf("  Element: %s (Z=%2d), atoms/molecule =%12.5E\n", LASYMB[IZZ-1], IZZ, STF[M][I]);
      if(STF[M][I] <= 0.0){ ErrorFunction(1109); return;}
      if(I > 0)
      {
        for(int K = 0; K <= I-1; K++)
        {
          if(IZZ == IZ[M][K]){ ErrorFunction(1110); return;}
        }
      }
      ZT[M] = ZT[M]+IZZ*STF[M][I];
      AT[M] = AT[M]+ATW[IZZ-1]*STF[M][I];
    }
    printf("  Density = %12.5E g/cm**3\n", RHO[M]);
    printf("  Number of electrons per molecule = %12.5E\n", ZT[M]);
    if(RHO[M] <= 0.0){ ErrorFunction(1111); return;}
    printf("  Mean excitation energy = %12.5E eV\n", EXPOT[M]);
    VMOL[M] = AVOG*RHO[M]/AT[M];
//  }

  //  ************  Atomic configuration.

  for(int I = 0; I < 99; I++)
  {
    NSHT[I] = 0;
    for(int J = 0; J < 30; J++)
    {
      EB[I][J] = 0.0;
      CP0[I][J] = 0.0;
      IFI[I][J] = 0;
      IKS[I][J] = 0;
    }
  }
  for(int I = 0; I < NELEM[M]; I++)
  {
    IZZ = IZ[M][I];
    //  ****  Loads element data only once. NSHT(IZZ) is used as a status
    //        indicator.
    if(NSHT[IZZ-1] == 0)
    {
      FILE* pdatconf = fopen("./pdfiles/pdatconf.p14", "r");
      for(int J = 0; J < 22; J++)
      {
        char BUFFER[256];
        fgets(BUFFER,256,pdatconf);
      }
      NS = 0;
      int IZZT = 0;
      for(int J = 0; J < 150000; J++)
      {
        int IS, IIZ;
        char CH2[3];
        char CH5[6];
        double EIE, CCP, GA1, GA2;
        if(EOF==fscanf(pdatconf,"%3d %4d %2s %5s %3d %lf %lf %lf %lf%*[^\n]",&IIZ,&IS,CH2,CH5,&IIF,&EIE,&CCP,&GA1,&GA2)){break;}
        getc(pdatconf);
        if(IIZ == IZZ)
        {
          NS = NS+1;
          if(NS > 30)
          {
            printf(" NS =%4d\n", NS);
            ErrorFunction(1112); return;
          }
          if(IS < 1 || IS > 30)
          {
            printf(" IS =%4d\n", IS);
            ErrorFunction(1113); return;
          }
          IZZT = IZZT+IIF;
          EB[IZZ-1][IS-1] = EIE;
          if(GA2 > 0.0)
          {
            ALW[IZZ-1][IS-1] = GA2;
          }
          else if(GA1 > 0.0)
          {
            ALW[IZZ-1][IS-1] = GA1;
          }
          else
          {
            ALW[IZZ-1][IS-1] = 0.0;
          }
          CP0[IZZ-1][IS-1] = CCP;
          IFI[IZZ-1][IS-1] = IIF;
          IKS[IZZ-1][NS-1] = IS;
        }
      }
      NSHT[IZZ-1] = NS;
      if(IZZ != IZZT){ ErrorFunction(1114); return;}
      fclose(pdatconf);
    }
  }

  //  ************  E/P inelastic scattering model parameters.

  //  ****  Set the oscillator table (analogous to ICRU37).

  for(int I = 0; I < NO; I++)
  {
    F[M][I] = 0.0;
    UI[M][I] = 0.0;
    WRI[M][I] = 0.0;
    KZ[M][I] = 0;
    KS[M][I] = 0;
  }
  
  for(int I = 0; I < NOM; I++)
  {
    FF[I] = 0.0;
    UUI[I] = 0.0;
    WWRI[I] = 0.0;
    FFJ0[I] = 0.0;
    KKZ[I] = 0;
    KKS[I] = 0;
  }
  double FT = 0.0;
  //  ****  The 1st oscillator corresponds to the conduction band, which is
  //  tentatively assumed to consist of valence electrons (each atom con-
  //  tributes a number of electrons equal to its lowest chemical valence).
  int NOS = 1;
  FF[0] = 0.0;
  UUI[0] = 0.0;
  FFJ0[0] = 0.0;
  KKZ[0] = 0;
  KKS[0] = 30;
  for(int I = 0; I < NELEM[M]; I++)
  {
    IZZ = IZ[M][I];
    FF[0] = FF[0]+fabs(VAL[IZZ-1])*STF[M][I];
    for(int K = 0; K < 30; K++)
    {
      int JS = IKS[IZZ-1][K];
      if(JS > 0)
      {
        NOS = NOS+1;
        if(NOS > NOM){ ErrorFunction(1115); return;}
        FF[NOS-1] = IFI[IZZ-1][JS-1]*STF[M][I];
        UUI[NOS-1] = EB[IZZ-1][JS-1];
        FFJ0[NOS-1] = CP0[IZZ-1][JS-1];
        KKZ[NOS-1] = IZZ;
        if(IZZ > 2 && JS < 17)
        {
          KKS[NOS-1] = JS;
        }
        else
        {
          KKS[NOS-1] = 30;
        }
        FT = FT+FF[NOS-1];
      }
    }
  }

  if(fabs(FT-ZT[M]) > 1.0E-10*ZT[M]){ ErrorFunction(1116); return;}
  //  ****  Oscillators are sorted by increasing ionisation energies.
  for(int I = 0; I < NOS-1; I++)
  {
    for(int J = I+1; J < NOS; J++)
    {
      if(UUI[I] >= UUI[J])
      {
        double SAVE = UUI[I];
        UUI[I] = UUI[J];
        UUI[J] = SAVE;
        SAVE = FF[I];
        FF[I] = FF[J];
        FF[J] = SAVE;
        SAVE = FFJ0[I];
        FFJ0[I] = FFJ0[J];
        FFJ0[J] = SAVE;
        int ISAVE = KKZ[I];
        KKZ[I] = KKZ[J];
        KKZ[J] = ISAVE;
        ISAVE = KKS[I];
        KKS[I] = KKS[J];
        KKS[J] = ISAVE;
      }
    }
  }

  //  ************  Plasma energy and conduction band excitations.

  OP2[M] = FOURPI*ZT[M]*VMOL[M]*(A0B*A0B*A0B)*(HREV*HREV);
  double OMEGA = sqrt(OP2[M]);
  double EPP = OMEGA*sqrt(FF[0]/ZT[M]);
  printf("\n  Estimated oscillator strength and energy of the plasmon:\n  Fcb = %12.5E,  Wcb = %12.5E, eV,\n  (for insulators, these quantities should be set equal to zero)\n", FF[0], EPP);

//MA  printf("\n  Do you wish to change the Fcb and Wcb values?   (1=yes,2=no)\n  (type 2 if you are not sure...)\n");
//MA  int IPLOSP;
  int IFCB = 0;
//MA  scanf("%d", &IPLOSP);

  double FP, EP;
//MA  if(IPLOSP == 1)
//MA  {
//MA    printf("\n");
//MA    printf(" Enter the oscillator strength Fcb and energy Wcb (in eV) of the plasmon ...\n");
//MA    scanf("%lf %lf", &FP, &EP);
//MA    if(FP < 0.5)
//MA    {
//MA      FP = 0.0;
//MA      EP = 0.0;
//MA    }
//MA    else if(EP < 0.1)
//MA    {
//MA      EP = OMEGA*sqrt(FP/ZT[M]);
//MA    }
//MA  }
//MA  else
//MA  {
    EP = EPP;
    FP = FF[0];
//MA  }
  printf("\n  Fcb =%12.5E, Wcb =%12.5E, eV\n", FP, EP);
  if(FP > ZT[M]+1.0E-13){ ErrorFunction(1117); return;}
  if(EP < 1.0 || FP < 0.5)
  {
    //  ****  Insulator. There is no conduction band.
    for(int J = 0; J < NOS-1; J++)
    {
      FF[J] = FF[J+1];
      UUI[J] = UUI[J+1];
      FFJ0[J] = FFJ0[J+1];
      KKZ[J] = KKZ[J+1];
      KKS[J] = KKS[J+1];
    }
    NOS = NOS-1;
  }
  else
  {
    //  ****  Conductor. Outer shells are 'moved' to the c.b.
    IFCB = 1;
    int IDEAD = 0;
    double FPP = FP;
    int I = 0;
    brkIt = false;
    while(!brkIt)
    {
      brkIt = true;
      I = I+1;
      if(FF[I] < FPP)
      {
        FPP = FPP-FF[I];
        FF[I] = 0.0;
        IDEAD = IDEAD+1;
        brkIt = false;
      }
      else
      {
        FF[I] = FF[I]-FPP;
        if(fabs(FF[I]) < 1.0E-12)
        {
          FP = FP+FF[I];
          FF[I] = 0.0;
          IDEAD = IDEAD+1;
        }
      }
    }
    FF[0] = FP;
    UUI[0] = 0.0;
    WWRI[0] = EP;
    FFJ0[0] = 0.75/sqrt(3.0*PI*PI*VMOL[M]*(A0B*A0B*A0B)*FP);
    KKZ[0] = 0;
    KKS[0] = 30;
    if(IDEAD > 0)
    {
      for(int J = 1; J < NOS-IDEAD; J++)
      {
        FF[J] = FF[J+IDEAD];
        UUI[J] = UUI[J+IDEAD];
        FFJ0[J] = FFJ0[J+IDEAD];
        KKZ[J] = KKZ[J+IDEAD];
        KKS[J] = KKS[J+IDEAD];
      }
      NOS = NOS-IDEAD;
    }
  }
  //  ****  Check f-sum rule.
  double SUM = 0.0;
  double FACT;
  for(int J = 0; J < NOS; J++)
  {
    SUM = SUM+FF[J];
  }
  if(fabs(SUM-ZT[M]) > 1.0E-6*ZT[M]){ ErrorFunction(1118); return;}
  if(fabs(SUM-ZT[M]) > 1.0E-12*ZT[M])
  {
    FACT = ZT[M]/SUM;
    for(int J = 0; J < NOS; J++)
    {
      FF[J] = FACT*FF[J];
    }
  }

  //  ****  Initial parameters for Compton scattering (before grouping).

  int NOSTC = NOS;
  double CSUMT = 0.0;
  for(int I = 0; I < NOSTC; I++)
  {
    FC[I] = FF[I];
    UIC[I] = UUI[I];
    FJ0C[I] = FFJ0[I];
    KZC[I] = KKZ[I];
    KSC[I] = KKS[I];
    CSUMT = CSUMT+FC[I]*FJ0C[I];
  }

  //  ************  Sternheimer's adjustment factor.

  double TST, AAL, AAU, AA = 0.0;
  if(NOS > 1)
  {
    TST = ZT[M]*log(EXPOT[M]);
    AAL = 0.5;
    AAU = 10.0;
    brkIt = false;
    while(!brkIt)
    {
      brkIt = true;
      AA = 0.5*(AAL+AAU);
      SUM = 0.0;
      for(int I = 0; I < NOS; I++)
      {
        if(I == 0 && IFCB == 1)
        {
          SUM = SUM+FF[0]*log(WWRI[0]);
        }
        else
        {
          double WI2 = ((AA*UUI[I])*(AA*UUI[I]))+0.666666666666666*(FF[I]/ZT[M])*(OMEGA*OMEGA);
          WWRI[I] = sqrt(WI2);
          SUM = SUM+FF[I]*log(WWRI[I]);
        }
      }
      if(SUM < TST)
      {
        AAL = AA;
      }
      else
      {
        AAU = AA;
      }
      if(AAU-AAL > 1.0E-14*AA) { brkIt = false;}
    }
  }
  else
  {
    UUI[0] = fabs(UUI[0]);
    WWRI[0] = EXPOT[M];
  }
  
  printf("\n  Sternheimer adjustment factor = %12.5E\n", AA);
  //  ****  Verification.
  double EXPT;
  EXPT = FF[0]*log(WWRI[0]);
  TST = FF[0];
  if(NOS > 1)
  {
    for(int I = 1; I < NOS; I++)
    {
      EXPT = EXPT+FF[I]*log(WWRI[I]);
      TST = TST+FF[I];
    }
  }
  
  if(fabs(TST-ZT[M]) > 1.0E-8*ZT[M])
  {
    printf(" TST-ZT[M] = %12.5E\n", TST-ZT[M]);
    ErrorFunction(1119); return;
  }
  EXPT = exp(EXPT/ZT[M]);
  if(fabs(EXPT-EXPOT[M]) > 1.0E-8*EXPOT[M])
  {
    printf("EXPT-EXPOT(M) =%12.5E\n", EXPT-EXPOT[M]);
    printf("Error in the calculation of the Sternheimer factor.\n");
    printf("\n  Number of oscillators  = %3d\n", NOS);
    for(int I = 0; I < NOS; I++)
    {
      printf("%4d %13.5E %13.5E %13.5E %13.5E %4d %4d", I, FF[I], UUI[I], WWRI[I], FFJ0[I], KKZ[I], KKS[I]);
    }
    ErrorFunction(1120); return;
  }
  
  //  ****  Selection of the lowest ionisation energy for inner shells.
  //  Only the K, L, M and N shells with ionisation energies greater than
  //  that of the N7 shell of the heaviest element in the material are
  //  considered as inner shells. As a result, the inner/outer character
  //  of an atomic shell depends on the composition of the material.
  
  int IZMAX = 0;
  for(int I = 0; I < NELEM[M]; I++)
  {
    if(IZMAX < IZ[M][I]){ IZMAX = IZ[M][I];}
  }
  int JBM = 0;
  for(int J = 0; J < 16; J++)
  {
    if(EB[IZMAX-1][J] > 50.0){ JBM = J;}
  }
  double WISCUT;
  if(50.0 > EB[IZMAX-1][JBM]-0.1){ WISCUT = 50.0;}
  else{ WISCUT = EB[IZMAX-1][JBM]-0.1;}
  printf("\n  Inner-shell lowest energy = %12.5E eV\n", WISCUT);
  WISCUT=(50.0 > 0.7*WISCUT ? 50.0 : 0.7*WISCUT);  // Allows K lines of low-Z atoms.

  int NOST;
  brkIt = false;
  while(!brkIt)
  {
    brkIt = true;
    NOST = NOS;
    for(int I = 0; I < NOST; I++)
    {
      FFT[I] = FF[I];
      UIT[I] = UUI[I];
      WRIT[I] = WWRI[I];
      KZT[I] = KKZ[I];
      KST[I] = KKS[I];
    }
    //  ****  Oscillators are sorted by increasing resonance energies.
    if(NOST > IFCB+1)
    {
      for(int I = IFCB; I < NOST-1; I++)
      {
        for(int J = I+1; J < NOST; J++)
        {
          if(WRIT[I] > WRIT[J])
          {
            double SAVE;
            SAVE = FFT[I];
            FFT[I] = FFT[J];
            FFT[J] = SAVE;
            SAVE = UIT[I];
            UIT[I] = UIT[J];
            UIT[J] = SAVE;
            SAVE = WRIT[I];
            WRIT[I] = WRIT[J];
            WRIT[J] = SAVE;
            int ISAVE = KZT[I];
            KZT[I] = KZT[J];
            KZT[J] = ISAVE;
            ISAVE = KST[I];
            KST[I] = KST[J];
            KST[J] = ISAVE;
          }
        }
      }
    }
      
    //  ****  Oscillators of outer shells with resonance energies differing
    //  by a factor less than RGROUP are grouped as a single oscillator.
      
    RGROUP = 1.05;
    bool brkIt2 = false;
    while(!brkIt2)
    {
      brkIt2 = true;
      int NINSH = 0;
      for(int J = 0; J < NOST; J++)
      {
        if(KST[J] < 17 && UIT[J] > WISCUT){ NINSH = NINSH+1;}
      }
      int IELIM = 0;
      if(NOST > IFCB+1)
      {
        for(int I = IFCB; I < NOST-1; I++)
        {
          if(KST[I] < 17 && UIT[I] > WISCUT){ continue;}
          if(KST[I+1] < 17 && UIT[I+1] > WISCUT){ continue;}
          if(WRIT[I] < 1.0 || WRIT[I+1] < 1.0){ continue;}
          if(WRIT[I+1] > RGROUP*WRIT[I]){ continue;}
          WRIT[I] = exp((FFT[I]*log(WRIT[I])+FFT[I+1]*log(WRIT[I+1]))/(FFT[I]+FFT[I+1]));
          UIT[I] = (FFT[I]*UIT[I]+FFT[I+1]*UIT[I+1])/(FFT[I]+FFT[I+1]);
          FFT[I] = FFT[I]+FFT[I+1];
          if(KZT[I] != KZT[I+1]){ KZT[I] = 0;}
          KST[I] = 30;
          if(I+1 < NOST-1)
          {
            for(int J = I+1; J < NOST-1; J++)
            {
              FFT[J] = FFT[J+1];
              UIT[J] = UIT[J+1];
              WRIT[J] = WRIT[J+1];
              KZT[J] = KZT[J+1];
              KST[J] = KST[J+1];
            }
          }
          IELIM = IELIM+1;
          FFT[NOST-1] = 0.0;
          UIT[NOST-1] = 0.0;
          WRIT[NOST-1] = 0.0;
          KZT[NOST-1] = 0;
          KST[NOST-1] = 0;
        }
      }
      if(IELIM > 0)
      {
        NOST = NOST-IELIM;
        if(NOST > IFCB+NINSH+1){ brkIt2 = false; continue;}
      }
    //  ****  E/P inelastic model parameters transferred to the final
    //        arrays.
      if(NOST < NO)
      {
        if(RGROUP < 1.5)
        {
          RGROUP = (RGROUP*RGROUP);
          if(NOST > IFCB+NINSH+4){ brkIt2 = false; continue;}
        }
        NOSC[M] = NOST;
        for(int I = 0; I < NOSC[M]; I++)
        {
          F[M][I] = FFT[I];
          UI[M][I] = UIT[I];
          WRI[M][I] = WRIT[I];
          if(UI[M][I] < 1.0E-3)
          {
            UI[M][I] = 0.0;
          }
          KZ[M][I] = KZT[I];
          if(UI[M][I] > WISCUT)
          {
            KS[M][I] = KST[I];
          }
          else
          {
            KS[M][I] = 30;
          }
        }
      }
      else
      {
        RGROUP = (RGROUP*RGROUP);
        if(RGROUP > 2.0)
        {
          WISCUT = 1.25*WISCUT;
          printf("  Inner-shell lowest energy = %12.5E eV\n", WISCUT);
          brkIt = false;
          break;
        }
        brkIt2 = false;
        continue;
      }
    }
  }
  
  printf("\n  E/P inel. grouping factor = %12.5E\n", RGROUP);
  
  //  ************  Compton (impulse approximation) parameters.
  
  //  ****  Shells are sorted by increasing ionisation energies.
  if(NOSTC > 1)
  {
    for(int I = 0; I < NOSTC-1; I++)
    {
      for(int J = I+1; J < NOSTC; J++)
      {
        if(UIC[I] > UIC[J])
        {
          double SAVE = FC[I];
          FC[I] = FC[J];
          FC[J] = SAVE;
          SAVE = UIC[I];
          UIC[I] = UIC[J];
          UIC[J] = SAVE;
          SAVE = FJ0C[I];
          FJ0C[I] = FJ0C[J];
          FJ0C[J] = SAVE;
          int ISAVE = KZC[I];
          KZC[I] = KZC[J];
          KZC[J] = ISAVE;
          ISAVE = KSC[I];
          KSC[I] = KSC[J];
          KSC[J] = ISAVE;
        }
      }
    }
  }
  
  //  ****  Outer shells with ionisation energies differing by a factor
  //  less than RGROUP are grouped as a single shell.
  
  RGROUP = 1.05;
  brkIt = false;
  while(!brkIt)
  {
    brkIt = true;
    int NINSH = 0;
    for(int J = 0; J < NOST; J++)
    {
      if(KSC[J] < 17 && UIC[J] > WISCUT){ NINSH = NINSH+1;}
    }
    int IELIM = 0;
    if(NOSTC > IFCB+1)
    {
      for(int I = IFCB; I < NOSTC-1; I++)
      {
        if(KSC[I] < 17 && UIC[I] > WISCUT){ continue;}
        if(KSC[I+1] < 17 && UIC[I+1] > WISCUT){ continue;}
        if(UIC[I] < 1.0 || UIC[I+1] < 1.0){ continue;}
        if(UIC[I+1] > RGROUP*UIC[I]){ continue;}
        UIC[I] = (FC[I]*UIC[I]+FC[I+1]*UIC[I+1])/(FC[I]+FC[I+1]);
        FJ0C[I] = (FC[I]*FJ0C[I]+FC[I+1]*FJ0C[I+1])/(FC[I]+FC[I+1]);
        FC[I] = FC[I]+FC[I+1];
        if(KZC[I] != KZC[I+1]){ KZC[I] = 0;}
        KSC[I] = 30;
        if(I+1 < NOSTC-1)
        {
          for(int J = I+1; J < NOSTC-1; J++)
          {
            FC[J] = FC[J+1];
            UIC[J] = UIC[J+1];
            FJ0C[J] = FJ0C[J+1];
            KZC[J] = KZC[J+1];
            KSC[J] = KSC[J+1];
          }
        }
        IELIM = IELIM+1;
        FC[NOSTC-1] = 0.0;
        UIC[NOSTC-1] = 0.0;
        FJ0C[NOSTC-1] = 0.0;
        KZC[NOSTC-1] = 0;
        KSC[NOSTC-1] = 0;
      }
    }
    if(IELIM > 0)
    {
      NOSTC = NOSTC-IELIM;
      if(NOSTC > IFCB+NINSH+1){ brkIt = false; continue;}
    }
    //  ****  Compton scattering model parameters transferred to the final
    //        arrays.
    if(NOSTC < NOCO)
    {
      if(RGROUP < 1.5)
      {
        RGROUP = (RGROUP*RGROUP);
        if(NOSTC > IFCB+NINSH+4){ brkIt = false; continue;}
      }
      NOSCCO[M] = NOSTC;
      for(int I = 0; I < NOSCCO[M]; I++)
      {
        FCO[M][I] = FC[I];
        UICO[M][I] = UIC[I];
        FJ0[M][I] = FJ0C[I];
        KZCO[M][I] = KZC[I];
        if(UICO[M][I] > WISCUT)
        {
          KSCO[M][I] = KSC[I];
        }
        else
        {
          KSCO[M][I] = 30;
        }
        CSUMT = CSUMT-FCO[M][I]*FJ0[M][I];
      }
      if(fabs(CSUMT) > 1.0E-9)
      {
        printf("  Residual sum = %12.5e\n", fabs(CSUMT));
        ErrorFunction(1121); return;
      }
    }
    else
    {
      RGROUP = (RGROUP*RGROUP);
      brkIt = false;
      continue;
    }
  }
  printf("    Compton grouping factor = %12.5E\n", RGROUP);
  
  //  ************  PENELOPE's input file.
  
  printf("\n PENELOPE''s material data file is being created.\n");
  FILE* AuxFile = IWR;
//MA  FILE* AuxFile;
//MA  if(IFILE == 1)
//MA  {
//MA    AuxFile = fopen(MFNAME, "w");
//MA  }
//MA  else
//MA  {
//MA    printf(" Enter path+name for this file (up to 80 characters) ...\n");
//MA    scanf("%80s", PFILE);
//MA    AuxFile = fopen(PFILE, "w");
//MA  }
  fprintf(AuxFile, " PENELOPE (v. 2018)  Material data file ...............\n");
  fprintf(AuxFile, " Material: %-62s\n", NAME);
  fprintf(AuxFile, " Mass density =%15.8E g/cm**3\n", RHO[M]);
  fprintf(AuxFile, " Number of elements in the molecule = %2d\n", NELEM[M]);
  for(int I = 0; I < NELEM[M]; I++)
  {
    fprintf(AuxFile, "   atomic number =%3d,  atoms/molecule =%15.8E\n", IZ[M][I], STF[M][I]);
  }
  fprintf(AuxFile, " Mean excitation energy =%15.8E eV\n", EXPOT[M]);
  fprintf(AuxFile, " Number of oscillators =%3d (E/P inelastic model)\n", NOSC[M]);
  for(int I = 0; I < NOSC[M]; I++)
  {
    fprintf(AuxFile, "%4d%16.8E%16.8E%16.8E%4d%4d\n", I+1, F[M][I], UI[M][I], WRI[M][I], KZ[M][I], KS[M][I]); 
  }
      
  fprintf(AuxFile, " Number of shells =%3d (Compton IA model)\n", NOSCCO[M]);
  for(int I = 0; I < NOSCCO[M]; I++)
  {
    fprintf(AuxFile, "%4d%16.8E%16.8E%16.8E%4d%4d\n", I+1, FCO[M][I], UICO[M][I], FJ0[M][I], KZCO[M][I], KSCO[M][I]);
    FJ0[M][I] = FJ0[M][I]*SL;
  }
      
  //  ****  Atomic relaxation data.
      
  for(int I = 0; I < NELEM[M]; I++)
  {
    IZZ = IZ[M][I];
    RELAXW(IZZ, AuxFile);
    if(IRETRN != 0){ return;}
  }
  //  ****  Energy grid (standard).
      
  int NES = 0;
  int IGRID = 0;
  double FGRID = 1.0;
  brkIt = false;
  while(!brkIt)
  {
    brkIt = true;
    IGRID = IGRID+1;
    double EV=EGRT[IGRID-1]*FGRID;
    if(IGRID == 17)
    {
      IGRID = 1;
      FGRID = 10.0*FGRID;
    }
    if(EV < 49.0){ brkIt = false; continue;}
    NES = NES+1;
    ES[NES-1] = EV;
    if(EV < 1.0E9){ brkIt = false; continue;}
  }
  //MACG This follows penelope source.
  //MACG Absolute grid must be used to define material file.
  //MACG This method is invoked prior PEINIT, so this is ok.
  CEGRID_.EMIN= 50.0;
  double EMAX = 1.0e9;
  EGRID(CEGRID_.EMIN,EMAX);
      
  double WCRM = 10.0;
  double WCCM = 0.0;
      
  //  **** Electron and positron inner-shell ionisation x-sections.
      
  ESIaW(M, AuxFile);
  if(IRETRN != 0){ return;}
  PSIaW(M, AuxFile);
  if(IRETRN != 0){ return;}
     
  //  ****  Bremsstrahlung emission,
      
  EBRaW(M, AuxFile);
  if(IRETRN != 0){ return;}
  double ZEQ = sqrt(ZBR2[M]);
  BRaAW(ZEQ, AuxFile);
  if(IRETRN != 0){ return;}
     
  fprintf(AuxFile, " *** Stopping powers for electrons and positrons,  NDATA =%4d\n", NES);
  fflush(AuxFile);
  double EE,ESTP,ERSTP,PSTP,PRSTP;
  double XH0,XH1,XH2,XS0,XS1,XS2,XT1,XT2,DELTA;
  for(int IE = 0; IE < NES; IE++)
  {
    EE = ES[IE];
    EINaT(EE,WCCM,XH0,XH1,XH2,XS0,XS1,XS2,XT1,XT2,DELTA,M+1);
    ESTP = (XS1+XH1)*VMOL[M]*1.0E-6/RHO[M];
    EBRaT(EE,WCRM,XH0,XH1,XH2,XS1,XS2,M+1);
    if(IRETRN != 0){ return;}
    ERSTP = (XS1+XH1)*VMOL[M]*1.0E-6/RHO[M];
    PINaT(EE,WCCM,XH0,XH1,XH2,XS0,XS1,XS2,XT1,XT2,DELTA,M+1);
    PSTP = (XS1+XH1)*VMOL[M]*1.0E-6/RHO[M];
    PBRaT(EE,WCRM,XH0,XH1,XH2,XS1,XS2,M+1);
    if(IRETRN != 0){ return;}
    PRSTP = (XS1+XH1)*VMOL[M]*1.0E-6/RHO[M];
    fprintf(AuxFile, "%10.3E%12.5E%12.5E%12.5E%12.5E\n", EE, ESTP, ERSTP, PSTP, PRSTP);
  }
      
  //  **** Electron and positron elastic x-sections.
  EELaW(M,AuxFile);
  if(IRETRN != 0){ return;}
  EELdW(M,AuxFile);  // Uses the ELSEPA database.
  if(IRETRN != 0){ return;}
      
  //  ****  Photon x-sections.
      
  GRAaW(M,AuxFile);
  if(IRETRN != 0){ return;}

  int NPTAB;
  GPPaW(EIT,XGP0,XGT0,NPTAB,M);
  if(IRETRN != 0){ return;}
  for(int I = 0; I < NES; I++)
  {
    PPE[I] = 0.0;
  }
  int NESS;
  MERGE2(ES,PPE,EIT,XGP0,ESS,PPT,NES,NPTAB,NESS);
  if(IRETRN != 0){ return;}
  for(int I = 0; I < NESS; I++)
  {
    XGP0[I] = PPT[I];
  }
  for(int I = 0; I < NES; I++)
  {
    PPE[I] = 0.0;
  }
  MERGE2(ES,PPE,EIT,XGT0,ESS,PPT,NES,NPTAB,NESS);
  if(IRETRN != 0){ return;}
  for(int I = 0; I < NESS; I++)
  {
    XGT0[I] = PPT[I];
  }
      
  fprintf(AuxFile, " *** Compton and pair-production cross sections,  NDATA =%4d\n", NESS);
  for(int IE = 0; IE < NESS; IE++)
  {
    EE = ESS[IE];
    double CSC;
    GCOaT(EE,CSC,M);
    if(CSC < 1.0E-35){ CSC=0.0;}
    if(EE < TREV+5.0){ XGP0[IE] = 0.0;}
    if(EE < 2.0*TREV+10.0){ XGT0[IE] = 0.0;}
    fprintf(AuxFile, "%10.3E%12.5E%12.5E%12.5E\n", EE, CSC, XGP0[IE], XGT0[IE]);
  }
      
  GPHaW(M,AuxFile);
  if(IRETRN != 0){ return;}
      
  fprintf(AuxFile, " PENELOPE (v. 2018)  End of material data file ........\n");
//MA  fclose(AuxFile);
}

