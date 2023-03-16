
//  *********************************************************************
//                       SUBROUTINE GPHaR
//  *********************************************************************
void PenInterface::GPHaR(int M, FILE* IRD, FILE* IWR, int INFO)
{
  //  This subroutine reads photoelectric cross sections of the elements in
  //  material M and prepares simulation tables.
  
  //  NOTE: The array SGPH(M,IE) defines a piecewise constant function that
  //  is larger than the actual photoelectric cross section. SGPH(M,IE) is
  //  defined as the largest value of the photoelectric x-section in the
  //  energy interval from ET(IE) to ET(IE+1). The photon mean free path is
  //  sampled by using this 'augmented' cross section and, to compensate
  //  for this, the photon survives (i.e., it is not absorbed) with a prob-
  //  ability such that the 'exact' photoelectric attenuation coefficient
  //  is reproduced. This trick allows subroutine JUMP to disregard the
  //  existence of absorption edges in the photoelectric x-section and to
  //  perform the tracking of photons faster.
  
  using namespace PENELOPE_mod;
  using namespace PENERROR_mod;

//////////////  using namespace CEGRID;
  using namespace COMPOS;
  using namespace CECUTR;
  using namespace CADATA;
  using namespace CGIMFP;
  using namespace CGPH00;
  using namespace CGPH01;
  
  char CS5[17][6];

  double XGPHR[NDIM][17], X1[NDIM], Y1[NDIM], X2[NDIM], Y2[NDIM];
  int ISH[17], IZZ, NSHR, NDATA;

  strcpy(CS5[0], "total");
  strcpy(CS5[1], "CS-K ");
  strcpy(CS5[2], "CS-L1");
  strcpy(CS5[3], "CS-L2");
  strcpy(CS5[4], "CS-L3");
  strcpy(CS5[5], "CS-M1");
  strcpy(CS5[6], "CS-M2");
  strcpy(CS5[7], "CS-M3");
  strcpy(CS5[8], "CS-M4");
  strcpy(CS5[9], "CS-M5");
  strcpy(CS5[10], "CS-N1");
  strcpy(CS5[11], "CS-N2");
  strcpy(CS5[12], "CS-N3");
  strcpy(CS5[13], "CS-N4");
  strcpy(CS5[14], "CS-N5");
  strcpy(CS5[15], "CS-N6");
  strcpy(CS5[16], "CS-N7");

  //  ************  Read element x-section tables

  for(int IEL = 0; IEL < NELEM[M-1]; IEL++)
  {
    fscanf(IRD, "%*40c%d%*11c%d%*10c%d%*[^\n]", &IZZ, &NSHR, &NDATA);
    getc(IRD);
    
    if(INFO >= 2){ fprintf(IWR, "\n *** Photoelectric cross sections,  IZ =%3i,  NSHELL =%3i,  NDATA = %4i\n", IZZ, NSHR, NDATA);}

    if(IZZ != IZ[M-1][IEL]){ ErrorFunction(1332); return;}
    if(NDATA > NDIM){ ErrorFunction(1333); return;}
    if(NSHR > 16){ ErrorFunction(1334); return;}

    for(int IS = 0; IS < NSHR+1; IS++)
    {
      fscanf(IRD, "%d", &ISH[IS]);
    }
    fscanf(IRD,"%*[^\n]");
    getc(IRD);
    
    for(int IE = 0; IE < NDATA; IE++)
    {
      fscanf(IRD, "%lf", &ER[IE]);
      for(int IS = 0; IS < NSHR+1; IS++)
      {
        fscanf(IRD, "%lf", &XGPHR[IE][IS]);
      }
      fscanf(IRD,"%*[^\n]");
      getc(IRD);
    }

  //  ****  Remove shells with ionisation energies less than 50 eV.

    int NSHA;
    if(NSHR > 1)
    {
      bool brkIt = false;
      NSHA = NSHR;
      for(int IS = NSHA+1-1; IS >= 1; IS--)
      {
        if(EB[IZZ-1][ISH[IS]-1] < 50.0)
        {
          NSHR = NSHR-1;
        }
        else
        {
          brkIt = true;
          break;
        }
      }
      if(NSHR < 0 && !brkIt){ NSHR = 0;}
    }
  
    if(INFO >= 2)
    {
      fprintf(IWR, "\n   Energy       ");
      for(int IS = 0; IS < NSHR+1; IS++)
      {
        if(IS<NSHR+1-1){fprintf(IWR, "%s       ", CS5[ISH[IS]+1-1]);}
        else{fprintf(IWR, "%s\n", CS5[ISH[IS]+1-1]);}
      }
      for(int IE = 0; IE < NDATA; IE++)
      {
        fprintf(IWR, "%12.5E", ER[IE]);
        for(int IS = 0; IS < NSHR+1; IS++)
        {
          fprintf(IWR, "%12.5E", XGPHR[IE][IS]);
        }
        fprintf(IWR, "\n");
      }
    }

    int IC;
    if(NPHS[IZZ-1] == 0)
    {
      IPHF[IZZ-1] = NCUR+1;
      if(NCUR+NDATA > NTP)
      {
        fprintf(IWR, "Insufficient memory storage in GPHaR.\n");
        fprintf(IWR, "Increase the value of the parameter NTP to %d\n", NCUR+NDATA);
        ErrorFunction(1335); return;
      }
      for(int IE = 0; IE < NDATA; IE++)
      {
        IC = NCUR+IE+1;
        EPH[IC-1] = log(ER[IE]);
        for(int IS = 0; IS < NSHR+1; IS++)
        {
          if( XGPHR[IE][ISH[IS]+1-1] > 1.0E-35){ XPH[IC-1][ISH[IS]+1-1] = log(XGPHR[IE][ISH[IS]+1-1]);}
          else{ XPH[IC-1][ISH[IS]+1-1] = log(1.0E-35);}
        }
      }
      NCUR = NCUR+NDATA;
      IPHL[IZZ-1] = NCUR;
      NPHS[IZZ-1] = NSHR;
    }
  }

  //  ****  Total photoelectric attenuation coefficient.

  IZZ = IZ[M-1][0];
  int IST = IPHF[IZZ-1];
  int LST = IPHL[IZZ-1];
  NPHD = 0;
  for(int I = IST-1; I < LST; I++)
  {
    NPHD = NPHD+1;
    ER[NPHD-1] = exp(EPH[I]);
    XSR[NPHD-1] = STF[M-1][0]*exp(XPH[I][0]);
  }
  int N1, N2;
  if(NELEM[M-1] > 1)
  {
    for(int IEL = 1; IEL < NELEM[M-1]; IEL++)
    {
      N1 = NPHD;
      for(int I = 0; I < N1; I++)
      {
        X1[I] = ER[I];
        Y1[I] = XSR[I];
      }
      IZZ = IZ[M-1][IEL];
      IST = IPHF[IZZ-1];
      LST = IPHL[IZZ-1];
      N2 = 0;
      for(int I = IST-1; I < LST; I++)
      {
        N2 = N2+1;
        X2[N2-1] = exp(EPH[I]);
        Y2[N2-1] = STF[M-1][IEL]*exp(XPH[I][0]);
      }
      MERGE2(X1,Y1,X2,Y2,ER,XSR,N1,N2,NPHD);
      if(IRETRN != 0){ return;}
    }
  }

  //  ****  Total photoelectric cross section at the simulation grid points
  //        (slightly increased to simplify the interpolation).

  for(int I = 0; I < NPHD; I++)
  {
    X1[I] = log(ER[I]);
    Y1[I] = log(XSR[I]*VMOL[M-1]);
  }

  double EG1, EG2, DX, F1, F2;
  for(int IE = 0; IE < NEGP-1; IE++)
  {
    int I1, I2;
    EG1 = CEGRID_.DLEMP[IE];
    FINDI(X1,EG1,NPHD,I1);
    if(I1 == NPHD){ I1 = NPHD-1;}
    DX = X1[I1+1-1]-X1[I1-1];
    if(DX > 1.0E-15)
    {
      F1 = Y1[I1-1]+(Y1[I1+1-1]-Y1[I1-1])*(EG1-X1[I1-1])/DX;
    }
    else
    {
      F1 = Y1[I1-1];
    }
    EG2 = CEGRID_.DLEMP[IE+1];
    FINDI(X1,EG2,NPHD,I2);
    if(I2 == NPHD){ I2 = NPHD-1;}
    DX = X1[I2+1-1]-X1[I2-1];
    if(DX > 1.0E-15)
    {
      F2 = Y1[I2-1]+(Y1[I2+1-1]-Y1[I2-1])*(EG2-X1[I2-1])/DX;
    }
    else
    {
      F2 = Y1[I2-1];
    }
      //  ****  To avoid interpolating the attenuation coefficient tables, we
      //        replace the photoelectric inverse mean free path (imfp) in each
      //        energy interval by its upper bound. The increase of the imfp
      //        is interpreted as the imfp of delta interactions.
      
    double FM;

    if(F1 > F2){ FM = F1;}
    else{ FM = F2;}
    if(I1+1 < I2)
    {
      for(int I = I1+1-1; I < I2; I++)
      {
        if(FM < Y1[I]){ FM = Y1[I];}        
      }
    }
    SGPH[M-1][IE] = exp(FM);
  }
  SGPH[M-1][NEGP-1] = SGPH[M-1][NEGP-1-1];
}

