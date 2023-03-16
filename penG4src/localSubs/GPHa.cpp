
//  *********************************************************************
//                       SUBROUTINE GPHa
//  *********************************************************************
void PenPhys::GPHa(double &ES, int &IZZ, int &ISH)
{
//  Simulation of photoelectric absorption in material M.

//  Output arguments:
//    ES .... kinetic energy of the photoelectron.
//    IZZ ... atomic number of the atom where absorption has occurred.
//    ISH ... atomic electron shell that has been ionized.

//  NOTE: JUMP uses a photoelectric cross section that is slightly larger
//  than its 'true' value. To correct for this, the photon is allowed to
//  'survive' a photoelectric event. Survival of the photon is flagged by
//  setting IZZ=0, ISH=0, ES=0.0D0 (the energy E of the photon is kept
//  unaltered.

  using namespace PENELOPE_mod;
/////////////  using namespace TRACK_mod;

  using namespace CECUTR;
/////////////  using namespace CEGRID;
  using namespace CADATA;
  using namespace COMPOS;
  using namespace CGIMFP;
  using namespace CGPH00;
  
  double ACP[35], IP[35];

  //  ****  Partial attenuation coefficients.

  bool brkIt;
  int I, IU, IT;
  double DEE, PCSL;
  double PTOT = 0.0;
  for(int IEL = 0; IEL < NELEM[MAT-1]; IEL++)
  {
    IZZ = IZ[MAT-1][IEL];
    //  ****  Binary search.
    I = IPHF[IZZ-1];
    IU = IPHL[IZZ-1];
    brkIt = false;
    while(!brkIt)
    {
      brkIt = true;
      IT=(I+IU)/2;
      if(CEGRID_.XEL > EPH[IT-1])
      {
        I = IT;
      }
      else
      {
        IU = IT;
      }
      if(IU-I > 1){ brkIt = false; continue;}
    }
    
    IP[IEL] = I;
    DEE = EPH[I+1-1]-EPH[I-1];
    if(DEE > 1.0E-15)
    {
      PCSL = XPH[I-1][0]+(XPH[I+1-1][0]-XPH[I-1][0])*(CEGRID_.XEL-EPH[I-1])/DEE;
    }
    else
    {
      PCSL = XPH[I-1][0];
    }
    PTOT = PTOT+STF[MAT-1][IEL]*exp(PCSL);
    ACP[IEL] = PTOT;
  }
  if(PTOT*VMOL[MAT-1] > SGPH[MAT-1][CEGRID_.KE-1])
  {
    printf("WARNING: SGPH is less than the actual mac.\n");
  }

  //  ****  Sample the active element.

  int IELAC;
  double TST = RAND(1.0)*SGPH[MAT-1][CEGRID_.KE-1]/VMOL[MAT-1];
  bool Return = true;
  for(int IEL = 0; IEL < NELEM[MAT-1]; IEL++)
  {
    if(ACP[IEL] > TST)
    {
      IELAC = IEL+1;
      IZZ = IZ[MAT-1][IEL];
      Return = false;
      break;
    }
  }

  if(Return)
  {
    //  ****  Delta interaction. Introduced to correct for the use of an
    //        upper bound of the photoelectric attenuation coefficient.
    IZZ = 0;  // Flags delta interactions.
    ISH = 0;
    ES = 0.0;
    return;
  }
  
  //  ****  Selection of the active shell.
  I = IP[IELAC-1];
  DEE = EPH[I+1-1]-EPH[I-1];
  double PIS = 0.0;
  brkIt = false;
  int J;
  if(DEE > 1.0E-15)
  {
    PTOT = exp(XPH[I-1][0]+(XPH[I+1-1][0]-XPH[I-1][0])*(CEGRID_.XEL-EPH[I-1])/DEE);
    TST = RAND(2.0)*PTOT;
    for(int IS = 0; IS < NPHS[IZZ-1]; IS++)
    {
      J = (IS+1)+1;
      PCSL = XPH[I-1][J-1]+(XPH[I+1-1][J-1]-XPH[I-1][J-1])*(CEGRID_.XEL-EPH[I-1])/DEE;
      PIS = PIS+exp(PCSL);
      if(PIS > TST)
      {
        ISH = IS+1;
        brkIt = true;
        break;
      }
    }
  }
  else
  {
    PTOT = exp(XPH[I-1][0]);
    TST = RAND(2.0)*PTOT;
    for(int IS = 0; IS < NPHS[IZZ-1]; IS++)
    {
      PIS = PIS+exp(XPH[I-1][IS+1]);
      if(PIS > TST)
      {
        ISH = IS+1;
        brkIt = true;
        break;
      }
    }
  }
  
  if(!brkIt){ ISH = 17;}

  //  ****  Photoelectron emission.

  double EBB;
  if(ISH < 17)
  {
    EBB = EB[IZZ-1][ISH-1];
    if(EBB > ECUTR[MAT-1])
    {
      ES = E-EBB;
    }
    else
    {
      ES = E;
      ISH = 17;
    }
  }
  else
  {
    ES = E;
  }
}

