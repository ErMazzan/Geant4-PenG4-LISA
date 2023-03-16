
//  *********************************************************************
//                       SUBROUTINE STORES
//  *********************************************************************
void PenPhys::STORES(double EI, double XI, double YI, double ZI, double UI, double VI, double WI, double WGHTI, int KPARI, int ILBI[5], int IPOLI)
{
  //  This subroutine stores the initial state of a new secondary particle
  //  in the secondary stack. The input values are:
  //     EI ........... initial energy.
  //     XI, YI, ZI ... initial position coordinates.
  //     UI, VI, WI ... initial direction cosines.
  //     WGHTI ........ weight (=1 in analogue simulation).
  //     KPARI ........ type of particle (1: electron, 2: photon,
  //                    3: positron).
  //     ILBI(5) ...... particle labels.
  //     IPOLI ........ polarisation flag.

  //  The parameter NMS fixes the size of the secondary stack (i.e. the
  //  maximum number of particles that can be stored). If this number is
  //  exceeded, a warning message is printed on unit 26. When the memory
  //  storage is exhausted, each new secondary particle is stored on the
  //  position of the less energetic secondary electron or photon already
  //  produced, which is thus discarded.

  using namespace PENELOPE_mod;
  using namespace PENERROR_mod;
/////////////////////  using namespace TRACK_mod;
/////////////////////  using namespace SECST;

  int IS, IE, IG;
  double EME, EMG;
  if(SECST_.NSEC < NMS)
  {
    SECST_.NSEC = SECST_.NSEC+1;
    IS = SECST_.NSEC;
  }
  else
  {
    if(IERSEC == 0)
    {
      printf("\n   *** WARNING: (STORES) not enough storage for secondaries.\n                EABS(KPAR,MAT) or the parameter NMS should be increased.\n");
      IERSEC = 1;
    }
    SECST_.NSEC = NMS;
    EME = 1.0E35;
    EMG = 1.0E35;
    IE = 0;
    IG = 0;
    for(int I = 0; I < NMS; I++)
    {
      if(SECST_.KPS[I] == 1)
      {
        if(SECST_.ES[I] < EME)
        {
          EME = SECST_.ES[I];
          IE = I+1;
        }
      }
      else if(SECST_.KPS[I] == 2)
      {
        if(SECST_.ES[I] < EMG)
        {
          EMG = SECST_.ES[I];
          IG = I+1;
        }
      }
    }

    if(IE > 0)
    {
      IS = IE;
    }
    else if(IG > 0)
    {
      IS = IG;
    }
    else
    {
      printf("\n   *** Not enough storage for secondary positrons.\n       JOB ABORTED.\n");
      IS = 0;
      ErrorFunction(1203); return;
    }
  }

  SECST_.ES[IS-1] = EI;
  SECST_.XS[IS-1] = XI;
  SECST_.YS[IS-1] = YI;
  SECST_.ZS[IS-1] = ZI;
  SECST_.US[IS-1] = UI;
  SECST_.VS[IS-1] = VI;
  SECST_.WS[IS-1] = WI;
  SECST_.WGHTS[IS-1] = WGHTI;
  SECST_.KPS[IS-1] = KPARI;
  SECST_.IBODYS[IS-1] = IBODY;
  SECST_.MS[IS-1] = MAT;
  SECST_.ILBS[0][IS-1] = ILBI[0];
  SECST_.ILBS[1][IS-1] = ILBI[1];
  SECST_.ILBS[2][IS-1] = ILBI[2];
  SECST_.ILBS[3][IS-1] = ILBI[3];
  SECST_.ILBS[4][IS-1] = ILBI[4];
  if(IPOLI == 1)
  {
    SECST_.SP1S[IS-1] = SP1;
    SECST_.SP2S[IS-1] = SP2;
    SECST_.SP3S[IS-1] = SP3;
    SECST_.IPOLS[IS-1] = IPOLI;
  }
  else
  {
    SECST_.SP1S[IS-1] = 0.0;
    SECST_.SP2S[IS-1] = 0.0;
    SECST_.SP3S[IS-1] = 0.0;
    SECST_.IPOLS[IS-1] = 0;
  }
  SECST_.PAGES[IS-1] = PAGE;
}

