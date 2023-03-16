
//  *********************************************************************
//                       SUBROUTINE RELAX
//  *********************************************************************
void PenPhys::RELAX(int &IZ, int &IS)
{
  //  This subroutine simulates the relaxation of a singly ionized atom of
  //  the element IZ with a vacancy in the IS shell (the K shell or an L
  //  or M subshell). This initial vacancy is filled by electrons from
  //  outer shells through radiative and non-radiative transitions, which
  //  may produce additional vacancies.
  
  //  We use the following notation to designate the possible transitions:
  //  *  Radiative: IS0-IS1 (an electron from the IS1 shell fills the
  //     vacancy in the IS0 shell, leaving a hole in the IS1 shell).
  //  *  Non-radiative: IS0-IS1-IS2 (an electron from the IS1 shell fills
  //     the vacancy in the IS0 shell, and the released energy is taken
  //     away by an electron in the IS2 shell; this process leaves two
  //     holes, in the IS1 and IS2 shells).
  //  The de-excitation cascade (i.e. the set of transitions that occur for
  //  a given initial vacancy) is sampled from the transition probabilities
  //  contained in the Livermore Evaluated Atomic Data Library (EADL). The
  //  energy of the radiation emitted in each transition is read from the
  //  PENELOPE database.

  //  The simulation of the de-excitation cascade is discontinued either
  //  when the K, L and M shells have been filled up or when there is not
  //  enough energy to produce 'active' radiation (with energy larger than
  //  EABS). The excitation energy of the residual ion is assumed to be
  //  deposited locally. We disregard the emission and transport of soft
  //  x-rays and slow electrons, whose energies are less than the binding
  //  energy of the N1 shell of the heavier element in the medium. This
  //  sets a lower limit for the energy interval that can be covered by the
  //  simulation program in a consistent way.
  
  //  De-excitation data for the loaded elements are stored in the common
  //  block /CRELAX/, in a form designed to minimize the amount of memory
  //  and to facilitate the random sampling. The quantities in the common
  //  block are the following:
  //  IFIRST(99,16) ... de-excitation data for a vacancy in the shell IS of
  //     the element IZ start at the position K=IFIRST(IZ,IS) in the
  //     storage arrays. The allowed values for IS are 1 to 16 (K shell
  //     and L, M and N subshells).
  //  ILAST(99,16) ... the de-excitation data for a vacancy in the shell
  //     IS of the element IZ end at the position K=ILAST(IZ,IS) in the
  //     storage arrays.
  //  IS1(K), IS2(K) ... shells that are active in the transition (see the
  //     shell label code below). For radiative transitions, IS2(K)=0.
  //  P(K) ... relative probability for the transition IS-IS1(K)-IS2(K).
  //  ET(K) ... energy of the secondary particle emitted in the transition.
  //  F(K), IAL(K) ... cutoff and alias values (Walker's sampling method).
  
  //  ---------------------------------------------------------------------
  //  Label code IS for electron shells:
  //      1 = K  (1s1/2),     11 = N2 (4p1/2),     21 = O5 (5d5/2),
  //      2 = L1 (2s1/2),     12 = N3 (4p3/2),     22 = O6 (5f5/2),
  //      3 = L2 (2p1/2),     13 = N4 (4d3/2),     23 = O7 (5f7/2),
  //      4 = L3 (2p3/2),     14 = N5 (4d5/2),     24 = P1 (6s1/2),
  //      5 = M1 (3s1/2),     15 = N6 (4f5/2),     25 = P2 (6p1/2),
  //      6 = M2 (3p1/2),     16 = N7 (4f7/2),     26 = P3 (6p3/2),
  //      7 = M3 (3p3/2),     17 = O1 (5s1/2),     27 = P4 (6d3/2),
  //      8 = M4 (3d3/2),     18 = O2 (5p1/2),     28 = P5 (6d5/2),
  //      9 = M5 (3d5/2),     19 = O3 (5p3/2),     29 = Q1 (7s1/2),
  //     10 = N1 (4s1/2),     20 = O4 (5d3/2),     30 = outer shells.
  //  ---------------------------------------------------------------------
  
  using namespace PENELOPE_mod;
  using namespace PENERROR_mod;
////////////////////  using namespace TRACK_mod;

  using namespace CADATA;
  using namespace CECUTR;
  using namespace CRELAX;
////////////////////  using namespace CHIST;

  const double TWOPI = PI+PI;
  double PTIM[256];
  int ISV[256];
  int KPARS;

  //  ****  Initialisation.

  if(IZ < 3 || IS > 16){ return;}
  //  ****  If the shell ionisation energy is less than ECUTR, the cascade
  //        is not followed.
  if(EB[IZ-1][IS-1] < ECUTR[MAT-1]){ return;}
  
  int NV = 1;
  ISV[0] = IS;
  double PAGE0 = PAGE;
  PTIM[0] = PAGE;

  //  ****  Next transition.

  bool brkIt = false;
  while(!brkIt)
  {
    brkIt = true;
    int ISP = ISV[NV-1];
    PAGE = PTIM[NV-1];
    int KF = IFIRST[IZ-1][ISP-1];
    int KL = ILAST[IZ-1][ISP-1];
    NV = NV-1;
      
    double RN, TST;
    int K1;
    if(KL > KF)
    {
    //  ****  Walker's sampling algorithm.
      RN = RAND(1.0)*(double)(KL-KF+1);
      K1 = (int)RN;
      TST = RN-(double)K1;
      if(TST > F[KF+K1-1])
      {
        KS=IAL[KF+K1-1];
      }
      else
      {
        KS = KF+K1;
      }
    }
    else
    {
      KS = KF;
    }
    //  ****  If MODER=0, control is returned to the calling program after
    //  determining the first transition, KS. Useful for testing the random
    //  sampling. For normal operation, we can comment out the following
    //  statement.
    if(MODER == 0){ return;}
    //  ****  If LAGE=.TRUE., the particle age is recorded.
    int NVIS;
    if(LAGE)
    {
      NVIS = 1;
      if(NV > 1)  // Multiple vacancies in the active shell?
      {
        for(int ISC = 0; ISC < NV; ISC++)
        {
          if(ISV[ISC] == ISP){ NVIS = NVIS+1;}
        }
      }
      PAGE = PAGE-(ALW[IZ-1][ISP-1]/double(NVIS))*log(RAND(2.0));
    }
      
      //  ****  Fluorescence radiation.
      
    int IS1K = IS1[KS-1];
    int IS2K = IS2[KS-1];
    if(IS2K == 0)
    {
      KPARS = 2;
      if(IS1K < 17)
      {
        if(EB[IZ-1][IS1K-1] > ECUTR[MAT-1])
        {
          NV = NV+1;
          ISV[NV-1] = IS1K;
          PTIM[NV-1] = PAGE;
        }
      }
    }
    else
    {
      KPARS = 1;
      if(IS1K < 17)
      {
        if(EB[IZ-1][IS1K-1] > ECUTR[MAT-1])
        {
          NV = NV+1;
          ISV[NV-1] = IS1K;
          PTIM[NV-1] = PAGE;
        }
      }
      if(IS2K < 17)
      {
        if(EB[IZ-1][IS2K-1] > ECUTR[MAT-1])
        {
          NV = NV+1;
          ISV[NV-1] = IS2K;
          PTIM[NV-1] = PAGE;
        }
      }
    }
      
      //  ****  The emitted particle is stored in the secondary stack when
      //        its energy ET(K) is greater than EABS.
      
    double WS, SDTS, DF, US, VS;
    if(ET[KS-1] > EABS[KPARS-1][MAT-1])
    {
      //  ****  Initial direction (isotropic).
      WS = -1.0+2.0*RAND(2.0);
      SDTS = sqrt(1.0-WS*WS);
      DF = TWOPI*RAND(3.0);
      US = cos(DF)*SDTS;
      VS = sin(DF)*SDTS;
      CHIST_.ILBA[0] = ILB[0]+1;
      CHIST_.ILBA[1] = KPAR;
      CHIST_.ILBA[3] = IZ*1000000+ISP*10000+IS1K*100+IS2K;

      CHIST_.ILBA[4] = ILB[4];
      STORES(ET[KS-1],X,Y,Z,US,VS,WS,WGHT,KPARS,CHIST_.ILBA,0);
      if(IRETRN != 0){ return;}
    }
      
      //  ****  Are there any unfilled vacancies in inner shells?
      
    if(NV > 0){ brkIt = false; continue;}
  }
  PAGE = PAGE0;
  
}

