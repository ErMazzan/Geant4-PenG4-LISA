
//  *********************************************************************
//                       SUBROUTINE SECPAR
//  *********************************************************************
void PenPhys::SECPAR(int &LEFT)
{
  //  This subroutine delivers the initial state of a secondary particle
  //  produced during the previous simulation of the shower. This particle
  //  is removed from the secondary stack, so that it will be lost if a new
  //  call to SECPAR is performed before simulating its trajectory up to
  //  the end.

  //  LEFT is the number of particles in the secondary stack at the calling
  //  time. When LEFT=0, the simulation of the shower has been completed.

  using namespace PENELOPE_mod;
///////////////////  using namespace TRACK_mod;

///////////////////  using namespace SECST;

  if(SECST_.NSEC > 0)
  {
    LEFT = SECST_.NSEC;
    E = SECST_.ES[SECST_.NSEC-1];
    X = SECST_.XS[SECST_.NSEC-1];
    Y = SECST_.YS[SECST_.NSEC-1];
    Z = SECST_.ZS[SECST_.NSEC-1];
    U = SECST_.US[SECST_.NSEC-1];
    V = SECST_.VS[SECST_.NSEC-1];
    W = SECST_.WS[SECST_.NSEC-1];
    WGHT = SECST_.WGHTS[SECST_.NSEC-1];
    KPAR = SECST_.KPS[SECST_.NSEC-1];
    IBODY = SECST_.IBODYS[SECST_.NSEC-1];
    MAT = SECST_.MS[SECST_.NSEC-1];
    IPOL = SECST_.IPOLS[SECST_.NSEC-1];
    SP1 = SECST_.SP1S[SECST_.NSEC-1];
    SP2 = SECST_.SP2S[SECST_.NSEC-1];
    SP3 = SECST_.SP3S[SECST_.NSEC-1];
    PAGE = SECST_.PAGES[SECST_.NSEC-1];
    for(int I = 0; I < 5; I++)
    {
      ILB[I] = SECST_.ILBS[I][SECST_.NSEC-1];
    }
    SECST_.NSEC = SECST_.NSEC-1;
  }
  else
  {
    LEFT = 0;
  }
}

