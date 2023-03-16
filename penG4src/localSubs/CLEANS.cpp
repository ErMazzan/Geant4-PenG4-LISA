
//  *********************************************************************
//                       SUBROUTINE CLEANS
//  *********************************************************************
void PenPhys::CLEANS()
{
  //  This subroutine initializes the secondary stack. It must be called
  //  before starting the simulation of each primary track.

  using namespace PENELOPE_mod;
//////////////  using namespace SECST;

  SECST_.NSEC = 0;
}

