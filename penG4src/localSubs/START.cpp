
//  *********************************************************************
//                       SUBROUTINE START
//  *********************************************************************
void PenPhys::START()
{
  //  This subroutine forces the next event to be an artificial soft event.
  //  It must be called when a new (primary or secondary) particle track is
  //  is started, and when it crosses an interface.

  using namespace PENELOPE_mod;
  using namespace PENERROR_mod;
//////////////  using namespace TRACK_mod;
  
//////////////  using namespace CEGRID;
//////////////  using namespace CJUMP1;

  if(E < CEGRID_.EMIN || E > 0.99999999*CEGRID_.EU)
  {
    printf("\n   *** Energy out of range. KPAR = %2d,  E = %12.5E eV\n       ILB s = %4d %4d %4d %11d %11d\n       EMIN = %12.5E eV,  EMAX = %12.5E eV.\n       Check the values of EABS(KPAR,M) and EMAX.\n", KPAR, E, ILB[0], ILB[1], ILB[2], ILB[3], ILB[4], CEGRID_.EL, CEGRID_.EU);
      
    ErrorFunction(1201); return;
  }
  MHINGE = 0;
  CJUMP1_.ELAST1 = E+1.0E30;
  CJUMP1_.ELAST2 = CJUMP1_.ELAST1;
}

