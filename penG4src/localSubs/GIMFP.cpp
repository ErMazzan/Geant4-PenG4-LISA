
//  *********************************************************************
//                       SUBROUTINE GIMFP
//  *********************************************************************
void PenPhys::GIMFP()
{
  //  This subroutine computes the inverse mean free paths for interactions
  //  of photons with the current energy in material M.

  using namespace PENELOPE_mod;
//////////////  using namespace TRACK_mod;

//////////////  using namespace CEGRID;
  using namespace CGIMFP;
//////////////  using namespace CJUMP0;

  CJUMP0_.P[0] = SGRA[MAT-1][CEGRID_.KE-1];
  CJUMP0_.P[1] = exp(SGCO[MAT-1][CEGRID_.KE-1]+(SGCO[MAT-1][CEGRID_.KE+1-1]-SGCO[MAT-1][CEGRID_.KE-1])*CEGRID_.XEK);
  CJUMP0_.P[2] = SGPH[MAT-1][CEGRID_.KE-1];
  if(E < 1.023E6)
  {
    CJUMP0_.P[3] = 0.0;
  }
  else
  {
    CJUMP0_.P[3] = exp(SGPP[MAT-1][CEGRID_.KE-1]+(SGPP[MAT-1][CEGRID_.KE+1-1]-SGPP[MAT-1][CEGRID_.KE-1])*CEGRID_.XEK);
  }
  CJUMP0_.P[7] = 0.0;
}

