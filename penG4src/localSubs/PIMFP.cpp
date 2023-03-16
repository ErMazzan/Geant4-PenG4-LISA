
//  *********************************************************************
//                       SUBROUTINE PIMFP
//  *********************************************************************
void PenPhys::PIMFP(int IEND)
{
  //  This subroutine computes the inverse mean free paths for hard inter-
  //  actions of positrons with the current energy in material M.

  using namespace PENELOPE_mod;
////////////  using namespace TRACK_mod;

////////////  using namespace CEGRID;
  using namespace CPIMFP;
////////////  using namespace CJUMP0;

  CJUMP0_.P[1] = exp(SPHEL[MAT-1][CEGRID_.KE-1]+(SPHEL[MAT-1][CEGRID_.KE+1-1]-SPHEL[MAT-1][CEGRID_.KE-1])*CEGRID_.XEK);
  CJUMP0_.P[2] = exp(SPHIN[MAT-1][CEGRID_.KE-1]+(SPHIN[MAT-1][CEGRID_.KE+1-1]-SPHIN[MAT-1][CEGRID_.KE-1])*CEGRID_.XEK);
  CJUMP0_.P[3] = exp(SPHBR[MAT-1][CEGRID_.KE-1]+(SPHBR[MAT-1][CEGRID_.KE+1-1]-SPHBR[MAT-1][CEGRID_.KE-1])*CEGRID_.XEK);
  CJUMP0_.P[4] = exp(SPISI[MAT-1][CEGRID_.KE-1]+(SPISI[MAT-1][CEGRID_.KE+1-1]-SPISI[MAT-1][CEGRID_.KE-1])*CEGRID_.XEK);
  CJUMP0_.P[5] = exp(SPAN[MAT-1][CEGRID_.KE-1]+(SPAN[MAT-1][CEGRID_.KE+1-1]-SPAN[MAT-1][CEGRID_.KE-1])*CEGRID_.XEK);
  CJUMP0_.P[7] = 0.0;
  if(IEND == 1){ return;}
  if(W1P[MAT-1][CEGRID_.KE+1-1] > -78.3)
  {
    CJUMP0_.W1 = exp(W1P[MAT-1][CEGRID_.KE-1]+(W1P[MAT-1][CEGRID_.KE+1-1]-W1P[MAT-1][CEGRID_.KE-1])*CEGRID_.XEK);
    CJUMP0_.W2 = exp(W2P[MAT-1][CEGRID_.KE-1]+(W2P[MAT-1][CEGRID_.KE+1-1]-W2P[MAT-1][CEGRID_.KE-1])*CEGRID_.XEK);
  }
  else
  {
    CJUMP0_.W1 = 0.0;
    CJUMP0_.W2 = 0.0;
  }
  if(T1P[MAT-1][CEGRID_.KE+1-1] > -78.3)
  {
    CJUMP0_.T1 = exp(T1P[MAT-1][CEGRID_.KE-1]+(T1P[MAT-1][CEGRID_.KE+1-1]-T1P[MAT-1][CEGRID_.KE-1])*CEGRID_.XEK);
    CJUMP0_.T2 = exp(T2P[MAT-1][CEGRID_.KE-1]+(T2P[MAT-1][CEGRID_.KE+1-1]-T2P[MAT-1][CEGRID_.KE-1])*CEGRID_.XEK);
  }
  else
  {
    CJUMP0_.T1 = 0.0;
    CJUMP0_.T2 = 0.0;
  }
}

