
//  *********************************************************************
//                       SUBROUTINE EIMFP
//  *********************************************************************
void PenPhys::EIMFP(int IEND)
{
  //  This subroutine computes the inverse mean free paths for hard inter-
  //  actions of electrons with the current energy in material M.

  using namespace PENELOPE_mod;
///////////  using namespace TRACK_mod;

///////////  using namespace CEGRID;
  using namespace CEIMFP;
///////////  using namespace CJUMP0;

  CJUMP0_.P[1] = exp(SEHEL[MAT-1][CEGRID_.KE-1]+(SEHEL[MAT-1][CEGRID_.KE+1-1]-SEHEL[MAT-1][CEGRID_.KE-1])*CEGRID_.XEK);
  CJUMP0_.P[2] = exp(SEHIN[MAT-1][CEGRID_.KE-1]+(SEHIN[MAT-1][CEGRID_.KE+1-1]-SEHIN[MAT-1][CEGRID_.KE-1])*CEGRID_.XEK);
  CJUMP0_.P[3] = exp(SEHBR[MAT-1][CEGRID_.KE-1]+(SEHBR[MAT-1][CEGRID_.KE+1-1]-SEHBR[MAT-1][CEGRID_.KE-1])*CEGRID_.XEK);
  CJUMP0_.P[4] = exp(SEISI[MAT-1][CEGRID_.KE-1]+(SEISI[MAT-1][CEGRID_.KE+1-1]-SEISI[MAT-1][CEGRID_.KE-1])*CEGRID_.XEK);
  CJUMP0_.P[7] = 0.0;
  if(IEND == 1){ return;}
  if(W1E[MAT-1][CEGRID_.KE+1-1] > -78.3)
  {
    CJUMP0_.W1 = exp(W1E[MAT-1][CEGRID_.KE-1]+(W1E[MAT-1][CEGRID_.KE+1-1]-W1E[MAT-1][CEGRID_.KE-1])*CEGRID_.XEK);
    CJUMP0_.W2 = exp(W2E[MAT-1][CEGRID_.KE-1]+(W2E[MAT-1][CEGRID_.KE+1-1]-W2E[MAT-1][CEGRID_.KE-1])*CEGRID_.XEK);
  }
  else
  {
    CJUMP0_.W1 = 0.0;
    CJUMP0_.W2 = 0.0;
  }
  if(T1E[MAT-1][CEGRID_.KE+1-1] > -78.3)
  {
    CJUMP0_.T1 = exp(T1E[MAT-1][CEGRID_.KE-1]+(T1E[MAT-1][CEGRID_.KE+1-1]-T1E[MAT-1][CEGRID_.KE-1])*CEGRID_.XEK);
    CJUMP0_.T2 = exp(T2E[MAT-1][CEGRID_.KE-1]+(T2E[MAT-1][CEGRID_.KE+1-1]-T2E[MAT-1][CEGRID_.KE-1])*CEGRID_.XEK);
  }
  else
  {
    CJUMP0_.T1 = 0.0;
    CJUMP0_.T2 = 0.0;
  }
}

