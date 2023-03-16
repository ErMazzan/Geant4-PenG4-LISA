
//  *********************************************************************
//                       SUBROUTINE EGRID
//  *********************************************************************
void PenInterface::EGRID(double EMINu, double EMAXu)
{
  
  //This subroutine sets the energy grid where transport functions are
  //tabulated. The grid is logarithmically spaced and we assume that it
  //is dense enough to permit accurate linear log-log interpolation of
  //he tabulated functions.

  using namespace PENELOPE_mod;
  using namespace PENERROR_mod;

//////////////  using namespace CEGRID;

  //  ****  Consistency of the interval end-points.

  if(EMINu < (double)MINEGRID){ EMINu = (double)MINEGRID;}
  if(EMINu > EMAXu - 1.0) //Mirem si la distancia entre EMINu i EMAXu, es d'1 eV com a minim
  {
  	printf("   EMIN =%11.4E eV, EMAX =%11.4E eV\n", EMINu, EMAXu);
	  ErrorFunction(1007);return;
  }

  //  ****  Energy grid points.

  CEGRID_.EMIN = EMINu;
  CEGRID_.EL = 0.99999*EMINu;
  CEGRID_.EU = 1.00001*EMAXu;
  CEGRID_.DLFC = log(CEGRID_.EU/CEGRID_.EL)/double(NEGP-1);
  CEGRID_.DLEMP1 = log(CEGRID_.EL);
  CEGRID_.DLEMP[0] = CEGRID_.DLEMP1;
  CEGRID_.ET[0] = CEGRID_.EL;
  for(int I = 1; I < NEGP; I++)
  {
    CEGRID_.DLEMP[I] = CEGRID_.DLEMP[I-1]+CEGRID_.DLFC;
    CEGRID_.ET[I] = exp(CEGRID_.DLEMP[I]);
  }
  CEGRID_.DLFC = (double)1.0/CEGRID_.DLFC;

  //  NOTE: To determine the interval KE where the energy E is located, we
  //  do the following,
  //     XEL=LOG(E)
  //     XE=1.0D0+(XEL-DLEMP1)*DLFC
  //     KE=XE
  //     XEK=XE-KE  ! 'fractional' part of XE (used for interpolation).
  //
}

