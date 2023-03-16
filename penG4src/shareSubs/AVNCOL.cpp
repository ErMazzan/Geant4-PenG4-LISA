
//  *********************************************************************
//                  FUNCTION AVNCOL
//  *********************************************************************
double PenInterface::AVNCOL(double E, int KPAR, int M, int ICOL)
{
  //  This function computes the mean number of interactions of type ICOL
  //  that particles of type KPAR with initial energy E undergo along their
  //  complete tracks in material M. If ICOL does not correspond to a hard
  //  interaction type, or E is outside the allowed range, the result is
  //  set equal to 0.0D0. For photons (KPAR=2) the result is the average
  //  number of interactions in a mean free path.

  using namespace PENELOPE_mod;
  using namespace PENERROR_mod;
  
/////////  using namespace CEGRID;
  using namespace COMPOS;
  using namespace CEIMFP;
  using namespace CPIMFP;
  using namespace CGIMFP;
  using namespace CGRA01;
  using namespace CGPH00;
  
  double HMF[4];
  double AUX[NEGP];

  double AVNCOL_RETURN = 0.0;
  if(E <= CEGRID_.EL || E >= CEGRID_.EU){ return AVNCOL_RETURN;}

  if(KPAR == 1)
  {
    //  ************  Electrons.
      
    if(ICOL == 2)
    {
      for(int _KE = 0; _KE < NEGP; _KE++)
      {
        AUX[_KE] = exp(SEHEL[M-1][_KE])/(CSTPE[M-1][_KE]+RSTPE[M-1][_KE]);
      }
      AVNCOL_RETURN = RMOMX(CEGRID_.ET,AUX,CEGRID_.EL,E,NEGP,0);
      if(IRETRN != 0){ return AVNCOL_RETURN;}
    }
    else if(ICOL == 3)
    {
      for(int _KE = 0; _KE < NEGP; _KE++)
      {
        AUX[_KE] = (exp(SEHIN[M-1][_KE])+exp(SEISI[M-1][_KE]))/(CSTPE[M-1][_KE]+RSTPE[M-1][_KE]);
      }
      AVNCOL_RETURN = RMOMX(CEGRID_.ET,AUX,CEGRID_.EL,E,NEGP,0);
      if(IRETRN != 0){ return AVNCOL_RETURN;}
    }
    else if(ICOL == 4)
    {
      for(int _KE = 0; _KE < NEGP; _KE++)
      {
        AUX[_KE] = exp(SEHBR[M-1][_KE])/(CSTPE[M-1][_KE]+RSTPE[M-1][_KE]);
      }
      AVNCOL_RETURN = RMOMX(CEGRID_.ET,AUX,CEGRID_.EL,E,NEGP,0);
      if(IRETRN != 0){ return AVNCOL_RETURN;}
    }
    else if(ICOL == 5)
    {
      for(int _KE = 0; _KE < NEGP; _KE++)
      {
        AUX[_KE] = exp(SEISI[M-1][_KE])/(CSTPE[M-1][_KE]+RSTPE[M-1][_KE]);
      }
      AVNCOL_RETURN = RMOMX(CEGRID_.ET,AUX,CEGRID_.EL,E,NEGP,0);
    }
    else if(ICOL == 8)
    {
      for(int _KE = 0; _KE < NEGP; _KE++)
      {
        AUX[_KE] = exp(SEAUX[M-1][_KE])/(CSTPE[M-1][_KE]+RSTPE[M-1][_KE]);
      }
      AVNCOL_RETURN = RMOMX(CEGRID_.ET,AUX,CEGRID_.EL,E,NEGP,0);
      if(IRETRN != 0){ return AVNCOL_RETURN;}
    }
    else if(ICOL == 0)
    {
      for(int _KE = 0; _KE < NEGP; _KE++)
      {
        AUX[_KE] = exp(SETOT[M-1][_KE])/(CSTPE[M-1][_KE]+RSTPE[M-1][_KE]);
      }
      AVNCOL_RETURN = RMOMX(CEGRID_.ET,AUX,CEGRID_.EL,E,NEGP,0);
      if(IRETRN != 0){ return AVNCOL_RETURN;}
    }
  }
  else if(KPAR == 2) 
  {    
    //  ************  Photons.
      
    CEGRID_.XEL = log(E);
    CEGRID_.XE = 1.0+(CEGRID_.XEL-CEGRID_.DLEMP1)*CEGRID_.DLFC;
    CEGRID_.KE = (int)CEGRID_.XE;
    if(CEGRID_.KE < 1){ CEGRID_.KE=1;}
    if(CEGRID_.KE >= NEGP){ CEGRID_.KE=NEGP-1;}
    CEGRID_.XEK = CEGRID_.XE-(double)CEGRID_.KE;
      
    int II = IED[CEGRID_.KE-1];
    int IU = IEU[CEGRID_.KE-1];
    bool brkIt = false;
    while(!brkIt)
    {
      brkIt = true;
      int IT = (II+IU)/2;
      if(CEGRID_.XEL > ERA[IT-1])
      {
        II = IT;
      }
      else
      {
        IU = IT;
      }
      if(IU-II > 1){ brkIt = false; continue;}
    }
    HMF[0] = exp(XSRA[M-1][II-1]+(XSRA[M-1][II+1-1]-XSRA[M-1][II-1])*(CEGRID_.XEL-ERA[II-1])/(ERA[II+1-1]-ERA[II-1]));
    HMF[1] = exp(SGCO[M-1][CEGRID_.KE-1]+(SGCO[M-1][CEGRID_.KE+1-1]-SGCO[M-1][CEGRID_.KE-1])*CEGRID_.XEK);
    double EE, XS;
    GPHaT(EE,XS,M);
    HMF[2] = XS*VMOL[M-1];
    if(E > 1.022E6)
    {
      HMF[3] = exp(SGPP[M-1][CEGRID_.KE-1]+(SGPP[M-1][CEGRID_.KE+1-1]-SGPP[M-1][CEGRID_.KE-1])*CEGRID_.XEK);
    }
    else
    {
      HMF[3] = 0.0;
    }
    double HMF_TOTAL = HMF[0]+HMF[1]+HMF[2]+HMF[3];
    if(HMF_TOTAL > 1.0E-35){ AVNCOL_RETURN = HMF[ICOL-1]/HMF_TOTAL;}
    else{ AVNCOL_RETURN = HMF[ICOL]/1.0E-35;}
  }
  else if(KPAR == 3)
  {      
    //  ************  Positrons.
      
    if(ICOL == 2)
    {
      for(int _KE = 0; _KE < NEGP; _KE++)
      {
        AUX[_KE] = exp(SPHEL[M-1][_KE])/(CSTPP[M-1][_KE]+RSTPP[M-1][_KE]);
      }
      AVNCOL_RETURN = RMOMX(CEGRID_.ET,AUX,CEGRID_.EL,E,NEGP,0);
      if(IRETRN != 0){ return AVNCOL_RETURN;}
    }
    else if(ICOL == 3)
    {
      for(int _KE = 0; _KE < NEGP; _KE++)
      {
        AUX[_KE] = (exp(SPHIN[M-1][_KE])+exp(SPISI[M-1][_KE]))/(CSTPP[M-1][_KE]+RSTPP[M-1][_KE]);
      }
      AVNCOL_RETURN = RMOMX(CEGRID_.ET,AUX,CEGRID_.EL,E,NEGP,0);
      if(IRETRN != 0){ return AVNCOL_RETURN;}
    }
    else if(ICOL == 4)
    {
      for(int _KE = 0; _KE < NEGP; _KE++)
      {
        AUX[_KE] = exp(SPHBR[M-1][_KE])/(CSTPP[M-1][_KE]+RSTPP[M-1][_KE]);
      }
      AVNCOL_RETURN = RMOMX(CEGRID_.ET,AUX,CEGRID_.EL,E,NEGP,0);
      if(IRETRN != 0){ return AVNCOL_RETURN;}
    }
    else if(ICOL == 5)
    {
      for(int _KE = 0; _KE < NEGP; _KE++)
      {
        AUX[_KE] = exp(SPISI[M-1][_KE])/(CSTPP[M-1][_KE]+RSTPP[M-1][_KE]);
      }
      AVNCOL_RETURN = RMOMX(CEGRID_.ET,AUX,CEGRID_.EL,E,NEGP,0);
      if(IRETRN != 0){ return AVNCOL_RETURN;}
    }
    else if(ICOL == 6)
    {
      for(int _KE = 0; _KE < NEGP; _KE++)
      {
        AUX[_KE] = exp(SPAN[M-1][_KE])/(CSTPP[M-1][_KE]+RSTPP[M-1][_KE]);
      }
      AVNCOL_RETURN = RMOMX(CEGRID_.ET,AUX,CEGRID_.EL,E,NEGP,0);
      if(IRETRN != 0){ return AVNCOL_RETURN;}
    }
    else if(ICOL == 8)
    {
      for(int _KE = 0; _KE < NEGP; _KE++)
      {
        AUX[_KE] = exp(SPAUX[M-1][_KE])/(CSTPP[M-1][_KE]+RSTPP[M-1][_KE]);
      }
      AVNCOL_RETURN = RMOMX(CEGRID_.ET,AUX,CEGRID_.EL,E,NEGP,0);
      if(IRETRN != 0){ return AVNCOL_RETURN;}
    }
    else if(ICOL == 0)
    {
      for(int _KE = 0; _KE < NEGP; _KE++)
      {
        AUX[_KE] = exp(SPTOT[M-1][_KE])/(CSTPP[M-1][_KE]+RSTPP[M-1][_KE]);
      }
      AVNCOL_RETURN = RMOMX(CEGRID_.ET,AUX,CEGRID_.EL,E,NEGP,0);
      if(IRETRN != 0){ return AVNCOL_RETURN;}
    }
  }
  
  return AVNCOL_RETURN;
}

