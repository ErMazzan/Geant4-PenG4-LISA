
//  *********************************************************************
//                  FUNCTION PHMFP
//  *********************************************************************
double PenInterface::PHMFP(double E, int KPAR, int M, int ICOL)
{
  //  This function computes the mean free path (in cm) of particles of
  //  type KPAR and energy E between hard interactions of kind ICOL in
  //  material M. If ICOL does not correspond to a hard interaction type,
  //  or E is outside the allowed range, the result is set equal to 1.0D35.

  using namespace PENELOPE_mod;

/////////////////  using namespace CEGRID;
  using namespace COMPOS;
  using namespace CEIMFP;
  using namespace CPIMFP;
  using namespace CGIMFP;
  using namespace CGRA01;
  using namespace CGPH00;
  
  double PHMFP_RETURN; //Substitueix a PGMFP en fortran i es el valor que retorna la funcio PHMFP
  if(E <= CEGRID_.EL || E >= CEGRID_.EU)
  {
    PHMFP_RETURN = 1.0E35;
    return PHMFP_RETURN;
  }
  CEGRID_.XEL = log(E);
  CEGRID_.XE = 1.0+(CEGRID_.XEL-CEGRID_.DLEMP1)*CEGRID_.DLFC;
  CEGRID_.KE = (int)CEGRID_.XE;
  if(CEGRID_.KE < 1){ CEGRID_.KE = 1;}
  if(CEGRID_.KE >= NEGP){ CEGRID_.KE = NEGP-1;}
  CEGRID_.XEK = CEGRID_.XE-(double)CEGRID_.KE;
  
  double HMFP = 1.0E-35;
  if(KPAR == 1) // 1 -> electrons
  {   
      //  ************  Electrons.
      
    if(ICOL == 2)
    {
      HMFP = exp(SEHEL[M-1][CEGRID_.KE-1]+(SEHEL[M-1][CEGRID_.KE+1-1]-SEHEL[M-1][CEGRID_.KE-1])*CEGRID_.XEK);
    }
    else if(ICOL == 3)
    {
      HMFP = exp(SEHIN[M-1][CEGRID_.KE-1]+(SEHIN[M-1][CEGRID_.KE+1-1]-SEHIN[M-1][CEGRID_.KE-1])*CEGRID_.XEK);
    }
    else if(ICOL == 4)
    {
      HMFP = exp(SEHBR[M-1][CEGRID_.KE-1]+(SEHBR[M-1][CEGRID_.KE+1-1]-SEHBR[M-1][CEGRID_.KE-1])*CEGRID_.XEK);
    }
    else if(ICOL == 5)
    {
      HMFP = exp(SEISI[M-1][CEGRID_.KE-1]+(SEISI[M-1][CEGRID_.KE+1-1]-SEISI[M-1][CEGRID_.KE-1])*CEGRID_.XEK);
    }
    else if(ICOL == 8)
    {
      HMFP = exp(SEAUX[M-1][CEGRID_.KE-1]+(SEAUX[M-1][CEGRID_.KE+1-1]-SEAUX[M-1][CEGRID_.KE-1])*CEGRID_.XEK);
    }
  }
  else if(KPAR == 2) // 2 -> fotons
  {
      //  ************  Photons.

    if(ICOL == 1)
    {
      int I = IED[CEGRID_.KE-1];
      int IU = IEU[CEGRID_.KE-1];
      int IT;
      bool brkIt = false;
      while(!brkIt)
      {
        brkIt = true;
        IT = (I+IU)/2;
        if(CEGRID_.XEL > ERA[IT-1])
        {
          I = IT;
        }
        else
        {
          IU = IT;
        }
        if(IU-I > 1){ brkIt = false; continue;}
      }
      HMFP = exp(XSRA[M-1][I-1]+(XSRA[M-1][I+1-1]-XSRA[M-1][I-1])*(CEGRID_.XEL-ERA[I-1])/(ERA[I+1-1]-ERA[I-1]));
    }
    else if(ICOL == 2)
    {
      HMFP = exp(SGCO[M-1][CEGRID_.KE-1]+(SGCO[M-1][CEGRID_.KE+1-1]-SGCO[M-1][CEGRID_.KE-1])*CEGRID_.XEK);
    }
    else if(ICOL == 3)
    {
      double XS;
      GPHaT(E,XS,M);
      HMFP = XS*VMOL[M-1];
    }
    else if(ICOL == 4 && E > 1.022E6)
    {
      HMFP = exp(SGPP[M-1][CEGRID_.KE-1]+(SGPP[M-1][CEGRID_.KE+1-1]-SGPP[M-1][CEGRID_.KE-1])*CEGRID_.XEK);
    }
    else if(ICOL == 8)
    {
      HMFP = exp(SGAUX[M-1][CEGRID_.KE-1]+(SGAUX[M-1][CEGRID_.KE+1-1]-SGAUX[M-1][CEGRID_.KE-1])*CEGRID_.XEK);
    }
  }
  else if(KPAR == 3) // 3 -> positrons
  {
      //  ************  Positrons.

    if(ICOL == 2)
    {
      HMFP = exp(SPHEL[M-1][CEGRID_.KE-1]+(SPHEL[M-1][CEGRID_.KE+1-1]-SPHEL[M-1][CEGRID_.KE-1])*CEGRID_.XEK);
    }
    else if(ICOL == 3)
    {
      HMFP = exp(SPHIN[M-1][CEGRID_.KE-1]+(SPHIN[M-1][CEGRID_.KE+1-1]-SPHIN[M-1][CEGRID_.KE-1])*CEGRID_.XEK);
    }
    else if(ICOL == 4)
    {
      HMFP = exp(SPHBR[M-1][CEGRID_.KE-1]+(SPHBR[M-1][CEGRID_.KE+1-1]-SPHBR[M-1][CEGRID_.KE-1])*CEGRID_.XEK);
    }
    else if(ICOL == 5)
    {
      HMFP = exp(SPISI[M-1][CEGRID_.KE-1]+(SPISI[M-1][CEGRID_.KE+1-1]-SPISI[M-1][CEGRID_.KE-1])*CEGRID_.XEK);
    }
    else if(ICOL == 6)
    {
      HMFP = exp(SPAN[M-1][CEGRID_.KE-1]+(SPAN[M-1][CEGRID_.KE+1-1]-SPAN[M-1][CEGRID_.KE-1])*CEGRID_.XEK);
    }
    else if(ICOL == 8)
    {
      HMFP = exp(SPAUX[M-1][CEGRID_.KE-1]+(SPAUX[M-1][CEGRID_.KE+1-1]-SPAUX[M-1][CEGRID_.KE-1])*CEGRID_.XEK);
    }
  }
  if(HMFP > 1.0E-35){ PHMFP_RETURN = 1.0/HMFP;}
  else{ PHMFP_RETURN = 1.0/1.0E-35;}
  return PHMFP_RETURN;
}

