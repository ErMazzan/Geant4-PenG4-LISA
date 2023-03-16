
//  *********************************************************************
//                       SUBROUTINE DCSEL0
//  *********************************************************************
void PenInterface::DCSEL0(double &E, int IELEC)
{
  //     This subroutine computes a table of the molecular elastic DCSs for
  //  electrons (IELEC=-1) or positrons (IELEC=+1) with kinetic energy  E
  //  (in eV) by log-log cubic spline interpolation in E of the DCS table
  //  prepared by subroutine ELINIT.

  using namespace PENERROR_mod;
  using namespace CDCSEP;
  
  const int NE=96;
  double Y[NE], A[NE], B[NE], C[NE], D[NE];

  if(E < 49.999 || E > 1.0001E8)
  {
    printf(" Error in DCSEL0: Energy out of range.\n");
    ErrorFunction(1354); return;
  }
  
  double EL = log(E);
  int JE;
  FINDI(ETL,EL,NE,JE);

  for(int IA = 0; IA < NA; IA++)
  {
    for(int IE = 0; IE < NE; IE++)
    {
      if(IELEC == -1)
    {
      Y[IE] = log(EDCS[IE][IA]);
    }
      else
    {
      Y[IE] = log(PDCS[IE][IA]);
    }
    }
    SPLINE(ETL,Y,A,B,C,D,0.0,0.0,NE);
    if(IRETRN != 0){ return;}
    DCSIL[IA] = A[JE-1]+EL*(B[JE-1]+EL*(C[JE-1]+EL*D[JE-1]));
    DCSI[IA] = exp(DCSIL[IA]);
  }

  for(int IE = 0; IE < NE; IE++)
  {
    if(IELEC == -1)
    {
      Y[IE] = log(ECS[IE]);
  }
    else
  {
      Y[IE] = log(PCS[IE]);
    }
  }
  SPLINE(ETL,Y,A,B,C,D,0.0,0.0,NE);
  if(IRETRN != 0){ return;}
  CSI = exp(A[JE-1]+EL*(B[JE-1]+EL*(C[JE-1]+EL*D[JE-1])));

  for(int IE = 0; IE < NE; IE++)
  {
    if(IELEC == -1)
  {
      Y[IE] = log(ETCS1[IE]);
  }
    else
  {
      Y[IE] = log(PTCS1[IE]);
    }
  }
  SPLINE(ETL,Y,A,B,C,D,0.0,0.0,NE);
  if(IRETRN != 0){ return;}
  TCS1I = exp(A[JE-1]+EL*(B[JE-1]+EL*(C[JE-1]+EL*D[JE-1])));

  for(int IE = 0; IE < NE; IE++)
  {
    if(IELEC == -1)
  {
      Y[IE] = log(ETCS2[IE]);
  }
    else
  {
      Y[IE] = log(PTCS2[IE]);
    }
  }
  SPLINE(ETL,Y,A,B,C,D,0.0,0.0,NE);
  if(IRETRN != 0){ return;}
  TCS2I = exp(A[JE-1]+EL*(B[JE-1]+EL*(C[JE-1]+EL*D[JE-1])));
}

