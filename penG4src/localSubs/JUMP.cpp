//  *********************************************************************
//                       SUBROUTINE JUMP
//  *********************************************************************
void PenPhys::JUMP(double &DSMAX, double &DS)
{
  //  Calculation of the free path from the starting point to the position
  //  of the next event and of the probabilities of occurrence of different
  //  events.

  //  Arguments:
  //    DSMAX .... maximum allowed step length (input),
  //    DS ....... segment length (output).

  //  Output, through module PENELOPE_mod:
  //    E0STEP ... energy at the beginning of the segment,
  //    DESOFT ... energy loss due to soft interactions along the step,
  //    SSOFT .... stopping power due to soft interactions,
  //               = DESOFT/step_length.

  using namespace PENELOPE_mod;
////////////  using namespace TRACK_mod;

////////////  using namespace CEGRID;
  using namespace CEIMFP;
  using namespace CPIMFP;
////////////  using namespace CJUMP0;
////////////  using namespace CJUMP1;

  if(KPAR == 1)
  {

      //  ************  Electrons (KPAR=1).

    if(MHINGE == 1)
    {
      if(E < CJUMP1_.ELAST1)
      {
        CEGRID_.XEL = log(E);
        CEGRID_.XE = 1.0+(CEGRID_.XEL-CEGRID_.DLEMP1)*CEGRID_.DLFC;
        CEGRID_.KE = (int)CEGRID_.XE;
        CEGRID_.XEK = CEGRID_.XE-(double)CEGRID_.KE;
        EIMFP(1);
        CJUMP1_.ELAST1 = E;
      }
      DS = CJUMP0_.DSR;
      return;
    }

    E0STEP = E;
    if(E < CJUMP1_.ELAST2)
    {
      CEGRID_.XEL = log(E);
      CEGRID_.XE = 1.0+(CEGRID_.XEL-CEGRID_.DLEMP1)*CEGRID_.DLFC;
      CEGRID_.KE = (int)CEGRID_.XE;
      CEGRID_.XEK = CEGRID_.XE-(double)CEGRID_.KE;
      EIMFP(2);
      CJUMP1_.ELAST2 = E;
      CJUMP1_.ELAST1 = E;
    }

      //  ****  Inverse hard mean free path (interaction probability per unit
      //        path length).

    CJUMP0_.ST = CJUMP0_.P[1]+CJUMP0_.P[2]+CJUMP0_.P[3]+CJUMP0_.P[4]+CJUMP0_.P[7];
    double DSMAXP = DSMAX;
    //  ****  The value of DSMAXP is randomised to eliminate dose artifacts
    //        at the end of the first step.
    DSMAXP = (0.75+RAND(1.0)*0.5)*DSMAXP;

      //  ****  Soft stopping interactions.
      //        KSOFTI=1, soft stopping is active,
      //        KSOFTI=0, soft stopping is not active.

    if(CJUMP0_.W1 > 1.0E-20)
    {
      CJUMP1_.KSOFTI = 1;
    //  ****  The maximum step length, DSMAXP, is determined in terms of the
    //        input DSMAX value (which is specified by the user) and the mean
    //        free path for hard interactions (1/ST).
      double DSMC = 4.0/CJUMP0_.ST;
      if(DSMAXP > DSMC)
      {
        DSMAXP = DSMC;
      }
      else if(DSMAXP < 1.0E-8)
      {
        DSMAXP = DSMC;
      }

    //  ****  Upper bound for the interaction probability along the step
    //        (including soft energy straggling).

      double EDE0 = CJUMP0_.W1*DSMAXP;
      double VDE0 = CJUMP0_.W2*DSMAXP;

      double FSEDE = 1.0-DW1EL[MAT-1][CEGRID_.KE-1]*EDE0;
      if(FSEDE < 0.75){ FSEDE = 0.75;}

      double FSVDE = 1.0-DW2EL[MAT-1][CEGRID_.KE-1]*EDE0;
      if(FSVDE < 0.75){ FSVDE = 0.75;}
    
      double EDEM = EDE0*FSEDE;
      double VDEM = VDE0*FSVDE;
      double W21 = VDEM/EDEM;
      double ELOWER;
      if(EDEM > 9.0*W21)
      {
        ELOWER = E-(EDEM+3.0*sqrt(VDEM));
        if(ELOWER < CEGRID_.EMIN){ ELOWER = CEGRID_.EMIN;}
      }
      else if(EDEM > 3.0*W21)
      {
        ELOWER = E-(EDEM+sqrt(3.0*VDEM));
        if(ELOWER < CEGRID_.EMIN){ ELOWER = CEGRID_.EMIN;}
      }
      else
      {
        ELOWER = E-1.5*(EDEM+W21);
        if(ELOWER < CEGRID_.EMIN){ ELOWER = CEGRID_.EMIN;}
      }
      double XE1 = 1.0+(log(ELOWER)-CEGRID_.DLEMP1)*CEGRID_.DLFC;
      int KE1 = (int)XE1;
      double XEK1 = XE1-(double)KE1;
      double STLWR = exp(SETOT[MAT-1][KE1-1]+(SETOT[MAT-1][KE1+1-1]-SETOT[MAT-1][KE1-1])*XEK1);
      if(CJUMP0_.ST < STLWR){ CJUMP0_.ST = STLWR;}
    }
    else
    {
      CJUMP1_.KSOFTI = 0;
      DESOFT = 0.0;
      SSOFT = 0.0;
    }

      //  ****  Soft elastic scattering.
      //        KSOFTE=1, soft scattering is active,
      //        KSOFTE=0, soft scattering is not active.

    if(CJUMP0_.T1 > 1.0E-20)
    {
      CJUMP1_.KSOFTE = 1;
    }
    else
    {
      CJUMP1_.KSOFTE = 0;
    }

      //  ****  Delta interactions.
      //        KDELTA=0, a hard interaction follows,
      //        KDELTA=1, a delta interaction follows.

    CJUMP0_.DST = -log(RAND(2.0))/CJUMP0_.ST;
    if(CJUMP0_.DST < DSMAXP)
    {
      CJUMP1_.KDELTA = 0;
    }
    else
    {
      CJUMP0_.DST = DSMAXP;
      CJUMP1_.KDELTA = 1;
    }

    if(CJUMP1_.KSOFTE+CJUMP1_.KSOFTI == 0)
    {
      MHINGE = 1;
      DS = CJUMP0_.DST;
    }
    else
    {
      DS = CJUMP0_.DST*RAND(3.0);
      CJUMP0_.DSR = CJUMP0_.DST-DS;
      if(CJUMP1_.KSOFTI == 1)
      {
        if(CJUMP0_.DST < 1.0E-8)
        {
          SSOFT = CJUMP0_.W1;
          DESOFT = SSOFT*CJUMP0_.DST;
        }
        else
        {
          double EDE0 = CJUMP0_.W1*CJUMP0_.DST;
          double VDE0 = CJUMP0_.W2*CJUMP0_.DST;
          double FSEDE = 1.0-DW1EL[MAT-1][CEGRID_.KE-1]*EDE0;
          if(FSEDE < 0.75){ FSEDE = 0.75;}

          double FSVDE = 1.0-DW2EL[MAT-1][CEGRID_.KE-1]*EDE0;
          if(FSVDE < 0.75){ FSVDE = 0.75;}
          double EDE = EDE0*FSEDE;
          double VDE = VDE0*FSVDE;
      //  ****  Generation of random values DE with mean EDE and variance VDE.
          double SIGMA = sqrt(VDE);
          if(SIGMA < 0.333333333*EDE)
          {
          //  ****  Truncated Gaussian distribution.
            DESOFT = EDE+RNDG3()*SIGMA;
          }
          else
          {
            double RU = RAND(4.0);
            double EDE2 = EDE*EDE;
            double VDE3 = 3.0*VDE;
            if(EDE2 < VDE3)
            {
              double PNULL = (VDE3-EDE2)/(VDE3+3.0*EDE2);
              if(RU < PNULL)
              {
                DESOFT = 0.0;
                SSOFT = 0.0;
                if(CJUMP1_.KSOFTE == 0)
                {
                  MHINGE = 1;
                  DS = CJUMP0_.DST;
                }
                else
                {
                  CJUMP1_.KSOFTI = 0;
                }
                return;
              }
              else
              {
            //  ****  Uniform distribution.
                DESOFT = 1.5*(EDE+VDE/EDE)*(RU-PNULL)/(1.0-PNULL);
              }
            }
            else
            {
              DESOFT = EDE+(2.0*RU-1.0)*sqrt(VDE3);
            }
          }
          SSOFT = DESOFT/CJUMP0_.DST;
        }
      }
    }
    return;
  }
  else if(KPAR == 3)
  {

      //  ************  Positrons (KPAR=3).

    if(MHINGE == 1)
    {
      if(E < CJUMP1_.ELAST1)
      {
        CEGRID_.XEL = log(E);
        CEGRID_.XE = 1.0+(CEGRID_.XEL-CEGRID_.DLEMP1)*CEGRID_.DLFC;
        CEGRID_.KE = (int)CEGRID_.XE;
        CEGRID_.XEK = CEGRID_.XE-(double)CEGRID_.KE;
        PIMFP(1);
        CJUMP1_.ELAST1 = E;
      }
      DS = CJUMP0_.DSR;
      return;
    }

    E0STEP = E;
    if(E < CJUMP1_.ELAST2)
    {
      CEGRID_.XEL = log(E);
      CEGRID_.XE = 1.0+(CEGRID_.XEL-CEGRID_.DLEMP1)*CEGRID_.DLFC;
      CEGRID_.KE = (int)CEGRID_.XE;
      CEGRID_.XEK = CEGRID_.XE-(double)CEGRID_.KE;
      PIMFP(2);
      CJUMP1_.ELAST2 = E;
      CJUMP1_.ELAST1 = E;
    }

      //  ****  Inverse hard mean free path (interaction probability per unit
      //        path length).

    CJUMP0_.ST = CJUMP0_.P[1]+CJUMP0_.P[2]+CJUMP0_.P[3]+CJUMP0_.P[4]+CJUMP0_.P[5]+CJUMP0_.P[7];
    double DSMAXP = DSMAX;
    //  ****  The value of DSMAXP is randomised to eliminate dose artifacts
    //        at the end of the first step.
    DSMAXP = (0.75+RAND(1.0)*0.5)*DSMAXP;

      //  ****  Soft stopping interactions.
      //        KSOFTI=1, soft stopping is active,
      //        KSOFTI=0, soft stopping is not active.

    if(CJUMP0_.W1 > 1.0E-20)
    {
      CJUMP1_.KSOFTI = 1;
    //  ****  The maximum step length, DSMAXP, is determined in terms of the
    //        input DSMAX value (which is specified by the user) and the mean
    //        free path for hard interactions (1/ST).
      double DSMC = 4.0/CJUMP0_.ST;
      if(DSMAXP > DSMC)
      {
        DSMAXP = DSMC;
      }
      else if(DSMAXP < 1.0E-8)
      {
        DSMAXP = DSMC;
      }

    //
    //  ****  Upper bound for the interaction probability along the step
    //        (including soft energy straggling).

      double EDE0 = CJUMP0_.W1*DSMAXP;
      double VDE0 = CJUMP0_.W2*DSMAXP;

      double FSEDE = 1.0-DW1PL[MAT-1][CEGRID_.KE-1]*EDE0;
      if(FSEDE < 0.75){ FSEDE = 0.75;}

      double FSVDE = 1.0-DW2PL[MAT-1][CEGRID_.KE-1]*EDE0;
      if(FSVDE < 0.75){ FSVDE = 0.75;}

      double EDEM = EDE0*FSEDE;
      double VDEM = VDE0*FSVDE;
      double W21 = VDEM/EDEM;
      double ELOWER;
      if(EDEM > 9.0*W21)
      {
        ELOWER = E-(EDEM+3.0*sqrt(VDEM));
        if(ELOWER < CEGRID_.EMIN){ ELOWER = CEGRID_.EMIN;}
      }
      else if(EDEM > 3.0*W21)
      {
        ELOWER = E-(EDEM+sqrt(3.0*VDEM));
        if(ELOWER < CEGRID_.EMIN){ ELOWER = CEGRID_.EMIN;}
      }
      else
      {
        ELOWER = E-1.5*(EDEM+W21);
        if( ELOWER < CEGRID_.EMIN){ ELOWER = CEGRID_.EMIN;}
      }
      double XE1 = 1.0+(log(ELOWER)-CEGRID_.DLEMP1)*CEGRID_.DLFC;
      int KE1 = (int)XE1;
      double XEK1 = XE1-(double)KE1;
      double STLWR = exp(SPTOT[MAT-1][KE1-1]+(SPTOT[MAT-1][KE1+1-1]-SPTOT[MAT-1][KE1-1])*XEK1);
      if(CJUMP0_.ST < STLWR){ CJUMP0_.ST = STLWR;}
    }
    else
    {
      CJUMP1_.KSOFTI = 0;
      DESOFT = 0.0;
      SSOFT = 0.0;
    }

      //  ****  Soft elastic scattering.
      //        KSOFTE=1, soft scattering is active,
      //        KSOFTE=0, soft scattering is not active.
      
    if(CJUMP0_.T1 > 1.0E-20)
    {
      CJUMP1_.KSOFTE = 1;
    }
    else
    {
      CJUMP1_.KSOFTE = 0;
    }
      
      //  ****  Delta interactions.
      //        KDELTA=0, a hard interaction follows,
      //        KDELTA=1, a delta interaction follows.
      
    CJUMP0_.DST = -log(RAND(2.0))/CJUMP0_.ST;
    if(CJUMP0_.DST < DSMAXP)
    {
      CJUMP1_.KDELTA = 0;
    }
    else
    {
      CJUMP0_.DST = DSMAXP;
      CJUMP1_.KDELTA = 1;
    }

    if(CJUMP1_.KSOFTE+CJUMP1_.KSOFTI == 0)
    {
      MHINGE = 1;
      DS = CJUMP0_.DST;
    }
    else
    {
      DS = CJUMP0_.DST*RAND(3.0);
      CJUMP0_.DSR = CJUMP0_.DST-DS;
      if(CJUMP1_.KSOFTI == 1)
      {
        if(CJUMP0_.DST < 1.0E-8)
        {
          SSOFT = CJUMP0_.W1;
          DESOFT = SSOFT*CJUMP0_.DST;
        }
        else
        {
          double EDE0 = CJUMP0_.W1*CJUMP0_.DST;
          double VDE0 = CJUMP0_.W2*CJUMP0_.DST;
          double FSEDE = 1.0-DW1PL[MAT-1][CEGRID_.KE-1]*EDE0;
          if(FSEDE < 0.75){ FSEDE = 0.75;}

          double FSVDE = 1.0-DW2PL[MAT-1][CEGRID_.KE-1]*EDE0;
          if(FSVDE < 0.75){ FSVDE = 0.75;}
      
          double EDE = EDE0*FSEDE;
          double VDE = VDE0*FSVDE;
      //  ****  Generation of random values DE with mean EDE and variance VDE.
          double SIGMA = sqrt(VDE);
          if(SIGMA < 0.333333333E0*EDE)
          {
          //  ****  Truncated Gaussian distribution.
            DESOFT = EDE+RNDG3()*SIGMA;
          }
          else
          {
            double RU = RAND(4.0);
            double EDE2 = EDE*EDE;
            double VDE3 = 3.0*VDE;
            if(EDE2 < VDE3)
            {
              double PNULL = (VDE3-EDE2)/(VDE3+3.0*EDE2);
              if(RU < PNULL)
              {
                DESOFT = 0.0;
                SSOFT = 0.0;
                if(CJUMP1_.KSOFTE == 0)
                {
                  MHINGE = 1;
                  DS = CJUMP0_.DST;
                }
                else
                {
                  CJUMP1_.KSOFTI = 0;
                }
                return;
              }
              else
              {
            //  ****  Uniform distribution.
                DESOFT = 1.5*(EDE+VDE/EDE)*(RU-PNULL)/(1.0-PNULL);
              }
            }
            else
            {
              DESOFT = EDE+(2.0*RU-1.0)*sqrt(VDE3);
            }
          }
          SSOFT = DESOFT/CJUMP0_.DST;
        }
      }
    }
    return;
  }
  else
  {
      //
      //  ************  Photons (KPAR=2).
      //
    if(E < CJUMP1_.ELAST1)
    {
      CEGRID_.XEL = log(E);
      CEGRID_.XE = 1.0+(CEGRID_.XEL-CEGRID_.DLEMP1)*CEGRID_.DLFC;
      CEGRID_.KE = (int)CEGRID_.XE;
      CEGRID_.XEK = CEGRID_.XE-(double)CEGRID_.KE;
      GIMFP();
      CJUMP1_.ELAST1 = E;
      CJUMP0_.ST = CJUMP0_.P[0]+CJUMP0_.P[1]+CJUMP0_.P[2]+CJUMP0_.P[3]+CJUMP0_.P[7];
    }

    DS = -log(RAND(1.0))/CJUMP0_.ST;
  }
  return;
}
