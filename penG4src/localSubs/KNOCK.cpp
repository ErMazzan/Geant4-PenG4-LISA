//  *********************************************************************
//                       SUBROUTINE KNOCK
//  *********************************************************************
void PenPhys::KNOCK(double &DE, int &ICOL)
{
  //  Simulation of random hinges and hard interaction events.

  //  Output arguments:
  //    DE ..... energy deposited by the particle in the material. It is
  //             usually equal to the difference between the energies
  //             before and after the interaction.
  //    ICOL ... kind of interaction suffered by the particle.

  using namespace PENELOPE_mod;
  using namespace PENERROR_mod;
///////////// using namespace TRACK_mod;

/////////////  using namespace CHIST;
  using namespace COMPOS;
/////////////  using namespace CEGRID;
  using namespace CADATA;
  using namespace CEIN;
  using namespace CGCO;
  using namespace CEBR;
  using namespace CEIMFP;
  using namespace CPIMFP;
  using namespace CELSEP;
/////////////  using namespace CJUMP0;
/////////////  using namespace CJUMP1;
  
  const double TWOPI = PI+PI;
  const double RREV = 1.0/REV;
  const double TREV = 2.0*REV;
  
  if(KPAR == 1)
  {

    //  ************  Electrons (KPAR=1).
    if(MHINGE == 1){}
    else
    {
    
      //  ****  Hinge, artificial soft event (ICOL=1).
    
      ICOL = 1;
      MHINGE = 1;
      //  ****  Energy loss.

      if(CJUMP1_.KSOFTI == 1)
      {
        DE = DESOFT;
        E = E-DE;
        if(E < EABS[0][MAT-1])
        {
          DE = E0STEP;
          E = 0.0;
          return;
        }
        E0STEP = E0STEP-SSOFT*(CJUMP0_.DST-CJUMP0_.DSR);
        if(CJUMP1_.KSOFTE == 0){ return;}
        CEGRID_.XEL = log(E0STEP);
        CEGRID_.XE = 1.0+(CEGRID_.XEL-CEGRID_.DLEMP1)*CEGRID_.DLFC;
        CEGRID_.KE = (int)CEGRID_.XE;
        CEGRID_.XEK = CEGRID_.XE-(double)CEGRID_.KE;
      }
      else
      {
        DE = 0.0;
      }

      //  ****  Angular deflection.

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
      if(CJUMP0_.T1 < 1.0E-20){ return;}
      //  ****  1st and 2nd moments of the angular distribution.
      double EMU1 = 0.5*(1.0-exp(-CJUMP0_.DST*CJUMP0_.T1));
      double EMU2 = EMU1-(1.0-exp(-CJUMP0_.DST*CJUMP0_.T2))/6.0;
      //  ****  Sampling from a two-bar histogram with these moments.
      double PNUM = 2.0*EMU1-3.0*EMU2;
      double PDEN = 1.0-2.0*EMU1;
      double PMU0 = PNUM/PDEN;
      double PA = PDEN+PMU0;
      double RND = RAND(2.0);
      double CDT;
      if(RND < PA)
      {
        CDT = 1.0-2.0*PMU0*(RND/PA);
      }
      else
      {
        CDT = 1.0-2.0*(PMU0+(1.0-PMU0)*((RND-PA)/(1.0-PA)));
      }
      double DF = TWOPI*RAND(3.0);
      DIRECT(CDT,DF,U,V,W);
      return;
    }

    //  ************  Hard event.

    MHINGE = 0;
    //  ****  A delta interaction (ICOL=7) occurs when the maximum
    //        allowed step length is exceeded.
    if(CJUMP1_.KDELTA == 1)
    {
      ICOL = 7;
      DE = 0.0;
      return;
    }
    //  ****  Random sampling of the interaction type.
    double STNOW = CJUMP0_.P[1]+CJUMP0_.P[2]+CJUMP0_.P[3]+CJUMP0_.P[4]+CJUMP0_.P[7];
    double STS = (STNOW > CJUMP0_.ST ? STNOW : CJUMP0_.ST)*RAND(4.0);
    double SS = CJUMP0_.P[1];
    if(SS > STS)   //GOTO 1200
    {       
      //  ****  Hard elastic collision (ICOL=2).
      double RMU;  
      ICOL = 2;
      if(E >= EELMAX[MAT-1])
      {
        double TRNDC = RNDCE[MAT-1][CEGRID_.KE-1]+(RNDCE[MAT-1][CEGRID_.KE+1-1]-RNDCE[MAT-1][CEGRID_.KE-1])*CEGRID_.XEK;
        double TA = exp(AE[MAT-1][CEGRID_.KE-1]+(AE[MAT-1][CEGRID_.KE+1-1]-AE[MAT-1][CEGRID_.KE-1])*CEGRID_.XEK);
        double TB = BE[MAT-1][CEGRID_.KE-1]+(BE[MAT-1][CEGRID_.KE+1-1]-BE[MAT-1][CEGRID_.KE-1])*CEGRID_.XEK;
        EELa(TA,TB,TRNDC,RMU);
      }
      else
      {

        double TRNDC = RNDCEd[MAT-1][CEGRID_.KE-1]+(RNDCEd[MAT-1][CEGRID_.KE+1-1]-RNDCEd[MAT-1][CEGRID_.KE-1])*CEGRID_.XEK;
        EELd(TRNDC,RMU);  // Uses the ELSEPA database.
      }
      double CDT = 1.0-(RMU+RMU);
      double DF = TWOPI*RAND(5.0);
      DIRECT(CDT,DF,U,V,W);
      DE = 0.0;
      return;   
    }
    SS = SS+CJUMP0_.P[2];   //GOTO 1300   
    if(SS > STS)
    {

      //  ****  Hard inelastic collision (ICOL=3).
      double EP,CDT,ES,CDTS;
      int IOSC;
      ICOL = 3;
      double DELTA = DEL[MAT-1][CEGRID_.KE-1]+(DEL[MAT-1][CEGRID_.KE+1-1]-DEL[MAT-1][CEGRID_.KE-1])*CEGRID_.XEK;
      EINa(E,DELTA,DE,EP,CDT,ES,CDTS,MAT,IOSC);
      //  ****  Scattering angles (primary electron).
      double DF = TWOPI*RAND(6.0);
      //  ****  Delta ray.
      if(ES > EABS[0][MAT-1])
      {
        double DFS = DF+PI;
        double US = U;
        double VS = V;
        double WS = W;
        DIRECT(CDTS,DFS,US,VS,WS);
        CHIST_.ILBA[0] = ILB[0]+1;
        CHIST_.ILBA[1] = KPAR;
        CHIST_.ILBA[2] = ICOL;
        CHIST_.ILBA[3] = 0;
        CHIST_.ILBA[4] = ILB[4];
        STORES(ES,X,Y,Z,US,VS,WS,WGHT,1,CHIST_.ILBA,0);
        if(IRETRN != 0){ return;}        
      }
      //  ****  New energy and direction.
      if(EP > EABS[0][MAT-1])
      {
        E = EP;
        DIRECT(CDT,DF,U,V,W);
      }
      else
      {
        DE = E;
        E = 0.0;
      }
      return;   
    }
    SS = SS+CJUMP0_.P[3];   //GOTO 1400
    if(SS > STS)
    {   
      //  ****  Hard bremsstrahlung emission (ICOL=4).
    
      ICOL = 4;
      EBRa(E,DE,MAT);
      //  ****  Bremsstrahlung photon.
      if(DE > EABS[1][MAT-1])
      {
        double CDTS;
        EBRaA(E,DE,CDTS,MAT);
        double DFS = TWOPI*RAND(7.0);
        double US = U;
        double VS = V;
        double WS = W;
        DIRECT(CDTS,DFS,US,VS,WS);
        CHIST_.ILBA[0] = ILB[0]+1;
        CHIST_.ILBA[1] = KPAR;
        CHIST_.ILBA[2] = ICOL;
        CHIST_.ILBA[3] = 0;
        CHIST_.ILBA[4] = ILB[4];
        STORES(DE,X,Y,Z,US,VS,WS,WGHT,2,CHIST_.ILBA,0);
        if(IRETRN != 0){ return;}
      }
      //  ****  New energy.
      E = E-DE;
      if(E < EABS[0][MAT-1])
      {
        DE = E+DE;
        E = 0.0;
      }
      return;
    }
    SS = SS+CJUMP0_.P[4];   //GOTO 1500
    if(SS > STS)
    {
     
      //  ****  Ionisation of an inner shell (ICOL=5).
      double EP,CDT,ES,CDTS;
      int IZA,ISA;
      ICOL = 5;
      double DELTA = DEL[MAT-1][CEGRID_.KE-1]+(DEL[MAT-1][CEGRID_.KE+1-1]-DEL[MAT-1][CEGRID_.KE-1])*CEGRID_.XEK;
      ESIa(E,DELTA,DE,EP,CDT,ES,CDTS,MAT,IZA,ISA);
      //  ****  Atomic relaxation.
      if(IZA > 2)
      {
        CHIST_.ILBA[2] = ICOL;
        RELAX(IZA,ISA);
        if(IRETRN != 0){ return;}
      }
      //  ****  Scattering angles (primary electron).
      double DF = TWOPI*RAND(8.0);
      //  ****  Delta ray.
      if(ES > EABS[0][MAT-1])
      {
        double DFS = DF+PI;
        double US = U;
        double VS = V;
        double WS = W;
        DIRECT(CDTS,DFS,US,VS,WS);
        CHIST_.ILBA[0] = ILB[0]+1;
        CHIST_.ILBA[1] = KPAR;
        CHIST_.ILBA[2] = ICOL;
        CHIST_.ILBA[3] = 0;
        CHIST_.ILBA[4] = ILB[4];
        STORES(ES,X,Y,Z,US,VS,WS,WGHT,1,CHIST_.ILBA,0);
        if(IRETRN != 0){ return;}
      }
      //  ****  New energy and direction.
      if(EP > EABS[0][MAT-1])
      {
        E = EP;
        DIRECT(CDT,DF,U,V,W);
      }
      else
      {
        DE = E;
        E = 0.0;
      }
      return;
    }
    SS = SS+CJUMP0_.P[7];     //GOTO 1800
    if(SS > STS)
    {   
      //  ****  Auxiliary fictitious mechanism (ICOL=8).
      
      ICOL = 8;
      DE = 0.0;
      EAUX();
      return;
    }
    //  ****  A delta interaction (ICOL=7) may occur when the total
    //        interaction probability per unit path length, ST, is
    //        larger than STNOW.
    ICOL = 7;
    DE = 0.0;
    return;
  }
  else if(KPAR == 2)
  {

    //  ************  Photons (KPAR=2).

    double STS = CJUMP0_.ST*RAND(1.0);
    double SS = CJUMP0_.P[0];
    if(SS > STS)   // GOTO 2100
    {
    //  ****  Rayleigh scattering (ICOL=1).
      double CDT;
      int IEFF; 
      DE = 0.0;
      GRAa(E,CDT,IEFF,MAT);
      //  ****  Delta interaction. Introduced to correct for the use of an
      //        upper bound of the Rayleigh attenuation coefficient.
      if(IEFF == 0)
      {
        ICOL = 7;
        return;
      }
      ICOL = 1;
      double DF;
      if(IPOL == 1)
      {
        DIRPOL(CDT,DF,0.0,SP1,SP2,SP3,U,V,W);
      }
      else
      {
        DF = TWOPI*RAND(2.0);
        DIRECT(CDT,DF,U,V,W);
      }
      ILB[0] = ILB[0]+1;
      ILB[1] = KPAR;
      ILB[2] = ICOL;
      return;
    }
    SS = SS+CJUMP0_.P[1];
    if(SS > STS)   //  GOTO 2200
    {
      //  ****  Compton scattering (ICOL=2).
    
      ICOL = 2;
      double EP,CDT,ES,CDTS;
      int IZA,ISA;
      GCOa(E,DE,EP,CDT,ES,CDTS,MAT,IZA,ISA);
      double US = U;
      double VS = V;
      double WS = W;
      double DF = -1.0;
      if(IZA > 0 && ISA < 17)
      {
        CHIST_.ILBA[2] = ICOL;
        RELAX(IZA,ISA);
        if(IRETRN != 0){ return;}
      }
      //  ****  New direction and energy.
      if(EP > EABS[1][MAT-1])
      {
        if(IPOL == 1)
        {
          double ECDT = E*RREV*(1.0-CDT);
          double CONS = ECDT*ECDT/(1.0+ECDT);
          DIRPOL(CDT,DF,CONS,SP1,SP2,SP3,U,V,W);
        }
        else
        {
          DF = TWOPI*RAND(3.0);
          DIRECT(CDT,DF,U,V,W);
        }
        E = EP;
      }
      else
      {
        DE = E;
        E = 0.0;
      }
      //  ****  Compton electron.
      if(ES > EABS[0][MAT-1])
      {
        if(DF < -0.5){ DF = TWOPI*RAND(4.0);}
        double DFS = DF+PI;
        DIRECT(CDTS,DFS,US,VS,WS);
        CHIST_.ILBA[0] = ILB[0]+1;
        CHIST_.ILBA[1] = KPAR;
        CHIST_.ILBA[2] = ICOL;
        CHIST_.ILBA[3] = 0;
        CHIST_.ILBA[4] = ILB[4];
        STORES(ES,X,Y,Z,US,VS,WS,WGHT,1,CHIST_.ILBA,0);
        if(IRETRN != 0){ return;}
      }
      ILB[0] = ILB[0]+1;
      ILB[1] = KPAR;
      ILB[2] = ICOL;
      return;
    }
    SS = SS+CJUMP0_.P[2];
    if(SS > STS)  // GOTO 2300
    {

      //  ****  Photoelectric absorption (ICOL=3).

      ICOL = 3;
      double ES;
      int IZA, ISA;
      GPHa(ES,IZA,ISA);
    //  ****  Delta interaction. Introduced to correct for the use of an
    //        upper bound of the photoelectric attenuation coefficient.
      if(IZA == 0)
      {
        ICOL = 7;
        DE = 0.0;
        return;
      }
    
      if(ES > EABS[0][MAT-1])
      {
        double CDTS;
        SAUTER(ES,CDTS);
        double DFS = TWOPI*RAND(5.0);
        double US = U;
        double VS = V;
        double WS = W;
        DIRECT(CDTS,DFS,US,VS,WS);
        CHIST_.ILBA[0] = ILB[0]+1;
        CHIST_.ILBA[1] = KPAR;
        CHIST_.ILBA[2] = ICOL;
        CHIST_.ILBA[3] = 0;
        CHIST_.ILBA[4] = ILB[4];
        STORES(ES,X,Y,Z,US,VS,WS,WGHT,1,CHIST_.ILBA,0);
        if(IRETRN != 0){ return;}
      }
      if(ISA < 17)
      {
        CHIST_.ILBA[2] = ICOL;
        RELAX(IZA,ISA);
        if(IRETRN != 0){ return;}
      }
      DE = E;
      E = 0.0;
      return;
    }
    SS = SS+CJUMP0_.P[3];
    if(SS > STS)  //  GOTO 2400
    {

      //  ****  Electron-positron pair production (ICOL=4).

      ICOL = 4;
      double EE,CDTE,EP,CDTP;
      int IZA,ISA;
      GPPa(EE,CDTE,EP,CDTP,IZA,ISA);
      DE = E;
      //  ****  Electron.
      if(EE > EABS[0][MAT-1])
      {
        double DF = TWOPI*RAND(6.0);
        double US = U;
        double VS = V;
        double WS = W;
        DIRECT(CDTE,DF,US,VS,WS);
        CHIST_.ILBA[0] = ILB[0]+1;
        CHIST_.ILBA[1] = KPAR;
        CHIST_.ILBA[2] = ICOL;
        CHIST_.ILBA[3] = 0;
        CHIST_.ILBA[4] = ILB[4];
        STORES(EE,X,Y,Z,US,VS,WS,WGHT,1,CHIST_.ILBA,0);
        if(IRETRN != 0){ return;}
      }
      //  ****  Positron.
      if(EP > EABS[2][MAT-1])
      {
        double DF = TWOPI*RAND(7.0);
        double US = U;
        double VS = V;
        double WS = W;
        DIRECT(CDTP,DF,US,VS,WS);
        CHIST_.ILBA[0] = ILB[0]+1;
        CHIST_.ILBA[1] = KPAR;
        CHIST_.ILBA[2] = ICOL;
        CHIST_.ILBA[3] = 0;
        CHIST_.ILBA[4] = ILB[4];
        STORES(EP,X,Y,Z,US,VS,WS,WGHT,3,CHIST_.ILBA,0);
        if(IRETRN != 0){ return;}
        //  ****  The positron carries a 'latent' energy of 1022 keV.
        DE = DE-TREV;
      }
      else
      {
        PANaR(EABS[1][MAT-1]);
        if(IRETRN != 0){ return;}
      }
      E = 0.0;
      //  ****  Atomic relaxation after triplet production.
      if(ISA < 17)
      {
        CHIST_.ILBA[2] = ICOL;
        RELAX(IZA,ISA);
        if(IRETRN != 0){ return;}
      }
      return;
    }
    SS = SS+CJUMP0_.P[7];
    if(SS > STS)
    {

      //  ****  Auxiliary fictitious mechanism (ICOL=8).

      ICOL = 8;
      DE = 0.0;
      GAUX();
      return;
    }
  }
  else if(KPAR == 3)
  {

    //  ************  Positrons (KPAR=3).
    if(MHINGE == 1){}
    else
    {
  
      //  ****  Hinge, artificial soft event (ICOL=1).
  
      ICOL = 1;
      MHINGE = 1;
      //  ****  Energy loss.
  
      if(CJUMP1_.KSOFTI == 1)
      {
        DE = DESOFT;
        E = E-DE;
        if(E < EABS[2][MAT-1])
        {
          PANaR(EABS[1][MAT-1]);  // Annihilation at rest.
          if(IRETRN != 0){ return;}
          DE = E0STEP+TREV;
          E = 0.0;
          return;
        }
        E0STEP = E0STEP-SSOFT*(CJUMP0_.DST-CJUMP0_.DSR);
        if(CJUMP1_.KSOFTE == 0){ return;}
        CEGRID_.XEL = log(E0STEP);
        CEGRID_.XE = 1.0+(CEGRID_.XEL-CEGRID_.DLEMP1)*CEGRID_.DLFC;
        CEGRID_.KE = (int)CEGRID_.XE;
        CEGRID_.XEK = CEGRID_.XE-(double)CEGRID_.KE;
      }
      else
      {
        DE = 0.0;
      }

      //  ****  Angular deflection.

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
      if(CJUMP0_.T1 < 1.0E-20){ return;}
      //  ****  1st and 2nd moments of the angular distribution.
      double EMU1 = 0.5*(1.0-exp(-CJUMP0_.DST*CJUMP0_.T1));
      double EMU2 = EMU1-(1.0-exp(-CJUMP0_.DST*CJUMP0_.T2))/6.0;
      //  ****  Sampling from a two-bar histogram with these moments.
      double PNUM = 2.0*EMU1-3.0*EMU2;
      double PDEN = 1.0-2.0*EMU1;
      double PMU0 = PNUM/PDEN;
      double PA = PDEN+PMU0;
      double RND = RAND(2.0);
      double CDT;
      if(RND < PA)
      {
        CDT = 1.0-2.0*PMU0*(RND/PA);
      }
      else
      {
        CDT = 1.0-2.0*(PMU0+(1.0-PMU0)*((RND-PA)/(1.0-PA)));
      }
      double DF = TWOPI*RAND(3.0);
      DIRECT(CDT,DF,U,V,W);
      return;
    }

    //  ************  Hard event.

    MHINGE = 0;
    //  ****  A delta interaction (ICOL=7) occurs when the maximum
    //        allowed step length is exceeded.
    if(CJUMP1_.KDELTA == 1)
    {
      ICOL = 7;
      DE = 0.0;
      return;
    }
    //  ****  Random sampling of the interaction type.
    double STNOW = CJUMP0_.P[1]+CJUMP0_.P[2]+CJUMP0_.P[3]+CJUMP0_.P[4]+CJUMP0_.P[5]+CJUMP0_.P[7];
    double STS = (STNOW > CJUMP0_.ST ? STNOW : CJUMP0_.ST)*RAND(4.0);
    double SS = CJUMP0_.P[1];
    if(SS > STS)  // GOTO 3200
    {

      //  ****  Hard elastic collision (ICOL=2).
      double RMU;
      ICOL = 2;
      if(E >= PELMAX[MAT-1])
      {
        double TRNDC = RNDCP[MAT-1][CEGRID_.KE-1]+(RNDCP[MAT-1][CEGRID_.KE+1-1]-RNDCP[MAT-1][CEGRID_.KE-1])*CEGRID_.XEK;
        double TA = exp(AP[MAT-1][CEGRID_.KE-1]+(AP[MAT-1][CEGRID_.KE+1-1]-AP[MAT-1][CEGRID_.KE-1])*CEGRID_.XEK);
        double TB = BP[MAT-1][CEGRID_.KE-1]+(BP[MAT-1][CEGRID_.KE+1-1]-BP[MAT-1][CEGRID_.KE-1])*CEGRID_.XEK;
        EELa(TA,TB,TRNDC,RMU);
      }
      else
      {
        double TRNDC = RNDCPd[MAT-1][CEGRID_.KE-1]+(RNDCPd[MAT-1][CEGRID_.KE+1-1]-RNDCPd[MAT-1][CEGRID_.KE-1])*CEGRID_.XEK;
        PELd(TRNDC,RMU);  // Uses the ELSEPA database.
      }
      double CDT = 1.0-(RMU+RMU);
      double DF = TWOPI*RAND(5.0);
      DIRECT(CDT,DF,U,V,W);
      DE = 0.0;
      return;
    }
    SS = SS+CJUMP0_.P[2];
    if(SS > STS)  // GOTO 3300
    {

      //  ****  Hard inelastic collision (ICOL=3).
      double EP, CDT, ES, CDTS;
      int IOSC;
      ICOL = 3;
      double DELTA = DEL[MAT-1][CEGRID_.KE-1]+(DEL[MAT-1][CEGRID_.KE+1-1]-DEL[MAT-1][CEGRID_.KE-1])*CEGRID_.XEK;
      PINa(E,DELTA,DE,EP,CDT,ES,CDTS,MAT,IOSC);
      //  ****  Scattering angles (primary positron).
      double DF = TWOPI*RAND(6.0);
      //  ****  Delta ray.
      if(ES > EABS[0][MAT-1])
      {
        double DFS = DF+PI;
        double US = U;
        double VS = V;
        double WS = W;
        DIRECT(CDTS,DFS,US,VS,WS);
        CHIST_.ILBA[0] = ILB[0]+1;
        CHIST_.ILBA[1] = KPAR;
        CHIST_.ILBA[2] = ICOL;
        CHIST_.ILBA[3] = 0;
        CHIST_.ILBA[4] = ILB[4];
        STORES(ES,X,Y,Z,US,VS,WS,WGHT,1,CHIST_.ILBA,0);
        if(IRETRN != 0){ return;}
      }
      //  ****  New energy and direction.
      if(EP > EABS[2][MAT-1])
      {
        E = EP;
        DIRECT(CDT,DF,U,V,W);
      }
      else
      {
        PANaR(EABS[1][MAT-1]);  // Annihilation at rest.
        if(IRETRN != 0){ return;}
        DE = E+TREV;
        E = 0.0;
      }
      return;
    }
    SS = SS+CJUMP0_.P[3];
    if(SS > STS)    // GOTO 3400
    {

      //  ****  Hard bremsstrahlung emission (ICOL=4).

      ICOL = 4;
      EBRa(E,DE,MAT);
      //  ****  Bremsstrahlung photon.
      if(DE > EABS[1][MAT-1])
      {
        double CDTS;
        EBRaA(E,DE,CDTS,MAT);
        double DFS = TWOPI*RAND(7.0);
        double US = U;
        double VS = V;
        double WS = W;
        DIRECT(CDTS,DFS,US,VS,WS);
        CHIST_.ILBA[0] = ILB[0]+1;
        CHIST_.ILBA[1] = KPAR;
        CHIST_.ILBA[2] = ICOL;
        CHIST_.ILBA[3] = 0;
        CHIST_.ILBA[4] = ILB[4];
        STORES(DE,X,Y,Z,US,VS,WS,WGHT,2,CHIST_.ILBA,0);
        if(IRETRN != 0){ return;}
      }
      //  ****  New energy.
      E = E-DE;
      if(E < EABS[2][MAT-1])
      {
        PANaR(EABS[1][MAT-1]);  // Annihilation at rest.
        if(IRETRN != 0){ return;}
        DE = E+DE+TREV;
        E = 0.0;
      }
      return;
    }
    SS = SS+CJUMP0_.P[4];
    if(SS > STS)    // GOTO 3500
    {

      //  ****  Ionisation of an inner shell (ICOL=5).
      double EP,CDT,ES,CDTS;
      int IZA,ISA;
      ICOL = 5;
      double DELTA = DEL[MAT-1][CEGRID_.KE-1]+(DEL[MAT-1][CEGRID_.KE+1-1]-DEL[MAT-1][CEGRID_.KE-1])*CEGRID_.XEK;
      PSIa(E,DELTA,DE,EP,CDT,ES,CDTS,MAT,IZA,ISA);
      //  ****  Atomic relaxation.
      if(IZA > 2)
      {
        CHIST_.ILBA[2] = ICOL;
        RELAX(IZA,ISA);
        if(IRETRN != 0){ return;}
      }
      //  ****  Scattering angles (primary electron).
      double DF = TWOPI*RAND(8.0);
      //  ****  Delta ray.
      if(ES > EABS[0][MAT-1])
      {
        double DFS = DF+PI;
        double US = U;
        double VS = V;
        double WS = W;
        DIRECT(CDTS,DFS,US,VS,WS);
        CHIST_.ILBA[0] = ILB[0]+1;
        CHIST_.ILBA[1] = KPAR;
        CHIST_.ILBA[2] = ICOL;
        CHIST_.ILBA[3] = 0;
        CHIST_.ILBA[4] = ILB[4];
        STORES(ES,X,Y,Z,US,VS,WS,WGHT,1,CHIST_.ILBA,0);
        if(IRETRN != 0){ return;}
      }
      //  ****  New energy and direction.
      if(EP > EABS[2][MAT-1])
      {
        E = EP;
        DIRECT(CDT,DF,U,V,W);
      }
      else
      {
        PANaR(EABS[1][MAT-1]);  // Annihilation at rest.
        if(IRETRN != 0){ return;}
        DE = E+TREV;
        E = 0.0;
      }
      return;
    }
    SS = SS+CJUMP0_.P[5];
    if(SS > STS)  // GOTO 3600
    {
      //  ****  Positron annihilation in flight (ICOL=6).
      double E1,CDT1,E2,CDT2;
      ICOL = 6;
      PANa(E,E1,CDT1,E2,CDT2,MAT);
      double DF = TWOPI*RAND(9.0);
      if(E1 > EABS[1][MAT-1])
      {
        double US = U;
        double VS = V;
        double WS = W;
        DIRECT(CDT1,DF,US,VS,WS);
        CHIST_.ILBA[0] = ILB[0]+1;
        CHIST_.ILBA[1] = KPAR;
        CHIST_.ILBA[2] = ICOL;
        CHIST_.ILBA[3] = 0;
        CHIST_.ILBA[4] = ILB[4];
        STORES(E1,X,Y,Z,US,VS,WS,WGHT,2,CHIST_.ILBA,0);
        if(IRETRN != 0){ return;}
      }
      if(E2 > EABS[1][MAT-1])
      {
        DF = DF+PI;
        double US = U;
        double VS = V;
        double WS = W;
        DIRECT(CDT2,DF,US,VS,WS);
        CHIST_.ILBA[0] = ILB[0]+1;
        CHIST_.ILBA[1] = KPAR;
        CHIST_.ILBA[2] = ICOL;
        CHIST_.ILBA[3] = 0;
        CHIST_.ILBA[4] = ILB[4];
        STORES(E2,X,Y,Z,US,VS,WS,WGHT,2,CHIST_.ILBA,0);
        if(IRETRN != 0){ return;}
      }
      DE = E+TREV;
      E = 0.0;
      return;
    }
    SS = SS+CJUMP0_.P[7];
    if(SS > STS)  // GOTO 3800
    {
      //  ****  Auxiliary fictitious mechanism (ICOL=8).
    
      ICOL = 8;
      DE = 0.0;
      PAUX();
      return;
    }
      //  ****  A delta interaction (ICOL=7) may occur when the total
      //        interaction probability per unit path length, ST, is
      //        larger than STNOW.
    ICOL = 7;
    DE = 0.0;
    return;
  }
  else
  {
    ErrorFunction(1202); return;
  }

}
