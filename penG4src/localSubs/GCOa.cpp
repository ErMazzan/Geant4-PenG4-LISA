
//  *********************************************************************
//                       SUBROUTINE GCOa
//  *********************************************************************
void PenPhys::GCOa(double E_, double &DE, double &EP, double &CDT, double &ES, double &CDTS, int M, int &IZZ,int &ISH)
{
  //  Random sampling of incoherent (Compton) scattering of photons. Relat-
  //  ivistic impulse approximation with analytical one-electron Compton
  //  profiles.

  //  Input arguments:
  //    E_ ..... incident photon energy (eV).
  //    M ..... material where photons propagate.
  //  Output argument:
  //    DE .... energy loss (eV).
  //    EP .... energy of the scattered photon (eV).
  //    CDT ... cosine of the polar scattering angle.
  //    ES .... energy of the emitted electron (eV).
  //    CDTS .. polar cosine of direction of the electron.
  //    IZZ ... atomic number of the atom where scattering has occurred.
  //    ISH ... atomic electron shell that has been ionized.

  using namespace PENELOPE_mod;

  using namespace COMPOS;
  using namespace CECUTR;
  using namespace CGCO;

  const double RREV = 1.0/REV;
  const double D2 = 1.4142135623731;
  const double D1 = 1.0/D2;
  const double D12 = 0.5;
  
  double RN[NOCO], PAC[NOCO];
  double TAU, TST;
  int    ISHELL;
  
  double EK = E_*RREV;
  double EK2 = EK+EK+1.0;
  double EKS = EK*EK;
  double EK1 = EKS-EK2-1.0;
  double TAUMIN = 1.0/EK2;
  double TAUM2 = TAUMIN*TAUMIN;
  double A1 = log(EK2);
  double A2 = A1+2.0*EK*(1.0+EK)*TAUM2;
  bool Egt5MeV = false;
  if(E_ > 5.0E6){ Egt5MeV = true;}
  else
  {

    //  ****  Incoherent scattering function for theta=PI.
      
    double S0 = 0.0;
    double AUX, PZOMC, RNI;
    for(int I = 0; I < NOSCCO[M-1]; I++)
    {
      if(UICO[M-1][I] < E_)
      {
        AUX = E_*(E_-UICO[M-1][I])*2.0;
        PZOMC = FJ0[M-1][I]*(AUX-REV*UICO[M-1][I])/(REV*sqrt(AUX+AUX+(UICO[M-1][I]*UICO[M-1][I])));
        if(PZOMC > 0.0)
        {
          RNI = 1.0-0.5*exp(D12-((D1+D2*PZOMC)*(D1+D2*PZOMC)));
        }
        else
        {
          RNI = 0.5*exp(D12-((D1-D2*PZOMC)*(D1-D2*PZOMC)));
        }
        S0 = S0+FCO[M-1][I]*RNI;
      }
    }
      
    //  ****  Sampling tau.
      
    double CDT1, S, A;
    bool brkIt = false;
    while(!brkIt)
    {
      brkIt = true;
      if(RAND(1.0)*A2 < A1)
      {
        TAU = pow(TAUMIN,RAND(2.0));
      }
      else
      {
        TAU = sqrt(1.0+RAND(3.0)*(TAUM2-1.0));
      }
      CDT1 = (1.0-TAU)/(EK*TAU);
      //  ****  Incoherent scattering function.
      S = 0.0;
      for(int I = 0; I < NOSCCO[M-1]; I++)
      {
        if(UICO[M-1][I] < E_)
        {
          AUX = E_*(E_-UICO[M-1][I])*CDT1;
          PZOMC = FJ0[M-1][I]*(AUX-REV*UICO[M-1][I])/(REV*sqrt(AUX+AUX+(UICO[M-1][I]*UICO[M-1][I])));
          if(PZOMC > 0.0)
          {
            RN[I] = 1.0-0.5*exp(D12-((D1+D2*PZOMC)*(D1+D2*PZOMC)));
          }
          else
          {
            RN[I] = 0.5*exp(D12-((D1-D2*PZOMC)*(D1-D2*PZOMC)));
          }
          S = S+FCO[M-1][I]*RN[I];
          PAC[I] = S;
        }
        else
        {
          PAC[I] = S;
        }
      }
      //  ****  Rejection function.
      TST = S*(1.0+TAU*(EK1+TAU*(EK2+TAU*EKS)))/(EKS*TAU*(1.0+TAU*TAU));
      if(RAND(4.0)*S0 > TST){ brkIt = false; continue;}
    }
    CDT = 1.0-CDT1;
      
    //  ****  Target electron shell.
      
    int JO;
    double XQC, FPZMAX, AF, FPZ, T, B1, B2;
    brkIt = false;
    while(!brkIt)
    {
      brkIt = true;
      TST = S*RAND(5.0);
      //  ****  Binary search.
      if(TST < PAC[0])
      {
        ISHELL = 1;
      }
      else
      {
        ISHELL = 1;
        JO = NOSCCO[M-1]+1;
        int I;
        bool brkIt2 = false;
        while(!brkIt2)
        {
          brkIt2 = true;
          I = (ISHELL+JO)/2;
          if(TST > PAC[I-1])
          {
            ISHELL = I;
          }
          else
          {
            JO = I;
          }
          if(JO-ISHELL > 1){ brkIt2 = false; continue;}
        }
        ISHELL = ISHELL+1;
      }
    
      //  ****  Projected momentum of the target electron.
    
      A = RAND(6.0)*RN[ISHELL-1];
      if(A < 0.5)
      {
        PZOMC = (D1-sqrt(D12-log(A+A)))/(D2*FJ0[M-1][ISHELL-1]);
      }
      else
      {
        PZOMC = (sqrt(D12-log(2.0-A-A))-D1)/(D2*FJ0[M-1][ISHELL-1]);
      }
      if(PZOMC < -1.0){ brkIt = false; continue;}
    
      //  ****  F(EP) rejection.
    
      XQC = 1.0+TAU*(TAU-2.0*CDT);
      AF = sqrt(XQC)*(1.0+TAU*(TAU-CDT)/XQC);
      if(AF > 0.0)
      {
        FPZMAX = 1.0+AF*0.2;
      }
      else
      {
        FPZMAX = 1.0-AF*0.2;
      }
    
      if(PZOMC > 0.2){ FPZ = 0.2;}
      else{ FPZ = PZOMC;}
    
      if(FPZ < -0.2){ FPZ = -0.2;}
    
      FPZ = 1.0+AF*FPZ;
      if(RAND(7.0)*FPZMAX > FPZ){ brkIt = false; continue;}
    }
    //  ****  Energy of the scattered photon.
      
    T = PZOMC*PZOMC;
    B1 = 1.0-T*TAU*TAU;
    B2 = 1.0-T*TAU*CDT;
    if(PZOMC > 0.0)
    {
      EP = E_*(TAU/B1)*(B2+sqrt(fabs(B2*B2-B1*(1.0-T))));
    }
    else
    {
      EP = E_*(TAU/B1)*(B2-sqrt(fabs(B2*B2-B1*(1.0-T))));
    }
  }

  if(Egt5MeV)
  {
    //  ****  No Doppler broadening for E_ greater than 5 MeV.
    bool brkIt = false;
    while(!brkIt)
    {
      brkIt = true;
      if(RAND(8.0)*A2 < A1)
      {
        TAU = pow(TAUMIN,RAND(9.0));
      }
      else
      {
        TAU = sqrt(1.0+RAND(10.0)*(TAUM2-1.0));
      }
      //  ****  Rejection function.
      TST = (1.0+TAU*(EK1+TAU*(EK2+TAU*EKS)))/(EKS*TAU*(1.0+TAU*TAU));
      if(RAND(11.0) > TST){ brkIt = false; continue;}
      EP = TAU*E_;
      CDT = 1.0-(1.0-TAU)/(EK*TAU);

      //  ****  Target electron shell.

      TST = RAND(12.0);
      //  ****  Binary search.
      if(TST < PTRSH[M-1][0])
      {
        ISHELL = 1;
      }
      else
      {
        ISHELL = 1;
        int JO = NOSCCO[M-1]+1;
        bool brkIt2 = false;
        while(!brkIt2)
        {
          brkIt2 = true;
          int I = (ISHELL+JO)/2;
          if(TST > PTRSH[M-1][I-1])
          {
            ISHELL = I;
          }
          else
          {
            JO = I;
          }
          if(JO-ISHELL > 1){ brkIt2 = false; continue;}
        }
        ISHELL = ISHELL+1;
      }
      if(EP > E_-UICO[M-1][ISHELL-1]){ brkIt = false; continue;}
    }
  }
  DE = E_-EP;
  if(KSCO[M-1][ISHELL-1] < 17)
  {
    if(UICO[M-1][ISHELL-1] > ECUTR[M-1])
    {
      ES = DE-UICO[M-1][ISHELL-1];
    }
    else
    {
      ES = DE;
    }
  }
  else
  {
    ES = DE;
  }
  
  double Q2 = E_*E_+EP*(EP-2.0*E_*CDT);
  if(Q2 > 1.0E-12)
  {
    CDTS = (E_-EP*CDT)/sqrt(Q2);
  }
  else
  {
    CDTS = 1.0;
  }

  IZZ = KZCO[M-1][ISHELL-1];
  ISH = KSCO[M-1][ISHELL-1];
}

