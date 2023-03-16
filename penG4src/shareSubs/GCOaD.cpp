
//  *********************************************************************
//                        FUNCTION GCOaD
//  *********************************************************************
double GCOaD(double CDT)
{
  //  Single differential cross section for photon Compton scattering by
  //  electrons in the IO-th shell, differential in the direction of the
  //  scattered photon only. Evaluated from the incoherent scattering
  //  function.

  //  The energy E of the primary photon is entered through common CGCO00.
  //  The output value GCOaD is the DCS per electron in units of PIELR2.

  using namespace PENELOPE_mod;

  using namespace CGCO;
  using namespace CGCO00;

  const double RREV = 1.0/REV;
  const double D2 = 1.4142135623731;
  const double D1 = 1.0/D2;
  const double D12 = 0.5;

  double GCOaD_RETURN;
  if(EE < UICO[MM][IOSC-1])
  {
    GCOaD_RETURN = 0.0;
    return GCOaD_RETURN;
  }
  //  ****  Energy of the Compton line.
  double CDT1 = 1.0-CDT;
  double EOEC = 1.0+(EE*RREV)*CDT1;
  double ECOE = 1.0/EOEC;
  //  ****  Klein-Nishina X-factor.
  double XKN = EOEC+ECOE-1.0+CDT*CDT;
  //  ****  Incoherent scattering function (analytical profile).
  double AUX = EE*(EE-UICO[MM][IOSC-1])*CDT1;
  double PIMAX = (AUX-REV*UICO[MM][IOSC-1])/(REV*sqrt(AUX+AUX+(UICO[MM][IOSC-1]*UICO[MM][IOSC-1])));
  double SIA;
  if(PIMAX > 0.0)
  {
    SIA = 1.0-0.5*exp(D12-pow(D1+D2*FJ0[MM][IOSC-1]*PIMAX,2));
  }
  else
  {
    SIA = 0.5*exp(D12-pow(D1-D2*FJ0[MM][IOSC-1]*PIMAX,2));
  }
  //  ****  1st order correction, integral of Pz times the Compton profile.
  //        Calculated approximately using a free-electron gas profile.
  double PF = 3.0/(4.0*FJ0[MM][IOSC-1]);
  double QCOE2, P2, DSPZ;
  if(fabs(PIMAX) < PF)
  {
    QCOE2 = 1.0+ECOE*ECOE-2.0*ECOE*CDT;
    P2 = PIMAX*PIMAX;
    DSPZ = sqrt(QCOE2)*(1.0+ECOE*(ECOE-CDT)/QCOE2)*FJ0[MM][IOSC-1]*0.25*(2*P2-(P2*P2)/(PF*PF)-(PF*PF));
    if(DSPZ > -SIA){ SIA = SIA+DSPZ;}
    else{ SIA = 0.0;}      
  }
  //  ****  Differential cross section (per electron, in units of PIELR2).
  GCOaD_RETURN = (ECOE*ECOE)*XKN*SIA;
  return GCOaD_RETURN;
}

