
//  *********************************************************************
//                       SUBROUTINE DIRPOL
//  *********************************************************************
void PenPhys::DIRPOL(double CDT, double &DF, double CONS, double &_SP1, double &_SP2, double &_SP3, double &_U, double &_V, double &_W)
{
  //     This subroutine computes the direction cosines _and_ the Stokes
  //  parameters of a polarised photon after scattering with a given polar
  //  angle.

  //  Input:  _U,_V,_W ... initial direction cosines.
  //          _SP1,_SP2,_SP3 ... initial Stokes parameters.
  //          CDT ..... cosine of the polar scattering angle.
  //          CONS .... constant in the PDF of the azimuthal angle.
  //  Output: _U,_V,_W ... new direction cosines.
  //          _SP1,_SP2,_SP3 ... new Stokes parameters.
  //          DF ...... azimuthal scattering angle.
  //          CDT and CONS remain unchanged.

  const double TWOPI = 2.0*PI;
  
  //  ****  Sampling the azimuthal scattering angle.

  double CDT2 = CDT*CDT;
  double CDT21 = CDT2+1.0;
  double PHA = CDT21+CONS;
  double PHB = 1.0-CDT2;
  double SP0MAX = PHA+PHB*sqrt(_SP1*_SP1+_SP3*_SP3+1.0E-35);
  double SDT, SDF, CDF, S2DF, C2DF, SP3P, SP0P;
  bool brkIt = false;
  while(!brkIt)
  {
    brkIt = true;
    DF = RAND(1.0)*TWOPI;
    SDF = sin(DF);
    CDF = cos(DF);
    S2DF = 2.0*SDF*CDF;
    C2DF = CDF*CDF-SDF*SDF;
    SP3P = S2DF*_SP1+C2DF*_SP3;  // Stokes parameter with new zero-azimuth.
    SP0P = PHA-PHB*SP3P;
    if(RAND(2.0)*SP0MAX > SP0P){ brkIt = false; continue;}
  }

  //  ****  Calculate new Stokes parameters.

  double SP1P = C2DF*_SP1-S2DF*_SP3;  // Stokes parameter with new zero-azimuth.
  double RSP0 = 1.0/SP0P;

  _SP1 = 2.0*CDT*SP1P*RSP0;
  _SP2 = (2.0+CONS)*CDT*_SP2*RSP0;
  _SP3 = (CDT21*SP3P-PHB)*RSP0;

  //  ****  Ensure normalisation.

  double UV = _U*_U+_V*_V;
  double UVW = UV+_W*_W;
  if(fabs(UVW-1.0) > 1.0E-13)
  {
    double FNORM = 1.0/sqrt(UVW);
    _U = FNORM*_U;
    _V = FNORM*_V;
    _W = FNORM*_W;
    UV = _U*_U+_V*_V;
  }

  //  ****  Calculate new direction.

  if(1.0-fabs(CDT) > 1.0E-8)
  {
    SDT = sqrt(PHB);
  }
  else
  {
    SDT = sqrt(2.0*(1.0-fabs(CDT)));
  }

  if(SDT < 1.0E-13)
  {
    if(CDT < 0.0)
    {
      _U = -_U;
      _V = -_V;
      _W = -_W;
    }
  }
  else
  {
    double SDTSDF = SDT*SDF;
    double SDTCDF = SDT*CDF;
    if(UV > 1.0E-26)
    {
      double SUV = sqrt(UV);
      double UN = _U/SUV;
      double VN = _V/SUV;
      _U = _U*CDT+(UN*_W*SDTCDF-VN*SDTSDF);
      _V = _V*CDT+(VN*_W*SDTCDF+UN*SDTSDF);
      _W = _W*CDT-SUV*SDTCDF;
    }
    else
    {
      if(_W > 0.0)
      {
        _U = SDTCDF;
        _V = SDTSDF;
        _W = CDT;
      }
      else
      {
        _U = -SDTCDF;
        _V = -SDTSDF;
        _W = -CDT;
      }
    }
  }
}

