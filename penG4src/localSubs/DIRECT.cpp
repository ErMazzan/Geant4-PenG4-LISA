
//  *********************************************************************
//                       SUBROUTINE DIRECT
//  *********************************************************************
void PenPhys::DIRECT(double &CDT, double &DF, double &U_, double &V_, double &W_)
{
  //  This subroutine computes the new direction cosines of the particle
  //  velocity after a collision with given polar and azimuthal scattering
  //  angles.

  //  Input:  U_,V_,W_ ... initial direction cosines.
  //          CDT ..... cosine of the polar scattering angle.
  //          DF ...... azimuthal scattering angle (rad).

  //  Output: U_,V_,W_ ... new direction cosines.
  //          CDT and DF remain unchanged.

  //  ****  Ensure normalisation.
  double FNORM;
  double UV = U_*U_+V_*V_;
  double UVW = UV+W_*W_;
  if(fabs(UVW-1.0) > 1.0E-13)
  {
    FNORM = 1.0/sqrt(UVW);
    U_ = FNORM*U_;
    V_ = FNORM*V_;
    W_ = FNORM*W_;
    UV = U_*U_+V_*V_;
  }

  //  ****  Calculate new direction.

  double SDT;
  if(1.0-fabs(CDT) > 1.0E-8)
  {
    SDT = sqrt(1.0-CDT*CDT);
  }
  else
  {
    SDT = sqrt(2.0*(1.0-fabs(CDT)));
  }

  if(SDT < 1.0E-13)
  {
    if(CDT < 0.0)
    {
      U_ = -U_;
      V_ = -V_;
      W_ = -W_;
    }
  }
  else
  {
    double SDTSDF = SDT*sin(DF);
    double SDTCDF = SDT*cos(DF);
    if(UV > 1.0E-26)
    {
      double SUV = sqrt(UV);
      double UN = U_/SUV;
      double VN = V_/SUV;
      U_ = U_*CDT+(UN*W_*SDTCDF-VN*SDTSDF);
      V_ = V_*CDT+(VN*W_*SDTCDF+UN*SDTSDF);
      W_ = W_*CDT-SUV*SDTCDF;
    }
    else
    {
      if(W_ > 0.0)
      {
        U_ = SDTCDF;
        V_ = SDTSDF;
        W_ = CDT;
      }
      else
      {
        U_ = -SDTCDF;
        V_ = -SDTSDF;
        W_ = -CDT;
      }
    }
  }
}

