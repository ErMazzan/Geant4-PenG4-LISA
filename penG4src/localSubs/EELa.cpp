
//  *********************************************************************
//                       SUBROUTINE EELa
//  *********************************************************************
void PenPhys::EELa(double A, double B, double RNDC, double &RMU)
{
  //  Simulation of hard elastic events. Modified-Wentzel model.
  //
  //  Input arguments:
  //    A, B ... angular distribution parameters.
  //    RNDC ... cutoff probability.
  //  Output values:
  //    RMU .... angular deflection, =(1-CDT)/2.

  double A1 = A+1.0;
  double B1, BB;
  double RND0, RND, RNDMB, RNDRC;
  if(B >= 0.0)
  {
      
    //  ****  Case I.
      
    double RMUAV = A*A1*log(A1/A)-A;
    B1 = 1.0-B;
    RND0 = B1*A1*RMUAV/(A+RMUAV);
    RND = RNDC+RAND(1.0)*(1.0-RNDC);
    if(RND < RND0)
    {
      RMU = RND*A/(B1*A1-RND);
    }
    else if(RND > RND0+B)
    {
      RNDMB = RND-B;
      RMU = RNDMB*A/(B1*A1-RNDMB);
    }
    else
    {
      RMU = RMUAV;
    }
  }
  else
  {

    //  ****  Case II.

    BB = -B;
    B1 = 1.0-BB;
    double RMUC = RNDC*A/(B1*A1-RNDC);
    double PW = B1*A*(1.0-RMUC)/(A+RMUC);
    if(RAND(2.0)*(BB+PW) < BB)
    {
      RMU = 0.5*(1.0+sqrt(RAND(3.0)));
    }
    else
    {
      RNDRC = RAND(3.0)*(1.0-RMUC);
      RMU = (A*RNDRC+A1*RMUC)/(A1-RNDRC);
    }
  }
  
}

