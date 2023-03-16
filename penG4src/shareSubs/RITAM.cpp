
//  *********************************************************************
//                       SUBROUTINE RITAM
//  *********************************************************************
void PenInterface::RITAM(double XD, double XU, double &XM0, double &XM1, double &XM2)
{
  //  Calculation of (restricted) momenta of a pdf, PDF(X), obtained from
  //  its RITA approximation.
  //
  //     XD, XU ... limits of the integration interval.
  //     XM0 ...... 0th order moment.
  //     XM1 ...... 1st order moment.
  //     XM2 ...... 2nd order moment.

  using namespace CRITA;

  const int NIP = 51;
  double XS[NIP], YS[NIP];

  XM0 = 0.0;
  XM1 = 0.0;
  XM2 = 0.0;
  for(int I = 0; I < NPM1; I++)
  {
    if(XT[I+1] >= XD && XT[I] <= XU)
    {
      double X1 = (XT[I] > XD ? XT[I] : XD);
      double X2 = (XT[I+1] < XU ? XT[I+1] : XU);
      double DX = (X2-X1)/double(NIP-1);

      for(int K = 0; K < NIP; K++)
      {
        XS[K] = X1+double(K+1-1)*DX;
        //  ****  Value of the RITA rational pdf at the point XS(K).
        double TAU = (XS[K]-XT[I])/(XT[I+1]-XT[I]);
        double CON1 = 2.0*B[I]*TAU;
        double CON2 = 1.0+B[I]+A[I]*(1.0-TAU);
        double ETA;
        if(fabs(CON1) > 1.0E-10*fabs(CON2))
        {
          ETA = CON2*(1.0-sqrt(1.0-2.0*TAU*CON1/(CON2*CON2)))/CON1;
        }
        else
        {
          ETA = TAU/CON2;
        }
        YS[K] = DPAC[I]*((1.0+(A[I]+B[I]*ETA)*ETA)*(1.0+(A[I]+B[I]*ETA)*ETA))/((1.0-B[I]*ETA*ETA)*(1.0+A[I]+B[I])*(XT[I+1]-XT[I]));
      }

    //  ****  Simpson's integration.

      double CONS = DX*3.3333333333333333E-1;
      double SUM = 0.0;
      for(int LLL = 2; LLL < NIP; LLL += 2)
      {
        SUM = SUM+YS[LLL-2]+4.0*YS[LLL-1]+YS[LLL];
      }
      XM0 = XM0+SUM*CONS;

      for(int K = 0; K < NIP; K++)
      {
        YS[K] = YS[K]*XS[K];
      }
      SUM = 0.0;
      for(int LLL = 2; LLL < NIP; LLL +=2)
      {
        SUM = SUM+YS[LLL-2]+4.0*YS[LLL-1]+YS[LLL];
      }
      XM1 = XM1+SUM*CONS;

      for(int K = 0; K < NIP; K++)
      {
        YS[K] = YS[K]*XS[K];
      }
      SUM = 0.0;
      for(int LLL = 2; LLL < NIP; LLL += 2)
      {
        SUM = SUM+YS[LLL-2]+4.0*YS[LLL-1]+YS[LLL];
      }
      XM2 = XM2+SUM*CONS;
    }
  }     
}
