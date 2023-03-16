
//  *********************************************************************
//                       SUBROUTINE EBRa
//  *********************************************************************
void PenPhys::EBRa( double &E_, double &W_, int &M)
{
  //  Simulation of bremsstrahlung emission by electrons or positrons in
  //  material M.

  using namespace PENELOPE_mod;

///////////  using namespace CEGRID;
  using namespace CEBR;

  if(WCR[M-1] > E_)
  {
    W_ = 0.0;
    return;
  }

  //  ****  Selection of the energy grid point.
  int IE;
  if(RAND(1.0) < CEGRID_.XEK)
  {
    IE = CEGRID_.KE+1;
  }
  else
  {
    IE = CEGRID_.KE;
  }
  //  ****  Pointer.
  double PT, W1, W2, DW, B, A, PMAX;
  int I, J, K;
  bool brkIt = false;
  while(!brkIt)
  {
    brkIt = true;
    PT = PBCUT[M-1][IE-1]+RAND(2.0)*(PACB[M-1][IE-1][NBW-1]-PBCUT[M-1][IE-1]);
      //  ****  Binary search of the W-interval.
    I = 1;
    J = NBW;
    bool brkIt2 = false;
    while(!brkIt2)
    {
      brkIt2 = true;
      K = (I+J)/2;
      if(PT > PACB[M-1][IE-1][K-1])
      {
        I = K;
      }
      else
      {
        J = K;
      }
      if(J-I > 1){ brkIt2 = false; continue;}
    }
      //  ****  Sampling the photon energy (rejection method).
    W1 = WB[I-1];
    W2 = WB[I+1-1];
    DW = W2-W1;
    B = DPDFB[M-1][IE-1][I-1]/DW;
    A = PDFB[M-1][IE-1][I-1]-B*W1;
    if(W1 < WBCUT[M-1][IE-1]){ W1 = WBCUT[M-1][IE-1];}
    if(W2 < W1)
    {
      printf(" **** WARNING: EBR. Conflicting end-point values.\n");
      W_ = W1;
      return;
    }

    PMAX = A+B*W1;
    if(PMAX < A+B*W2){ PMAX = A+B*W2;}

    brkIt2 = false;
    while(!brkIt2)
    {
      brkIt2 = true;
      W_ = W1*pow(W2/W1,RAND(3.0));
      if(RAND(4.0)*PMAX > A+B*W_){ brkIt2 = false; continue;}
    }
    W_ = W_*E_;
    if(W_ < WCR[M-1]){brkIt = false; continue;}
  }
     
}

