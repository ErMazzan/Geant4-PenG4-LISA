
//  *********************************************************************
//                       SUBROUTINE PEILB4
//  *********************************************************************
void PenInterface::PEILB4(int &ILB4, int &IZZ, int &IS1, int &IS2, int &IS3)
{
  //  This subroutine parses the value of the label ILB(4) and returns the
  //  atomic number (IZZ) and the active electron shells (IS1,IS2,IS3).
  
  IZZ = ILB4/1000000;
  int ITMP = ILB4-IZZ*1000000;
  IS1 = ITMP/10000;
  ITMP = ITMP-IS1*10000;
  IS2 = ITMP/100;
  IS3 = ITMP-IS2*100;
}

