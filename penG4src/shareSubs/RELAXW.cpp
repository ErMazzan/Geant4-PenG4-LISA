
//  *********************************************************************
//                       SUBROUTINE RELAXW
//  *********************************************************************
void PenInterface::RELAXW(int &IZ, FILE* IWR)
{
  //  This subroutine produces a table of atomic relaxation data for the
  //  element IZ, and prints it on unit IWR. The output table is part of
  //  PENELOPE's material definition file.
  
  //  Data are read from file 'pdrelax.p11', which contains data pertaining
  //  to singly ionized atoms with the initial vacancy in one of the K, L,
  //  M and N shells. This file was prepared from the Livermore Evaluated
  //  Atomic Data Library (EADL). The energies of x-ray lines were replaced
  //  by more accurate experimental and theoretical values given by
  //  Deslattes et al. (2003) -K and L shells- and by Burr (1967) -M
  //  shells.
  
  //  NOTE: The transition probabilities and emission energies can be
  //  modified by editing the material data file. For each initial vacancy,
  //  the sum of transition probabilities _must_ be equal to unity.

  using namespace PENERROR_mod;
  using namespace CADATA;
  
  char CSH5[30][6] = {"1s1/2","2s1/2","2p1/2","2p3/2","3s1/2","3p1/2","3p3/2","3d3/2","3d5/2","4s1/2","4p1/2","4p3/2","4d3/2","4d5/2","4f5/2","4f7/2","5s1/2","5p1/2","5p3/2","5d3/2","5d5/2","5f5/2","5f7/2","6s1/2","6p1/2","6p3/2","6d3/2","6d5/2","7s1/2"," free"};
  
  const int NM = 2500;
  int IS0[NM], IS1[NM], IS2[NM];
  double P[NM], EI[NM], EE[99];
      
  int NT = 0;
  double ET;

  if(NSHT[IZ-1] <= 0){ ErrorFunction(1344); return;}
  FILE* pdrelax = fopen("./pdfiles/pdrelax.p11", "r");
  int IZR, IS0R;
  int Nombre_elements_llegits = fscanf(pdrelax, "%d %d%*[^\n]", &IZR, &IS0R); // Ignores the data.
  getc(pdrelax);
  if( Nombre_elements_llegits == 2)
  {
    NT = 0;
    int IS1R, IS2R;
    double PR, EIN;
    for(int I = 0; I < 150000; I++)
    {
      Nombre_elements_llegits = fscanf(pdrelax, "%d %d %d %d %lf %lf%*[^\n]", &IZR, &IS0R, &IS1R, &IS2R, &PR, &EIN);
      getc(pdrelax);
      if(Nombre_elements_llegits != 6){ break;}

      if(IZR == IZ)
      {
        NT = NT+1;
        if(NT > NM){ ErrorFunction(1345); return;}
        IS0[NT-1] = IS0R;
        IS1[NT-1] = IS1R;
        IS2[NT-1] = IS2R;
        P[NT-1] = PR;
        EI[NT-1] = EIN;
      }
    }
  }
  fclose(pdrelax);

  fprintf(IWR, " *** RELAX:  Z =%3d,  no. of shells =%3d,  no. of transitions =%5d\n", IZ, NSHT[IZ-1], NT);

  for(int I = 0; I < 99; I++)
  {
    EE[I] = 0.0;
  }
  int KS;
  for(int J = 0; J < 30; J++)
  {
    KS = IKS[IZ-1][J];
    if(KS > 0)
    {
      if(IFI[IZ-1][KS-1] != 0)
      {
        fprintf(IWR, " %3d %5s %1d %12.5E %12.5E %12.5E\n", KS, CSH5[KS-1], IFI[IZ-1][KS-1], EB[IZ-1][KS-1], ALW[IZ-1][KS-1], CP0[IZ-1][KS-1]);
        EE[KS-1] = EB[IZ-1][KS-1];
      }
    }
  }
  
  if(NT > 0)
  {
    for(int I = 0; I < NT; I++)
    {
      if(IS2[I] == 0)
      {
        if(EI[I] < 1.0)
        {
          ET = EE[IS0[I]-1]-EE[IS1[I]-1];
        }
        else
        {
          ET = EI[I];
        }
      }
      else
      {
        if(EI[I] < 1.0)
        {
          ET = EE[IS0[I]-1]-EE[IS1[I]-1]-EE[IS2[I]-1];
        }
        else
        {
          ET = EI[I];
        }
      }
      if(ET < 1.0){ ET = 1.0;}
      fprintf(IWR, " %3d%3d%3d %12.5E %12.5E\n", IS0[I], IS1[I], IS2[I], P[I], ET);
    }
  }     
}

