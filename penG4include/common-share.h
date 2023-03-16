 
// En aquest fitxer guardarem les variables dels commons i moduls amb el mateix nom que el namespace

#ifndef COMMON_SHARE_H
#define COMMON_SHARE_H

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "PenelopeDefines.hh"

namespace PENELOPE_mod
{
  //  ****  Simulation parameters (must be defined before calling the
  //        initialisation subroutine PEINIT).

  //  ----  Absorption energies, EABS(KPAR,MAT).
  double EABS[3][MAXMAT];

  //  ----  Electron/positron transport parameters.
  double C1[MAXMAT];
  double C2[MAXMAT];
  double WCC[MAXMAT];
  double WCR[MAXMAT];
  
  //  ****  Global information on the material system (defined by
  //        subroutine PEINIT).
  //  ----  Number of materials present.
  int NMAT;
  //  ----  Material densities and its reciprocals.
  double DEN[MAXMAT];
  double RDEN[MAXMAT];

}

namespace CECUTR
{
  //  ****  Simulation parameters.
  double ECUTR[MAXMAT];
}
namespace CSGAWR
{
  int ISGAW; //Controls warning messages from SUMGA.
}
namespace COMPOS
{
  //  ****  Composition data.
  double STF[MAXMAT][30], ZT[MAXMAT], AT[MAXMAT], RHO[MAXMAT], VMOL[MAXMAT];
  int IZ[MAXMAT][30], NELEM[MAXMAT];
  
}
namespace CADATA
{
  //  ****  Element data.
  double EB[99][30], ALW[99][30], CP0[99][30];
  int IFI[99][30], IKS[99][30], NSHT[99];

//  *********************************************************************
//                       BLOCK DATA PENDAT
//  *********************************************************************
//      BLOCK DATA PENDAT

//  Physical data for the elements Z=1-99.
  
//  ************  Chemical symbols of the elements.
// COMPTE ! LASYMB[2] definida en moltes subrutines està mal. Ve del fortran char*2 LASYMB. ELIMINAR !

  char LASYMB[99][3] = {"H ","He","Li","Be","B ","C ","N ","O ","F ","Ne","Na","Mg","Al","Si","P ","S ","Cl","Ar","K ","Ca","Sc","Ti","V ","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y ","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I ","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W ","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U ","Np","Pu","Am","Cm","Bk","Cf","Es"};

  //  ************  Atomic weights (mean relative atomic masses).

  double ATW[99] = {1.0079, 4.0026, 6.9410, 9.0122, 1.0811E1, 1.2011E1, 1.4007E1, 1.5999E1, 1.8998E1, 2.0179E1, 2.2990E1, 2.4305E1, 2.6982E1, 2.8086E1, 3.0974E1, 3.2066E1, 3.5453E1, 3.9948E1, 3.9098E1, 4.0078E1, 4.4956E1, 4.7880E1, 5.0942E1, 5.1996E1, 5.4938E1, 5.5847E1, 5.8933E1, 5.8690E1, 6.3546E1, 6.5390E1, 6.9723E1, 7.2610E1, 7.4922E1, 7.8960E1, 7.9904E1, 8.3800E1, 8.5468E1, 8.7620E1, 8.8906E1, 9.1224E1, 9.2906E1, 9.5940E1, 9.7907E1, 1.0107E2, 1.0291E2, 1.0642E2, 1.0787E2, 1.1241E2, 1.1482E2, 1.1871E2, 1.2175E2, 1.2760E2, 1.2690E2, 1.3129E2, 1.3291E2, 1.3733E2, 1.3891E2, 1.4012E2, 1.4091E2, 1.4424E2, 1.4491E2, 1.5036E2, 1.5196E2, 1.5725E2, 1.5893E2, 1.6250E2, 1.6493E2, 1.6726E2, 1.6893E2, 1.7304E2, 1.7497E2, 1.7849E2, 1.8095E2, 1.8385E2, 1.8621E2, 1.9020E2, 1.9222E2, 1.9508E2, 1.9697E2, 2.0059E2, 2.0438E2, 2.0720E2, 2.0898E2, 2.0898E2, 2.0999E2, 2.2202E2, 2.2302E2, 2.2603E2, 2.2703E2, 2.3204E2, 2.3104E2, 2.3803E2, 2.3705E2, 2.3905E2, 2.4306E2, 2.4707E2, 2.4707E2, 2.5108E2, 2.5208E2};

  //  ************  Mean excitation energies of the elements (eV).

  double EPX[99] ={ 19.2, 41.8, 40.0, 63.7, 76.0, 81.0, 82.0, 95.0, 115.0, 137.0, 149.0, 156.0, 166.0, 173.0,173.0, 180.0, 174.0, 188.0, 190.0, 191.0, 216.0, 233.0, 245.0, 257.0, 272.0, 286.0, 297.0, 311.0, 322.0, 330.0, 334.0, 350.0, 347.0, 348.0, 343.0, 352.0, 363.0, 366.0, 379.0, 393.0, 417.0, 424.0, 428.0, 441.0, 449.0, 470.0, 470.0, 469.0, 488.0, 488.0, 487.0, 485.0, 491.0, 482.0, 488.0, 491.0, 501.0, 523.0, 535.0, 546.0, 560.0, 574.0, 580.0, 591.0, 614.0, 628.0, 650.0, 658.0, 674.0, 684.0, 694.0, 705.0, 718.0, 727.0, 736.0, 746.0, 757.0, 790.0, 790.0, 800.0, 810.0, 823.0, 823.0, 830.0, 825.0, 794.0, 827.0, 826.0, 841.0, 847.0, 878.0, 890.0, 902.0, 921.0, 934.0, 939.0, 952.0, 966.0, 980.0};

  //  ************  Pair-production cross section parameters.

  //  ****  Screening parameter (R mc/hbar).
  double RSCR[99] = {1.2281E2, 7.3167E1, 6.9228E1, 6.7301E1, 6.4696E1, 6.1228E1, 5.7524E1, 5.4033E1, 5.0787E1, 4.7851E1, 4.6373E1, 4.5401E1, 4.4503E1, 4.3815E1, 4.3074E1, 4.2321E1, 4.1586E1, 4.0953E1, 4.0524E1, 4.0256E1, 3.9756E1, 3.9144E1, 3.8462E1, 3.7778E1, 3.7174E1, 3.6663E1, 3.5986E1, 3.5317E1, 3.4688E1, 3.4197E1, 3.3786E1, 3.3422E1, 3.3068E1, 3.2740E1, 3.2438E1, 3.2143E1, 3.1884E1, 3.1622E1, 3.1438E1, 3.1142E1, 3.0950E1, 3.0758E1, 3.0561E1, 3.0285E1, 3.0097E1, 2.9832E1, 2.9581E1, 2.9411E1, 2.9247E1, 2.9085E1, 2.8930E1, 2.8721E1, 2.8580E1, 2.8442E1, 2.8312E1, 2.8139E1, 2.7973E1, 2.7819E1, 2.7675E1, 2.7496E1, 2.7285E1, 2.7093E1, 2.6911E1, 2.6705E1, 2.6516E1, 2.6304E1, 2.6108E1, 2.5929E1, 2.5730E1, 2.5577E1, 2.5403E1, 2.5245E1, 2.5100E1, 2.4941E1, 2.4790E1, 2.4655E1, 2.4506E1, 2.4391E1, 2.4262E1, 2.4145E1, 2.4039E1, 2.3922E1, 2.3813E1, 2.3712E1, 2.3621E1, 2.3523E1, 2.3430E1, 2.3331E1, 2.3238E1, 2.3139E1, 2.3048E1, 2.2967E1, 2.2833E1, 2.2694E1, 2.2624E1, 2.2545E1, 2.2446E1, 2.2358E1, 2.2264E1};
  //  ****  Asymptotic triplet contribution (eta).
  double ETA[99] = {1.1570, 1.1690, 1.2190, 1.2010, 1.1890, 1.1740, 1.1760, 1.1690, 1.1630, 1.1570, 1.1740, 1.1830, 1.1860, 1.1840, 1.1800, 1.1780, 1.1750, 1.1700, 1.1800, 1.1870, 1.1840, 1.1800, 1.1770, 1.1660, 1.1690, 1.1660, 1.1640, 1.1620, 1.1540, 1.1560, 1.1570, 1.1580, 1.1570, 1.1580, 1.1580, 1.1580, 1.1660, 1.1730, 1.1740, 1.1750, 1.1700, 1.1690, 1.1720, 1.1690, 1.1680, 1.1640, 1.1670, 1.1700, 1.1720, 1.1740, 1.1750, 1.1780, 1.1790, 1.1800, 1.1870, 1.1940, 1.1970, 1.1960, 1.1940, 1.1940, 1.1940, 1.1940, 1.1940, 1.1960, 1.1970, 1.1960, 1.1970, 1.1970, 1.1980, 1.1980, 1.2000, 1.2010, 1.2020, 1.2040, 1.2050, 1.2060, 1.2080, 1.2070, 1.2080, 1.2120, 1.2150, 1.2180, 1.2210, 1.2240, 1.2270, 1.2300, 1.2370, 1.2430, 1.2470, 1.2500, 1.2510, 1.2520, 1.2550, 1.2560, 1.2570, 1.2590, 1.2620, 1.2620, 1.2650};

}

namespace CRANGE
{
  //  ****  Energy grid and interpolation constants for the current energy.
  double RANGE[3][MAXMAT][NEGP], RANGEL[3][MAXMAT][NEGP];  
}
namespace CEIN
{
  //  ****  E/P inelastic collisions.
  double EXPOT[MAXMAT], OP2[MAXMAT], F[MAXMAT][NO], UI[MAXMAT][NO], WRI[MAXMAT][NO];
  int KZ[MAXMAT][NO], KS[MAXMAT][NO], NOSC[MAXMAT];  
}
namespace CEINTF
{
  //  ****  E/P inelastic collisions.
  double T1EI[NEGP], T2EI[NEGP], T1PI[NEGP], T2PI[NEGP];
}
namespace CEIN00
{
  //  ****  Partial cross sections of individual shells/oscillators.
  // Associat a les routines PEMATR i EINaT solament.
  double SEH0[NO], SEH1[NO], SEH2[NO], SES0[NO], SES1[NO], SES2[NO], SET0[NO], SET1[NO], SET2[NO];  
}
namespace CPIN00
{
  //  ****  Partial cross sections of individual shells/oscillators.
  // Associat a les routines PEMATR i PINaT solament. COMPTE !! Les variables del COMMON tene noms differents en les dues.
  double SPH0[NO], SPH1[NO], SPH2[NO], SPS0[NO], SPS1[NO], SPS2[NO], SPT0[NO], SPT1[NO], SPT2[NO];  
}
namespace CESI0
{
  //  ****  Inner-shell ionisation by electron and positron impact.
  double XESI[NRP][16];
  int IESIF[99], IESIL[99], NSESI[99], NCURE;
}
namespace CPSI0
{
  //  ****  Inner-shell ionisation by electron and positron impact.
  double XPSI[NRP][16];
  int IPSIF[99], IPSIL[99], NSPSI[99], NCURP;  
}
namespace CEINAC
{
  //  ****  Electron inelastic coll. and inner-shell ionisation tables.
  double EINAC[MAXMAT][NEGP][NO];
  int IEIN[MAXMAT][NO], NEIN[MAXMAT];  
}
namespace CESIAC
{
  //  ****  Electron inelastic coll. and inner-shell ionisation tables.
  double ESIAC[MAXMAT][NEGP][NO];
  int IESI[MAXMAT][NO], NESI[MAXMAT];
}
namespace CESIN
{
  //  ****  Electron inelastic coll. and inner-shell ionisation tables.
  double XSEIN[NEGP][NO], XSESI[NEGP][NO];
  int ISIE[NO];
}
namespace CPINAC
{
  //  ****  Positron inelastic coll. and inner-shell ionisation tables.
  double PINAC[MAXMAT][NEGP][NO];
  int IPIN[MAXMAT][NO], NPIN[MAXMAT];  
}
namespace CPSIAC
{
  //  ****  Positron inelastic coll. and inner-shell ionisation tables.
  double PSIAC[MAXMAT][NEGP][NO];
  int IPSI[MAXMAT][NO], NPSI[MAXMAT];  
}
namespace CPSIN
{
  //  ****  Positron inelastic coll. and inner-shell ionisation tables.
  double XSPIN[NEGP][NO], XSPSI[NEGP][NO];
  int ISIP[NO];
}
namespace CGCO
{
  //  ****  Compton scattering.
  double FCO[MAXMAT][NOCO], UICO[MAXMAT][NOCO], FJ0[MAXMAT][NOCO], PTRSH[MAXMAT][NOCO];
  int KZCO[MAXMAT][NOCO], KSCO[MAXMAT][NOCO], NOSCCO[MAXMAT];  
}
namespace CEIMFP
{
  //  ****  Electron simulation tables.
  double SEHEL[MAXMAT][NEGP], SEHIN[MAXMAT][NEGP], SEISI[MAXMAT][NEGP], SEHBR[MAXMAT][NEGP], SEAUX[MAXMAT][NEGP], SETOT[MAXMAT][NEGP], CSTPE[MAXMAT][NEGP], RSTPE[MAXMAT][NEGP], DEL[MAXMAT][NEGP], W1E[MAXMAT][NEGP], W2E[MAXMAT][NEGP], DW1EL[MAXMAT][NEGP], DW2EL[MAXMAT][NEGP], RNDCE[MAXMAT][NEGP], AE[MAXMAT][NEGP], BE[MAXMAT][NEGP], T1E[MAXMAT][NEGP], T2E[MAXMAT][NEGP];  
}
namespace CLAS1E
{
  //  ****  Electron simulation tables.
  double TSTPE[MAXMAT][NEGP], TSTRE[MAXMAT][NEGP], TRL1E[MAXMAT][NEGP], TRL2E[MAXMAT][NEGP];
}
namespace CPIMFP
{
  //  ****  Positron simulation tables.
  double SPHEL[MAXMAT][NEGP], SPHIN[MAXMAT][NEGP], SPISI[MAXMAT][NEGP], SPHBR[MAXMAT][NEGP], SPAN[MAXMAT][NEGP], SPAUX[MAXMAT][NEGP], SPTOT[MAXMAT][NEGP], CSTPP[MAXMAT][NEGP], RSTPP[MAXMAT][NEGP], W1P[MAXMAT][NEGP], W2P[MAXMAT][NEGP], DW1PL[MAXMAT][NEGP], DW2PL[MAXMAT][NEGP], RNDCP[MAXMAT][NEGP], AP[MAXMAT][NEGP], BP[MAXMAT][NEGP], T1P[MAXMAT][NEGP], T2P[MAXMAT][NEGP];
}
namespace CLAS1P
{
  //  ****  Positron simulation tables.
  double TSTPP[MAXMAT][NEGP], TSTRP[MAXMAT][NEGP], TRL1P[MAXMAT][NEGP], TRL2P[MAXMAT][NEGP];  
}
namespace CEEL00
{
  //  ****  Elastic scattering of electrons and positrons.
  double EJT[NEGP], XE0[NEGP], XE1[NEGP], XE2[NEGP], XP0[NEGP], XP1[NEGP], XP2[NEGP], T1E0[NEGP], T2E0[NEGP], T1P0[NEGP], T2P0[NEGP], EJTL[NEGP], FJL[NEGP], A[NEGP], B[NEGP], C[NEGP], D[NEGP];  
}
namespace CBRYLD
{
  //  ****  Electron and positron radiative yields.
  double EBRY[MAXMAT][NEGP], PBRY[MAXMAT][NEGP];  
}
namespace CGIMFP
{
  //  ****  Photon simulation tables.
  double SGRA[MAXMAT][NEGP], SGCO[MAXMAT][NEGP], SGPH[MAXMAT][NEGP], SGPP[MAXMAT][NEGP], SGAUX[MAXMAT][NEGP];  
}
namespace CGPH01
{
  //  ****  Photon simulation tables.
  double ER[NDIM], XSR[NDIM];
  int NPHD;  
}
namespace CGPP01
{
  //  ****  Photon simulation tables.
  double TRIP[MAXMAT][NEGP];
}
namespace CEBR
{
  //  ****  Bremsstrahlung emission.
  double WB[NBW], PBCUT[MAXMAT][NEGP], WBCUT[MAXMAT][NEGP], PDFB[MAXMAT][NEGP][NBW], DPDFB[MAXMAT][NEGP][NBW], PACB[MAXMAT][NEGP][NBW], ZBR2[MAXMAT];  
}
namespace CGRA01
{
  //  ****  Rayleigh scattering.
  double FF[MAXMAT][NQ], ERA[NEX], XSRA[MAXMAT][NEX];
  int IED[NEGP], IEU[NEGP], NE;
}
namespace CGPH00
{
  //  ****  Photoelectric cross sections.
  double EPH[NTP], XPH[NTP][17];
  int IPHF[99], IPHL[99], NPHS[99], NCUR;  
}
namespace CEIN01
{
  double EI, EE, CPS, AMOL, MOM;
}
namespace CPIN01
{
  double EI,CPS,BHA1,BHA2,BHA3,BHA4;
  int MOM;
}
namespace CEBR01
{
  double EBT[NBE], XS[NBE][NBW], TXS[NBE], X[NBE], Y[NBE];
}
namespace CEBR02
{
  double P0[MAXMAT][NEGP][NBW];
}
namespace CBRANG
{
  const int NET = 7;
  const int NKT = 41;
  double BET[NET], BK[NKT], BP1[MAXMAT][NET][NKT][4], BP2[MAXMAT][NET][NKT][4], ZBEQ[MAXMAT];
}
namespace CGRA03
{
  double QRA[NP][MAXMAT], PRA[NP][MAXMAT], DPRA[NP][MAXMAT], ARA[NP][MAXMAT], BRA[NP][MAXMAT], PMAX[NEGP][MAXMAT];
  int ITLRA[NP][MAXMAT], ITURA[NP][MAXMAT];  
}
namespace CGRA00
{
  double FACTE, Q2MAX;
  int MM, MOM;
}
namespace CGRA02
{
  double QQ[NQ], AR[MAXMAT][NQ], BR[MAXMAT][NQ], CR[MAXMAT][NQ], DR[MAXMAT][NQ], FF0[MAXMAT], QQM;  
}
// en penelope.f els noms dels arrays són diferents
//namespace CRITA
//{
//  const unsigned int NM=512;
//  double QTI[NM], PACI[NM], DPACI[NM], AI[NM], BI[NM];
//  int ITLI[NM], ITUI[NM], NPM1I;  
//}
// en rita.f els noms dels arrays són diferents
//Podriem resoldre aquest problema emprant punters
namespace CRITA
{
  const int NM=512;
  double XT[NM], PAC[NM], DPAC[NM], A[NM], B[NM];
  int IL[NM], IU[NM], NPM1;  
  double* XTI=XT;
  double* QTI=XT;
  double* PACI=PAC;
  double* DPACI=DPAC;
  double* AI=A;
  double* BI=B;
  int* ITLI=IL;
  int* ITUI=IU;
  int* NPM1I=&NPM1;
}
namespace CGCO00
{
  // GCOaT, GCOaD (canvi de E --> EE, M --> MM, IO --> IOSC)
  double EE;
  int MM, IOSC;
}
namespace CGPP00
{
  double ZEQPP[MAXMAT], F0[MAXMAT][2], BCB[MAXMAT];
}
namespace CRELAX
{
  double P[NRX], ET[NRX], F[NRX];
  int IS0[NRX], IS1[NRX], IS2[NRX], IFIRST[99][16], ILAST[99][16], NCUR, KS, MODER;
  int* IAL=IS0;
}
namespace CDCSEP
{
  #define NE 96
  double ETS[NE], ETL[NE], TH[NA], THR[NA], XMU[NA], XMUL[NA], ECS[NE], ETCS1[NE], ETCS2[NE], EDCS[NE][NA], PCS[NE], PTCS1[NE], PTCS2[NE], PDCS[NE][NA], DCSI[NA], DCSIL[NA], CSI, TCS1I, TCS2I;  
  #undef NE
}
namespace CEEL01
{
  double EJT[NEGP], XE0[NEGP], XE1[NEGP], XE2[NEGP], XP0[NEGP], XP1[NEGP], XP2[NEGP], T1E0[NEGP], T2E0[NEGP], T1P0[NEGP], T2P0[NEGP], EJTL[NEGP], FL[NEGP], A[NEGP], B[NEGP], C[NEGP], D[NEGP];  
}
namespace CEELDB
{
  #define NP_P 128
  double XSE[NP_P][NEGP][MAXMAT], PSE[NP_P][NEGP][MAXMAT], ASE[NP_P][NEGP][MAXMAT], BSE[NP_P][NEGP][MAXMAT];
  int ITLE[NP_P][NEGP][MAXMAT], ITUE[NP_P][NEGP][MAXMAT];  
  #undef NP_P
}
namespace CPELDB
{
  #define NP_P 128
  double XSP[NP_P][NEGP][MAXMAT], PSP[NP_P][NEGP][MAXMAT], ASP[NP_P][NEGP][MAXMAT], BSP[NP_P][NEGP][MAXMAT];
  int ITLP[NP_P][NEGP][MAXMAT], ITUP[NP_P][NEGP][MAXMAT];
  #undef NP_P 
}
namespace CELSEP
{
  double EELMAX[MAXMAT], PELMAX[MAXMAT], RNDCEd[MAXMAT][NEGP], RNDCPd[MAXMAT][NEGP];  
}
namespace CSUMGA
{
  int IERGA, NCALL;  // Error code, no. of function calls.
}
namespace CRNDG3
{
  double X_[NR], A[NR], B[NR], F[NR];
  int KA[NR], NPM1;  
}
namespace CRITAA
{
  const unsigned int NM=512;
  double XA[NM], AA[NM], BA[NM], FA[NM];
  int IA[NM], NPM1A; 
}

namespace CRITAN
{
  double CNORM; 
}

//************************** PENGEON *************************

namespace PENGEOM_mod
{

  //  ****  Geometry definition parameters and I/O quantities.

  //  ----  Geometry array sizes.
  //        Maximum numbers of surfaces, bodies, and limiting elements.
  const int NS = 10000;
  const int NB = 5000;
  const int NXG = 250;
  //        Number of bodies in the material system (given by PENGEOM).
  int NBODY;

  inline int initialize_chararray(char DUMMY[NB][5])
  {
    for(unsigned int i = 0; i < NB; i++)
    {
      strcpy(DUMMY[i], "    ");
    }
    return 0;
  }
  //  ----  Body aliases (user labels).
  char BALIAS[NB][5];
  int dummyint=initialize_chararray(BALIAS);
  //char** BALIAS = new char*[NB];
  //for(unsigned int i = 0; i < NB; i++)
  //{
  //  BALIAS[i] = new char[5] = {' ', ' ', ' ', ' ', '\0'};
  //}

  //  ----  Body materials. MATER(KB) is the material in body KB.
  inline int initialize_intarray(int DUMMY[NB])
  {
    for(unsigned int i = 0; i < NB; i++)
    {
      DUMMY[i]=0;
    }
    return 0;
  }
  int MATER[NB];
  int dummyint2=initialize_intarray(MATER);
  
  //int* MATER = new int[NB];
  //for(unsigned int i = 0; i < NB; i++)
  //{
  //  MATER[i] = 0;
  //}

  //  ----  Detector definition.
  //        KDET(KB)=ID if body KB is part of detector ID.
  int KDET[NB];
  int dummyint3=initialize_intarray(KDET);

  //int* KDET = new int[NB];
  //for(unsigned int i = 0; i < NB; i++)
  //{
  //  KDET[i] = 0;
  //}

  //  ****  Warning messages for accidental undershots or round-off errors
  //        are issued when LVERB=.TRUE.

  bool LVERB = false;

  //  ****  Last step features (output from subroutine STEP).
  //
  //  ----  Travelled path length, including segments in void volumes.
  double DSTOT;
  //  ----  Label of the last interface crossed by the particle before
  //        entering a material body, or when leaving a material body to
  //        enter a void volume (defined only when NCROSS /= 0).
  int KSLAST;
  
}

namespace QSURF
{
  double AXX[PENGEOM_mod::NS];
  double AXY[PENGEOM_mod::NS];
  double AXZ[PENGEOM_mod::NS];
  double AYY[PENGEOM_mod::NS];
  double AYZ[PENGEOM_mod::NS];
  double AZZ[PENGEOM_mod::NS];
  double AX[PENGEOM_mod::NS];
  double AY[PENGEOM_mod::NS];
  double AZ[PENGEOM_mod::NS];
  double A0[PENGEOM_mod::NS];
  int NSURF;
  int KPLANE[PENGEOM_mod::NS];  
}
namespace QBODY
{
  int KBODY[PENGEOM_mod::NB][PENGEOM_mod::NXG];
  int KBOMO[PENGEOM_mod::NB];
}
namespace QTREE
{
  int NBODYS;
  int KMOTH[PENGEOM_mod::NB];
  int KDGHT[PENGEOM_mod::NB][PENGEOM_mod::NXG];
  int KSURF[PENGEOM_mod::NB][PENGEOM_mod::NXG];
  int KFLAG[PENGEOM_mod::NB][PENGEOM_mod::NXG];

  int KSP[PENGEOM_mod::NS];
  int NWARN;  
}
//////////////////
//PENEASY/////////
//////////////////
namespace BIGSmod
{
//*******************************************************************
//*    Vars for the BIGS source.                                    *
//*******************************************************************
  bool active,isPosGauss;
  int parsrc,matsrc,nspc,xipol;

  const unsigned int nemax=32000;

  double espc[nemax],pspc[nemax],despc[nemax],esigma;
  double rotv[3][3],rotbox[3][3],cossrc0,dcossrc,phi0,dphi;
  double boxside[3],ctrans[3],xsigma,ysigma,xsp1,xsp2,xsp3;
  const double kfact = 6.0;   // Maximum num. of sigmas in Gaussian dispersion
}
namespace dataTypesMod
{
//*******************************************************************
//*    Defines the size in bytes of the data types employed by      *
//*    penEasy. (Superseeded in Fortran2008 by means of the         *
//*    ISO_FORTRAN_ENV module.)                                     *
//*******************************************************************

  int sizeOfChar  = (int)sizeof(char);
  int sizeOfInt   = (int)sizeof(int);  // CAUTION: may not apply to some (older) systems
  int sizeOfInt4  = (int)sizeof(int);
  int sizeOfInt8  = (int)sizeof(long int);
  int sizeOfReal  = (int)sizeof(float);
  int sizeOfReal8 = (int)sizeof(double);
}
namespace ctrsimMod
{

//*******************************************************************
//*    Vars for simulation control.                                 *
//*******************************************************************
  double nhist, nhistmax, time0, atime;
  double lastComm, lastFresh, refresh, lastDump, refreshDump; 
}
namespace dsmaxMod
{

//*******************************************************************
//*    Vars for step length control.                                *
//*******************************************************************
      double DSMAX[MAXMAT];      
}
namespace dumpMod
{
//*******************************************************************
//*    Vars for dump file.                                          *
//*******************************************************************

  char dumpfilen[80];
  FILE* restartfile;
  FILE* dumpfile;
  int DumpFile;
  int RestartFile;
}
namespace EDPmod
{

  bool active;
  int matdet;
  double edptmp[MAXMAT],edep[MAXMAT],edep2[MAXMAT];
  double unclimit;
  
}
namespace SDDmod
{

  bool active;
  int  prtxyz;
  int  nx,ny,nz,nxny,nbin;
  double xmin,ymin,zmin,dx,dy,dz,idx,idy,idz,unclimit;
  const double nbinmax=1.024E9; // Max num bins (do not exceed int*4)
  double memSDD=0.0; //Memory used by arrays (bytes)
  //ALLOCATABLE
  double* edptmp;
  double* edep;
  double* edep2;
  double* idens;
  double* nlast;
  
}
namespace CDDmod
{

  bool active;
  int  prtxyz;
  int  nr,nz;
  const int nbinmax=100000;
  double edptmp[nbinmax];
  double edep[nbinmax];
  double edep2[nbinmax];
  double idens[nbinmax];
  double nlast[nbinmax];
  double dr,dz,idr,idz,unclimit,rmin,zmin,r2min,r2max; 
}
//Moduls del PENEASY per a accedir als commons del PENELOPE. Farem que aquests moduls canvien el seu nom i agafen el del common del PENELOPE corresponent
#define RSEEDcommonMod RSEED
#define SECSTcommonMod SECST
#define COMPOScommonMod COMPOS
#define PENGEOMcommonsMod PENGEOM_mod
//#define CSIMPHcommonMod  //Per a protons. No he trobat aquest common

namespace tim001
{
  double lasthour,tshift,rtime0;
  float  utime0;
}
namespace timrestart
{
  double drtime,dutime;
}

namespace PENVARED_mod
{
  //  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  //  ****  Variance-reduction parameters.

  //  ----  Parameter values defined for each body. NBV is the maximum
  //        number of bodies. When using PENGEOM, NBV must have the same
  //        value as the parameter NB in module PENGEOM_mod.
  const int NBV = 5000;
  //  ----  Forcing factors, FORCE(IBODY,KPAR,ICOL).
  inline int initialize_doublearray(double DUMMY[NBV][3][8])
  {
    for(unsigned int i = 0; i < NBV; i++)
    {
      for(unsigned int j = 0; j < 3; j++)
      {
        for(unsigned int k = 0; k < 8; k++)
        {
          DUMMY[i][j][k]=1.0;
        }
      }
    }
    return 0;
  }
  double FORCE[NBV][3][8];
  int dummydouble=initialize_doublearray(FORCE);
  double WFORCE;
        //  ----  Bremsstrahlung splitting numbers, IBRSPL(IBODY).
  inline int initialize_intarray(int DUMMY[NBV])
  {
    for(unsigned int i = 0; i < NBV; i++)
    {
      DUMMY[i]=1;
    }
    return 0;
  }
  int IBRSPL[NBV];
  int dummyint=initialize_intarray(IBRSPL);
        //  ----  Energy deposited in the last event (unweighted).
  double DEA;

  //  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
}

namespace forcingMod
{
  const int maxmat=MAXMAT;
  bool   analog[maxmat][4], isforcing, active;
  double minwght;
}
namespace splitMod
{
  bool   active;
  int    signsplit,mode,nsplit;
  int    matsplit;
  double rotS[3][3],rot[3][3],rotz[3],rorx[3][3],rory[3][3],rorxy[3][3];
  double dr[3],phi0,dphi,idphi,minwght;
}
namespace russiaMod
{
  bool   active;
  int    matrussia;
  double maxwght,psurv;
}

namespace CFORCF
{
  double POR[8];
  double P0[8];
  int    IBR;
  bool   LFORC[8];
}

namespace CWOODC
{
  double STMV[NEGP];
  double STMC[NEGP];
  double STMAX;
}

//PENVOX

namespace geovoxMod
{
//*******************************************************************
//*    PENVOX vars                                                  *
//*******************************************************************
  char voxfilen[81];
  bool isQuad,isOverlay;  // Indicate: .true. if quadrics defined; .true. if quadrics overlay voxels
  int granul;
  const int maxGranul=100; // Max num. points in voxel side to compute mass
  int nx,ny,nz,nxy,nvox,bodymask,matmask;
  double dx,dy,dz,volvox,idx,idy,idz,vbb[3];
  double memvox=0.0;
  //ALLOCATABLE
  int* matvox;
  float* densvox;
  float* massvox;
  double* massvoxd;  // Temporary array to perform integration in DPrec
}

namespace partvoxMod
{
//*******************************************************************
//*    Particle state vars in the voxelized geo.                    *
//*******************************************************************
  int xvox,yvox,zvox,absvox;
  double uold,vold,ivx,ivy,ivz,sdx,sdy,sdz;
}

namespace srcpsf
{
  double es,xs,ys,zs,us,vs,ws,wghts;
  int kpars,ilbs[5],dns;
//  long int kpars,ilbs[5],dns;
}

namespace srcps1
{
  double rot[3][3],xshift,yshift,zshift;
  int split;
  FILE* ufile;
  bool formatiaea,active;
}

namespace srcps2
{
//  long long int npart,restartnpart,usednpart;
//  long long int nline,restartnline,usednline;
//  long int npart,restartnpart,usednpart;
//  long int nline,restartnline,usednline;
  int npart,restartnpart,usednpart;
  int nline,restartnline,usednline;
}

namespace formatVer
{
  bool format2008;
}

namespace srciaeaps
{
//  long long int npiaeaMax,nTopiaea;
//  long int npiaeaMax,nTopiaea;
//  long int sourceRead;
  int npiaeaMax,nTopiaea;
  int sourceRead;
  int xtraPos[7];
  bool isXtraStored[7];
}

//*******************************************************************
//*  Usernames fo the PSF tally.                                      *
//*******************************************************************

namespace PSFmod
{
  bool active;
  int  detmat,formatiaea;
  FILE* psfunit;
//  long long int nhisti,nhlast,npar[4];
  long int nhisti,nhlast,npar[4];
}

namespace PSFendinoMod
{
  int gkpar,gilb[5];
  double   gen,gx,gy,gz,gu,gv,gw;
}

namespace PSFiaeaMod
{
  int sourceWrite;
}

//*******************************************************************
//*  Usernames fo the PSF tally.                                    *
//*******************************************************************

namespace VDDmod
{
  bool active;
  int prtxyz,prtdens;
  int xvoxmin,xvoxmax,yvoxmin,yvoxmax,zvoxmin,zvoxmax;
  double unclimit,memVDD;//=0.0; // Memory used by arrays (bytes)
  //ALLOCATABLE
  double* nlast;
  double* edptmp;
  double* edep;
  double* edep2;
}

//*******************************************************************
//*  Usernames fo the FTL tally.                                    *
//*******************************************************************

namespace FTLmod
{
  bool active,isLinScale;
  int  nbin,matdet;
  const int nbinmax=32000; 
  double flutmp[4][nbinmax],flu[4][nbinmax],flu2[4][nbinmax];
  double egrid[nbinmax],ebingrd[nbinmax];
  double emin,ebin,iebin,eratio,unclimit;
}

//*******************************************************************
//*  Usernames fo the PCS tally.                                    *
//*******************************************************************

namespace PCSmod
{
  bool active;
  int  matdet,nbin;
  const int nbinmax=32000; 
  double countmp[nbinmax][4],counts[nbinmax][4],count2[nbinmax][4];
  double tcount2[4],edeptmp[4],edep[4],edep2[4];
  double tedep2,emin,ebin,iebin,unclimit;
}

//*******************************************************************
//*  Usernames fo the PTS tally.                                    *
//*******************************************************************

namespace PTSmod
{
  bool active;
  FILE* ptsunit;
  long int ntrack,ntrackmax;
  double xlast,ylast,zlast;
}

//*******************************************************************
//*  Usernames fo the PTS tally.                                    *
//*******************************************************************

namespace SPDmod
{

  bool active;
  int  prtxyz,nbin;
  const int nbinmax=32000;
  double edptmp[nbinmax],edep[nbinmax],edep2[nbinmax];
  double idens[nbinmax],nlast[nbinmax];
  double rmin,r2min,r2max,dr,idr,unclimit;
}

//*******************************************************************
//*  Usernames fo the PHS tally.                                    *
//*******************************************************************

namespace PHSmod
{

  bool active,isGaussConvo;
  int  matdet,nbin;
  const int nbinmax=32000;
  double counts[nbinmax],count2[nbinmax];
  double edptmp,emin,ebin,iebin,unclimit,a0,a1;
}

//*******************************************************************
//*  Usernames fo the PHS tally.                                    *
//*******************************************************************

namespace PIDmod
{

  bool active,isLab,isEspread,isEinteg,isPCM;
  int detmode,repformat,collided,selectCollision;
  long int bodydet,matdet,nx,ny,nebin,ijmax,ijkmax;
  double lab2det[2][3],r0[3],dx,dy,idx,idy,c0,c1;
  double unclimit,ethr,emin,emax,ebin,iebin;
  const double fwhm2sig2=1.0/(8.0*log(2.0E0));
  const double npixmax=1.024E9;    // Max num pixels=32k squared (do not exceed int*4)
  double memPID=0.0;      // Memory used by arrays (bytes)
  double* image;
  double* image2;
  double* edeptmp;
  double* nlast;
  double* wghtlast;
}

//*******************************************************************
//*  Username for the LOCATE and STEPSI of PENGEOM                  *
//*******************************************************************

namespace FUZZYGEO
{
  const double FUZZL = 1.0E-12;  
}

//*******************************************************************
//*  Commons for psampler                                           *
//*******************************************************************

namespace CPARIN
{
  double ZTOTM,EXPOTM;
  double FM[NO],UIM[NO],WRIM[NO],QRIM[NO],WTH[NO],FCIN[NO];
  int NOSCM;

/*  double& ZTOT=ZTOTM;
  double& EXPOT=EXPOTM;
  double  (&F)[NO]=FM;
  double  (&UI)[NO]=UIM;
  double  (&WRI)[NO]=WRIM;
  double  (&QRI)[NO]=QRIM;
  int&    NOSC=NOSCM; */
}

namespace CDCSIN
{
  double EM,DELTA,WCCM;
  int IELEC,MOM;
}

namespace DOICSC
{
  double EI,EE,BETA2,RMU1,RMU2,RMU3,AMOL,BHA1,BHA2,BHA3,BHA4,CDIST,C1,C2,CCLO;
  int    JELEC;
}

namespace CADATA_PSAMPLER
{
  double ATW[99],EPX[99],RSCR[99],ETA[99],EB[99][30],ALW[99][30],CP0[99][30];
  int IFI[99][30],IKS[99][30],NSHT[99],LASYMB[99];
}

namespace CEBRA0
{
  double EE, DEE;
  int MM,MOM;
}

namespace CPAN01
{
  double GAM;
  int MOM;
}

namespace CVERIF
{
  const int NVF=12000;
  double XC1[NVF],XC2[NVF],XC3[NVF],XC4[NVF];
  int KC1[NVF],KC2[NVF],KC3[NVF],KC4[NVF],NX1,NX2,NX3,NX4,NCHAN,NK2,NK3,NK4;
}

namespace CGCO01
{
  double EE,EPP;
  int MM;
}

namespace CGPHVS
{
  double GAM,BETA;
  int MOM;
}

namespace CGPP02
{
  double GAM,FACT,G0,BCBT;
  int MOM;
}

namespace PENERROR_mod    //! Initialization errors.
{
//       SAVE   ! Saves all items in the module.
  char REASON[128];   // Warning/error message.
  int IRETRN = 0;         // Return code.
  int IERSEC = 0;     // > 0 => Secondary stack overflow.
}

#endif
