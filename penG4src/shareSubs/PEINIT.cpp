
//  *********************************************************************
//                       SUBROUTINE PEINIT
//  *********************************************************************
void PenInterface::PEINIT(double EMAX, unsigned int NMATER, FILE *IWR, int INFO,
            char PMFILE[MAXMAT][81])
{
  //  Input of material data and initialisation of simulation routines.
  //
  //  Each material is defined through an input file, which is created by
  //  the program MATERIAL using information contained in the database.
  //  This file can be modified by the user if more accurate interaction
  //  data are available.
  //
  //  Input arguments:
  //    EMAX ... maximum particle energy (kinetic energy for electrons and
  //             positrons) used in the simulation. Note: Positrons with
  //             energy E may produce photons with energy E+1.022E6.
  //    NMATER .... number of materials in the geometry.
  //    PMFILE .... array of MAXMAT 20-character strings. The first NMATER
  //             elements are the filenames of the material data files.
  //             The file PMFILE(M) contains radiation interaction data
  //             for material M (that is, the order is important!).
  //    IWR .... output unit.
  //    INFO ... determines the amount of information that is written on
  //             the output file,
  //               INFO=1 (or less), minimal (composition data only).
  //               INFO=2, medium (same information as in the material
  //                 definition data file, useful to check that the struc-
  //                 ture of the latter is correct).
  //               INFO=3 or larger, full information, including tables of
  //                 interaction properties used in the simulation.
  //
  //  For the preliminary computations, PEINIT needs to know the absorption
  //  energies EABS(KPAR,M) and the simulation parameters C1(M), C2(M),
  //  WCC(M) and WCR(M). This information is introduced through the module
  //  PENELOPE_mod, which has to be loaded before calling PEINIT.
  //

  using namespace PENELOPE_mod;
  using namespace PENERROR_mod;
  
  using namespace CECUTR;
  using namespace CSGAWR;
  
  char LIT[3];
  //  ****  Simulation parameters.
  double EABS0[3][MAXMAT];
  
  fprintf(IWR ,"\n **********************************\n **   PENELOPE  (version 2018)   **");
  fprintf(IWR ,"\n **********************************\n");
  ISGAW=0;

  NMAT = NMATER;
  for(int M = 0; M < NMAT; M++)
  {
    if(EABS[0][M] < 49.999)
    {
	    fprintf(IWR, "\nMaterial %3d, EABS(1,%2d) = %11.4E eV\n *** ERROR: electron absorption energy cannot be less than %.5E eV\n", M, M, EABS[0][M], double(MINEABSE));
	    ErrorFunction(1001);return;
	  }
    EABS0[0][M] = EABS[0][M];
    if(EABS[1][M] < 49.999)
	  {
	    fprintf(IWR, "\nMaterial %3d, EABS(2,%2d) = %11.4E eV\n *** ERROR: photon absorption energy cannot be less than %.5E eV\n", M, M, EABS[1][M], double(MINEABSPh));
	    ErrorFunction(1002);return;
	  }
    EABS0[1][M] = EABS[1][M]; 
    if(EABS[2][M] < 49.999)
	  {
	    fprintf(IWR, "\nMaterial %3d, EABS(3,%2d) = %11.4E eV\n *** ERROR: positron absorption energy cannot be less than %.5E eV\n", M, M, EABS[2][M], double(MINEABSPo));
	    ErrorFunction(1003);return;
	  }
    EABS0[2][M] = EABS[2][M];
  }


  //  ****  Lower limit of the energy grid.

  double EMIN = 1.0E35;
  for(int M = 0; M < NMAT; M++)
  {
	  if(EABS[0][M] >= EMAX)
  	{
	    EABS[0][M] = EMAX*(double)0.9999;
	    fprintf(IWR, "\n WARNING: EABS(1,%2d) has been modified because is greater than EMAX\n", M);
	  }
	  if(EABS[1][M] >= EMAX)
	  {
	    EABS[1][M] = EMAX*(double)0.9999;
	    fprintf(IWR, "\n WARNING: EABS(2,%2d) has been modified because is greater than EMAX\n", M);
	  }
	  if(EABS[2][M] >= EMAX)
	  {
      EABS[2][M] = EMAX*(double)0.9999;
	    fprintf(IWR, "\n WARNING: EABS(3,%2d) has been modified because is greater than EMAX\n", M);
	  }
	  if(EABS[0][M] < EMIN)
	  {
	    EMIN = EABS[0][M];
	  }
	  if(EABS[1][M] < EMIN)
	  {
	    EMIN = EABS[1][M];
	  }
	  if(EABS[2][M] < EMIN)
	  {
	    EMIN = EABS[2][M];
	  }
  }
      
  if(EMIN < (double)MINEGRID){ EMIN = (double)MINEGRID;}

  fprintf(IWR, "\n EMIN =%11.4E eV,  EMAX =%11.4E eV\n", EMIN, EMAX);

  if(EMAX < EMIN+(double)MINDIMGRID)
  {
    fprintf(IWR, "\n ERROR: The energy interval between EMIN and EMAX is less than %.5E eV\n", MINDIMGRID);
	  ErrorFunction(1004);return;
  }
      
  if(NMAT > MAXMAT)
  {
	  fprintf(IWR, "\n ERROR: You are using %d materials and the maximum material number is %d\n", NMAT, MAXMAT);
	  ErrorFunction(1005);return;
  }
      
  if(INFO > 2){ fprintf(IWR, "\n NOTE: 1 mtu = 1 g/cm**2\n");}
		     
  EGRID(EMIN,EMAX);  // Defines the simulation energy grid.
  ESIa0();  // Initialises electron impact ionisation routines.
  PSIa0();  // Initialises positron impact ionisation routines.
  GPHa0();  // Initialises photoelectric routines.
  RELAX0();  // Initialises atomic relaxation routines.
  RNDG30();  // Initialises the Gaussian sampling routine.
  if(IRETRN != 0){ return;}

  for(int M = 0; M < NMAT; M++)
  {
  	if(M == 0){ LIT[0]='s'; LIT[1]='t';}
	  if(M == 1){ LIT[0]='n'; LIT[1]='d';}
	  if(M == 2){ LIT[0]='r'; LIT[1]='d';}
	  if(M > 2){ LIT[0]='t'; LIT[1]='h';}
	  LIT[2]='\0';
	  fprintf(IWR, "\n\n **********************\n **  %2d%s material   **\n **********************\n", M+1, LIT);
	  
	  fprintf(IWR, "\n Material data file: %-20s\n", PMFILE[M]);
	  FILE *IRD = fopen(PMFILE[M],"r");
	  if(IRD == NULL)
	  {
	    fprintf(IWR, "Error: Material file %s could not be opened\n", PMFILE[M]);
	    ErrorFunction(1006);return;
	  }
	  
	  //  ****  Energy limits and thresholds.
	  
	  fprintf(IWR, "\n *** Simulation parameters:\n");
	  fprintf(IWR, "     Electron absorption energy =%11.4E eV\n", EABS[0][M]);
	  fprintf(IWR, "       Photon absorption energy =%11.4E eV\n", EABS[1][M]);
	  fprintf(IWR, "     Positron absorption energy =%11.4E eV\n", EABS[2][M]);

	  if(fabs(C1[M]) < 0.2){C1[M] = fabs(C1[M]);}
	  else{C1[M] = 0.2;}
	  
	  if(fabs(C2[M]) < 0.2){C2[M] = fabs(C2[M]);}
	  else{C2[M] = 0.2;}
	  
	  if(fabs(WCC[M]) < EMAX){WCC[M] = fabs(WCC[M]);}
	  else{WCC[M] = EMAX;}

	  if(WCC[M] > EABS[0][M]){ WCC[M] = EABS[0][M];}
	  if(WCR[M] > EABS[1][M]){ WCR[M] = EABS[1][M];}
	  if(WCR[M] < 0.0){fprintf(IWR, "*** Warning: soft radiative losses are switched off in material number %3d\n", M);}

	  fprintf(IWR, "      C1 =%11.4E,       C2 =%11.4E\n     WCC =%11.4E eV,   WCR =%11.4E eV\n\n", C1[M], C2[M], WCC[M], (WCR[M]>10.0 ? WCR[M]:10.0));

	  if(EABS[0][M] < EABS[1][M]){ ECUTR[M] = EABS[0][M];}
	  else{ ECUTR[M] = EABS[1][M];}
     
	  PEMATR(M+1,IRD,IWR,INFO);
	  fclose(IRD);
  }
      
      //  ****  Restore the user values of EABS(.), to avoid inconsistencies
      //        in the main program.

  for(int M=0; M < NMAT; M++)
  {
  	EABS[0][M]=EABS0[0][M];
	  EABS[1][M]=EABS0[1][M];
	  EABS[2][M]=EABS0[2][M];
  }
}

