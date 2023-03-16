 
// En aquest fitxer guardarem les variables dels commons i moduls amb el mateix nom que el namespace

#ifndef COMMON_DEFS_H
#define COMMON_DEFS_H 1

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

//Llistat de parametres i variables globals
#define MAXMAT  20 //  ****  Maximum number of materials in the geometry.
#define NMS 1000  //  ****  Size of the secondary stack (maximum number of particles that can be stored).
#define NEGP  200 //  ****  Energy interpolation, number of grid points.
#define MINEABSE 50.0 // **** Minima energia d'absorcio en electrons (eV)
#define MINEABSPh 50.0 // **** Minima energia d'absorcio en fotons (eV)
#define MINEABSPo 50.0 // **** Minima energia d'absorcio en positrons (eV)
#define MINEGRID 50.0 // **** Minima energia del "grid" (eV). Deuria ser la minima de les energies d'absorcio.
#define MINDIMGRID 10.0 // **** Minima distancia entre el valor maxim i minim del "grid" (eV)

#define A0B 5.2917721092E-9  // Bohr radius (cm)
#define HREV 27.21138505  // Hartree energy (eV)
#define AVOG 6.02214129E23  // Avogadro's number
#define SL 137.035999074  // Speed of light (1/alpha)
#define PI 3.1415926535897932 // Valor de PI
#define REV 5.10998928E5  // Electron rest energy (eV)
#define SL 137.035999074  // Speed of light (1/alpha)

//Aquests son parametres dins de commons:
#define NO 512 //****  E/P inelastic collisions.
#define NRP 8000 //  ****  Inner-shell ionisation by electron and positron impact.
#define NOCO 512 //  ****  Compton scattering.
#define NDIM 12000 //  ****  Photon simulation tables.
#define NBW 32 //  ****  Bremsstrahlung emission.
#define NBE 57
#define NEGP1 12000 //  ****  'Standard' energy grid.

//  ****  Rayleigh scattering.
#define NQ 250
#define NEX 1024
#define NTP 12000 //  ****  Photoelectric cross sections.

#define NP 150
#define NP2 300 //(NP2 = NP+NP)
//#define NPM1 149 //(NPM1 = NP-1)
#define NPMINUS1 149 //(NPM1 = NP-1)

#define NOM 1000
#define ELRAD 2.8179403267E-13  // Class. electron radius (cm)

//#define NM 512 //MULTIDEFINED !!!!!
#define NEM 10000

#define NRX 60000

#define HBAR 6.58211928E-16

//define NE 96
#define NA 606

#define NR 128

//SUMGA necessita passar una funció com a paràmetre. Aquest és el typedef que cal.
 typedef double (* FCT_SUMGA)(double);
//RITA0 necessita passar una funció com a paràmetre. Aquest és el typedef que cal.
 typedef double (* PDF_RITA)(double);

struct CEGRID
{
  // ****  Energy grid and interpolation constants for the current energy.
  double EMIN, EL, EU, ET[NEGP], DLEMP[NEGP], DLEMP1, DLFC, XEL, XE, XEK;
  int KE;
};

//static functions
double PINaDS(double);
double EINaDS(double);
double GRAaF2(double);
double DCSEL(double);
double RNDG3F(double);
double GCOaD(double);
double GRAaD(double);

void FINDI(double*, double, int, int&);
void ErrorFunction(int);

#endif
