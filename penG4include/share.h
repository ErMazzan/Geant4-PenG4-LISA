/*
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C    PPPPP   EEEEEE  N    N  EEEEEE  L        OOOO   PPPPP   EEEEEE    C
C    P    P  E       NN   N  E       L       O    O  P    P  E         C
C    P    P  E       N N  N  E       L       O    O  P    P  E         C
C    PPPPP   EEEE    N  N N  EEEE    L       O    O  PPPPP   EEEE      C
C    P       E       N   NN  E       L       O    O  P       E         C
C    P       EEEEEE  N    N  EEEEEE  LLLLLL   OOOO   P       EEEEEE    C
C                                                                      C
C                                                   (version 2018).    C
C                                                                      C
C  Subroutine package for Monte Carlo simulation of coupled electron-  C
C  photon transport in homogeneous media.                              C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  PENELOPE/PENGEOM (version 2018)                                     C
C  Copyright (c) 2001-2018                                             C
C  Universitat de Barcelona                                            C
C                                                                      C
C  Permission to use, copy, modify, distribute and sell this software  C
C  and its documentation for any purpose is hereby granted without     C
C  fee, provided that the above copyright notice appears in all        C
C  copies and that both that copyright notice and this permission      C
C  notice appear in all supporting documentation. The Universitat de   C
C  Barcelona makes no representations about the suitability of this    C
C  software for any purpose. It is provided "as is" without express    C
C  or implied warranty.                                                C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  The convention used to name the interaction simulation subroutines is
C  the following:
C  - The first letter indicates the particle (E for electrons, P for
C    positrons, G for photons).
C  - The second and third letters denote the interaction mechanism
C    (EL for elastic, IN for inelastic, SI for inner-shell ionisation,
C    BR for bremsstrahlung, AN for annihilation, RA for Rayleigh, CO for
C    Compton, PH for photoelectric and PP for pair production).
C  - The fourth (lowercase) letter indicates the theoretical model
C    used to describe the interactions. This serves to distinguish the
C    default model (denoted by the letter 'a') from alternative models.
C  - The random sampling routines have four-letter names. Auxiliary
C    routines, which perform specific calculations, have longer names,
C    with the fifth and subsequent letters and/or numbers indicating
C    the kind of calculation (T for total x-section, D for differen-
C    tial x-section) or action (W for write data in a file, R for read
C    data from a file, I for initialisation of simulation algorithm).
C
C  The present subroutines may print warning and error messages in
C  the I/O unit 26. This is the default output unit in the example main
C  programs PENCYL and PENMAIN.
C
C  Subroutine PEMATW connects files to the I/O units 3 (input) and 7
C  (output). However, this does not conflict with the main program,
C  because PEMATW is not invoked during simulation. This subroutine is
C  used only by the program MATERIAL, to generate material data files.
C
C  Subroutine PEINIT connects material definition files to the input
C  unit 3; this unit is closed before returning to the main program.
C
*/

//  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#ifndef PENELOPE_SHARE_H
#define PENELOPE_SHARE_H

#include <stdio.h>
#include <math.h>

void PEINIT(double, unsigned int, FILE*, int, char PMFILE[MAXMAT][81]);
void EGRID(double, double);
void PEMATR(unsigned int, FILE*, FILE*, int);
void  PEMATW(int, char MFNAME[20]);
double PRANGE(double, int, int);
double AVNCOL(double, int, int, int);
double PHMFP(double, int, int, int);
double RYIELD(double, int, int);
void EELaR(int, FILE*, FILE*, int);
void EELa0(double&, double&, double&, double&, double&, double&, double&, double&, double&);
void EELaW(int&, FILE*);
void EINaT(double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, int);
void EINaT1(double&, double&, double&, double, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&);
//double EINaDS(double);
void PINaT(double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, int);
void PINaT1(double&, double&, double&, double, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&);
//double PINaDS(double);
//void ESIb(int&, int&, int&);
void ESIaR(int, FILE*, FILE*, int);
void ESIa0();
void ESIaW(int&, FILE*);
//void PSIb(int&, int&, int&);
void PSIaR(int, FILE*, FILE*, int);
void PSIa0();
void PSIaW(int&, FILE*);
void EBRaR(double&, int, FILE*, FILE*, int);
void EBRaW(int&, FILE*);
void EBRaT(double&, double&, double&, double&, double&, double&, double&, int);
void PBRaT(double&, double&, double&, double&, double&, double&, double&, int);
double RLMOM(double*, double*, double, int, int);
void RLPAC(double*, double*, double*, int);
void BRaAR(int, FILE*, FILE*, int);
void BRaAW(double&, FILE*);
void PANaT(double&, double&);
void GRAbT(double&, double&, int&);
//double GRAaD(double);
//double GRAaF2(double);
void GRAaTI(double&, double&, int);
void GRAaR(int, FILE*, FILE*, int&);
void GRAaW(int&, FILE*);
void GCOaT(double&, double&, int&);
//double GCOaD(double);
void GPHaT(double&, double&, int&);
void GPHaR(int, FILE*, FILE*, int);
void GPHa0();
void GPHaW(int&, FILE*);
void GPPa0(int);
void GPPaW(double*, double*, double*, int&, int&);
void PEILB4(int&, int&, int&, int&, int&);
void RELAXE(int&, int&, int&, int&, double&, double&);
void RELAX0();
void RELAXR(FILE*, FILE*, int);
void RELAXW(int&, FILE*);
void EELdW(int&, FILE*);
void EELdR(int, FILE*, FILE*, int&);
void ELINIT(int*, double*, int&);
void DCSEL0(double&, int);
//double DCSEL(double);
void MERGE2(double*, double*, double*, double*, double*, double*, int&, int&, int&);
void SORT2(double*, double*, int&);
void SPLINE(double*, double*, double*, double*, double*, double*, double, double, int);
//void FINDI(double*, double, int, int&);
void SINTEG(double*, double*, double*, double*, double*, double&, double&, double&, int);
void SLAG6(double&, double*, double*, int);
double SUMGA(FCT_SUMGA, double, double, double);
double RMOMX(double*, double*, double, double, int, int);
void RNDG30();
//double RNDG3F(double);
int IRND(double*, int*, int&);
double RITA();
double RITAI();
void IRND0(double*, double*, int*, int&);
void RITA0(PDF_RITA, double, double, int&, int&, double&, int, int&);
void RITAI0(PDF_RITA, double, double, int, int, double&, int, int&);
void RITAV(double&, double&, double&);
void RITAM(double, double, double&, double&, double&);

void PEMATZ(int&,FILE*);
void PEMATC(char*,int,int*,double*,double,FILE*);

#endif
