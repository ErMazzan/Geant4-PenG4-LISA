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

#ifndef PENELOPE_LOCAL_H
#define PENELOPE_LOCAL_H

#include <stdio.h>
#include <math.h>

void JUMP(double&, double&);
void KNOCK(double&, int&);
void EELa(double, double, double, double&);
void EINa(double, double, double&, double&, double&, double&, double&, int, int&);
void PINa(double, double, double&, double&, double&, double&, double&, int, int&);
void ESIa(double, double, double&, double&, double&, double&, double&, int, int&, int&);
void PSIa(double, double, double&, double&, double&, double&, double&, int, int&, int&);
void EBRa(double&, double&, int&);
void EBRaA(double&, double&, double&, int&);
void PANaR(double&);
void PANa(double&, double&, double&, double&, double&, int&);
void GRAa(double&, double&, int&, int&);
void GCOa(double, double&, double&, double&, double&, double&, int, int&, int&);
void GPHa(double&, int&, int&);
void SAUTER(double&, double&);
void GPPa(double&, double&, double&, double&, int&, int&);
void SCHIFF(double&, double&, double&);
void RELAX(int&, int&);
void PELd(double&, double&);
double RNDG3();
void CLEANS();
void START();
void DIRECT(double&, double&, double&, double&, double&);
void DIRPOL(double, double&, double, double&, double&, double&, double&, double&, double&);
void STORES(double, double, double, double, double, double, double, double, int, int ILBI[5], int);
void SECPAR(int&);
void EIMFP(int);
void PIMFP(int);
void GIMFP();
void EAUX();
void PAUX();
void GAUX();
void EELd(double&, double&);
void RAND0(int);
double RAND(double);

#endif
