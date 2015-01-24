#define RNDFILE "example/water/random.in"
#line 228 "/home/jyy/Dev/splash2/codes/null_macros/c.m4.null.POSIX"

#line 1 "bndry.C"
/*************************************************************************/
/*                                                                       */
/*  Copyright (c) 1994 Stanford University                               */
/*                                                                       */
/*  All rights reserved.                                                 */
/*                                                                       */
/*  Permission is given to use, copy, and modify this software for any   */
/*  non-commercial purpose as long as this copyright notice is not       */
/*  removed.  All other uses, including redistribution in whole or in    */
/*  part, are forbidden without prior written permission.                */
/*                                                                       */
/*  This software is provided with absolutely no warranty and no         */
/*  support.                                                             */
/*                                                                       */
/*************************************************************************/


#line 17
#include <pthread.h>
#line 17
#include <sys/time.h>
#line 17
#include <unistd.h>
#line 17
#include <stdlib.h>
#line 17
extern pthread_t PThreadTable[];
#line 17

#include "mdvar.h"
#include "parameters.h"
#include "mddata.h"
#include "split.h"
#include "global.h"

/* this routine puts the molecules back inside the box if they are out */
void BNDRY(long ProcID)
{
    long mol, dir;
    double *extra_p;

    /* for each molecule */
    for (mol = StartMol[ProcID]; mol < StartMol[ProcID+1]; mol++) {
        /* for each direction */
        for ( dir = XDIR; dir <= ZDIR; dir++ ) {
            extra_p = VAR[mol].F[DISP][dir];
            /* if the oxygen atom is out of the box */
            if (extra_p[O] > BOXL) {
                /* move all three atoms back in the box */
                extra_p[H1] -= BOXL;
                extra_p[O]  -= BOXL;
                extra_p[H2] -= BOXL;
            }
            else if (extra_p[O] < 0.00) {
                extra_p[H1] += BOXL;
                extra_p[O]  += BOXL;
                extra_p[H2] += BOXL;
            }
        } /* for dir */
    } /* for mol */
} /* end of subroutine BNDRY */
#line 228 "/home/jyy/Dev/splash2/codes/null_macros/c.m4.null.POSIX"

#line 1 "cnstnt.C"
/*************************************************************************/
/*                                                                       */
/*  Copyright (c) 1994 Stanford University                               */
/*                                                                       */
/*  All rights reserved.                                                 */
/*                                                                       */
/*  Permission is given to use, copy, and modify this software for any   */
/*  non-commercial purpose as long as this copyright notice is not       */
/*  removed.  All other uses, including redistribution in whole or in    */
/*  part, are forbidden without prior written permission.                */
/*                                                                       */
/*  This software is provided with absolutely no warranty and no         */
/*  support.                                                             */
/*                                                                       */
/*************************************************************************/


#line 17
#include <pthread.h>
#line 17
#include <sys/time.h>
#line 17
#include <unistd.h>
#line 17
#include <stdlib.h>
#line 17
extern pthread_t PThreadTable[];
#line 17

#include <stdio.h>
#include "water.h"
#include "wwpot.h"
#include "cnst.h"
#include "frcnst.h"
#include "fileio.h"
#include "parameters.h"
#include "global.h"

/* set up some constants
 * N : NORDER + 1 = 7 for a sixth-order method
 * C : DIMENSION C(N,N)
 */
void CNSTNT(long N, double *C)
{
    long NN,N1,K1;
    double TN,TK,CM;

    /* molecular constants for water in angstrom, radian, and a.m.u. */

    NATOMS = 3;
    ROH = 0.9572;
    ROHI = ONE/ROH;
    ROHI2 = ROHI*ROHI;
    ANGLE = 1.824218;
    OMAS = 15.99945;
    HMAS = 1.007825;
    WTMOL = OMAS+TWO*HMAS;

    /* units used to scale variables (in c.g.s.) */

    UNITT = 1.0e-15;
    UNITL = 1.0e-8;
    UNITM = 1.6605655e-24;
    BOLTZ = 1.380662e-16;
    AVGNO = 6.022045e23;

    /* force constants scaled (divided) by (UNITM/UNITT**2) */

    FC11 =  0.512596;
    FC33 =  0.048098;
    FC12 = -0.005823;
    FC13 =  0.016452;
    FC111 = -0.57191;
    FC333 = -0.007636;
    FC112 = -0.001867;
    FC113 = -0.002047;
    FC123 = -0.03083;
    FC133 = -0.0094245;
    FC1111 =  0.8431;
    FC3333 = -0.00193;
    FC1112 = -0.0030;
    FC1122 =  0.0036;
    FC1113 = -0.012;
    FC1123 =  0.0060;
    FC1133 = -0.0048;
    FC1233 =  0.0211;
    FC1333 =  0.006263;

    /* water-water interaction parameters */

    QQ = 0.07152158;
    A1 = 455.313100;
    B1 = 5.15271070;
    A2 = 0.27879839;
    B2 = 2.76084370;
    A3 = 0.60895706;
    B3 = 2.96189550;
    A4 = 0.11447336;
    B4 = 2.23326410;
    CM = 0.45682590;
    AB1 = A1*B1;
    AB2 = A2*B2;
    AB3 = A3*B3;
    AB4 = A4*B4;
    C1 = ONE-CM;
    C2 = 0.50*CM;
    QQ2 = 2.00*QQ;
    QQ4 = 2.00*QQ2;

    /*  calculate the coefficients of taylor series expansion */
    /*     for F(X), F"(X), F""(X), ...... (with DELTAT**N/N] included) */
    /*     in C(1,1),..... C(1,2),..... C(1,3),....... */

    C[1] = ONE;
    for (N1=2;N1<=N;N1++) {
        NN = N1-1;
        TN = NN;
        C[N1] = ONE;
        TK = ONE;
        for (K1=2;K1<=N1;K1++) {
            C[(K1-1)*N+NN] = C[(K1-2)*N+NN+1]*TN/TK;
            NN = NN-1;
            TN = TN-ONE;
            TK = TK+ONE;
        }
    }


    /* predictor-corrector constants for 2nd order differential equation */

    PCC[2] = ONE;
    N1 = N-1;
    switch(N1) {
    case 1:
    case 2:
        fprintf(six,"***** ERROR: THE ORDER HAS TO BE GREATER THAN 2 ****");
        break;
    case 3:
        PCC[0] = ONE/6.00;
        PCC[1] = FIVE/6.00;
        PCC[3] = ONE/3.00;
        break;
    case 4:
        PCC[0] = (double) 19.00/120.00;
        PCC[1] = (double) 3.00/4.00;
        PCC[3] = ONE/2.00;
        PCC[4] = ONE/12.00;
        break;
    case 5:
        PCC[0] = (double) 3.00/20.00;
        PCC[1] = (double) 251.00/360.00;
        PCC[3] = (double) 11.00/18.00;
        PCC[4] = ONE/6.00;
        PCC[5] = ONE/60.00;
        break;
    case 6:
        PCC[0] = (double) 863.00/6048.00;
        PCC[1] = (double) 665.00/1008.00;
        PCC[3] = (double) 25.00/36.00;
        PCC[4] = (double) 35.00/144.00;
        PCC[5] = ONE/24.00;
        PCC[6] = ONE/360.00;
        break;
    case 7:
        PCC[0] = (double) 275.00/2016.00;
        PCC[1] = (double) 19087.00/30240.00;
        PCC[3] = (double) 137.00/180.00;
        PCC[4] = FIVE/16.00;
        PCC[5] = (double) 17.00/240.00;
        PCC[6] = ONE/120.00;
        PCC[7] = ONE/2520.00;
        break;
    default:
        break;
    }
}           /* end of subroutine CNSTNT */
#line 228 "/home/jyy/Dev/splash2/codes/null_macros/c.m4.null.POSIX"

#line 1 "cshift.C"
/*************************************************************************/
/*                                                                       */
/*  Copyright (c) 1994 Stanford University                               */
/*                                                                       */
/*  All rights reserved.                                                 */
/*                                                                       */
/*  Permission is given to use, copy, and modify this software for any   */
/*  non-commercial purpose as long as this copyright notice is not       */
/*  removed.  All other uses, including redistribution in whole or in    */
/*  part, are forbidden without prior written permission.                */
/*                                                                       */
/*  This software is provided with absolutely no warranty and no         */
/*  support.                                                             */
/*                                                                       */
/*************************************************************************/


#line 17
#include <pthread.h>
#line 17
#include <sys/time.h>
#line 17
#include <unistd.h>
#line 17
#include <stdlib.h>
#line 17
extern pthread_t PThreadTable[];
#line 17

#include <math.h>
#include "global.h"

/* return the value of a with the same sign as b */

#define sign(a,b)  (b < 0 ) ? ( (a < 0) ? a : -a) : ( (a < 0) ? -a : a)

  /* compute some relevant distances between the two input molecules to
     this routine. if they are greater than the cutoff radius, compute
     these distances as if one of the particles were at its mirror image
     (periodic boundary conditions).
     used by the intermolecular interactions routines */
void CSHIFT(double *XA, double *XB, double XMA, double XMB, double *XL, double BOXH, double BOXL)
{

    long I;

    XL[0] = XMA-XMB;
    XL[1] = XMA-XB[0];
    XL[2] = XMA-XB[2];
    XL[3] = XA[0]-XMB;
    XL[4] = XA[2]-XMB;
    XL[5] = XA[0]-XB[0];
    XL[6] = XA[0]-XB[2];
    XL[7] = XA[2]-XB[0];
    XL[8] = XA[2]-XB[2];
    XL[9] = XA[1]-XB[1];
    XL[10] = XA[1]-XB[0];
    XL[11] = XA[1]-XB[2];
    XL[12] = XA[0]-XB[1];
    XL[13] = XA[2]-XB[1];

    /* go through all 14 distances computed */
    for (I = 0; I <  14; I++) {
        /* if the value is greater than the cutoff radius */
        if (fabs(XL[I]) > BOXH) {
            XL[I]  =  XL[I] -(sign(BOXL,XL[I]));
        }
    } /* for */
} /* end of subroutine CSHIFT */

#line 228 "/home/jyy/Dev/splash2/codes/null_macros/c.m4.null.POSIX"

#line 1 "initia.C"
/*************************************************************************/
/*                                                                       */
/*  Copyright (c) 1994 Stanford University                               */
/*                                                                       */
/*  All rights reserved.                                                 */
/*                                                                       */
/*  Permission is given to use, copy, and modify this software for any   */
/*  non-commercial purpose as long as this copyright notice is not       */
/*  removed.  All other uses, including redistribution in whole or in    */
/*  part, are forbidden without prior written permission.                */
/*                                                                       */
/*  This software is provided with absolutely no warranty and no         */
/*  support.                                                             */
/*                                                                       */
/*************************************************************************/


#line 17
#include <pthread.h>
#line 17
#include <sys/time.h>
#line 17
#include <unistd.h>
#line 17
#include <stdlib.h>
#line 17
extern pthread_t PThreadTable[];
#line 17

#include "math.h"
#include "stdio.h"
#include "mdvar.h"
#include "water.h"
#include "cnst.h"
#include "fileio.h"
#include "parameters.h"
#include "mddata.h"
#include "split.h"
#include "global.h"

void INITIA()
{
    /*   this routine initializes the positions of the molecules along
         a regular cubical lattice, and randomizes the initial velocities of
         the atoms.  The random numbers used in the initialization of velocities
         are read from the file random.in, which must be in the current working
         directory in which the program is run  */

    FILE *random_numbers;       /* points to input file containing
                                   pseudo-random numbers for initializing
                                   velocities */
    double XMAS[4], XT[4], YT[4], Z;
    double SUX, SUY, SUZ, SUMX, SUMY, SUMZ, FAC;
    long mol=0;
    long atom=0;
    long deriv;

    random_numbers = fopen(RNDFILE,"r");
    if (random_numbers == NULL) {
        fprintf(stderr,"Error in opening file random.in\n");
        fflush(stderr);
        exit(-1);
    }

    XMAS[1]=sqrt(OMAS*HMAS);
    XMAS[0]=HMAS;
    XMAS[2]=HMAS;

    /* .....assign positions */
    {
        double NS = pow((double) NMOL, 1.0/3.0) - 0.00001;
        double XS = BOXL/NS;
        double ZERO = XS * 0.50;
        double WCOS = ROH * cos(ANGLE * 0.5);
        double WSIN = ROH * sin(ANGLE * 0.5);
        long i,j,k;

        printf("\nNS = %.16f\n",NS);
        printf("BOXL = %10f\n",BOXL);
        printf("CUTOFF = %10f\n",CUTOFF);
        printf("XS = %10f\n",XS);
        printf("ZERO = %g\n",ZERO);
        printf("WCOS = %f\n",WCOS);
        printf("WSIN = %f\n",WSIN);
        fflush(stdout);

#ifdef RANDOM
        /* if we want to initialize to a random distribution of displacements
           for the molecules, rather than a distribution along a regular lattice
           spaced according to intermolecular distances in water */
        srand(1023);
        for (i = 0; i < NMOL; i++) {
            VAR[mol].F[DISP][XDIR][O] = xrand(0, BOXL);
            VAR[mol].F[DISP][XDIR][H1] = VAR[mol].F[DISP][XDIR][O] + WCOS;
            VAR[mol].F[DISP][XDIR][H2] = VAR[mol].F[DISP][XDIR][H1];
            VAR[mol].F[DISP][YDIR][O] = xrand(0, BOXL);
            VAR[mol].F[DISP][YDIR][H1] = VAR[mol].F[DISP][YDIR][O] + WSIN;
            VAR[mol].F[DISP][YDIR][H2] = VAR[mol].F[DISP][YDIR][O] - WSIN;
            VAR[mol].F[DISP][ZDIR][O] = xrand(0, BOXL);
            VAR[mol].F[DISP][ZDIR][H1] = VAR[mol].F[DISP][ZDIR][O];
            VAR[mol].F[DISP][ZDIR][H2] = VAR[mol].F[DISP][ZDIR][O];
        }
#else
        /* not random initial placement, but rather along a regular
           lattice.  This is the default and the prefered initialization
           since random does not necessarily make sense from the viewpoint
           of preserving bond distances */

        fprintf(six, "***** NEW RUN STARTING FROM REGULAR LATTICE *****\n");
        fflush(six);
        XT[2] = ZERO;
        mol = 0;
        for (i=0; i < NS; i+=1) {
            XT[1]=XT[2]+WCOS;
            XT[3]=XT[1];
            YT[2]=ZERO;
            for (j=0; j < NS; j+=1) {
                YT[1]=YT[2]+WSIN;
                YT[3]=YT[2]-WSIN;
                Z=ZERO;
                for (k = 0; k < NS; k++) {
                    for (atom = 0; atom < NATOMS; atom +=1) {
                        VAR[mol].F[DISP][XDIR][atom] = XT[atom+1];
                        VAR[mol].F[DISP][YDIR][atom] = YT[atom+1];
                        VAR[mol].F[DISP][ZDIR][atom] = Z;
                    }
                    mol += 1;
                    Z=Z+XS;
                }
                YT[2]=YT[2]+XS;
            }
            XT[2]=XT[2]+XS;
        }

        if (NMOL != mol) {
            printf("Lattice init error: total mol %ld != NMOL %ld\n", mol, NMOL);
            exit(-1);
        }
#endif
    }

    /* ASSIGN RANDOM MOMENTA */
    fscanf(random_numbers,"%lf",&SUX);

    SUMX=0.0;
    SUMY=0.0;
    SUMZ=0.0;
    /*   read pseudo-random numbers from input file random.in */
    for (mol = 0; mol < NMOL; mol++) {
        for (atom = 0; atom < NATOMS; atom++) {
            fscanf(random_numbers,"%lf",&VAR[mol].F[VEL][XDIR][atom]);
            fscanf(random_numbers,"%lf",&VAR[mol].F[VEL][YDIR][atom]);
            fscanf(random_numbers,"%lf",&VAR[mol].F[VEL][ZDIR][atom]);
            SUMX = SUMX + VAR[mol].F[VEL][XDIR][atom];
            SUMY = SUMY + VAR[mol].F[VEL][YDIR][atom];
            SUMZ = SUMZ + VAR[mol].F[VEL][ZDIR][atom];
            for (deriv = ACC; deriv < MAXODR; deriv++) {
                VAR[mol].F[deriv][XDIR][atom] = 0.0;
                VAR[mol].F[deriv][YDIR][atom] = 0.0;
                VAR[mol].F[deriv][ZDIR][atom] = 0.0;
            }
        } /* atoms */
    } /* molecules */

    /* find average momenta per atom */
    SUMX=SUMX/(NATOMS*NMOL);
    SUMY=SUMY/(NATOMS*NMOL);
    SUMZ=SUMZ/(NATOMS*NMOL);

    /*  find normalization factor so that <k.e.>=KT/2  */
    SUX=0.0;
    SUY=0.0;
    SUZ=0.0;
    for (mol = 0; mol < NMOL; mol++) {
        SUX = SUX + (pow( (VAR[mol].F[VEL][XDIR][H1] - SUMX),2.0)
                     +pow( (VAR[mol].F[VEL][XDIR][H2] - SUMX),2.0))/HMAS
                         +pow( (VAR[mol].F[VEL][XDIR][O]  - SUMX),2.0)/OMAS;

        SUY = SUY + (pow( (VAR[mol].F[VEL][YDIR][H1] - SUMY),2.0)
                     +pow( (VAR[mol].F[VEL][YDIR][H2] - SUMY),2.0))/HMAS
                         +pow( (VAR[mol].F[VEL][YDIR][O]  - SUMY),2.0)/OMAS;

        SUZ = SUZ + (pow( (VAR[mol].F[VEL][ZDIR][H1] - SUMZ),2.0)
                     +pow( (VAR[mol].F[VEL][ZDIR][H2] - SUMZ),2.0))/HMAS
                         +pow( (VAR[mol].F[VEL][ZDIR][O]  - SUMZ),2.0)/OMAS;
    }
    FAC=BOLTZ*TEMP*NATMO/UNITM * pow((UNITT*TSTEP/UNITL),2.0);
    SUX=sqrt(FAC/SUX);
    SUY=sqrt(FAC/SUY);
    SUZ=sqrt(FAC/SUZ);

    /* normalize individual velocities so that there are no bulk
       momenta  */
    XMAS[1]=OMAS;
    for (mol = 0; mol < NMOL; mol++) {
        for (atom = 0; atom < NATOMS; atom++) {
            VAR[mol].F[VEL][XDIR][atom] = ( VAR[mol].F[VEL][XDIR][atom] -
                                           SUMX) * SUX/XMAS[atom];
            VAR[mol].F[VEL][YDIR][atom] = ( VAR[mol].F[VEL][YDIR][atom] -
                                           SUMY) * SUY/XMAS[atom];
            VAR[mol].F[VEL][ZDIR][atom] = ( VAR[mol].F[VEL][ZDIR][atom] -
                                           SUMZ) * SUZ/XMAS[atom];
        } /* for atom */
    } /* for mol */

    fclose(random_numbers);

} /* end of subroutine INITIA */

/*
 * XRAND: generate floating-point random number.
 */

double xrand(double xl, double xh)
{
    double x;

    x=(xl + (xh - xl) * ((double) rand()) / 2147483647.0);
    return (x);
}
#line 228 "/home/jyy/Dev/splash2/codes/null_macros/c.m4.null.POSIX"

#line 1 "interf.C"
/*************************************************************************/
/*                                                                       */
/*  Copyright (c) 1994 Stanford University                               */
/*                                                                       */
/*  All rights reserved.                                                 */
/*                                                                       */
/*  Permission is given to use, copy, and modify this software for any   */
/*  non-commercial purpose as long as this copyright notice is not       */
/*  removed.  All other uses, including redistribution in whole or in    */
/*  part, are forbidden without prior written permission.                */
/*                                                                       */
/*  This software is provided with absolutely no warranty and no         */
/*  support.                                                             */
/*                                                                       */
/*************************************************************************/


#line 17
#include <pthread.h>
#line 17
#include <sys/time.h>
#line 17
#include <unistd.h>
#line 17
#include <stdlib.h>
#line 17
extern pthread_t PThreadTable[];
#line 17

#include "math.h"
#include "mdvar.h"
#include "water.h"
#include "wwpot.h"
#include "cnst.h"
#include "parameters.h"
#include "mddata.h"
#include "split.h"
#include "global.h"

double ****PFORCES;

/* in this version of interf, a private force array is maintained */
/* for every process.  A process computes interactions into its   */
/* private force array, and later updates the shared destination  */
/* array with locks, but only updates those locations that it     */
/* computed something for                                         */

void INTERF(long DEST, double *VIR, long ProcID)
{
    /* This routine gets called both from main() and from mdmain().
       When called from main(), it is used to estimate the initial
       accelerations by computing intermolecular forces.  When called
       from mdmain(), it is used to compute intermolecular forces.
       The parameter DEST specifies whether results go into the
       accelerations or the forces. Uses routine UPDATE_FORCES in this
       file, and routine CSHIFT in file cshift.U */
    /*
      .....this routine calculates inter-molecular interaction forces
      the distances are arranged in the order  M-M, M-H1, M-H3, H1-M,
      H3-M, H1-H3, H1-H1, H3-H1, H3-H3, O-O, O-H1, O-H3, H1-O, H3-O,
      where the M are "centers" of the molecules.
      */

    long mol, comp, dir, icomp;
    long comp_last, half_mol;
    long    KC, K;
    double YL[15], XL[15], ZL[15], RS[15], FF[15], RL[15]; /* per-
                                                              interaction arrays that hold some computed distances */
    double  FTEMP;
    double LVIR = 0.0;
    double *temp_p;

    { /* initialize PFORCES array */

        long ct1,ct2,ct3;

        for (ct1 = 0; ct1<NMOL; ct1++)
            for (ct2 = 0; ct2<NDIR; ct2++)
                for (ct3 = 0; ct3<NATOM; ct3++)
                    PFORCES[ProcID][ct1][ct2][ct3] = 0;

    }

    half_mol = NMOL/2;
    for (mol = StartMol[ProcID]; mol < StartMol[ProcID+1]; mol++) {
        comp_last = mol + half_mol;
        if (NMOL%2 == 0) {
            if ((half_mol <= mol) && (mol%2 == 0)) {
                comp_last--;
            }
            if ((mol < half_mol) && (comp_last%2 == 1)) {
                comp_last--;
            }
        }
        for (icomp = mol+1; icomp <= comp_last; icomp++) {
            comp = icomp;
            if (comp > NMOL1) comp = comp%NMOL;

            /*  compute some intermolecular distances */

            CSHIFT(VAR[mol].F[DISP][XDIR],VAR[comp].F[DISP][XDIR],
                   VAR[mol].VM[XDIR],VAR[comp].VM[XDIR],XL,BOXH,BOXL);
            CSHIFT(VAR[mol].F[DISP][YDIR],VAR[comp].F[DISP][YDIR],
                   VAR[mol].VM[YDIR],VAR[comp].VM[YDIR],YL,BOXH,BOXL);
            CSHIFT(VAR[mol].F[DISP][ZDIR],VAR[comp].F[DISP][ZDIR],
                   VAR[mol].VM[ZDIR],VAR[comp].VM[ZDIR],ZL,BOXH,BOXL);

            KC=0;
            for (K = 0; K < 9; K++) {
                RS[K]=XL[K]*XL[K]+YL[K]*YL[K]+ZL[K]*ZL[K];
                if (RS[K] > CUT2)
                    KC++;
            } /* for K */

            if (KC != 9) {
                for (K = 0; K < 14; K++)
                    FF[K]=0.0;
                if (RS[0] < CUT2) {
                    FF[0]=QQ4/(RS[0]*sqrt(RS[0]))+REF4;
                    LVIR = LVIR + FF[0]*RS[0];
                } /* if */
                for (K = 1; K < 5; K++) {
                    if (RS[K] < CUT2) {
                        FF[K]= -QQ2/(RS[K]*sqrt(RS[K]))-REF2;
                        LVIR = LVIR + FF[K]*RS[K];
                    } /* if */
                    if (RS[K+4] <= CUT2) {
                        RL[K+4]=sqrt(RS[K+4]);
                        FF[K+4]=QQ/(RS[K+4]*RL[K+4])+REF1;
                        LVIR = LVIR + FF[K+4]*RS[K+4];
                    } /* if */
                } /* for K */
                if (KC == 0) {
                    RS[9]=XL[9]*XL[9]+YL[9]*YL[9]+ZL[9]*ZL[9];
                    RL[9]=sqrt(RS[9]);
                    FF[9]=AB1*exp(-B1*RL[9])/RL[9];
                    LVIR = LVIR + FF[9]*RS[9];
                    for (K = 10; K < 14; K++) {
                        FTEMP=AB2*exp(-B2*RL[K-5])/RL[K-5];
                        FF[K-5]=FF[K-5]+FTEMP;
                        LVIR= LVIR+FTEMP*RS[K-5];
                        RS[K]=XL[K]*XL[K]+YL[K]*YL[K]+ZL[K]*ZL[K];
                        RL[K]=sqrt(RS[K]);
                        FF[K]=(AB3*exp(-B3*RL[K])-AB4*exp(-B4*RL[K]))/RL[K];
                        LVIR = LVIR + FF[K]*RS[K];
                    } /* for K */
                } /* if KC == 0 */

                UPDATE_FORCES(mol, comp, XL, YL, ZL, FF, ProcID);

            }  /* if KC != 9 */
        } /* for comp */
    } /* for mol */

    /*  accumulate the running sum from private
        per-interaction partial sums   */
    {pthread_mutex_lock(&(gl->InterfVirLock));};
    *VIR = *VIR + LVIR;
    {pthread_mutex_unlock(&(gl->InterfVirLock));};

    /* at the end of the above force-computation, comp_last */
    /* contains the number of the last molecule (no modulo) */
    /* that this process touched                            */

    if (comp_last > NMOL1) {
        for (mol = StartMol[ProcID]; mol < NMOL; mol++) {
            {pthread_mutex_lock(&gl->MolLock[mol % MAXLCKS]);};
            for ( dir = XDIR; dir  <= ZDIR; dir++) {
                temp_p = VAR[mol].F[DEST][dir];
                temp_p[H1] += PFORCES[ProcID][mol][dir][H1];
                temp_p[O]  += PFORCES[ProcID][mol][dir][O];
                temp_p[H2] += PFORCES[ProcID][mol][dir][H2];
            }
            {pthread_mutex_unlock(&gl->MolLock[mol % MAXLCKS]);};
        }
        comp = comp_last % NMOL;
        for (mol = 0; ((mol <= comp) && (mol < StartMol[ProcID])); mol++) {
            {pthread_mutex_lock(&gl->MolLock[mol % MAXLCKS]);};
            for ( dir = XDIR; dir  <= ZDIR; dir++) {
                temp_p = VAR[mol].F[DEST][dir];
                temp_p[H1] += PFORCES[ProcID][mol][dir][H1];
                temp_p[O]  += PFORCES[ProcID][mol][dir][O];
                temp_p[H2] += PFORCES[ProcID][mol][dir][H2];
            }
            {pthread_mutex_unlock(&gl->MolLock[mol % MAXLCKS]);};
        }
    }
    else{
        for (mol = StartMol[ProcID]; mol <= comp_last; mol++) {
            {pthread_mutex_lock(&gl->MolLock[mol % MAXLCKS]);};
            for ( dir = XDIR; dir  <= ZDIR; dir++) {
                temp_p = VAR[mol].F[DEST][dir];
                temp_p[H1] += PFORCES[ProcID][mol][dir][H1];
                temp_p[O]  += PFORCES[ProcID][mol][dir][O];
                temp_p[H2] += PFORCES[ProcID][mol][dir][H2];
            }
            {pthread_mutex_unlock(&gl->MolLock[mol % MAXLCKS]);};
        }
    }

    /* wait till all forces are updated */

    {
#line 191
	unsigned long	Error, Cycle;
#line 191
	long		Cancel, Temp;
#line 191

#line 191
	Error = pthread_mutex_lock(&(gl->InterfBar).mutex);
#line 191
	if (Error != 0) {
#line 191
		printf("Error while trying to get lock in barrier.\n");
#line 191
		exit(-1);
#line 191
	}
#line 191

#line 191
	Cycle = (gl->InterfBar).cycle;
#line 191
	if (++(gl->InterfBar).counter != (NumProcs)) {
#line 191
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 191
		while (Cycle == (gl->InterfBar).cycle) {
#line 191
			Error = pthread_cond_wait(&(gl->InterfBar).cv, &(gl->InterfBar).mutex);
#line 191
			if (Error != 0) {
#line 191
				break;
#line 191
			}
#line 191
		}
#line 191
		pthread_setcancelstate(Cancel, &Temp);
#line 191
	} else {
#line 191
		(gl->InterfBar).cycle = !(gl->InterfBar).cycle;
#line 191
		(gl->InterfBar).counter = 0;
#line 191
		Error = pthread_cond_broadcast(&(gl->InterfBar).cv);
#line 191
	}
#line 191
	pthread_mutex_unlock(&(gl->InterfBar).mutex);
#line 191
};

    /* divide final forces by masses */

    for (mol = StartMol[ProcID]; mol < StartMol[ProcID+1]; mol++) {
        for ( dir = XDIR; dir  <= ZDIR; dir++) {
            temp_p = VAR[mol].F[DEST][dir];
            temp_p[H1] = temp_p[H1] * FHM;
            temp_p[O]  = temp_p[O] * FOM;
            temp_p[H2] = temp_p[H2] * FHM;
        } /* for dir */
    } /* for mol */

}/* end of subroutine INTERF */

  /* from the computed distances etc., compute the
     intermolecular forces and update the force (or
     acceleration) locations */
void UPDATE_FORCES(long mol, long comp, double *XL, double *YL, double *ZL, double *FF, long ProcID)
{
    long K;
    double G110[3], G23[3], G45[3], TT1[3], TT[3], TT2[3];
    double GG[15][3];
    double *tx_p, *ty_p, *tz_p;

    /*   CALCULATE X-COMPONENT FORCES */
    for (K = 0; K < 14; K++)  {
        GG[K+1][XDIR] = FF[K]*XL[K];
        GG[K+1][YDIR] = FF[K]*YL[K];
        GG[K+1][ZDIR] = FF[K]*ZL[K];
    }

    G110[XDIR] = GG[10][XDIR]+GG[1][XDIR]*C1;
    G110[YDIR] = GG[10][YDIR]+GG[1][YDIR]*C1;
    G110[ZDIR] = GG[10][ZDIR]+GG[1][ZDIR]*C1;
    G23[XDIR] = GG[2][XDIR]+GG[3][XDIR];
    G23[YDIR] = GG[2][YDIR]+GG[3][YDIR];
    G23[ZDIR] = GG[2][ZDIR]+GG[3][ZDIR];
    G45[XDIR]=GG[4][XDIR]+GG[5][XDIR];
    G45[YDIR]=GG[4][YDIR]+GG[5][YDIR];
    G45[ZDIR]=GG[4][ZDIR]+GG[5][ZDIR];
    TT1[XDIR] =GG[1][XDIR]*C2;
    TT1[YDIR] =GG[1][YDIR]*C2;
    TT1[ZDIR] =GG[1][ZDIR]*C2;
    TT[XDIR] =G23[XDIR]*C2+TT1[XDIR];
    TT[YDIR] =G23[YDIR]*C2+TT1[YDIR];
    TT[ZDIR] =G23[ZDIR]*C2+TT1[ZDIR];
    TT2[XDIR]=G45[XDIR]*C2+TT1[XDIR];
    TT2[YDIR]=G45[YDIR]*C2+TT1[YDIR];
    TT2[ZDIR]=G45[ZDIR]*C2+TT1[ZDIR];
    /* lock locations for the molecule to be updated */
    tx_p = PFORCES[ProcID][mol][XDIR];
    ty_p = PFORCES[ProcID][mol][YDIR];
    tz_p = PFORCES[ProcID][mol][ZDIR];
    tx_p[H1] +=
        GG[6][XDIR]+GG[7][XDIR]+GG[13][XDIR]+TT[XDIR]+GG[4][XDIR];
    tx_p[O] +=
        G110[XDIR] + GG[11][XDIR] +GG[12][XDIR]+C1*G23[XDIR];
    tx_p[H2] +=
        GG[8][XDIR]+GG[9][XDIR]+GG[14][XDIR]+TT[XDIR]+GG[5][XDIR];
    ty_p[H1] +=
        GG[6][YDIR]+GG[7][YDIR]+GG[13][YDIR]+TT[YDIR]+GG[4][YDIR];
    ty_p[O]  +=
        G110[YDIR]+GG[11][YDIR]+GG[12][YDIR]+C1*G23[YDIR];
    ty_p[H2] +=
        GG[8][YDIR]+GG[9][YDIR]+GG[14][YDIR]+TT[YDIR]+GG[5][YDIR];
    tz_p[H1] +=
        GG[6][ZDIR]+GG[7][ZDIR]+GG[13][ZDIR]+TT[ZDIR]+GG[4][ZDIR];
    tz_p[O]  +=
        G110[ZDIR]+GG[11][ZDIR]+GG[12][ZDIR]+C1*G23[ZDIR];
    tz_p[H2] +=
        GG[8][ZDIR]+GG[9][ZDIR]+GG[14][ZDIR]+TT[ZDIR]+GG[5][ZDIR];

    tx_p = PFORCES[ProcID][comp][XDIR];
    ty_p = PFORCES[ProcID][comp][YDIR];
    tz_p = PFORCES[ProcID][comp][ZDIR];
    tx_p[H1] +=
        -GG[6][XDIR]-GG[8][XDIR]-GG[11][XDIR]-TT2[XDIR]-GG[2][XDIR];
    tx_p[O] +=
        -G110[XDIR]-GG[13][XDIR]-GG[14][XDIR]-C1*G45[XDIR];
    tx_p[H2] +=
        -GG[7][XDIR]-GG[9][XDIR]-GG[12][XDIR]-TT2[XDIR]-GG[3][XDIR];
    ty_p[H1] +=
        -GG[6][YDIR]-GG[8][YDIR]-GG[11][YDIR]-TT2[YDIR]-GG[2][YDIR];
    ty_p[O] +=
        -G110[YDIR]-GG[13][YDIR]-GG[14][YDIR]-C1*G45[YDIR];
    ty_p[H2] +=
        -GG[7][YDIR]-GG[9][YDIR]-GG[12][YDIR]-TT2[YDIR]-GG[3][YDIR];
    tz_p[H1] +=
        -GG[6][ZDIR]-GG[8][ZDIR]-GG[11][ZDIR]-TT2[ZDIR]-GG[2][ZDIR];
    tz_p[O] +=
        -G110[ZDIR]-GG[13][ZDIR]-GG[14][ZDIR]-C1*G45[ZDIR];
    tz_p[H2] +=
        -GG[7][ZDIR]-GG[9][ZDIR]-GG[12][ZDIR]-TT2[ZDIR]-GG[3][ZDIR];
}           /* end of subroutine UPDATE_FORCES */
#line 228 "/home/jyy/Dev/splash2/codes/null_macros/c.m4.null.POSIX"

#line 1 "intraf.C"
/*************************************************************************/
/*                                                                       */
/*  Copyright (c) 1994 Stanford University                               */
/*                                                                       */
/*  All rights reserved.                                                 */
/*                                                                       */
/*  Permission is given to use, copy, and modify this software for any   */
/*  non-commercial purpose as long as this copyright notice is not       */
/*  removed.  All other uses, including redistribution in whole or in    */
/*  part, are forbidden without prior written permission.                */
/*                                                                       */
/*  This software is provided with absolutely no warranty and no         */
/*  support.                                                             */
/*                                                                       */
/*************************************************************************/


#line 17
#include <pthread.h>
#line 17
#include <sys/time.h>
#line 17
#include <unistd.h>
#line 17
#include <stdlib.h>
#line 17
extern pthread_t PThreadTable[];
#line 17

#include "math.h"
#include "frcnst.h"
#include "mdvar.h"
#include "water.h"
#include "wwpot.h"
#include "parameters.h"
#include "mddata.h"
#include "split.h"
#include "global.h"

void INTRAF(double *VIR, long ProcID)
{
    /*
      .....this routine calculates the intra-molecular force/mass acting on
      each atom.
      FC11, FC12, FC13, AND FC33 are the quardratic force constants
      FC111, FC112, ....... ETC. are the cubic      force constants
      FC1111, FC1112 ...... ETC. are the quartic    force constants
      */

    double SUM, R1, R2, VR1[4], VR2[4], COS, SIN;
    double DT, DTS, DR1, DR1S, DR2, DR2S, R1S, R2S, DR11[4], DR23[4];
    double DT1[4], DT3[4], F1, F2, F3, T1, T2;
    long mol, dir, atom;
    double LVIR;  /* process keeps a local copy of the sum,
                     to reduce synchronized updates*/
    double *temp_p;

    /* loop through the molecules */
    for (mol = StartMol[ProcID]; mol < StartMol[ProcID+1];  mol++) {
        SUM=0.0;
        R1=0.0;
        R2=0.0;
        /* loop through the three directions */
        for (dir = XDIR; dir <= ZDIR; dir++) {
            temp_p = VAR[mol].F[DISP][dir];
            VAR[mol].VM[dir] = C1 * temp_p[O]
                + C2 * (temp_p[H1] +
                        temp_p[H2] );
            VR1[dir] = temp_p[O] - temp_p[H1];
            R1 += VR1[dir] * VR1[dir] ;
            VR2[dir] = temp_p[O] - temp_p[H2];
            R2 += VR2[dir] * VR2[dir];
            SUM += VR1[dir] * VR2[dir];
        } /* for dir */

        R1=sqrt(R1);
        R2=sqrt(R2);

        /* calculate cos(THETA), sin(THETA), delta(R1),
           delta(R2), and delta(THETA) */
        COS=SUM/(R1*R2);
        SIN=sqrt(ONE-COS*COS);
        DT=(acos(COS)-ANGLE)*ROH;
        DTS=DT*DT;
        DR1=R1-ROH;
        DR1S=DR1*DR1;
        DR2=R2-ROH;
        DR2S=DR2*DR2;

        /* calculate derivatives of R1/X1, R2/X3, THETA/X1, and THETA/X3 */

        R1S=ROH/(R1*SIN);
        R2S=ROH/(R2*SIN);
        for (dir = XDIR; dir <= ZDIR; dir++) {
            DR11[dir]=VR1[dir]/R1;
            DR23[dir]=VR2[dir]/R2;
            DT1[dir]=(-DR23[dir]+DR11[dir]*COS)*R1S;
            DT3[dir]=(-DR11[dir]+DR23[dir]*COS)*R2S;
        } /* for dir */

        /* calculate forces */
        F1=FC11*DR1+FC12*DR2+FC13*DT;
        F2=FC33*DT +FC13*(DR1+DR2);
        F3=FC11*DR2+FC12*DR1+FC13*DT;
        F1=F1+(3.0*FC111*DR1S+FC112*(2.0*DR1+DR2)*DR2
               +2.0*FC113*DR1*DT+FC123*DR2*DT+FC133*DTS)*ROHI;
        F2=F2+(3.0*FC333*DTS+FC113*(DR1S+DR2S)
               +FC123*DR1*DR2+2.0*FC133*(DR1+DR2)*DT)*ROHI;
        F3=F3+(3.0*FC111*DR2S+FC112*(2.0*DR2+DR1)*DR1
               +2.0*FC113*DR2*DT+FC123*DR1*DT+FC133*DTS)*ROHI;
        F1=F1+(4.0*FC1111*DR1S*DR1+FC1112*(3.0*DR1S+DR2S)
               *DR2+2.0*FC1122*DR1*DR2S+3.0*FC1113*DR1S*DT
               +FC1123*(2.0*DR1+DR2)*DR2*DT+(2.0*FC1133*DR1
                                             +FC1233*DR2+FC1333*DT)*DTS)*ROHI2;
        F2=F2+(4.0*FC3333*DTS*DT+FC1113*(DR1S*DR1+DR2S*DR2)
               +FC1123*(DR1+DR2)*DR1*DR2+2.0*FC1133*(DR1S+DR2S)
               *DT+2.0*FC1233*DR1*DR2*DT+3.0*FC1333*(DR1+DR2)*DTS)
            *ROHI2;
        F3=F3+(4.0*FC1111*DR2S*DR2+FC1112*(3.0*DR2S+DR1S)
               *DR1+2.0*FC1122*DR1S*DR2+3.0*FC1113*DR2S*DT
               +FC1123*(2.0*DR2+DR1)*DR1*DT+(2.0*FC1133*DR2
                                             +FC1233*DR1+FC1333*DT)*DTS)*ROHI2;

        for (dir = XDIR; dir <= ZDIR; dir++) {
            temp_p = VAR[mol].F[FORCES][dir];
            T1=F1*DR11[dir]+F2*DT1[dir];
            temp_p[H1] = T1;
            T2=F3*DR23[dir]+F2*DT3[dir];
            temp_p[H2] = T2;
            temp_p[O] = -(T1+T2);
        } /* for dir */
    } /* for mol */

    /* calculate summation of the product of the displacement and computed
       force for every molecule, direction, and atom */

    LVIR=0.0;
    for (mol = StartMol[ProcID]; mol < StartMol[ProcID+1];  mol++)
        for ( dir = XDIR; dir <= ZDIR; dir++)
            for (atom = 0; atom < NATOM; atom++)
                LVIR += VAR[mol].F[DISP][dir][atom] *
                    VAR[mol].F[FORCES][dir][atom];

    {pthread_mutex_lock(&(gl->IntrafVirLock));};
    *VIR =  *VIR + LVIR;
    {pthread_mutex_unlock(&(gl->IntrafVirLock));};
} /* end of subroutine INTRAF */
#line 228 "/home/jyy/Dev/splash2/codes/null_macros/c.m4.null.POSIX"

#line 1 "kineti.C"
/*************************************************************************/
/*                                                                       */
/*  Copyright (c) 1994 Stanford University                               */
/*                                                                       */
/*  All rights reserved.                                                 */
/*                                                                       */
/*  Permission is given to use, copy, and modify this software for any   */
/*  non-commercial purpose as long as this copyright notice is not       */
/*  removed.  All other uses, including redistribution in whole or in    */
/*  part, are forbidden without prior written permission.                */
/*                                                                       */
/*  This software is provided with absolutely no warranty and no         */
/*  support.                                                             */
/*                                                                       */
/*************************************************************************/


#line 17
#include <pthread.h>
#line 17
#include <sys/time.h>
#line 17
#include <unistd.h>
#line 17
#include <stdlib.h>
#line 17
extern pthread_t PThreadTable[];
#line 17

#include "math.h"
#include "mdvar.h"
#include "parameters.h"
#include "mddata.h"
#include "split.h"
#include "global.h"

  /* this routine computes kinetic energy in each of the three spatial
     dimensions, and puts the computed values in the SUM array */
void KINETI(double *SUM, double HMAS, double OMAS, long ProcID)
{
    long dir, mol;
    double S;

    /* loop over the three directions */
    for (dir = XDIR; dir <= ZDIR; dir++) {
        S=0.0;
        /* loop over the molecules */
        for (mol = StartMol[ProcID]; mol < StartMol[ProcID+1]; mol++) {
            double *tempptr = VAR[mol].F[VEL][dir];
            S += ( tempptr[H1] * tempptr[H1] +
                  tempptr[H2] * tempptr[H2] ) * HMAS
                      + (tempptr[O] * tempptr[O]) * OMAS;
        }
        {pthread_mutex_lock(&(gl->KinetiSumLock));};
        SUM[dir]+=S;
        {pthread_mutex_unlock(&(gl->KinetiSumLock));};
    } /* for */
} /* end of subroutine KINETI */

#line 228 "/home/jyy/Dev/splash2/codes/null_macros/c.m4.null.POSIX"

#line 1 "mdmain.C"
/*************************************************************************/
/*                                                                       */
/*  Copyright (c) 1994 Stanford University                               */
/*                                                                       */
/*  All rights reserved.                                                 */
/*                                                                       */
/*  Permission is given to use, copy, and modify this software for any   */
/*  non-commercial purpose as long as this copyright notice is not       */
/*  removed.  All other uses, including redistribution in whole or in    */
/*  part, are forbidden without prior written permission.                */
/*                                                                       */
/*  This software is provided with absolutely no warranty and no         */
/*  support.                                                             */
/*                                                                       */
/*************************************************************************/


#line 17
#include <pthread.h>
#line 17
#include <sys/time.h>
#line 17
#include <unistd.h>
#line 17
#include <stdlib.h>
#line 17
extern pthread_t PThreadTable[];
#line 17

#include "stdio.h"
#include "parameters.h"
#include "mdvar.h"
#include "water.h"
#include "wwpot.h"
#include "cnst.h"
#include "mddata.h"
#include "fileio.h"
#include "split.h"
#include "global.h"

/************************************************************************/

/* routine that implements the time-steps. Called by main routine and calls others */
double MDMAIN(long NSTEP, long NPRINT, long NSAVE, long NORD1, long ProcID)
{
    double XTT;
    long i;
    double POTA,POTR,POTRF;
    double XVIR,AVGT,TEN;
    double TTMV = 0.0, TKIN = 0.0, TVIR = 0.0;

    /*.......ESTIMATE ACCELERATION FROM F/M */
    INTRAF(&gl->VIR,ProcID);

    {
#line 43
	unsigned long	Error, Cycle;
#line 43
	long		Cancel, Temp;
#line 43

#line 43
	Error = pthread_mutex_lock(&(gl->start).mutex);
#line 43
	if (Error != 0) {
#line 43
		printf("Error while trying to get lock in barrier.\n");
#line 43
		exit(-1);
#line 43
	}
#line 43

#line 43
	Cycle = (gl->start).cycle;
#line 43
	if (++(gl->start).counter != (NumProcs)) {
#line 43
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 43
		while (Cycle == (gl->start).cycle) {
#line 43
			Error = pthread_cond_wait(&(gl->start).cv, &(gl->start).mutex);
#line 43
			if (Error != 0) {
#line 43
				break;
#line 43
			}
#line 43
		}
#line 43
		pthread_setcancelstate(Cancel, &Temp);
#line 43
	} else {
#line 43
		(gl->start).cycle = !(gl->start).cycle;
#line 43
		(gl->start).counter = 0;
#line 43
		Error = pthread_cond_broadcast(&(gl->start).cv);
#line 43
	}
#line 43
	pthread_mutex_unlock(&(gl->start).mutex);
#line 43
};

    INTERF(ACC,&gl->VIR,ProcID);

    {
#line 47
	unsigned long	Error, Cycle;
#line 47
	long		Cancel, Temp;
#line 47

#line 47
	Error = pthread_mutex_lock(&(gl->start).mutex);
#line 47
	if (Error != 0) {
#line 47
		printf("Error while trying to get lock in barrier.\n");
#line 47
		exit(-1);
#line 47
	}
#line 47

#line 47
	Cycle = (gl->start).cycle;
#line 47
	if (++(gl->start).counter != (NumProcs)) {
#line 47
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 47
		while (Cycle == (gl->start).cycle) {
#line 47
			Error = pthread_cond_wait(&(gl->start).cv, &(gl->start).mutex);
#line 47
			if (Error != 0) {
#line 47
				break;
#line 47
			}
#line 47
		}
#line 47
		pthread_setcancelstate(Cancel, &Temp);
#line 47
	} else {
#line 47
		(gl->start).cycle = !(gl->start).cycle;
#line 47
		(gl->start).counter = 0;
#line 47
		Error = pthread_cond_broadcast(&(gl->start).cv);
#line 47
	}
#line 47
	pthread_mutex_unlock(&(gl->start).mutex);
#line 47
};

    /* MOLECULAR DYNAMICS LOOP OVER ALL TIME-STEPS */

    for (i=1;i <= NSTEP; i++) {
        TTMV=TTMV+1.00;

        /* reset simulator stats at beginning of second time-step */

        /* POSSIBLE ENHANCEMENT:  Here's where one start measurements to avoid
           cold-start effects.  Recommended to do this at the beginning of the
           second timestep; i.e. if (i == 2).
           */

        /* initialize various shared sums */
        if (ProcID == 0) {
            long dir;
            if (i >= 2) {
                {
#line 65
	struct timeval	FullTime;
#line 65

#line 65
	gettimeofday(&FullTime, NULL);
#line 65
	(gl->trackstart) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 65
};
            }
            gl->VIR = 0.0;
            gl->POTA = 0.0;
            gl->POTR = 0.0;
            gl->POTRF = 0.0;
            for (dir = XDIR; dir <= ZDIR; dir++)
                gl->SUM[dir] = 0.0;
        }

        if ((ProcID == 0) && (i >= 2)) {
            {
#line 76
	struct timeval	FullTime;
#line 76

#line 76
	gettimeofday(&FullTime, NULL);
#line 76
	(gl->intrastart) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 76
};
        }

        {
#line 79
	unsigned long	Error, Cycle;
#line 79
	long		Cancel, Temp;
#line 79

#line 79
	Error = pthread_mutex_lock(&(gl->start).mutex);
#line 79
	if (Error != 0) {
#line 79
		printf("Error while trying to get lock in barrier.\n");
#line 79
		exit(-1);
#line 79
	}
#line 79

#line 79
	Cycle = (gl->start).cycle;
#line 79
	if (++(gl->start).counter != (NumProcs)) {
#line 79
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 79
		while (Cycle == (gl->start).cycle) {
#line 79
			Error = pthread_cond_wait(&(gl->start).cv, &(gl->start).mutex);
#line 79
			if (Error != 0) {
#line 79
				break;
#line 79
			}
#line 79
		}
#line 79
		pthread_setcancelstate(Cancel, &Temp);
#line 79
	} else {
#line 79
		(gl->start).cycle = !(gl->start).cycle;
#line 79
		(gl->start).counter = 0;
#line 79
		Error = pthread_cond_broadcast(&(gl->start).cv);
#line 79
	}
#line 79
	pthread_mutex_unlock(&(gl->start).mutex);
#line 79
};
        PREDIC(TLC,NORD1,ProcID);
        INTRAF(&gl->VIR,ProcID);
        {
#line 82
	unsigned long	Error, Cycle;
#line 82
	long		Cancel, Temp;
#line 82

#line 82
	Error = pthread_mutex_lock(&(gl->start).mutex);
#line 82
	if (Error != 0) {
#line 82
		printf("Error while trying to get lock in barrier.\n");
#line 82
		exit(-1);
#line 82
	}
#line 82

#line 82
	Cycle = (gl->start).cycle;
#line 82
	if (++(gl->start).counter != (NumProcs)) {
#line 82
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 82
		while (Cycle == (gl->start).cycle) {
#line 82
			Error = pthread_cond_wait(&(gl->start).cv, &(gl->start).mutex);
#line 82
			if (Error != 0) {
#line 82
				break;
#line 82
			}
#line 82
		}
#line 82
		pthread_setcancelstate(Cancel, &Temp);
#line 82
	} else {
#line 82
		(gl->start).cycle = !(gl->start).cycle;
#line 82
		(gl->start).counter = 0;
#line 82
		Error = pthread_cond_broadcast(&(gl->start).cv);
#line 82
	}
#line 82
	pthread_mutex_unlock(&(gl->start).mutex);
#line 82
};

        if ((ProcID == 0) && (i >= 2)) {
            {
#line 85
	struct timeval	FullTime;
#line 85

#line 85
	gettimeofday(&FullTime, NULL);
#line 85
	(gl->intraend) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 85
};
            gl->intratime += gl->intraend - gl->intrastart;
        }


        if ((ProcID == 0) && (i >= 2)) {
            {
#line 91
	struct timeval	FullTime;
#line 91

#line 91
	gettimeofday(&FullTime, NULL);
#line 91
	(gl->interstart) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 91
};
        }

        INTERF(FORCES,&gl->VIR,ProcID);

        if ((ProcID == 0) && (i >= 2)) {
            {
#line 97
	struct timeval	FullTime;
#line 97

#line 97
	gettimeofday(&FullTime, NULL);
#line 97
	(gl->interend) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 97
};
            gl->intertime += gl->interend - gl->interstart;
        }

        if ((ProcID == 0) && (i >= 2)) {
            {
#line 102
	struct timeval	FullTime;
#line 102

#line 102
	gettimeofday(&FullTime, NULL);
#line 102
	(gl->intrastart) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 102
};
        }

        CORREC(PCC,NORD1,ProcID);

        BNDRY(ProcID);

        KINETI(gl->SUM,HMAS,OMAS,ProcID);

        {
#line 111
	unsigned long	Error, Cycle;
#line 111
	long		Cancel, Temp;
#line 111

#line 111
	Error = pthread_mutex_lock(&(gl->start).mutex);
#line 111
	if (Error != 0) {
#line 111
		printf("Error while trying to get lock in barrier.\n");
#line 111
		exit(-1);
#line 111
	}
#line 111

#line 111
	Cycle = (gl->start).cycle;
#line 111
	if (++(gl->start).counter != (NumProcs)) {
#line 111
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 111
		while (Cycle == (gl->start).cycle) {
#line 111
			Error = pthread_cond_wait(&(gl->start).cv, &(gl->start).mutex);
#line 111
			if (Error != 0) {
#line 111
				break;
#line 111
			}
#line 111
		}
#line 111
		pthread_setcancelstate(Cancel, &Temp);
#line 111
	} else {
#line 111
		(gl->start).cycle = !(gl->start).cycle;
#line 111
		(gl->start).counter = 0;
#line 111
		Error = pthread_cond_broadcast(&(gl->start).cv);
#line 111
	}
#line 111
	pthread_mutex_unlock(&(gl->start).mutex);
#line 111
};

        if ((ProcID == 0) && (i >= 2)) {
            {
#line 114
	struct timeval	FullTime;
#line 114

#line 114
	gettimeofday(&FullTime, NULL);
#line 114
	(gl->intraend) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 114
};
            gl->intratime += gl->intraend - gl->intrastart;
        }

        TKIN=TKIN+gl->SUM[0]+gl->SUM[1]+gl->SUM[2];
        TVIR=TVIR-gl->VIR;

        /*  check if potential energy is to be computed, and if
            printing and/or saving is to be done, this time step.
            Note that potential energy is computed once every NPRINT
            time-steps */

        if (((i % NPRINT) == 0) || ( (NSAVE > 0) && ((i % NSAVE) == 0))){

            if ((ProcID == 0) && (i >= 2)) {
                {
#line 129
	struct timeval	FullTime;
#line 129

#line 129
	gettimeofday(&FullTime, NULL);
#line 129
	(gl->interstart) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 129
};
            }

            /*  call potential energy computing routine */
            POTENG(&gl->POTA,&gl->POTR,&gl->POTRF,ProcID);
            {
#line 134
	unsigned long	Error, Cycle;
#line 134
	long		Cancel, Temp;
#line 134

#line 134
	Error = pthread_mutex_lock(&(gl->start).mutex);
#line 134
	if (Error != 0) {
#line 134
		printf("Error while trying to get lock in barrier.\n");
#line 134
		exit(-1);
#line 134
	}
#line 134

#line 134
	Cycle = (gl->start).cycle;
#line 134
	if (++(gl->start).counter != (NumProcs)) {
#line 134
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 134
		while (Cycle == (gl->start).cycle) {
#line 134
			Error = pthread_cond_wait(&(gl->start).cv, &(gl->start).mutex);
#line 134
			if (Error != 0) {
#line 134
				break;
#line 134
			}
#line 134
		}
#line 134
		pthread_setcancelstate(Cancel, &Temp);
#line 134
	} else {
#line 134
		(gl->start).cycle = !(gl->start).cycle;
#line 134
		(gl->start).counter = 0;
#line 134
		Error = pthread_cond_broadcast(&(gl->start).cv);
#line 134
	}
#line 134
	pthread_mutex_unlock(&(gl->start).mutex);
#line 134
};

            if ((ProcID == 0) && (i >= 2)) {
                {
#line 137
	struct timeval	FullTime;
#line 137

#line 137
	gettimeofday(&FullTime, NULL);
#line 137
	(gl->interend) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 137
};
                gl->intertime += gl->interend - gl->interstart;
            }

            POTA=gl->POTA*FPOT;
            POTR=gl->POTR*FPOT;
            POTRF=gl->POTRF*FPOT;

            /* compute some values to print */
            XVIR=TVIR*FPOT*0.50/TTMV;
            AVGT=TKIN*FKIN*TEMP*2.00/(3.00*TTMV);
            TEN=(gl->SUM[0]+gl->SUM[1]+gl->SUM[2])*FKIN;
            XTT=POTA+POTR+POTRF+TEN;

            if ((i % NPRINT) == 0 && ProcID == 0) {
                fprintf(six,"     %5ld %14.5lf %12.5lf %12.5lf  \
                %12.5lf\n %16.3lf %16.5lf %16.5lf\n",
                        i,TEN,POTA,POTR,POTRF,XTT,AVGT,XVIR);
            }
        }

        /* wait for everyone to finish time-step */
        {
#line 159
	unsigned long	Error, Cycle;
#line 159
	long		Cancel, Temp;
#line 159

#line 159
	Error = pthread_mutex_lock(&(gl->start).mutex);
#line 159
	if (Error != 0) {
#line 159
		printf("Error while trying to get lock in barrier.\n");
#line 159
		exit(-1);
#line 159
	}
#line 159

#line 159
	Cycle = (gl->start).cycle;
#line 159
	if (++(gl->start).counter != (NumProcs)) {
#line 159
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 159
		while (Cycle == (gl->start).cycle) {
#line 159
			Error = pthread_cond_wait(&(gl->start).cv, &(gl->start).mutex);
#line 159
			if (Error != 0) {
#line 159
				break;
#line 159
			}
#line 159
		}
#line 159
		pthread_setcancelstate(Cancel, &Temp);
#line 159
	} else {
#line 159
		(gl->start).cycle = !(gl->start).cycle;
#line 159
		(gl->start).counter = 0;
#line 159
		Error = pthread_cond_broadcast(&(gl->start).cv);
#line 159
	}
#line 159
	pthread_mutex_unlock(&(gl->start).mutex);
#line 159
};

        if ((ProcID == 0) && (i >= 2)) {
            {
#line 162
	struct timeval	FullTime;
#line 162

#line 162
	gettimeofday(&FullTime, NULL);
#line 162
	(gl->trackend) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 162
};
            gl->tracktime += gl->trackend - gl->trackstart;
        }
    } /* for i */

    return(XTT);

} /* end of subroutine MDMAIN */
#line 228 "/home/jyy/Dev/splash2/codes/null_macros/c.m4.null.POSIX"

#line 1 "poteng.C"
/*************************************************************************/
/*                                                                       */
/*  Copyright (c) 1994 Stanford University                               */
/*                                                                       */
/*  All rights reserved.                                                 */
/*                                                                       */
/*  Permission is given to use, copy, and modify this software for any   */
/*  non-commercial purpose as long as this copyright notice is not       */
/*  removed.  All other uses, including redistribution in whole or in    */
/*  part, are forbidden without prior written permission.                */
/*                                                                       */
/*  This software is provided with absolutely no warranty and no         */
/*  support.                                                             */
/*                                                                       */
/*************************************************************************/


#line 17
#include <pthread.h>
#line 17
#include <sys/time.h>
#line 17
#include <unistd.h>
#line 17
#include <stdlib.h>
#line 17
extern pthread_t PThreadTable[];
#line 17

#include "mdvar.h"
#include "frcnst.h"
#include "water.h"
#include "wwpot.h"
#include "math.h"
#include "parameters.h"
#include "mddata.h"
#include "split.h"
#include "global.h"

void POTENG(double *POTA, double *POTR, double *PTRF, long ProcID)
{

    /*
      this routine calculates the potential energy of the system.
      FC11 ,FC12, FC13, and FC33 are the quardratic force constants
      */

    long mol,comp;
    long half_mol;
    long    KC, K;
    double R1, R2, RX, COS, DT, DR1, DR2, DR1S, DR2S, DRP;
    double XL[15], YL[15], ZL[15], RS[15], RL[15];
    double DTS;
    double LPOTA, LPOTR, LPTRF;
    double *tx_p, *ty_p, *tz_p;

    /*  compute intra-molecular potential energy */
    LPOTA=0.0;
    for (mol = StartMol[ProcID]; mol < StartMol[ProcID+1]; mol++) {
        double dx1, dy1, dz1, dx2, dy2, dz2;
        tx_p = VAR[mol].F[DISP][XDIR];
        ty_p = VAR[mol].F[DISP][YDIR];
        tz_p = VAR[mol].F[DISP][ZDIR];
        VAR[mol].VM[XDIR] = C1 * tx_p[ O] +
            C2 * (tx_p[H1] +
                  tx_p[H2] );
        VAR[mol].VM[YDIR] = C1*ty_p[O] +
            C2*(ty_p[H1] +
                ty_p[H2] );
        VAR[mol].VM[ZDIR] = C1*tz_p[O] +
            C2*(tz_p[H1] +
                tz_p[H2] );
        dx1 = tx_p[O]-tx_p[H1];
        dy1 = ty_p[O]-ty_p[H1];
        dz1 = tz_p[O]-tz_p[H1];

        dx2 = tx_p[O]-tx_p[H2];
        dy2 = ty_p[O]-ty_p[H2];
        dz2 = tz_p[O]-tz_p[H2];
        R1 = dx1*dx1 + dy1*dy1 + dz1*dz1;
        R2 = dx2*dx2 + dy2*dy2 + dz2*dz2;
        RX = dx1*dx2 + dy1*dy2 + dz1*dz2;

        R1=sqrt(R1);
        R2=sqrt(R2);
        COS=RX/(R1*R2);
        DT=(acos(COS)-ANGLE)*ROH;
        DR1=R1-ROH;
        DR2=R2-ROH;
        DR1S=DR1*DR1;
        DR2S=DR2*DR2;
        DRP=DR1+DR2;
        DTS=DT*DT;
        LPOTA += (FC11*(DR1S+DR2S)+FC33*DTS)*0.5
            +FC12*DR1*DR2+FC13*DRP*DT
                +(FC111*(DR1S*DR1+DR2S*DR2)+FC333*DTS*DT+FC112*DRP*DR1*DR2+
                  FC113*(DR1S+DR2S)*DT+FC123*DR1*DR2*DT+FC133*DRP*DTS)*ROHI;
        LPOTA += (FC1111*(DR1S*DR1S+DR2S*DR2S)+FC3333*DTS*DTS+
                  FC1112*(DR1S+DR2S)*DR1*DR2+FC1122*DR1S*DR2S+
                  FC1113*(DR1S*DR1+DR2S*DR2)*DT+FC1123*DRP*DR1*DR2*DT+
                  FC1133*(DR1S+DR2S)*DTS+FC1233*DR1*DR2*DTS+
                  FC1333*DRP*DTS*DT)*ROHI2;
    } /* for mol */

    {
#line 93
	unsigned long	Error, Cycle;
#line 93
	long		Cancel, Temp;
#line 93

#line 93
	Error = pthread_mutex_lock(&(gl->PotengBar).mutex);
#line 93
	if (Error != 0) {
#line 93
		printf("Error while trying to get lock in barrier.\n");
#line 93
		exit(-1);
#line 93
	}
#line 93

#line 93
	Cycle = (gl->PotengBar).cycle;
#line 93
	if (++(gl->PotengBar).counter != (NumProcs)) {
#line 93
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 93
		while (Cycle == (gl->PotengBar).cycle) {
#line 93
			Error = pthread_cond_wait(&(gl->PotengBar).cv, &(gl->PotengBar).mutex);
#line 93
			if (Error != 0) {
#line 93
				break;
#line 93
			}
#line 93
		}
#line 93
		pthread_setcancelstate(Cancel, &Temp);
#line 93
	} else {
#line 93
		(gl->PotengBar).cycle = !(gl->PotengBar).cycle;
#line 93
		(gl->PotengBar).counter = 0;
#line 93
		Error = pthread_cond_broadcast(&(gl->PotengBar).cv);
#line 93
	}
#line 93
	pthread_mutex_unlock(&(gl->PotengBar).mutex);
#line 93
};

    /*  compute inter-molecular potential energy */
    LPOTR=0.0;
    LPTRF=0.0;
    half_mol = NMOL/2;
    for (mol = StartMol[ProcID]; mol < StartMol[ProcID+1]; mol++) {
        long comp_last = mol + half_mol;
        long icomp;
        if (NMOL%2 == 0) {
            if ((half_mol <= mol) && (mol%2 == 0)) {
                comp_last--;
            }
            if ((mol < half_mol) && (comp_last%2 == 1)) {
                comp_last--;
            }
        }
        for (icomp = mol+1; icomp <= comp_last; icomp++) {
            comp = icomp;
            if (comp > NMOL1) comp = comp%NMOL;
            CSHIFT(VAR[mol].F[DISP][XDIR],VAR[comp].F[DISP][XDIR],
                   VAR[mol].VM[XDIR], VAR[comp].VM[XDIR],XL,BOXH,BOXL);
            CSHIFT(VAR[mol].F[DISP][YDIR],VAR[comp].F[DISP][YDIR],
                   VAR[mol].VM[YDIR], VAR[comp].VM[YDIR],YL,BOXH,BOXL);
            CSHIFT(VAR[mol].F[DISP][ZDIR],VAR[comp].F[DISP][ZDIR],
                   VAR[mol].VM[ZDIR], VAR[comp].VM[ZDIR],ZL,BOXH,BOXL);
            KC=0;
            for (K = 0; K < 9; K++) {
                RS[K]=XL[K]*XL[K]+YL[K]*YL[K]+ZL[K]*ZL[K];
                if (RS[K] > CUT2)
                    KC=KC+1;
            } /* for k */
            if (KC != 9) {
                for (K = 0; K < 9; K++) {
                    if (RS[K] <= CUT2) {
                        RL[K]=sqrt(RS[K]);
                    }
                    else {
                        RL[K]=CUTOFF;
                        RS[K]=CUT2;
                    } /* else */
                } /* for K */
                LPOTR= LPOTR-QQ2/RL[1]-QQ2/RL[2]-QQ2/RL[3]-QQ2/RL[4]
                    +QQ /RL[5]+QQ /RL[6]+QQ /RL[7]+QQ /RL[8]
                        +QQ4/RL[0];
                LPTRF= LPTRF-REF2*RS[0]-REF1*((RS[5]+RS[6]+RS[7]+RS[8])*0.5
                                              -RS[1]-RS[2]-RS[3]-RS[4]);

                if (KC <= 0) {
                    for (K = 9; K <  14; K++) {
                        RL[K]=sqrt(XL[K]*XL[K]+YL[K]*YL[K]+ZL[K]*ZL[K]);
                    }
                    LPOTR= LPOTR+A1* exp(-B1*RL[9])
                        +A2*(exp(-B2*RL[ 5])+exp(-B2*RL[ 6])
                             +exp(-B2*RL[ 7])+exp(-B2*RL[ 8]))
                            +A3*(exp(-B3*RL[10])+exp(-B3*RL[11])
                                 +exp(-B3*RL[12])+exp(-B3*RL[13]))
                                -A4*(exp(-B4*RL[10])+exp(-B4*RL[11])
                                     +exp(-B4*RL[12])+exp(-B4*RL[13]));
                } /* if KC <= 0 */
            } /* if KC != 9 */
        } /* for comp */
    } /* for mol */

    /* update shared sums from computed  private sums */
    {pthread_mutex_lock(&(gl->PotengSumLock));};
    *POTA = *POTA + LPOTA;
    *POTR = *POTR + LPOTR;
    *PTRF = *PTRF + LPTRF;
    {pthread_mutex_unlock(&(gl->PotengSumLock));};
} /* end of subroutine POTENG */
#line 228 "/home/jyy/Dev/splash2/codes/null_macros/c.m4.null.POSIX"

#line 1 "predcor.C"
/*************************************************************************/
/*                                                                       */
/*  Copyright (c) 1994 Stanford University                               */
/*                                                                       */
/*  All rights reserved.                                                 */
/*                                                                       */
/*  Permission is given to use, copy, and modify this software for any   */
/*  non-commercial purpose as long as this copyright notice is not       */
/*  removed.  All other uses, including redistribution in whole or in    */
/*  part, are forbidden without prior written permission.                */
/*                                                                       */
/*  This software is provided with absolutely no warranty and no         */
/*  support.                                                             */
/*                                                                       */
/*************************************************************************/


#line 17
#include <pthread.h>
#line 17
#include <sys/time.h>
#line 17
#include <unistd.h>
#line 17
#include <stdlib.h>
#line 17
extern pthread_t PThreadTable[];
#line 17

#include "mdvar.h"
#include "parameters.h"
#include "mddata.h"
#include "split.h"
#include "global.h"

/* predicts new values for displacement and its five derivatives
 *
 * NOR1 : NOR1 = NORDER + 1 = 7 (for a sixth-order method)
 */
void PREDIC(double *C, long NOR1, long ProcID)
{
    /*   this routine calculates predicted F(X), F'(X), F''(X), ... */

    long JIZ;
    long  JI;
    long  L;
    double S;
    long func, mol, dir, atom;

    JIZ=2;
    /* .....loop over F(X), F'(X), F''(X), ..... */
    for (func = 0; func < NORDER; func++) {
        for (mol = StartMol[ProcID]; mol < StartMol[ProcID+1]; mol++)
            for ( dir = 0; dir < NDIR; dir++)
                for ( atom = 0; atom < NATOM; atom++ ) {
                    JI = JIZ;
                    /* sum over Taylor Series */
                    S = 0.0;
                    for ( L = func; L < NORDER; L++) {
                        S += C[JI] * VAR[mol].F[L+1][dir][atom];
                        JI++;
                    } /* for */
                    VAR[mol].F[func][dir][atom] += S;
                } /* for atom */
        JIZ += NOR1;
    } /* for func */
} /* end of subroutine PREDIC */

/* corrects the predicted values, based on forces etc. computed in the interim
 *
 * PCC : the predictor-corrector constants
 * NOR1: NORDER + 1 = 7 for a sixth-order method)
 */

void CORREC(double *PCC, long NOR1, long ProcID)
{
    /*
      .....this routine calculates corrected F(X), F'(X), F"(X), ....
      from corrected F(X) = predicted F(X) + PCC(1)*(FR-SD)
      where SD is predicted accl. F"(X) and FR is computed
      accl. (force/mass) at predicted position
      */

    double Y;
    long mol, dir, atom, func;

    for (mol = StartMol[ProcID]; mol < StartMol[ProcID+1]; mol++) {
        for (dir = 0; dir < NDIR; dir++) {
            for (atom = 0; atom < NATOM; atom++) {
                Y = VAR[mol].F[FORCES][dir][atom] - VAR[mol].F[ACC][dir][atom];
                for ( func = 0; func < NOR1; func++)
                    VAR[mol].F[func][dir][atom] += PCC[func] * Y;
            } /* for atom */
        } /* for dir */
    } /* for mol */

} /* end of subroutine CORREC */
#line 228 "/home/jyy/Dev/splash2/codes/null_macros/c.m4.null.POSIX"

#line 1 "syscons.C"
/*************************************************************************/
/*                                                                       */
/*  Copyright (c) 1994 Stanford University                               */
/*                                                                       */
/*  All rights reserved.                                                 */
/*                                                                       */
/*  Permission is given to use, copy, and modify this software for any   */
/*  non-commercial purpose as long as this copyright notice is not       */
/*  removed.  All other uses, including redistribution in whole or in    */
/*  part, are forbidden without prior written permission.                */
/*                                                                       */
/*  This software is provided with absolutely no warranty and no         */
/*  support.                                                             */
/*                                                                       */
/*************************************************************************/


#line 17
#include <pthread.h>
#line 17
#include <sys/time.h>
#line 17
#include <unistd.h>
#line 17
#include <stdlib.h>
#line 17
extern pthread_t PThreadTable[];
#line 17

#include "stdio.h"
#include <math.h>

#include "parameters.h"
#include "mdvar.h"
#include "water.h"
#include "wwpot.h"
#include "cnst.h"
#include "mddata.h"
#include "global.h"

void SYSCNS()                    /* sets up some system constants */
{
    TSTEP=TSTEP/UNITT;        /* time between steps */
    NATMO=NATOMS*NMOL;        /* total number of atoms in system */
    NATMO3=NATMO*3; /* number of atoms * number of spatial dimensions */
    FPOT= UNITM * pow((UNITL/UNITT),2.0) / (BOLTZ*TEMP*NATMO);
    FKIN=FPOT*0.50/(TSTEP*TSTEP);
    BOXL= pow( (NMOL*WTMOL*UNITM/RHO),(1.00/3.00));  /* computed
                                                        length of the cubical "box".  Note that box size is
                                                        computed as being large enough to handle the input
                                                        number of water molecules */

    BOXL=BOXL/UNITL;    /* normalized length of computational box */

    BOXH=BOXL*0.50; /* half the box length, used in
                       computing cutoff radius */

    if (CUTOFF == 0.0) {
        CUTOFF=max(BOXH,CUTOFF);    /* cutoff radius is max of BOXH
                                       and default (= 0); i.e. CUTOFF
                                       radius is set to half the normalized
                                       box length */
    }
    if (CUTOFF > 11.0) CUTOFF = 11.0; /* cutoff never greater than 11
					  Angstrom*/

    REF1= -QQ/(CUTOFF*CUTOFF*CUTOFF);
    REF2=2.00*REF1;
    REF4=2.00*REF2;
    CUT2=CUTOFF*CUTOFF;       /* square of cutoff radius,  used
                                 to actually decide whether an
                                 interaction should be computed in
                                 INTERF and POTENG */
    FHM=(TSTEP*TSTEP*0.50)/HMAS;
    FOM=(TSTEP*TSTEP*0.50)/OMAS;
    NMOL1=NMOL-1;
}       /* end of subroutine SYSCNS */
#line 228 "/home/jyy/Dev/splash2/codes/null_macros/c.m4.null.POSIX"

#line 1 "water.C"
/*************************************************************************/
/*                                                                       */
/*  Copyright (c) 1994 Stanford University                               */
/*                                                                       */
/*  All rights reserved.                                                 */
/*                                                                       */
/*  Permission is given to use, copy, and modify this software for any   */
/*  non-commercial purpose as long as this copyright notice is not       */
/*  removed.  All other uses, including redistribution in whole or in    */
/*  part, are forbidden without prior written permission.                */
/*                                                                       */
/*  This software is provided with absolutely no warranty and no         */
/*  support.                                                             */
/*                                                                       */
/*************************************************************************/

/*  Usage:   water < infile,
    where infile has 10 fields which can be described in order as
    follows:

    TSTEP:   the physical time interval (in sec) between timesteps.
    Good default is 1e-15.
    NMOL:    the number of molecules to be simulated.
    NSTEP:   the number of timesteps to be simulated.
    NORDER:  the order of the predictor-corrector method to be used.
    set this to 6.
    NSAVE:   the frequency with which to save data in data collection.
    Set to 0 always.
    NRST:    the frequency with which to write RST file: set to 0 always (not used).
    NPRINT:  the frequency with which to compute potential energy.
    i.e. the routine POTENG is called every NPRINT timesteps.
    It also computes intermolecular as well as intramolecular
    interactions, and hence is very expensive.
    NFMC:    Not used (historical artifact).  Set to anything, say 0.
    NumProcs: the number of processors to be used.
    CUTOFF:  the cutoff radius to be used (in Angstrom,
    floating-point).  In a real simulation, this
    will be set to 0 here in which case the program will
    compute it itself (and set it to about 11 Angstrom.
    It can be set by the user if they want
    to use an artificially small cutoff radius, for example
    to control the number of boxes created for small problems
    (and not have fewer boxes than processors).
    */


#line 46
#include <pthread.h>
#line 46
#include <sys/time.h>
#line 46
#include <unistd.h>
#line 46
#include <stdlib.h>
#line 46
#define MAX_THREADS 32
#line 46
pthread_t PThreadTable[MAX_THREADS];
#line 46

#include <stdio.h>
#include <string.h>
#include "split.h"

/*  include files for declarations  */
#include "parameters.h"
#include "mdvar.h"
#include "water.h"
#include "wwpot.h"
#include "cnst.h"
#include "mddata.h"
#include "fileio.h"
#include "frcnst.h"
#include "global.h"

long NMOL,NORDER,NATMO,NATMO3,NMOL1;
long NATOMS;
long I2;

double TLC[100], FPOT, FKIN;
double TEMP,RHO,TSTEP,BOXL,BOXH,CUTOFF,CUT2;
double R3[128],R1;
double UNITT,UNITL,UNITM,BOLTZ,AVGNO,PCC[11];
double OMAS,HMAS,WTMOL,ROH,ANGLE,FHM,FOM,ROHI,ROHI2;
double QQ,A1,B1,A2,B2,A3,B3,A4,B4,AB1,AB2,AB3,AB4,C1,C2,QQ2,QQ4,REF1,REF2,REF4;
double FC11,FC12,FC13,FC33,FC111,FC333,FC112,FC113,FC123,FC133,FC1111,FC3333,FC1112,FC1122,FC1113,FC1123,FC1133,FC1233,FC1333;

FILE *six;

molecule_type *VAR;

struct GlobalMemory *gl;        /* pointer to the Global Memory
                                   structure, which contains the lock,
                                   barrier, and some scalar variables */


long NSTEP, NSAVE, NRST, NPRINT,NFMC;
long NORD1;
long II;                         /*  variables explained in common.h */
long i;
long NDATA;
long NFRST=11;
long NFSV=10;
long LKT=0;

long StartMol[MAXPROCS+1];       /* number of the first molecule
                                   to be handled by this process; used
                                   for static scheduling     */
long MolsPerProc;                /* number of mols per processor */
long NumProcs;                   /* number of processors being used;
                                   run-time input           */
double XTT;

int main(int argc, char **argv)
{
    /* default values for the control parameters of the driver */
    /* are in parameters.h */

    if ((argc == 2) &&((strncmp(argv[1],"-h",strlen("-h")) == 0) || (strncmp(argv[1],"-H",strlen("-H")) == 0))) {
        printf("Usage:  WATER-NSQUARED < infile, where the contents of infile can be\nobtained from the comments at the top of water.C and the first scanf \nin main() in water.C\n\n");
        exit(0);
    }

    /*  POSSIBLE ENHANCEMENT:  Here's where one might bind the main process
        (process 0) to a processor if one wanted to. Others can be bound in
        the WorkStart routine.
        */

    six = stdout;   /* output file */

    TEMP  =298.0;
    RHO   =0.9980;
    CUTOFF=0.0;

    /* read input */

    freopen("infile", "r", stdin);
    if (scanf("%lf%ld%ld%ld%ld%ld%ld%ld%ld%lf",&TSTEP, &NMOL, &NSTEP, &NORDER,
              &NSAVE, &NRST, &NPRINT, &NFMC,&NumProcs, &CUTOFF) != 10)
        fprintf(stderr,"ERROR: Usage: water < infile, which must have 10 fields, see SPLASH documentation or comment at top of water.C\n");

    if (NMOL > MAXLCKS) {
        fprintf(stderr, "Just so you know ... Lock array in global.H has size %ld < %ld (NMOL)\n code will still run correctly but there may be lock contention\n\n", MAXLCKS, NMOL);
    }
    if (argc > 1) NumProcs = atoi(argv[1]);

    printf("Using %ld procs on %ld steps of %ld mols\n", NumProcs, NSTEP, NMOL);
    printf("Other parameters:\n\tTSTEP = %8.2e\n\tNORDER = %ld\n\tNSAVE = %ld\n",TSTEP,NORDER,NSAVE);
    printf("\tNRST = %ld\n\tNPRINT = %ld\n\tNFMC = %ld\n\tCUTOFF = %lf\n\n",NRST,NPRINT,NFMC,CUTOFF);


    /* SET UP SCALING FACTORS AND CONSTANTS */

    NORD1=NORDER+1;

    CNSTNT(NORD1,TLC);  /* sub. call to set up constants */


    { /* Do memory initializations */
        long pid;
        long mol_size = sizeof(molecule_type) * NMOL;
        long gmem_size = sizeof(struct GlobalMemory);

        /*  POSSIBLE ENHANCEMENT:  One might bind the first process to
            a processor here, even before the other (child) processes are
            bound later in mdmain().
            */

        {;};  /* macro call to initialize
                                      shared memory etc. */

        /* allocate space for main (VAR) data structure as well as
           synchronization variables */

        /*  POSSIBLE ENHANCEMENT: One might want to allocate a process's
            portion of the VAR array and what it points to in its local
            memory */

        VAR = (molecule_type *) valloc(mol_size);;
        gl = (struct GlobalMemory *) valloc(gmem_size);;

        /*  POSSIBLE ENHANCEMENT: One might want to allocate  process i's
            PFORCES[i] array in its local memory */

        PFORCES = (double ****) valloc(NumProcs * sizeof (double ***));;
        { long i,j,k;

          for (i = 0; i < NumProcs; i++) {
              PFORCES[i] = (double ***) valloc(NMOL * sizeof (double **));;
              for (j = 0; j < NMOL; j++) {
                  PFORCES[i][j] = (double **) valloc(NDIR * sizeof (double *));;
                  for (k = 0; k < NDIR; k++) {
                      PFORCES[i][j][k] = (double *) valloc(NATOM * sizeof (double));;
                  }
              }
          }
      }
        /* macro calls to initialize synch varibles  */

        {
#line 184
	unsigned long	Error;
#line 184

#line 184
	Error = pthread_mutex_init(&(gl->start).mutex, NULL);
#line 184
	if (Error != 0) {
#line 184
		printf("Error while initializing barrier.\n");
#line 184
		exit(-1);
#line 184
	}
#line 184

#line 184
	Error = pthread_cond_init(&(gl->start).cv, NULL);
#line 184
	if (Error != 0) {
#line 184
		printf("Error while initializing barrier.\n");
#line 184
		pthread_mutex_destroy(&(gl->start).mutex);
#line 184
		exit(-1);
#line 184
	}
#line 184

#line 184
	(gl->start).counter = 0;
#line 184
	(gl->start).cycle = 0;
#line 184
};
	{
#line 185
	unsigned long	Error;
#line 185

#line 185
	Error = pthread_mutex_init(&(gl->InterfBar).mutex, NULL);
#line 185
	if (Error != 0) {
#line 185
		printf("Error while initializing barrier.\n");
#line 185
		exit(-1);
#line 185
	}
#line 185

#line 185
	Error = pthread_cond_init(&(gl->InterfBar).cv, NULL);
#line 185
	if (Error != 0) {
#line 185
		printf("Error while initializing barrier.\n");
#line 185
		pthread_mutex_destroy(&(gl->InterfBar).mutex);
#line 185
		exit(-1);
#line 185
	}
#line 185

#line 185
	(gl->InterfBar).counter = 0;
#line 185
	(gl->InterfBar).cycle = 0;
#line 185
};
	{
#line 186
	unsigned long	Error;
#line 186

#line 186
	Error = pthread_mutex_init(&(gl->PotengBar).mutex, NULL);
#line 186
	if (Error != 0) {
#line 186
		printf("Error while initializing barrier.\n");
#line 186
		exit(-1);
#line 186
	}
#line 186

#line 186
	Error = pthread_cond_init(&(gl->PotengBar).cv, NULL);
#line 186
	if (Error != 0) {
#line 186
		printf("Error while initializing barrier.\n");
#line 186
		pthread_mutex_destroy(&(gl->PotengBar).mutex);
#line 186
		exit(-1);
#line 186
	}
#line 186

#line 186
	(gl->PotengBar).counter = 0;
#line 186
	(gl->PotengBar).cycle = 0;
#line 186
};
        {pthread_mutex_init(&(gl->IOLock), NULL);};
        {pthread_mutex_init(&(gl->IndexLock), NULL);};
        {pthread_mutex_init(&(gl->IntrafVirLock), NULL);};
        {pthread_mutex_init(&(gl->InterfVirLock), NULL);};
        {pthread_mutex_init(&(gl->FXLock), NULL);};
        {pthread_mutex_init(&(gl->FYLock), NULL);};
        {pthread_mutex_init(&(gl->FZLock), NULL);};
        if (NMOL < MAXLCKS) {
            {
#line 195
	unsigned long	i, Error;
#line 195

#line 195
	for (i = 0; i < NMOL; i++) {
#line 195
		Error = pthread_mutex_init(&gl->MolLock[i], NULL);
#line 195
		if (Error != 0) {
#line 195
			printf("Error while initializing array of locks.\n");
#line 195
			exit(-1);
#line 195
		}
#line 195
	}
#line 195
};
        }
        else {
            {
#line 198
	unsigned long	i, Error;
#line 198

#line 198
	for (i = 0; i < MAXLCKS; i++) {
#line 198
		Error = pthread_mutex_init(&gl->MolLock[i], NULL);
#line 198
		if (Error != 0) {
#line 198
			printf("Error while initializing array of locks.\n");
#line 198
			exit(-1);
#line 198
		}
#line 198
	}
#line 198
};
        }
        {pthread_mutex_init(&(gl->KinetiSumLock), NULL);};
        {pthread_mutex_init(&(gl->PotengSumLock), NULL);};

        /* set up control for static scheduling */

        MolsPerProc = NMOL/NumProcs;
        StartMol[0] = 0;
        for (pid = 1; pid < NumProcs; pid += 1) {
            StartMol[pid] = StartMol[pid-1] + MolsPerProc;
        }
        StartMol[NumProcs] = NMOL;
    }

    SYSCNS();    /* sub. call to initialize system constants  */

    fprintf(six,"\nTEMPERATURE                = %8.2f K\n",TEMP);
    fprintf(six,"DENSITY                    = %8.5f G/C.C.\n",RHO);
    fprintf(six,"NUMBER OF MOLECULES        = %8ld\n",NMOL);
    fprintf(six,"NUMBER OF PROCESSORS       = %8ld\n",NumProcs);
    fprintf(six,"TIME STEP                  = %8.2e SEC\n",TSTEP);
    fprintf(six,"ORDER USED TO SOLVE F=MA   = %8ld \n",NORDER);
    fprintf(six,"NO. OF TIME STEPS          = %8ld \n",NSTEP);
    fprintf(six,"FREQUENCY OF DATA SAVING   = %8ld \n",NSAVE);
    fprintf(six,"FREQUENCY TO WRITE RST FILE= %8ld \n",NRST);
    fprintf(six,"SPHERICAL CUTOFF RADIUS    = %8.4f ANGSTROM\n",CUTOFF);
    fflush(six);


    /* initialization routine; also reads displacements and
       sets up random velocities*/
    INITIA();

    /*.....start molecular dynamic loop */

    gl->tracktime = 0;
    gl->intratime = 0;
    gl->intertime = 0;

    /* initialize Index to 1 so that the first created child gets
       id 1, not 0 */

    gl->Index = 1;

    if (NSAVE > 0)  /* not true for input decks provided */
	fprintf(six,"COLLECTING X AND V DATA AT EVERY %4ld TIME STEPS \n",NSAVE);

    /* spawn helper processes, each getting its unique process id */
    {
#line 247
	struct timeval	FullTime;
#line 247

#line 247
	gettimeofday(&FullTime, NULL);
#line 247
	(gl->computestart) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 247
};
    {
#line 248
	long	i, Error;
#line 248

#line 248
	for (i = 0; i < (NumProcs) - 1; i++) {
#line 248
		Error = pthread_create(&PThreadTable[i], NULL, (void * (*)(void *))(WorkStart), NULL);
#line 248
		if (Error != 0) {
#line 248
			printf("Error in pthread_create().\n");
#line 248
			exit(-1);
#line 248
		}
#line 248
	}
#line 248

#line 248
	WorkStart();
#line 248
};

    /* macro to make main process wait for all others to finish */
    {
#line 251
	unsigned long	i, Error;
#line 251
	for (i = 0; i < (NumProcs) - 1; i++) {
#line 251
		Error = pthread_join(PThreadTable[i], NULL);
#line 251
		if (Error != 0) {
#line 251
			printf("Error in pthread_join().\n");
#line 251
			exit(-1);
#line 251
		}
#line 251
	}
#line 251
};
    {
#line 252
	struct timeval	FullTime;
#line 252

#line 252
	gettimeofday(&FullTime, NULL);
#line 252
	(gl->computeend) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 252
};

    printf("COMPUTESTART (after initialization) = %lu\n",gl->computestart);
    printf("COMPUTEEND = %lu\n",gl->computeend);
    printf("COMPUTETIME (after initialization) = %lu\n",gl->computeend-gl->computestart);
    printf("Measured Time (2nd timestep onward) = %lu\n",gl->tracktime);
    printf("Intramolecular time only (2nd timestep onward) = %lu\n",gl->intratime);
    printf("Intermolecular time only (2nd timestep onward) = %lu\n",gl->intertime);
    printf("Other time (2nd timestep onward) = %lu\n",gl->tracktime - gl->intratime - gl->intertime);

    printf("\nExited Happily with XTT = %g (note: XTT value is garbage if NPRINT > NSTEP)\n", XTT);

    {exit(0);};
} /* main.c */

void WorkStart() /* routine that each created process starts at;
                    it simply calls the timestep routine */
{
    long ProcID;
    double LocalXTT;

    {pthread_mutex_lock(&(gl->IndexLock));};
    ProcID = gl->Index++;
    {pthread_mutex_unlock(&(gl->IndexLock));};

    {;};
    {;};
    {;};

    ProcID = ProcID % NumProcs;

    /*  POSSIBLE ENHANCEMENT:  Here's where one might bind processes
        to processors if one wanted to.
        */

    LocalXTT = MDMAIN(NSTEP,NPRINT,NSAVE,NORD1,ProcID);
    if (ProcID == 0) {
	    XTT = LocalXTT;
    }
}
