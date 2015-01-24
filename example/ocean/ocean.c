#line 228 "/home/jyy/Dev/splash2/codes/null_macros/c.m4.null.POSIX"

#line 1 "jacobcalc2.C"
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

/* Does the arakawa jacobian calculation (of the x and y matrices,
   putting the results in the z matrix) for a subblock. */


#line 20
#include <pthread.h>
#line 20
#include <sys/time.h>
#line 20
#include <unistd.h>
#line 20
#include <stdlib.h>
#line 20
extern pthread_t PThreadTable[];
#line 20


#include <stdio.h>
#include <math.h>
#include <time.h>
#include "decs.h"

void jacobcalc2(double ****x, double ****y, double ****z, long psiindex, long pid, long firstrow, long lastrow, long firstcol, long lastcol)
{
   double f1;
   double f2;
   double f3;
   double f4;
   double f5;
   double f6;
   double f7;
   double f8;
   long iindex;
   long indexp1;
   long indexm1;
   long im1;
   long ip1;
   long i;
   long j;
   long jj;
   double **t2a;
   double **t2b;
   double **t2c;
   double *t1a;
   double *t1b;
   double *t1c;
   double *t1d;
   double *t1e;
   double *t1f;
   double *t1g;

   t2a = z[pid][psiindex];
   if ((gp[pid].neighbors[UP] == -1) && (gp[pid].neighbors[LEFT] == -1)) {
     t2a[0][0]=0.0;
   }
   if ((gp[pid].neighbors[DOWN] == -1) && (gp[pid].neighbors[LEFT] == -1)) {
     t2a[im-1][0]=0.0;
   }
   if ((gp[pid].neighbors[UP] == -1) && (gp[pid].neighbors[RIGHT] == -1)) {
     t2a[0][jm-1]=0.0;
   }
   if ((gp[pid].neighbors[DOWN] == -1) && (gp[pid].neighbors[RIGHT] == -1)) {
     t2a[im-1][jm-1]=0.0;
   }

   t2a = x[pid][psiindex];
   jj = gp[pid].neighbors[UPLEFT];
   if (jj != -1) {
     t2a[0][0]=x[jj][psiindex][im-2][jm-2];
   }
   jj = gp[pid].neighbors[UPRIGHT];
   if (jj != -1) {
     t2a[0][jm-1]=x[jj][psiindex][im-2][1];
   }
   jj = gp[pid].neighbors[DOWNLEFT];
   if (jj != -1) {
     t2a[im-1][0]=x[jj][psiindex][1][jm-2];
   }
   jj = gp[pid].neighbors[DOWNRIGHT];
   if (jj != -1) {
     t2a[im-1][jm-1]=x[jj][psiindex][1][1];
   }

   t2a = y[pid][psiindex];
   jj = gp[pid].neighbors[UPLEFT];
   if (jj != -1) {
     t2a[0][0]=y[jj][psiindex][im-2][jm-2];
   }
   jj = gp[pid].neighbors[UPRIGHT];
   if (jj != -1) {
     t2a[0][jm-1]=y[jj][psiindex][im-2][1];
   }
   jj = gp[pid].neighbors[DOWNLEFT];
   if (jj != -1) {
     t2a[im-1][0]=y[jj][psiindex][1][jm-2];
   }
   jj = gp[pid].neighbors[DOWNRIGHT];
   if (jj != -1) {
     t2a[im-1][jm-1]=y[jj][psiindex][1][1];
   }

   t2a = x[pid][psiindex];
   if (gp[pid].neighbors[UP] == -1) {
     jj = gp[pid].neighbors[LEFT];
     if (jj != -1) {
       t2a[0][0] = x[jj][psiindex][0][jm-2];
     } else {
       jj = gp[pid].neighbors[DOWN];
       if (jj != -1) {
         t2a[im-1][0] = x[jj][psiindex][1][0];
       }
     }
     jj = gp[pid].neighbors[RIGHT];
     if (jj != -1) {
       t2a[0][jm-1] = x[jj][psiindex][0][1];
     } else {
       jj = gp[pid].neighbors[DOWN];
       if (jj != -1) {
         t2a[im-1][jm-1] = x[jj][psiindex][1][jm-1];
       }
     }
   } else if (gp[pid].neighbors[DOWN] == -1) {
     jj = gp[pid].neighbors[LEFT];
     if (jj != -1) {
       t2a[im-1][0] = x[jj][psiindex][im-1][jm-2];
     } else {
       jj = gp[pid].neighbors[UP];
       if (jj != -1) {
         t2a[0][0] = x[jj][psiindex][im-2][0];
       }
     }
     jj = gp[pid].neighbors[RIGHT];
     if (jj != -1) {
       t2a[im-1][jm-1] = x[jj][psiindex][im-1][1];
     } else {
       jj = gp[pid].neighbors[UP];
       if (jj != -1) {
         t2a[0][jm-1] = x[jj][psiindex][im-2][jm-1];
       }
     }
   } else if (gp[pid].neighbors[LEFT] == -1) {
     jj = gp[pid].neighbors[UP];
     if (jj != -1) {
       t2a[0][0] = x[jj][psiindex][im-2][0];
     }
     jj = gp[pid].neighbors[DOWN];
     if (jj != -1) {
       t2a[im-1][0] = x[jj][psiindex][1][0];
     }
   } else if (gp[pid].neighbors[RIGHT] == -1) {
     jj = gp[pid].neighbors[UP];
     if (jj != -1) {
       t2a[0][jm-1] = x[jj][psiindex][im-2][jm-1];
     }
     jj = gp[pid].neighbors[DOWN];
     if (jj != -1) {
       t2a[im-1][jm-1] = x[jj][psiindex][1][jm-1];
     }
   }

   t2a = y[pid][psiindex];
   if (gp[pid].neighbors[UP] == -1) {
     jj = gp[pid].neighbors[LEFT];
     if (jj != -1) {
       t2a[0][0] = y[jj][psiindex][0][jm-2];
     } else {
       jj = gp[pid].neighbors[DOWN];
       if (jj != -1) {
         t2a[im-1][0] = y[jj][psiindex][1][0];
       }
     }
     jj = gp[pid].neighbors[RIGHT];
     if (jj != -1) {
       t2a[0][jm-1] = y[jj][psiindex][0][1];
     } else {
       jj = gp[pid].neighbors[DOWN];
       if (jj != -1) {
         t2a[im-1][jm-1] = y[jj][psiindex][1][jm-1];
       }
     }
   } else if (gp[pid].neighbors[DOWN] == -1) {
     jj = gp[pid].neighbors[LEFT];
     if (jj != -1) {
       t2a[im-1][0] = y[jj][psiindex][im-1][jm-2];
     } else {
       jj = gp[pid].neighbors[UP];
       if (jj != -1) {
         t2a[0][0] = y[jj][psiindex][im-2][0];
       }
     }
     jj = gp[pid].neighbors[RIGHT];
     if (jj != -1) {
       t2a[im-1][jm-1] = y[jj][psiindex][im-1][1];
     } else {
       jj = gp[pid].neighbors[UP];
       if (jj != -1) {
         t2a[0][jm-1] = y[jj][psiindex][im-2][jm-1];
       }
     }
   } else if (gp[pid].neighbors[LEFT] == -1) {
     jj = gp[pid].neighbors[UP];
     if (jj != -1) {
       t2a[0][0] = y[jj][psiindex][im-2][0];
     }
     jj = gp[pid].neighbors[DOWN];
     if (jj != -1) {
       t2a[im-1][0] = y[jj][psiindex][1][0];
     }
   } else if (gp[pid].neighbors[RIGHT] == -1) {
     jj = gp[pid].neighbors[UP];
     if (jj != -1) {
       t2a[0][jm-1] = y[jj][psiindex][im-2][jm-1];
     }
     jj = gp[pid].neighbors[DOWN];
     if (jj != -1) {
       t2a[im-1][jm-1] = y[jj][psiindex][1][jm-1];
     }
   }

   t2a = y[pid][psiindex];
   j = gp[pid].neighbors[UP];
   if (j != -1) {
     t1a = (double *) t2a[0];
     t1b = (double *) y[j][psiindex][im-2];
     for (i=1;i<=lastcol;i++) {
       t1a[i] = t1b[i];
     }
   }
   j = gp[pid].neighbors[DOWN];
   if (j != -1) {
     t1a = (double *) t2a[im-1];
     t1b = (double *) y[j][psiindex][1];
     for (i=1;i<=lastcol;i++) {
       t1a[i] = t1b[i];
     }
   }
   j = gp[pid].neighbors[LEFT];
   if (j != -1) {
     t2b = y[j][psiindex];
     for (i=1;i<=lastrow;i++) {
       t2a[i][0] = t2b[i][jm-2];
     }
   }
   j = gp[pid].neighbors[RIGHT];
   if (j != -1) {
     t2b = y[j][psiindex];
     for (i=1;i<=lastrow;i++) {
       t2a[i][jm-1] = t2b[i][1];
     }
   }

   t2a = x[pid][psiindex];
   j = gp[pid].neighbors[UP];
   if (j != -1) {
     t1a = (double *) t2a[0];
     t1b = (double *) x[j][psiindex][im-2];
     for (i=1;i<=lastcol;i++) {
       t1a[i] = t1b[i];
     }
   }
   j = gp[pid].neighbors[DOWN];
   if (j != -1) {
     t1a = (double *) t2a[im-1];
     t1b = (double *) x[j][psiindex][1];
     for (i=1;i<=lastcol;i++) {
       t1a[i] = t1b[i];
     }
   }
   j = gp[pid].neighbors[LEFT];
   if (j != -1) {
     t2b = x[j][psiindex];
     for (i=1;i<=lastrow;i++) {
       t2a[i][0] = t2b[i][jm-2];
     }
   }
   j = gp[pid].neighbors[RIGHT];
   if (j != -1) {
     t2b = x[j][psiindex];
     for (i=1;i<=lastrow;i++) {
       t2a[i][jm-1] = t2b[i][1];
     }
   }

   t2a = x[pid][psiindex];
   t2b = y[pid][psiindex];
   t2c = z[pid][psiindex];
   for (i=firstrow;i<=lastrow;i++) {
     ip1 = i+1;
     im1 = i-1;
     t1a = (double *) t2a[i];
     t1b = (double *) t2b[i];
     t1c = (double *) t2c[i];
     t1d = (double *) t2b[ip1];
     t1e = (double *) t2b[im1];
     t1f = (double *) t2a[ip1];
     t1g = (double *) t2a[im1];
     for (iindex=firstcol;iindex<=lastcol;iindex++) {
       indexp1 = iindex+1;
       indexm1 = iindex-1;
       f1 = (t1b[indexm1]+t1d[indexm1]-
             t1b[indexp1]-t1d[indexp1])*
            (t1f[iindex]-t1a[iindex]);
       f2 = (t1e[indexm1]+t1b[indexm1]-
             t1e[indexp1]-t1b[indexp1])*
            (t1a[iindex]-t1g[iindex]);
       f3 = (t1d[iindex]+t1d[indexp1]-
             t1e[iindex]-t1e[indexp1])*
            (t1a[indexp1]-t1a[iindex]);
       f4 = (t1d[indexm1]+t1d[iindex]-
             t1e[indexm1]-t1e[iindex])*
            (t1a[iindex]-t1a[indexm1]);
       f5 = (t1d[iindex]-t1b[indexp1])*
            (t1f[indexp1]-t1a[iindex]);
       f6 = (t1b[indexm1]-t1e[iindex])*
            (t1a[iindex]-t1g[indexm1]);
       f7 = (t1b[indexp1]-t1e[iindex])*
            (t1g[indexp1]-t1a[iindex]);
       f8 = (t1d[iindex]-t1b[indexm1])*
            (t1a[iindex]-t1f[indexm1]);

       t1c[iindex] = factjacob*(f1+f2+f3+f4+f5+f6+f7+f8);

     }
   }

   if (gp[pid].neighbors[UP] == -1) {
     t1c = (double *) t2c[0];
     for (j=firstcol;j<=lastcol;j++) {
       t1c[j] = 0.0;
     }
   }
   if (gp[pid].neighbors[DOWN] == -1) {
     t1c = (double *) t2c[im-1];
     for (j=firstcol;j<=lastcol;j++) {
       t1c[j] = 0.0;
     }
   }
   if (gp[pid].neighbors[LEFT] == -1) {
     for (j=firstrow;j<=lastrow;j++) {
       t2c[j][0] = 0.0;
     }
   }
   if (gp[pid].neighbors[RIGHT] == -1) {
     for (j=firstrow;j<=lastrow;j++) {
       t2c[j][jm-1] = 0.0;
     }
   }

}

#line 228 "/home/jyy/Dev/splash2/codes/null_macros/c.m4.null.POSIX"

#line 1 "jacobcalc.C"
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

/* Does the arakawa jacobian calculation (of the x and y matrices,
   putting the results in the z matrix) for a subblock.  */


#line 20
#include <pthread.h>
#line 20
#include <sys/time.h>
#line 20
#include <unistd.h>
#line 20
#include <stdlib.h>
#line 20
extern pthread_t PThreadTable[];
#line 20


#include <stdio.h>
#include <math.h>
#include <time.h>
#include "decs.h"

void jacobcalc(double ***x, double ***y, double ***z, long pid, long firstrow, long lastrow, long firstcol, long lastcol)
{
   double f1;
   double f2;
   double f3;
   double f4;
   double f5;
   double f6;
   double f7;
   double f8;
   long iindex;
   long indexp1;
   long indexm1;
   long im1;
   long ip1;
   long i;
   long j;
   long jj;
   double **t2a;
   double **t2b;
   double **t2c;
   double *t1a;
   double *t1b;
   double *t1c;
   double *t1d;
   double *t1e;
   double *t1f;
   double *t1g;

   t2a = (double **) z[pid];
   if ((gp[pid].neighbors[UP] == -1) && (gp[pid].neighbors[LEFT] == -1)) {
     t2a[0][0]=0.0;
   }
   if ((gp[pid].neighbors[DOWN] == -1) && (gp[pid].neighbors[LEFT] == -1)) {
     t2a[im-1][0]=0.0;
   }
   if ((gp[pid].neighbors[UP] == -1) && (gp[pid].neighbors[RIGHT] == -1)) {
     t2a[0][jm-1]=0.0;
   }
   if ((gp[pid].neighbors[DOWN] == -1) && (gp[pid].neighbors[RIGHT] == -1)) {
     t2a[im-1][jm-1]=0.0;
   }

   t2a = (double **) x[pid];
   jj = gp[pid].neighbors[UPLEFT];
   if (jj != -1) {
     t2a[0][0]=x[jj][im-2][jm-2];
   }
   jj = gp[pid].neighbors[UPRIGHT];
   if (jj != -1) {
     t2a[0][jm-1]=x[jj][im-2][1];
   }
   jj = gp[pid].neighbors[DOWNLEFT];
   if (jj != -1) {
     t2a[im-1][0]=x[jj][1][jm-2];
   }
   jj = gp[pid].neighbors[DOWNRIGHT];
   if (jj != -1) {
     t2a[im-1][jm-1]=x[jj][1][1];
   }

   t2a = (double **) y[pid];
   jj = gp[pid].neighbors[UPLEFT];
   if (jj != -1) {
     t2a[0][0]=y[jj][im-2][jm-2];
   }
   jj = gp[pid].neighbors[UPRIGHT];
   if (jj != -1) {
     t2a[0][jm-1]=y[jj][im-2][1];
   }
   jj = gp[pid].neighbors[DOWNLEFT];
   if (jj != -1) {
     t2a[im-1][0]=y[jj][1][jm-2];
   }
   jj = gp[pid].neighbors[DOWNRIGHT];
   if (jj != -1) {
     t2a[im-1][jm-1]=y[jj][1][1];
   }

   t2a = (double **) x[pid];
   if (gp[pid].neighbors[UP] == -1) {
     jj = gp[pid].neighbors[LEFT];
     if (jj != -1) {
       t2a[0][0] = x[jj][0][jm-2];
     } else {
       jj = gp[pid].neighbors[DOWN];
       if (jj != -1) {
         t2a[im-1][0] = x[jj][1][0];
       }
     }
     jj = gp[pid].neighbors[RIGHT];
     if (jj != -1) {
       t2a[0][jm-1] = x[jj][0][1];
     } else {
       jj = gp[pid].neighbors[DOWN];
       if (jj != -1) {
         t2a[im-1][jm-1] = x[jj][1][jm-1];
       }
     }
   } else if (gp[pid].neighbors[DOWN] == -1) {
     jj = gp[pid].neighbors[LEFT];
     if (jj != -1) {
       t2a[im-1][0] = x[jj][im-1][jm-2];
     } else {
       jj = gp[pid].neighbors[UP];
       if (jj != -1) {
         t2a[0][0] = x[jj][im-2][0];
       }
     }
     jj = gp[pid].neighbors[RIGHT];
     if (jj != -1) {
       t2a[im-1][jm-1] = x[jj][im-1][1];
     } else {
       jj = gp[pid].neighbors[UP];
       if (jj != -1) {
         t2a[0][jm-1] = x[jj][im-2][jm-1];
       }
     }
   } else if (gp[pid].neighbors[LEFT] == -1) {
     jj = gp[pid].neighbors[UP];
     if (jj != -1) {
       t2a[0][0] = x[jj][im-2][0];
     }
     jj = gp[pid].neighbors[DOWN];
     if (jj != -1) {
       t2a[im-1][0] = x[jj][1][0];
     }
   } else if (gp[pid].neighbors[RIGHT] == -1) {
     jj = gp[pid].neighbors[UP];
     if (jj != -1) {
       t2a[0][jm-1] = x[jj][im-2][jm-1];
     }
     jj = gp[pid].neighbors[DOWN];
     if (jj != -1) {
       t2a[im-1][jm-1] = x[jj][1][jm-1];
     }
   }

   t2a = (double **) y[pid];
   if (gp[pid].neighbors[UP] == -1) {
     jj = gp[pid].neighbors[LEFT];
     if (jj != -1) {
       t2a[0][0] = y[jj][0][jm-2];
     } else {
       jj = gp[pid].neighbors[DOWN];
       if (jj != -1) {
         t2a[im-1][0] = y[jj][1][0];
       }
     }
     jj = gp[pid].neighbors[RIGHT];
     if (jj != -1) {
       t2a[0][jm-1] = y[jj][0][1];
     } else {
       jj = gp[pid].neighbors[DOWN];
       if (jj != -1) {
         t2a[im-1][jm-1] = y[jj][1][jm-1];
       }
     }
   } else if (gp[pid].neighbors[DOWN] == -1) {
     jj = gp[pid].neighbors[LEFT];
     if (jj != -1) {
       t2a[im-1][0] = y[jj][im-1][jm-2];
     } else {
       jj = gp[pid].neighbors[UP];
       if (jj != -1) {
         t2a[0][0] = y[jj][im-2][0];
       }
     }
     jj = gp[pid].neighbors[RIGHT];
     if (jj != -1) {
       t2a[im-1][jm-1] = y[jj][im-1][1];
     } else {
       jj = gp[pid].neighbors[UP];
       if (jj != -1) {
         t2a[0][jm-1] = y[jj][im-2][jm-1];
       }
     }
   } else if (gp[pid].neighbors[LEFT] == -1) {
     jj = gp[pid].neighbors[UP];
     if (jj != -1) {
       t2a[0][0] = y[jj][im-2][0];
     }
     jj = gp[pid].neighbors[DOWN];
     if (jj != -1) {
       t2a[im-1][0] = y[jj][1][0];
     }
   } else if (gp[pid].neighbors[RIGHT] == -1) {
     jj = gp[pid].neighbors[UP];
     if (jj != -1) {
       t2a[0][jm-1] = y[jj][im-2][jm-1];
     }
     jj = gp[pid].neighbors[DOWN];
     if (jj != -1) {
       t2a[im-1][jm-1] = y[jj][1][jm-1];
     }
   }

   j = gp[pid].neighbors[UP];
   if (j != -1) {
     t1a = (double *) t2a[0];
     t1b = (double *) y[j][im-2];
     for (i=1;i<=lastcol;i++) {
       t1a[i] = t1b[i];
     }
   }
   j = gp[pid].neighbors[DOWN];
   if (j != -1) {
     t1a = (double *) t2a[im-1];
     t1b = (double *) y[j][1];
     for (i=1;i<=lastcol;i++) {
       t1a[i] = t1b[i];
     }
   }
   j = gp[pid].neighbors[LEFT];
   if (j != -1) {
     t2b = (double **) y[j];
     for (i=1;i<=lastrow;i++) {
       t2a[i][0] = t2b[i][jm-2];
     }
   }
   j = gp[pid].neighbors[RIGHT];
   if (j != -1) {
     t2b = (double **) y[j];
     for (i=1;i<=lastrow;i++) {
       t2a[i][jm-1] = t2b[i][1];
     }
   }

   t2a = (double **) x[pid];
   j = gp[pid].neighbors[UP];
   if (j != -1) {
     t1a = (double *) t2a[0];
     t1b = (double *) x[j][im-2];
     for (i=1;i<=lastcol;i++) {
       t1a[i] = t1b[i];
     }
   }
   j = gp[pid].neighbors[DOWN];
   if (j != -1) {
     t1a = (double *) t2a[im-1];
     t1b = (double *) x[j][1];
     for (i=1;i<=lastcol;i++) {
       t1a[i] = t1b[i];
     }
   }
   j = gp[pid].neighbors[LEFT];
   if (j != -1) {
     t2b = (double **) x[j];
     for (i=1;i<=lastrow;i++) {
       t2a[i][0] = t2b[i][jm-2];
     }
   }
   j = gp[pid].neighbors[RIGHT];
   if (j != -1) {
     t2b = (double **) x[j];
     for (i=1;i<=lastrow;i++) {
       t2a[i][jm-1] = t2b[i][1];
     }
   }

   t2a = (double **) x[pid];
   t2b = (double **) y[pid];
   t2c = (double **) z[pid];
   for (i=firstrow;i<=lastrow;i++) {
     ip1 = i+1;
     im1 = i-1;
     t1a = (double *) t2a[i];
     t1b = (double *) t2b[i];
     t1c = (double *) t2c[i];
     t1d = (double *) t2b[ip1];
     t1e = (double *) t2b[im1];
     t1f = (double *) t2a[ip1];
     t1g = (double *) t2a[im1];
     for (iindex=firstcol;iindex<=lastcol;iindex++) {
       indexp1 = iindex+1;
       indexm1 = iindex-1;
       f1 = (t1b[indexm1]+t1d[indexm1]-
             t1b[indexp1]-t1d[indexp1])*
            (t1f[iindex]-t1a[iindex]);
       f2 = (t1e[indexm1]+t1b[indexm1]-
             t1e[indexp1]-t1b[indexp1])*
            (t1a[iindex]-t1g[iindex]);
       f3 = (t1d[iindex]+t1d[indexp1]-
             t1e[iindex]-t1e[indexp1])*
            (t1a[indexp1]-t1a[iindex]);
       f4 = (t1d[indexm1]+t1d[iindex]-
             t1e[indexm1]-t1e[iindex])*
            (t1a[iindex]-t1a[indexm1]);
       f5 = (t1d[iindex]-t1b[indexp1])*
            (t1f[indexp1]-t1a[iindex]);
       f6 = (t1b[indexm1]-t1e[iindex])*
            (t1a[iindex]-t1g[indexm1]);
       f7 = (t1b[indexp1]-t1e[iindex])*
            (t1g[indexp1]-t1a[iindex]);
       f8 = (t1d[iindex]-t1b[indexm1])*
            (t1a[iindex]-t1f[indexm1]);

       t1c[iindex] = factjacob*(f1+f2+f3+f4+f5+f6+f7+f8);
     }
   }

   if (gp[pid].neighbors[UP] == -1) {
     t1c = (double *) t2c[0];
     for (j=firstcol;j<=lastcol;j++) {
       t1c[j] = 0.0;
     }
   }
   if (gp[pid].neighbors[DOWN] == -1) {
     t1c = (double *) t2c[im-1];
     for (j=firstcol;j<=lastcol;j++) {
       t1c[j] = 0.0;
     }
   }
   if (gp[pid].neighbors[LEFT] == -1) {
     for (j=firstrow;j<=lastrow;j++) {
       t2c[j][0] = 0.0;
     }
   }
   if (gp[pid].neighbors[RIGHT] == -1) {
     for (j=firstrow;j<=lastrow;j++) {
       t2c[j][jm-1] = 0.0;
     }
   }

}

#line 228 "/home/jyy/Dev/splash2/codes/null_macros/c.m4.null.POSIX"

#line 1 "laplacalc.C"
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

/* Performs the laplacian calculation for a subblock */


#line 19
#include <pthread.h>
#line 19
#include <sys/time.h>
#line 19
#include <unistd.h>
#line 19
#include <stdlib.h>
#line 19
extern pthread_t PThreadTable[];
#line 19


#include <stdio.h>
#include <math.h>
#include <time.h>
#include "decs.h"

void laplacalc(long procid, double ****x, double ****z, long psiindex, long firstrow, long lastrow, long firstcol, long lastcol)
{
   long iindex;
   long indexp1;
   long indexm1;
   long ip1;
   long im1;
   long i;
   long j;
   double **t2a;
   double **t2b;
   double *t1a;
   double *t1b;
   double *t1c;
   double *t1d;

   t2a = (double **) x[procid][psiindex];
   j = gp[procid].neighbors[UP];
   if (j != -1) {
     t1a = (double *) t2a[0];
     t1b = (double *) x[j][psiindex][im-2];
     for (i=1;i<=lastcol;i++) {
       t1a[i] = t1b[i];
     }
   }
   j = gp[procid].neighbors[DOWN];
   if (j != -1) {
     t1a = (double *) t2a[im-1];
     t1b = (double *) x[j][psiindex][1];
     for (i=1;i<=lastcol;i++) {
       t1a[i] = t1b[i];
     }
   }
   j = gp[procid].neighbors[LEFT];
   if (j != -1) {
     t2b = (double **) x[j][psiindex];
     for (i=1;i<=lastrow;i++) {
       t2a[i][0] = t2b[i][jm-2];
     }
   }
   j = gp[procid].neighbors[RIGHT];
   if (j != -1) {
     t2b = (double **) x[j][psiindex];
     for (i=1;i<=lastrow;i++) {
       t2a[i][jm-1] = t2b[i][1];
     }
   }

   t2a = (double **) x[procid][psiindex];
   t2b = (double **) z[procid][psiindex];
   for (i=firstrow;i<=lastrow;i++) {
     ip1 = i+1;
     im1 = i-1;
     t1a = (double *) t2a[i];
     t1b = (double *) t2b[i];
     t1c = (double *) t2a[ip1];
     t1d = (double *) t2a[im1];
     for (iindex=firstcol;iindex<=lastcol;iindex++) {
       indexp1 = iindex+1;
       indexm1 = iindex-1;
       t1b[iindex] = factlap*(t1c[iindex]+
			      t1d[iindex]+t1a[indexp1]+
                              t1a[indexm1]-4.*t1a[iindex]);
     }
   }

   if (gp[procid].neighbors[UP] == -1) {
     t1b = (double *) t2b[0];
     for (j=firstcol;j<=lastcol;j++) {
       t1b[j] = 0.0;
     }
   }
   if (gp[procid].neighbors[DOWN] == -1) {
     t1b = (double *) t2b[im-1];
     for (j=firstcol;j<=lastcol;j++) {
       t1b[j] = 0.0;
     }
   }
   if (gp[procid].neighbors[LEFT] == -1) {
     for (j=firstrow;j<=lastrow;j++) {
       t2b[j][0] = 0.0;
     }
   }
   if (gp[procid].neighbors[RIGHT] == -1) {
     for (j=firstrow;j<=lastrow;j++) {
       t2b[j][jm-1] = 0.0;
     }
   }

}
#line 228 "/home/jyy/Dev/splash2/codes/null_macros/c.m4.null.POSIX"

#line 1 "linkup.C"
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

/* Set all the pointers to the proper locations for the q_multi and
   rhs_multi data structures */


#line 20
#include <pthread.h>
#line 20
#include <sys/time.h>
#line 20
#include <unistd.h>
#line 20
#include <stdlib.h>
#line 20
extern pthread_t PThreadTable[];
#line 20


#include "decs.h"

void link_all()
{
  long i;
  long j;

  for (j=0;j<nprocs;j++) {
    linkup(psium[j]);
    linkup(psilm[j]);
    linkup(psib[j]);
    linkup(ga[j]);
    linkup(gb[j]);
    linkup(work2[j]);
    linkup(work3[j]);
    linkup(work6[j]);
    linkup(tauz[j]);
    linkup(oldga[j]);
    linkup(oldgb[j]);
    for (i=0;i<=1;i++) {
      linkup(psi[j][i]);
      linkup(psim[j][i]);
      linkup(work1[j][i]);
      linkup(work4[j][i]);
      linkup(work5[j][i]);
      linkup(work7[j][i]);
      linkup(temparray[j][i]);
    }
  }
  link_multi();
}

void linkup(double **row_ptr)
{
  long i;
  double *a;
  double **row;
  double **y;
  long x_part;
  long y_part;

  x_part = (jm-2)/xprocs + 2;
  y_part = (im-2)/yprocs + 2;
  row = row_ptr;
  y = row + y_part;
  a = (double *) y;
  for (i=0;i<y_part;i++) {
    *row = (double *) a;
    row++;
    a += x_part;
  }
}

void link_multi()
{
  long i;
  long j;
  long l;
  double *a;
  double **row;
  double **y;
  unsigned long z;
  unsigned long zz;
  long x_part;
  long y_part;
  unsigned long d_size;

  z = ((unsigned long) q_multi + nprocs*sizeof(double ***));

  if (nprocs%2 == 1) {         /* To make sure that the actual data
                                  starts double word aligned, add an extra
                                  pointer */
    z += sizeof(double ***);
  }

  d_size = numlev*sizeof(double **);
  if (numlev%2 == 1) {       /* To make sure that the actual data
                                starts double word aligned, add an extra
                                pointer */
    d_size += sizeof(double **);
  }
  for (i=0;i<numlev;i++) {
    d_size += ((imx[i]-2)/yprocs+2)*((jmx[i]-2)/xprocs+2)*sizeof(double)+
             ((imx[i]-2)/yprocs+2)*sizeof(double *);
  }
  for (i=0;i<nprocs;i++) {
    q_multi[i] = (double ***) z;
    z += d_size;
  }
  for (j=0;j<nprocs;j++) {
    zz = (unsigned long) q_multi[j];
    zz += numlev*sizeof(double **);
    if (numlev%2 == 1) {       /* To make sure that the actual data
                                  starts double word aligned, add an extra
                                  pointer */
      zz += sizeof(double **);
    }
    for (i=0;i<numlev;i++) {
      d_size = ((imx[i]-2)/yprocs+2)*((jmx[i]-2)/xprocs+2)*sizeof(double)+
               ((imx[i]-2)/yprocs+2)*sizeof(double *);
      q_multi[j][i] = (double **) zz;
      zz += d_size;
    }
  }

  for (l=0;l<numlev;l++) {
    x_part = (jmx[l]-2)/xprocs + 2;
    y_part = (imx[l]-2)/yprocs + 2;
    for (j=0;j<nprocs;j++) {
      row = q_multi[j][l];
      y = row + y_part;
      a = (double *) y;
      for (i=0;i<y_part;i++) {
        *row = (double *) a;
        row++;
        a += x_part;
      }
    }
  }

  z = ((unsigned long) rhs_multi + nprocs*sizeof(double ***));
  if (nprocs%2 == 1) {         /* To make sure that the actual data
                                  starts double word aligned, add an extra
                                  pointer */
    z += sizeof(double ***);
  }

  d_size = numlev*sizeof(double **);
  if (numlev%2 == 1) {       /* To make sure that the actual data
                                starts double word aligned, add an extra
                                pointer */
    d_size += sizeof(double **);
  }
  for (i=0;i<numlev;i++) {
    d_size += ((imx[i]-2)/yprocs+2)*((jmx[i]-2)/xprocs+2)*sizeof(double)+
             ((imx[i]-2)/yprocs+2)*sizeof(double *);
  }
  for (i=0;i<nprocs;i++) {
    rhs_multi[i] = (double ***) z;
    z += d_size;
  }
  for (j=0;j<nprocs;j++) {
    zz = (unsigned long) rhs_multi[j];
    zz += numlev*sizeof(double **);
    if (numlev%2 == 1) {       /* To make sure that the actual data
                                  starts double word aligned, add an extra
                                  pointer */
      zz += sizeof(double **);
    }
    for (i=0;i<numlev;i++) {
      d_size = ((imx[i]-2)/yprocs+2)*((jmx[i]-2)/xprocs+2)*sizeof(double)+
               ((imx[i]-2)/yprocs+2)*sizeof(double *);
      rhs_multi[j][i] = (double **) zz;
      zz += d_size;
    }
  }

  for (l=0;l<numlev;l++) {
    x_part = (jmx[l]-2)/xprocs + 2;
    y_part = (imx[l]-2)/yprocs + 2;
    for (j=0;j<nprocs;j++) {
      row = rhs_multi[j][l];
      y = row + y_part;
      a = (double *) y;
      for (i=0;i<y_part;i++) {
        *row = (double *) a;
        row++;
        a += x_part;
      }
    }
  }

}


#line 228 "/home/jyy/Dev/splash2/codes/null_macros/c.m4.null.POSIX"

#line 1 "main.C"
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

/*************************************************************************/
/*                                                                       */
/*  SPLASH Ocean Code                                                    */
/*                                                                       */
/*  This application studies the role of eddy and boundary currents in   */
/*  influencing large-scale ocean movements.  This implementation uses   */
/*  dynamically allocated four-dimensional arrays for grid data storage. */
/*                                                                       */
/*  Command line options:                                                */
/*                                                                       */
/*     -nN : Simulate NxN ocean.  N must be (power of 2)+2.              */
/*     -pP : P = number of processors.  P must be power of 2.            */
/*     -eE : E = error tolerance for iterative relaxation.               */
/*     -rR : R = distance between grid points in meters.                 */
/*     -tT : T = timestep in seconds.                                    */
/*     -s  : Print timing statistics.                                    */
/*     -o  : Print out relaxation residual values.                       */
/*     -h  : Print out command line options.                             */
/*                                                                       */
/*  Default: OCEAN -n130 -p1 -e1e-7 -r20000.0 -t28800.0                  */
/*                                                                       */
/*  NOTE: This code works under both the FORK and SPROC models.          */
/*                                                                       */
/*************************************************************************/


#line 42
#include <pthread.h>
#line 42
#include <sys/time.h>
#line 42
#include <unistd.h>
#line 42
#include <stdlib.h>
#line 42
#define MAX_THREADS 32
#line 42
pthread_t PThreadTable[MAX_THREADS];
#line 42


#define DEFAULT_N      258
#define DEFAULT_P        1
#define DEFAULT_E        1e-7
#define DEFAULT_T    28800.0
#define DEFAULT_R    20000.0
#define UP               0
#define DOWN             1
#define LEFT             2
#define RIGHT            3
#define UPLEFT           4
#define UPRIGHT          5
#define DOWNLEFT         6
#define DOWNRIGHT        7
#define PAGE_SIZE     4096

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "decs.h"

struct multi_struct *multi;
struct global_struct *global;
struct locks_struct *locks;
struct bars_struct *bars;

double ****psi;
double ****psim;
double ***psium;
double ***psilm;
double ***psib;
double ***ga;
double ***gb;
double ****work1;
double ***work2;
double ***work3;
double ****work4;
double ****work5;
double ***work6;
double ****work7;
double ****temparray;
double ***tauz;
double ***oldga;
double ***oldgb;
double *f;
double ****q_multi;
double ****rhs_multi;

long nprocs = DEFAULT_P;
double h1 = 1000.0;
double h3 = 4000.0;
double h = 5000.0;
double lf = -5.12e11;
double res = DEFAULT_R;
double dtau = DEFAULT_T;
double f0 = 8.3e-5;
double beta = 2.0e-11;
double gpr = 0.02;
long im = DEFAULT_N;
long jm;
double tolerance = DEFAULT_E;
double eig2;
double ysca;
long jmm1;
double pi;
double t0 = 0.5e-4 ;
double outday0 = 1.0;
double outday1 = 2.0;
double outday2 = 2.0;
double outday3 = 2.0;
double factjacob;
double factlap;
long numlev;
long *imx;
long *jmx;
double *lev_res;
double *lev_tol;
double maxwork = 10000.0;

struct Global_Private *gp;

double *i_int_coeff;
double *j_int_coeff;
long xprocs;
long yprocs;
long *xpts_per_proc;
long *ypts_per_proc;
long minlevel;
long do_stats = 0;
long do_output = 0;

int main(int argc, char *argv[])
{
   long i;
   long j;
   long k;
   long x_part;
   long y_part;
   long d_size;
   long itemp;
   long jtemp;
   double procsqrt;
   long temp = 0;
   double min_total;
   double max_total;
   double avg_total;
   double min_multi;
   double max_multi;
   double avg_multi;
   double min_frac;
   double max_frac;
   double avg_frac;
   long ch;
   extern char *optarg;
   unsigned long computeend;
   unsigned long start;

   {
#line 161
	struct timeval	FullTime;
#line 161

#line 161
	gettimeofday(&FullTime, NULL);
#line 161
	(start) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 161
}

   while ((ch = getopt(argc, argv, "n:p:e:r:t:soh")) != -1) {
     switch(ch) {
     case 'n': im = atoi(optarg);
               if (log_2(im-2) == -1) {
                 printerr("Grid must be ((power of 2)+2) in each dimension\n");
                 exit(-1);
               }
               break;
     case 'p': nprocs = atoi(optarg);
               if (nprocs < 1) {
                 printerr("P must be >= 1\n");
                 exit(-1);
               }
               if (log_2(nprocs) == -1) {
                 printerr("P must be a power of 2\n");
                 exit(-1);
               }
               break;
     case 'e': tolerance = atof(optarg); break;
     case 'r': res = atof(optarg); break;
     case 't': dtau = atof(optarg); break;
     case 's': do_stats = !do_stats; break;
     case 'o': do_output = !do_output; break;
     case 'h': printf("Usage: OCEAN <options>\n\n");
               printf("options:\n");
               printf("  -nN : Simulate NxN ocean.  N must be (power of 2)+2.\n");
               printf("  -pP : P = number of processors.  P must be power of 2.\n");
               printf("  -eE : E = error tolerance for iterative relaxation.\n");
               printf("  -rR : R = distance between grid points in meters.\n");
               printf("  -tT : T = timestep in seconds.\n");
               printf("  -s  : Print timing statistics.\n");
               printf("  -o  : Print out relaxation residual values.\n");
               printf("  -h  : Print out command line options.\n\n");
               printf("Default: OCEAN -n%1d -p%1d -e%1g -r%1g -t%1g\n",
                       DEFAULT_N,DEFAULT_P,DEFAULT_E,DEFAULT_R,DEFAULT_T);
               exit(0);
               break;
     }
   }

   {;}

   jm = im;
   printf("\n");
   printf("Ocean simulation with W-cycle multigrid solver\n");
   printf("    Processors                         : %1ld\n",nprocs);
   printf("    Grid size                          : %1ld x %1ld\n",im,jm);
   printf("    Grid resolution (meters)           : %0.2f\n",res);
   printf("    Time between relaxations (seconds) : %0.0f\n",dtau);
   printf("    Error tolerance                    : %0.7g\n",tolerance);
   printf("\n");

   xprocs = 0;
   yprocs = 0;
   procsqrt = sqrt((double) nprocs);
   j = (long) procsqrt;
   while ((xprocs == 0) && (j > 0)) {
     k = nprocs / j;
     if (k * j == nprocs) {
       if (k > j) {
         xprocs = j;
         yprocs = k;
       } else {
         xprocs = k;
         yprocs = j;
       }
     }
     j--;
   }
   if (xprocs == 0) {
     printerr("Could not find factors for subblocking\n");
     exit(-1);
   }

   minlevel = 0;
   itemp = 1;
   jtemp = 1;
   numlev = 0;
   minlevel = 0;
   while (itemp < (im-2)) {
     itemp = itemp*2;
     jtemp = jtemp*2;
     if ((itemp/yprocs > 1) && (jtemp/xprocs > 1)) {
       numlev++;
     }
   }

   if (numlev == 0) {
     printerr("Must have at least 2 grid points per processor in each dimension\n");
     exit(-1);
   }

   imx = (long *) valloc(numlev*sizeof(long));;
   jmx = (long *) valloc(numlev*sizeof(long));;
   lev_res = (double *) valloc(numlev*sizeof(double));;
   lev_tol = (double *) valloc(numlev*sizeof(double));;
   i_int_coeff = (double *) valloc(numlev*sizeof(double));;
   j_int_coeff = (double *) valloc(numlev*sizeof(double));;
   xpts_per_proc = (long *) valloc(numlev*sizeof(long));;
   ypts_per_proc = (long *) valloc(numlev*sizeof(long));;

   imx[numlev-1] = im;
   jmx[numlev-1] = jm;
   lev_res[numlev-1] = res;
   lev_tol[numlev-1] = tolerance;

   for (i=numlev-2;i>=0;i--) {
     imx[i] = ((imx[i+1] - 2) / 2) + 2;
     jmx[i] = ((jmx[i+1] - 2) / 2) + 2;
     lev_res[i] = lev_res[i+1] * 2;
   }

   for (i=0;i<numlev;i++) {
     xpts_per_proc[i] = (jmx[i]-2) / xprocs;
     ypts_per_proc[i] = (imx[i]-2) / yprocs;
   }
   for (i=numlev-1;i>=0;i--) {
     if ((xpts_per_proc[i] < 2) || (ypts_per_proc[i] < 2)) {
       minlevel = i+1;
       break;
     }
   }

   for (i=0;i<numlev;i++) {
     temp += imx[i];
   }
   temp = 0;
   j = 0;
   for (k=0;k<numlev;k++) {
     for (i=0;i<imx[k];i++) {
       j++;
       temp += jmx[k];
     }
   }

   d_size = nprocs*sizeof(double ***);
   psi = (double ****) valloc(d_size);;
   psim = (double ****) valloc(d_size);;
   work1 = (double ****) valloc(d_size);;
   work4 = (double ****) valloc(d_size);;
   work5 = (double ****) valloc(d_size);;
   work7 = (double ****) valloc(d_size);;
   temparray = (double ****) valloc(d_size);;

   d_size = 2*sizeof(double **);
   for (i=0;i<nprocs;i++) {
     psi[i] = (double ***) valloc(d_size);;
     psim[i] = (double ***) valloc(d_size);;
     work1[i] = (double ***) valloc(d_size);;
     work4[i] = (double ***) valloc(d_size);;
     work5[i] = (double ***) valloc(d_size);;
     work7[i] = (double ***) valloc(d_size);;
     temparray[i] = (double ***) valloc(d_size);;
   }

   d_size = nprocs*sizeof(double **);
   psium = (double ***) valloc(d_size);;
   psilm = (double ***) valloc(d_size);;
   psib = (double ***) valloc(d_size);;
   ga = (double ***) valloc(d_size);;
   gb = (double ***) valloc(d_size);;
   work2 = (double ***) valloc(d_size);;
   work3 = (double ***) valloc(d_size);;
   work6 = (double ***) valloc(d_size);;
   tauz = (double ***) valloc(d_size);;
   oldga = (double ***) valloc(d_size);;
   oldgb = (double ***) valloc(d_size);;

   gp = (struct Global_Private *) valloc((nprocs+1)*sizeof(struct Global_Private));;
   for (i=0;i<nprocs;i++) {
     gp[i].rel_num_x = (long *) valloc(numlev*sizeof(long));;
     gp[i].rel_num_y = (long *) valloc(numlev*sizeof(long));;
     gp[i].eist = (long *) valloc(numlev*sizeof(long));;
     gp[i].ejst = (long *) valloc(numlev*sizeof(long));;
     gp[i].oist = (long *) valloc(numlev*sizeof(long));;
     gp[i].ojst = (long *) valloc(numlev*sizeof(long));;
     gp[i].rlist = (long *) valloc(numlev*sizeof(long));;
     gp[i].rljst = (long *) valloc(numlev*sizeof(long));;
     gp[i].rlien = (long *) valloc(numlev*sizeof(long));;
     gp[i].rljen = (long *) valloc(numlev*sizeof(long));;
     gp[i].multi_time = 0;
     gp[i].total_time = 0;
   }

   subblock();

   x_part = (jm - 2)/xprocs + 2;
   y_part = (im - 2)/yprocs + 2;

   d_size = x_part*y_part*sizeof(double) + y_part*sizeof(double *);

   global = (struct global_struct *) valloc(sizeof(struct global_struct));;
   for (i=0;i<nprocs;i++) {
     psi[i][0] = (double **) valloc(d_size);;
     psi[i][1] = (double **) valloc(d_size);;
     psim[i][0] = (double **) valloc(d_size);;
     psim[i][1] = (double **) valloc(d_size);;
     psium[i] = (double **) valloc(d_size);;
     psilm[i] = (double **) valloc(d_size);;
     psib[i] = (double **) valloc(d_size);;
     ga[i] = (double **) valloc(d_size);;
     gb[i] = (double **) valloc(d_size);;
     work1[i][0] = (double **) valloc(d_size);;
     work1[i][1] = (double **) valloc(d_size);;
     work2[i] = (double **) valloc(d_size);;
     work3[i] = (double **) valloc(d_size);;
     work4[i][0] = (double **) valloc(d_size);;
     work4[i][1] = (double **) valloc(d_size);;
     work5[i][0] = (double **) valloc(d_size);;
     work5[i][1] = (double **) valloc(d_size);;
     work6[i] = (double **) valloc(d_size);;
     work7[i][0] = (double **) valloc(d_size);;
     work7[i][1] = (double **) valloc(d_size);;
     temparray[i][0] = (double **) valloc(d_size);;
     temparray[i][1] = (double **) valloc(d_size);;
     tauz[i] = (double **) valloc(d_size);;
     oldga[i] = (double **) valloc(d_size);;
     oldgb[i] = (double **) valloc(d_size);;
   }
   f = (double *) valloc(im*sizeof(double));;

   multi = (struct multi_struct *) valloc(sizeof(struct multi_struct));;

   d_size = numlev*sizeof(double **);
   if (numlev%2 == 1) {         /* To make sure that the actual data
                                   starts double word aligned, add an extra
                                   pointer */
     d_size += sizeof(double **);
   }
   for (i=0;i<numlev;i++) {
     d_size += ((imx[i]-2)/yprocs+2)*((jmx[i]-2)/xprocs+2)*sizeof(double)+
              ((imx[i]-2)/yprocs+2)*sizeof(double *);
   }

   d_size *= nprocs;

   if (nprocs%2 == 1) {         /* To make sure that the actual data
                                   starts double word aligned, add an extra
                                   pointer */
     d_size += sizeof(double ***);
   }

   d_size += nprocs*sizeof(double ***);
   q_multi = (double ****) valloc(d_size);;
   rhs_multi = (double ****) valloc(d_size);;

   locks = (struct locks_struct *) valloc(sizeof(struct locks_struct));;
   bars = (struct bars_struct *) valloc(sizeof(struct bars_struct));;

   {pthread_mutex_init(&(locks->idlock), NULL);}
   {pthread_mutex_init(&(locks->psiailock), NULL);}
   {pthread_mutex_init(&(locks->psibilock), NULL);}
   {pthread_mutex_init(&(locks->donelock), NULL);}
   {pthread_mutex_init(&(locks->error_lock), NULL);}
   {pthread_mutex_init(&(locks->bar_lock), NULL);}

#if defined(MULTIPLE_BARRIERS)
   {
#line 420
	unsigned long	Error;
#line 420

#line 420
	Error = pthread_mutex_init(&(bars->iteration).mutex, NULL);
#line 420
	if (Error != 0) {
#line 420
		printf("Error while initializing barrier.\n");
#line 420
		exit(-1);
#line 420
	}
#line 420

#line 420
	Error = pthread_cond_init(&(bars->iteration).cv, NULL);
#line 420
	if (Error != 0) {
#line 420
		printf("Error while initializing barrier.\n");
#line 420
		pthread_mutex_destroy(&(bars->iteration).mutex);
#line 420
		exit(-1);
#line 420
	}
#line 420

#line 420
	(bars->iteration).counter = 0;
#line 420
	(bars->iteration).cycle = 0;
#line 420
}
   {
#line 421
	unsigned long	Error;
#line 421

#line 421
	Error = pthread_mutex_init(&(bars->gsudn).mutex, NULL);
#line 421
	if (Error != 0) {
#line 421
		printf("Error while initializing barrier.\n");
#line 421
		exit(-1);
#line 421
	}
#line 421

#line 421
	Error = pthread_cond_init(&(bars->gsudn).cv, NULL);
#line 421
	if (Error != 0) {
#line 421
		printf("Error while initializing barrier.\n");
#line 421
		pthread_mutex_destroy(&(bars->gsudn).mutex);
#line 421
		exit(-1);
#line 421
	}
#line 421

#line 421
	(bars->gsudn).counter = 0;
#line 421
	(bars->gsudn).cycle = 0;
#line 421
}
   {
#line 422
	unsigned long	Error;
#line 422

#line 422
	Error = pthread_mutex_init(&(bars->p_setup).mutex, NULL);
#line 422
	if (Error != 0) {
#line 422
		printf("Error while initializing barrier.\n");
#line 422
		exit(-1);
#line 422
	}
#line 422

#line 422
	Error = pthread_cond_init(&(bars->p_setup).cv, NULL);
#line 422
	if (Error != 0) {
#line 422
		printf("Error while initializing barrier.\n");
#line 422
		pthread_mutex_destroy(&(bars->p_setup).mutex);
#line 422
		exit(-1);
#line 422
	}
#line 422

#line 422
	(bars->p_setup).counter = 0;
#line 422
	(bars->p_setup).cycle = 0;
#line 422
}
   {
#line 423
	unsigned long	Error;
#line 423

#line 423
	Error = pthread_mutex_init(&(bars->p_redph).mutex, NULL);
#line 423
	if (Error != 0) {
#line 423
		printf("Error while initializing barrier.\n");
#line 423
		exit(-1);
#line 423
	}
#line 423

#line 423
	Error = pthread_cond_init(&(bars->p_redph).cv, NULL);
#line 423
	if (Error != 0) {
#line 423
		printf("Error while initializing barrier.\n");
#line 423
		pthread_mutex_destroy(&(bars->p_redph).mutex);
#line 423
		exit(-1);
#line 423
	}
#line 423

#line 423
	(bars->p_redph).counter = 0;
#line 423
	(bars->p_redph).cycle = 0;
#line 423
}
   {
#line 424
	unsigned long	Error;
#line 424

#line 424
	Error = pthread_mutex_init(&(bars->p_soln).mutex, NULL);
#line 424
	if (Error != 0) {
#line 424
		printf("Error while initializing barrier.\n");
#line 424
		exit(-1);
#line 424
	}
#line 424

#line 424
	Error = pthread_cond_init(&(bars->p_soln).cv, NULL);
#line 424
	if (Error != 0) {
#line 424
		printf("Error while initializing barrier.\n");
#line 424
		pthread_mutex_destroy(&(bars->p_soln).mutex);
#line 424
		exit(-1);
#line 424
	}
#line 424

#line 424
	(bars->p_soln).counter = 0;
#line 424
	(bars->p_soln).cycle = 0;
#line 424
}
   {
#line 425
	unsigned long	Error;
#line 425

#line 425
	Error = pthread_mutex_init(&(bars->p_subph).mutex, NULL);
#line 425
	if (Error != 0) {
#line 425
		printf("Error while initializing barrier.\n");
#line 425
		exit(-1);
#line 425
	}
#line 425

#line 425
	Error = pthread_cond_init(&(bars->p_subph).cv, NULL);
#line 425
	if (Error != 0) {
#line 425
		printf("Error while initializing barrier.\n");
#line 425
		pthread_mutex_destroy(&(bars->p_subph).mutex);
#line 425
		exit(-1);
#line 425
	}
#line 425

#line 425
	(bars->p_subph).counter = 0;
#line 425
	(bars->p_subph).cycle = 0;
#line 425
}
   {
#line 426
	unsigned long	Error;
#line 426

#line 426
	Error = pthread_mutex_init(&(bars->sl_prini).mutex, NULL);
#line 426
	if (Error != 0) {
#line 426
		printf("Error while initializing barrier.\n");
#line 426
		exit(-1);
#line 426
	}
#line 426

#line 426
	Error = pthread_cond_init(&(bars->sl_prini).cv, NULL);
#line 426
	if (Error != 0) {
#line 426
		printf("Error while initializing barrier.\n");
#line 426
		pthread_mutex_destroy(&(bars->sl_prini).mutex);
#line 426
		exit(-1);
#line 426
	}
#line 426

#line 426
	(bars->sl_prini).counter = 0;
#line 426
	(bars->sl_prini).cycle = 0;
#line 426
}
   {
#line 427
	unsigned long	Error;
#line 427

#line 427
	Error = pthread_mutex_init(&(bars->sl_psini).mutex, NULL);
#line 427
	if (Error != 0) {
#line 427
		printf("Error while initializing barrier.\n");
#line 427
		exit(-1);
#line 427
	}
#line 427

#line 427
	Error = pthread_cond_init(&(bars->sl_psini).cv, NULL);
#line 427
	if (Error != 0) {
#line 427
		printf("Error while initializing barrier.\n");
#line 427
		pthread_mutex_destroy(&(bars->sl_psini).mutex);
#line 427
		exit(-1);
#line 427
	}
#line 427

#line 427
	(bars->sl_psini).counter = 0;
#line 427
	(bars->sl_psini).cycle = 0;
#line 427
}
   {
#line 428
	unsigned long	Error;
#line 428

#line 428
	Error = pthread_mutex_init(&(bars->sl_onetime).mutex, NULL);
#line 428
	if (Error != 0) {
#line 428
		printf("Error while initializing barrier.\n");
#line 428
		exit(-1);
#line 428
	}
#line 428

#line 428
	Error = pthread_cond_init(&(bars->sl_onetime).cv, NULL);
#line 428
	if (Error != 0) {
#line 428
		printf("Error while initializing barrier.\n");
#line 428
		pthread_mutex_destroy(&(bars->sl_onetime).mutex);
#line 428
		exit(-1);
#line 428
	}
#line 428

#line 428
	(bars->sl_onetime).counter = 0;
#line 428
	(bars->sl_onetime).cycle = 0;
#line 428
}
   {
#line 429
	unsigned long	Error;
#line 429

#line 429
	Error = pthread_mutex_init(&(bars->sl_phase_1).mutex, NULL);
#line 429
	if (Error != 0) {
#line 429
		printf("Error while initializing barrier.\n");
#line 429
		exit(-1);
#line 429
	}
#line 429

#line 429
	Error = pthread_cond_init(&(bars->sl_phase_1).cv, NULL);
#line 429
	if (Error != 0) {
#line 429
		printf("Error while initializing barrier.\n");
#line 429
		pthread_mutex_destroy(&(bars->sl_phase_1).mutex);
#line 429
		exit(-1);
#line 429
	}
#line 429

#line 429
	(bars->sl_phase_1).counter = 0;
#line 429
	(bars->sl_phase_1).cycle = 0;
#line 429
}
   {
#line 430
	unsigned long	Error;
#line 430

#line 430
	Error = pthread_mutex_init(&(bars->sl_phase_2).mutex, NULL);
#line 430
	if (Error != 0) {
#line 430
		printf("Error while initializing barrier.\n");
#line 430
		exit(-1);
#line 430
	}
#line 430

#line 430
	Error = pthread_cond_init(&(bars->sl_phase_2).cv, NULL);
#line 430
	if (Error != 0) {
#line 430
		printf("Error while initializing barrier.\n");
#line 430
		pthread_mutex_destroy(&(bars->sl_phase_2).mutex);
#line 430
		exit(-1);
#line 430
	}
#line 430

#line 430
	(bars->sl_phase_2).counter = 0;
#line 430
	(bars->sl_phase_2).cycle = 0;
#line 430
}
   {
#line 431
	unsigned long	Error;
#line 431

#line 431
	Error = pthread_mutex_init(&(bars->sl_phase_3).mutex, NULL);
#line 431
	if (Error != 0) {
#line 431
		printf("Error while initializing barrier.\n");
#line 431
		exit(-1);
#line 431
	}
#line 431

#line 431
	Error = pthread_cond_init(&(bars->sl_phase_3).cv, NULL);
#line 431
	if (Error != 0) {
#line 431
		printf("Error while initializing barrier.\n");
#line 431
		pthread_mutex_destroy(&(bars->sl_phase_3).mutex);
#line 431
		exit(-1);
#line 431
	}
#line 431

#line 431
	(bars->sl_phase_3).counter = 0;
#line 431
	(bars->sl_phase_3).cycle = 0;
#line 431
}
   {
#line 432
	unsigned long	Error;
#line 432

#line 432
	Error = pthread_mutex_init(&(bars->sl_phase_4).mutex, NULL);
#line 432
	if (Error != 0) {
#line 432
		printf("Error while initializing barrier.\n");
#line 432
		exit(-1);
#line 432
	}
#line 432

#line 432
	Error = pthread_cond_init(&(bars->sl_phase_4).cv, NULL);
#line 432
	if (Error != 0) {
#line 432
		printf("Error while initializing barrier.\n");
#line 432
		pthread_mutex_destroy(&(bars->sl_phase_4).mutex);
#line 432
		exit(-1);
#line 432
	}
#line 432

#line 432
	(bars->sl_phase_4).counter = 0;
#line 432
	(bars->sl_phase_4).cycle = 0;
#line 432
}
   {
#line 433
	unsigned long	Error;
#line 433

#line 433
	Error = pthread_mutex_init(&(bars->sl_phase_5).mutex, NULL);
#line 433
	if (Error != 0) {
#line 433
		printf("Error while initializing barrier.\n");
#line 433
		exit(-1);
#line 433
	}
#line 433

#line 433
	Error = pthread_cond_init(&(bars->sl_phase_5).cv, NULL);
#line 433
	if (Error != 0) {
#line 433
		printf("Error while initializing barrier.\n");
#line 433
		pthread_mutex_destroy(&(bars->sl_phase_5).mutex);
#line 433
		exit(-1);
#line 433
	}
#line 433

#line 433
	(bars->sl_phase_5).counter = 0;
#line 433
	(bars->sl_phase_5).cycle = 0;
#line 433
}
   {
#line 434
	unsigned long	Error;
#line 434

#line 434
	Error = pthread_mutex_init(&(bars->sl_phase_6).mutex, NULL);
#line 434
	if (Error != 0) {
#line 434
		printf("Error while initializing barrier.\n");
#line 434
		exit(-1);
#line 434
	}
#line 434

#line 434
	Error = pthread_cond_init(&(bars->sl_phase_6).cv, NULL);
#line 434
	if (Error != 0) {
#line 434
		printf("Error while initializing barrier.\n");
#line 434
		pthread_mutex_destroy(&(bars->sl_phase_6).mutex);
#line 434
		exit(-1);
#line 434
	}
#line 434

#line 434
	(bars->sl_phase_6).counter = 0;
#line 434
	(bars->sl_phase_6).cycle = 0;
#line 434
}
   {
#line 435
	unsigned long	Error;
#line 435

#line 435
	Error = pthread_mutex_init(&(bars->sl_phase_7).mutex, NULL);
#line 435
	if (Error != 0) {
#line 435
		printf("Error while initializing barrier.\n");
#line 435
		exit(-1);
#line 435
	}
#line 435

#line 435
	Error = pthread_cond_init(&(bars->sl_phase_7).cv, NULL);
#line 435
	if (Error != 0) {
#line 435
		printf("Error while initializing barrier.\n");
#line 435
		pthread_mutex_destroy(&(bars->sl_phase_7).mutex);
#line 435
		exit(-1);
#line 435
	}
#line 435

#line 435
	(bars->sl_phase_7).counter = 0;
#line 435
	(bars->sl_phase_7).cycle = 0;
#line 435
}
   {
#line 436
	unsigned long	Error;
#line 436

#line 436
	Error = pthread_mutex_init(&(bars->sl_phase_8).mutex, NULL);
#line 436
	if (Error != 0) {
#line 436
		printf("Error while initializing barrier.\n");
#line 436
		exit(-1);
#line 436
	}
#line 436

#line 436
	Error = pthread_cond_init(&(bars->sl_phase_8).cv, NULL);
#line 436
	if (Error != 0) {
#line 436
		printf("Error while initializing barrier.\n");
#line 436
		pthread_mutex_destroy(&(bars->sl_phase_8).mutex);
#line 436
		exit(-1);
#line 436
	}
#line 436

#line 436
	(bars->sl_phase_8).counter = 0;
#line 436
	(bars->sl_phase_8).cycle = 0;
#line 436
}
   {
#line 437
	unsigned long	Error;
#line 437

#line 437
	Error = pthread_mutex_init(&(bars->sl_phase_9).mutex, NULL);
#line 437
	if (Error != 0) {
#line 437
		printf("Error while initializing barrier.\n");
#line 437
		exit(-1);
#line 437
	}
#line 437

#line 437
	Error = pthread_cond_init(&(bars->sl_phase_9).cv, NULL);
#line 437
	if (Error != 0) {
#line 437
		printf("Error while initializing barrier.\n");
#line 437
		pthread_mutex_destroy(&(bars->sl_phase_9).mutex);
#line 437
		exit(-1);
#line 437
	}
#line 437

#line 437
	(bars->sl_phase_9).counter = 0;
#line 437
	(bars->sl_phase_9).cycle = 0;
#line 437
}
   {
#line 438
	unsigned long	Error;
#line 438

#line 438
	Error = pthread_mutex_init(&(bars->sl_phase_10).mutex, NULL);
#line 438
	if (Error != 0) {
#line 438
		printf("Error while initializing barrier.\n");
#line 438
		exit(-1);
#line 438
	}
#line 438

#line 438
	Error = pthread_cond_init(&(bars->sl_phase_10).cv, NULL);
#line 438
	if (Error != 0) {
#line 438
		printf("Error while initializing barrier.\n");
#line 438
		pthread_mutex_destroy(&(bars->sl_phase_10).mutex);
#line 438
		exit(-1);
#line 438
	}
#line 438

#line 438
	(bars->sl_phase_10).counter = 0;
#line 438
	(bars->sl_phase_10).cycle = 0;
#line 438
}
   {
#line 439
	unsigned long	Error;
#line 439

#line 439
	Error = pthread_mutex_init(&(bars->error_barrier).mutex, NULL);
#line 439
	if (Error != 0) {
#line 439
		printf("Error while initializing barrier.\n");
#line 439
		exit(-1);
#line 439
	}
#line 439

#line 439
	Error = pthread_cond_init(&(bars->error_barrier).cv, NULL);
#line 439
	if (Error != 0) {
#line 439
		printf("Error while initializing barrier.\n");
#line 439
		pthread_mutex_destroy(&(bars->error_barrier).mutex);
#line 439
		exit(-1);
#line 439
	}
#line 439

#line 439
	(bars->error_barrier).counter = 0;
#line 439
	(bars->error_barrier).cycle = 0;
#line 439
}
#else
   {
#line 441
	unsigned long	Error;
#line 441

#line 441
	Error = pthread_mutex_init(&(bars->barrier).mutex, NULL);
#line 441
	if (Error != 0) {
#line 441
		printf("Error while initializing barrier.\n");
#line 441
		exit(-1);
#line 441
	}
#line 441

#line 441
	Error = pthread_cond_init(&(bars->barrier).cv, NULL);
#line 441
	if (Error != 0) {
#line 441
		printf("Error while initializing barrier.\n");
#line 441
		pthread_mutex_destroy(&(bars->barrier).mutex);
#line 441
		exit(-1);
#line 441
	}
#line 441

#line 441
	(bars->barrier).counter = 0;
#line 441
	(bars->barrier).cycle = 0;
#line 441
}
#endif

   link_all();

   multi->err_multi = 0.0;
   i_int_coeff[0] = 0.0;
   j_int_coeff[0] = 0.0;
   for (i=0;i<numlev;i++) {
     i_int_coeff[i] = 1.0/(imx[i]-1);
     j_int_coeff[i] = 1.0/(jmx[i]-1);
   }

/* initialize constants and variables

   id is a global shared variable that has fetch-and-add operations
   performed on it by processes to obtain their pids.   */

   global->id = 0;
   global->psibi = 0.0;
   pi = atan(1.0);
   pi = 4.*pi;

   factjacob = -1./(12.*res*res);
   factlap = 1./(res*res);
   eig2 = -h*f0*f0/(h1*h3*gpr);

   jmm1 = jm-1 ;
   ysca = ((double) jmm1)*res ;

   im = (imx[numlev-1]-2)/yprocs + 2;
   jm = (jmx[numlev-1]-2)/xprocs + 2;

   if (do_output) {
     printf("                       MULTIGRID OUTPUTS\n");
   }

   {
#line 478
	long	i, Error;
#line 478

#line 478
	for (i = 0; i < (nprocs) - 1; i++) {
#line 478
		Error = pthread_create(&PThreadTable[i], NULL, (void * (*)(void *))(slave), NULL);
#line 478
		if (Error != 0) {
#line 478
			printf("Error in pthread_create().\n");
#line 478
			exit(-1);
#line 478
		}
#line 478
	}
#line 478

#line 478
	slave();
#line 478
};
   {
#line 479
	unsigned long	i, Error;
#line 479
	for (i = 0; i < (nprocs) - 1; i++) {
#line 479
		Error = pthread_join(PThreadTable[i], NULL);
#line 479
		if (Error != 0) {
#line 479
			printf("Error in pthread_join().\n");
#line 479
			exit(-1);
#line 479
		}
#line 479
	}
#line 479
};
   {
#line 480
	struct timeval	FullTime;
#line 480

#line 480
	gettimeofday(&FullTime, NULL);
#line 480
	(computeend) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 480
}

   printf("\n");
   printf("                       PROCESS STATISTICS\n");
   printf("                  Total          Multigrid         Multigrid\n");
   printf(" Proc             Time             Time            Fraction\n");
   printf("    0   %15.0f    %15.0f        %10.3f\n", gp[0].total_time,gp[0].multi_time, gp[0].multi_time/gp[0].total_time);

   if (do_stats) {
     min_total = max_total = avg_total = gp[0].total_time;
     min_multi = max_multi = avg_multi = gp[0].multi_time;
     min_frac = max_frac = avg_frac = gp[0].multi_time/gp[0].total_time;
     for (i=1;i<nprocs;i++) {
       if (gp[i].total_time > max_total) {
         max_total = gp[i].total_time;
       }
       if (gp[i].total_time < min_total) {
         min_total = gp[i].total_time;
       }
       if (gp[i].multi_time > max_multi) {
         max_multi = gp[i].multi_time;
       }
       if (gp[i].multi_time < min_multi) {
         min_multi = gp[i].multi_time;
       }
       if (gp[i].multi_time/gp[i].total_time > max_frac) {
         max_frac = gp[i].multi_time/gp[i].total_time;
       }
       if (gp[i].multi_time/gp[i].total_time < min_frac) {
         min_frac = gp[i].multi_time/gp[i].total_time;
       }
       avg_total += gp[i].total_time;
       avg_multi += gp[i].multi_time;
       avg_frac += gp[i].multi_time/gp[i].total_time;
     }
     avg_total = avg_total / nprocs;
     avg_multi = avg_multi / nprocs;
     avg_frac = avg_frac / nprocs;
     for (i=1;i<nprocs;i++) {
       printf("  %3ld   %15.0f    %15.0f        %10.3f\n", i,gp[i].total_time,gp[i].multi_time, gp[i].multi_time/gp[i].total_time);
     }
     printf("  Avg   %15.0f    %15.0f        %10.3f\n", avg_total,avg_multi,avg_frac);
     printf("  Min   %15.0f    %15.0f        %10.3f\n", min_total,min_multi,min_frac);
     printf("  Max   %15.0f    %15.0f        %10.3f\n", max_total,max_multi,max_frac);
   }
   printf("\n");

   global->starttime = start;
   printf("                       TIMING INFORMATION\n");
   printf("Start time                        : %16lu\n", global->starttime);
   printf("Initialization finish time        : %16lu\n", global->trackstart);
   printf("Overall finish time               : %16lu\n", computeend);
   printf("Total time with initialization    : %16lu\n", computeend-global->starttime);
   printf("Total time without initialization : %16lu\n", computeend-global->trackstart);
   printf("    (excludes first timestep)\n");
   printf("\n");

   {exit(0);}
}

long log_2(long number)
{
  long cumulative = 1;
  long out = 0;
  long done = 0;

  while ((cumulative < number) && (!done) && (out < 50)) {
    if (cumulative == number) {
      done = 1;
    } else {
      cumulative = cumulative * 2;
      out ++;
    }
  }

  if (cumulative == number) {
    return(out);
  } else {
    return(-1);
  }
}

void printerr(char *s)
{
  fprintf(stderr,"ERROR: %s\n",s);
}

#line 228 "/home/jyy/Dev/splash2/codes/null_macros/c.m4.null.POSIX"

#line 1 "multi.C"
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

/* Shared memory implementation of the multigrid method
   Implementation uses red-black gauss-seidel relaxation
   iterations, w cycles, and the method of half-injection for
   residual computation. */


#line 22
#include <pthread.h>
#line 22
#include <sys/time.h>
#line 22
#include <unistd.h>
#line 22
#include <stdlib.h>
#line 22
extern pthread_t PThreadTable[];
#line 22


#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "decs.h"

/* perform multigrid (w cycles)                                     */
void multig(long my_id)
{
   long iter;
   double wu;
   double errp;
   long m;
   long flag1;
   long flag2;
   long k;
   long my_num;
   double wmax;
   double local_err;
   double red_local_err;
   double black_local_err;
   double g_error;

   flag1 = 0;
   flag2 = 0;
   iter = 0;
   m = numlev-1;
   wmax = maxwork;
   my_num = my_id;
   wu = 0.0;

   k = m;
   g_error = 1.0e30;
   while ((!flag1) && (!flag2)) {
     errp = g_error;
     iter++;
     if (my_num == MASTER) {
       multi->err_multi = 0.0;
     }

/* barrier to make sure all procs have finished intadd or rescal   */
/* before proceeding with relaxation                               */
#if defined(MULTIPLE_BARRIERS)
     {
#line 67
	unsigned long	Error, Cycle;
#line 67
	long		Cancel, Temp;
#line 67

#line 67
	Error = pthread_mutex_lock(&(bars->error_barrier).mutex);
#line 67
	if (Error != 0) {
#line 67
		printf("Error while trying to get lock in barrier.\n");
#line 67
		exit(-1);
#line 67
	}
#line 67

#line 67
	Cycle = (bars->error_barrier).cycle;
#line 67
	if (++(bars->error_barrier).counter != (nprocs)) {
#line 67
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 67
		while (Cycle == (bars->error_barrier).cycle) {
#line 67
			Error = pthread_cond_wait(&(bars->error_barrier).cv, &(bars->error_barrier).mutex);
#line 67
			if (Error != 0) {
#line 67
				break;
#line 67
			}
#line 67
		}
#line 67
		pthread_setcancelstate(Cancel, &Temp);
#line 67
	} else {
#line 67
		(bars->error_barrier).cycle = !(bars->error_barrier).cycle;
#line 67
		(bars->error_barrier).counter = 0;
#line 67
		Error = pthread_cond_broadcast(&(bars->error_barrier).cv);
#line 67
	}
#line 67
	pthread_mutex_unlock(&(bars->error_barrier).mutex);
#line 67
}
#else
     {
#line 69
	unsigned long	Error, Cycle;
#line 69
	long		Cancel, Temp;
#line 69

#line 69
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 69
	if (Error != 0) {
#line 69
		printf("Error while trying to get lock in barrier.\n");
#line 69
		exit(-1);
#line 69
	}
#line 69

#line 69
	Cycle = (bars->barrier).cycle;
#line 69
	if (++(bars->barrier).counter != (nprocs)) {
#line 69
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 69
		while (Cycle == (bars->barrier).cycle) {
#line 69
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 69
			if (Error != 0) {
#line 69
				break;
#line 69
			}
#line 69
		}
#line 69
		pthread_setcancelstate(Cancel, &Temp);
#line 69
	} else {
#line 69
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 69
		(bars->barrier).counter = 0;
#line 69
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 69
	}
#line 69
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 69
}
#endif
     copy_black(k,my_num);

     relax(k,&red_local_err,RED_ITER,my_num);

/* barrier to make sure all red computations have been performed   */
#if defined(MULTIPLE_BARRIERS)
     {
#line 77
	unsigned long	Error, Cycle;
#line 77
	long		Cancel, Temp;
#line 77

#line 77
	Error = pthread_mutex_lock(&(bars->error_barrier).mutex);
#line 77
	if (Error != 0) {
#line 77
		printf("Error while trying to get lock in barrier.\n");
#line 77
		exit(-1);
#line 77
	}
#line 77

#line 77
	Cycle = (bars->error_barrier).cycle;
#line 77
	if (++(bars->error_barrier).counter != (nprocs)) {
#line 77
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 77
		while (Cycle == (bars->error_barrier).cycle) {
#line 77
			Error = pthread_cond_wait(&(bars->error_barrier).cv, &(bars->error_barrier).mutex);
#line 77
			if (Error != 0) {
#line 77
				break;
#line 77
			}
#line 77
		}
#line 77
		pthread_setcancelstate(Cancel, &Temp);
#line 77
	} else {
#line 77
		(bars->error_barrier).cycle = !(bars->error_barrier).cycle;
#line 77
		(bars->error_barrier).counter = 0;
#line 77
		Error = pthread_cond_broadcast(&(bars->error_barrier).cv);
#line 77
	}
#line 77
	pthread_mutex_unlock(&(bars->error_barrier).mutex);
#line 77
}
#else
     {
#line 79
	unsigned long	Error, Cycle;
#line 79
	long		Cancel, Temp;
#line 79

#line 79
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
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
	Cycle = (bars->barrier).cycle;
#line 79
	if (++(bars->barrier).counter != (nprocs)) {
#line 79
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 79
		while (Cycle == (bars->barrier).cycle) {
#line 79
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
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
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 79
		(bars->barrier).counter = 0;
#line 79
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 79
	}
#line 79
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 79
}
#endif
     copy_red(k,my_num);

     relax(k,&black_local_err,BLACK_ITER,my_num);

/* compute max local error from red_local_err and black_local_err  */

     if (red_local_err > black_local_err) {
       local_err = red_local_err;
     } else {
       local_err = black_local_err;
     }

/* update the global error if necessary                         */

     {pthread_mutex_lock(&(locks->error_lock));}
     if (local_err > multi->err_multi) {
       multi->err_multi = local_err;
     }
     {pthread_mutex_unlock(&(locks->error_lock));}

/* a single relaxation sweep at the finest level is one unit of    */
/* work                                                            */

     wu+=pow((double)4.0,(double)k-m);

/* barrier to make sure all processors have checked local error    */
#if defined(MULTIPLE_BARRIERS)
     {
#line 108
	unsigned long	Error, Cycle;
#line 108
	long		Cancel, Temp;
#line 108

#line 108
	Error = pthread_mutex_lock(&(bars->error_barrier).mutex);
#line 108
	if (Error != 0) {
#line 108
		printf("Error while trying to get lock in barrier.\n");
#line 108
		exit(-1);
#line 108
	}
#line 108

#line 108
	Cycle = (bars->error_barrier).cycle;
#line 108
	if (++(bars->error_barrier).counter != (nprocs)) {
#line 108
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 108
		while (Cycle == (bars->error_barrier).cycle) {
#line 108
			Error = pthread_cond_wait(&(bars->error_barrier).cv, &(bars->error_barrier).mutex);
#line 108
			if (Error != 0) {
#line 108
				break;
#line 108
			}
#line 108
		}
#line 108
		pthread_setcancelstate(Cancel, &Temp);
#line 108
	} else {
#line 108
		(bars->error_barrier).cycle = !(bars->error_barrier).cycle;
#line 108
		(bars->error_barrier).counter = 0;
#line 108
		Error = pthread_cond_broadcast(&(bars->error_barrier).cv);
#line 108
	}
#line 108
	pthread_mutex_unlock(&(bars->error_barrier).mutex);
#line 108
}
#else
     {
#line 110
	unsigned long	Error, Cycle;
#line 110
	long		Cancel, Temp;
#line 110

#line 110
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 110
	if (Error != 0) {
#line 110
		printf("Error while trying to get lock in barrier.\n");
#line 110
		exit(-1);
#line 110
	}
#line 110

#line 110
	Cycle = (bars->barrier).cycle;
#line 110
	if (++(bars->barrier).counter != (nprocs)) {
#line 110
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 110
		while (Cycle == (bars->barrier).cycle) {
#line 110
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 110
			if (Error != 0) {
#line 110
				break;
#line 110
			}
#line 110
		}
#line 110
		pthread_setcancelstate(Cancel, &Temp);
#line 110
	} else {
#line 110
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 110
		(bars->barrier).counter = 0;
#line 110
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 110
	}
#line 110
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 110
}
#endif
     g_error = multi->err_multi;

/* barrier to make sure master does not cycle back to top of loop  */
/* and reset global->err before we read it and decide what to do   */
#if defined(MULTIPLE_BARRIERS)
     {
#line 117
	unsigned long	Error, Cycle;
#line 117
	long		Cancel, Temp;
#line 117

#line 117
	Error = pthread_mutex_lock(&(bars->error_barrier).mutex);
#line 117
	if (Error != 0) {
#line 117
		printf("Error while trying to get lock in barrier.\n");
#line 117
		exit(-1);
#line 117
	}
#line 117

#line 117
	Cycle = (bars->error_barrier).cycle;
#line 117
	if (++(bars->error_barrier).counter != (nprocs)) {
#line 117
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 117
		while (Cycle == (bars->error_barrier).cycle) {
#line 117
			Error = pthread_cond_wait(&(bars->error_barrier).cv, &(bars->error_barrier).mutex);
#line 117
			if (Error != 0) {
#line 117
				break;
#line 117
			}
#line 117
		}
#line 117
		pthread_setcancelstate(Cancel, &Temp);
#line 117
	} else {
#line 117
		(bars->error_barrier).cycle = !(bars->error_barrier).cycle;
#line 117
		(bars->error_barrier).counter = 0;
#line 117
		Error = pthread_cond_broadcast(&(bars->error_barrier).cv);
#line 117
	}
#line 117
	pthread_mutex_unlock(&(bars->error_barrier).mutex);
#line 117
}
#else
     {
#line 119
	unsigned long	Error, Cycle;
#line 119
	long		Cancel, Temp;
#line 119

#line 119
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 119
	if (Error != 0) {
#line 119
		printf("Error while trying to get lock in barrier.\n");
#line 119
		exit(-1);
#line 119
	}
#line 119

#line 119
	Cycle = (bars->barrier).cycle;
#line 119
	if (++(bars->barrier).counter != (nprocs)) {
#line 119
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 119
		while (Cycle == (bars->barrier).cycle) {
#line 119
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 119
			if (Error != 0) {
#line 119
				break;
#line 119
			}
#line 119
		}
#line 119
		pthread_setcancelstate(Cancel, &Temp);
#line 119
	} else {
#line 119
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 119
		(bars->barrier).counter = 0;
#line 119
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 119
	}
#line 119
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 119
}
#endif

     if (g_error >= lev_tol[k]) {
       if (wu > wmax) {
/* max work exceeded                                               */
         flag1 = 1;
	 fprintf(stderr,"ERROR: Maximum work limit %0.5f exceeded\n",wmax);
	 exit(-1);
       } else {
/* if we have not converged                                        */
         if ((k != 0) && (g_error/errp >= 0.6) &&
	     (k > minlevel)) {
/* if need to go to coarser grid                                   */

           copy_borders(k,my_num);
           copy_rhs_borders(k,my_num);

/* This bar is needed because the routine rescal uses the neighbor's
   border points to compute s4.  We must ensure that the neighbor's
   border points have been written before we try computing the new
   rescal values                                                   */

#if defined(MULTIPLE_BARRIERS)
	   {
#line 143
	unsigned long	Error, Cycle;
#line 143
	long		Cancel, Temp;
#line 143

#line 143
	Error = pthread_mutex_lock(&(bars->error_barrier).mutex);
#line 143
	if (Error != 0) {
#line 143
		printf("Error while trying to get lock in barrier.\n");
#line 143
		exit(-1);
#line 143
	}
#line 143

#line 143
	Cycle = (bars->error_barrier).cycle;
#line 143
	if (++(bars->error_barrier).counter != (nprocs)) {
#line 143
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 143
		while (Cycle == (bars->error_barrier).cycle) {
#line 143
			Error = pthread_cond_wait(&(bars->error_barrier).cv, &(bars->error_barrier).mutex);
#line 143
			if (Error != 0) {
#line 143
				break;
#line 143
			}
#line 143
		}
#line 143
		pthread_setcancelstate(Cancel, &Temp);
#line 143
	} else {
#line 143
		(bars->error_barrier).cycle = !(bars->error_barrier).cycle;
#line 143
		(bars->error_barrier).counter = 0;
#line 143
		Error = pthread_cond_broadcast(&(bars->error_barrier).cv);
#line 143
	}
#line 143
	pthread_mutex_unlock(&(bars->error_barrier).mutex);
#line 143
}
#else
	   {
#line 145
	unsigned long	Error, Cycle;
#line 145
	long		Cancel, Temp;
#line 145

#line 145
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 145
	if (Error != 0) {
#line 145
		printf("Error while trying to get lock in barrier.\n");
#line 145
		exit(-1);
#line 145
	}
#line 145

#line 145
	Cycle = (bars->barrier).cycle;
#line 145
	if (++(bars->barrier).counter != (nprocs)) {
#line 145
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 145
		while (Cycle == (bars->barrier).cycle) {
#line 145
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 145
			if (Error != 0) {
#line 145
				break;
#line 145
			}
#line 145
		}
#line 145
		pthread_setcancelstate(Cancel, &Temp);
#line 145
	} else {
#line 145
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 145
		(bars->barrier).counter = 0;
#line 145
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 145
	}
#line 145
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 145
}
#endif

           rescal(k,my_num);

/* transfer residual to rhs of coarser grid                        */
           lev_tol[k-1] = 0.3 * g_error;
           k = k-1;
           putz(k,my_num);
/* make initial guess on coarser grid zero                         */
           g_error = 1.0e30;
         }
       }
     } else {
/* if we have converged at this level                              */
       if (k == m) {
/* if finest grid, we are done                                     */
         flag2 = 1;
       } else {
/* else go to next finest grid                                     */

         copy_borders(k,my_num);

         intadd(k,my_num);
/* changes the grid values at the finer level.  rhs at finer level */
/* remains what it already is                                      */
         k++;
         g_error = 1.0e30;
       }
     }
   }
   if (do_output) {
     if (my_num == MASTER) {
       printf("iter %ld, level %ld, residual norm %12.8e, work = %7.3f\n", iter,k,multi->err_multi,wu);
     }
   }
}

/* perform red or black iteration (not both)                    */
void relax(long k, double *err, long color, long my_num)
{
   long i;
   long j;
   long iend;
   long jend;
   long oddistart;
   long oddjstart;
   long evenistart;
   long evenjstart;
   double a;
   double h;
   double factor;
   double maxerr;
   double newerr;
   double oldval;
   double newval;
   double **t2a;
   double **t2b;
   double *t1a;
   double *t1b;
   double *t1c;
   double *t1d;

   i = 0;
   j = 0;

   *err = 0.0;
   h = lev_res[k];

/* points whose sum of row and col index is even do a red iteration, */
/* others do a black				                     */

   evenistart = gp[my_num].eist[k];
   evenjstart = gp[my_num].ejst[k];
   oddistart = gp[my_num].oist[k];
   oddjstart = gp[my_num].ojst[k];

   iend = gp[my_num].rlien[k];
   jend = gp[my_num].rljen[k];

   factor = 4.0 - eig2 * h * h ;
   maxerr = 0.0;
   t2a = (double **) q_multi[my_num][k];
   t2b = (double **) rhs_multi[my_num][k];
   if (color == RED_ITER) {
     for (i=evenistart;i<iend;i+=2) {
       t1a = (double *) t2a[i];
       t1b = (double *) t2b[i];
       t1c = (double *) t2a[i-1];
       t1d = (double *) t2a[i+1];
       for (j=evenjstart;j<jend;j+=2) {
         a = t1a[j+1] + t1a[j-1] +
	     t1c[j] + t1d[j] -
	     t1b[j] ;
         oldval = t1a[j];
         newval = a / factor;
         newerr = oldval - newval;
         t1a[j] = newval;
         if (fabs(newerr) > maxerr) {
           maxerr = fabs(newerr);
         }
       }
     }
     for (i=oddistart;i<iend;i+=2) {
       t1a = (double *) t2a[i];
       t1b = (double *) t2b[i];
       t1c = (double *) t2a[i-1];
       t1d = (double *) t2a[i+1];
       for (j=oddjstart;j<jend;j+=2) {
         a = t1a[j+1] + t1a[j-1] +
	     t1c[j] + t1d[j] -
	     t1b[j] ;
         oldval = t1a[j];
         newval = a / factor;
         newerr = oldval - newval;
         t1a[j] = newval;
         if (fabs(newerr) > maxerr) {
           maxerr = fabs(newerr);
         }
       }
     }
   } else if (color == BLACK_ITER) {
     for (i=evenistart;i<iend;i+=2) {
       t1a = (double *) t2a[i];
       t1b = (double *) t2b[i];
       t1c = (double *) t2a[i-1];
       t1d = (double *) t2a[i+1];
       for (j=oddjstart;j<jend;j+=2) {
         a = t1a[j+1] + t1a[j-1] +
	     t1c[j] + t1d[j] -
	     t1b[j] ;
         oldval = t1a[j];
         newval = a / factor;
         newerr = oldval - newval;
         t1a[j] = newval;
         if (fabs(newerr) > maxerr) {
           maxerr = fabs(newerr);
         }
       }
     }
     for (i=oddistart;i<iend;i+=2) {
       t1a = (double *) t2a[i];
       t1b = (double *) t2b[i];
       t1c = (double *) t2a[i-1];
       t1d = (double *) t2a[i+1];
       for (j=evenjstart;j<jend;j+=2) {
         a = t1a[j+1] + t1a[j-1] +
	     t1c[j] + t1d[j] -
	     t1b[j] ;
         oldval = t1a[j];
         newval = a / factor;
         newerr = oldval - newval;
         t1a[j] = newval;
         if (fabs(newerr) > maxerr) {
           maxerr = fabs(newerr);
         }
       }
     }
   }
   *err = maxerr;
}

/* perform half-injection to next coarsest level                */
void rescal(long kf, long my_num)
{
   long ic;
   long if17;
   long jf;
   long jc;
   long krc;
   long istart;
   long iend;
   long jstart;
   long jend;
   double hf;
   double hc;
   double s;
   double s1;
   double s2;
   double s3;
   double s4;
   double factor;
   double int1;
   double int2;
   double i_int_factor;
   double j_int_factor;
   double int_val;
   long i_off;
   long j_off;
   long up_proc;
   long left_proc;
   long im;
   long jm;
   double temp;
   double temp2;
   double **t2a;
   double **t2b;
   double **t2c;
   double *t1a;
   double *t1b;
   double *t1c;
   double *t1d;
   double *t1e;
   double *t1f;
   double *t1g;
   double *t1h;

   krc = kf - 1;
   hc = lev_res[krc];
   hf = lev_res[kf];
   i_off = gp[my_num].rownum*ypts_per_proc[krc];
   j_off = gp[my_num].colnum*xpts_per_proc[krc];
   up_proc = gp[my_num].neighbors[UP];
   left_proc = gp[my_num].neighbors[LEFT];
   im = (imx[kf]-2)/yprocs;
   jm = (jmx[kf]-2)/xprocs;

   istart = gp[my_num].rlist[krc];
   jstart = gp[my_num].rljst[krc];
   iend = gp[my_num].rlien[krc] - 1;
   jend = gp[my_num].rljen[krc] - 1;

   factor = 4.0 - eig2 * hf * hf;

   t2a = (double **) q_multi[my_num][kf];
   t2b = (double **) rhs_multi[my_num][kf];
   t2c = (double **) rhs_multi[my_num][krc];
   if17=2*(istart-1);
   for(ic=istart;ic<=iend;ic++) {
     if17+=2;
     i_int_factor = (ic+i_off) * i_int_coeff[krc] * 0.5;
     jf = 2 * (jstart - 1);
     t1a = (double *) t2a[if17];
     t1b = (double *) t2b[if17];
     t1c = (double *) t2c[ic];
     t1d = (double *) t2a[if17-1];
     t1e = (double *) t2a[if17+1];
     t1f = (double *) t2a[if17-2];
     t1g = (double *) t2a[if17-3];
     t1h = (double *) t2b[if17-2];
     for(jc=jstart;jc<=jend;jc++) {
       jf+=2;
       j_int_factor = (jc+j_off)*j_int_coeff[krc] * 0.5;

/*             method of half-injection uses 2.0 instead of 4.0 */

/* do bilinear interpolation */
       s = t1a[jf+1] + t1a[jf-1] + t1d[jf] + t1e[jf];
       s1 = 2.0 * (t1b[jf] - s + factor * t1a[jf]);
       if (((if17 == 2) && (gp[my_num].neighbors[UP] == -1)) ||
	   ((jf == 2) && (gp[my_num].neighbors[LEFT] == -1))) {
          s2 = 0;
          s3 = 0;
          s4 = 0;
       } else if ((if17 == 2) || (jf == 2)) {
	  if (jf == 2) {
	    temp = q_multi[left_proc][kf][if17][jm-1];
          } else {
            temp = t1a[jf-3];
          }
          s = t1a[jf-1] + temp + t1d[jf-2] + t1e[jf-2];
          s2 = 2.0 * (t1b[jf-2] - s + factor * t1a[jf-2]);
	  if (if17 == 2) {
	    temp = q_multi[up_proc][kf][im-1][jf];
          } else {
            temp = t1g[jf];
          }
          s = t1f[jf+1]+ t1f[jf-1]+ temp + t1d[jf];
          s3 = 2.0 * (t1h[jf] - s + factor * t1f[jf]);
	  if (jf == 2) {
	    temp = q_multi[left_proc][kf][if17-2][jm-1];
          } else {
            temp = t1f[jf-3];
          }
	  if (if17 == 2) {
	    temp2 = q_multi[up_proc][kf][im-1][jf-2];
          } else {
            temp2 = t1g[jf-2];
          }
          s = t1f[jf-1]+ temp + temp2 + t1d[jf-2];
          s4 = 2.0 * (t1h[jf-2] - s + factor * t1f[jf-2]);
       } else {
          s = t1a[jf-1] + t1a[jf-3] + t1d[jf-2] + t1e[jf-2];
          s2 = 2.0 * (t1b[jf-2] - s + factor * t1a[jf-2]);
          s = t1f[jf+1]+ t1f[jf-1]+ t1g[jf] +   t1d[jf];
          s3 = 2.0 * (t1h[jf] - s + factor * t1f[jf]);
          s = t1f[jf-1]+ t1f[jf-3]+ t1g[jf-2]+ t1d[jf-2];
          s4 = 2.0 * (t1h[jf-2] - s + factor * t1f[jf-2]);
       }
       int1 = j_int_factor*s4 + (1.0-j_int_factor)*s3;
       int2 = j_int_factor*s2 + (1.0-j_int_factor)*s1;
       int_val = i_int_factor*int1+(1.0-i_int_factor)*int2;
       t1c[jc] = i_int_factor*int1+(1.0-i_int_factor)*int2;
     }
   }
}

/* perform interpolation and addition to next finest grid       */
void intadd(long kc, long my_num)
{
   long ic;
   long if17;
   long jf;
   long jc;
   long kf;
   long istart;
   long jstart;
   long iend;
   long jend;
   double hc;
   double hf;
   double int1;
   double int2;
   double i_int_factor1;
   double j_int_factor1;
   double i_int_factor2;
   double j_int_factor2;
   long i_off;
   long j_off;
   double **t2a;
   double **t2b;
   double *t1a;
   double *t1b;
   double *t1c;
   double *t1d;
   double *t1e;

   kf = kc + 1;
   hc = lev_res[kc];
   hf = lev_res[kf];

   istart = gp[my_num].rlist[kc];
   jstart = gp[my_num].rljst[kc];
   iend = gp[my_num].rlien[kc] - 1;
   jend = gp[my_num].rljen[kc] - 1;
   i_off = gp[my_num].rownum*ypts_per_proc[kc];
   j_off = gp[my_num].colnum*xpts_per_proc[kc];

   t2a = (double **) q_multi[my_num][kc];
   t2b = (double **) q_multi[my_num][kf];
   if17 = 2*(istart-1);
   for(ic=istart;ic<=iend;ic++) {
     if17+=2;
     i_int_factor1= ((imx[kc]-2)-(ic+i_off-1)) * (i_int_coeff[kf]);
     i_int_factor2= (ic+i_off) * i_int_coeff[kf];
     jf = 2*(jstart-1);

     t1a = (double *) t2a[ic];
     t1b = (double *) t2a[ic-1];
     t1c = (double *) t2a[ic+1];
     t1d = (double *) t2b[if17];
     t1e = (double *) t2b[if17-1];
     for(jc=jstart;jc<=jend;jc++) {
       jf+=2;
       j_int_factor1= ((jmx[kc]-2)-(jc+j_off-1)) * (j_int_coeff[kf]);
       j_int_factor2= (jc+j_off) * j_int_coeff[kf];

       int1 = j_int_factor1*t1a[jc-1] + (1.0-j_int_factor1)*t1a[jc];
       int2 = j_int_factor1*t1b[jc-1] + (1.0-j_int_factor1)*t1b[jc];
       t1e[jf-1] += i_int_factor1*int2 + (1.0-i_int_factor1)*int1;
       int2 = j_int_factor1*t1c[jc-1] + (1.0-j_int_factor1)*t1c[jc];
       t1d[jf-1] += i_int_factor2*int2 + (1.0-i_int_factor2)*int1;
       int1 = j_int_factor2*t1a[jc+1] + (1.0-j_int_factor2)*t1a[jc];
       int2 = j_int_factor2*t1b[jc+1] + (1.0-j_int_factor2)*t1b[jc];
       t1e[jf] += i_int_factor1*int2 + (1.0-i_int_factor1)*int1;
       int2 = j_int_factor2*t1c[jc+1] + (1.0-j_int_factor2)*t1c[jc];
       t1d[jf] += i_int_factor2*int2 + (1.0-i_int_factor2)*int1;
     }
   }
}

/* initialize a grid to zero in parallel                        */
void putz(long k, long my_num)
{
   long i;
   long j;
   long istart;
   long jstart;
   long iend;
   long jend;
   double **t2a;
   double *t1a;

   istart = gp[my_num].rlist[k];
   jstart = gp[my_num].rljst[k];
   iend = gp[my_num].rlien[k];
   jend = gp[my_num].rljen[k];

   t2a = (double **) q_multi[my_num][k];
   for (i=istart;i<=iend;i++) {
     t1a = (double *) t2a[i];
     for (j=jstart;j<=jend;j++) {
       t1a[j] = 0.0;
     }
   }
}

void copy_borders(long k, long pid)
{
  long i;
  long j;
  long jj;
  long im;
  long jm;
  long lastrow;
  long lastcol;
  double **t2a;
  double **t2b;
  double *t1a;
  double *t1b;

   im = (imx[k]-2)/yprocs + 2;
   jm = (jmx[k]-2)/xprocs + 2;
   lastrow = (imx[k]-2)/yprocs;
   lastcol = (jmx[k]-2)/xprocs;

   t2a = (double **) q_multi[pid][k];
   jj = gp[pid].neighbors[UPLEFT];
   if (jj != -1) {
     t2a[0][0]=q_multi[jj][k][im-2][jm-2];
   }
   jj = gp[pid].neighbors[UPRIGHT];
   if (jj != -1) {
     t2a[0][jm-1]=q_multi[jj][k][im-2][1];
   }
   jj = gp[pid].neighbors[DOWNLEFT];
   if (jj != -1) {
     t2a[im-1][0]=q_multi[jj][k][1][jm-2];
   }
   jj = gp[pid].neighbors[DOWNRIGHT];
   if (jj != -1) {
     t2a[im-1][jm-1]=q_multi[jj][k][1][1];
   }

   if (gp[pid].neighbors[UP] == -1) {
     jj = gp[pid].neighbors[LEFT];
     if (jj != -1) {
       t2a[0][0] = q_multi[jj][k][0][jm-2];
     } else {
       jj = gp[pid].neighbors[DOWN];
       if (jj != -1) {
         t2a[im-1][0] = q_multi[jj][k][1][0];
       }
     }
     jj = gp[pid].neighbors[RIGHT];
     if (jj != -1) {
       t2a[0][jm-1] = q_multi[jj][k][0][1];
     } else {
       jj = gp[pid].neighbors[DOWN];
       if (jj != -1) {
         t2a[im-1][jm-1] = q_multi[jj][k][1][jm-1];
       }
     }
   } else if (gp[pid].neighbors[DOWN] == -1) {
     jj = gp[pid].neighbors[LEFT];
     if (jj != -1) {
       t2a[im-1][0] = q_multi[jj][k][im-1][jm-2];
     } else {
       jj = gp[pid].neighbors[UP];
       if (jj != -1) {
         t2a[0][0] = q_multi[jj][k][im-2][0];
       }
     }
     jj = gp[pid].neighbors[RIGHT];
     if (jj != -1) {
       t2a[im-1][jm-1] = q_multi[jj][k][im-1][1];
     } else {
       jj = gp[pid].neighbors[UP];
       if (jj != -1) {
         t2a[0][jm-1] = q_multi[jj][k][im-2][jm-1];
       }
     }
   } else if (gp[pid].neighbors[LEFT] == -1) {
     jj = gp[pid].neighbors[UP];
     if (jj != -1) {
       t2a[0][0] = q_multi[jj][k][im-2][0];
     }
     jj = gp[pid].neighbors[DOWN];
     if (jj != -1) {
       t2a[im-1][0] = q_multi[jj][k][1][0];
     }
   } else if (gp[pid].neighbors[RIGHT] == -1) {
     jj = gp[pid].neighbors[UP];
     if (jj != -1) {
       t2a[0][jm-1] = q_multi[jj][k][im-2][jm-1];
     }
     jj = gp[pid].neighbors[DOWN];
     if (jj != -1) {
       t2a[im-1][jm-1] = q_multi[jj][k][1][jm-1];
     }
   }

   j = gp[pid].neighbors[UP];
   if (j != -1) {
     t1a = (double *) t2a[0];
     t1b = (double *) q_multi[j][k][im-2];
     for (i=1;i<=lastcol;i++) {
       t1a[i] = t1b[i];
     }
   }
   j = gp[pid].neighbors[DOWN];
   if (j != -1) {
     t1a = (double *) t2a[im-1];
     t1b = (double *) q_multi[j][k][1];
     for (i=1;i<=lastcol;i++) {
       t1a[i] = t1b[i];
     }
   }
   j = gp[pid].neighbors[LEFT];
   if (j != -1) {
     t2b = (double **) q_multi[j][k];
     for (i=1;i<=lastrow;i++) {
       t2a[i][0] = t2b[i][jm-2];
     }
   }
   j = gp[pid].neighbors[RIGHT];
   if (j != -1) {
     t2b = (double **) q_multi[j][k];
     for (i=1;i<=lastrow;i++) {
       t2a[i][jm-1] = t2b[i][1];
     }
   }

}

void copy_rhs_borders(long k, long procid)
{
   long i;
   long j;
   long im;
   long jm;
   long lastrow;
   long lastcol;
   double **t2a;
   double **t2b;
   double *t1a;
   double *t1b;

   im = (imx[k]-2)/yprocs+2;
   jm = (jmx[k]-2)/xprocs+2;
   lastrow = (imx[k]-2)/yprocs;
   lastcol = (jmx[k]-2)/xprocs;

   t2a = (double **) rhs_multi[procid][k];
   if (gp[procid].neighbors[UPLEFT] != -1) {
     j = gp[procid].neighbors[UPLEFT];
     t2a[0][0] = rhs_multi[j][k][im-2][jm-2];
   }

   if (gp[procid].neighbors[UP] != -1) {
     j = gp[procid].neighbors[UP];
     if (j != -1) {
       t1a = (double *) t2a[0];
       t1b = (double *) rhs_multi[j][k][im-2];
       for (i=2;i<=lastcol;i+=2) {
         t1a[i] = t1b[i];
       }
     }
   }
   if (gp[procid].neighbors[LEFT] != -1) {
     j = gp[procid].neighbors[LEFT];
     if (j != -1) {
       t2b = (double **) rhs_multi[j][k];
       for (i=2;i<=lastrow;i+=2) {
         t2a[i][0] = t2b[i][jm-2];
       }
     }
   }
}

void copy_red(long k, long procid)
{
   long i;
   long j;
   long im;
   long jm;
   long lastrow;
   long lastcol;
   double **t2a;
   double **t2b;
   double *t1a;
   double *t1b;

   im = (imx[k]-2)/yprocs+2;
   jm = (jmx[k]-2)/xprocs+2;
   lastrow = (imx[k]-2)/yprocs;
   lastcol = (jmx[k]-2)/xprocs;

   t2a = (double **) q_multi[procid][k];
   j = gp[procid].neighbors[UP];
   if (j != -1) {
     t1a = (double *) t2a[0];
     t1b = (double *) q_multi[j][k][im-2];
     for (i=2;i<=lastcol;i+=2) {
       t1a[i] = t1b[i];
     }
   }
   j = gp[procid].neighbors[DOWN];
   if (j != -1) {
     t1a = (double *) t2a[im-1];
     t1b = (double *) q_multi[j][k][1];
     for (i=1;i<=lastcol;i+=2) {
       t1a[i] = t1b[i];
     }
   }
   j = gp[procid].neighbors[LEFT];
   if (j != -1) {
     t2b = (double **) q_multi[j][k];
     for (i=2;i<=lastrow;i+=2) {
       t2a[i][0] = t2b[i][jm-2];
     }
   }
   j = gp[procid].neighbors[RIGHT];
   if (j != -1) {
     t2b = (double **) q_multi[j][k];
     for (i=1;i<=lastrow;i+=2) {
       t2a[i][jm-1] = t2b[i][1];
     }
   }
}

void copy_black(long k, long procid)
{
   long i;
   long j;
   long im;
   long jm;
   long lastrow;
   long lastcol;
   double **t2a;
   double **t2b;
   double *t1a;
   double *t1b;

   im = (imx[k]-2)/yprocs+2;
   jm = (jmx[k]-2)/xprocs+2;
   lastrow = (imx[k]-2)/yprocs;
   lastcol = (jmx[k]-2)/xprocs;

   t2a = (double **) q_multi[procid][k];
   j = gp[procid].neighbors[UP];
   if (j != -1) {
     t1a = (double *) t2a[0];
     t1b = (double *) q_multi[j][k][im-2];
     for (i=1;i<=lastcol;i+=2) {
       t1a[i] = t1b[i];
     }
   }
   j = gp[procid].neighbors[DOWN];
   if (j != -1) {
     t1a = (double *) t2a[im-1];
     t1b = (double *) q_multi[j][k][1];
     for (i=2;i<=lastcol;i+=2) {
       t1a[i] = t1b[i];
     }
   }
   j = gp[procid].neighbors[LEFT];
   if (j != -1) {
     t2b = (double **) q_multi[j][k];
     for (i=1;i<=lastrow;i+=2) {
       t2a[i][0] = t2b[i][jm-2];
     }
   }
   j = gp[procid].neighbors[RIGHT];
   if (j != -1) {
     t2b = (double **) q_multi[j][k];
     for (i=2;i<=lastrow;i+=2) {
       t2a[i][jm-1] = t2b[i][1];
     }
   }
}

#line 228 "/home/jyy/Dev/splash2/codes/null_macros/c.m4.null.POSIX"

#line 1 "slave1.C"
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

/*    ****************
      subroutine slave
      ****************  */


#line 21
#include <pthread.h>
#line 21
#include <sys/time.h>
#line 21
#include <unistd.h>
#line 21
#include <stdlib.h>
#line 21
extern pthread_t PThreadTable[];
#line 21


#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "decs.h"

void slave()
{
   long i;
   long j;
   long nstep;
   long iindex;
   long iday;
   double ysca1;
   double y;
   double factor;
   double sintemp;
   double curlt;
   double ressqr;
   long istart;
   long iend;
   long jstart;
   long jend;
   long ist;
   long ien;
   long jst;
   long jen;
   double fac;
   long dayflag=0;
   long dhourflag=0;
   long endflag=0;
   long firstrow;
   long lastrow;
   long numrows;
   long firstcol;
   long lastcol;
   long numcols;
   long psiindex;
   double psibipriv;
   double ttime;
   double dhour;
   double day;
   long procid;
   long psinum;
   long j_off = 0;
   unsigned long t1;
   double **t2a;
   double **t2b;
   double *t1a;
   double *t1b;
   double *t1c;
   double *t1d;

   ressqr = lev_res[numlev-1] * lev_res[numlev-1];

   {pthread_mutex_lock(&(locks->idlock));}
     procid = global->id;
     global->id = global->id+1;
   {pthread_mutex_unlock(&(locks->idlock));}

#if defined(MULTIPLE_BARRIERS)
   {
#line 84
	unsigned long	Error, Cycle;
#line 84
	long		Cancel, Temp;
#line 84

#line 84
	Error = pthread_mutex_lock(&(bars->sl_prini).mutex);
#line 84
	if (Error != 0) {
#line 84
		printf("Error while trying to get lock in barrier.\n");
#line 84
		exit(-1);
#line 84
	}
#line 84

#line 84
	Cycle = (bars->sl_prini).cycle;
#line 84
	if (++(bars->sl_prini).counter != (nprocs)) {
#line 84
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 84
		while (Cycle == (bars->sl_prini).cycle) {
#line 84
			Error = pthread_cond_wait(&(bars->sl_prini).cv, &(bars->sl_prini).mutex);
#line 84
			if (Error != 0) {
#line 84
				break;
#line 84
			}
#line 84
		}
#line 84
		pthread_setcancelstate(Cancel, &Temp);
#line 84
	} else {
#line 84
		(bars->sl_prini).cycle = !(bars->sl_prini).cycle;
#line 84
		(bars->sl_prini).counter = 0;
#line 84
		Error = pthread_cond_broadcast(&(bars->sl_prini).cv);
#line 84
	}
#line 84
	pthread_mutex_unlock(&(bars->sl_prini).mutex);
#line 84
}
#else
   {
#line 86
	unsigned long	Error, Cycle;
#line 86
	long		Cancel, Temp;
#line 86

#line 86
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 86
	if (Error != 0) {
#line 86
		printf("Error while trying to get lock in barrier.\n");
#line 86
		exit(-1);
#line 86
	}
#line 86

#line 86
	Cycle = (bars->barrier).cycle;
#line 86
	if (++(bars->barrier).counter != (nprocs)) {
#line 86
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 86
		while (Cycle == (bars->barrier).cycle) {
#line 86
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 86
			if (Error != 0) {
#line 86
				break;
#line 86
			}
#line 86
		}
#line 86
		pthread_setcancelstate(Cancel, &Temp);
#line 86
	} else {
#line 86
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 86
		(bars->barrier).counter = 0;
#line 86
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 86
	}
#line 86
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 86
}
#endif
/* POSSIBLE ENHANCEMENT:  Here is where one might pin processes to
   processors to avoid migration. */

/* POSSIBLE ENHANCEMENT:  Here is where one might distribute
   data structures across physically distributed memories as
   desired.

   One way to do this is as follows.  The function allocate(START,SIZE,I)
   is assumed to place all addresses x such that
   (START <= x < START+SIZE) on node I.

   long d_size;
   unsigned long g_size;
   unsigned long mg_size;

   if (procid == MASTER) {
     g_size = ((jmx[numlev-1]-2)/xprocs+2)*((imx[numlev-1]-2)/yprocs+2)*siz
eof(double) +
              ((imx[numlev-1]-2)/yprocs+2)*sizeof(double *);

     mg_size = numlev*sizeof(double **);
     for (i=0;i<numlev;i++) {
       mg_size+=((imx[i]-2)/yprocs+2)*((jmx[i]-2)/xprocs+2)*sizeof(double)+
                ((imx[i]-2)/yprocs+2)*sizeof(double *);
     }
     for (i= 0;i<nprocs;i++) {
       d_size = 2*sizeof(double **);
       allocate((unsigned long) psi[i],d_size,i);
       allocate((unsigned long) psim[i],d_size,i);
       allocate((unsigned long) work1[i],d_size,i);
       allocate((unsigned long) work4[i],d_size,i);
       allocate((unsigned long) work5[i],d_size,i);
       allocate((unsigned long) work7[i],d_size,i);
       allocate((unsigned long) temparray[i],d_size,i);
       allocate((unsigned long) psi[i][0],g_size,i);
       allocate((unsigned long) psi[i][1],g_size,i);
       allocate((unsigned long) psim[i][0],g_size,i);
       allocate((unsigned long) psim[i][1],g_size,i);
       allocate((unsigned long) psium[i],g_size,i);
       allocate((unsigned long) psilm[i],g_size,i);
       allocate((unsigned long) psib[i],g_size,i);
       allocate((unsigned long) ga[i],g_size,i);
       allocate((unsigned long) gb[i],g_size,i);
       allocate((unsigned long) work1[i][0],g_size,i);
       allocate((unsigned long) work1[i][1],g_size,i);
       allocate((unsigned long) work2[i],g_size,i);
       allocate((unsigned long) work3[i],g_size,i);
       allocate((unsigned long) work4[i][0],g_size,i);
       allocate((unsigned long) work4[i][1],g_size,i);
       allocate((unsigned long) work5[i][0],g_size,i);
       allocate((unsigned long) work5[i][1],g_size,i);
       allocate((unsigned long) work6[i],g_size,i);
       allocate((unsigned long) work7[i][0],g_size,i);
       allocate((unsigned long) work7[i][1],g_size,i);
       allocate((unsigned long) temparray[i][0],g_size,i);
       allocate((unsigned long) temparray[i][1],g_size,i);
       allocate((unsigned long) tauz[i],g_size,i);
       allocate((unsigned long) oldga[i],g_size,i);
       allocate((unsigned long) oldgb[i],g_size,i);
       d_size = numlev * sizeof(long);
       allocate((unsigned long) gp[i].rel_num_x,d_size,i);
       allocate((unsigned long) gp[i].rel_num_y,d_size,i);
       allocate((unsigned long) gp[i].eist,d_size,i);
       allocate((unsigned long) gp[i].ejst,d_size,i);
       allocate((unsigned long) gp[i].oist,d_size,i);
       allocate((unsigned long) gp[i].ojst,d_size,i);
       allocate((unsigned long) gp[i].rlist,d_size,i);
       allocate((unsigned long) gp[i].rljst,d_size,i);
       allocate((unsigned long) gp[i].rlien,d_size,i);
       allocate((unsigned long) gp[i].rljen,d_size,i);

       allocate((unsigned long) q_multi[i],mg_size,i);
       allocate((unsigned long) rhs_multi[i],mg_size,i);
       allocate((unsigned long) &(gp[i]),sizeof(struct Global_Private),i);
     }
   }

*/

   t2a = (double **) oldga[procid];
   t2b = (double **) oldgb[procid];
   for (i=0;i<im;i++) {
     t1a = (double *) t2a[i];
     t1b = (double *) t2b[i];
     for (j=0;j<jm;j++) {
        t1a[j] = 0.0;
        t1b[j] = 0.0;
     }
   }

   firstcol = 1;
   lastcol = firstcol + gp[procid].rel_num_x[numlev-1] - 1;
   firstrow = 1;
   lastrow = firstrow + gp[procid].rel_num_y[numlev-1] - 1;
   numcols = gp[procid].rel_num_x[numlev-1];
   numrows = gp[procid].rel_num_y[numlev-1];
   j_off = gp[procid].colnum*numcols;

   if (procid > nprocs/2) {
      psinum = 2;
   } else {
      psinum = 1;
   }

/* every process gets its own copy of the timing variables to avoid
   contention at shared memory locations.  here, these variables
   are initialized.  */

   ttime = 0.0;
   dhour = 0.0;
   nstep = 0 ;
   day = 0.0;

   ysca1 = 0.5*ysca;
   if (procid == MASTER) {
     t1a = (double *) f;
     for (iindex = 0;iindex<=jmx[numlev-1]-1;iindex++) {
       y = ((double) iindex)*res;
       t1a[iindex] = f0+beta*(y-ysca1);
     }
   }

   t2a = (double **) psium[procid];
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[0][0]=0.0;
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[im-1][0]=0.0;
   }
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[0][jm-1]=0.0;
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[im-1][jm-1]=0.0;
   }
   if (gp[procid].neighbors[UP] == -1) {
     t1a = (double *) t2a[0];
     for(j=firstcol;j<=lastcol;j++) {
       t1a[j] = 0.0;
     }
   }
   if (gp[procid].neighbors[DOWN] == -1) {
     t1a = (double *) t2a[im-1];
     for(j=firstcol;j<=lastcol;j++) {
       t1a[j] = 0.0;
     }
   }
   if (gp[procid].neighbors[LEFT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][0] = 0.0;
     }
   }
   if (gp[procid].neighbors[RIGHT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][jm-1] = 0.0;
     }
   }

   for(i=firstrow;i<=lastrow;i++) {
     t1a = (double *) t2a[i];
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
       t1a[iindex] = 0.0;
     }
   }
   t2a = (double **) psilm[procid];
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[0][0]=0.0;
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[im-1][0]=0.0;
   }
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[0][jm-1]=0.0;
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[im-1][jm-1]=0.0;
   }
   if (gp[procid].neighbors[UP] == -1) {
     t1a = (double *) t2a[0];
     for(j=firstcol;j<=lastcol;j++) {
       t1a[j] = 0.0;
     }
   }
   if (gp[procid].neighbors[DOWN] == -1) {
     t1a = (double *) t2a[im-1];
     for(j=firstcol;j<=lastcol;j++) {
       t1a[j] = 0.0;
     }
   }
   if (gp[procid].neighbors[LEFT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][0] = 0.0;
     }
   }
   if (gp[procid].neighbors[RIGHT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][jm-1] = 0.0;
     }
   }
   for(i=firstrow;i<=lastrow;i++) {
     t1a = (double *) t2a[i];
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
       t1a[iindex] = 0.0;
     }
   }

   t2a = (double **) psib[procid];
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[0][0]=1.0;
   }
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[0][jm-1]=1.0;
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[im-1][0]=1.0;
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[im-1][jm-1]=1.0;
   }
   if (gp[procid].neighbors[UP] == -1) {
     t1a = (double *) t2a[0];
     for(j=firstcol;j<=lastcol;j++) {
       t1a[j] = 1.0;
     }
   }
   if (gp[procid].neighbors[DOWN] == -1) {
     t1a = (double *) t2a[im-1];
     for(j=firstcol;j<=lastcol;j++) {
       t1a[j] = 1.0;
     }
   }
   if (gp[procid].neighbors[LEFT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][0] = 1.0;
     }
   }
   if (gp[procid].neighbors[RIGHT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][jm-1] = 1.0;
     }
   }
   for(i=firstrow;i<=lastrow;i++) {
     t1a = (double *) t2a[i];
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
       t1a[iindex] = 0.0;
     }
   }

/* wait until all processes have completed the above initialization  */
#if defined(MULTIPLE_BARRIERS)
   {
#line 338
	unsigned long	Error, Cycle;
#line 338
	long		Cancel, Temp;
#line 338

#line 338
	Error = pthread_mutex_lock(&(bars->sl_prini).mutex);
#line 338
	if (Error != 0) {
#line 338
		printf("Error while trying to get lock in barrier.\n");
#line 338
		exit(-1);
#line 338
	}
#line 338

#line 338
	Cycle = (bars->sl_prini).cycle;
#line 338
	if (++(bars->sl_prini).counter != (nprocs)) {
#line 338
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 338
		while (Cycle == (bars->sl_prini).cycle) {
#line 338
			Error = pthread_cond_wait(&(bars->sl_prini).cv, &(bars->sl_prini).mutex);
#line 338
			if (Error != 0) {
#line 338
				break;
#line 338
			}
#line 338
		}
#line 338
		pthread_setcancelstate(Cancel, &Temp);
#line 338
	} else {
#line 338
		(bars->sl_prini).cycle = !(bars->sl_prini).cycle;
#line 338
		(bars->sl_prini).counter = 0;
#line 338
		Error = pthread_cond_broadcast(&(bars->sl_prini).cv);
#line 338
	}
#line 338
	pthread_mutex_unlock(&(bars->sl_prini).mutex);
#line 338
}
#else
   {
#line 340
	unsigned long	Error, Cycle;
#line 340
	long		Cancel, Temp;
#line 340

#line 340
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 340
	if (Error != 0) {
#line 340
		printf("Error while trying to get lock in barrier.\n");
#line 340
		exit(-1);
#line 340
	}
#line 340

#line 340
	Cycle = (bars->barrier).cycle;
#line 340
	if (++(bars->barrier).counter != (nprocs)) {
#line 340
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 340
		while (Cycle == (bars->barrier).cycle) {
#line 340
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 340
			if (Error != 0) {
#line 340
				break;
#line 340
			}
#line 340
		}
#line 340
		pthread_setcancelstate(Cancel, &Temp);
#line 340
	} else {
#line 340
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 340
		(bars->barrier).counter = 0;
#line 340
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 340
	}
#line 340
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 340
}
#endif
/* compute psib array (one-time computation) and integrate into psibi */

   istart = 1;
   iend = istart + gp[procid].rel_num_y[numlev-1] - 1;
   jstart = 1;
   jend = jstart + gp[procid].rel_num_x[numlev-1] - 1;
   ist = istart;
   ien = iend;
   jst = jstart;
   jen = jend;

   if (gp[procid].neighbors[UP] == -1) {
     istart = 0;
   }
   if (gp[procid].neighbors[LEFT] == -1) {
     jstart = 0;
   }
   if (gp[procid].neighbors[DOWN] == -1) {
     iend = im-1;
   }
   if (gp[procid].neighbors[RIGHT] == -1) {
     jend = jm-1;
   }

   t2a = (double **) rhs_multi[procid][numlev-1];
   t2b = (double **) psib[procid];
   for(i=istart;i<=iend;i++) {
     t1a = (double *) t2a[i];
     t1b = (double *) t2b[i];
     for(j=jstart;j<=jend;j++) {
       t1a[j] = t1b[j] * ressqr;
     }
   }
   t2a = (double **) q_multi[procid][numlev-1];
   if (gp[procid].neighbors[UP] == -1) {
     t1a = (double *) t2a[0];
     t1b = (double *) t2b[0];
     for(j=jstart;j<=jend;j++) {
       t1a[j] = t1b[j];
     }
   }
   if (gp[procid].neighbors[DOWN] == -1) {
     t1a = (double *) t2a[im-1];
     t1b = (double *) t2b[im-1];
     for(j=jstart;j<=jend;j++) {
       t1a[j] = t1b[j];
     }
   }
   if (gp[procid].neighbors[LEFT] == -1) {
     for(i=istart;i<=iend;i++) {
       t2a[i][0] = t2b[i][0];
     }
   }
   if (gp[procid].neighbors[RIGHT] == -1) {
     for(i=istart;i<=iend;i++) {
       t2a[i][jm-1] = t2b[i][jm-1];
     }
   }
#if defined(MULTIPLE_BARRIERS)
   {
#line 401
	unsigned long	Error, Cycle;
#line 401
	long		Cancel, Temp;
#line 401

#line 401
	Error = pthread_mutex_lock(&(bars->sl_psini).mutex);
#line 401
	if (Error != 0) {
#line 401
		printf("Error while trying to get lock in barrier.\n");
#line 401
		exit(-1);
#line 401
	}
#line 401

#line 401
	Cycle = (bars->sl_psini).cycle;
#line 401
	if (++(bars->sl_psini).counter != (nprocs)) {
#line 401
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 401
		while (Cycle == (bars->sl_psini).cycle) {
#line 401
			Error = pthread_cond_wait(&(bars->sl_psini).cv, &(bars->sl_psini).mutex);
#line 401
			if (Error != 0) {
#line 401
				break;
#line 401
			}
#line 401
		}
#line 401
		pthread_setcancelstate(Cancel, &Temp);
#line 401
	} else {
#line 401
		(bars->sl_psini).cycle = !(bars->sl_psini).cycle;
#line 401
		(bars->sl_psini).counter = 0;
#line 401
		Error = pthread_cond_broadcast(&(bars->sl_psini).cv);
#line 401
	}
#line 401
	pthread_mutex_unlock(&(bars->sl_psini).mutex);
#line 401
}
#else
   {
#line 403
	unsigned long	Error, Cycle;
#line 403
	long		Cancel, Temp;
#line 403

#line 403
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 403
	if (Error != 0) {
#line 403
		printf("Error while trying to get lock in barrier.\n");
#line 403
		exit(-1);
#line 403
	}
#line 403

#line 403
	Cycle = (bars->barrier).cycle;
#line 403
	if (++(bars->barrier).counter != (nprocs)) {
#line 403
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 403
		while (Cycle == (bars->barrier).cycle) {
#line 403
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 403
			if (Error != 0) {
#line 403
				break;
#line 403
			}
#line 403
		}
#line 403
		pthread_setcancelstate(Cancel, &Temp);
#line 403
	} else {
#line 403
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 403
		(bars->barrier).counter = 0;
#line 403
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 403
	}
#line 403
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 403
}
#endif
   t2a = (double **) psib[procid];
   j = gp[procid].neighbors[UP];
   if (j != -1) {
     t1a = (double *) t2a[0];
     t1b = (double *) psib[j][im-2];
     for (i=1;i<jm-1;i++) {
       t1a[i] = t1b[i];
     }
   }
   j = gp[procid].neighbors[DOWN];
   if (j != -1) {
     t1a = (double *) t2a[im-1];
     t1b = (double *) psib[j][1];
     for (i=1;i<jm-1;i++) {
       t1a[i] = t1b[i];
     }
   }
   j = gp[procid].neighbors[LEFT];
   if (j != -1) {
     t2b = (double **) psib[j];
     for (i=1;i<im-1;i++) {
       t2a[i][0] = t2b[i][jm-2];
     }
   }
   j = gp[procid].neighbors[RIGHT];
   if (j != -1) {
     t2b = (double **) psib[j];
     for (i=1;i<im-1;i++) {
       t2a[i][jm-1] = t2b[i][1];
     }
   }

   t2a = (double **) q_multi[procid][numlev-1];
   t2b = (double **) psib[procid];
   fac = 1.0 / (4.0 - ressqr*eig2);
   for(i=ist;i<=ien;i++) {
     t1a = (double *) t2a[i];
     t1b = (double *) t2b[i];
     t1c = (double *) t2b[i-1];
     t1d = (double *) t2b[i+1];
     for(j=jst;j<=jen;j++) {
       t1a[j] = fac * (t1d[j]+t1c[j]+t1b[j+1]+t1b[j-1] -
                   ressqr*t1b[j]);
     }
   }

   multig(procid);

   for(i=istart;i<=iend;i++) {
     t1a = (double *) t2a[i];
     t1b = (double *) t2b[i];
     for(j=jstart;j<=jend;j++) {
       t1b[j] = t1a[j];
     }
   }
#if defined(MULTIPLE_BARRIERS)
   {
#line 461
	unsigned long	Error, Cycle;
#line 461
	long		Cancel, Temp;
#line 461

#line 461
	Error = pthread_mutex_lock(&(bars->sl_prini).mutex);
#line 461
	if (Error != 0) {
#line 461
		printf("Error while trying to get lock in barrier.\n");
#line 461
		exit(-1);
#line 461
	}
#line 461

#line 461
	Cycle = (bars->sl_prini).cycle;
#line 461
	if (++(bars->sl_prini).counter != (nprocs)) {
#line 461
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 461
		while (Cycle == (bars->sl_prini).cycle) {
#line 461
			Error = pthread_cond_wait(&(bars->sl_prini).cv, &(bars->sl_prini).mutex);
#line 461
			if (Error != 0) {
#line 461
				break;
#line 461
			}
#line 461
		}
#line 461
		pthread_setcancelstate(Cancel, &Temp);
#line 461
	} else {
#line 461
		(bars->sl_prini).cycle = !(bars->sl_prini).cycle;
#line 461
		(bars->sl_prini).counter = 0;
#line 461
		Error = pthread_cond_broadcast(&(bars->sl_prini).cv);
#line 461
	}
#line 461
	pthread_mutex_unlock(&(bars->sl_prini).mutex);
#line 461
}
#else
   {
#line 463
	unsigned long	Error, Cycle;
#line 463
	long		Cancel, Temp;
#line 463

#line 463
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 463
	if (Error != 0) {
#line 463
		printf("Error while trying to get lock in barrier.\n");
#line 463
		exit(-1);
#line 463
	}
#line 463

#line 463
	Cycle = (bars->barrier).cycle;
#line 463
	if (++(bars->barrier).counter != (nprocs)) {
#line 463
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 463
		while (Cycle == (bars->barrier).cycle) {
#line 463
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 463
			if (Error != 0) {
#line 463
				break;
#line 463
			}
#line 463
		}
#line 463
		pthread_setcancelstate(Cancel, &Temp);
#line 463
	} else {
#line 463
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 463
		(bars->barrier).counter = 0;
#line 463
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 463
	}
#line 463
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 463
}
#endif
/* update the local running sum psibipriv by summing all the resulting
   values in that process's share of the psib matrix   */

   t2a = (double **) psib[procid];
   psibipriv=0.0;
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     psibipriv = psibipriv + 0.25*(t2a[0][0]);
   }
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     psibipriv = psibipriv + 0.25*(t2a[0][jm-1]);
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     psibipriv=psibipriv+0.25*(t2a[im-1][0]);
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     psibipriv=psibipriv+0.25*(t2a[im-1][jm-1]);
   }
   if (gp[procid].neighbors[UP] == -1) {
     t1a = (double *) t2a[0];
     for(j=firstcol;j<=lastcol;j++) {
       psibipriv = psibipriv + 0.5*t1a[j];
     }
   }
   if (gp[procid].neighbors[DOWN] == -1) {
     t1a = (double *) t2a[im-1];
     for(j=firstcol;j<=lastcol;j++) {
       psibipriv = psibipriv + 0.5*t1a[j];
     }
   }
   if (gp[procid].neighbors[LEFT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       psibipriv = psibipriv + 0.5*t2a[j][0];
     }
   }
   if (gp[procid].neighbors[RIGHT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       psibipriv = psibipriv + 0.5*t2a[j][jm-1];
     }
   }
   for(i=firstrow;i<=lastrow;i++) {
     t1a = (double *) t2a[i];
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
       psibipriv = psibipriv + t1a[iindex];
     }
   }

/* update the shared variable psibi by summing all the psibiprivs
   of the individual processes into it.  note that this combined
   private and shared sum method avoids accessing the shared
   variable psibi once for every element of the matrix.  */

   {pthread_mutex_lock(&(locks->psibilock));}
     global->psibi = global->psibi + psibipriv;
   {pthread_mutex_unlock(&(locks->psibilock));}

/* initialize psim matrices

   if there is more than one process, then split the processes
   between the two psim matrices; otherwise, let the single process
   work on one first and then the other   */

   for(psiindex=0;psiindex<=1;psiindex++) {
     t2a = (double **) psim[procid][psiindex];
     if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
       t2a[0][0] = 0.0;
     }
     if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
       t2a[im-1][0] = 0.0;
     }
     if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
       t2a[0][jm-1] = 0.0;
     }
     if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
       t2a[im-1][jm-1] = 0.0;
     }
     if (gp[procid].neighbors[UP] == -1) {
       t1a = (double *) t2a[0];
       for(j=firstcol;j<=lastcol;j++) {
         t1a[j] = 0.0;
       }
     }
     if (gp[procid].neighbors[DOWN] == -1) {
       t1a = (double *) t2a[im-1];
       for(j=firstcol;j<=lastcol;j++) {
         t1a[j] = 0.0;
       }
     }
     if (gp[procid].neighbors[LEFT] == -1) {
       for(j=firstrow;j<=lastrow;j++) {
         t2a[j][0] = 0.0;
       }
     }
     if (gp[procid].neighbors[RIGHT] == -1) {
       for(j=firstrow;j<=lastrow;j++) {
         t2a[j][jm-1] = 0.0;
       }
     }
     for(i=firstrow;i<=lastrow;i++) {
       t1a = (double *) t2a[i];
       for(iindex=firstcol;iindex<=lastcol;iindex++) {
         t1a[iindex] = 0.0;
       }
     }
   }

/* initialize psi matrices the same way  */

   for(psiindex=0;psiindex<=1;psiindex++) {
     t2a = (double **) psi[procid][psiindex];
     if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
       t2a[0][0] = 0.0;
     }
     if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
       t2a[0][jm-1] = 0.0;
     }
     if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
       t2a[im-1][0] = 0.0;
     }
     if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
       t2a[im-1][jm-1] = 0.0;
     }
     if (gp[procid].neighbors[UP] == -1) {
       t1a = (double *) t2a[0];
       for(j=firstcol;j<=lastcol;j++) {
         t1a[j] = 0.0;
       }
     }
     if (gp[procid].neighbors[DOWN] == -1) {
       t1a = (double *) t2a[im-1];
       for(j=firstcol;j<=lastcol;j++) {
         t1a[j] = 0.0;
       }
     }
     if (gp[procid].neighbors[LEFT] == -1) {
       for(j=firstrow;j<=lastrow;j++) {
         t2a[j][0] = 0.0;
       }
     }
     if (gp[procid].neighbors[RIGHT] == -1) {
       for(j=firstrow;j<=lastrow;j++) {
         t2a[j][jm-1] = 0.0;
       }
     }
     for(i=firstrow;i<=lastrow;i++) {
       t1a = (double *) t2a[i];
       for(iindex=firstcol;iindex<=lastcol;iindex++) {
         t1a[iindex] = 0.0;
       }
     }
   }

/* compute input curl of wind stress */

   t2a = (double **) tauz[procid];
   ysca1 = .5*ysca;
   factor= -t0*pi/ysca1;
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[0][0] = 0.0;
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[im-1][0] = 0.0;
   }
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     sintemp = pi*((double) jm-1+j_off)*res/ysca1;
     sintemp = sin(sintemp);
     t2a[0][jm-1] = factor*sintemp;
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     sintemp = pi*((double) jm-1+j_off)*res/ysca1;
     sintemp = sin(sintemp);
     t2a[im-1][jm-1] = factor*sintemp;
   }
   if (gp[procid].neighbors[UP] == -1) {
     t1a = (double *) t2a[0];
     for(j=firstcol;j<=lastcol;j++) {
       sintemp = pi*((double) j+j_off)*res/ysca1;
       sintemp = sin(sintemp);
       curlt = factor*sintemp;
       t1a[j] = curlt;
     }
   }
   if (gp[procid].neighbors[DOWN] == -1) {
     t1a = (double *) t2a[im-1];
     for(j=firstcol;j<=lastcol;j++) {
       sintemp = pi*((double) j+j_off)*res/ysca1;
       sintemp = sin(sintemp);
       curlt = factor*sintemp;
       t1a[j] = curlt;
     }
   }
   if (gp[procid].neighbors[LEFT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][0] = 0.0;
     }
   }
   if (gp[procid].neighbors[RIGHT] == -1) {
     sintemp = pi*((double) jm-1+j_off)*res/ysca1;
     sintemp = sin(sintemp);
     curlt = factor*sintemp;
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][jm-1] = curlt;
     }
   }
   for(i=firstrow;i<=lastrow;i++) {
     t1a = (double *) t2a[i];
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
       sintemp = pi*((double) iindex+j_off)*res/ysca1;
       sintemp = sin(sintemp);
       curlt = factor*sintemp;
       t1a[iindex] = curlt;
     }
   }
#if defined(MULTIPLE_BARRIERS)
   {
#line 678
	unsigned long	Error, Cycle;
#line 678
	long		Cancel, Temp;
#line 678

#line 678
	Error = pthread_mutex_lock(&(bars->sl_onetime).mutex);
#line 678
	if (Error != 0) {
#line 678
		printf("Error while trying to get lock in barrier.\n");
#line 678
		exit(-1);
#line 678
	}
#line 678

#line 678
	Cycle = (bars->sl_onetime).cycle;
#line 678
	if (++(bars->sl_onetime).counter != (nprocs)) {
#line 678
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 678
		while (Cycle == (bars->sl_onetime).cycle) {
#line 678
			Error = pthread_cond_wait(&(bars->sl_onetime).cv, &(bars->sl_onetime).mutex);
#line 678
			if (Error != 0) {
#line 678
				break;
#line 678
			}
#line 678
		}
#line 678
		pthread_setcancelstate(Cancel, &Temp);
#line 678
	} else {
#line 678
		(bars->sl_onetime).cycle = !(bars->sl_onetime).cycle;
#line 678
		(bars->sl_onetime).counter = 0;
#line 678
		Error = pthread_cond_broadcast(&(bars->sl_onetime).cv);
#line 678
	}
#line 678
	pthread_mutex_unlock(&(bars->sl_onetime).mutex);
#line 678
}
#else
   {
#line 680
	unsigned long	Error, Cycle;
#line 680
	long		Cancel, Temp;
#line 680

#line 680
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 680
	if (Error != 0) {
#line 680
		printf("Error while trying to get lock in barrier.\n");
#line 680
		exit(-1);
#line 680
	}
#line 680

#line 680
	Cycle = (bars->barrier).cycle;
#line 680
	if (++(bars->barrier).counter != (nprocs)) {
#line 680
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 680
		while (Cycle == (bars->barrier).cycle) {
#line 680
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 680
			if (Error != 0) {
#line 680
				break;
#line 680
			}
#line 680
		}
#line 680
		pthread_setcancelstate(Cancel, &Temp);
#line 680
	} else {
#line 680
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 680
		(bars->barrier).counter = 0;
#line 680
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 680
	}
#line 680
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 680
}
#endif

/***************************************************************
 one-time stuff over at this point
 ***************************************************************/

   while (!endflag) {
     while ((!dayflag) || (!dhourflag)) {
       dayflag = 0;
       dhourflag = 0;
       if (nstep == 1) {
         if (procid == MASTER) {
            {
#line 693
	struct timeval	FullTime;
#line 693

#line 693
	gettimeofday(&FullTime, NULL);
#line 693
	(global->trackstart) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 693
}
         }
	 if ((procid == MASTER) || (do_stats)) {
	   {
#line 696
	struct timeval	FullTime;
#line 696

#line 696
	gettimeofday(&FullTime, NULL);
#line 696
	(t1) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 696
};
           gp[procid].total_time = t1;
           gp[procid].multi_time = 0;
	 }
/* POSSIBLE ENHANCEMENT:  Here is where one might reset the
   statistics that one is measuring about the parallel execution */
       }

       slave2(procid,firstrow,lastrow,numrows,firstcol,lastcol,numcols);

/* update time and step number
   note that these time and step variables are private i.e. every
   process has its own copy and keeps track of its own time  */

       ttime = ttime + dtau;
       nstep = nstep + 1;
       day = ttime/86400.0;

       if (day > ((double) outday0)) {
         dayflag = 1;
         iday = (long) day;
         dhour = dhour+dtau;
         if (dhour >= 86400.0) {
	   dhourflag = 1;
         }
       }
     }
     dhour = 0.0;

     t2a = (double **) psium[procid];
     t2b = (double **) psim[procid][0];
     if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
       t2a[0][0] = t2a[0][0]+t2b[0][0];
     }
     if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
       t2a[im-1][0] = t2a[im-1][0]+t2b[im-1][0];
     }
     if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
       t2a[0][jm-1] = t2a[0][jm-1]+t2b[0][jm-1];
     }
     if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
       t2a[im-1][jm-1] = t2a[im-1][jm-1] +
				   t2b[im-1][jm-1];
     }
     if (gp[procid].neighbors[UP] == -1) {
       t1a = (double *) t2a[0];
       t1b = (double *) t2b[0];
       for(j=firstcol;j<=lastcol;j++) {
         t1a[j] = t1a[j]+t1b[j];
       }
     }
     if (gp[procid].neighbors[DOWN] == -1) {
       t1a = (double *) t2a[im-1];
       t1b = (double *) t2b[im-1];
       for(j=firstcol;j<=lastcol;j++) {
         t1a[j] = t1a[j] + t1b[j];
       }
     }
     if (gp[procid].neighbors[LEFT] == -1) {
       for(j=firstrow;j<=lastrow;j++) {
         t2a[j][0] = t2a[j][0]+t2b[j][0];
       }
     }
     if (gp[procid].neighbors[RIGHT] == -1) {
       for(j=firstrow;j<=lastrow;j++) {
         t2a[j][jm-1] = t2a[j][jm-1] +
				  t2b[j][jm-1];
       }
     }
     for(i=firstrow;i<=lastrow;i++) {
       t1a = (double *) t2a[i];
       t1b = (double *) t2b[i];
       for(iindex=firstcol;iindex<=lastcol;iindex++) {
         t1a[iindex] = t1a[iindex] + t1b[iindex];
       }
     }

/* update values of psilm array to psilm + psim[2]  */

     t2a = (double **) psilm[procid];
     t2b = (double **) psim[procid][1];
     if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
       t2a[0][0] = t2a[0][0]+t2b[0][0];
     }
     if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
       t2a[im-1][0] = t2a[im-1][0]+t2b[im-1][0];
     }
     if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
       t2a[0][jm-1] = t2a[0][jm-1]+t2b[0][jm-1];
     }
     if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
       t2a[im-1][jm-1] = t2a[im-1][jm-1] +
				   t2b[im-1][jm-1];
     }
     if (gp[procid].neighbors[UP] == -1) {
       t1a = (double *) t2a[0];
       t1b = (double *) t2b[0];
       for(j=firstcol;j<=lastcol;j++) {
         t1a[j] = t1a[j]+t1b[j];
       }
     }
     if (gp[procid].neighbors[DOWN] == -1) {
       t1a = (double *) t2a[im-1];
       t1b = (double *) t2b[im-1];
       for(j=firstcol;j<=lastcol;j++) {
         t1a[j] = t1a[j]+t1b[j];
       }
     }
     if (gp[procid].neighbors[LEFT] == -1) {
       for(j=firstrow;j<=lastrow;j++) {
         t2a[j][0] = t2a[j][0]+t2b[j][0];
       }
     }
     if (gp[procid].neighbors[RIGHT] == -1) {
       for(j=firstrow;j<=lastrow;j++) {
         t2a[j][jm-1] = t2a[j][jm-1] + t2b[j][jm-1];
       }
     }
     for(i=firstrow;i<=lastrow;i++) {
       t1a = (double *) t2a[i];
       t1b = (double *) t2b[i];
       for(iindex=firstcol;iindex<=lastcol;iindex++) {
         t1a[iindex] = t1a[iindex] + t1b[iindex];
       }
     }
     if (iday >= (long) outday3) {
       endflag = 1;
     }
  }
  if ((procid == MASTER) || (do_stats)) {
    {
#line 826
	struct timeval	FullTime;
#line 826

#line 826
	gettimeofday(&FullTime, NULL);
#line 826
	(t1) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 826
};
    gp[procid].total_time = t1-gp[procid].total_time;
  }
}

#line 228 "/home/jyy/Dev/splash2/codes/null_macros/c.m4.null.POSIX"

#line 1 "slave2.C"
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

/*    ****************
      subroutine slave2
      ****************  */


#line 21
#include <pthread.h>
#line 21
#include <sys/time.h>
#line 21
#include <unistd.h>
#line 21
#include <stdlib.h>
#line 21
extern pthread_t PThreadTable[];
#line 21


#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "decs.h"

void slave2(long procid, long firstrow, long lastrow, long numrows, long firstcol, long lastcol, long numcols)
{
   long i;
   long j;
   long iindex;
   double hh1;
   double hh3;
   double hinv;
   double h1inv;
   long istart;
   long iend;
   long jstart;
   long jend;
   long ist;
   long ien;
   long jst;
   long jen;
   double fac;
   double ressqr;
   double psiaipriv;
   double f4;
   double timst;
   long psiindex;
   long i_off;
   long j_off;
   long multi_start;
   long multi_end;
   double **t2a;
   double **t2b;
   double **t2c;
   double **t2d;
   double **t2e;
   double **t2f;
   double **t2g;
   double **t2h;
   double *t1a;
   double *t1b;
   double *t1c;
   double *t1d;
   double *t1e;
   double *t1f;
   double *t1g;
   double *t1h;

   ressqr = lev_res[numlev-1] * lev_res[numlev-1];
   i_off = gp[procid].rownum*numrows;
   j_off = gp[procid].colnum*numcols;

/*   ***************************************************************

          f i r s t     p h a s e   (of timestep calculation)

     ***************************************************************/

   t2a = (double **) ga[procid];
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[0][0]=0.0;
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[im-1][0]=0.0;
   }
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[0][jm-1]=0.0;
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[im-1][jm-1]=0.0;
   }
   if (gp[procid].neighbors[UP] == -1) {
     t1a = (double *) t2a[0];
     for(j=firstcol;j<=lastcol;j++) {
       t1a[j] = 0.0;
     }
   }
   if (gp[procid].neighbors[DOWN] == -1) {
     t1a = (double *) t2a[im-1];
     for(j=firstcol;j<=lastcol;j++) {
       t1a[j] = 0.0;
     }
   }
   if (gp[procid].neighbors[LEFT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][0] = 0.0;
     }
   }
   if (gp[procid].neighbors[RIGHT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][jm-1] = 0.0;
     }
   }
   for(i=firstrow;i<=lastrow;i++) {
     t1a = (double *) t2a[i];
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
       t1a[iindex] = 0.0;
     }
   }

   t2a = (double **) gb[procid];
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[0][0]=0.0;
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[im-1][0]=0.0;
   }
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[0][jm-1]=0.0;
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[im-1][jm-1]=0.0;
   }
   if (gp[procid].neighbors[UP] == -1) {
     t1a = (double *) t2a[0];
     for(j=firstcol;j<=lastcol;j++) {
       t1a[j] = 0.0;
     }
   }
   if (gp[procid].neighbors[DOWN] == -1) {
     t1a = (double *) t2a[im-1];
     for(j=firstcol;j<=lastcol;j++) {
       t1a[j] = 0.0;
     }
   }
   if (gp[procid].neighbors[LEFT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][0] = 0.0;
     }
   }
   if (gp[procid].neighbors[RIGHT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][jm-1] = 0.0;
     }
   }
   for(i=firstrow;i<=lastrow;i++) {
     t1a = (double *) t2a[i];
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
       t1a[iindex] = 0.0;
     }
   }

/* put the laplacian of psi{1,3} in work1{1,2}
   note that psi(i,j,2) represents the psi3 array in
   the original equations  */

   for(psiindex=0;psiindex<=1;psiindex++) {
     t2a = (double **) work1[procid][psiindex];
     if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
       t2a[0][0] = 0;
     }
     if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
       t2a[im-1][0] = 0;
     }
     if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
       t2a[0][jm-1] = 0;
     }
     if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
       t2a[im-1][jm-1] = 0;
     }
     laplacalc(procid,psi,work1,psiindex,
	       firstrow,lastrow,firstcol,lastcol);
   }

/* set values of work2 array to psi1 - psi3   */

   t2a = (double **) work2[procid];
   t2b = (double **) psi[procid][0];
   t2c = (double **) psi[procid][1];
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[0][0] = t2b[0][0]-t2c[0][0];
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[im-1][0] = t2b[im-1][0]-t2c[im-1][0];
   }
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[0][jm-1] = t2b[0][jm-1]-t2c[0][jm-1];
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[im-1][jm-1] = t2b[im-1][jm-1] -
				 t2c[im-1][jm-1];
   }
   if (gp[procid].neighbors[UP] == -1) {
     t1a = (double *) t2a[0];
     t1b = (double *) t2b[0];
     t1c = (double *) t2c[0];
     for(j=firstcol;j<=lastcol;j++) {
       t1a[j] = t1b[j]-t1c[j];
     }
   }
   if (gp[procid].neighbors[DOWN] == -1) {
     t1a = (double *) t2a[im-1];
     t1b = (double *) t2b[im-1];
     t1c = (double *) t2c[im-1];
     for(j=firstcol;j<=lastcol;j++) {
       t1a[j] = t1b[j]-t1c[j];
     }
   }
   if (gp[procid].neighbors[LEFT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][0] = t2b[j][0]-t2c[j][0];
     }
   }
   if (gp[procid].neighbors[RIGHT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][jm-1] = t2b[j][jm-1]-t2c[j][jm-1];
     }
   }
   for(i=firstrow;i<=lastrow;i++) {
     t1a = (double *) t2a[i];
     t1b = (double *) t2b[i];
     t1c = (double *) t2c[i];
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
         t1a[iindex] = t1b[iindex] - t1c[iindex];
     }
   }

/* set values of work3 array to h3/h * psi1 + h1/h * psi3  */

   t2a = (double **) work3[procid];
   hh3 = h3/h;
   hh1 = h1/h;
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[0][0] = hh3*t2a[0][0]+hh1*t2c[0][0];
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[im-1][0] = hh3*t2a[im-1][0] +
			      hh1*t2c[im-1][0];
   }
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[0][jm-1] = hh3*t2a[0][jm-1] +
			      hh1*t2c[0][jm-1];
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[im-1][jm-1] = hh3*t2a[im-1][jm-1] +
				 hh1*t2c[im-1][jm-1];
   }
   if (gp[procid].neighbors[UP] == -1) {
     for(j=firstcol;j<=lastcol;j++) {
       t2a[0][j] = hh3*t2a[0][j]+hh1*t2c[0][j];
     }
   }
   if (gp[procid].neighbors[DOWN] == -1) {
     for(j=firstcol;j<=lastcol;j++) {
       t2a[im-1][j] = hh3*t2a[im-1][j] +
				hh1*t2c[im-1][j];
     }
   }
   if (gp[procid].neighbors[LEFT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][0] = hh3*t2a[j][0]+hh1*t2c[j][0];
     }
   }
   if (gp[procid].neighbors[RIGHT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][jm-1] = hh3*t2a[j][jm-1] +
				hh1*t2c[j][jm-1];
     }
   }
   for(i=firstrow;i<=lastrow;i++) {
     t1a = (double *) t2a[i];
     t1c = (double *) t2c[i];
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
        t1a[iindex] = hh3*t1a[iindex] + hh1*t1c[iindex];
     }
   }

/* set values of temparray{1,3} to psim{1,3}  */

   for(psiindex=0;psiindex<=1;psiindex++) {
     t2a = (double **) temparray[procid][psiindex];
     t2b = (double **) psi[procid][psiindex];
     if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
       t2a[0][0] = t2b[0][0];
     }
     if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
       t2a[im-1][0] = t2b[im-1][0];
     }
     if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
       t2a[0][jm-1] = t2b[0][jm-1];
     }
     if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
       t2a[im-1][jm-1] = t2b[im-1][jm-1];
     }
     if (gp[procid].neighbors[UP] == -1) {
       for(j=firstcol;j<=lastcol;j++) {
         t2a[0][j] = t2b[0][j];
       }
     }
     if (gp[procid].neighbors[DOWN] == -1) {
       for(j=firstcol;j<=lastcol;j++) {
         t2a[im-1][j] = t2b[im-1][j];
       }
     }
     if (gp[procid].neighbors[LEFT] == -1) {
       for(j=firstrow;j<=lastrow;j++) {
         t2a[j][0] = t2b[j][0];
       }
     }
     if (gp[procid].neighbors[RIGHT] == -1) {
       for(j=firstrow;j<=lastrow;j++) {
         t2a[j][jm-1] = t2b[j][jm-1];
       }
     }

     for(i=firstrow;i<=lastrow;i++) {
       t1a = (double *) t2a[i];
       t1b = (double *) t2b[i];
       for(iindex=firstcol;iindex<=lastcol;iindex++) {
         t1a[iindex] = t1b[iindex];
       }
     }
   }
#if defined(MULTIPLE_BARRIERS)
   {
#line 339
	unsigned long	Error, Cycle;
#line 339
	long		Cancel, Temp;
#line 339

#line 339
	Error = pthread_mutex_lock(&(bars->sl_phase_1).mutex);
#line 339
	if (Error != 0) {
#line 339
		printf("Error while trying to get lock in barrier.\n");
#line 339
		exit(-1);
#line 339
	}
#line 339

#line 339
	Cycle = (bars->sl_phase_1).cycle;
#line 339
	if (++(bars->sl_phase_1).counter != (nprocs)) {
#line 339
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 339
		while (Cycle == (bars->sl_phase_1).cycle) {
#line 339
			Error = pthread_cond_wait(&(bars->sl_phase_1).cv, &(bars->sl_phase_1).mutex);
#line 339
			if (Error != 0) {
#line 339
				break;
#line 339
			}
#line 339
		}
#line 339
		pthread_setcancelstate(Cancel, &Temp);
#line 339
	} else {
#line 339
		(bars->sl_phase_1).cycle = !(bars->sl_phase_1).cycle;
#line 339
		(bars->sl_phase_1).counter = 0;
#line 339
		Error = pthread_cond_broadcast(&(bars->sl_phase_1).cv);
#line 339
	}
#line 339
	pthread_mutex_unlock(&(bars->sl_phase_1).mutex);
#line 339
}
#else
   {
#line 341
	unsigned long	Error, Cycle;
#line 341
	long		Cancel, Temp;
#line 341

#line 341
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 341
	if (Error != 0) {
#line 341
		printf("Error while trying to get lock in barrier.\n");
#line 341
		exit(-1);
#line 341
	}
#line 341

#line 341
	Cycle = (bars->barrier).cycle;
#line 341
	if (++(bars->barrier).counter != (nprocs)) {
#line 341
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 341
		while (Cycle == (bars->barrier).cycle) {
#line 341
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 341
			if (Error != 0) {
#line 341
				break;
#line 341
			}
#line 341
		}
#line 341
		pthread_setcancelstate(Cancel, &Temp);
#line 341
	} else {
#line 341
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 341
		(bars->barrier).counter = 0;
#line 341
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 341
	}
#line 341
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 341
}
#endif
/*     *******************************************************

              s e c o n d   p h a s e

       *******************************************************

   set values of psi{1,3} to psim{1,3}   */

   for(psiindex=0;psiindex<=1;psiindex++) {
     t2a = (double **) psi[procid][psiindex];
     t2b = (double **) psim[procid][psiindex];
     if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
       t2a[0][0] = t2b[0][0];
     }
     if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
       t2a[0][jm-1] = t2b[0][jm-1];
     }
     if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
       t2a[im-1][0] = t2b[im-1][0];
     }
     if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
       t2a[im-1][jm-1] = t2b[im-1][jm-1];
     }
     if (gp[procid].neighbors[UP] == -1) {
       for(j=firstcol;j<=lastcol;j++) {
         t2a[0][j] = t2b[0][j];
       }
     }
     if (gp[procid].neighbors[DOWN] == -1) {
       for(j=firstcol;j<=lastcol;j++) {
         t2a[im-1][j] = t2b[im-1][j];
       }
     }
     if (gp[procid].neighbors[LEFT] == -1) {
       for(j=firstrow;j<=lastrow;j++) {
         t2a[j][0] = t2b[j][0];
       }
     }
     if (gp[procid].neighbors[RIGHT] == -1) {
       for(j=firstrow;j<=lastrow;j++) {
         t2a[j][jm-1] = t2b[j][jm-1];
       }
     }

     for(i=firstrow;i<=lastrow;i++) {
       t1a = (double *) t2a[i];
       t1b = (double *) t2b[i];
       for(iindex=firstcol;iindex<=lastcol;iindex++) {
         t1a[iindex] = t1b[iindex];
       }
     }
   }

/* put the laplacian of the psim array
   into the work7 array; first part of a three-laplacian
   calculation to compute the friction terms  */

   for(psiindex=0;psiindex<=1;psiindex++) {
     t2a = (double **) work7[procid][psiindex];
     if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
       t2a[0][0] = 0;
     }
     if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
       t2a[im-1][0] = 0;
     }
     if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
       t2a[0][jm-1] = 0;
     }
     if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
       t2a[im-1][jm-1] = 0;
     }
     laplacalc(procid,psim,work7,psiindex,
	       firstrow,lastrow,firstcol,lastcol);
   }

/* to the values of the work1{1,2} arrays obtained from the
   laplacians of psi{1,2} in the previous phase, add to the
   elements of every column the corresponding value in the
   one-dimenional f array  */

   for(psiindex=0;psiindex<=1;psiindex++) {
     t2a = (double **) work1[procid][psiindex];
     if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
       t2a[0][0] = t2a[0][0] + f[0];
     }
     if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
       t2a[im-1][0] = t2a[im-1][0] + f[0];
     }
     if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
       t2a[0][jm-1] = t2a[0][jm-1] + f[jmx[numlev-1]-1];
     }
     if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
       t2a[im-1][jm-1]=t2a[im-1][jm-1] + f[jmx[numlev-1]-1];
     }
     if (gp[procid].neighbors[UP] == -1) {
       for(j=firstcol;j<=lastcol;j++) {
         t2a[0][j] = t2a[0][j] + f[j+j_off];
       }
     }
     if (gp[procid].neighbors[DOWN] == -1) {
       for(j=firstcol;j<=lastcol;j++) {
         t2a[im-1][j] = t2a[im-1][j] + f[j+j_off];
       }
     }
     if (gp[procid].neighbors[LEFT] == -1) {
       for(j=firstrow;j<=lastrow;j++) {
         t2a[j][0] = t2a[j][0] + f[j+i_off];
       }
     }
     if (gp[procid].neighbors[RIGHT] == -1) {
       for(j=firstrow;j<=lastrow;j++) {
         t2a[j][jm-1] = t2a[j][jm-1] + f[j+i_off];
       }
     }
     for(i=firstrow;i<=lastrow;i++) {
       t1a = (double *) t2a[i];
       for(iindex=firstcol;iindex<=lastcol;iindex++) {
         t1a[iindex]=t1a[iindex] + f[iindex+j_off];
       }
     }
   }
#if defined(MULTIPLE_BARRIERS)
   {
#line 465
	unsigned long	Error, Cycle;
#line 465
	long		Cancel, Temp;
#line 465

#line 465
	Error = pthread_mutex_lock(&(bars->sl_phase_2).mutex);
#line 465
	if (Error != 0) {
#line 465
		printf("Error while trying to get lock in barrier.\n");
#line 465
		exit(-1);
#line 465
	}
#line 465

#line 465
	Cycle = (bars->sl_phase_2).cycle;
#line 465
	if (++(bars->sl_phase_2).counter != (nprocs)) {
#line 465
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 465
		while (Cycle == (bars->sl_phase_2).cycle) {
#line 465
			Error = pthread_cond_wait(&(bars->sl_phase_2).cv, &(bars->sl_phase_2).mutex);
#line 465
			if (Error != 0) {
#line 465
				break;
#line 465
			}
#line 465
		}
#line 465
		pthread_setcancelstate(Cancel, &Temp);
#line 465
	} else {
#line 465
		(bars->sl_phase_2).cycle = !(bars->sl_phase_2).cycle;
#line 465
		(bars->sl_phase_2).counter = 0;
#line 465
		Error = pthread_cond_broadcast(&(bars->sl_phase_2).cv);
#line 465
	}
#line 465
	pthread_mutex_unlock(&(bars->sl_phase_2).mutex);
#line 465
}
#else
   {
#line 467
	unsigned long	Error, Cycle;
#line 467
	long		Cancel, Temp;
#line 467

#line 467
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 467
	if (Error != 0) {
#line 467
		printf("Error while trying to get lock in barrier.\n");
#line 467
		exit(-1);
#line 467
	}
#line 467

#line 467
	Cycle = (bars->barrier).cycle;
#line 467
	if (++(bars->barrier).counter != (nprocs)) {
#line 467
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 467
		while (Cycle == (bars->barrier).cycle) {
#line 467
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 467
			if (Error != 0) {
#line 467
				break;
#line 467
			}
#line 467
		}
#line 467
		pthread_setcancelstate(Cancel, &Temp);
#line 467
	} else {
#line 467
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 467
		(bars->barrier).counter = 0;
#line 467
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 467
	}
#line 467
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 467
}
#endif
/* 	*******************************************************

                 t h i r d   p h a s e

 	*******************************************************

   put the jacobian of the work1{1,2} and psi{1,3} arrays
   (the latter currently in temparray) in the work5{1,2} arrays  */

   for(psiindex=0;psiindex<=1;psiindex++) {
     jacobcalc2(work1,temparray,work5,psiindex,procid,firstrow,lastrow,
	       firstcol,lastcol);
   }

/* set values of psim{1,3} to temparray{1,3}  */

   for(psiindex=0;psiindex<=1;psiindex++) {
     t2a = (double **) psim[procid][psiindex];
     t2b = (double **) temparray[procid][psiindex];
     if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
       t2a[0][0] = t2b[0][0];
     }
     if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
       t2a[im-1][0] = t2b[im-1][0];
     }
     if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
       t2a[0][jm-1] = t2b[0][jm-1];
     }
     if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
       t2a[im-1][jm-1] = t2b[im-1][jm-1];
     }
     if (gp[procid].neighbors[UP] == -1) {
       t1a = (double *) t2a[0];
       t1b = (double *) t2b[0];
       for(j=firstcol;j<=lastcol;j++) {
         t1a[j] = t1b[j];
       }
     }
     if (gp[procid].neighbors[DOWN] == -1) {
       t1a = (double *) t2a[im-1];
       t1b = (double *) t2b[im-1];
       for(j=firstcol;j<=lastcol;j++) {
         t1a[j] = t1b[j];
       }
     }
     if (gp[procid].neighbors[LEFT] == -1) {
       for(j=firstrow;j<=lastrow;j++) {
         t2a[j][0] = t2b[j][0];
       }
     }
     if (gp[procid].neighbors[RIGHT] == -1) {
       for(j=firstrow;j<=lastrow;j++) {
         t2a[j][jm-1] = t2b[j][jm-1];
       }
     }
     for(i=firstrow;i<=lastrow;i++) {
       t1a = (double *) t2a[i];
       t1b = (double *) t2b[i];
       for(iindex=firstcol;iindex<=lastcol;iindex++) {
         t1a[iindex] = t1b[iindex];
       }
     }
   }

/* put the laplacian of the work7{1,2} arrays in the work4{1,2}
   arrays; second step in the three-laplacian friction calculation  */

   for(psiindex=0;psiindex<=1;psiindex++) {
     laplacalc(procid,work7,work4,psiindex,
	       firstrow,lastrow,firstcol,lastcol);
   }
#if defined(MULTIPLE_BARRIERS)
   {
#line 541
	unsigned long	Error, Cycle;
#line 541
	long		Cancel, Temp;
#line 541

#line 541
	Error = pthread_mutex_lock(&(bars->sl_phase_3).mutex);
#line 541
	if (Error != 0) {
#line 541
		printf("Error while trying to get lock in barrier.\n");
#line 541
		exit(-1);
#line 541
	}
#line 541

#line 541
	Cycle = (bars->sl_phase_3).cycle;
#line 541
	if (++(bars->sl_phase_3).counter != (nprocs)) {
#line 541
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 541
		while (Cycle == (bars->sl_phase_3).cycle) {
#line 541
			Error = pthread_cond_wait(&(bars->sl_phase_3).cv, &(bars->sl_phase_3).mutex);
#line 541
			if (Error != 0) {
#line 541
				break;
#line 541
			}
#line 541
		}
#line 541
		pthread_setcancelstate(Cancel, &Temp);
#line 541
	} else {
#line 541
		(bars->sl_phase_3).cycle = !(bars->sl_phase_3).cycle;
#line 541
		(bars->sl_phase_3).counter = 0;
#line 541
		Error = pthread_cond_broadcast(&(bars->sl_phase_3).cv);
#line 541
	}
#line 541
	pthread_mutex_unlock(&(bars->sl_phase_3).mutex);
#line 541
}
#else
   {
#line 543
	unsigned long	Error, Cycle;
#line 543
	long		Cancel, Temp;
#line 543

#line 543
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 543
	if (Error != 0) {
#line 543
		printf("Error while trying to get lock in barrier.\n");
#line 543
		exit(-1);
#line 543
	}
#line 543

#line 543
	Cycle = (bars->barrier).cycle;
#line 543
	if (++(bars->barrier).counter != (nprocs)) {
#line 543
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 543
		while (Cycle == (bars->barrier).cycle) {
#line 543
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 543
			if (Error != 0) {
#line 543
				break;
#line 543
			}
#line 543
		}
#line 543
		pthread_setcancelstate(Cancel, &Temp);
#line 543
	} else {
#line 543
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 543
		(bars->barrier).counter = 0;
#line 543
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 543
	}
#line 543
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 543
}
#endif
/*     *******************************************************

                f o u r t h   p h a s e

       *******************************************************

   put the jacobian of the work2 and work3 arrays in the work6
   array  */

   jacobcalc(work2,work3,work6,procid,firstrow,lastrow,firstcol,lastcol);

/* put the laplacian of the work4{1,2} arrays in the work7{1,2}
   arrays; third step in the three-laplacian friction calculation  */

   for(psiindex=0;psiindex<=1;psiindex++) {
     laplacalc(procid,work4,work7,psiindex,
	       firstrow,lastrow,firstcol,lastcol);
   }
#if defined(MULTIPLE_BARRIERS)
   {
#line 564
	unsigned long	Error, Cycle;
#line 564
	long		Cancel, Temp;
#line 564

#line 564
	Error = pthread_mutex_lock(&(bars->sl_phase_4).mutex);
#line 564
	if (Error != 0) {
#line 564
		printf("Error while trying to get lock in barrier.\n");
#line 564
		exit(-1);
#line 564
	}
#line 564

#line 564
	Cycle = (bars->sl_phase_4).cycle;
#line 564
	if (++(bars->sl_phase_4).counter != (nprocs)) {
#line 564
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 564
		while (Cycle == (bars->sl_phase_4).cycle) {
#line 564
			Error = pthread_cond_wait(&(bars->sl_phase_4).cv, &(bars->sl_phase_4).mutex);
#line 564
			if (Error != 0) {
#line 564
				break;
#line 564
			}
#line 564
		}
#line 564
		pthread_setcancelstate(Cancel, &Temp);
#line 564
	} else {
#line 564
		(bars->sl_phase_4).cycle = !(bars->sl_phase_4).cycle;
#line 564
		(bars->sl_phase_4).counter = 0;
#line 564
		Error = pthread_cond_broadcast(&(bars->sl_phase_4).cv);
#line 564
	}
#line 564
	pthread_mutex_unlock(&(bars->sl_phase_4).mutex);
#line 564
}
#else
   {
#line 566
	unsigned long	Error, Cycle;
#line 566
	long		Cancel, Temp;
#line 566

#line 566
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 566
	if (Error != 0) {
#line 566
		printf("Error while trying to get lock in barrier.\n");
#line 566
		exit(-1);
#line 566
	}
#line 566

#line 566
	Cycle = (bars->barrier).cycle;
#line 566
	if (++(bars->barrier).counter != (nprocs)) {
#line 566
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 566
		while (Cycle == (bars->barrier).cycle) {
#line 566
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 566
			if (Error != 0) {
#line 566
				break;
#line 566
			}
#line 566
		}
#line 566
		pthread_setcancelstate(Cancel, &Temp);
#line 566
	} else {
#line 566
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 566
		(bars->barrier).counter = 0;
#line 566
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 566
	}
#line 566
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 566
}
#endif
/*     *******************************************************

                f i f t h   p h a s e

       *******************************************************

   use the values of the work5, work6 and work7 arrays
   computed in the previous time-steps to compute the
   ga and gb arrays   */

   hinv = 1.0/h;
   h1inv = 1.0/h1;

   t2a = (double **) ga[procid];
   t2b = (double **) gb[procid];
   t2c = (double **) work5[procid][0];
   t2d = (double **) work5[procid][1];
   t2e = (double **) work7[procid][0];
   t2f = (double **) work7[procid][1];
   t2g = (double **) work6[procid];
   t2h = (double **) tauz[procid];
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[0][0] = t2c[0][0]-t2d[0][0] +
			eig2*t2g[0][0]+h1inv*t2h[0][0] +
			lf*t2e[0][0]-lf*t2f[0][0];
     t2b[0][0] = hh1*t2c[0][0]+hh3*t2d[0][0] +
			hinv*t2h[0][0]+lf*hh1*t2e[0][0] +
		        lf*hh3*t2f[0][0];
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[im-1][0] = t2c[im-1][0]-t2d[im-1][0] +
	   eig2*t2g[im-1][0] + h1inv*t2h[im-1][0] +
	   lf*t2e[im-1][0] - lf*t2f[im-1][0];
     t2b[im-1][0] = hh1*t2c[im-1][0] +
	   hh3*t2d[im-1][0] + hinv*t2h[im-1][0] +
	   lf*hh1*t2e[im-1][0] + lf*hh3*t2f[im-1][0];
   }
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[0][jm-1] = t2c[0][jm-1]-t2d[0][jm-1]+
	   eig2*t2g[0][jm-1]+h1inv*t2h[0][jm-1] +
	   lf*t2e[0][jm-1]-lf*t2f[0][jm-1];
     t2b[0][jm-1] = hh1*t2c[0][jm-1] +
	   hh3*t2d[0][jm-1]+hinv*t2h[0][jm-1] +
	   lf*hh1*t2e[0][jm-1]+lf*hh3*t2f[0][jm-1];
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[im-1][jm-1] = t2c[im-1][jm-1] -
	   t2d[im-1][jm-1]+eig2*t2g[im-1][jm-1] +
	   h1inv*t2h[im-1][jm-1]+lf*t2e[im-1][jm-1] -
	   lf*t2f[im-1][jm-1];
     t2b[im-1][jm-1] = hh1*t2c[im-1][jm-1] +
	   hh3*t2d[im-1][jm-1]+hinv*t2h[im-1][jm-1] +
	   lf*hh1*t2e[im-1][jm-1] +
	   lf*hh3*t2f[im-1][jm-1];
   }
   if (gp[procid].neighbors[UP] == -1) {
     t1a = (double *) t2a[0];
     t1b = (double *) t2b[0];
     t1c = (double *) t2c[0];
     t1d = (double *) t2d[0];
     t1e = (double *) t2e[0];
     t1f = (double *) t2f[0];
     t1g = (double *) t2g[0];
     t1h = (double *) t2h[0];
     for(j=firstcol;j<=lastcol;j++) {
       t1a[j] = t1c[j]-t1d[j] +
	   eig2*t1g[j]+h1inv*t1h[j] +
	   lf*t1e[j]-lf*t1f[j];
       t1b[j] = hh1*t1c[j] +
	   hh3*t1d[j]+hinv*t1h[j] +
	   lf*hh1*t1e[j]+lf*hh3*t1f[j];
     }
   }
   if (gp[procid].neighbors[DOWN] == -1) {
     t1a = (double *) t2a[im-1];
     t1b = (double *) t2b[im-1];
     t1c = (double *) t2c[im-1];
     t1d = (double *) t2d[im-1];
     t1e = (double *) t2e[im-1];
     t1f = (double *) t2f[im-1];
     t1g = (double *) t2g[im-1];
     t1h = (double *) t2h[im-1];
     for(j=firstcol;j<=lastcol;j++) {
       t1a[j] = t1c[j] -
	   t1d[j]+eig2*t1g[j] +
	   h1inv*t1h[j]+lf*t1e[j] -
	   lf*t1f[j];
       t1b[j] = hh1*t1c[j] +
	   hh3*t1d[j]+hinv*t1h[j] +
	   lf*hh1*t1e[j]+lf*hh3*t1f[j];
     }
   }
   if (gp[procid].neighbors[LEFT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][0] = t2c[j][0]-t2d[j][0] +
	   eig2*t2g[j][0]+h1inv*t2h[j][0] +
	   lf*t2e[j][0]-lf*t2f[j][0];
       t2b[j][0] = hh1*t2c[j][0] +
	   hh3*t2d[j][0]+hinv*t2h[j][0] +
	   lf*hh1*t2e[j][0]+lf*hh3*t2f[j][0];
     }
   }
   if (gp[procid].neighbors[RIGHT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][jm-1] = t2c[j][jm-1] -
	   t2d[j][jm-1]+eig2*t2g[j][jm-1] +
	   h1inv*t2h[j][jm-1]+lf*t2e[j][jm-1] -
	   lf*t2f[j][jm-1];
       t2b[j][jm-1] = hh1*t2c[j][jm-1] +
	   hh3*t2d[j][jm-1]+hinv*t2h[j][jm-1] +
	   lf*hh1*t2e[j][jm-1]+lf*hh3*t2f[j][jm-1];
     }
   }

   for(i=firstrow;i<=lastrow;i++) {
     t1a = (double *) t2a[i];
     t1b = (double *) t2b[i];
     t1c = (double *) t2c[i];
     t1d = (double *) t2d[i];
     t1e = (double *) t2e[i];
     t1f = (double *) t2f[i];
     t1g = (double *) t2g[i];
     t1h = (double *) t2h[i];
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
       t1a[iindex] = t1c[iindex] -
	   t1d[iindex]+eig2*t1g[iindex] +
	   h1inv*t1h[iindex]+lf*t1e[iindex] -
	   lf*t1f[iindex];
       t1b[iindex] = hh1*t1c[iindex] +
	   hh3*t1d[iindex]+hinv*t1h[iindex] +
	   lf*hh1*t1e[iindex] +
	   lf*hh3*t1f[iindex];
     }
   }
#if defined(MULTIPLE_BARRIERS)
   {
#line 703
	unsigned long	Error, Cycle;
#line 703
	long		Cancel, Temp;
#line 703

#line 703
	Error = pthread_mutex_lock(&(bars->sl_phase_5).mutex);
#line 703
	if (Error != 0) {
#line 703
		printf("Error while trying to get lock in barrier.\n");
#line 703
		exit(-1);
#line 703
	}
#line 703

#line 703
	Cycle = (bars->sl_phase_5).cycle;
#line 703
	if (++(bars->sl_phase_5).counter != (nprocs)) {
#line 703
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 703
		while (Cycle == (bars->sl_phase_5).cycle) {
#line 703
			Error = pthread_cond_wait(&(bars->sl_phase_5).cv, &(bars->sl_phase_5).mutex);
#line 703
			if (Error != 0) {
#line 703
				break;
#line 703
			}
#line 703
		}
#line 703
		pthread_setcancelstate(Cancel, &Temp);
#line 703
	} else {
#line 703
		(bars->sl_phase_5).cycle = !(bars->sl_phase_5).cycle;
#line 703
		(bars->sl_phase_5).counter = 0;
#line 703
		Error = pthread_cond_broadcast(&(bars->sl_phase_5).cv);
#line 703
	}
#line 703
	pthread_mutex_unlock(&(bars->sl_phase_5).mutex);
#line 703
}
#else
   {
#line 705
	unsigned long	Error, Cycle;
#line 705
	long		Cancel, Temp;
#line 705

#line 705
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 705
	if (Error != 0) {
#line 705
		printf("Error while trying to get lock in barrier.\n");
#line 705
		exit(-1);
#line 705
	}
#line 705

#line 705
	Cycle = (bars->barrier).cycle;
#line 705
	if (++(bars->barrier).counter != (nprocs)) {
#line 705
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 705
		while (Cycle == (bars->barrier).cycle) {
#line 705
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 705
			if (Error != 0) {
#line 705
				break;
#line 705
			}
#line 705
		}
#line 705
		pthread_setcancelstate(Cancel, &Temp);
#line 705
	} else {
#line 705
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 705
		(bars->barrier).counter = 0;
#line 705
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 705
	}
#line 705
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 705
}
#endif
/*     *******************************************************

               s i x t h   p h a s e

       *******************************************************  */

   istart = 1;
   iend = istart + gp[procid].rel_num_y[numlev-1] - 1;
   jstart = 1;
   jend = jstart + gp[procid].rel_num_x[numlev-1] - 1;
   ist = istart;
   ien = iend;
   jst = jstart;
   jen = jend;

   if (gp[procid].neighbors[UP] == -1) {
     istart = 0;
   }
   if (gp[procid].neighbors[LEFT] == -1) {
     jstart = 0;
   }
   if (gp[procid].neighbors[DOWN] == -1) {
     iend = im-1;
   }
   if (gp[procid].neighbors[RIGHT] == -1) {
     jend = jm-1;
   }
   t2a = (double **) rhs_multi[procid][numlev-1];
   t2b = (double **) ga[procid];
   t2c = (double **) oldga[procid];
   t2d = (double **) q_multi[procid][numlev-1];
   for(i=istart;i<=iend;i++) {
     t1a = (double *) t2a[i];
     t1b = (double *) t2b[i];
     for(j=jstart;j<=jend;j++) {
       t1a[j] = t1b[j] * ressqr;
     }
   }

   if (gp[procid].neighbors[UP] == -1) {
     t1d = (double *) t2d[0];
     t1b = (double *) t2b[0];
     for(j=jstart;j<=jend;j++) {
       t1d[j] = t1b[j];
     }
   }
   if (gp[procid].neighbors[DOWN] == -1) {
     t1d = (double *) t2d[im-1];
     t1b = (double *) t2b[im-1];
     for(j=jstart;j<=jend;j++) {
       t1d[j] = t1b[j];
     }
   }
   if (gp[procid].neighbors[LEFT] == -1) {
     for(i=istart;i<=iend;i++) {
       t2d[i][0] = t2b[i][0];
     }
   }
   if (gp[procid].neighbors[RIGHT] == -1) {
     for(i=istart;i<=iend;i++) {
       t2d[i][jm-1] = t2b[i][jm-1];
     }
   }

   fac = 1.0 / (4.0 - ressqr*eig2);
   for(i=ist;i<=ien;i++) {
     t1d = (double *) t2d[i];
     t1c = (double *) t2c[i];
     for(j=jst;j<=jen;j++) {
       t1d[j] = t1c[j];
     }
   }

   if ((procid == MASTER) || (do_stats)) {
     {
#line 781
	struct timeval	FullTime;
#line 781

#line 781
	gettimeofday(&FullTime, NULL);
#line 781
	(multi_start) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 781
};
   }

   multig(procid);

   if ((procid == MASTER) || (do_stats)) {
     {
#line 787
	struct timeval	FullTime;
#line 787

#line 787
	gettimeofday(&FullTime, NULL);
#line 787
	(multi_end) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 787
};
     gp[procid].multi_time += (multi_end - multi_start);
   }

/* the shared sum variable psiai is initialized to 0 at
   every time-step  */

   if (procid == MASTER) {
     global->psiai=0.0;
   }

/*  copy the solution for use as initial guess in next time-step  */

   for(i=istart;i<=iend;i++) {
     t1b = (double *) t2b[i];
     t1c = (double *) t2c[i];
     t1d = (double *) t2d[i];
     for(j=jstart;j<=jend;j++) {
       t1b[j] = t1d[j];
       t1c[j] = t1d[j];
     }
   }
#if defined(MULTIPLE_BARRIERS)
   {
#line 810
	unsigned long	Error, Cycle;
#line 810
	long		Cancel, Temp;
#line 810

#line 810
	Error = pthread_mutex_lock(&(bars->sl_phase_6).mutex);
#line 810
	if (Error != 0) {
#line 810
		printf("Error while trying to get lock in barrier.\n");
#line 810
		exit(-1);
#line 810
	}
#line 810

#line 810
	Cycle = (bars->sl_phase_6).cycle;
#line 810
	if (++(bars->sl_phase_6).counter != (nprocs)) {
#line 810
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 810
		while (Cycle == (bars->sl_phase_6).cycle) {
#line 810
			Error = pthread_cond_wait(&(bars->sl_phase_6).cv, &(bars->sl_phase_6).mutex);
#line 810
			if (Error != 0) {
#line 810
				break;
#line 810
			}
#line 810
		}
#line 810
		pthread_setcancelstate(Cancel, &Temp);
#line 810
	} else {
#line 810
		(bars->sl_phase_6).cycle = !(bars->sl_phase_6).cycle;
#line 810
		(bars->sl_phase_6).counter = 0;
#line 810
		Error = pthread_cond_broadcast(&(bars->sl_phase_6).cv);
#line 810
	}
#line 810
	pthread_mutex_unlock(&(bars->sl_phase_6).mutex);
#line 810
}
#else
   {
#line 812
	unsigned long	Error, Cycle;
#line 812
	long		Cancel, Temp;
#line 812

#line 812
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 812
	if (Error != 0) {
#line 812
		printf("Error while trying to get lock in barrier.\n");
#line 812
		exit(-1);
#line 812
	}
#line 812

#line 812
	Cycle = (bars->barrier).cycle;
#line 812
	if (++(bars->barrier).counter != (nprocs)) {
#line 812
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 812
		while (Cycle == (bars->barrier).cycle) {
#line 812
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 812
			if (Error != 0) {
#line 812
				break;
#line 812
			}
#line 812
		}
#line 812
		pthread_setcancelstate(Cancel, &Temp);
#line 812
	} else {
#line 812
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 812
		(bars->barrier).counter = 0;
#line 812
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 812
	}
#line 812
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 812
}
#endif
/*     *******************************************************

                s e v e n t h   p h a s e

       *******************************************************

   every process computes the running sum for its assigned portion
   in a private variable psiaipriv   */

   psiaipriv=0.0;
   t2a = (double **) ga[procid];
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     psiaipriv = psiaipriv + 0.25*(t2a[0][0]);
   }
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     psiaipriv = psiaipriv + 0.25*(t2a[0][jm-1]);
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     psiaipriv=psiaipriv+0.25*(t2a[im-1][0]);
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     psiaipriv=psiaipriv+0.25*(t2a[im-1][jm-1]);
   }
   if (gp[procid].neighbors[UP] == -1) {
     t1a = (double *) t2a[0];
     for(j=firstcol;j<=lastcol;j++) {
       psiaipriv = psiaipriv + 0.5*t1a[j];
     }
   }
   if (gp[procid].neighbors[DOWN] == -1) {
     t1a = (double *) t2a[im-1];
     for(j=firstcol;j<=lastcol;j++) {
       psiaipriv = psiaipriv + 0.5*t1a[j];
     }
   }
   if (gp[procid].neighbors[LEFT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       psiaipriv = psiaipriv + 0.5*t2a[j][0];
     }
   }
   if (gp[procid].neighbors[RIGHT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       psiaipriv = psiaipriv + 0.5*t2a[j][jm-1];
     }
   }
   for(i=firstrow;i<=lastrow;i++) {
     t1a = (double *) t2a[i];
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
       psiaipriv = psiaipriv + t1a[iindex];
     }
   }

/* after computing its private sum, every process adds that to the
   shared running sum psiai  */

   {pthread_mutex_lock(&(locks->psiailock));}
   global->psiai = global->psiai + psiaipriv;
   {pthread_mutex_unlock(&(locks->psiailock));}
#if defined(MULTIPLE_BARRIERS)
   {
#line 873
	unsigned long	Error, Cycle;
#line 873
	long		Cancel, Temp;
#line 873

#line 873
	Error = pthread_mutex_lock(&(bars->sl_phase_7).mutex);
#line 873
	if (Error != 0) {
#line 873
		printf("Error while trying to get lock in barrier.\n");
#line 873
		exit(-1);
#line 873
	}
#line 873

#line 873
	Cycle = (bars->sl_phase_7).cycle;
#line 873
	if (++(bars->sl_phase_7).counter != (nprocs)) {
#line 873
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 873
		while (Cycle == (bars->sl_phase_7).cycle) {
#line 873
			Error = pthread_cond_wait(&(bars->sl_phase_7).cv, &(bars->sl_phase_7).mutex);
#line 873
			if (Error != 0) {
#line 873
				break;
#line 873
			}
#line 873
		}
#line 873
		pthread_setcancelstate(Cancel, &Temp);
#line 873
	} else {
#line 873
		(bars->sl_phase_7).cycle = !(bars->sl_phase_7).cycle;
#line 873
		(bars->sl_phase_7).counter = 0;
#line 873
		Error = pthread_cond_broadcast(&(bars->sl_phase_7).cv);
#line 873
	}
#line 873
	pthread_mutex_unlock(&(bars->sl_phase_7).mutex);
#line 873
}
#else
   {
#line 875
	unsigned long	Error, Cycle;
#line 875
	long		Cancel, Temp;
#line 875

#line 875
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 875
	if (Error != 0) {
#line 875
		printf("Error while trying to get lock in barrier.\n");
#line 875
		exit(-1);
#line 875
	}
#line 875

#line 875
	Cycle = (bars->barrier).cycle;
#line 875
	if (++(bars->barrier).counter != (nprocs)) {
#line 875
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 875
		while (Cycle == (bars->barrier).cycle) {
#line 875
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 875
			if (Error != 0) {
#line 875
				break;
#line 875
			}
#line 875
		}
#line 875
		pthread_setcancelstate(Cancel, &Temp);
#line 875
	} else {
#line 875
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 875
		(bars->barrier).counter = 0;
#line 875
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 875
	}
#line 875
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 875
}
#endif
/*      *******************************************************

                e i g h t h   p h a s e

        *******************************************************

   augment ga(i,j) with [-psiai/psibi]*psib(i,j) */

   f4 = (-global->psiai)/(global->psibi);

   t2a = (double **) ga[procid];
   t2b = (double **) psib[procid];
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[0][0] = t2a[0][0]+f4*t2b[0][0];
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[im-1][0] = t2a[im-1][0]+f4*t2b[im-1][0];
   }
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[0][jm-1] = t2a[0][jm-1]+f4*t2b[0][jm-1];
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[im-1][jm-1] = t2a[im-1][jm-1] +
			      f4*t2b[im-1][jm-1];
   }
   if (gp[procid].neighbors[UP] == -1) {
     t1a = (double *) t2a[0];
     t1b = (double *) t2b[0];
     for(j=firstcol;j<=lastcol;j++) {
       t1a[j] = t1a[j]+f4*t1b[j];
     }
   }
   if (gp[procid].neighbors[DOWN] == -1) {
     t1a = (double *) t2a[im-1];
     t1b = (double *) t2b[im-1];
     for(j=firstcol;j<=lastcol;j++) {
       t1a[j] = t1a[j]+f4*t1b[j];
     }
   }
   if (gp[procid].neighbors[LEFT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][0] = t2a[j][0]+f4*t2b[j][0];
     }
   }
   if (gp[procid].neighbors[RIGHT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][jm-1] = t2a[j][jm-1]+f4*t2b[j][jm-1];
     }
   }
   for(i=firstrow;i<=lastrow;i++) {
     t1a = (double *) t2a[i];
     t1b = (double *) t2b[i];
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
       t1a[iindex] = t1a[iindex]+f4*t1b[iindex];
     }
   }

   t2a = (double **) rhs_multi[procid][numlev-1];
   t2b = (double **) gb[procid];
   t2c = (double **) oldgb[procid];
   t2d = (double **) q_multi[procid][numlev-1];
   for(i=istart;i<=iend;i++) {
     t1a = (double *) t2a[i];
     t1b = (double *) t2b[i];
     for(j=jstart;j<=jend;j++) {
       t1a[j] = t1b[j] * ressqr;
     }
   }
   if (gp[procid].neighbors[UP] == -1) {
     t1d = (double *) t2d[0];
     t1b = (double *) t2b[0];
     for(j=jstart;j<=jend;j++) {
       t1d[j] = t1b[j];
     }
   }
   if (gp[procid].neighbors[DOWN] == -1) {
     t1d = (double *) t2d[im-1];
     t1b = (double *) t2b[im-1];
     for(j=jstart;j<=jend;j++) {
       t1d[j] = t1b[j];
     }
   }
   if (gp[procid].neighbors[LEFT] == -1) {
     for(i=istart;i<=iend;i++) {
       t2d[i][0] = t2b[i][0];
     }
   }
   if (gp[procid].neighbors[RIGHT] == -1) {
     for(i=istart;i<=iend;i++) {
       t2d[i][jm-1] = t2b[i][jm-1];
     }
   }

   fac = 1.0 / (4.0 - ressqr*eig2);
   for(i=ist;i<=ien;i++) {
     t1d = (double *) t2d[i];
     t1c = (double *) t2c[i];
     for(j=jst;j<=jen;j++) {
       t1d[j] = t1c[j];
     }
   }

   if ((procid == MASTER) || (do_stats)) {
     {
#line 980
	struct timeval	FullTime;
#line 980

#line 980
	gettimeofday(&FullTime, NULL);
#line 980
	(multi_start) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 980
};
   }

   multig(procid);

   if ((procid == MASTER) || (do_stats)) {
     {
#line 986
	struct timeval	FullTime;
#line 986

#line 986
	gettimeofday(&FullTime, NULL);
#line 986
	(multi_end) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 986
};
     gp[procid].multi_time += (multi_end - multi_start);
   }

   for(i=istart;i<=iend;i++) {
     t1b = (double *) t2b[i];
     t1c = (double *) t2c[i];
     t1d = (double *) t2d[i];
     for(j=jstart;j<=jend;j++) {
       t1b[j] = t1d[j];
       t1c[j] = t1d[j];
     }
   }
#if defined(MULTIPLE_BARRIERS)
   {
#line 1000
	unsigned long	Error, Cycle;
#line 1000
	long		Cancel, Temp;
#line 1000

#line 1000
	Error = pthread_mutex_lock(&(bars->sl_phase_8).mutex);
#line 1000
	if (Error != 0) {
#line 1000
		printf("Error while trying to get lock in barrier.\n");
#line 1000
		exit(-1);
#line 1000
	}
#line 1000

#line 1000
	Cycle = (bars->sl_phase_8).cycle;
#line 1000
	if (++(bars->sl_phase_8).counter != (nprocs)) {
#line 1000
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 1000
		while (Cycle == (bars->sl_phase_8).cycle) {
#line 1000
			Error = pthread_cond_wait(&(bars->sl_phase_8).cv, &(bars->sl_phase_8).mutex);
#line 1000
			if (Error != 0) {
#line 1000
				break;
#line 1000
			}
#line 1000
		}
#line 1000
		pthread_setcancelstate(Cancel, &Temp);
#line 1000
	} else {
#line 1000
		(bars->sl_phase_8).cycle = !(bars->sl_phase_8).cycle;
#line 1000
		(bars->sl_phase_8).counter = 0;
#line 1000
		Error = pthread_cond_broadcast(&(bars->sl_phase_8).cv);
#line 1000
	}
#line 1000
	pthread_mutex_unlock(&(bars->sl_phase_8).mutex);
#line 1000
}
#else
   {
#line 1002
	unsigned long	Error, Cycle;
#line 1002
	long		Cancel, Temp;
#line 1002

#line 1002
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 1002
	if (Error != 0) {
#line 1002
		printf("Error while trying to get lock in barrier.\n");
#line 1002
		exit(-1);
#line 1002
	}
#line 1002

#line 1002
	Cycle = (bars->barrier).cycle;
#line 1002
	if (++(bars->barrier).counter != (nprocs)) {
#line 1002
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 1002
		while (Cycle == (bars->barrier).cycle) {
#line 1002
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 1002
			if (Error != 0) {
#line 1002
				break;
#line 1002
			}
#line 1002
		}
#line 1002
		pthread_setcancelstate(Cancel, &Temp);
#line 1002
	} else {
#line 1002
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 1002
		(bars->barrier).counter = 0;
#line 1002
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 1002
	}
#line 1002
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 1002
}
#endif
/*      *******************************************************

                n i n t h   p h a s e

        *******************************************************

   put appropriate linear combinations of ga and gb in work2 and work3;
   note that here (as in most cases) the constant multipliers are made
   private variables; the specific order in which things are done is
   chosen in order to hopefully reuse things brought into the cache

   note that here again we choose to have all processes share the work
   on both matrices despite the fact that the work done per element
   is the same, because the operand matrices are the same in both cases */

   t2a = (double **) ga[procid];
   t2b = (double **) gb[procid];
   t2c = (double **) work2[procid];
   t2d = (double **) work3[procid];
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2c[0][0] = t2b[0][0]-hh1*t2a[0][0];
     t2d[0][0] = t2b[0][0]+hh3*t2a[0][0];
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2c[im-1][0] = t2b[im-1][0]-hh1*t2a[im-1][0];
     t2d[im-1][0] = t2b[im-1][0]+hh3*t2a[im-1][0];
   }
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2c[0][jm-1] = t2b[0][jm-1]-hh1*t2a[0][jm-1];
     t2d[0][jm-1] = t2b[0][jm-1]+hh3*t2a[0][jm-1];
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2c[im-1][jm-1] = t2b[im-1][jm-1] -
				 hh1*t2a[im-1][jm-1];
     t2d[im-1][jm-1] = t2b[im-1][jm-1] +
				 hh3*t2a[im-1][jm-1];
   }
   if (gp[procid].neighbors[UP] == -1) {
     t1a = (double *) t2a[0];
     t1b = (double *) t2b[0];
     t1c = (double *) t2c[0];
     t1d = (double *) t2d[0];
     for(j=firstcol;j<=lastcol;j++) {
       t1d[j] = t1b[j]+hh3*t1a[j];
       t1c[j] = t1b[j]-hh1*t1a[j];
     }
   }
   if (gp[procid].neighbors[DOWN] == -1) {
     t1a = (double *) t2a[im-1];
     t1b = (double *) t2b[im-1];
     t1c = (double *) t2c[im-1];
     t1d = (double *) t2d[im-1];
     for(j=firstcol;j<=lastcol;j++) {
       t1d[j] = t1b[j]+hh3*t1a[j];
       t1c[j] = t1b[j]-hh1*t1a[j];
     }
   }
   if (gp[procid].neighbors[LEFT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2d[j][0] = t2b[j][0]+hh3*t2a[j][0];
       t2c[j][0] = t2b[j][0]-hh1*t2a[j][0];
     }
   }
   if (gp[procid].neighbors[RIGHT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2d[j][jm-1] = t2b[j][jm-1]+hh3*t2a[j][jm-1];
       t2c[j][jm-1] = t2b[j][jm-1]-hh1*t2a[j][jm-1];
     }
   }

   for(i=firstrow;i<=lastrow;i++) {
     t1a = (double *) t2a[i];
     t1b = (double *) t2b[i];
     t1c = (double *) t2c[i];
     t1d = (double *) t2d[i];
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
       t1d[iindex] = t1b[iindex] + hh3*t1a[iindex];
       t1c[iindex] = t1b[iindex] - hh1*t1a[iindex];
     }
   }
#if defined(MULTIPLE_BARRIERS)
   {
#line 1085
	unsigned long	Error, Cycle;
#line 1085
	long		Cancel, Temp;
#line 1085

#line 1085
	Error = pthread_mutex_lock(&(bars->sl_phase_9).mutex);
#line 1085
	if (Error != 0) {
#line 1085
		printf("Error while trying to get lock in barrier.\n");
#line 1085
		exit(-1);
#line 1085
	}
#line 1085

#line 1085
	Cycle = (bars->sl_phase_9).cycle;
#line 1085
	if (++(bars->sl_phase_9).counter != (nprocs)) {
#line 1085
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 1085
		while (Cycle == (bars->sl_phase_9).cycle) {
#line 1085
			Error = pthread_cond_wait(&(bars->sl_phase_9).cv, &(bars->sl_phase_9).mutex);
#line 1085
			if (Error != 0) {
#line 1085
				break;
#line 1085
			}
#line 1085
		}
#line 1085
		pthread_setcancelstate(Cancel, &Temp);
#line 1085
	} else {
#line 1085
		(bars->sl_phase_9).cycle = !(bars->sl_phase_9).cycle;
#line 1085
		(bars->sl_phase_9).counter = 0;
#line 1085
		Error = pthread_cond_broadcast(&(bars->sl_phase_9).cv);
#line 1085
	}
#line 1085
	pthread_mutex_unlock(&(bars->sl_phase_9).mutex);
#line 1085
}
#else
   {
#line 1087
	unsigned long	Error, Cycle;
#line 1087
	long		Cancel, Temp;
#line 1087

#line 1087
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 1087
	if (Error != 0) {
#line 1087
		printf("Error while trying to get lock in barrier.\n");
#line 1087
		exit(-1);
#line 1087
	}
#line 1087

#line 1087
	Cycle = (bars->barrier).cycle;
#line 1087
	if (++(bars->barrier).counter != (nprocs)) {
#line 1087
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 1087
		while (Cycle == (bars->barrier).cycle) {
#line 1087
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 1087
			if (Error != 0) {
#line 1087
				break;
#line 1087
			}
#line 1087
		}
#line 1087
		pthread_setcancelstate(Cancel, &Temp);
#line 1087
	} else {
#line 1087
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 1087
		(bars->barrier).counter = 0;
#line 1087
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 1087
	}
#line 1087
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 1087
}
#endif
/*      *******************************************************

                t e n t h    p h a s e

        *******************************************************/


   timst = 2*dtau;

/* update the psi{1,3} matrices by adding 2*dtau*work3 to each */

   t2a = (double **) psi[procid][0];
   t2b = (double **) work3[procid];
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[0][0] = t2a[0][0] + timst*t2b[0][0];
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[im-1][0] = t2a[im-1][0] +
			       timst*t2b[im-1][0];
   }
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[0][jm-1] = t2a[0][jm-1] +
			       timst*t2b[0][jm-1];
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[im-1][jm-1] = t2a[im-1][jm-1] +
				  timst*t2b[im-1][jm-1];
   }
   if (gp[procid].neighbors[UP] == -1) {
     t1a = (double *) t2a[0];
     t1b = (double *) t2b[0];
     for(j=firstcol;j<=lastcol;j++) {
       t1a[j] = t1a[j] + timst*t1b[j];
     }
   }
   if (gp[procid].neighbors[DOWN] == -1) {
     t1a = (double *) t2a[im-1];
     t1b = (double *) t2b[im-1];
     for(j=firstcol;j<=lastcol;j++) {
       t1a[j] = t1a[j] + timst*t1b[j];
     }
   }
   if (gp[procid].neighbors[LEFT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][0] = t2a[j][0] + timst*t2b[j][0];
     }
   }
   if (gp[procid].neighbors[RIGHT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][jm-1] = t2a[j][jm-1] +
				 timst*t2b[j][jm-1];
     }
   }
   for(i=firstrow;i<=lastrow;i++) {
     t1a = (double *) t2a[i];
     t1b = (double *) t2b[i];
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
         t1a[iindex] = t1a[iindex] + timst*t1b[iindex];
     }
   }

   t2a = (double **) psi[procid][1];
   t2b = (double **) work2[procid];
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[0][0] = t2a[0][0] + timst*t2b[0][0];
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[im-1][0] = t2a[im-1][0] +
			       timst*t2b[im-1][0];
   }
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[0][jm-1] = t2a[0][jm-1] +
			       timst*t2b[0][jm-1];
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[im-1][jm-1] = t2a[im-1][jm-1] +
				  timst*t2b[im-1][jm-1];
   }
   if (gp[procid].neighbors[UP] == -1) {
     t1a = (double *) t2a[0];
     t1b = (double *) t2b[0];
     for(j=firstcol;j<=lastcol;j++) {
       t1a[j] = t1a[j] + timst*t1b[j];
     }
   }
   if (gp[procid].neighbors[DOWN] == -1) {
     t1a = (double *) t2a[im-1];
     t1b = (double *) t2b[im-1];
     for(j=firstcol;j<=lastcol;j++) {
       t1a[j] = t1a[j] + timst*t1b[j];
     }
   }
   if (gp[procid].neighbors[LEFT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][0] = t2a[j][0] + timst*t2b[j][0];
     }
   }
   if (gp[procid].neighbors[RIGHT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][jm-1] = t2a[j][jm-1] +
				 timst*t2b[j][jm-1];
     }
   }

   for(i=firstrow;i<=lastrow;i++) {
     t1a = (double *) t2a[i];
     t1b = (double *) t2b[i];
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
         t1a[iindex] = t1a[iindex] + timst*t1b[iindex];
     }
   }
#if defined(MULTIPLE_BARRIERS)
   {
#line 1201
	unsigned long	Error, Cycle;
#line 1201
	long		Cancel, Temp;
#line 1201

#line 1201
	Error = pthread_mutex_lock(&(bars->sl_phase_10).mutex);
#line 1201
	if (Error != 0) {
#line 1201
		printf("Error while trying to get lock in barrier.\n");
#line 1201
		exit(-1);
#line 1201
	}
#line 1201

#line 1201
	Cycle = (bars->sl_phase_10).cycle;
#line 1201
	if (++(bars->sl_phase_10).counter != (nprocs)) {
#line 1201
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 1201
		while (Cycle == (bars->sl_phase_10).cycle) {
#line 1201
			Error = pthread_cond_wait(&(bars->sl_phase_10).cv, &(bars->sl_phase_10).mutex);
#line 1201
			if (Error != 0) {
#line 1201
				break;
#line 1201
			}
#line 1201
		}
#line 1201
		pthread_setcancelstate(Cancel, &Temp);
#line 1201
	} else {
#line 1201
		(bars->sl_phase_10).cycle = !(bars->sl_phase_10).cycle;
#line 1201
		(bars->sl_phase_10).counter = 0;
#line 1201
		Error = pthread_cond_broadcast(&(bars->sl_phase_10).cv);
#line 1201
	}
#line 1201
	pthread_mutex_unlock(&(bars->sl_phase_10).mutex);
#line 1201
}
#else
   {
#line 1203
	unsigned long	Error, Cycle;
#line 1203
	long		Cancel, Temp;
#line 1203

#line 1203
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 1203
	if (Error != 0) {
#line 1203
		printf("Error while trying to get lock in barrier.\n");
#line 1203
		exit(-1);
#line 1203
	}
#line 1203

#line 1203
	Cycle = (bars->barrier).cycle;
#line 1203
	if (++(bars->barrier).counter != (nprocs)) {
#line 1203
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 1203
		while (Cycle == (bars->barrier).cycle) {
#line 1203
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 1203
			if (Error != 0) {
#line 1203
				break;
#line 1203
			}
#line 1203
		}
#line 1203
		pthread_setcancelstate(Cancel, &Temp);
#line 1203
	} else {
#line 1203
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 1203
		(bars->barrier).counter = 0;
#line 1203
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 1203
	}
#line 1203
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 1203
}
#endif
}
#line 228 "/home/jyy/Dev/splash2/codes/null_macros/c.m4.null.POSIX"

#line 1 "subblock.C"
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
#include <math.h>
#include "decs.h"

void subblock()

{
   long i;
   long j;
   long k;
   long xportion;
   long xextra;
   long yportion;
   long yextra;
   long my_num;

/* Determine starting coord and number of points to process in     */
/* each direction                                                  */

   for (i=0;i<numlev;i++) {
     xportion = (jmx[i] - 2) / xprocs;
     xextra = (jmx[i] - 2) % xprocs;
     for (j=0;j<xprocs;j++) {
       for (k=0;k<yprocs;k++) {
         gp[k*xprocs+j].rel_num_x[i] = xportion;
       }
     }
     yportion = (imx[i] - 2) / yprocs;
     yextra = (imx[i] - 2) % yprocs;
     for (j=0;j<yprocs;j++) {
       for (k=0;k<xprocs;k++) {
         gp[j*xprocs+k].rel_num_y[i] = yportion;
       }
     }
   }

   for (my_num=0;my_num<nprocs;my_num++) {
     for (i=0;i<numlev;i++) {
       gp[my_num].rlist[i] = 1;
       gp[my_num].rljst[i] = 1;
       gp[my_num].rlien[i] = gp[my_num].rlist[i] + gp[my_num].rel_num_y[i];
       gp[my_num].rljen[i] = gp[my_num].rljst[i] + gp[my_num].rel_num_x[i];
       gp[my_num].eist[i] = gp[my_num].rlist[i] + 1;
       gp[my_num].oist[i] = gp[my_num].rlist[i];
       gp[my_num].ejst[i] = gp[my_num].rljst[i] + 1;
       gp[my_num].ojst[i] = gp[my_num].rljst[i];
     }
   }
  for (i=0;i<nprocs;i++) {
    gp[i].neighbors[LEFT] = -1;
    gp[i].neighbors[RIGHT] = -1;
    gp[i].neighbors[UP] = -1;
    gp[i].neighbors[DOWN] = -1;
    gp[i].neighbors[UPLEFT] = -1;
    gp[i].neighbors[UPRIGHT] = -1;
    gp[i].neighbors[DOWNLEFT] = -1;
    gp[i].neighbors[DOWNRIGHT] = -1;
    if (i >= xprocs) {
      gp[i].neighbors[UP] = i-xprocs;
    }
    if (i < nprocs-xprocs) {
      gp[i].neighbors[DOWN] = i+xprocs;
    }
    if ((i % xprocs) > 0) {
      gp[i].neighbors[LEFT] = i-1;
    }
    if ((i % xprocs) < (xprocs-1)) {
      gp[i].neighbors[RIGHT] = i+1;
    }
    j = gp[i].neighbors[UP];
    if (j != -1) {
      if ((j % xprocs) > 0) {
        gp[i].neighbors[UPLEFT] = j-1;
      }
      if ((j % xprocs) < (xprocs-1)) {
        gp[i].neighbors[UPRIGHT] = j+1;
      }
    }
    j = gp[i].neighbors[DOWN];
    if (j != -1) {
      if ((j % xprocs) > 0) {
        gp[i].neighbors[DOWNLEFT] = j-1;
      }
      if ((j % xprocs) < (xprocs-1)) {
        gp[i].neighbors[DOWNRIGHT] = j+1;
      }
    }
  }
  for (i=0;i<nprocs;i++) {
    gp[i].rownum = i/xprocs;
    gp[i].colnum = i%xprocs;
  }
}

