#ifndef __DEC_H__
#define __DEC_H__

#line 228 "/home/jyy/Dev/splash2/codes/null_macros/c.m4.null.POSIX"

#line 1 "decs.H"
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

#define MASTER            0
#define RED_ITER          0
#define BLACK_ITER        1
#define UP                0
#define DOWN              1
#define LEFT              2
#define RIGHT             3
#define UPLEFT            4
#define UPRIGHT           5
#define DOWNLEFT          6
#define DOWNRIGHT         7
#define PAGE_SIZE      4096

struct multi_struct {
   double err_multi;
};

extern struct multi_struct *multi;

struct global_struct {
   long id;
   long starttime;
   long trackstart;
   double psiai;
   double psibi;
};

extern struct global_struct *global;

extern double eig2;
extern double ysca;
extern long jmm1;
extern double pi;
extern double t0;

extern double ****psi;
extern double ****psim;
extern double ***psium;
extern double ***psilm;
extern double ***psib;
extern double ***ga;
extern double ***gb;
extern double ****work1;
extern double ***work2;
extern double ***work3;
extern double ****work4;
extern double ****work5;
extern double ***work6;
extern double ****work7;
extern double ****temparray;
extern double ***tauz;
extern double ***oldga;
extern double ***oldgb;
extern double *f;
extern double ****q_multi;
extern double ****rhs_multi;

struct locks_struct {
   pthread_mutex_t (idlock);
   pthread_mutex_t (psiailock);
   pthread_mutex_t (psibilock);
   pthread_mutex_t (donelock);
   pthread_mutex_t (error_lock);
   pthread_mutex_t (bar_lock);
};

extern struct locks_struct *locks;

struct bars_struct {
#if defined(MULTIPLE_BARRIERS)
   
#line 87
struct {
#line 87
	pthread_mutex_t	mutex;
#line 87
	pthread_cond_t	cv;
#line 87
	unsigned long	counter;
#line 87
	unsigned long	cycle;
#line 87
} (iteration);
#line 87

   
#line 88
struct {
#line 88
	pthread_mutex_t	mutex;
#line 88
	pthread_cond_t	cv;
#line 88
	unsigned long	counter;
#line 88
	unsigned long	cycle;
#line 88
} (gsudn);
#line 88

   
#line 89
struct {
#line 89
	pthread_mutex_t	mutex;
#line 89
	pthread_cond_t	cv;
#line 89
	unsigned long	counter;
#line 89
	unsigned long	cycle;
#line 89
} (p_setup);
#line 89

   
#line 90
struct {
#line 90
	pthread_mutex_t	mutex;
#line 90
	pthread_cond_t	cv;
#line 90
	unsigned long	counter;
#line 90
	unsigned long	cycle;
#line 90
} (p_redph);
#line 90

   
#line 91
struct {
#line 91
	pthread_mutex_t	mutex;
#line 91
	pthread_cond_t	cv;
#line 91
	unsigned long	counter;
#line 91
	unsigned long	cycle;
#line 91
} (p_soln);
#line 91

   
#line 92
struct {
#line 92
	pthread_mutex_t	mutex;
#line 92
	pthread_cond_t	cv;
#line 92
	unsigned long	counter;
#line 92
	unsigned long	cycle;
#line 92
} (p_subph);
#line 92

   
#line 93
struct {
#line 93
	pthread_mutex_t	mutex;
#line 93
	pthread_cond_t	cv;
#line 93
	unsigned long	counter;
#line 93
	unsigned long	cycle;
#line 93
} (sl_prini);
#line 93

   
#line 94
struct {
#line 94
	pthread_mutex_t	mutex;
#line 94
	pthread_cond_t	cv;
#line 94
	unsigned long	counter;
#line 94
	unsigned long	cycle;
#line 94
} (sl_psini);
#line 94

   
#line 95
struct {
#line 95
	pthread_mutex_t	mutex;
#line 95
	pthread_cond_t	cv;
#line 95
	unsigned long	counter;
#line 95
	unsigned long	cycle;
#line 95
} (sl_onetime);
#line 95

   
#line 96
struct {
#line 96
	pthread_mutex_t	mutex;
#line 96
	pthread_cond_t	cv;
#line 96
	unsigned long	counter;
#line 96
	unsigned long	cycle;
#line 96
} (sl_phase_1);
#line 96

   
#line 97
struct {
#line 97
	pthread_mutex_t	mutex;
#line 97
	pthread_cond_t	cv;
#line 97
	unsigned long	counter;
#line 97
	unsigned long	cycle;
#line 97
} (sl_phase_2);
#line 97

   
#line 98
struct {
#line 98
	pthread_mutex_t	mutex;
#line 98
	pthread_cond_t	cv;
#line 98
	unsigned long	counter;
#line 98
	unsigned long	cycle;
#line 98
} (sl_phase_3);
#line 98

   
#line 99
struct {
#line 99
	pthread_mutex_t	mutex;
#line 99
	pthread_cond_t	cv;
#line 99
	unsigned long	counter;
#line 99
	unsigned long	cycle;
#line 99
} (sl_phase_4);
#line 99

   
#line 100
struct {
#line 100
	pthread_mutex_t	mutex;
#line 100
	pthread_cond_t	cv;
#line 100
	unsigned long	counter;
#line 100
	unsigned long	cycle;
#line 100
} (sl_phase_5);
#line 100

   
#line 101
struct {
#line 101
	pthread_mutex_t	mutex;
#line 101
	pthread_cond_t	cv;
#line 101
	unsigned long	counter;
#line 101
	unsigned long	cycle;
#line 101
} (sl_phase_6);
#line 101

   
#line 102
struct {
#line 102
	pthread_mutex_t	mutex;
#line 102
	pthread_cond_t	cv;
#line 102
	unsigned long	counter;
#line 102
	unsigned long	cycle;
#line 102
} (sl_phase_7);
#line 102

   
#line 103
struct {
#line 103
	pthread_mutex_t	mutex;
#line 103
	pthread_cond_t	cv;
#line 103
	unsigned long	counter;
#line 103
	unsigned long	cycle;
#line 103
} (sl_phase_8);
#line 103

   
#line 104
struct {
#line 104
	pthread_mutex_t	mutex;
#line 104
	pthread_cond_t	cv;
#line 104
	unsigned long	counter;
#line 104
	unsigned long	cycle;
#line 104
} (sl_phase_9);
#line 104

   
#line 105
struct {
#line 105
	pthread_mutex_t	mutex;
#line 105
	pthread_cond_t	cv;
#line 105
	unsigned long	counter;
#line 105
	unsigned long	cycle;
#line 105
} (sl_phase_10);
#line 105

   
#line 106
struct {
#line 106
	pthread_mutex_t	mutex;
#line 106
	pthread_cond_t	cv;
#line 106
	unsigned long	counter;
#line 106
	unsigned long	cycle;
#line 106
} (error_barrier);
#line 106

#else
   
#line 108
struct {
#line 108
	pthread_mutex_t	mutex;
#line 108
	pthread_cond_t	cv;
#line 108
	unsigned long	counter;
#line 108
	unsigned long	cycle;
#line 108
} (barrier);
#line 108

#endif
};

extern struct bars_struct *bars;

extern double factjacob;
extern double factlap;

struct Global_Private {
  char pad[PAGE_SIZE];
  long *rel_num_x;
  long *rel_num_y;
  long *eist;
  long *ejst;
  long *oist;
  long *ojst;
  long *rlist;
  long *rljst;
  long *rlien;
  long *rljen;
  long rownum;
  long colnum;
  long neighbors[8];
  double multi_time;
  double total_time;
};

extern struct Global_Private *gp;

extern double *i_int_coeff;
extern double *j_int_coeff;
extern long xprocs;
extern long yprocs;

extern long numlev;
extern long *imx;
extern long *jmx;
extern double *lev_res;
extern double *lev_tol;
extern double maxwork;
extern long *xpts_per_proc;
extern long *ypts_per_proc;
extern long minlevel;
extern double outday0;
extern double outday1;
extern double outday2;
extern double outday3;

extern long nprocs;
extern double h1;
extern double h3;
extern double h;
extern double lf;
extern double res;
extern double dtau;
extern double f0;
extern double beta;
extern double gpr;
extern long im;
extern long jm;
extern long do_stats;
extern long do_output;
extern long *multi_times;
extern long *total_times;

/*
 * jacobcalc.C
 */
void jacobcalc(double ***x, double ***y, double ***z, long pid, long firstrow, long lastrow, long firstcol, long lastcol);

/*
 * jacobcalc2.C
 */
void jacobcalc2(double ****x, double ****y, double ****z, long psiindex, long pid, long firstrow, long lastrow, long firstcol, long lastcol);

/*
 * laplacalc.C
 */
void laplacalc(long procid, double ****x, double ****z, long psiindex, long firstrow, long lastrow, long firstcol, long lastcol);

/*
 * linkup.C
 */
void link_all(void);
void linkup(double **row_ptr);
void link_multi(void);

/*
 * main.C
 */
long log_2(long number);
void printerr(char *s);

/*
 * multi.C
 */
void multig(long my_id);
void relax(long k, double *err, long color, long my_num);
void rescal(long kf, long my_num);
void intadd(long kc, long my_num);
void putz(long k, long my_num);
void copy_borders(long k, long pid);
void copy_rhs_borders(long k, long procid);
void copy_red(long k, long procid);
void copy_black(long k, long procid);

/*
 * slave1.C
 */
void slave(void);

/*
 * slave2.C
 */
void slave2(long procid, long firstrow, long lastrow, long numrows, long firstcol, long lastcol, long numcols);

/*
 * subblock.C
 */
void subblock(void);

#endif
