#ifndef _PS_
#define _PS_

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>
#include "../Parameter_files/COSMOLOGY.H"
#include "../Parameter_files/INIT_PARAMS.H"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#ifdef INHOMO_FEEDBACK
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#endif
#include "cosmo_progs.c"
#include "misc.c"

/* New in v1.1 */
#define ERFC_NPTS (int) 75
#define ERFC_PARAM_DELTA (float) 0.1
static double log_erfc_table[ERFC_NPTS], erfc_params[ERFC_NPTS];
static gsl_interp_accel *erfc_acc;
static gsl_spline *erfc_spline;

#define NR_END 1
#define FREE_ARG char*

#define MM 7
#define NSTACK 50

#define FUNC(x,y,z,xx,yy) ((*func)(x,y,z,xx,yy))
#define FUNC2(x1,x2,x3,x4,x5,x6) ((*func)(x1,x2,x3,x4,x5,x6))
#define EPS2 3.0e-11

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

#define NGaussLegendre 40  //defines the number of points in the Gauss-Legendre quadrature integration

#define SPLINE_NPTS (int) 250
#define NGLhigh 100
#define NGLlow 100

#define Nhigh 200
#define Nlow 100
#define NMass 2000

/* New in v2 - part 1 of 4 start */
#define NSFR_high 200
#define NSFR_low 250
#define NGL_SFR 100 
#ifdef INHOMO_FEEDBACK
#define NMTURN 100
#endif //INHOMO_FEEDBACK
/* Number of interpolation points for the interpolation table for z'' 
                This is the same parameter in 21CMMC */
#define zpp_interp_points (int) (300) 
/* New in v2 - part 1 of 4 end */

static gsl_interp_accel *Fcoll_spline_acc;
static gsl_spline *Fcoll_spline;

/* New in v2 - part 2 of 4: begin */
static double log10_overdense_spline_SFR[NSFR_low];
static gsl_interp_accel *NionLow_spline_acc;
#ifdef INHOMO_FEEDBACK
static gsl_interp_accel *NionLow_spline_accm;
static double log10_Mturn_spline_SFR[NMTURN];
static double log10_Mturn_spline_SFRm[NMTURN];
static float log10_Mturn_spline_SFR_float[NMTURN];
static float log10_Mturn_spline_SFRm_float[NMTURN];
static gsl_interp_accel *NionLow_spline_acc_Mturn;
static gsl_interp_accel *NionLow_spline_accm_Mturn;
static gsl_spline2d *NionLow_spline;
static gsl_spline2d *NionLow_splinem;
static double log10_Nion_spline[NSFR_low*NMTURN];
static double log10_Nion_splinem[NSFR_low*NMTURN];
#else //INHOMO_FEEDBACK
static gsl_spline *NionLow_spline;
static double log10_Nion_spline[NSFR_low];
#ifdef MINI_HALO
static double log10_Nion_splinem[NSFR_low];
static gsl_interp_accel *NionLow_spline_accm;
static gsl_spline *NionLow_splinem;
#endif
#endif //INHOMO_FEEDBACK

void initialiseGL_Nion(int n, float M_Min, float M_Max);
#ifdef INHOMO_FEEDBACK
void Nion_Spline_density(float Overdensity, float log10_M_MINa, float *splined_value);
void initialise_Nion_spline(float z, float Mmax, float Mmin, float MMinTurnover, float Alpha_star, float Alpha_esc, float Fstar10, float Fesc10, float Mlim_Fstar, float Mlim_Fesc);
#else //INHOMO_FEEDBACK
void Nion_Spline_density(float Overdensity, float *splined_value);
void initialise_Nion_spline(float z, float Mmax, float Mmin, float MassTurnover, float Alpha_star, float Alpha_esc, float Fstar10, float Fesc10, float Mlim_Fstar, float Mlim_Fesc);
#endif //INHOMO_FEEDBACK
float Mass_limit (float logM, float PL, float FRAC);
void bisection(float *x, float xlow, float xup, int *iter);
float Mass_limit_bisection(float Mmin, float Mmax, float PL, float FRAC);

static double Nion_z_val[zpp_interp_points]; // For Ts.c
static double SFRD_val[zpp_interp_points]; 
//static double z_val[zpp_interp_points],Nion_z_val[zpp_interp_points]; // For Ts.c
//static double z_X_val[zpp_interp_points],SFRD_val[zpp_interp_points]; 
static gsl_interp_accel *Nion_z_spline_acc;
static gsl_spline *Nion_z_spline;
static gsl_interp_accel *SFRD_ST_z_spline_acc;
static gsl_spline *SFRD_ST_z_spline;
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
static double Nion_z_valm[zpp_interp_points*NMTURN],SFRD_valm[zpp_interp_points*NMTURN]; // For Ts.c
#else
static double Nion_z_valm[zpp_interp_points],SFRD_valm[zpp_interp_points]; // For Ts.c
#endif
static gsl_interp_accel *Nion_z_spline_accm;
#ifdef INHOMO_FEEDBACK
static gsl_interp_accel *Nion_z_spline_accm_Mturn;
static gsl_spline2d *Nion_z_splinem;
#else
static gsl_spline *Nion_z_splinem;
#endif
static gsl_interp_accel *SFRD_ST_z_spline_accm;
#ifdef INHOMO_FEEDBACK
static gsl_interp_accel *SFRD_ST_z_spline_accm_Mturn;
static gsl_spline2d *SFRD_ST_z_splinem;
#else
static gsl_spline *SFRD_ST_z_splinem;
#endif
#endif
void initialise_Nion_ST_spline(int Nbin, double z_val[], float Mmin, double M_MINa_interp_table[], float Alpha_star, float Alpha_esc, float Fstar10, float Fesc10);
void Nion_ST_z(float z, float *splined_value);
void initialise_SFRD_ST_spline(int Nbin, double z_val[], float Mmin, double M_MINa_interp_table[], float Alpha_star, float Fstar10);
void SFRD_ST_z(float z, float *splined_value);
float *zpp_table;
double *log10_overdense_low_table;
float *Overdense_high_table;
#ifdef INHOMO_FEEDBACK
double *log10_overdense_low_table_Mturn;
float *Overdense_high_table_Mturn;
#endif
float **SFRD_z_high_table;
double **log10_SFRD_z_low_table;
#ifdef MINI_HALO
double **log10_SFRD_z_low_tablem;
float **SFRD_z_high_tablem;
#endif

#ifdef INHOMO_FEEDBACK
void initialise_SFRD_Conditional_table(int Nsteps_zp, int Nfilter, float z[], double R[], float Mmin, float Alpha_star, float Fstar10);
#else
#ifdef REION_SM
void initialise_SFRD_Conditional_table(int Nsteps_zp, int Nfilter, float z[], double R[], float Mmin, float Alpha_star, float Fstar10, double REION_SM13_Z_RE, double REION_SM13_DELTA_Z_RE, double REION_SM13_DELTA_Z_SC);
#else
void initialise_SFRD_Conditional_table(int Nsteps_zp, int Nfilter, float z[], double R[], float Mmin, float MassTurnover, float Alpha_star, float Fstar10);
#endif
#endif
void initialise_Xray_Fcollz_SFR_Conditional(int R_ct, int zp_int1, int zp_int2);


struct parameters_gsl_SFR_int_{
    double z_obs;
#ifdef CONTEMPORANEOUS_DUTYCYCLE
    double zp_obs;
#endif
    double Mdrop;
    double pl_star;
    double pl_esc;
    double frac_star;
    double frac_esc;
    double LimitMass_Fstar;
    double LimitMass_Fesc;
};

#ifdef MINI_HALO
struct parameters_gsl_SFR_intm_{
    double z_obs;
#ifdef CONTEMPORANEOUS_DUTYCYCLE
    double zp_obs;
#endif
    double Mdrop;
    double Mdrop_2;
    double pl_star;
    double frac_star;
    double LimitMass_Fstar;
};
#endif

struct parameters_gsl_SFR_con_int_{
    double z_obs;
    double Mval;
    double delta1;
    double delta2;
    double Mdrop;
    double pl_star;
    double pl_esc;
    double frac_star;
    double frac_esc;
    double LimitMass_Fstar;
    double LimitMass_Fesc;
};
/* New in v2 - part 2 of 4: end */
#ifdef MINI_HALO
struct parameters_gsl_SFR_con_intm_{
    double z_obs;
    double Mval;
    double delta1;
    double delta2;
    double pl_star;
    double Mdropm1;
    double Mdropm2;
    double frac_starm;
    double LimitMass_Fstarm;
};
#endif

unsigned long *lvector(long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);

float *vector(long nl, long nh);
void free_vector(float *v, long nl, long nh);

void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);

void gauleg(float x1, float x2, float x[], float w[], int n);

double FgtrlnM_general(double lnM, void *params);
double FgtrM_general(float z, float M1, float M_Max, float M2, float MFeedback, float alpha, float delta1, float delta2);

float FgtrConditionalM_second(float z, float M1, float M2, float MFeedback, float alpha, float delta1, float delta2);
float dNdM_conditional_second(float z, float M1, float M2, float delta1, float delta2);

float *Mass_Spline, *Sigma_Spline, *dSigmadm_Spline, *second_derivs_sigma, *second_derivs_dsigma;

void initialiseSplinedSigmaM(float M_Min, float M_Max);

/* New in v2 - part 3 of 4: start */
float *xi_SFR,*wi_SFR;
float *Overdense_spline_SFR;
float *Nion_spline;
#ifdef INHOMO_FEEDBACK
float *second_derivs_Nion[3];
float *Nion_splinem;
float *second_derivs_Nionm[3];
#else
float *second_derivs_Nion;
#ifdef MINI_HALO
float *Nion_splinem;
float *second_derivs_Nionm;
#endif
#endif

float Nion_ConditionallnM_GL(float M, struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con);
float GaussLegendreQuad_Nion(int n, float z, float M2, float delta1, float delta2, float MassTurnover, float Alpha_star, float Alpha_esc, float Fstar10, float Fesc10, float Mlim_Fstar, float Mlim_Fesc);
double dNion_ConditionallnM(double lnM, void *params);
double Nion_ConditionalM(double z, double M1, double M2, double delta1, double delta2, double MassTurnover, double Alpha_star, double Alpha_esc, double Fstar10, double Fesc10, double Mlim_Fstar, double Mlim_Fesc);

double dNion_ST(double lnM, void *params);
double Nion_ST(double z, double M_Min, double MassTurnover, double Alpha_star, double Alpha_esc, double Fstar10, double Fesc10, double Mlim_Fstar, double Mlim_Fesc);
/* New in v2 - part 3 of 4: end */

double sigma_norm, theta_cmb, omhh, z_equality, y_d, sound_horizon, alpha_nu, f_nu, f_baryon, beta_c, R_CUTOFF;


/*****     FUNCTION PROTOTYPES     *****/
double init_ps(); /* initialize global variables, MUST CALL THIS FIRST!!! returns R_CUTOFF */
void free_ps(); /* deallocates the gsl structures from init_ps */
float mean_SFRD(double z); // returns the mean star formation rate density at z in M_sun yr^-1 Mpc^-3
double splined_erfc(double); /* returns erfc for x>=0, using cubic spline in logy-x space */
double deltolindel(float del, float z); /* converts a non-linear overdensity, del, at z to a linear overdensity at z=0 */
double lindeltodel(float lindel, float z); /* converts a linear overdensity, del, at z=0 to a non-linear overdensity at redshift z */
double power_in_k(double k); /* Returns the value of the linear power spectrum density (i.e. <|delta_k|^2>/V) at a given k mode at z=0 */
double RtoM(double); /* R in Mpc, M in Msun */
double MtoR(double); /* R in Mpc, M in Msun */
double M_J_WDM(); /* returns the "effective Jeans mass" corresponding to the gas analog of WDM ; eq. 10 in BHO 2001 */
double sheth_delc(double del, double sig);
double dNdM_st(double z, double M);
double dNdM(double z, double M);
double dnbiasdM(double M, float z, double M_o, float del_o); /* dnbiasdM */
double FgtrM(double z, double M);  //calculates the fraction of mass contained in haloes with mass > M at redshift z
double FgtrM_st(double z, double M);  //calculates the fraction of mass contained in haloes with mass > M at redshift z, with Sheth-Tormen correction
double FgtrM_bias(double z, double M, double del_bias, double sig_bias);  //calculates the fraction of mass contained in haloes with mass > M at redshift z, in regions with a linear overdensity of del_bias, and standard deviation sig_bias
double sigmaparam_FgtrM_bias(float z, float sigsmallR, float del_bias, float sig_bias);/* Uses sigma parameters instead of Mass for scale */

double FgtrM_bias_BL08(double z, double M, double del_bias, double sig_bias); // as above, but this version uses the hybrid perscription of Barkana & Loeb 2004 (specifically the separate integral version of eq. 2 in Barkana & Loeb 2008)
double dicke(double z); //calculates the dicke growth function at redshift z
double ddickedz(double z); /* Redshift derivative of the growth function at z */
double ddickedt(double z); /* Time derivative of the growth function at z */
double sigma_z0(double M); //calculates sigma at z=0 (no dicke)
double dsigmasqdm_z0(double M); //calculates d(sigma^2)/dm at z=0 (i.e. does not include dicke growth)
double TFmdm(double k); //Eisenstien & Hu power spectrum transfer function
void TFset_parameters();
float get_R_c();  // returns R_CUTOFF
double get_M_min_ion(float z);
/***************************************/



/* Returns the minimum source mass for ionizing sources, according to user specifications */
double get_M_min_ion(float z){
  double MMIN;

  MMIN = M_TURNOVER;

  // check for WDM
#ifdef P_CUTOFF
  if ( MMIN < M_J_WDM())
    MMIN = M_J_WDM();
#endif //P_CUTOFF
  return MMIN;
}

/* Returns the minimum source mass for x-ray sources, according to user specifications */
double get_M_min_xray(float z){
  return get_M_min_ion(z);
}


/* returns the "effective Jeans mass" in Msun
   corresponding to the gas analog of WDM ; eq. 10 in Barkana+ 2001 */
double M_J_WDM(){
  double z_eq, fudge=60;
#ifndef P_CUTOFF 
    return 0;
#else //P_CUTOFF
  z_eq = 3600*(OMm-OMb)*hlittle*hlittle/0.15;
  return fudge*3.06e8 * (1.5/g_x) * sqrt((OMm-OMb)*hlittle*hlittle/0.15) * pow(M_WDM, -4) * pow(z_eq/3000.0, 1.5);
#endif //P_CUTOFF
}


/* converts a non-linear overdensity, del, at z to a linear overdensity at z=0 */
double deltolindel(float del, float z){
  float onepdel = 1.0+del;
  return ( 1.68647 - 1.35*pow(onepdel,-2/3.0) + 0.78785*pow(onepdel,-0.58661) - 1.12431*pow(onepdel,-0.5) )/dicke(z);
}

/* converts a linear overdensity, del, at z=0 to a non-linear overdensity at redshift z */
double lindeltodel(float lindel, float z){
  float prev_lindelguess, delcrit, delguess;
  float lindelguess, delmin, delmax, epsilon = 1.0e-7;

  // set the critical density corresponding to virialization
  // this will be maximum allowed del
  delcrit = Deltac_nonlinear(z)*rho_critz(z)/(OMm*RHOcrit*pow(1+z, 3)) - 1;

  delmin = -1;
  delmax = 500;
  prev_lindelguess = -1e10;
  while (1){
    delguess = 0.5*(delmax+delmin);
    lindelguess = deltolindel(delguess, z);
    //fprintf(stderr, "%e\t%e\n", delmin, delmax);
          //    fprintf(stderr, "%e\t%e\t%e\n\n", delguess, lindelguess, lindel);
    if ((fabs((lindelguess-lindel)/lindel) < epsilon ) ||
    (fabs(lindelguess-lindel) < epsilon ) || 
    (fabs(prev_lindelguess - lindelguess) < TINY ))// close enough, or resolution loop
      return delguess;

    if (lindelguess > lindel)
      delmax = delguess;
    else
      delmin = delguess;

    // check if we are above delcrit (see above)
    if (delmin > delcrit){
      //      printf("exced max at lindel=%e\n", lindel);
      return delcrit;
    }

    prev_lindelguess = lindelguess;
  }
}


/* R in Mpc, M in Msun */
double RtoM(double R){
  // set M according to M<->R conversion defined by the filter type in ../Parameter_files/COSMOLOGY.H
  if (FILTER == 0) //top hat M = (4/3) PI <rho> R^3
    return (4.0/3.0)*PI*pow(R,3)*(OMm*RHOcrit);
  else if (FILTER == 1) //gaussian: M = (2PI)^1.5 <rho> R^3
    return pow(2*PI, 1.5) * OMm*RHOcrit * pow(R, 3);
  else // filter not defined
    fprintf(stderr, "No such filter = %i.\nResults are bogus.\n", FILTER);
  return -1;
}

/* R in Mpc, M in Msun */
double MtoR(double M){
  // set R according to M<->R conversion defined by the filter type in ../Parameter_files/COSMOLOGY.H
  if (FILTER == 0) //top hat M = (4/3) PI <rho> R^3
    return pow(3*M/(4*PI*OMm*RHOcrit), 1.0/3.0);
  else if (FILTER == 1) //gaussian: M = (2PI)^1.5 <rho> R^3
    return pow( M/(pow(2*PI, 1.5) * OMm * RHOcrit), 1.0/3.0 );
  else // filter not defined
    fprintf(stderr, "No such filter = %i.\nResults are bogus.\n", FILTER);
  return -1;
}


/* equation (5) from jenkis et al. (2001) */
double f_jenkins(float del, double sigsq){
  if (del < 0){  fprintf(stderr, "ERROR:  In function f_jenkins del_o must be less than del_1 = del_crit/dicke(z)!\nAborting...\n"); return 0; }

  //  fprintf(stderr, "%f\t%f\n", del, sqrt(sigsq));
  return sqrt(2/PI) * del/sqrt(sigsq) * pow(E, -0.5*del*del/sigsq);
}


float get_R_c(){
  return R_CUTOFF;
}

/* sheth correction to delta crit */
double sheth_delc(double del, double sig){
  return sqrt(SHETH_a)*del*(1 + SHETH_b*pow(sig*sig/(SHETH_a*del*del), SHETH_c));
}


/* dnbiasdM */
double dnbiasdM(double M, float z, double M_o, float del_o){
  double sigsq, del, sig_one, sig_o;

  if ((M_o-M) < TINY){
    fprintf(stderr, "WARNING:  In function dnbiasdM: M must be less than M_o!\nAborting...\n");
    return -1;
  }
  del = Deltac/dicke(z) - del_o;
  if (del < 0){  fprintf(stderr, "ERROR:  In function dnbiasdM: del_o must be less than del_1 = del_crit/dicke(z)!\nAborting...\n"); return 0; }
  sig_o = sigma_z0(M_o);
  sig_one = sigma_z0(M);
  sigsq = sig_one*sig_one - sig_o*sig_o;
  return -(RHOcrit*OMm)/M /sqrt(2*PI) *del*pow(sigsq,-1.5)*pow(E, -0.5*del*del/sigsq)*dsigmasqdm_z0(M);
}


/*
  FUNCTION dNdM(z, M)
  Computes the Press_schechter mass function with Sheth-Torman correction for ellipsoidal collapse at 
  redshift z, and dark matter halo mass M (in solar masses).

  The return value is the number density per unit mass of halos in the mass range M to M+dM in units of:
  comoving Mpc^-3 Msun^-1

  Reference: Sheth, Mo, Torman 2001
*/
double dNdM_st(double z, double M){
  double sigma, dsigmadm, nuhat, dicke_growth;

  dicke_growth = dicke(z);
  sigma = sigma_z0(M) * dicke_growth;
  dsigmadm = dsigmasqdm_z0(M) * dicke_growth*dicke_growth/(2.0*sigma);
  nuhat = sqrt(SHETH_a) * Deltac / sigma;
  
  return (-OMm*RHOcrit/M) * (dsigmadm/sigma) * sqrt(2/PI)*SHETH_A * (1+ pow(nuhat, -2*SHETH_p)) * nuhat * pow(E, -nuhat*nuhat/2.0);
}



/*
  FUNCTION dNdM(z, M)
  Computes the Press_schechter mass function at 
  redshift z, and dark matter halo mass M (in solar masses).

  The return value is the number density per unit mass of halos in the mass range M to M+dM in units of:
  comoving Mpc^-3 Msun^-1

  Reference: Padmanabhan, pg. 214
*/
double dNdM(double z, double M){
  double sigma, dsigmadm, dicke_growth;

  dicke_growth = dicke(z);
  sigma = sigma_z0(M) * dicke_growth;
  dsigmadm = dsigmasqdm_z0(M) * (dicke_growth*dicke_growth/(2*sigma));

  return (-OMm*RHOcrit/M) * sqrt(2/PI) * (Deltac/(sigma*sigma)) * dsigmadm * pow(E, -(Deltac*Deltac)/(2*sigma*sigma));
}


/*
  FUNCTION FgtrM_st(z, M)
  Computes the fraction of mass contained in haloes with mass > M at redshift z
  Uses Sheth-Torman correction
*/
double dFdlnM_st (double lnM, void *params){
  double z = *(double *)params;
  double M = exp(lnM);
  return dNdM_st(z, M) * M * M;
}
double FgtrM_st(double z, double M){
  double result, error, lower_limit, upper_limit;
  gsl_function F;
  double rel_tol  = 0.001; //<- relative tolerance
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1000);

  F.function = &dFdlnM_st;
  F.params = &z;
  lower_limit = log(M);
  upper_limit = log(FMAX(1e16, M*100));
               
  gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
               1000, GSL_INTEG_GAUSS61, w, &result, &error); 
  gsl_integration_workspace_free (w);

  return result / (OMm*RHOcrit);
}


/*
  FUNCTION FgtrM(z, M)
  Computes the fraction of mass contained in haloes with mass > M at redshift z
*/
double FgtrM(double z, double M){
  double del, sig;

  del = Deltac/dicke(z); //regular spherical collapse delta
  sig = sigma_z0(M);

  return splined_erfc(del / (sqrt(2)*sig));
}


/*
  calculates the fraction of mass contained in haloes with mass > M at redshift z, in regions with a linear overdensity of del_bias, and standard deviation sig_bias
 */
double FgtrM_bias(double z, double M, double del_bias, double sig_bias){
  double del, sig, sigsmallR;

  sigsmallR = sigma_z0(M);

  if (!(sig_bias < sigsmallR)){ // biased region is smaller that halo!
    fprintf(stderr, "FgtrM_bias: Biased region is smaller than halo!\nResult is bogus.\n");
    return 0;
  }

  del = Deltac/dicke(z) - del_bias;
  sig = sqrt(sigsmallR*sigsmallR - sig_bias*sig_bias);

   return splined_erfc(del / (sqrt(2)*sig));
}


/* Uses sigma parameters instead of Mass for scale */
double sigmaparam_FgtrM_bias(float z, float sigsmallR, float del_bias, float sig_bias){
  double del, sig;

  if (!(sig_bias < sigsmallR)){ // biased region is smaller that halo!
    fprintf(stderr, "local_FgtrM_bias: Biased region is smaller than halo!\nResult is bogus.\n");
    fprintf(stderr, "In function: zpp = %.4f, sigma_Tmin = %.4e, sigma_atR = %.4e, delta = %.5f\n",z,sigsmallR,sig_bias,del_bias);
    return 0;
  }

  del = Deltac/dicke(z) - del_bias;
  sig = sqrt(sigsmallR*sigsmallR - sig_bias*sig_bias);

  return splined_erfc(del / (sqrt(2)*sig));
}


/*
  Calculates the fraction of mass contained in haloes with mass > M at redshift z, in regions with a linear overdensity of del_bias, and standard deviation sig_bias.
  This version uses the hybrid perscription of Barkana & Loeb 2004 (specifically the separate
  integral version of eq. 2 in Barkana & Loeb 2008)
 */
double FgtrM_bias_BL08(double z, double M, double del_bias, double sig_bias){
  return FgtrM_st(z, M) / FgtrM(z, M) * FgtrM_bias(z, M, del_bias, sig_bias);
}


/*
  FUNCTION dicke(z)
  Computes the dicke growth function at redshift z, i.e. the z dependance part of sigma

  References: Peebles, "Large-Scale...", pg.53 (eq. 11.16). Includes omega<=1
  Nonzero Lambda case from Liddle et al, astro-ph/9512102, eqs. 6-8.
  and quintessence case from Wang et al, astro-ph/9804015

  Normalized to dicke(z=0)=1
*/
double dicke(double z){
  double omegaM_z, dick_z, dick_0, x, x_0;
  double tiny = 1e-4;

  if (fabs(OMm-1.0) < tiny){ //OMm = 1 (Einstein de-Sitter)
    return 1.0/(1.0+z);
  }
  else if ( (OMl > (-tiny)) && (fabs(OMl+OMm+OMr-1.0) < 0.01) && (fabs(wl+1.0) < tiny) ){
    //this is a flat, cosmological CONSTANT universe, with only lambda, matter and radiation
    //it is taken from liddle et al.
    omegaM_z = OMm*pow(1+z,3) / ( OMl + OMm*pow(1+z,3) + OMr*pow(1+z,4) );
    dick_z = 2.5*omegaM_z / ( 1.0/70.0 + omegaM_z*(209-omegaM_z)/140.0 + pow(omegaM_z, 4.0/7.0) );
    dick_0 = 2.5*OMm / ( 1.0/70.0 + OMm*(209-OMm)/140.0 + pow(OMm, 4.0/7.0) );
    return dick_z / (dick_0 * (1.0+z));
  }
  else if ( (OMtot < (1+tiny)) && (fabs(OMl) < tiny) ){ //open, zero lambda case (peebles, pg. 53)
    x_0 = 1.0/(OMm+0.0) - 1.0;
    dick_0 = 1 + 3.0/x_0 + 3*log(sqrt(1+x_0)-sqrt(x_0))*sqrt(1+x_0)/pow(x_0,1.5);
    x = fabs(1.0/(OMm+0.0) - 1.0) / (1+z);
    dick_z = 1 + 3.0/x + 3*log(sqrt(1+x)-sqrt(x))*sqrt(1+x)/pow(x,1.5);
    return dick_z/dick_0;
  }
  else if ( (OMl > (-tiny)) && (fabs(OMtot-1.0) < tiny) && (fabs(wl+1) > tiny) ){
    fprintf(stderr, "IN WANG\n");
    return -1;
  }

  fprintf(stderr, "No growth function!!! Output will be fucked up.");
  return -1;
}


/* redshift derivative of the growth function at z */
double ddicke_dz(double z){
  float dz = 1e-10;
  double omegaM_z, ddickdz, dick_0, x, x_0, domegaMdz;

   return (dicke(z+dz)-dicke(z))/dz;
}

/* Time derivative of the growth function at z */
double ddickedt(double z){
  float dz = 1e-10;
  double omegaM_z, ddickdz, dick_0, x, x_0, domegaMdz;
  double tiny = 1e-4;

   return (dicke(z+dz)-dicke(z))/dz/dtdz(z); // lazy non-analytic form getting

  if (fabs(OMm-1.0) < tiny){ //OMm = 1 (Einstein de-Sitter)
    return -pow(1+z,-2)/dtdz(z);
  }
  else if ( (OMl > (-tiny)) && (fabs(OMl+OMm+OMr-1.0) < 0.01) && (fabs(wl+1.0) < tiny) ){
    //this is a flat, cosmological CONSTANT universe, with only lambda, matter and radiation
    //it is taken from liddle et al.
    omegaM_z = OMm*pow(1+z,3) / ( OMl + OMm*pow(1+z,3) + OMr*pow(1+z,4) );
    domegaMdz = omegaM_z*3/(1+z) - OMm*pow(1+z,3)*pow(OMl + OMm*pow(1+z,3) + OMr*pow(1+z,4), -2) * (3*OMm*(1+z)*(1+z) + 4*OMr*pow(1+z,3));
    dick_0 = OMm / ( 1.0/70.0 + OMm*(209-OMm)/140.0 + pow(OMm, 4.0/7.0) );

    ddickdz = (domegaMdz/(1+z)) * (1.0/70.0*pow(omegaM_z,-2) + 1.0/140.0 + 3.0/7.0*pow(omegaM_z, -10.0/3.0)) * pow(1.0/70.0/omegaM_z + (209.0-omegaM_z)/140.0 + pow(omegaM_z, -3.0/7.0) , -2);
    ddickdz -= pow(1+z,-2)/(1.0/70.0/omegaM_z + (209.0-omegaM_z)/140.0 + pow(omegaM_z, -3.0/7.0));

    return ddickdz / dick_0 / dtdz(z);
  }

  fprintf(stderr, "No growth function!!! Output will be fucked up.");
  return -1;
}

/*
  FUNCTION sigma_z0(M)
  Returns the standard deviation of the normalized, density excess (delta(x)) field,
  smoothed on the comoving scale of M (see filter definitions for M<->R conversion).
  The sigma is evaluated at z=0, with the time evolution contained in the dicke(z) factor,
  i.e. sigma(M,z) = sigma_z0(m) * dicke(z)

  normalized so that sigma_z0(M->8/h Mpc) = SIGMA8 in ../Parameter_files/COSMOLOGY.H

  NOTE: volume is normalized to = 1, so this is equvalent to the mass standard deviation

  M is in solar masses

  References: Padmanabhan, pg. 210, eq. 5.107
*/
double dsigma_dk(double k, void *params){
  double p, w, T, gamma, q, aa, bb, cc, kR;
  double R = *(double *)params;

  // get the power spectrum.. choice of 5:
  if (POWER_SPECTRUM == 0){ // Eisenstein & Hu
    T = TFmdm(k);
    // check if we should cuttoff power spectrum according to Bode et al. 2000 transfer function
#ifdef P_CUTOFF 
    T *= pow(1. + pow(BODE_e*k*R_CUTOFF, 2*BODE_v), -BODE_n/BODE_v);
#endif //P_CUTOFF
    p = pow(k, POWER_INDEX) * T * T;
  }
  else if (POWER_SPECTRUM == 1){ // BBKS
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    q = k / (hlittle*gamma);
    T = (log(1.0+2.34*q)/(2.34*q)) * 
      pow( 1.0+3.89*q + pow(16.1*q, 2) + pow( 5.46*q, 3) + pow(6.71*q, 4), -0.25);
     p = pow(k, POWER_INDEX) * T * T;
  }
  else if (POWER_SPECTRUM == 2){ // Efstathiou,G., Bond,J.R., and White,S.D.M., MNRAS,258,1P (1992)
    gamma = 0.25;
    aa = 6.4/(hlittle*gamma);
    bb = 3.0/(hlittle*gamma);
    cc = 1.7/(hlittle*gamma);
    p = pow(k, POWER_INDEX) / pow( 1+pow( aa*k + pow(bb*k, 1.5) + pow(cc*k,2), 1.13), 2.0/1.13 );
  }
  else if (POWER_SPECTRUM == 3){ // Peebles, pg. 626
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    aa = 8.0 / (hlittle*gamma);
    bb = 4.7 / pow(hlittle*gamma, 2);
    p = pow(k, POWER_INDEX) / pow(1 + aa*k + bb*k*k, 2);
  }
  else if (POWER_SPECTRUM == 4){ // White, SDM and Frenk, CS, 1991, 379, 52
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    aa = 1.7/(hlittle*gamma);
    bb = 9.0/pow(hlittle*gamma, 1.5);
    cc = 1.0/pow(hlittle*gamma, 2);
    p = pow(k, POWER_INDEX) * 19400.0 / pow(1 + aa*k + bb*pow(k, 1.5) + cc*k*k, 2);
  }
  else{
    fprintf(stderr, "No such power spectrum defined: %i\nOutput is bogus.\n", POWER_SPECTRUM);
    p = 0;
  }


  // now get the value of the window function
  // NOTE: only use top hat for SIGMA8 normalization
  kR = k*R;
  if ( (FILTER == 0) || (sigma_norm < 0) ){ // top hat
    if ( (kR) < 1.0e-4 ){ w = 1.0;} // w converges to 1 as (kR) -> 0
    else { w = 3.0 * (sin(kR)/pow(kR, 3) - cos(kR)/pow(kR, 2));}
  }
  else if (FILTER == 1){ // gaussian of width 1/R
    w = pow(E, -kR*kR/2.0);
  }
  else {
    fprintf(stderr, "No such filter: %i\nOutput is bogus.\n", FILTER);
    w=0;
  }

  return k*k*p*w*w;
}

double sigma_z0(double M){
  double result, error, lower_limit, upper_limit;
  gsl_function F;
  double rel_tol  = FRACT_FLOAT_ERR*10; //<- relative tolerance
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1000);
  double kstart, kend;
  double R = MtoR(M);

  // now lets do the integral for sigma and scale it with sigma_norm
  kstart = 1.0e-99/R;
  kend = 350.0/R;
  lower_limit = kstart;//log(kstart);
  upper_limit = kend;//log(kend);

  F.function = &dsigma_dk;
  F.params   = &R;
  gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
               1000, GSL_INTEG_GAUSS61, w, &result, &error); 
  gsl_integration_workspace_free (w);
  
  return sigma_norm * sqrt(result);
}


/*
  Returns the value of the linear power spectrum DENSITY (i.e. <|delta_k|^2>/V)
  at a given k mode linearly extrapolated to z=0
*/
double power_in_k(double k){
  double p, T, gamma, q, aa, bb, cc;

  // get the power spectrum.. choice of 5:
  if (POWER_SPECTRUM == 0){ // Eisenstein & Hu
    T = TFmdm(k);
    // check if we should cuttoff power spectrum according to Bode et al. 2000 transfer function
#ifdef P_CUTOFF 
    T *= pow(1 + pow(BODE_e*k*R_CUTOFF, 2*BODE_v), -BODE_n/BODE_v);
#endif //P_CUTOFF
    p = pow(k, POWER_INDEX) * T * T;
    //p = pow(k, POWER_INDEX - 0.05*log(k/0.05)) * T * T; //running, alpha=0.05
  }
  else if (POWER_SPECTRUM == 1){ // BBKS
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    q = k / (hlittle*gamma);
    T = (log(1.0+2.34*q)/(2.34*q)) * 
      pow( 1.0+3.89*q + pow(16.1*q, 2) + pow( 5.46*q, 3) + pow(6.71*q, 4), -0.25);
     p = pow(k, POWER_INDEX) * T * T;
  }
  else if (POWER_SPECTRUM == 2){ // Efstathiou,G., Bond,J.R., and White,S.D.M., MNRAS,258,1P (1992)
    gamma = 0.25;
    aa = 6.4/(hlittle*gamma);
    bb = 3.0/(hlittle*gamma);
    cc = 1.7/(hlittle*gamma);
    p = pow(k, POWER_INDEX) / pow( 1+pow( aa*k + pow(bb*k, 1.5) + pow(cc*k,2), 1.13), 2.0/1.13 );
  }
  else if (POWER_SPECTRUM == 3){ // Peebles, pg. 626
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    aa = 8.0 / (hlittle*gamma);
    bb = 4.7 / pow(hlittle*gamma, 2);
    p = pow(k, POWER_INDEX) / pow(1 + aa*k + bb*k*k, 2);
  }
  else if (POWER_SPECTRUM == 4){ // White, SDM and Frenk, CS, 1991, 379, 52
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    aa = 1.7/(hlittle*gamma);
    bb = 9.0/pow(hlittle*gamma, 1.5);
    cc = 1.0/pow(hlittle*gamma, 2);
    p = pow(k, POWER_INDEX) * 19400.0 / pow(1 + aa*k + bb*pow(k, 1.5) + cc*k*k, 2);
  }
  else{
    fprintf(stderr, "No such power spectrum defined: %i\nOutput is bogus.\n", POWER_SPECTRUM);
    p = 0;
  }


  return p*TWOPI*PI*sigma_norm*sigma_norm;
}


/*
  FUNCTION dsigmasqdm_z0(M)
  returns  d/dm (sigma^2) (see function sigma), in units of Msun^-1
*/
double dsigmasq_dm(double k, void *params){
  double p, w, T, gamma, q, aa, bb, cc, dwdr, drdm, kR;
  double R = *(double *)params;

  // get the power spectrum.. choice of 5:
  if (POWER_SPECTRUM == 0){ // Eisenstein & Hu ApJ, 1999, 511, 5
    T = TFmdm(k);
    // check if we should cuttoff power spectrum according to Bode et al. 2000 transfer function
#ifdef P_CUTOFF
    T *= pow(1 + pow(BODE_e*k*R_CUTOFF, 2*BODE_v), -BODE_n/BODE_v);
#endif //P_CUTOFF
    p = pow(k, POWER_INDEX) * T * T;
    //p = pow(k, POWER_INDEX - 0.05*log(k/0.05)) * T * T; //running, alpha=0.05
  }
  else if (POWER_SPECTRUM == 1){ // BBKS
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    q = k / (hlittle*gamma);
    T = (log(1.0+2.34*q)/(2.34*q)) * 
      pow( 1.0+3.89*q + pow(16.1*q, 2) + pow( 5.46*q, 3) + pow(6.71*q, 4), -0.25);
    p = pow(k, POWER_INDEX) * T * T;
  }
  else if (POWER_SPECTRUM == 2){ // Efstathiou,G., Bond,J.R., and White,S.D.M., MNRAS,258,1P (1992)
    gamma = 0.25;
    aa = 6.4/(hlittle*gamma);
    bb = 3.0/(hlittle*gamma);
    cc = 1.7/(hlittle*gamma);
    p = pow(k, POWER_INDEX) / pow( 1+pow( aa*k + pow(bb*k, 1.5) + pow(cc*k,2), 1.13), 2.0/1.13 );
  }
  else if (POWER_SPECTRUM == 3){ // Peebles, pg. 626
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    aa = 8.0 / (hlittle*gamma);
    bb = 4.7 / (hlittle*gamma);
    p = pow(k, POWER_INDEX) / pow(1 + aa*k + bb*k*k, 2);
  }
  else if (POWER_SPECTRUM == 4){ // White, SDM and Frenk, CS, 1991, 379, 52
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    aa = 1.7/(hlittle*gamma);
    bb = 9.0/pow(hlittle*gamma, 1.5);
    cc = 1.0/pow(hlittle*gamma, 2);
    p = pow(k, POWER_INDEX) * 19400.0 / pow(1 + aa*k + pow(bb*k, 1.5) + cc*k*k, 2);
  }
  else{
    fprintf(stderr, "No such power spectrum defined: %i\nOutput is bogus.\n", POWER_SPECTRUM);
    p = 0;
  }


  // now get the value of the window function
  kR = k * R;
  if (FILTER == 0){ // top hat
    if ( (kR) < 1.0e-4 ){ w = 1.0; }// w converges to 1 as (kR) -> 0
    else { w = 3.0 * (sin(kR)/pow(kR, 3) - cos(kR)/pow(kR, 2));}

    // now do d(w^2)/dm = 2 w dw/dr dr/dm
    if ( (kR) < 1.0e-10 ){  dwdr = 0;}
    else{ dwdr = 9*cos(kR)*k/pow(kR,3) + 3*sin(kR)*(1 - 3/(kR*kR))/(kR*R);}
        //3*k*( 3*cos(kR)/pow(kR,3) + sin(kR)*(-3*pow(kR, -4) + 1/(kR*kR)) );}
      //     dwdr = -1e8 * k / (R*1e3);
    drdm = 1.0 / (4.0*PI * OMm*RHOcrit * R*R);
  }
  else if (FILTER == 1){ // gaussian of width 1/R
    w = pow(E, -kR*kR/2.0);
    dwdr = - k*kR * w;
    drdm = 1.0 / (pow(2*PI, 1.5) * OMm*RHOcrit * 3*R*R);
  }
  else {
    fprintf(stderr, "No such filter: %i\nOutput is bogus.\n", FILTER);
    w=0;
  }

  //  printf("%e\t%e\t%e\t%e\t%e\t%e\t%e\n", k, R, p, w, dwdr, drdm, dsigmadk[1]);
  return k*k*p*2*w*dwdr*drdm;
}

double dsigmasqdm_z0(double M){
  double result, error, lower_limit, upper_limit;
  gsl_function F;
  double rel_tol  = FRACT_FLOAT_ERR*10; //<- relative tolerance
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1000);
  double kstart, kend;
  double R = MtoR(M);

  // now lets do the integral for sigma and scale it with sigma_norm
  kstart = 1.0e-99/R;
  kend = 350.0/R;
  lower_limit = kstart;//log(kstart);
  upper_limit = kend;//log(kend);

  F.function = &dsigmasq_dm;
  F.params = &R;
  gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
               1000, GSL_INTEG_GAUSS61, w, &result, &error); 
  gsl_integration_workspace_free (w);
  
  return sigma_norm * sigma_norm * result;
}



/*
  FUNCTION TFmdm is the power spectrum transfer function from Eisenstein & Hu ApJ, 1999, 511, 5
*/
double TFmdm(double k){
  double q, gamma_eff, q_eff, TF_m, q_nu;

  q = k*pow(theta_cmb,2)/omhh;
  gamma_eff=sqrt(alpha_nu) + (1.0-sqrt(alpha_nu))/(1.0+pow(0.43*k*sound_horizon, 4));     
  q_eff = q/gamma_eff;
  TF_m= log(E+1.84*beta_c*sqrt(alpha_nu)*q_eff);
  TF_m /= TF_m + pow(q_eff,2) * (14.4 + 325.0/(1.0+60.5*pow(q_eff,1.11)));
  q_nu = 3.92*q/sqrt(f_nu/N_nu);
  TF_m *= 1.0 + (1.2*pow(f_nu,0.64)*pow(N_nu,0.3+0.6*f_nu)) / 
    (pow(q_nu,-1.6)+pow(q_nu,0.8));

  //   printf("%f  %e  %f  %f  %f  %f\n",omhh,f_nu,f_baryon,N_nu,y_d,alpha_nu);
  // printf("%f  %f  %f  %f\n", beta_c,sound_horizon,theta_cmb,z_equality);
  //printf("%f  %e  %f  %f  %f\n\n",q, k, gamma_eff, q_nu, TF_m);


  return TF_m;
}


void TFset_parameters(){
  double z_drag, R_drag, R_equality, p_c, p_cb, f_c, f_cb, f_nub, k_equality;

  z_equality = 25000*omhh*pow(theta_cmb, -4) - 1.0;
  k_equality = 0.0746*omhh/(theta_cmb*theta_cmb);

  z_drag = 0.313*pow(omhh,-0.419) * (1 + 0.607*pow(omhh, 0.674));
  z_drag = 1 + z_drag*pow(OMb*hlittle*hlittle, 0.238*pow(omhh, 0.223));
  z_drag *= 1291 * pow(omhh, 0.251) / (1 + 0.659*pow(omhh, 0.828));

  y_d = (1 + z_equality) / (1.0 + z_drag);

  R_drag = 31.5 * OMb*hlittle*hlittle * pow(theta_cmb, -4) * 1000 / (1.0 + z_drag);
  R_equality = 31.5 * OMb*hlittle*hlittle * pow(theta_cmb, -4) * 1000 / (1.0 + z_equality);

  sound_horizon = 2.0/3.0/k_equality * sqrt(6.0/R_equality) * 
    log( (sqrt(1+R_drag) + sqrt(R_drag+R_equality)) / (1.0 + sqrt(R_equality)) );

  p_c = -(5 - sqrt(1 + 24*(1 - f_nu-f_baryon)))/4.0;
  p_cb = -(5 - sqrt(1 + 24*(1 - f_nu)))/4.0;
  f_c = 1 - f_nu - f_baryon;
  f_cb = 1 - f_nu;
  f_nub = f_nu+f_baryon;

  alpha_nu = (f_c/f_cb) * (2*(p_c+p_cb)+5)/(4*p_cb+5.0);
  alpha_nu *= 1 - 0.553*f_nub+0.126*pow(f_nub,3);
  alpha_nu /= 1-0.193*sqrt(f_nu)+0.169*f_nu;
  alpha_nu *= pow(1+y_d, p_c-p_cb);
  alpha_nu *= 1+ (p_cb-p_c)/2.0 * (1.0+1.0/(4.0*p_c+3.0)/(4.0*p_cb+7.0))/(1.0+y_d);
  beta_c = 1.0/(1.0-0.949*f_nub);
}




double init_ps(){
  double result, error, lower_limit, upper_limit;
  gsl_function F;
  double rel_tol  = FRACT_FLOAT_ERR*10; //<- relative tolerance
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1000);
  double kstart, kend;
  int i;
  double x;

  // Set cuttoff scale for WDM (eq. 4 in Barkana et al. 2001) in comoving Mpc
  R_CUTOFF = 0.201*pow((OMm-OMb)*hlittle*hlittle/0.15, 0.15)*pow(g_x/1.5, -0.29)*pow(M_WDM, -1.15);

#ifdef P_CUTOFF
  fprintf(stderr, "For M_DM = %.2e keV, R_CUTOFF is: %.2e comoving Mpc\n", M_WDM, R_CUTOFF);
#endif //P_CUTOFF

  omhh = OMm*hlittle*hlittle;
  theta_cmb = T_cmb / 2.7;

  // Translate Parameters into forms GLOBALVARIABLES form
  f_nu = OMn/OMm;
  f_baryon = OMb/OMm;
  if (f_nu < TINY) f_nu = 1e-10;
  if (f_baryon < TINY) f_baryon = 1e-10;


  TFset_parameters();

  sigma_norm = -1;

  double R = 8.0/hlittle;
  kstart = 1.0e-99/R;
  kend = 350.0/R;
  lower_limit = kstart;//log(kstart);
  upper_limit = kend;//log(kend);

  F.function = &dsigma_dk;
  F.params   = &R;

  gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
               1000, GSL_INTEG_GAUSS61, w, &result, &error); 
  gsl_integration_workspace_free (w);

  sigma_norm = SIGMA8/sqrt(result); //takes care of volume factor


  /* initialize the lookup table for erfc */
  /*
  for (i=0; i<=ERFC_NPTS; i++){
    erfc_params[i] = i*ERFC_PARAM_DELTA;
    log_erfc_table[i] = log(erfcc(erfc_params[i]));
  }
  // Set up spline table
  erfc_acc   = gsl_interp_accel_alloc ();
  erfc_spline  = gsl_spline_alloc (gsl_interp_cspline, ERFC_NPTS);
  gsl_spline_init(erfc_spline, erfc_params, log_erfc_table, ERFC_NPTS);
  */

  return R_CUTOFF;
}

void free_ps(){
  /*    gsl_spline_free (erfc_spline);
    gsl_interp_accel_free(erfc_acc);
  */
  return;
}

double splined_erfc(double x){
  if (x < 0){
    //    fprintf(stderr, "WARNING: Negative value %e passed to splined_erfc. Returning 1\n", x);
    return 1;
  }
  return erfcc(x); // the interpolation below doesn't seem to be stable in Ts.c
  if (x > ERFC_PARAM_DELTA*(ERFC_NPTS-1))
    return erfcc(x);
  else
    return exp(gsl_spline_eval(erfc_spline, x, erfc_acc));
}

float FgtrConditionalM_second(float z, float M1, float M2, float MFeedback, float alpha, float delta1, float delta2) {
    
    return exp(M1)*pow(exp(M1)/MFeedback,alpha)*dNdM_conditional_second(z,M1,M2,delta1,delta2)/sqrt(2.*PI);
}

float dNdM_conditional_second(float z, float M1, float M2, float delta1, float delta2){
    
    float sigma1, sigma2, dsigmadm, dicke_growth,dsigma_val;
    
    M1 = exp(M1);
    M2 = exp(M2);
    
    dicke_growth = dicke(z);
    
    splint(Mass_Spline-1,Sigma_Spline-1,second_derivs_sigma-1,(int)NMass,M1,&(sigma1));
    splint(Mass_Spline-1,Sigma_Spline-1,second_derivs_sigma-1,(int)NMass,M2,&(sigma2));
    
    sigma1 = sigma1*sigma1;
    sigma2 = sigma2*sigma2;

    //printf("sigma_min^2 = %.4e sigma^2 = %.4e\n",sigma1,sigma2);
    
    splint(Mass_Spline-1,dSigmadm_Spline-1,second_derivs_dsigma-1,(int)NMass,M1,&(dsigma_val));
    
    dsigmadm = -pow(10.,dsigma_val)/(2.0*sigma1); // This is actually sigma1^{2} as calculated above, however, it should just be sigma1. It cancels with the same factor below. Why I have decided to write it like that I don't know!
    //if(z > 33 && z < 35) 
    //printf("dNdM_conditional: z = %.4f, sigma_min^2 = %.4e sigma^2 = %.4e, delta1 = %.4e, delta2 = %.4e \n",z,sigma1,sigma2,delta1,delta2);
    if((sigma1 > sigma2)) {
        
        return -(( delta1 - delta2 )/dicke_growth)*( 2.*sigma1*dsigmadm )*( exp( - ( delta1 - delta2 )*( delta1 - delta2 )/( 2.*dicke_growth*dicke_growth*( sigma1 - sigma2 ) ) ) )/(pow( sigma1 - sigma2, 1.5));
    }
    else if(sigma1==sigma2) {
        
        return -(( delta1 - delta2 )/dicke_growth)*( 2.*sigma1*dsigmadm )*( exp( - ( delta1 - delta2 )*( delta1 - delta2 )/( 2.*dicke_growth*dicke_growth*( 1.e-6 ) ) ) )/(pow( 1.e-6, 1.5));
        
    }
    else {
        return 0.;
    }
}

void gauleg(float x1, float x2, float x[], float w[], int n)
//Given the lower and upper limits of integration x1 and x2, and given n, this routine returns arrays x[1..n] and w[1..n] of length n,
//containing the abscissas and weights of the Gauss- Legendre n-point quadrature formula.
{
    
    int m,j,i;
    double z1,z,xm,xl,pp,p3,p2,p1;
    
    m=(n+1)/2;
    xm=0.5*(x2+x1);
    xl=0.5*(x2-x1);
    for (i=1;i<=m;i++) {
        //High precision is a good idea for this routine.
        //The roots are symmetric in the interval, so we only have to find half of them.
        //Loop over the desired roots.
        
        z=cos(3.141592654*(i-0.25)/(n+0.5));
        
        //Starting with the above approximation to the ith root, we enter the main loop of refinement by Newton’s method.
        do {
            p1=1.0;
            p2=0.0;
            for (j=1;j<=n;j++) {
                //Loop up the recurrence relation to get the Legendre polynomial evaluated at z.
                p3=p2;
                p2=p1;
                p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
            }
            //p1 is now the desired Legendre polynomial. We next compute pp, its derivative, by a standard relation involving also p2,
            //the polynomial of one lower order.
            pp=n*(z*p1-p2)/(z*z-1.0);
            z1=z;
            z=z1-p1/pp;
        } while (fabs(z-z1) > EPS2);
        x[i]=xm-xl*z;
        x[n+1-i]=xm+xl*z;
        w[i]=2.0*xl/((1.0-z*z)*pp*pp);
        w[n+1-i]=w[i];
    }
}

void nrerror(char error_text[])
{
    fprintf(stderr,"Numerical Recipes run-time error...\n");
    fprintf(stderr,"%s\n",error_text);
    fprintf(stderr,"...now exiting to system...\n");
    exit(1);
}

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
    float *v;
    v = (float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
    if(!v) nrerror("allocation failure in vector()");
    return v - nl + NR_END;
}

void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
    free((FREE_ARG) (v+nl-NR_END));
}

void spline(float x[], float y[], int n, float yp1, float ypn, float y2[])
/*Given arrays x[1..n] and y[1..n] containing a tabulated function, i.e., yi = f(xi), with
 x1 <x2 < :: : < xN, and given values yp1 and ypn for the first derivative of the interpolating
 function at points 1 and n, respectively, this routine returns an array y2[1..n] that contains
 the second derivatives of the interpolating function at the tabulated points xi. If yp1 and/or
 ypn are equal to 1e30 or larger, the routine is signaled to set the corresponding boundary
 condition for a natural spline, with zero second derivative on that boundary.*/
{
    int i,k;
    float p,qn,sig,un,*u;
    int na,nb,check;
    u=vector(1,n-1);
    if (yp1 > 0.99e30)                     // The lower boundary condition is set either to be "natural"
        y2[1]=u[1]=0.0;
    else {                                 // or else to have a specified first derivative.
        y2[1] = -0.5;
        u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
    }
    for (i=2;i<=n-1;i++) {                              //This is the decomposition loop of the tridiagonal algorithm.
        sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);                //y2 and u are used for temporary
        na = 1;
        nb = 1;
        check = 0;
        while(((float)(x[i+na*1]-x[i-nb*1])==(float)0.0)) {
            check = check + 1;
            if(check%2==0) {
                na = na + 1;
            }
            else {
                nb = nb + 1;
            }
            sig=(x[i]-x[i-1])/(x[i+na*1]-x[i-nb*1]);
        }
        p=sig*y2[i-1]+2.0;                                //storage of the decomposed
        y2[i]=(sig-1.0)/p;                                //  factors.
        u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
        u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
        
        if(((float)(x[i+1]-x[i])==(float)0.0) || ((float)(x[i]-x[i-1])==(float)0.0)) {
            na = 0;
            nb = 0;
            check = 0;
            while((float)(x[i+na*1]-x[i-nb])==(float)(0.0) || ((float)(x[i+na]-x[i-nb*1])==(float)0.0)) {
                check = check + 1;
                if(check%2==0) {
                    na = na + 1;
                }
                else {
                    nb = nb + 1;
                }
            }
            u[i]=(y[i+1]-y[i])/(x[i+na*1]-x[i-nb]) - (y[i]-y[i-1])/(x[i+na]-x[i-nb*1]);
            
            u[i]=(6.0*u[i]/(x[i+na*1]-x[i-nb*1])-sig*u[i-1])/p;
            
        }
    }
    if (ypn > 0.99e30)                        //The upper boundary condition is set either to be "natural"
        qn=un=0.0;
    else {                                    //or else to have a specified first derivative.
        qn=0.5;
        un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
    }
    y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
    
    for (k=n-1;k>=1;k--) {                      //This is the backsubstitution loop of the tridiagonal
        y2[k]=y2[k]*y2[k+1]+u[k];               //algorithm.
    }
    free_vector(u,1,n-1);
}


void splint(float xa[], float ya[], float y2a[], int n, float x, float *y)
/*Given the arrays xa[1..n] and ya[1..n], which tabulate a function (with the xai's in order),
 and given the array y2a[1..n], which is the output from spline above, and given a value of
 x, this routine returns a cubic-spline interpolated value y.*/
{
    void nrerror(char error_text[]);
    int klo,khi,k;
    float h,b,a;
    klo=1;                                                  // We will find the right place in the table by means of
    khi=n;                                                  //bisection. This is optimal if sequential calls to this
    while (khi-klo > 1) {                                   //routine are at random values of x. If sequential calls
        k=(khi+klo) >> 1;                                     //are in order, and closely spaced, one would do better
        if (xa[k] > x) khi=k;                                 //to store previous values of klo and khi and test if
        else klo=k;                                           //they remain appropriate on the next call.
    }                                                           // klo and khi now bracket the input value of x.
    h=xa[khi]-xa[klo];
    if (h == 0.0) nrerror("Bad xa input to routine splint");    //The xa's must be distinct.
    a=(xa[khi]-x)/h;
    b=(x-xa[klo])/h;                                            //Cubic spline polynomial is now evaluated.
    *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

#ifdef INHOMO_FEEDBACK
void spline2d(float xi[], float yi[], float zi[], int nx, int ny, float *z2[])
/* Given arrays x[1..n],  y[1..m], and z[11,12,...1n, 21,22,..2n,..,m1,m2,..,mn] containing a tabulated 
 * function, i.e. zij = f(xi, yj), with x1 < x2 <... < xn; y1 < y2 < ... < ym, and assuming the firsti
 * derivative of the interpolating function at the edges z1*, zm*, z*1 and z*n are 0, this routine 
 * returns an array z2[11,12,...1n, 21,22,..2n,..,m1,m2,..,mn] that contains the second derivatives of
 * the interpolating function at the tabulated points xi and yj. I don't do the unit-offset following
 * rwstrickland/multicubic-spline */
{
  void nrerror(char error_text[]);
  int nw, i, j;
  float *w, *s;
  float *x, *y, *z;
  float *sx, *sy, *sxy;
  x = xi;
  y = yi;
  z = zi;
  sx = z2[0];
  sy = z2[1];
  sxy = z2[2];

  nw = (nx > ny) ? nx : ny;
  s = (float *) malloc(nw * sizeof(float));
  w = (float *) malloc(nw * sizeof(float));
  if (s == NULL || w == NULL)
    nrerror("malloc failed");

  for (j = 0; j < ny; j++) {
    for (i = 0; i < nx; i++)
      w[i] = z[i+nx*j];
    spline(x-1, w-1, nx, 0, 0, s-1);
    for (i = 0; i < nx; i++)
      sx[i+nx*j] = s[i];
  }
  
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++)
      w[j] = z[i+nx*j];
    spline(y-1, w-1, ny, 0, 0, s-1);
    for (j = 0; j < ny; j++)
      sy[i+nx*j] = s[j];

    for (j = 0; j < ny; j++)
      w[j] = sx[i+nx*j];
    spline(y-1, w-1, ny, 0, 0, s-1);
    for (j = 0; j < ny; j++)
      sxy[i+nx*j] = s[j];
  }
  free(s);
  free(w);
}

void splint2d(float xi[], float yi[], float zi[], float *z2[], int nx, int ny, float xx, float yy, float *zz)
{
  void nrerror(char error_text[]);
  double a, b, g, c, d;
  int i0, i1, j0, j1, tmp;
  int k00, k01, k10, k11;
  double dx, dy;
  double f0, f1, s0, s1;
  float *x, *y, *z;
  float *sx, *sy, *sxy;

  x = xi;
  y = yi;
  z = zi;
  sx = z2[0];
  sy = z2[1];
  sxy = z2[2];

  /* Linear search for interval containing each independent variable */
  i0=0;
  i1=nx-1;
  while (i1-i0 > 1) {
    tmp=(i1+i0) >> 1;
    if (x[tmp] > xx) i1=tmp;
    else i0=tmp;
  }

  j0=0;
  j1=ny-1;
  while (j1-j0 > 1) {
    tmp=(j1+j0) >> 1;
    if (y[tmp] > yy) j1=tmp;
    else j0=tmp;
  }

  dx = x[i1] - x[i0];
  if (dx <= 0.)
    nrerror("Bad xa input to routine splint");    //The xa's must be distinct.

  dy = y[j1] - y[j0];
  if (dy <= 0.)
    nrerror("Bad xa input to routine splint");    //The xa's must be distinct.

  k00 = i0 + nx*j0;
  k01 = i0 + nx*j1;
  k10 = i1 + nx*j0;
  k11 = i1 + nx*j1;

  /* Interpolation coefficients for y direction */

  b = (yy - y[j0])/dy;
  a = 1.-b;
  g = dy*dy/6.;
  c = g*(a*a*a - a);
  d = g*(b*b*b - b);

  /* Interpolate in y direction for f(x[], y) to get f at points U and V */

  f0 = a*z[k00] + b*z[k01] + c*sy[k00] + d*sy[k01];
  f1 = a*z[k10] + b*z[k11] + c*sy[k10] + d*sy[k11];

  /* Interpolate in y direction for 2nd derivative of f(x[], y) at U and V */

  s0 = a*sx[k00] + b*sx[k01] + c*sxy[k00] + d*sxy[k01];
  s1 = a*sx[k10] + b*sx[k11] + c*sxy[k10] + d*sxy[k11];

  /* Interpolate in x direction for value of function at our point P */

  b = (xx - x[i0])/dx;
  a = 1.-b;
  g = dx*dx/6.;
  c = g*(a*a*a - a);
  d = g*(b*b*b - b);

  *zz = a*f0 + b*f1 + c*s0 + d*s1;

}
#endif //INHOMO_FEEDBACK

unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
    unsigned long *v;
    v = (unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
    if(!v) nrerror("allocation failure in lvector()");
    return v - nl + NR_END;
}

void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
    free((FREE_ARG) (v+nl-NR_END));
}


void initialiseSplinedSigmaM(float M_Min, float M_Max)
{
    int i;
    float Mass;
    
    Mass_Spline = calloc(NMass,sizeof(float));
    Sigma_Spline = calloc(NMass,sizeof(float));
    dSigmadm_Spline = calloc(NMass,sizeof(float));
    second_derivs_sigma = calloc(NMass,sizeof(float));
    second_derivs_dsigma = calloc(NMass,sizeof(float));
    
#pragma omp parallel shared(Mass_Spline, Sigma_Spline, dSigmadm_Spline, M_Min, M_Max) private(i)
{
#pragma omp for
    for(i=0;i<NMass;i++) {
        Mass_Spline[i] = pow(10., log10(M_Min) + (float)i/(NMass-1)*( log10(M_Max) - log10(M_Min) ) );
        Sigma_Spline[i] = sigma_z0(Mass_Spline[i]);
        dSigmadm_Spline[i] = log10(-dsigmasqdm_z0(Mass_Spline[i]));
    }
}
    spline(Mass_Spline-1,Sigma_Spline-1,NMass,0,0,second_derivs_sigma-1);
    spline(Mass_Spline-1,dSigmadm_Spline-1,NMass,0,0,second_derivs_dsigma-1);
}


/* New in v2 - part 4 of 4 start */

/*
 FUNCTION Mass_limit_bisection(M_min, M_max, power-law, normalization)
 Compute a threshold of halo mass above/below which f_{\ast}(M)/f_{esc}(M) exceeds unity/zero.
*/
float Mass_limit (float logM, float PL, float FRAC) {
    return FRAC*pow(pow(10.,logM)/1e10,PL);
}
void bisection(float *x, float xlow, float xup, int *iter){
    *x=(xlow + xup)/2.;
    ++(*iter);
}

float Mass_limit_bisection(float Mmin, float Mmax, float PL, float FRAC){
  int i, iter, max_iter=200;
  float rel_tol=0.001;
  float logMlow, logMupper, x, x1;
  iter = 0;
  logMlow = log10(Mmin);
  logMupper = log10(Mmax);
  
  if (PL < 0.) {
    if (Mass_limit(logMlow,PL,FRAC) < 1.) {
      return Mmin;
    }
  }
  else if (PL > 0.) {
    if (Mass_limit(logMupper,PL,FRAC) < 1.) {
      return Mmax;
    }
  }
  else
    return 0;
  bisection(&x, logMlow, logMupper, &iter);
  do {
    if((Mass_limit(logMlow,PL,FRAC)-1.)*(Mass_limit(x,PL,FRAC)-1.) < 0.) 
      logMupper = x;
    else 
      logMlow = x;
    bisection(&x1, logMlow, logMupper, &iter);
    if(fabs(x1-x) < rel_tol) {
      return pow(10.,x1);        
    }
    x = x1;
  }
  while(iter < max_iter);
  fprintf(stderr, "\n Failed to find a mass limit to regulate stellar fraction/escape fraction is between 0 and 1.\n");
  fprintf(stderr, " The solution does not converge or iterations are not sufficient\n");
  return -1;
}

/*
 FUNCTION Nion_ST(z, Mturn)
 Computes the mean number of IGM ionizing photon per baryon at redshift z.
 Uses Sheth-Torman collapse fraction.
 see eq. (17) in Park et al. 2018
 Note this function DO NOT include ION_EFF_FACTOR, defined as f_{\ast,10}*f_{esc,10}*N_{\gamma/b}. 
 The actual Nion in eq. (17) is the value multiplied by ION_EFF_FACTOR.
 */
double dNion_ST(double lnM, void *params){
    struct parameters_gsl_SFR_int_ vals = *(struct parameters_gsl_SFR_int_ *)params;

    double M = exp(lnM);
    float z = vals.z_obs;
    double MassTurnover = vals.Mdrop;
    double Alpha_star = vals.pl_star;
    double Alpha_esc = vals.pl_esc;
    double Fstar10 = vals.frac_star;
    double Fesc10 = vals.frac_esc;
    double Mlim_Fstar = vals.LimitMass_Fstar;
    double Mlim_Fesc = vals.LimitMass_Fesc;

    double Fstar, Fesc;

    if (Alpha_star > 0. && M > Mlim_Fstar)
        Fstar = 1./Fstar10;
    else if (Alpha_star < 0. && M < Mlim_Fstar)
        Fstar = 1/Fstar10;
    else
        Fstar = pow(M/1e10,Alpha_star);

    if (Alpha_esc > 0. && M > Mlim_Fesc)
        Fesc = 1./Fesc10;
    else if (Alpha_esc < 0. && M < Mlim_Fesc)
        Fesc = 1./Fesc10;
    else
        Fesc = pow(M/1e10,Alpha_esc);
    return dNdM_st(z,M) * M * M * exp(-MassTurnover/M) * Fstar * Fesc;
}

#ifdef MINI_HALO
double dNion_STm(double lnM, void *params){
    struct parameters_gsl_SFR_intm_ vals = *(struct parameters_gsl_SFR_intm_ *)params;

    double M = exp(lnM);
    float z = vals.z_obs;
    double MassTurnover = vals.Mdrop;
    double MassTurnover2 = vals.Mdrop_2;
    double Alpha_star = vals.pl_star;
    double Fstar10 = vals.frac_star;
    double Mlim_Fstar = vals.LimitMass_Fstar;

    double Fstar;

    if (Alpha_star > 0. && M > Mlim_Fstar)
        Fstar = 1./Fstar10;
    else if (Alpha_star < 0. && M < Mlim_Fstar)
        Fstar = 1/Fstar10;
    else
        Fstar = pow(M/1e10,Alpha_star);

    return dNdM_st(z,M) * M * M * exp(-MassTurnover/M) * exp(-M/MassTurnover2) * Fstar;
}
#endif

double Nion_ST(double z, double M_Min, double MassTurnover, double Alpha_star, double Alpha_esc, double Fstar10, double Fesc10, double Mlim_Fstar, double Mlim_Fesc){
    double result, error, lower_limit, upper_limit;
    gsl_function F;
    double rel_tol = 0.001; //<- relative tolerance
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

    struct parameters_gsl_SFR_int_ parameters_gsl_SFR = {
        .z_obs = z,
        .Mdrop = MassTurnover,
        .pl_star = Alpha_star,
        .pl_esc = Alpha_esc,
        .frac_star = Fstar10,
        .frac_esc = Fesc10,
        .LimitMass_Fstar = Mlim_Fstar,
        .LimitMass_Fesc = Mlim_Fesc,
    };

    F.function = &dNion_ST;
    F.params = &parameters_gsl_SFR;
    lower_limit = log(M_Min);
    upper_limit = log(FMAX(1e16, M_Min*100));

    gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                        1000, GSL_INTEG_GAUSS61, w, &result, &error);
    gsl_integration_workspace_free (w);
    return result / (OMm*RHOcrit);

}
#ifdef MINI_HALO
double Nion_STm(double z, double M_Min, double MassTurnover_lower, double MassTurnover_upper, double Alpha_star, double Fstar10, double Mlim_Fstar){
    double result, error, lower_limit, upper_limit;
    gsl_function F;
    double rel_tol = 0.001; //<- relative tolerance
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

    struct parameters_gsl_SFR_intm_ parameters_gsl_SFR = {
        .z_obs = z,
        .Mdrop = MassTurnover_lower,
        .Mdrop_2 = MassTurnover_upper,
        .pl_star = Alpha_star,
        .frac_star = Fstar10,
        .LimitMass_Fstar = Mlim_Fstar,
    };

    F.function = &dNion_STm;
    F.params = &parameters_gsl_SFR;
    lower_limit = log(M_Min);
    upper_limit = log(FMAX(1e16, M_Min*100));

    gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                        1000, GSL_INTEG_GAUSS61, w, &result, &error);
    gsl_integration_workspace_free (w);

    return result / (OMm*RHOcrit);

}
#endif

void initialiseGL_Nion(int n, float M_Min, float M_Max){
    //calculates the weightings and the positions for Gauss-Legendre quadrature.
    gauleg(log(M_Min),log(M_Max),xi_SFR,wi_SFR,n);

}

float Nion_ConditionallnM_GL(float lnM, struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con){
    float M = exp(lnM);
    float z = parameters_gsl_SFR_con.z_obs;
    float M2 = parameters_gsl_SFR_con.Mval;
    float del1 = parameters_gsl_SFR_con.delta1;
    float del2 = parameters_gsl_SFR_con.delta2;
    float MassTurnover = parameters_gsl_SFR_con.Mdrop;
    float Alpha_star = parameters_gsl_SFR_con.pl_star;
    float Alpha_esc = parameters_gsl_SFR_con.pl_esc;
    float Fstar10 = parameters_gsl_SFR_con.frac_star;
    float Fesc10 = parameters_gsl_SFR_con.frac_esc;
    float Mlim_Fstar = parameters_gsl_SFR_con.LimitMass_Fstar;
    float Mlim_Fesc = parameters_gsl_SFR_con.LimitMass_Fesc;

    float Fstar,Fesc;

    if (Alpha_star > 0. && M > Mlim_Fstar)
        Fstar = 1./Fstar10;
    else if (Alpha_star < 0. && M < Mlim_Fstar)
        Fstar = 1./Fstar10;
    else
        Fstar = pow(M/1e10,Alpha_star);

    if (Alpha_esc > 0. && M > Mlim_Fesc)
        Fesc = 1./Fesc10;
    else if (Alpha_esc < 0. && M < Mlim_Fesc)
        Fesc = 1./Fesc10;
    else
        Fesc = pow(M/1e10,Alpha_esc);

    return M*exp(-MassTurnover/M)*Fstar*Fesc*dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI);

}
#ifdef MINI_HALO
float Nion_ConditionallnM_GLm(float lnM, struct parameters_gsl_SFR_con_intm_ parameters_gsl_SFR_con){
    float M = exp(lnM);
    float z = parameters_gsl_SFR_con.z_obs;
    float M2 = parameters_gsl_SFR_con.Mval;
    float del1 = parameters_gsl_SFR_con.delta1;
    float del2 = parameters_gsl_SFR_con.delta2;
    float Alpha_star = parameters_gsl_SFR_con.pl_star;
    float MassTurnoverm = parameters_gsl_SFR_con.Mdropm1;
    float Mcrit_atom = parameters_gsl_SFR_con.Mdropm2;
    float Fstar10m = parameters_gsl_SFR_con.frac_starm;
    float Mlim_Fstarm = parameters_gsl_SFR_con.LimitMass_Fstarm;

    float Fstarm;

    if (Alpha_star > 0. && M > Mlim_Fstarm)
        Fstarm = 1./Fstar10m;
    else if (Alpha_star < 0. && M < Mlim_Fstarm)
        Fstarm = 1./Fstar10m;
    else
        Fstarm = pow(M/1e10,Alpha_star);

    return M* exp(-MassTurnoverm / M) * exp(-M / Mcrit_atom) * Fstarm*dNdM_conditional_second(z, log(M), M2, del1, del2) / sqrt(2. * PI);
}
#endif

float GaussLegendreQuad_Nion(int n, float z, float M2, float delta1, float delta2, float MassTurnover, float Alpha_star, float Alpha_esc, float Fstar10, float Fesc10, float Mlim_Fstar, float Mlim_Fesc) {
    //Performs the Gauss-Legendre quadrature.
    int i;

    float integrand, x;
    integrand = 0.;

    struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
        .z_obs = z,
        .Mval = M2,
        .delta1 = delta1,
        .delta2 = delta2,
        .Mdrop = MassTurnover,
        .pl_star = Alpha_star,
        .pl_esc = Alpha_esc,
        .frac_star = Fstar10,
        .frac_esc = Fesc10,
        .LimitMass_Fstar = Mlim_Fstar,
        .LimitMass_Fesc = Mlim_Fesc
    };

    if(delta2 > delta1){
        return 1.;
    }
    else{
        for(i=1; i<(n+1); i++){
            x = xi_SFR[i];
            integrand += wi_SFR[i]*Nion_ConditionallnM_GL(x,parameters_gsl_SFR_con);
        }
        return integrand;
    }

}
#ifdef MINI_HALO
float GaussLegendreQuad_Nionm(int n, float z, float M2, float delta1, float delta2, float Alpha_star, float MassTurnoverm, float Mcrit_atom, float Fstar10m, float Mlim_Fstarm) {
    //Performs the Gauss-Legendre quadrature.
    int i;

    float integrand, x;
    integrand = 0.;

    struct parameters_gsl_SFR_con_intm_ parameters_gsl_SFR_con = {
        .z_obs = z,
        .Mval = M2,
        .delta1 = delta1,
        .delta2 = delta2,
        .pl_star = Alpha_star,
        .Mdropm1 = MassTurnoverm,
        .Mdropm2 = Mcrit_atom,
        .frac_starm = Fstar10m,
        .LimitMass_Fstarm = Mlim_Fstarm
    };

    if(delta2 > delta1){
        return 1.;
    }
    else{
        for(i=1; i<(n+1); i++){
            x = xi_SFR[i];
            integrand += wi_SFR[i]*Nion_ConditionallnM_GLm(x,parameters_gsl_SFR_con);
        }
        return integrand;
    }

}
#endif


double dNion_ConditionallnM(double lnM, void *params) {
    struct parameters_gsl_SFR_con_int_ vals = *(struct parameters_gsl_SFR_con_int_ *)params;
    double M = exp(lnM); // linear scale
    double z = vals.z_obs;
    double M2 = vals.Mval; // natural log scale
    double del1 = vals.delta1;
    double del2 = vals.delta2;
    double MassTurnover = vals.Mdrop;
    double Alpha_star = vals.pl_star;
    double Alpha_esc = vals.pl_esc;
    double Fstar10 = vals.frac_star;
    double Fesc10 = vals.frac_esc;
    double Mlim_Fstar = vals.LimitMass_Fstar;
    double Mlim_Fesc = vals.LimitMass_Fesc;

    double Fstar,Fesc;

    if (Alpha_star > 0. && M > Mlim_Fstar)
        Fstar = 1./Fstar10;
    else if (Alpha_star < 0. && M < Mlim_Fstar)
        Fstar = 1./Fstar10;
    else
        Fstar = pow(M/1e10,Alpha_star);

    if (Alpha_esc > 0. && M > Mlim_Fesc)
        Fesc = 1./Fesc10;
    else if (Alpha_esc < 0. && M < Mlim_Fesc)
        Fesc = 1./Fesc10;
    else 
        Fesc = pow(M/1e10,Alpha_esc);

    return M*exp(-MassTurnover/M)*Fstar*Fesc*dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI);
}

#ifdef MINI_HALO
double dNion_ConditionallnMm(double lnM, void *params) {
    struct parameters_gsl_SFR_con_intm_ vals = *(struct parameters_gsl_SFR_con_intm_ *)params;
    double M = exp(lnM); // linear scale
    double z = vals.z_obs;
    double M2 = vals.Mval; // natural log scale
    double del1 = vals.delta1;
    double del2 = vals.delta2;
    double Alpha_star = vals.pl_star;
    double MassTurnoverm = vals.Mdropm1;
    double Mcrit_atom = vals.Mdropm2;
    double Fstar10m = vals.frac_starm;
    double Mlim_Fstarm = vals.LimitMass_Fstarm;

    double Fstarm;

    if (Alpha_star > 0. && M > Mlim_Fstarm)
        Fstarm = 1./Fstar10m;
    else if (Alpha_star < 0. && M < Mlim_Fstarm)
        Fstarm = 1./Fstar10m;
    else
        Fstarm = pow(M/1e10,Alpha_star);

    return M*exp(-MassTurnoverm / M) * exp(-M / Mcrit_atom) * Fstarm *
           dNdM_conditional_second(z, log(M), M2, del1, del2) / sqrt(2. * PI);
}
#endif

double Nion_ConditionalM(double z, double M1, double M2, double delta1, double delta2, double MassTurnover, double Alpha_star, double Alpha_esc, double Fstar10, double Fesc10, double Mlim_Fstar, double Mlim_Fesc) {
    double result, error, lower_limit, upper_limit;
    gsl_function F;
    double rel_tol = 0.005; //<- relative tolerance
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

    struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
        .z_obs = z,
        .Mval = M2,
        .delta1 = delta1,
        .delta2 = delta2,
        .Mdrop = MassTurnover,
        .pl_star = Alpha_star,
        .pl_esc = Alpha_esc,
        .frac_star = Fstar10,
        .frac_esc = Fesc10,
        .LimitMass_Fstar = Mlim_Fstar,
        .LimitMass_Fesc = Mlim_Fesc
    };

    F.function = &dNion_ConditionallnM;
    F.params = &parameters_gsl_SFR_con;
    lower_limit = M1;
    upper_limit = M2;

    if(delta2 >= delta1)
      result = 1.;
    else {
      gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                          1000, GSL_INTEG_GAUSS61, w, &result, &error);
      gsl_integration_workspace_free (w);
    }
    return result;

}

#ifdef MINI_HALO
double Nion_ConditionalMm(double z, double M1, double M2, double delta1, double delta2, double Alpha_star, double MassTurnoverm, double Mcrit_atom, double Fstar10m, double Mlim_Fstarm) {
    double result, error, lower_limit, upper_limit;
    gsl_function F;
    double rel_tol = 0.005; //<- relative tolerance
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

    struct parameters_gsl_SFR_con_intm_ parameters_gsl_SFR_con = {
        .z_obs = z,
        .Mval = M2,
        .delta1 = delta1,
        .delta2 = delta2,
        .pl_star = Alpha_star,
        .Mdropm1 = MassTurnoverm,
        .Mdropm2 = Mcrit_atom,
        .frac_starm = Fstar10m,
        .LimitMass_Fstarm = Mlim_Fstarm
    };

    F.function = &dNion_ConditionallnMm;
    F.params = &parameters_gsl_SFR_con;
    lower_limit = M1;
    upper_limit = M2;

    if(delta2 >= delta1)
      result = 1.;
    else {
      gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                          1000, GSL_INTEG_GAUSS61, w, &result, &error);
      gsl_integration_workspace_free (w);
    }
    return result;

}
#endif

// Set up interpolation table for the number of IGM ionizing photons per baryon and initialise interploation.
// This function includes conditional mass function over the density field at a given redshift.
// Split density range into low and high density, since to increase accuray in high density region.
// see eq. (17) in Park et al. 2018
// Note this function DO NOT include ION_EFF_FACTOR, defined as f_{\ast,10}*f_{esc,10}*N_{\gamma/b}. 
// The actual Nion in eq. (17) is the value multiplied by ION_EFF_FACTOR.
#ifdef INHOMO_FEEDBACK
void initialise_Nion_spline(float z, float Mmax, float Mmin, float MMinTurnover, float Alpha_star, float Alpha_esc, float Fstar10, float Fesc10, float Mlim_Fstar, float Mlim_Fesc)
#else //INHOMO_FEEDBACK
void initialise_Nion_spline(float z, float Mmax, float Mmin, float MassTurnover, float Alpha_star, float Alpha_esc, float Fstar10, float Fesc10, float Mlim_Fstar, float Mlim_Fesc)
#endif //INHOMO_FEEDBACK
{
    double overdense_val;
    double overdense_large_high = Deltac, overdense_large_low = 1.5;
    double overdense_small_high = 1.5, overdense_small_low = -1. + 9e-8;

    int i;
    NionLow_spline_acc = gsl_interp_accel_alloc ();
#ifdef INHOMO_FEEDBACK
    double MassTurnover;
    double LogMassTurnover_high = 10. + 9e-8;
    double LogMassTurnover_low = log10(MMinTurnover) - 9e-8;

    int j;
    NionLow_spline_acc_Mturn = gsl_interp_accel_alloc ();
    NionLow_spline = gsl_spline2d_alloc (gsl_interp2d_bicubic, NSFR_low, NMTURN);
#pragma omp parallel shared(overdense_small_low, overdense_small_high, log10_overdense_spline_SFR, log10_Mturn_spline_SFR, log10_Mturn_spline_SFR_float, LogMassTurnover_low, LogMassTurnover_high, log10_Nion_spline, z,Mmax, Alpha_star,Alpha_esc,Fstar10,Fesc10,Mlim_Fstar,Mlim_Fesc) private(i, overdense_val, j, MassTurnover)
#else //INHOMO_FEEDBACK
    NionLow_spline = gsl_spline_alloc (gsl_interp_cspline, NSFR_low);
#pragma omp parallel shared(overdense_small_low, overdense_small_high, log10_overdense_spline_SFR, log10_Nion_spline, z,Mmax, MassTurnover,Alpha_star,Alpha_esc,Fstar10,Fesc10,Mlim_Fstar,Mlim_Fesc) private(i, overdense_val)
#endif //INHOMO_FEEDBACK
{
#pragma omp for
    for (i=0; i<NSFR_low; i++){
        overdense_val = log10(1. + overdense_small_low) + (double)i/((double)NSFR_low-1.)*(log10(1.+overdense_small_high)-log10(1.+overdense_small_low));
        log10_overdense_spline_SFR[i] = overdense_val;

#ifdef INHOMO_FEEDBACK
        for (j=0; j<NMTURN; j++){
          log10_Mturn_spline_SFR[j] = LogMassTurnover_low + (double)j/((double)NMTURN-1.)*(LogMassTurnover_high-LogMassTurnover_low);
          log10_Mturn_spline_SFR_float[j] = (float)(log10_Mturn_spline_SFR[j]);
          MassTurnover =  pow(10., log10_Mturn_spline_SFR[j]);
          log10_Nion_spline[i+j*NSFR_low] = log10(GaussLegendreQuad_Nion(NGL_SFR,z,log(Mmax),Deltac,pow(10.,log10_overdense_spline_SFR[i])-1.,MassTurnover,Alpha_star,Alpha_esc,Fstar10,Fesc10,Mlim_Fstar,Mlim_Fesc));
          if(log10_Nion_spline[i+j*NSFR_low] < -40.)
            log10_Nion_spline[i+j*NSFR_low] = -40.;
        }
#else //INHOMO_FEEDBACK
        log10_Nion_spline[i] = log10(GaussLegendreQuad_Nion(NGL_SFR,z,log(Mmax),Deltac,pow(10.,log10_overdense_spline_SFR[i])-1.,MassTurnover,Alpha_star,Alpha_esc,Fstar10,Fesc10,Mlim_Fstar,Mlim_Fesc));
        if(log10_Nion_spline[i] < -40.)
            log10_Nion_spline[i] = -40.;
#endif
    }
}
#ifdef INHOMO_FEEDBACK
    gsl_spline2d_init(NionLow_spline, log10_overdense_spline_SFR, log10_Mturn_spline_SFR, log10_Nion_spline, NSFR_low, NMTURN);
#else //INHOMO_FEEDBACK
    gsl_spline_init(NionLow_spline, log10_overdense_spline_SFR, log10_Nion_spline, NSFR_low);
#endif //INHOMO_FEEDBACK

#ifdef INHOMO_FEEDBACK
#pragma omp parallel shared(Overdense_spline_SFR, overdense_large_low, overdense_large_high, log10_Mturn_spline_SFR, Nion_spline, z,Mmin,Mmax, Alpha_star,Alpha_esc,Fstar10,Fesc10,Mlim_Fstar,Mlim_Fesc) private(i, j, MassTurnover)
#else
#pragma omp parallel shared(Overdense_spline_SFR, overdense_large_low, overdense_large_high, Nion_spline, z,Mmin,Mmax, MassTurnover,Alpha_star,Alpha_esc,Fstar10,Fesc10,Mlim_Fstar,Mlim_Fesc) private(i)
#endif
{
#pragma omp for
    for(i=0;i<NSFR_high;i++) {
        Overdense_spline_SFR[i] = overdense_large_low + (float)i/((float)NSFR_high-1.)*(overdense_large_high - overdense_large_low);

#ifdef INHOMO_FEEDBACK
        for (j=0; j<NMTURN; j++){
          //log10_Mturn_spline_SFR[j] = LogMassTurnover_low + (double)j/((double)NMTURN-1.)*(LogMassTurnover_high-LogMassTurnover_low);
          //log10_Mturn_spline_SFR_float[j] = (float)(log10_Mturn_spline_SFR[j]);
          MassTurnover =  pow(10., log10_Mturn_spline_SFR[j]);
          Nion_spline[i+j*NSFR_high] = Nion_ConditionalM(z,log(Mmin),log(Mmax),Deltac,Overdense_spline_SFR[i],MassTurnover,Alpha_star,Alpha_esc,Fstar10,Fesc10,Mlim_Fstar,Mlim_Fesc);
          if(Nion_spline[i+j*NSFR_high]<0.)
            Nion_spline[i+j*NSFR_high]=pow(10.,-40.0);
        }
#else //INHOMO_FEEDBACK
        Nion_spline[i] = Nion_ConditionalM(z,log(Mmin),log(Mmax),Deltac,Overdense_spline_SFR[i],MassTurnover,Alpha_star,Alpha_esc,Fstar10,Fesc10,Mlim_Fstar,Mlim_Fesc);
        if(Nion_spline[i]<0.)
          Nion_spline[i]=pow(10.,-40.0);
#endif //INHOMO_FEEDBACK
    }
}
#ifdef INHOMO_FEEDBACK
    spline2d(Overdense_spline_SFR,log10_Mturn_spline_SFR_float, Nion_spline,NSFR_high, NMTURN,second_derivs_Nion);
#else //INHOMO_FEEDBACK
    spline(Overdense_spline_SFR-1,Nion_spline-1,NSFR_high,0,0,second_derivs_Nion-1);
#endif //INHOMO_FEEDBACK
}

#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
void initialise_Nion_splinem(float z, float Mmax, float Mmin, float Alpha_star, float MminTurnoverm, float Mcrit_atom, float Fstar10m, float Mlim_Fstarm)
#else //INHOMO_FEEDBACK
void initialise_Nion_splinem(float z, float Mmax, float Mmin, float Alpha_star, float MassTurnoverm, float Mcrit_atom, float Fstar10m, float Mlim_Fstarm)
#endif //INHOMO_FEEDBACK
{
    //double overdense_val;
    //double overdense_large_high = Deltac, overdense_large_low = 1.5;
   // double overdense_small_high = 1.5, overdense_small_low = -1. + 9e-8;
    int i;
    NionLow_spline_accm = gsl_interp_accel_alloc ();
#ifdef INHOMO_FEEDBACK
    double MassTurnoverm;
    double LogMassTurnover_high = 10. + 9e-8;
    double LogMassTurnover_low = log10(MminTurnoverm) - 9e-8;

    int j;
    NionLow_spline_accm_Mturn = gsl_interp_accel_alloc ();
    NionLow_splinem = gsl_spline2d_alloc (gsl_interp2d_bicubic, NSFR_low, NMTURN);
#else //INHOMO_FEEDBACK
    NionLow_splinem = gsl_spline_alloc (gsl_interp_cspline, NSFR_low);
#endif //INHOMO_FEEDBACK

#ifdef INHOMO_FEEDBACK
#pragma omp parallel shared(log10_Mturn_spline_SFRm, log10_Mturn_spline_SFRm_float, LogMassTurnover_low, LogMassTurnover_high, log10_Nion_splinem, z,Mmax, log10_overdense_spline_SFR, Alpha_star, Mcrit_atom,Fstar10m,Mlim_Fstarm) private(i, j, MassTurnoverm)
#else
#pragma omp parallel shared(log10_Nion_splinem, z,Mmax, log10_overdense_spline_SFR, Alpha_star, MassTurnoverm, Mcrit_atom,Fstar10m,Mlim_Fstarm) private(i)
#endif
{
#pragma omp for
    for (i=0; i<NSFR_low; i++){
        //overdense_val = log10(1. + overdense_small_low) + (double)i/((double)NSFR_low-1.)*(log10(1.+overdense_small_high)-log10(1.+overdense_small_low));
        //log10_overdense_spline_SFR[i] = overdense_val;
        
#ifdef INHOMO_FEEDBACK
        for (j=0; j<NMTURN; j++){
          log10_Mturn_spline_SFRm[j] = LogMassTurnover_low + (double)j/((double)NMTURN-1.)*(LogMassTurnover_high-LogMassTurnover_low);
          log10_Mturn_spline_SFRm_float[j] = (float)(log10_Mturn_spline_SFRm[j]);
          MassTurnoverm =  pow(10., log10_Mturn_spline_SFRm[j]);
          log10_Nion_splinem[i+j*NSFR_low] = log10(GaussLegendreQuad_Nionm(NGL_SFR,z,log(Mmax),Deltac,pow(10.,log10_overdense_spline_SFR[i])-1.,Alpha_star,MassTurnoverm,Mcrit_atom,Fstar10m,Mlim_Fstarm));
          if(log10_Nion_splinem[i+j*NSFR_low] < -40.)
            log10_Nion_splinem[i+j*NSFR_low] = -40.;
        }
#else //INHOMO_FEEDBACK
        log10_Nion_splinem[i] = log10(GaussLegendreQuad_Nionm(NGL_SFR,z,log(Mmax),Deltac,pow(10.,log10_overdense_spline_SFR[i])-1.,Alpha_star,MassTurnoverm,Mcrit_atom,Fstar10m,Mlim_Fstarm));
        if(log10_Nion_splinem[i] < -40.)
            log10_Nion_splinem[i] = -40.;
#endif
    }
}
#ifdef INHOMO_FEEDBACK
    gsl_spline2d_init(NionLow_splinem, log10_overdense_spline_SFR, log10_Mturn_spline_SFRm, log10_Nion_splinem, NSFR_low, NMTURN);
#else //INHOMO_FEEDBACK
    gsl_spline_init(NionLow_splinem, log10_overdense_spline_SFR, log10_Nion_splinem, NSFR_low);
#endif //INHOMO_FEEDBACK

#ifdef INHOMO_FEEDBACK
#pragma omp parallel shared(log10_Mturn_spline_SFRm, Nion_splinem, z, Mmin, Mmax, Overdense_spline_SFR, Alpha_star, Mcrit_atom,Fstar10m,Mlim_Fstarm) private(i, j, MassTurnoverm)
#else
#pragma omp parallel shared(Nion_splinem, z, Mmin, Mmax, Overdense_spline_SFR, Alpha_star,MassTurnoverm, Mcrit_atom,Fstar10m,Mlim_Fstarm) private(i)
#endif
{
#pragma omp for
    for(i=0;i<NSFR_high;i++) {
        //Overdense_spline_SFR[i] = overdense_large_low + (float)i/((float)NSFR_high-1.)*(overdense_large_high - overdense_large_low);

#ifdef INHOMO_FEEDBACK
        for (j=0; j<NMTURN; j++){
          //log10_Mturn_spline_SFRm[j] = LogMassTurnover_low + (double)j/((double)NMTURN-1.)*(LogMassTurnover_high-LogMassTurnover_low);
          //log10_Mturn_spline_SFRm_float[j] = (float)(log10_Mturn_spline_SFRm[j]);
          MassTurnoverm =  pow(10., log10_Mturn_spline_SFRm[j]);
          Nion_splinem[i+j*NSFR_high] = Nion_ConditionalMm(z,log(Mmin),log(Mmax),Deltac,Overdense_spline_SFR[i],Alpha_star,MassTurnoverm,Mcrit_atom,Fstar10m,Mlim_Fstarm);
          if(Nion_splinem[i+j*NSFR_high]<0.)
            Nion_splinem[i+j*NSFR_high]=pow(10.,-40.0);
        }
#else //INHOMO_FEEDBACK
        Nion_splinem[i] = Nion_ConditionalMm(z,log(Mmin),log(Mmax),Deltac,Overdense_spline_SFR[i],Alpha_star,MassTurnoverm,Mcrit_atom,Fstar10m,Mlim_Fstarm);
        if(Nion_splinem[i]<0.)
          Nion_splinem[i]=pow(10.,-40.0);
#endif //INHOMO_FEEDBACK
    }
}
#ifdef INHOMO_FEEDBACK
    spline2d(Overdense_spline_SFR,log10_Mturn_spline_SFRm_float, Nion_splinem,NSFR_high, NMTURN,second_derivs_Nionm);
#else //INHOMO_FEEDBACK
    spline(Overdense_spline_SFR-1,Nion_splinem-1,NSFR_high,0,0,second_derivs_Nionm-1);
#endif //INHOMO_FEEDBACK
}
#endif

// Find the number of IGM ionizing photons per baryon at a given overdensity using interpolation.
#ifdef INHOMO_FEEDBACK
void Nion_Spline_density(float Overdensity, float log10_M_MINa, float *splined_value)
#else //INHOMO_FEEDBACK
void Nion_Spline_density(float Overdensity, float *splined_value)
#endif //INHOMO_FEEDBACK
{
    int i;
    float returned_value;

    if(Overdensity<1.5) {
        if(Overdensity<-1.) {
            returned_value = 0;
        }
        else {
#ifdef INHOMO_FEEDBACK
            if (log10_M_MINa>10)
              returned_value = gsl_spline2d_eval(NionLow_spline, log10(Overdensity+1.), 10, NionLow_spline_acc, NionLow_spline_acc_Mturn);

            else
              returned_value = gsl_spline2d_eval(NionLow_spline, log10(Overdensity+1.), log10_M_MINa, NionLow_spline_acc, NionLow_spline_acc_Mturn);
#else //INHOMO_FEEDBACK
            returned_value = gsl_spline_eval(NionLow_spline, log10(Overdensity+1.), NionLow_spline_acc);
#endif //INHOMO_FEEDBACK
            returned_value = pow(10.,returned_value);
        }
    }
    else {
        if(Overdensity<0.99*Deltac) {
        // Usage of 0.99*Deltac arises due to the fact that close to the critical density, the collapsed fraction becomes a little unstable
        // However, such densities should always be collapsed, so just set f_coll to unity. 
        // Additionally, the fraction of points in this regime relative to the entire simulation volume is extremely small.
#ifdef INHOMO_FEEDBACK
          if (log10_M_MINa>10) 
            splint2d(Overdense_spline_SFR,log10_Mturn_spline_SFR_float,Nion_spline,second_derivs_Nion,NSFR_high,NMTURN,Overdensity,10,&(returned_value));
          else
            splint2d(Overdense_spline_SFR,log10_Mturn_spline_SFR_float,Nion_spline,second_derivs_Nion,NSFR_high,NMTURN,Overdensity,log10_M_MINa,&(returned_value));
#else //INHOMO_FEEDBACK
          splint(Overdense_spline_SFR-1,Nion_spline-1,second_derivs_Nion-1,NSFR_high,Overdensity,&(returned_value));
#endif //INHOMO_FEEDBACK
        }
        else {
            returned_value = 1.;
        }
    }
    if(returned_value > 1.) returned_value = 1.;
    *splined_value = returned_value > 0. ? returned_value : 1e-40;
}
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
void Nion_Spline_densitym(float Overdensity, float log10_M_MINm, float *splined_value)
#else //INHOMO_FEEDBACK
void Nion_Spline_densitym(float Overdensity, float *splined_value)
#endif //INHOMO_FEEDBACK
{
    int i;
    float returned_value;

    if(Overdensity<1.5) {
        if(Overdensity<-1.) {
            returned_value = 0;
        }
        else {
#ifdef INHOMO_FEEDBACK
            if (log10_M_MINm>10)
              returned_value = gsl_spline2d_eval(NionLow_splinem, log10(Overdensity+1.), 10, NionLow_spline_accm, NionLow_spline_accm_Mturn);
            else
              returned_value = gsl_spline2d_eval(NionLow_splinem, log10(Overdensity+1.), log10_M_MINm, NionLow_spline_accm, NionLow_spline_accm_Mturn);
#else //INHOMO_FEEDBACK
            returned_value = gsl_spline_eval(NionLow_splinem, log10(Overdensity+1.), NionLow_spline_accm);
#endif //INHOMO_FEEDBACK
            returned_value = pow(10.,returned_value);
        }
    }
    else {
        if(Overdensity<0.99*Deltac) {
        // Usage of 0.99*Deltac arises due to the fact that close to the critical density, the collapsed fraction becomes a little unstable
        // However, such densities should always be collapsed, so just set f_coll to unity. 
        // Additionally, the fraction of points in this regime relative to the entire simulation volume is extremely small.
#ifdef INHOMO_FEEDBACK
            if (log10_M_MINm>10)
              splint2d(Overdense_spline_SFR,log10_Mturn_spline_SFRm_float,Nion_splinem,second_derivs_Nionm,NSFR_high,NMTURN,Overdensity,10,&(returned_value));
            else
              splint2d(Overdense_spline_SFR,log10_Mturn_spline_SFRm_float,Nion_splinem,second_derivs_Nionm,NSFR_high,NMTURN,Overdensity,log10_M_MINm,&(returned_value));
#else //INHOMO_FEEDBACK
            splint(Overdense_spline_SFR-1,Nion_splinem-1,second_derivs_Nionm-1,NSFR_high,Overdensity,&(returned_value));
#endif //INHOMO_FEEDBACK
        }
        else {
            returned_value = 1.;
        }
    }
    if(returned_value > 1.) returned_value = 1.;
    *splined_value = returned_value > 0. ? returned_value : 1e-40;
}
#endif //MINI_HALO

#ifdef CONTEMPORANEOUS_DUTYCYCLE
double dDeltaNion_ST(double lnM, void *params){
    struct parameters_gsl_SFR_int_ vals = *(struct parameters_gsl_SFR_int_ *)params;

    double M = exp(lnM);
    float z = vals.z_obs;
    float zp = vals.zp_obs;
    double MassTurnover = vals.Mdrop;
    double Alpha_star = vals.pl_star;
    double Alpha_esc = vals.pl_esc;
    double Fstar10 = vals.frac_star;
    double Fesc10 = vals.frac_esc;
    double Mlim_Fstar = vals.LimitMass_Fstar;
    double Mlim_Fesc = vals.LimitMass_Fesc;

    double Fstar, Fesc;


    if (Alpha_star > 0. && M > Mlim_Fstar)
        Fstar = 1./Fstar10;
    else if (Alpha_star < 0. && M < Mlim_Fstar)
        Fstar = 1/Fstar10;
    else
        Fstar = pow(M/1e10,Alpha_star);

    if (Alpha_esc > 0. && M > Mlim_Fesc)
        Fesc = 1./Fesc10;
    else if (Alpha_esc < 0. && M < Mlim_Fesc)
        Fesc = 1./Fesc10;
    else
        Fesc = pow(M/1e10,Alpha_esc);
    return (dNdM_st(z,M) - dNdM_st(zp,M)) * M * M * exp(-MassTurnover/M) * Fstar * Fesc;
}

#ifdef MINI_HALO
double dDeltaNion_STm(double lnM, void *params){
    struct parameters_gsl_SFR_intm_ vals = *(struct parameters_gsl_SFR_intm_ *)params;

    double M = exp(lnM);
    float z = vals.z_obs;
    float zp = vals.zp_obs;
    double MassTurnover = vals.Mdrop;
    double MassTurnover2 = vals.Mdrop_2;
    double Alpha_star = vals.pl_star;
    double Fstar10 = vals.frac_star;
    double Mlim_Fstar = vals.LimitMass_Fstar;

    double Fstar;

    if (Alpha_star > 0. && M > Mlim_Fstar)
        Fstar = 1./Fstar10;
    else if (Alpha_star < 0. && M < Mlim_Fstar)
        Fstar = 1/Fstar10;
    else
        Fstar = pow(M/1e10,Alpha_star);

    return (dNdM_st(z,M) - dNdM_st(zp,M)) * M * M * exp(-MassTurnover/M) * exp(-M/MassTurnover2) * Fstar;
}
#endif

double DeltaNion_ST(double z, double zp, double M_Min, double MassTurnover, double Alpha_star, double Alpha_esc, double Fstar10, double Fesc10, double Mlim_Fstar, double Mlim_Fesc){
    double result, error, lower_limit, upper_limit;
    gsl_function F;
    double rel_tol = 0.001; //<- relative tolerance
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

    struct parameters_gsl_SFR_int_ parameters_gsl_SFR = {
        .z_obs = z,
        .zp_obs = zp,
        .Mdrop = MassTurnover,
        .pl_star = Alpha_star,
        .pl_esc = Alpha_esc,
        .frac_star = Fstar10,
        .frac_esc = Fesc10,
        .LimitMass_Fstar = Mlim_Fstar,
        .LimitMass_Fesc = Mlim_Fesc,
    };

    F.function = &dDeltaNion_ST;
    F.params = &parameters_gsl_SFR;
    lower_limit = log(M_Min);
    upper_limit = log(FMAX(1e16, M_Min*100));

    gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                        1000, GSL_INTEG_GAUSS61, w, &result, &error);
    gsl_integration_workspace_free (w);

    return result / (OMm*RHOcrit);

}
#ifdef MINI_HALO
double DeltaNion_STm(double z, double zp, double M_Min, double MassTurnover_lower, double MassTurnover_upper, double Alpha_star, double Fstar10, double Mlim_Fstar){
    double result, error, lower_limit, upper_limit;
    gsl_function F;
    double rel_tol = 0.001; //<- relative tolerance
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

    struct parameters_gsl_SFR_intm_ parameters_gsl_SFR = {
        .z_obs = z,
        .zp_obs = zp,
        .Mdrop = MassTurnover_lower,
        .Mdrop_2 = MassTurnover_upper,
        .pl_star = Alpha_star,
        .frac_star = Fstar10,
        .LimitMass_Fstar = Mlim_Fstar,
    };

    F.function = &dDeltaNion_STm;
    F.params = &parameters_gsl_SFR;
    lower_limit = log(M_Min);
    upper_limit = log(FMAX(1e16, M_Min*100));

    gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                        1000, GSL_INTEG_GAUSS61, w, &result, &error);
    gsl_integration_workspace_free (w);

    return result / (OMm*RHOcrit);

}
#endif

float GaussLegendreQuad_DeltaNion(int n, float z, float zp, float M2, float delta1, float delta2, float MassTurnover, float Alpha_star, float Alpha_esc, float Fstar10, float Fesc10, float Mlim_Fstar, float Mlim_Fesc) {
  if(delta2 > delta1)
    return 0.;

  struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
      .z_obs = z,
      .Mval = M2,
      .delta1 = delta1,
      .delta2 = delta2,
      .Mdrop = MassTurnover,
      .pl_star = Alpha_star,
      .pl_esc = Alpha_esc,
      .frac_star = Fstar10,
      .frac_esc = Fesc10,
      .LimitMass_Fstar = Mlim_Fstar,
      .LimitMass_Fesc = Mlim_Fesc
  };

  struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_conp = {
      .z_obs = zp,
      .Mval = M2,
      .delta1 = delta1,
      .delta2 = delta2,
      .Mdrop = MassTurnover,
      .pl_star = Alpha_star,
      .pl_esc = Alpha_esc,
      .frac_star = Fstar10,
      .frac_esc = Fesc10,
      .LimitMass_Fstar = Mlim_Fstar,
      .LimitMass_Fesc = Mlim_Fesc
  };

  //Performs the Gauss-Legendre quadrature.
  int i;
  float integrand = 0.;
  for(i=1; i<(n+1); i++)
    integrand += wi_SFR[i]*(Nion_ConditionallnM_GL(xi_SFR[i],parameters_gsl_SFR_con) - Nion_ConditionallnM_GL(xi_SFR[i],parameters_gsl_SFR_conp));

    if (integrand<=0)
      integrand = 1e-40;
  return integrand;

}
#ifdef MINI_HALO
float GaussLegendreQuad_DeltaNionm(int n, float z, float zp, float M2, float delta1, float delta2, float Alpha_star, float MassTurnoverm, float Mcrit_atom, float Fstar10m, float Mlim_Fstarm) {
    if(delta2 > delta1)
      return 0.;

    struct parameters_gsl_SFR_con_intm_ parameters_gsl_SFR_con = {
        .z_obs = z,
        .Mval = M2,
        .delta1 = delta1,
        .delta2 = delta2,
        .pl_star = Alpha_star,
        .Mdropm1 = MassTurnoverm,
        .Mdropm2 = Mcrit_atom,
        .frac_starm = Fstar10m,
        .LimitMass_Fstarm = Mlim_Fstarm
    };

    struct parameters_gsl_SFR_con_intm_ parameters_gsl_SFR_conp = {
        .z_obs = zp,
        .Mval = M2,
        .delta1 = delta1,
        .delta2 = delta2,
        .pl_star = Alpha_star,
        .Mdropm1 = MassTurnoverm,
        .Mdropm2 = Mcrit_atom,
        .frac_starm = Fstar10m,
        .LimitMass_Fstarm = Mlim_Fstarm
    };

    //Performs the Gauss-Legendre quadrature.
    int i;
    float integrand = 0.;
    for(i=1; i<(n+1); i++)
      integrand += wi_SFR[i]*(Nion_ConditionallnM_GLm(xi_SFR[i],parameters_gsl_SFR_con) - Nion_ConditionallnM_GLm(xi_SFR[i],parameters_gsl_SFR_conp));

    // TODO: Sometimes this happens, need to get back and check what is going on here
    if (integrand<=0)
      integrand = 1e-40;
    return integrand;

}
#endif

double DeltaNion_ConditionalM(double z, double zp, double M1, double M2, double delta1, double delta2, double MassTurnover, double Alpha_star, double Alpha_esc, double Fstar10, double Fesc10, double Mlim_Fstar, double Mlim_Fesc) {
    if(delta2 >= delta1)
      return 0.;

    double result, resultp, error, lower_limit, upper_limit;
    gsl_function F;
    double rel_tol = 0.005; //<- relative tolerance
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

    struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
        .z_obs = z,
        .Mval = M2,
        .delta1 = delta1,
        .delta2 = delta2,
        .Mdrop = MassTurnover,
        .pl_star = Alpha_star,
        .pl_esc = Alpha_esc,
        .frac_star = Fstar10,
        .frac_esc = Fesc10,
        .LimitMass_Fstar = Mlim_Fstar,
        .LimitMass_Fesc = Mlim_Fesc
    };

    struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_conp = {
        .z_obs = zp,
        .Mval = M2,
        .delta1 = delta1,
        .delta2 = delta2,
        .Mdrop = MassTurnover,
        .pl_star = Alpha_star,
        .pl_esc = Alpha_esc,
        .frac_star = Fstar10,
        .frac_esc = Fesc10,
        .LimitMass_Fstar = Mlim_Fstar,
        .LimitMass_Fesc = Mlim_Fesc
    };

    lower_limit = M1;
    upper_limit = M2;
    F.function = &dNion_ConditionallnM;

    F.params = &parameters_gsl_SFR_con;
    gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                        1000, GSL_INTEG_GAUSS61, w, &result, &error);
    F.params = &parameters_gsl_SFR_conp;
    gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                        1000, GSL_INTEG_GAUSS61, w, &resultp, &error);
    gsl_integration_workspace_free (w);
    
    return result - resultp;

}

#ifdef MINI_HALO
double DeltaNion_ConditionalMm(double z, double zp, double M1, double M2, double delta1, double delta2, double Alpha_star, double MassTurnoverm, double Mcrit_atom, double Fstar10m, double Mlim_Fstarm) {
    if(delta2 >= delta1)
      return 0.;

    double result, resultp, error, lower_limit, upper_limit;
    gsl_function F;
    double rel_tol = 0.005; //<- relative tolerance
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

    struct parameters_gsl_SFR_con_intm_ parameters_gsl_SFR_con = {
        .z_obs = z,
        .Mval = M2,
        .delta1 = delta1,
        .delta2 = delta2,
        .pl_star = Alpha_star,
        .Mdropm1 = MassTurnoverm,
        .Mdropm2 = Mcrit_atom,
        .frac_starm = Fstar10m,
        .LimitMass_Fstarm = Mlim_Fstarm
    };

    struct parameters_gsl_SFR_con_intm_ parameters_gsl_SFR_conp = {
        .z_obs = zp,
        .Mval = M2,
        .delta1 = delta1,
        .delta2 = delta2,
        .pl_star = Alpha_star,
        .Mdropm1 = MassTurnoverm,
        .Mdropm2 = Mcrit_atom,
        .frac_starm = Fstar10m,
        .LimitMass_Fstarm = Mlim_Fstarm
    };

    lower_limit = M1;
    upper_limit = M2;
    F.function = &dNion_ConditionallnMm;

    F.params = &parameters_gsl_SFR_con;
    gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                        1000, GSL_INTEG_GAUSS61, w, &result, &error);
    F.params = &parameters_gsl_SFR_conp;
    gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                        1000, GSL_INTEG_GAUSS61, w, &resultp, &error);
    gsl_integration_workspace_free (w);

    return result - resultp;

}
#endif

#ifdef INHOMO_FEEDBACK
void initialise_DeltaNion_spline(float z, float zp, float Mmax, float Mmin, float MminTurnover, float Alpha_star, float Alpha_esc, float Fstar10, float Fesc10, float Mlim_Fstar, float Mlim_Fesc)
#else //INHOMO_FEEDBACK
void initialise_DeltaNion_spline(float z, float zp, float Mmax, float Mmin, float MassTurnover, float Alpha_star, float Alpha_esc, float Fstar10, float Fesc10, float Mlim_Fstar, float Mlim_Fesc)
#endif //INHOMO_FEEDBACK
{
    double overdense_val;
    double overdense_large_high = Deltac, overdense_large_low = 1.5;
    double overdense_small_high = 1.5, overdense_small_low = -1. + 9e-8;

    int i;
    NionLow_spline_acc = gsl_interp_accel_alloc ();
#ifdef INHOMO_FEEDBACK
    double MassTurnover;
    double LogMassTurnover_high = 10. + 9e-8;
    double LogMassTurnover_low = log10(MminTurnover) - 9e-8;

    int j;
    NionLow_spline_acc_Mturn = gsl_interp_accel_alloc ();
    NionLow_spline = gsl_spline2d_alloc (gsl_interp2d_bicubic, NSFR_low, NMTURN);
#else //INHOMO_FEEDBACK
    NionLow_spline = gsl_spline_alloc (gsl_interp_cspline, NSFR_low);
#endif //INHOMO_FEEDBACK

#ifdef INHOMO_FEEDBACK
#pragma omp parallel shared(overdense_small_low, overdense_small_high, log10_overdense_spline_SFR, log10_Mturn_spline_SFR, log10_Mturn_spline_SFR_float, LogMassTurnover_low, LogMassTurnover_high, log10_Nion_spline, z,zp,Mmax, Alpha_star,Alpha_esc,Fstar10,Fesc10,Mlim_Fstar,Mlim_Fesc) private(i, overdense_val, j, MassTurnover)
#else
#pragma omp parallel shared(overdense_small_low, overdense_small_high, log10_overdense_spline_SFR, log10_Nion_spline, z,zp,Mmax,MassTurnover, Alpha_star,Alpha_esc,Fstar10,Fesc10,Mlim_Fstar,Mlim_Fesc) private(i, overdense_val)
#endif
{
#pragma omp for
    for (i=0; i<NSFR_low; i++){
        overdense_val = log10(1. + overdense_small_low) + (double)i/((double)NSFR_low-1.)*(log10(1.+overdense_small_high)-log10(1.+overdense_small_low));
        log10_overdense_spline_SFR[i] = overdense_val;

#ifdef INHOMO_FEEDBACK
        for (j=0; j<NMTURN; j++){
          log10_Mturn_spline_SFR[j] = LogMassTurnover_low + (double)j/((double)NMTURN-1.)*(LogMassTurnover_high-LogMassTurnover_low);
          log10_Mturn_spline_SFR_float[j] = (float)(log10_Mturn_spline_SFR[j]);
          MassTurnover =  pow(10., log10_Mturn_spline_SFR[j]);
          log10_Nion_spline[i+j*NSFR_low] = log10(GaussLegendreQuad_DeltaNion(NGL_SFR,z,zp,log(Mmax),Deltac,pow(10.,log10_overdense_spline_SFR[i])-1.,MassTurnover,Alpha_star,Alpha_esc,Fstar10,Fesc10,Mlim_Fstar,Mlim_Fesc));
          if(log10_Nion_spline[i+j*NSFR_low] < -40.)
              log10_Nion_spline[i+j*NSFR_low] = -40.;
        }
#else //INHOMO_FEEDBACK
        log10_Nion_spline[i] = log10(GaussLegendreQuad_DeltaNion(NGL_SFR,z,zp,log(Mmax),Deltac,pow(10.,log10_overdense_spline_SFR[i])-1.,MassTurnover,Alpha_star,Alpha_esc,Fstar10,Fesc10,Mlim_Fstar,Mlim_Fesc));
        if(log10_Nion_spline[i] < -40.)
            log10_Nion_spline[i] = -40.;
#endif //INHOMO_FEEDBACK
    }
}
#ifdef INHOMO_FEEDBACK
    gsl_spline2d_init(NionLow_spline, log10_overdense_spline_SFR, log10_Mturn_spline_SFR, log10_Nion_spline, NSFR_low, NMTURN);
#else //INHOMO_FEEDBACK
    gsl_spline_init(NionLow_spline, log10_overdense_spline_SFR, log10_Nion_spline, NSFR_low);
#endif //INHOMO_FEEDBACK

#ifdef INHOMO_FEEDBACK
#pragma omp parallel shared(Overdense_spline_SFR, overdense_large_low, overdense_large_high, log10_Mturn_spline_SFR, Nion_spline, z,zp,Mmin,Mmax,Alpha_star,Alpha_esc,Fstar10,Fesc10,Mlim_Fstar,Mlim_Fesc) private(i, j, MassTurnover)
#else
#pragma omp parallel shared(Overdense_spline_SFR, overdense_large_low, overdense_large_high, Nion_spline, z,zp,Mmin,Mmax,MassTurnover,Alpha_star,Alpha_esc,Fstar10,Fesc10,Mlim_Fstar,Mlim_Fesc) private(i)
#endif
{
#pragma omp for
    for(i=0;i<NSFR_high;i++) {
        Overdense_spline_SFR[i] = overdense_large_low + (float)i/((float)NSFR_high-1.)*(overdense_large_high - overdense_large_low);

#ifdef INHOMO_FEEDBACK
        for (j=0; j<NMTURN; j++){
            //log10_Mturn_spline_SFR[j] = LogMassTurnover_low + (double)j/((double)NMTURN-1.)*(LogMassTurnover_high-LogMassTurnover_low);
            //log10_Mturn_spline_SFR_float[j] = (float)(log10_Mturn_spline_SFR[j]);
            MassTurnover =  pow(10., log10_Mturn_spline_SFR[j]);
            Nion_spline[i+j*NSFR_high] = DeltaNion_ConditionalM(z,zp,log(Mmin),log(Mmax),Deltac,Overdense_spline_SFR[i],MassTurnover,Alpha_star,Alpha_esc,Fstar10,Fesc10,Mlim_Fstar,Mlim_Fesc);
            if(Nion_spline[i+j*NSFR_high]<0.)
                Nion_spline[i+j*NSFR_high]=pow(10.,-40.0);
        }
#else //INHOMO_FEEDBACK
        Nion_spline[i] = DeltaNion_ConditionalM(z,zp,log(Mmin),log(Mmax),Deltac,Overdense_spline_SFR[i],MassTurnover,Alpha_star,Alpha_esc,Fstar10,Fesc10,Mlim_Fstar,Mlim_Fesc);
        if(Nion_spline[i]<0.)
            Nion_spline[i]=pow(10.,-40.0);
#endif //INHOMO_FEEDBACK
    }
}
#ifdef INHOMO_FEEDBACK
    spline2d(Overdense_spline_SFR,log10_Mturn_spline_SFR_float, Nion_spline,NSFR_high, NMTURN,second_derivs_Nion);
#else //INHOMO_FEEDBACK
    spline(Overdense_spline_SFR-1,Nion_spline-1,NSFR_high,0,0,second_derivs_Nion-1);
#endif //INHOMO_FEEDBACK
}

#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
void initialise_DeltaNion_splinem(float z, float zp, float Mmax, float Mmin, float Alpha_star, float MminTurnoverm, float Mcrit_atom, float Fstar10m, float Mlim_Fstarm)
#else //INHOMO_FEEDBACK    
void initialise_DeltaNion_splinem(float z, float zp, float Mmax, float Mmin, float Alpha_star, float MassTurnoverm, float Mcrit_atom, float Fstar10m, float Mlim_Fstarm)
#endif //INHOMO_FEEDBACK
{
    //double overdense_val;
    //double overdense_large_high = Deltac, overdense_large_low = 1.5;
   // double overdense_small_high = 1.5, overdense_small_low = -1. + 9e-8;
    int i;
    NionLow_spline_accm = gsl_interp_accel_alloc ();
#ifdef INHOMO_FEEDBACK
    double MassTurnoverm;
    double LogMassTurnover_high = 10. + 9e-8;
    double LogMassTurnover_low = log10(MminTurnoverm) - 9e-8;

    int j;
    NionLow_spline_accm_Mturn = gsl_interp_accel_alloc ();
    NionLow_splinem = gsl_spline2d_alloc (gsl_interp2d_bicubic, NSFR_low, NMTURN);
#else //INHOMO_FEEDBACK
    NionLow_splinem = gsl_spline_alloc (gsl_interp_cspline, NSFR_low);
#endif //INHOMO_FEEDBACK

#ifdef INHOMO_FEEDBACK
#pragma omp parallel shared(log10_Mturn_spline_SFRm, log10_Mturn_spline_SFRm_float, LogMassTurnover_low, LogMassTurnover_high, log10_Nion_splinem, z,zp,Mmax, log10_overdense_spline_SFR, Alpha_star, Mcrit_atom,Fstar10m,Mlim_Fstarm) private(i, j, MassTurnoverm)
#else
#pragma omp parallel shared(log10_Nion_splinem, z,zp,Mmax,  log10_overdense_spline_SFR, Alpha_star,MassTurnoverm, Mcrit_atom,Fstar10m,Mlim_Fstarm) private(i)
#endif
{
#pragma omp for
    for (i=0; i<NSFR_low; i++){
        //overdense_val = log10(1. + overdense_small_low) + (double)i/((double)NSFR_low-1.)*(log10(1.+overdense_small_high)-log10(1.+overdense_small_low));
        //log10_overdense_spline_SFR[i] = overdense_val;

#ifdef INHOMO_FEEDBACK
        for (j=0; j<NMTURN; j++){
            log10_Mturn_spline_SFRm[j] = LogMassTurnover_low + (double)j/((double)NMTURN-1.)*(LogMassTurnover_high-LogMassTurnover_low);
            log10_Mturn_spline_SFRm_float[j] = (float)(log10_Mturn_spline_SFRm[j]);
            MassTurnoverm =  pow(10., log10_Mturn_spline_SFRm[j]);
            log10_Nion_splinem[i+j*NSFR_low] = log10(GaussLegendreQuad_DeltaNionm(NGL_SFR,z,zp,log(Mmax),Deltac,pow(10.,log10_overdense_spline_SFR[i])-1.,Alpha_star,MassTurnoverm,Mcrit_atom,Fstar10m,Mlim_Fstarm));
            if(log10_Nion_splinem[i+j*NSFR_low] < -40.)
              log10_Nion_splinem[i+j*NSFR_low] = -40.;
        }
#else //INHOMO_FEEDBACK
        log10_Nion_splinem[i] = log10(GaussLegendreQuad_DeltaNionm(NGL_SFR,z,zp,log(Mmax),Deltac,pow(10.,log10_overdense_spline_SFR[i])-1.,Alpha_star,MassTurnoverm,Mcrit_atom,Fstar10m,Mlim_Fstarm));
        if(log10_Nion_splinem[i] < -40.)
            log10_Nion_splinem[i] = -40.;
#endif
    }
}
#ifdef INHOMO_FEEDBACK
    gsl_spline2d_init(NionLow_splinem, log10_overdense_spline_SFR, log10_Mturn_spline_SFRm, log10_Nion_splinem, NSFR_low, NMTURN);   
#else //INHOMO_FEEDBACK
    gsl_spline_init(NionLow_splinem, log10_overdense_spline_SFR, log10_Nion_splinem, NSFR_low);
#endif //INHOMO_FEEDBACK


#ifdef INHOMO_FEEDBACK
#pragma omp parallel shared(Overdense_spline_SFR, log10_Mturn_spline_SFRm, Nion_spline, z,zp,Mmin,Mmax,Alpha_star,Mcrit_atom,Fstar10m,Mlim_Fstarm) private(i, j, MassTurnoverm)
#else
#pragma omp parallel shared(Overdense_spline_SFR, Nion_spline, z,zp,Mmin,Mmax,MassTurnoverm,Alpha_star,Mcrit_atom,Fstar10m,Mlim_Fstarm) private(i)
#endif
{
#pragma omp for
    for(i=0;i<NSFR_high;i++) {
        //Overdense_spline_SFR[i] = overdense_large_low + (float)i/((float)NSFR_high-1.)*(overdense_large_high - overdense_large_low);

#ifdef INHOMO_FEEDBACK
        for (j=0; j<NMTURN; j++){
          //log10_Mturn_spline_SFRm[j] = LogMassTurnover_low + (double)j/((double)NMTURN-1.)*(LogMassTurnover_high-LogMassTurnover_low);
          //log10_Mturn_spline_SFRm_float[j] = (float)(log10_Mturn_spline_SFRm[j]);
          MassTurnoverm =  pow(10., log10_Mturn_spline_SFRm[j]);
          Nion_splinem[i+j*NSFR_high] = DeltaNion_ConditionalMm(z,zp,log(Mmin),log(Mmax),Deltac,Overdense_spline_SFR[i],Alpha_star,MassTurnoverm,Mcrit_atom,Fstar10m,Mlim_Fstarm);
          if(Nion_splinem[i+j*NSFR_high]<0.)
            Nion_splinem[i+j*NSFR_high]=pow(10.,-40.0);
        }
#else //INHOMO_FEEDBACK
        Nion_splinem[i] = DeltaNion_ConditionalMm(z,zp,log(Mmin),log(Mmax),Deltac,Overdense_spline_SFR[i],Alpha_star,MassTurnoverm,Mcrit_atom,Fstar10m,Mlim_Fstarm);
        if(Nion_splinem[i]<0.)
            Nion_splinem[i]=pow(10.,-40.0);
#endif //INHOMO_FEEDBACK
    }
}
#ifdef INHOMO_FEEDBACK
    spline2d(Overdense_spline_SFR,log10_Mturn_spline_SFRm_float, Nion_splinem,NSFR_high, NMTURN,second_derivs_Nionm);
#else //INHOMO_FEEDBACK
    spline(Overdense_spline_SFR-1,Nion_splinem-1,NSFR_high,0,0,second_derivs_Nionm-1);
#endif //INHOMO_FEEDBACK
}
#endif //MINI_HALO

// Find the number of IGM ionizing photons per baryon at a given overdensity using interpolation.
#ifdef INHOMO_FEEDBACK
void DeltaNion_Spline_density(float Overdensity, float log10_M_MINa, float *splined_value)
#else //INHOMO_FEEDBACK
void DeltaNion_Spline_density(float Overdensity, float *splined_value)
#endif //INHOMO_FEEDBACK
{
    float returned_value;

    if(Overdensity<1.5) {
        if(Overdensity<-1.) {
            returned_value = 0;
        }
        else {
#ifdef INHOMO_FEEDBACK
            if (log10_M_MINa>10)
                returned_value = gsl_spline2d_eval(NionLow_spline, log10(Overdensity+1), 10, NionLow_spline_acc, NionLow_spline_acc_Mturn);
            else
                returned_value = gsl_spline2d_eval(NionLow_spline, log10(Overdensity+1), log10_M_MINa, NionLow_spline_acc, NionLow_spline_acc_Mturn);
#else //INHOMO_FEEDBACK
            returned_value = gsl_spline_eval(NionLow_spline, log10(Overdensity+1.), NionLow_spline_acc);
#endif //INHOMO_FEEDBACK
            returned_value = pow(10.,returned_value);
        }
    }
    else {
        if(Overdensity<0.99*Deltac) {
        // Usage of 0.99*Deltac arises due to the fact that close to the critical density, the collapsed fraction becomes a little unstable
        // However, such densities should always be collapsed, so just set f_coll to unity. 
        // Additionally, the fraction of points in this regime relative to the entire simulation volume is extremely small.
#ifdef INHOMO_FEEDBACK
            if (log10_M_MINa>10)
                splint2d(Overdense_spline_SFR,log10_Mturn_spline_SFR_float,Nion_spline,second_derivs_Nion,NSFR_high,NMTURN,Overdensity,10,&(returned_value));
            else
                splint2d(Overdense_spline_SFR,log10_Mturn_spline_SFR_float,Nion_spline,second_derivs_Nion,NSFR_high,NMTURN,Overdensity,log10_M_MINa,&(returned_value));
#else //INHOMO_FEEDBACK
            splint(Overdense_spline_SFR-1,Nion_spline-1,second_derivs_Nion-1,NSFR_high,Overdensity,&(returned_value));
#endif //INHOMO_FEEDBACK
        }
        else {
            returned_value = 0.;
        }
    }
    *splined_value = returned_value > 0. ? returned_value : 1e-40;
}

#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
void DeltaNion_Spline_densitym(float Overdensity, float log10_M_MINm, float *splined_value)
#else //INHOMO_FEEDBACK
void DeltaNion_Spline_densitym(float Overdensity, float *splined_value)
#endif //INHOMO_FEEDBACK
{
    float returned_value;

    if(Overdensity<1.5) {
        if(Overdensity<-1.) {
            returned_value = 0;
        }
        else {
#ifdef INHOMO_FEEDBACK
            if (log10_M_MINm>10)
                returned_value = gsl_spline2d_eval(NionLow_splinem, log10(Overdensity+1.), 10, NionLow_spline_accm, NionLow_spline_accm_Mturn);
            else
                returned_value = gsl_spline2d_eval(NionLow_splinem, log10(Overdensity+1.), log10_M_MINm, NionLow_spline_accm, NionLow_spline_accm_Mturn);
#else //INHOMO_FEEDBACK
            returned_value = gsl_spline_eval(NionLow_splinem, log10(Overdensity+1.), NionLow_spline_accm);
#endif //INHOMO_FEEDBACK
            returned_value = pow(10.,returned_value);
        }
    }
    else {
        if(Overdensity<0.99*Deltac) {
        // Usage of 0.99*Deltac arises due to the fact that close to the critical density, the collapsed fraction becomes a little unstable
        // However, such densities should always be collapsed, so just set f_coll to unity. 
        // Additionally, the fraction of points in this regime relative to the entire simulation volume is extremely small.
#ifdef INHOMO_FEEDBACK
            if (log10_M_MINm>10)
                splint2d(Overdense_spline_SFR,log10_Mturn_spline_SFRm_float,Nion_splinem,second_derivs_Nionm,NSFR_high,NMTURN,Overdensity,10,&(returned_value));
            else
                splint2d(Overdense_spline_SFR,log10_Mturn_spline_SFRm_float,Nion_splinem,second_derivs_Nionm,NSFR_high,NMTURN,Overdensity,log10_M_MINm,&(returned_value));
#else //INHOMO_FEEDBACK
            splint(Overdense_spline_SFR-1,Nion_splinem-1,second_derivs_Nionm-1,NSFR_high,Overdensity,&(returned_value));
#endif //INHOMO_FEEDBACK
        }
        else {
            returned_value = 0.;
        }
    }
    *splined_value = returned_value > 0. ? returned_value : 1e-40;
}
#endif //MINI_HALO
#endif //CONTEMPORANEOUS_DUTYCYCLE

// Set up interpolation table for the mean number of IGM ionizing photons per baryon and initialise interploation.
// compute 'Nion_ST' corresponding to an array of redshift.
void initialise_Nion_ST_spline(int Nbin, double z_val[], float Mmin, double M_MINa_interp_table[], float Alpha_star, float Alpha_esc, float Fstar10, float Fesc10){
    int i;

    float Mlim_Fstar = Mass_limit_bisection(Mmin, 1e16, Alpha_star, Fstar10);
    float Mlim_Fesc  = Mass_limit_bisection(Mmin, 1e16, Alpha_esc, Fesc10);

    Nion_z_spline_acc = gsl_interp_accel_alloc ();
    Nion_z_spline = gsl_spline_alloc (gsl_interp_cspline, Nbin);
#pragma omp parallel shared(Nbin, Nion_z_val, z_val, Mmin, M_MINa_interp_table,  Alpha_star, Alpha_esc, Fstar10, Fesc10, Mlim_Fstar, Mlim_Fesc) private(i)
{
#pragma omp for    
    for (i=0; i<Nbin; i++){
      Nion_z_val[i] = Nion_ST(z_val[i], Mmin, M_MINa_interp_table[i], Alpha_star, Alpha_esc, Fstar10, Fesc10, Mlim_Fstar, Mlim_Fesc);
	}
}
    gsl_spline_init(Nion_z_spline, z_val, Nion_z_val, Nbin);
}

#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
void initialise_Nion_ST_splinem(int Nbin, double z_val[], float Mmin, double log10_Mturn[], double Mcrit_atom_interp_table[], float Alpha_star, float Fstar10m)
#else
void initialise_Nion_ST_splinem(int Nbin, double z_val[], float Mmin, double M_MINm_interp_table[], double Mcrit_atom_interp_table[], float Alpha_star, float Fstar10m)
#endif
{
    int i;

    float Mlim_Fstarm = Mass_limit_bisection(Mmin, 1e16, Alpha_star, Fstar10m);

    Nion_z_spline_accm = gsl_interp_accel_alloc ();
#ifdef INHOMO_FEEDBACK
    int j;
    double MassTurnover;

    Nion_z_spline_accm_Mturn = gsl_interp_accel_alloc ();
    Nion_z_splinem = gsl_spline2d_alloc (gsl_interp2d_bicubic, Nbin, NMTURN);
#else
    Nion_z_splinem = gsl_spline_alloc (gsl_interp_cspline, Nbin);
#endif
#ifdef INHOMO_FEEDBACK
#pragma omp parallel shared(Nbin, Nion_z_valm, z_val, Mmin, Mcrit_atom_interp_table, Alpha_star, Fstar10m, Mlim_Fstarm) private(i, j, MassTurnover)
#else
#pragma omp parallel shared(Nbin, Nion_z_valm, z_val, Mmin, M_MINm_interp_table, Mcrit_atom_interp_table, Alpha_star, Fstar10m, Mlim_Fstarm) private(i)
#endif
{
#pragma omp for
    for (i=0; i<Nbin; i++){
#ifdef INHOMO_FEEDBACK
        for (j=0; j<NMTURN; j++){
            MassTurnover =  pow(10., log10_Mturn[j]);
            Nion_z_valm[i+j*Nbin] = Nion_STm(z_val[i], Mmin, MassTurnover, Mcrit_atom_interp_table[i], Alpha_star, Fstar10m, Mlim_Fstarm);
        }
#else
        Nion_z_valm[i] = Nion_STm(z_val[i], Mmin, M_MINm_interp_table[i], Mcrit_atom_interp_table[i], Alpha_star, Fstar10m, Mlim_Fstarm);
#endif
    }
}
#ifdef INHOMO_FEEDBACK
    gsl_spline2d_init(Nion_z_splinem, z_val,log10_Mturn, Nion_z_valm, Nbin, NMTURN);
#else
    gsl_spline_init(Nion_z_splinem, z_val, Nion_z_valm, Nbin);
#endif
}
#endif

void Nion_ST_z(float z, float *splined_value){
    float returned_value;

    returned_value = gsl_spline_eval(Nion_z_spline, z, Nion_z_spline_acc);
    *splined_value = returned_value > 0. ? returned_value : 1e-40;
}

#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
void Nion_ST_zm(float z, float logM_MINm_ave, float *splined_value)
#else
void Nion_ST_zm(float z, float *splined_value)
#endif
{
    float returned_value;

#ifdef INHOMO_FEEDBACK
    if (logM_MINm_ave<5)
      returned_value = gsl_spline2d_eval(Nion_z_splinem, z, 5, Nion_z_spline_accm, Nion_z_spline_accm_Mturn);
    else if (logM_MINm_ave>10)
      returned_value = gsl_spline2d_eval(Nion_z_splinem, z, 10, Nion_z_spline_accm, Nion_z_spline_accm_Mturn);
    else
      returned_value = gsl_spline2d_eval(Nion_z_splinem, z, logM_MINm_ave, Nion_z_spline_accm, Nion_z_spline_accm_Mturn);
#else
    returned_value = gsl_spline_eval(Nion_z_splinem, z, Nion_z_spline_accm);
#endif
    *splined_value = returned_value > 0. ? returned_value : 1e-40;
}
#endif

void initialise_SFRD_ST_spline(int Nbin, double z_val[], float Mmin, double M_MINa_interp_table[], float Alpha_star, float Fstar10){
    int i;
    
    float Mlim_Fstar = Mass_limit_bisection(Mmin, 1e16, Alpha_star, Fstar10);

    SFRD_ST_z_spline_acc = gsl_interp_accel_alloc ();
    SFRD_ST_z_spline = gsl_spline_alloc (gsl_interp_cspline, Nbin);
#pragma omp parallel shared(Nbin, SFRD_val, z_val, Mmin, M_MINa_interp_table, Alpha_star, Fstar10, Mlim_Fstar) private(i)
{
#pragma omp for
    for (i=0; i<Nbin; i++)
        SFRD_val[i] = Nion_ST(z_val[i], Mmin, M_MINa_interp_table[i], Alpha_star, 0., Fstar10, 1.,Mlim_Fstar,0.);
}
    gsl_spline_init(SFRD_ST_z_spline, z_val, SFRD_val, Nbin);
}

#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
void initialise_SFRD_ST_splinem(int Nbin, double z_val[], float Mmin, double log10_Mturn[], double Mcrit_atom_interp_table[], float Alpha_star, float Fstar10m)
#else
void initialise_SFRD_ST_splinem(int Nbin, double z_val[], float Mmin, double M_MINm_interp_table[], double Mcrit_atom_interp_table[], float Alpha_star, float Fstar10m)
#endif
{
    int i;

    float Mlim_Fstarm = Mass_limit_bisection(Mmin, 1e16, Alpha_star, Fstar10m);

    SFRD_ST_z_spline_accm = gsl_interp_accel_alloc ();
#ifdef INHOMO_FEEDBACK
    int j;
    double MassTurnover;

    SFRD_ST_z_spline_accm_Mturn = gsl_interp_accel_alloc ();
    SFRD_ST_z_splinem = gsl_spline2d_alloc (gsl_interp2d_bicubic, Nbin, NMTURN);
#else
    SFRD_ST_z_splinem = gsl_spline_alloc (gsl_interp_cspline, Nbin);
#endif
#ifdef INHOMO_FEEDBACK
#pragma omp parallel shared(Nbin, log10_Mturn, SFRD_valm, z_val, Mmin, Mcrit_atom_interp_table, Alpha_star, Fstar10m, Mlim_Fstarm) private(i, MassTurnover, j)
#else
#pragma omp parallel shared(Nbin, SFRD_valm, z_val, Mmin, M_MINm_interp_table, Mcrit_atom_interp_table, Alpha_star, Fstar10m, Mlim_Fstarm) private(i)
#endif
{
//#pragma omp for
    for (i=0; i<Nbin; i++){
#ifdef INHOMO_FEEDBACK
        for (j=0; j<NMTURN; j++){
            MassTurnover =  pow(10., log10_Mturn[j]);
            SFRD_valm[i+j*Nbin] = Nion_STm(z_val[i], Mmin, MassTurnover, Mcrit_atom_interp_table[i], Alpha_star, Fstar10m, Mlim_Fstarm);
        }
#else
        SFRD_valm[i] = Nion_STm(z_val[i], Mmin, M_MINm_interp_table[i], Mcrit_atom_interp_table[i], Alpha_star, Fstar10m, Mlim_Fstarm);
#endif
    }
}
#ifdef INHOMO_FEEDBACK
    gsl_spline2d_init(SFRD_ST_z_splinem, z_val, log10_Mturn, SFRD_valm, Nbin, NMTURN);
#else
    gsl_spline_init(SFRD_ST_z_splinem, z_val, SFRD_valm, Nbin);
#endif
}
#endif

void SFRD_ST_z(float z, float *splined_value){
    float returned_value;

    returned_value = gsl_spline_eval(SFRD_ST_z_spline, z, SFRD_ST_z_spline_acc);
    *splined_value = returned_value > 0. ? returned_value : 1e-40;
}

#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
void SFRD_ST_zm(float z, float logM_MINm_ave, float *splined_value)
#else
void SFRD_ST_zm(float z, float *splined_value)
#endif
{
    float returned_value;

#ifdef INHOMO_FEEDBACK
    if (logM_MINm_ave<5)
      returned_value = gsl_spline2d_eval(SFRD_ST_z_splinem, z, 5, SFRD_ST_z_spline_accm, SFRD_ST_z_spline_accm_Mturn);
    else if (logM_MINm_ave>10)
      returned_value = gsl_spline2d_eval(SFRD_ST_z_splinem, z, 10, SFRD_ST_z_spline_accm, SFRD_ST_z_spline_accm_Mturn);
    else
      returned_value = gsl_spline2d_eval(SFRD_ST_z_splinem, z, logM_MINm_ave, SFRD_ST_z_spline_accm, SFRD_ST_z_spline_accm_Mturn);
#else
    returned_value = gsl_spline_eval(SFRD_ST_z_splinem, z, SFRD_ST_z_spline_accm);
#endif
    *splined_value = returned_value > 0. ? returned_value : 1e-40;
}
#endif

#ifdef INHOMO_FEEDBACK
void initialise_SFRD_Conditional_table(int Nsteps_zp, int Nfilter, float z[], double R[], float Mmin, float Alpha_star, float Fstar10)
#else
#ifdef REION_SM
void initialise_SFRD_Conditional_table(int Nsteps_zp, int Nfilter, float z[], double R[], float Mmin, float Alpha_star, float Fstar10, double REION_SM13_Z_RE, double REION_SM13_DELTA_Z_RE, double REION_SM13_DELTA_Z_SC)
#else
void initialise_SFRD_Conditional_table(int Nsteps_zp, int Nfilter, float z[], double R[], float Mmin, float MassTurnover, float Alpha_star, float Fstar10)
#endif
#endif
{
    double overdense_val;
    double overdense_large_high = Deltac, overdense_large_low = 1.5;
    double overdense_small_high = 1.5, overdense_small_low = -1. + 9e-8;
    double overdense_low_table[NSFR_low];
    float Mmax,Mlim_Fstar;
    int i,j,k,i_tot;
    float MassTurn;
#ifdef MINI_HALO
    double Mcrit_atom, Mcrit_RE;
#endif

    Mmax = RtoM(R[Nfilter-1]);
    Mlim_Fstar = Mass_limit_bisection(Mmin, Mmax, Alpha_star, Fstar10);
    initialiseSplinedSigmaM(Mmin,Mmax);
    fprintf(stderr, "In initialise_Fcollz_SFR_Conditional_table: Rmin = %6.4f, Rmax = %6.4f, Mmin = %.4e, Mmax = %.4e\n",
                                                                R[0],R[Nfilter-1],RtoM(R[0]),RtoM(R[Nfilter-1]));
#pragma omp parallel shared(overdense_small_low, overdense_small_high, log10_overdense_low_table, overdense_low_table) private(overdense_val, i)
{
#pragma omp for
    for (i=0; i<NSFR_low; i++) {
      overdense_val = log10(1. + overdense_small_low) + (double)i/((double)NSFR_low-1.)*(log10(1.+overdense_small_high)-log10(1.+overdense_small_low));
      log10_overdense_low_table[i] = overdense_val;
      overdense_low_table[i] = pow(10.,log10_overdense_low_table[i]) - 1.;
    }
}

#pragma omp parallel shared(Overdense_high_table, overdense_large_low, overdense_large_high) private(i)
{
#pragma omp for
    for (i=0; i<NSFR_high;i++) {
      Overdense_high_table[i] = overdense_large_low + (float)i/((float)NSFR_high-1.)*(overdense_large_high - overdense_large_low);
    }
}

    for (k=0; k < Nsteps_zp; k++) {
      i_tot = Nfilter*k;
      for (j=0; j < Nfilter; j++) {
        Mmax = RtoM(R[j]);
        initialiseGL_Nion(NGL_SFR, Mmin, Mmax);
#ifdef MINI_HALO
        Mcrit_atom = atomic_cooling_threshold(z[i_tot+j]);
#ifdef INHOMO_FEEDBACK
		MassTurn   = Mcrit_atom;
#else //INHOMO_FEEDBACK
#ifdef REION_SM
        Mcrit_RE   = reionization_feedback(z[i_tot+j], REION_SM13_Z_RE, REION_SM13_DELTA_Z_RE, REION_SM13_DELTA_Z_SC);
#else //REION_SM
        Mcrit_RE   = MassTurnover;
#endif //REION_SM
        MassTurn   = Mcrit_RE > Mcrit_atom ? Mcrit_RE : Mcrit_atom;
#endif //INHOMO_FEEDBACK
#else // MINI_HALO
        MassTurn   = MassTurnover;
#endif //MINI_HALO

#pragma omp parallel shared(log10_SFRD_z_low_table, i_tot, j, z, Mmax, overdense_low_table, MassTurn,Alpha_star, Fstar10, Mlim_Fstar) private(i)
{
#pragma omp for
        for (i=0; i<NSFR_low; i++){
            log10_SFRD_z_low_table[i_tot+j][i] = log10(GaussLegendreQuad_Nion(NGL_SFR,z[i_tot+j],log(Mmax),Deltac,overdense_low_table[i],MassTurn,Alpha_star,0.,Fstar10,1.,Mlim_Fstar,0.));
            if(log10_SFRD_z_low_table[i_tot+j][i] < -40.) log10_SFRD_z_low_table[i_tot+j][i] = -40.;
        }
}

#pragma omp parallel shared(SFRD_z_high_table, i_tot, j, z, Mmin, Mmax, Overdense_high_table, MassTurn,Alpha_star,Fstar10, Mlim_Fstar) private(i)
{
#pragma omp for
        for(i=0;i<NSFR_high;i++) {
            SFRD_z_high_table[i_tot+j][i] = Nion_ConditionalM(z[i_tot+j],log(Mmin),log(Mmax),Deltac,Overdense_high_table[i],MassTurn,Alpha_star,0.,Fstar10,1.,Mlim_Fstar,0.);
            if(SFRD_z_high_table[i_tot+j][i]<0.) SFRD_z_high_table[i_tot+j][i]=pow(10.,-40.0);
        }
}
      }
    }
}
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
void initialise_SFRD_Conditional_tablem(int Nsteps_zp, int Nfilter, float z[], double R[], float Mmin, double logMturn[], float Alpha_star, float Fstar10m)
#else
#ifdef REION_SM
void initialise_SFRD_Conditional_tablem(int Nsteps_zp, int Nfilter, float z[], double R[], float Mmin, float Alpha_star, float Fstar10m, double REION_SM13_Z_RE, double REION_SM13_DELTA_Z_RE, double REION_SM13_DELTA_Z_SC)
#else
void initialise_SFRD_Conditional_tablem(int Nsteps_zp, int Nfilter, float z[], double R[], float Mmin, float MassTurnover, float Alpha_star, float Fstar10m)
#endif
#endif
{
//    double overdense_val;
//    double overdense_large_high = Deltac, overdense_large_low = 1.5;
//    double overdense_small_high = 1.5, overdense_small_low = -1. + 9e-8;
    double overdense_low_table[NSFR_low];
    float Mmax,Mlim_Fstarm;
    int i,j,k,i_tot;
#ifdef INHOMO_FEEDBACK
    int q;
#endif
    double Mcrit_atom, Mcrit_LW, Mcrit_RE;
    float MassTurn;

    Mmax = RtoM(R[Nfilter-1]);
    Mlim_Fstarm = Mass_limit_bisection(Mmin, Mmax, Alpha_star, Fstar10m);
//    initialiseSplinedSigmaM(Mmin,Mmax);
//    fprintf(stderr, "In initialise_Fcollz_SFR_Conditional_table: Rmin = %6.4f, Rmax = %6.4f, Mmin = %.4e, Mmax = %.4e\n",
//                                                                R[0],R[Nfilter-1],RtoM(R[0]),RtoM(R[Nfilter-1]));
#pragma omp parallel shared(log10_overdense_low_table, overdense_low_table) private(i)
{
#pragma omp for
    for (i=0; i<NSFR_low; i++) {
      //overdense_val = log10(1. + overdense_small_low) + (double)i/((double)NSFR_low-1.)*(log10(1.+overdense_small_high)-log10(1.+overdense_small_low));
      //log10_overdense_low_table[i] = overdense_val;
      overdense_low_table[i] = pow(10.,log10_overdense_low_table[i]) - 1.;
    }
}
    //for (i=0; i<NSFR_high;i++) {
    //  Overdense_high_table[i] = overdense_large_low + (float)i/((float)NSFR_high-1.)*(overdense_large_high - overdense_large_low);
    //}
#ifdef INHOMO_FEEDBACK
#pragma omp parallel shared(logMturn, log10_overdense_low_table_Mturn, Overdense_high_table_Mturn) private(i)
{
#pragma omp for
    for (i=0; i<NMTURN;i++) {
      log10_overdense_low_table_Mturn[i] = logMturn[i];
      Overdense_high_table_Mturn[i] = (float)logMturn[i];
    }
}
#endif
    for (k=0; k < Nsteps_zp; k++) {
      i_tot = Nfilter*k;
      for (j=0; j < Nfilter; j++) {
        Mmax = RtoM(R[j]);
        initialiseGL_Nion(NGL_SFR, Mmin, Mmax);
        Mcrit_atom = atomic_cooling_threshold(z[i_tot+j]);
#ifdef INHOMO_FEEDBACK
        for (q=0; q<NMTURN; q++){
          MassTurn = pow(10, logMturn[q]);
#pragma omp parallel shared(log10_SFRD_z_low_tablem, i_tot, j, q, z, Mmax, overdense_low_table, Alpha_star,MassTurn,Mcrit_atom,Fstar10m,Mlim_Fstarm) private(i)
{
#pragma omp for
          for (i=0; i<NSFR_low; i++){
            log10_SFRD_z_low_tablem[i_tot+j][i+q*NSFR_low] = log10(GaussLegendreQuad_Nionm(NGL_SFR,z[i_tot+j],log(Mmax),Deltac,overdense_low_table[i],Alpha_star,MassTurn,Mcrit_atom,Fstar10m,Mlim_Fstarm));
            if(log10_SFRD_z_low_tablem[i_tot+j][i+q*NSFR_low] < -40.) log10_SFRD_z_low_tablem[i_tot+j][i+q*NSFR_low] = -40.;
          }
}

#pragma omp parallel shared(SFRD_z_high_tablem, i_tot, j, q, z, Mmin, Mmax, Overdense_high_table, Alpha_star,MassTurn,Mcrit_atom,Fstar10m,Mlim_Fstarm) private(i)
{
#pragma omp for
          for(i=0;i<NSFR_high;i++) {
            SFRD_z_high_tablem[i_tot+j][i+q*NSFR_high] = Nion_ConditionalMm(z[i_tot+j],log(Mmin),log(Mmax),Deltac,Overdense_high_table[i],Alpha_star,MassTurn,Mcrit_atom,Fstar10m,Mlim_Fstarm);
            if(SFRD_z_high_tablem[i_tot+j][i+q*NSFR_high]<0.) SFRD_z_high_tablem[i_tot+j][i+q*NSFR_high]=pow(10.,-40.0);
          }
}
        }
#else
        Mcrit_LW   = atomic_cooling_threshold(z[i_tot+j]);
#ifdef REION_SM
        Mcrit_RE   = reionization_feedback(z[i_tot+j], REION_SM13_Z_RE, REION_SM13_DELTA_Z_RE, REION_SM13_DELTA_Z_SC);
#else //REION_SM
        Mcrit_RE   = MassTurnover;
#endif //REION_SM
        MassTurn   = Mcrit_RE > Mcrit_LW ? Mcrit_RE : Mcrit_LW;
#pragma omp parallel shared(log10_SFRD_z_low_tablem, i_tot, j, z, Mmax, overdense_low_table, Alpha_star,MassTurn,Mcrit_atom,Fstar10m,Mlim_Fstarm) private(i)
{
#pragma omp for
        for (i=0; i<NSFR_low; i++){
            log10_SFRD_z_low_tablem[i_tot+j][i] = log10(GaussLegendreQuad_Nionm(NGL_SFR,z[i_tot+j],log(Mmax),Deltac,overdense_low_table[i],Alpha_star,MassTurn,Mcrit_atom,Fstar10m,Mlim_Fstarm));
            if(log10_SFRD_z_low_tablem[i_tot+j][i] < -40.) log10_SFRD_z_low_tablem[i_tot+j][i] = -40.;
        }
}

#pragma omp parallel shared(SFRD_z_high_tablem, i_tot, j, z, Mmin, Mmax, Overdense_high_table, Alpha_star,MassTurn,Mcrit_atom,Fstar10m,Mlim_Fstarm) private(i)
{
#pragma omp for
        for(i=0;i<NSFR_high;i++) {
            SFRD_z_high_tablem[i_tot+j][i] = Nion_ConditionalMm(z[i_tot+j],log(Mmin),log(Mmax),Deltac,Overdense_high_table[i],Alpha_star,MassTurn,Mcrit_atom,Fstar10m,Mlim_Fstarm);
            if(SFRD_z_high_tablem[i_tot+j][i]<0.) SFRD_z_high_tablem[i_tot+j][i]=pow(10.,-40.0);
        }
}
#endif
      }
    }
}
#endif

/* returns the mean star formation rate density at z in M_sun yr^-1 Mpc^-3 */
double mean_SFRD_dlnMhalo(double lnM, void *params){
  double z = *(double *)params;
  double M = exp(lnM);
  double f_ast = STELLAR_BARYON_FRAC * pow(M/1.0e10, STELLAR_BARYON_PL);
  double dndM = dNdM_st(z, M);
  if (f_ast > 1)
    f_ast = 1;

  return dndM * f_ast * exp(-M_TURNOVER/M) * M * M * OMb/OMm; //extra M for the dlnM
}

float mean_SFRD(double z){
  double result, error, lower_limit, upper_limit;
  gsl_function F;
  double timescale, rel_tol  = 0.001; //<- relative tolerance
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1000);

  F.function = &mean_SFRD_dlnMhalo;
  F.params = &z;
  lower_limit = log(M_TURNOVER/50.0);
  upper_limit = log(FMAX(1e16, M_TURNOVER*100));

  gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
               1000, GSL_INTEG_GAUSS61, w, &result, &error); 
  gsl_integration_workspace_free (w);

  timescale = t_STAR/hubble(z)/SperYR;

  return result / timescale;
}

#ifdef REION_SM
#define TOL_RE (double) 0.01

struct Param_REION_SM{
  double Z_RE;
  double DELTA_Z_RE;
  double DELTA_Z_SC;
};

void estimating_reionization(double ION_EFF_FACTOR, double ION_EFF_FACTOR_MINI, double ALPHA_STAR, double F_STAR10, double ALPHA_ESC, double F_ESC10, double F_STAR10m, double *REION_SM13_Z_RE, double *REION_SM13_DELTA_Z_RE, double *REION_SM13_DELTA_Z_SC){
  double zmax = 35.;
  double zmin = 4.5;
  double Q0, tt, zz; 
  //10Myr should be enough, the integration is slow, you don't want to use a too high cadence.
  double dt = 1e7*SperYR; 
  double Mcrit_RE, Mcrit_LW, Mcrit_atom;
  double M_MINa, M_MINm, M_MIN;
  double Mlim_Fstar, Mlim_Fesc, Mlim_Fstarm;
  double *ans, *redshifts, *redshifts_prev;
  double *Nion0, *Nion1;
  struct Param_REION_SM params = {
      .Z_RE = 9999.,
      .DELTA_Z_RE = 9999.,
      .DELTA_Z_SC = 9999.
  };
  double REION_SM13_Z_RE_updated = 9.3;
  double REION_SM13_DELTA_Z_RE_updated = 1.;
  int Nzz, i;
  FILE *F;
  gsl_interp_accel *zQ_spline_acc;
  gsl_spline *zQ_spline;

  // first calculate the number of redshift that we will be scrolling
  tt  = 0;
  Nzz = 0;
  zz  = zmax;
  while (zz > zmin){
      tt += dt;
      zz = ttoz(tt, zmax);
      Nzz += 1;
  }

  // allocate the memory
  ans            = (double *)malloc(sizeof(double) * Nzz);
  redshifts      = (double *)malloc(sizeof(double) * Nzz);
  redshifts_prev = (double *)malloc(sizeof(double) * Nzz);
  Nion0          = (double *)malloc(sizeof(double) * Nzz);
  Nion1          = (double *)malloc(sizeof(double) * Nzz);
  zQ_spline_acc  = gsl_interp_accel_alloc ();
  zQ_spline      = gsl_spline_alloc (gsl_interp_cspline, Nzz);

  // let's record the redshifts and previous redshifts required for Nion1
  tt = 0.;
  i  = 0 ;
  zz = zmax;
  while (zz > zmin){
      redshifts[i] = zz;
      redshifts_prev[i] = ttoz(dt,zz);
      tt += dt;
      zz  = ttoz(tt, zmax);
      i += 1;
  }

  fprintf(stderr, "You choose to use Sobacchi & Mesinger 2013 uniform photonheating background,\n\
I am trying to estimate the redshift and duration of reionization...\n");
  int Ntry = 1;
  // let's find the results!
  while ( (fabs(REION_SM13_Z_RE_updated - params.Z_RE) > TOL_RE) ||
          (fabs(REION_SM13_DELTA_Z_RE_updated - params.DELTA_Z_RE) > TOL_RE) ){
    params.Z_RE       = REION_SM13_Z_RE_updated;
    params.DELTA_Z_RE = REION_SM13_DELTA_Z_RE_updated;
    params.DELTA_Z_SC = soundcrossing_timescale_inredshift(params.Z_RE);
    fprintf(stderr, "Guessing %03d...z_re = %g, Delta z_re = %g, Delta z_sc = %g\n", Ntry, params.Z_RE, params.DELTA_Z_RE, params.DELTA_Z_SC);

//#pragma omp parallel shared(Nzz, redshifts, params,  ALPHA_STAR, F_STAR10, ALPHA_ESC, F_ESC10, F_STAR10m, ION_EFF_FACTOR, ION_EFF_FACTOR_MINI, Nion0, Nion1) private(i, Mcrit_RE, Mcrit_LW, Mcrit_atom, M_MINa, M_MINm, M_MIN, Mlim_Fstar, Mlim_Fesc, Mlim_Fstarm)
{
//#pragma omp for

    for (i=0; i<Nzz; i++){
      Mcrit_RE   = reionization_feedback(redshifts[i], params.Z_RE, params.DELTA_Z_RE, params.DELTA_Z_SC);
      Mcrit_LW   = lyman_werner_threshold(redshifts[i]);
      Mcrit_atom = atomic_cooling_threshold(redshifts[i]);
  
      M_MINa = Mcrit_RE > Mcrit_atom ? Mcrit_RE : Mcrit_atom;
      M_MINm = Mcrit_RE > Mcrit_LW   ? Mcrit_RE : Mcrit_LW;
      M_MIN  = M_MINa > M_MINm       ? M_MINm : M_MINa;
      M_MIN /= 50;
  
      Mlim_Fstar  = Mass_limit_bisection(M_MIN, 1e16, ALPHA_STAR, F_STAR10);
      Mlim_Fesc   = Mass_limit_bisection(M_MIN, 1e16, ALPHA_ESC, F_ESC10);
      Mlim_Fstarm = Mass_limit_bisection(M_MIN, 1e16, ALPHA_STAR, F_STAR10m);

      Nion0[i]  = ION_EFF_FACTOR      * Nion_ST(redshifts[i], M_MIN, M_MINa, ALPHA_STAR, ALPHA_ESC, F_STAR10, F_ESC10, Mlim_Fstar, Mlim_Fesc);
      Nion0[i] += ION_EFF_FACTOR_MINI * Nion_STm(redshifts[i], M_MIN, M_MINm, Mcrit_atom, ALPHA_STAR, F_STAR10m, Mlim_Fstarm);
  
      Nion1[i]  = ION_EFF_FACTOR      * Nion_ST(redshifts_prev[i], M_MIN, M_MINa, ALPHA_STAR, ALPHA_ESC, F_STAR10, F_ESC10, Mlim_Fstar, Mlim_Fesc);
      Nion1[i] += ION_EFF_FACTOR_MINI * Nion_STm(redshifts_prev[i], M_MIN, M_MINm, Mcrit_atom, ALPHA_STAR, F_STAR10m, Mlim_Fstarm);
    }
} 

    Q0 = 0;
    for (i=0; i<Nzz;i++){
      //Trec0 = 0.93 * 1e9 * SperYR * pow(C_HII/3.,-1) * pow(T_0/2e4,0.7) * pow((1.+redshifts[i])/7.,-3);
      //Q1    = Q0 + Nion1 - Nion0 - Q0 / Trec0 * dt;
      //assuming C_HII = 3, T_0 = 2e4
      Q0    += Nion1[i] - Nion0[i] - Q0 / 9.3e8 / SperYR * pow((1.+redshifts[i])/7., 3.) * dt;
      ans[i] = Q0;
    }
    
    // in case we cannot find the result
    if (ans[Nzz-1] < 0.62)
      ans[Nzz-1] = 1.0;

    // interpolate z vs ans
    gsl_spline_init(zQ_spline, ans, redshifts, Nzz);
    
    // find Z_RE and DELTA_Z_RE
    REION_SM13_Z_RE_updated       = gsl_spline_eval(zQ_spline, 0.5, zQ_spline_acc);
    REION_SM13_DELTA_Z_RE_updated = gsl_spline_eval(zQ_spline, 0.4, zQ_spline_acc) - gsl_spline_eval(zQ_spline, 0.6, zQ_spline_acc);
    Ntry += 1;
  }
  fprintf(stderr, "converged!\n");

  *REION_SM13_Z_RE       = params.Z_RE;
  *REION_SM13_DELTA_Z_RE = params.DELTA_Z_RE;
  *REION_SM13_DELTA_Z_SC = params.DELTA_Z_SC;
 
  // also save it so I don't need to do it again...
  if(F = fopen("../Parameter_files/REION_SM.H", "w")){
    if(fwrite(&params, sizeof(struct Param_REION_SM), 1, F) !=1)
      fprintf(stderr,  "reionization_helper_progs.c: Error writing reionization parameter file, ../Parameter_files/REION_SM.H\n");
    else
      fclose(F);
  }
  else
    fprintf(stderr,  "reionization_helper_progs.c: Error creating reionization parameter file, ../Parameter_files/REION_SM.H");

  free(ans);
  free(redshifts);
  free(redshifts_prev);
  free(Nion0);
  free(Nion1);
  gsl_spline_free (zQ_spline);
  gsl_interp_accel_free (zQ_spline_acc);
}

void reading_reionization_SM13parameters(double *REION_SM13_Z_RE, double *REION_SM13_DELTA_Z_RE, double *REION_SM13_DELTA_Z_SC, FILE *F){
  struct Param_REION_SM params;
  if(fread(&params, sizeof(struct Param_REION_SM), 1, F) !=1)
    fprintf(stderr,  "reionization_helper_progs.c: Error reading reionization parameter file, ../Parameter_files/REION_SM.H\nAborting...\n");
  else
    fclose(F);
  *REION_SM13_Z_RE       = params.Z_RE;
  *REION_SM13_DELTA_Z_RE = params.DELTA_Z_RE;
  *REION_SM13_DELTA_Z_SC = params.DELTA_Z_SC;
  fprintf(stderr, "Reading from ../Parameter_files/REION_SM.H...z_re = %g, Delta z_re = %g, Delta z_sc = %g\n", params.Z_RE, params.DELTA_Z_RE, params.DELTA_Z_SC);
}
#endif //REION_SM

#endif
