#ifdef REION_SM
#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"
#include "../Parameter_files/HEAT_PARAMS.H"
#include "../Cosmo_c_files/cosmo_progs.c"
#include "../Cosmo_c_files/ps.c"

#define TOL_RE (double) 0.01

struct Param_REION_SM{
    double REION_SM13_Z_RE;
    double REION_SM13_DELTA_Z_RE;
    double REION_SM13_DELTA_Z_SC;
} Param_REION_SM;

void estimating_reionization(){
  double zmax = 35.;
  double zmin = 0.;
  double Q0 = 0.;
  double tt = 0;
  double dt = 1e6*SperYR;
  double Mcrit_RE, Mcrit_LW, Mcrit_atom;
  double M_MINa, M_MINm, M_MIN;
  double Mlim_Fstar, Mlim_Fesc, Mlim_Fstarm;
  Param_REION_SM *params;
  params->REION_SM13_Z_RE = 9999.
  params->REION_SM13_DELTA_Z_RE = 9999.;
  params->REION_SM13_DELTA_Z_SC = 9999.; 
  double REION_SM13_Z_RE_updated = 9.3;
  double REION_SM13_DELTA_Z_RE_updated = 1.;
  double zz; 
  double *ans;
  double *redshifts, *redshifts_prev;
  int Nzz, i;
  gsl_interp_accel *zQ_spline_acc;
  gsl_spline *zQ_spline;

  double ION_EFF_FACTOR      = N_GAMMA_UV      * F_STAR10  * F_ESC10;
  double ION_EFF_FACTOR_MINI = N_GAMMA_UV_MINI * F_STAR10m * F_ESC10m;

  // first calculate the number of redshift that we will be scrolling
  Nzz = 0;
  zz  = zmax;
  while (zz > zmin){
      tt += dt;
      zz = ttoz(tt, zmax);
      Nzz += 1;
  }

  // allocate the memory
  ans            = (double *) malloc(sizeof(double) * Nzz);
  redshifts      = (double *) malloc(sizeof(double) * Nzz);
  redshifts_prev = (double *) malloc(sizeof(double) * Nzz);

  // let's record the redshifts and previous redshifts required for Nion1
  tt = 0.;
  i  = 0;
  zz = zmax;
  while (zz > zmin){
      redshifts[i] = zz;
      redshifts_prev[i] = ttoz(dt,zz);
      tt += dt;
      zz  = ttoz(tt, zmax);
      i += 1;
  }

  // let's find the results!
  while ( (fabs(REION_SM13_Z_RE_updated - params->REION_SM13_Z_RE) > TOL_RE) ||
          (fabs(REION_SM13_DELTA_Z_RE_updated - params->REION_SM13_DELTA_Z_RE) > TOL_RE) ){
    params->REION_SM13_Z_RE       = REION_SM13_Z_RE_updated;
    params->REION_SM13_DELTA_Z_RE = REION_SM13_DELTA_Z_RE_updated;
    params->REION_SM13_DELTA_Z_SC = soundcrossing_timescale(REION_SM13_Z_RE);

    zQ_spline_acc = gsl_interp_accel_alloc ();
    zQ_spline     = gsl_spline_alloc (gsl_interp_cspline, Nzz);

    Q0  = 0.;
    for (i=0; i<Nzz; i++){

      Mcrit_RE   = reionization_feedback(redshifts[i], params->REION_SM13_Z_RE, params->REION_SM13_DELTA_Z_RE, params->REION_SM13_DELTA_Z_SC);
      Mcrit_LW   = lyman_werner_threshold(redshifts[i]);
      Mcrit_atom = atomic_cooling_threshold(redshifts[i]);
  
      M_MINa = Mcrit_RE > Mcrit_atom ? Mcrit_RE : Mcrit_atom;
      M_MINm = Mcrit_RE > Mcrit_LW   ? Mcrit_RE : Mcrit_LW;
      M_MIN  = M_MINa > M_MINm       ? M_MINm : M_MINa;
      M_MIN /= 50;
  
      Mlim_Fstar  = Mass_limit_bisection(M_MIN, 1e16, ALPHA_STAR, F_STAR10);
      Mlim_Fesc   = Mass_limit_bisection(M_MIN, 1e16, ALPHA_ESC, F_ESC10);
      Mlim_Fstarm = Mass_limit_bisection(M_MIN, 1e16, ALPHA_STAR, F_STAR10m);
  
      Nion0  = ION_EFF_FACTOR      * Nion_ST(redshifts[i], M_MIN, M_MINa, ALPHA_STAR, ALPHA_ESC, F_STAR10, F_ESC10, Mlim_Fstar, Mlim_Fesc);
      Nion0 += ION_EFF_FACTOR_MINI * Nion_STm(redshifts[i], M_MIN, M_MINm, Mcrit_atom, ALPHA_STAR, F_STAR10m, Mlim_Fstarm);
  
      Nion1  = ION_EFF_FACTOR      * Nion_ST(redshifts_prev[i], M_MIN, M_MINa, ALPHA_STAR, ALPHA_ESC, F_STAR10, F_ESC10, Mlim_Fstar, Mlim_Fesc);
      Nion1 += ION_EFF_FACTOR_MINI * Nion_STm(redshifts_prev[i], M_MIN, M_MINm, Mcrit_atom, ALPHA_STAR, F_STAR10m, Mlim_Fstarm);
  
      //Trec0 = 0.93 * 1e9 * SperYR * pow(C_HII/3.,-1) * pow(T_0/2e4,0.7) * pow((1.+redshifts[i])/7.,-3);
      //Q1    = Q0 + Nion1 - Nion0 - Q0 / Trec0 * dt;
      //assuming C_HII = 3, T_0 = 2e4
      Q1     = Q0 + Nion1 - Nion0 - Q0 / 9.3e8 / SperYR * pow((1.+redshifts[i])/7., 3.) * dt;
      ans[i] = Q1 > 1 ? 1 : Q1;
      Q0     = ans[i];
    } 

    // interpolate z vs ans
    gsl_spline_init(zQ_spline, ans, redshifts, Nzz);
    
    // find Z_RE and DELTA_Z_RE
    REION_SM13_Z_RE_updated       = gsl_spline_eval(Q_spline, 1.0, Q_spline_acc);
    REION_SM13_DELTA_Z_RE_updated = gsl_spline_eval(Q_spline, 0.4, Q_spline_acc) - gsl_spline_eval(Q_spline, 0.6, Q_spline_acc);
  }
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);

  if(F = fopen("../Parameter_files/REION_SM.H", "w")){
    if(fwrite(params, sizeof(struct Param_REION_SM), 1, F) !=1)
      fprintf(stderr,  "reionization_helper_progs.c: Error writing reionization parameter file, ../Parameter_files/REION_SM.H");
  }
  else
    fprintf(stderr,  "reionization_helper_progs.c: Error creating reionization parameter file, ../Parameter_files/REION_SM.H");
}

void reading_reionization_SM13parameters(double *REION_SM13_Z_RE, double *REION_SM13_DELTA_Z_RE, double *REION_SM13_DELTA_Z_SC){
  Param_REION_SM *params;
  if(F = fopen("../Parameter_files/REION_SM.H", "r")){
    if(fread(params, sizeof(struct Param_REION_SM), 1, F) !=1)
      fprintf(stderr,  "reionization_helper_progs.c: Error reading reionization parameter file, ../Parameter_files/REION_SM.H");
  }
  else
      fprintf(stderr,  "reionization_helper_progs.c: Error openning reionization parameter file, ../Parameter_files/REION_SM.H");
  *REION_SM13_Z_RE       = params->REION_SM13_Z_RE;
  *REION_SM13_DELTA_Z_RE = params->REION_SM13_DELTA_Z_RE;
  *REION_SM13_DELTA_Z_SC = params->REION_SM13_DELTA_Z_SC;
}
#endif
