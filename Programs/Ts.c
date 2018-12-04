#include "heating_helper_progs.c"

/* 
  This is completed version.
*/

/*
  Program Ts calculates the spin temperature field, according to the perscription outlined in
  Mesinger, Furlanetto, Cen (2010).  The fluctuating component is sourced by the mean EPS collapsed fraction in spherical annyli surrounding each cell.

  Usage: Ts <REDSHIFT> [reload zp redshift]  [<stellar fraction for 10^10 Msun halos> <power law index for stellar fraction halo mass scaling> 
       <escape fraction for 10^10 Msun halos> <power law index for escape fraction halo mass scaling>
       <turn-over scale for the duty cycle of galaxies, in units of halo mass> <Soft band X-ray luminosity>]

  The last optional argument is the z' redshift output from which to reload intermediate
  evolution files in ../Boxes/Ts_evolution/

  Memory usage (in floats)~ (<NUMBER OF FILTER STEPS> + 3) x HII_DIM^3

  Author: Andrei Mesinger
  Date: 9.2.2009
*/


// New in v2
void init_21cmMC_arrays() { 

  int i;

  xi_SFR = calloc((NGL_SFR+1),sizeof(float));
  wi_SFR = calloc((NGL_SFR+1),sizeof(float));

  for (i=0; i < NUM_FILTER_STEPS_FOR_Ts; i++){
    SFRDLow_zpp_spline_acc[i] = gsl_interp_accel_alloc ();
    SFRDLow_zpp_spline[i] = gsl_spline_alloc (gsl_interp_cspline, NSFR_low);
    second_derivs_Nion_zpp[i] = calloc(NSFR_high,sizeof(float));

#ifdef MINI_HALO
    SFRDLow_zpp_spline_accm[i] = gsl_interp_accel_alloc ();
#ifdef INHOMO_FEEDBACK
    SFRDLow_zpp_spline_accm_Mturn[i] = gsl_interp_accel_alloc ();
    SFRDLow_zpp_splinem[i] = gsl_spline2d_alloc (gsl_interp2d_bicubic, NSFR_low, NMTURN);
    second_derivs_Nion_zppm[i][0] = (float *) malloc(NSFR_high*NMTURN*sizeof(float));
    second_derivs_Nion_zppm[i][1] = (float *) malloc(NSFR_high*NMTURN*sizeof(float));
    second_derivs_Nion_zppm[i][2] = (float *) malloc(NSFR_high*NMTURN*sizeof(float));
#else
    SFRDLow_zpp_splinem[i] = gsl_spline_alloc (gsl_interp_cspline, NSFR_low);
    second_derivs_Nion_zppm[i] = calloc(NSFR_high,sizeof(float));
#endif
#endif
  }

  redshift_interp_table = calloc(NUM_FILTER_STEPS_FOR_Ts*Nsteps_zp, sizeof(float)); // New

  log10_overdense_low_table = calloc(NSFR_low,sizeof(double));
#ifdef INHOMO_FEEDBACK
  log10_overdense_low_table_Mturn = calloc(NMTURN,sizeof(double));
#endif

  log10_SFRD_z_low_table = (double **)calloc(NUM_FILTER_STEPS_FOR_Ts*Nsteps_zp,sizeof(double *)); //New
#ifdef MINI_HALO
  log10_SFRD_z_low_tablem = (double **)calloc(NUM_FILTER_STEPS_FOR_Ts*Nsteps_zp,sizeof(double *)); //New
#endif
  for(i=0;i<NUM_FILTER_STEPS_FOR_Ts*Nsteps_zp;i++){  // New
#ifdef INHOMO_FEEDBACK
    log10_SFRD_z_low_table[i] = (double *)calloc(NSFR_low*NMTURN,sizeof(double));
    log10_SFRD_z_low_tablem[i] = (double *)calloc(NSFR_low*NMTURN,sizeof(double));
#else
    log10_SFRD_z_low_table[i] = (double *)calloc(NSFR_low,sizeof(double));
#ifdef MINI_HALO
    log10_SFRD_z_low_tablem[i] = (double *)calloc(NSFR_low,sizeof(double));
#endif
#endif
  }


  Overdense_high_table = calloc(NSFR_high,sizeof(float));
#ifdef INHOMO_FEEDBACK
  Overdense_high_table_Mturn = calloc(NMTURN,sizeof(float));
#endif

  SFRD_z_high_table = (float **)calloc(NUM_FILTER_STEPS_FOR_Ts*Nsteps_zp,sizeof(float *)); //New
#ifdef MINI_HALO
  SFRD_z_high_tablem = (float **)calloc(NUM_FILTER_STEPS_FOR_Ts*Nsteps_zp,sizeof(float *)); //New
#endif
  for(i=0;i<NUM_FILTER_STEPS_FOR_Ts*Nsteps_zp;i++){  // New
#ifdef INHOMO_FEEDBACK
    SFRD_z_high_table[i] = (float *)calloc(NSFR_high*NMTURN,sizeof(float));
    SFRD_z_high_tablem[i] = (float *)calloc(NSFR_high*NMTURN,sizeof(float));
#else
    SFRD_z_high_table[i] = (float *)calloc(NSFR_high,sizeof(float));
#ifdef MINI_HALO
    SFRD_z_high_tablem[i] = (float *)calloc(NSFR_high,sizeof(float));
#endif
#endif
  }
}

void destroy_21cmMC_arrays() {

  int i;

  free(xi_SFR);
  free(wi_SFR);
  for (i=0; i < NUM_FILTER_STEPS_FOR_Ts; i++){
    gsl_interp_accel_free (SFRDLow_zpp_spline_acc[i]);
    gsl_spline_free (SFRDLow_zpp_spline[i]);
    free(second_derivs_Nion_zpp[i]);

#ifdef MINI_HALO
    gsl_interp_accel_free (SFRDLow_zpp_spline_accm[i]);
#ifdef INHOMO_FEEDBACK
    gsl_interp_accel_free (SFRDLow_zpp_spline_accm_Mturn[i]);
    gsl_spline2d_free (SFRDLow_zpp_splinem[i]);
    free(second_derivs_Nion_zppm[i][0]);
    free(second_derivs_Nion_zppm[i][1]);
    free(second_derivs_Nion_zppm[i][2]);
#else
    gsl_spline_free (SFRDLow_zpp_splinem[i]);
    free(second_derivs_Nion_zppm[i]);
#endif
#endif
  }

  free(redshift_interp_table);

  free(log10_overdense_low_table);
#ifdef INHOMO_FEEDBACK
  free(log10_overdense_low_table_Mturn);
#endif

  for(i=0;i<NUM_FILTER_STEPS_FOR_Ts*Nsteps_zp;i++) {
    free(log10_SFRD_z_low_table[i]);
#ifdef MINI_HALO
    free(log10_SFRD_z_low_tablem[i]);
#endif
  }
  free(log10_SFRD_z_low_table);
#ifdef MINI_HALO
  free(log10_SFRD_z_low_tablem);
#endif

  free(Overdense_high_table);
#ifdef INHOMO_FEEDBACK
  free(Overdense_high_table_Mturn);
#endif

  for(i=0;i<NUM_FILTER_STEPS_FOR_Ts*Nsteps_zp;i++) {
      free(SFRD_z_high_table[i]);
#ifdef MINI_HALO
      free(SFRD_z_high_tablem[i]);
#endif
  }
  free(SFRD_z_high_table);
#ifdef MINI_HALO
  free(SFRD_z_high_tablem);
#endif

  gsl_spline_free (SFRD_ST_z_spline);
  gsl_interp_accel_free (SFRD_ST_z_spline_acc);
  gsl_spline_free (Nion_z_spline);
  gsl_interp_accel_free (Nion_z_spline_acc);
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
  gsl_spline2d_free (SFRD_ST_z_splinem);
  gsl_interp_accel_free (SFRD_ST_z_spline_accm_Mturn);
  gsl_spline2d_free (Nion_z_splinem);
  gsl_interp_accel_free (Nion_z_spline_accm_Mturn);
#else
  gsl_spline_free (SFRD_ST_z_splinem);
#endif
  gsl_interp_accel_free (SFRD_ST_z_spline_accm);
  gsl_interp_accel_free (Nion_z_spline_accm);
#endif
}

/* Maximum allowed value for the kinetic temperature.
   Useful to set to avoid some spurious behaviour 
   when the code is run with redshift poor resolution 
   and very high X-ray heating efficiency */
#define MAX_TK (float) 5e4 
#define box_ct_increment (int) (HII_TOT_NUM_PIXELS/1e5+1)


int main(int argc, char ** argv){
  fftwf_complex *box, *unfiltered_box;
  fftwf_plan plan;
  unsigned long long ct, sample_ct;
  int R_ct,i,j,k, COMPUTE_Ts, x_e_ct;
  float REDSHIFT, growth_factor_z, R, R_factor, zp, mu_for_Ts, filling_factor_of_HI_zp;
  int ithread;
  float *Tk_box, *x_e_box, *Ts, J_star_Lya, dzp, prev_zp, zpp, prev_zpp, prev_R;
  FILE *F, *GLOBAL_EVOL, *OUT;
  char filename[500];
  float dz, zeta_ion_eff, Tk_BC, xe_BC, nu, zprev, zcurr, curr_delNL0[NUM_FILTER_STEPS_FOR_Ts];
  double *evolve_ans, ans[2], Tk_ave, J_alpha_ave, xalpha_ave, J_alpha_tot, Xheat_ave,Xion_ave;
#ifdef INHOMO_FEEDBACK
  double dansdz[6], J_LW_ave, nu_nplus1;
#else
  double dansdz[5];
#endif
  double freq_int_heat_tbl[x_int_NXHII][NUM_FILTER_STEPS_FOR_Ts], freq_int_ion_tbl[x_int_NXHII][NUM_FILTER_STEPS_FOR_Ts], freq_int_lya_tbl[x_int_NXHII][NUM_FILTER_STEPS_FOR_Ts];
#ifdef MINI_HALO
  double freq_int_heat_tblm[x_int_NXHII][NUM_FILTER_STEPS_FOR_Ts], freq_int_ion_tblm[x_int_NXHII][NUM_FILTER_STEPS_FOR_Ts], freq_int_lya_tblm[x_int_NXHII][NUM_FILTER_STEPS_FOR_Ts];
#endif
  int goodSteps,badSteps;
  int m_xHII_low, m_xHII_high, n_ct, zp_ct;
  double freq_int_heat[NUM_FILTER_STEPS_FOR_Ts], freq_int_ion[NUM_FILTER_STEPS_FOR_Ts], freq_int_lya[NUM_FILTER_STEPS_FOR_Ts];
#ifdef MINI_HALO
  double freq_int_heatm[NUM_FILTER_STEPS_FOR_Ts], freq_int_ionm[NUM_FILTER_STEPS_FOR_Ts], freq_int_lyam[NUM_FILTER_STEPS_FOR_Ts];
#endif
  double nuprime, fcoll_R, Ts_ave;
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
  float *J_21_LW=NULL, *log10_Mcrit_LW, log10_Mcrit_mol;
  fftwf_complex *log10_Mcrit_LW_unfiltered=NULL, *log10_Mcrit_LW_filtered=NULL;
#endif
  double fcollm_R;
#ifdef REION_SM
  double REION_SM13_Z_RE, REION_SM13_DELTA_Z_RE, REION_SM13_DELTA_Z_SC;
#endif
#endif
  float *delNL0[NUM_FILTER_STEPS_FOR_Ts], delNL_zpp, xHII_call, curr_xalpha;
  float z, Jalpha, TK, TS, xe, deltax;
  time_t start_time, curr_time;
  double J_alpha_threads[NUMCORES], xalpha_threads[NUMCORES], Xheat_threads[NUMCORES],
   Xion_threads[NUMCORES], lower_int_limit;
#ifdef INHOMO_FEEDBACK
  double J_LW_threads[NUMCORES];
#endif
  float Splined_Nion_ST_zp, Splined_SFRD_ST_zpp,ION_EFF_FACTOR,fcoll; // New in v2
#ifdef MINI_HALO
  float Splined_Nion_ST_zpm, Splined_SFRD_ST_zppm,ION_EFF_FACTOR_MINI,fcollm; // New in v2.1
#endif
  float zp_table; //New in v2
  int counter,arr_num; // New in v2
  double Luminosity_conversion_factor;
#ifdef MINI_HALO
  double Luminosity_conversion_factorm;
#endif
  float prev_zp_temp, zp_temp;
  int RESTART = 0;

 
  /**********  BEGIN INITIALIZATION   **************************************/
  //New in v2
#ifdef SHARP_CUTOFF
    if (argc == 3){
      RESTART = 1;
      zp = atof(argv[2]);
    }
    else if (argc != 2){
      fprintf(stderr, "Usage: Ts <REDSHIFT>  [reload zp redshift]\nAborting...\n");
      return -1;
    }
    X_LUMINOSITY = pow(10.,L_X);
    F_STAR10 = STELLAR_BARYON_FRAC;
    M_MIN = M_TURNOVER;
#else //SHARP_CUTOFF
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
    M_MIN = 1e5;
    if (argc  == 12) {
      RESTART = 1;
      zp = atof(argv[2]);
      F_STAR10 = atof(argv[3]);
      ALPHA_STAR = atof(argv[4]);
      F_ESC10 = atof(argv[5]);
      ALPHA_ESC = atof(argv[6]);
      T_AST = atof(argv[7]);
      X_LUMINOSITY = pow(10.,atof(argv[8]));
      F_STAR10m = atof(argv[9]);
      F_ESC10m = atof(argv[10]);
      X_LUMINOSITYm = pow(10.,atof(argv[11]));
    }
    else if (argc == 11) {
      F_STAR10 = atof(argv[2]);
      ALPHA_STAR = atof(argv[3]);
      F_ESC10 = atof(argv[4]);
      ALPHA_ESC = atof(argv[5]);
      T_AST = atof(argv[6]);
      X_LUMINOSITY = pow(10.,atof(argv[7]));
      F_STAR10m = atof(argv[8]);
      F_ESC10m = atof(argv[9]);
      X_LUMINOSITYm = pow(10.,atof(argv[10]));
    }
    else if (argc == 3) {
      RESTART = 1;
      zp = atof(argv[2]);
      F_STAR10 = STELLAR_BARYON_FRAC;
      ALPHA_STAR = STELLAR_BARYON_PL;
      F_ESC10 = ESC_FRAC;
      ALPHA_ESC = ESC_PL;
      T_AST = t_STAR;
      X_LUMINOSITY = pow(10.,L_X);
      F_STAR10m = STELLAR_BARYON_FRAC_MINI;
      F_ESC10m = ESC_FRAC_MINI;
      X_LUMINOSITYm = pow(10.,L_X_MINI);
    }
    else if (argc == 2) {
      F_STAR10 = STELLAR_BARYON_FRAC;
      ALPHA_STAR = STELLAR_BARYON_PL;
      F_ESC10 = ESC_FRAC;
      ALPHA_ESC = ESC_PL;
      T_AST = t_STAR;
      X_LUMINOSITY = pow(10.,L_X);
      F_STAR10m = STELLAR_BARYON_FRAC_MINI;
      F_ESC10m = ESC_FRAC_MINI;
      X_LUMINOSITYm = pow(10.,L_X_MINI);
    }
    else {
      fprintf(stderr, "Usage: Ts <REDSHIFT> [reload zp redshift] [<f_star10> <alpha_star> <f_esc10> <alpha_esc> <M_turn> <t_star> <X_luminosity>] \nAborting...\n");
      return -1;
    }
#else //INHOMO_FEEDBACK
    M_MIN = 1e16; //calculated later
#ifdef REION_SM
    if (argc  == 12) {
      RESTART = 1;
      zp = atof(argv[2]);
      F_STAR10 = atof(argv[3]);
      ALPHA_STAR = atof(argv[4]);
      F_ESC10 = atof(argv[5]);
      ALPHA_ESC = atof(argv[6]);
      T_AST = atof(argv[7]);
      X_LUMINOSITY = pow(10.,atof(argv[8]));
      F_STAR10m = atof(argv[9]);
      F_ESC10m = atof(argv[10]);
      X_LUMINOSITYm = pow(10.,atof(argv[11]));
    }
    else if (argc == 11) {
      F_STAR10 = atof(argv[2]);
      ALPHA_STAR = atof(argv[3]);
      F_ESC10 = atof(argv[4]);
      ALPHA_ESC = atof(argv[5]);
      T_AST = atof(argv[6]);
      X_LUMINOSITY = pow(10.,atof(argv[7]));
      F_STAR10m = atof(argv[8]);
      F_ESC10m = atof(argv[9]);
      X_LUMINOSITYm = pow(10.,atof(argv[10]));
    }
    else if (argc == 3) {
      RESTART = 1;
      zp = atof(argv[2]);
      F_STAR10 = STELLAR_BARYON_FRAC;
      ALPHA_STAR = STELLAR_BARYON_PL;
      F_ESC10 = ESC_FRAC;
      ALPHA_ESC = ESC_PL;
      T_AST = t_STAR;
      X_LUMINOSITY = pow(10.,L_X);
      F_STAR10m = STELLAR_BARYON_FRAC_MINI;
      F_ESC10m = ESC_FRAC_MINI;
      X_LUMINOSITYm = pow(10.,L_X_MINI);
    }
    else if (argc == 2) {
      F_STAR10 = STELLAR_BARYON_FRAC;
      ALPHA_STAR = STELLAR_BARYON_PL;
      F_ESC10 = ESC_FRAC;
      ALPHA_ESC = ESC_PL;
      T_AST = t_STAR;
      X_LUMINOSITY = pow(10.,L_X);
      F_STAR10m = STELLAR_BARYON_FRAC_MINI;
      F_ESC10m = ESC_FRAC_MINI;
      X_LUMINOSITYm = pow(10.,L_X_MINI);
    }
    else {
      fprintf(stderr, "Usage: Ts <REDSHIFT> [reload zp redshift] [<f_star10> <alpha_star> <f_esc10> <alpha_esc> <M_turn> <t_star> <X_luminosity>] \nAborting...\n");
      return -1;
    }
#else //REION_SM
    if (argc  == 13) {
      RESTART = 1;
      zp = atof(argv[2]);
      F_STAR10 = atof(argv[3]);
      ALPHA_STAR = atof(argv[4]);
      F_ESC10 = atof(argv[5]);
      ALPHA_ESC = atof(argv[6]);
      M_TURN = atof(argv[7]); 
      T_AST = atof(argv[8]);
      X_LUMINOSITY = pow(10.,atof(argv[9]));
      F_STAR10m = atof(argv[10]);
      F_ESC10m = atof(argv[11]);
      X_LUMINOSITYm = pow(10.,atof(argv[12]));
    }
    else if (argc == 12) {
      F_STAR10 = atof(argv[2]);
      ALPHA_STAR = atof(argv[3]);
      F_ESC10 = atof(argv[4]);
      ALPHA_ESC = atof(argv[5]);
      M_TURN = atof(argv[6]);
      T_AST = atof(argv[7]);
      X_LUMINOSITY = pow(10.,atof(argv[8]));
      F_STAR10m = atof(argv[9]);
      F_ESC10m = atof(argv[10]);
      X_LUMINOSITYm = pow(10.,atof(argv[11]));
    }
    else if (argc == 3) {
      RESTART = 1;
      zp = atof(argv[2]);
      F_STAR10 = STELLAR_BARYON_FRAC;
      ALPHA_STAR = STELLAR_BARYON_PL;
      F_ESC10 = ESC_FRAC;
      ALPHA_ESC = ESC_PL;
      M_TURN = M_TURNOVER; 
      T_AST = t_STAR;
      X_LUMINOSITY = pow(10.,L_X);
      F_STAR10m = STELLAR_BARYON_FRAC_MINI;
      F_ESC10m = ESC_FRAC_MINI;
      X_LUMINOSITYm = pow(10.,L_X_MINI);
    }
    else if (argc == 2) {
      F_STAR10 = STELLAR_BARYON_FRAC;
      ALPHA_STAR = STELLAR_BARYON_PL;
      F_ESC10 = ESC_FRAC;
      ALPHA_ESC = ESC_PL;
      M_TURN = M_TURNOVER; 
      T_AST = t_STAR;
      X_LUMINOSITY = pow(10.,L_X);
      F_STAR10m = STELLAR_BARYON_FRAC_MINI;
      F_ESC10m = ESC_FRAC_MINI;
      X_LUMINOSITYm = pow(10.,L_X_MINI);
    }
    else {
      fprintf(stderr, "Usage: Ts <REDSHIFT> [reload zp redshift] [<f_star10> <alpha_star> <f_esc10> <alpha_esc> <M_turn> <t_star> <X_luminosity>] \nAborting...\n");
      return -1;
    }
#endif //REION_SM
#endif//INHOMO_FEEDBACK
    ION_EFF_FACTOR      = N_GAMMA_UV      * F_STAR10  * F_ESC10;
    ION_EFF_FACTOR_MINI = N_GAMMA_UV_MINI * F_STAR10m * F_ESC10m;
#else //MINI_HALO
    if (argc  == 10) {
      RESTART = 1;
      zp = atof(argv[2]);
      F_STAR10 = atof(argv[3]);
      ALPHA_STAR = atof(argv[4]);
      F_ESC10 = atof(argv[5]);
      ALPHA_ESC = atof(argv[6]);
      M_TURN = atof(argv[7]); 
      T_AST = atof(argv[8]);
      X_LUMINOSITY = pow(10.,atof(argv[9]));
    }
    else if (argc == 9) {
      F_STAR10 = atof(argv[2]);
      ALPHA_STAR = atof(argv[3]);
      F_ESC10 = atof(argv[4]);
      ALPHA_ESC = atof(argv[5]);
      M_TURN = atof(argv[6]);
      T_AST = atof(argv[7]);
      X_LUMINOSITY = pow(10.,atof(argv[8]));
    }
    else if (argc == 3) {
      RESTART = 1;
      zp = atof(argv[2]);
      F_STAR10 = STELLAR_BARYON_FRAC;
      ALPHA_STAR = STELLAR_BARYON_PL;
      F_ESC10 = ESC_FRAC;
      ALPHA_ESC = ESC_PL;
      M_TURN = M_TURNOVER; 
      T_AST = t_STAR;
      X_LUMINOSITY = pow(10.,L_X);
    }
    else if (argc == 2) {
      F_STAR10 = STELLAR_BARYON_FRAC;
      ALPHA_STAR = STELLAR_BARYON_PL;
      F_ESC10 = ESC_FRAC;
      ALPHA_ESC = ESC_PL;
      M_TURN = M_TURNOVER; 
      T_AST = t_STAR;
      X_LUMINOSITY = pow(10.,L_X);
    }
    else {
      fprintf(stderr, "Usage: Ts <REDSHIFT> [reload zp redshift] [<f_star10> <alpha_star> <f_esc10> <alpha_esc> <M_turn> <t_star> <X_luminosity>] \nAborting...\n");
      return -1;
    }
    ION_EFF_FACTOR = N_GAMMA_UV * F_STAR10 * F_ESC10;
    M_MIN = M_TURNOVER;
#endif //MINI_HALO
#endif //SHARP_CUTOFF
  REDSHIFT = atof(argv[1]);
  system("mkdir -p ../Log_files");
  system("mkdir -p ../Output_files");
  system("mkdir -p ../Boxes/Ts_evolution/");
  system("mkdir -p ../Output_files/Ts_outs/");
  system("cp ../Parameter_files/* ../Output_files/Ts_outs/");
  system("cp ../Parameter_files/* ../Boxes/Ts_evolution/");
  init_ps();
  omp_set_num_threads(NUMCORES);
  growth_factor_z = dicke(REDSHIFT);
   
  // open log file
  if (!(LOG = fopen("../Log_files/Ts_log", "w") ) ){
    fprintf(stderr, "Unable to open log file for writting\nAborting...\n");
    return -1;
  }
 
  // Initialize some interpolation tables
  if (init_heat() < 0){
    fclose(LOG); fclose(GLOBAL_EVOL);
    return -1;
  }
 
  // check if we are in the really high z regime before the first stars; if so, simple
  if (REDSHIFT > Z_HEAT_MAX){
    xe = xion_RECFAST(REDSHIFT,0);
    TK = T_RECFAST(REDSHIFT,0);
    
    // open input
    sprintf(filename, "../Boxes/updated_smoothed_deltax_z%06.2f_%i_%.0fMpc", 
        REDSHIFT, HII_DIM, BOX_LEN);
    F = fopen(filename, "rb");
    if ( !(F = fopen(filename, "rb") ) ){
      fprintf(stderr, "Error opening file %s for reading.\nAborting...\n", filename);
      fprintf(LOG, "Error opening file %s for reading.\nAborting...\n", filename);
      destruct_heat(); return -1;
    }
    fprintf(stderr, "Opened density file %s for reading\n", filename);
    fprintf(LOG, "Opened density file %s for reading\n", filename);
 
    // open output
    // New in v2
#ifdef SHARP_CUTOFF
    sprintf(filename, "../Boxes/Ts_z%06.2f_L_X%.1e_alphaX%.1f_MminX%.1e_zetaIon%.2f_Pop%i_%i_%.0fMpc", REDSHIFT, X_LUMINOSITY, X_RAY_SPEC_INDEX, M_MIN, HII_EFF_FACTOR, Pop, HII_DIM, BOX_LEN); 
#else //SHARP_CUTOFF
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
    sprintf(filename, "../Boxes/Ts_z%06.2f_L_X%.1e_alphaX%.1f_f_star%06.4f_alpha_star%06.4f_t_star%06.4f_f_star10m%06.4f_f_esc10m%06.4f_L_Xm%.1e_alphaXm%.1f_%i_%.0fMpc", REDSHIFT, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, T_AST, F_STAR10m, F_ESC10m, X_LUMINOSITYm, X_RAY_SPEC_INDEX_MINI, HII_DIM, BOX_LEN); 
#else //INHOMO_FEEDBACK
#ifdef REION_SM
    sprintf(filename, "../Boxes/Ts_z%06.2f_L_X%.1e_alphaX%.1f_f_star%06.4f_alpha_star%06.4f_Msm_t_star%06.4f_f_star10m%06.4f_f_esc10m%06.4f_L_Xm%.1e_alphaXm%.1f_%i_%.0fMpc", REDSHIFT, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, T_AST, F_STAR10m, F_ESC10m, X_LUMINOSITYm, X_RAY_SPEC_INDEX_MINI, HII_DIM, BOX_LEN); 
#else //REION_SM
    sprintf(filename, "../Boxes/Ts_z%06.2f_L_X%.1e_alphaX%.1f_f_star%06.4f_alpha_star%06.4f_Mturn%.1e_t_star%06.4f_f_star10m%06.4f_f_esc10m%06.4f_L_Xm%.1e_alphaXm%.1f_%i_%.0fMpc", REDSHIFT, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, M_TURN, T_AST, F_STAR10m, F_ESC10m, X_LUMINOSITYm, X_RAY_SPEC_INDEX_MINI, HII_DIM, BOX_LEN); 
#endif //REION_SM
#endif //INHOMO_FEEDBACK
#else //MINI_HALO
    sprintf(filename, "../Boxes/Ts_z%06.2f_L_X%.1e_alphaX%.1f_f_star%06.4f_alpha_star%06.4f_Mturn%.1e_t_star%06.4f_Pop%i_%i_%.0fMpc", REDSHIFT, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, M_TURN, T_AST, Pop, HII_DIM, BOX_LEN); 
#endif //MINI_HALO
#endif //SHARP_CUTOFF
    if (!(OUT=fopen(filename, "wb"))){
      fprintf(stderr, "Ts.c: WARNING: Unable to open output file %s\n", filename);
      fprintf(LOG, "Ts.c: WARNING: Unable to open output file %s\n", filename);
      destruct_heat(); return -1;
    }
    fprintf(stderr, "Opened TS file %s for writting\n", filename);
    fprintf(LOG, "Opened TS file %s for writting\n", filename);
 
    // read file
    for (i=0; i<HII_DIM; i++){
      for (j=0; j<HII_DIM; j++){
        for (k=0; k<HII_DIM; k++){
      if (fread(&deltax, sizeof(float), 1, F)!=1){
       fprintf(stderr, "Error reading-in binary density file\nAborting...\n");
       fprintf(LOG, "Error reading-in binary density file\nAborting...\n");
       destruct_heat(); return -1;
      }
 
      // compute the spin temperature
     TS = get_Ts(REDSHIFT, deltax, TK, xe, 0, &curr_xalpha);
 
     // and print it out
     if (fwrite(&TS, sizeof(float), 1, OUT)!=1){
       fprintf(stderr, "Ts.c: Write error occured while writting Tk box.\n");
       fprintf(LOG, "Ts.c: Write error occured while writting Tk box.\n");
       destruct_heat(); return -1;
      }
 
        }
      }
    }
 
    destruct_heat(); fclose(F); fclose(OUT);
    return 0;
  }
 
 
  // open global evolution output file
  // New in v2
#ifdef SHARP_CUTOFF
  sprintf(filename, "../Output_files/Ts_outs/global_evolution_zetaIon%.2f_Nsteps%i_zprimestepfactor%.3f_L_X%.1e_alphaX%.1f_TvirminX%.1e_Pop%i_%i_%.0fMpc", HII_EFF_FACTOR, NUM_FILTER_STEPS_FOR_Ts, ZPRIME_STEP_FACTOR, X_LUMINOSITY, X_RAY_SPEC_INDEX, M_TURN, Pop, HII_DIM, BOX_LEN);
    if (argc > 2) // restarting
      GLOBAL_EVOL = fopen(filename, "a");
    else
      GLOBAL_EVOL = fopen(filename, "w");
#else //SHARP_CUTOFF
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
  sprintf(filename, "../Output_files/Ts_outs/global_evolution_Nsteps%i_zprimestepfactor%.3f_L_X%.1e_alphaX%.1f_f_star10%06.4f_alpha_star%06.4f_f_esc10%06.4f_alpha_esc%06.4f_t_star%06.4f_f_star10m%06.4f_f_esc10m%06.4f_L_Xm%.1e_alphaXm%.1f_%i_%.0fMpc", NUM_FILTER_STEPS_FOR_Ts, ZPRIME_STEP_FACTOR, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, T_AST, F_STAR10m, F_ESC10m, X_LUMINOSITYm, X_RAY_SPEC_INDEX_MINI, HII_DIM, BOX_LEN);
#else //INHOMO_FEEDBACK
#ifdef REION_SM
  sprintf(filename, "../Output_files/Ts_outs/global_evolution_Nsteps%i_zprimestepfactor%.3f_L_X%.1e_alphaX%.1f_f_star10%06.4f_alpha_star%06.4f_f_esc10%06.4f_alpha_esc%06.4f_Msm_t_star%06.4f_f_star10m%06.4f_f_esc10m%06.4f_L_Xm%.1e_alphaXm%.1f_%i_%.0fMpc", NUM_FILTER_STEPS_FOR_Ts, ZPRIME_STEP_FACTOR, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, T_AST, F_STAR10m, F_ESC10m, X_LUMINOSITYm, X_RAY_SPEC_INDEX_MINI, HII_DIM, BOX_LEN);
#else //REION_SM
  sprintf(filename, "../Output_files/Ts_outs/global_evolution_Nsteps%i_zprimestepfactor%.3f_L_X%.1e_alphaX%.1f_f_star10%06.4f_alpha_star%06.4f_f_esc10%06.4f_alpha_esc%06.4f_Mturn%.1e_t_star%06.4f_f_star10m%06.4f_f_esc10m%06.4f_L_Xm%.1e_alphaXm%.1f_%i_%.0fMpc", NUM_FILTER_STEPS_FOR_Ts, ZPRIME_STEP_FACTOR, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, T_AST, F_STAR10m, F_ESC10m, X_LUMINOSITYm, X_RAY_SPEC_INDEX_MINI, HII_DIM, BOX_LEN);
#endif //REION_SM
#endif //INHOMO_FEEDBACK
#else //MINI_HALO
  sprintf(filename, "../Output_files/Ts_outs/global_evolution_Nsteps%i_zprimestepfactor%.3f_L_X%.1e_alphaX%.1f_f_star10%06.4f_alpha_star%06.4f_f_esc10%06.4f_alpha_esc%06.4f_Mturn%.1e_t_star%06.4f_Pop%i_%i_%.0fMpc", NUM_FILTER_STEPS_FOR_Ts, ZPRIME_STEP_FACTOR, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, T_AST, Pop, HII_DIM, BOX_LEN);
#endif //MINI_HALO
    if (argc == 3 || argc == 9) // restarting
      GLOBAL_EVOL = fopen(filename, "a");
    else
      GLOBAL_EVOL = fopen(filename, "w");
#endif //SHARP_CUTOFF
  if (!GLOBAL_EVOL){
    fprintf(stderr, "Unable to open global evolution file at %s\nAborting...\n",
        filename);
    fprintf(LOG, "Unable to open global evolution file at %s\nAborting...\n",
        filename);
    fclose(LOG);
    return -1;
  }
 
  // set boundary conditions for the evolution equations->  values of Tk and x_e at Z_HEAT_MAX
  if (XION_at_Z_HEAT_MAX > 0) // user has opted to use his/her own value
    xe_BC = XION_at_Z_HEAT_MAX;
  else// will use the results obtained from recfast
    xe_BC = xion_RECFAST(Z_HEAT_MAX,0);
  if (TK_at_Z_HEAT_MAX > 0)
    Tk_BC = TK_at_Z_HEAT_MAX;
  else
    Tk_BC = T_RECFAST(Z_HEAT_MAX,0);


  /******  Now allocate large arrays  ******/

  // allocate memory for the nonlinear density field and open file
  sprintf(filename, "../Boxes/updated_smoothed_deltax_z%06.2f_%i_%.0fMpc", 
      REDSHIFT, HII_DIM, BOX_LEN);
  F = fopen(filename, "rb");
  if ( !(F = fopen(filename, "rb") ) ){
    fprintf(stderr, "Error opening file %s for reading.\nAborting...\n", filename);
    fprintf(LOG, "Error opening file %s for reading.\nAborting...\n", filename);
    fclose(LOG); fclose(GLOBAL_EVOL);
    destruct_heat();
    return -1;
  }
  if (!(box = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS))){
    fprintf(stderr, "Error in memory allocation for %s\nAborting...\n", filename);
    fprintf(LOG, "Error in memory allocation for %s\nAborting...\n", filename);
    fclose(LOG); fclose(GLOBAL_EVOL); fclose(F);   destruct_heat();
    return -1;
  }
  if (!(unfiltered_box = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS))){
    fprintf(stderr, "Error in memory allocation for %s\nAborting...\n", filename);
    fprintf(LOG, "Error in memory allocation for %s\nAborting...\n", filename);
    fclose(LOG); fclose(GLOBAL_EVOL);fclose(F);   destruct_heat(); fftwf_free(box);
    return -1;
  }
  fprintf(stderr, "Reading in deltax box\n");
  fprintf(LOG, "Reading in deltax box\n");
  for (i=0; i<HII_DIM; i++){
    for (j=0; j<HII_DIM; j++){
      for (k=0; k<HII_DIM; k++){
    if (fread((float *)unfiltered_box + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, F)!=1){
      fprintf(stderr, "Error reading-in binary file %s\nAborting...\n", filename);
      fprintf(LOG, "Error reading-in binary file %s\nAborting...\n", filename);
      fftwf_free(box); fclose(GLOBAL_EVOL); fclose(F); fclose(LOG); fftwf_free(unfiltered_box);
      destruct_heat();
      return -1;
    }
      }
    }
  }
  fclose(F);

#ifdef INHOMO_FEEDBACK
  J_21_LW  = (float *) fftwf_malloc(sizeof(float)*HII_TOT_NUM_PIXELS);
  log10_Mcrit_LW = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS);
  log10_Mcrit_LW_unfiltered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
  log10_Mcrit_LW_filtered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
  if (!log10_Mcrit_LW_unfiltered || !log10_Mcrit_LW_filtered || !J_21_LW || !log10_Mcrit_LW){
    fprintf(stderr, "Ts.c: Error allocating memory for feedback boxes\nAborting...\n");
    return -1;
  }
#endif

  /*** Transform unfiltered box to k-space to prepare for filtering ***/
  fprintf(stderr, "begin initial ffts, time=%06.2f min\n", (double)clock()/CLOCKS_PER_SEC/60.0);
  fprintf(LOG, "begin initial ffts, time=%06.2f min\n", (double)clock()/CLOCKS_PER_SEC/60.0);
  plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)unfiltered_box, (fftwf_complex *)unfiltered_box, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);
  fftwf_cleanup();
  // remember to add the factor of VOLUME/TOT_NUM_PIXELS when converting from real space to k-space
  // Note: we will leave off factor of VOLUME, in anticipation of the inverse FFT below
  for (ct=0; ct<HII_KSPACE_NUM_PIXELS; ct++){
    unfiltered_box[ct] /= (float)HII_TOT_NUM_PIXELS;
  }
  fprintf(stderr, "end initial ffts, time=%06.2f min\n", (double)clock()/CLOCKS_PER_SEC/60.0);
  fprintf(LOG, "end initial ffts, time=%06.2f min\n", (double)clock()/CLOCKS_PER_SEC/60.0);


  /*** Create the z=0 non-linear density fields smoothed on scale R to be used in computing fcoll ***/
  R = L_FACTOR*BOX_LEN/(float)HII_DIM;
  R_factor = pow(R_XLy_MAX/R, 1/(float)NUM_FILTER_STEPS_FOR_Ts);
  //  R_factor = pow(E, log(HII_DIM)/(float)NUM_FILTER_STEPS_FOR_Ts);
  for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
    R_values[R_ct] = R;
    sigma_atR[R_ct] = sigma_z0(RtoM(R));
    fprintf(stderr, "Processing scale R= %06.2fMpc, time=%06.2f min\n", R, 
        (double)clock()/CLOCKS_PER_SEC/60.0);
    fprintf(LOG, "Processing scale R= %06.2fMpc, time=%06.2f min\n", R, 
        (double)clock()/CLOCKS_PER_SEC/60.0);
    if (! (delNL0[R_ct] = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS))){
      fprintf(stderr, "Error in memory allocation\nAborting...\n");
      fprintf(LOG, "Error in memory allocation\nAborting...\n");
      fclose(LOG); fclose(GLOBAL_EVOL);fftwf_free(box);  fftwf_free(unfiltered_box);
      for(ct=0; ct<R_ct; ct++)
    free(delNL0[ct]);
      destruct_heat();
      return -1;
    }

    // copy over unfiltered box
    memcpy(box, unfiltered_box, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    if (R_ct > 0){ // don't filter on cell size
      HII_filter(box, HEAT_FILTER, R);
    }

    // now fft back to real space
    plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)box, (float *)box, FFTW_ESTIMATE);
    fftwf_execute(plan);

    // copy over the values
#pragma omp parallel shared(delNL0, R_ct, box, growth_factor_z) private(i, j, k)
{
#pragma omp for
    for (i=0; i<HII_DIM; i++){
      for (j=0; j<HII_DIM; j++){
    for (k=0; k<HII_DIM; k++){
      delNL0[R_ct][HII_R_INDEX(i,j,k)] = *((float *) box + HII_R_FFT_INDEX(i,j,k));
      if (delNL0[R_ct][HII_R_INDEX(i,j,k)] < -1){ // correct for alliasing in the filtering step
        delNL0[R_ct][HII_R_INDEX(i,j,k)] = -1+FRACT_FLOAT_ERR;
      }
      // and linearly extrapolate to z=0
      delNL0[R_ct][HII_R_INDEX(i,j,k)] /= growth_factor_z; 
    }
      }
    }
}

    R *= R_factor;
  } //end for loop through the filter scales R
  
  fftwf_destroy_plan(plan);
  fftwf_cleanup();
  fftwf_free(box); fftwf_free(unfiltered_box);// we don't need this anymore

  // now lets allocate memory for our kinetic temperature and residual neutral fraction boxes
  if (!(Tk_box = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS))){
    fprintf(stderr, "Error in memory allocation for Tk box\nAborting...\n");
    fprintf(LOG, "Error in memory allocation for Tk box\nAborting...\n");
    fclose(LOG);fclose(GLOBAL_EVOL);
    for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
      free(delNL0[R_ct]);
    }
    destruct_heat();
    return -1;
  }
  if (!(x_e_box = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS))){
    fprintf(stderr, "Error in memory allocation for xe box\nAborting...\n");
    fprintf(LOG, "Error in memory allocation for xe box\nAborting...\n");
    fclose(LOG);  free(Tk_box);fclose(GLOBAL_EVOL);
    for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
      free(delNL0[R_ct]);
    }
    destruct_heat();
    return -1;
  }

  // and finally allocate memory for the spin temperature box
  if (!(Ts = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS))){
    fprintf(stderr, "Error in memory allocation for Ts box\nAborting...\n");
    fprintf(LOG, "Error in memory allocation for Ts box\nAborting...\n");
    fclose(LOG);  fclose(GLOBAL_EVOL);free(Tk_box); free(x_e_box);
    for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
      free(delNL0[R_ct]);
    }
    destruct_heat();
    return -1;
  }


  // and initialize to the boundary values at Z_HEAT_END
  if (!RESTART){ // we are not restarting
#pragma omp parallel shared(Tk_box, x_e_box, Tk_BC, xe_BC) private(ct)
{
#pragma omp for
    for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){
      Tk_box[ct] = Tk_BC;
      x_e_box[ct] = xe_BC;
    }
}
    x_e_ave = xe_BC;
    Tk_ave = Tk_BC;

    fprintf(stderr, "Starting at at z_max=%f, Tk=%f, x_e=%e\n", Z_HEAT_MAX, Tk_ave, x_e_ave);
  }
  else{ // we need to load the evolution files from the intermediate output
    // first Tk
#ifdef SHARP_CUTOFF
    sprintf(filename, "../Boxes/Ts_evolution/Tk_zprime%06.2f_L_X%.1e_alphaX%.1f_Mmin%.1e_zetaIon%.2f_Pop%i_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, M_MIN, HII_EFF_FACTOR, Pop, HII_DIM, BOX_LEN);
#else //SHARP_CUTOFF
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
    sprintf(filename, "../Boxes/Ts_evolution/Tk_zprime%06.2f_L_X%.1e_alphaX%.1f_f_star10_%06.4f_alpha_star%06.4f_f_esc10_%06.4f_alpha_esc%06.4f_t_star%06.4f_f_star10m%06.4f_f_esc10m%06.4f_L_Xm%.1e_alphaXm%.1f_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, T_AST, F_STAR10m, F_ESC10m, X_LUMINOSITYm, X_RAY_SPEC_INDEX_MINI, HII_DIM, BOX_LEN);
#else //INHOMO_FEEDBACK
#ifdef REION_SM
    sprintf(filename, "../Boxes/Ts_evolution/Tk_zprime%06.2f_L_X%.1e_alphaX%.1f_f_star10_%06.4f_alpha_star%06.4f_f_esc10_%06.4f_alpha_esc%06.4f_Msm_t_star%06.4f_f_star10m%06.4f_f_esc10m%06.4f_L_Xm%.1e_alphaXm%.1f_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, T_AST, F_STAR10m, F_ESC10m, X_LUMINOSITYm, X_RAY_SPEC_INDEX_MINI, HII_DIM, BOX_LEN);
#else //REION_SM
    sprintf(filename, "../Boxes/Ts_evolution/Tk_zprime%06.2f_L_X%.1e_alphaX%.1f_f_star10_%06.4f_alpha_star%06.4f_f_esc10_%06.4f_alpha_esc%06.4f_Mturn%.1e_t_star%06.4f_f_star10m%06.4f_f_esc10m%06.4f_L_Xm%.1e_alphaXm%.1f_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, T_AST, F_STAR10m, F_ESC10m, X_LUMINOSITYm, X_RAY_SPEC_INDEX_MINI, HII_DIM, BOX_LEN);
#endif //REION_SM
#endif //INHOMO_FEEDBACK
#else //MINI_HALO
    sprintf(filename, "../Boxes/Ts_evolution/Tk_zprime%06.2f_L_X%.1e_alphaX%.1f_f_star10_%06.4f_alpha_star%06.4f_f_esc10_%06.4f_alpha_esc%06.4f_Mturn%.1e_t_star%06.4f_Pop%i_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, T_AST, Pop, HII_DIM, BOX_LEN);
#endif //MINI_HALO
#endif //SHARP_CUTOFF
    if (!(F=fopen(filename, "rb"))){
      fprintf(stderr, "Ts.c: WARNING: Unable to open input file %s\nAborting\n", filename);
      fprintf(LOG, "Ts.c: WARNING: Unable to open input file %s\nAborting\n", filename);
      fclose(LOG); fclose(GLOBAL_EVOL); free(Tk_box); free(x_e_box); free(Ts);
      for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
        free(delNL0[R_ct]);
      }
      destruct_heat();
      return -1;
    }
    else{
      if (mod_fread(Tk_box, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
        fprintf(stderr, "Ts.c: Write error occured while reading Tk box.\nAborting\n");
        fprintf(LOG, "Ts.c: Write error occured while reading Tk box.\nAborting\n");
        fclose(LOG);  free(Tk_box); free(x_e_box); free(Ts);
        for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
          free(delNL0[R_ct]);
        }
        destruct_heat();
      }
      fclose(F);
    }
    // then xe_neutral
    // New in v2
#ifdef SHARP_CUTOFF
      sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_L_X%.1e_alphaX%.1f_Mmin%.1e_zetaIon%.2f_Pop%i_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, M_MIN, HII_EFF_FACTOR, Pop, HII_DIM, BOX_LEN);
#else //SHARP_CUTOFF
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
      sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_L_X%.1e_alphaX%.1f_f_star10_%06.4f_alpha_star%06.4f_f_esc10_%06.4f_alpha_esc%06.4f_t_star%06.4f_f_star10m%06.4f_f_esc10m%06.4f_L_Xm%.1e_alphaXm%.1f_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, T_AST, F_STAR10m, F_ESC10m, X_LUMINOSITYm, X_RAY_SPEC_INDEX_MINI, HII_DIM, BOX_LEN); 
#else //INHOMO_FEEDBACK
#ifdef REION_SM
      sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_L_X%.1e_alphaX%.1f_f_star10_%06.4f_alpha_star%06.4f_f_esc10_%06.4f_alpha_esc%06.4f_Msm_t_star%06.4f_f_star10m%06.4f_f_esc10m%06.4f_L_Xm%.1e_alphaXm%.1f_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, T_AST, F_STAR10m, F_ESC10m, X_LUMINOSITYm, X_RAY_SPEC_INDEX_MINI, HII_DIM, BOX_LEN); 
#else //REION_SM
      sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_L_X%.1e_alphaX%.1f_f_star10_%06.4f_alpha_star%06.4f_f_esc10_%06.4f_alpha_esc%06.4f_Mturn%.1e_t_star%06.4f_f_star10m%06.4f_f_esc10m%06.4f_L_Xm%.1e_alphaXm%.1f_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, T_AST, F_STAR10m, F_ESC10m, X_LUMINOSITYm, X_RAY_SPEC_INDEX_MINI, HII_DIM, BOX_LEN); 
#endif //REION_SM
#endif //INHOMO_FEEDBACK
#else //MINI_HALO
      sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_L_X%.1e_alphaX%.1f_f_star10_%06.4f_alpha_star%06.4f_f_esc10_%06.4f_alpha_esc%06.4f_Mturn%.1e_t_star%06.4f_Pop%i_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, T_AST, Pop, HII_DIM, BOX_LEN); 
#endif //MINI_HALO
#endif //SHARP_CUTOFF
    if (!(F=fopen(filename, "rb"))){
      fprintf(stderr, "Ts.c: WARNING: Unable to open output file %s\nAborting\n", filename);
      fprintf(LOG, "Ts.c: WARNING: Unable to open output file %s\nAborting\n", filename);
      fclose(LOG);  free(Tk_box); free(x_e_box); free(Ts);
      for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
        free(delNL0[R_ct]);
      }
      destruct_heat();
    }
    else{
      if (mod_fread(x_e_box, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
        fprintf(stderr, "Ts.c: Write error occured while reading xe box.\n");
        fprintf(LOG, "Ts.c: Write error occured while reading xe box.\n");
        fclose(LOG);  free(Tk_box); free(x_e_box); free(Ts);
        for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
          free(delNL0[R_ct]);
        }
        destruct_heat();
      }
      fclose(F);
    }
    Tk_ave = x_e_ave = 0;
    for (box_ct=0; box_ct<HII_TOT_NUM_PIXELS; box_ct++){
      Tk_ave += Tk_box[box_ct];
      x_e_ave += x_e_box[box_ct];
    }
    Tk_ave /= (double) HII_TOT_NUM_PIXELS;
    x_e_ave /= (double) HII_TOT_NUM_PIXELS;
    fprintf(stderr, "Rebooting from z'=%f output. <Tk> = %f. <xe> = %e\n", zp, Tk_ave, x_e_ave);
    
  }

  /***************    END INITIALIZATION   *********************************/

  /*********  FOR DEBUGGING, set IGM to be homogeneous for testing purposes 
    for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++)
    for (box_ct=0; box_ct<HII_TOT_NUM_PIXELS;box_ct++)
      delNL0[R_ct][box_ct] = 0;
    /*  *********/

  // main trapezoidal integral over z' (see eq. ? in Mesinger et al. 2009)
  if (!RESTART){
    Nsteps_zp = 0;
    zp = REDSHIFT*1.0001; //higher for rounding
    while (zp < Z_HEAT_MAX) { 
      Nsteps_zp += 1;
      zp = ((1+zp)*ZPRIME_STEP_FACTOR - 1);
    }
    prev_zp = Z_HEAT_MAX;
  }
  else{
    prev_zp_temp = zp;
    zp_temp = zp;
    Nsteps_zp = 0;
    zp = REDSHIFT*1.0001; //higher for rounding
    while (zp < Z_HEAT_MAX) { 
      Nsteps_zp += 1;
      zp = ((1+zp)*ZPRIME_STEP_FACTOR - 1);
    }
    prev_zp = Z_HEAT_MAX;
  }
  
  zp = ((1+zp)/ ZPRIME_STEP_FACTOR - 1);
  dzp = zp - prev_zp;
  if (RESTART == 1) {
    zp_temp = ((1+zp_temp)/ ZPRIME_STEP_FACTOR - 1);
    dzp = zp_temp - prev_zp_temp;
  }
  zp_ct=0;
  COMPUTE_Ts = 0;
  // New in v2
#ifndef SHARP_CUTOFF
  init_21cmMC_arrays();
  // Find the highest and lowest redshfit to initialise interpolation of the mean number of IGM ionizing photon per baryon
  determine_zpp_min = REDSHIFT*0.999;
  for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
      if (R_ct==0){
          prev_zpp = zp;
          prev_R = 0; 
      }    
      else{
          prev_zpp = zpp_edge[R_ct-1];
          prev_R = R_values[R_ct-1];
      }    
      zpp_edge[R_ct] = prev_zpp - (R_values[R_ct] - prev_R)*CMperMPC / drdz(prev_zpp); // cell size
      zpp = (zpp_edge[R_ct]+prev_zpp)*0.5; // average redshift value of shell: z'' + 0.5 * dz'' 
  }    
  determine_zpp_max = zpp*1.001;

  for (i=0; i<zpp_interp_points;i++) {
    zpp_interp_table[i] = determine_zpp_min + (determine_zpp_max - determine_zpp_min)*(float)i/((float)zpp_interp_points-1.0);
#ifdef MINI_HALO
    Mcrit_atom_interp_table[i] = atomic_cooling_threshold(zpp_interp_table[i]);
#ifdef INHOMO_FEEDBACK
    // NOTE: in Ts.c reionization feedback is ignored ifdef INHOMO_FEEDBACK
    Mcrit_RE_interp_table[i]   = 0.;
    M_MINa_interp_table[i]     = Mcrit_atom_interp_table[i];
#else //INHOMO_FEEDBACK
#ifdef REION_SM
    if(F = fopen("../Parameter_files/REION_SM.H", "r"))
      reading_reionization_SM13parameters(&REION_SM13_Z_RE, &REION_SM13_DELTA_Z_RE, &REION_SM13_DELTA_Z_SC, F);
    else
      estimating_reionization(ION_EFF_FACTOR, ION_EFF_FACTOR_MINI, ALPHA_STAR, F_STAR10, ALPHA_ESC, F_ESC10, F_STAR10m,
                              &REION_SM13_Z_RE, &REION_SM13_DELTA_Z_RE, &REION_SM13_DELTA_Z_SC);
    Mcrit_RE_interp_table[i]   = reionization_feedback(zpp_interp_table[i], REION_SM13_Z_RE, REION_SM13_DELTA_Z_RE, REION_SM13_DELTA_Z_SC);
#else //REION_SM
    Mcrit_RE_interp_table[i]   = M_TURN;
#endif //REION_SM
    M_MINa_interp_table[i]     = Mcrit_RE_interp_table[i] > Mcrit_atom_interp_table[i] ? Mcrit_RE_interp_table[i] : Mcrit_atom_interp_table[i];

    Mcrit_LW_interp_table[i]   = lyman_werner_threshold(zpp_interp_table[i]);
    M_MINm_interp_table[i]     = Mcrit_RE_interp_table[i] > Mcrit_LW_interp_table[i]   ? Mcrit_RE_interp_table[i] : Mcrit_LW_interp_table[i];

    if(M_MIN > M_MINa_interp_table[i])
      M_MIN = M_MINa_interp_table[i];
    if(M_MIN > M_MINm_interp_table[i])
      M_MIN = M_MINm_interp_table[i];
    fprintf(stderr, "zpp=%f, Mcrit_RE=%g, Mcrit_atom=%g, Mcrit_LW=%g, Mmin_a=%g, Mmin_m=%g\n",
            zpp_interp_table[i], Mcrit_RE_interp_table[i], Mcrit_atom_interp_table[i],
            Mcrit_LW_interp_table[i] ,M_MINa_interp_table[i],M_MINm_interp_table[i]);
#endif // INHOMO_FEEDBACK
#else //MINI_HALO
    M_MINa_interp_table[i]     = M_TURN;
#endif //MINI_HALO
  }

#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
#pragma omp parallel shared(log10_Mturn_interp_table) private(i)
{
#pragma omp for
  for (i=0; i<NMTURN; i++)
    log10_Mturn_interp_table[i] = 5. - 9e-8 + (double)i/((double)NMTURN-1.)*(5.+ 1.8e-7);
}
  M_MIN  = 1e5;
#else //INHOMO_FEEDBACK
  M_MIN /= 50;
#endif //INHOMO_FEEDBACK
#else //MINI_HALO
  M_MIN  = M_TURN/50;
#endif //MINI_HALO

  Mlim_Fstar  = Mass_limit_bisection(M_MIN, 1e16, ALPHA_STAR, F_STAR10);
  Mlim_Fesc   = Mass_limit_bisection(M_MIN, 1e16, ALPHA_ESC, F_ESC10);
#ifdef MINI_HALO
  Mlim_Fstarm = Mass_limit_bisection(M_MIN, 1e16, ALPHA_STAR, F_STAR10m);
#endif //MINI_HALO
  fprintf(stderr, "setting minimum mass for integral at %g\n",M_MIN);

  /* initialise interpolation of the mean number of IGM ionizing photon per baryon for global reionization.
     compute 'Nion_ST' corresponding to an array of redshift. */
  initialise_Nion_ST_spline(zpp_interp_points, zpp_interp_table, M_MIN, M_MINa_interp_table, ALPHA_STAR, ALPHA_ESC, F_STAR10, F_ESC10);
  fprintf(stderr, "\n Completed initialise Nion_ST, Time = %06.2f min \n",(double)clock()/CLOCKS_PER_SEC/60.0);
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
  initialise_Nion_ST_splinem(zpp_interp_points, zpp_interp_table, M_MIN, log10_Mturn_interp_table, Mcrit_atom_interp_table, ALPHA_STAR, F_STAR10m);
#else
  initialise_Nion_ST_splinem(zpp_interp_points, zpp_interp_table, M_MIN, M_MINm_interp_table, Mcrit_atom_interp_table, ALPHA_STAR, F_STAR10m);
#endif
  fprintf(stderr, "\n Completed initialise Nion_STm for mini halos, Time = %06.2f min \n",(double)clock()/CLOCKS_PER_SEC/60.0);
#endif //MINI_HALO
  
  /* initialise interpolation of the mean SFRD.
     compute 'Nion_ST' corresponding to an array of redshift, but assume f_{esc10} = 1 and \alpha_{esc} = 0. */
  initialise_SFRD_ST_spline(zpp_interp_points, zpp_interp_table, M_MIN, M_MINa_interp_table, ALPHA_STAR, F_STAR10);
  fprintf(stderr, "\n Completed initialise SFRD using Sheth-Tormen halo mass function, Time = %06.2f min \n",(double)clock()/CLOCKS_PER_SEC/60.0);
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
  initialise_SFRD_ST_splinem(zpp_interp_points, zpp_interp_table, M_MIN, log10_Mturn_interp_table, Mcrit_atom_interp_table, ALPHA_STAR, F_STAR10m);
#else
  initialise_SFRD_ST_splinem(zpp_interp_points, zpp_interp_table, M_MIN, M_MINm_interp_table, Mcrit_atom_interp_table, ALPHA_STAR, F_STAR10m);
#endif
  fprintf(stderr, "\n Completed initialise SFRD for mini halos using Sheth-Tormen halo mass function, Time = %06.2f min \n",(double)clock()/CLOCKS_PER_SEC/60.0);
#endif //MINI_HALO
  
  // initialise redshift table corresponding to all the redshifts to initialise interpolation for the conditional mass function.
  zp_table = zp;
  counter = 0;
  for (i=0; i<Nsteps_zp; i++) {
    for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
        if (R_ct==0){
            prev_zpp = zp_table;
            prev_R = 0; 
        }    
        else{
            prev_zpp = zpp_edge[R_ct-1];
            prev_R = R_values[R_ct-1];
        }    
        zpp_edge[R_ct] = prev_zpp - (R_values[R_ct] - prev_R)*CMperMPC / drdz(prev_zpp); // cell size
        zpp = (zpp_edge[R_ct]+prev_zpp)*0.5; // average redshift value of shell: z'' + 0.5 * dz''
        redshift_interp_table[counter] = zpp;
        counter += 1;
    }    
    prev_zp = zp_table;
    zp_table = ((1+prev_zp) / ZPRIME_STEP_FACTOR - 1);
  } 

  /* generate a table for interpolation of the SFRD using the conditional mass function, as functions of 
  filtering scale, redshift and overdensity.
     See eq. (8) in Park et al. (2018)
     Note that at a given zp, zpp values depends on the filtering scale R. */
  // Note from YQ: due to the initial configuration left from v2, I will construct M_MINa_interp_table and M_MINm_interp_table inside these two
  // functions (which is fine because we are only going to use them once, although there are unnecessary duplicated calculations...)
#ifdef INHOMO_FEEDBACK
  initialise_SFRD_Conditional_table(Nsteps_zp,NUM_FILTER_STEPS_FOR_Ts,redshift_interp_table,R_values, M_MIN, ALPHA_STAR, F_STAR10);
#else
#ifdef REION_SM
  initialise_SFRD_Conditional_table(Nsteps_zp,NUM_FILTER_STEPS_FOR_Ts,redshift_interp_table,R_values, M_MIN, ALPHA_STAR, F_STAR10, REION_SM13_Z_RE, REION_SM13_DELTA_Z_RE, REION_SM13_DELTA_Z_SC);
#else
  initialise_SFRD_Conditional_table(Nsteps_zp,NUM_FILTER_STEPS_FOR_Ts,redshift_interp_table,R_values, M_MIN, M_TURN, ALPHA_STAR, F_STAR10);
#endif
#endif
  fprintf(stderr, "\n Generated the table of SFRD using conditional mass function = %06.2f min \n",(double)clock()/CLOCKS_PER_SEC/60.0);

#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
  initialise_SFRD_Conditional_tablem(Nsteps_zp,NUM_FILTER_STEPS_FOR_Ts,redshift_interp_table,R_values, M_MIN, log10_Mturn_interp_table, ALPHA_STAR, F_STAR10m);
#else
#ifdef REION_SM
  initialise_SFRD_Conditional_tablem(Nsteps_zp,NUM_FILTER_STEPS_FOR_Ts,redshift_interp_table,R_values, M_MIN, ALPHA_STAR, F_STAR10m, REION_SM13_Z_RE, REION_SM13_DELTA_Z_RE, REION_SM13_DELTA_Z_SC);
#else
  initialise_SFRD_Conditional_tablem(Nsteps_zp,NUM_FILTER_STEPS_FOR_Ts,redshift_interp_table,R_values, M_MIN, M_TURN, ALPHA_STAR, F_STAR10m);
#endif
#endif

#endif

  fprintf(stderr, "\n Generated the table of SFRD using conditional mass function for mini halos = %06.2f min \n",(double)clock()/CLOCKS_PER_SEC/60.0);
#endif //SHARP_CUTOFF

  if (RESTART == 1){
    zp = zp_temp;
    prev_zp = prev_zp_temp;
  }

  counter = 0;
  while (zp > REDSHIFT){

#ifndef SHARP_CUTOFF
    // New in v2: initialise interpolation of SFRD over zpp and overdensity.
    arr_num = NUM_FILTER_STEPS_FOR_Ts*counter; // New
    fprintf(stderr, "\n [z=%.3f] constructing the interpolation table...", zp);
/*#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
#pragma omp parallel shared(log10_overdense_low_table,log10_SFRD_z_low_table,arr_num,Overdense_high_table,SFRD_z_high_table,second_derivs_Nion_zpp, log10_overdense_low_table_Mturn,Overdense_high_table_Mturn,log10_SFRD_z_low_tablem, SFRD_z_high_tablem, second_derivs_Nion_zppm, SFRDLow_zpp_spline, SFRDLow_zpp_splinem) private(i)
#else
#pragma omp parallel shared(log10_overdense_low_table,log10_SFRD_z_low_table,arr_num,Overdense_high_table,SFRD_z_high_table,second_derivs_Nion_zpp, log10_SFRD_z_low_tablem, SFRD_z_high_tablem, second_derivs_Nion_zppm, SFRDLow_zpp_spline, SFRDLow_zpp_splinem) private(i)
#endif
#pragma omp parallel shared(log10_overdense_low_table,log10_SFRD_z_low_table,arr_num,Overdense_high_table,SFRD_z_high_table,second_derivs_Nion_zpp, SFRDLow_zpp_spline) private(i)
#endif*/
{
//#pragma omp for
//NOTE THAT this part seems to be not thread safe!!
    for (i=0; i<NUM_FILTER_STEPS_FOR_Ts; i++) {
      gsl_spline_init(SFRDLow_zpp_spline[i], log10_overdense_low_table, log10_SFRD_z_low_table[arr_num + i], NSFR_low);
      spline(Overdense_high_table-1,SFRD_z_high_table[arr_num + i]-1,NSFR_high,0,0,second_derivs_Nion_zpp[i]-1); 
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
      gsl_spline2d_init(SFRDLow_zpp_splinem[i], log10_overdense_low_table, log10_overdense_low_table_Mturn, log10_SFRD_z_low_tablem[arr_num + i], NSFR_low, NMTURN);
      spline2d(Overdense_high_table,Overdense_high_table_Mturn,SFRD_z_high_tablem[arr_num + i],NSFR_high,NMTURN,second_derivs_Nion_zppm[i]); 
#else
      gsl_spline_init(SFRDLow_zpp_splinem[i], log10_overdense_low_table, log10_SFRD_z_low_tablem[arr_num + i], NSFR_low);
      spline(Overdense_high_table-1,SFRD_z_high_tablem[arr_num + i]-1,NSFR_high,0,0,second_derivs_Nion_zppm[i]-1); 
#endif
#endif
    }
}
    fprintf(stderr, "done = %06.2f min \n",(double)clock()/CLOCKS_PER_SEC/60.0);
#endif

    // check if we will next compute the spin temperature (i.e. if this is the final zp step)
#ifndef Ts_verbose
    if (((1+zp) / ZPRIME_STEP_FACTOR) < (REDSHIFT+1))
#endif
    {COMPUTE_Ts = 1;}

    // check if we are in the really high z regime before the first stars..
#ifdef SHARP_CUTOFF
    if (FgtrM(zp, M_MIN) < 1e-15 )
      NO_LIGHT = 1;
    else
      NO_LIGHT = 0;
#else //SHARP_CUTOFF
    Nion_ST_z(zp,&(Splined_Nion_ST_zp));
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
    log10_Mcrit_mol = log10(lyman_werner_threshold(zp, 0.));
    sprintf(filename, "../Boxes/J_21_LW_z%06.2f_L_X%.1e_alphaX%.1f_f_star%06.4f_alpha_star%06.4f_f_esc%06.4f_alpha_esc%06.4f_t_star%06.4f_f_star10m%06.4f_f_esc10m%06.4f_L_Xm%.1e_alphaXm%.1f_%i_%.0fMpc", prev_zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, T_AST, F_STAR10m, F_ESC10m, X_LUMINOSITYm, X_RAY_SPEC_INDEX_MINI, HII_DIM, BOX_LEN); 
    if (F=fopen(filename, "rb")){
      if (mod_fread(J_21_LW, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
        fprintf(stderr, "Ts.c: Read error occured while reading J_21_LW box!\n");
        return -1;
      } 
      fclose(F);
    }
    else{
      for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++)
        J_21_LW[ct] = 0.0;
    }
    log10_Mcrit_LW_ave = 0;
#pragma omp parallel shared(zp, log10_Mcrit_LW_unfiltered, J_21_LW) private(i,j,k) reduction(+:log10_Mcrit_LW_ave)
{
#pragma omp for
  for (i=0; i<HII_DIM; i++){
    for (j=0; j<HII_DIM; j++){
      for (k=0; k<HII_DIM; k++){
        *((float *)log10_Mcrit_LW_unfiltered + HII_R_FFT_INDEX(i,j,k)) = log10(lyman_werner_threshold(zp, J_21_LW[HII_R_INDEX(i,j,k)]));
        log10_Mcrit_LW_ave += *((float *)log10_Mcrit_LW_unfiltered + HII_R_FFT_INDEX(i,j,k));
      }
    }
  }
}
    log10_Mcrit_LW_ave /= HII_TOT_NUM_PIXELS;
    Mcrit_atom_glob  = atomic_cooling_threshold(zp);
    Nion_ST_zm(zp,log10_Mcrit_LW_ave,&(Splined_Nion_ST_zpm));

    // NEED TO FILTER Mcrit_LW!!!
    /*** Transform unfiltered box to k-space to prepare for filtering ***/
    plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)log10_Mcrit_LW_unfiltered, (fftwf_complex *)log10_Mcrit_LW_unfiltered, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    fftwf_cleanup();
    // remember to add the factor of VOLUME/TOT_NUM_PIXELS when converting from real space to k-space
    // Note: we will leave off factor of VOLUME, in anticipation of the inverse FFT below
#pragma omp parallel shared(log10_Mcrit_LW_unfiltered) private(ct)
{
#pragma omp for
    for (ct=0; ct<HII_KSPACE_NUM_PIXELS; ct++)
      log10_Mcrit_LW_unfiltered[ct] /= (float)HII_TOT_NUM_PIXELS;
}

#else //INHOMO_FEEDBACK
    Nion_ST_zm(zp,&(Splined_Nion_ST_zpm));
#endif //INHOMO_FEEDBACK
    if ( Splined_Nion_ST_zp + Splined_Nion_ST_zpm < 1e-15 )
#else
    if ( Splined_Nion_ST_zp < 1e-15 )
#endif
      NO_LIGHT = 1;
    else
      NO_LIGHT = 0;
#endif //SHARP_CUTOFF

    //New in v2
#ifdef SHARP_CUTOFF
    filling_factor_of_HI_zp = 1 - HII_EFF_FACTOR * FgtrM_st(zp, M_MIN) / (1.0 - x_e_ave);
#else
    // let's initialize an array of redshifts (z'') corresponding to the 
    // far edge of the dz'' filtering shells
    // and the corresponding minimum halo scale, sigma_Tmin, 
    // as well as an array of the frequency integrals
    fprintf(stderr, "Initializing look-up tables. Time=%06.2f min\n", (double)clock()/CLOCKS_PER_SEC/60.0);
    fprintf(LOG, "Initializing look-up tables. Time=%06.2f min\n", (double)clock()/CLOCKS_PER_SEC/60.0);
    time(&start_time);
#ifdef MINI_HALO
    filling_factor_of_HI_zp = 1 - (ION_EFF_FACTOR * Splined_Nion_ST_zp + ION_EFF_FACTOR_MINI * Splined_Nion_ST_zpm) / (1.0 - x_e_ave); // fcoll including f_esc
#else
    filling_factor_of_HI_zp = 1 - ION_EFF_FACTOR * Splined_Nion_ST_zp / (1.0 - x_e_ave); // fcoll including f_esc
#endif
#endif
    if (filling_factor_of_HI_zp > 1) filling_factor_of_HI_zp=1;
    if (filling_factor_of_HI_zp < 0) filling_factor_of_HI_zp=0;

    for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
      if (R_ct==0){
        prev_zpp = zp;
        prev_R = 0;
      }
      else{
        prev_zpp = zpp_edge[R_ct-1];
        prev_R = R_values[R_ct-1];
      }
    
      zpp_edge[R_ct] = prev_zpp - (R_values[R_ct] - prev_R)*CMperMPC / drdz(prev_zpp); // cell size
      zpp = (zpp_edge[R_ct]+prev_zpp)*0.5; // average redshift value of shell: z'' + 0.5 * dz''
      if (zpp - redshift_interp_table[arr_num+R_ct] > 1e-3) fprintf(stderr, "zpp = %.4f, zpp_array = %.4f\n", zpp, redshift_interp_table[arr_num+R_ct]);
#ifdef SHARP_CUTOFF 
      sigma_Tmin[R_ct] =  sigma_z0(M_MIN); // In v2 sigma_Tmin doesn't need to be an array, just a constant.
#endif

#ifdef INHOMO_FEEDBACK
      // copy over unfiltered box
      memcpy(log10_Mcrit_LW_filtered, log10_Mcrit_LW_unfiltered, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
      if (R_ct > 0){// don't filter on cell size
        HII_filter(log10_Mcrit_LW_filtered, HEAT_FILTER, R_values[R_ct]);
	  }

      // now fft back to real space
      plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)log10_Mcrit_LW_filtered, (float *)log10_Mcrit_LW_filtered, FFTW_ESTIMATE);
      fftwf_execute(plan);
      fftwf_destroy_plan(plan);
      fftwf_cleanup();

      // I don't know how box_ct_increment works... so just do the same copying thing...
#pragma omp parallel shared(log10_Mcrit_LW, log10_Mcrit_LW_filtered, log10_Mcrit_mol) private(i, j, k)
{
#pragma omp for
      for (i=0; i<HII_DIM; i++){
        for (j=0; j<HII_DIM; j++){
          for (k=0; k<HII_DIM; k++){
            log10_Mcrit_LW[HII_R_INDEX(i,j,k)] = *((float *) log10_Mcrit_LW_filtered + HII_R_FFT_INDEX(i,j,k));
            if(log10_Mcrit_LW[HII_R_INDEX(i,j,k)] < log10_Mcrit_mol)
              log10_Mcrit_LW[HII_R_INDEX(i,j,k)] = log10_Mcrit_mol;
            if (log10_Mcrit_LW[HII_R_INDEX(i,j,k)] > 10)
              log10_Mcrit_LW[HII_R_INDEX(i,j,k)] = 10;
          }
        }
      }
}
#endif

      // let's now normalize the total collapse fraction so that the mean is the
      // Sheth-Torman collapse fraction
      fcoll_R = 0;
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
	  log10_Mcrit_LW_ave = 0;
#endif
      fcollm_R = 0;
#endif
      sample_ct=0;
#ifdef SHARP_CUTOFF
#pragma omp parallel shared(zpp, sigma_Tmin, R_ct, delNL0, sigma_atR) private(box_ct) reduction(+:sample_ct, fcoll_R)
#else //SHARP_CUTOFF
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
#pragma omp parallel shared(zpp,arr_num, R_ct, delNL0, Overdense_high_table, SFRD_z_high_table, SFRD_z_high_tablem, second_derivs_Nion_zpp, second_derivs_Nion_zppm, log10_Mcrit_LW, Overdense_high_table_Mturn, SFRDLow_zpp_spline, SFRDLow_zpp_spline_acc, SFRDLow_zpp_splinem, SFRDLow_zpp_spline_accm, SFRDLow_zpp_spline_accm_Mturn) private(box_ct, delNL_zpp, fcoll, Splined_Fcoll, fcollm, Splined_Fcollm) reduction(+:sample_ct, fcoll_R, fcollm_R, log10_Mcrit_LW_ave)
#else //INHOMO_FEEDBACK
#pragma omp parallel shared(zpp,arr_num, R_ct, delNL0, Overdense_high_table, SFRD_z_high_table, SFRD_z_high_tablem, second_derivs_Nion_zpp, second_derivs_Nion_zppm, SFRDLow_zpp_spline, SFRDLow_zpp_spline_acc, SFRDLow_zpp_splinem, SFRDLow_zpp_spline_accm) private(box_ct, delNL_zpp, fcoll, Splined_Fcoll, fcollm, Splined_Fcollm) reduction(+:sample_ct, fcoll_R, fcollm_R)
#endif //INHOMO_FEEDBACK
#else //MINI_HALO
#pragma omp parallel shared(zpp,arr_num, R_ct, delNL0, Overdense_high_table, SFRD_z_high_table, second_derivs_Nion_zpp, SFRDLow_zpp_spline, SFRDLow_zpp_spline_acc) private(box_ct, delNL_zpp, fcoll, Splined_Fcoll) reduction(+:sample_ct, fcoll_R)
#endif //MINI_HALO
#endif //SHARP_CUTOFF
{
#pragma omp for
      for (box_ct=0; box_ct<HII_TOT_NUM_PIXELS; box_ct+=box_ct_increment){
        sample_ct++;
        // New in v2
#ifdef SHARP_CUTOFF
        fcoll_R += sigmaparam_FgtrM_bias(zpp, sigma_Tmin[R_ct], delNL0[R_ct][box_ct], sigma_atR[R_ct]);
#else //SHARP_CUTOFF
        delNL_zpp    = delNL0[R_ct][box_ct]*dicke(zpp);
        //---------- interpolation for fcoll starts ----------
        // Here 'fcoll' is not the collpased fraction, but leave this name as is to simplify the variable name.
#ifdef INHOMO_FEEDBACK
		log10_Mcrit_LW_ave += log10_Mcrit_LW[box_ct];
#endif
        if (delNL_zpp < 1.5){
          if (delNL_zpp < -1.) {
            fcoll = 0;
#ifdef MINI_HALO
            fcollm = 0;
#endif
          }    
          else {
            fcoll = gsl_spline_eval(SFRDLow_zpp_spline[R_ct], log10(delNL_zpp+1.), SFRDLow_zpp_spline_acc[R_ct]);
            fcoll = pow(10., fcoll);
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
            fcollm = gsl_spline2d_eval(SFRDLow_zpp_splinem[R_ct], log10(delNL_zpp+1.), log10_Mcrit_LW[box_ct], SFRDLow_zpp_spline_accm[R_ct], SFRDLow_zpp_spline_accm_Mturn[R_ct]);
#else 
            fcollm = gsl_spline_eval(SFRDLow_zpp_splinem[R_ct], log10(delNL_zpp+1.), SFRDLow_zpp_spline_accm[R_ct]);
#endif //INHOMO_FEEDBACK
            fcollm = pow(10., fcollm);
#endif //MINI_HALO
          }
        }
        else {
          if (delNL_zpp < 0.99*Deltac) {
            // Usage of 0.99*Deltac arises due to the fact that close to the critical density, the collapsed fraction becomes a little unstable
            // However, such densities should always be collapsed, so just set f_coll to unity. 
            // Additionally, the fraction of points in this regime relative to the entire simulation volume is extremely small.
            splint(Overdense_high_table-1,SFRD_z_high_table[arr_num+R_ct]-1,second_derivs_Nion_zpp[R_ct]-1,NSFR_high,delNL_zpp,&(fcoll));
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
            splint2d(Overdense_high_table,Overdense_high_table_Mturn,SFRD_z_high_tablem[arr_num+R_ct],second_derivs_Nion_zppm[R_ct],NSFR_high,NMTURN,delNL_zpp,log10_Mcrit_LW[box_ct],&(fcollm));
#else //INHOMO_FEEDBACK
            splint(Overdense_high_table-1,SFRD_z_high_tablem[arr_num+R_ct]-1,second_derivs_Nion_zppm[R_ct]-1,NSFR_high,delNL_zpp,&(fcollm));
#endif //INHOMO_FEEDBACK
#endif //MINI_HALO
          }
          else {
            fcoll = 1.;
#ifdef MINI_HALO
            fcollm = 1.;
#endif
          }    
        }
		if (fcoll > 1.)
			fcoll = 1.;
        Splined_Fcoll = fcoll > 0 ? fcoll : 1e-40;
#ifdef MINI_HALO
		if (fcollm > 1.)
			fcollm = 1.;
        Splined_Fcollm = fcollm > 0 ? fcollm : 1e-40;
#endif
        //---------- interpolation for fcoll is done ----------
        fcoll_R += Splined_Fcoll;
#ifdef MINI_HALO
        fcollm_R += Splined_Fcollm;
#endif
#endif //SHARP_CUTOFF
      }
}

      fcoll_R /= (double) sample_ct;
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
	  log10_Mcrit_LW_ave /= (double) sample_ct;
#endif
      fcollm_R /= (double) sample_ct;
#endif

#ifdef SHARP_CUTOFF
      ST_over_PS[R_ct] = FgtrM_st(zpp, M_MIN) / fcoll_R;
#else //SHARP_CUTOFF
      SFRD_ST_z(zpp,&(Splined_SFRD_ST_zpp));
      ST_over_PS[R_ct] = Splined_SFRD_ST_zpp / fcoll_R; 
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
      // TODO: now for mini halos the ST/PS normalization factor is the mean of ST when feedback from RE is ignored
      // and LW is considered homogeneous (averaging Mcrit_LW over cells) to the mean of PS when inhomogeneous LW feedback is considered.
      // This means that it might not be as good as in previous versions, because now the fluctuation in fcoll is also 
      // due to feedbacks in M_MINm. e.g. if the distribution of Mcrit_LW is too large
      // this is not interpolated value
      SFRD_ST_zm(zpp,log10_Mcrit_LW_ave,&(Splined_SFRD_ST_zppm));
#else
      SFRD_ST_zm(zpp,&(Splined_SFRD_ST_zppm));
#endif //INHOMO_FEEDBACK
      ST_over_PSm[R_ct] = Splined_SFRD_ST_zppm / fcollm_R; 
#endif //MINI_HALO
#endif //SHARP_CUTOFF

//#ifdef DEBUG_ON
#ifdef SHARP_CUTOFF
      fprintf(LOG, "R=%06.2fMpc, ST/PS=%g, mean_ST=%g, mean_ps=%g, ratios of mean=%g\n",
         R_values[R_ct],
         ST_over_PS[R_ct], 
         FgtrM_st(zpp, M_MIN), 
         FgtrM(zpp, M_MIN),
         FgtrM_st(zpp, M_MIN)/FgtrM(zpp, M_MIN)
         );
#else //SHARP_CUTOFF
#ifdef MINI_HALO
      fprintf(LOG, "scale R        = %06.2fMpc\n \
ST/PS          = (atomic:%g, molecular:%g)\n \
mean_ST        = (atomic:%g, molecular:%g)\n \
mean_PS        = %g\n \
ratios of mean = (atomic:%g, molecular:%g)\n", 
         R_values[R_ct],
         ST_over_PS[R_ct],ST_over_PSm[R_ct], 
         Splined_SFRD_ST_zpp,Splined_SFRD_ST_zppm,
         FgtrM(zpp, M_MIN),
         Splined_SFRD_ST_zpp/FgtrM(zpp, M_MIN),Splined_SFRD_ST_zppm/FgtrM(zpp, M_MIN)
         );
#else //MINI_HALO
      fprintf(LOG, "R=%06.2fMpc, ST/PS=%g, mean_ST=%g, mean_ps=%g, ratios of mean=%g\n", 
         R_values[R_ct],
         ST_over_PS[R_ct], 
         Splined_SFRD_ST_zpp,
         FgtrM(zpp, M_MIN),
         Splined_SFRD_ST_zpp/FgtrM(zpp, M_MIN)
         );
#endif //MINI_HALO
#endif //SHARP_CUTOFF
//#endif //DEBUG_ON

#ifdef SHARP_CUTOFF
      lower_int_limit = FMAX(nu_tau_one(zp, zpp, x_e_ave, filling_factor_of_HI_zp), NU_X_THRESH);
#else //SHARP_CUTOFF
#ifdef MINI_HALO
      lower_int_limit = FMAX(nu_tau_one(zp, zpp, x_e_ave, filling_factor_of_HI_zp, ION_EFF_FACTOR, ION_EFF_FACTOR_MINI), NU_X_THRESH);
#else //MINI_HALO
      lower_int_limit = FMAX(nu_tau_one(zp, zpp, x_e_ave, filling_factor_of_HI_zp, ION_EFF_FACTOR), NU_X_THRESH);
#endif //MINI_HALO
#endif //SHARP_CUTOFF

/***************  PARALLELIZED LOOP ******************************************************************/
      // set up frequency integral table for later interpolation for the cell's x_e value
#ifdef MINI_HALO
#pragma omp parallel shared(freq_int_heat_tbl, freq_int_ion_tbl, COMPUTE_Ts, freq_int_lya_tbl, freq_int_heat_tblm, freq_int_ion_tblm, freq_int_lya_tblm, zp, R_ct, x_e_ave, x_int_XHII, x_int_Energy, x_int_fheat, x_int_n_Lya, x_int_nion_HI, x_int_nion_HeI, x_int_nion_HeII, lower_int_limit) private(x_e_ct)
#else
#pragma omp parallel shared(freq_int_heat_tbl, freq_int_ion_tbl, COMPUTE_Ts, freq_int_lya_tbl, zp, R_ct, x_e_ave, x_int_XHII, x_int_Energy, x_int_fheat, x_int_n_Lya, x_int_nion_HI, x_int_nion_HeI, x_int_nion_HeII, lower_int_limit) private(x_e_ct)
#endif
{
#pragma omp for
  for (x_e_ct = 0; x_e_ct < x_int_NXHII; x_e_ct++){
    freq_int_heat_tbl[x_e_ct][R_ct] = integrate_over_nu(zp, x_int_XHII[x_e_ct], lower_int_limit, 0);
    freq_int_ion_tbl[x_e_ct][R_ct] = integrate_over_nu(zp, x_int_XHII[x_e_ct], lower_int_limit, 1);
    if (COMPUTE_Ts)
      freq_int_lya_tbl[x_e_ct][R_ct] = integrate_over_nu(zp, x_int_XHII[x_e_ct], lower_int_limit, 2);
#ifdef MINI_HALO
    freq_int_heat_tblm[x_e_ct][R_ct] = integrate_over_num(zp, x_int_XHII[x_e_ct], lower_int_limit, 0);
    freq_int_ion_tblm[x_e_ct][R_ct] = integrate_over_num(zp, x_int_XHII[x_e_ct], lower_int_limit, 1);
    if (COMPUTE_Ts)
      freq_int_lya_tblm[x_e_ct][R_ct] = integrate_over_num(zp, x_int_XHII[x_e_ct], lower_int_limit, 2);
#endif
  }
} // end omp declaration
/***************  END PARALLELIZED LOOP ******************************************************************/

  // and create the sum over Lya transitions from direct Lyn flux
      sum_lyn[R_ct] = 0;
#ifdef MINI_HALO
      sum_lynm[R_ct] = 0;
#ifdef INHOMO_FEEDBACK
      sum_lyLWn[R_ct] = 0;
      sum_lyLWnm[R_ct] = 0;
#endif
#endif
      for (n_ct=NSPEC_MAX; n_ct>=2; n_ct--){
        if (zpp > zmax(zp, n_ct))
          continue;

        nuprime = nu_n(n_ct)*(1+zpp)/(1.0+zp);
#ifdef MINI_HALO
        sum_lyn[R_ct]  += frecycle(n_ct) * spectral_emissivity(nuprime, 0, 2);
        sum_lynm[R_ct] += frecycle(n_ct) * spectral_emissivity(nuprime, 0, 3);
#ifdef INHOMO_FEEDBACK
        nu_nplus1 = nu_n(n_ct + 1);
        if (nuprime < NU_LW_THRESH)
          nuprime = NU_LW_THRESH;
        sum_lyLWn[R_ct]  += spectral_emissivity(nuprime, 3, 2);
        sum_lyLWnm[R_ct] += spectral_emissivity(nuprime, 3, 3);
#endif
#else
        sum_lyn[R_ct] += frecycle(n_ct) * spectral_emissivity(nuprime, 0);
#endif
      }
    } // end loop over R_ct filter steps

    time(&curr_time);
    fprintf(stderr, "Finishing initializing look-up tables.  It took %06.2f min on the main thread. Time elapsed (total for all threads)=%06.2f\n", difftime(curr_time, start_time)/60.0, (double)clock()/CLOCKS_PER_SEC/60.0);
    fprintf(LOG, "Finishing initializing look-up tables.  It took %06.2f min on the main thread. Total time elapsed (total for all threads)=%06.2f\n", difftime(curr_time, start_time)/60.0, (double)clock()/CLOCKS_PER_SEC/60.0);
    fflush(NULL);

    // scroll through each cell and update the temperature and residual ionization fraction
    growth_factor_zp = dicke(zp);
    dgrowth_factor_dzp = ddicke_dz(zp);
    dt_dzp = dtdz(zp);
    // New in v2
    // Conversion of the input bolometric luminosity to a ZETA_X, as used to be used in Ts.c
    // Conversion here means the code otherwise remains the same as the original Ts.c
    if(fabs(X_RAY_SPEC_INDEX - 1.0) < 0.000001) {
      Luminosity_conversion_factor = NU_X_THRESH * log( NU_X_BAND_MAX/NU_X_THRESH );
      Luminosity_conversion_factor = 1./Luminosity_conversion_factor;
    }    
    else {
      Luminosity_conversion_factor = pow( NU_X_BAND_MAX , 1. - X_RAY_SPEC_INDEX ) - pow( NU_X_THRESH , 1. - X_RAY_SPEC_INDEX ) ;
      Luminosity_conversion_factor = 1./Luminosity_conversion_factor;
      Luminosity_conversion_factor *= pow( NU_X_THRESH, - X_RAY_SPEC_INDEX )*(1 - X_RAY_SPEC_INDEX);
    }    
    // Finally, convert to the correct units. NU_over_EV*hplank as only want to divide by eV -> erg (owing to the definition of Luminosity)
    Luminosity_conversion_factor *= (3.1556226e7)/(hplank);
    const_zp_prefactor = ( X_LUMINOSITY * Luminosity_conversion_factor ) / NU_X_THRESH * C 
             * F_STAR10 * OMb * RHOcrit * pow(CMperMPC, -3) * pow(1+zp, X_RAY_SPEC_INDEX+3);
#ifdef MINI_HALO
    if(fabs(X_RAY_SPEC_INDEX_MINI - 1.0) < 0.000001) {
      Luminosity_conversion_factorm = NU_X_THRESH * log( NU_X_BAND_MAX/NU_X_THRESH );
      Luminosity_conversion_factorm = 1./Luminosity_conversion_factorm;
    }    
    else {
      Luminosity_conversion_factorm = pow( NU_X_BAND_MAX , 1. - X_RAY_SPEC_INDEX_MINI ) - pow( NU_X_THRESH , 1. - X_RAY_SPEC_INDEX_MINI ) ;
      Luminosity_conversion_factorm = 1./Luminosity_conversion_factorm;
      Luminosity_conversion_factorm *= pow( NU_X_THRESH, - X_RAY_SPEC_INDEX_MINI )*(1 - X_RAY_SPEC_INDEX_MINI);
    }
    // Finally, convert to the correct units. NU_over_EV*hplank as only want to divide by eV -> erg (owing to the definition of Luminosity)
    Luminosity_conversion_factorm *= (3.1556226e7)/(hplank);
    const_zp_prefactorm = ( X_LUMINOSITYm * Luminosity_conversion_factorm ) / NU_X_THRESH * C 
             * F_STAR10m * OMb * RHOcrit * pow(CMperMPC, -3) * pow(1+zp, X_RAY_SPEC_INDEX_MINI+3);
#endif


    /********  LOOP THROUGH BOX *************/
    fprintf(stderr, "Looping through box at z'=%f, time elapsed  (total for all threads)= %06.2f min\n", zp, (double)clock()/CLOCKS_PER_SEC/60.0);
    fprintf(LOG, "Looping through box at z'=%f, time elapsed  (total for all threads)= %06.2f min\n", zp, (double)clock()/CLOCKS_PER_SEC/60.0);
    fflush(NULL);
    time(&start_time);
    for (ct=0; ct<NUMCORES; ct++){
      J_alpha_threads[ct] = xalpha_threads[ct] = Xheat_threads[ct] = Xion_threads[ct] = 0;
#ifdef INHOMO_FEEDBACK
	  J_LW_threads[ct] = 0;
#endif
	}
    /***************  PARALLELIZED LOOP ******************************************************************/
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
#pragma omp parallel shared(COMPUTE_Ts, Tk_box, x_e_box, x_e_ave, delNL0, freq_int_heat_tbl, freq_int_ion_tbl, freq_int_lya_tbl,freq_int_heat_tblm, freq_int_ion_tblm, freq_int_lya_tblm, zp, dzp, Ts, x_int_XHII, x_int_Energy, x_int_fheat, x_int_n_Lya, x_int_nion_HI, x_int_nion_HeI, x_int_nion_HeII, growth_factor_zp, dgrowth_factor_dzp, NO_LIGHT, zpp_edge, sigma_atR, sigma_Tmin, ST_over_PS, ST_over_PSm, sum_lyn,sum_lynm,sum_lyLWn, sum_lyLWnm, const_zp_prefactor, const_zp_prefactorm, M_MIN_at_z, M_MIN_at_zp, dt_dzp, J_alpha_threads, J_LW_threads, xalpha_threads, Xheat_threads, Xion_threads, M_MIN, R_values, Mcrit_atom_glob,log10_Mcrit_LW_ave,J_21_LW,arr_num) private(box_ct, ans, xHII_call, R_ct, curr_delNL0, m_xHII_low, m_xHII_high, freq_int_heat, freq_int_ion, freq_int_lya, freq_int_heatm, freq_int_ionm, freq_int_lyam, dansdz, J_alpha_tot, curr_xalpha)
#else
#pragma omp parallel shared(COMPUTE_Ts, Tk_box, x_e_box, x_e_ave, delNL0, freq_int_heat_tbl, freq_int_ion_tbl, freq_int_lya_tbl,freq_int_heat_tblm, freq_int_ion_tblm, freq_int_lya_tblm, zp, dzp, Ts, x_int_XHII, x_int_Energy, x_int_fheat, x_int_n_Lya, x_int_nion_HI, x_int_nion_HeI, x_int_nion_HeII, growth_factor_zp, dgrowth_factor_dzp, NO_LIGHT, zpp_edge, sigma_atR, sigma_Tmin, ST_over_PS, ST_over_PSm, sum_lyn,sum_lynm, const_zp_prefactor, const_zp_prefactorm, M_MIN_at_z, M_MIN_at_zp, dt_dzp, J_alpha_threads, xalpha_threads, Xheat_threads, Xion_threads,arr_num) private(box_ct, ans, xHII_call, R_ct, curr_delNL0, m_xHII_low, m_xHII_high, freq_int_heat, freq_int_ion, freq_int_lya, freq_int_heatm, freq_int_ionm, freq_int_lyam, dansdz, J_alpha_tot, curr_xalpha)
#endif
#else
#pragma omp parallel shared(COMPUTE_Ts, Tk_box, x_e_box, x_e_ave, delNL0, freq_int_heat_tbl, freq_int_ion_tbl, freq_int_lya_tbl, zp, dzp, Ts, x_int_XHII, x_int_Energy, x_int_fheat, x_int_n_Lya, x_int_nion_HI, x_int_nion_HeI, x_int_nion_HeII, growth_factor_zp, dgrowth_factor_dzp, NO_LIGHT, zpp_edge, sigma_atR, sigma_Tmin, ST_over_PS, sum_lyn, const_zp_prefactor, M_MIN_at_z, M_MIN_at_zp, dt_dzp, J_alpha_threads, xalpha_threads, Xheat_threads, Xion_threads,arr_num) private(box_ct, ans, xHII_call, R_ct, curr_delNL0, m_xHII_low, m_xHII_high, freq_int_heat, freq_int_ion, freq_int_lya, dansdz, J_alpha_tot, curr_xalpha)
#endif
    {
#pragma omp for
      
      for (box_ct=0; box_ct<HII_TOT_NUM_PIXELS; box_ct++){
        if (!COMPUTE_Ts && (Tk_box[box_ct] > MAX_TK)) //just leave it alone and go to next value
        continue;

        // set to current values before updating
        ans[0] = x_e_box[box_ct];
        ans[1] = Tk_box[box_ct];

        /*
        if (DEBUG_ON){
          if (isnan(ans[0]))
            fprintf(stderr, "Problem at cell %llu, x_e=%e\n", box_ct, ans[0]);
          if (isnan(ans[1]))
            fprintf(stderr, "Problem at cell %llu, Tk=%e\n", box_ct, ans[1]);
        }
        */

        xHII_call = x_e_box[box_ct];

        // Check if ionized fraction is within boundaries; if not, adjust to be within
        if (xHII_call > x_int_XHII[x_int_NXHII-1]*0.999) 
          xHII_call = x_int_XHII[x_int_NXHII-1]*0.999;
        else if (xHII_call < x_int_XHII[0])
          xHII_call = 1.001*x_int_XHII[0];
 
        //interpolate to correct nu integral value based on the cell's ionization state
        for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
          curr_delNL0[R_ct] = delNL0[R_ct][box_ct];
          m_xHII_low = locate_xHII_index(xHII_call);
          m_xHII_high = m_xHII_low + 1;

          // heat
          freq_int_heat[R_ct] = (freq_int_heat_tbl[m_xHII_high][R_ct] - 
                                freq_int_heat_tbl[m_xHII_low][R_ct]) / 
                                (x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
          freq_int_heat[R_ct] *= (xHII_call - x_int_XHII[m_xHII_low]);
          freq_int_heat[R_ct] += freq_int_heat_tbl[m_xHII_low][R_ct];

          // ionization
          freq_int_ion[R_ct] = (freq_int_ion_tbl[m_xHII_high][R_ct] - 
                               freq_int_ion_tbl[m_xHII_low][R_ct]) / 
                               (x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
          freq_int_ion[R_ct] *= (xHII_call - x_int_XHII[m_xHII_low]);
          freq_int_ion[R_ct] += freq_int_ion_tbl[m_xHII_low][R_ct];

#ifdef MINI_HALO
          // heat
          freq_int_heatm[R_ct] = (freq_int_heat_tblm[m_xHII_high][R_ct] - 
                                 freq_int_heat_tblm[m_xHII_low][R_ct]) / 
                                 (x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
          freq_int_heatm[R_ct] *= (xHII_call - x_int_XHII[m_xHII_low]);
          freq_int_heatm[R_ct] += freq_int_heat_tblm[m_xHII_low][R_ct];

          // ionization
          freq_int_ionm[R_ct] = (freq_int_ion_tblm[m_xHII_high][R_ct] - 
                                freq_int_ion_tblm[m_xHII_low][R_ct]) / 
                                (x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
          freq_int_ionm[R_ct] *= (xHII_call - x_int_XHII[m_xHII_low]);
          freq_int_ionm[R_ct] += freq_int_ion_tblm[m_xHII_low][R_ct];
#endif
    
          // lya
          if (COMPUTE_Ts){
            freq_int_lya[R_ct] = (freq_int_lya_tbl[m_xHII_high][R_ct] - 
                                 freq_int_lya_tbl[m_xHII_low][R_ct]) / 
                                 (x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
            freq_int_lya[R_ct] *= (xHII_call - x_int_XHII[m_xHII_low]);
            freq_int_lya[R_ct] += freq_int_lya_tbl[m_xHII_low][R_ct];
#ifdef MINI_HALO
            freq_int_lyam[R_ct] = (freq_int_lya_tblm[m_xHII_high][R_ct] - 
                                  freq_int_lya_tblm[m_xHII_low][R_ct]) / 
                                  (x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
            freq_int_lyam[R_ct] *= (xHII_call - x_int_XHII[m_xHII_low]);
            freq_int_lyam[R_ct] += freq_int_lya_tblm[m_xHII_low][R_ct];
          }
#endif
        }

        /********  finally compute the redshift derivatives *************/
#ifdef MINI_HALO
        evolveInt(zp, arr_num, curr_delNL0, freq_int_heat, freq_int_ion, freq_int_lya,
                  freq_int_heatm, freq_int_ionm, freq_int_lyam,
                  COMPUTE_Ts, ans, dansdz);//, M_TURN,ALPHA_STAR,F_STAR10,T_AST);
#else
        evolveInt(zp, arr_num, curr_delNL0, freq_int_heat, freq_int_ion, freq_int_lya,
                  COMPUTE_Ts, ans, dansdz);//, M_TURN,ALPHA_STAR,F_STAR10,T_AST);
#endif
 
        //update quantities
        x_e_box[box_ct] += dansdz[0] * dzp; // remember dzp is negative
        if (x_e_box[box_ct] > 1) // can do this late in evolution if dzp is too large
          x_e_box[box_ct] = 1 - FRACT_FLOAT_ERR;
        else if (x_e_box[box_ct] < 0)
          x_e_box[box_ct] = 0;
        if (Tk_box[box_ct] < MAX_TK)
          Tk_box[box_ct] += dansdz[1] * dzp;
        if (Tk_box[box_ct]<0) // spurious bahaviour of the trapazoidalintegrator. generally overcooling in underdensities
          Tk_box[box_ct] = T_cmb*(1+zp);
       
        if (COMPUTE_Ts){
          J_alpha_tot = dansdz[2]; //not really d/dz, but the lya flux
          Ts[box_ct] = get_Ts(zp, curr_delNL0[0]*growth_factor_zp,
                       Tk_box[box_ct], x_e_box[box_ct], J_alpha_tot, &curr_xalpha);
          J_alpha_threads[omp_get_thread_num()] += J_alpha_tot;
          xalpha_threads[omp_get_thread_num()] += curr_xalpha;
          Xheat_threads[omp_get_thread_num()] += dansdz[3];
          Xion_threads[omp_get_thread_num()] += dansdz[4];
#ifdef INHOMO_FEEDBACK
          J_21_LW[box_ct] = dansdz[5];
          J_LW_threads[omp_get_thread_num()] += dansdz[5];
#endif
        }
      }
    } // end parallelization pragma

/***************  END PARALLELIZED LOOP ******************************************************************/
    time(&curr_time);
    fprintf(stderr, "End scrolling through the box, which took %06.2f min\n", difftime(curr_time, start_time)/60.0);
    fprintf(LOG, "End scrolling through the box, which took %06.2f min\n", difftime(curr_time, start_time)/60.0);
    fflush(NULL);

    // compute new average values
    x_e_ave = 0; Tk_ave = 0; Ts_ave = 0; J_alpha_ave = 0; xalpha_ave = 0; Xheat_ave=0; Xion_ave=0;
#ifdef INHOMO_FEEDBACK
    J_LW_ave = 0;
#endif
    for (box_ct=0; box_ct<HII_TOT_NUM_PIXELS; box_ct++){
      x_e_ave += x_e_box[box_ct];
      Tk_ave += Tk_box[box_ct];
      if (COMPUTE_Ts)
    Ts_ave += Ts[box_ct];
    }
    for (ct=0; ct<NUMCORES; ct++){
      J_alpha_ave += J_alpha_threads[ct];
      xalpha_ave += xalpha_threads[ct];
      Xheat_ave += Xheat_threads[ct];
      Xion_ave += Xion_threads[ct];
#ifdef INHOMO_FEEDBACK
      J_LW_ave += J_LW_threads[ct];
#endif
    }
    Ts_ave /= (double)HII_TOT_NUM_PIXELS;
    x_e_ave /= (double)HII_TOT_NUM_PIXELS;
    Tk_ave /= (double)HII_TOT_NUM_PIXELS;
    J_alpha_ave /= (double)HII_TOT_NUM_PIXELS;
    xalpha_ave /= (double)HII_TOT_NUM_PIXELS;
    Xheat_ave /= (double)HII_TOT_NUM_PIXELS;
    Xion_ave /= (double)HII_TOT_NUM_PIXELS;
#ifdef INHOMO_FEEDBACK
    J_LW_ave /= (double)HII_TOT_NUM_PIXELS;
#endif
    // write to global evolution file
#ifdef INHOMO_FEEDBACK
    fprintf(GLOBAL_EVOL, "%f\t%f\t%f\t%e\t%f\t%f\t%e\t%e\t%e\t%e\t%e\n", zp, filling_factor_of_HI_zp, Tk_ave, x_e_ave, Ts_ave, T_cmb*(1+zp), J_alpha_ave, xalpha_ave, Xheat_ave, Xion_ave, J_LW_ave * 1e-21);
#else
    fprintf(GLOBAL_EVOL, "%f\t%f\t%f\t%e\t%f\t%f\t%e\t%e\t%e\t%e\n", zp, filling_factor_of_HI_zp, Tk_ave, x_e_ave, Ts_ave, T_cmb*(1+zp), J_alpha_ave, xalpha_ave, Xheat_ave, Xion_ave);
#endif
    fflush(NULL);

    // output these intermediate boxes
#ifndef Ts_verbose
    if (++zp_ct >= 10) // print every 10th z' evolution step, in case we need to restart
#endif //Ts_verbose
    {
      zp_ct=0;
      fprintf(stderr, "Writting the intermediate output at zp = %.4f, <Tk>=%f, <x_e>=%e\n", zp, Tk_ave, x_e_ave);
      fprintf(LOG, "Writting the intermediate output at zp = %.4f, <Tk>=%f, <x_e>=%e\n", zp, Tk_ave, x_e_ave);
      fflush(NULL);

      // first Tk
        // New v2
#ifdef SHARP_CUTOFF
      sprintf(filename, "../Boxes/Ts_evolution/Tk_zprime%06.2f_L_X%.1e_alphaX%.1f_Mmin%.1e_zetaIon%.2f_Pop%i_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, M_MIN, HII_EFF_FACTOR, Pop, HII_DIM, BOX_LEN);
#else //SHARP_CUTOFF
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
      sprintf(filename, "../Boxes/Ts_evolution/Tk_zprime%06.2f_L_X%.1e_alphaX%.1f_f_star10_%06.4f_alpha_star%06.4f_f_esc10_%06.4f_alpha_esc%06.4f_t_star%06.4f_f_star10m%06.4f_f_esc10m%06.4f_L_Xm%.1e_alphaXm%.1f_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, T_AST, F_STAR10m, F_ESC10m, X_LUMINOSITYm, X_RAY_SPEC_INDEX_MINI, HII_DIM, BOX_LEN);
#else //INHOMO_FEEDBACK
#ifdef REION_SM
      sprintf(filename, "../Boxes/Ts_evolution/Tk_zprime%06.2f_L_X%.1e_alphaX%.1f_f_star10_%06.4f_alpha_star%06.4f_f_esc10_%06.4f_alpha_esc%06.4f_Msm_t_star%06.4f_f_star10m%06.4f_f_esc10m%06.4f_L_Xm%.1e_alphaXm%.1f_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, T_AST, F_STAR10m, F_ESC10m, X_LUMINOSITYm, X_RAY_SPEC_INDEX_MINI, HII_DIM, BOX_LEN);
#else //REION_SM
      sprintf(filename, "../Boxes/Ts_evolution/Tk_zprime%06.2f_L_X%.1e_alphaX%.1f_f_star10_%06.4f_alpha_star%06.4f_f_esc10_%06.4f_alpha_esc%06.4f_Mturn%.1e_t_star%06.4f_f_star10m%06.4f_f_esc10m%06.4f_L_Xm%.1e_alphaXm%.1f_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, T_AST, F_STAR10m, F_ESC10m, X_LUMINOSITYm, X_RAY_SPEC_INDEX_MINI, HII_DIM, BOX_LEN);
#endif //REION_SM
#endif //INHOMO_FEEDBACK
#else //MINI_HALO
      sprintf(filename, "../Boxes/Ts_evolution/Tk_zprime%06.2f_L_X%.1e_alphaX%.1f_f_star10_%06.4f_alpha_star%06.4f_f_esc10_%06.4f_alpha_esc%06.4f_Mturn%.1e_t_star%06.4f_Pop%i_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, T_AST, Pop, HII_DIM, BOX_LEN);
#endif //MINI_HALO
#endif //SHARP_CUTOFF
      if (!(F=fopen(filename, "wb"))){
        fprintf(stderr, "Ts.c: WARNING: Unable to open output file %s\n", filename);
        fprintf(LOG, "Ts.c: WARNING: Unable to open output file %s\n", filename);
      }
      else{
        if (mod_fwrite(Tk_box, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
          fprintf(stderr, "Ts.c: Write error occured while writting Tk box.\n");
          fprintf(LOG, "Ts.c: Write error occured while writting Tk box.\n");
        }
        fclose(F);
      }
      // then xe_neutral
        // New in v2
#ifdef SHARP_CUTOFF
      sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_L_X%.1e_alphaX%.1f_Mmin%.1e_zetaIon%.2f_Pop%i_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, M_MIN, HII_EFF_FACTOR, Pop, HII_DIM, BOX_LEN);
#else //SHARP_CUTOFF
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
      sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_L_X%.1e_alphaX%.1f_f_star10_%06.4f_alpha_star%06.4f_f_esc10_%06.4f_alpha_esc%06.4f_t_star%06.4f_f_star10m%06.4f_f_esc10m%06.4f_L_Xm%.1e_alphaXm%.1f_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, T_AST, F_STAR10m, F_ESC10m, X_LUMINOSITYm, X_RAY_SPEC_INDEX_MINI, HII_DIM, BOX_LEN);
#else //INHOMO_FEEDBACK
#ifdef REION_SM
      sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_L_X%.1e_alphaX%.1f_f_star10_%06.4f_alpha_star%06.4f_f_esc10_%06.4f_alpha_esc%06.4f_Msm_t_star%06.4f_f_star10m%06.4f_f_esc10m%06.4f_L_Xm%.1e_alphaXm%.1f_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, T_AST, F_STAR10m, F_ESC10m, X_LUMINOSITYm, X_RAY_SPEC_INDEX_MINI, HII_DIM, BOX_LEN);
#else //REION_SM
      sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_L_X%.1e_alphaX%.1f_f_star10_%06.4f_alpha_star%06.4f_f_esc10_%06.4f_alpha_esc%06.4f_Mturn%.1e_t_star%06.4f_f_star10m%06.4f_f_esc10m%06.4f_L_Xm%.1e_alphaXm%.1f_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, T_AST, F_STAR10m, F_ESC10m, X_LUMINOSITYm, X_RAY_SPEC_INDEX_MINI, HII_DIM, BOX_LEN);
#endif //REION_SM
#endif //INHOMO_FEEDBACK
#else //MINI_HALO
      sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_L_X%.1e_alphaX%.1f_f_star10_%06.4f_alpha_star%06.4f_f_esc10_%06.4f_alpha_esc%06.4f_Mturn%.1e_t_star%06.4f_Pop%i_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, T_AST, Pop, HII_DIM, BOX_LEN);
#endif //MINI_HALO
#endif //SHARP_CUTOFF
      if (!(F=fopen(filename, "wb"))){
        fprintf(stderr, "Ts.c: WARNING: Unable to open output file %s\n", filename);
        fprintf(LOG, "Ts.c: WARNING: Unable to open output file %s\n", filename);
      }
      else{
        if (mod_fwrite(x_e_box, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
          fprintf(stderr, "Ts.c: Write error occured while writting Tk box.\n");
          fprintf(LOG, "Ts.c: Write error occured while writting Tk box.\n");
        }
        fclose(F);
      }
    }

    // and the spin temperature if desired
    if ( COMPUTE_Ts ){
      // New in v2
#ifdef SHARP_CUTOFF
      sprintf(filename, "../Boxes/Ts_z%06.2f_L_X%.1e_alphaX%.1f_TvirminX%.1e_zetaIon%.2f_Pop%i_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, M_TURN, HII_EFF_FACTOR, Pop, HII_DIM, BOX_LEN); 
#else //SHARP_CUTOFF
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
      sprintf(filename, "../Boxes/Ts_z%06.2f_L_X%.1e_alphaX%.1f_f_star%06.4f_alpha_star%06.4f_f_esc%06.4f_alpha_esc%06.4f_t_star%06.4f_f_star10m%06.4f_f_esc10m%06.4f_L_Xm%.1e_alphaXm%.1f_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, T_AST, F_STAR10m, F_ESC10m, X_LUMINOSITYm, X_RAY_SPEC_INDEX_MINI, HII_DIM, BOX_LEN); 
#else //INHOMO_FEEDBACK
#ifdef REION_SM
      sprintf(filename, "../Boxes/Ts_z%06.2f_L_X%.1e_alphaX%.1f_f_star%06.4f_alpha_star%06.4f_f_esc%06.4f_alpha_esc%06.4f_Msm_t_star%06.4f_f_star10m%06.4f_f_esc10m%06.4f_L_Xm%.1e_alphaXm%.1f_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, T_AST, F_STAR10m, F_ESC10m, X_LUMINOSITYm, X_RAY_SPEC_INDEX_MINI, HII_DIM, BOX_LEN); 
#else //REION_SM
      sprintf(filename, "../Boxes/Ts_z%06.2f_L_X%.1e_alphaX%.1f_f_star%06.4f_alpha_star%06.4f_f_esc%06.4f_alpha_esc%06.4f_Mturn%.1e_t_star%06.4f_f_star10m%06.4f_f_esc10m%06.4f_L_Xm%.1e_alphaXm%.1f_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, T_AST, F_STAR10m, F_ESC10m, X_LUMINOSITYm, X_RAY_SPEC_INDEX_MINI, HII_DIM, BOX_LEN); 
#endif //REION_SM
#endif //INHOMO_FEEDBACK
#else //MINI_HALO
      sprintf(filename, "../Boxes/Ts_z%06.2f_L_X%.1e_alphaX%.1f_f_star%06.4f_alpha_star%06.4f_f_esc%06.4f_alpha_esc%06.4f_Mturn%.1e_t_star%06.4f_Pop%i_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, T_AST, Pop, HII_DIM, BOX_LEN); 
#endif //MINI_HALO
#endif //SHARP_CUTOFF
      if (!(F=fopen(filename, "wb"))){
        fprintf(stderr, "Ts.c: WARNING: Unable to open output file %s\n", filename);
        fprintf(LOG, "Ts.c: WARNING: Unable to open output file %s\n", filename);
      }
      else{
        if (mod_fwrite(Ts, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
          fprintf(stderr, "Ts.c: Write error occured while writting Tk box.\n");
          fprintf(LOG, "Ts.c: Write error occured while writting Tk box.\n");
        }
        fclose(F);
      }
#ifdef INHOMO_FEEDBACK
      sprintf(filename, "../Boxes/J_21_LW_z%06.2f_L_X%.1e_alphaX%.1f_f_star%06.4f_alpha_star%06.4f_f_esc%06.4f_alpha_esc%06.4f_t_star%06.4f_f_star10m%06.4f_f_esc10m%06.4f_L_Xm%.1e_alphaXm%.1f_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, T_AST, F_STAR10m, F_ESC10m, X_LUMINOSITYm, X_RAY_SPEC_INDEX_MINI, HII_DIM, BOX_LEN); 
      if (!(F=fopen(filename, "wb"))){
        fprintf(stderr, "Ts.c: WARNING: Unable to open output file %s\n", filename);
        fprintf(LOG, "Ts.c: WARNING: Unable to open output file %s\n", filename);
      }
      else{
        if (mod_fwrite(J_21_LW, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
          fprintf(stderr, "Ts.c: Write error occured while writting Tk box.\n");
          fprintf(LOG, "Ts.c: Write error occured while writting Tk box.\n");
        }
        fclose(F);
      }
#endif
    }

    prev_zp = zp;
    zp = ((1+prev_zp) / ZPRIME_STEP_FACTOR - 1);
    dzp = zp - prev_zp;
    counter += 1;
  } // end main integral loop over z'

  //deallocate
  fclose(LOG); fclose(GLOBAL_EVOL); free(Tk_box); free(x_e_box); free(Ts);
#ifndef SHARP_CUTOFF
  destroy_21cmMC_arrays();
#ifdef INHOMO_FEEDBACK
  fftwf_free(log10_Mcrit_LW_unfiltered);
  fftwf_free(log10_Mcrit_LW_filtered);
  free(log10_Mcrit_LW);
  free(J_21_LW);
#endif
#endif
  for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++)
    free(delNL0[R_ct]);
  destruct_heat();
  return 0;
}
