#include "bubble_helper_progs.c"
#include "heating_helper_progs.c"

/*
  USAGE: find_HII_bubbles [-p <num of processors>] <redshift> [<previous redshift>]
  [<stellar fraction for 10^10 Msun halos> <power law index for stellar fraction halo mass scaling> 
   <escape fraction for 10^10 Msun halos> <power law index for escape fraction halo mass scaling>
   <turn-over scale for the duty cycle of galaxies, in units of halo mass>] [<Soft band X-ray luminosity>]

   final optional parameter is MFP (depricated in v2)

  Program FIND_HII_BUBBLES generates the ionization field in Boxes/

  The minimum command line parameter is the redshift and, if the INHOMO_RECO flag is set, also the previous redshift.
  (the previous redshift is only used if computing recombinations)

  If optional parameters are not passed, then the 
  values in ANAL_PARAMS.H are used.

  NOTE: the optional argument of thread number including the -p flag, MUST
  be the first two arguments.  If these are omitted, num_threads defaults
  to NUMCORES in INIT_PARAMS.H

  See ANAL_PARAMS.H for updated parameters and algorithms!

  Author: Andrei Mesinger
  Date: 01/10/07

  Support for Alpha power law scaling of zeta with halo mass added by Bradley Greig, 2016
*/

float *Fcoll;
#ifdef MINI_HALO
float *Fcollm;
#endif //MINI_HALO

void init_21cmMC_arrays() { // defined in Cosmo_c_files/ps.c
    
    Overdense_spline_SFR = calloc(NSFR_high,sizeof(float)); // New in v2
    Nion_spline = calloc(NSFR_high,sizeof(float));
    second_derivs_Nion = calloc(NSFR_high,sizeof(float));
#ifdef MINI_HALO
    Nion_splinem = calloc(NSFR_high,sizeof(float));
    second_derivs_Nionm = calloc(NSFR_high,sizeof(float));
#endif //MINI_HALO
    xi_SFR = calloc((NGL_SFR+1),sizeof(float));
    wi_SFR = calloc((NGL_SFR+1),sizeof(float));
}

void destroy_21cmMC_arrays() {
    
    free(Mass_Spline);
    free(Sigma_Spline);
    free(dSigmadm_Spline);
    free(second_derivs_sigma);
    free(second_derivs_dsigma);

    free(Overdense_spline_SFR); // New in v2
    free(Nion_spline);
    free(second_derivs_Nion);
#ifdef MINI_HALO
    free(Nion_splinem);
    free(second_derivs_Nionm);
#endif //MINI_HALO
    free(xi_SFR);
    free(wi_SFR);
    gsl_spline_free(NionLow_spline);
    gsl_interp_accel_free(NionLow_spline_acc);
#ifdef MINI_HALO
    gsl_spline_free(NionLow_splinem);
    gsl_interp_accel_free(NionLow_spline_accm);
#endif //MINI_HALO
}


FILE *LOG;
unsigned long long SAMPLING_INTERVAL = (((unsigned long long)(HII_TOT_NUM_PIXELS/1.0e6)) + 1); //used to sample density field to compute mean collapsed fraction

int parse_arguments(int argc, char ** argv, int * num_th, int * arg_offset, float * F_STAR10, 
#ifdef SHARP_CUTOFF
                    float * ALPHA_STAR, float * F_ESC10, float * ALPHA_ESC, float * M_TURN, float * T_AST, double * X_LUMINOSITY,
#else //SHARP_CUTOFF
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
                    float * ALPHA_STAR, float * F_ESC10, float * ALPHA_ESC, float * T_AST, double * X_LUMINOSITY,
#else  //INHOMO_FEEDBACK
#ifdef REION_SM
                    float * ALPHA_STAR, float * F_ESC10, float * ALPHA_ESC, float * T_AST, double * X_LUMINOSITY,
#else //REION_SM
                    float * ALPHA_STAR, float * F_ESC10, float * ALPHA_ESC, float * M_TURN, float * T_AST, double * X_LUMINOSITY,
#endif //REION_SM
#endif //INHOMO_FEEDBACK
                    float * F_STAR10m, float * F_ESC10m, double * X_LUMINOSITYm,
#else //MINI_HALO
                    float * ALPHA_STAR, float * F_ESC10, float * ALPHA_ESC, float * M_TURN, float * T_AST, double * X_LUMINOSITY,
#endif //MINI_HALO
#endif //SHARP_CUTOFF
                    float * MFP, float * REDSHIFT, float * PREV_REDSHIFT){

  int min_argc = 2;

#ifdef INHOMO_RECO
  min_argc++;
#endif //INHOMO_RECO
  
  // check if user wants to specify the multi threading
  if ((argc > min_argc) && (argv[1][0]=='-') && ((argv[1][1]=='p') || (argv[1][1]=='P'))){
    // user specified num proc
    *num_th = atoi(argv[2]);
    *arg_offset = 2;
  }
  else{
    *num_th = NUMCORES;
    *arg_offset = 0;
  }

  // parse the remaining arguments
#ifdef SHARP_CUTOFF
  if (argc == (*arg_offset + min_argc)){
    *MFP = R_BUBBLE_MAX;
    *F_STAR10 = STELLAR_BARYON_FRAC;     
    *ALPHA_STAR = STELLAR_BARYON_PL;
    *F_ESC10 = ESC_FRAC;
    *ALPHA_ESC = ESC_PL;
    *M_TURN = M_TURNOVER;
    *T_AST = t_STAR;
    *X_LUMINOSITY = 0;
  }
  else{ return 0;} // format is not allowed
#else //SHARP_CUTOFF
#ifdef USE_TS_IN_21CM
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
  if (argc == (*arg_offset + min_argc+9)){
    *F_STAR10 = atof(argv[*arg_offset + min_argc]);
    *ALPHA_STAR = atof(argv[*arg_offset + min_argc+1]);
    *F_ESC10 = atof(argv[*arg_offset + min_argc+2]);
    *ALPHA_ESC = atof(argv[*arg_offset + min_argc+3]);
    *T_AST = atof(argv[*arg_offset + min_argc+4]);
    *X_LUMINOSITY = pow(10.,atof(argv[*arg_offset + min_argc+5]));
    *F_STAR10m = atof(argv[*arg_offset + min_argc+6]);
    *F_ESC10m = atof(argv[*arg_offset + min_argc+7]);
    *X_LUMINOSITYm = pow(10.,atof(argv[*arg_offset + min_argc+8]));
  }
  else if (argc == (*arg_offset + min_argc)){ //These parameters give the result which is the same with the default model.
    *F_STAR10 = STELLAR_BARYON_FRAC;       
    *ALPHA_STAR = STELLAR_BARYON_PL;
    *F_ESC10 = ESC_FRAC;
    *ALPHA_ESC = ESC_PL;
    *T_AST = t_STAR;
    *X_LUMINOSITY = pow(10.,L_X);
    *F_STAR10m = STELLAR_BARYON_FRAC_MINI;
    *F_ESC10m = ESC_FRAC_MINI;
    *X_LUMINOSITYm = pow(10.,L_X_MINI);
  }
  else{ return 0;} // format is not allowed
#else //INHOMO_FEEDBACK
#ifdef REION_SM
  if (argc == (*arg_offset + min_argc+9)){
    *F_STAR10 = atof(argv[*arg_offset + min_argc]);
    *ALPHA_STAR = atof(argv[*arg_offset + min_argc+1]);
    *F_ESC10 = atof(argv[*arg_offset + min_argc+2]);
    *ALPHA_ESC = atof(argv[*arg_offset + min_argc+3]);
    *T_AST = atof(argv[*arg_offset + min_argc+4]);
    *X_LUMINOSITY = pow(10.,atof(argv[*arg_offset + min_argc+5]));
    *F_STAR10m = atof(argv[*arg_offset + min_argc+6]);
    *F_ESC10m = atof(argv[*arg_offset + min_argc+7]);
    *X_LUMINOSITYm = pow(10.,atof(argv[*arg_offset + min_argc+8]));
  }
  else if (argc == (*arg_offset + min_argc)){ //These parameters give the result which is the same with the default model.
    *F_STAR10 = STELLAR_BARYON_FRAC;       
    *ALPHA_STAR = STELLAR_BARYON_PL;
    *F_ESC10 = ESC_FRAC;
    *ALPHA_ESC = ESC_PL;
    *T_AST = t_STAR;
    *X_LUMINOSITY = pow(10.,L_X);
    *F_STAR10m = STELLAR_BARYON_FRAC_MINI;
    *F_ESC10m = ESC_FRAC_MINI;
    *X_LUMINOSITYm = pow(10.,L_X_MINI);
  }
  else{ return 0;} // format is not allowed
#else //REION_SM
  if (argc == (*arg_offset + min_argc+10)){
    *F_STAR10 = atof(argv[*arg_offset + min_argc]);
    *ALPHA_STAR = atof(argv[*arg_offset + min_argc+1]);
    *F_ESC10 = atof(argv[*arg_offset + min_argc+2]);
    *ALPHA_ESC = atof(argv[*arg_offset + min_argc+3]);
    *M_TURN = atof(argv[*arg_offset + min_argc+4]); 
    *T_AST = atof(argv[*arg_offset + min_argc+5]);
    *X_LUMINOSITY = pow(10.,atof(argv[*arg_offset + min_argc+6]));
    *F_STAR10m = atof(argv[*arg_offset + min_argc+7]);
    *F_ESC10m = atof(argv[*arg_offset + min_argc+8]);
    *X_LUMINOSITYm = pow(10.,atof(argv[*arg_offset + min_argc+9]));
  }
  else if (argc == (*arg_offset + min_argc)){ //These parameters give the result which is the same with the default model.
    *F_STAR10 = STELLAR_BARYON_FRAC;       
    *ALPHA_STAR = STELLAR_BARYON_PL;
    *F_ESC10 = ESC_FRAC;
    *ALPHA_ESC = ESC_PL;
    *M_TURN = M_TURNOVER;
    *T_AST = t_STAR;
    *X_LUMINOSITY = pow(10.,L_X);
    *F_STAR10m = STELLAR_BARYON_FRAC_MINI;
    *F_ESC10m = ESC_FRAC_MINI;
    *X_LUMINOSITYm = pow(10.,L_X_MINI);
  }
  else{ return 0;} // format is not allowed
#endif //REION_SM
#endif //INHOMO_FEEDBACK
#else // MINI_HALO
  if (argc == (*arg_offset + min_argc+7)){
    *F_STAR10 = atof(argv[*arg_offset + min_argc]);
    *ALPHA_STAR = atof(argv[*arg_offset + min_argc+1]);
    *F_ESC10 = atof(argv[*arg_offset + min_argc+2]);
    *ALPHA_ESC = atof(argv[*arg_offset + min_argc+3]);
    *M_TURN = atof(argv[*arg_offset + min_argc+4]); 
    *T_AST = atof(argv[*arg_offset + min_argc+5]);
    *X_LUMINOSITY = pow(10.,atof(argv[*arg_offset + min_argc+6]));
  }
  else if (argc == (*arg_offset + min_argc)){ //These parameters give the result which is the same with the default model.
    *F_STAR10 = STELLAR_BARYON_FRAC;       
    *ALPHA_STAR = STELLAR_BARYON_PL;
    *F_ESC10 = ESC_FRAC;
    *ALPHA_ESC = ESC_PL;
    *M_TURN = M_TURNOVER;
    *T_AST = t_STAR;
    *X_LUMINOSITY = pow(10.,L_X);
  }
  else{ return 0;} // format is not allowed
#endif //MINI_HALO
#else //USE_TS_IN_21CM
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
  if (argc == (*arg_offset + min_argc+6)){
    *F_STAR10 = atof(argv[*arg_offset + min_argc]);
    *ALPHA_STAR = atof(argv[*arg_offset + min_argc+1]);
    *F_ESC10 = atof(argv[*arg_offset + min_argc+2]);
    *ALPHA_ESC = atof(argv[*arg_offset + min_argc+3]);
    *T_AST = t_STAR;
    *X_LUMINOSITY = 0;
    *F_STAR10m = atof(argv[*arg_offset + min_argc+4]);
    *F_ESC10m = atof(argv[*arg_offset + min_argc+5]);
    *X_LUMINOSITYm = 0;
  }
  else if (argc == (*arg_offset + min_argc)){ //These parameters give the result which is the same with the default model.
    *F_STAR10 = STELLAR_BARYON_FRAC;        
    *ALPHA_STAR = STELLAR_BARYON_PL;
    *F_ESC10 = ESC_FRAC;
    *ALPHA_ESC = ESC_PL;
    *T_AST = t_STAR;
    *X_LUMINOSITY = 0;
    *F_STAR10m = STELLAR_BARYON_FRAC_MINI;
    *F_ESC10m = ESC_FRAC_MINI;
    *X_LUMINOSITYm = 0;
  }
  else{ return 0;} // format is not allowed
#else //INHOMO_FEEDBACK
#ifdef REION_SM
  if (argc == (*arg_offset + min_argc+6)){
    *F_STAR10 = atof(argv[*arg_offset + min_argc]);
    *ALPHA_STAR = atof(argv[*arg_offset + min_argc+1]);
    *F_ESC10 = atof(argv[*arg_offset + min_argc+2]);
    *ALPHA_ESC = atof(argv[*arg_offset + min_argc+3]);
    *T_AST = t_STAR;
    *X_LUMINOSITY = 0;
    *F_STAR10m = atof(argv[*arg_offset + min_argc+4]);
    *F_ESC10m = atof(argv[*arg_offset + min_argc+5]);
    *X_LUMINOSITYm = 0;
  }
  else if (argc == (*arg_offset + min_argc)){ //These parameters give the result which is the same with the default model.
    *F_STAR10 = STELLAR_BARYON_FRAC;        
    *ALPHA_STAR = STELLAR_BARYON_PL;
    *F_ESC10 = ESC_FRAC;
    *ALPHA_ESC = ESC_PL;
    *T_AST = t_STAR;
    *X_LUMINOSITY = 0;
    *F_STAR10m = STELLAR_BARYON_FRAC_MINI;
    *F_ESC10m = ESC_FRAC_MINI;
    *X_LUMINOSITYm = 0;
  }
  else{ return 0;} // format is not allowed
#else //REION_SM
  if (argc == (*arg_offset + min_argc+7)){
    *F_STAR10 = atof(argv[*arg_offset + min_argc]);
    *ALPHA_STAR = atof(argv[*arg_offset + min_argc+1]);
    *F_ESC10 = atof(argv[*arg_offset + min_argc+2]);
    *ALPHA_ESC = atof(argv[*arg_offset + min_argc+3]);
    *M_TURN = atof(argv[*arg_offset + min_argc+4]); // Input value M_TURN
    *T_AST = t_STAR;
    *X_LUMINOSITY = 0;
    *F_STAR10m = atof(argv[*arg_offset + min_argc+5]);
    *F_ESC10m = atof(argv[*arg_offset + min_argc+6]);
    *X_LUMINOSITYm = 0;
  }
  else if (argc == (*arg_offset + min_argc)){ //These parameters give the result which is the same with the default model.
    *F_STAR10 = STELLAR_BARYON_FRAC;        
    *ALPHA_STAR = STELLAR_BARYON_PL;
    *F_ESC10 = ESC_FRAC;
    *ALPHA_ESC = ESC_PL;
    *M_TURN = M_TURNOVER;
    *T_AST = t_STAR;
    *X_LUMINOSITY = 0;
    *F_STAR10m = STELLAR_BARYON_FRAC_MINI;
    *F_ESC10m = ESC_FRAC_MINI;
    *X_LUMINOSITYm = 0;
  }
  else{ return 0;} // format is not allowed
#endif //REION_SM
#endif //INHOMO_FEEDBACK
#else //MINI_HALO
  if (argc == (*arg_offset + min_argc+5)){
    *F_STAR10 = atof(argv[*arg_offset + min_argc]);
    *ALPHA_STAR = atof(argv[*arg_offset + min_argc+1]);
    *F_ESC10 = atof(argv[*arg_offset + min_argc+2]);
    *ALPHA_ESC = atof(argv[*arg_offset + min_argc+3]);
    *M_TURN = atof(argv[*arg_offset + min_argc+4]); // Input value M_TURN
    *T_AST = t_STAR;
    *X_LUMINOSITY = 0;
  }
  else if (argc == (*arg_offset + min_argc)){ //These parameters give the result which is the same with the default model.
    *F_STAR10 = STELLAR_BARYON_FRAC;        
    *ALPHA_STAR = STELLAR_BARYON_PL;
    *F_ESC10 = ESC_FRAC;
    *ALPHA_ESC = ESC_PL;
    *M_TURN = M_TURNOVER;
    *T_AST = t_STAR;
    *X_LUMINOSITY = 0;
  }
  else{ return 0;} // format is not allowed
#endif //MINI_HALO
#endif //USE_TS_IN_21CM
  *MFP = R_BUBBLE_MAX;
#endif // SHARP_CUTOFF 

  *REDSHIFT = atof(argv[*arg_offset+1]);
#ifdef INHOMO_RECO
  *PREV_REDSHIFT = atof(argv[*arg_offset+2]);
  if (*PREV_REDSHIFT <= *REDSHIFT ){
    fprintf(stderr, "find_HII_bubbles: previous redshift must be larger than current redshift!!!\nAborting...\n");
    return -1;
  }
#else //INHOMO_RECO
  *PREV_REDSHIFT = *REDSHIFT+0.2; // dummy variable which is not used
#endif //INHOMO_RECO

#ifdef SHARP_CUTOFF
   fprintf(stderr, "SHARP_CUTOFF (ON)\n");
   fprintf(stderr, "find_HII_bubbles: command line parameters are as follows\nnum threads=%i, f_star10=%g, alpha_star=%g, f_esc10=%g, alpha_esc=%g, Mturn=%g, t_star=%g, L_X=%g, z=%g, prev z=%g\n",
      *num_th, *F_STAR10, *ALPHA_STAR, *F_ESC10, *ALPHA_ESC, *M_TURN, *T_AST, *X_LUMINOSITY, *REDSHIFT, *PREV_REDSHIFT);
#else //SHARP_CUTOFF
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
   fprintf(stderr, "SHARP_CUTOFF (OFF) MINI_HALO (ON) INHOMO_FEEDBACK (ON)\n");
   fprintf(stderr, "find_HII_bubbles: command line parameters are as follows\nnum threads=%i, f_star10=%g, alpha_star=%g, f_esc10=%g, alpha_esc=%g, t_star=%g, L_X=%g, f_star10m=%g, f_esc10m=%g, L_Xm=%g, z=%g, prev z=%g\n",
      *num_th, *F_STAR10, *ALPHA_STAR, *F_ESC10, *ALPHA_ESC, *T_AST, *X_LUMINOSITY, *F_STAR10m, *F_ESC10m, *X_LUMINOSITYm, *REDSHIFT, *PREV_REDSHIFT);
#else //INHOMO_FEEDBACK
#ifdef REION_SM
   fprintf(stderr, "SHARP_CUTOFF (OFF) MINI_HALO (ON) INHOMO_FEEDBACK (OFF) REION_SM (ON)\n");
   fprintf(stderr, "find_HII_bubbles: command line parameters are as follows\nnum threads=%i, f_star10=%g, alpha_star=%g, f_esc10=%g, alpha_esc=%g, t_star=%g, L_X=%g, f_star10m=%g, f_esc10m=%g, L_Xm=%g, z=%g, prev z=%g\n",
      *num_th, *F_STAR10, *ALPHA_STAR, *F_ESC10, *ALPHA_ESC, *T_AST, *X_LUMINOSITY, *F_STAR10m, *F_ESC10m, *X_LUMINOSITYm, *REDSHIFT, *PREV_REDSHIFT);
#else //REION_SM
   fprintf(stderr, "SHARP_CUTOFF (OFF) MINI_HALO (ON) INHOMO_FEEDBACK (OFF) REION_SM (OFF)\n");
   fprintf(stderr, "find_HII_bubbles: command line parameters are as follows\nnum threads=%i, f_star10=%g, alpha_star=%g, f_esc10=%g, alpha_esc=%g, Mturn=%g, t_star=%g, L_X=%g, f_star10m=%g, f_esc10m=%g, L_Xm=%g, z=%g, prev z=%g\n",
      *num_th, *F_STAR10, *ALPHA_STAR, *F_ESC10, *ALPHA_ESC, *M_TURN, *T_AST, *X_LUMINOSITY, *F_STAR10m, *F_ESC10m, *X_LUMINOSITYm, *REDSHIFT, *PREV_REDSHIFT);
#endif //REION_SM
#endif //INHOMO_FEEDBACK
#else //MINI_HALO
   fprintf(stderr, "SHARP_CUTOFF (OFF) MINI_HALO (OFF)\n");
   fprintf(stderr, "find_HII_bubbles: command line parameters are as follows\nnum threads=%i, f_star10=%g, alpha_star=%g, f_esc10=%g, alpha_esc=%g, Mturn=%g, t_star=%g, L_X=%g, z=%g, prev z=%g\n",
      *num_th, *F_STAR10, *ALPHA_STAR, *F_ESC10, *ALPHA_ESC, *M_TURN, *T_AST, *X_LUMINOSITY, *REDSHIFT, *PREV_REDSHIFT);
#endif //MINI_HALO
#endif //SHARP_CUTOFF

  return 1;
}



/********** MAIN PROGRAM **********/
int main(int argc, char ** argv){
  char filename[1000], error_message[1000];
  FILE *F = NULL, *pPipe = NULL;
  float REDSHIFT, PREV_REDSHIFT, mass, R, xf, yf, zf, growth_factor, pixel_mass, cell_length_factor, massofscaleR;
  float ave_M_coll_cell, ave_N_min_cell, ION_EFF_FACTOR, M_MIN;
  int x,y,z, N_halos_in_cell, N_min_cell, LAST_FILTER_STEP, num_th, arg_offset, i,j,k;
  unsigned long long ct, ion_ct, sample_ct;
  float f_coll_crit, pixel_volume,  density_over_mean, erfc_num, erfc_denom, erfc_denom_cell, res_xH, Splined_Fcoll;
  float *xH=NULL, TVIR_MIN, MFP, xHI_from_xrays, std_xrays, *z_re=NULL, *Gamma12=NULL, *mfp=NULL;
  fftwf_complex *M_coll_unfiltered=NULL, *M_coll_filtered=NULL, *deltax_unfiltered=NULL, *deltax_filtered=NULL, *xe_unfiltered=NULL, *xe_filtered=NULL;
  fftwf_complex *N_rec_unfiltered=NULL, *N_rec_filtered=NULL;
  fftwf_plan plan;
  double global_xH=0, ave_xHI_xrays, ave_den, ST_over_PS, mean_f_coll_st, f_coll, ave_fcoll, dNrec;
  const gsl_rng_type * T=NULL;
  gsl_rng * r=NULL;
  double t_ast, dfcolldt, Gamma_R_prefactor, rec;
  float nua, dnua, temparg, Gamma_R, z_eff;
  float F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, Mlim_Fstar, Mlim_Fesc; //New in v2
#ifdef MINI_HALO
  float F_STAR10m, F_ESC10m, Mlim_Fstarm, ION_EFF_FACTOR_MINI,M_MINm, M_MINa, Splined_Fcollm, dfcolldtm,Gamma_R_prefactorm,ST_over_PSm,f_collm, Mcrit_atom, Mcrit_LW, Mcrit_RE, mean_f_coll_stm; //New in v2.1
  double X_LUMINOSITYm;
#else //MINI_HALO
  float M_MINa;
#endif //MINI_HALO
  double X_LUMINOSITY;
  float fabs_dtdz, ZSTEP;
  const float dz = 0.01;
  *error_message = '\0';

  double aveR = 0;
  unsigned long long Rct = 0;

    
  /*************************************************************************************/  
  /******** BEGIN INITIALIZATION ********/
  /*************************************************************************************/  

  // PARSE COMMAND LINE ARGUMENTS
#ifdef SHARP_CUTOFF
  if( !parse_arguments(argc, argv, &num_th, &arg_offset, &F_STAR10, &ALPHA_STAR, &F_ESC10,
             &ALPHA_ESC, &M_TURN, &T_AST, &X_LUMINOSITY, &MFP, &REDSHIFT, &PREV_REDSHIFT))
  {
      fprintf(stderr, "find_HII_bubbles <redshift> [<previous redshift>] \n \
          Aborting...\n                               \
      Check that your inclusion (or not) of [<previous redshift>] is consistent with the INHOMO_RECO flag in ../Parameter_files/ANAL_PARAMS.H\nAborting...\n");
      return -1;
  }
#else //SHARP_CUTOFF
#ifdef MINI_HALO
  Mcrit_atom = atomic_cooling_threshold(REDSHIFT);
#ifdef INHOMO_FEEDBACK
  Mcrit_LW   = lyman_werner_threshold(REDSHIFT, 0); // no LW suppression, for checking Dark Age
  if( !parse_arguments(argc, argv, &num_th, &arg_offset, &F_STAR10, &ALPHA_STAR, &F_ESC10,
             &ALPHA_ESC, &T_AST, &X_LUMINOSITY, &F_STAR10m, &F_ESC10m, &X_LUMINOSITYm, &MFP, &REDSHIFT, &PREV_REDSHIFT)){
      fprintf(stderr, "find_HII_bubbles <redshift> [<previous redshift>] \n \
      additional optional arguments: <f_star10> <alpha,star> <f_esc10> <alpha,esc> [<t_star>] [<f_star10,m> <f_esc10,m>]\n \
      Check that your inclusion (or not) of [<previous redshift>] is consistent with the INHOMO_RECO flag in ../Parameter_files/ANAL_PARAMS.H\n \
   Also check that your inclusion (or not) of [<t_star>] is consistent with the USE_TS_IN_21CM flag in ../Parameter_files/HEAT_PARAMS.H\nAborting...\n");
      return -1;
  }
#else //INHOMO_FEEDBACK
  Mcrit_LW   = lyman_werner_threshold(REDSHIFT);
#ifdef REION_SM
  if( !parse_arguments(argc, argv, &num_th, &arg_offset, &F_STAR10, &ALPHA_STAR, &F_ESC10,
             &ALPHA_ESC, &T_AST, &X_LUMINOSITY, &F_STAR10m, &F_ESC10m, &X_LUMINOSITYm, &MFP, &REDSHIFT, &PREV_REDSHIFT)){
      fprintf(stderr, "find_HII_bubbles <redshift> [<previous redshift>] \n \
      additional optional arguments: <f_star10> <alpha,star> <f_esc10> <alpha,esc> [<t_star>] [<f_star10,m> <f_esc10,m>]\n \
      Check that your inclusion (or not) of [<previous redshift>] is consistent with the INHOMO_RECO flag in ../Parameter_files/ANAL_PARAMS.H\n \
   Also check that your inclusion (or not) of [<t_star>] is consistent with the USE_TS_IN_21CM flag in ../Parameter_files/HEAT_PARAMS.H\nAborting...\n");
      return -1;
  }
  double REION_SM13_Z_RE, REION_SM13_DELTA_Z_RE, REION_SM13_DELTA_Z_SC;
  reading_reionization_SM13parameters(&REION_SM13_Z_RE, &REION_SM13_DELTA_Z_RE, &REION_SM13_DELTA_Z_SC);
  Mcrit_RE   = reionization_feedback(REDSHIFT, REION_SM13_Z_RE, REION_SM13_DELTA_Z_RE, REION_SM13_DELTA_Z_SC);
#else //REION_SM
  if( !parse_arguments(argc, argv, &num_th, &arg_offset, &F_STAR10, &ALPHA_STAR, &F_ESC10,
             &ALPHA_ESC, &M_TURN, &T_AST, &X_LUMINOSITY, &F_STAR10m, &F_ESC10m, &X_LUMINOSITYm, &MFP, &REDSHIFT, &PREV_REDSHIFT)){
      fprintf(stderr, "find_HII_bubbles <redshift> [<previous redshift>] \n \
      additional optional arguments: <f_star10> <alpha,star> <f_esc10> <alpha,esc> <M_TURNOVER>] [<t_star>] [<f_star10,m> <f_esc10,m>]\n \
      Check that your inclusion (or not) of [<previous redshift>] is consistent with the INHOMO_RECO flag in ../Parameter_files/ANAL_PARAMS.H\n \
   Also check that your inclusion (or not) of [<t_star>] is consistent with the USE_TS_IN_21CM flag in ../Parameter_files/HEAT_PARAMS.H\nAborting...\n");
      return -1;
  }
  Mcrit_RE   = M_TURN;
#endif //REION_SM
#endif //INHOMO_FEEDBACK
#else //MINI_HALO
  if( !parse_arguments(argc, argv, &num_th, &arg_offset, &F_STAR10, &ALPHA_STAR, &F_ESC10,
             &ALPHA_ESC, &M_TURN, &T_AST, &X_LUMINOSITY, &MFP, &REDSHIFT, &PREV_REDSHIFT)){
      fprintf(stderr, "find_HII_bubbles <redshift> [<previous redshift>] \n \
      additional optional arguments: <f_star10> <alpha,star> <f_esc10> <alpha,esc> <M_TURNOVER>] [<t_star>]\n \
      Check that your inclusion (or not) of [<previous redshift>] is consistent with the INHOMO_RECO flag in ../Parameter_files/ANAL_PARAMS.H\n \
   Also check that your inclusion (or not) of [<t_star>] is consistent with the USE_TS_IN_21CM flag in ../Parameter_files/HEAT_PARAMS.H\nAborting...\n");
      return -1;
  }
#endif //MINI_HALO
#endif //SHARP_CUTOFF

  ZSTEP = PREV_REDSHIFT - REDSHIFT;
  fabs_dtdz = fabs(dtdz(REDSHIFT));
  t_ast = T_AST * t_hubble(REDSHIFT);
  growth_factor = dicke(REDSHIFT);
  pixel_volume = pow(BOX_LEN/(float)HII_DIM, 3);
  pixel_mass = RtoM(L_FACTOR*BOX_LEN/(float)HII_DIM);
  // this parameter choice is sensitive to noise on the cell size, at least for the typical
  // cell sizes in RT simulations.  it probably doesn't matter for larger cell sizes.
  cell_length_factor = L_FACTOR;

#ifdef USE_HALO_FIELD
  if ((FIND_BUBBLE_ALGORITHM==2)&& ((BOX_LEN/(float)HII_DIM) < 1)) // fairly arbitrary length based on 2 runs i did
    cell_length_factor = 1;
#endif //USE_HALO_FIELD
  init_ps();
  Fcoll = (float *) malloc(sizeof(float)*HII_TOT_FFT_NUM_PIXELS);
#ifdef MINI_HALO
  Fcollm = (float *) malloc(sizeof(float)*HII_TOT_FFT_NUM_PIXELS);
#endif //MINI_HALO
#ifdef INHOMO_RECO 
  init_MHR();
#endif //INHOMO_RECO 
#ifndef SHARP_CUTOFF
    // Set the minimum halo mass hosting ionizing source mass.
    // For constant ionizing efficiency parameter M_MIN is set to be M_TURN which is a sharp cut-off.
    // For the new parametrization the number of halos hosting active galaxies (i.e. the duty cycle) is assumed to
    // exponentially decrease below M_TURNOVER Msun, : fduty \propto e^(- M_TURNOVER / M)
    // In this case, we define M_MIN = M_TURN/50.
    // v2.1 split halos in to atomic and molecular populations, with different fduty see ANAL_PARAMS.H
  init_21cmMC_arrays();
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
  M_MINa = Mcrit_atom;
  M_MINm = Mcrit_LW; // basically it's molecular cooling threshold
  M_MIN  = 1e5;
#else //INHOMO_FEEDBACK
  //TODO: here is inconsistent with the M_MIN set in Ts, which is the minimum of all scrolled redshifts.
  //Shouldn't matter too much, but better change it later.
  M_MINa = Mcrit_RE > Mcrit_atom ? Mcrit_RE : Mcrit_atom;
  M_MINm = Mcrit_RE > Mcrit_LW   ? Mcrit_RE : Mcrit_LW;
  M_MIN  = M_MINa > M_MINm       ? M_MINm : M_MINa;
  M_MIN /= 50;
#endif //INHOMO_FEEDBACK
  Mlim_Fstarm         = Mass_limit_bisection(M_MIN, 1e16, ALPHA_STAR, F_STAR10m);
  ION_EFF_FACTOR_MINI = N_GAMMA_UV_MINI * F_STAR10m * F_ESC10m;
#else //MINI_HALO
  M_MINa = M_TURN;
  M_MIN  = M_TURN/50;
#endif //MINI_HALO
  Mlim_Fstar     = Mass_limit_bisection(M_MIN, 1e16, ALPHA_STAR, F_STAR10);
  Mlim_Fesc      = Mass_limit_bisection(M_MIN, 1e16, ALPHA_ESC, F_ESC10);
  ION_EFF_FACTOR = N_GAMMA_UV      * F_STAR10  * F_ESC10;
#else //SHARP_CUTOFF
  M_MIN          = M_TURNOVER;
  ION_EFF_FACTOR = N_GAMMA_UV * STELLAR_BARYON_FRAC * ESC_FRAC; // Constant ionizing efficiency parameter.
#endif //SHARP_CUTOFF

  // check for WDM
#ifdef P_CUTOFF 
  if ( M_MIN < M_J_WDM()){
    fprintf(stderr, "The default Jeans mass of %e Msun is smaller than the scale supressed by the effective pressure of WDM.\n", M_MIN);
    M_MIN = M_J_WDM();
    fprintf(stderr, "Setting a new effective Jeans mass from WDM pressure supression of %e Msun\n", M_MIN);
  }
#endif //P_CUTOFF
  initialiseSplinedSigmaM(M_MIN,5e16); // set up interpolation table for computing sigma_M
  
  // INITIALIZE THREADS
  if (fftwf_init_threads()==0){
    fprintf(stderr, "find_HII_bubbles: ERROR: problemin itializing fftwf threads\nAborting\n.");
    return -1;
  }
  omp_set_num_threads(num_th);
  fftwf_plan_with_nthreads(num_th);
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  
  // OPEN LOG FILE
  system("mkdir -p ../Log_files");
  sprintf(filename, "../Log_files/HII_bubble_log_file_%d", getpid());
  LOG = fopen(filename, "w");
  if (!LOG){
    fprintf(stderr, "find_HII_bubbles.c: Error opening log file\n");
  }

  // allocate memory for the neutral fraction box
  xH = (float *) fftwf_malloc(sizeof(float)*HII_TOT_NUM_PIXELS);
  if (!xH){
    strcpy(error_message, "find_HII_bubbles.c: Error allocating memory for xH box\nAborting...\n");
    goto CLEANUP;
  }
  for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){    xH[ct] = 1;  }

  // compute the mean collpased fraction at this redshift
#ifdef SHARP_CUTOFF
  mean_f_coll_st = FgtrM_st(REDSHIFT, M_MIN);
#else //SHARP_CUTOFF
  // Here 'mean_f_coll_st' is not the mean collpased fraction, but leave this name as is to simplify the variable name.
  // Nion_ST * ION_EFF_FACTOR = the mean number of IGM ionizing photon per baryon
  // see eq. (17) in Park et al. 2018
  //
  // Note from YQ, ifdef IHHOMO_FEEDBACK, this calculate the mean collpased fraction when there is no LW or reionization
  mean_f_coll_st = Nion_ST(REDSHIFT, M_MIN, M_MINa, ALPHA_STAR, ALPHA_ESC, F_STAR10, F_ESC10, Mlim_Fstar, Mlim_Fesc);
#ifdef MINI_HALO
  // Nion_ST * ION_EFF_FACTOR + Nion_STm * ION_EFF_FACTOR_MINI = the mean number of IGM ionizing photon per baryon
  mean_f_coll_stm = Nion_STm(REDSHIFT, M_MIN, M_MINm, Mcrit_atom, ALPHA_STAR, F_STAR10m, Mlim_Fstarm);
#endif //MINI_HALO
#endif //SHARP_CUTOFF

  /**********  CHECK IF WE ARE IN THE DARK AGES ******************************/
  // lets check if we are going to bother with computing the inhmogeneous field at all...
#ifdef MINI_HALO
  if ((mean_f_coll_st*ION_EFF_FACTOR + mean_f_coll_stm*ION_EFF_FACTOR_MINI < HII_ROUND_ERR)) // way too small to ionize anything...//New in v2.1
#else //MINI_HALO
  if ((mean_f_coll_st*ION_EFF_FACTOR < HII_ROUND_ERR)) // way too small to ionize anything...//New in v2
#endif //MINI_HALO
  {
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
    fprintf(stderr, "INHOMO_FEEDBACK is ON!\n \
                     For the purpose of checking if we are in the dark ages, the calculation here assumes no LW or RE background\n \
                     The ST mean collapse fraction is %e (atomic) and %e (molecular),\n \
                     which is much smaller than the effective critical collapse fraction of %e (atomic) and %e (molecular)\n \
                     I will just declare everything to be neutral\n", mean_f_coll_st, mean_f_coll_stm, 1./ION_EFF_FACTOR, 1./ION_EFF_FACTOR_MINI);
    fprintf(LOG, "INHOMO_FEEDBACK is ON!\n \
                  For the purpose of checking if we are in the dark ages, the calculation here assumes no LW or RE background\n \
                  The ST mean collapse fraction is %e (atomic) and %e (molecular),\n \
                  which is much smaller than the effective critical collapse fraction of %e (atomic) and %e (molecular)\n \
                  I will just declare everything to be neutral\n", mean_f_coll_st, mean_f_coll_stm, 1./ION_EFF_FACTOR, 1./ION_EFF_FACTOR_MINI);
#else //INHOMO_FEEDBACK
    fprintf(stderr, "The ST mean collapse fraction is %e (atomic) and %e (molecular),\n \
                     which is much smaller than the effective critical collapse fraction of %e (atomic) and %e (molecular)\n \
                     I will just declare everything to be neutral\n", mean_f_coll_st, mean_f_coll_stm, 1./ION_EFF_FACTOR, 1./ION_EFF_FACTOR_MINI);
    fprintf(LOG, "The ST mean collapse fraction is %e (atomic) and %e (molecular),\n \
                  which is much smaller than the effective critical collapse fraction of %e (atomic) and %e (molecular)\n \
                  I will just declare everything to be neutral\n", mean_f_coll_st, mean_f_coll_stm, 1./ION_EFF_FACTOR, 1./ION_EFF_FACTOR_MINI);
#endif //INHOMO_FEEDBACK
#else //MINI_HALO
    fprintf(stderr, "The ST mean collapse fraction is %e, which is much smaller than the effective critical collapse fraction of %e\n \
                     I will just declare everything to be neutral\n", mean_f_coll_st, 1./ION_EFF_FACTOR);
    fprintf(LOG, "The ST mean collapse fraction is %e, which is much smaller than the effective critical collapse fraction of %e\n \
                  I will just declare everything to be neutral\n", mean_f_coll_st, 1./ION_EFF_FACTOR);
#endif //MINI_HALO

#ifdef USE_TS_IN_21CM
#ifndef SHARP_CUTOFF
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
    sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_L_X%.1e_alphaX%.1f_f_star10_%06.4f_alpha_star%06.4f_f_esc10_%06.4f_alpha_esc%06.4f_t_star%06.4f_f_star10m%06.4f_f_esc10m%06.4f_L_Xm%.1e_alphaXm%.1f_%i_%.0fMpc", REDSHIFT, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, T_AST, F_STAR10m, F_ESC10m, X_LUMINOSITYm, X_RAY_SPEC_INDEX_MINI, HII_DIM, BOX_LEN); 
    fprintf(stderr, "filename: %s\n",filename);
#else //INHOMO_FEEDBACK
#ifdef REION_SM
    sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_L_X%.1e_alphaX%.1f_f_star10_%06.4f_alpha_star%06.4f_f_esc10_%06.4f_alpha_esc%06.4f_Msm_t_star%06.4f_f_star10m%06.4f_f_esc10m%06.4f_L_Xm%.1e_alphaXm%.1f_%i_%.0fMpc", REDSHIFT, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, T_AST, F_STAR10m, F_ESC10m, X_LUMINOSITYm, X_RAY_SPEC_INDEX_MINI, HII_DIM, BOX_LEN); 
    fprintf(stderr, "filename: %s\n",filename);
#else //REION_SM
    sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_L_X%.1e_alphaX%.1f_f_star10_%06.4f_alpha_star%06.4f_f_esc10_%06.4f_alpha_esc%06.4f_Mturn%.1e_t_star%06.4f_f_star10m%06.4f_f_esc10m%06.4f_L_Xm%.1e_alphaXm%.1f_%i_%.0fMpc", REDSHIFT, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, T_AST, F_STAR10m, F_ESC10m, X_LUMINOSITYm, X_RAY_SPEC_INDEX_MINI, HII_DIM, BOX_LEN); 
    fprintf(stderr, "filename: %s\n",filename);
#endif //REION_SM
#endif //INHOMO_FEEDBACK
#else //MINI_HALO
    sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_L_X%.1e_alphaX%.1f_f_star10_%06.4f_alpha_star%06.4f_f_esc10_%06.4f_alpha_esc%06.4f_Mturn%.1e_t_star%06.4f_Pop%i_%i_%.0fMpc", REDSHIFT, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, T_AST, Pop, HII_DIM, BOX_LEN); 
    fprintf(stderr, "filename: %s\n",filename);
#endif //MINI_HALO
#else //SHARP_CUTOFF
    sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_L_X%.1e_alphaX%.1f_Mmin%.1e_zetaIon%.2f_Pop%i_%i_%.0fMpc", REDSHIFT, X_LUMINOSITY, X_RAY_SPEC_INDEX, M_MIN, HII_EFF_FACTOR, Pop, HII_DIM, BOX_LEN); 
#endif //SHARP_CUTOFF
    if (!(F = fopen(filename, "rb"))){
      fprintf(stderr, "find_HII_bubbles: Unable to open x_e file at %s\nAborting...\n", filename);
      fprintf(LOG, "find_HII_bubbles: Unable to open x_e file at %s\nAborting...\n", filename);
      fclose(LOG); fftwf_free(xH); fftwf_cleanup_threads();
      free(Fcoll);
#ifdef MINI_HALO
      free(Fcollm);
#endif
      free_ps();  
#ifndef SHARP_CUTOFF
      destroy_21cmMC_arrays();
#endif
      return -1;
    }
    for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){
      if (fread(&xH[ct], sizeof(float), 1, F)!=1){
        strcpy(error_message, "find_HII_bubbles.c: Read error occured while reading xe box.\nAborting...\n");
        goto CLEANUP;
      }
      xH[ct] = 1-xH[ct]; // convert from x_e to xH
      if (xH[ct]<0) xH[ct] = 0; //  should not happen....
      global_xH += xH[ct];
    }
    fclose(F);
    F = NULL;
    global_xH /= (double)HII_TOT_NUM_PIXELS;
#else //USE_TS_IN_21CM
    // find the neutral fraction
    init_heat();
    global_xH = 1 - xion_RECFAST(REDSHIFT, 0);;
    destruct_heat();
    for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){
      xH[ct] = global_xH;
    }
#endif //USE_TS_IN_21CM
  
    // print out the xH box
    global_xH_m = global_xH;
    switch(FIND_BUBBLE_ALGORITHM){
      case 2:
#ifndef SHARP_CUTOFF
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
#ifdef USE_HALO_FIELD
        sprintf(filename, "../Boxes/xH_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_f_star10m%.4f_f_esc10m%.4f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, F_STAR10m, F_ESC10m, HII_FILTER, MFP, HII_DIM, BOX_LEN);
#else //USE_HALO_FIELD
        sprintf(filename, "../Boxes/xH_nohalos_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_f_star10m%.4f_f_esc10m%.4f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, F_STAR10m, F_ESC10m, HII_FILTER, MFP, HII_DIM, BOX_LEN);
#endif //USE_HALO_FIELD
#else //INHOMO_FEEDBACK
#ifdef REION_SM
#ifdef USE_HALO_FIELD
        sprintf(filename, "../Boxes/xH_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_Msm_f_star10m%.4f_f_esc10m%.4f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, F_STAR10m, F_ESC10m, HII_FILTER, MFP, HII_DIM, BOX_LEN);
#else //USE_HALO_FIELD
        sprintf(filename, "../Boxes/xH_nohalos_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_Msm_f_star10m%.4f_f_esc10m%.4f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, F_STAR10m, F_ESC10m, HII_FILTER, MFP, HII_DIM, BOX_LEN);
#endif //USE_HALO_FIELD
#else //REION_SM
#ifdef USE_HALO_FIELD
        sprintf(filename, "../Boxes/xH_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_Mturn%.2e_f_star10m%.4f_f_esc10m%.4f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, F_STAR10m, F_ESC10m, HII_FILTER, MFP, HII_DIM, BOX_LEN);
#else //USE_HALO_FIELD
        sprintf(filename, "../Boxes/xH_nohalos_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_Mturn%.2e_f_star10m%.4f_f_esc10m%.4f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, F_STAR10m, F_ESC10m, HII_FILTER, MFP, HII_DIM, BOX_LEN);
#endif //USE_HALO_FIELD
#endif //REION_SM
#endif //INHOMO_FEEDBACK
#else //MINI_HALO
#ifdef USE_HALO_FIELD
        sprintf(filename, "../Boxes/xH_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_Mturn%.2e_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, HII_FILTER, MFP, HII_DIM, BOX_LEN);
#else //USE_HALO_FIELD
        sprintf(filename, "../Boxes/xH_nohalos_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_Mturn%.2e_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, HII_FILTER, MFP, HII_DIM, BOX_LEN);
#endif //USE_HALO_FIELD
#endif //MINI_HALO
#else //SHARP_CUTOFF
#ifdef USE_HALO_FIELD
        sprintf(filename, "../Boxes/xH_z%06.2f_nf%f_eff%.1f_effPLindex0_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, ION_EFF_FACTOR, HII_FILTER, M_MIN, MFP, HII_DIM, BOX_LEN);
#else //USE_HALO_FIELD
        sprintf(filename, "../Boxes/xH_nohalos_z%06.2f_nf%f_eff%.1f_effPLindex0_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, ION_EFF_FACTOR, HII_FILTER, M_MIN, MFP, HII_DIM, BOX_LEN);
#endif //USE_HALO_FIELD
#endif //SHARP_CUTOFF
        break;
      default:
#ifndef SHARP_CUTOFF
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
#ifdef USE_HALO_FIELD
        sprintf(filename, "../Boxes/sphere_xH_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_f_star10m%.4f_f_esc10m%.4f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, F_STAR10m, F_ESC10m, HII_FILTER, MFP, HII_DIM, BOX_LEN);
#else //USE_HALO_FIELD
        sprintf(filename, "../Boxes/sphere_xH_nohalos_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_f_star10m%.4f_f_esc10m%.4f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, F_STAR10m, F_ESC10m, HII_FILTER, MFP, HII_DIM, BOX_LEN);
#endif //USE_HALO_FIELD
#else //INHOMO_FEEDBACK
#ifdef REION_SM
#ifdef USE_HALO_FIELD
        sprintf(filename, "../Boxes/sphere_xH_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_Msm_f_star10m%.4f_f_esc10m%.4f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, F_STAR10m, F_ESC10m, HII_FILTER, MFP, HII_DIM, BOX_LEN);
#else //USE_HALO_FIELD
        sprintf(filename, "../Boxes/sphere_xH_nohalos_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_Msm_f_star10m%.4f_f_esc10m%.4f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, F_STAR10m, F_ESC10m, HII_FILTER, MFP, HII_DIM, BOX_LEN);
#endif //USE_HALO_FIELD
#else //REION_SM
#ifdef USE_HALO_FIELD
        sprintf(filename, "../Boxes/sphere_xH_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_Mturn%.2e_f_star10m%.4f_f_esc10m%.4f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, F_STAR10m, F_ESC10m, HII_FILTER, MFP, HII_DIM, BOX_LEN);
#else //USE_HALO_FIELD
        sprintf(filename, "../Boxes/sphere_xH_nohalos_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_Mturn%.2e_f_star10m%.4f_f_esc10m%.4f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, F_STAR10m, F_ESC10m, HII_FILTER, MFP, HII_DIM, BOX_LEN);
#endif //USE_HALO_FIELD
#endif //REION_SM
#endif //INHOMO_FEEDBACK
#else //MINI_HALO
#ifdef USE_HALO_FIELD
        sprintf(filename, "../Boxes/sphere_xH_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_Mturn%.2e_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, HII_FILTER, MFP, HII_DIM, BOX_LEN);
#else //USE_HALO_FIELD
        sprintf(filename, "../Boxes/sphere_xH_nohalos_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_Mturn%.2e_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, HII_FILTER, MFP, HII_DIM, BOX_LEN);
#endif //USE_HALO_FIELD
#endif //MINI_HALO
#else //SHARP_CUTOFF
#ifdef USE_HALO_FIELD
        sprintf(filename, "../Boxes/sphere_xH_z%06.2f_nf%f_eff%.1f_effPLindex0_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, ION_EFF_FACTOR, HII_FILTER, M_MIN, MFP, HII_DIM, BOX_LEN);
#else //USE_HALO_FIELD
        sprintf(filename, "../Boxes/sphere_xH_nohalos_z%06.2f_nf%f_eff%.1f_effPLindex0_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, ION_EFF_FACTOR, HII_FILTER, M_MIN, MFP, HII_DIM, BOX_LEN);
#endif //USE_HALO_FIELD
#endif //SHARP_CUTOFF
    }
    F = fopen(filename, "wb");
    fprintf(LOG, "Neutral fraction is %f\nNow writting xH box at %s\n", global_xH, filename);
    fprintf(stderr, "Neutral fraction is %f\nNow writting xH box at %s\n", global_xH, filename);
    if (mod_fwrite(xH, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
      fprintf(stderr, "find_HII_bubbles.c: Write error occured while writting xH box.\n");
      fprintf(LOG, "find_HII_bubbles.c: Write error occured while writting xH box.\n");
    }
    free_ps(); fclose(F); F = NULL; fclose(LOG); fftwf_free(xH); fftwf_cleanup_threads();
    free(Fcoll);
#ifdef MINI_HALO
    free(Fcollm);
#endif
#ifndef SHARP_CUTOFF
     destroy_21cmMC_arrays();
#endif
     return (int) (global_xH * 100);
  }
  /*************   END CHECK TO SEE IF WE ARE STILL IN THE DARK AGES *************/

  // ARE WE INCLUDING PRE-IONIZATION FROM X-RAYS??
#ifdef USE_TS_IN_21CM
  xe_unfiltered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
  xe_filtered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
  if (!xe_filtered || !xe_unfiltered){
    strcpy(error_message, "find_HII_bubbles.c: Error allocating memory for xe boxes\nAborting...\n");
    goto CLEANUP;
  }

  // and read-in
#ifndef SHARP_CUTOFF
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
  sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_L_X%.1e_alphaX%.1f_f_star10_%06.4f_alpha_star%06.4f_f_esc10_%06.4f_alpha_esc%06.4f_t_star%06.4f_f_star10m%.4f_f_esc10m%.4f_L_Xm%.1e_alphaXm%.1f_%i_%.0fMpc", REDSHIFT, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, T_AST, F_STAR10m, F_ESC10m, X_LUMINOSITYm, X_RAY_SPEC_INDEX_MINI, HII_DIM, BOX_LEN);
#else //INHOMO_FEEDBACK
#ifdef REION_SM
  sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_L_X%.1e_alphaX%.1f_f_star10_%06.4f_alpha_star%06.4f_f_esc10_%06.4f_alpha_esc%06.4f_Msm_t_star%06.4f_f_star10m%.4f_f_esc10m%.4f_L_Xm%.1e_alphaXm%.1f_%i_%.0fMpc", REDSHIFT, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, T_AST, F_STAR10m, F_ESC10m, X_LUMINOSITYm, X_RAY_SPEC_INDEX_MINI, HII_DIM, BOX_LEN);
#else //REION_SM
  sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_L_X%.1e_alphaX%.1f_f_star10_%06.4f_alpha_star%06.4f_f_esc10_%06.4f_alpha_esc%06.4f_Mturn%.1e_t_star%06.4f_f_star10m%.4f_f_esc10m%.4f_L_Xm%.1e_alphaXm%.1f_%i_%.0fMpc", REDSHIFT, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, T_AST, F_STAR10m, F_ESC10m, X_LUMINOSITYm, X_RAY_SPEC_INDEX_MINI, HII_DIM, BOX_LEN);
#endif //REION_SM
#endif //INHOMO_FEEDBACK
#else //MINI_HALO
  sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_L_X%.1e_alphaX%.1f_f_star10_%06.4f_alpha_star%06.4f_f_esc10_%06.4f_alpha_esc%06.4f_Mturn%.1e_t_star%06.4f_Pop%i_%i_%.0fMpc", REDSHIFT, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, T_AST, Pop, HII_DIM, BOX_LEN);
#endif //MINI_HALO
#else //SHARP_CUTOFF
  sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_L_X%.1e_alphaX%.1f_Mmin%.1e_zetaIon%.2f_Pop%i_%i_%.0fMpc", REDSHIFT, X_LUMINOSITY, X_RAY_SPEC_INDEX, M_MIN, HII_EFF_FACTOR, Pop, HII_DIM, BOX_LEN); 
#endif //SHARP_CUTOFF
  if (!(F = fopen(filename, "rb"))){
    strcpy(error_message, "find_HII_bubbles.c: Unable to open x_e file at ");
    strcat(error_message, filename);
    strcat(error_message, "\nAborting...\n");
    goto CLEANUP;
  }
  for (i=0; i<HII_DIM; i++){
    for (j=0; j<HII_DIM; j++){
      for (k=0; k<HII_DIM; k++){
        if (fread((float *)xe_unfiltered + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, F)!=1){
          strcpy(error_message, "find_HII_bubbles.c: Read error occured while reading xe box.\nAborting...\n");
          goto CLEANUP;
        }
      }
    }
  }
  fclose(F);
  F = NULL;
#endif //USE_TS_IN_21CM


  // ARE WE USING A DISCRETE HALO FIELD (identified in the ICs with find_halos.c and evolved  with update_halos.c)
#ifdef USE_HALO_FIELD
  // allocate memory for the smoothed halo field
  M_coll_unfiltered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
  M_coll_filtered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
  if (!M_coll_unfiltered || !M_coll_filtered){
    strcpy(error_message, "find_HII_bubbles.c: Error allocating memory for M_coll boxes\nAborting...\n");
    goto CLEANUP;
  }
  for (ct=0; ct<HII_TOT_FFT_NUM_PIXELS; ct++){    *((float *)M_coll_unfiltered + ct) = 0;  }
    
  // read in the halo list
  sprintf(filename, "../Output_files/Halo_lists/updated_halos_z%06.2f_%i_%.0fMpc", REDSHIFT, DIM, BOX_LEN);
  F = fopen(filename, "r");
  if (!F){
    strcpy(error_message, "find_HII_bubbles.c: Unable to open halo list file: ");
    strcat(error_message, filename);
    strcat(error_message, "\nAborting...\n");
    goto CLEANUP;
  }
  // now read in all halos above our threshold into the smoothed halo field
  fscanf(F, "%e %f %f %f", &mass, &xf, &yf, &zf);
  while (!feof(F) && (mass>=M_MIN)){
    x = xf*HII_DIM;
    y = yf*HII_DIM;
    z = zf*HII_DIM;
    *((float *)M_coll_unfiltered + HII_R_FFT_INDEX(x, y, z)) += mass;
    fscanf(F, "%e %f %f %f", &mass, &xf, &yf, &zf);
  }
  fclose(F);
  F = NULL;
#endif //USE_HALO_FIELD

    
  // ALLOCATE AND READ-IN THE EVOLVED DENSITY FIELD
  deltax_unfiltered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
  deltax_filtered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
  if (!deltax_unfiltered || !deltax_filtered){
    strcpy(error_message, "find_HII_bubbles.c: Error allocating memory for deltax boxes\nAborting...\n");
    goto CLEANUP;
  }
  fprintf(stderr, "Reading in deltax box\n");
  fprintf(LOG, "Reading in deltax box\n");

  sprintf(filename, "../Boxes/updated_smoothed_deltax_z%06.2f_%i_%.0fMpc", REDSHIFT, HII_DIM, BOX_LEN);
  F = fopen(filename, "rb");
  if (!F){
    strcpy(error_message, "find_HII_bubbles.c: Unable to open file: ");
    strcat(error_message, filename);
    strcat(error_message, "\nAborting...\n");
    goto CLEANUP;
  }
  for (i=0; i<HII_DIM; i++){
    for (j=0; j<HII_DIM; j++){
      for (k=0; k<HII_DIM; k++){
        if (fread((float *)deltax_unfiltered + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, F)!=1){
          strcpy(error_message, "find_HII_bubbles.c: Read error occured while reading deltax box.\n");
          goto CLEANUP;
        }
      }
    }
  }
  fclose(F);
  F = NULL;

  // ALLOCATE AND INITIALIZE ADDITIONAL BOXES NEEDED TO KEEP TRACK OF RECOMBINATIONS (Sobacchi & Mesinger 2014; NEW IN v1.3)
#ifdef INHOMO_RECO //  flag in ANAL_PARAMS.H to determine whether to compute recombinations or not
  z_re = (float *) fftwf_malloc(sizeof(float)*HII_TOT_NUM_PIXELS); // the redshift at which the cell is ionized
  Gamma12 = (float *) fftwf_malloc(sizeof(float)*HII_TOT_NUM_PIXELS);  // stores the ionizing backgroud
  N_rec_unfiltered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS); // cumulative number of recombinations
  N_rec_filtered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
  if (!z_re || !N_rec_filtered || !N_rec_unfiltered || !Gamma12){
    strcpy(error_message, "find_HII_bubbles.c: Error allocating memory for recombination boxes boxes\nAborting...\n");
    goto CLEANUP;
  }
  sprintf(filename, "../Boxes/z_first_ionization_z%06.2f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc",  PREV_REDSHIFT, HII_FILTER, MFP, HII_DIM, BOX_LEN);
  if (F=fopen(filename, "rb")){  // this is the first call for this run, i.e. at the highest redshift
    //check if some read error occurs
    if (mod_fread(z_re, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
    strcpy(error_message, "find_HII_bubbles.c: Read error occured while reading z_re box!\n");
    goto CLEANUP;
    }
  }
  else{ // this is the first call (highest redshift)
    for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++)
      z_re[ct] = -1.0;
  }
  sprintf(filename, "../Boxes/Nrec_z%06.2f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", PREV_REDSHIFT, HII_FILTER, MFP, HII_DIM, BOX_LEN);
  if (F=fopen(filename, "rb")){ // we had prvious boxes
    //check if some read error occurs
    for (i=0; i<HII_DIM; i++){
      for (j=0; j<HII_DIM; j++){
        for (k=0; k<HII_DIM; k++){
          if (fread((float *)N_rec_unfiltered + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, F)!=1){
            strcpy(error_message, "find_HII_bubbles.c: Read error occured while reading N_rec box!\n");
            goto CLEANUP;
          }
        }
      }
    }
  }
  else{
    fprintf(stderr, "Earliest redshift call, initializing Nrec to zeros\n");
    // initialize N_rec
    for (i=0; i<HII_DIM; i++){
      for (j=0; j<HII_DIM; j++){
        for (k=0; k<HII_DIM; k++){
          *((float *)N_rec_unfiltered + HII_R_FFT_INDEX(i,j,k)) = 0.0;
        }
      }
    }
  }
#endif //INHOMO_RECO

  // do the fft to get the k-space M_coll field and deltax field
  fprintf(LOG, "begin initial ffts, clock=%06.2f\n", (double)clock()/CLOCKS_PER_SEC);
  fflush(LOG);
#ifdef USE_HALO_FIELD
  plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)M_coll_unfiltered, (fftwf_complex *)M_coll_unfiltered, FFTW_ESTIMATE);
  fftwf_execute(plan);
#endif //USE_HALO_FIELD
#ifdef USE_TS_IN_21CM
  plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)xe_unfiltered, (fftwf_complex *)xe_unfiltered, FFTW_ESTIMATE);
  fftwf_execute(plan);
#endif //USE_TS_IN_21CM 
#ifdef INHOMO_RECO
  plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)N_rec_unfiltered, (fftwf_complex *)N_rec_unfiltered, FFTW_ESTIMATE);
  fftwf_execute(plan);
#endif //INHOMO_RECO
  plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)deltax_unfiltered, (fftwf_complex *)deltax_unfiltered, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);
  fftwf_cleanup();
  // remember to add the factor of VOLUME/TOT_NUM_PIXELS when converting from
  //  real space to k-space
  // Note: we will leave off factor of VOLUME, in anticipation of the inverse FFT below
  for (ct=0; ct<HII_KSPACE_NUM_PIXELS; ct++){
#ifdef USE_HALO_FIELD
    M_coll_unfiltered[ct] /= (double)HII_TOT_NUM_PIXELS;
#endif //USE_HALO_FIELD
#ifdef USE_TS_IN_21CM
    xe_unfiltered[ct] /= (double)HII_TOT_NUM_PIXELS;
#endif //USE_TS_IN_21CM
#ifdef INHOMO_RECO
    N_rec_unfiltered[ct] /= (double)HII_TOT_NUM_PIXELS;
#endif //INHOMO_RECO
    deltax_unfiltered[ct] /= (HII_TOT_NUM_PIXELS+0.0);
  }
  fprintf(LOG, "end initial ffts, clock=%06.2f\n", (double)clock()/CLOCKS_PER_SEC);
  fflush(LOG);


 /*********************** DESTRUCTOR ************************************************/
  CLEANUP: if (*error_message != '\0'){
      fprintf(stderr, error_message);
      fprintf(LOG, error_message);
      fftwf_free(xH);
      fftwf_free(z_re);
      fftwf_free(Gamma12);
      fclose(LOG);
      fftwf_free(deltax_unfiltered);
      fftwf_free(deltax_filtered);
      fftwf_free(M_coll_unfiltered);
      fftwf_free(M_coll_filtered);
      fftwf_cleanup_threads();
      free_ps();
#ifdef INHOMO_RECO
      free_MHR();
#endif //INHOMO_RECO
      fftwf_free(xe_filtered);
      fftwf_free(xe_unfiltered);
      fftwf_free(N_rec_unfiltered);
      fftwf_free(N_rec_filtered);
      free(Fcoll);
#ifdef MINI_HALO
      free(Fcollm);
#endif
#ifndef SHARP_CUTOFF
      destroy_21cmMC_arrays();
#endif
      fftwf_destroy_plan(plan);
      fftwf_cleanup();
      
      return -1;
  }
    

  /*************************************************************************************/
  /*****************  END OF INITIALIZATION ********************************************/
  /*************************************************************************************/

  
  /*************************************************************************************/
  /***************** LOOP THROUGH THE FILTER RADII (in Mpc)  ***************************/
  /*************************************************************************************/
  // set the max radius we will use, making sure we are always sampling the same values of radius
  // (this avoids aliasing differences w redshift)
  R=fmax(R_BUBBLE_MIN, (cell_length_factor*BOX_LEN/(float)HII_DIM));
  while (R < fmin(MFP, L_FACTOR*BOX_LEN)) { R*= DELTA_R_HII_FACTOR;} 
  R /= DELTA_R_HII_FACTOR;
  LAST_FILTER_STEP = 0;
  while (!LAST_FILTER_STEP && (M_MIN < RtoM(R)) ){

    if ( ((R/DELTA_R_HII_FACTOR) <= (cell_length_factor*BOX_LEN/(float)HII_DIM)) || ((R/DELTA_R_HII_FACTOR) <= R_BUBBLE_MIN) ){
      LAST_FILTER_STEP = 1;
      R = fmax(cell_length_factor*BOX_LEN/(double)HII_DIM, R_BUBBLE_MIN);
    }

    fprintf(LOG, "before memcpy, clock=%06.2f\n", (double)clock()/CLOCKS_PER_SEC);
    fflush(LOG);
#ifdef USE_HALO_FIELD
    memcpy(M_coll_filtered, M_coll_unfiltered, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
#endif //USE_HALO_FIELD
#ifdef USE_TS_IN_21CM
    memcpy(xe_filtered, xe_unfiltered, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
#endif //USE_TS_IN_21CM
#ifdef INHOMO_RECO
    memcpy(N_rec_filtered, N_rec_unfiltered, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
#endif //INHOMO_RECO
    memcpy(deltax_filtered, deltax_unfiltered, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);


    // if this scale is not the size of the cells, we need to filter the fields
    if (!LAST_FILTER_STEP || (R > cell_length_factor*BOX_LEN/(double)HII_DIM) ){
      fprintf(LOG, "begin filter, clock=%06.2f\n", (double)clock()/CLOCKS_PER_SEC);
      fflush(LOG);
#ifdef USE_HALO_FIELD
      HII_filter(M_coll_filtered, HII_FILTER, R);
#endif //USE_HALO_FIELD
#ifdef USE_TS_IN_21CM
      HII_filter(xe_filtered, HII_FILTER, R);
#endif //USE_TS_IN_21CM
#ifdef INHOMO_RECO
      HII_filter(N_rec_filtered, HII_FILTER, R);
#endif //INHOMO_RECO
      HII_filter(deltax_filtered, HII_FILTER, R);
      fprintf(LOG, "end filter, clock=%06.2f\n", (double)clock()/CLOCKS_PER_SEC);
      fflush(LOG);
    }
    
    // do the FFT to get back to real space
    fprintf(LOG, "begin fft with R=%f, clock=%06.2f\n", R, (double)clock()/CLOCKS_PER_SEC);
    fflush(LOG);
#ifdef USE_HALO_FIELD
    plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)M_coll_filtered, (float *)M_coll_filtered, FFTW_ESTIMATE);
    fftwf_execute(plan);
#endif //USE_HALO_FIELD
#ifdef USE_TS_IN_21CM
    plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)xe_filtered, (float *)xe_filtered, FFTW_ESTIMATE);
    fftwf_execute(plan);
#endif //USE_TS_IN_21CM
#ifdef INHOMO_RECO
    plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)N_rec_filtered, (float *)N_rec_filtered, FFTW_ESTIMATE);
    fftwf_execute(plan);
#endif //INHOMO_RECO
    plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)deltax_filtered, (float *)deltax_filtered, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fprintf(LOG, "end fft with R=%f, clock=%06.2f\n", R, (double)clock()/CLOCKS_PER_SEC);
    fflush(LOG);


    // perform sanity checks to account for aliasing effects
    for (x=0; x<HII_DIM; x++){
      for (y=0; y<HII_DIM; y++){
        for (z=0; z<HII_DIM; z++){
      
          // delta cannot be less than -1
          *((float *)deltax_filtered + HII_R_FFT_INDEX(x,y,z)) = 
            FMAX(*((float *)deltax_filtered + HII_R_FFT_INDEX(x,y,z)) , -1+FRACT_FLOAT_ERR);

          // <N_rec> cannot be less than zero
#ifdef INHOMO_RECO
          *((float *)N_rec_filtered + HII_R_FFT_INDEX(x,y,z)) = 
          FMAX(*((float *)N_rec_filtered + HII_R_FFT_INDEX(x,y,z)) , 0.0);
#endif //INHOMO_RECO
        
          // collapsed mass cannot be less than zero
#ifdef USE_HALO_FIELD
          *((float *)M_coll_filtered + HII_R_FFT_INDEX(x,y,z)) = 
            FMAX(*((float *)M_coll_filtered + HII_R_FFT_INDEX(x,y,z)) , 0.0);
#endif //USE_HALO_FIELD

          // x_e has to be between zero and unity
#ifdef USE_TS_IN_21CM
          *((float *)xe_filtered + HII_R_FFT_INDEX(x,y,z)) = FMAX(*((float *)xe_filtered + HII_R_FFT_INDEX(x,y,z)) , 0);
          *((float *)xe_filtered + HII_R_FFT_INDEX(x,y,z)) = FMIN(*((float *)xe_filtered + HII_R_FFT_INDEX(x,y,z)) , 0.999);
#endif //USE_TS_IN_21CM
        }
      }
    } //  end sanity check in box
    
    // normalize the analytic collapse fractions if we operating on the density field
    f_coll = 0.0;
#ifdef MINI_HALO
    f_collm = 0.0;
#endif
    massofscaleR = RtoM(R);
#ifndef USE_HALO_FIELD
    fprintf(LOG, "begin f_coll normalization, clock=%06.2f\n", (double)clock()/CLOCKS_PER_SEC);
    fflush(LOG);
    temparg =  2*(pow(sigma_z0(M_MIN), 2) - pow(sigma_z0(massofscaleR), 2) );
    erfc_denom = sqrt(temparg);
        
#ifndef SHARP_CUTOFF
    initialiseGL_Nion(NGL_SFR, M_MIN, massofscaleR);
    initialise_Nion_spline(REDSHIFT, massofscaleR,M_MIN,M_MINa,ALPHA_STAR,ALPHA_ESC,F_STAR10,F_ESC10,Mlim_Fstar,Mlim_Fesc);
#ifdef MINI_HALO
    initialise_Nion_splinem(REDSHIFT, massofscaleR,M_MIN,ALPHA_STAR,M_MINm,Mcrit_atom,F_STAR10m,Mlim_Fstarm);
#endif //SHARP_CUTOFF
#endif //MINI_HALO
#endif //USE_HALO_FIELD

    for (x=0; x<HII_DIM; x++){
      for (y=0; y<HII_DIM; y++){
        for (z=0; z<HII_DIM; z++){
#ifdef USE_HALO_FIELD
          Splined_Fcoll = *((float *)M_coll_filtered + HII_R_FFT_INDEX(x,y,z)) / (massofscaleR*density_over_mean);
          Splined_Fcoll *= (4/3.0)*PI*pow(R,3) / pixel_volume;
#else //USE_HALO_FIELD
          density_over_mean = 1.0 + *((float *)deltax_filtered + HII_R_FFT_INDEX(x,y,z));
          if ( (density_over_mean - 1) < Deltac){ // we are not resolving collapsed structures
#ifndef SHARP_CUTOFF
            // Here again, 'Splined_Fcoll' and 'f_coll' are not the collpased fraction, but leave this name as is to simplify the variable name.
            // f_coll * ION_EFF_FACTOR = the number of IGM ionizing photon per baryon at a given overdensity.
            // see eq. (17) in Park et al. 2018
            Nion_Spline_density(density_over_mean - 1,&(Splined_Fcoll));
#ifdef MINI_HALO
            Nion_Spline_densitym(density_over_mean - 1,&(Splined_Fcollm));
#endif //MINI_HALO
#else //SHARP_CUTOFF
            erfc_num = (Deltac - (density_over_mean-1)) /  growth_factor;
            Splined_Fcoll = splined_erfc(erfc_num/erfc_denom);
#endif //SHARP_CUTOFF
          }
          else { // the entrire cell belongs to a collpased halo...  this is rare...
            Splined_Fcoll =  1.0;
#ifdef MINI_HALO
            Splined_Fcollm =  1.0;
#endif //MINI_HALO
          }
#endif //USE_HALO_FIELD

          // save the value of the collasped fraction into the Fcoll array
          Fcoll[HII_R_FFT_INDEX(x,y,z)] = Splined_Fcoll;
          f_coll += Splined_Fcoll;        
#ifdef MINI_HALO
          Fcollm[HII_R_FFT_INDEX(x,y,z)] = Splined_Fcollm;
          f_collm += Splined_Fcollm;
#endif
        }
      }
    } //  end loop through Fcoll box

    f_coll  /= (double) HII_TOT_NUM_PIXELS; // ave PS fcoll for this filter scale
    ST_over_PS = mean_f_coll_st/f_coll; // normalization ratio used to adjust the PS conditional collapsed fraction
#ifdef MINI_HALO
    f_collm /= (double) HII_TOT_NUM_PIXELS; // ave PS fcoll for this filter scale
    ST_over_PSm = mean_f_coll_stm/f_collm; // normalization ratio used to adjust the PS conditional collapsed fraction
#endif
    fprintf(LOG, "end f_coll normalization if, clock=%06.2f\n", (double)clock()/CLOCKS_PER_SEC);
    fflush(LOG);
    
    //     fprintf(stderr, "Last filter %i, R_filter=%f, fcoll=%f, ST_over_PS=%f, mean normalized fcoll=%f\n", LAST_FILTER_STEP, R, f_coll, ST_over_PS, f_coll*ST_over_PS);

    /****************************************************************************/
    /************  MAIN LOOP THROUGH THE BOX FOR THIS FILTER SCALE **************/
    /****************************************************************************/
    fprintf(LOG, "Start of the main loop through the box for this filter scale, clock=%06.2f\n", (double)clock()/CLOCKS_PER_SEC);
    fflush(LOG);
    // now lets scroll through the filtered box
    rec = ave_xHI_xrays = ave_den = ave_fcoll = std_xrays = 0;
    ion_ct=0;
    xHI_from_xrays = 1;
    Gamma_R_prefactor  = pow(1+REDSHIFT, 2) * (R*CMperMPC) * SIGMA_HI * ALPHA_UVB / (ALPHA_UVB+2.75) * N_b0 * ION_EFF_FACTOR / 1.0e-12;
#ifdef MINI_HALO
    Gamma_R_prefactorm = Gamma_R_prefactor / ION_EFF_FACTOR * ION_EFF_FACTOR_MINI;
#endif //MINI_HALO

    for (x=0; x<HII_DIM; x++){
      for (y=0; y<HII_DIM; y++){
        for (z=0; z<HII_DIM; z++){

          density_over_mean = 1.0 + *((float *)deltax_filtered + HII_R_FFT_INDEX(x,y,z));

          f_coll  = ST_over_PS * Fcoll[HII_R_FFT_INDEX(x,y,z)];
#ifdef MINI_HALO
          f_collm = ST_over_PSm * Fcollm[HII_R_FFT_INDEX(x,y,z)];
#endif //MINI_HALO
      
          // if this is the last filter step, prepare to account for poisson fluctuations in the sub grid halo number...
          // this is very approximate as it doesn't sample the halo mass function but merely samples a number of halos of a characterisic mass
          if (LAST_FILTER_STEP){
            ave_M_coll_cell = f_coll * pixel_mass * density_over_mean;
            ave_N_min_cell = ave_M_coll_cell / M_MIN; // ave # of M_MIN halos in cell
            N_halos_in_cell = (int) gsl_ran_poisson(r, N_POISSON);
          }

#ifdef INHOMO_RECO
          dfcolldt  = f_coll  / t_ast;
          Gamma_R = Gamma_R_prefactor * dfcolldt;
#ifdef MINI_HALO
          dfcolldtm = f_collm / t_ast;
          Gamma_R  += Gamma_R_prefactorm * dfcolldtm;
#endif //MINI_HALO
          rec = (*((float *)N_rec_filtered + HII_R_FFT_INDEX(x,y,z))); // number of recombinations per mean baryon
          rec /= density_over_mean; // number of recombinations per baryon inside <R>
#endif //INHOMO_RECO
      
          // adjust the denominator of the collapse fraction for the residual electron fraction in the neutral medium
#ifdef USE_TS_IN_21CM
          xHI_from_xrays =  1 - (*((float *)xe_filtered + HII_R_FFT_INDEX(x,y,z)));
#endif //USE_TS_IN_21CM

          // check if fully ionized!
#ifdef MINI_HALO
          if ( (f_coll*ION_EFF_FACTOR + f_collm*ION_EFF_FACTOR_MINI > xHI_from_xrays*(1.0+rec)) )
#else //MINI_HALO
          if ( (f_coll*ION_EFF_FACTOR > xHI_from_xrays*(1.0+rec)) )
#endif //MINI_HALO
          {
        
            // if this is the first crossing of the ionization barrier for this cell (largest R), record the gamma
            // this assumes photon-starved growth of HII regions...  breaks down post EoR
#ifdef INHOMO_RECO
            if (xH[HII_R_INDEX(x, y, z)] > FRACT_FLOAT_ERR){
              Gamma12[HII_R_INDEX(x, y, z)] = Gamma_R;
              aveR += R;
              Rct++;
            }

            // keep track of the first time this cell is ionized (earliest time)
            if (z_re[HII_R_INDEX(x, y, z)] < 0)
              z_re[HII_R_INDEX(x, y, z)] = REDSHIFT;          
#endif //INHOMO_RECO

        
            // FLAG CELL(S) AS IONIZED
            if (FIND_BUBBLE_ALGORITHM == 2) // center method
              xH[HII_R_INDEX(x, y, z)] = 0;
            else if (FIND_BUBBLE_ALGORITHM == 1) // sphere method
              update_in_sphere(xH, HII_DIM, R/BOX_LEN, x/(HII_DIM+0.0), y/(HII_DIM+0.0), z/(HII_DIM+0.0));
            else{
              fprintf(stderr, "Incorrect choice of find bubble algorithm set in ANAL_PARAMS.H.\nSetting center method...");
              fprintf(LOG, "Incorrect choice of find bubble algorithm set in ANAL_PARAMS.H.\nSetting center method...");
              xH[HII_R_INDEX(x, y, z)] = 0;
            }
        
          } // end ionized
        
          // If not fully ionized, then assign partial ionizations 
          else if (LAST_FILTER_STEP && (xH[HII_R_INDEX(x, y, z)] > TINY)){

#ifndef USE_HALO_FIELD
            if (ave_N_min_cell < N_POISSON){ // add poissonian fluctuations to the nalo number
#ifdef MINI_HALO
              // TODO: check with andrei
              f_coll  = N_halos_in_cell * (ave_M_coll_cell / (float) N_POISSON ) / (pixel_mass*density_over_mean) * ((f_coll / (f_coll + f_collm)));
              f_collm = N_halos_in_cell * (ave_M_coll_cell / (float) N_POISSON ) / (pixel_mass*density_over_mean) - f_coll;
#else //MINI_HALO
              f_coll = N_halos_in_cell * (ave_M_coll_cell / (float) N_POISSON ) / (pixel_mass*density_over_mean);
#endif //MINI_HALO
              if (ave_M_coll_cell  < (M_MIN/5.0)){
                f_coll = 0;
#ifdef MINI_HALO
                f_collm = 0;
#endif //MINI_HALO
              }
            }
#endif //USE_HALO_FIELD
          
#ifdef MINI_HALO
            res_xH = xHI_from_xrays - f_coll * ION_EFF_FACTOR - f_collm * ION_EFF_FACTOR_MINI;
#else //MINI_HALO
            res_xH = xHI_from_xrays - f_coll * ION_EFF_FACTOR;
#endif //MINI_HALO
            // and make sure fraction doesn't blow up for underdense pixels
            if (res_xH < 0)
              res_xH = 0;
            else if (res_xH > 1)
              res_xH = 1;
          
            xH[HII_R_INDEX(x, y, z)] = res_xH;
          } // end partial ionizations at last filtering step

          /*
          if ((x==0) && (y==19) && (z==107)){ //DEBUGGING
            fprintf(stderr, "R=%g Mpc, <Delta>_R=%g, <fcoll>_R=%g, <rec>_V=%g, xH=%g, Gamma12=%g, z_re=%g\n",
                R, density_over_mean, f_coll, rec, xH[HII_R_INDEX(x, y, z)], Gamma12[HII_R_INDEX(x, y, z)], z_re[HII_R_INDEX(x, y, z)]);
            fprintf(LOG, "R=%g Mpc, <Delta>_R=%g, <fcoll>_R=%g, <rec>_V=%g, xH=%g, Gamma12=%g, z_re=%g\n",
                R, density_over_mean, f_coll, rec, xH[HII_R_INDEX(x, y, z)], Gamma12[HII_R_INDEX(x, y, z)], z_re[HII_R_INDEX(x, y, z)]);
            fflush(NULL);
          }
          */
        
        } // z
      } // y
    } // x
    R /= DELTA_R_HII_FACTOR;
  } // END OF LOOP THROUGH FILTER RADII

  // find the neutral fraction
  global_xH = 0;
  for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++)
    global_xH += xH[ct];
  global_xH /= (float)HII_TOT_NUM_PIXELS;


  // update the N_rec field
#ifdef INHOMO_RECO
  //fft to get the real N_rec  and delta fields
  plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)N_rec_unfiltered, (float *)N_rec_unfiltered, FFTW_ESTIMATE);
  fftwf_execute(plan);
  plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)deltax_unfiltered, (float *)deltax_unfiltered, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);
  fftwf_cleanup();
  for (x=0; x<HII_DIM; x++){
    for (y=0; y<HII_DIM; y++){
      for (z=0; z<HII_DIM; z++){

        density_over_mean = 1.0 + (*((float *)deltax_unfiltered + HII_R_FFT_INDEX(x,y,z)));
        z_eff = (1+REDSHIFT) * pow(density_over_mean, 1.0/3.0) - 1;
        dNrec = splined_recombination_rate(z_eff, Gamma12[HII_R_INDEX(x,y,z)]) * fabs_dtdz * ZSTEP * (1 - xH[HII_R_INDEX(x,y,z)]);
          *((float *)N_rec_unfiltered + HII_R_FFT_INDEX(x,y,z)) += dNrec; 

        /*
        if ((x==0) && (y==19) && (z==107)){ //DEBUGGING
          fprintf(stderr, "Delta=%g, z_eff=%g, dNrec=%g, Nrec=%g\n\n", density_over_mean, z_eff, dNrec, *((float *)N_rec_unfiltered + HII_R_FFT_INDEX(x,y,z)));
          fprintf(LOG, "Delta=%g, z_eff=%g, dNrec=%g, Nrec=%g\n\n", density_over_mean, z_eff, dNrec, *((float *)N_rec_unfiltered + HII_R_FFT_INDEX(x,y,z)));
          fflush(NULL);
        }
        */
    
      }
    }
  }
#endif  //INHOMO_RECO

  /*********************************************************************************************/
  /************  END OF CALCULATION FOR THIS REDSHIFT.  NOW WE WILL PRINT OUT THE BOXES ********/
  /*********************************************************************************************/


  fprintf(stderr, "ave R used to compute Gamma in ionized regions = %g Mpc\n", aveR / (double)Rct);
  fprintf(LOG, "ave R used to compute Gamma in ionized regions = %g Mpc\n", aveR / (double)Rct);

  
#ifdef INHOMO_RECO
  // N_rec box
  sprintf(filename, "../Boxes/Nrec_z%06.2f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT,HII_FILTER, MFP, HII_DIM, BOX_LEN);
  if (!(F = fopen(filename, "wb"))){
    sprintf(error_message, "find_HII_bubbles: ERROR: unable to open file for writting Nrec box!\n");
    goto CLEANUP;
  }
  else{
    for (i=0; i<HII_DIM; i++){
      for (j=0; j<HII_DIM; j++){
        for (k=0; k<HII_DIM; k++){
          if(fwrite((float *)N_rec_unfiltered + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, F)!=1){
            sprintf(error_message, "find_HII_bubbles.c: Write error occured while writting N_rec box.\n");
            goto CLEANUP;
          }
        }
      }
    }
    fclose(F);
    F = NULL;
  }

  // Write z_re in the box
  sprintf(filename, "../Boxes/z_first_ionization_z%06.2f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, HII_FILTER, MFP, HII_DIM, BOX_LEN);
  if (!(F = fopen(filename, "wb"))){
    sprintf(error_message, "find_HII_bubbles: ERROR: unable to open file for writting z_re box!\n");
    goto CLEANUP;
  }
  else{
    if (mod_fwrite(z_re, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
      sprintf(error_message, "find_HII_bubbles: ERROR: unable to open file for writting z_re box!\n");
      goto CLEANUP;
    }
  fclose(F);
  F = NULL;
  }

  // Gamma12 box
  sprintf(filename, "../Boxes/Gamma12aveHII_z%06.2f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, HII_FILTER, MFP, HII_DIM, BOX_LEN);
  if (!(F = fopen(filename, "wb"))){
    sprintf(error_message, "find_HII_bubbles: ERROR: unable to open file for writting gamma box!\n");
    goto CLEANUP;
  }
  else{
    if (mod_fwrite(Gamma12, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
      sprintf(error_message, "find_HII_bubbles.c: Write error occured while writting gamma box.\n");
      goto CLEANUP;
    }
    fclose(F);
    F = NULL;
  }
#endif //INHOMO_RECO
    
  // print out the xH box
  switch(FIND_BUBBLE_ALGORITHM){
    case 2:
#ifndef SHARP_CUTOFF
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
#ifdef USE_HALO_FIELD
      sprintf(filename, "../Boxes/xH_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_f_star10m%.4f_f_esc10m%.4f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, F_STAR10m, F_ESC10m, HII_FILTER, MFP, HII_DIM, BOX_LEN);
#else //USE_HALO_FIELD
      sprintf(filename, "../Boxes/xH_nohalos_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_f_star10m%.4f_f_esc10m%.4f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, F_STAR10m, F_ESC10m, HII_FILTER, MFP, HII_DIM, BOX_LEN);
#endif //USE_HALO_FIELD
#else //INHOMO_FEEDBACK
#ifdef REION_SM
#ifdef USE_HALO_FIELD
      sprintf(filename, "../Boxes/xH_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_Msm_f_star10m%.4f_f_esc10m%.4f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, F_STAR10m, F_ESC10m, HII_FILTER, MFP, HII_DIM, BOX_LEN);
#else //USE_HALO_FIELD
      sprintf(filename, "../Boxes/xH_nohalos_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_Msm_f_star10m%.4f_f_esc10m%.4f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, F_STAR10m, F_ESC10m, HII_FILTER, MFP, HII_DIM, BOX_LEN);
#endif //USE_HALO_FIELD
#else //REION_SM
#ifdef USE_HALO_FIELD
      sprintf(filename, "../Boxes/xH_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_Mturn%.2e_f_star10m%.4f_f_esc10m%.4f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, F_STAR10m, F_ESC10m, HII_FILTER, MFP, HII_DIM, BOX_LEN);
#else //USE_HALO_FIELD
      sprintf(filename, "../Boxes/xH_nohalos_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_Mturn%.2e_f_star10m%.4f_f_esc10m%.4f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, F_STAR10m, F_ESC10m, HII_FILTER, MFP, HII_DIM, BOX_LEN);
#endif //USE_HALO_FIELD
#endif //REION_SM
#endif //INHOMO_FEEDBACK
#else  //MINI_HALO
#ifdef USE_HALO_FIELD
      sprintf(filename, "../Boxes/xH_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_Mturn%.2e_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, HII_FILTER, MFP, HII_DIM, BOX_LEN);
#else //USE_HALO_FIELD
      sprintf(filename, "../Boxes/xH_nohalos_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_Mturn%.2e_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, HII_FILTER, MFP, HII_DIM, BOX_LEN);
#endif //USE_HALO_FIELD
#endif //MINI_HALO
#else //SHARP_CUTOFF
#ifdef USE_HALO_FIELD
      sprintf(filename, "../Boxes/xH_z%06.2f_nf%f_eff%.1f_effPLindex0_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, ION_EFF_FACTOR, HII_FILTER, M_MIN, MFP, HII_DIM, BOX_LEN);
#else //USE_HALO_FIELD
      sprintf(filename, "../Boxes/xH_nohalos_z%06.2f_nf%f_eff%.1f_effPLindex0_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, ION_EFF_FACTOR, HII_FILTER, M_MIN, MFP, HII_DIM, BOX_LEN);
#endif //USE_HALO_FIELD
#endif //SHARP_CUTOFF
      break;
    default:
#ifndef SHARP_CUTOFF
#ifdef MINI_HALO
#ifdef INHOMO_FEEDBACK
#ifdef USE_HALO_FIELD
      sprintf(filename, "../Boxes/sphere_xH_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_f_star10m%.4f_f_esc10m%.4f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, F_STAR10m, F_ESC10m, HII_FILTER, MFP, HII_DIM, BOX_LEN);
#else //USE_HALO_FIELD
      sprintf(filename, "../Boxes/sphere_xH_nohalos_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_f_star10m%.4f_f_esc10m%.4f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, F_STAR10m, F_ESC10m, HII_FILTER, MFP, HII_DIM, BOX_LEN);
#endif //USE_HALO_FIELD
#else //INHOMO_FEEDBACK
#ifdef REION_SM
#ifdef USE_HALO_FIELD
      sprintf(filename, "../Boxes/sphere_xH_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_Msm_f_star10m%.4f_f_esc10m%.4f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, F_STAR10m, F_ESC10m, HII_FILTER, MFP, HII_DIM, BOX_LEN);
#else //USE_HALO_FIELD
      sprintf(filename, "../Boxes/sphere_xH_nohalos_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_Msm_f_star10m%.4f_f_esc10m%.4f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, F_STAR10m, F_ESC10m, HII_FILTER, MFP, HII_DIM, BOX_LEN);
#endif //USE_HALO_FIELD
#else //REION_SM
#ifdef USE_HALO_FIELD
      sprintf(filename, "../Boxes/sphere_xH_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_Mturn%.2e_f_star10m%.4f_f_esc10m%.4f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, F_STAR10m, F_ESC10m, HII_FILTER, MFP, HII_DIM, BOX_LEN);
#else //USE_HALO_FIELD
      sprintf(filename, "../Boxes/sphere_xH_nohalos_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_Mturn%.2e_f_star10m%.4f_f_esc10m%.4f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, F_STAR10m, F_ESC10m, HII_FILTER, MFP, HII_DIM, BOX_LEN);
#endif //USE_HALO_FIELD
#endif //REION_SM
#endif //INHOMO_FEEDBACK
#else //MINI_HALO
#ifdef USE_HALO_FIELD
      sprintf(filename, "../Boxes/sphere_xH_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_Mturn%.2e_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, HII_FILTER, MFP, HII_DIM, BOX_LEN);
#else //USE_HALO_FIELD
      sprintf(filename, "../Boxes/sphere_xH_nohalos_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_Mturn%.2e_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, HII_FILTER, MFP, HII_DIM, BOX_LEN);
#endif //USE_HALO_FIELD
#endif //MINI_HALO
#else //SHARP_CUTOFF
#ifdef USE_HALO_FIELD
      sprintf(filename, "../Boxes/sphere_xH_z%06.2f_nf%f_eff%.1f_effPLindex0_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, ION_EFF_FACTOR, HII_FILTER, M_MIN, MFP, HII_DIM, BOX_LEN);
#else //USE_HALO_FIELD
      sprintf(filename, "../Boxes/sphere_xH_nohalos_z%06.2f_nf%f_eff%.1f_effPLindex0_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, ION_EFF_FACTOR, HII_FILTER, M_MIN, MFP, HII_DIM, BOX_LEN);
#endif //USE_HALO_FIELD
#endif //SHARP_CUTOFF
  }
  if (!(F = fopen(filename, "wb"))){
    fprintf(stderr, "find_HII_bubbles: ERROR: unable to open file %s for writting!\n", filename);
    fprintf(LOG, "find_HII_bubbles: ERROR: unable to open file %s for writting!\n", filename);
    global_xH = -1;
  }
  else {
    fprintf(LOG, "Neutral fraction is %f\nNow writting xH box at %s\n", global_xH, filename);
    fprintf(stderr, "Neutral fraction is %f\nNow writting xH box at %s\n", global_xH, filename);
    fflush(LOG);
    if (mod_fwrite(xH, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
      fprintf(stderr, "find_HII_bubbles.c: Write error occured while writting xH box.\n");
      fprintf(LOG, "find_HII_bubbles.c: Write error occured while writting xH box.\n");
      global_xH = -1;
    }
    fclose(F);
    F = NULL;
  }

  fprintf(stderr, error_message);
  fprintf(LOG, error_message);
  fftwf_free(xH);
  fftwf_free(z_re);
  fftwf_free(Gamma12);
  fclose(LOG);
  fftwf_free(deltax_unfiltered);
  fftwf_free(deltax_filtered);
  fftwf_free(M_coll_unfiltered);
  fftwf_free(M_coll_filtered);
  fftwf_cleanup_threads();
  if (F){ fclose(F);}
  free_ps();
  free_MHR();
  fftwf_free(xe_filtered);
  fftwf_free(xe_unfiltered);
  fftwf_free(N_rec_unfiltered);
  fftwf_free(N_rec_filtered);
#ifndef SHARP_CUTOFF
  destroy_21cmMC_arrays();
#endif
      
  return 0;
}
