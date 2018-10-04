#ifndef _ANAL_PARAM_H_
#define _ANAL_PARAM_H_

/**************************************************************************************************************/
/**************** BEGIN USER CHANGABLE DEFINITIONS   **********************************************************/

/***********************************************
  GALAXY PARAMETERS
  NEW in v2:

  The ionizing efficiency zeta is allowed to scale with the halo mass:
  \zeta = N_{\gammaUV} x f_{esc}(M_{h}) x f_{\ast}(M_{h}),

  where:
  M_{h} is a halo mass,
  f_{\ast} = f_{\ast,10} x (M_{h}/10^10 Msun)^\alpha_{\ast}
           = STELLAR_BARYON_FRAC x (M_{h}/10^10 Msun)^STELLAR_BARYON_PL,
  and
  f_{esc} = f_{esc,10} x (M_{h}/10^10 Msun)^\alpha_{esc}
          = ESC_FRAC x (M_{h}/10^10 Msun)^ESC_PL.

  The number of halos hosting active star-forming galaxies (i.e. the duty cycle) is assumed to
  exponentially decrease below M_TURNOVER Msun, : fduty \propto e^(- M_TURNOVER / M)

  See eq.(16) and equations in Section 2 of Park et al. (2018): https://arxiv.org/pdf/1809.08995.pdf

  This code still support the constant ionizing efficiency with sharp cut-off for minimum halo mass.
  See 'SHARP_CUTOFF' in this file.

***********************************************/
#define N_GAMMA_UV (float) (5000) // number of ionizing photons per stellar baryon
#define STELLAR_BARYON_FRAC (double) (0.05) // defined for 10^10 Msun halos
#define STELLAR_BARYON_PL (double) (0.5) //power law index w/ halo mass
#define ESC_FRAC (double) (0.1) // ionizing escape fraction, defined for 10^10 Msun halos
#define ESC_PL (double) (-0.5) //power law index w/ halo mass
//#define M_TURNOVER (double) (5e8) // Halo mass threshold for efficient star formation (in Msun). 

/***********************************************
  MINI HALO PARAMETERS
  NEW in v2.1:

  The ionizing efficiency zeta is allowed to scale with the halo mass
  in two populations -- atomic and molecular cooling halos:
  \zeta_ato = N_{\gammaUV}     x f_{esc}     x f_{\ast} 
  \zeta_mol = N_{\gammaUV,mol} x f_{esc,mol} x f_{\ast,mol}
  \zeta     = \zeta_ato + \zeta_mol

  where:
  f_{\ast}     = STELLAR_BARYON_FRAC x (M_{halo}/10^10 M_{sol})^STELLAR_BARYON_PL
  f_{esc}      = ESC_FRAC            x (M_{halo}/10^10 M_{sol})^ESC_PL

  f_{\ast,mol} = STELLAR_BARYON_FRAC_MINI x (M_{halo}/10^10 M_{sol})^STELLAR_BARYON_PL
  f_{esc, mol} = ESC_FRAC_MINI

  Note that in the current version
  1. the scaling of f_{\ast,mol} is the same as f_\ast
  2. f_{esc, mol} is a constant
  3. t_STAR is the same in the two populations
  4. f_{\ast}, f_{esc}, f_{\ast,mol}, f_{esc, mol} need to be smaller than 1

  -----------
  DUTY CYCLES
  -----------
  1. atomic cooling halos:
  f_\mathrm{duty,atom} = \exp[-M_\mathrm{turn,atom}(z,\mathbf{r})/M_\mathrm{h}]
  with 
  M_\mathrm{turn,atom}(z,\mathbf{r}) = \max[M_\mathrm{crit,RE}(z,\mathbf{r}),M_\mathrm{crit,atom}(z)]
  
  2. molecular cooling halos:
  f_\mathrm{duty,mol}  = \exp[-M_\mathrm{turn,mol,low}(z,\mathbf{r})/M_\mathrm{h}] \times 
                         \exp[-M_\mathrm{h}/M_\mathrm{turn,mol,high}(z)]
  with
  M_\mathrm{turn,mol,high}(z)           = M_\mathrm{crit,atom}(z)
  M_\mathrm{turn,mol,low}(z,\mathbf{r}) = \max[M_\mathrm{crit,RE}(z,\mathbf{r}),M_\mathrm{crit,LW}(z,\mathbf{r})]


  ------------------------
  FEEDBACK AND SUPPRESSION
  ------------------------
  1. atomic cooling threshold:
  M_\mathrm{crit,atom} corresponds to the atomic cooling threshold (T~1e4K), which essentially split halos
                       into atomic and molecular cooling halos.
  M_\mathrm{crit,atom} = 1e8/h * ((1e4 / 1.98e4) / (0.59/0.6) / ((1.+z)/10))**(1.5)
                          ((Omega_m / Omega_m_z(z)) * (Delta_c(z)/18/pi**2))**(1./3)
  with 
  Omega_m_z(z) = 1./(1.+(Omega_L/Omega_m/(1.+z)**3))
  Delta_c(z)   = 18 * pi**2 + 82 * (Omega_m_z(z) - 1.0) -39 * (Omega_m_z(z) - 1.0) ** 2

  2. Reionization feedback (following Sobacchi & Mesinger 2013a,b; 2014)
  M_\mathrm{crit,RE}(z,\mathbf{r}) = 3e9 \mathrm{M}_\odot [\Gamma_{\mathrm{halo, HII}}(z,\mathbf{r})/1e-12 s^{-1}]^0.17
                                     \times ((1+z)/10)^-2.1 \times {1-[(1+z)/1+z_\mathrm{IN}(z,\mathbf{r})]^2.0}^2.5

  3. Lyman-Werner (LW) suppression (following Visbal, Haiman & Bryan 2017)
  M_\mathrm{crit,LW}(z,\mathbf{r}) = 2.5e5 \mathrm{M}_\odot \times {1+6.96[4PI J_\mathrm{LW,21}(z,\mathbf{r})]^0.47}
                                     \times ((1+z)/26)^-1.5

***********************************************/
/*
  This version still keeps support for the constant ionizing efficiency with sharp cut-off for minimum halo mass.
  If SHARP_CUTOFF = 1, use the constant ionizing efficency parameter \zeta = N_{\gammaUV} x f_{esc} x f_{\ast}
  with sharp cut-off at M_TURNOVER.
  New in v2.1: Now SHARP_CUTOFF works as a flag, commit it if do not want it
*/
//#define SHARP_CUTOFF 1
#ifndef SHARP_CUTOFF
  #define MINI_HALO 1             // operate as a flag, comment this line if MINI_HALO version is not required
  #ifdef MINI_HALO
    #define N_GAMMA_UV_MINI (float) (50000) // number of ionizing photons per stellar baryon for mini halos
    #define STELLAR_BARYON_FRAC_MINI (double) (0.01) // defined for 10^10 Msun mini halos
    #define ESC_FRAC_MINI (double) (0.5)    // ionizing escape fraction for mini halos
  
    #define INHOMO_FEEDBACK 1 // operate as a flag, comment this line if homogeneous feedback is adopted
    #ifdef INHOMO_FEEDBACK 
      // LW suppression follows ?
      // --
      // use Sobacchi & Mesinger 2013a; 2014 prescription for reionization
      #define REION_SM13_M0 (double) 3e9
      #define REION_SM13_A  (double) 0.17
      #define REION_SM13_B  (double) -2.1
      #define REION_SM13_C  (double) 2.0
      #define REION_SM13_D  (double) 2.5
    #else  // let's do homogeneous feedback
      // use Wise+12 prescription for LW suppression
      #define LW_WISE12_A (double) -2.356
      #define LW_WISE12_B (double) 0.4562
    #define LW_WISE12_C (double) -0.02680
    #define LW_WISE12_D (double) 5.882e-4
    #define LW_WISE12_E (double) -5.056e-6
    // --
    #define REION_SM 1 // operate as a flag, ifdef use obacchi & Mesinger 2013b for homogeneous reionization 
    #ifndef REION_SM   // otherwise use the free parameter M_TURNOVER for homogeneous reionization
      #define M_TURNOVER (double) (5e8) 
    #endif  // end of ifndef REION_SM
  #endif  // end of ifdef INHOMO_FEEDBACK
#endif  // end of ifdef MINI_HALO
#endif // end of ifdef SHARP_CUTOFF

/**********************************************
  RECOMBINATION PARAMETERS
***********************************************/
/*
  INHOMO_RECO = 1 compute inhomogeneous recombinations, ala Sobacchi & Mesinger (2014).
  This *requires* stepping down in redshift, e.g. using drive_logZscroll_Ts.c or  drive_zscroll_noTs.c
*/
#define INHOMO_RECO (int) (1)

/*
  Typical star-formation timescale, in units of the Hubble time.  This is used to compute the UVB if inhomo_reco is on.
  This is also used when comparing againts observed UV LF w. 21CMMC, 
  as the UV luminosity is proportional to the stellar mass / (t*H^-1)
 */
#define t_STAR (float) (0.5)

/*
  Power law index of the UVB during the EoR.  This is only used if INHOMO_RECO is on, in order to compute the
  local mean free path inside the cosmic HII regions
 */
#define ALPHA_UVB (float) (5)


/*****************   END OF MAIN PHYSICAL PARAMETERS.  ADDITIONAL OPTIONS ARE PRESENT BELOW THIS LINE *********/
/**************************************************************************************************************/
/**************************************************************************************************************/


/*
  Maximum horizon for ionizing photons, in comoving Mpc.
  The bubble filtering procedure uses the minimum
  value of L_FACTOR*BOX_LEN and R_BUBBLE_MAX as the maximum allowed
  radius of a single, cohesive ionized region. (in practice, the ionization
  field is extreemly insensitive to this choice).

  If setting IHOMO_RECO below to 1, R_BUBBLE_MAX 
  physically corresponds to the mfp at the end of the EoR.
  In this case, it should be at least ~20Mpc, with results being extremely insensitive.
*/
#define R_BUBBLE_MAX (float) (50)


/*
  EVOLVE_DENSITY_LINEARLY = 1, evolve the density field with linear theory.
                If choosing this option, make sure that your cell size is
        in the linear regime at the redshift of interest

  EVOLVE_DENSITY_LINEARLY = 0, evolve the density field with 1st order perturbation theory.
                If choosing this option, make sure that you resolve small 
        enough scales, roughly we find BOX_LEN/DIM should be < 1Mpc
 */
#define EVOLVE_DENSITY_LINEARLY (int) (0)


/*
  If set to 1, the ZA density field is additionally smoothed (asside from the implicit
  boxcar smoothing performed when re-binning the ICs from DIM to HII_DIM) with a Gaussian
  filter of width R_smooth_density*BOX_LEN/HII_DIM.  The implicit boxcar smoothing in perturb_field.c
  bins the density field on scale DIM/HII_DIM, similar to what Lagrangian codes do 
  when constructing Eulerian grids. In other words, the density field, \delta,
  is quantized into (DIM/HII_DIM)^3 values. If your usage requires smooth density fields,
  it is recommended to set SMOOTH_EVOLVED_FIELD to 1.  This also decreases the shot noise
  present in all grid based codes, though it overcompensates by an effective loss in resolution.
  New in v1.1
 */
#define SMOOTH_EVOLVED_DENSITY_FIELD (int) (1)
#define R_smooth_density (float) (0.2)

/*
  Use second-order Lagrangian perturbation theory (2LPT).
  Set this to 1 if the density field or the halo positions are extrapolated to low redshifts.
  The current implementation is very naive and add a factor ~6 to the memory requirements.
  Reference: Scoccimarro R., 1998, MNRAS, 299, 1097-1118 Appendix D  
*/
#define SECOND_ORDER_LPT_CORRECTIONS (int) (1)

/*
  Allows one to set a flag allowing find_HII_bubbles to skip constructing the 
  ionization field if it is estimated that the mean neutral fraction, <xH>, is
  within HII_ROUND_ERR of 1. In other words, if <xH> > 1-HII_ROUND_ERR, 
  then find_HII_bubbles just prints a homogeneous xHI field
  of  1's.

  This is a new option in v1.1. Previous versions had a hardcoded value of 1e-15.
*/
#define HII_ROUND_ERR (float) (1e-3)


/*
  Choice of:
  1 - Mesinger & Furlanetto 2007 method of overlaping spheres:
  paint an ionized sphere with radius R, centered on pixel
  where R is filter radius
  This method, while somewhat more accurate, is slower than (2) especially
  in mostly ionized unverses, so only use for lower resolution boxes (HII_DIM<~400)
  2 - Center pixel only method (Zahn et al. 2007). this is faster.
*/
#define FIND_BUBBLE_ALGORITHM (int) (2)

/* 
   Minimum radius of an HII region in cMpc.  One can set this to 0, but should be careful with
   shot noise if the find_HII_bubble algorithm is run on a fine, non-linear density grid.
 */
#define R_BUBBLE_MIN (float) (L_FACTOR*1)


/*
  1 = Use the filtered and adjusted halo field to construct the ionization and 21-cm fields
  0 = Use the mean collapse fraction (eq. 2 in Barkana & Loeb 2008) to construct ...

  Note: option (1) is more accurate as the halo field captures the Poisson and clustering signature
  which dominates the collapse fraction field on small scales (i.e. where the mean halo number 
  is small).  However, option (0) is faster as it bypasses the halo finder, and it becomes
  increasingly accurate as the typical bubble size increases.
*/
#define USE_HALO_FIELD (int) (0)


/*
  If not using the halo field to generate HII regions, we provide the option of
  including Poisson scatter in the number of sources obtained through the conditional
  collapse fraction (which only gives the *mean* collapse fraction on a particular
  scale.  If the predicted mean collapse fraction is < N_POISSON * M_MIN,
  then Poisson scatter is added to mimic discrete halos on the subgrid scale
  (see Zahn+ 2010).

  NOTE: If you are interested in snapshots of the same realization at several redshifts,
  it is recommended to turn off this feature, as halos can stocastically
  "pop in and out of" existance from one redshift to the next...
*/
#define N_POISSON (int) (-1)


/*                                                                                                    
  Parameter choice of whether to use velocity corrections in 21-cm fields
  1=use velocities in delta_T; 0=do not use velocities

and 

  Maximum velocity gradient along the line of sight in units of the hubble parameter at z.
  This is only used in computing the 21cm fields.
  Note, setting this too high can add spurious 21cm power in the early stages, due to the
  1-e^-tau ~ tau approximation (see my 21cm intro paper and mao+2011).  However, this is still
  a good approximation at the <~10% level. 
  redshift space distortions.
  NOTE: since v1.4, the default is to ignore this approximation, using SUBCELL_RSD below
*/
#define T_USE_VELOCITIES (int) (1)
#define MAX_DVDR (float) (0.2)


/*
 Whether to include line-of-sight redshift space distortions. 
 
 This perturbs the 21cm brightness temperature into subcells (number set by N_RSD_STEPS), perturbs the 
 cell boundaries according to the peculiar velocity and re-grids the perturbed cells back to the
 original resolution.
 
 Uses the approach outlined in Jensen et al. (2013).
*/
#define SUBCELL_RSD (int) (1)

/*
 Defines the number of subcells to be used for calculating the redshift space distortions.
 
 Found 20 to be able to provide sufficient accuracy, Greig & Mesinger, 2018
*/
#define N_RSD_STEPS (int) (20)


/*
  0 = plot 21cm temperature power spectrum in non-dimensional units
  1 = plot 21cm...  in mK^2
*/
#define DIMENSIONAL_T_POWER_SPEC (int) (1)

#define DELTA_R_FACTOR (float) (1.1) // factor by which to scroll through filter radius for halos

#define DELTA_R_HII_FACTOR (float) (1.1) // factor by which to scroll through filter radius for bubbles

/*
Factor of the halo's radius, R, so that the effective radius
is R_eff = R_OVERLAP_FACTOR * R.  Halos whose centers are less than
R_eff away from another halo are not allowed.
R_OVERLAP_FACTOR = 1 is fully disjoint 
R_OVERLAP_FACTOR = 0 means that centers are allowed to luy on the edges of
neighboring halos
*/
#define R_OVERLAP_FACTOR (float) (1.0)

/*
0 = delta_crit is constant 1.68
1 = delta_crit is the sheth tormen ellipsoidal collapse
    correction to delta_crit (see ps.c)
*/
#define DELTA_CRIT_MODE (int) (1)


/*
Filter for the density field used to generate the halo field with EPS
0 = use real space top hat filter to smooth density field
1 = use k-space top hat filter
2 = use gaussian filter
*/
#define HALO_FILTER (int) (0)

/*
Filter for the Halo or density field used to generate ionization field
0 = use real space top hat filter
1 = use k-space top hat filter
2 = use gaussian filter
*/
#define HII_FILTER (int) (1)

/*
  0 = don't run the optimization code for zscroll (only useful at low-z)
  1 = do run ...
*/
#define OPTIMIZE (int) (0)
/*
  Minimum mass for which the optimization algorithm will be used
*/
#define OPTIMIZE_MIN_MASS (float) (1e11)

#define SIZE_RANDOM_SEED (-23456789) // seed for the size dist random number generator
#define LOS_RANDOM_SEED (-123456789) // seed for the extract LOS random number generator

#define INITIAL_REDSHIFT (float) 300 // used to perturb field

/************  END USER CHANGABLE DEFINITIONS  **************/
/************************************************************************************************/
/************************************************************************************************/

#define HII_EFF_FACTOR (float) ( N_GAMMA_UV * STELLAR_BARYON_FRAC * ESC_FRAC) 
#ifdef MINI_HALO
#define HII_EFF_FACTOR_MINI (float) ( N_GAMMA_UV_MINI * STELLAR_BARYON_FRAC_MINI * ESC_FRAC_MINI) 
#endif

/*
  Component of the velocity to be used in 21-cm temperature maps
  3 = z
*/
#define VELOCITY_COMPONENT (int) (3)

#define L_FACTOR (float) (0.620350491) // factor relating cube length to filter radius = (4PI/3)^(-1/3)

#define HII_D (unsigned long long) (HII_DIM)
#define HII_MIDDLE (HII_DIM/2)
#define HII_MID ((unsigned long long)HII_MIDDLE)

#define HII_TOT_NUM_PIXELS (unsigned long long)(HII_D*HII_D*HII_D)
#define HII_TOT_FFT_NUM_PIXELS ((unsigned long long)(HII_D*HII_D*2llu*(HII_MID+1llu)))
#define HII_KSPACE_NUM_PIXELS ((unsigned long long)(HII_D*HII_D*(HII_MID+1llu)))

/* INDEXING MACROS */
#define HII_C_INDEX(x,y,z)((unsigned long long)((z)+(HII_MID+1llu)*((y)+HII_D*(x))))// for 3D complex array
#define HII_R_FFT_INDEX(x,y,z)((unsigned long long)((z)+2llu*(HII_MID+1llu)*((y)+HII_D*(x)))) // for 3D real array with the FFT padding
#define HII_R_INDEX(x,y,z)((unsigned long long)((z)+HII_D*((y)+HII_D*(x)))) // for 3D real array with no padding

#endif
