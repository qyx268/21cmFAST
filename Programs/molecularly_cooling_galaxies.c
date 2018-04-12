#include "../Parameter_files/ANAL_PARAMS.H"

/*
 functions that needs for calculate molecularly-cooled galaxies
 Yuxiang Qin (Yuxiang.L.Qin@gmail.com)
*/

float M_TURNOVER (float log10_M_feedback, float z_mol, float DELTA_z_mol, float z)
{
    float S = exp( (z - z_mol) / DELTA_z_mol );
    return pow(10, (log10_M_feedback + LOG10_M_MOL * S) / (1. + S));
}

