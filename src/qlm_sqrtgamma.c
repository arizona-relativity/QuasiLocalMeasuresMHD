#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <math.h>

void qlm_compute_sqrt_gamma(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /* TODO: This operation is massively inefficient!!!! It is performed at every
   iteration, but it is needed only when the horizon properties are computed.
   The optimal solution would be to call this routine only when it is needed. */

#pragma omp parallel for
  for(int index = 0; index < cctk_lsh[2] * cctk_lsh[1] * cctk_lsh[0]; index++){
        CCTK_REAL gxxL = gxx[index];
        CCTK_REAL gxyL = gxy[index];
        CCTK_REAL gxzL = gxz[index];
        CCTK_REAL gyyL = gyy[index];
        CCTK_REAL gyzL = gyz[index];
        CCTK_REAL gzzL = gzz[index];

        spatial_vol_in_qlm[index] = sqrt(-gxzL*gxzL*gyyL + 2*gxyL*gxzL*gyzL - gxxL*gyzL*gyzL - gxyL*gxyL*gzzL + gxxL*gyyL*gzzL);
      }
}
