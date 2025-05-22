#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"


subroutine qlm_symmetries( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT ierr

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "QuasiLocalMeasuresMHD::spatial_vol_in_qlm" )

end subroutine qlm_symmetries
