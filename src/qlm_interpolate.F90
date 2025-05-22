#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



! TODO:
! instead of interpolating to the symmetry points, copy them



! A convenient shortcut - changed the name of the pointer trying to fix a call error
#define QLM_PTR(x) CCTK_PointerTo(x)



subroutine qlm_interpolate (CCTK_ARGUMENTS, hn)
  use cctk
  use qlm_variables
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: hn

  CCTK_INT,  parameter :: izero = 0
  integer,   parameter :: ik = kind(izero)
  integer,   parameter :: sk = kind(interpolator)
  CCTK_REAL, parameter :: one = 1

  CCTK_REAL, parameter :: poison_value = -42

  integer      :: len_coordsystem
  integer      :: len_interpolator
  integer      :: len_interpolator_options

  character    :: fort_coordsystem*100
  character    :: fort_interpolator*100
  character    :: fort_interpolator_options*1000

  integer      :: nvars

  integer      :: coord_handle
  integer      :: interp_handle
  integer      :: options_table

  integer      :: ninputs
  integer      :: noutputs

  CCTK_REAL, allocatable :: xcoord(:,:)
  CCTK_REAL, allocatable :: ycoord(:,:)
  CCTK_REAL, allocatable :: zcoord(:,:)

  integer      :: ind_gxx, ind_gxy, ind_gxz, ind_gyy, ind_gyz, ind_gzz
  integer      :: ind_kxx, ind_kxy, ind_kxz, ind_kyy, ind_kyz, ind_kzz
  integer      :: ind_alpha
  integer      :: ind_spatial_vol_in_qlm
  integer      :: ind_betax, ind_betay, ind_betaz
  integer      :: ind_Bx, ind_By, ind_Bz
  integer      :: ind_ttt
  integer      :: ind_ttx, ind_tty, ind_ttz
  integer      :: ind_txx, ind_txy, ind_txz, ind_tyy, ind_tyz, ind_tzz

  integer      :: coord_type
  CCTK_POINTER :: coords(3)
  CCTK_INT     :: inputs(30)
  CCTK_INT     :: output_types(102)
  CCTK_POINTER :: outputs(102)
  CCTK_INT     :: operand_indices(102)
  CCTK_INT     :: operation_codes(102)
  integer      :: npoints

  character    :: msg*1000

  integer      :: ni, nj

  integer      :: ierr



  if (veryverbose/=0) then
     call CCTK_INFO ("Interpolating 3d grid functions")
  end if



  if (shift_state==0) then
     call CCTK_WARN (0, "The shift must have storage")
  end if

!!$  if (stress_energy_state==0) then
!!$     call CCTK_WARN (0, "The stress-energy tensor must have storage")
!!$  end if



  ! Get coordinate system
  call CCTK_FortranString &
       (len_coordsystem, int(coordsystem,sk), fort_coordsystem)
  call CCTK_CoordSystemHandle (coord_handle, fort_coordsystem)
  if (coord_handle<0) then
     write (msg, '("The coordinate system """, a, """ does not exist")') &
          trim(fort_coordsystem)
     call CCTK_WARN (0, msg)
  end if

  ! Get interpolator
  call CCTK_FortranString &
       (len_interpolator, int(interpolator,sk), fort_interpolator)
  call CCTK_InterpHandle (interp_handle, fort_interpolator)
  if (interp_handle<0) then
     write (msg, '("The interpolator """,a,""" does not exist")') &
          trim(fort_interpolator)
     call CCTK_WARN (0, msg)
  end if

  ! Get interpolator options
  call CCTK_FortranString &
       (len_interpolator_options, int(interpolator_options,sk), &
       fort_interpolator_options)
  call Util_TableCreateFromString (options_table, fort_interpolator_options)
  if (options_table<0) then
     write (msg, '("The interpolator_options """,a,""" have a wrong syntax")') &
          trim(fort_interpolator_options)
     call CCTK_WARN (0, msg)
  end if



  if (hn > 0) then

     ni = qlm_ntheta(hn)
     nj = qlm_nphi(hn)

     allocate (xcoord(ni,nj))
     allocate (ycoord(ni,nj))
     allocate (zcoord(ni,nj))

     xcoord(:,:) = qlm_x(:ni,:nj,hn)
     ycoord(:,:) = qlm_y(:ni,:nj,hn)
     zcoord(:,:) = qlm_z(:ni,:nj,hn)

  end if



  ! TODO: check the excision mask

  ! Get variable indices
  call CCTK_VarIndex (ind_gxx  , "ADMBase::gxx"   )
  call CCTK_VarIndex (ind_gxy  , "ADMBase::gxy"   )
  call CCTK_VarIndex (ind_gxz  , "ADMBase::gxz"   )
  call CCTK_VarIndex (ind_gyy  , "ADMBase::gyy"   )
  call CCTK_VarIndex (ind_gyz  , "ADMBase::gyz"   )
  call CCTK_VarIndex (ind_gzz  , "ADMBase::gzz"   )
  call CCTK_VarIndex (ind_kxx  , "ADMBase::kxx"   )
  call CCTK_VarIndex (ind_kxy  , "ADMBase::kxy"   )
  call CCTK_VarIndex (ind_kxz  , "ADMBase::kxz"   )
  call CCTK_VarIndex (ind_kyy  , "ADMBase::kyy"   )
  call CCTK_VarIndex (ind_kyz  , "ADMBase::kyz"   )
  call CCTK_VarIndex (ind_kzz  , "ADMBase::kzz"   )
  call CCTK_VarIndex (ind_alpha, "ADMBase::alp"   )
  call CCTK_VarIndex (ind_spatial_vol_in_qlm, "QuasiLocalMeasuresMHD::spatial_vol_in_qlm")
  call CCTK_VarIndex (ind_betax, "ADMBase::betax" )
  call CCTK_VarIndex (ind_betay, "ADMBase::betay" )
  call CCTK_VarIndex (ind_betaz, "ADMBase::betaz" )
  call CCTK_VarIndex (ind_Bx, "IllinoisGRMHD::Bx" )
  call CCTK_VarIndex (ind_By, "IllinoisGRMHD::By" )
  call CCTK_VarIndex (ind_Bz, "IllinoisGRMHD::Bz" )
  if (stress_energy_state /= 0) then
     call CCTK_VarIndex (ind_ttt  , "TmunuBase::eTtt")
     call CCTK_VarIndex (ind_ttx  , "TmunuBase::eTtx")
     call CCTK_VarIndex (ind_tty  , "TmunuBase::eTty")
     call CCTK_VarIndex (ind_ttz  , "TmunuBase::eTtz")
     call CCTK_VarIndex (ind_txx  , "TmunuBase::eTxx")
     call CCTK_VarIndex (ind_txy  , "TmunuBase::eTxy")
     call CCTK_VarIndex (ind_txz  , "TmunuBase::eTxz")
     call CCTK_VarIndex (ind_tyy  , "TmunuBase::eTyy")
     call CCTK_VarIndex (ind_tyz  , "TmunuBase::eTyz")
     call CCTK_VarIndex (ind_tzz  , "TmunuBase::eTzz")
  else
     ind_ttt = -1
     ind_ttx = -1
     ind_tty = -1
     ind_ttz = -1
     ind_txx = -1
     ind_txy = -1
     ind_txz = -1
     ind_tyy = -1
     ind_tyz = -1
     ind_tzz = -1
  end if



  ! Set up the interpolator arguments
  coord_type = CCTK_VARIABLE_REAL
  if (hn > 0) then
     npoints = ni * nj
     coords(:) = (/ QLM_PTR(xcoord), QLM_PTR(ycoord), QLM_PTR(zcoord) /)
  else
     npoints = 0
     coords(:) = CCTK_NullPointer()
  end if

  inputs = (/ &
       ind_gxx, ind_gxy, ind_gxz, ind_gyy, ind_gyz, ind_gzz, &
       ind_kxx, ind_kxy, ind_kxz, ind_kyy, ind_kyz, ind_kzz, &
       ind_alpha, &
       ind_spatial_vol_in_qlm, &
       ind_betax, ind_betay, ind_betaz, &
       ind_Bx, ind_By, ind_Bz, &
       ind_ttt, &
       ind_ttx, ind_tty, ind_ttz, &
       ind_txx, ind_txy, ind_txz, ind_tyy, ind_tyz, ind_tzz /)

  call CCTK_NumVars (nvars)
  if (nvars < 0) call CCTK_WARN (0, "internal error")
  if (any(inputs /= -1 .and. (inputs < 0 .or. inputs >= nvars))) then
     call CCTK_WARN (0, "internal error")
  end if

  operand_indices = (/ &
       00, 01, 02, 03, 04, 05, & ! g_ij
       00, 01, 02, 03, 04, 05, & ! g_ij,k
       00, 01, 02, 03, 04, 05, &
       00, 01, 02, 03, 04, 05, &
       00, 01, 02, 03, 04, 05, & ! g_ij,kl
       00, 01, 02, 03, 04, 05, &
       00, 01, 02, 03, 04, 05, &
       00, 01, 02, 03, 04, 05, &
       00, 01, 02, 03, 04, 05, &
       00, 01, 02, 03, 04, 05, &
       06, 07, 08, 09, 10, 11, & ! K_ij
       06, 07, 08, 09, 10, 11, & ! K_ij,k
       06, 07, 08, 09, 10, 11, &
       06, 07, 08, 09, 10, 11, &
       12, &                    ! alp
       13, &                    ! spatial_vol_in_qlm
       14, 15, 16, &            ! beta^i
       17, 18, 19, &            ! B^i
       20, &                    ! T_tt
       21, 22, 23, &            ! T_ti
       24, 25, 26, 27, 28, 29 /) ! T_ij

  operation_codes = (/ &
       0, 0, 0, 0, 0, 0, &      ! g_ij
       1, 1, 1, 1, 1, 1, &      ! g_ij,k
       2, 2, 2, 2, 2, 2, &
       3, 3, 3, 3, 3, 3, &
       11, 11, 11, 11, 11, 11, & ! g_ij,kl
       12, 12, 12, 12, 12, 12, &
       13, 13, 13, 13, 13, 13, &
       22, 22, 22, 22, 22, 22, &
       23, 23, 23, 23, 23, 23, &
       33, 33, 33, 33, 33, 33, &
       0, 0, 0, 0, 0, 0, &      ! K_ij
       1, 1, 1, 1, 1, 1, &      ! K_ij,k
       2, 2, 2, 2, 2, 2, &
       3, 3, 3, 3, 3, 3, &
       0, &                     ! alp
       0, &                     ! spatial_vol_in_qlm
       0, 0, 0, &               ! beta^i
       0, 0, 0, &               ! B^i
       0, &                     ! T_tt
       0, 0, 0, &               ! T_ti
       0, 0, 0, 0, 0, 0 /)      ! T_ij

  output_types(:) = CCTK_VARIABLE_REAL
  if (hn > 0) then
     outputs = (/ &
          QLM_PTR(qlm_gxx), QLM_PTR(qlm_gxy), QLM_PTR(qlm_gxz), QLM_PTR(qlm_gyy), QLM_PTR(qlm_gyz), QLM_PTR(qlm_gzz), &
          QLM_PTR(qlm_dgxxx), QLM_PTR(qlm_dgxyx), QLM_PTR(qlm_dgxzx), QLM_PTR(qlm_dgyyx), QLM_PTR(qlm_dgyzx), QLM_PTR(qlm_dgzzx), &
          QLM_PTR(qlm_dgxxy), QLM_PTR(qlm_dgxyy), QLM_PTR(qlm_dgxzy), QLM_PTR(qlm_dgyyy), QLM_PTR(qlm_dgyzy), QLM_PTR(qlm_dgzzy), &
          QLM_PTR(qlm_dgxxz), QLM_PTR(qlm_dgxyz), QLM_PTR(qlm_dgxzz), QLM_PTR(qlm_dgyyz), QLM_PTR(qlm_dgyzz), QLM_PTR(qlm_dgzzz), &
          QLM_PTR(qlm_ddgxxxx), QLM_PTR(qlm_ddgxyxx), QLM_PTR(qlm_ddgxzxx), QLM_PTR(qlm_ddgyyxx), QLM_PTR(qlm_ddgyzxx), QLM_PTR(qlm_ddgzzxx), &
          QLM_PTR(qlm_ddgxxxy), QLM_PTR(qlm_ddgxyxy), QLM_PTR(qlm_ddgxzxy), QLM_PTR(qlm_ddgyyxy), QLM_PTR(qlm_ddgyzxy), QLM_PTR(qlm_ddgzzxy), &
          QLM_PTR(qlm_ddgxxxz), QLM_PTR(qlm_ddgxyxz), QLM_PTR(qlm_ddgxzxz), QLM_PTR(qlm_ddgyyxz), QLM_PTR(qlm_ddgyzxz), QLM_PTR(qlm_ddgzzxz), &
          QLM_PTR(qlm_ddgxxyy), QLM_PTR(qlm_ddgxyyy), QLM_PTR(qlm_ddgxzyy), QLM_PTR(qlm_ddgyyyy), QLM_PTR(qlm_ddgyzyy), QLM_PTR(qlm_ddgzzyy), &
          QLM_PTR(qlm_ddgxxyz), QLM_PTR(qlm_ddgxyyz), QLM_PTR(qlm_ddgxzyz), QLM_PTR(qlm_ddgyyyz), QLM_PTR(qlm_ddgyzyz), QLM_PTR(qlm_ddgzzyz), &
          QLM_PTR(qlm_ddgxxzz), QLM_PTR(qlm_ddgxyzz), QLM_PTR(qlm_ddgxzzz), QLM_PTR(qlm_ddgyyzz), QLM_PTR(qlm_ddgyzzz), QLM_PTR(qlm_ddgzzzz), &
          QLM_PTR(qlm_kxx), QLM_PTR(qlm_kxy), QLM_PTR(qlm_kxz), QLM_PTR(qlm_kyy), QLM_PTR(qlm_kyz), QLM_PTR(qlm_kzz), &
          QLM_PTR(qlm_dkxxx), QLM_PTR(qlm_dkxyx), QLM_PTR(qlm_dkxzx), QLM_PTR(qlm_dkyyx), QLM_PTR(qlm_dkyzx), QLM_PTR(qlm_dkzzx), &
          QLM_PTR(qlm_dkxxy), QLM_PTR(qlm_dkxyy), QLM_PTR(qlm_dkxzy), QLM_PTR(qlm_dkyyy), QLM_PTR(qlm_dkyzy), QLM_PTR(qlm_dkzzy), &
          QLM_PTR(qlm_dkxxz), QLM_PTR(qlm_dkxyz), QLM_PTR(qlm_dkxzz), QLM_PTR(qlm_dkyyz), QLM_PTR(qlm_dkyzz), QLM_PTR(qlm_dkzzz), &
          QLM_PTR(qlm_alpha), &
          QLM_PTR(qlm_spatial_vol_in_qlm), &
          QLM_PTR(qlm_betax), QLM_PTR(qlm_betay), QLM_PTR(qlm_betaz), &
          QLM_PTR(qlm_Bx), QLM_PTR(qlm_By), QLM_PTR(qlm_Bz), &
          QLM_PTR(qlm_ttt), &
          QLM_PTR(qlm_ttx), QLM_PTR(qlm_tty), QLM_PTR(qlm_ttz), &
          QLM_PTR(qlm_txx), QLM_PTR(qlm_txy), QLM_PTR(qlm_txz), QLM_PTR(qlm_tyy), QLM_PTR(qlm_tyz), QLM_PTR(qlm_tzz) /)
  else
     outputs(:) = CCTK_NullPointer()
  end if



  ninputs = size(inputs)
  noutputs = size(outputs)



#if 0
  ! Poison the output variables
  call poison (qlm_gxx    )
  call poison (qlm_gxy    )
  call poison (qlm_gxz    )
  call poison (qlm_gyy    )
  call poison (qlm_gyz    )
  call poison (qlm_gzz    )
  call poison (qlm_dgxxx  )
  call poison (qlm_dgxyx  )
  call poison (qlm_dgxzx  )
  call poison (qlm_dgyyx  )
  call poison (qlm_dgyzx  )
  call poison (qlm_dgzzx  )
  call poison (qlm_dgxxy  )
  call poison (qlm_dgxyy  )
  call poison (qlm_dgxzy  )
  call poison (qlm_dgyyy  )
  call poison (qlm_dgyzy  )
  call poison (qlm_dgzzy  )
  call poison (qlm_dgxxz  )
  call poison (qlm_dgxyz  )
  call poison (qlm_dgxzz  )
  call poison (qlm_dgyyz  )
  call poison (qlm_dgyzz  )
  call poison (qlm_dgzzz  )
  call poison (qlm_ddgxxxx)
  call poison (qlm_ddgxyxx)
  call poison (qlm_ddgxzxx)
  call poison (qlm_ddgyyxx)
  call poison (qlm_ddgyzxx)
  call poison (qlm_ddgzzxx)
  call poison (qlm_ddgxxxy)
  call poison (qlm_ddgxyxy)
  call poison (qlm_ddgxzxy)
  call poison (qlm_ddgyyxy)
  call poison (qlm_ddgyzxy)
  call poison (qlm_ddgzzxy)
  call poison (qlm_ddgxxxz)
  call poison (qlm_ddgxyxz)
  call poison (qlm_ddgxzxz)
  call poison (qlm_ddgyyxz)
  call poison (qlm_ddgyzxz)
  call poison (qlm_ddgzzxz)
  call poison (qlm_ddgxxyy)
  call poison (qlm_ddgxyyy)
  call poison (qlm_ddgxzyy)
  call poison (qlm_ddgyyyy)
  call poison (qlm_ddgyzyy)
  call poison (qlm_ddgzzyy)
  call poison (qlm_ddgxxyz)
  call poison (qlm_ddgxyyz)
  call poison (qlm_ddgxzyz)
  call poison (qlm_ddgyyyz)
  call poison (qlm_ddgyzyz)
  call poison (qlm_ddgzzyz)
  call poison (qlm_ddgxxzz)
  call poison (qlm_ddgxyzz)
  call poison (qlm_ddgxzzz)
  call poison (qlm_ddgyyzz)
  call poison (qlm_ddgyzzz)
  call poison (qlm_ddgzzzz)
  call poison (qlm_kxx    )
  call poison (qlm_kxy    )
  call poison (qlm_kxz    )
  call poison (qlm_kyy    )
  call poison (qlm_kyz    )
  call poison (qlm_kzz    )
  call poison (qlm_dkxxx  )
  call poison (qlm_dkxyx  )
  call poison (qlm_dkxzx  )
  call poison (qlm_dkyyx  )
  call poison (qlm_dkyzx  )
  call poison (qlm_dkzzx  )
  call poison (qlm_dkxxy  )
  call poison (qlm_dkxyy  )
  call poison (qlm_dkxzy  )
  call poison (qlm_dkyyy  )
  call poison (qlm_dkyzy  )
  call poison (qlm_dkzzy  )
  call poison (qlm_dkxxz  )
  call poison (qlm_dkxyz  )
  call poison (qlm_dkxzz  )
  call poison (qlm_dkyyz  )
  call poison (qlm_dkyzz  )
  call poison (qlm_dkzzz  )
  call poison (qlm_alpha  )
  call poison (qlm_spatial_vol_in_qlm )
  call poison (qlm_betax  )
  call poison (qlm_betay  )
  call poison (qlm_betaz  )
  call poison (qlm_Bx  )
  call poison (qlm_By  )
  call poison (qlm_Bz  )
  call poison (qlm_ttt    )
  call poison (qlm_ttx    )
  call poison (qlm_tty    )
  call poison (qlm_ttz    )
  call poison (qlm_txx    )
  call poison (qlm_txy    )
  call poison (qlm_txz    )
  call poison (qlm_tyy    )
  call poison (qlm_tyz    )
  call poison (qlm_tzz    )
#endif



  ! Call the interpolator
  call Util_TableSetIntArray &
       (ierr, options_table, noutputs, &
       operand_indices, "operand_indices")
  if (ierr /= 0) call CCTK_WARN (0, "internal error")
  call Util_TableSetIntArray &
       (ierr, options_table, noutputs, &
       operation_codes, "operation_codes")
  if (ierr /= 0) call CCTK_WARN (0, "internal error")

  call CCTK_InterpGridArrays &
       (ierr, cctkGH, 3, &
       interp_handle, options_table, coord_handle, &
       npoints, coord_type, coords, &
       ninputs, inputs, &
       noutputs, output_types, outputs)

  if (ierr /= 0) then
     if (hn > 0) then
        qlm_calc_error(hn) = 1
     end if
     call CCTK_WARN (1, "Interpolator failed")
     return
  end if



  ! Unpack the variables
  if (hn > 0) then

     call unpack (qlm_gxx    , ni, nj)
     call unpack (qlm_gxy    , ni, nj)
     call unpack (qlm_gxz    , ni, nj)
     call unpack (qlm_gyy    , ni, nj)
     call unpack (qlm_gyz    , ni, nj)
     call unpack (qlm_gzz    , ni, nj)
     call unpack (qlm_dgxxx  , ni, nj)
     call unpack (qlm_dgxyx  , ni, nj)
     call unpack (qlm_dgxzx  , ni, nj)
     call unpack (qlm_dgyyx  , ni, nj)
     call unpack (qlm_dgyzx  , ni, nj)
     call unpack (qlm_dgzzx  , ni, nj)
     call unpack (qlm_dgxxy  , ni, nj)
     call unpack (qlm_dgxyy  , ni, nj)
     call unpack (qlm_dgxzy  , ni, nj)
     call unpack (qlm_dgyyy  , ni, nj)
     call unpack (qlm_dgyzy  , ni, nj)
     call unpack (qlm_dgzzy  , ni, nj)
     call unpack (qlm_dgxxz  , ni, nj)
     call unpack (qlm_dgxyz  , ni, nj)
     call unpack (qlm_dgxzz  , ni, nj)
     call unpack (qlm_dgyyz  , ni, nj)
     call unpack (qlm_dgyzz  , ni, nj)
     call unpack (qlm_dgzzz  , ni, nj)
     call unpack (qlm_ddgxxxx, ni, nj)
     call unpack (qlm_ddgxyxx, ni, nj)
     call unpack (qlm_ddgxzxx, ni, nj)
     call unpack (qlm_ddgyyxx, ni, nj)
     call unpack (qlm_ddgyzxx, ni, nj)
     call unpack (qlm_ddgzzxx, ni, nj)
     call unpack (qlm_ddgxxxy, ni, nj)
     call unpack (qlm_ddgxyxy, ni, nj)
     call unpack (qlm_ddgxzxy, ni, nj)
     call unpack (qlm_ddgyyxy, ni, nj)
     call unpack (qlm_ddgyzxy, ni, nj)
     call unpack (qlm_ddgzzxy, ni, nj)
     call unpack (qlm_ddgxxxz, ni, nj)
     call unpack (qlm_ddgxyxz, ni, nj)
     call unpack (qlm_ddgxzxz, ni, nj)
     call unpack (qlm_ddgyyxz, ni, nj)
     call unpack (qlm_ddgyzxz, ni, nj)
     call unpack (qlm_ddgzzxz, ni, nj)
     call unpack (qlm_ddgxxyy, ni, nj)
     call unpack (qlm_ddgxyyy, ni, nj)
     call unpack (qlm_ddgxzyy, ni, nj)
     call unpack (qlm_ddgyyyy, ni, nj)
     call unpack (qlm_ddgyzyy, ni, nj)
     call unpack (qlm_ddgzzyy, ni, nj)
     call unpack (qlm_ddgxxyz, ni, nj)
     call unpack (qlm_ddgxyyz, ni, nj)
     call unpack (qlm_ddgxzyz, ni, nj)
     call unpack (qlm_ddgyyyz, ni, nj)
     call unpack (qlm_ddgyzyz, ni, nj)
     call unpack (qlm_ddgzzyz, ni, nj)
     call unpack (qlm_ddgxxzz, ni, nj)
     call unpack (qlm_ddgxyzz, ni, nj)
     call unpack (qlm_ddgxzzz, ni, nj)
     call unpack (qlm_ddgyyzz, ni, nj)
     call unpack (qlm_ddgyzzz, ni, nj)
     call unpack (qlm_ddgzzzz, ni, nj)
     call unpack (qlm_kxx    , ni, nj)
     call unpack (qlm_kxy    , ni, nj)
     call unpack (qlm_kxz    , ni, nj)
     call unpack (qlm_kyy    , ni, nj)
     call unpack (qlm_kyz    , ni, nj)
     call unpack (qlm_kzz    , ni, nj)
     call unpack (qlm_dkxxx  , ni, nj)
     call unpack (qlm_dkxyx  , ni, nj)
     call unpack (qlm_dkxzx  , ni, nj)
     call unpack (qlm_dkyyx  , ni, nj)
     call unpack (qlm_dkyzx  , ni, nj)
     call unpack (qlm_dkzzx  , ni, nj)
     call unpack (qlm_dkxxy  , ni, nj)
     call unpack (qlm_dkxyy  , ni, nj)
     call unpack (qlm_dkxzy  , ni, nj)
     call unpack (qlm_dkyyy  , ni, nj)
     call unpack (qlm_dkyzy  , ni, nj)
     call unpack (qlm_dkzzy  , ni, nj)
     call unpack (qlm_dkxxz  , ni, nj)
     call unpack (qlm_dkxyz  , ni, nj)
     call unpack (qlm_dkxzz  , ni, nj)
     call unpack (qlm_dkyyz  , ni, nj)
     call unpack (qlm_dkyzz  , ni, nj)
     call unpack (qlm_dkzzz  , ni, nj)
     call unpack (qlm_alpha  , ni, nj)
     call unpack (qlm_spatial_vol_in_qlm , ni, nj)
     call unpack (qlm_betax  , ni, nj)
     call unpack (qlm_betay  , ni, nj)
     call unpack (qlm_betaz  , ni, nj)
     call unpack (qlm_Bx  , ni, nj)
     call unpack (qlm_By  , ni, nj)
     call unpack (qlm_Bz  , ni, nj)
     if (stress_energy_state /= 0) then
        call unpack (qlm_ttt    , ni, nj)
        call unpack (qlm_ttx    , ni, nj)
        call unpack (qlm_tty    , ni, nj)
        call unpack (qlm_ttz    , ni, nj)
        call unpack (qlm_txx    , ni, nj)
        call unpack (qlm_txy    , ni, nj)
        call unpack (qlm_txz    , ni, nj)
        call unpack (qlm_tyy    , ni, nj)
        call unpack (qlm_tyz    , ni, nj)
        call unpack (qlm_tzz    , ni, nj)
     else
        qlm_ttt = 0
        qlm_ttx = 0
        qlm_tty = 0
        qlm_ttz = 0
        qlm_txx = 0
        qlm_txy = 0
        qlm_txz = 0
        qlm_tyy = 0
        qlm_tyz = 0
        qlm_tzz = 0
     end if



#if 0
     ! Check for poison
     call poison_check (qlm_gxx    , "qlm_gxx    ")
     call poison_check (qlm_gxy    , "qlm_gxy    ")
     call poison_check (qlm_gxz    , "qlm_gxz    ")
     call poison_check (qlm_gyy    , "qlm_gyy    ")
     call poison_check (qlm_gyz    , "qlm_gyz    ")
     call poison_check (qlm_gzz    , "qlm_gzz    ")
     call poison_check (qlm_dgxxx  , "qlm_dgxxx  ")
     call poison_check (qlm_dgxyx  , "qlm_dgxyx  ")
     call poison_check (qlm_dgxzx  , "qlm_dgxzx  ")
     call poison_check (qlm_dgyyx  , "qlm_dgyyx  ")
     call poison_check (qlm_dgyzx  , "qlm_dgyzx  ")
     call poison_check (qlm_dgzzx  , "qlm_dgzzx  ")
     call poison_check (qlm_dgxxy  , "qlm_dgxxy  ")
     call poison_check (qlm_dgxyy  , "qlm_dgxyy  ")
     call poison_check (qlm_dgxzy  , "qlm_dgxzy  ")
     call poison_check (qlm_dgyyy  , "qlm_dgyyy  ")
     call poison_check (qlm_dgyzy  , "qlm_dgyzy  ")
     call poison_check (qlm_dgzzy  , "qlm_dgzzy  ")
     call poison_check (qlm_dgxxz  , "qlm_dgxxz  ")
     call poison_check (qlm_dgxyz  , "qlm_dgxyz  ")
     call poison_check (qlm_dgxzz  , "qlm_dgxzz  ")
     call poison_check (qlm_dgyyz  , "qlm_dgyyz  ")
     call poison_check (qlm_dgyzz  , "qlm_dgyzz  ")
     call poison_check (qlm_dgzzz  , "qlm_dgzzz  ")
     call poison_check (qlm_ddgxxxx, "qlm_ddgxxxx")
     call poison_check (qlm_ddgxyxx, "qlm_ddgxyxx")
     call poison_check (qlm_ddgxzxx, "qlm_ddgxzxx")
     call poison_check (qlm_ddgyyxx, "qlm_ddgyyxx")
     call poison_check (qlm_ddgyzxx, "qlm_ddgyzxx")
     call poison_check (qlm_ddgzzxx, "qlm_ddgzzxx")
     call poison_check (qlm_ddgxxxy, "qlm_ddgxxxy")
     call poison_check (qlm_ddgxyxy, "qlm_ddgxyxy")
     call poison_check (qlm_ddgxzxy, "qlm_ddgxzxy")
     call poison_check (qlm_ddgyyxy, "qlm_ddgyyxy")
     call poison_check (qlm_ddgyzxy, "qlm_ddgyzxy")
     call poison_check (qlm_ddgzzxy, "qlm_ddgzzxy")
     call poison_check (qlm_ddgxxxz, "qlm_ddgxxxz")
     call poison_check (qlm_ddgxyxz, "qlm_ddgxyxz")
     call poison_check (qlm_ddgxzxz, "qlm_ddgxzxz")
     call poison_check (qlm_ddgyyxz, "qlm_ddgyyxz")
     call poison_check (qlm_ddgyzxz, "qlm_ddgyzxz")
     call poison_check (qlm_ddgzzxz, "qlm_ddgzzxz")
     call poison_check (qlm_ddgxxyy, "qlm_ddgxxyy")
     call poison_check (qlm_ddgxyyy, "qlm_ddgxyyy")
     call poison_check (qlm_ddgxzyy, "qlm_ddgxzyy")
     call poison_check (qlm_ddgyyyy, "qlm_ddgyyyy")
     call poison_check (qlm_ddgyzyy, "qlm_ddgyzyy")
     call poison_check (qlm_ddgzzyy, "qlm_ddgzzyy")
     call poison_check (qlm_ddgxxyz, "qlm_ddgxxyz")
     call poison_check (qlm_ddgxyyz, "qlm_ddgxyyz")
     call poison_check (qlm_ddgxzyz, "qlm_ddgxzyz")
     call poison_check (qlm_ddgyyyz, "qlm_ddgyyyz")
     call poison_check (qlm_ddgyzyz, "qlm_ddgyzyz")
     call poison_check (qlm_ddgzzyz, "qlm_ddgzzyz")
     call poison_check (qlm_ddgxxzz, "qlm_ddgxxzz")
     call poison_check (qlm_ddgxyzz, "qlm_ddgxyzz")
     call poison_check (qlm_ddgxzzz, "qlm_ddgxzzz")
     call poison_check (qlm_ddgyyzz, "qlm_ddgyyzz")
     call poison_check (qlm_ddgyzzz, "qlm_ddgyzzz")
     call poison_check (qlm_ddgzzzz, "qlm_ddgzzzz")
     call poison_check (qlm_kxx    , "qlm_kxx    ")
     call poison_check (qlm_kxy    , "qlm_kxy    ")
     call poison_check (qlm_kxz    , "qlm_kxz    ")
     call poison_check (qlm_kyy    , "qlm_kyy    ")
     call poison_check (qlm_kyz    , "qlm_kyz    ")
     call poison_check (qlm_kzz    , "qlm_kzz    ")
     call poison_check (qlm_dkxxx  , "qlm_dkxxx  ")
     call poison_check (qlm_dkxyx  , "qlm_dkxyx  ")
     call poison_check (qlm_dkxzx  , "qlm_dkxzx  ")
     call poison_check (qlm_dkyyx  , "qlm_dkyyx  ")
     call poison_check (qlm_dkyzx  , "qlm_dkyzx  ")
     call poison_check (qlm_dkzzx  , "qlm_dkzzx  ")
     call poison_check (qlm_dkxxy  , "qlm_dkxxy  ")
     call poison_check (qlm_dkxyy  , "qlm_dkxyy  ")
     call poison_check (qlm_dkxzy  , "qlm_dkxzy  ")
     call poison_check (qlm_dkyyy  , "qlm_dkyyy  ")
     call poison_check (qlm_dkyzy  , "qlm_dkyzy  ")
     call poison_check (qlm_dkzzy  , "qlm_dkzzy  ")
     call poison_check (qlm_dkxxz  , "qlm_dkxxz  ")
     call poison_check (qlm_dkxyz  , "qlm_dkxyz  ")
     call poison_check (qlm_dkxzz  , "qlm_dkxzz  ")
     call poison_check (qlm_dkyyz  , "qlm_dkyyz  ")
     call poison_check (qlm_dkyzz  , "qlm_dkyzz  ")
     call poison_check (qlm_dkzzz  , "qlm_dkzzz  ")
     call poison_check (qlm_alpha  , "qlm_alpha  ")
     call poison_check (qlm_spatial_vol_in_qlm , "qlm_spatial_vol_in_qlm ")
     call poison_check (qlm_betax  , "qlm_betax  ")
     call poison_check (qlm_betay  , "qlm_betay  ")
     call poison_check (qlm_betaz  , "qlm_betaz  ")
     call poison_check (qlm_Bx     , "qlm_Bx     ")
     call poison_check (qlm_By     , "qlm_By     ")
     call poison_check (qlm_Bz     , "qlm_Bz     ")
     call poison_check (qlm_ttt    , "qlm_ttt    ")
     call poison_check (qlm_ttx    , "qlm_ttx    ")
     call poison_check (qlm_tty    , "qlm_tty    ")
     call poison_check (qlm_ttz    , "qlm_ttz    ")
     call poison_check (qlm_txx    , "qlm_txx    ")
     call poison_check (qlm_txy    , "qlm_txy    ")
     call poison_check (qlm_txz    , "qlm_txz    ")
     call poison_check (qlm_tyy    , "qlm_tyy    ")
     call poison_check (qlm_tyz    , "qlm_tyz    ")
     call poison_check (qlm_tzz    , "qlm_tzz    ")
#endif

  end if



  ! Free interpolator options
  call Util_TableDestroy (ierr, options_table)



  if (hn > 0) then

     qlm_have_valid_data(hn) = 1

     deallocate (xcoord)
     deallocate (ycoord)
     deallocate (zcoord)

  end if



contains

  subroutine pack (arr, ni, nj)
    integer,   intent(in)    :: ni, nj
    CCTK_REAL, intent(inout) :: arr(:,:)
    CCTK_REAL :: tmp(ni,nj)
    tmp(:,:) = arr(:ni, :nj)
    call copy (arr, tmp, size(tmp))
  end subroutine pack

  subroutine unpack (arr, ni, nj)
    integer,   intent(in)    :: ni, nj
    CCTK_REAL, intent(inout) :: arr(:,:)
    CCTK_REAL :: tmp(ni,nj)
    call copy (tmp, arr, size(tmp))
    arr(:ni, :nj) = tmp(:,:)
    arr(ni+1:, :nj) = 0
    arr(:, nj+1:) = 0
  end subroutine unpack

  subroutine copy (a, b, n)
    integer,   intent(in)  :: n
    CCTK_REAL, intent(out) :: a(n)
    CCTK_REAL, intent(in)  :: b(n)
    a = b
  end subroutine copy

  subroutine poison (arr)
    CCTK_REAL, intent(out) :: arr(:,:)
    arr = poison_value
  end subroutine poison

  subroutine poison_check (arr, name)
    CCTK_REAL,    intent(in) :: arr(:,:)
    character(*), intent(in) :: name
    character*1000 :: msg
!!$    integer        :: i, j
    if (any(arr==poison_value)) then
       write (msg, '("Poison found in ",a)') trim(name)
       call CCTK_WARN (CCTK_WARN_ALERT, msg)
!!$       do j=1,size(arr,2)
!!$          do i=1,size(arr,1)
!!$             print '(2i6)', i,j
!!$          end do
!!$       end do
    end if
  end subroutine poison_check

end subroutine qlm_interpolate
