#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

! We want to compute the surface integral of *F on the horizon.
! What we want to do is: Q = ∫_S *F
! *F is a 2-form: *F = *F_xy dx ∧ dy + *F_yz dy ∧ dz + *F_xz dx ∧ dz
! So we have to compute integrals of the form ∫_S *F_xy dx ∧ dy
!
! Let s be the shape of the horizon,
! S is parametrized with:
! S(θ, φ) = (x_0 + s(θ, φ) sin θ cos φ,
!            y_0 + s(θ, φ) sin θ sin φ,
!            z_0 + s(θ, φ) cos θ)
!
! The integral is ∫_θ ∫_φ *F_xy(θ, φ) det J[xy <-> θφ] dθ dφ
! With J[xy <-> θφ] Jacobian of the transformation
!
! What we have is (Aphi, A) in Cartesian coordinates, so the steps are:
! 1. Compute F in Cartesian coordinates F = dA, F_ij = A_[i,j]
! 2. Compute *F in Cartesian coordinates *F_ij = 1/2 eps_ijkl F^kl
! 3. Split *F in the various components
! --------------------------------
! That's not true anymore! With FaradayBase we have directly Ex, Ey, Ez
! 4. Compute the derivatives of  s(θ, φ) with respect to θ and φ
! 5. Compute the determinant of the transormations
!    (qlm_compute_jacobian)
! 6. Compute the integral

subroutine qlm_compute_jacobian (CCTK_ARGUMENTS, hn)
  use cctk
  use qlm_derivs
  use qlm_variables
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: hn

  CCTK_REAL :: theta, phi_an
  CCTK_REAL :: sin_theta, cos_theta, sin_phi, cos_phi
  CCTK_REAL :: sin_theta_cos_phi, sin_theta_sin_phi
  CCTK_REAL :: cos_theta_cos_phi, cos_theta_sin_phi
  CCTK_REAL :: deriv_shape(1,2)
  CCTK_REAL :: s, s_d_theta, s_d_phi ! s = shape
  CCTK_REAL :: Jxy_11, Jxy_12, Jxy_21, Jxy_22
  CCTK_REAL :: Jyz_11, Jyz_12, Jyz_21, Jyz_22
  CCTK_REAL :: Jxz_11, Jxz_12, Jxz_21, Jxz_22

  CCTK_REAL :: delta_space(2)

  integer :: i, j

  delta_space(:) = (/ qlm_delta_theta(hn), qlm_delta_phi(hn) /)

  if (verbose/=0) then
     call CCTK_INFO ("Computing the determinant of the transormations")
  end if

  ! Calculate the coordinates
  do j = 1, qlm_nphi(hn)
     do i = 1, qlm_ntheta(hn)
        theta  = qlm_origin_theta(hn) + (i-1)*qlm_delta_theta(hn)
        phi_an = qlm_origin_phi(hn)   + (j-1)*qlm_delta_phi(hn)

        sin_theta = sin(theta)
        cos_theta = cos(theta)
        sin_phi   = sin(phi_an)
        cos_phi   = cos(phi_an)

        sin_theta_cos_phi = sin_theta * cos_phi
        sin_theta_sin_phi = sin_theta * sin_phi
        cos_theta_cos_phi = cos_theta * cos_phi
        cos_theta_sin_phi = cos_theta * sin_phi

        s = qlm_shape(i,j,hn)

        deriv_shape(1,1:2) = deriv (qlm_shape(:,:,hn), i, j, delta_space)

        s_d_theta = deriv_shape(1,1)
        s_d_phi   = deriv_shape(1,2)

        Jxy_11 = s_d_theta * sin_theta_cos_phi + s * cos_theta_cos_phi
        Jxy_12 = s_d_phi   * sin_theta_cos_phi - s * sin_theta_sin_phi
        Jxy_21 = s_d_theta * sin_theta_sin_phi + s * cos_theta_sin_phi
        Jxy_22 = s_d_phi   * sin_theta_sin_phi + s * sin_theta_cos_phi

        Jyz_11 = s_d_theta * sin_theta_sin_phi + s * cos_theta_sin_phi
        Jyz_12 = s_d_phi   * sin_theta_sin_phi + s * sin_theta_cos_phi
        Jyz_21 = s_d_theta * cos_theta         - s * sin_theta
        Jyz_22 = s_d_phi   * cos_theta

        Jxz_11 = s_d_theta * sin_theta_cos_phi + s * cos_theta_cos_phi
        Jxz_12 = s_d_phi   * sin_theta_cos_phi - s * sin_theta_sin_phi
        Jxz_21 = s_d_theta * cos_theta         - s * sin_theta
        Jxz_22 = s_d_phi   * cos_theta

        qlm_det_Jxy(i,j,hn) = Jxy_11 * Jxy_22 - Jxy_12 * Jxy_21
        qlm_det_Jyz(i,j,hn) = Jyz_11 * Jyz_22 - Jyz_12 * Jyz_21
        qlm_det_Jxz(i,j,hn) = Jxz_11 * Jxz_22 - Jxz_12 * Jxz_21

     end do
  end do

end subroutine qlm_compute_jacobian

