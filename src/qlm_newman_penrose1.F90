#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine qlm_calc_newman_penrose1 (CCTK_ARGUMENTS, hn)
  use adm_metric
  use cctk
  use classify
  use constants
  use qlm_derivs
  use qlm_variables
  use ricci4
  use tensor
  use tensor4
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: hn

  CCTK_REAL, parameter :: zero=0, one=1, two=2, half=1/two
  CCTK_REAL    :: ll(0:3), nn(0:3)
  CCTK_COMPLEX :: mmm(0:3)
  CCTK_REAL    :: nabla_ll(0:3,0:3), nabla_nn(0:3,0:3)
  CCTK_COMPLEX :: nabla_mm(0:3,0:3)

  CCTK_REAL    :: t0, t1, t2
  logical      :: ce0, ce1, ce2
  CCTK_REAL    :: delta_space(2)

  integer      :: i, j
  integer      :: aa, bb
  CCTK_REAL    :: theta, phi_an

  if (veryverbose/=0) then
     call CCTK_INFO ("Calculating Newman-Penrose quantities")
  end if

  t0 = qlm_time(hn)
  t1 = qlm_time_p(hn)
  t2 = qlm_time_p_p(hn)

  ce0 = qlm_have_valid_data(hn) == 0
  ce1 = qlm_have_valid_data_p(hn) == 0
  ce2 = qlm_have_valid_data_p_p(hn) == 0

  delta_space(:) = (/ qlm_delta_theta(hn), qlm_delta_phi(hn) /)

  if (qlm_nghoststheta(hn)<2 .or. qlm_nghostsphi(hn)<2) call CCTK_WARN (0, "internal error")

  ! Calculate the coordinates
  do j = 1+qlm_nghostsphi(hn), qlm_nphi(hn)-qlm_nghostsphi(hn)
     do i = 1+qlm_nghoststheta(hn), qlm_ntheta(hn)-qlm_nghoststheta(hn)
        theta = qlm_origin_theta(hn) + (i-1)*qlm_delta_theta(hn)
        phi_an   = qlm_origin_phi(hn)   + (j-1)*qlm_delta_phi(hn)

        ! Get the stuff from the arrays
        ll(0) = qlm_l0(i,j,hn)
        ll(1) = qlm_l1(i,j,hn)
        ll(2) = qlm_l2(i,j,hn)
        ll(3) = qlm_l3(i,j,hn)

        nn(0) = qlm_n0(i,j,hn)
        nn(1) = qlm_n1(i,j,hn)
        nn(2) = qlm_n2(i,j,hn)
        nn(3) = qlm_n3(i,j,hn)

        mmm(0) = qlm_m0(i,j,hn)
        mmm(1) = qlm_m1(i,j,hn)
        mmm(2) = qlm_m2(i,j,hn)
        mmm(3) = qlm_m3(i,j,hn)

        nabla_ll(:,:) = qlm_tetrad_derivs(i,j)%nabla_ll(:,:)
        nabla_nn(:,:) = qlm_tetrad_derivs(i,j)%nabla_nn(:,:)
        nabla_mm(:,:) = qlm_tetrad_derivs(i,j)%nabla_mm(:,:)



        ! kappa   = - m^a  l^b D_b l_a
        ! tau     = - m^a  n^b D_b l_a
        ! sigma   = - m^a  m^b D_b l_a
        ! rho     = - m^a ~m^b D_b l_a   [ ~m = congj(m) ]
        qlm_npkappa  (i,j,hn) = 0
        qlm_nptau    (i,j,hn) = 0
        qlm_npsigma  (i,j,hn) = 0
        qlm_nprho    (i,j,hn) = 0
        ! epsilon = - 1/2 ( n^a  l^b D_b l_a - ~m^a  l^b D_b m_a )
        ! gamma   = - 1/2 ( n^a  n^b D_b l_a - ~m^a  n^b D_b m_a )
        ! beta    = - 1/2 ( n^a  m^b D_b l_a - ~m^a  m^b D_b m_a )
        ! alpha   = - 1/2 ( n^a ~m^b D_b l_a - ~m^a ~m^b D_b m_a )
        qlm_npepsilon(i,j,hn) = 0
        qlm_npgamma  (i,j,hn) = 0
        qlm_npbeta   (i,j,hn) = 0
        qlm_npalpha  (i,j,hn) = 0
        ! pi      = ~m^a  l^b D_b n_a
        ! nu      = ~m^a  n^b D_b n_a
        ! mu      = ~m^a  m^b D_b n_a
        ! lambda  = ~m^a ~m^b D_b n_a
        qlm_nppi     (i,j,hn) = 0
        qlm_npnu     (i,j,hn) = 0
        qlm_npmu     (i,j,hn) = 0
        qlm_nplambda (i,j,hn) = 0

!!$        qlm_lie_l_npsigma(i,j,hn) = 0
!!$        qlm_lie_n_npsigma(i,j,hn) = 0

!!$        ! Theta is the expansion
!!$        ! Theta_(l) is (- 2 Re rho)
!!$        ! Theta_(n) is (+ 2 Re mu)
!!$        ! Theta_(T) is [Theta_(l) + Theta_(n)] / sqrt(2)
!!$        ! Theta_(S) is [Theta_(l) - Theta_(n)] / sqrt(2)

!!$        qlm_lie_l_theta_l(i,j,hn) = 0
!!$        qlm_lie_l_theta_n(i,j,hn) = 0
!!$        qlm_lie_n_theta_l(i,j,hn) = 0
!!$        qlm_lie_n_theta_n(i,j,hn) = 0

        do aa=0,3
           do bb=0,3
              qlm_npkappa  (i,j,hn) = qlm_npkappa  (i,j,hn) - mmm(aa) *       ll(bb)  * nabla_ll(aa,bb)
              qlm_nptau    (i,j,hn) = qlm_nptau    (i,j,hn) - mmm(aa) *       nn(bb)  * nabla_ll(aa,bb)
              qlm_npsigma  (i,j,hn) = qlm_npsigma  (i,j,hn) - mmm(aa) *       mmm(bb)  * nabla_ll(aa,bb)
              qlm_nprho    (i,j,hn) = qlm_nprho    (i,j,hn) - mmm(aa) * conjg(mmm(bb)) * nabla_ll(aa,bb)
              qlm_npepsilon(i,j,hn) = qlm_npepsilon(i,j,hn) - half * (nn(aa) *       ll(bb)  * nabla_ll(aa,bb) - conjg(mmm(aa)) *       ll(bb)  * nabla_mm(aa,bb))
              qlm_npgamma  (i,j,hn) = qlm_npgamma  (i,j,hn) - half * (nn(aa) *       nn(bb)  * nabla_ll(aa,bb) - conjg(mmm(aa)) *       nn(bb)  * nabla_mm(aa,bb))
              qlm_npbeta   (i,j,hn) = qlm_npbeta   (i,j,hn) - half * (nn(aa) *       mmm(bb)  * nabla_ll(aa,bb) - conjg(mmm(aa)) *       mmm(bb)  * nabla_mm(aa,bb))
              qlm_npalpha  (i,j,hn) = qlm_npalpha  (i,j,hn) - half * (nn(aa) * conjg(mmm(bb)) * nabla_ll(aa,bb) - conjg(mmm(aa)) * conjg(mmm(bb)) * nabla_mm(aa,bb))
              qlm_nppi     (i,j,hn) = qlm_nppi     (i,j,hn) + conjg(mmm(aa)) *       ll(bb)  * nabla_nn(aa,bb)
              qlm_npnu     (i,j,hn) = qlm_npnu     (i,j,hn) + conjg(mmm(aa)) *       nn(bb)  * nabla_nn(aa,bb)
              qlm_npmu     (i,j,hn) = qlm_npmu     (i,j,hn) + conjg(mmm(aa)) *       mmm(bb)  * nabla_nn(aa,bb)
              qlm_nplambda (i,j,hn) = qlm_nplambda (i,j,hn) + conjg(mmm(aa)) * conjg(mmm(bb)) * nabla_nn(aa,bb)

           end do
        end do

     end do
  end do
end subroutine qlm_calc_newman_penrose1
