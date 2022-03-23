module MultiFluidDE
  use DarkEnergyInterface
  use results
  use constants
  use classes
  implicit none

  private
  ! JVR - can we use these quantities from results.f90 without needing to bring them to this file via the Init subroutine?
  real(dl) :: grhocrit
  real(dl) :: grhode_today
  real(dl) :: Omega_de

  type, extends(TDarkEnergyModel) :: TMultiFluidDE
    real(dl), dimension(max_num_of_fluids * max_num_of_params) :: de_params
    integer, dimension(max_num_of_fluids) :: models
    ! JVR - TODO go back to a parameter array
    real(dl) :: w0, wa
    real(dl) :: zc, fde_zc, wn
    real(dl) :: grhonode_zc
  contains
  procedure :: Init => TMultiFluidDE_Init
  procedure :: w_de => TMultiFluidDE_w_de
  procedure :: grho_de => TMultiFluidDE_grho_de
  procedure :: BackgroundDensityAndPressure => TMultiFluidDE_density
  procedure :: Effective_w_wa => TMultiFluidDE_Effective_w_wa
  procedure :: PerturbedStressEnergy => TMultiFluidDE_PerturbedStressEnergy
  procedure :: PerturbationEvolve => TMultiFluidDE_PerturbationEvolve
  procedure :: diff_rhopi_Add_Term => TMultiFluidDE_diff_rhopi_Add_Term
  procedure, nopass :: PythonClass => TMultiFluidDE_PythonClass
  procedure, nopass :: SelfPointer => TMultiFluidDE_SelfPointer
  procedure :: ReadParams => TMultiFluidDE_ReadParams
  procedure :: w_derivative
  end type TMultiFluidDE

  public TMultiFluidDE

  contains

  ! --------------------- Background Equations ---------------------

  function TMultiFluidDE_w_de(this, a) result (w_de)
    class(TMultiFluidDE) :: this
    real(dl), dimension(max_num_of_fluids) :: w_de
    real(dl), intent(in) :: a
    integer :: i
    real(dl) :: zc, ac, fde_zc, wn, z

    ! i - counts the component: 1 is the late component, 2 is the early component.
    w_de = 0
    do i = 1, this%num_of_components
      if (i == 1) then
        if (this%models(i) == 1) then
          ! w_constant model
          w_de(i) = this%w0
        else if (this%models(i) == 2) then
          ! w0wa with PPF - allows phantom divide crossing
          w_de(i) = this%w0 + this%wa * (1-a)
        else
          stop "[Multifluid DE] Invalid dark energy model for fluid 1"
        end if

      else if (i == 2) then
        if (this%models(i) == 1) then
          ! Axion effective fluid from arXiv:1806.10608
          zc = this%zc
          ac = 1._dl/(1._dl + zc)
          fde_zc = this%fde_zc
          wn = this%wn
          w_de(i) = (1._dl + wn)/(1._dl + (ac/a)**(3._dl*(1._dl + wn))) - 1._dl
        else
          stop "[Multifluid DE] Invalid dark energy model for fluid 2"
        end if
      else 
        stop "[Multifluid DE] Invalid dark energy fluid"
      end if
    end do
  end function TMultiFluidDE_w_de

  function TMultiFluidDE_grho_de(this, a) result(grho_de)
    class(TMultiFluidDE) :: this
    real(dl), dimension(max_num_of_fluids) :: grho_de
    real(dl), intent(in) :: a
    integer :: i
    real(dl) :: w0, wa, zc, ac, fde_zc, oma_zc, wn, z

    ! Returns 8*pi*G*rho_de, no scale factors
    ! Assuming EDE is negligible today
    ! grhode_today = grhocrit * (1 - Omega_baryons - Omega_cdm - Omega_neutrinos - Omega_photons)
    grho_de = 0
    do i = 1, this%num_of_components
      if (i == 1) then
        if (this%models(i) == 1) then
          ! w_constant model
          !w0 = this%de_params(max_num_of_params*(i-1) + 1)
          w0 = this%w0
          grho_de(i) = grhode_today * a**(-3*(1 + w0))
        else if (this%models(i) == 2) then
          w0 = this%w0
          wa = this%wa
          grho_de(i) = grhode_today * a**(-3*(1 + w0 + wa)) * exp(-3 * wa * (1-a))
        else
          stop "[Multifluid DE] Invalid dark energy model for fluid 1"
        end if

      else if (i == 2) then
        if (this%models(i) == 1) then
          ! Axion effective fluid from arXiv:1806.10608
          z = 1/a - 1
          !zc = this%de_params(max_num_of_params*(i-1) + 1)
          !fde_zc = this%de_params(max_num_of_params*(i-1) + 2)
          !wn = this%de_params(max_num_of_params*(i-1) + 3)
          zc = this%zc
          fde_zc = this%fde_zc
          wn = this%wn
          ac = 1._dl/(1._dl + zc)
          oma_zc = fde_zc * this%grhonode_zc / (grhocrit * (1 - fde_zc))
          grho_de(i) = grhocrit * (2._dl * oma_zc)/(1 + ( (1+zc)/(1+z) )**(3*(1+wn)))
        else
          stop "[Multifluid DE] Invalid dark energy model for fluid 2"
        end if
      else 
        stop "[Multifluid DE] Invalid dark energy fluid"
      end if
    end do
  end function TMultiFluidDE_grho_de

  subroutine TMultiFluidDE_density(this, grhov, a, grhov_t, w)
    ! Get grhov_t = 8*pi*G*rho_de*a**2 and (optionally) equation of state at scale factor a
    class(TMultiFluidDE), intent(inout) :: this
    real(dl), intent(in) :: grhov, a
    real(dl), dimension(max_num_of_fluids), intent(out) :: grhov_t
    real(dl), optional, dimension(max_num_of_fluids), intent(out) :: w
    integer :: i

    if (a > 1e-10) then
      grhov_t = this%grho_de(a) * a**2
    else
      grhov_t = 0
    end if
    if (present(w)) then
      w = this%w_de(a)
    end if
  end subroutine TMultiFluidDE_density

  ! --------------------- Perturbation Equations ---------------------

  subroutine TMultiFluidDE_PerturbedStressEnergy(this, dgrhoe, dgqe, &
    a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
    class(TMultiFluidDE), intent(inout) :: this
    real(dl), dimension(max_num_of_fluids), intent(out) :: dgrhoe, dgqe
    real(dl), dimension(max_num_of_fluids), intent(in) :: grhov_t
    real(dl), dimension(max_num_of_fluids), intent(in) :: w
    real(dl), intent(in) ::  a, dgq, dgrho, grho, gpres_noDE, etak, adotoa, k, kf1
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: w_ix
    integer :: i, delta_index

    ! TODO clean variables
    dgrhoe = 0
    dgqe = 0
    do i=1, this%num_of_components
      delta_index = w_ix + 2*(i - 1)
      if (i == 1 .and. this%models(1) == 2) then ! w0wa requires PPF
        call PPF_Perturbations(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, &
        etak, adotoa, k, kf1, ay, ayprime, w_ix)
      else
        if (this%models(1) == 2) then
          delta_index = delta_index - 1 ! PPF has only one perturbation equation
        end if
        ! Usual fluid definitions
        dgrhoe(i) = ay(delta_index) * grhov_t(i)
        dgqe(i) = ay(delta_index + 1) * grhov_t(i) * (1 + w(i))
      end if
    end do
  end subroutine TMultiFluidDE_PerturbedStressEnergy

  function w_derivative(this, a, i) result(dwda)
    ! Returns dw/da, necessary for calculating the adiabatic sound speed
    ! Check equations 22, 23, 24 in https://arxiv.org/pdf/1806.10608.pdf
    class(TMultiFluidDE), intent(in) :: this
    integer, intent(in) :: i
    real(dl), intent(in) :: a
    real(dl) :: dwda
    real(dl) :: ac, zc, fde_zc
    real(dl) :: wn

    select case (i)
      case(1)
        if (this%models(i) == 1) then
          ! w_constant model
          dwda = 0
        else
          stop "[Multifluid DE] Invalid dark energy model for fluid 1"
        end if

      case(2)
        if (this%models(i) == 1) then
          ! Axion effective fluid from arXiv:1806.10608
          zc = this%zc
          wn = this%wn
          ac = 1._dl/(1._dl + zc)
          fde_zc = this%fde_zc
          dwda = 3 * ac * (1+wn)**2 * (ac/a)**(-1 + 3*(1+wn)) &
          / (a**2 * (1 + (ac/a)**(3*(1+wn)) )**2) 
        else
          stop "[Multifluid DE] Invalid dark energy model for fluid 2"
        end if
      case default
        stop "[Multifluid DE] Invalid dark energy fluid"
    end select
  end function w_derivative

  subroutine TMultiFluidDE_PerturbationEvolve(this, ayprime, w, w_ix, &
    a, adotoa, k, z, y)
    class(TMultiFluidDE), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: y(:)
    real(dl), dimension(max_num_of_fluids), intent(in) :: w
    real(dl), intent(in) :: a, adotoa, k, z
    integer, intent(in) :: w_ix
    real(dl) :: Hv3_over_k, v, delta
    real(dl), parameter :: cs2 = 1._dl
    integer :: i, delta_index

    ! Fluid equations in synchronous gauge https://arxiv.org/pdf/1806.10608.pdf
    ! Assuming cs_2 = 1 for all fluids
    ! PPF has only one equation for Gamma and it is initialized in the PerturbedStressEnergy subroutine
    do i=1, this%num_of_components
      delta_index = w_ix + 2*(i-1)
      if (this%models(1) == 2) then
        delta_index = delta_index - 1
        if (i == 1) then
          cycle
        end if
      end if
      delta = y(delta_index)
      v = y(delta_index + 1)
      Hv3_over_k =  3 * adotoa * v / k
      ! Density perturbation evolution
      ayprime(delta_index) = - 3 * adotoa * (cs2 - w(i)) * delta &
                                - 3 * adotoa * (1 + w(i)) * (cs2 - w(i)) * Hv3_over_k &
                                - (1 + w(i)) * k * v - (1 + w(i)) * k * z
      ! The remaining term is the one involving the adiabatic sound speed
      ayprime(delta_index) = ayprime(delta_index) - adotoa*Hv3_over_k*a*this%w_derivative(a,i)
      ! JVR - checked the term involving 1/3 * (d(ln(1+w))/lna)
      ! Velocity
      if (abs(w(i)+1) > 1e-6) then
          ayprime(delta_index + 1) = - adotoa * (1 - 3 * cs2) * v + &
              k * cs2 * delta / (1 + w(i))
      else
          ayprime(delta_index + 1) = 0
      end if
    end do
  end subroutine TMultiFluidDE_PerturbationEvolve

  subroutine TMultiFluidDE_Effective_w_wa(this, w, wa)
    class(TMultiFluidDE), intent(inout) :: this
    real(dl), intent(out) :: w, wa
    ! JVR - For w_const, effective w_wa = w
    ! Don't need effective wwa for EDE
    if (this%models(1) == 1) then
      !w = this%de_params(1)
      w = this%w0
      wa = 0
    else if (this%models(1) == 2) then
      w = this%w0
      wa = this%wa
    else
      stop "[Multifluid DE] Effective w0 wa not implemented for this DE model"
    end if
  end subroutine TMultiFluidDE_Effective_w_wa

  ! --------------------- PPF specific subroutines ---------------------

  subroutine PPF_Perturbations(this, dgrhoe, dgqe, &
    a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
    class(TMultiFluidDE), intent(inout) :: this
    real(dl), dimension(max_num_of_fluids), intent(out) :: dgrhoe, dgqe
    real(dl), dimension(max_num_of_fluids), intent(in) :: grhov_t
    real(dl), dimension(max_num_of_fluids), intent(in) :: w
    real(dl), intent(in) ::  a, dgq, dgrho, grho, gpres_noDE, etak, adotoa, k, kf1
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: w_ix
    integer :: i, delta_index
    real(dl) :: grhoT, vT, k2, sigma, S_Gamma, ckH, Gamma, Gammadot, Fa, c_Gamma_ppf

    k2=k**2
    grhoT = grho - grhov_t(1)
    vT = dgq / (grhoT + gpres_noDE)
    Gamma = ay(w_ix)
    c_Gamma_ppf = 0.4_dl

    !sigma for ppf
    sigma = (etak + (dgrho + 3 * adotoa / k * dgq) / 2._dl / k) / kf1 - &
        k * Gamma
    sigma = sigma / adotoa

    S_Gamma = grhov_t(1) * (1 + w(1)) * (vT + sigma) * k / adotoa / 2._dl / k2
    ckH = c_Gamma_ppf * k / adotoa

    if (ckH * ckH .gt. 3.d1) then ! ckH^2 > 30 ?????????
        Gamma = 0
        Gammadot = 0.d0
    else
        Gammadot = S_Gamma / (1 + ckH * ckH) - Gamma - ckH * ckH * Gamma
        Gammadot = Gammadot * adotoa
    endif
    ayprime(w_ix) = Gammadot !Set this here, and don't use PerturbationEvolve

    Fa = 1 + 3 * (grhoT + gpres_noDE) / 2._dl / k2 / kf1
    dgqe(1) = S_Gamma - Gammadot / adotoa - Gamma
    dgqe(1) = -dgqe(1) / Fa * 2._dl * k * adotoa + vT * grhov_t(1) * (1 + w(1))
    dgrhoe(1) = -2 * k2 * kf1 * Gamma - 3 / k * adotoa * dgqe(1)
  end subroutine PPF_Perturbations

  function TMultiFluidDE_diff_rhopi_Add_Term(this, dgrhoe, dgqe, grho, gpres, w,  grhok, adotoa, &
        Kf1, k, grhov_t, z, k2, yprime, y, w_ix) result(ppiedot)
    ! [PPF] Get derivative of anisotropic stress
    class(TMultifluidDE), intent(in) :: this
    real(dl), intent(in) :: dgrhoe(max_num_of_fluids), dgqe(max_num_of_fluids), grho, gpres, w(max_num_of_fluids), grhok, adotoa, &
        k, grhov_t(max_num_of_fluids), z, k2, yprime(:), y(:), Kf1
    integer, intent(in) :: w_ix
    real(dl) :: ppiedot, hdotoh

    if (this%models(1) /= 2) then
      ppiedot = 0
    else
      hdotoh = (-3._dl * grho - 3._dl * gpres - 2._dl * grhok) / 6._dl / adotoa
      ppiedot = 3._dl * dgrhoe(1) + dgqe(1) * &
            (12._dl / k * adotoa + k / adotoa - 3._dl / k * (adotoa + hdotoh)) + &
            grhov_t(1) * (1 + w(1)) * k * z / adotoa - 2._dl * k2 * Kf1 * &
            (yprime(w_ix) / adotoa - 2._dl * y(w_ix))
      ppiedot = ppiedot * adotoa / Kf1
    end if
    
  end function TMultiFluidDE_diff_rhopi_Add_Term

 ! --------------------- Technical Subroutines ---------------------

  subroutine TMultiFluidDE_Init(this, State)
    use classes
    use results
    class(TMultiFluidDE), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State
    real(dl) :: ac

    ! JVR - If we need any information from other components, Init is the place to get it

    this%is_cosmological_constant = .false.
    this%num_perturb_equations = 2 * this%num_of_components
    if (this%models(1) == 2) then
      this%num_perturb_equations = this%num_perturb_equations - 1
    end if

    ac = 1._dl/(1+this%zc)

    select type (State)
    type is (CAMBdata)
      Omega_de = State%Omega_de
      grhode_today = State%grhov
      grhocrit = State%grhocrit
      this%grhonode_zc = State%grho_no_de(ac) / ac**4
    end select
  end subroutine TMultiFluidDE_Init

  subroutine TMultiFluidDE_ReadParams(this, Ini)
    use IniObjects
    class(TMultiFluidDE) :: this
    class(TIniFile), intent(in) :: Ini
    integer :: i, j
    character(len=30) :: read_model

    do i = 1, this%num_of_components
      write(read_model, "(A12,I1)") "model_fluid_", i
      this%models(i) = Ini%Read_Int(trim(read_model), 1)
    end do

    do i = 1, this%num_of_components
      do j = 1, max_num_of_params
        write(read_model, "(A12,I1)") "de_fluid_", i, "_param_", j
        this%models(i) = Ini%Read_Int(trim(read_model), 0)
      end do
    end do
  end subroutine TMultiFluidDE_ReadParams

  function TMultiFluidDE_PythonClass()
    character(LEN=:), allocatable :: TMultiFluidDE_PythonClass
    TMultiFluidDE_PythonClass = 'MultiFluidDE'
  end function TMultiFluidDE_PythonClass

  subroutine TMultiFluidDE_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TMultiFluidDE), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType
  end subroutine TMultiFluidDE_SelfPointer

end module MultiFluidDE