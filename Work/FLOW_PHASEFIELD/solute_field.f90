!
! the solute module, aiming to deal with the solute
!
module solute_field_module

  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use multifab_fill_ghost_module
  !use ml_restriction_module
  use pf_utility_module
  use ml_restrict_fill_module

  implicit none

  private

  public :: advance_solute_field

contains

subroutine advance_solute_field(mla,flw,phi_old,phi_new,dx,dt,the_bc_tower,pfpara)
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: phi_old(:),phi_new(:)
    type(multifab) , intent(in   ) :: flw(:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: dt
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(pf_para),   intent(in   ) :: pfpara
 
    select case(pfpara%do_solve_flow)
    case(0)
      call advance_solute_field_den(mla,phi_old,phi_new,dx,dt,the_bc_tower,pfpara)
    case(1)
      call advance_solute_field_flw(mla,flw,phi_old,phi_new,dx,dt,the_bc_tower,pfpara)
    end select

end subroutine advance_solute_field

subroutine advance_solute_field_den(mla,phi_old,phi_new,dx,dt,the_bc_tower,pfpara)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: phi_old(:),phi_new(:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: dt
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(pf_para),   intent(in   ) :: pfpara

    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: dm, ng_p, i,n,nlevs
    real(kind=dp_t) :: var_h, var_l,var_dt
    real(kind=dp_t) :: rkm_con1, rkm_con2, rkm_con3

    real(kind=dp_t), pointer ::  pfo(:,:,:,:)
    real(kind=dp_t), pointer ::  pfb(:,:,:,:)
    real(kind=dp_t), pointer ::  pfn(:,:,:,:)
    real(kind=dp_t), pointer ::  rkm(:,:,:,:)

    rkm_con1 = 2.d0/3.d0
    rkm_con2 = 0.25d0
    rkm_con3 = 0.75d0

    dm   = mla%dim
    nlevs = mla%nlevel
    ng_p = phi_old(1)%ng

    var_h = 0.d0
    var_l = -1.d0
    var_dt = rkm_con1 * dt

    ! calculate the RK matrices
    do n=1,nlevs
    do i=1,nfabs(phi_old(n))
       pfo  => dataptr(phi_old(n),i)
       pfn  => dataptr(phi_new(n),i)
       lo = lwb(get_box(phi_old(n),i))
       hi = upb(get_box(phi_old(n),i))

       select case(dm)
       case (2)
          call cal_rkmat_uc_2d(pfn(:,:,1,2),pfo(:,:,1,1),pfn(:,:,1,1),pfo(:,:,1,2),ng_p, &
                               lo,hi,dx(n),dt,pfpara%thers,pfpara%DM,pfpara%sk,&
                               -1.d0,0.d0)

       case (3)
          call cal_rkmat_uc_3d(pfn(:,:,:,2),pfo(:,:,:,1),pfn(:,:,:,1),pfo(:,:,:,2),ng_p, &
                               lo,hi,dx(n),dt,pfpara%thers,pfpara%DM,pfpara%sk,&
                               -1.d0,0.d0)

       end select
    end do
    end do

    ! update the new phase field and orientaion field
    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary_c(phi_new(nlevs),2,1)

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(phi_new(nlevs),2,1,1,the_bc_tower%bc_tower_array(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1
          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction_c(phi_new(n-1),2,phi_new(n),2,mla%mba%rr(n-1,:),1)

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(phi_new(n),phi_new(n-1),phi_new(n)%ng,mla%mba%rr(n-1,:), &
                                         the_bc_tower%bc_tower_array(n-1), &
                                         the_bc_tower%bc_tower_array(n), &
                                         2,1,1)
       end do

    end if

end subroutine advance_solute_field_den

! subroutine calculate_rkmatrices_uc_2d(...) aims to calculate the RK matrices for phase field
! input variable:
!            phi    :    solute field
!            ad     :    predefined orientation, should be calculated
!            uc     :    the solute field, matrix
!            th     :    the temperature field, matrix
!            noise  :    the noise matrix
!            fluxx  :    flux along x direction
!            fluxy  :    flux along y direction
!            ng_p   :    number of ghost cells used for phase field
!            ng_f   :    number of ghost cells used for flux
!            lo     :    lower corner of the current box
!            hi     :    higher corner of the current box
!            dx     :    spatial step
!            thers  :    a small value indicating the closeness to zero
!            kk     :    hamonic symmetric factor, could be 4 for cubic or 6 for hcp
!            Le     :    the Lewis number, i.e. ratio between thermal and solutal diffusivities
!            Mc00   :    the solute scaled factor, |m|*c00*(1-k)/(L/Cp)
!            p_k    :    the partition coefficient
!            lamda  :    the scaling parameter
! output variables:
!              K_R    :    RK matrices

subroutine cal_rkmat_uc_2d(K_R,pf_old,pf_new,uc_old,ng_p,lo,hi,dx,dt,thers,DM,p_k,var_l,var_h)

  integer           :: lo(2),hi(2),ng_p
  double precision  :: K_R(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: pf_old(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: pf_new(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: uc_old(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: dx,dt,DM,p_k,thers,var_l,var_h

  ! local variables
  integer           :: i,j,x_c,y_c
  double precision  :: pA,pE,pft,JR,JL,JT,JB

  pA  = 0.d0
  pE  = 0.d0
  pft = 0.d0
  JR  = 0.d0
  JL  = 0.d0
  JT  = 0.d0
  JB  = 0.d0

  x_c = 1
  y_c = 2

  ! do the loops
  !$omp parallel do private(j,i,pft,JR,JL,JT,JB,pA,pE)
  do j=lo(2),hi(2)
    do i=lo(1),hi(1)

      pft = (pf_new(i,j) - pf_old(i,j)) / dt

      call cal_sur_sol_flux_2d(x_c,JL,i,j,pf_old,pf_new,uc_old,ng_p,lo,hi,dx,dt,thers,DM,p_k)
      call cal_sur_sol_flux_2d(x_c,JR,i+1,j,pf_old,pf_new,uc_old,ng_p,lo,hi,dx,dt,thers,DM,p_k)
      call cal_sur_sol_flux_2d(y_c,JB,i,j,pf_old,pf_new,uc_old,ng_p,lo,hi,dx,dt,thers,DM,p_k)
      call cal_sur_sol_flux_2d(y_c,JT,i,j+1,pf_old,pf_new,uc_old,ng_p,lo,hi,dx,dt,thers,DM,p_k)

      pA = 0.5d0 * ( 1.d0 + (1.d0 - p_k) * uc_old(i,j) ) * pft
      pE = (1.d0 + p_k) / 2.d0 - ((1.d0 - p_k)/2.d0) * pf_old(i,j)

      K_R(i,j) = ( JR - JL + JT - JB ) / dx + pA
      K_R(i,j) = K_R(i,j) / pE

      K_R(i,j) = uc_old(i,j) + dt * K_R(i,j)

      ! cut off the bouds
      !if (K_R(i,j) > var_h) then
      !   K_R(i,j) = var_h
      !else if(K_R(i,j) < var_l) then
      !   K_R(i,j) = var_l
      !end if

    end do
  end do
  !$omp end parallel do

end subroutine cal_rkmat_uc_2d

! subroutine calculate_rkmatrices_uc_3d(...) aims to calculate the RK matrices for phase field
! input variable:
!            phi    :    solute field
!            ad     :    predefined orientation, should be calculated
!            uc     :    the solute field, matrix
!            th     :    the temperature field, matrix
!            noise  :    the noise matrix
!            fluxx  :    flux along x direction
!            fluxy  :    flux along y direction
!            ng_p   :    number of ghost cells used for phase field
!            ng_f   :    number of ghost cells used for flux
!            lo     :    lower corner of the current box
!            hi     :    higher corner of the current box
!            dx     :    spatial step
!            thers  :    a small value indicating the closeness to zero
!            kk     :    hamonic symmetric factor, could be 4 for cubic or 6 for hcp
!            Le     :    the Lewis number, i.e. ratio between thermal and solutal diffusivities
!            Mc00   :    the solute scaled factor, |m|*c00*(1-k)/(L/Cp)
!            p_k    :    the partition coefficient
!            lamda  :    the scaling parameter
! output variables:
!              K_R    :    RK matrices

subroutine cal_rkmat_uc_3d(K_R,pf_old,pf_new,uc_old,ng_p,lo,hi,dx,dt,thers,DM,p_k,var_l,var_h)

  integer           :: lo(3),hi(3),ng_p
  double precision  :: K_R(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: pf_old(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: pf_new(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: uc_old(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: dx,dt,DM,p_k,thers,var_l,var_h

  ! local variables
  integer           :: i,j,k,x_c,y_c,z_c
  double precision  :: pA,pE,pft,JR,JL,JT,JB,JF,JK

  pA  = 0.d0
  pE  = 0.d0
  pft = 0.d0
  JR  = 0.d0
  JL  = 0.d0
  JT  = 0.d0
  JB  = 0.d0
  JF  = 0.d0
  JK  = 0.d0

  x_c = 1
  y_c = 2
  z_c = 3

  ! do the loops
  !$omp parallel do private(k,j,i,pft,JL,JR,JB,JT,JK,JF,pA,pE)
  do k=lo(3),hi(3)
    do j=lo(2),hi(2)
      do i=lo(1),hi(1)

        pft = (pf_new(i,j,k) - pf_old(i,j,k)) / dt

        call cal_sur_sol_flux_3d(x_c,JL,i,  j,  k,  pf_old,pf_new,uc_old,ng_p,lo,hi,dx,dt,thers,DM,p_k)
        call cal_sur_sol_flux_3d(x_c,JR,i+1,j,  k,  pf_old,pf_new,uc_old,ng_p,lo,hi,dx,dt,thers,DM,p_k)
        call cal_sur_sol_flux_3d(y_c,JB,i,  j,  k,  pf_old,pf_new,uc_old,ng_p,lo,hi,dx,dt,thers,DM,p_k)
        call cal_sur_sol_flux_3d(y_c,JT,i,  j+1,k,  pf_old,pf_new,uc_old,ng_p,lo,hi,dx,dt,thers,DM,p_k)
        call cal_sur_sol_flux_3d(z_c,JK,i,  j,  k,  pf_old,pf_new,uc_old,ng_p,lo,hi,dx,dt,thers,DM,p_k)
        call cal_sur_sol_flux_3d(z_c,JF,i,  j,  k+1,pf_old,pf_new,uc_old,ng_p,lo,hi,dx,dt,thers,DM,p_k)

        pA = 0.5d0 * ( 1.d0 + (1.d0 - p_k) * uc_old(i,j,k) ) * pft
        pE = (1.d0 + p_k) / 2.d0 - ((1.d0 - p_k)/2.d0) * pf_old(i,j,k)

        K_R(i,j,k) = ( JR - JL + JT - JB + JF - JK ) / dx + pA
        K_R(i,j,k) = K_R(i,j,k) / pE

        K_R(i,j,k) = uc_old(i,j,k) + dt * K_R(i,j,k)

        ! cut off the bouds
        !if (K_R(i,j,k) > var_h) then
        !   K_R(i,j,k) = var_h
        !else if(K_R(i,j,k) < var_l) then
        !   K_R(i,j,k) = var_l
        !end if

      end do
    end do
  end do
  !$omp end parallel do

end subroutine cal_rkmat_uc_3d

! calculate the flux along x direction (s_case = 1), remember to use
!                 (i,j)    to calculate flux at left and
!                 (i+1,j)  to calculate flux at right
! calculate the flux along y direction (s_case = 2), remember to use
!                 (i,j)    to calculate flux at bottom and
!                 (i,j+1)  to calculate flux at top

subroutine cal_sur_sol_flux_2d(s_case,flux,i,j,pf_old,pf_new,uc_old,ng_p,lo,hi,dx,dt,thers,DM,p_k)

  integer           :: lo(2),hi(2),ng_p,s_case,i,j
  double precision  :: pf_old(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: pf_new(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: uc_old(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: dx,dt,thers,DM,p_k,flux

  ! local variables
  double precision  :: dfan_dx,dfan_dy,dfan_dt,fan_c,fan_p,Urel,mod_c,Js,Ja

  flux    = 0.d0
  dfan_dx = 0.d0
  dfan_dy = 0.d0
  dfan_dt = 0.d0
  fan_c   = 0.d0
  fan_p   = 0.d0
  Urel    = 0.d0
  mod_c   = 0.d0
  Js      = 0.d0
  Ja      = 0.d0

  select case (s_case)
    case (1)
      dfan_dx = ( pf_old(i,j) - pf_old(i-1,j) ) / dx
      fan_p   = 0.5d0 * (pf_old(i,j) + pf_old(i-1,j))
      Js      = DM * 0.5 * (1.d0 - fan_p) * (uc_old(i,j) - uc_old(i-1,j)) / dx

      if (abs(dfan_dx) <= thers) then
          Ja = 0.d0
      else
          dfan_dy = ( pf_old(i,j+1) + pf_old(i-1,j+1) - pf_old(i,j-1) - pf_old(i-1,j-1) ) / (4.d0 * dx)
          fan_c   = 0.5d0 * (pf_new(i,j) + pf_new(i-1,j))
          Urel    = 0.5d0 * (uc_old(i,j) + uc_old(i-1,j))
          mod_c = sqrt(dfan_dx*dfan_dx+dfan_dy*dfan_dy)
          dfan_dt = (fan_c - fan_p) / dt
          Ja = sqrt(2.d0) / 4.d0 * ( 1.d0 + ( 1.d0 - p_k) *Urel) * dfan_dt * dfan_dx / mod_c
      end if

      flux = Js + Ja
      !flux = Js

  case (2)

    dfan_dy = ( pf_old(i,j) - pf_old(i,j-1) ) / dx
    fan_p   = 0.5d0 * (pf_old(i,j) + pf_old(i,j-1))
    Js      = DM * 0.5 * (1.d0 - fan_p) * (uc_old(i,j) - uc_old(i,j-1)) / dx

    if (abs(dfan_dy) <= thers) then
        Ja = 0.d0
    else
        dfan_dx = ( pf_old(i+1,j) + pf_old(i+1,j-1) - pf_old(i-1,j) - pf_old(i-1,j-1) ) / (4.d0 * dx)
        fan_c   = 0.5d0 * (pf_new(i,j) + pf_new(i,j-1))
        Urel    = 0.5d0 * (uc_old(i,j) + uc_old(i,j-1))
        mod_c = sqrt(dfan_dx*dfan_dx+dfan_dy*dfan_dy)
        dfan_dt = (fan_c - fan_p) / dt
        Ja = sqrt(2.d0) / 4.d0 * ( 1.d0 + ( 1.d0 - p_k) *Urel) * dfan_dt * dfan_dy / mod_c
    end if

    flux = Js + Ja
    !flux = Js

  end select

end subroutine cal_sur_sol_flux_2d

subroutine cal_sur_sol_flux_3d(s_case,flux,i,j,k,pf_old,pf_new,uc_old,ng_p,lo,hi,dx,dt,thers,DM,p_k)

  integer           :: lo(3),hi(3),ng_p,s_case,i,j,k
  double precision  :: pf_old(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: pf_new(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: uc_old(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: dx,dt,thers,DM,p_k,flux

  ! local variables
  double precision  :: dfan_dx,dfan_dy,dfan_dz,dfan_dt
  double precision  :: fan_c,fan_p,Urel,mod_c,Js,Ja

  flux    = 0.d0
  dfan_dx = 0.d0
  dfan_dy = 0.d0
  dfan_dz =0.d0
  dfan_dt = 0.d0
  fan_c   = 0.d0
  fan_p   = 0.d0
  Urel    = 0.d0
  mod_c   = 0.d0
  Js      = 0.d0
  Ja      = 0.d0

  select case (s_case)
    case (1)
      dfan_dx = ( pf_old(i,j,k) - pf_old(i-1,j,k) ) / dx
      fan_p   = 0.5d0 * (pf_old(i,j,k) + pf_old(i-1,j,k))
      Js      = DM * 0.5 * (1.d0 - fan_p) * (uc_old(i,j,k) - uc_old(i-1,j,k)) / dx

      if (abs(dfan_dx) <= thers) then
          Ja = 0.d0
      else
          ! calculate related terms
          dfan_dy = ( pf_old(i,j+1,k) + pf_old(i-1,j+1,k) - pf_old(i,j-1,k) - pf_old(i-1,j-1,k) ) / (4.d0 * dx)
          dfan_dz = ( pf_old(i,j,k+1) + pf_old(i-1,j,k+1) - pf_old(i,j,k-1) - pf_old(i-1,j,k-1) ) / (4.d0 * dx)
          fan_c   = 0.5d0 * (pf_new(i,j,k) + pf_new(i-1,j,k))
          Urel    = 0.5d0 * (uc_old(i,j,k) + uc_old(i-1,j,k))
          mod_c = sqrt(dfan_dx**2 + dfan_dy**2 + dfan_dz**2)
          dfan_dt = (fan_c - fan_p) / dt

          ! calculate the antitrapping current
          Ja = sqrt(2.d0) / 4.d0 * ( 1.d0 + ( 1.d0 - p_k) * Urel) * dfan_dt * dfan_dx / mod_c

      end if

      flux = Js + Ja
      !flux = Js

  case (2)
    dfan_dy = ( pf_old(i,j,k) - pf_old(i,j-1,k) ) / dx
    fan_p   = 0.5d0 * (pf_old(i,j,k) + pf_old(i,j-1,k))
    Js      = DM * 0.5 * (1.d0 - fan_p) * (uc_old(i,j,k) - uc_old(i,j-1,k)) / dx

    if (abs(dfan_dy) <= thers) then
        Ja = 0.d0
    else
        ! calculate related terms
        dfan_dx = ( pf_old(i+1,j,k) + pf_old(i+1,j-1,k) - pf_old(i-1,j,k) - pf_old(i-1,j-1,k) ) / (4.d0 * dx)
        dfan_dz = ( pf_old(i,j,k+1) + pf_old(i,j-1,k+1) - pf_old(i,j,k-1) - pf_old(i,j-1,k-1) ) / (4.d0 * dx)
        fan_c   = 0.5d0 * (pf_new(i,j,k) + pf_new(i,j-1,k))
        Urel    = 0.5d0 * (uc_old(i,j,k) + uc_old(i,j-1,k))
        mod_c = sqrt(dfan_dx**2 + dfan_dy**2 + dfan_dz**2)
        dfan_dt = (fan_c - fan_p) / dt

        ! calculate the antitrapping current
        Ja = sqrt(2.d0) / 4.d0 * ( 1.d0 + ( 1.d0 - p_k) *Urel) * dfan_dt * dfan_dy / mod_c

    end if

    flux = Js + Ja
    !flux = Js

  case (3)
    dfan_dz = ( pf_old(i,j,k) - pf_old(i,j,k-1) ) / dx
    fan_p   = 0.5d0 * (pf_old(i,j,k) + pf_old(i,j,k-1))
    Js      = DM * 0.5 * (1.d0 - fan_p) * (uc_old(i,j,k) - uc_old(i,j,k-1)) / dx

    if (abs(dfan_dz) <= thers) then
        Ja = 0.d0
    else
        !calculate related terms
        dfan_dx = ( pf_old(i+1,j,k) + pf_old(i+1,j,k-1) - pf_old(i-1,j,k) - pf_old(i-1,j,k-1) ) / (4.d0 * dx)
        dfan_dy = ( pf_old(i,j+1,k) + pf_old(i,j+1,k-1) - pf_old(i,j-1,k) - pf_old(i,j-1,k-1) ) / (4.d0 * dx)
        fan_c   = 0.5d0 * (pf_new(i,j,k) + pf_new(i,j,k-1))
        Urel    = 0.5d0 * (uc_old(i,j,k) + uc_old(i,j,k-1))
        mod_c = sqrt(dfan_dx**2 + dfan_dy**2 + dfan_dz**2)
        dfan_dt = (fan_c - fan_p) / dt

        ! calculate the antitrapping current
        Ja = sqrt(2.d0) / 4.d0 * ( 1.d0 + ( 1.d0 - p_k) *Urel) * dfan_dt * dfan_dz / mod_c

    end if

    flux = Js + Ja
    !flux = Js

  end select

end subroutine cal_sur_sol_flux_3d

subroutine advance_solute_field_flw(mla,flw,phi_old,phi_new,dx,dt,the_bc_tower,pfpara)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: phi_old(:),phi_new(:)
    type(multifab) , intent(in   ) :: flw(:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: dt
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(pf_para),   intent(in   ) :: pfpara

    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: dm, ng_p, i,n,nlevs
    real(kind=dp_t) :: var_h, var_l,var_dt
    real(kind=dp_t) :: rkm_con1, rkm_con2, rkm_con3

    real(kind=dp_t), pointer ::  pfo(:,:,:,:)
    real(kind=dp_t), pointer ::  pfb(:,:,:,:)
    real(kind=dp_t), pointer ::  pfn(:,:,:,:)
    real(kind=dp_t), pointer ::  rkm(:,:,:,:)
    real(kind=dp_t), pointer ::  pfl(:,:,:,:)

    if(pfpara%solve_flow_only == 1) return

    rkm_con1 = 2.d0/3.d0
    rkm_con2 = 0.25d0
    rkm_con3 = 0.75d0

    dm   = mla%dim
    nlevs = mla%nlevel
    ng_p = phi_old(1)%ng

    var_h = 0.d0
    var_l = -1.d0
    var_dt = rkm_con1 * dt

    ! calculate the RK matrices
    do n=1,nlevs
    do i=1,nfabs(phi_old(n))
       pfo  => dataptr(phi_old(n),i)
       pfn  => dataptr(phi_new(n),i)
       pfl  => dataptr(flw(n),i)
       lo = lwb(get_box(phi_old(n),i))
       hi = upb(get_box(phi_old(n),i))

       select case(dm)
       case (2)
          call cal_rkmat_uc_2d_flw(pfn(:,:,1,2),pfl(:,:,1,1),pfl(:,:,1,2),pfo(:,:,1,1),pfn(:,:,1,1),&
                                   pfo(:,:,1,2),ng_p,lo,hi,dx(n),dt,pfpara%thers,pfpara%DM,pfpara%sk,&
                                   pfpara%flw_amp,-1.d0,0.d0)

       case (3)
          call cal_rkmat_uc_3d_flw(pfn(:,:,:,2),pfl(:,:,:,1),pfl(:,:,:,2),pfl(:,:,:,3),pfo(:,:,:,1),&
                                   pfn(:,:,:,1),pfo(:,:,:,2),ng_p,lo,hi,dx(n),dt,pfpara%thers,&
                                   pfpara%DM,pfpara%sk,pfpara%flw_amp,-1.d0,0.d0)

       end select
    end do
    end do

    ! update the new phase field and orientaion field
    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary_c(phi_new(nlevs),2,1)

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(phi_new(nlevs),2,1,1,the_bc_tower%bc_tower_array(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1
          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction_c(phi_new(n-1),2,phi_new(n),2,mla%mba%rr(n-1,:),1)

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(phi_new(n),phi_new(n-1),phi_new(n)%ng,mla%mba%rr(n-1,:), &
                                         the_bc_tower%bc_tower_array(n-1), &
                                         the_bc_tower%bc_tower_array(n), &
                                         2,1,1)
       end do

    end if

end subroutine advance_solute_field_flw

! subroutine calculate_rkmatrices_uc_2d(...) aims to calculate the RK matrices for phase field
! input variable:
!            phi    :    solute field
!            ad     :    predefined orientation, should be calculated
!            uc     :    the solute field, matrix
!            th     :    the temperature field, matrix
!            noise  :    the noise matrix
!            fluxx  :    flux along x direction
!            fluxy  :    flux along y direction
!            ng_p   :    number of ghost cells used for phase field
!            ng_f   :    number of ghost cells used for flux
!            lo     :    lower corner of the current box
!            hi     :    higher corner of the current box
!            dx     :    spatial step
!            thers  :    a small value indicating the closeness to zero
!            kk     :    hamonic symmetric factor, could be 4 for cubic or 6 for hcp
!            Le     :    the Lewis number, i.e. ratio between thermal and solutal diffusivities
!            Mc00   :    the solute scaled factor, |m|*c00*(1-k)/(L/Cp)
!            p_k    :    the partition coefficient
!            lamda  :    the scaling parameter
! output variables:
!              K_R    :    RK matrices

subroutine cal_rkmat_uc_2d_flw(K_R,u_x,u_y,pf_old,pf_new,uc_old,ng_p,lo,hi,dx,dt,thers,DM,&
                               p_k,flw_amp,var_l,var_h)

  integer           :: lo(2),hi(2),ng_p
  double precision  :: K_R(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: u_x(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: u_y(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: pf_old(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: pf_new(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: uc_old(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: dx,dt,DM,p_k,thers,var_l,var_h,flw_amp

  ! local variables
  integer           :: i,j,x_c,y_c
  double precision  :: pA,pE,pft,JR,JL,JT,JB
  double precision  :: p_flw, JL_F, JR_F,JT_F, JB_F
  double precision  :: fan_term, u_term

  pA  = 0.d0; pE  = 0.d0; pft = 0.d0
  JR  = 0.d0; JL  = 0.d0; JT  = 0.d0; JB  = 0.d0
  p_flw = 0.d0; JL_F  = 0.d0; JR_F  = 0.d0; JT_F  = 0.d0; JB_F  = 0.d0
  fan_term= 0.d0; u_term = 0.d0;

  x_c = 1
  y_c = 2

  ! do the loops
  !$omp parallel do private(j,i,pft,JR,JL,JT,JB,pA,pE)
  do j=lo(2),hi(2)
    do i=lo(1),hi(1)

      pft = (pf_new(i,j) - pf_old(i,j)) / dt

      call cal_sur_sol_flux_2d_flw(x_c,JL_F,JL,i,j,  pf_old,pf_new,uc_old,ng_p,lo,hi,dx,dt,thers,DM,p_k)
      call cal_sur_sol_flux_2d_flw(x_c,JR_F,JR,i+1,j,pf_old,pf_new,uc_old,ng_p,lo,hi,dx,dt,thers,DM,p_k)
      call cal_sur_sol_flux_2d_flw(y_c,JB_F,JB,i,j,  pf_old,pf_new,uc_old,ng_p,lo,hi,dx,dt,thers,DM,p_k)
      call cal_sur_sol_flux_2d_flw(y_c,JT_F,JT,i,j+1,pf_old,pf_new,uc_old,ng_p,lo,hi,dx,dt,thers,DM,p_k)

      pA = 0.5d0 * ( 1.d0 + (1.d0 - p_k) * uc_old(i,j) ) * pft
      pE = (1.d0 + p_k) / 2.d0 - ((1.d0 - p_k)/2.d0) * pf_old(i,j)

      ! let's calculate the convection term explicitly
      fan_term = 1.d0 + p_k - (1.d0 - p_k)*pf_new(i,j)
      u_term   = 1.d0 + (1.d0 - p_k)*uc_old(i,j)
      p_flw = u_x(i,j)*( fan_term *(uc_old(i+1,j) -uc_old(i-1,j)) -u_term*( pf_new(i+1,j) - pf_new(i-1,j)) ) + &
              u_y(i,j)*( fan_term *(uc_old(i,j+1) -uc_old(i,j-1)) -u_term*( pf_new(i,j+1) - pf_new(i,j-1)) )
      p_flw = p_flw*flw_amp*(1.d0 - pf_new(i,j)) / (8.d0*dx)
      !p_flw = flw_amp*0.5d0*(1.d0 - pf_new(i,j))*( u_x(i,j)*(JR_F-JL_F) + u_y(i,j)*(JT_F-JB_F) )/(2.d0*dx*(1.d0 - p_k))

      K_R(i,j) = ( JR - JL + JT - JB ) / dx + pA - p_flw
      K_R(i,j) = K_R(i,j) / pE

      K_R(i,j) = uc_old(i,j) + dt * K_R(i,j)

      ! cut off the bouds
      !if (K_R(i,j) > var_h) then
      !   K_R(i,j) = var_h
      !else if(K_R(i,j) < var_l) then
      !   K_R(i,j) = var_l
      !end if

    end do
  end do
  !$omp end parallel do

end subroutine cal_rkmat_uc_2d_flw

subroutine cal_rkmat_uc_3d_flw(K_R,u_x,u_y,u_z,pf_old,pf_new,uc_old,ng_p,lo,hi,dx,dt,&
                               thers,DM,p_k,flw_amp,var_l,var_h)

  integer           :: lo(3),hi(3),ng_p
  double precision  :: K_R(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: u_x(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: u_y(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: u_z(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: pf_old(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: pf_new(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: uc_old(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: dx,dt,DM,p_k,thers,var_l,var_h,flw_amp

  ! local variables
  integer           :: i,j,k,x_c,y_c,z_c
  double precision  :: pA,pE,pft,JR,JL,JT,JB,JF,JK
  double precision  :: p_flw, JL_F, JR_F,JT_F, JB_F,JF_F,JK_F
  double precision  :: fan_term, u_term

  pA    = 0.d0; pE    = 0.d0; pft   = 0.d0
  JR    = 0.d0; JL    = 0.d0; JT    = 0.d0
  JB    = 0.d0; JF    = 0.d0; JK    = 0.d0
  p_flw = 0.d0; JL_F  = 0.d0; JR_F  = 0.d0; JT_F  = 0.d0
  JB_F  = 0.d0; JF_F  = 0.d0; JK_F  = 0.d0
  x_c   = 1; y_c   = 2; z_c   = 3; fan_term= 0.d0; u_term = 0.d0;

  ! do the loops
  !$omp parallel do private(k,j,i,pft,JL,JR,JB,JT,JK,JF,pA,pE)
  do k=lo(3),hi(3)
    do j=lo(2),hi(2)
      do i=lo(1),hi(1)

        pft = (pf_new(i,j,k) - pf_old(i,j,k)) / dt

        call cal_sur_sol_flux_3d_flw(x_c,JL_F,JL,i,  j,  k,  pf_old,pf_new,uc_old,ng_p,lo,hi,dx,dt,thers,DM,p_k)
        call cal_sur_sol_flux_3d_flw(x_c,JR_F,JR,i+1,j,  k,  pf_old,pf_new,uc_old,ng_p,lo,hi,dx,dt,thers,DM,p_k)
        call cal_sur_sol_flux_3d_flw(y_c,JB_F,JB,i,  j,  k,  pf_old,pf_new,uc_old,ng_p,lo,hi,dx,dt,thers,DM,p_k)
        call cal_sur_sol_flux_3d_flw(y_c,JT_F,JT,i,  j+1,k,  pf_old,pf_new,uc_old,ng_p,lo,hi,dx,dt,thers,DM,p_k)
        call cal_sur_sol_flux_3d_flw(z_c,JK_F,JK,i,  j,  k,  pf_old,pf_new,uc_old,ng_p,lo,hi,dx,dt,thers,DM,p_k)
        call cal_sur_sol_flux_3d_flw(z_c,JF_F,JF,i,  j,  k+1,pf_old,pf_new,uc_old,ng_p,lo,hi,dx,dt,thers,DM,p_k)

        pA = 0.5d0 * ( 1.d0 + (1.d0 - p_k) * uc_old(i,j,k) ) * pft
        pE = (1.d0 + p_k) / 2.d0 - ((1.d0 - p_k)/2.d0) * pf_old(i,j,k)

        fan_term = 1.d0 + p_k - (1.d0 - p_k)*pf_new(i,j,k)
        u_term   = 1.d0 + (1.d0 - p_k)*uc_old(i,j,k)
        p_flw = u_x(i,j,k)*( fan_term *(uc_old(i+1,j,k) -uc_old(i-1,j,k)) -u_term*( pf_new(i+1,j,k) - pf_new(i-1,j,k)) ) + &
                u_y(i,j,k)*( fan_term *(uc_old(i,j+1,k) -uc_old(i,j-1,k)) -u_term*( pf_new(i,j+1,k) - pf_new(i,j-1,k)) ) + &
                u_z(i,j,k)*( fan_term *(uc_old(i,j,k+1) -uc_old(i,j,k-1)) -u_term*( pf_new(i,j,k+1) - pf_new(i,j,k-1)) )
        p_flw = p_flw*flw_amp*(1.d0 - pf_new(i,j,k)) / (8.d0*dx)

        !p_flw = flw_amp*0.5d0*(1.d0 - pf_old(i,j,k))*(u_x(i,j,k)*(JR_F-JL_F)+u_y(i,j,k)*(JT_F-JB_F)+u_z(i,j,k)*&
        !        (JF_F-JK_F))/(2.d0*dx*(1.d0-p_k))

        K_R(i,j,k) = ( JR - JL + JT - JB + JF - JK ) / dx + pA - p_flw
        K_R(i,j,k) = K_R(i,j,k) / pE

        K_R(i,j,k) = uc_old(i,j,k) + dt * K_R(i,j,k)

        ! cut off the bouds
        !if (K_R(i,j,k) > var_h) then
        !   K_R(i,j,k) = var_h
        !else if(K_R(i,j,k) < var_l) then
        !   K_R(i,j,k) = var_l
        !end if

      end do
    end do
  end do
  !$omp end parallel do

end subroutine cal_rkmat_uc_3d_flw



! calculate the flux along x direction (s_case = 1), remember to use
!                 (i,j)    to calculate flux at left and
!                 (i+1,j)  to calculate flux at right
! calculate the flux along y direction (s_case = 2), remember to use
!                 (i,j)    to calculate flux at bottom and
!                 (i,j+1)  to calculate flux at top

subroutine cal_sur_sol_flux_2d_flw(s_case,flux_flw,flux,i,j,pf_old,pf_new,uc_old,ng_p,lo,hi,dx,dt,thers,DM,p_k)

  integer           :: lo(2),hi(2),ng_p,s_case,i,j
  !double precision  :: ux(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  !double precision  :: uy(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: pf_old(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: pf_new(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: uc_old(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: dx,dt,thers,DM,p_k,flux,flux_flw

  ! local variables
  double precision  :: dfan_dx,dfan_dy,dfan_dt,fan_c,fan_p
  double precision  :: Urel,mod_c,Js,Ja,fan_p2

  dfan_dx  = 0.d0; dfan_dy  = 0.d0; dfan_dt  = 0.d0
  fan_c    = 0.d0; fan_p    = 0.d0; Urel     = 0.d0 ;mod_c    = 0.d0
  Js       = 0.d0; Ja       = 0.d0; fan_p2   = 0.d0

  !---------------------------------------The flow part------------------------------------------------
  !double precision  :: fan_flw = 0.d0
  !----------------------------------------------------------------------------------------------------

  flux = 0.d0
  flux_flw = 0.d0

  select case (s_case)
    case (1)
      dfan_dx = ( pf_old(i,j) - pf_old(i-1,j) ) / dx
      fan_p   = 0.5d0 * (pf_old(i,j) + pf_old(i-1,j))
      Urel    = 0.5d0 * (uc_old(i,j) + uc_old(i-1,j))
      Js      = DM * 0.5 * (1.d0 - fan_p) * (uc_old(i,j) - uc_old(i-1,j)) / dx

      !---------------------------------------The flow part------------------------------------------------
      !fan_flw = 0.5d0 * (ux(i,j) + ux(i-1,j)) * fan_p * flw_amp
      !Js      = Js - 0.5d0 * fan_flw * ( 1.d0+p_k-(1.d0-p_k)*fan_p ) * (uc_old(i,j) - uc_old(i-1,j))
      
      !fan_p2   = 0.5d0 * (pf_new(i,j) + pf_new(i-1,j))
      !flux_flw = (1.d0+p_k-(1.d0-p_k)*fan_p2)*(1.d0+(1.d0-p_k)*Urel)
      !----------------------------------------------------------------------------------------------------

      if (abs(dfan_dx) <= thers) then
          Ja = 0.d0
      else
          dfan_dy = ( pf_old(i,j+1) + pf_old(i-1,j+1) - pf_old(i,j-1) - pf_old(i-1,j-1) ) / (4.d0 * dx)
          fan_c   = 0.5d0 * (pf_new(i,j) + pf_new(i-1,j))
          !Urel    = 0.5d0 * (uc_old(i,j) + uc_old(i-1,j))
          mod_c = sqrt(dfan_dx*dfan_dx+dfan_dy*dfan_dy)
          dfan_dt = (fan_c - fan_p) / dt
          Ja = sqrt(2.d0) / 4.d0 * ( 1.d0 + ( 1.d0 - p_k) *Urel) * dfan_dt * dfan_dx / mod_c

          !---------------------------------------The flow part------------------------------------------------
          !Ja = Ja + 0.5d0 * fan_flw * (1.d0+(1.d0-p_k)*Urel)* ( pf_old(i,j) - pf_old(i-1,j) )
          !----------------------------------------------------------------------------------------------------
      end if

      flux = Js + Ja
      !flux = Js

  case (2)

    dfan_dy = ( pf_old(i,j) - pf_old(i,j-1) ) / dx
    fan_p   = 0.5d0 * (pf_old(i,j) + pf_old(i,j-1))
    Urel    = 0.5d0 * (uc_old(i,j) + uc_old(i,j-1))
    Js      = DM * 0.5 * (1.d0 - fan_p) * (uc_old(i,j) - uc_old(i,j-1)) / dx

    !---------------------------------------The flow part------------------------------------------------
    !fan_flw = 0.5d0 * (uy(i,j) + uy(i,j-1)) * fan_p * flw_amp
    !Js      = Js - 0.5d0 * fan_flw * ( 1.d0+p_k-(1.d0-p_k)*fan_p ) * (uc_old(i,j) - uc_old(i,j-1))
    !fan_p2   = 0.5d0 * (pf_new(i,j) + pf_new(i,j-1))
    !flux_flw = (1.d0+p_k-(1.d0-p_k)*fan_p2)*(1.d0+(1.d0-p_k)*Urel)
    !----------------------------------------------------------------------------------------------------

    if (abs(dfan_dy) <= thers) then
        Ja = 0.d0
    else
        dfan_dx = ( pf_old(i+1,j) + pf_old(i+1,j-1) - pf_old(i-1,j) - pf_old(i-1,j-1) ) / (4.d0 * dx)
        fan_c   = 0.5d0 * (pf_new(i,j) + pf_new(i,j-1))
        !Urel    = 0.5d0 * (uc_old(i,j) + uc_old(i,j-1))
        mod_c = sqrt(dfan_dx*dfan_dx+dfan_dy*dfan_dy)
        dfan_dt = (fan_c - fan_p) / dt
        Ja = sqrt(2.d0) / 4.d0 * ( 1.d0 + ( 1.d0 - p_k) *Urel) * dfan_dt * dfan_dy / mod_c

        !---------------------------------------The flow part------------------------------------------------
        !Ja = Ja + 0.5d0 * fan_flw * (1.d0+(1.d0-p_k)*Urel)* ( pf_old(i,j) - pf_old(i,j-1) )
        !----------------------------------------------------------------------------------------------------
    end if

    flux = Js + Ja
    !flux = Js

  end select

end subroutine cal_sur_sol_flux_2d_flw


subroutine cal_sur_sol_flux_3d_flw(s_case,flux_flw,flux,i,j,k,pf_old,pf_new,uc_old,ng_p,lo,hi,dx,dt,thers,DM,p_k)

  integer           :: lo(3),hi(3),ng_p,s_case,i,j,k
  !double precision  :: ux(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  !double precision  :: uy(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  !double precision  :: uz(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: pf_old(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: pf_new(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: uc_old(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: dx,dt,thers,DM,p_k,flux,flux_flw

  ! local variables
  double precision  :: dfan_dx,dfan_dy,dfan_dz,dfan_dt
  double precision  :: fan_c,fan_p,Urel,mod_c,Js,Ja,fan_p2
  !double precision  :: fan_flw = 0.d0

  dfan_dx   = 0.d0
  dfan_dy   = 0.d0
  dfan_dz   = 0.d0
  dfan_dt   = 0.d0
  fan_c     = 0.d0
  fan_p     = 0.d0
  Urel      = 0.d0
  mod_c     = 0.d0
  Js        = 0.d0
  Ja        = 0.d0
  fan_p2    = 0.d0

  flux = 0.d0
  flux_flw = 0.d0

  select case (s_case)
    case (1)
      dfan_dx = ( pf_old(i,j,k) - pf_old(i-1,j,k) ) / dx
      fan_p   = 0.5d0 * (pf_old(i,j,k) + pf_old(i-1,j,k))
      Urel    = 0.5d0 * (uc_old(i,j,k) + uc_old(i-1,j,k))
      Js      = DM * 0.5 * (1.d0 - fan_p) * (uc_old(i,j,k) - uc_old(i-1,j,k)) / dx

      !---------------------------------------The flow part------------------------------------------------
      !fan_flw = 0.5d0 * (ux(i,j,k) + ux(i-1,j,k)) * fan_p * flw_amp
      !Js      = Js - 0.5d0 * fan_flw * ( 1.d0+p_k-(1.d0-p_k)*fan_p ) * (uc_old(i,j,k) - uc_old(i-1,j,k))
       
       !fan_p2   = 0.5d0 * (pf_new(i,j,k) + pf_new(i-1,j,k))
       !flux_flw = (1.d0+p_k-(1.d0-p_k)*fan_p2)*(1.d0+(1.d0-p_k)*Urel)
      !----------------------------------------------------------------------------------------------------

      if (abs(dfan_dx) <= thers) then
          Ja = 0.d0
      else
          ! calculate related terms
          dfan_dy = ( pf_old(i,j+1,k) + pf_old(i-1,j+1,k) - pf_old(i,j-1,k) - pf_old(i-1,j-1,k) ) / (4.d0 * dx)
          dfan_dz = ( pf_old(i,j,k+1) + pf_old(i-1,j,k+1) - pf_old(i,j,k-1) - pf_old(i-1,j,k-1) ) / (4.d0 * dx)
          fan_c   = 0.5d0 * (pf_new(i,j,k) + pf_new(i-1,j,k))
          !Urel    = 0.5d0 * (uc_old(i,j,k) + uc_old(i-1,j,k))
          mod_c = sqrt(dfan_dx**2 + dfan_dy**2 + dfan_dz**2)
          dfan_dt = (fan_c - fan_p) / dt

          ! calculate the antitrapping current
          Ja = sqrt(2.d0) / 4.d0 * ( 1.d0 + ( 1.d0 - p_k) * Urel) * dfan_dt * dfan_dx / mod_c

          !---------------------------------------The flow part------------------------------------------------
          !Ja = Ja + 0.5d0 * fan_flw * (1.d0+(1.d0-p_k)*Urel)* ( pf_old(i,j,k) - pf_old(i-1,j,k) )
          !----------------------------------------------------------------------------------------------------

      end if

      flux = Js + Ja
      !flux = Js

  case (2)
    dfan_dy = ( pf_old(i,j,k) - pf_old(i,j-1,k) ) / dx
    fan_p   = 0.5d0 * (pf_old(i,j,k) + pf_old(i,j-1,k))
    Urel    = 0.5d0 * (uc_old(i,j,k) + uc_old(i,j-1,k))
    Js      = DM * 0.5 * (1.d0 - fan_p) * (uc_old(i,j,k) - uc_old(i,j-1,k)) / dx

    !---------------------------------------The flow part------------------------------------------------
    !fan_flw = 0.5d0 * (uy(i,j,k) + uy(i,j-1,k)) * fan_p *flw_amp
    !Js      = Js - 0.5d0 * fan_flw * ( 1.d0+p_k-(1.d0-p_k)*fan_p ) * (uc_old(i,j,k) - uc_old(i,j-1,k))
     !fan_p2   = 0.5d0 * (pf_new(i,j,k) + pf_new(i,j-1,k))
     !flux_flw = (1.d0+p_k-(1.d0-p_k)*fan_p2)*(1.d0+(1.d0-p_k)*Urel)
    !----------------------------------------------------------------------------------------------------

    if (abs(dfan_dy) <= thers) then
        Ja = 0.d0
    else
        ! calculate related terms
        dfan_dx = ( pf_old(i+1,j,k) + pf_old(i+1,j-1,k) - pf_old(i-1,j,k) - pf_old(i-1,j-1,k) ) / (4.d0 * dx)
        dfan_dz = ( pf_old(i,j,k+1) + pf_old(i,j-1,k+1) - pf_old(i,j,k-1) - pf_old(i,j-1,k-1) ) / (4.d0 * dx)
        fan_c   = 0.5d0 * (pf_new(i,j,k) + pf_new(i,j-1,k))
        !Urel    = 0.5d0 * (uc_old(i,j,k) + uc_old(i,j-1,k))
        mod_c = sqrt(dfan_dx**2 + dfan_dy**2 + dfan_dz**2)
        dfan_dt = (fan_c - fan_p) / dt

        ! calculate the antitrapping current
        Ja = sqrt(2.d0) / 4.d0 * ( 1.d0 + ( 1.d0 - p_k) *Urel) * dfan_dt * dfan_dy / mod_c

        !---------------------------------------The flow part------------------------------------------------
        !Ja = Ja + 0.5d0 * fan_flw * (1.d0+(1.d0-p_k)*Urel)* ( pf_old(i,j,k) - pf_old(i,j-1,k) )
        !----------------------------------------------------------------------------------------------------

    end if

    flux = Js + Ja
    !flux = Js

  case (3)
    dfan_dz = ( pf_old(i,j,k) - pf_old(i,j,k-1) ) / dx
    fan_p   = 0.5d0 * (pf_old(i,j,k) + pf_old(i,j,k-1))
    Urel    = 0.5d0 * (uc_old(i,j,k) + uc_old(i,j,k-1))
    Js      = DM * 0.5 * (1.d0 - fan_p) * (uc_old(i,j,k) - uc_old(i,j,k-1)) / dx

    !---------------------------------------The flow part------------------------------------------------
    !fan_flw = 0.5d0 * (uz(i,j,k) + uz(i,j,k-1)) * fan_p * flw_amp
    !Js      = Js - 0.5d0 * fan_flw * ( 1.d0+p_k-(1.d0-p_k)*fan_p ) * (uc_old(i,j,k) - uc_old(i,j,k-1))
     !fan_p2   = 0.5d0 * (pf_new(i,j,k) + pf_new(i,j,k-1))
     !flux_flw = (1.d0+p_k-(1.d0-p_k)*fan_p2)*(1.d0+(1.d0-p_k)*Urel)
    !----------------------------------------------------------------------------------------------------

    if (abs(dfan_dz) <= thers) then
        Ja = 0.d0
    else
        !calculate related terms
        dfan_dx = ( pf_old(i+1,j,k) + pf_old(i+1,j,k-1) - pf_old(i-1,j,k) - pf_old(i-1,j,k-1) ) / (4.d0 * dx)
        dfan_dy = ( pf_old(i,j+1,k) + pf_old(i,j+1,k-1) - pf_old(i,j-1,k) - pf_old(i,j-1,k-1) ) / (4.d0 * dx)
        fan_c   = 0.5d0 * (pf_new(i,j,k) + pf_new(i,j,k-1))
        !Urel    = 0.5d0 * (uc_old(i,j,k) + uc_old(i,j,k-1))
        mod_c = sqrt(dfan_dx**2 + dfan_dy**2 + dfan_dz**2)
        dfan_dt = (fan_c - fan_p) / dt

        ! calculate the antitrapping current
        Ja = sqrt(2.d0) / 4.d0 * ( 1.d0 + ( 1.d0 - p_k) *Urel) * dfan_dt * dfan_dz / mod_c

        !---------------------------------------The flow part------------------------------------------------
        !Ja = Ja + 0.5d0 * fan_flw * (1.d0+(1.d0-p_k)*Urel)* ( pf_old(i,j,k) - pf_old(i,j,k-1) )
        !----------------------------------------------------------------------------------------------------

    end if

    flux = Js + Ja
    !flux = Js

  end select

end subroutine cal_sur_sol_flux_3d_flw

end module solute_field_module
