!
! This module is designed to work with the phase field
!
module phase_field_module

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

  public :: advance_phase_field

contains

  subroutine advance_phase_field(mla,ad_ori,phi_old,phi_new,dx,dt,the_bc_tower,pfpara,noise)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: ad_ori(:),phi_old(:),phi_new(:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: dt
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(pf_para),   intent(in   ) :: pfpara
    type(multifab) , intent(in   ) :: noise(:)

    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: dm, ng_p, i, n,nlevs
    real(kind=dp_t) :: var_h, var_l,var_dt
    real(kind=dp_t) :: rkm_con1, rkm_con2, rkm_con3

    real(kind=dp_t), pointer ::  pfo(:,:,:,:)
    real(kind=dp_t), pointer ::  pfn(:,:,:,:)
    real(kind=dp_t), pointer ::  noi(:,:,:,:)
    real(kind=dp_t), pointer ::  por(:,:,:,:)
    real(kind=dp_t), pointer ::  pfb(:,:,:,:)

    rkm_con1 = 2.d0/3.d0; rkm_con2 = 0.25d0; rkm_con3 = 0.75d0;
    dm   = mla%dim
    nlevs = mla%nlevel
    ng_p = phi_old(1)%ng

    if(pfpara%do_solve_flow == 1 .and. pfpara%solve_flow_only == 1) then
       do n=1,nlevs
           call multifab_copy(phi_new(n), phi_old(n), 1)
       end do
       
       return
    end if

    var_h = 1.d0
    var_l = -1.d0
    var_dt = rkm_con1 * dt

    ! calculate the RK matrices
    do n=1,nlevs
       do i=1,nfabs(phi_old(n))
          por  => dataptr(ad_ori(n),i)
          pfo  => dataptr(phi_old(n),i)
          pfn  => dataptr(phi_new(n),i)
          noi  => dataptr(noise(n), i)
          lo = lwb(get_box(phi_old(n),i))
          hi = upb(get_box(phi_old(n),i))

          select case(dm)
          case (2)
          ! calculate the K1
          call cal_rkmat_pf_2d(pfn(:,:,1,1),pfo(:,:,1,1),por(:,:,1,1),pfo(:,:,1,2),&
                                  pfo(:,:,1,3),noi(:,:,1,1),ng_p,lo,hi,dx(n),pfpara%thers,&
                                  pfpara%kk,pfpara%Le,pfpara%Mc00,pfpara%sk,pfpara%lamda,&
                                  pfpara%anis,pfpara%R_M,dt,-1.d0,1.d0 )

          case (3)
          call cal_rkmat_pf_3d(pfn(:,:,:,1),pfo(:,:,:,1),por(:,:,:,1),pfo(:,:,:,2),&
                                  pfo(:,:,:,3),noi(:,:,:,1),ng_p,lo,hi,dx(n),pfpara%thers,&
                                  pfpara%kk,pfpara%Le,pfpara%Mc00,pfpara%sk,pfpara%lamda,&
                                  pfpara%anis,pfpara%anis_0,pfpara%R_M,dt,-1.d0,1.d0 )

          end select

       end do
    end do

    call ml_restrict_and_fill(nlevs, phi_new, mla%mba%rr, the_bc_tower%bc_tower_array,1,1,1)
    call ml_restrict_and_fill(nlevs, ad_ori, mla%mba%rr, the_bc_tower%bc_tower_array,1,1,1)
    
    ! update the new phase field and orientaion field
    !if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
    !   call multifab_fill_boundary_c(phi_new(nlevs),1,1)

       ! fill non-periodic domain boundary ghost cells
    !   call multifab_physbc(phi_new(nlevs),1,1,1,the_bc_tower%bc_tower_array(nlevs))

       ! Do the same for the orientation field
    !   call multifab_fill_boundary_c(ad_ori(nlevs),1,1)
    !   call multifab_physbc(ad_ori(nlevs),1,1,1,the_bc_tower%bc_tower_array(nlevs))

    !else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
    !   do n=nlevs,2,-1
          ! set level n-1 data to be the average of the level n data covering it
    !      call ml_cc_restriction_c(phi_new(n-1),1,phi_new(n),1,mla%mba%rr(n-1,:),1)

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
    !      call multifab_fill_ghost_cells(phi_new(n),phi_new(n-1),phi_new(n)%ng,mla%mba%rr(n-1,:), &
    !                                     the_bc_tower%bc_tower_array(n-1), &
    !                                     the_bc_tower%bc_tower_array(n), &
    !                                     1,1,1)

          ! Do the same for orientation field
    !      call ml_cc_restriction(ad_ori(n-1),ad_ori(n),mla%mba%rr(n-1,:))
    !      call multifab_fill_ghost_cells(ad_ori(n),ad_ori(n-1),ad_ori(n)%ng,mla%mba%rr(n-1,:), &
    !                                     the_bc_tower%bc_tower_array(n-1), &
    !                                     the_bc_tower%bc_tower_array(n), &
    !                                     1,1,1)

    !   end do

    !end if


  end subroutine advance_phase_field

! subroutine calculate_rkmatrices_pf_2d(...) aims to calculate the RK matrices for phase field in 3-D
! input variable:
!            phi    :    phase field
!            ad     :    predefined orientation, should be calculated
!            uc     :    the solute field, matrix
!            th     :    the temperature field, matrix
!            noise  :    the noise matrix
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
!            anis   :    the anisotropy strength
! output variables:
!              K_R    :    RK matrices

subroutine cal_rkmat_pf_2d(K_R,phi,ad,uc,th,noise,ng_p,lo,hi,dx,thers,kk,Le,Mc00,&
                           p_k,lamda,anis,R_M,dt,var_l,var_h)

  integer           :: lo(2),hi(2),ng_p
  double precision  :: phi(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: ad(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: uc(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: th(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: noise(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: K_R(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: dx,dt,var_l,var_h,thers,kk,Le,Mc00,p_k,lamda,anis,R_M(:,:,:)

  ! local variables
  integer           :: i,j, x_c,y_c,sd_id
  double precision  :: de1X,de1Y,A0,ebar,mod_all,A_F0,fA,BT,pfs,JR,JL,JT,JB

  ! set the values for these commonly used
  ! Currently the values are all for cubic geometry
  A0 = 1.d0 - 3.d0 * anis
  ebar = 4.d0 * anis / A0
  x_c = 1;y_c = 2;sd_id = 1;

  ! do the loops
  !$omp parallel do private(j,i,de1X,de1Y,fA,BT,JR,JL,JT,JB,mod_all,pfs)
  do j=lo(2),hi(2)
    do i=lo(1),hi(1)

      de1X = (phi(i+1,j) - phi(i-1,j))/(2.d0 * dx)
      de1Y = (phi(i,j+1) - phi(i,j-1))/(2.d0 * dx)

      !fA = 1.d0/Le + Mc00*(1.d0+(1.d0-p_k)*uc(i,j))
      fA = Mc00* ( 1.d0+(1.d0-p_k)*uc(i,j) )
      BT = 1.d0

      !Set the fluxes to be zero
      JR = 0.d0
      JL = 0.d0
      JT = 0.d0
      JB = 0.d0

      mod_all = sqrt(de1X**2 + de1Y**2)

      if (mod_all > thers) then
          ! Retrieve the orientation of the transition cell, i,j
          call retriveMax_2d_R(sd_id,phi,ad,ng_p,i,j,lo,hi)

          ! Calculate the anisotropy factor
          call calculate_rotate_roi_2d_R2(A_F0,de1X,de1Y,R_M(sd_id,:,:),A0,ebar)

          ! Calculate the BT terms
          BT = A_F0 * A_F0

          ! Calculate the surrounding fluxes
          call cal_sur_pf_flux_2d_R(x_c,JL,phi,i,  j,   ng_p,lo,hi,kk,R_M(sd_id,:,:),dx,thers,A0,ebar)
          call cal_sur_pf_flux_2d_R(x_c,JR,phi,i+1,j,   ng_p,lo,hi,kk,R_M(sd_id,:,:),dx,thers,A0,ebar)
          call cal_sur_pf_flux_2d_R(y_c,JB,phi,i,  j,   ng_p,lo,hi,kk,R_M(sd_id,:,:),dx,thers,A0,ebar)
          call cal_sur_pf_flux_2d_R(y_c,JT,phi,i,  j+1, ng_p,lo,hi,kk,R_M(sd_id,:,:),dx,thers,A0,ebar)

      end if

      !pfs = lamda * (th(i,j) *(1.d0 + noise(i,j)) + Mc00 * uc(i,j))
      pfs = lamda * ( Mc00* ( th(i,j) +  uc(i,j) )  + noise(i,j)  )

      K_R(i,j) = ( JR - JL + JT - JB ) / dx + &
               phi(i,j) * ( 1.d0-phi(i,j) * phi(i,j) )-( ( 1.d0 - phi(i,j) * phi(i,j) )**2 ) * pfs
      K_R(i,j) = K_R(i,j)/(fA * BT)

      K_R(i,j) = phi(i,j) + dt * K_R(i,j)

      ! cut off the bouds
      !if (K_R(i,j) > var_h) then
      !   K_R(i,j) = var_h
      !else if(K_R(i,j) < var_l) then
      !   K_R(i,j) = var_l
      !end if

    end do
  end do
  !$omp end parallel do

end subroutine cal_rkmat_pf_2d

! subroutine calculate_rkmatrices_pf_3d(...) aims to calculate the RK matrices for phase field in 3-D
! input variable:
!            phi    :    phase field
!            ad     :    predefined orientation, should be calculated
!            uc     :    the solute field, matrix
!            th     :    the temperature field, matrix
!            noise  :    the noise matrix
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
!            anis   :    the anisotropy strength
! output variables:
!              K_R    :    RK matrices

subroutine cal_rkmat_pf_3d(K_R,phi,ad,uc,th,noise,ng_p,lo,hi,dx,thers,kk,Le,Mc00,&
                           p_k,lamda,anis,anis_0,R_M,dt,var_l,var_h)

  integer           :: lo(3),hi(3),ng_p
  double precision  :: phi(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: ad(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: uc(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: th(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: noise(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: K_R(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: dx,dt,var_l,var_h,thers,kk,Le,Mc00,p_k,lamda,anis,anis_0,R_M(:,:,:)

  ! local variables
  integer           :: i,j,k,x_c,y_c,z_c,sd_id
  double precision  :: de1X,de1Y,de1Z,A_F0,fA,BT,pfs,JR,JL,JT,JB,JF,JK
  double precision  :: A0,ebar,mod_all

  ! set the values for these commonly used
  ! Currently the values are all for cubic geometry
  !A0 = 1.d0 - 3.d0 * anis
  !ebar = 4.d0 * anis / A0
  A0 = anis
  ebar = anis_0
  x_c = 1;y_c = 2;z_c = 3;sd_id = 1;

  ! do the loops
  !$omp parallel do private(k,j,i,de1X,de1Y,de1Z,fA,BT,JF,JK,JR,JL,JT,JB,mod_all,pfs)
  do k=lo(3),hi(3)
    do j=lo(2),hi(2)
      do i=lo(1),hi(1)

        de1X = (phi(i+1,j,k) - phi(i-1,j,k))/(2.d0 * dx)
        de1Y = (phi(i,j+1,k) - phi(i,j-1,k))/(2.d0 * dx)
        de1Z = (phi(i,j,k+1) - phi(i,j,k-1))/(2.d0 * dx)

        !fA = 1.d0/Le + Mc00*(1.d0+(1.d0-p_k)*uc(i,j))
        fA = Mc00* ( 1.d0+(1.d0-p_k)*uc(i,j,k) )
        BT = 1.d0

        !Set the fluxes to be zero
        JF = 0.d0
        JK = 0.d0
        JL = 0.d0
        JR = 0.d0
        JT = 0.d0
        JB = 0.d0

        mod_all = sqrt(de1X**2 + de1Y**2 + de1Z**2)

        if (mod_all > thers) then

          ! Retrieve the seed_index of the transition cell, i,j
          call retriveMax_3d_R(sd_id,phi,ad,ng_p,i,j,k,lo,hi)

          ! Calculate the anisotropy factor
          call calculate_rotate_roi_3d_R2(A_F0,de1X,de1Y,de1Z,R_M(sd_id,:,:),A0,ebar)

          ! Calculate the BT terms
          BT = A_F0 * A_F0

          ! Calculate the surrounding fluxes
          call cal_sur_pf_flux_3d_R(x_c,JL,phi,i,  j,  k,  ng_p,lo,hi,kk,R_M(sd_id,:,:),dx,thers,A0,ebar)
          call cal_sur_pf_flux_3d_R(x_c,JR,phi,i+1,j,  k,  ng_p,lo,hi,kk,R_M(sd_id,:,:),dx,thers,A0,ebar)
          call cal_sur_pf_flux_3d_R(y_c,JB,phi,i,  j,  k,  ng_p,lo,hi,kk,R_M(sd_id,:,:),dx,thers,A0,ebar)
          call cal_sur_pf_flux_3d_R(y_c,JT,phi,i,  j+1,k,  ng_p,lo,hi,kk,R_M(sd_id,:,:),dx,thers,A0,ebar)
          call cal_sur_pf_flux_3d_R(z_c,JK,phi,i,  j,  k,  ng_p,lo,hi,kk,R_M(sd_id,:,:),dx,thers,A0,ebar)
          call cal_sur_pf_flux_3d_R(z_c,JF,phi,i  ,j,  k+1,ng_p,lo,hi,kk,R_M(sd_id,:,:),dx,thers,A0,ebar)

        end if

        !pfs = lamda * (th(i,j) *(1.d0 + noise(i,j)) + Mc00 * uc(i,j))
        pfs = lamda * ( Mc00* ( th(i,j,k)+  uc(i,j,k) )  + noise(i,j,k)  )

        K_R(i,j,k) = ( JR - JL + JT - JB + JF - JK ) / dx + &
                     phi(i,j,k) * ( 1.d0-phi(i,j,k) * phi(i,j,k) )-( ( 1.d0 - phi(i,j,k) * phi(i,j,k) )**2 ) * pfs
        K_R(i,j,k) = K_R(i,j,k)/(fA * BT)

        K_R(i,j,k) = phi(i,j,k) + dt * K_R(i,j,k)

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

end subroutine cal_rkmat_pf_3d

! this routine could only be called when de1X and de1Y are not equal to zero at the same time
subroutine calculate_rotate_roi_2d_R(A_F0,A_F1_x,A_F1_y,de1X,de1Y,R_M,A0,ebar,thers)

  double precision :: A_F0,A_F1_x,A_F1_y,de1X,de1Y,A0,ebar,thers
  double precision :: R_M(:,:)

  ! local variables
  double precision :: sq_mod, de1X_p, de1Y_p, anis_E

  sq_mod = de1X**2 + de1Y**2

  if(sqrt(sq_mod) <= thers) then
     A_F0 = 1.d0
     A_F1_x = 0.d0
     A_F1_y = 0.d0

     return

  end if

  ! calculate the derivative with prime
  de1X_p = R_M(1,1) * de1X + R_M(1,2) * de1Y
  de1Y_p = R_M(2,1) * de1X + R_M(2,2) * de1Y

  anis_E = (de1X_p**4 + de1Y_p**4 )/(sq_mod**2)

  ! calculate the A_F0
  A_F0 = A0*( 1.d0 + ebar* anis_E)

  ! calculate the A_F1_x, A_F1_y, A_F1_z
  A_F1_x = 4.d0 * A0 * ebar * de1X_p * ( de1X_p**2/sq_mod - anis_E)
  A_F1_y = 4.d0 * A0 * ebar * de1Y_p * ( de1Y_p**2/sq_mod - anis_E)

end subroutine calculate_rotate_roi_2d_R

subroutine calculate_rotate_roi_2d_R2(A_F0,de1X,de1Y,R_M,A0,ebar)

  double precision :: A_F0,de1X,de1Y,A0,ebar
  double precision :: R_M(:,:)

  ! local variables
  double precision :: sq_mod, de1X_p, de1Y_p,  anis_E

  sq_mod = de1X**2 + de1Y**2

  ! calculate the derivative with prime
  de1X_p = R_M(1,1) * de1X + R_M(1,2) * de1Y
  de1Y_p = R_M(2,1) * de1X + R_M(2,2) * de1Y

  anis_E = (de1X_p**4 + de1Y_p**4)/(sq_mod**2)

  ! calculate the A_F0
  A_F0 = A0*( 1.d0 + ebar* anis_E)

end subroutine calculate_rotate_roi_2d_R2

subroutine cal_sur_pf_flux_2d_R(s_case,flux,phi,i,j,ng_p,lo,hi,kk,R_M,dx,thers,A0,ebar)

  integer             :: lo(2),hi(2),ng_p,i,j,s_case
  double precision    :: flux,kk,dx,thers,A0,ebar,R_M(2,2)
  double precision    :: phi(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)

  ! local variables
  double precision    :: de1X,de1Y,de1X_p,de1Y_p
  double precision    :: A_F0,A_F1_x,A_F1_y

  ! setup initial values
  flux = 0.d0

  select case (s_case)
    case (1)

      ! calculate the first derivatives
      de1X = ( phi(i,j) - phi(i-1,j) ) / dx
      de1Y = ( phi(i,j+1) + phi(i-1,j+1) - phi(i,j-1) - phi(i-1,j-1) ) / (4.d0 * dx)

      ! calculate the A_F0, A_F1 terms
      call calculate_rotate_roi_2d_R(A_F0,A_F1_x,A_F1_y,de1X,de1Y,R_M,A0,ebar,thers)

      ! calculate the flux
      flux = A_F0*( A_F0*de1X + R_M(1,1)*A_F1_x + R_M(2,1)*A_F1_y )

    case (2)

      ! calculate the first derivatives
      de1X = ( phi(i+1,j) + phi(i+1,j-1) - phi(i-1,j) - phi(i-1,j-1) ) / (4.d0 * dx)
      de1Y = ( phi(i,j) - phi(i,j-1) ) / dx

      ! calculate the A_F0, A_F1 terms
      call calculate_rotate_roi_2d_R(A_F0,A_F1_x,A_F1_y,de1X,de1Y,R_M,A0,ebar,thers)

      ! calculate the flux
      flux = A_F0*(A_F0*de1Y + R_M(1,2)*A_F1_x + R_M(2,2)*A_F1_y )

  end select

end subroutine cal_sur_pf_flux_2d_R

subroutine retriveMax_2d_R(seed_index,phi,ad,ng_p,i_d,j_d,lo,hi)

integer           :: i_d,j_d,ng_p,lo(2),hi(2),seed_index
double precision  :: phi(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
double precision  :: ad(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)

!local variables
integer           :: i,j
double precision  :: pfmax

! set the pfmax variable to be itself
pfmax = phi(i_d,j_d)

! do the loops
do j=j_d-1,j_d+1
   do i=i_d-1,i_d+1
      if(phi(i,j) > pfmax) then
         pfmax = phi(i,j)
         ad(i_d,j_d) = ad(i,j)
      end if
   end do
end do

! set the related values
seed_index = int( ad(i_d,j_d) )

if(seed_index <= 0.d0) seed_index = 1

end subroutine retriveMax_2d_R

! this routine could only be called when de1X and de1Y are not equal to zero at the same time
subroutine calculate_rotate_roi_3d_R_BK(A_F0,A_F1_x,A_F1_y,A_F1_z,de1X,de1Y,de1Z,R_M,A0,ebar,thers)

  double precision :: A_F0,A_F1_x,A_F1_y,A_F1_z,de1X,de1Y,de1Z,A0,ebar,thers
  double precision :: R_M(:,:)

  ! local variables
  double precision :: sq_mod, de1X_p, de1Y_p, de1Z_p, anis_E

  sq_mod = de1X**2 + de1Y**2 + de1Z**2

  if(sqrt(sq_mod) <= thers) then
     A_F0 = 1.d0
     A_F1_x = 0.d0
     A_F1_y = 0.d0
     A_F1_z = 0.d0

     return

  end if

  ! calculate the derivative with prime
  de1X_p = R_M(1,1) * de1X + R_M(1,2) * de1Y + R_M(1,3) * de1Z
  de1Y_p = R_M(2,1) * de1X + R_M(2,2) * de1Y + R_M(2,3) * de1Z
  de1Z_p = R_M(3,1) * de1X + R_M(3,2) * de1Y + R_M(3,3) * de1Z

  anis_E = (de1X_p**4 + de1Y_p**4 + de1Z_p**4)/(sq_mod**2)

  ! calculate the A_F0
  A_F0 = A0*( 1.d0 + ebar* anis_E)

  ! calculate the A_F1_x, A_F1_y, A_F1_z
  A_F1_x = 4.d0 * A0 * ebar * de1X_p * ( de1X_p**2/sq_mod - anis_E)
  A_F1_y = 4.d0 * A0 * ebar * de1Y_p * ( de1Y_p**2/sq_mod - anis_E)
  A_F1_z = 4.d0 * A0 * ebar * de1Z_p * ( de1Z_p**2/sq_mod - anis_E)

end subroutine calculate_rotate_roi_3d_R_BK

subroutine calculate_rotate_roi_3d_R2_BK(A_F0,de1X,de1Y,de1Z,R_M,A0,ebar)

  double precision :: A_F0,de1X,de1Y,de1Z,A0,ebar
  double precision :: R_M(:,:)

  ! local variables
  double precision :: sq_mod, de1X_p, de1Y_p, de1Z_p, anis_E

  sq_mod = de1X**2 + de1Y**2 + de1Z**2

  ! calculate the derivative with prime
  de1X_p = R_M(1,1) * de1X + R_M(1,2) * de1Y + R_M(1,3) * de1Z
  de1Y_p = R_M(2,1) * de1X + R_M(2,2) * de1Y + R_M(2,3) * de1Z
  de1Z_p = R_M(3,1) * de1X + R_M(3,2) * de1Y + R_M(3,3) * de1Z

  anis_E = (de1X_p**4 + de1Y_p**4 + de1Z_p**4)/(sq_mod**2)

  ! calculate the A_F0
  A_F0 = A0*( 1.d0 + ebar* anis_E)

end subroutine calculate_rotate_roi_3d_R2_BK

! this routine could only be called when de1X and de1Y are not equal to zero at the same time
subroutine calculate_rotate_roi_3d_R(A_F0,A_F1_x,A_F1_y,A_F1_z,de1X,de1Y,de1Z,R_M,anis_1,anis_2,thers)

  double precision :: A_F0,A_F1_x,A_F1_y,A_F1_z,de1X,de1Y,de1Z,anis_1,anis_2,thers
  double precision :: R_M(:,:)

  ! local variables
  double precision :: sq_mod, de1X_p, de1Y_p, de1Z_p, anis_E, anis_S, anis_E6
  double precision :: A0,ebar_1,ebar_2, sq_d1X, sq_d1Y, sq_d1Z,tri_mod,sq_sqmod
  double precision :: p_Q6, p_Qz, p_Qx, p_Qy

  ebar_1 = anis_1
  ebar_2 = anis_2

  !A0  = 1.d0 - 3.d0 * anis_1
  !ebar_1 = 4.d0 * anis_1 / A0
  !ebar_2 = anis_2 / A0

  sq_mod = de1X**2 + de1Y**2 + de1Z**2

  if(sqrt(sq_mod) <= thers) then
     A_F0 = 1.d0
     A_F1_x = 0.d0
     A_F1_y = 0.d0
     A_F1_z = 0.d0

     return

  end if

  ! calculate the derivative with prime
  de1X_p = R_M(1,1) * de1X + R_M(1,2) * de1Y + R_M(1,3) * de1Z
  de1Y_p = R_M(2,1) * de1X + R_M(2,2) * de1Y + R_M(2,3) * de1Z
  de1Z_p = R_M(3,1) * de1X + R_M(3,2) * de1Y + R_M(3,3) * de1Z

  sq_d1X = de1X_p**2
  sq_d1Y = de1Y_p**2
  sq_d1Z = de1Z_p**2
  tri_mod = sq_mod**3
  sq_sqmod = sq_mod**2

  anis_E = (de1X_p**4 + de1Y_p**4 + de1Z_p**4)/sq_sqmod
  anis_S = sq_d1X * sq_d1Y * sq_d1Z / tri_mod
  !anis_E6 = (de1X_p**6 + de1Y_p**6 + de1Z_p**6)/tri_mod

  ! calculate the A_F0
  ! anisotropy functional - type 1
  A_F0 = 1.d0 + ebar_1 * (anis_E - 0.6d0) + &
         ebar_2*(3.d0 * anis_E + 66.d0 * anis_S - 17.d0/7.d0)


  ! calculate the A_F1_x, A_F1_y, A_F1_z
  A_F1_x = (ebar_1 + 3.d0 * ebar_2) * 4.d0 *de1X_p *( sq_d1X/sq_mod - anis_E) &
           + 132.d0 * ebar_2 * de1X_p * sq_d1Y * sq_d1Z * &
           (sq_mod - 3.d0* sq_d1X) / tri_mod

  A_F1_y = (ebar_1 + 3.d0 * ebar_2) * 4.d0 *de1Y_p *( sq_d1Y/sq_mod - anis_E) &
           + 132.d0 * ebar_2 * de1Y_p * sq_d1X * sq_d1Z * &
           (sq_mod - 3.d0* sq_d1Y) / tri_mod

  A_F1_z = (ebar_1 + 3.d0 * ebar_2) * 4.d0 *de1Z_p *( sq_d1Z/sq_mod - anis_E) &
           + 132.d0 * ebar_2 * de1Z_p * sq_d1X * sq_d1Y * &
           (sq_mod - 3.d0* sq_d1Z) / tri_mod

  ! The anisotropy functional - type 2,
  !A_F0 = A0*(1.d0 + ebar_1*anis_E+ebar_2*( anis_E6 + 30.d0*anis_S) )

  ! The three derivatives
  !A_F1_x = 4.d0 * A0 * ebar_1 * de1X_p * ( sq_d1X/sq_mod - anis_E) + &
  !         6.d0 * A0 * ebar_2 * de1X_p * ( de1X_p**4/sq_sqmod - anis_E6) + &
  !         60.d0* A0 * ebar_2 * de1X_p * sq_d1Y * sq_d1Z * &
  !         (sq_mod - 3.d0* sq_d1X) / tri_mod

  !A_F1_y = 4.d0 * A0 * ebar_1 * de1Y_p * ( sq_d1Y/sq_mod - anis_E) + &
  !         6.d0 * A0 * ebar_2 * de1Y_p * ( de1Y_p**4/sq_sqmod - anis_E6) + &
  !         60.d0* A0 * ebar_2 * de1Y_p * sq_d1X * sq_d1Z * &
  !         (sq_mod - 3.d0* sq_d1Y) / tri_mod

  !A_F1_z = 4.d0 * A0 * ebar_1 * de1Z_p * ( sq_d1Z/sq_mod - anis_E) + &
  !         6.d0 * A0 * ebar_2 * de1Z_p * ( de1Z_p**4/sq_sqmod - anis_E6) + &
  !         60.d0* A0 * ebar_2 * de1Z_p * sq_d1X * sq_d1Y * &
  !         (sq_mod - 3.d0* sq_d1Z) / tri_mod

end subroutine calculate_rotate_roi_3d_R

! this routine could only be called when de1X and de1Y are not equal to zero at the same time
subroutine calculate_rotate_roi_3d_R2(A_F0,de1X,de1Y,de1Z,R_M,anis_1,anis_2)

  double precision :: A_F0,de1X,de1Y,de1Z,anis_1,anis_2
  double precision :: R_M(:,:)

  ! local variables
  double precision :: sq_mod, de1X_p, de1Y_p, de1Z_p, anis_E, anis_S,anis_E6
  double precision :: A0,ebar_1,ebar_2

  ebar_1 = anis_1
  ebar_2 = anis_2

  !A0  = 1.d0 - 3.d0 * anis_1
  !ebar_1 = 4.d0 * anis_1 / A0
  !ebar_2 = anis_2 / A0

  sq_mod = de1X**2 + de1Y**2 + de1Z**2

  ! calculate the derivative with prime
  de1X_p = R_M(1,1) * de1X + R_M(1,2) * de1Y + R_M(1,3) * de1Z
  de1Y_p = R_M(2,1) * de1X + R_M(2,2) * de1Y + R_M(2,3) * de1Z
  de1Z_p = R_M(3,1) * de1X + R_M(3,2) * de1Y + R_M(3,3) * de1Z

  anis_E = (de1X_p**4 + de1Y_p**4 + de1Z_p**4)/(sq_mod**2)
  anis_S = (de1X_p**2) * (de1Y_p**2) * (de1Z_p**2) / (sq_mod**3)
  anis_E6 = (de1X_p**6 + de1Y_p**6 + de1Z_p**6)/(sq_mod**3)

  ! calculate the A_F0
  A_F0 = 1.d0 + ebar_1 * (anis_E - 0.6d0) + &
         ebar_2*(3.d0 * anis_E + 66.d0 * anis_S - 17.d0/7.d0)

  !A_F0 = A0*(1.d0 + ebar_1*anis_E+ebar_2*( anis_E6 + 30.d0*anis_S) )

end subroutine calculate_rotate_roi_3d_R2

! this routine could only be called when de1X and de1Y are not equal to zero at the same time
subroutine calculate_rotate_roi_3d_R_USE(A_F0,A_F1_x,A_F1_y,A_F1_z,de1X,de1Y,de1Z,R_M,ebar_1,ebar_2,thers)

  double precision :: A_F0,A_F1_x,A_F1_y,A_F1_z,de1X,de1Y,de1Z,ebar_1,ebar_2,thers
  double precision :: R_M(:,:)

  ! local variables
  double precision :: sq_mod, sq2_mod,de1X_p, de1Y_p, de1Z_p,sq_de1X,sq_de1Y,sq_de1Z
  double precision :: p_Q6, p_Qx, p_Qy, p_Qz, p_Qxyz

  sq_mod = de1X**2 + de1Y**2 + de1Z**2

  if(sqrt(sq_mod) <= thers) then
     A_F0 = 1.d0
     A_F1_x = 0.d0
     A_F1_y = 0.d0
     A_F1_z = 0.d0

     return

  end if

  ! calculate the derivative with prime
  de1X_p = R_M(1,1) * de1X + R_M(1,2) * de1Y + R_M(1,3) * de1Z
  de1Y_p = R_M(2,1) * de1X + R_M(2,2) * de1Y + R_M(2,3) * de1Z
  de1Z_p = R_M(3,1) * de1X + R_M(3,2) * de1Y + R_M(3,3) * de1Z

  ! Squares of the first derivatives
  sq_de1X = de1X_p**2
  sq_de1Y = de1Y_p**2
  sq_de1Z = de1Z_p**2

  ! A temporary variable, equals to |fan|^4
  sq2_mod = sq_mod**2

  p_Qxyz = ( sq_de1X**2+sq_de1Y**2-10.d0*sq_de1X*sq_de1Y )/ sq2_mod

  p_Q6 = ( sq_de1X**3-sq_de1Y**3-15.d0*sq_de1X*sq_de1Y*(sq_de1X - sq_de1Y) )/( sq_mod**3 )
  p_Qz = sq_de1Z / sq_mod
  p_Qx = p_Qxyz + 4.d0 * ( sq_de1Y**2 )/sq2_mod
  p_Qy = p_Qxyz + 4.d0 * ( sq_de1X**2 )/sq2_mod

  ! calculate the A_F0
  ! anisotropy functional - type 1
  A_F0 = 1.d0 + ebar_1 * p_Q6 + ebar_2 * p_Qz

  ! calculate the A_F1_x, A_F1_y, A_F1_z
  A_F1_x = 2.d0 * de1X_p * (  3.d0* ebar_1*( p_Qx-p_Q6) - ebar_2*p_Qz )

  A_F1_y = 2.d0 * de1Y_p * ( -3.d0* ebar_1*( p_Qy+p_Q6) - ebar_2*p_Qz )

  A_F1_z = 2.d0 * de1Z_p * ( -3.d0* ebar_1*p_Q6 + ebar_2* (1.d0 - p_Qz) )


end subroutine calculate_rotate_roi_3d_R_USE

! this routine could only be called when de1X and de1Y are not equal to zero at the same time
subroutine calculate_rotate_roi_3d_R2_USE(A_F0,de1X,de1Y,de1Z,R_M,ebar_1,ebar_2)

  double precision :: A_F0,de1X,de1Y,de1Z,ebar_1,ebar_2
  double precision :: R_M(:,:)

  ! local variables
  double precision :: sq_mod,de1X_p,de1Y_p,de1Z_p,sq_de1X,sq_de1Y,sq_de1Z
  double precision :: p_Q6,p_Qz

  sq_mod = de1X**2 + de1Y**2 + de1Z**2

  ! calculate the derivative with prime
  de1X_p = R_M(1,1) * de1X + R_M(1,2) * de1Y + R_M(1,3) * de1Z
  de1Y_p = R_M(2,1) * de1X + R_M(2,2) * de1Y + R_M(2,3) * de1Z
  de1Z_p = R_M(3,1) * de1X + R_M(3,2) * de1Y + R_M(3,3) * de1Z

  sq_de1X = de1X_p**2
  sq_de1Y = de1Y_p**2
  sq_de1Z = de1Z_p**2

  p_Q6 = ( sq_de1X**3 - sq_de1Y**3 - 15.d0*sq_de1X*sq_de1Y*(sq_de1X - sq_de1Y) )/( sq_mod**3 )
  p_Qz = sq_de1Z / sq_mod

  ! calculate the A_F0
  A_F0 = 1.d0 + ebar_1 * p_Q6 + ebar_2 * p_Qz

end subroutine calculate_rotate_roi_3d_R2_USE

subroutine cal_sur_pf_flux_3d_R(s_case,flux,phi,i,j,k,ng_p,lo,hi,kk,R_M,dx,thers,A0,ebar)

  integer             :: lo(3),hi(3),ng_p,i,j,k,s_case
  double precision    :: flux,kk,dx,thers,A0,ebar,R_M(3,3)
  double precision    :: phi(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)

  ! local variables
  double precision    :: de1X,de1Y,de1Z,de1X_p,de1Y_p,de1Z_p
  double precision    :: A_F0,A_F1_x,A_F1_y,A_F1_z,mod_all_sq

  ! setup initial values
  flux = 0.d0

  select case (s_case)
    case (1)

      ! calculate the first derivatives
      de1X = ( phi(i,j,k) - phi(i-1,j,k) ) / dx
      de1Y = ( phi(i,j+1,k) + phi(i-1,j+1,k) - phi(i,j-1,k) - phi(i-1,j-1,k) ) / (4.d0 * dx)
      de1Z = ( phi(i,j,k+1) + phi(i-1,j,k+1) - phi(i,j,k-1) - phi(i-1,j,k-1) ) / (4.d0 * dx)

      ! calculate the A_F0, A_F1 terms
      call calculate_rotate_roi_3d_R(A_F0,A_F1_x,A_F1_y,A_F1_z,de1X,de1Y,de1Z,R_M,A0,ebar,thers)

      ! calculate the flux
      flux = A_F0*( A_F0*de1X + R_M(1,1)*A_F1_x + R_M(2,1)*A_F1_y + R_M(3,1)*A_F1_z )

    case (2)

      ! calculate the first derivatives
      de1X = ( phi(i+1,j,k) + phi(i+1,j-1,k) - phi(i-1,j,k) - phi(i-1,j-1,k) ) / (4.d0 * dx)
      de1Y = ( phi(i,j,k) - phi(i,j-1,k) ) / dx
      de1Z = ( phi(i,j,k+1) + phi(i,j-1,k+1) - phi(i,j,k-1) - phi(i,j-1,k-1) ) / (4.d0 * dx)

      ! calculate the A_F0, A_F1 terms
      call calculate_rotate_roi_3d_R(A_F0,A_F1_x,A_F1_y,A_F1_z,de1X,de1Y,de1Z,R_M,A0,ebar,thers)

      ! calculate the flux
      flux = A_F0*(A_F0*de1Y + R_M(1,2)*A_F1_x + R_M(2,2)*A_F1_y + R_M(3,2)*A_F1_z )

    case (3)
      ! calculate the first derivatives
      de1X = ( phi(i+1,j,k) + phi(i+1,j,k-1) - phi(i-1,j,k) - phi(i-1,j,k-1) ) / (4.d0 * dx)
      de1Y = ( phi(i,j+1,k) + phi(i,j+1,k-1) - phi(i,j-1,k) - phi(i,j-1,k-1) ) / (4.d0 * dx)
      de1Z = ( phi(i,j,k) - phi(i,j,k-1) ) / dx

      ! calculate the A_F0, A_F1 terms
      call calculate_rotate_roi_3d_R(A_F0,A_F1_x,A_F1_y,A_F1_z,de1X,de1Y,de1Z,R_M,A0,ebar,thers)

      flux = A_F0*(A_F0*de1Z + R_M(1,3)*A_F1_x + R_M(2,3)*A_F1_y + R_M(3,3)*A_F1_z )

  end select

end subroutine cal_sur_pf_flux_3d_R

subroutine retriveMax_3d_R(seed_index,phi,ad,ng_p,i_d,j_d,k_d,lo,hi)

integer           :: i_d,j_d,k_d,ng_p,lo(3),hi(3),seed_index
double precision  :: phi(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
double precision  :: ad(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)

!local variables
integer           :: i,j,k
double precision  :: pfmax

! set the pfmax variable to be itself
pfmax = phi(i_d,j_d,k_d)

! do the loops
do k=k_d-1,k_d+1
   do j=j_d-1,j_d+1
     do i=i_d-1,i_d+1
         if(phi(i,j,k) > pfmax) then
            pfmax = phi(i,j,k)
            ad(i_d,j_d,k_d) = ad(i,j,k)
         end if
      end do
   end do
end do

! set the related values
seed_index = int( ad(i_d,j_d,k_d) )

if(seed_index <= 0.d0) seed_index = 1

end subroutine retriveMax_3d_R

end module phase_field_module




