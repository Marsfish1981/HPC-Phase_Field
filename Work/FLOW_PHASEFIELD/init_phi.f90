! the initilization module, used to initialize the phase field data
! Created by Z. Guo on Jan 17, 2014
! Last modified ib Jan 22, 2014

module init_phi_module

  use multifab_module
  use ml_layout_module
  use define_bc_module
  use multifab_physbc_module
  use multifab_fill_ghost_module
  !use ml_restriction_module
  use pf_utility_module
  use ml_restrict_fill_module

  implicit none

  private

  public :: init_phi_on_level, init_phi,trans_result_load
  public :: init_flw_on_level, init_flw,trans_result_load_flw


contains

  subroutine init_phi_on_level(phi,ad,dx,seed_pos,prob_lo,prob_hi,the_bc_level,pfpara)

    type(multifab) , intent(inout) :: phi,ad
    type(pf_para),   intent(in   ) :: pfpara
    real(kind=dp_t), intent(in   ) :: dx
    real(kind=dp_t), intent(inout) :: seed_pos(:,:)
    real(kind=dp_t), intent(in   ) :: prob_lo(phi%dim)
    real(kind=dp_t), intent(in   ) :: prob_hi(phi%dim)
    type(bc_level) , intent(in   ) :: the_bc_level

    ! local
    integer i,ng,dm
    integer :: lo(phi%dim), hi(phi%dim)

    real(kind=dp_t), pointer :: dp(:,:,:,:)
    real(kind=dp_t), pointer :: da(:,:,:,:)

    ng = phi%ng
    dm = phi%dim

    do i=1,nfabs(phi)
       dp => dataptr(phi,i)
       da => dataptr(ad,i)
       lo = lwb(get_box(phi,i))
       hi = upb(get_box(phi,i))

       select case(dm)
       case (2)
          call init_phi_2d(da(:,:,1,1),dp(:,:,1,1),dp(:,:,1,2),dp(:,:,1,3),ng,lo,hi,prob_lo,prob_hi,dx, &
                         pfpara%ori_def,pfpara%temp_h_l,pfpara%temp_h_r,pfpara%temp_l,pfpara%seed_radius,seed_pos,&
                         pfpara%seed_num, pfpara%seed_type)
       case (3)
          call init_phi_3d(da(:,:,:,1),dp(:,:,:,1),dp(:,:,:,2),dp(:,:,:,3),ng,lo,hi,prob_lo,prob_hi,dx, &
                         pfpara%ori_def,pfpara%temp_h_l,pfpara%temp_h_r,pfpara%temp_l,pfpara%seed_radius,seed_pos,&
                         pfpara%seed_num, pfpara%seed_type)
       end select
    end do

    call multifab_fill_boundary(phi)
    call multifab_fill_boundary(ad)

    call multifab_physbc(phi,1,1,3,the_bc_level)
    call multifab_physbc(ad,1,1,1,the_bc_level)

  end subroutine init_phi_on_level

  subroutine init_phi(mla,phi,ad,dx,seed_pos,prob_lo,prob_hi,the_bc_tower,pfpara)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: phi(:),ad(:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(inout) :: seed_pos(:,:)
    real(kind=dp_t), intent(in   ) :: prob_lo(mla%dim)
    real(kind=dp_t), intent(in   ) :: prob_hi(mla%dim)
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(pf_para),   intent(in   ) :: pfpara

    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: nlevs, dm, ng, i, n

    real(kind=dp_t), pointer :: dp(:,:,:,:)
    real(kind=dp_t), pointer :: da(:,:,:,:)

    ng = phi(1)%ng
    dm = mla%dim
    nlevs = mla%nlevel

    do n=1,nlevs
       do i=1,nfabs(phi(n))
          dp => dataptr(phi(n),i)
          da => dataptr(ad(n),i)
          lo = lwb(get_box(phi(n),i))
          hi = upb(get_box(phi(n),i))

          select case(dm)
          case (2)
             call init_phi_2d(da(:,:,1,1),dp(:,:,1,1),dp(:,:,1,2),dp(:,:,1,3),ng,lo,hi,prob_lo,prob_hi,dx(n), &
                              pfpara%ori_def,pfpara%temp_h_l,pfpara%temp_h_r,pfpara%temp_l,pfpara%seed_radius,seed_pos,&
                              pfpara%seed_num, pfpara%seed_type)
          case (3)
             call init_phi_3d(da(:,:,:,1),dp(:,:,:,1),dp(:,:,:,2),dp(:,:,:,3),ng,lo,hi,prob_lo,prob_hi,dx(n), &
                              pfpara%ori_def,pfpara%temp_h_l,pfpara%temp_h_r,pfpara%temp_l,pfpara%seed_radius,seed_pos,&
                              pfpara%seed_num, pfpara%seed_type)
          end select
       end do

    end do

    call ml_restrict_and_fill(nlevs, phi, mla%mba%rr, the_bc_tower%bc_tower_array,1,1,3)
    call ml_restrict_and_fill(nlevs, ad, mla%mba%rr, the_bc_tower%bc_tower_array,1,1,1)

  end subroutine init_phi

  subroutine init_phi_2d (phi_ad,phi_pf,phi_uc,phi_th,ng,lo,hi,prob_lo,prob_hi,dx, &
                          ori,t_h_l,t_h_r,t_l,s_r,seed_pos,seed_num,seed_type)

  integer          :: lo(2), hi(2), ng,seed_num,seed_type
  double precision :: phi_ad(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
  double precision :: phi_pf(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
  double precision :: phi_uc(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
  double precision :: phi_th(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
  double precision :: prob_lo(2),prob_hi(2)
  double precision :: dx,ori,t_h_l,t_h_r,t_l,s_r
  real(kind=dp_t)  :: seed_pos(seed_num,4)

  ! local variables
  integer          :: i,j,s_i
  double precision :: x,x_l,y,r2,y_l,t_h,d_c= sqrt(2.d0)
  !double precision :: dis_sep(seed_num, 2), dis_unit

  ! The domain top to bottom distance
  y_l = prob_hi(2) - prob_lo(2)
  x_l = prob_hi(1) - prob_lo(1)

  select case (seed_type)
  case(1)
    ! columnar dendrite growth, started from the bottom of y direction
    !$omp parallel do private(i,j,y)
    do j = lo(2), hi(2)
       y = prob_lo(2) + (dble(j)+0.5d0) * dx
       do i = lo(1), hi(1)
        x = prob_lo(1) + (dble(i)+0.5d0) * dx

          phi_pf(i,j) = -1.d0
          phi_uc(i,j) = -1.d0  ! initialized as zero, model A, no need further change
          phi_ad(i,j) = -1.d0  ! initialized as zero, i.e. along axis

          t_h = t_h_l + (t_h_r - t_h_l) * x / x_l

          phi_th(i,j) =  t_l + (t_h - t_l) * y / y_l  ! setup a temperature gradient

          if (y < s_r) then
             phi_pf(i,j)     = 1.d0
             phi_ad(i,j)     = atan(ori)
          end if
       end do
    end do
    !$omp end parallel do

  case(2)
    ! columnar dendrite growth, started from the bottom of y direction
    !$omp parallel do private(i,j,y)
    do j = lo(2), hi(2)
       y = prob_lo(2) + (dble(j)+0.5d0) * dx
       do i = lo(1), hi(1)
        x = prob_lo(1) + (dble(i)+0.5d0) * dx

          phi_pf(i,j) = -1.d0
          phi_uc(i,j) = -1.d0  ! initialized as zero, model A, no need further change
          phi_ad(i,j) = -1.d0  ! initialized as zero, i.e. along axis

          t_h = t_h_l + (t_h_r - t_h_l) * x / x_l

          phi_th(i,j) =  t_l + (t_h - t_l) * y / y_l  ! setup a temperature gradient

          if (x <= s_r .or. y <= s_r .or. x >= (x_l-s_r) .or. y>= (y_l-s_r) ) then
             phi_pf(i,j)     = 1.d0
             phi_ad(i,j)     = atan(ori)
          end if
       end do
    end do
    !$omp end parallel do

  case default
    ! equiaxed dendrite growth
    !$omp parallel do private(i,j,y)
    do j = lo(2), hi(2)
     y = prob_lo(2) + (dble(j)+0.5d0) * dx
     do i = lo(1), hi(1)
        x = prob_lo(1) + (dble(i)+0.5d0) * dx

        phi_pf(i,j) = -1.d0
        phi_uc(i,j) = -1.d0  ! initialized as zero, model A, no need further change
        phi_ad(i,j) = -1.d0  ! initialized as zero, i.e. along axis

        t_h = t_h_l + (t_h_r - t_h_l) * x / x_l

        phi_th(i,j) =  t_l + (t_h - t_l) * y / y_l  ! setup a temperature gradient

        do s_i = 1,seed_num

           r2 = sqrt( (x - seed_pos(s_i,1))**2 + (y - seed_pos(s_i,2))**2)

           !phi_pf(i,j)     = -tanh(r2-20.d0)

           if (r2 < s_r) then
              phi_pf(i,j)     = 1.d0
              phi_ad(i,j)     = dble(s_i)
           end if

        end do

     end do
    end do
    !$omp end parallel do

  end select

end subroutine init_phi_2d

subroutine init_phi_3d (phi_ad,phi_pf,phi_uc,phi_th,ng,lo,hi,prob_lo,prob_hi,dx,&
                        ori,t_h_l,t_h_r,t_l,s_r,seed_pos,seed_num,seed_type)

  integer          :: lo(3), hi(3), ng,seed_num,seed_type
  double precision :: phi_ad(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
  double precision :: phi_pf(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
  double precision :: phi_uc(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
  double precision :: phi_th(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
  double precision :: prob_lo(3),prob_hi(3)
  double precision :: dx, ori,t_h_l,t_h_r,t_l,s_r
  real(kind=dp_t)  :: seed_pos(seed_num,4)

  ! local variables
  integer          :: i,j,k,s_i
  double precision :: x,x_l,y,z,r2,y_l,t_h,d_c=sqrt(2.d0)
  double precision :: x_left, x_right, y_top, y_bottom

  ! The domain top to bottom distance
  y_l = prob_hi(2) - prob_lo(2)
  x_l = prob_hi(1) - prob_lo(1)

  select case (seed_type)
  case(1)
    !$omp parallel do private(i,j,k,y)
     do k = lo(3), hi(3)
       do j = lo(2), hi(2)
         y = prob_lo(2) + (dble(j)+0.5d0) * dx
         do i = lo(1), hi(1)
            x = prob_lo(1) + (dble(i)+0.5d0) * dx

            phi_pf(i,j,k) = -1.d0
            phi_uc(i,j,k) = -1.d0  ! initialized as zero, model A, no need further change
            phi_ad(i,j,k) = -1.d0  ! initialized as zero, i.e. along axis

            t_h = t_h_l + (t_h_r - t_h_l) * x / x_l

            phi_th(i,j,k) =  t_l + (t_h - t_l) * y / y_l  ! setup a temperature gradient


            if (y < s_r) then
                phi_pf(i,j,k)     = 1.d0
                phi_ad(i,j,k)     = atan( ori )
            end if

          end do  !i
        end do  !j
     end do  !k
     !$omp end parallel do

  case (2)
    !$omp parallel do private(i,j,k,z,y,x,s_i,r2)
     do k = lo(3), hi(3)
       !z = prob_lo(3) + (dble(k)+0.5d0) * dx
       do j = lo(2), hi(2)
         y = prob_lo(2) + (dble(j)+0.5d0) * dx
         do i = lo(1), hi(1)
            x = prob_lo(1) + (dble(i)+0.5d0) * dx

            phi_pf(i,j,k) = -1.d0
            phi_uc(i,j,k) = -1.d0  ! initialized as zero, model A, no need further change
            phi_ad(i,j,k) = -1.d0  ! initialized as zero, i.e. along axis

            t_h = t_h_l + (t_h_r - t_h_l) * x / x_l

            phi_th(i,j,k) =  t_l + (t_h - t_l) * y / y_l  ! setup a temperature gradient

            do s_i = 1,seed_num

               r2 = sqrt( (x - seed_pos(s_i,1))**2 + (y - seed_pos(s_i,2))**2  )

               !phi_pf(i,j,k)     = -tanh(r2-20.d0)

               if (r2 < s_r) then
                  phi_pf(i,j,k)     = 1.d0
                  phi_ad(i,j,k)     = dble(s_i)
               end if
            end do

          end do  !i
        end do  !j
     end do  !k
     !$omp end parallel do

  case (3)
    !$omp parallel do private(i,j,k,z,y,x,s_i,r2)
     do k = lo(3), hi(3)
       !z = prob_lo(3) + (dble(k)+0.5d0) * dx
       do j = lo(2), hi(2)
         y = prob_lo(2) + (dble(j)+0.5d0) * dx
         do i = lo(1), hi(1)
            x = prob_lo(1) + (dble(i)+0.5d0) * dx

            phi_pf(i,j,k) = -1.d0
            phi_uc(i,j,k) = -1.d0  ! initialized as zero, model A, no need further change
            phi_ad(i,j,k) = -1.d0  ! initialized as zero, i.e. along axis

            t_h = t_h_l + (t_h_r - t_h_l) * x / x_l

            phi_th(i,j,k) =  t_l + (t_h - t_l) * y / y_l  ! setup a temperature gradient

            do s_i = 1,seed_num

               x_left   = seed_pos(s_i,1) - s_r; x_right = seed_pos(s_i,1) + s_r
               y_bottom = seed_pos(s_i,2) - s_r; y_top   = seed_pos(s_i,2) + s_r

               if ( (x >= x_left   .and. x <= x_right) .and. &
                    (y >= y_bottom .and. y <= y_top  ) ) then
                  phi_pf(i,j,k)     = 1.d0
                  phi_ad(i,j,k)     = dble(s_i)
               end if
            end do

          end do  !i
        end do  !j
     end do  !k
     !$omp end parallel do

  case default
    !$omp parallel do private(i,j,k,z,y,x,s_i,r2)
     do k = lo(3), hi(3)
       z = prob_lo(3) + (dble(k)+0.5d0) * dx
       do j = lo(2), hi(2)
         y = prob_lo(2) + (dble(j)+0.5d0) * dx
         do i = lo(1), hi(1)
            x = prob_lo(1) + (dble(i)+0.5d0) * dx

            phi_pf(i,j,k) = -1.d0
            phi_uc(i,j,k) = -1.d0  ! initialized as zero, model A, no need further change
            phi_ad(i,j,k) = -1.d0  ! initialized as zero, i.e. along axis

            t_h = t_h_l + (t_h_r - t_h_l) * x / x_l

            phi_th(i,j,k) =  t_l + (t_h - t_l) * y / y_l  ! setup a temperature gradient

            do s_i = 1,seed_num

               r2 = sqrt( (x - seed_pos(s_i,1))**2 + (y - seed_pos(s_i,2))**2 + (z - seed_pos(s_i,3))**2 )

               !phi_pf(i,j,k)     = -tanh(r2-20.d0)

               if (r2 < s_r) then
                  phi_pf(i,j,k)     = 1.d0
                  phi_ad(i,j,k)     = dble(s_i)
               end if
            end do

          end do  !i
        end do  !j
     end do  !k
     !$omp end parallel do

  end select

end subroutine init_phi_3d

subroutine trans_result_load(mla,phi_old,ad_ori,phi_load,dx,the_bc_tower,time,&
                               prob_lo,prob_hi,pfpara)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: phi_old(:),ad_ori(:)
    type(multifab) , pointer       :: phi_load(:)
    real(kind=dp_t), intent(in   ) :: dx(:),time
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(pf_para),   intent(in   ) :: pfpara
    real(kind=dp_t), intent(in   ) :: prob_lo(mla%dim)
    real(kind=dp_t), intent(in   ) :: prob_hi(mla%dim)

    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: nlevs, dm, ng, i, n

    real(kind=dp_t), pointer :: dp(:,:,:,:)
    real(kind=dp_t), pointer :: da(:,:,:,:)
    real(kind=dp_t), pointer :: dl(:,:,:,:)

    ng = phi_old(1)%ng
    dm = mla%dim
    nlevs = mla%nlevel

    if(pfpara%plot_mode .eq. 3) then

    do n=1,nlevs
       do i=1,nfabs(phi_old(n))
          dp => dataptr(phi_old(n),i)
          da => dataptr(ad_ori(n),i)
          dl => dataptr(phi_load(n),i)
          lo = lwb(get_box(phi_old(n),i))
          hi = upb(get_box(phi_old(n),i))

          select case(dm)
          case (2)
             call load_phi_2d(da(:,:,1,1),dp(:,:,1,1),dp(:,:,1,2),dp(:,:,1,3),&
                              dl(:,:,1,1),dl(:,:,1,2),dl(:,:,1,3), time, &
                              ng,lo,hi,prob_lo,prob_hi,dx(n),pfpara%temp_h_l,pfpara%temp_h_r,pfpara%temp_l,&
                              pfpara%sk,pfpara%Rate_Cooling,pfpara%temp_lowcut)

          case (3)
             call load_phi_3d(da(:,:,:,1),dp(:,:,:,1),dp(:,:,:,2),dp(:,:,:,3),&
                              dl(:,:,:,1),dl(:,:,:,2),dl(:,:,:,3), time, &
                              ng,lo,hi,prob_lo,prob_hi,dx(n),pfpara%temp_h_l,pfpara%temp_h_r,pfpara%temp_l,&
                              pfpara%sk,pfpara%Rate_Cooling,pfpara%temp_lowcut)
          end select
       end do

    end do

    else if(pfpara%plot_mode .eq. 0) then

    do n=1,nlevs
       do i=1,nfabs(phi_old(n))
          dp => dataptr(phi_old(n),i)
          da => dataptr(ad_ori(n),i)
          dl => dataptr(phi_load(n),i)
          lo = lwb(get_box(phi_old(n),i))
          hi = upb(get_box(phi_old(n),i))

          select case(dm)
          case (2)
             call load_phi_2d_0(da(:,:,1,1),dp(:,:,1,1),dp(:,:,1,2),dp(:,:,1,3),&
                              dl(:,:,1,1),dl(:,:,1,2),dl(:,:,1,3),dl(:,:,1,4), time, &
                              ng,lo,hi,prob_lo,prob_hi,dx(n),pfpara%temp_h_l,pfpara%temp_h_r,pfpara%temp_l,&
                              pfpara%sk,pfpara%Rate_Cooling,pfpara%temp_lowcut,pfpara%coupled_mode,&
                              pfpara%cal_tem_mode)

          case (3)
             call load_phi_3d_0(da(:,:,:,1),dp(:,:,:,1),dp(:,:,:,2),dp(:,:,:,3),&
                              dl(:,:,:,1),dl(:,:,:,2),dl(:,:,:,3),dl(:,:,1,4), time, &
                              ng,lo,hi,prob_lo,prob_hi,dx(n),pfpara%temp_h_l,pfpara%temp_h_r,pfpara%temp_l,&
                              pfpara%sk,pfpara%Rate_Cooling,pfpara%temp_lowcut,pfpara%coupled_mode,&
                              pfpara%cal_tem_mode)
          end select
       end do

    end do

    end if

    call ml_restrict_and_fill(nlevs, phi_old, mla%mba%rr, the_bc_tower%bc_tower_array,1,1,3)
    call ml_restrict_and_fill(nlevs, ad_ori, mla%mba%rr, the_bc_tower%bc_tower_array,1,1,1)

end subroutine trans_result_load

subroutine load_phi_2d (phi_ad,phi_pf,phi_uc,phi_th,lod_ad,lod_pf,lod_trc,time,&
                        ng,lo,hi,prob_lo,prob_hi,dx,t_h_l,t_h_r,t_l,p_k,r_cool,t_cut)

  integer          :: lo(2), hi(2), ng
  double precision :: phi_ad(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
  double precision :: phi_pf(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
  double precision :: phi_uc(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
  double precision :: phi_th(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
  double precision :: lod_ad(lo(1):hi(1),lo(2):hi(2))
  double precision :: lod_pf(lo(1):hi(1),lo(2):hi(2))
  double precision :: lod_trc(lo(1):hi(1),lo(2):hi(2))
  double precision :: prob_lo(2),prob_hi(2)
  double precision :: dx,t_h_l,t_h_r,t_l,p_k,time,r_cool,t_cut

  ! local variables
  integer          :: i,j
  double precision :: y,y_l,x,x_l,t_h

  ! The domain top to bottom distance
  y_l = prob_hi(2) - prob_lo(2)
  x_l = prob_hi(1) - prob_lo(1)

  ! columnar dendrite growth, started from the bottom of y direction
  !$omp parallel do private(i,j,y)
  do j = lo(2), hi(2)
    y = prob_lo(2) + (dble(j)+0.5d0) * dx
    do i = lo(1), hi(1)
        x = prob_lo(1) + (dble(i)+0.5d0) * dx

       phi_ad(i,j) = lod_ad(i,j)
       phi_pf(i,j) = lod_pf(i,j)
       phi_uc(i,j) = 2.d0*p_k*lod_trc(i,j)/(1.d0+p_k-(1.d0-p_k)*lod_pf(i,j)) - 1.d0
       phi_uc(i,j) = phi_uc(i,j) / (1.d0-p_k)

       t_h = t_h_l + (t_h_r - t_h_l) * x / x_l

       phi_th(i,j) =  t_l + (t_h - t_l) * y / y_l

       if(phi_th(i,j) < t_cut) then
           phi_th(i,j) = t_cut
       end if

    end do
  end do
  !$omp end parallel do

end subroutine load_phi_2d

subroutine load_phi_3d (phi_ad,phi_pf,phi_uc,phi_th,lod_ad,lod_pf,lod_trc,time,&
                        ng,lo,hi,prob_lo,prob_hi,dx,t_h_l,t_h_r,t_l,p_k,r_cool,t_cut)

  integer          :: lo(3), hi(3), ng
  double precision :: phi_ad(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
  double precision :: phi_pf(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
  double precision :: phi_uc(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
  double precision :: phi_th(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
  double precision :: lod_ad(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  double precision :: lod_pf(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  double precision :: lod_trc(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  double precision :: prob_lo(3),prob_hi(3)
  double precision :: dx, t_h_l,t_h_r,t_l,p_k,time,r_cool,t_cut

  ! local variables
  integer          :: i,j,k
  double precision :: y,y_l,x,x_l,t_h

  ! The domain top to bottom distance
  y_l = prob_hi(2) - prob_lo(2)
  x_l = prob_hi(1) - prob_lo(1)

  !$omp parallel do private(i,j,k,y)
  do k = lo(3), hi(3)
    do j = lo(2), hi(2)
       y = prob_lo(2) + (dble(j)+0.5d0) * dx
       do i = lo(1), hi(1)
        x = prob_lo(1) + (dble(i)+0.5d0) * dx

         phi_ad(i,j,k) = lod_ad(i,j,k)
         phi_pf(i,j,k) = lod_pf(i,j,k)
         phi_uc(i,j,k) = 2.d0*p_k*lod_trc(i,j,k)/(1.d0+p_k-(1.d0-p_k)*lod_pf(i,j,k)) - 1.d0
         phi_uc(i,j,k) = phi_uc(i,j,k) / (1.d0-p_k)

         t_h = t_h_l + (t_h_r - t_h_l) * x / x_l

         phi_th(i,j,k) =  t_l + (t_h - t_l) * y / y_l

         if(phi_th(i,j,k) < t_cut) then
            phi_th(i,j,k) = t_cut
         end if

       end do  !i
    end do  !j
  end do  !k
  !$omp end parallel do

end subroutine load_phi_3d

subroutine load_phi_2d_0 (phi_ad,phi_pf,phi_uc,phi_th,lod_ad,lod_pf,lod_uc,lod_th,time,&
                          ng,lo,hi,prob_lo,prob_hi,dx,t_h_l,t_h_r,t_l,p_k,r_cool,t_cut,&
                          coupled_mode,cal_tem_mode)

  integer          :: lo(2), hi(2), ng, coupled_mode,cal_tem_mode
  double precision :: phi_ad(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
  double precision :: phi_pf(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
  double precision :: phi_uc(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
  double precision :: phi_th(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
  double precision :: lod_ad(lo(1):hi(1),lo(2):hi(2))
  double precision :: lod_pf(lo(1):hi(1),lo(2):hi(2))
  double precision :: lod_uc(lo(1):hi(1),lo(2):hi(2))
  double precision :: lod_th(lo(1):hi(1),lo(2):hi(2))
  double precision :: prob_lo(2),prob_hi(2)
  double precision :: dx,t_h_l,t_h_r,t_l,p_k,time,r_cool,t_cut

  ! local variables
  integer          :: i,j
  double precision :: y,y_l,x,x_l,t_h

  ! The domain top to bottom distance
  y_l = prob_hi(2) - prob_lo(2)
  x_l = prob_hi(1) - prob_lo(1)

  ! columnar dendrite growth, started from the bottom of y direction
  do j = lo(2), hi(2)
    y = prob_lo(2) + (dble(j)+0.5d0) * dx
    do i = lo(1), hi(1)
      x = prob_lo(1) + (dble(i)+0.5d0) * dx

       phi_ad(i,j) = lod_ad(i,j)
       phi_pf(i,j) = lod_pf(i,j)
       phi_uc(i,j) = lod_uc(i,j)

       if(coupled_mode==1 .and. cal_tem_mode==1) then
          phi_th(i,j) = lod_th(i,j)

       else

          t_h = t_h_l + (t_h_r - t_h_l) * x / x_l

          phi_th(i,j) =  t_l + (t_h - t_l) * y / y_l

          if(phi_th(i,j) < t_cut) then
             phi_th(i,j) = t_cut
          end if

       end if

    end do
  end do

end subroutine load_phi_2d_0

subroutine load_phi_3d_0 (phi_ad,phi_pf,phi_uc,phi_th,lod_ad,lod_pf,lod_uc,lod_th,time,&
                          ng,lo,hi,prob_lo,prob_hi,dx,t_h_l,t_h_r,t_l,p_k,r_cool,t_cut, &
                          coupled_mode,cal_tem_mode)

  integer          :: lo(3), hi(3), ng, coupled_mode,cal_tem_mode
  double precision :: phi_ad(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
  double precision :: phi_pf(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
  double precision :: phi_uc(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
  double precision :: phi_th(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
  double precision :: lod_ad(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  double precision :: lod_pf(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  double precision :: lod_uc(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  double precision :: lod_th(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  double precision :: prob_lo(3),prob_hi(3)
  double precision :: dx, t_h_l,t_h_r,t_l,p_k,time,r_cool,t_cut

  ! local variables
  integer          :: i,j,k
  double precision :: y,y_l,x,x_l,t_h

  ! The domain top to bottom distance
  y_l = prob_hi(2) - prob_lo(2)
  x_l = prob_hi(1) - prob_lo(1)


  do k = lo(3), hi(3)
    do j = lo(2), hi(2)
       y = prob_lo(2) + (dble(j)+0.5d0) * dx
       do i = lo(1), hi(1)
        x = prob_lo(1) + (dble(i)+0.5d0) * dx

         phi_ad(i,j,k) = lod_ad(i,j,k)
         phi_pf(i,j,k) = lod_pf(i,j,k)
         phi_uc(i,j,k) = lod_uc(i,j,k)

         if(coupled_mode==1 .and. cal_tem_mode==1) then
            phi_th(i,j,k) = lod_th(i,j,k)
         else

            t_h = t_h_l + (t_h_r - t_h_l) * x / x_l

            phi_th(i,j,k) =  t_l + (t_h - t_l) * y / y_l

            if(phi_th(i,j,k) < t_cut) then
               phi_th(i,j,k) = t_cut
            end if

         end if

       end do  !i
    end do  !j
  end do  !k

end subroutine load_phi_3d_0

subroutine init_flw_on_level(flw,phi,the_bc_level,pfpara)

    type(multifab) , intent(inout) :: flw
    type(multifab) , intent(inout) :: phi
    type(pf_para),   intent(in   ) :: pfpara
    type(bc_level) , intent(in   ) :: the_bc_level

    ! local
    integer i,ng,dm
    integer :: lo(phi%dim), hi(phi%dim)
    integer :: iter_st(3),iter_end(3),icom_start,ncomponent

    real(kind=dp_t), pointer :: pfo(:,:,:,:)
    real(kind=dp_t), pointer :: u_vec(:,:,:,:)
    real(kind=dp_t), pointer :: f_mat(:,:,:,:)
    real(kind=dp_t), pointer :: k_mat(:,:,:,:)

    ng = phi%ng
    dm = phi%dim

    select case(dm)
    case (2)
      icom_start = 1
      ncomponent = 21
      iter_st  = (/1,4,13/)
      iter_end = (/3,9,9/)
    case (3)
      icom_start = 1
      ncomponent = 42
      iter_st  = (/1,5,24/)
      iter_end = (/4,19,19/)
    end select

    do i=1,nfabs(phi)
       pfo => dataptr(phi,i)
       u_vec => dataptr(flw,i,iter_st(1),iter_end(1))
       f_mat => dataptr(flw,i,iter_st(2),iter_end(2))
       k_mat => dataptr(flw,i,iter_st(3),iter_end(3))
       lo = lwb(get_box(phi,i))
       hi = upb(get_box(phi,i))

       select case(dm)
       case (2)
          call init_flw_2d(pfo(:,:,1,1),u_vec(:,:,1,:),f_mat(:,:,1,:),k_mat(:,:,1,:),ng,lo,hi,&
                           pfpara%flw_rho)
       case (3)
          call init_flw_3d(pfo(:,:,:,1),u_vec(:,:,:,:),f_mat(:,:,:,:),k_mat(:,:,:,:),ng,lo,hi,&
                           pfpara%flw_rho)
       end select
    end do

    call multifab_fill_boundary_c(flw,icom_start,ncomponent,ng)
    call multifab_physbc(flw,icom_start,1,ncomponent,the_bc_level)

    if(pfpara%kill_phi_flag .eq. 1) call multifab_setval_c(phi,-1.d0,1,1)
    
  end subroutine init_flw_on_level

  subroutine init_flw(mla,flw,phi,the_bc_tower,pfpara)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: flw(:)
    type(multifab) , intent(inout) :: phi(:)
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(pf_para),   intent(in   ) :: pfpara

    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: nlevs, dm, ng, i, n
    integer :: iter_st(3),iter_end(3),icom_start,ncomponent

    real(kind=dp_t), pointer :: pfo(:,:,:,:)
    real(kind=dp_t), pointer :: u_vec(:,:,:,:)
    real(kind=dp_t), pointer :: f_mat(:,:,:,:)
    real(kind=dp_t), pointer :: k_mat(:,:,:,:)

    ng = phi(1)%ng
    dm = mla%dim
    nlevs = mla%nlevel

    select case(dm)
    case (2)
      icom_start = 1
      ncomponent = 21
      iter_st  = (/1,4,13/)
      iter_end = (/3,9,9/)
    case (3)
      icom_start = 1
      ncomponent = 42
      iter_st  = (/1,5,24/)
      iter_end = (/4,19,19/)
    end select

    do n=1,nlevs
       do i=1,nfabs(phi(n))
          pfo => dataptr(phi(n),i)
          u_vec => dataptr(flw(n),i,iter_st(1),iter_end(1))
          f_mat => dataptr(flw(n),i,iter_st(2),iter_end(2))
          k_mat => dataptr(flw(n),i,iter_st(3),iter_end(3))
          lo = lwb(get_box(phi(n),i))
          hi = upb(get_box(phi(n),i))

          select case(dm)
          case (2)
             call init_flw_2d(pfo(:,:,1,1),u_vec(:,:,1,:),f_mat(:,:,1,:),k_mat(:,:,1,:),ng,lo,hi,&
                              pfpara%flw_rho)
          case (3)
             call init_flw_3d(pfo(:,:,:,1),u_vec(:,:,:,:),f_mat(:,:,:,:),k_mat(:,:,:,:),ng,lo,hi,&
                              pfpara%flw_rho)
          end select

       end do

       if(pfpara%kill_phi_flag .eq. 1) call multifab_setval_c(phi(n),-1.d0,1,1)
    end do

    ! restrict the multi-level data, and
    ! fill all boundaries: same-level, coarse-fine, periodic, and domain boundaries
    call ml_restrict_and_fill(nlevs, flw, mla%mba%rr, the_bc_tower%bc_tower_array,icom_start,1,ncomponent)

  end subroutine init_flw

  subroutine init_flw_2d (phi,u_vec,f_mat,k_mat,ng,lo,hi,rho)

  integer          :: lo(2), hi(2), ng,kill_phi_flag
  double precision :: phi(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
  double precision :: u_vec(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,1:3)
  double precision :: f_mat(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,1:9)
  double precision :: k_mat(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,1:9)
  double precision :: rho

  ! local variables
  integer          :: i,j,iter

  ! local varables
  integer,parameter:: q=9
  double precision, parameter :: p0_2d = 4.d0/9.d0, p1_2d = 1.d0/9.d0, p2_2d=1.d0/36.d0
  integer          :: ex(q)=(/ 0,1,0,-1,0,1,-1,-1,1 /),ey(q)=(/ 0,0,1,0,-1,1,1,-1,-1 /)
  double precision :: w(q)=(/ p0_2d,p1_2d,p1_2d,p1_2d,p1_2d,p2_2d,p2_2d,p2_2d,p2_2d /) 

  ! equiaxed dendrite growth
    do j = lo(2), hi(2)
     do i = lo(1), hi(1)

        u_vec(i,j,1) =0.d0
        u_vec(i,j,2) =0.d0
        u_vec(i,j,3) =rho

        do iter=1,q
            f_mat(i,j,iter)=w(iter)*u_vec(i,j,3)

        end do

        k_mat(i,j,:)=0.d0

     end do
    end do

end subroutine init_flw_2d

subroutine init_flw_3d (phi,u_vec,f_mat,k_mat,ng,lo,hi,rho)

  integer          :: lo(3), hi(3), ng
  double precision :: phi(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
  double precision :: u_vec(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,1:4)
  double precision :: f_mat(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,1:19)
  double precision :: k_mat(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,1:19)
  double precision :: rho

  ! local variables
  integer          :: i,j,k,iter
  integer,parameter:: q=19
  double precision, parameter :: p0_3d = 1.d0/3.d0, p1_3d = 1.d0/18.d0, p2_3d=1.d0/36.d0
  integer          :: ex(q)=(/ 0,1,-1,0,0,0,0,1,1,-1,-1,1,-1,1,-1,0,0,0,0 /)
  integer          :: ey(q)=(/ 0,0,0,1,-1,0,0,1,-1,1,-1,0,0,0,0,1,1,-1,-1 /)
  integer          :: ez(q)=(/ 0,0,0,0,0,1,-1,0,0,0,0,1,1,-1,-1,1,-1,1,-1 /)
  double precision :: w(q)=(/ p0_3d,p1_3d,p1_3d,p1_3d,p1_3d,p1_3d,p1_3d,p2_3d,p2_3d,p2_3d,&
                           p2_3d,p2_3d,p2_3d,p2_3d,p2_3d,p2_3d,p2_3d,p2_3d,p2_3d /)  

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
         do i = lo(1), hi(1)

            u_vec(i,j,k,1:3) =0.d0
            u_vec(i,j,k,4) =rho

            do iter=1,q
               f_mat(i,j,k,iter)=w(iter)*u_vec(i,j,k,4)

            end do

            k_mat(i,j,k,:) = 0.d0

          end do  !i
        end do  !j
     end do  !k

end subroutine init_flw_3d

subroutine trans_result_load_flw(mla,flw,phi_old,ad_ori,phi_load,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: flw(:),phi_old(:),ad_ori(:)
    type(multifab) , pointer       :: phi_load(:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: nlevs, dm, ng, i, n,nc_flw_0

    real(kind=dp_t), pointer :: dp(:,:,:,:)
    real(kind=dp_t), pointer :: da(:,:,:,:)
    real(kind=dp_t), pointer :: dl(:,:,:,:)
    real(kind=dp_t), pointer :: df(:,:,:,:)
    real(kind=dp_t), pointer :: dll(:,:,:,:)

    ng = phi_old(1)%ng
    dm = mla%dim
    nlevs = mla%nlevel

    if(dm .eq. 2) then
       nc_flw_0 = 12
    else if(dm .eq. 3) then
       nc_flw_0 = 23
    end if

    do n=1,nlevs
       do i=1,nfabs(phi_old(n))
          dp => dataptr(phi_old(n),i)
          da => dataptr(ad_ori(n),i)
          dl => dataptr(phi_load(n),i,1,4)
          dll => dataptr(phi_load(n),i,5,nc_flw_0)
          df => dataptr(flw(n),i)
          lo = lwb(get_box(phi_old(n),i))
          hi = upb(get_box(phi_old(n),i))

          select case(dm)
          case (2)
             call load_flw_2d(da(:,:,1,1),dp(:,:,1,1),dp(:,:,1,2),dp(:,:,1,3),df(:,:,1,:),&
                              dl(:,:,1,:),dll(:,:,1,:),ng,lo,hi)

          case (3)
             call load_flw_3d(da(:,:,:,1),dp(:,:,:,1),dp(:,:,:,2),dp(:,:,:,3),df(:,:,:,:),&
                              dl(:,:,:,:),dll(:,:,:,:),ng,lo,hi)
          end select
       end do

    end do

    !call ml_restrict_and_fill(nlevs, phi_old, mla%mba%rr, the_bc_tower%bc_tower_array,1,1,3)
    !call ml_restrict_and_fill(nlevs, ad_ori, mla%mba%rr, the_bc_tower%bc_tower_array,1,1,1)
    !call ml_restrict_and_fill(nlevs, flw, mla%mba%rr, the_bc_tower%bc_tower_array,1,1,nc_flw_0)

    do n=1,nlevs
      call multifab_fill_boundary_c(ad_ori(n), 1, 1, ng)
      call multifab_fill_boundary_c(phi_old(n), 1, 3, ng)
      call multifab_fill_boundary_c(flw(n), 1, nc_flw_0, ng)

      call multifab_physbc(ad_ori(n),1,1,1,the_bc_tower%bc_tower_array(n))
      call multifab_physbc(phi_old(n),1,1,3,the_bc_tower%bc_tower_array(n))
      call multifab_physbc(flw(n),1,1,nc_flw_0,the_bc_tower%bc_tower_array(n))
    end do

end subroutine trans_result_load_flw

subroutine load_flw_2d (phi_ad,phi_pf,phi_uc,phi_th,flw,phi_load,flw_load,ng,lo,hi)

  integer          :: lo(2), hi(2), ng
  double precision :: phi_ad(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
  double precision :: phi_pf(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
  double precision :: phi_uc(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
  double precision :: phi_th(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
  double precision :: flw(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,1:21)
  double precision :: phi_load(lo(1):hi(1),lo(2):hi(2),1:4)
  double precision :: flw_load(lo(1):hi(1),lo(2):hi(2),1:12)

  ! local variables
  integer          :: i,j,q

  ! columnar dendrite growth, started from the bottom of y direction
  !$omp parallel do private(i,j,y)
  do j = lo(2), hi(2)
    do i = lo(1), hi(1)

       phi_ad(i,j) = phi_load(i,j,1)
       phi_pf(i,j) = phi_load(i,j,2)
       phi_uc(i,j) = phi_load(i,j,3)
       phi_th(i,j) = phi_load(i,j,4)

       flw(i,j,1:12)  = flw_load(i,j,1:12)
       flw(i,j,13:21) = 0.d0
       
    end do
  end do
  !$omp end parallel do

end subroutine load_flw_2d

subroutine load_flw_3d (phi_ad,phi_pf,phi_uc,phi_th,flw,phi_load,flw_load,ng,lo,hi)

  integer          :: lo(3), hi(3), ng
  double precision :: phi_ad(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
  double precision :: phi_pf(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
  double precision :: phi_uc(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
  double precision :: phi_th(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
  double precision :: flw(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,1:42)
  double precision :: phi_load(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:4)
  double precision :: flw_load(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:23)


  ! local variables
  integer          :: i,j,k,q

  !$omp parallel do private(i,j,k,y)
  do k = lo(3), hi(3)
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

         phi_ad(i,j,k) = phi_load(i,j,k,1)
         phi_pf(i,j,k) = phi_load(i,j,k,2)
         phi_uc(i,j,k) = phi_load(i,j,k,3)
         phi_th(i,j,k) = phi_load(i,j,k,4)
         flw(i,j,k,1:23)  = flw_load(i,j,k,1:23)
         flw(i,j,k,24:42) = 0.d0

       end do  !i
    end do  !j
  end do  !k
  !$omp end parallel do

end subroutine load_flw_3d

end module init_phi_module
