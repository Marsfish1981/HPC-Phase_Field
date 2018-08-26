module pf_utility_module

  use multifab_module
  use layout_module
  use ml_layout_module
  use define_bc_module
  use bc_module
  !use ml_restriction_module
  use multifab_physbc_module
  use multifab_fill_ghost_module
  use ml_restrict_fill_module

  implicit none

  type pf_para

     integer          ::   dim       ! the dimension of the problem
     integer          ::   no_noise  ! switch noise if set this value to 1
     integer          ::   plot_mode  !0 = AD,PF,TH,UC,TRC, 1 = AD,PF,TH,UC,TRC,GRA, 2 = PF,TRC
     integer          ::   coupled_mode ! 0 = temperature is specified, 1 = temperature is calculated,
     integer          ::   cal_tem_mode ! 0 = explicit mode, 1 = implicit mode
                                        ! The calculation is defined using pairs of [coupled_mode,cal_tem_mode]
                                        ! [0,0] : T is copy from previous time step
                                        ! [0,1] : T is the sum of T0 and Rate_Cooling * time
                                        ! [1,0] : T is calculated using explicit algorithm
                                        ! [1,1] : T is calculated using MG algorithm
     integer          ::   cooling_mode

     Real(kind=dp_t)  ::   lamda    ! the scaling parameter
     Real(kind=dp_t)  ::   anis     ! the strength of anisotropy
     Real(kind=dp_t)  ::   anis_0     ! the strength of anisotropy
     Real(kind=dp_t)  ::   Le       ! the Lewis number
     Real(kind=dp_t)  ::   kk       ! the hamonic symmetric factor
     Real(kind=dp_t)  ::   DM       ! the solute diffusivity
     Real(kind=dp_t)  ::   sk      ! the solute partition parameter
     Real(kind=dp_t)  ::   Mc00     ! the solute scaling parameter
     Real(kind=dp_t)  ::   a1     ! constant
     Real(kind=dp_t)  ::   a2     ! constant
     Real(kind=dp_t)  ::   thers    ! a small parameter
     Real(kind=dp_t)  ::   noiamp  ! the amplitude of noise
     Real(kind=dp_t)  ::   temp_h_l  ! the parameter used for tag related to UC
     Real(kind=dp_t)  ::   temp_h_r  ! the parameter used for tag related to UC
     Real(kind=dp_t)  ::   temp_l  ! the parameter used for tag related to UC
     Real(kind=dp_t)  ::   t_DM    ! the diffusivity of temperature
     Real(kind=dp_t)  ::   r_hyper  ! the ratio between L/Cp and deltaT0
     Real(kind=dp_t)  ::   Rate_Cooling    ! the scaling parameter
     Real(kind=dp_t)  ::   temp_lowcut    ! the lowest temperature allowed

     integer          ::   period_type
     integer          ::   period_NT
     Real(kind=dp_t)  ::   period_p0

     Real(kind=dp_t)  ::   dx       ! the spatial step
     Real(kind=dp_t)  ::   para_uc  ! the parameter used for tag related to UC
     Real(kind=dp_t)  ::   para_th  ! the parameter used for tag related to TH
     Real(kind=dp_t)  ::   amr_thers
     Real(kind=dp_t)  ::   safe_n_init  ! the number of cells be tagged initially around the seed center

     Real(kind=dp_t)  ::   seed_radius  ! the parameter used for tag related to UC
     integer          ::   seed_num  ! the parameter used for tag related to UC
     integer          ::   seed_type  ! 0:equiaxed, 1:columnar
     Real(kind=dp_t)  ::   ori_def  ! the parameter used for tag related to UC

     ! For solving flow
     integer          ::   do_solve_flow
     integer          ::   solve_flow_only
     integer          ::   kill_phi_flag
     integer          ::   kill_body_force
     integer          ::   kill_lorentz_force
     Real(kind=dp_t)  ::   flw_amp   ! Amplified factor
     integer          ::   flw_skip_level_n
     Real(kind=dp_t)  ::   gacc  ! acceleration of gravity
     Real(kind=dp_t)  ::   belta_c ! solute-density coefficient = beltac*c00/(2k)
     Real(kind=dp_t)  ::   dt_flw
     Real(kind=dp_t)  ::   multi_add_tau
     Real(kind=dp_t)  ::   flw_constrain_limit

     Real(kind=dp_t)  ::   lorentz_omega
     Real(kind=dp_t)  ::   lorentz_amp_x   
     Real(kind=dp_t)  ::   lorentz_amp_y 
     Real(kind=dp_t)  ::   lorentz_amp_z 

     ! For solving flow
     integer          ::   flw_calculate_mode
     Real(kind=dp_t)  ::   u_inlet_def(3,2)
     Real(kind=dp_t)  ::   u_inlet(3,2)   ! the inlet velocity
     Real(kind=dp_t)  ::   tau_LBM   ! the relaxation time
     Real(kind=dp_t)  ::   lbm_viscosity  ! the viscosity of the LBM method
     Real(kind=dp_t)  ::   lbm_scale_dx    ! the length used to scale the length
     Real(kind=dp_t)  ::   lbm_scale_dt    ! the time used to scale the time
     Real(kind=dp_t)  ::   lbm_ratio_dx_W0 ! the ratio between dx and W0
     Real(kind=dp_t)  ::   lbm_drag_hc     ! the coefficient from Beckmann's
     Real(kind=dp_t)  ::   thers_solid  ! the solid fraction to cut off the flow
     Real(kind=dp_t)  ::   para_u  ! the solid fraction to cut off the flow
     Real(kind=dp_t)  ::   flw_rho  ! the solid fraction to cut off the flow
     Real(kind=dp_t)  ::   thers_tau
     integer          ::   cal_flw_stf_lev
	   integer          ::   cal_boundary_spec
     integer          ::   dim_lbm
     integer          ::   nsteps_for_upper
     Real(kind=dp_t)  ::   uinlet_set_init
     integer          ::   flag_set_uinlet
     Real(kind=dp_t)  ::   uinlet_set_step
     
     integer, pointer, dimension(:,:)             :: bc_flag => Null()
     real(kind = dp_t), pointer, dimension(:)     :: tau_r => Null()
     real(kind = dp_t), pointer, dimension(:,:,:) :: R_M => Null()
     real(kind = dp_t), pointer, dimension(:)     :: T_D => Null()

     integer, pointer, dimension(:,:)             :: lbm_e => Null()
     real(kind = dp_t), pointer, dimension(:)     :: lbm_w => Null()

  end type pf_para

  private

  public :: pf_para, gen_random, update_rkpf_2d, update_rkpf_3d, update_rkth_2d, update_rkth_3d
  public :: pf_para_destroy,pf_para_set_seed_info,follow_tip_info,write_tkfile,write_status_file
  public :: recalculate_regrid_number,retrieve_max_velocity,write_status_file2,pf_para_set_temp_segment

  public :: pf_para_set_flw_tau,pf_para_set_lbm
  public :: pf_set_uinlet,retrieve_error_msg,recalculate_time_step

contains

subroutine pf_para_set_temp_segment(pfpara, temp_vec, num_seg, num_tstep)

    type (pf_para), intent(inout) :: pfpara
    real(kind=dp_t),intent(in   ) :: temp_vec(:,:)
    integer,        intent(in   ) :: num_seg, num_tstep

    ! local variables
    integer :: i,j,nc = 2,id_start, id_end
    double precision :: slope

    ! allocate the memory
    allocate( pfpara%T_D(num_tstep) )

    do i= 1,num_tstep

        pfpara%T_D(i) = temp_vec(1,2)

    end do

    if(num_seg .ge. 2) then

      ! setup the temperature vector
      do j=1,num_seg-1

        ! linear intepolation
        slope = (temp_vec(j+1,2) - temp_vec(j,2))/( temp_vec(j+1,1) - temp_vec(j,1) )
        id_start = int( temp_vec(j,1)   )
        id_end   = int( temp_vec(j+1,1) )

        do i= id_start,id_end

          pfpara%T_D(i) = temp_vec(j,2) + slope*( ( dble(i) - temp_vec(j,1) ) )

        end do
      end do

    end if



end subroutine pf_para_set_temp_segment

subroutine pf_para_destroy(pfpara)

   type (pf_para), intent(inout) :: pfpara

   if ( associated(pfpara%R_M) ) then
       deallocate(pfpara%R_M)
   end if

   if ( associated(pfpara%T_D) ) then
       deallocate(pfpara%T_D)
   end if

   if ( associated(pfpara%bc_flag) ) then
       deallocate(pfpara%bc_flag)
   end if

   if ( associated(pfpara%lbm_e) ) then
       deallocate(pfpara%lbm_e)
   end if

   if ( associated(pfpara%lbm_w) ) then
       deallocate(pfpara%lbm_w)
   end if

   if ( associated(pfpara%tau_r) ) then
       deallocate(pfpara%tau_r)
   end if

end subroutine pf_para_destroy

subroutine pf_para_set_seed_info(pfpara,seed_pos,seed_num,r_type)

    type (pf_para), intent(inout) :: pfpara
    real(kind=dp_t),intent(in   ) :: seed_pos(:,:)
    integer,        intent(in   ) :: seed_num,r_type

    ! local variables
    integer i,j,dim
    dim = pfpara%dim

    if(seed_num .le. 0) return

    allocate(pfpara%R_M(seed_num,dim,dim))

    ! calculate the R_M
    select case(dim)
    case(2)

        do i=1,seed_num
           call cal_rot_matrix_2d( pfpara%R_M(i,:,:),atan( seed_pos(i,4) ) )
        end do

    case(3)

       if(r_type .eq. 0) then

          do i=1,seed_num
              call cal_rot_matrix_3d(pfpara%R_M(i,:,:), &
                                     atan( seed_pos(i,4) ), &
                                     seed_pos(i,5),seed_pos(i,6),seed_pos(i,7) )
          end do
       else if(r_type .eq. 1) then
          do i=1,seed_num
              call cal_rot_matrix_angle_3d(pfpara%R_M(i,:,:), &
                                           int(seed_pos(i,4)), &
                                           atan(seed_pos(i,5)),&
                                           atan(seed_pos(i,6)),&
                                           atan(seed_pos(i,7)) )
          end do
       end if
    end select

end subroutine pf_para_set_seed_info

subroutine cal_rot_matrix_2d(R_M,ori_z)

   double precision :: R_M(2,2), ori_z

   ! local variables
   double precision :: sin_z, cos_z

   ! calculate the sin, cos things
   sin_z = sin(ori_z)
   cos_z = cos(ori_z)

   ! calculate the R_M terms
   R_M(1,1) = cos_z
   R_M(1,2) = sin_z

   R_M(2,1) = -sin_z
   R_M(2,2) = cos_z

end subroutine cal_rot_matrix_2d

subroutine cal_rot_matrix_3d(R_M,ori,n_x,n_y,n_z)

   double precision :: R_M(3,3), ori,n_x, n_y, n_z

   ! local variables
   double precision :: sin_g, cos_g,nom_x,nom_y,nom_z,mod_c

   mod_c = sqrt(n_x**2 + n_y**2 + n_z**2)
   nom_x = n_x / mod_c
   nom_y = n_y / mod_c
   nom_z = n_z / mod_c

   sin_g = sin( ori )
   cos_g = cos( ori )

   ! calculate the R_M terms
   R_M(1,1) = cos_g + (1.d0 - cos_g) * nom_x**2
   R_M(1,2) = (1.d0-cos_g)*nom_x*nom_y-sin_g*nom_z
   R_M(1,3) = (1.d0-cos_g)*nom_x*nom_z+sin_g*nom_y

   R_M(2,1) = (1.d0-cos_g)*nom_x*nom_y+sin_g*nom_z
   R_M(2,2) = cos_g + (1.d0 - cos_g) * nom_y**2
   R_M(2,3) = (1.d0-cos_g)*nom_y*nom_z-sin_g*nom_x

   R_M(3,1) = (1.d0-cos_g)*nom_x*nom_z-sin_g*nom_y
   R_M(3,2) = (1.d0-cos_g)*nom_z*nom_y+sin_g*nom_x
   R_M(3,3) = cos_g + (1.d0 - cos_g) * nom_z**2

end subroutine cal_rot_matrix_3d

subroutine cal_rot_matrix_angle_3d(R_M,r_order,ang_x,ang_y,ang_z)

    double precision :: R_M(3,3), ang_x,ang_y,ang_z
    integer          :: r_order

    ! local variables
    double precision :: R_X(3,3), R_Y(3,3), R_Z(3,3),R_T(3,3)

    ! calculate the three rotation matrices
    call cal_rot_matrix_3d(R_X,ang_x,1.d0,0.d0,0.d0)
    call cal_rot_matrix_3d(R_Y,ang_y,0.d0,1.d0,0.d0)
    call cal_rot_matrix_3d(R_Z,ang_z,0.d0,0.d0,1.d0)

    ! calculate the final rotation matrix
    select case(r_order)
    case(1)
        call matrix_multi_matrix(R_T,R_Y,R_X,3)
        call matrix_multi_matrix(R_M,R_Z,R_T,3)
    case(2)
        call matrix_multi_matrix(R_T,R_X,R_Y,3)
        call matrix_multi_matrix(R_M,R_Z,R_T,3)
    case(3)
        call matrix_multi_matrix(R_T,R_Z,R_X,3)
        call matrix_multi_matrix(R_M,R_Y,R_T,3)
    case(4)
        call matrix_multi_matrix(R_T,R_X,R_Z,3)
        call matrix_multi_matrix(R_M,R_Y,R_T,3)
    case(5)
        call matrix_multi_matrix(R_T,R_Y,R_Z,3)
        call matrix_multi_matrix(R_M,R_X,R_T,3)
    case(6)
        call matrix_multi_matrix(R_T,R_Z,R_Y,3)
        call matrix_multi_matrix(R_M,R_X,R_T,3)
    end select

end subroutine cal_rot_matrix_angle_3d

subroutine matrix_multi_matrix(mat_0, mat_1, mat_2,dm)

    double precision :: mat_0(:,:), mat_1(:,:), mat_2(:,:)
    integer          :: dm

    ! local variables
    integer i,j,k

    do i=1,dm
        do j=1,dm
            mat_0(i,j) = 0.d0
        end do
    end do


    do i=1,dm
        do j=1,dm
            do k=1,dm
               mat_0(i,j) = mat_0(i,j) + mat_1(i,k) * mat_2(k,j)
            end do
        end do
    end do


end subroutine matrix_multi_matrix

subroutine gen_random(mla,noise,pfpara,time_step,the_bc_tower)

  type(ml_layout), intent(in   ) :: mla
  type(multifab) , intent(inout) :: noise(:)
  type(pf_para),   intent(in   ) :: pfpara
  integer,         intent(in   ) :: time_step
  type(bc_tower) , intent(in   ) :: the_bc_tower

  ! local variables
  integer :: lo(mla%dim), hi(mla%dim)
  integer :: dm, ng_p, i,n,nlevs

  real(kind=dp_t), pointer ::  noi(:,:,:,:)

  if(pfpara%solve_flow_only == 1) return

  dm   = mla%dim
  nlevs = mla%nlevel
  ng_p = noise(1)%ng

  do n=1,nlevs
  do i=1,nfabs(noise(n))
    noi  => dataptr(noise(n),i)
    lo = lwb(get_box(noise(n),i))
    hi = upb(get_box(noise(n),i))

    select case(dm)
    case (2)
      call gen_norm_random_2d(noi(:,:,1,1),ng_p,lo,hi,i,pfpara%noiamp,pfpara%no_noise,&
                              time_step,pfpara%period_type, pfpara%period_p0, pfpara%period_NT)

    case (3)
      call gen_norm_random_3d(noi(:,:,:,1),ng_p,lo,hi,i,pfpara%noiamp,pfpara%no_noise,&
                              time_step,pfpara%period_type, pfpara%period_p0, pfpara%period_NT)

    end select
  end do
  end do

  call ml_restrict_and_fill(nlevs, noise, mla%mba%rr, the_bc_tower%bc_tower_array,1,1,1)
  
  !if (nlevs .eq. 1) then

    ! Do the same for the orientation field
  !  call multifab_fill_boundary_c(noise(nlevs),1,1)
  !  call multifab_physbc(noise(nlevs),1,1,1,the_bc_tower%bc_tower_array(nlevs))

  !else

    ! the loop over nlevs must count backwards to make sure the finer grids are done first
  !  do n=nlevs,2,-1
      ! set level n-1 data to be the average of the level n data covering it
  !    call ml_cc_restriction_c(noise(n-1),1,noise(n),1,mla%mba%rr(n-1,:),1)
  !    call multifab_fill_ghost_cells(noise(n),noise(n-1),noise(n)%ng,mla%mba%rr(n-1,:), &
  !                                       the_bc_tower%bc_tower_array(n-1), &
  !                                       the_bc_tower%bc_tower_array(n), &
  !                                       1,1,1)

  !  end do

  !end if

end subroutine gen_random

! This routine is used to generate the random numbers (Gaussian) for the uncertanty
subroutine gen_norm_random_2d(noise,ng,lo,hi,nbox,noi_amp,no_noise, time_step,period_type, period_p0, period_NT)

  integer           :: lo(2),hi(2),ng,nbox,no_noise,period_type,time_step, period_NT
  double precision  :: noise(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
  double precision  :: noi_amp,period_p0
  real ( kind = 8 ) :: r8_normal_01

  ! local variables
  integer            :: i,j
  integer (kind = 4) :: seed
  integer            :: count_0,count_rate, count_max ! for system clock use

  double precision   :: d_PI = 3.14159265358979d0, var_flag = 0.d0

  select case(no_noise)
  case (1)

    select case(period_type)
    case (1) ! sinowave
        var_flag = period_p0 * sin(2.d0*d_PI* dble(time_step)/dble(period_NT) )

    case (2)
        if( mod(floor(dble(time_step)/dble(period_NT)),2) .eq. 0) then
            var_flag = -period_p0
        else
            var_flag = period_p0
        end if
    end select

    ! do the loop and set all noise to be zero
    do j=lo(2),hi(2)
      do i=lo(1),hi(1)
        noise(i,j) = var_flag
      end do
    end do

  case (2)

    ! Get the system clock information, used to relate to the random number generator
    call system_clock(count_0, count_rate, count_max)

    ! Calculate the seed
    seed = nbox + count_0 * 1000.d0 / count_rate

    select case(period_type)
    case (1) ! sinowave
        var_flag = period_p0 * sin(2.d0*d_PI* dble(time_step)/dble(period_NT) )

    case (2)
        if( mod(floor(dble(time_step)/dble(period_NT)),2) .eq. 0) then
            var_flag = -period_p0
        else
            var_flag = period_p0
        end if
    end select

    ! do the loop and set all noise to be zero
    do j=lo(2),hi(2)
      do i=lo(1),hi(1)
        noise(i,j) = var_flag + r8_normal_01 ( seed ) * noi_amp
      end do
    end do

  case default
    ! Get the system clock information, used to relate to the random number generator
    call system_clock(count_0, count_rate, count_max)

    ! Calculate the seed
    seed = nbox + count_0 * 1000.d0 / count_rate

    ! do the loop and set noise to be Gaussian
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
         noise(i,j) = r8_normal_01 ( seed ) * noi_amp
       end do
    end do

  end select

end subroutine gen_norm_random_2d

subroutine gen_norm_random_3d(noise,ng,lo,hi,nbox,noi_amp,no_noise, time_step,period_type, period_p0, period_NT)

  integer           :: lo(3),hi(3),ng,nbox,no_noise,period_type,time_step, period_NT
  double precision  :: noise(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
  double precision  :: noi_amp,period_p0
  real ( kind = 8 ) :: r8_normal_01

  ! local variables
  integer            :: i,j,k
  integer (kind = 4) :: seed
  integer            :: count_0,count_rate, count_max ! for system clock use

  double precision   :: d_PI = 3.14159265358979d0, var_flag

  select case(no_noise)
  case (1)

    select case(period_type)
    case (1) ! sinowave
        var_flag = period_p0 * sin(2.d0*d_PI* dble(time_step)/dble(period_NT) )

    case (2)
        if( mod(floor(dble(time_step)/dble(period_NT)),2) .eq. 0) then
            var_flag = -period_p0
        else
            var_flag = period_p0
        end if
    end select

    ! do the loop
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          noise(i,j,k) = var_flag
        end do
      end do
    end do

  case (2)

    ! Get the system clock information, used to relate to the random number generator
    call system_clock(count_0, count_rate, count_max)

    ! Calculate the seed
    seed = nbox + count_0 * 1000.d0 / count_rate

    select case(period_type)
    case (1) ! sinowave
        var_flag = period_p0 * sin(2.d0*d_PI* dble(time_step)/dble(period_NT) )

    case (2)
        if( mod(floor(dble(time_step)/dble(period_NT)),2) .eq. 0) then
            var_flag = -period_p0
        else
            var_flag = period_p0
        end if
    end select

    ! do the loop
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          noise(i,j,k) = var_flag + r8_normal_01 ( seed ) * noi_amp
        end do
      end do
    end do

  case default
    ! Get the system clock information, used to relate to the random number generator
    call system_clock(count_0, count_rate, count_max)

    ! Calculate the seed
    seed = nbox + count_0 * 1000.d0 / count_rate

    ! do the loop
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          noise(i,j,k) = r8_normal_01 ( seed ) * noi_amp
        end do
      end do
    end do

  end select

end subroutine gen_norm_random_3d

! This subroutine is used to update the variable specified, since currently
subroutine update_rkpf_2d(pf_old,pf_new,rkmat,ng_p,lo,hi,dt,var_h,var_l)

integer           :: lo(2),hi(2),ng_p
double precision  :: pf_old(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
double precision  :: pf_new(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
double precision  :: rkmat(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
double precision  :: dt,var_h,var_l

! local variables
integer           :: i,j

! do the loops
do j=lo(2),hi(2)
  do i=lo(1),hi(1)

    pf_new(i,j) = pf_old(i,j) + dt * rkmat(i,j)

    ! cut off the bouds
    if (pf_new(i,j) > var_h) then
       pf_new(i,j) = var_h
    else if(pf_new(i,j) < var_l) then
       pf_new(i,j) = var_l
    end if

  end do
end do

end subroutine update_rkpf_2d

subroutine update_rkpf_3d(pf_old,pf_new,rkmat,ng_p,lo,hi,dt,var_h,var_l)

  integer           :: lo(3),hi(3),ng_p
  double precision  :: pf_old(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: pf_new(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: rkmat(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: dt,var_h,var_l

  ! local variables
  integer           :: i,j,k

  ! do the loops
  do k=lo(3),hi(3)
    do j=lo(2),hi(2)
      do i=lo(1),hi(1)

        pf_new(i,j,k) = pf_old(i,j,k) + dt * rkmat(i,j,k)

        ! cut off the bouds
        if (pf_new(i,j,k) > var_h) then
           pf_new(i,j,k) = var_h
        else if(pf_new(i,j,k) < var_l) then
           pf_new(i,j,k) = var_l
        end if

      end do
    end do
  end do

end subroutine update_rkpf_3d

subroutine update_rkth_2d(pf_old,pf_new,rkmat,ng_p,lo,hi,dt,var_h,var_l)

integer           :: lo(2),hi(2),ng_p
double precision  :: pf_old(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
double precision  :: pf_new(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
double precision  :: rkmat(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
double precision  :: dt,var_h,var_l

! local variables
integer           :: i,j

! do the loops
do j=lo(2),hi(2)
  do i=lo(1),hi(1)
    pf_new(i,j) = pf_old(i,j)
  end do
end do

end subroutine update_rkth_2d

subroutine update_rkth_3d(pf_old,pf_new,rkmat,ng_p,lo,hi,dt,var_h,var_l)

integer           :: lo(3),hi(3),ng_p
double precision  :: pf_old(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
double precision  :: pf_new(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
double precision  :: rkmat(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
double precision  :: dt,var_h,var_l

! local variables
integer           :: i,j,k

! do the loops
do k=lo(3),hi(3)
  do j=lo(2),hi(2)
    do i=lo(1),hi(1)
      pf_new(i,j,k) = pf_old(i,j,k)
    end do
  end do
end do

end subroutine update_rkth_3d

subroutine follow_tip_info(un,time, mla,pf,ml_dx)

  type(ml_layout), intent(in   ) :: mla
  type(multifab) , intent(in   ) :: pf(:)
  real(kind=dp_t), intent(in   ) :: ml_dx(:), time
  integer,         intent(in   ) :: un

  ! local variables
  integer :: lo(mla%dim), hi(mla%dim),n
  integer :: dm, ng_p, i,j_tip = 0
  logical :: box_marked = .false.

  double precision  :: dx,tip_info(3)
  double precision, allocatable :: root_recv(:)
  real(kind=dp_t), pointer ::  pd(:,:,:,:)
  real(kind=dp_t) :: tip_position

  if ( parallel_IOProcessor() ) then
     allocate( root_recv( 3*parallel_nprocs() ) )
  else
     allocate( root_recv(3) )
  end if

  tip_info(1) = -1.d0
  tip_info(2) = 0.d0
  tip_info(3) = 0.d0

  tip_position = 0.d0

  dm   = mla%dim
  ng_p = pf(1)%ng
  n    = mla%nlevel
  dx = ml_dx(n)

  do i=1,nfabs(pf(n))
    pd  => dataptr(pf(n),i)
    lo = lwb(get_box(pf(n),i))
    hi = upb(get_box(pf(n),i))

    select case(dm)
    case (2)
      call follow_tip_2d(tip_position,box_marked,j_tip,pd(:,:,1,1),ng_p,lo,hi,dx)

    case (3)
      call follow_tip_3d(tip_position,box_marked,j_tip,pd(:,:,:,1),ng_p,lo,hi,dx)

    end select

    if(box_marked) then
       tip_info(1) = 1.d0
       tip_info(2) = dble(j_tip)
       tip_info(3) = tip_position
       box_marked = .false.
    end if

  end do

  !gather information to the root process
  call parallel_gather(tip_info, root_recv, 3)

  if ( parallel_IOProcessor() ) then
      do i = 1,parallel_nprocs()
         if( root_recv((i-1)*3+1) .gt. 0 ) then
             j_tip = int(root_recv((i-1)*3+2))
             tip_position = root_recv(i*3)
         end if
      end do

      write(unit=un, fmt='(es20.5e3, es20.5e3)') time,tip_position
  end if

  ! free the memory
  deallocate(root_recv)

end subroutine follow_tip_info

subroutine follow_tip_info2(istep, dirname,tip_position,cal_meancurv,mla,pf,ml_dx,N_range,&
                           plot_para_num,write_para_file,ratio_select)

  real(kind=dp_t), intent(  out) :: tip_position,cal_meancurv
  type(ml_layout), intent(in   ) :: mla
  type(multifab) , intent(in   ) :: pf(:)
  real(kind=dp_t), intent(in   ) :: ml_dx(:),ratio_select
  integer,         intent(in   ) :: N_range,istep,plot_para_num,write_para_file
  character(len=*), intent(in)   :: dirname

  ! local variables
  integer :: lo(mla%dim), hi(mla%dim),n
  integer :: dm, ng_p, i,j,j_tip = 0,j_min_ck,pb_count=0,n_ptnum = 0,nrows = 0,n_para_count=0
  logical :: box_marked = .false.
  double precision  :: dx,tip_info(3),pb_volume
  double precision, allocatable :: root_recv(:),parabolic_vec(:)!,collect_para_vec(:)
  !double precision, allocatable :: parab_pts(:,:)
  real(kind=dp_t), pointer ::  pd(:,:,:,:)

  ! pb_count is the number of points in each parabolic line
  ! The offset 2 is the flag and number of points
  pb_count = 2*N_range + 2

  if ( parallel_IOProcessor() ) then
     allocate( root_recv( 3*parallel_nprocs() ) )
     !allocate( collect_para_vec( pb_count * parallel_nprocs() ))
  else
     allocate( root_recv(3) )
     !allocate( collect_para_vec( 3 ) )
     !allocate( parab_pts(3, 2) )
  end if
  allocate( parabolic_vec( pb_count ) )
  !allocate( parab_pts(3, 2) )

  tip_info(1) = -1.d0
  tip_info(2) = 0.d0
  tip_info(3) = 0.d0

  tip_position = 0.d0
  cal_meancurv = 0.d0

  dm   = mla%dim
  ng_p = pf(1)%ng
  n    = mla%nlevel
  dx = ml_dx(n)

  do i=1,nfabs(pf(n))
    pd  => dataptr(pf(n),i)
    lo = lwb(get_box(pf(n),i))
    hi = upb(get_box(pf(n),i))

    select case(dm)
    case (2)
      call follow_tip_2d(tip_position,box_marked,j_tip,pd(:,:,1,1),ng_p,lo,hi,dx)

    case (3)
      call follow_tip_3d(tip_position,box_marked,j_tip,pd(:,:,:,1),ng_p,lo,hi,dx)

    end select

    if(box_marked) then
       tip_info(1) = 1.d0
       tip_info(2) = dble(j_tip)
       tip_info(3) = tip_position
       box_marked = .false.
    end if

  end do

  !gather information to the root process
  call parallel_gather(tip_info, root_recv, 3)

  if ( parallel_IOProcessor() ) then
      do i = 1,parallel_nprocs()
         if( root_recv((i-1)*3+1) .gt. 0 ) then
             j_tip = int(root_recv((i-1)*3+2))
             tip_position = root_recv(i*3)
         end if
      end do
  end if

  ! Now broadcast the j index of the tip to every process
  call parallel_bcast(j_tip)

  ! We use the j_min and N_setup to pinpoint the range that needs checking
  ! The checking loop starts from max(lo(2),j_min-N_setup) and min(j_min, hi(2))
  ! for each j in [st,ed], check if pf(i,j)>0 and pf(i+1,j)<0 or visvers, if so, record i
  ! Add a storage vector with the i and count++
  ! The final format of the vector in each process is [flag, count, i1, j1, i2, j2 ...]
  ! Use parallel gather to retrieve info from every process to the root
  ! Retrive the vector info at the root and store the data to the file
  j_min_ck = j_tip - N_range + 1

  n_para_count = 0
  pb_volume = 0.d0

  parabolic_vec(1) = -1.d0
  parabolic_vec(2) = 0.d0

  do i=1,nfabs(pf(n))
    pd  => dataptr(pf(n),i)
    lo = lwb(get_box(pf(n),i))
    hi = upb(get_box(pf(n),i))

    select case(dm)
    case (2)
      !call rt_parabolic_2d(parabolic_vec,j_min_ck,j_tip,pd(:,:,1,1),ng_p,lo,hi,dx,n_para_count)
      call rt_parabolic_volume_2d(pb_volume,j_min_ck,j_tip,pd(:,:,1,1),ng_p,lo,hi,dx)
    case (3)
      !call rt_parabolic_3d(parabolic_vec,j_min_ck,j_tip,pd(:,:,:,1),ng_p,lo,hi,dx,n_para_count)
      call rt_parabolic_volume_3d(pb_volume,j_min_ck,j_tip,pd(:,:,:,1),ng_p,lo,hi,dx)

    end select

  end do

  ! collect run_time from each processor and store the maximum
  call parallel_reduce(cal_meancurv, pb_volume, MPI_SUM, &
                         proc = parallel_IOProcessorNode())

  !gather information to the root process
  !call parallel_gather(parabolic_vec, collect_para_vec, pb_count)

  ! Retrive info from the collected data
  !if ( parallel_IOProcessor() ) then

  !    n_ptnum = 0
  !    nrows = 0

  !    do i = 1,parallel_nprocs()
  !       if( collect_para_vec( (i-1)*pb_count + 1 ) .gt. 0 ) then
  !           n_ptnum = n_ptnum + collect_para_vec( (i-1)*pb_count + 2 )
  !       end if
  !    end do

  !    allocate( parab_pts(n_ptnum, 2) )

  !    do i = 1,parallel_nprocs()
  !       if( collect_para_vec( (i-1)*pb_count + 1 ) .gt. 0 ) then
  !           n_ptnum = collect_para_vec( (i-1)*pb_count + 2 )
  !           do j =1,n_ptnum
  !              parab_pts(j+nrows,1) = collect_para_vec( (i-1)*pb_count + 2 + (j-1)*2 + 1 )
  !              parab_pts(j+nrows,2) = collect_para_vec( (i-1)*pb_count + 2 + (j-1)*2 + 2 )
  !           end do
  !           nrows = nrows + n_ptnum
  !       end if
  !    end do

      ! calculate the curvature
  !    call cal_mean_curvature(cal_meancurv, parab_pts,nrows,0.5d0*dx,tip_position,ratio_select,dm)

  !    if ( (write_para_file .eq. 1) .and. ( mod(istep,plot_para_num) .eq. 0 ) ) then
         ! write out the parabolic files in the PRB folder
  !       call write_parafile(dirname,parab_pts,0.5d0*dx,tip_position, nrows,istep)

  !    end if

  !else
  !  allocate( parab_pts(2, 2) )

  !end if

  ! free the memory
  deallocate(root_recv,parabolic_vec)!,collect_para_vec,parab_pts)

end subroutine follow_tip_info2

subroutine rt_parabolic_2d(pb_vec,j_min_def,j_max_def,pf,ng,lo,hi,dx,n_para_count)

integer           :: lo(2),hi(2),ng, j_min_def,j_max_def,n_para_count
double precision  :: pf(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
double precision  :: dx,pb_vec(:)

! local variables
integer           :: i,j,j_min,j_max
double precision  :: tip_position_x = 0.d0, tip_position_y = 0.d0

! set up the range
j_min = max(j_min_def, lo(2))
j_max = min(j_max_def, hi(2))

! do the loops
do j=j_min,j_max
   do i=lo(1),hi(1)
      ! if the sign of the phase field swaps that means we track down the tip position
      if( ( pf(i,j) .ge. 0.d0 .and. pf(i+1,j) .le. 0.d0 ) .or. ( pf(i,j) .le. 0.d0 .and. pf(i+1,j) .ge. 0.d0 ) ) then
         tip_position_x = (dble(i)+0.5d0) * dx - pf(i,j) * dx/( pf(i+1,j) - pf(i,j) )
         tip_position_y = (dble(j)+0.5d0) * dx
         n_para_count = n_para_count + 1
         pb_vec(1) = 1.d0
         pb_vec(2) = dble(n_para_count)
         pb_vec(2+2*n_para_count-1) = tip_position_x
         pb_vec(2+2*n_para_count)   = tip_position_y
      end if
   end do
end do

end subroutine rt_parabolic_2d

subroutine rt_parabolic_3d(pb_vec,j_min_def,j_max_def,pf,ng,lo,hi,dx,n_para_count)

integer           :: lo(3),hi(3),ng, j_min_def,j_max_def,n_para_count
double precision  :: pf(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
double precision  :: dx,pb_vec(:)

! local variables
integer           :: i,j,k,j_min,j_max
double precision  :: tip_position_x = 0.d0, tip_position_y = 0.d0

! set up the range
j_min = max(j_min_def, lo(2))
j_max = min(j_max_def, hi(2))

if( lo(3) .eq. 0 ) then
    k = lo(3)
    ! do the loops
    do j=j_min,j_max
      do i=lo(1),hi(1)
         ! if the sign of the phase field swaps that means we track down the tip position
         if( ( pf(i,j,k) .ge. 0.d0 .and. pf(i+1,j,k) .le. 0.d0 ) .or. ( pf(i,j,k) .le. 0.d0 .and. pf(i+1,j,k) .ge. 0.d0 ) ) then
            tip_position_x = (dble(i)+0.5d0) * dx - pf(i,j,k) * dx/( pf(i+1,j,k) - pf(i,j,k) )
            tip_position_y = (dble(j)+0.5d0) * dx
            n_para_count = n_para_count + 1
            pb_vec(1) = 1.d0
            pb_vec(2) = dble(n_para_count)
            pb_vec(2+2*n_para_count-1) = tip_position_x
            pb_vec(2+2*n_para_count)   = tip_position_y
         end if
      end do
   end do

end if

end subroutine rt_parabolic_3d

subroutine follow_tip_2d(pos_tip,box_marked,j_min_tip,pf,ng,lo,hi,dx)

integer           :: lo(2),hi(2),ng, j_min_tip
double precision  :: pf(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
logical           :: box_marked
double precision  :: dx,pos_tip

! local variables
integer           :: i,j

if(lo(1) .eq. 0) then
   ! the first line
   i = lo(1)

   ! do the loops
   do j=lo(2),hi(2)
      ! if the sign of the phase field swaps that means we track down the tip position
      if( ( pf(i,j) .ge. 0.d0 .and. pf(i,j+1) .le. 0.d0 ) .or. ( pf(i,j) .le. 0.d0 .and. pf(i,j+1) .ge. 0.d0 ) ) then
         pos_tip = (dble(j)+0.5d0) * dx - pf(i,j) * dx/( pf(i,j+1) - pf(i,j) )
         !cal_meancurv = 0.d0
         j_min_tip = j
         box_marked = .true.
      end if
   end do

end if

end subroutine follow_tip_2d

subroutine follow_tip_3d(pos_tip,box_marked,j_min_tip,pf,ng,lo,hi,dx)

integer           :: lo(3),hi(3),ng,j_min_tip
double precision  :: pf(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
logical           :: box_marked
double precision  :: dx,pos_tip

! local variables
integer           :: i,j,k

if( ( lo(1) .eq. 0 ) .and. ( lo(3) .eq. 0 ) ) then
   ! the first line
   i = lo(1)
   k = lo(3)

   ! do the loops
   do j=lo(2),hi(2)
      ! if the sign of the phase field swaps that means we track down the tip position
      if( ( pf(i,j,k) .ge. 0.d0 .and. pf(i,j+1,k) .le. 0.d0 ) .or. ( pf(i,j,k) .le. 0.d0 .and. pf(i,j+1,k) .ge. 0.d0 ) ) then
         pos_tip = (dble(j)+0.5d0) * dx - pf(i,j,k) * dx/( pf(i,j+1,k) - pf(i,j,k) )
         !cal_meancurv = 0.d0
         j_min_tip = j
         box_marked = .true.
      end if
   end do

end if


end subroutine follow_tip_3d

subroutine write_tkfile(dirname,tk_info, nlength)
    real(dp_t)    , intent(in   ) :: tk_info(:,:)
    integer       , intent(in   ) :: nlength
    character(len=*), intent(in)   :: dirname

    ! local variables
    integer       :: un = 1022, i
    character(len=128) :: tkfilename = "TK_INFO.dat"

    ! open a file name TK_INFO.dat, replace the file everytime write the tip information
    open(unit=un, file = trim(dirname) // "/" // trim(tkfilename), &
         form = "formatted", access = "sequential", &
         status = "replace", action = "write")

    ! write all the information down
    do i=1,nlength
       write(unit=un, fmt='(3es27.17e3)') tk_info(i,1), tk_info(i,2), tk_info(i,3)
    end do

    ! close the file after processing
    close(unit=un)


end subroutine write_tkfile

subroutine write_parafile(dirname,para_vec,tip_x,tip_y, nlength,filenum)
    real(dp_t)    , intent(in   )  :: para_vec(:,:),tip_x,tip_y
    integer       , intent(in   )  :: nlength,filenum
    character(len=*), intent(in)   :: dirname

    ! local variables
    integer       :: un = 1023, i
    character(len=128) :: parafilename

    ! define the name of the plotfile that will be written
    write(unit=parafilename,fmt='("PARA_",i6.6)') filenum

    ! open a file name TK_INFO.dat, replace the file everytime write the tip information
    open(unit=un, file = trim(dirname) // "/" // trim(parafilename), &
         form = "formatted", access = "sequential", &
         status = "replace", action = "write")

    ! write all the information down
    write(unit=un, fmt='(3es27.17e3)') 0.d0, tip_y, 0.d0

    ! write others
    do i=1,nlength
       write(unit=un, fmt='(3es27.17e3)') para_vec(i,1)-tip_x, para_vec(i,2), para_vec(i,2)-tip_y
    end do

    ! close the file after processing
    close(unit=un)


end subroutine write_parafile

subroutine cal_mean_curvature(var_curve, para_vec,nlength,tip_x,tip_y,ratio_select,dm)

    double precision :: var_curve, para_vec(:,:)
    integer          :: nlength,dm
    double precision :: tip_x, tip_y,ratio_select

    ! local variables
    integer   ::  i,nlength_select
    double precision :: vec_mod(nlength,2)
    double precision :: var_x2 = 0.d0, var_x3 = 0.d0, var_x4 = 0.d0
    double precision :: var_yx2 = 0.d0, var_yx1 = 0.d0, p_a, p_b

    ! Select the point for fitting
    nlength_select = ratio_select*nlength

    ! reset the points
    do i =1,nlength_select
        vec_mod(i,1) = para_vec(i,1) - tip_x
        vec_mod(i,2) = para_vec(i,2) - tip_y
    end do

    ! calculate rhe curvature
    do i = 1,nlength_select
        var_x2 = var_x2 + vec_mod(i,1)**2
        var_x3 = var_x3 + vec_mod(i,1)**3
        var_x4 = var_x4 + vec_mod(i,1)**4
        var_yx2 = var_yx2 + vec_mod(i,2)* vec_mod(i,1)**2
        var_yx1 = var_yx1 + vec_mod(i,2)* vec_mod(i,1)
    end do

    p_b = (var_yx2/var_x4 - var_yx1/var_x3)/( var_x3/var_x4 - var_x2/var_x3 )
    p_a = (var_yx1 - var_x2*p_b)/var_x3

    var_curve = 4.d0 * p_a**2 /( 1.d0 + p_b**2 )**3

    select case(dm)
    case (2)
       ! calculate the curvature
       var_curve = sqrt(var_curve)
    case(3)
        var_curve = 2.d0*sqrt(var_curve)
    end select

end subroutine cal_mean_curvature

subroutine write_status_file(un,dt,dx,regrid_time, start_time, istep, nstep,ada_regriding, &
                             regrid_int_num,v_mag_max, Rate_Cooling, min_temp)

    !character(len=*), intent(in)   :: dirname
    real(dp_t)    , intent(in   )  :: regrid_time, start_time
    real(dp_t)    , intent(in   )  :: dt,dx,v_mag_max,Rate_Cooling, min_temp
    integer       , intent(in   )  :: istep, nstep,ada_regriding,un,regrid_int_num

    !integer          :: un, i
    integer          :: i
    real(dp_t)       :: elapsed_time,regrid_time_IOproc, elapsed_time_IOproc
    real(dp_t)       :: temp_cal
    !character(len=128) :: status_filename

    elapsed_time = parallel_wtime() - start_time

    call parallel_barrier()

    ! collect run_time from each processor and store the maximum
    call parallel_reduce(regrid_time_IOproc, regrid_time, MPI_MAX, &
                         proc = parallel_IOProcessorNode())

    ! collect run_time from each processor and store the maximum
    call parallel_reduce(elapsed_time_IOproc, elapsed_time, MPI_MAX, &
                         proc = parallel_IOProcessorNode())

    if ( parallel_IOProcessor() ) then

        temp_cal = min_temp - dt * dble(istep) * Rate_Cooling

        !un = 1019

        ! define the name of the plotfile that will be written
        !write(unit=status_filename,fmt='("CAL_STATUS.dat")' )

        ! open a file name TK_INFO.dat, replace the file everytime write the tip information
        !open(unit=un, file = trim(dirname) // "/" // trim(status_filename), &
        !     form = "formatted", access = "sequential", &
        !     status = "replace", action = "write")

        !write(unit=un, fmt='("----------------Calculation status description-----------------")')
        !write(unit=un, fmt='("Finished steps             :", i8.7)') istep
        if(ada_regriding .eq. 1) then
            write(unit=un, fmt='(i8.7, es10.2e3, es10.2e3, i15.3, es15.2e3, es15.2e3)') &
                              istep, elapsed_time_IOproc, regrid_time_IOproc, &
                              regrid_int_num,v_mag_max, temp_cal
        else
           write(unit=un, fmt='(i8.7, es10.2e3, es10.2e3, es10.2e3)') &
                              istep, elapsed_time_IOproc, regrid_time_IOproc, temp_cal
        end if
        !write(unit=un, fmt='("Total steps                :", i8.7)') nstep
        !write(unit=un, fmt='("Elapsed time (s)           :", es10.2e3)') elapsed_time_IOproc
        !write(unit=un, fmt='("Regrid  time (s)           :", es10.2e3)') regrid_time_IOproc
        !write(unit=un, fmt='("Elapsed time per step (s)  :", es10.2e3)') elapsed_time_IOproc / dble(istep)
        !write(unit=un, fmt='("Regrid  time per step (s)  :", es10.2e3)') regrid_time_IOproc / dble(istep)
        !write(unit=un, fmt='("---------------------------------------------------------------")')

        !if(ada_regriding .eq. 1) then
        !   write(unit=un, fmt='("----------------Regriding status description-----------------")')
        !   write(unit=un, fmt='("Time step                  :", es10.2e3)') dt
        !   write(unit=un, fmt='("Reference velocity         :", es10.2e3)') dx/dt
        !   write(unit=un, fmt='("regriding time        regriding number        max_v_magnitude")')

        !  do i=1,regrid_count
        !      write(unit=un, fmt='(es15.2e3, i15.3, es15.2e3 )') &
        !                       regrid_int_temp(i,1),int( regrid_int_temp(i,2) ),regrid_int_temp(i,3)
        !   end do

        !   write(unit=un, fmt='("---------------------------------------------------------------")')

        !end if

        ! close the file after processing
        !close(unit=un)

    end if


end subroutine write_status_file

subroutine write_status_file2(un,start_time,rigrid_time, istep)

   real(dp_t)    , intent(in   )  :: start_time,rigrid_time
   integer       , intent(in   )  :: istep, un

   real(dp_t)       :: elapsed_time
   real(dp_t)       :: elapsed_time_IOproc, regrid_time_IOproc

   elapsed_time = parallel_wtime() - start_time

   call parallel_barrier()

   ! collect run_time from each processor and store the maximum
   !call parallel_reduce(regrid_time_IOproc, rigrid_time, MPI_MAX, &
   !                      proc = parallel_IOProcessorNode())

   ! collect run_time from each processor and store the maximum
   !call parallel_reduce(elapsed_time_IOproc, elapsed_time, MPI_MAX, &
   !                     proc = parallel_IOProcessorNode())

!   if ( parallel_IOProcessor() ) then

      write(un, '(i8.7, es10.2e3, es10.2e3)') &
                           istep, elapsed_time, rigrid_time
                           !istep, elapsed_time_IOproc, regrid_time_IOproc

 !  end if


end subroutine write_status_file2

subroutine retrieve_max_velocity(istep,v_mag_max,mla,phi_old,phi_new,dt,dx,thers_gradient)
    real(kind=dp_t), intent(inout) :: v_mag_max
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: phi_old(:),phi_new(:)
    real(kind=dp_t), intent(in   ) :: dt,dx(:),thers_gradient
    integer        , intent(in   ) :: istep

    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: dm, ng_p, i,j, n,nlevs
    double precision :: v_mag_max_sub,v_mag_store

    ! debug
    !double precision :: locate_pos(mla%dim + 1)
    !double precision :: recv_locate_pos( parallel_nprocs()*(mla%dim + 1) )

    real(kind=dp_t), pointer ::  pfo(:,:,:,:)
    real(kind=dp_t), pointer ::  pfn(:,:,:,:)

    dm   = mla%dim
    nlevs = mla%nlevel
    ng_p = phi_old(1)%ng

    v_mag_max_sub = v_mag_max
    v_mag_store = v_mag_max

    n = nlevs

    ! calculate the RK matrices
    do i=1,nfabs(phi_old(n))
       pfo  => dataptr(phi_old(n),i)
       pfn  => dataptr(phi_new(n),i)
       lo = lwb(get_box(phi_old(n),i))
       hi = upb(get_box(phi_old(n),i))

       select case(dm)
       case (2)
       ! calculate the K1
       call cal_vmax_2d(v_mag_max_sub,pfn(:,:,1,1),pfo(:,:,1,1),ng_p,lo,hi,dx(n),dt,thers_gradient)

       case (3)
       call cal_vmax_3d(v_mag_max_sub,pfn(:,:,:,1),pfo(:,:,:,1),ng_p,lo,hi,dx(n),dt,thers_gradient)

       end select

    end do

    ! collect run_time from each processor and store the maximum
    call parallel_reduce(v_mag_max, v_mag_max_sub, MPI_MAX, &
                         proc = parallel_IOProcessorNode())

    !debug
    !gather information to the root process
    !call parallel_gather(locate_pos, recv_locate_pos, dm+1)

    !if ( parallel_IOProcessor() ) then
    !    v_mag_max_sub = v_mag_store
    !    do i=1,parallel_nprocs()
    !        if( recv_locate_pos( (i-1)*(dm+1) + 1 ) .gt. v_mag_max_sub) then
    !            v_mag_max_sub = recv_locate_pos( (i-1)*(dm+1) + 1 )
    !            do j=1,dm
    !               location_pos(j) = recv_locate_pos( (i-1)*(dm+1) + 1 + j )
    !            end do
    !            location_pos(4) = dble(istep)
    !        end if
    !    end do

    !end if

end subroutine retrieve_max_velocity

subroutine cal_vmax_2d(v_max,pfn,pfo,ng_p,lo,hi,dx,dt,thers_gradient)
  integer           :: lo(2),hi(2),ng_p
  double precision  :: pfn(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: pfo(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: v_max,dx,dt,thers_gradient

  ! local variables
  integer           :: i,j
  double precision  :: de1_x, de1_y, de1_t,v_mag=0.d0,mod_1,var_target

  ! do the loops
  do j=lo(2),hi(2)
    do i=lo(1),hi(1)

      de1_x = (pfo(i+1,j) - pfo(i-1,j))/(2.d0 * dx)
      de1_y = (pfo(i,j+1) - pfo(i,j-1))/(2.d0 * dx)

      mod_1 = sqrt(de1_x**2 + de1_y**2)
      var_target = mod_1*dx

      if (var_target > thers_gradient) then

         de1_t = (pfn(i,j) - pfo(i,j)) / dt
         v_mag = abs(de1_t)/mod_1

         if(v_mag .gt. v_max) then
             v_max = v_mag
             !locate_pos(1) = v_max
             !locate_pos(2) = dble(i)
             !locate_pos(3) = dble(j)
         end if

      end if

    end do
  end do

end subroutine cal_vmax_2d

subroutine cal_vmax_3d(v_max,pfn,pfo,ng_p,lo,hi,dx,dt,thers_gradient)
  integer           :: lo(3),hi(3),ng_p
  double precision  :: pfn(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: pfo(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: v_max,dx,dt,thers_gradient

  ! local variables
  integer           :: i,j,k
  double precision  :: de1_x, de1_y, de1_z, de1_t,v_mag=0.d0,mod_1,var_target

  ! do the loops
  do k=lo(3),hi(3)
    do j=lo(2),hi(2)
      do i=lo(1),hi(1)

        de1_x = (pfo(i+1,j,k) - pfo(i-1,j,k))/(2.d0 * dx)
        de1_y = (pfo(i,j+1,k) - pfo(i,j-1,k))/(2.d0 * dx)
        de1_z = (pfo(i,j,k+1) - pfo(i,j,k-1))/(2.d0 * dx)

        mod_1 = sqrt(de1_x**2 + de1_y**2 + de1_z**2)
        var_target = mod_1*dx

        if (var_target > thers_gradient) then

           de1_t = (pfn(i,j,k) - pfo(i,j,k)) / dt
           v_mag = abs(de1_t)/mod_1

           if(v_mag .gt. v_max) then
               v_max = v_mag
               !locate_pos(1) = v_max
               !locate_pos(2) = dble(i)
               !locate_pos(3) = dble(j)
               !locate_pos(4) = dble(k)
           end if

        end if

      end do
    end do
  end do

end subroutine cal_vmax_3d

subroutine recalculate_regrid_number(mla,regrid_int,v_mag_max,dt,dx_ml,regrid_amp_ref,rg_max,rg_min)

   type(ml_layout)   , intent(in   ) :: mla
   real(kind=dp_t)   , intent(in   ) :: dt,dx_ml(:)
   integer           , intent(  out) :: regrid_int
   real(kind=dp_t)   , intent(inout) :: v_mag_max
   integer           , intent(in   ) :: regrid_amp_ref,rg_max,rg_min

   ! local variables
   integer          :: nlevs
   double precision :: v_tip_ref, small_thers = 1.0E-8,dx

   nlevs = mla%nlevel
   dx = dx_ml(nlevs)

   if( parallel_IOProcessor() ) then

      v_tip_ref = dx/dt

      v_mag_max = max(v_mag_max, small_thers)

      regrid_int = regrid_amp_ref * int(v_tip_ref/v_mag_max)

      regrid_int = max(regrid_int, rg_min)
      regrid_int = min(regrid_int, rg_max)

   end if

   ! Now broadcast the j index of the tip to every process
   call parallel_bcast(regrid_int)

end subroutine recalculate_regrid_number

subroutine retrieve_curvature(mla,curv,phi,dx,thers_gradient)
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(  out) :: curv(:)
    type(multifab) , intent(in   ) :: phi(:)
    real(kind=dp_t), intent(in   ) :: dx(:),thers_gradient

    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: dm, ng_p, n,nlevs,i

    real(kind=dp_t), pointer ::  pf(:,:,:,:)
    real(kind=dp_t), pointer ::  pcv(:,:,:,:)

    dm   = mla%dim
    nlevs = mla%nlevel
    ng_p = phi(1)%ng

    n = nlevs

    ! calculate the RK matrices
    do i=1,nfabs(phi(n))
       pcv  => dataptr(curv(n),i)
       pf   => dataptr(phi(n),i)
       lo   = lwb(get_box(phi(n),i))
       hi   = upb(get_box(phi(n),i))

       select case(dm)
       case (2)
       ! calculate the K1
       call cal_curvature_2d(pcv(:,:,1,1),pf(:,:,1,1),ng_p,lo,hi,dx(n),thers_gradient)

       case (3)
       call cal_curvature_3d(pcv(:,:,:,1),pf(:,:,:,1),ng_p,lo,hi,dx(n),thers_gradient)

       end select

    end do

    ! collect run_time from each processor and store the maximum
    !call parallel_reduce(v_mag_max, v_mag_max_sub, MPI_MAX, &
    !                     proc = parallel_IOProcessorNode())


end subroutine retrieve_curvature

subroutine cal_curvature_2d(curv,pf,ng_p,lo,hi,dx,thers_gradient)
  integer           :: lo(2),hi(2),ng_p
  double precision  :: pf(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: curv(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: dx,thers_gradient

  ! local variables
  integer           :: i,j
  double precision  :: var_target, mod_1
  double precision  :: de1_x, de1_y, de2_xx, de2_yy, de2_xy

  ! do the loops
  do j=lo(2),hi(2)
    do i=lo(1),hi(1)

       curv(i,j) = 0.d0

       de1_x  = (pf(i+1,j) - pf(i-1,j))/(2.d0 * dx)
       de1_y  = (pf(i,j+1) - pf(i,j-1))/(2.d0 * dx)

       mod_1 = sqrt(de1_x**2 + de1_y**2)
       var_target = mod_1*dx

      if (var_target > thers_gradient) then

         de2_xx = (pf(i+1,j) - 2.d0*pf(i,j) + pf(i-1,j))/(dx * dx)
         de2_yy = (pf(i,j+1) - 2.d0*pf(i,j) + pf(i,j-1))/(dx * dx)
         de2_xy = (pf(i+1,j+1) + pf(i-1,j-1) - pf(i+1,j-1) - pf(i-1,j+1) )/(dx * dx)

         curv(i,j) = de2_xx * de1_y**2 + de2_yy * de1_x**2 - 2.d0*de1_x*de1_y*de2_xy
         curv(i,j) = -0.5*curv(i,j) / mod_1**3

      end if

    end do
  end do

end subroutine cal_curvature_2d

subroutine cal_curvature_3d(curv,pf,ng_p,lo,hi,dx,thers_gradient)
  integer           :: lo(3),hi(3),ng_p
  double precision  :: pf(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: curv(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: dx,thers_gradient

  ! local variables
  integer           :: i,j,k
  double precision  :: var_target, mod_1
  double precision  :: de1_x,de1_y,de1_z,de2_xx,de2_yy,de2_zz,de2_xy,de2_yz,de2_xz

  ! do the loops
  do k=lo(3),hi(3)
    do j=lo(2),hi(2)
      do i=lo(1),hi(1)

        curv(i,j,k) = 0.d0

        de1_x = (pf(i+1,j,k) - pf(i-1,j,k))/(2.d0 * dx)
        de1_y = (pf(i,j+1,k) - pf(i,j-1,k))/(2.d0 * dx)
        de1_z = (pf(i,j,k+1) - pf(i,j,k-1))/(2.d0 * dx)

        mod_1 = sqrt(de1_x**2 + de1_y**2 + de1_z**2)
        var_target = mod_1*dx

        if (var_target > thers_gradient) then

           de2_xx = (pf(i+1,j,k) - 2.d0*pf(i,j,k) + pf(i-1,j,k))/(dx * dx)
           de2_yy = (pf(i,j+1,k) - 2.d0*pf(i,j,k) + pf(i,j-1,k))/(dx * dx)
           de2_zz = (pf(i,j,k+1) - 2.d0*pf(i,j,k) + pf(i,j,k-1))/(dx * dx)

           de2_xy = (pf(i+1,j+1,k) + pf(i-1,j-1,k) - pf(i-1,j+1,k)- pf(i+1,j-1,k) )/(dx * dx)
           de2_xz = (pf(i+1,j,k+1) + pf(i-1,j,k-1) - pf(i-1,j,k+1)- pf(i+1,j,k-1) )/(dx * dx)
           de2_yz = (pf(i,j+1,k+1) + pf(i,j-1,k-1) - pf(i,j-1,k+1)- pf(i,j+1,k-1) )/(dx * dx)

           curv(i,j,k) = de2_xx * (de1_y**2 + de1_z**2) + de2_yy * (de1_x**2 + de1_z**2) + de2_zz * (de1_x**2 + de1_y**2)
           curv(i,j,k) = curv(i,j,k) - 2.d0*(de1_x*de1_y*de2_xy + de1_x*de1_z*de2_xz + de1_y*de1_z*de2_yz)
           curv(i,j,k) = -0.5*curv(i,j,k) / mod_1**3

        end if

      end do
    end do
  end do

end subroutine cal_curvature_3d

subroutine rt_parabolic_volume_2d(pb_volume,j_min_def,j_max_def,pf,ng,lo,hi,dx)

integer           :: lo(2),hi(2),ng, j_min_def,j_max_def
double precision  :: pf(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
double precision  :: dx,pb_volume

! local variables
integer           :: i,j,j_min,j_max

! set up the range
j_min = max(j_min_def, lo(2))
j_max = min(j_max_def, hi(2))

! do the loops
do j=j_min,j_max
   do i=lo(1),hi(1)
      ! if the sign of the phase field swaps that means we track down the tip position
      if( pf(i,j) .ge. 0.d0  ) then
         pb_volume = pb_volume + dx**2

      end if
   end do
end do

end subroutine rt_parabolic_volume_2d

subroutine rt_parabolic_volume_3d(pb_volume,j_min_def,j_max_def,pf,ng,lo,hi,dx)

integer           :: lo(3),hi(3),ng, j_min_def,j_max_def
double precision  :: pf(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
double precision  :: dx,pb_volume

! local variables
integer           :: i,j,k,j_min,j_max

! set up the range
j_min = max(j_min_def, lo(2))
j_max = min(j_max_def, hi(2))


! do the loops
do k = lo(3),hi(3)
    do j=j_min,j_max
      do i=lo(1),hi(1)
         ! if the sign of the phase field swaps that means we track down the tip position
         if(  pf(i,j,k) .ge. 0.d0 ) then
            pb_volume = pb_volume + dx**3
         end if
      end do
   end do
end do

end subroutine rt_parabolic_volume_3d

subroutine pf_para_set_flw_tau(tau, tau0, max_levels)

   real(kind=dp_t),intent(inout) :: tau(:)
   real(kind=dp_t),intent(in   ) :: tau0
   integer,        intent(in   ) :: max_levels

   integer :: n

   do n=1,max_levels
      tau(n) = (tau0 - 0.5d0)*(2.d0**(n-1))+0.5d0
   end do

end subroutine pf_para_set_flw_tau

subroutine pf_para_set_lbm(pfpara,dim,dx)
   type (pf_para),  intent(inout) :: pfpara
   integer,         intent(in   ) :: dim
   real(kind=dp_t), intent(in   ) :: dx  ! this is the dimensionless spatial step

   integer :: q
   double precision, parameter :: p0_2d = 4.d0/9.d0, p1_2d = 1.d0/9.d0, p2_2d=1.d0/36.d0
   double precision, parameter :: p0_3d = 1.d0/3.d0, p1_3d = 1.d0/18.d0, p2_3d=1.d0/36.d0
   double precision, parameter :: lbm_drag_hc = 2.757d0

   select case(dim)
    case (2)
      pfpara%dim_lbm = 9
    case (3)
      pfpara%dim_lbm = 19
   end select

   q = pfpara%dim_lbm

   allocate(pfpara%lbm_e(dim,q))
   allocate(pfpara%lbm_w(q))

   select case(dim)
    case (2)
      pfpara%lbm_e(1,1:q)=(/ 0, 1, 0,-1, 0, 1,-1,-1, 1 /)
      pfpara%lbm_e(2,1:q)=(/ 0, 0, 1, 0,-1, 1, 1,-1,-1 /)
      pfpara%lbm_w(1:q)=(/ p0_2d,p1_2d,p1_2d,p1_2d,p1_2d,p2_2d,p2_2d,p2_2d,p2_2d /) 
    case (3)
      pfpara%lbm_e(1,1:q)=(/ 0,1,-1,0, 0,0, 0,1, 1,-1,-1, 1,-1, 1,-1, 0, 0, 0, 0 /)
      pfpara%lbm_e(2,1:q)=(/ 0,0, 0,1,-1,0, 0,1,-1, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1 /)
      pfpara%lbm_e(3,1:q)=(/ 0,0, 0,0, 0,1,-1,0, 0, 0, 0, 1, 1,-1,-1, 1,-1, 1,-1 /)
      pfpara%lbm_w(1:q)=(/ p0_3d,p1_3d,p1_3d,p1_3d,p1_3d,p1_3d,p1_3d,p2_3d,p2_3d,p2_3d,&
                           p2_3d,p2_3d,p2_3d,p2_3d,p2_3d,p2_3d,p2_3d,p2_3d,p2_3d /)  
   end select

   pfpara%lbm_viscosity    = 1.d0/3.d0 * (pfpara%tau_LBM - 0.5d0)
   pfpara%lbm_ratio_dx_W0  = dx
   pfpara%lbm_drag_hc      = lbm_drag_hc

end subroutine pf_para_set_lbm

subroutine pf_set_uinlet(pfpara,istep)
  type (pf_para), intent(inout) :: pfpara
  integer       , intent(in  )  :: istep

  ! local variables
  integer          ::  nsteps_for_upper, dm, i
  double precision ::  ratio_climb

  select case(pfpara%flag_set_uinlet)
  
  case (0)

    !do dm=1,pfpara%dim 
    !  do i=1,2  
    !    if(pfpara%bc_flag(dm,i) == INLET) then
    !      pfpara%u_inlet(dm,i) = pfpara%u_inlet_def(dm,i)
    !    end if
    !  end do
    !end do

  case (1)

    nsteps_for_upper = pfpara%nsteps_for_upper
    ratio_climb = ( dble(istep) / dble(nsteps_for_upper) )

    do dm=1,pfpara%dim 
      do i=1,2  
        if(pfpara%bc_flag(dm,i) == INLET) then
        
          ! do the calculation
          pfpara%u_inlet(dm,i) = pfpara%uinlet_set_init + &
              ratio_climb * (pfpara%u_inlet_def(dm,i) - pfpara%uinlet_set_init)
        
          ! cut off the values for lower and upper
          if(pfpara%u_inlet(dm,i) .lt. pfpara%uinlet_set_init) then
             pfpara%u_inlet(dm,i) = pfpara%uinlet_set_init
          else if(pfpara%u_inlet(dm,i) .gt. pfpara%u_inlet_def(dm,i)) then
             pfpara%u_inlet(dm,i) = pfpara%u_inlet_def(dm,i)
          end if

        end if
      end do
    end do

  case (2)

    nsteps_for_upper = pfpara%nsteps_for_upper
    ratio_climb = dble(istep / nsteps_for_upper)

    do dm=1,pfpara%dim 
      do i=1,2  
        if(pfpara%bc_flag(dm,i) == INLET) then
        
          ! do the calculation
          pfpara%u_inlet(dm,i) = pfpara%uinlet_set_init + &
              ratio_climb * pfpara%uinlet_set_step
        
          ! cut off the values for lower and upper
          if(pfpara%u_inlet(dm,i) .lt. pfpara%uinlet_set_init) then
             pfpara%u_inlet(dm,i) = pfpara%uinlet_set_init
          else if(pfpara%u_inlet(dm,i) .gt. pfpara%u_inlet_def(dm,i)) then
             pfpara%u_inlet(dm,i) = pfpara%u_inlet_def(dm,i)
          end if

        end if
      end do
    end do

  end select

end subroutine pf_set_uinlet

subroutine retrieve_error_msg(mla,flw,thers_var,msg)
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: flw(:)
    real(kind=dp_t), intent(in   ) :: thers_var
    logical,         intent(inout) :: msg

    double precision :: nm_inf

    nm_inf = multifab_norm_inf_c(flw(1), 1, 2)

    nm_inf = abs(nm_inf)

    if(thers_var .le. nm_inf) then
      msg = .true.
    else
      msg = .false.
    end if

end subroutine retrieve_error_msg

subroutine recalculate_time_step(mla,dt,n_cal,dt0,flw,thres_safe,dx_f,dx_c,istep,time,un)
  real(kind=dp_t), intent(inout) :: dt
  integer,         intent(inout) :: n_cal
  type(ml_layout), intent(in   ) :: mla
  type(multifab) , intent(in   ) :: flw(:)
  real(kind=dp_t), intent(in   ) :: thres_safe,dx_f,dx_c,dt0,time
  integer,         intent(in   ) :: un,istep

  ! local variables
  double precision :: u_mag_max, u_mag_ref,data_bc(2)

  n_cal = 1
  u_mag_max = 0.d0
  u_mag_ref = thres_safe * dx_f / dx_c

  call retrieve_max_velocity_flw(mla,u_mag_max,flw)

  if ( parallel_IOProcessor() ) then
    data_bc(1) = u_mag_max/u_mag_ref + 1.d0
    data_bc(2) = dt0 / data_bc(1)
  end if

  call parallel_bcast(data_bc)
  
  n_cal = int(data_bc(1))
  dt    = data_bc(2)

  if ( parallel_IOProcessor() ) then  
     write(un, '(i8.7,es10.2e3,es10.2e3,es10.2e3,es10.2e3)') istep,time,u_mag_max,dt,dble(n_cal)
  end if

end subroutine recalculate_time_step 

subroutine retrieve_max_velocity_flw(mla,u_mag_max,flw)
   real(kind=dp_t), intent(inout) :: u_mag_max
   type(ml_layout), intent(in   ) :: mla
   type(multifab) , intent(in   ) :: flw(:)

   integer ::  ng_p,i,n,dm,nlevs
   integer :: lo(mla%dim), hi(mla%dim)
   real(kind=dp_t), pointer :: fl(:,:,:,:)
   real(kind=dp_t) :: u_mag_max_sep

   ng_p = flw(1)%ng
   dm   = mla%dim
   nlevs = mla%nlevel

   u_mag_max     = 0.d0
   u_mag_max_sep = 0.d0

   do n=1,nlevs
   do i=1,nfabs(flw(n))
       fl  => dataptr(flw(n),i)
       lo = lwb(get_box(flw(n),i))
       hi = upb(get_box(flw(n),i))

       select case(dm)
       case (2)
         call retrieve_max_velocity_flw_2d(u_mag_max_sep,fl(:,:,1,1),fl(:,:,1,2),ng_p,lo,hi)

       case (3)
         call retrieve_max_velocity_flw_3d(u_mag_max_sep,fl(:,:,:,1),fl(:,:,:,2),fl(:,:,:,3),ng_p,lo,hi)
       end select

    end do
    end do

    call parallel_reduce(u_mag_max, u_mag_max_sep, MPI_MAX, proc = parallel_IOProcessorNode())

end subroutine retrieve_max_velocity_flw

subroutine retrieve_max_velocity_flw_2d(u_mag_max,ux,uy,ng_p,lo,hi)
  integer           :: ng_p,lo(2),hi(2)
  double precision  :: ux(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: uy(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: u_mag_max

  integer :: i,j
  real(dp_t) :: var_temp

  do j=lo(2)-ng_p,hi(2)+ng_p
    do i=lo(1)-ng_p,hi(1)+ng_p
      var_temp = sqrt(ux(i,j)**2 + uy(i,j)**2)
      u_mag_max = max(u_mag_max,var_temp)
    end do
  end do

end subroutine retrieve_max_velocity_flw_2d

subroutine retrieve_max_velocity_flw_3d(u_mag_max,ux,uy,uz,ng_p,lo,hi)
  integer           :: ng_p,lo(3),hi(3)
  double precision  :: ux(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: uy(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: uz(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: u_mag_max

  integer :: i,j,k
  real(dp_t) :: var_temp

  do k=lo(3)-ng_p,hi(3)+ng_p
    do j=lo(2)-ng_p,hi(2)+ng_p
      do i=lo(1)-ng_p,hi(1)+ng_p
        var_temp = sqrt(ux(i,j,k)**2 + uy(i,j,k)**2+uz(i,j,k)**2)
        u_mag_max = max(u_mag_max,var_temp)
      end do
    end do
  end do

end subroutine retrieve_max_velocity_flw_3d

end module pf_utility_module




