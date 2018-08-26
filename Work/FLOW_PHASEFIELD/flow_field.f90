! ----------------------------------------------------------------------------------------------------------------
! the advance_flow_field_module, solving the N-S equations using LBM
! Originally created by ZPG, on Aug 18, 2016
! Modified by ZPG, on December 27, 2016
! use call advance_flow_field(mla,flw,phi_new,pfpara,the_bc_tower ):
! where mla     is the ml_layout structure,
!       flw     is the ml_multifab data for flow field
!       phi_new is the ml_multifab data for phase field
!       pfpara  is the pf_utilities data for parameter
!       the_bc_tower is the boundary condition tower data
!-----------------------------------------------------------------------------------------------------------------
module advance_flow_field_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use multifab_fill_ghost_module
  use pf_utility_module
  use ml_restrict_fill_module

  implicit none

  public  :: advance_flow_field

contains

subroutine advance_flow_field(mla,flw,flw_old,flw_new,surf,phi_new,pfpara,the_bc_tower,dt,istep,dx,&
                              prob_lo,prob_hi,need_regrid,time,ufile)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: flw(:),flw_old(:),flw_new(:),surf(:)
    type(multifab) , intent(inout) :: phi_new(:)
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(pf_para),   intent(inout) :: pfpara
    Real(kind=dp_t), intent(in   ) :: dt,time
    integer,         intent(in   ) :: istep
    real(kind=dp_t), intent(in   ) :: dx(:), prob_lo(mla%dim),prob_hi(mla%dim)
    logical,         intent(in   ) :: need_regrid
    integer,intent(in  ), optional :: ufile

    if( (pfpara%do_solve_flow == 0) ) return;

    select case(pfpara%flw_calculate_mode)
    ! case 0 is designed for phase field coupled LBM, which proceeds from finer grid to coarser ones
    case (0)
      if(present(ufile)) then
        call advance_flow_field_tb(mla,flw,phi_new,surf,pfpara,the_bc_tower,dt,istep,dx,&
                                   prob_lo,prob_hi,need_regrid,time,ufile)
      else
        call advance_flow_field_tb(mla,flw,phi_new,surf,pfpara,the_bc_tower,dt,istep,dx,&
                                   prob_lo,prob_hi,need_regrid,time)
      end if
    ! case 1 is designed for general case, which proceeds from coarser grid to finer ones
    case (1)
      if(present(ufile)) then
        call advance_flow_field_bt(mla,flw,flw_old,flw_new,surf,phi_new,pfpara,the_bc_tower,dt,istep,dx,&
                                   prob_lo,prob_hi,need_regrid,time,ufile)
      else
        call advance_flow_field_bt(mla,flw,flw_old,flw_new,surf,phi_new,pfpara,the_bc_tower,dt,istep,dx,&
                                   prob_lo,prob_hi,need_regrid,time)
      end if
    end select

end subroutine advance_flow_field

subroutine advance_flow_field_bt(mla,flw,flw_old,flw_new,surf,phi_new,pfpara,the_bc_tower,dt,istep,dx,&
                                 prob_lo,prob_hi,need_regrid,time,ufile)
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: flw(:),flw_old(:),flw_new(:)
    type(multifab) , intent(inout) :: phi_new(:),surf(:)
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(pf_para),   intent(inout) :: pfpara
    Real(kind=dp_t), intent(in   ) :: dt,time
    integer,         intent(in   ) :: istep
    real(kind=dp_t), intent(in   ) :: dx(:), prob_lo(mla%dim),prob_hi(mla%dim)
    logical,         intent(in   ) :: need_regrid
    integer,intent(in  ), optional :: ufile

    ! local variables
    integer          :: n,nlevels,n_max_levels,loops,n_sub_loops,nc_p,ng_p
    integer          :: dm,icom_start,ncomponent
    double precision :: alpha
    logical          :: need_cal_solid_flag

    !  time related routine for speed up
    !integer, parameter :: num_time = 5
    !integer :: ii
    !double precision :: cal_time_sub(mla%nlevel,num_time), start_sub(mla%nlevel,num_time)
    !double precision :: cal_time_flw(mla%nlevel,4)

    need_cal_solid_flag = .true.
    n_max_levels = mla%nlevel
    nlevels = mla%nlevel - pfpara%flw_skip_level_n
    dm = mla%dim 
    ng_p = flw(1)%ng
    nc_p = flw(1)%nc

    select case(dm)
    case (2)
      icom_start = 1
      ncomponent = 12
    case (3)
      icom_start = 1
      ncomponent = 23
    end select

    do n=1,nlevels
      call multifab_copy_c(flw_old(n),icom_start, flw(n),icom_start,ncomponent, ng_p)
      !cal_time_sub(n,:) = 0.d0
      !cal_time_flw(n,:) = 0.d0
    end do

    need_cal_solid_flag = need_regrid .or. (  pfpara%solve_flow_only /= 1) .or. &
                          ( (pfpara%solve_flow_only == 1)  .and. (istep == 1) ) 
    
    if(need_cal_solid_flag) then
      do n=1,n_max_levels
         !call cal_interface_flag_on_level(phi_new(n),pfpara,3)
         call cal_interface_flag_on_level(surf(n),phi_new(n),pfpara,1)
      end do
      !call ml_restrict_and_fill(n_max_levels, surf, mla%mba%rr, the_bc_tower%bc_tower_array,1,1,1)
    end if

    ! advance the first level, which spans the whole domain
    !start_sub(1,1) = parallel_wtime()

    !start_sub(1,4) = parallel_wtime()
    call advance_flow_field_on_level(1,flw(1),phi_new(1),surf(1),pfpara,&
                           the_bc_tower,dt,dx(1),prob_lo, prob_hi,time)!,cal_time_flw(1,:))
    !cal_time_sub(1,4) = cal_time_sub(1,4) + parallel_wtime() - start_sub(1,4)

    !start_sub(1,3) = parallel_wtime()
    !call multifab_fill_boundary_c(flw(1),icom_start ,ncomponent, ng_p)
    !cal_time_sub(1,3) = cal_time_sub(1,3) + parallel_wtime() - start_sub(1,3)

    !cal_time_sub(1,1) = cal_time_sub(1,1) + parallel_wtime() - start_sub(1,1)

    ! now climb upwardly for different levels
    do n = 2,nlevels
       !start_sub(n,1) = parallel_wtime()

       ! The subcycle number equals to 2 to the n-1
       ! For each subcycle, the inner boundary condition, i.e. ghost cells on the fine level
       ! must be updated using the intepolation value from the coarse grid level at the specific time
       n_sub_loops = 2**(n-1)

       do loops = 1,n_sub_loops
           !start_sub(n,2) = parallel_wtime()

           ! we do not need to do this for the last step of updating
           !if(loops .lt. n_sub_loops) then
              
           ! calculate the flw_old values using flw, i.e. flw_old = flw_old + (flw-flw_old)*loops/n_sub_loops
           call transit_multifab_flw_c_f(alpha,flw_new(n-1),flw_old(n-1),flw(n-1),surf(n-1),n-1,loops,n_sub_loops,&
                                         pfpara%tau_LBM,icom_start,ncomponent)

           !cal_time_sub(n,2) = cal_time_sub(n,2) + parallel_wtime() - start_sub(n,2)
           !start_sub(n,3) = parallel_wtime()
                
           ! calculate the flw_old values using flw, i.e. flw_old = flw_old + (flw-flw_old)*loops/n_sub_loops
           call multifab_fill_ghost_cells_flw(flw(n),flw_new(n-1),alpha,flw(n)%ng,mla%mba%rr(n-1,:), &
                                      the_bc_tower%bc_tower_array(n-1), &
                                      the_bc_tower%bc_tower_array(n), &
                                      icom_start,1,ncomponent)

           !call multifab_physbc(flw(n),icom_start,1,ncomponent,the_bc_tower%bc_tower_array(n))

           !cal_time_sub(n,3) = cal_time_sub(n,3) + parallel_wtime() - start_sub(n,3)
           !start_sub(n,4) = parallel_wtime()
           !end if

           ! first call this to update for the first time, the boundary is from time ZERO
           call advance_flow_field_on_level(n,flw(n),phi_new(n),surf(n),pfpara,the_bc_tower,dt,&
                                            dx(n),prob_lo, prob_hi,time)
           !cal_time_sub(n,4) = cal_time_sub(n,4) + parallel_wtime() - start_sub(n,4)

       end do
       !cal_time_sub(n,1) = cal_time_sub(n,1) + parallel_wtime() - start_sub(n,1)
    end do

    do n=nlevels,2,-1
       !start_sub(n,1) = parallel_wtime()
       !start_sub(n,5) = parallel_wtime()

       ! calculate the alpha
       call cal_tau_alpha(alpha, n-1,pfpara%tau_LBM,1)

       ! Do the same for orientation field
       call ml_cc_restriction_c_flw(flw(n-1),icom_start,flw(n),icom_start,&
                                    mla%mba%rr(n-1,:),ncomponent,1.d0/alpha)

       ! calculate the flw_old values using flw, i.e. flw_old = flw_old + (flw-flw_old)*loops/n_sub_loops
       ! this is a linear intepolation, might as well using others high order ones
       !call multifab_fill_ghost_cells_flw(flw(n),flw(n-1),alpha,flw(n)%ng,mla%mba%rr(n-1,:), &
       !                                  the_bc_tower%bc_tower_array(n-1), &
       !                                  the_bc_tower%bc_tower_array(n), &
       !                                  icom_start,1,ncomponent)
       !cal_time_sub(n,5) = cal_time_sub(n,5) + parallel_wtime() - start_sub(n,5)
       !cal_time_sub(n,1) = cal_time_sub(n,1) + parallel_wtime() - start_sub(n,1)

    end do

    ! the following action is for velocity and density, since they are all conserved variables
    ! we use normal routines to do the restriction and intepolatio
    ! pass back all the updated values
    !call ml_restrict_and_fill(nlevels, flw, mla%mba%rr, the_bc_tower%bc_tower_array,icom_start,1,ncomponent)

    ! do the intepolation stuff, possibly for saving computing time
    if(pfpara%flw_skip_level_n .gt. 0) then

       do n=nlevels,n_max_levels-1

          call fillpatch_flw(flw(n+1),flw(n),flw(n+1)%ng,mla%mba%rr(n,:), &
                                  the_bc_tower%bc_tower_array(n), &
                                  the_bc_tower%bc_tower_array(n+1), &
                                  1,1,1,dm+1,n,pfpara%tau_LBM)

          call check_solid_velocity_on_level(flw(n+1),surf(n+1))

       end do

    end if

    !if(parallel_IOProcessor() .and. present(ufile) .and. (mod(istep,20) .eq. 0) ) then
    !   write(ufile, '("---------------------------------------------------------------------------------")')
       
    !   do n=1,nlevels
    !      write(ufile, '("Grid Level = :", i8.7, i8.7, es10.2e3, es10.2e3, es10.2e3, es10.2e3, es10.2e3)') &
    !                    istep, n, cal_time_sub(n,1), cal_time_sub(n,2), cal_time_sub(n,3),cal_time_sub(n,4), &
    !                    cal_time_sub(n,5)
    !      write(ufile, '("In Flow    = :", es10.2e3, es10.2e3, es10.2e3, es10.2e3)') &
    !                    cal_time_flw(n,1),cal_time_flw(n,2),cal_time_flw(n,3),cal_time_flw(n,4)
    !   end do
    !   write(ufile, '("--------------------------------------------------------------------------------")')
    !end if

end subroutine advance_flow_field_bt

subroutine advance_flow_field_tb(mla,flw,phi_new,surf,pfpara,the_bc_tower,dt,istep,dx,&
                                 prob_lo,prob_hi,need_regrid,time,ufile)
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: flw(:)
    type(multifab) , intent(inout) :: phi_new(:),surf(:)
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(pf_para),   intent(inout) :: pfpara
    Real(kind=dp_t), intent(in   ) :: dt,time
    integer,         intent(in   ) :: istep
    real(kind=dp_t), intent(in   ) :: dx(:), prob_lo(mla%dim),prob_hi(mla%dim)
    logical,         intent(in   ) :: need_regrid
    integer,intent(in  ), optional :: ufile

    ! local variables
    integer :: n,nlevels,n_max_levels,loops,n_sub_loops,nc_p,ng_p
    integer :: dm,icom_start,ncomponent

    double precision :: alpha
    logical          :: need_cal_solid_flag

    !  to extract calculation time to Trackinfo.dat
    !double precision :: cal_time_sub(mla%nlevel,5), start_sub(mla%nlevel,5)
    !double precision :: cal_time_flw(mla%nlevel,4)

    need_cal_solid_flag = .true.

    ! The maximum number of grid levels
    n_max_levels = mla%nlevel

    ! The levels that requires calculation
    nlevels = mla%nlevel - pfpara%flw_skip_level_n

    dm = mla%dim 
    ng_p = flw(1)%ng
    nc_p = flw(1)%nc

    select case(dm)
    case (2)
      icom_start = 1 
      ncomponent = 12
    case (3)
      icom_start = 1 
      ncomponent = 23
    end select

    !do n=1,nlevels
    !  cal_time_sub(n,:) = 0.d0
    !  cal_time_flw(n,:) = 0.d0
    !end do

    ! If the phase field is calculated, this will be called for each time step
    ! If only the flow is calculated, this will be called for only the first time step
    need_cal_solid_flag = need_regrid .or. (  pfpara%solve_flow_only /= 1) .or. &
                          ( (pfpara%solve_flow_only == 1)  .and. (istep == 1) ) 
    
    ! Calculate the surface flag sf, 0 for interface, 1 for solid, -1 for liquid
    if(need_cal_solid_flag) then
      do n=1,n_max_levels
         ! calculate the interface flag, we use the temperature data to store this, 
         ! However, it will be replaced by the temperature old values each time the function
         ! is called
         !call cal_interface_flag_on_level(phi_new(n),pfpara,3)
         call cal_interface_flag_on_level(surf(n),phi_new(n),pfpara,1)
      end do
    end if

    ! advance the top level for each time step
    ! We do not need to call multifab_physbc because for flow, the boundary condition at domain 
    ! boundaries is taken care of explicitly
    !start_sub(nlevels,1) = parallel_wtime()

    ! calculate the flow
    !start_sub(nlevels,4) = parallel_wtime()
    call advance_flow_field_on_level(nlevels,flw(nlevels),phi_new(nlevels),surf(nlevels),pfpara,&
                                     the_bc_tower,dt,dx(nlevels),prob_lo, prob_hi,time)!,cal_time_flw(nlevels,:))
    !cal_time_sub(nlevels,4) = cal_time_sub(nlevels,4) + parallel_wtime() - start_sub(nlevels,4)

    ! exchange the ghost cells on the top level
    !start_sub(nlevels,3) = parallel_wtime()
    !call multifab_fill_boundary_c(flw(nlevels),icom_start ,ncomponent, ng_p)
    !cal_time_sub(nlevels,3) = cal_time_sub(nlevels,3) + parallel_wtime() - start_sub(nlevels,3)

    !cal_time_sub(nlevels,1) = cal_time_sub(nlevels,1) + parallel_wtime() - start_sub(nlevels,1)

    ! now climb downwardly for different levels
    do n = nlevels-1,1,-1

       !start_sub(n,1) = parallel_wtime()

       ! The subcycle number equals to 2 to the N-n
       n_sub_loops = 2**(n_max_levels - n)

       if( mod(istep, n_sub_loops) == 0 ) then
           
           !start_sub(n,4) = parallel_wtime()

           ! when the time step meets the requirement, do the flow calculation
           call advance_flow_field_on_level(n,flw(n),phi_new(n),surf(n),pfpara,the_bc_tower,dt,&
                                            dx(n),prob_lo, prob_hi,time)!,cal_time_flw(n,:))

           !cal_time_sub(n,4) = cal_time_sub(n,4) + parallel_wtime() - start_sub(n,4)

       end if

       !cal_time_sub(n,1) = cal_time_sub(n,1) + parallel_wtime() - start_sub(n,1)

    end do

    do n=nlevels,2,-1
       !start_sub(n,1) = parallel_wtime()

       ! subcycle number equals to 2 to the N-(n-1)
       ! synchronization
       n_sub_loops = 2**(n_max_levels - (n-1) )

       if( mod(istep, n_sub_loops) == 0 ) then

          !start_sub(n,5) = parallel_wtime()
          
          ! calculate the alpha
          ! this routine is in fillpatch.f90
          call cal_tau_alpha(alpha, n-1,pfpara%tau_LBM,1)

          ! Restrict the fine grid values to the coarse ones using special routines
          ! f(n-1) = 1/alfa*feq(n) + (1 - 1/alfa)*f(n)
          call ml_cc_restriction_c_flw(flw(n-1),icom_start,flw(n),icom_start,&
                                    mla%mba%rr(n-1,:),ncomponent,1.d0/alpha)
          !cal_time_sub(n,5) = cal_time_sub(n,5) + parallel_wtime() - start_sub(n,5)

          !start_sub(n,3) = parallel_wtime()
     
          ! Interpolate the coarse grid values to the fine ones using special routines
          ! f(n) = alfa*feq(n) + (1 - alfa)*f(n)
          call multifab_fill_ghost_cells_flw(flw(n),flw(n-1),alpha,flw(n)%ng,mla%mba%rr(n-1,:), &
                                      the_bc_tower%bc_tower_array(n-1), &
                                      the_bc_tower%bc_tower_array(n), &
                                      icom_start,1,ncomponent)

          !cal_time_sub(n,3) = cal_time_sub(n,3) + parallel_wtime() - start_sub(n,3)

       end if

       !cal_time_sub(n,1) = cal_time_sub(n,1) + parallel_wtime() - start_sub(n,1)

    end do

    ! do the intepolation stuff, possibly for saving computing time
    if(pfpara%flw_skip_level_n .gt. 0) then

       do n=nlevels,n_max_levels-1

          call fillpatch_flw(flw(n+1),flw(n),flw(n+1)%ng,mla%mba%rr(n,:), &
                                  the_bc_tower%bc_tower_array(n), &
                                  the_bc_tower%bc_tower_array(n+1), &
                                  1,1,1,dm+1,n,pfpara%tau_LBM)

          call check_solid_velocity_on_level(flw(n+1),surf(n+1))

       end do

    end if

    !if(parallel_IOProcessor() .and. present(ufile) .and. (mod(istep,20) .eq. 0) ) then
    !   write(ufile, '("---------------------------------------------------------------------------------")')
       
    !   do n=1,nlevels
    !      write(ufile, '("Grid Level = :", i8.7, i8.7, es10.2e3, es10.2e3, es10.2e3, es10.2e3, es10.2e3)') &
    !                    istep, n, cal_time_sub(n,1), cal_time_sub(n,2), cal_time_sub(n,3),cal_time_sub(n,4), &
    !                    cal_time_sub(n,5)
    !      write(ufile, '("In Flow    = :", es10.2e3, es10.2e3, es10.2e3, es10.2e3)') &
    !                    cal_time_flw(n,1),cal_time_flw(n,2),cal_time_flw(n,3),cal_time_flw(n,4)
    !   end do
    !   write(ufile, '("--------------------------------------------------------------------------------")')
    !end if

end subroutine advance_flow_field_tb

subroutine transit_multifab_flw_c_f0(alpha,n,flw_old,flw,flw_new,loops,n_sub_loops,tau_0,icom_start,ncomponent)
   Real(kind=dp_t), intent(  out) :: alpha
   type(multifab) , intent(inout) :: flw_new
   type(multifab) , intent(in   ) :: flw_old
   type(multifab) , intent(in   ) :: flw
   integer,         intent(in   ) :: loops, n_sub_loops,n,icom_start,ncomponent
   Real(kind=dp_t), intent(in   ) :: tau_0

   !type(multifab)     :: flw_equ

   integer ::  ng_p,nc
   double precision ::  var_inter = 0.d0,tau_c,tau_f
   ng_p = flw%ng
   nc = flw%nc

   var_inter = dble(loops-1) / dble(n_sub_loops)

   call multifab_copy_c(flw_new,icom_start, flw_old,icom_start,ncomponent, ng_p)

   if(loops .gt. 1) then
      call multifab_mult_mult_s_c(flw_new, icom_start, (1.d0/var_inter-1.d0),ncomponent,ng_p)
      call multifab_plus_plus_c(flw_new, icom_start, flw, icom_start, ncomponent,ng_p)
      call multifab_mult_mult_s_c(flw_new, icom_start, var_inter,ncomponent,ng_p)
   end if

   call cal_tau_alpha(alpha, n,tau_0,1)

end subroutine transit_multifab_flw_c_f0

subroutine transit_multifab_flw_c_f(alpha,flw_new,flw_old,flw,phi,n,loops,n_sub_loops,tau_0,icom_start,ncomponent)
   Real(kind=dp_t), intent(  out) :: alpha
   type(multifab) , intent(inout) :: flw_new
   type(multifab) , intent(in   ) :: flw_old,flw,phi
   integer,         intent(in   ) :: loops, n_sub_loops,n,icom_start,ncomponent
   Real(kind=dp_t), intent(in   ) :: tau_0

   !type(multifab)     :: flw_equ

   integer ::  ng_p,i,dm
   double precision ::  var_inter = 0.d0
   integer :: lo(flw%dim), hi(flw%dim)

   real(kind=dp_t), pointer :: flo(:,:,:,:)
   real(kind=dp_t), pointer :: fln(:,:,:,:)
   real(kind=dp_t), pointer :: fl(:,:,:,:)
   real(kind=dp_t), pointer :: pf(:,:,:,:)

   ng_p = flw%ng
   dm   = flw%dim

   var_inter = dble(loops-1) / dble(n_sub_loops)

   do i=1,nfabs(flw)
       flo => dataptr(flw_old,i,icom_start,ncomponent)
       fln => dataptr(flw_new,i,icom_start,ncomponent)
       fl  => dataptr(flw,i,icom_start,ncomponent)
       pf  => dataptr(phi,i)
       
       lo = lwb(get_box(flw,i))
       hi = upb(get_box(flw,i))

       select case(dm)
       case (2)
         call transit_multifab_flw_2d(fln(:,:,1,:),flo(:,:,1,:),fl(:,:,1,:),pf(:,:,1,1),ng_p,lo,hi,var_inter)

       case (3)
         call transit_multifab_flw_3d(fln(:,:,:,:),flo(:,:,:,:),fl(:,:,:,:),pf(:,:,:,1),ng_p,lo,hi,var_inter)

       end select

    end do

   call cal_tau_alpha(alpha, n,tau_0,1)

end subroutine transit_multifab_flw_c_f

subroutine transit_multifab_flw_2d(flw_new,flw_old,flw,sf,ng_p,lo,hi,var_inter)
  integer           :: ng_p,lo(2),hi(2)
  double precision  :: flw_new(lo(1)-ng_p:,lo(2)-ng_p:,1:)
  double precision  :: flw_old(lo(1)-ng_p:,lo(2)-ng_p:,1:)
  double precision  :: flw(lo(1)-ng_p:,lo(2)-ng_p:,1:)
  double precision  :: sf(lo(1)-ng_p:,lo(2)-ng_p:)
  real(kind=dp_t), intent(in   ) :: var_inter

  integer :: i,j
  real(dp_t) :: var_temp

  var_temp = 1.d0 - var_inter

  do j=lo(2)-ng_p,hi(2)+ng_p
    do i=lo(1)-ng_p,hi(1)+ng_p

      !if(sf(i,j) <= -1.d0) then

         flw_new(i,j,:) = flw_old(i,j,:)*var_temp + var_inter * flw(i,j,:)

      !end if

    end do
  end do


end subroutine transit_multifab_flw_2d

subroutine transit_multifab_flw_3d(flw_new,flw_old,flw,sf,ng_p,lo,hi,var_inter)
  integer           :: ng_p,lo(3),hi(3)
  double precision  :: flw_new(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:,1:)
  double precision  :: flw_old(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:,1:)
  double precision  :: flw(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:,1:)
  double precision  :: sf(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
  real(kind=dp_t), intent(in   ) :: var_inter

  integer :: i,j,k
  real(dp_t) :: var_temp

  var_temp = 1.d0 - var_inter

  do k=lo(3)-ng_p,hi(3)+ng_p
    do j=lo(2)-ng_p,hi(2)+ng_p
      do i=lo(1)-ng_p,hi(1)+ng_p

        !if(sf(i,j,k) <= -1.d0) then
          flw_new(i,j,k,:) = flw_old(i,j,k,:)*var_temp + var_inter * flw(i,j,k,:)
        !end if

      end do
    end do
  end do

end subroutine transit_multifab_flw_3d

subroutine advance_flow_field_on_level(lvl,flw,phi_new,surf,pfpara,the_bc_tower,dt,dx,prob_lo,prob_hi,time)
    integer,         intent(in   ) :: lvl
    type(multifab) , intent(inout) :: flw
    type(multifab) , intent(in   ) :: phi_new,surf
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(pf_para),   intent(inout) :: pfpara
    Real(kind=dp_t), intent(in   ) :: dt,time
    real(kind=dp_t), intent(in   ) :: dx, prob_lo(flw%dim),prob_hi(flw%dim)

    ! put f and equilibrium values inside the vectors
    call advance_k_series_on_level(lvl,flw,the_bc_tower,pfpara,phi_new,surf,time)

    ! advance velocity and f series
    call advance_velocityandf_series_on_level(lvl,flw,phi_new,surf,pfpara,the_bc_tower)

end subroutine advance_flow_field_on_level

subroutine advance_k_series_on_level(n,flw,the_bc_tower,pfpara,phi_new,surf,time,c_case)
    integer,         intent(in   ) :: n
    type(multifab) , intent(inout) :: flw
    type(multifab) , intent(in   ) :: phi_new,surf
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(pf_para),   intent(in   ) :: pfpara
    real(dp_t),      intent(in   ) :: time
    integer, intent(in), optional  :: c_case

    ! local variables
    integer :: lo(flw%dim), hi(flw%dim)
    integer :: dm, ng_p, i,d_iter,l_case
    integer :: iter_st(3),iter_end(3)

    real(kind=dp_t), pointer :: u_vec(:,:,:,:)
    real(kind=dp_t), pointer :: f_mat(:,:,:,:)
    real(kind=dp_t), pointer :: k_mat(:,:,:,:)
    real(kind=dp_t), pointer :: pf_new(:,:,:,:)
    real(kind=dp_t), pointer :: sf(:,:,:,:)

    dm   = flw%dim
    ng_p = flw%ng
    l_case = 2; if(present(c_case)) l_case = c_case

    select case(dm)
    case (2)
      iter_st  = (/1,4,13/)
      iter_end = (/3,9,9/)
    case (3)
      iter_st  = (/1,5,24/)
      iter_end = (/4,19,19/)
    end select

    do i=1,nfabs(flw)
       u_vec => dataptr(flw,i,iter_st(1),iter_end(1))
       f_mat => dataptr(flw,i,iter_st(2),iter_end(2))
       k_mat => dataptr(flw,i,iter_st(3),iter_end(3))
       pf_new => dataptr(phi_new,i)
       sf    => dataptr(surf,i)
       lo = lwb(get_box(flw,i))
       hi = upb(get_box(flw,i))

       select case(dm)
       case (2)
         call cal_kseries_2d(k_mat(:,:,1,:),u_vec(:,:,1,:),f_mat(:,:,1,:),ng_p,lo,hi,pfpara%tau_r(n),&
                             pfpara%gacc,pfpara%belta_c,pf_new(:,:,1,1),pf_new(:,:,1,2),sf(:,:,1,1),&
                             pfpara%sk,pfpara%dt_flw,pfpara%dim_lbm,pfpara%lbm_e(1,:),pfpara%lbm_e(2,:),&
                             pfpara%lbm_w(:),pfpara%kill_body_force,pfpara%lbm_drag_hc,pfpara%lbm_viscosity,&
                             pfpara%lbm_ratio_dx_W0,l_case,time,pfpara%lorentz_omega,&
                             pfpara%kill_lorentz_force,pfpara%lorentz_amp_x,pfpara%lorentz_amp_y)

       case (3)
         call cal_kseries_3d(k_mat(:,:,:,:),u_vec(:,:,:,:),f_mat(:,:,:,:),ng_p,lo,hi,pfpara%tau_r(n),& 
                             pfpara%gacc,pfpara%belta_c,pf_new(:,:,:,1),pf_new(:,:,:,2),sf(:,:,:,1),&
                             pfpara%sk,pfpara%dt_flw,pfpara%dim_lbm,pfpara%lbm_e(1,:),pfpara%lbm_e(2,:),&
                             pfpara%lbm_e(3,:),pfpara%lbm_w(:),pfpara%kill_body_force,pfpara%lbm_drag_hc,&
                             pfpara%lbm_viscosity,pfpara%lbm_ratio_dx_W0,l_case,time,&
                             pfpara%lorentz_omega,pfpara%kill_lorentz_force,pfpara%lorentz_amp_x,&
                             pfpara%lorentz_amp_y,pfpara%lorentz_amp_z)

       end select

    end do

  end subroutine advance_k_series_on_level


! This subroutine updates both velocity and the f_mat variables
subroutine advance_velocityandf_series_on_level(n,flw,phi_new,surf,pfpara,the_bc_tower,c_case)
    integer,         intent(in   ) :: n
    type(multifab) , intent(in   ) :: phi_new,surf
    type(multifab) , intent(inout) :: flw
    type(pf_para),   intent(in   ) :: pfpara
    type(bc_tower) , intent(in   ) :: the_bc_tower
    integer, intent(in), optional  :: c_case

    ! local variables
    integer :: lo(flw%dim), hi(flw%dim),sw_flag(flw%dim,2)
    integer :: dm, ng_p, i,d_iter,l_case
    integer :: iter_st(3),iter_end(3),icom_start,ncomponent
    type(box) :: pd

    real(kind=dp_t), pointer :: u_vec(:,:,:,:)
    real(kind=dp_t), pointer :: f_mat(:,:,:,:)
    real(kind=dp_t), pointer :: k_mat(:,:,:,:)
    real(kind=dp_t), pointer :: pfn(:,:,:,:)

    dm   = flw%dim
    ng_p = flw%ng
    l_case = 1; if(present(c_case)) l_case = c_case

    select case(dm)
    case (2)
      icom_start = 1
      ncomponent = 12
      iter_st  = (/1,4,13/)
      iter_end = (/3,9,9/)
    case (3)
      icom_start = 1
      ncomponent = 23
      iter_st  = (/1,5,24/)
      iter_end = (/4,19,19/)
    end select

    ! calculate the RK matrices
    pd = layout_get_pd(flw%la)

    do i=1,nfabs(flw)
       u_vec => dataptr(flw,i,iter_st(1),iter_end(1))
       f_mat => dataptr(flw,i,iter_st(2),iter_end(2))
       k_mat => dataptr(flw,i,iter_st(3),iter_end(3))
       pfn  => dataptr(surf,i)
       lo = lwb(get_box(flw,i))
       hi = upb(get_box(flw,i))

       sw_flag(:,:) = 0
       do d_iter=1,dm
         if( lo(d_iter) == lwb(pd,d_iter) ) sw_flag(d_iter,1) = 1
         if( hi(d_iter) == upb(pd,d_iter) ) sw_flag(d_iter,2) = 1
       end do

       select case(dm)
       case (2)
          call cal_vandf_2d( u_vec(:,:,1,:),f_mat(:,:,1,:),k_mat(:,:,1,:),& 
                             pfn(:,:,1,1),ng_p, lo,hi,pfpara%dim_lbm,pfpara%lbm_e(1,:),&
                             pfpara%lbm_e(2,:),pfpara%lbm_w(:))
          
          call update_boundary_2d_X(u_vec(:,:,1,:),f_mat(:,:,1,:),pfn(:,:,1,1),ng_p,lo,hi,&
                              pfpara%u_inlet(1,1),pfpara%u_inlet(1,2),sw_flag,pfpara%bc_flag(1,:),&
								              pfpara%cal_boundary_spec,pfpara%dim_lbm,pfpara%lbm_e(1,:),&
                              pfpara%lbm_e(2,:),pfpara%lbm_w(:))
          
          call update_boundary_2d_Y(u_vec(:,:,1,:),f_mat(:,:,1,:),pfn(:,:,1,1),ng_p,lo,hi,&
                              pfpara%u_inlet(2,1),pfpara%u_inlet(2,2),sw_flag,pfpara%bc_flag(2,:),&
								              pfpara%cal_boundary_spec,pfpara%dim_lbm,pfpara%lbm_e(1,:),&
                              pfpara%lbm_e(2,:),pfpara%lbm_w(:))
          

       case (3)
          call cal_vandf_3d( u_vec(:,:,:,:),f_mat(:,:,:,:),k_mat(:,:,:,:),& 
                             pfn(:,:,:,1),ng_p, lo,hi,pfpara%dim_lbm,pfpara%lbm_e(1,:),&
                             pfpara%lbm_e(2,:),pfpara%lbm_e(3,:),pfpara%lbm_w(:) )

          call update_boundary_3d_X(u_vec(:,:,:,:),f_mat(:,:,:,:),pfn(:,:,:,1),ng_p,lo,hi,&
                              pfpara%u_inlet(1,1),pfpara%u_inlet(1,2),sw_flag,pfpara%bc_flag(1,:),&
								              pfpara%cal_boundary_spec,pfpara%dim_lbm,pfpara%lbm_e(1,:),&
                              pfpara%lbm_e(2,:),pfpara%lbm_e(3,:),pfpara%lbm_w(:))
          
          call update_boundary_3d_Y(u_vec(:,:,:,:),f_mat(:,:,:,:),pfn(:,:,:,1),ng_p,lo,hi,&
                              pfpara%u_inlet(2,1),pfpara%u_inlet(2,2),sw_flag,pfpara%bc_flag(2,:),&
								              pfpara%cal_boundary_spec,pfpara%dim_lbm,pfpara%lbm_e(1,:),&
                              pfpara%lbm_e(2,:),pfpara%lbm_e(3,:),pfpara%lbm_w(:))
          
          call update_boundary_3d_Z(u_vec(:,:,:,:),f_mat(:,:,:,:),pfn(:,:,:,1),ng_p,lo,hi,&
                              pfpara%u_inlet(3,1),pfpara%u_inlet(3,2),sw_flag,pfpara%bc_flag(3,:),&
								              pfpara%cal_boundary_spec,pfpara%dim_lbm,pfpara%lbm_e(1,:),&
                              pfpara%lbm_e(2,:),pfpara%lbm_e(3,:),pfpara%lbm_w(:))

       end select
    end do

    call multifab_fill_boundary_c(flw, icom_start, ncomponent, ng_p)
    call multifab_physbc(flw,icom_start,1,ncomponent,the_bc_tower%bc_tower_array(n))
    
end subroutine advance_velocityandf_series_on_level

! Let's follow the trasition that we always put output variables at first
subroutine cal_kseries_2d(k_mat,u_vec,f_mat,ng_p,lo,hi,tau,gacc,belta_c,pf,uc,sf,k_p,dt_flw,q,&
                          ex,ey,w,kill_body_force,l_hc,l_viscosity,l_dx_w0,l_case,time,&
                          lorentz_omega,kill_lorentz_force,lorentz_amp_x,lorentz_amp_y)
  integer           :: lo(2),hi(2),ng_p,q
  integer           :: ex(q),ey(q),kill_body_force,l_case,kill_lorentz_force
  double precision  :: k_mat(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,1:9)
  double precision  :: f_mat(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,1:9)
  double precision  :: u_vec(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,1:3)
  double precision  :: pf(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: uc(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: sf(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: tau,gacc,belta_c,k_p,dt_flw,w(q),lorentz_amp_x,lorentz_amp_y
  double precision  :: l_hc, l_viscosity, l_dx_w0,time,lorentz_omega

  ! local variables
  integer           :: i,j, iter,ip(q),jp(q)
  double precision  :: var_1,var_2, var_sub, var_temp, var_sin
  double precision  :: F_drag, F_body, F_lorentz, v_sp_1,v_sp_2,fl

  var_1   = 0.d0; var_2   = 0.d0; var_sub = 0.d0
  F_body  = 0.d0; F_drag  = 0.d0; F_lorentz = 0.d0; var_sin = 0.d0
  v_sp_1  = l_hc *l_viscosity*l_dx_w0**2
  var_sin = sin(lorentz_omega*time)

 select case(l_case)
 case(1)
  ! do the loops
  do j=lo(2),hi(2)
    do i=lo(1),hi(1)
  
      do iter=1,q
        ip(iter)=i-ex(iter)
        jp(iter)=j-ey(iter)
      end do

      if(sf(i,j) .lt. 0.d0) then

        ! This loop is rather large, it updates the k_mat based on the value from both u_vec and f_mat
        do iter=1,q

          if( sf(ip(iter),jp(iter)) .lt. 0.d0 ) then
  
               var_1 = ex(iter)*u_vec(ip(iter),jp(iter),1)+ey(iter)*u_vec(ip(iter),jp(iter),2)
               var_2 = u_vec(ip(iter),jp(iter),1)**2+u_vec(ip(iter),jp(iter),2)**2

               !fl = 0.5d0 * ( 1.d0 - pf(ip(iter),jp(iter)) )

               !F_drag = -6.d0*v_sp_1*u_vec(ip(iter),jp(iter),3)*u_vec(ip(iter),jp(iter),3)*fl*(1.d0-fl)**2
               !F_drag = F_drag*w(iter)*( var_1 - var_2 +3.d0*var_1*var_1)

               if(kill_body_force == 0) then
                  fl = 0.5d0 * ( 1.d0 - pf(ip(iter),jp(iter)) )
                  var_sub = (1.d0+(1.d0-k_p)*uc(ip(iter),jp(iter)))*&
                            (1.d0+k_p-(1.d0-k_p)*pf(ip(iter),jp(iter))) -2.d0*k_p

                  F_body = 3.d0*u_vec(ip(iter),jp(iter),3)**2*fl
                  F_body = F_body*w(iter)*belta_c*var_sub*gacc
                  F_body = F_body*( ey(iter)-u_vec(ip(iter),jp(iter),2) +3.d0*var_1*ey(iter) )
               end if

               k_mat(i,j,iter) = f_mat(ip(iter),jp(iter),iter) + &
               ( w(iter)*u_vec(ip(iter),jp(iter),3)*(1.0+3.0*var_1+4.5*var_1**2-1.5*var_2)-&
               f_mat(ip(iter),jp(iter),iter) )/tau + F_body ! + F_drag
            
            else
               
               k_mat(i,j,iter) = f_mat(ip(iter),jp(iter),iter)
               
            end if

        end do

      else

        ! This loop is rather large, it updates the k_mat based on the value from both u_vec and f_mat
        do iter=1,q
          k_mat(i,j,iter) = f_mat(ip(iter),jp(iter),iter)
        end do

     end if

    end do
  end do

 case(3)

  ! do the loops
  do j=lo(2),hi(2)
    do i=lo(1),hi(1)

      if(sf(i,j) .lt. 0.d0) cycle;
      
      do iter=1,q
        ip(iter)=i-ex(iter)
        jp(iter)=j-ey(iter)
      end do

      ! This loop is rather large, it updates the k_mat based on the value from both u_vec and f_mat
      do iter=1,q
        k_mat(i,j,iter) = f_mat(ip(iter),jp(iter),iter)
      end do

      var_temp = k_mat(i,j,2); k_mat(i,j,2)=k_mat(i,j,4); k_mat(i,j,4)=var_temp
      var_temp = k_mat(i,j,3); k_mat(i,j,3)=k_mat(i,j,5); k_mat(i,j,5)=var_temp
      var_temp = k_mat(i,j,6); k_mat(i,j,6)=k_mat(i,j,8); k_mat(i,j,8)=var_temp
      var_temp = k_mat(i,j,7); k_mat(i,j,7)=k_mat(i,j,9); k_mat(i,j,9)=var_temp

    end do
  end do

 case(2)

  ! do the loops
  do j=lo(2),hi(2)
    do i=lo(1),hi(1)

      if(sf(i,j) .gt. 0.d0) cycle;
  
      do iter=1,q
        ip(iter)=i-ex(iter)
        jp(iter)=j-ey(iter)
      end do

      !if(sf(i,j) .lt. 0.d0) then
      ! This loop is rather large, it updates the k_mat based on the value from both u_vec and f_mat
      do iter=1,q

        if( sf(ip(iter),jp(iter)) .lt. 0.d0 ) then
  
            var_1 = ex(iter)*u_vec(ip(iter),jp(iter),1)+ey(iter)*u_vec(ip(iter),jp(iter),2)
            var_2 = u_vec(ip(iter),jp(iter),1)**2+u_vec(ip(iter),jp(iter),2)**2

            fl = 0.5d0 * ( 1.d0 - pf(ip(iter),jp(iter)) )   
            F_drag = -6.d0*v_sp_1*u_vec(ip(iter),jp(iter),3)**2*(1.d0-fl)**2
            F_drag = F_drag*w(iter)*( var_1 - var_2 +3.d0*var_1**2)

            if(kill_body_force == 0) then
               var_sub = (1.d0+(1.d0-k_p)*uc(ip(iter),jp(iter)))*&
                         (1.d0+k_p-(1.d0-k_p)*pf(ip(iter),jp(iter))) -2.d0*k_p

               F_body = 3.d0*u_vec(ip(iter),jp(iter),3)**2*fl
               F_body = F_body*w(iter)*belta_c*var_sub*gacc
               F_body = F_body*( ey(iter)-u_vec(ip(iter),jp(iter),2) + 3.d0* var_1*ey(iter) )
            end if

            if(kill_lorentz_force == 0) then

               F_lorentz = 3.d0*u_vec(ip(iter),jp(iter),3)**2*w(iter)*fl*var_sin
               F_lorentz = F_lorentz * &
                         ( lorentz_amp_x *( ex(iter)-u_vec(ip(iter),jp(iter),1) +3.d0*var_1*ex(iter) ) + &
                           lorentz_amp_y *( ey(iter)-u_vec(ip(iter),jp(iter),2) +3.d0*var_1*ey(iter) ) )

            end if

            k_mat(i,j,iter) = f_mat(ip(iter),jp(iter),iter) + &
               ( w(iter)*u_vec(ip(iter),jp(iter),3)*(1.0+3.0*var_1+4.5*var_1**2-1.5*var_2)-&
               f_mat(ip(iter),jp(iter),iter) )/tau + F_body + F_lorentz + F_drag
            
        else
               
            k_mat(i,j,iter) = f_mat(ip(iter),jp(iter),iter)
               
        end if

      end do

      !else
      ! This loop is rather large, it updates the k_mat based on the value from both u_vec and f_mat

      if(sf(i,j) .eq. 0.d0) then

        var_temp = k_mat(i,j,2); k_mat(i,j,2)=k_mat(i,j,4); k_mat(i,j,4)=var_temp
        var_temp = k_mat(i,j,3); k_mat(i,j,3)=k_mat(i,j,5); k_mat(i,j,5)=var_temp
        var_temp = k_mat(i,j,6); k_mat(i,j,6)=k_mat(i,j,8); k_mat(i,j,8)=var_temp
        var_temp = k_mat(i,j,7); k_mat(i,j,7)=k_mat(i,j,9); k_mat(i,j,9)=var_temp

      end if

     !end if

    end do
  end do

  end select

end subroutine cal_kseries_2d

subroutine cal_kseries_3d(k_mat,u_vec,f_mat,ng_p,lo,hi,tau,gacc,belta_c,pf,uc,sf,k_p,dt_flw,&
                          q,ex,ey,ez,w,kill_body_force,l_hc,l_viscosity,l_dx_w0,l_case,time,&
                          lorentz_omega,kill_lorentz_force,lorentz_amp_x,lorentz_amp_y,lorentz_amp_z)
  integer          :: lo(3),hi(3),ng_p,q
  integer          :: ex(q),ey(q),ez(q),kill_body_force,l_case,kill_lorentz_force
  double precision :: pf(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision :: uc(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision :: sf(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision :: u_vec(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p,1:4)
  double precision :: k_mat(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p,1:19)
  double precision :: f_mat(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p,1:19)
  double precision :: tau,gacc,belta_c,k_p,dt_flw,w(q),lorentz_amp_x,lorentz_amp_y,lorentz_amp_z
  double precision :: l_hc, l_viscosity, l_dx_w0,time,lorentz_omega

  ! local variables
  integer           :: i,j,k,iter,ip(q),jp(q),kp(q)
  double precision  :: var_temp_1,var_temp_2, var_sub, var_temp,var_sin
  double precision  :: F_drag, F_body, v_sp_1,v_sp_2,fl,F_lorentz

  var_temp_1 = 0.d0; var_temp_2 = 0.d0; var_sub = 0.d0
  F_body  = 0.d0; F_drag  = 0.d0; F_lorentz = 0.d0
  v_sp_1  = l_hc *l_viscosity*l_dx_w0*l_dx_w0
  var_sin = sin(lorentz_omega*time)

 select case(l_case)
 case(1)

  ! do the loops
  do k=lo(3),hi(3)
    do j=lo(2),hi(2)
      do i=lo(1),hi(1)

        do iter=1,q
          ip(iter) = i-ex(iter)
          jp(iter) = j-ey(iter)
          kp(iter) = k-ez(iter)
        end do
        
        ! only calculate the liquid inside, i.e. excluding the these on the domain surface
        if(sf(i,j,k) .lt. 0.d0) then

           do iter=1,q

             ! For liquid cells, both streaming and collision are performed
             if(sf(ip(iter),jp(iter),kp(iter)) .lt. -1.d0) then

                 var_temp_1 = ex(iter)*u_vec(ip(iter),jp(iter),kp(iter),1)+&
                           ey(iter)*u_vec(ip(iter),jp(iter),kp(iter),2)+&
                           ez(iter)*u_vec(ip(iter),jp(iter),kp(iter),3)

                 var_temp_2 = u_vec(ip(iter),jp(iter),kp(iter),1)**2+&
                           u_vec(ip(iter),jp(iter),kp(iter),2)**2+&
                           u_vec(ip(iter),jp(iter),kp(iter),3)**2

                 !fl = 0.5d0 * ( 1.d0 - pf(ip(iter),jp(iter),kp(iter)) )

                 !F_drag = -6.d0*v_sp_1*u_vec(ip(iter),jp(iter),kp(iter),4)**2*fl*(1.d0-fl)**2
                 !F_drag = F_drag*w(iter)*( var_1 - var_2 +3.d0*var_1*var_1)

                 fl = 0.5d0 * ( 1.d0 - pf(ip(iter),jp(iter),kp(iter)) )

                 if(kill_body_force == 0) then

                    var_sub = (1.d0+(1.d0-k_p)*uc(ip(iter),jp(iter),kp(iter)))*&
                              (1.d0+k_p-(1.d0-k_p)*pf(ip(iter),jp(iter),kp(iter))) -2.d0*k_p

                    F_body = 3.d0*u_vec(ip(iter),jp(iter),kp(iter),4)**2*fl
                    F_body = F_body*w(iter)*belta_c*var_sub*gacc
                    F_body = F_body*( ey(iter)-u_vec(ip(iter),jp(iter),kp(iter),2) +3.d0*var_temp_1*ey(iter) )
                 end if

                 if(kill_lorentz_force == 0) then

                   F_lorentz = 3.d0*u_vec(ip(iter),jp(iter),kp(iter),4)*w(iter)*fl*var_sin
                   F_lorentz = F_lorentz * &
                    ( lorentz_amp_x *( ex(iter)-u_vec(ip(iter),jp(iter),kp(iter),1) +3.d0*var_temp_1*ex(iter) ) + &
                      lorentz_amp_y *( ey(iter)-u_vec(ip(iter),jp(iter),kp(iter),2) +3.d0*var_temp_1*ey(iter) ) + &
                      lorentz_amp_z *( ez(iter)-u_vec(ip(iter),jp(iter),kp(iter),3) +3.d0*var_temp_1*ez(iter) ) )

                 end if

                 k_mat(i,j,k,iter) = f_mat(ip(iter),jp(iter),kp(iter),iter) +&
                 ( w(iter)*u_vec(ip(iter),jp(iter),kp(iter),4)*(1.0+3.0*var_temp_1+4.5*var_temp_1**2-1.5*var_temp_2) - &
                                f_mat(ip(iter),jp(iter),kp(iter),iter) )/tau + F_body + F_lorentz!+ F_drag

             else

                 k_mat(i,j,k,iter) = f_mat(ip(iter),jp(iter),kp(iter),iter)

             end if

           end do

        else 
          
          do iter=1,q
            k_mat(i,j,k,iter) = f_mat(ip(iter),jp(iter),kp(iter),iter)
          end do

        end if

      end do
    end do
  end do

 case(3)

  ! do the loops
  do k=lo(3),hi(3)
    do j=lo(2),hi(2)
      do i=lo(1),hi(1)

        if(sf(i,j,k) .lt. 0.d0) cycle;

        do iter=1,q
          ip(iter) = i-ex(iter)
          jp(iter) = j-ey(iter)
          kp(iter) = k-ez(iter)
        end do
        
        do iter=1,q
          k_mat(i,j,k,iter) = f_mat(ip(iter),jp(iter),kp(iter),iter)
        end do

        ! Bouncing back
        var_temp = k_mat(i,j,k,2); k_mat(i,j,k,2)=k_mat(i,j,k,3); k_mat(i,j,k,3)=var_temp
        var_temp = k_mat(i,j,k,4); k_mat(i,j,k,4)=k_mat(i,j,k,5); k_mat(i,j,k,5)=var_temp
        var_temp = k_mat(i,j,k,6); k_mat(i,j,k,6)=k_mat(i,j,k,7); k_mat(i,j,k,7)=var_temp

        var_temp = k_mat(i,j,k,8); k_mat(i,j,k,8)=k_mat(i,j,k,11); k_mat(i,j,k,11)=var_temp
        var_temp = k_mat(i,j,k,9); k_mat(i,j,k,9)=k_mat(i,j,k,10); k_mat(i,j,k,10)=var_temp

        var_temp = k_mat(i,j,k,12); k_mat(i,j,k,12)=k_mat(i,j,k,15); k_mat(i,j,k,15)=var_temp
        var_temp = k_mat(i,j,k,13); k_mat(i,j,k,13)=k_mat(i,j,k,14); k_mat(i,j,k,14)=var_temp

        var_temp = k_mat(i,j,k,16); k_mat(i,j,k,16)=k_mat(i,j,k,19); k_mat(i,j,k,19)=var_temp
        var_temp = k_mat(i,j,k,17); k_mat(i,j,k,17)=k_mat(i,j,k,18); k_mat(i,j,k,18)=var_temp

      end do
    end do
  end do

 case(2)

  ! do the loops
  do k=lo(3),hi(3)
    do j=lo(2),hi(2)
      do i=lo(1),hi(1)

        if(sf(i,j,k) .gt. 0.d0) cycle;

        do iter=1,q
          ip(iter) = i-ex(iter)
          jp(iter) = j-ey(iter)
          kp(iter) = k-ez(iter)
        end do
        
        ! only calculate the liquid inside, i.e. excluding the these on the domain surface
        !if(sf(i,j,k) .lt. 0.d0) then
        do iter=1,q

          ! For liquid cells, both streaming and collision are performed
          if(sf(ip(iter),jp(iter),kp(iter)) .lt. 0.d0) then

              var_temp_1 = ex(iter)*u_vec(ip(iter),jp(iter),kp(iter),1)+&
                           ey(iter)*u_vec(ip(iter),jp(iter),kp(iter),2)+&
                           ez(iter)*u_vec(ip(iter),jp(iter),kp(iter),3)

              var_temp_2 = u_vec(ip(iter),jp(iter),kp(iter),1)**2+&
                           u_vec(ip(iter),jp(iter),kp(iter),2)**2+&
                           u_vec(ip(iter),jp(iter),kp(iter),3)**2

              fl = 0.5d0 * ( 1.d0 - pf(ip(iter),jp(iter),kp(iter)) )

              F_drag = -6.d0*v_sp_1*u_vec(ip(iter),jp(iter),kp(iter),4)**2*fl*(1.d0-fl)**2
              F_drag = F_drag*w(iter)*( var_temp_1 - var_temp_2 +3.d0*var_temp_1**2)

              if(kill_body_force == 0) then

                 var_sub = (1.d0+(1.d0-k_p)*uc(ip(iter),jp(iter),kp(iter)))*&
                           (1.d0+k_p-(1.d0-k_p)*pf(ip(iter),jp(iter),kp(iter))) -2.d0*k_p

                 F_body = 3.d0*u_vec(ip(iter),jp(iter),kp(iter),4)**2*fl
                 F_body = F_body*w(iter)*belta_c*var_sub*gacc
                 F_body = F_body*( ey(iter)-u_vec(ip(iter),jp(iter),kp(iter),2) +3.d0*var_temp_1*ey(iter) )

              end if

              if(kill_lorentz_force == 0) then

                  F_lorentz = 3.d0*u_vec(ip(iter),jp(iter),kp(iter),4)**2*w(iter)*fl*var_sin
                  F_lorentz = F_lorentz * &
                    ( lorentz_amp_x *( ex(iter)-u_vec(ip(iter),jp(iter),kp(iter),1) +3.d0*var_temp_1*ex(iter) ) + &
                      lorentz_amp_y *( ey(iter)-u_vec(ip(iter),jp(iter),kp(iter),2) +3.d0*var_temp_1*ey(iter) ) + &
                      lorentz_amp_z *( ez(iter)-u_vec(ip(iter),jp(iter),kp(iter),3) +3.d0*var_temp_1*ez(iter) ) )

              end if

              k_mat(i,j,k,iter) = f_mat(ip(iter),jp(iter),kp(iter),iter) +&
                 ( w(iter)*u_vec(ip(iter),jp(iter),kp(iter),4)*(1.0+3.0*var_temp_1+4.5*var_temp_1**2-1.5*var_temp_2) - &
                                  f_mat(ip(iter),jp(iter),kp(iter),iter) )/tau + F_body + F_lorentz + F_drag

          else

              k_mat(i,j,k,iter) = f_mat(ip(iter),jp(iter),kp(iter),iter)

          end if

        end do

        !else 
          
        if(sf(i,j,k) .eq. 0.d0) then

          ! Bouncing back
          var_temp = k_mat(i,j,k,2); k_mat(i,j,k,2)=k_mat(i,j,k,3); k_mat(i,j,k,3)=var_temp
          var_temp = k_mat(i,j,k,4); k_mat(i,j,k,4)=k_mat(i,j,k,5); k_mat(i,j,k,5)=var_temp
          var_temp = k_mat(i,j,k,6); k_mat(i,j,k,6)=k_mat(i,j,k,7); k_mat(i,j,k,7)=var_temp

          var_temp = k_mat(i,j,k,8); k_mat(i,j,k,8)=k_mat(i,j,k,11); k_mat(i,j,k,11)=var_temp
          var_temp = k_mat(i,j,k,9); k_mat(i,j,k,9)=k_mat(i,j,k,10); k_mat(i,j,k,10)=var_temp

          var_temp = k_mat(i,j,k,12); k_mat(i,j,k,12)=k_mat(i,j,k,15); k_mat(i,j,k,15)=var_temp
          var_temp = k_mat(i,j,k,13); k_mat(i,j,k,13)=k_mat(i,j,k,14); k_mat(i,j,k,14)=var_temp

          var_temp = k_mat(i,j,k,16); k_mat(i,j,k,16)=k_mat(i,j,k,19); k_mat(i,j,k,19)=var_temp
          var_temp = k_mat(i,j,k,17); k_mat(i,j,k,17)=k_mat(i,j,k,18); k_mat(i,j,k,18)=var_temp

        end if

      end do
    end do
  end do

 end select

end subroutine cal_kseries_3d

subroutine cal_vandf_2d(u_vec,f_mat,k_mat,sf,ng_p,lo,hi,q,ex,ey,w)
  integer           :: lo(2),hi(2),ng_p,q 
  integer           :: ex(q),ey(q)
  double precision  :: u_vec(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,1:3)
  double precision  :: f_mat(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,1:9)
  double precision  :: k_mat(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,1:9)
  double precision  :: sf(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: w(q)

  ! local variables
  integer           :: i,j,iter

  ! do the loops
  do j=lo(2),hi(2)
    do i=lo(1),hi(1)

     if(sf(i,j) .ge. 0.d0) then
       u_vec(i,j,1:2)=0.d0
       cycle
     end if
     
     do iter=1,q
       f_mat(i,j,iter)=k_mat(i,j,iter)
     end do

       u_vec(i,j,:)=0.d0

       do iter=1,q
         u_vec(i,j,3)=u_vec(i,j,3)+f_mat(i,j,iter)
	       u_vec(i,j,1)=u_vec(i,j,1)+ex(iter)*f_mat(i,j,iter)
	       u_vec(i,j,2)=u_vec(i,j,2)+ey(iter)*f_mat(i,j,iter)
       end do

       u_vec(i,j,1)=u_vec(i,j,1)/u_vec(i,j,3)
       u_vec(i,j,2)=u_vec(i,j,2)/u_vec(i,j,3)

    end do
  end do

end subroutine cal_vandf_2d

subroutine cal_vandf_3d(u_vec,f_mat,k_mat,sf,ng_p,lo,hi,q,ex,ey,ez,w)
  integer           :: lo(3),hi(3),ng_p,q 
  integer           :: ex(q),ey(q),ez(q)
  double precision  :: sf(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: u_vec(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p,1:4)
  double precision  :: f_mat(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p,1:19)
  double precision  :: k_mat(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p,1:19)
  double precision  :: w(q)

  ! local variables
  integer           :: i,j,k,iter

  ! do the loops
  do k=lo(3),hi(3)
    do j=lo(2),hi(2)
      do i=lo(1),hi(1)

        if(sf(i,j,k) .ge. 0.d0) then
          u_vec(i,j,k,1:3)=0.d0
          cycle
        end if

        do iter=1,q
          f_mat(i,j,k,iter) = k_mat(i,j,k,iter)
        end do

        u_vec(i,j,k,:)=0.d0

        do iter=1,q
          u_vec(i,j,k,4) = u_vec(i,j,k,4) + f_mat(i,j,k,iter)
          u_vec(i,j,k,1) = u_vec(i,j,k,1) + ex(iter)*f_mat(i,j,k,iter)
          u_vec(i,j,k,2) = u_vec(i,j,k,2) + ey(iter)*f_mat(i,j,k,iter)
          u_vec(i,j,k,3) = u_vec(i,j,k,3) + ez(iter)*f_mat(i,j,k,iter)
        end do

	      u_vec(i,j,k,1)=u_vec(i,j,k,1)/u_vec(i,j,k,4)
	      u_vec(i,j,k,2)=u_vec(i,j,k,2)/u_vec(i,j,k,4)
	      u_vec(i,j,k,3)=u_vec(i,j,k,3)/u_vec(i,j,k,4)

      end do
    end do
  end do

end subroutine cal_vandf_3d

subroutine cal_interface_flag_on_level(surf,phi,pfpara,ncomp)
    type(multifab) , intent(inout) :: surf
    type(multifab) , intent(in   ) :: phi
    type(pf_para),   intent(in   ) :: pfpara
    integer,         intent(in   ) :: ncomp

    ! local variables
    integer :: lo(phi%dim), hi(phi%dim),sw_flag(phi%dim,2),d_iter
    integer :: dm, ng_p, i

    real(kind=dp_t), pointer ::  pfo(:,:,:,:)
    real(kind=dp_t), pointer ::  sf(:,:,:,:)
    type(box) :: pd

    dm = phi%dim
    ng_p = phi%ng
    pd = layout_get_pd(phi%la)
    
    do i=1,nfabs(phi)
          
       pfo  => dataptr(phi,i)
       sf   => dataptr(surf,i)
        
       lo = lwb(get_box(phi,i))
       hi = upb(get_box(phi,i))

       sw_flag(:,:) = 0
       do d_iter=1,dm
         if( lo(d_iter) == lwb(pd,d_iter) ) sw_flag(d_iter,1) = 1
         if( hi(d_iter) == upb(pd,d_iter) ) sw_flag(d_iter,2) = 1
       end do

       select case(dm)
          case (2)
          ! calculate the K1
             call cal_interface_flag_2d(sf(:,:,1,ncomp),pfo(:,:,1,1),ng_p,lo,hi,&
                                        pfpara%thers_solid,sw_flag)

          case (3)
             call cal_interface_flag_3d(sf(:,:,:,ncomp),pfo(:,:,:,1),ng_p,lo,hi,&
                                        pfpara%thers_solid,sw_flag)

       end select

    end do

    call multifab_fill_boundary_c(surf,ncomp ,1, ng_p)
    
end subroutine cal_interface_flag_on_level

subroutine cal_interface_flag_2d(sf,pf,ng_p,lo,hi,thers_solid,sw_flag)
  integer           :: lo(2),hi(2),ng_p,sw_flag(2,2)
  double precision  :: pf(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: sf(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: thers_solid

  integer           :: i,j,ii,jj
  double precision  :: solid_f, liquid_f, num_p
  logical           :: find_liquid

  num_p = 9.d0
  find_liquid = .false.

  ! do the loops
  solid_f  = thers_solid
  liquid_f = -1.d0
  
  do j=lo(2),hi(2)
    do i=lo(1),hi(1)
    
       sf(i,j) = -1.d0

       if(pf(i,j) >= solid_f) then

       loop1: do jj = j-1,j+1
         do ii = i-1,i+1
            if(pf(ii,jj) < solid_f) then
               find_liquid = .true.
               exit loop1
            end if
         end do
       end do loop1

       if(find_liquid) then
          sf(i,j) = 0.d0
          find_liquid = .false.
       else
          sf(i,j) = 1.d0
       end if

       end if

    end do
  end do

  if(sw_flag(1,1) == 1) then
    i = lo(1)
    do j=lo(2),hi(2)
      if(sf(i,j) == -1.d0) sf(i,j) = -2.d0;
    end do
  end if

  if(sw_flag(1,2) == 1) then
    i = hi(1)
    do j=lo(2),hi(2)
      if(sf(i,j) == -1.d0) sf(i,j) = -2.d0;
    end do
  end if

  if(sw_flag(2,1) == 1) then
    j = lo(2)
    do i=lo(1),hi(1)
      if(sf(i,j) == -1.d0) sf(i,j) = -2.d0;
    end do
  end if

  if(sw_flag(2,2) == 1) then
    j = hi(2)
    do i=lo(1),hi(1)
      if(sf(i,j) == -1.d0) sf(i,j) = -2.d0;
    end do
  end if

end subroutine cal_interface_flag_2d

subroutine cal_interface_flag_3d(sf,pf,ng_p,lo,hi,thers_solid,sw_flag)

  integer           :: lo(3),hi(3),ng_p,sw_flag(3,2)
  double precision  :: pf(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: sf(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: thers_solid

  integer           :: i,j,k,ii,jj,kk
  double precision  :: solid_f, liquid_f, num_p
  logical           :: find_liquid

  find_liquid = .false.
  num_p = 27.d0
  
  ! do the loops
  solid_f  = thers_solid
  liquid_f = -1.d0
  
  do k=lo(3),hi(3)
  do j=lo(2),hi(2)
    do i=lo(1),hi(1)

       sf(i,j,k) = -1.d0

       if(pf(i,j,k) >= solid_f) then

       loop1: do kk = k-1,k+1
          do jj = j-1,j+1
            do ii = i-1,i+1
               if(pf(ii,jj,kk) < solid_f) then
                  find_liquid = .true.
                  exit loop1
               end if
            end do
          end do
       end do loop1

       if(find_liquid) then
          sf(i,j,k) = 0.d0
          find_liquid = .false.
       else
          sf(i,j,k) = 1.d0
       end if

       end if
       
    end do
  end do
  end do

  if(sw_flag(1,1) == 1) then
    i = lo(1)
     do k=lo(3),hi(3)
       do j=lo(2),hi(2)
         if(sf(i,j,k) == -1.d0) sf(i,j,k) = -2.d0;
       end do
     end do
  end if

  if(sw_flag(1,2) == 1) then
    i = hi(1)
     do k=lo(3),hi(3)
       do j=lo(2),hi(2)
         if(sf(i,j,k) == -1.d0) sf(i,j,k) = -2.d0;
       end do
     end do
  end if

  if(sw_flag(2,1) == 1) then
    j = lo(2)
     do k=lo(3),hi(3)
       do i=lo(1),hi(1)
         if(sf(i,j,k) == -1.d0) sf(i,j,k) = -2.d0;
       end do
     end do
  end if

  if(sw_flag(2,2) == 1) then
    j = hi(2)
     do k=lo(3),hi(3)
       do i=lo(1),hi(1)
         if(sf(i,j,k) == -1.d0) sf(i,j,k) = -2.d0;
       end do
     end do
  end if

  if(sw_flag(3,1) == 1) then
    k = lo(3)
     do j=lo(2),hi(2)
       do i=lo(1),hi(1)
         if(sf(i,j,k) == -1.d0) sf(i,j,k) = -2.d0;
       end do
     end do
  end if

  if(sw_flag(3,2) == 1) then
    k = hi(3)
     do j=lo(2),hi(2)
       do i=lo(1),hi(1)
         if(sf(i,j,k) == -1.d0) sf(i,j,k) = -2.d0;
       end do
     end do
  end if

end subroutine cal_interface_flag_3d

! Use three functions to specify the boundary conditions, using one appears to be two nested
! use u_inlet_lo and u_inlet_hi to specify the density/pressure boundary conditions 
subroutine update_boundary_2d_X(u_vec,f_mat,sf,ng_p,lo,hi,u_inlet_lo,u_inlet_hi,sw_flag,bc_flag,&
                                boundary_type,q,ex,ey,w)
  integer           :: lo(2),hi(2),ng_p,sw_flag(2,2),bc_flag(2),boundary_type,q
  integer           :: ex(q),ey(q)
  double precision  :: sf(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: u_vec(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,1:3)
  double precision  :: f_mat(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,1:9)
  double precision  :: u_inlet_lo, u_inlet_hi,w(q)

  ! local variables
  integer,parameter :: sq=3,sq_p=3 
  integer           :: i,j,iter,is,ie,js,je,ccc
  double precision  :: rho_u,v_specify,var_1,var_2
  double precision  :: var_temp_1,var_temp_2,var_temp_3,var_temp_4
  integer           :: dsi(sq),dso(sq),dsp(sq_p)

  dsi(1:sq)=(/2,6,9/)
  dso(1:sq)=(/4,8,7/)
  dsp(1:sq_p)=(/1,3,5/)
  
  select case (boundary_type)
  case (0) ! velocity is specified
  if(sw_flag(1,1) == 1 .and. bc_flag(1) /= PERIODIC) then
    ! the index of x
    i = lo(1)

    js = lo(2)
    je = hi(2)
	
	  if(bc_flag(1) == INLET) then
      v_specify = u_inlet_lo  

      do j=js,je
        if(sf(i,j) == -2.d0) u_vec(i,j,1) = v_specify;
      end do
    
    else if(bc_flag(1) == OUTLET) then

      do j=js,je
        if(sf(i,j) == -2.d0) u_vec(i,j,1) = u_vec(i+1,j,1);
      end do

    else if(bc_flag(1) == NO_SLIP_WALL) then
      v_specify = 0.d0 

      do j=js,je
        if(sf(i,j) == -2.d0) u_vec(i,j,1:2) = v_specify;
      end do

    else if(bc_flag(1) == SYMMETRY) then

      do j=js,je
        if(sf(i,j) == -2.d0) then
           u_vec(i,j,2) = u_vec(i+1,j,2)
           u_vec(i,j,1) = 0.d0
        end if
      end do

    end if
	
    ! do the loop for y and z
    do j=js,je
      if(sf(i,j) == -2.d0) then
        
        u_vec(i,j,3) = u_vec(i+1,j,3)	

        var_temp_2 = u_vec(i,j,1)**2+u_vec(i,j,2)**2
        var_temp_4 = u_vec(i+1,j,1)**2+u_vec(i+1,j,2)**2
        do iter=1,q
          var_temp_1 = ex(iter)*u_vec(i,j,1)+ey(iter)*u_vec(i,j,2)
          var_temp_3 = ex(iter)*u_vec(i+1,j,1)+ey(iter)*u_vec(i+1,j,2)
          f_mat(i,j,iter) = w(iter)*u_vec(i,  j,3)*(1.0+3.0*var_temp_1+4.5*var_temp_1**2-1.5*var_temp_2)+&
          f_mat(i+1,j,iter)-w(iter)*u_vec(i+1,j,3)*(1.0+3.0*var_temp_3+4.5*var_temp_3**2-1.5*var_temp_4)
        end do  

      end if
      
    end do

  end if

  if(sw_flag(1,2) == 1 .and. bc_flag(2) /= PERIODIC) then
    i = hi(1)

    js = lo(2) 
    je = hi(2)

    if(bc_flag(2) == OUTLET ) then

      do j=js,je
        if(sf(i,j) == -2.d0) u_vec(i,j,1) = u_vec(i-1,j,1);
      end do

    else if(bc_flag(2) == INLET) then
      v_specify = u_inlet_hi  

      do j=js,je
        if(sf(i,j) == -2.d0) u_vec(i,j,1) = v_specify;
      end do

    else if(bc_flag(2) == NO_SLIP_WALL) then
      v_specify = 0.d0  

      do j=js,je
        if(sf(i,j) == -2.d0) u_vec(i,j,1:2) = v_specify;
      end do

    else if(bc_flag(2) == SYMMETRY) then

      do j=js,je
        if(sf(i,j) == -2.d0) then
          u_vec(i,j,2) = u_vec(i-1,j,2)
          u_vec(i,j,1) = 0.d0
        end if
      end do

    end if

    do j=js,je
	    if(sf(i,j) == -2.d0) then
	      
        u_vec(i,j,3) = u_vec(i-1,j,3)
          
        var_temp_2 = u_vec(i,j,1)**2+u_vec(i,j,2)**2
        var_temp_4 = u_vec(i-1,j,1)**2+u_vec(i-1,j,2)**2
        do iter=1,q
          var_temp_1 = ex(iter)*u_vec(i,j,1)+ey(iter)*u_vec(i,j,2)
          var_temp_3 = ex(iter)*u_vec(i-1,j,1)+ey(iter)*u_vec(i-1,j,2)
          f_mat(i,j,iter) = w(iter)*u_vec(i,  j,3)*(1.0+3.0*var_temp_1+4.5*var_temp_1**2-1.5*var_temp_2)+&
          f_mat(i-1,j,iter)-w(iter)*u_vec(i-1,j,3)*(1.0+3.0*var_temp_3+4.5*var_temp_3**2-1.5*var_temp_4)
        end do

	    end if 
    end do

  end if
  
  case (1) ! pressure is specified
  if(sw_flag(1,1) == 1 .and. bc_flag(1) /= PERIODIC) then
    ! the index of x
    i = lo(1)

    js = lo(2)
    je = hi(2)
	
	  if(bc_flag(1) == INLET .or. bc_flag(1) == OUTLET) then
      v_specify = u_inlet_lo  

      do j=js,je
        if(sf(i,j) == -2.d0) then
          u_vec(i,j,3) = v_specify
          !u_vec(i,j,1) = u_vec(i+1,j,1)
          !u_vec(i,j,2) = 0.d0
        end if
      end do

    else if(bc_flag(1) == NO_SLIP_WALL) then

      do j=js,je
        if(sf(i,j) == -2.d0) then
          !u_vec(i,j,3) = u_vec(i+1,j,3)
          u_vec(i,j,1:2) = 0.d0
        end if
      end do

    end if
	
    ! do the loop for y and z
    do j=js,je
      if(sf(i,j) == -2.d0) then
        
        !u_vec(i,j,3) = u_vec(i+1,j,3)	
        
        var_temp_2 = u_vec(i,j,1)**2+u_vec(i,j,2)**2
        var_temp_4 = u_vec(i+1,j,1)**2+u_vec(i+1,j,2)**2
        do iter=1,q
          var_temp_1 = ex(iter)*u_vec(i,j,1)+ey(iter)*u_vec(i,j,2)
          var_temp_3 = ex(iter)*u_vec(i+1,j,1)+ey(iter)*u_vec(i+1,j,2)
          f_mat(i,j,iter) = w(iter)*u_vec(i,  j,3)*(1.0+3.0*var_temp_1+4.5*var_temp_1**2-1.5*var_temp_2)+&
          f_mat(i+1,j,iter)-w(iter)*u_vec(i+1,j,3)*(1.0+3.0*var_temp_3+4.5*var_temp_3**2-1.5*var_temp_4)
        end do  

      end if
      
    end do

  end if

  if(sw_flag(1,2) == 1 .and. bc_flag(2) /= PERIODIC) then
    i = hi(1)

    js = lo(2) 
    je = hi(2)

    if(bc_flag(2) == OUTLET .or. bc_flag(2) == INLET) then
      v_specify = u_inlet_hi

      do j=js,je
        if(sf(i,j) == -2.d0) then
          u_vec(i,j,3) = v_specify
          !u_vec(i,j,1) = u_vec(i-1,j,1)
          !u_vec(i,j,2) = 0.d0
        end if
      end do

    else if(bc_flag(2) == NO_SLIP_WALL) then

      do j=js,je
        if(sf(i,j) == -2.d0) then
          !u_vec(i,j,3) = u_vec(i-1,j,3)
          u_vec(i,j,1:2) = 0.d0
        end if
      end do

    end if

    do j=js,je
	    if(sf(i,j) == -2.d0) then
	      
        !u_vec(i,j,3) = u_vec(i-1,j,3)
          
        var_temp_2 = u_vec(i,j,1)**2+u_vec(i,j,2)**2
        var_temp_4 = u_vec(i-1,j,1)**2+u_vec(i-1,j,2)**2
        do iter=1,q
          var_temp_1 = ex(iter)*u_vec(i,j,1)+ey(iter)*u_vec(i,j,2)
          var_temp_3 = ex(iter)*u_vec(i-1,j,1)+ey(iter)*u_vec(i-1,j,2)
          f_mat(i,j,iter) = w(iter)*u_vec(i,  j,3)*(1.0+3.0*var_temp_1+4.5*var_temp_1**2-1.5*var_temp_2)+&
          f_mat(i-1,j,iter)-w(iter)*u_vec(i-1,j,3)*(1.0+3.0*var_temp_3+4.5*var_temp_3**2-1.5*var_temp_4)
        end do

	    end if 
    end do

  end if
  end select

end subroutine update_boundary_2d_X

! Use three functions to specify the boundary conditions, using one appears to be two nested
! use u_inlet_lo and u_inlet_hi to specify the density/pressure boundary conditions 
subroutine update_boundary_2d_Y(u_vec,f_mat,sf,ng_p,lo,hi,u_inlet_lo,u_inlet_hi,sw_flag,bc_flag,&
                                boundary_type,q,ex,ey,w)

  integer           :: lo(2),hi(2),ng_p,sw_flag(2,2),bc_flag(2),boundary_type,q
  integer           :: ex(q), ey(q)
  double precision  :: sf(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: u_vec(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,1:3)
  double precision  :: f_mat(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,1:9)
  double precision  :: u_inlet_lo, u_inlet_hi,w(q)

  ! local variables
  integer,parameter :: sq=3,sq_p=3 
  integer           :: i,j,iter,is,ie,js,je
  double precision  :: rho_u,v_specify,var_1,var_2,var_3,var_4
  double precision  :: var_temp_1,var_temp_2,var_temp_3,var_temp_4
  integer           :: dsi(sq),dso(sq),dsp(sq_p)

  dsi(1:sq)=(/3,6,7/)
  dso(1:sq)=(/5,8,9/)
  dsp(1:sq_p)=(/1,2,4/)
  
  select case (boundary_type)
  case (0) ! velocity is specified
  if(sw_flag(2,1) == 1 .and. bc_flag(1) /= PERIODIC) then
    ! the index of x
    j = lo(2)
	
    is = lo(1)
    ie = hi(1)

	  if(bc_flag(1) == INLET) then

      v_specify = u_inlet_lo

      do i=is,ie
	      if(sf(i,j) == -2.d0) u_vec(i,j,2) = v_specify;

      end do

    else if(bc_flag(1) == OUTLET ) then

      do i=is,ie
	      if(sf(i,j) == -2.d0) u_vec(i,j,2) = u_vec(i,j+1,2);
      end do

    else if(bc_flag(1) == NO_SLIP_WALL) then
      v_specify = 0.d0

      do i=is,ie
	      if(sf(i,j) == -2.d0) u_vec(i,j,1:2) = v_specify;

      end do

    else if(bc_flag(1) == SYMMETRY) then

      do i=is,ie
	      if(sf(i,j) == -2.d0) then
          u_vec(i,j,1) = u_vec(i,j+1,1)
          u_vec(i,j,2) = 0.d0
        end if
      end do

    end if

    do i=is,ie
	    if(sf(i,j) == -2.d0) then     
        u_vec(i,j,3) = u_vec(i,j+1,3)

        var_temp_2 = u_vec(i,j,1)**2+u_vec(i,j,2)**2
        var_temp_4 = u_vec(i,j+1,1)**2+u_vec(i,j+1,2)**2
        
        do iter=1,q
          var_temp_1 = ex(iter)*u_vec(i,j,1)+ey(iter)*u_vec(i,j,2)
          var_temp_3 = ex(iter)*u_vec(i,j+1,1)+ey(iter)*u_vec(i,j+1,2)

          f_mat(i,j,iter) = w(iter)*u_vec(i,j,  3)*(1.0+3.0*var_temp_1+4.5*var_temp_1**2-1.5*var_temp_2)+&
          f_mat(i,j+1,iter)-w(iter)*u_vec(i,j+1,3)*(1.0+3.0*var_temp_3+4.5*var_temp_3**2-1.5*var_temp_4)
        end do
	  
	    end if
    end do

  end if

  if(sw_flag(2,2) == 1 .and. bc_flag(2) /= PERIODIC) then
     j = hi(2)

     is = lo(1)
     ie = hi(1)

     if(bc_flag(2) == OUTLET) then 

       do i=is,ie
         if(sf(i,j) == -2.d0) u_vec(i,j,2) = u_vec(i,j-1,2);
       end do

     else if(bc_flag(2) == INLET) then
       v_specify = u_inlet_hi

      do i=is,ie
	      if(sf(i,j) == -2.d0) u_vec(i,j,2) = v_specify;
      end do

      else if(bc_flag(2) == NO_SLIP_WALL) then
        v_specify = 0.d0

        do i=is,ie
	        if(sf(i,j) == -2.d0) u_vec(i,j,1:2) = v_specify;
        end do

      else if(bc_flag(2) == SYMMETRY) then
        
        do i=is,ie
         if(sf(i,j) == -2.d0) then
           u_vec(i,j,1) = u_vec(i,j-1,1)
           u_vec(i,j,2) = 0.d0
         end if
        end do

     end if
	  
     do i=is,ie
	
	     if(sf(i,j) == -2.d0) then
	        
         u_vec(i,j,3) = u_vec(i,j-1,3)

         var_temp_2 = u_vec(i,j,1)**2+u_vec(i,j,2)**2
         var_temp_4 = u_vec(i,j-1,1)**2+u_vec(i,j-1,2)**2
         
         do iter=1,q
            var_temp_1 = ex(iter)*u_vec(i,j,1)+ey(iter)*u_vec(i,j,2)
            var_temp_3 = ex(iter)*u_vec(i,j-1,1)+ey(iter)*u_vec(i,j-1,2)
            f_mat(i,j,iter) = w(iter)*u_vec(i,j,  3)*(1.0+3.0*var_temp_1+4.5*var_temp_1**2-1.5*var_temp_2)+&
            f_mat(i,j-1,iter)-w(iter)*u_vec(i,j-1,3)*(1.0+3.0*var_temp_3+4.5*var_temp_3**2-1.5*var_temp_4)
         end do
	  
	     end if
      
     end do

  end if
  
  case (1) ! pressure is specified
  if(sw_flag(2,1) == 1 .and. bc_flag(1) /= PERIODIC) then
    ! the index of x
    j = lo(2)
	
    is = lo(1)
    ie = hi(1)

	  if(bc_flag(1) == INLET .or. bc_flag(1) == OUTLET) then

      v_specify = u_inlet_lo

      do i=is,ie
	      if(sf(i,j) == -2.d0) then
          u_vec(i,j,3) = v_specify
          !u_vec(i,j,2) = u_vec(i,j+1,2)
          !u_vec(i,j,1) = 0.d0
        end if

      end do

    else if(bc_flag(1) == NO_SLIP_WALL) then

       do i=is,ie
	      if(sf(i,j) == -2.d0) then
          !u_vec(i,j,3) = u_vec(i,j+1,3)
          u_vec(i,j,1:2) = 0.d0
        end if
      end do    

    end if

    do i=is,ie
	    if(sf(i,j) == -2.d0) then     
        !u_vec(i,j,3) = u_vec(i,j+1,3)

        var_temp_2 = u_vec(i,j,1)**2+u_vec(i,j,2)**2
        var_temp_4 = u_vec(i,j+1,1)**2+u_vec(i,j+1,2)**2
        
        do iter=1,q
          var_temp_1 = ex(iter)*u_vec(i,j,1)+ey(iter)*u_vec(i,j,2)
          var_temp_3 = ex(iter)*u_vec(i,j+1,1)+ey(iter)*u_vec(i,j+1,2)

          f_mat(i,j,iter) = w(iter)*u_vec(i,j,  3)*(1.0+3.0*var_temp_1+4.5*var_temp_1**2-1.5*var_temp_2)+&
          f_mat(i,j+1,iter)-w(iter)*u_vec(i,j+1,3)*(1.0+3.0*var_temp_3+4.5*var_temp_3**2-1.5*var_temp_4)
        end do
	  
	    end if
    end do

  end if

  if(sw_flag(2,2) == 1 .and. bc_flag(2) /= PERIODIC) then
     j = hi(2)

     is = lo(1)
     ie = hi(1)

     if(bc_flag(2) == OUTLET .or. bc_flag(2) == INLET) then 

       v_specify = u_inlet_hi

       do i=is,ie
         if(sf(i,j) == -2.d0) then
           u_vec(i,j,3) = v_specify
           !u_vec(i,j,2) = u_vec(i,j-1,2)
           !u_vec(i,j,1) = 0.d0
         end if
       end do

     else if(bc_flag(2) == NO_SLIP_WALL) then

       do i=is,ie
	      if(sf(i,j) == -2.d0) then
          !u_vec(i,j,3) = u_vec(i,j-1,3)
          u_vec(i,j,1:2) = 0.d0
        end if
      end do    

     end if
	  
     do i=is,ie
	
	     if(sf(i,j) == -2.d0) then
	       
         !u_vec(i,j,3) = u_vec(i,j-1,3)

          var_temp_2 = u_vec(i,j,1)**2+u_vec(i,j,2)**2
          var_temp_4 = u_vec(i,j-1,1)**2+u_vec(i,j-1,2)**2
          do iter=1,q
            var_temp_1 = ex(iter)*u_vec(i,j,1)+ey(iter)*u_vec(i,j,2)
            var_temp_3 = ex(iter)*u_vec(i,j-1,1)+ey(iter)*u_vec(i,j-1,2)
            f_mat(i,j,iter) = w(iter)*u_vec(i,j,  3)*(1.0+3.0*var_temp_1+4.5*var_temp_1**2-1.5*var_temp_2)+&
            f_mat(i,j-1,iter)-w(iter)*u_vec(i,j-1,3)*(1.0+3.0*var_temp_3+4.5*var_temp_3**2-1.5*var_temp_4)
          end do
	  
	      end if
      
      end do

    end if
  
  end select

end subroutine update_boundary_2d_Y

! The boundary condition seems quite complicated, it seems better if i split it into three functions.
! Two kinds of boundary conditions are supported, either velocity or pressure can be specified on the boundary.
! Beside INLET and OUTLET, NO_SLIP_WALL can also be specified.
! use u_inlet_lo and u_inlet_hi to specify the density/pressure boundary conditions 
subroutine update_boundary_3d_X(u_vec,f_mat,sf,ng_p,lo,hi,u_inlet_lo,u_inlet_hi,sw_flag,bc_flag,&
                                boundary_type,q,ex,ey,ez,w)

  integer           :: lo(3),hi(3),ng_p,sw_flag(3,2),bc_flag(2),boundary_type,q
  integer           :: ex(q),ey(q),ez(q)
  double precision  :: sf(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: u_vec(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p,1:4)
  double precision  :: f_mat(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p,1:19)
  double precision  :: u_inlet_lo, u_inlet_hi,w(q)

  ! local variables
  integer,parameter :: sq=5,sq_p=9 
  integer           :: i,j,k,iter
  double precision  :: rho_u,v_spec
  double precision  :: var_temp_1,var_temp_2,var_temp_3,var_temp_4
  integer           :: dsi(sq),dso(sq), dsp(sq_p)
  integer           :: is,ie,js,je,ks,ke

  dsi(1:sq)=(/2, 8, 9,12,14/)
  dso(1:sq)=(/3,11,10,15,13/)
  dsp(1:sq_p)=(/1,4,5,6,7,16,19,18,17/)
  
  select case (boundary_type)
  case (0) ! velocity is specified
  if(sw_flag(1,1) == 1 .and. bc_flag(1) /= PERIODIC) then
    ! the index of x
    i = lo(1)

    js = lo(2) 
    je = hi(2) 
    ks = lo(3)
    ke = hi(3)
  
    if(bc_flag(1) == INLET ) then
	    v_spec = u_inlet_lo
      
      do k=ks,ke
        do j=js,je
	        if(sf(i,j,k) == -2.d0) u_vec(i,j,k,1) = v_spec;
        end do
      end do

    else if(bc_flag(1) == OUTLET) then

      do k=ks,ke
        do j=js,je
	        if(sf(i,j,k) == -2.d0) u_vec(i,j,k,1) = u_vec(i+1,j,k,1);
        end do
      end do

    else if(bc_flag(1) == NO_SLIP_WALL) then
      v_spec = 0.d0
      
      do k=ks,ke
        do j=js,je
	        if(sf(i,j,k) == -2.d0) u_vec(i,j,k,1:3) = v_spec;
        end do
      end do

    else if(bc_flag(1) == SYMMETRY) then

      do k=ks,ke
        do j=js,je
	        if(sf(i,j,k) == -2.d0) then
            u_vec(i,j,k,2) = u_vec(i+1,j,k,2)
            u_vec(i,j,k,3) = u_vec(i+1,j,k,3)
            u_vec(i,j,k,1) = 0.d0
          end if
        end do
      end do

    end if
	  
      ! do the loop for y and z
      do k=ks,ke
      do j=js,je
	
	    if(sf(i,j,k) == -2.d0) then

        u_vec(i,j,k,4) = u_vec(i+1,j,k,4)

        var_temp_2 = u_vec(i,j,k,1)**2+u_vec(i,j,k,2)**2+u_vec(i,j,k,3)**2
        var_temp_4 = u_vec(i+1,j,k,1)**2+u_vec(i+1,j,k,2)**2+u_vec(i+1,j,k,3)**2
        do iter=1,q
          var_temp_1 = ex(iter)*u_vec(i,j,k,1)+ey(iter)*u_vec(i,j,k,2)+ez(iter)*u_vec(i,j,k,3)
          var_temp_3 = ex(iter)*u_vec(i+1,j,k,1)+ey(iter)*u_vec(i+1,j,k,2)+ez(iter)*u_vec(i+1,j,k,3)
          f_mat(i,j,k,iter) = w(iter)*u_vec(i,  j,k,4)*(1.0+3.0*var_temp_1+4.5*var_temp_1**2-1.5*var_temp_2)+&
          f_mat(i+1,j,k,iter)-w(iter)*u_vec(i+1,j,k,4)*(1.0+3.0*var_temp_3+4.5*var_temp_3**2-1.5*var_temp_4)
        end do
          
	    end if
      
      end do
      end do
  
  end if

  if(sw_flag(1,2) == 1 .and. bc_flag(2) /= PERIODIC) then
    i = hi(1)

    js = lo(2) 
    je = hi(2) 
    ks = lo(3) 
    ke = hi(3)
    
    if(bc_flag(2) == OUTLET) then

      do k=ks,ke
        do j=js,je
	        if(sf(i,j,k) == -2.d0) u_vec(i,j,k,1) = u_vec(i-1,j,k,1);
        end do
      end do

    else if(bc_flag(2) == INLET) then
      v_spec = u_inlet_hi
      
      do k=ks,ke
        do j=js,je
	        if(sf(i,j,k) == -2.d0) u_vec(i,j,k,1) = v_spec;
        end do
      end do

    else if(bc_flag(2) == NO_SLIP_WALL) then
      v_spec = 0.d0
      
      do k=ks,ke
        do j=js,je
	        if(sf(i,j,k) == -2.d0) u_vec(i,j,k,1:3) = v_spec;
        end do
      end do

    else if(bc_flag(2) == SYMMETRY) then

      do k=ks,ke
        do j=js,je
	        if(sf(i,j,k) == -2.d0) then
            u_vec(i,j,k,2) = u_vec(i-1,j,k,2)
            u_vec(i,j,k,3) = u_vec(i-1,j,k,3)
            u_vec(i,j,k,1) = 0.d0
          end if
        end do
      end do      

    end if
		
      do k=ks,ke
      do j=js,je
	  
	    if(sf(i,j,k) == -2.d0) then

        u_vec(i,j,k,4) = u_vec(i-1,j,k,4)

        var_temp_2 = u_vec(i,j,k,1)**2+u_vec(i,j,k,2)**2+u_vec(i,j,k,3)**2
        var_temp_4 = u_vec(i-1,j,k,1)**2+u_vec(i-1,j,k,2)**2+u_vec(i-1,j,k,3)**2
        do iter=1,q
          var_temp_1 = ex(iter)*u_vec(i,j,k,1)+ey(iter)*u_vec(i,j,k,2)+ez(iter)*u_vec(i,j,k,3)  
          var_temp_3 = ex(iter)*u_vec(i-1,j,k,1)+ey(iter)*u_vec(i-1,j,k,2)+ez(iter)*u_vec(i-1,j,k,3)    
          f_mat(i,j,k,iter) = w(iter)*u_vec(i,  j,k,4)*(1.0+3.0*var_temp_1+4.5*var_temp_1**2-1.5*var_temp_2)+&
          f_mat(i-1,j,k,iter)-w(iter)*u_vec(i-1,j,k,4)*(1.0+3.0*var_temp_3+4.5*var_temp_3**2-1.5*var_temp_4)
        end do
	  
	    end if
      
      end do
      end do

  end if
  
  case (1) ! density or pressure is sepcified
  if(sw_flag(1,1) == 1 .and. bc_flag(1) /= PERIODIC) then
    ! the index of x
    i = lo(1)

    js = lo(2) 
    je = hi(2) 
    ks = lo(3)
    ke = hi(3)
  
    if(bc_flag(1) == INLET .or. bc_flag(1) == OUTLET ) then
	    v_spec = u_inlet_lo
      
      do k=ks,ke
        do j=js,je
	        if(sf(i,j,k) == -2.d0) u_vec(i,j,k,4) = v_spec;
        end do
      end do

    else if(bc_flag(1) == NO_SLIP_WALL) then
      v_spec = 0.d0
      
      do k=ks,ke
        do j=js,je
	        if(sf(i,j,k) == -2.d0) u_vec(i,j,k,1:3) = v_spec;
        end do
      end do

    end if
	  
      ! do the loop for y and z
      do k=ks,ke
      do j=js,je
	
	    if(sf(i,j,k) == -2.d0) then

        !u_vec(i,j,k,4) = u_vec(i+1,j,k,4)

        var_temp_2 = u_vec(i,j,k,1)**2+u_vec(i,j,k,2)**2+u_vec(i,j,k,3)**2
        var_temp_4 = u_vec(i+1,j,k,1)**2+u_vec(i+1,j,k,2)**2+u_vec(i+1,j,k,3)**2
        do iter=1,q
          var_temp_1 = ex(iter)*u_vec(i,j,k,1)+ey(iter)*u_vec(i,j,k,2)+ez(iter)*u_vec(i,j,k,3)
          var_temp_3 = ex(iter)*u_vec(i+1,j,k,1)+ey(iter)*u_vec(i+1,j,k,2)+ez(iter)*u_vec(i+1,j,k,3)
          f_mat(i,j,k,iter) = w(iter)*u_vec(i,  j,k,4)*(1.0+3.0*var_temp_1+4.5*var_temp_1**2-1.5*var_temp_2)+&
          f_mat(i+1,j,k,iter)-w(iter)*u_vec(i+1,j,k,4)*(1.0+3.0*var_temp_3+4.5*var_temp_3**2-1.5*var_temp_4)
        end do
          
	    end if
      
      end do
      end do
  
  end if

  if(sw_flag(1,2) == 1 .and. bc_flag(2) /= PERIODIC) then
    i = hi(1)

    js = lo(2) 
    je = hi(2) 
    ks = lo(3) 
    ke = hi(3)
    
    if(bc_flag(2) == OUTLET .or. bc_flag(2) == INLET) then

      v_spec = u_inlet_hi

      do k=ks,ke
        do j=js,je
	        if(sf(i,j,k) == -2.d0) u_vec(i,j,k,4) = v_spec;
        end do
      end do

    else if(bc_flag(2) == NO_SLIP_WALL) then
      v_spec = 0.d0
      
      do k=ks,ke
        do j=js,je
	        if(sf(i,j,k) == -2.d0) u_vec(i,j,k,1:3) = v_spec;
        end do
      end do

    end if
		
      do k=ks,ke
      do j=js,je
	  
	    if(sf(i,j,k) == -2.d0) then

        !u_vec(i,j,k,4) = u_vec(i-1,j,k,4)

        var_temp_2 = u_vec(i,j,k,1)**2+u_vec(i,j,k,2)**2+u_vec(i,j,k,3)**2
        var_temp_4 = u_vec(i-1,j,k,1)**2+u_vec(i-1,j,k,2)**2+u_vec(i-1,j,k,3)**2
        do iter=1,q
          var_temp_1 = ex(iter)*u_vec(i,j,k,1)+ey(iter)*u_vec(i,j,k,2)+ez(iter)*u_vec(i,j,k,3)  
          var_temp_3 = ex(iter)*u_vec(i-1,j,k,1)+ey(iter)*u_vec(i-1,j,k,2)+ez(iter)*u_vec(i-1,j,k,3)    
          f_mat(i,j,k,iter) = w(iter)*u_vec(i,  j,k,4)*(1.0+3.0*var_temp_1+4.5*var_temp_1**2-1.5*var_temp_2)+&
          f_mat(i-1,j,k,iter)-w(iter)*u_vec(i-1,j,k,4)*(1.0+3.0*var_temp_3+4.5*var_temp_3**2-1.5*var_temp_4)
        end do
	  
	    end if
      
      end do
      end do

  end if
  
  end select

end subroutine update_boundary_3d_X

! Use three functions to specify the boundary conditions, using one appears to be two nested
subroutine update_boundary_3d_Y(u_vec,f_mat,sf,ng_p,lo,hi,u_inlet_lo,u_inlet_hi,sw_flag,bc_flag,&
                                boundary_type,q,ex,ey,ez,w)

  integer           :: lo(3),hi(3),ng_p,sw_flag(3,2),bc_flag(2),boundary_type,q
  integer           :: ex(q), ey(q), ez(q)
  double precision  :: sf(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: u_vec(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p,1:4)
  double precision  :: f_mat(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p,1:19)
  double precision  :: u_inlet_lo, u_inlet_hi, w(q)

  ! local variables
  integer,parameter :: sq=5,sq_p=9
  integer           :: i,j,k,iter
  double precision  :: rho_u,v_spec
  double precision  :: var_temp_1,var_temp_2,var_temp_3,var_temp_4
  integer           :: dsi(sq),dso(sq), dsp(sq_p)
  integer           :: is,ie,js,je,ks,ke

  dsi(1:sq)=(/4,8,10,16,17/)
  dso(1:sq)=(/5,11,9,19,18/) 
  dsp(1:sq_p)=(/1,2,3,6,7,13,14,12,15/)
  
  select case (boundary_type)
  case (0) ! velocity is specified
  if(sw_flag(2,1) == 1 .and. bc_flag(1) /= PERIODIC) then
    j = lo(2)

    is = lo(1) 
    ie = hi(1) 
    ks = lo(3) 
    ke = hi(3) 
  
    if(bc_flag(1) == INLET) then
    
	    v_spec = u_inlet_lo

      do k=ks,ke
        do i=is,ie
	        if(sf(i,j,k) == -2.d0) u_vec(i,j,k,2) = v_spec;
        end do
      end do

    else if(bc_flag(1) == OUTLET) then
      
      do k=ks,ke
        do i=is,ie
	        if(sf(i,j,k) == -2.d0) u_vec(i,j,k,2) = u_vec(i,j+1,k,2);
        end do
      end do

    else if(bc_flag(1) == NO_SLIP_WALL) then
      v_spec = 0.d0

      do k=ks,ke
        do i=is,ie
	        if(sf(i,j,k) == -2.d0) u_vec(i,j,k,1:3) = v_spec;
        end do
      end do

    else if(bc_flag(1) == SYMMETRY) then
      
      do k=ks,ke
        do i=is,ie
	        if(sf(i,j,k) == -2.d0) then
            u_vec(i,j,k,1) = u_vec(i,j+1,k,1)
            u_vec(i,j,k,3) = u_vec(i,j+1,k,3)
            u_vec(i,j,k,2) = 0.d0
          end if
        end do
      end do

	  end if  

      do k=ks,ke
      do i=is,ie
	
	    if(sf(i,j,k) == -2.d0) then

        u_vec(i,j,k,4) = u_vec(i,j+1,k,4)

        var_temp_2 = u_vec(i,j,k,1)**2+u_vec(i,j,k,2)**2+u_vec(i,j,k,3)**2
        var_temp_4 = u_vec(i,j+1,k,1)**2+u_vec(i,j+1,k,2)**2+u_vec(i,j+1,k,3)**2
        do iter=1,q
          var_temp_1 = ex(iter)*u_vec(i,j,k,1)+ey(iter)*u_vec(i,j,k,2)+ez(iter)*u_vec(i,j,k,3)
          var_temp_3 = ex(iter)*u_vec(i,j+1,k,1)+ey(iter)*u_vec(i,j+1,k,2)+ez(iter)*u_vec(i,j+1,k,3)
          f_mat(i,j,k,iter) = w(iter)*u_vec(i,j,  k,4)*(1.0+3.0*var_temp_1+4.5*var_temp_1**2-1.5*var_temp_2)+&
          f_mat(i,j+1,k,iter)-w(iter)*u_vec(i,j+1,k,4)*(1.0+3.0*var_temp_3+4.5*var_temp_3**2-1.5*var_temp_4)
        end do
	  
	    end if
      
      end do
      end do

  end if

  if(sw_flag(2,2) == 1 .and. bc_flag(2) /= PERIODIC) then
    j = hi(2)

    is = lo(1) 
    ie = hi(1) 
    ks = lo(3) 
    ke = hi(3) 
    
    if(bc_flag(2) == OUTLET) then

      do k=ks,ke
        do i=is,ie  
	        if(sf(i,j,k) == -2.d0) u_vec(i,j,k,2) = u_vec(i,j-1,k,2);
        end do
      end do

    else if(bc_flag(2) == INLET) then
      v_spec = u_inlet_hi
      
      do k=ks,ke
        do i=is,ie
	        if(sf(i,j,k) == -2.d0) u_vec(i,j,k,2) = v_spec;
        end do
      end do

    else if(bc_flag(2) == NO_SLIP_WALL) then
      v_spec = 0.d0

      do k=ks,ke
        do i=is,ie
	        if(sf(i,j,k) == -2.d0) u_vec(i,j,k,1:3) = v_spec;
        end do
      end do

    else if(bc_flag(2) == SYMMETRY) then

      do k=ks,ke
        do i=is,ie  
	        if(sf(i,j,k) == -2.d0) then
            u_vec(i,j,k,1) = u_vec(i,j-1,k,1)
            u_vec(i,j,k,3) = u_vec(i,j-1,k,3)
            u_vec(i,j,k,2) = 0.d0
          end if
        end do
      end do

    end if

      do k=ks,ke
      do i=is,ie
	  
	    if(sf(i,j,k) == -2.d0) then

        u_vec(i,j,k,4) = u_vec(i,j-1,k,4)

        var_temp_2 = u_vec(i,j,k,1)**2+u_vec(i,j,k,2)**2+u_vec(i,j,k,3)**2
        var_temp_4 = u_vec(i,j-1,k,1)**2+u_vec(i,j-1,k,2)**2+u_vec(i,j-1,k,3)**2
        do iter=1,q
           var_temp_1 = ex(iter)*u_vec(i,j,k,1)+ey(iter)*u_vec(i,j,k,2)+ez(iter)*u_vec(i,j,k,3) 
           var_temp_3 = ex(iter)*u_vec(i,j-1,k,1)+ey(iter)*u_vec(i,j-1,k,2)+ez(iter)*u_vec(i,j-1,k,3)
           f_mat(i,j,k,iter) = w(iter)*u_vec(i,j,  k,4)*(1.0+3.0*var_temp_1+4.5*var_temp_1**2-1.5*var_temp_2)+&
           f_mat(i,j-1,k,iter)-w(iter)*u_vec(i,j-1,k,4)*(1.0+3.0*var_temp_3+4.5*var_temp_3**2-1.5*var_temp_4)
        end do
	  
	    end if
      
      end do
      end do

  end if
  
  case (1) ! pressure is specified
  if(sw_flag(2,1) == 1 .and. bc_flag(1) /= PERIODIC) then
    j = lo(2)

    is = lo(1) 
    ie = hi(1) 
    ks = lo(3) 
    ke = hi(3) 
  
    if(bc_flag(1) == INLET .or. bc_flag(1) == OUTLET ) then
    
	    v_spec = u_inlet_lo

      do k=ks,ke
        do i=is,ie
	        if(sf(i,j,k) == -2.d0) u_vec(i,j,k,4) = v_spec;
        end do
      end do

    else if(bc_flag(1) == NO_SLIP_WALL) then
      v_spec = 0.d0

      do k=ks,ke
        do i=is,ie
	        if(sf(i,j,k) == -2.d0) u_vec(i,j,k,1:3) = v_spec;
        end do
      end do

	  end if  

      do k=ks,ke
      do i=is,ie
	
	    if(sf(i,j,k) == -2.d0) then

        !u_vec(i,j,k,4) = u_vec(i,j+1,k,4)

        var_temp_2 = u_vec(i,j,k,1)**2+u_vec(i,j,k,2)**2+u_vec(i,j,k,3)**2
        var_temp_4 = u_vec(i,j+1,k,1)**2+u_vec(i,j+1,k,2)**2+u_vec(i,j+1,k,3)**2
        do iter=1,q
          var_temp_1 = ex(iter)*u_vec(i,j,k,1)+ey(iter)*u_vec(i,j,k,2)+ez(iter)*u_vec(i,j,k,3)
          var_temp_3 = ex(iter)*u_vec(i,j+1,k,1)+ey(iter)*u_vec(i,j+1,k,2)+ez(iter)*u_vec(i,j+1,k,3)
          f_mat(i,j,k,iter) = w(iter)*u_vec(i,j,  k,4)*(1.0+3.0*var_temp_1+4.5*var_temp_1**2-1.5*var_temp_2)+&
          f_mat(i,j+1,k,iter)-w(iter)*u_vec(i,j+1,k,4)*(1.0+3.0*var_temp_3+4.5*var_temp_3**2-1.5*var_temp_4)
        end do
	  
	    end if
      
      end do
      end do

  end if

  if(sw_flag(2,2) == 1 .and. bc_flag(2) /= PERIODIC) then
    j = hi(2)

    is = lo(1) 
    ie = hi(1) 
    ks = lo(3) 
    ke = hi(3) 
    
    if(bc_flag(2) == OUTLET .or. bc_flag(2) == INLET ) then

       v_spec = u_inlet_hi

      do k=ks,ke
        do i=is,ie  
	        if(sf(i,j,k) == -2.d0) u_vec(i,j,k,4) = v_spec;
        end do
      end do

    else if(bc_flag(2) == NO_SLIP_WALL) then
      v_spec = 0.d0

      do k=ks,ke
        do i=is,ie
	        if(sf(i,j,k) == -2.d0) u_vec(i,j,k,1:3) = v_spec;
        end do
      end do

    end if

      do k=ks,ke
      do i=is,ie
	  
	    if(sf(i,j,k) == -2.d0) then

        !u_vec(i,j,k,4) = u_vec(i,j-1,k,4)

        var_temp_2 = u_vec(i,j,k,1)**2+u_vec(i,j,k,2)**2+u_vec(i,j,k,3)**2
        var_temp_4 = u_vec(i,j-1,k,1)**2+u_vec(i,j-1,k,2)**2+u_vec(i,j-1,k,3)**2
        do iter=1,q
           var_temp_1 = ex(iter)*u_vec(i,j,k,1)+ey(iter)*u_vec(i,j,k,2)+ez(iter)*u_vec(i,j,k,3) 
           var_temp_3 = ex(iter)*u_vec(i,j-1,k,1)+ey(iter)*u_vec(i,j-1,k,2)+ez(iter)*u_vec(i,j-1,k,3)
           f_mat(i,j,k,iter) = w(iter)*u_vec(i,j,  k,4)*(1.0+3.0*var_temp_1+4.5*var_temp_1**2-1.5*var_temp_2)+&
           f_mat(i,j-1,k,iter)-w(iter)*u_vec(i,j-1,k,4)*(1.0+3.0*var_temp_3+4.5*var_temp_3**2-1.5*var_temp_4)
        end do
	  
	    end if
      
      end do
      end do

  end if
  
  end select

end subroutine update_boundary_3d_Y

! Use three functions to specify the boundary conditions, using one appears to be two nested
subroutine update_boundary_3d_Z(u_vec,f_mat,sf,ng_p,lo,hi,u_inlet_lo,u_inlet_hi,sw_flag,bc_flag,&
                                boundary_type,q,ex,ey,ez,w)

  integer           :: lo(3),hi(3),ng_p,sw_flag(3,2),bc_flag(2),boundary_type,q
  integer           :: ex(q),ey(q),ez(q)
  double precision  :: sf(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: u_vec(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p,1:4)
  double precision  :: f_mat(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p,1:19)
  double precision  :: u_inlet_lo, u_inlet_hi,w(q)

  ! local variables
  integer,parameter :: sq=5,sq_p=9
  integer           :: i,j,k,iter
  double precision  :: rho_u,v_spec
  double precision  :: var_temp_1,var_temp_2,var_temp_3,var_temp_4
  integer           :: dsi(sq),dso(sq), dsp(sq_p)
  integer           :: is,ie,js,je,ks,ke
  
  dsi(1:sq)=(/6,12,13,16,18/)
  dso(1:sq)=(/7,15,14,19,17/)
  dsp(1:sq_p)=(/1,2,3,4,5,9,10,8,11/)
  
  select case (boundary_type)
  case (0) ! velocity is specified
  if(sw_flag(3,1) == 1 .and. bc_flag(1) /= PERIODIC) then
    
	  k = lo(3)

    is = lo(1)
    ie = hi(1)
    js = lo(2) 
    je = hi(2) 
	
    if(bc_flag(1) == INLET) then    

	    v_spec = u_inlet_lo
      do j=js,je
        do i=is,ie
	        if(sf(i,j,k) == -2.d0) u_vec(i,j,k,3) = v_spec;
        end do
      end do

    else if(bc_flag(1) == OUTLET) then
      
      do j=js,je
        do i=is,ie
	        if(sf(i,j,k) == -2.d0) u_vec(i,j,k,3) = u_vec(i,j,k+1,3);
        end do
      end do

    else if(bc_flag(1) == NO_SLIP_WALL) then    
	    v_spec = 0.d0

      do j=js,je
        do i=is,ie
	        if(sf(i,j,k) == -2.d0) u_vec(i,j,k,1:3) = v_spec;
        end do
      end do

    else if(bc_flag(1) == SYMMETRY) then
      
      do j=js,je
        do i=is,ie
	        if(sf(i,j,k) == -2.d0) then
            u_vec(i,j,k,1) = u_vec(i,j,k+1,1)
            u_vec(i,j,k,2) = u_vec(i,j,k+1,2)
            u_vec(i,j,k,3) = 0.d0
          end if
        end do
      end do

    end if

      do j=js,je
      do i=is,ie
	
	    if(sf(i,j,k) == -2.d0) then

        u_vec(i,j,k,4) = u_vec(i,j,k+1,4)
        var_temp_2 = u_vec(i,j,k,1)**2+u_vec(i,j,k,2)**2+u_vec(i,j,k,3)**2
        var_temp_4 = u_vec(i,j,k+1,1)**2+u_vec(i,j,k+1,2)**2+u_vec(i,j,k+1,3)**2
        do iter=1,q
          var_temp_1 = ex(iter)*u_vec(i,j,k,1)+ey(iter)*u_vec(i,j,k,2)+ez(iter)*u_vec(i,j,k,3)    
          var_temp_3 = ex(iter)*u_vec(i,j,k+1,1)+ey(iter)*u_vec(i,j,k+1,2)+ez(iter)*u_vec(i,j,k+1,3)
          f_mat(i,j,k,  iter)=w(iter)*u_vec(i,j,k,  4)*(1.0+3.0*var_temp_1+4.5*var_temp_1**2-1.5*var_temp_2)+&
          f_mat(i,j,k+1,iter)-w(iter)*u_vec(i,j,k+1,4)*(1.0+3.0*var_temp_3+4.5*var_temp_3**2-1.5*var_temp_4)
        end do
	  
	    end if
      
      end do
      end do
    
  end if

  if(sw_flag(3,2) == 1 .and. bc_flag(2) /= PERIODIC) then
    k = hi(3)

    is = lo(1) 
    ie = hi(1) 
    js = lo(2) 
    je = hi(2) 
    
    if(bc_flag(2) == OUTLET) then

      do j=js,je
        do i=is,ie
	        if(sf(i,j,k) == -2.d0) u_vec(i,j,k,3) = u_vec(i,j,k-1,3);
        end do
      end do

    else if(bc_flag(2) == INLET) then
      v_spec = u_inlet_hi
      
      do j=js,je
        do i=is,ie
	        if(sf(i,j,k) == -2.d0) u_vec(i,j,k,3) = v_spec;
        end do
      end do
    
    else if(bc_flag(2) == NO_SLIP_WALL) then    
	    v_spec = 0.d0

      do j=js,je
        do i=is,ie
	        if(sf(i,j,k) == -2.d0) u_vec(i,j,k,1:3) = v_spec;
        end do
      end do

    else if(bc_flag(2) == SYMMETRY) then

      do j=js,je
        do i=is,ie
	        if(sf(i,j,k) == -2.d0) then
            u_vec(i,j,k,1) = u_vec(i,j,k-1,1)
            u_vec(i,j,k,2) = u_vec(i,j,k-1,2)
            u_vec(i,j,k,3) = 0.d0
          end if
        end do
      end do

    end if

      do j=js,je
      do i=is,ie
	
	    if(sf(i,j,k) == -2.d0) then

        u_vec(i,j,k,4) = u_vec(i,j,k-1,4)

        var_temp_2 = u_vec(i,j,k,1)**2+u_vec(i,j,k,2)**2+u_vec(i,j,k,3)**2
        var_temp_4 = u_vec(i,j,k-1,1)**2+u_vec(i,j,k-1,2)**2+u_vec(i,j,k-1,3)**2
        do iter=1,q
           var_temp_1 = ex(iter)*u_vec(i,j,k,1)+ey(iter)*u_vec(i,j,k,2)+ez(iter)*u_vec(i,j,k,3)
           var_temp_3 = ex(iter)*u_vec(i,j,k-1,1)+ey(iter)*u_vec(i,j,k-1,2)+ez(iter)*u_vec(i,j,k-1,3)
           f_mat(i,j,k,  iter)=w(iter)*u_vec(i,j,k,  4)*(1.0+3.0*var_temp_1+4.5*var_temp_1**2-1.5*var_temp_2)+&
           f_mat(i,j,k-1,iter)-w(iter)*u_vec(i,j,k-1,4)*(1.0+3.0*var_temp_3+4.5*var_temp_3**2-1.5*var_temp_4)
        end do
	  
	    end if
      
      end do
      end do

  end if
  
  case (1) ! pressure is specified
  if(sw_flag(3,1) == 1 .and. bc_flag(1) /= PERIODIC) then
    
	  k = lo(3)

    is = lo(1)
    ie = hi(1)
    js = lo(2) 
    je = hi(2) 
	
    if(bc_flag(1) == INLET .or. bc_flag(1) == OUTLET) then    

	    v_spec = u_inlet_lo

      do j=js,je
        do i=is,ie
	        if(sf(i,j,k) == -2.d0) u_vec(i,j,k,4) = v_spec;
        end do
      end do

    else if(bc_flag(1) == NO_SLIP_WALL) then    
	    v_spec = 0.d0

      do j=js,je
        do i=is,ie
	        if(sf(i,j,k) == -2.d0) u_vec(i,j,k,1:3) = v_spec;
        end do
      end do

    end if

      do j=js,je
      do i=is,ie
	
	    if(sf(i,j,k) == -2.d0) then

        !u_vec(i,j,k,4) = u_vec(i,j,k+1,4)
        var_temp_2 = u_vec(i,j,k,1)**2+u_vec(i,j,k,2)**2+u_vec(i,j,k,3)**2
        var_temp_4 = u_vec(i,j,k+1,1)**2+u_vec(i,j,k+1,2)**2+u_vec(i,j,k+1,3)**2
        do iter=1,q
          var_temp_1 = ex(iter)*u_vec(i,j,k,1)+ey(iter)*u_vec(i,j,k,2)+ez(iter)*u_vec(i,j,k,3)    
          var_temp_3 = ex(iter)*u_vec(i,j,k+1,1)+ey(iter)*u_vec(i,j,k+1,2)+ez(iter)*u_vec(i,j,k+1,3)
          f_mat(i,j,k,  iter)=w(iter)*u_vec(i,j,k,  4)*(1.0+3.0*var_temp_1+4.5*var_temp_1**2-1.5*var_temp_2)+&
          f_mat(i,j,k+1,iter)-w(iter)*u_vec(i,j,k+1,4)*(1.0+3.0*var_temp_3+4.5*var_temp_3**2-1.5*var_temp_4)
        end do
	  
	    end if
      
      end do
      end do
    
  end if

  if(sw_flag(3,2) == 1 .and. bc_flag(2) /= PERIODIC) then
    k = hi(3)

    is = lo(1) 
    ie = hi(1) 
    js = lo(2) 
    je = hi(2) 
    
    if(bc_flag(2) == OUTLET .or. bc_flag(2) == INLET ) then

      v_spec = u_inlet_hi

      do j=js,je
        do i=is,ie
	        if(sf(i,j,k) == -2.d0) u_vec(i,j,k,4) = v_spec;
        end do
      end do

    else if(bc_flag(2) == NO_SLIP_WALL) then    
	    v_spec = 0.d0

      do j=js,je
        do i=is,ie
	        if(sf(i,j,k) == -2.d0) u_vec(i,j,k,1:3) = v_spec;
        end do
      end do
      
    end if

      do j=js,je
      do i=is,ie
	
	    if(sf(i,j,k) == -2.d0) then

        !u_vec(i,j,k,4) = u_vec(i,j,k-1,4)

        var_temp_2 = u_vec(i,j,k,1)**2+u_vec(i,j,k,2)**2+u_vec(i,j,k,3)**2
        var_temp_4 = u_vec(i,j,k-1,1)**2+u_vec(i,j,k-1,2)**2+u_vec(i,j,k-1,3)**2
        do iter=1,q
           var_temp_1 = ex(iter)*u_vec(i,j,k,1)+ey(iter)*u_vec(i,j,k,2)+ez(iter)*u_vec(i,j,k,3)
           var_temp_3 = ex(iter)*u_vec(i,j,k-1,1)+ey(iter)*u_vec(i,j,k-1,2)+ez(iter)*u_vec(i,j,k-1,3)
           f_mat(i,j,k,  iter)=w(iter)*u_vec(i,j,k,  4)*(1.0+3.0*var_temp_1+4.5*var_temp_1**2-1.5*var_temp_2)+&
           f_mat(i,j,k-1,iter)-w(iter)*u_vec(i,j,k-1,4)*(1.0+3.0*var_temp_3+4.5*var_temp_3**2-1.5*var_temp_4)
        end do
	  
	    end if
      
      end do
      end do

  end if
  
  end select

end subroutine update_boundary_3d_Z


function cal_inlet_velocity_2d(u_inlet,j,dx,y_pd) result(u_v)

    double precision :: u_v
    double precision :: dx, y_pd, u_inlet
    integer          :: j
    double precision :: y_pos = 0.d0

    u_v = u_inlet

    y_pos = dble(j)*dx
    u_v = 4.d0*u_inlet*y_pos*(y_pd-y_pos)/y_pd**2
    
end function cal_inlet_velocity_2d

function cal_inlet_velocity_3d(u_inlet,j,k,dx,y_pd,z_pd) result(u_v)

    double precision :: u_v
    double precision :: dx, y_pd, u_inlet,z_pd
    integer          :: j,k

    double precision :: y_pos = 0.d0,z_pos = 0.d0

    u_v = u_inlet

    y_pos = dble(j)*dx
    z_pos = dble(k)*dx
    u_v = 16.d0*u_inlet*y_pos*z_pos*(y_pd-y_pos)*(z_pd-z_pos)/(y_pd**2 * z_pd**2)
    
end function cal_inlet_velocity_3d

subroutine check_solid_velocity_on_level(flw,phi)
    type(multifab) , intent(inout) :: flw
    type(multifab) , intent(in   ) :: phi

    ! local variables
    integer :: lo(phi%dim), hi(phi%dim)
    integer :: dm, ng_p, i

    real(kind=dp_t), pointer ::  pfo(:,:,:,:)
    real(kind=dp_t), pointer ::  pfl(:,:,:,:)

    dm = flw%dim
    ng_p = flw%ng
    
    do i=1,nfabs(flw)
          
       pfl  => dataptr(flw,i)
       pfo  => dataptr(phi,i)
        
       lo = lwb(get_box(flw,i))
       hi = upb(get_box(flw,i))

       select case(dm)
          case (2)
             call check_solid_velocity_2d(pfl(:,:,1,1),pfl(:,:,1,2),pfo(:,:,1,1),ng_p,lo,hi)
          case (3)
             call check_solid_velocity_3d(pfl(:,:,:,1),pfl(:,:,:,2),pfl(:,:,:,3),pfo(:,:,:,1),ng_p,lo,hi)
       end select

    end do

    call multifab_fill_boundary_c(flw,1 ,dm, ng_p)
    
end subroutine check_solid_velocity_on_level

subroutine check_solid_velocity_2d(ux,uy,sf,ng_p,lo,hi)
  integer           :: lo(2),hi(2),ng_p
  double precision  :: sf(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: ux(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
  double precision  :: uy(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)

  integer           :: i,j

  ! do the loops
  do j=lo(2),hi(2)
    do i=lo(1),hi(1)
    
       if(sf(i,j) >= 0.d0) then
         ux(i,j) = 0.d0
         uy(i,j) = 0.d0
       end if

    end do
  end do

end subroutine check_solid_velocity_2d

subroutine check_solid_velocity_3d(ux,uy,uz,sf,ng_p,lo,hi)

  integer           :: lo(3),hi(3),ng_p
  double precision  :: sf(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: ux(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: uy(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
  double precision  :: uz(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)

  integer           :: i,j,k

  ! do the loops
  do k=lo(3),hi(3)
  do j=lo(2),hi(2)
    do i=lo(1),hi(1)

       if(sf(i,j,k) >= 0.d0) then
         ux(i,j,k) = 0.d0
         uy(i,j,k) = 0.d0
         uz(i,j,k) = 0.d0
       end if
       
    end do
  end do
  end do

end subroutine check_solid_velocity_3d

end module advance_flow_field_module
