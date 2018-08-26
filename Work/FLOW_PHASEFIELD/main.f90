program main

  use boxlib
  use multifab_module
  use bl_IO_module
  use ml_layout_module
  use init_phi_module
  use write_plotfile_module
  use define_bc_module
  use make_new_grids_module
  use regrid_module
  use pf_utility_module
  use phase_field_module
  use solute_field_module
  use temperature_field_module
  use advance_flow_field_module
  use multifab_fill_ghost_module

  implicit none

  ! stuff you can set with the inputs file (otherwise use default values below)
  integer    :: max_levs, dim, nsteps_0,plot_int_1,nsteps, plot_int,status_int, n_cell_x,n_cell_y,n_cell_z
  integer    :: amr_buf_width, cluster_minwidth, cluster_blocking_factor,ada_regriding, temp_seg_num
  integer    :: regrid_int,compute_mode,reset_plotnum,write_file_mode,N_para_Range,write_para_file
  integer    :: bc_x_lo, bc_x_hi, bc_y_lo, bc_y_hi, bc_z_lo, bc_z_hi,regrid_int_1, max_grid_size
  integer    :: regrid_int_max,regrid_int_min,regrid_amp_ref, rotation_type,period_type,period_NT
  integer    :: dt_ratio,seed_num,seed_type,no_noise,plot_mode,coupled_mode,cal_tem_mode,plot_para_num
  integer    :: cooling_mode,un, farg, narg,un1,un2,un3,un4,kill_lorentz_force
  integer    :: flw_calculate_mode,n_comp_flw_uf,cal_boundary_spec,do_solve_flow, solve_flow_only
  integer    :: fl_x_lo, fl_x_hi, fl_y_lo, fl_y_hi, fl_z_lo, fl_z_hi,n_comp_pf,kill_body_force
  integer    :: flw_skip_level_n,kill_phi_flag,i_pro,j_pro, n_pro,n_comp_flw,n_calpf,iter_pf
  integer    :: istep,i,j,n,nl,nlevs, ada_regriding_count,runtime_hours, runtime_minutes, runtime_seconds

  ! stuff you can set with the inputs file, related parameters for phase field
  real(dp_t) :: lamda, anis, anis_0,Le, kk, DM, s_k, Mc00, p_a1, p_a2, thers, noi_amp, dx_def,dt01,dt08
  real(dp_t) :: ori_def,temp_h_l, temp_h_r,temp_l,seed_radius,para_uc,para_th,amr_thers,safe_n_init,t_DM,r_hyper
  real(dp_t) :: Rate_Cooling,temp_lowcut,ratio_select,regrid_time, elapsed_time, start_regrid_time
  real(dp_t) :: cluster_min_eff,v_mag_max,thers_gradient, period_p0,lorentz_omega
  real(dp_t) :: start_flw_time, total_flw_time,gacc, belta_c,multi_add_tau,thers_tau
  real(dp_t) :: u_inlet_x_lo,u_inlet_x_hi,u_inlet_y_lo,u_inlet_y_hi,u_inlet_z_lo,u_inlet_z_hi
  real(dp_t) :: thers_solid,flw_rho,para_u,tau_LBM,flw_amp,dummy,flw_constrain_limit
  real(dp_t) :: dt,time,start_time,run_time,run_time_IOproc,lorentz_amp_x,lorentz_amp_y,lorentz_amp_z
  
  ! dummy indices using for reading in inputs file
  character(len=128) :: inputs_file_name,load_folder_name,parafile_name
  character(len=128) :: status_filename1,status_filename2,status_filename3
  character(len=32)  :: seed_arg_name
  character(len=10)  :: vchar,time_str
  character(len=256) :: myformat

  ! Other variables
  logical                     :: new_grid,need_inputs_file,found_inputs_file,need_regrid,err_msg
  type(box)                   :: bx
  type(ml_boxarray)           :: mba
  type(ml_layout)             :: mla
  type(pf_para)               :: pfpara
  type(bc_tower)              :: the_bc_tower

  ! allocatable variables
  integer       , allocatable :: lo(:), hi(:), phys_bc(:,:)
  logical       , allocatable :: is_periodic(:)
  real(dp_t)    , allocatable :: prob_lo(:), prob_hi(:), tk_tip(:,:),dx(:)
  real(dp_t)    , allocatable :: seed_pos(:,:)

  ! will be allocated with max_levs components
  type(multifab), allocatable :: ad_ori(:),phi_old(:),phi_new(:),noise(:)
  type(multifab), allocatable :: flw(:),flw_old(:),flw_new(:),surf(:)
  type(layout)  , allocatable :: la_array(:)
  type(multifab), pointer     :: phi_load(:) => Null()

  namelist /probin/ max_levs,dim,nsteps_0,nsteps,plot_int,plot_int_1,status_int,n_cell_x,n_cell_y,&
       n_cell_z,max_grid_size,amr_buf_width,cluster_minwidth, cluster_blocking_factor, cluster_min_eff,&
       ada_regriding,regrid_int, regrid_amp_ref, regrid_int_max,regrid_int_min,thers_gradient,plot_mode,compute_mode,&
       coupled_mode, cal_tem_mode,cooling_mode,temp_seg_num,write_file_mode,write_para_file, N_para_Range,&
       plot_para_num, ratio_select, bc_x_lo, bc_x_hi, bc_y_lo, bc_y_hi, bc_z_lo, bc_z_hi, no_noise,lamda,&
       anis,anis_0,Le,kk,s_k,Mc00,p_a1,p_a2,thers,noi_amp,ori_def,temp_h_l,temp_h_r,temp_l,seed_radius,seed_num,&
       seed_type,r_hyper,Rate_Cooling,temp_lowcut,period_type, period_p0, period_NT,rotation_type,dx_def, &
       dt01,dt08,dt_ratio,para_uc,para_th,amr_thers,safe_n_init,load_folder_name,reset_plotnum, do_solve_flow, &
       solve_flow_only,kill_phi_flag,kill_body_force,kill_lorentz_force,tau_LBM,flw_amp, thers_solid,flw_rho, para_u,&
       flw_skip_level_n,gacc, belta_c,multi_add_tau,flw_constrain_limit,lorentz_omega, lorentz_amp_x,lorentz_amp_y,&
       lorentz_amp_z, fl_x_lo,fl_x_hi,fl_y_lo,fl_y_hi,fl_z_lo,fl_z_hi,&
       flw_calculate_mode,tau_LBM,thers_tau, u_inlet_x_lo,u_inlet_x_hi,u_inlet_y_lo,u_inlet_y_hi,&
       u_inlet_z_lo,u_inlet_z_hi,cal_boundary_spec


  ! if running in parallel, this will print out the number of MPI
  call boxlib_initialize()

  ! the program began
  start_time = parallel_wtime()

  err_msg                   = .false.
  n_pro                     = 120
  max_levs                  = 3
  dim                       = 2
  nsteps_0                  = 1000
  nsteps                    = 1000
  plot_int                  = 100
  plot_int_1                = 5000
  status_int                = 10
  n_cell_x                  = 256
  n_cell_y                  = 256
  n_cell_z                  = 256
  max_grid_size             = 64

  n_calpf                   = 1

  amr_buf_width             = 2
  cluster_minwidth          = 16
  cluster_blocking_factor   = 8
  cluster_min_eff           = 0.7d0
  regrid_int                = 10
  compute_mode              = 0

  lamda         = 10.d0
  anis          = 0.15d0
  anis_0        = -0.02d0
  Le            = 1.d0
  kk            = 4.d0
  DM            = 3.1335d0
  s_k           = 0.15d0
  Mc00          = 0.15d0   
  p_a1          = 0.8839d0
  p_a2          = 0.6267d0
  thers         = 1.0E-6
  noi_amp       = 0.d0
  dx_def        = 0.8d0
  ori_def       = 0.d0
  temp_h_l      = 0.7d0
  temp_h_r      = 0.7d0
  temp_l        = 0.7d0
  seed_radius   = 3.d0
  seed_num      = 1
  seed_type     = 0
  r_hyper       = 2.d0
  no_noise      = 0
  plot_mode     = 0
  coupled_mode  = 0
  cal_tem_mode  = 0
  cooling_mode  = 0
  para_uc       = 0.5d0
  para_th       = 0.d0
  Rate_Cooling  = 0.d0
  temp_lowcut   = 0.7d0
  reset_plotnum = 0
  write_file_mode = 0
  N_para_Range    = 0
  plot_para_num   = 1
  write_para_file = 0
  ratio_select    = 0.5
  regrid_amp_ref  = 1
  regrid_int_max  = 200
  regrid_int_min  = 10
  ada_regriding   = 0
  thers_gradient  = 0.d0
  rotation_type   = 0
  period_type     = 1
  period_p0       = 0.d0
  period_NT       = 1

  ! For flow routine
  flw_calculate_mode  = 1
  do_solve_flow       = 0
  solve_flow_only     = 0
  tau_LBM             = 1.d0
  flw_amp             = 1.d0
  thers_solid         = 0.9d0
  flw_rho             = 0.d0
  para_u              = 0.d0
  flw_skip_level_n    = 0
  gacc                = 0.d0
  belta_c             = 0.d0
  multi_add_tau       = 0.0
  flw_constrain_limit = 0.99d0
  kill_phi_flag       = 0
  kill_body_force     = 1
  kill_lorentz_force  = 1
  cal_boundary_spec   = 1
  lorentz_omega       = 0.d0
  lorentz_amp_x      = 0.0
  lorentz_amp_y      = 0.0
  lorentz_amp_z      = 0.0
 
  fl_x_lo       = 15 
  fl_x_hi       = 15 
  fl_y_lo       = 15 
  fl_y_hi       = 15 
  fl_z_lo       = 15 
  fl_z_hi       = 15 

  u_inlet_x_lo  = 1.d0
  u_inlet_x_hi  = 1.d0
  u_inlet_y_lo  = 1.d0
  u_inlet_y_hi  = 1.d0 
  u_inlet_z_lo  = 1.d0
  u_inlet_z_hi  = 1.d0

  ! allowable options for this example are
  bc_x_lo       = -1 
  bc_x_hi       = -1
  bc_y_lo       = -1 
  bc_y_hi       = -1 
  bc_z_lo       = -1 
  bc_z_hi       = -1 

  ! read inputs file and overwrite any default values
  narg = command_argument_count()
  need_inputs_file  = .true.
  farg = 1
  if ( need_inputs_file .AND. narg >= 1 ) then
     call get_command_argument(farg, value = inputs_file_name)
     inquire(file = inputs_file_name, exist = found_inputs_file )
     if ( found_inputs_file ) then
        farg = farg + 1
        un = unit_new()
        open(unit=un, file = inputs_file_name, status = 'old', action = 'read')
        read(unit=un, nml = probin)
        close(unit=un)
        !need_inputs_file = .false.
     end if

     ! allocate the seed information based on the seed_num read from the input file
     ! the vector has dim + 1 dimensions, x, y, z and ori
     if(seed_num .ge. 1) allocate( seed_pos(seed_num,7) )

     call get_command_argument(farg, value = inputs_file_name)
     inquire(file = inputs_file_name, exist = found_inputs_file )
     if (found_inputs_file) then
        farg = farg + 1
        un = unit_new()
        open(unit=un, file = inputs_file_name, status = 'old', action = 'read')

        ! the length of the seed in the file must be longer than the seed_num
        do i=1,seed_num
            read(un,*) seed_pos(i,1), seed_pos(i,2), seed_pos(i,3), &
                       seed_pos(i,4), seed_pos(i,5), seed_pos(i,6), &
                       seed_pos(i,7)
        end do

        close(unit=un)
        need_inputs_file = .false.
     end if

  end if

  ! now that we have dim, we can allocate these
  allocate(lo(dim),hi(dim))
  allocate(is_periodic(dim))
  allocate(prob_lo(dim),prob_hi(dim))
  allocate(phys_bc(dim,2))

  ! now that we have max_levs, we can allocate these
  allocate(dx(max_levs))
  allocate(ad_ori(max_levs))
  allocate(phi_old(max_levs))
  allocate(phi_new(max_levs))
  allocate(noise(max_levs))
  allocate(la_array(max_levs))

  ! If the flow is calculated
  allocate(flw(max_levs))
  allocate(surf(max_levs))
  allocate(flw_old(max_levs))
  allocate(flw_new(max_levs))
  allocate(pfpara%bc_flag(dim,2))
  allocate(pfpara%tau_r(max_levs))

  DM   = lamda * p_a2
  t_DM = Le * DM

  ! Setup the parameters
  pfpara%dim          = dim
  pfpara%lamda        = lamda
  pfpara%anis         = anis
  pfpara%anis_0       = anis_0
  pfpara%Le           = Le
  pfpara%kk           = kk
  pfpara%DM           = DM
  pfpara%sk           = s_k
  pfpara%Mc00         = Mc00
  pfpara%a1           = p_a1
  pfpara%a2           = p_a2
  pfpara%thers        = thers
  pfpara%noiamp       = noi_amp
  pfpara%dx           = dx_def
  pfpara%para_uc      = para_uc
  pfpara%para_th      = para_th
  pfpara%amr_thers    = amr_thers
  pfpara%ori_def      = ori_def
  pfpara%temp_h_l     = temp_h_l
  pfpara%temp_h_r     = temp_h_r
  pfpara%temp_l       = temp_l
  pfpara%seed_radius  = seed_radius
  pfpara%seed_num     = seed_num
  pfpara%no_noise     = no_noise
  pfpara%safe_n_init  = safe_n_init
  pfpara%seed_type    = seed_type
  pfpara%plot_mode    = plot_mode
  pfpara%coupled_mode = coupled_mode
  pfpara%cal_tem_mode = cal_tem_mode
  pfpara%t_DM         = t_DM
  pfpara%r_hyper      = r_hyper
  pfpara%Rate_Cooling = Rate_Cooling
  pfpara%temp_lowcut  = temp_lowcut
  pfpara%cooling_mode = cooling_mode
  pfpara%period_type  = period_type
  pfpara%period_p0    = period_p0
  pfpara%period_NT    = period_NT

  ! For the flow routine
  pfpara%do_solve_flow    = do_solve_flow
  pfpara%solve_flow_only  = solve_flow_only
  pfpara%flw_rho          = flw_rho
  pfpara%para_u           = para_u
  pfpara%thers_solid      = thers_solid
  pfpara%flw_amp          = flw_amp
  pfpara%flw_skip_level_n = flw_skip_level_n
  pfpara%gacc             = gacc
  pfpara%belta_c          = belta_c
  pfpara%dt_flw           = 1.d0
  pfpara%tau_LBM          = tau_LBM
  pfpara%multi_add_tau    = multi_add_tau
  pfpara%kill_phi_flag    = kill_phi_flag
  pfpara%kill_body_force  = kill_body_force
  pfpara%kill_lorentz_force = kill_lorentz_force
  pfpara%lorentz_omega    = lorentz_omega
  pfpara%lorentz_amp_x    = lorentz_amp_x 
  pfpara%lorentz_amp_y    = lorentz_amp_y 
  pfpara%lorentz_amp_z    = lorentz_amp_z 
  pfpara%flw_constrain_limit = flw_constrain_limit

  pfpara%flw_calculate_mode   = flw_calculate_mode
  pfpara%thers_tau            = thers_tau
  pfpara%cal_boundary_spec    = cal_boundary_spec

  pfpara%u_inlet_def(1,1)     = u_inlet_x_lo
  pfpara%u_inlet_def(1,2)     = u_inlet_x_hi
  pfpara%u_inlet_def(2,1)     = u_inlet_y_lo
  pfpara%u_inlet_def(2,2)     = u_inlet_y_hi
  pfpara%u_inlet_def(3,1)     = u_inlet_z_lo
  pfpara%u_inlet_def(3,2)     = u_inlet_z_hi
  pfpara%u_inlet(:,:)         = pfpara%u_inlet_def(:,:)

  select case(dim)
    case (2)
      n_comp_flw    = 21
      n_comp_flw_uf = 12
    case (3)
      n_comp_flw    = 42
      n_comp_flw_uf = 23
  end select

  n_comp_pf = 2

  ! put all the domain boundary conditions into bc_flag
  pfpara%bc_flag(1,1) = fl_x_lo
  pfpara%bc_flag(1,2) = fl_x_hi
  pfpara%bc_flag(2,1) = fl_y_lo
  pfpara%bc_flag(2,2) = fl_y_hi
  if (dim .eq. 3) then
     pfpara%bc_flag(3,1) = fl_z_lo
     pfpara%bc_flag(3,2) = fl_z_hi
  end if

  ! make sure to switch off th tag if not coupled
  if(coupled_mode .eq. 0) then
      pfpara%para_th = 0.d0
  end if

  ! put all the domain boundary conditions into phys_bc
  phys_bc(1,1) = bc_x_lo
  phys_bc(1,2) = bc_x_hi
  phys_bc(2,1) = bc_y_lo
  phys_bc(2,2) = bc_y_hi
  if (dim .eq. 3) then
     phys_bc(3,1) = bc_z_lo
     phys_bc(3,2) = bc_z_hi
  end if

  ! build an array indicating periodicity in each direction
  is_periodic(:) = .false.
  do i=1,dim
     if (phys_bc(i,1) .eq. -1 .and. phys_bc(i,2) .ne. -1) then
        call bl_error("Invalid BC's - both lo and hi need to be periodic")
     end if
     if (phys_bc(i,2) .eq. -1 .and. phys_bc(i,1) .ne. -1) then
        call bl_error("Invalid BC's - both lo and hi need to be periodic")
     end if
     if (phys_bc(i,1) .eq. -1 .and. phys_bc(i,2) .eq. -1) then
        is_periodic(i) = .true.
     end if
  end do

  call cluster_set_minwidth(cluster_minwidth)
  call cluster_set_blocking_factor(cluster_blocking_factor)
  call cluster_set_min_eff(cluster_min_eff)

  prob_lo(:) = 0.d0
  prob_hi(1) = n_cell_x * dx_def
  prob_hi(2) = n_cell_y * dx_def
  if(dim .eq. 3) then
     prob_hi(3) = n_cell_z * dx_def
  end if

  call pf_para_set_flw_tau(pfpara%tau_r, pfpara%tau_LBM, max_levs)

  !initialize the seed inormation
  do i=1,seed_num
    do j=1,dim
       seed_pos(i,j) = prob_lo(j) + seed_pos(i,j) * ( prob_hi(j) - prob_lo(j) )
    end do
  end do

  ! initialize the seed information inside the pfpara
  call pf_para_set_seed_info(pfpara,seed_pos,seed_num,rotation_type)

  ! tell the_bc_tower about max_levs, dim, and phys_bc
  call bc_tower_init(the_bc_tower,max_levs,dim,phys_bc)

  ! define the name of the plotfile that will be written
  write(unit=parafile_name,fmt='("CAL_INFO")')

  if ( parallel_IOProcessor() ) then
      call fabio_mkdir(parafile_name)

      ! the calculate data to input into
      call fabio_mkdir("CAL_DATA")

      un3 = 829
      un4 = 839
      un1 = 849
      un2 = 819

      ! define the name of the plotfile that will be written
      write(unit=status_filename2,fmt='("calpara")' )
      write(unit=status_filename3,fmt='("calstat")' )
      write(unit=status_filename1,fmt='("tipoper")' )

      ! open a file name TK_INFO.dat, replace the file everytime write the tip information
      open(unit=un3, file = trim(parafile_name) // "/" // trim(status_filename2), &
              form = "formatted", access = "sequential", &
              status = "replace", action = "write")

      write(unit=un3, nml = probin)

      ! open a file name TK_INFO.dat, replace the file everytime write the tip information
      open(unit=un4, file = trim(parafile_name) // "/" // trim(status_filename3), &
              form = "formatted", access = "sequential", &
              status = "replace", action = "write")

      ! open a file name TK_INFO.dat, replace the file everytime write the tip information
      open(unit=un1, file = trim(parafile_name) // "/" // trim(status_filename1), &
              form = "formatted", access = "sequential", &
              status = "replace", action = "write")
      
      write(unit=status_filename1,fmt='("time_step")' )

      ! open a file name TK_INFO.dat, replace the file everytime write the tip information
      open(unit=un2, file = trim(parafile_name) // "/" // trim(status_filename1), &
              form = "formatted", access = "sequential", &
              status = "replace", action = "write")

  end if

  ! Decided by the computing mode, if it is 1 then the load mode starts
  if(compute_mode .ge. 1) then  ! the load mode

     if ( parallel_IOProcessor() ) then
         print*, '--> Loading data ... '
         write(un4, '("--> Loading data ...")')
     end if

     ! Read the ml_Boxarray mba from the input file
     call fabio_ml_boxarray_read(mba, time, load_folder_name)

     ! Read the fab data from the input file
     call fabio_ml_multifab_read_d2(mla,phi_load,mba,load_folder_name,is_periodic)

     ! get the levels
     nlevs = mla%nlevel

     ! set grid spacing at each level
     dx(1) = dx_def
     do n=2,max_levs
        dx(n) = dx(n-1) / mba%rr(n-1,1)
     end do

     if ( parallel_IOProcessor() ) then
         print*, '--> Data loading successfully ... '
         write(un4, '("--> Data loading successfully ...")')
     end if

  else

     ! tell mba about max_levs and dimensionality of problem
     call ml_boxarray_build_n(mba,max_levs,dim)

     ! tell mba about the ref_ratio between levels
     do n=2,max_levs
        mba%rr(n-1,:) = 2
     end do

     ! set grid spacing at each level
     dx(1) = dx_def 
     do n=2,max_levs
        dx(n) = dx(n-1) / mba%rr(n-1,1)
     end do

     ! create a box from (0,0) to (n_cell-1,n_cell-1)
     lo(:) = 0
     hi(1) = n_cell_x - 1
     hi(2) = n_cell_y - 1
     if(dim .eq. 3) then
        hi(3) = n_cell_z - 1
     end if

     ! calculation resumes
     bx = make_box(lo,hi)

     ! tell mba about the problem domain at each level
     mba%pd(1) = bx
     do n=2,max_levs
        mba%pd(n) = refine(mba%pd(n-1),mba%rr((n-1),:))
     end do

     ! initialize the boxarray at level 1 to be one single box
     call boxarray_build_bx(mba%bas(1),bx)

     ! overwrite the boxarray at level 1 to respect max_grid_size
     call boxarray_maxsize(mba%bas(1),max_grid_size)

     ! build the level 1 layout
     call layout_build_ba(la_array(1),mba%bas(1),mba%pd(1),is_periodic)

     ! build the level 1 multifab with 1 component and 1 ghost cell
     call multifab_build(ad_ori(1),la_array(1),1,1)
     call multifab_build(phi_old(1),la_array(1),3,1)
     call multifab_build(phi_new(1),la_array(1),3,1)
     call multifab_build(noise(1),la_array(1),1,1)

     ! For the flow routine
     if(do_solve_flow .eq. 1) then
        call multifab_build(flw(1),la_array(1),n_comp_flw,1)
     !else
     !   call multifab_build(flw(1),la_array(1),1,1)
     end if

     ! define level 1 of the_bc_tower
     call bc_tower_level_build(the_bc_tower,1,la_array(1))

     ! initialize phi on level 1
     call init_phi_on_level(phi_old(1),ad_ori(1),dx(1),seed_pos,prob_lo,prob_hi,&
                         the_bc_tower%bc_tower_array(1),pfpara)

     ! For the flow
     if(do_solve_flow .eq. 1) then
        call init_flw_on_level(flw(1),phi_old(1),the_bc_tower%bc_tower_array(1),pfpara)
     end if

     nl = 1
     new_grid = .true.

     do while ( (nl .lt. max_levs) .and. (new_grid) )

        call make_new_grids_init(new_grid,la_array(nl),la_array(nl+1),flw(nl),phi_old(nl),seed_pos,&
                              pfpara,dx(nl),amr_buf_width,mba%rr(nl,1),nl,max_grid_size,prob_lo,prob_hi)

        if (new_grid) then

           ! tell mba about the finer level boxarray
           call copy(mba%bas(nl+1),get_boxarray(la_array(nl+1)))

           ! Build the level nl+1 data
           call multifab_build(ad_ori(nl+1), la_array(nl+1),1,1)
           call multifab_build(phi_old(nl+1),la_array(nl+1),3,1)
           call multifab_build(phi_new(nl+1),la_array(nl+1),3,1)
           call multifab_build(noise(nl+1),  la_array(nl+1),1,1)

           if(do_solve_flow .eq. 1) then
             call multifab_build(flw(nl+1),la_array(nl+1),n_comp_flw,1)
           !else
           !  call multifab_build(flw(nl+1),la_array(nl+1),1,1)
           end if

           ! define level nl+1 of the_bc_tower
           call bc_tower_level_build(the_bc_tower,nl+1,la_array(nl+1))

           ! initialize phi on level nl+1
           call init_phi_on_level(phi_old(nl+1),ad_ori(nl+1),dx(nl+1),seed_pos,prob_lo,prob_hi, &
                               the_bc_tower%bc_tower_array(nl+1),pfpara)

           if(do_solve_flow .eq. 1) then
              call init_flw_on_level(flw(nl+1),phi_old(nl+1),the_bc_tower%bc_tower_array(nl+1),pfpara)
           end if

           nl = nl+1

        end if

     end do

     nlevs = nl

     ! destroy phi - we are going to build it again using the new multilevel
     do n=1,nlevs
        call destroy(ad_ori(n))
        call destroy(phi_old(n))
        call destroy(phi_new(n))
        call destroy(noise(n))
        if(do_solve_flow .eq. 1) then
          call destroy(flw(n))
        end if
     end do

     if (nlevs .ge. 3) then
        call enforce_proper_nesting(mba,la_array,max_grid_size)
     end if

     do n=1,nlevs
        call destroy(la_array(n))
     end do

     ! tell mla that there are nlevs levels, not max_levs
     call ml_layout_restricted_build(mla,mba,nlevs,is_periodic)

  end if ! compute_mode else

  ! this makes sure the boundary conditions are properly defined everywhere
  do n = 1,nlevs
     call bc_tower_level_build(the_bc_tower,n,mla%la(n))
  end do

  do n=1,nlevs
     call multifab_build(ad_ori(n),mla%la(n),1,1)
     call multifab_build(phi_old(n),mla%la(n),3,1)
     call multifab_build(phi_new(n),mla%la(n),3,1)
     call multifab_build(noise(n),mla%la(n),1,1)

     if(do_solve_flow .eq. 1) then
       call multifab_build(flw(n),mla%la(n),n_comp_flw,1)
       call multifab_build(flw_old(n),mla%la(n),n_comp_flw_uf,1)
       call multifab_build(flw_new(n),mla%la(n),n_comp_flw_uf,1)
       call multifab_build(surf(n),mla%la(n),1,1)
     !else
     !  call multifab_build(flw(n),mla%la(n),1,1)
     !  call multifab_build(flw_old(n),mla%la(n),1,1)
     !  call multifab_build(flw_new(n),mla%la(n),1,1)
     !  call multifab_build(surf(n),mla%la(n),1,1)
     end if

  end do

  ! setup the parameters for the LBM method
  call pf_para_set_lbm(pfpara,dim,dx(nlevs))

  call destroy(mba)

  if(compute_mode .ge. 1) then

     if ( parallel_IOProcessor() ) then
         print*, '--> Transfering and initializing loading data ... '
         write(un4, '("--> Transfering and initializing loading data ...")')
     end if

     if(do_solve_flow .eq. 1) then
       call trans_result_load_flw(mla,flw,phi_old,ad_ori,phi_load,the_bc_tower)
     else
       call trans_result_load(mla,phi_old,ad_ori,phi_load,dx,the_bc_tower,time,&
                            prob_lo,prob_hi,pfpara)
     end if

     do n=1,nlevs
        call destroy(phi_load(n))
     end do

  else

     call init_phi(mla,phi_old,ad_ori,dx,seed_pos,prob_lo,prob_hi,the_bc_tower,pfpara)

     ! copy the content from new multifab to old multifab
     do n=1,nlevs
         call multifab_copy(phi_new(n), phi_old(n), 1)
     end do

     if(do_solve_flow .eq. 1) then
        call init_flw(mla,flw,phi_old,the_bc_tower,pfpara)
     end if

     istep = 0
     time = 0.d0

  end if

  total_flw_time = 0

  dt01 = dt01 * dx(max_levs)**2 /(2.d0*dim*DM)
  dt08 = dt08 * dx(max_levs)**2 /(2.d0*dim*DM)

  if(compute_mode .eq. 0) then

    dt = dt01
    pfpara%period_p0 = 0.d0
    elapsed_time = parallel_wtime()

    if(pfpara%solve_flow_only /= 1) then

      if ( parallel_IOProcessor() ) then
        print *, '--> Pre-calculate the data for seed evolution'
      end if

      do istep=1,nsteps_0

        ! generate random numbers for uncertainty
        call gen_random(mla,noise,pfpara,istep,the_bc_tower)

        ! advance phase field
        call advance_phase_field(mla,ad_ori,phi_old,phi_new,dx,dt,the_bc_tower,pfpara,noise)

        ! For the flow routine
        start_flw_time = parallel_wtime()

        !call advance_flow_field(mla,flw,flw_old, flw_new, phi_new,pfpara,the_bc_tower,&
        !                    dt,istep,dx,prob_lo,prob_hi,un4)

        total_flw_time = total_flw_time + parallel_wtime()-start_flw_time

        ! need to change here when flow is considered
        call advance_solute_field(mla,flw,phi_old,phi_new,dx,dt,the_bc_tower,pfpara)

        ! advance temperature field
        call advance_temp_field(mla,phi_old,phi_new,dx,dt,reset_plotnum,the_bc_tower,pfpara,time, &
                               temp_l,temp_h_l,prob_lo,prob_hi,istep,0)

        ! copy the content from new multifab to old multifab
        do n=1,nlevs
           call multifab_copy(phi_old(n), phi_new(n), 1)
        end do

        if ( parallel_IOProcessor() ) then
            ! run time counting 
            run_time = parallel_wtime() - elapsed_time

            runtime_hours = int(run_time / 3600.0)
            runtime_minutes = int((run_time - real(runtime_hours) * 3600.0) / 60.0)
            runtime_seconds = int(run_time - real(runtime_hours) * 3600.0 - real(runtime_minutes) * 60.0)

            write(time_str, '(i4,a,i2,a,i2)') runtime_hours, ':', runtime_minutes, ':', runtime_seconds

            ! i_pro = int(100.0*real(istep)/real(nsteps_0))
            ! write(vchar,'(i3)') n_pro  ! n control the number of backspace char
            ! myformat='('//trim(adjustl(vchar))//'a,i10,i10,a,a15)'
            ! write(*,myformat,advance='no') (char(8),j_pro=1,n_pro), istep, i_pro, '%',time_str
            print*, 'Pre-simulation steps: ', istep, ' , elapsed time: ', time_str

        end if

        time = time + dt

      end do

      if ( parallel_IOProcessor() ) then
        print *, char(13)
      end if

    end if

    call write_plotfile(mla,flw,ad_ori,phi_old,surf,0,dx,time,prob_lo,prob_hi,pfpara)

  end if

  if( compute_mode .le. 1) then

     regrid_time = 0
     v_mag_max  = 0.d0
     elapsed_time = parallel_wtime()
     need_regrid = .false.

     call parallel_barrier()

     ada_regriding_count = 0
     dt = dt08
     pfpara%period_p0 = period_p0

     if ( parallel_IOProcessor() ) then
        print *, '--> Normal calculate the data for dendrite growth'
      end if

     do istep=1,nsteps

        if(ada_regriding .eq. 1) then
            ada_regriding_count = ada_regriding_count + 1

            need_regrid = (istep > 1 .and. max_levs > 1 .and. regrid_int > 0 .and. &
                           (mod(ada_regriding_count,regrid_int) .eq. 0) )
        else
            need_regrid = (istep > 1 .and. max_levs > 1 .and. regrid_int > 0 .and. &
                           (mod(istep-1,regrid_int) .eq. 0) )
        end if

        if(compute_mode .eq. 1) then
            need_regrid = need_regrid .or. (istep .eq. 1)
        end if

        if ( need_regrid ) then

           start_regrid_time = parallel_wtime()

           call regrid(mla,flw,flw_old,flw_new,surf,phi_old,ad_ori,phi_new,noise,nlevs,max_levs,dx, &
                       the_bc_tower,amr_buf_width,max_grid_size,pfpara,n_comp_flw_uf)

           if(ada_regriding .eq. 1 ) then

              ! reset the regriding count
              ada_regriding_count = 0

              ! recalculate the number regrid_int
              call recalculate_regrid_number(mla,regrid_int,v_mag_max,dt,dx,regrid_amp_ref, &
                                          regrid_int_max,regrid_int_min)

            end if

           regrid_time = regrid_time + parallel_wtime() - start_regrid_time

           ! reset the magnitude of the velocity
           v_mag_max = 0.d0

        end if

        !do iter_pf=1,n_calpf

        ! generate random numbers for uncertainty
        call gen_random(mla,noise,pfpara,istep,the_bc_tower)

        ! advance phase field
        call advance_phase_field(mla,ad_ori,phi_old,phi_new,dx,dt,the_bc_tower,pfpara,noise)

        !if(n_calpf == 1) then
          
        ! For the flow routine
        start_flw_time = parallel_wtime()     

        call advance_flow_field(mla,flw,flw_old, flw_new,surf,phi_new,pfpara,the_bc_tower,&
                            dt,istep,dx,prob_lo,prob_hi,need_regrid,time,un4)

        total_flw_time = total_flw_time + parallel_wtime()-start_flw_time

        ! this subroutine specifies whether the flow is too fast, if so the phase field must be
        ! proceeded in multiple steps
        !call recalculate_time_step(mla,dt,n_calpf,dt08,flw,0.9d0,dx(max_levs),dx(1),istep,time,un2)

        !end if

        ! change here for the advection induced by solute or temperature gradients
        call advance_solute_field(mla,flw,phi_old,phi_new,dx,dt,the_bc_tower,pfpara)

        ! advance temperature field
        call advance_temp_field(mla,phi_old,phi_new,dx,dt,reset_plotnum,the_bc_tower,pfpara,time, &
                               temp_l,temp_h_l,prob_lo,prob_hi,istep,1)

        if(ada_regriding .eq. 1 ) then

          ! Retrieve the maximum magnitude of the velocity
          call retrieve_max_velocity(istep,v_mag_max,mla,phi_old,phi_new,dt,dx,thers_gradient)

        end if

        ! copy the content from new multifab to old multifab
        do n=1,nlevs
            call multifab_copy(phi_old(n), phi_new(n), 1)
        end do

        time = time + dt

        !end do

        if(write_file_mode .eq. 1) then

           ! measure up the information
           call follow_tip_info(un1,time,mla,phi_new,dx)

        end if

        if ( mod(istep,plot_int) .eq. 0 .or. (istep .eq. nsteps) ) then
           call write_plotfile(mla,flw,ad_ori,phi_old,surf,istep+reset_plotnum,dx,time,prob_lo,prob_hi,pfpara)
        end if

        if ( (pfpara%do_solve_flow == 1) .and. ( pfpara%plot_mode .eq. 0) .and. ( mod(istep,plot_int_1) .eq. 0 ) ) then
           pfpara%plot_mode = 1
           call write_plotfile(mla,flw,ad_ori,phi_old,surf,istep+reset_plotnum+1,dx,time,prob_lo,prob_hi,pfpara)
           pfpara%plot_mode = 0
        end if

        ! Write out the calculation information, mainly the computatuon time
        if ( mod(istep,status_int) .eq. 0) then
           if ( parallel_IOProcessor() ) then
               !print*,istep, time, parallel_wtime()-elapsed_time, regrid_int, regrid_time, total_flw_time
               write(un4, '(i8.7, es10.2e3, es10.2e3, i8.7, es10.2e3, es10.2e3)') &
                           istep, time, parallel_wtime()-elapsed_time, regrid_int, regrid_time, total_flw_time
               
            end if
        end if

        if ( parallel_IOProcessor() ) then
            ! run time counting 
            run_time = parallel_wtime() - elapsed_time

            runtime_hours = int(run_time / 3600.0)
            runtime_minutes = int((run_time - real(runtime_hours) * 3600.0) / 60.0)
            runtime_seconds = int(run_time - real(runtime_hours) * 3600.0 - real(runtime_minutes) * 60.0)

            write(time_str, '(i4,a,i2,a,i2)') runtime_hours, ':', runtime_minutes, ':', runtime_seconds

            ! i_pro = int(100.0*real(istep)/real(nsteps))
            ! write(vchar,'(i3)') n_pro  ! n control the number of backspace char
            ! myformat='('//trim(adjustl(vchar))//'a,i10,i10,a,a15)'
            ! write(*,myformat,advance='no') (char(8),j_pro=1,n_pro), istep, i_pro, '%',time_str
            print*, 'Normal-simulation steps: ', istep, ' , elapsed time: ', time_str

        end if
        
        ! if(pfpara%do_solve_flow == 1) then
        !   call retrieve_error_msg(mla,flw,pfpara%flw_constrain_limit,err_msg)
        !   if(err_msg) exit;
        ! end if

     end do

     if ( parallel_IOProcessor() ) then
        close(un3)
        close(un4)
        close(un1)
        close(un2)
     end if

     ! modify the writing mode
     if(pfpara%do_solve_flow == 1 .and. err_msg ) then
       pfpara%plot_mode = 1  
       call write_plotfile(mla,flw,ad_ori,phi_old,surf,istep+reset_plotnum,dx,time,prob_lo,prob_hi,pfpara)
     end if

  end if 

  !if ( parallel_IOProcessor() ) print*, 'step 1';

  do n=1,nlevs
     call destroy(ad_ori(n))
     call destroy(phi_old(n))
     call destroy(phi_new(n))
     call destroy(noise(n))
     
     if(pfpara%do_solve_flow == 1) then
       call destroy(flw(n))
       call destroy(flw_old(n))
       call destroy(flw_new(n))
       call destroy(surf(n))
     end if
  end do

  !if ( parallel_IOProcessor() ) print*, 'step 2';

  call destroy(mla)
  call bc_tower_destroy(the_bc_tower)
  call pf_para_destroy(pfpara)

  !if ( parallel_IOProcessor() ) print*, 'step 3';

  deallocate(dx,lo,hi,is_periodic,prob_lo,prob_hi)
  if(seed_num .ge. 1) deallocate(seed_pos)

  ! deallocate temporary boxarrays and communication mappings
  call layout_flush_copyassoc_cache ()

  ! check for memory that should have been deallocated
  if ( parallel_IOProcessor() ) then
     print *, char(13)
     print*, 'MEMORY STATS AT END OF PROGRAM'
     print*, ' '
  end if

  call print(multifab_mem_stats(),    "    multifab")
  call print(fab_mem_stats(),         "         fab")
  call print(boxarray_mem_stats(),    "    boxarray")
  call print(layout_mem_stats(),      "      layout")
  call print(boxassoc_mem_stats(),    "    boxassoc")
  call print(fgassoc_mem_stats(),     "     fgassoc")
  call print(syncassoc_mem_stats(),   "   syncassoc")
  call print(copyassoc_mem_stats(),   "   copyassoc")
  call print(fluxassoc_mem_stats(),   "   fluxassoc")

  run_time = parallel_wtime() - start_time

  ! collect run_time from each processor and store the maximum
  call parallel_reduce(run_time_IOproc, run_time, MPI_MAX, &
                       proc = parallel_IOProcessorNode())

  if ( parallel_IOProcessor() ) then
     print*,"Run time (s) =",run_time_IOproc
  end if

  call boxlib_finalize()

end program main
