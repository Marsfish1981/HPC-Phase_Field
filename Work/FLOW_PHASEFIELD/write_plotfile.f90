module write_plotfile_module

  use ml_layout_module
  use multifab_module
  use fabio_module
  use pf_utility_module

  implicit none

contains

  subroutine write_plotfile(mla,flw,ad,phi,sf,istep,dx,time,prob_lo,prob_hi,pfpara,phi_new,dt)

    type(ml_layout), intent(in   )           :: mla
    type(multifab) , intent(in   )           :: phi(:),ad(:),flw(:),sf(:)
    integer        , intent(in   )           :: istep
    real(dp_t)     , intent(in   )           :: dx(:),time
    real(dp_t)     , intent(in   )           :: prob_lo(mla%dim), prob_hi(mla%dim)
    type(pf_para),   intent(in   )           :: pfpara
    type(multifab) , intent(in   ), optional :: phi_new(:)
    real(dp_t)     , intent(in   ), optional :: dt

    select case(pfpara%do_solve_flow)
    case(0)
      if(present(phi_new) .and. present(dt)) then
        call write_plotfile_den(mla,ad,phi,istep,dx,time,prob_lo,prob_hi,pfpara,phi_new,dt)
      else
        call write_plotfile_den(mla,ad,phi,istep,dx,time,prob_lo,prob_hi,pfpara)
      end if
    case(1)
      call write_plotfile_flw(mla,flw,ad,phi,sf,istep,dx,time,prob_lo,prob_hi,pfpara)
    end select

  end subroutine write_plotfile

  ! Define plot mode
  !plot_mode:
  !           0 = AD,PF,UC,TH
  !           1 = AD,PF,UC,TH,VM,GRA
  !           2 = PF,TRC
  !           3 = AD,PF,TRC
  !           4 = PF,GRA
  !           5 = AD,PF

  subroutine write_plotfile_den(mla,ad,phi,istep,dx,time,prob_lo,prob_hi,pfpara,phi_new,dt)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: phi(:)
    type(multifab) , intent(in   ) :: ad(:)
    integer        , intent(in   ) :: istep
    real(dp_t)     , intent(in   ) :: dx(:),time
    real(dp_t)     , intent(in   ) :: prob_lo(mla%dim), prob_hi(mla%dim)
    type(pf_para),   intent(in   ) :: pfpara
    type(multifab) , intent(in   ),optional :: phi_new(:)
    real(dp_t)     , intent(in   ),optional :: dt

    ! local variables
    integer           :: n, nlevs, nvarl = 4
    character(len=40)  :: plotfile_name
    character(len=20),allocatable :: variable_names(:)
    integer :: un

    type(multifab), allocatable :: plotdata(:)

    ! dimensioned as an array with size dm for fabio_ml_multifab_write_d
    real(dp_t) :: dx_vec(mla%dim)

    if(pfpara%plot_mode .eq. 1) then
        nvarl = 6
    else if(pfpara%plot_mode .eq. 2) then
        nvarl = 2
    else if(pfpara%plot_mode .eq. 3) then
        nvarl = 3
    else if(pfpara%plot_mode .eq. 4) then
        nvarl = 2
    else if(pfpara%plot_mode .eq. 5) then
        nvarl = 2
    else if(pfpara%plot_mode .eq. 6) then
        nvarl = 1
    else if(pfpara%plot_mode .eq. 7) then
        nvarl = 1
    end if

    allocate( variable_names(nvarl) )

    if(pfpara%plot_mode .eq. 2) then
       variable_names(1) = "phase_field"
       variable_names(2) = "Solute_trc"
    else if(pfpara%plot_mode .eq. 3) then
       variable_names(1) = "orient_field"
       variable_names(2) = "phase_field"
       variable_names(3) = "Solute_trc"
    else if(pfpara%plot_mode .eq. 4) then
       variable_names(1) = "orient_field"
       variable_names(2) = "phase_field"
    else if(pfpara%plot_mode .eq. 5) then
       variable_names(1) = "parallel_dis"
       variable_names(2) = "phase_field"
    else if(pfpara%plot_mode .eq. 6) then
       variable_names(1) = "phase_field"
    else if(pfpara%plot_mode .eq. 7) then
       variable_names(1) = "phase_field"
    else
       variable_names(1) = "orientation"
       variable_names(2) = "phase_field"
       variable_names(3) = "solute_field"
       variable_names(4) = "temp_field"

       if(pfpara%plot_mode .eq. 1) then
          variable_names(5) = "velocity_mag"
          variable_names(6) = "gra_criteor"
       end if

    end if

    nlevs = mla%nlevel

    dx_vec(:) = dx(1)

    allocate(plotdata(nlevs))

    do n=1,nlevs
       ! build plotdata with 1 component and 0 ghost cells
       call multifab_build(plotdata(n),mla%la(n),nvarl,0)

       if(pfpara%plot_mode .eq. 2) then
          call multifab_copy_c(plotdata(n),1,phi(n),1,1)

          ! calculate the trc term
          call cal_phi_tag(plotdata(n), 2,ad(n),phi(n), dx(n),pfpara%para_uc,pfpara%para_th,pfpara%sk,0)

       else if(pfpara%plot_mode .eq. 3) then

          call multifab_copy_c(plotdata(n),1,ad(n),1,1)
          call multifab_copy_c(plotdata(n),2,phi(n),1,1)

          ! calculate the trc term
          call cal_phi_tag(plotdata(n), 3,ad(n),phi(n), dx(n),pfpara%para_uc,pfpara%para_th,pfpara%sk,0)

          !call multifab_copy_c(plotdata(n),3,phi(n),3,1)

       else if(pfpara%plot_mode .eq. 4) then
          call multifab_copy_c(plotdata(n),1,ad(n),1,1)
          call multifab_copy_c(plotdata(n),2,phi(n),1,1)

          ! calculate the trc term
          !call cal_phi_tag(plotdata(n), 2,ad(n),phi(n), dx(n),pfpara%para_uc,pfpara%para_th,pfpara%sk,1)

       else if(pfpara%plot_mode .eq. 5) then
          !call multifab_copy_c(plotdata(n),1,ad(n),1,1)
          call multifab_setval_c(plotdata(n),dble(parallel_myproc()),1,1)
          call multifab_copy_c(plotdata(n),2,phi(n),1,1)

       else if(pfpara%plot_mode .eq. 6) then
          call multifab_copy_c(plotdata(n),1,phi(n),1,1)

          ! calculate the trc term
          !call cal_phi_tag(plotdata(n), 2,ad(n),phi(n), dx(n),pfpara%para_uc,pfpara%para_th,pfpara%sk,2)
       else if(pfpara%plot_mode .eq. 7) then
          !call multifab_copy_c(plotdata(n),1,ad(n),1,1)
          !call multifab_copy_c(plotdata(n),2,phi(n),1,1)

          ! calculate the trc term
          !call cal_phi_tag(plotdata(n), 3,ad(n),phi(n), dx(n),pfpara%para_uc,pfpara%para_th,pfpara%sk,1,phi_new(n),dt)
          call cal_phi_tag_test(plotdata(n), 1, dx(n),prob_lo,prob_hi)

       else
          ! copy the state into plotdata
          call multifab_copy_c(plotdata(n),1,ad(n),1,1)
          call multifab_copy_c(plotdata(n),2,phi(n),1,3)

          ! calculate the gradient stuff in component 5
          if(pfpara%plot_mode .eq. 1) then
              ! calculate the trc term
              call cal_phi_tag(plotdata(n), 5,ad(n),phi(n), dx(n),pfpara%para_uc,pfpara%para_th,pfpara%sk,1,phi_new(n),dt)
              call cal_phi_tag(plotdata(n), 6,ad(n),phi(n), dx(n),pfpara%para_uc,pfpara%para_th,pfpara%sk,1)
          end if
       end if

    end do

    ! define the name of the plotfile that will be written
    write(unit=plotfile_name,fmt='("CAL_DATA/sim",i7.7)') istep

    ! write the plotfile
    call fabio_ml_multifab_write_d(plotdata, mla%mba%rr(:,1), plotfile_name, variable_names, &
                                   mla%mba%pd(1), prob_lo, prob_hi, time, dx_vec)

    ! make sure to destroy the multifab or you'll leak memory
    do n=1,nlevs
       call multifab_destroy(plotdata(n))
    end do

  end subroutine write_plotfile_den

  !------------------------------FLOW--------------------------------------------------------------
  !plot_mode:
  !               0 = AD,PF,SF,UC,TH,UX,UY,UZ
  !               1 = AD,PF,UC,TH,UX,UY,UZ,RHO,F1-F9(2D),F1-F19(3D)
  subroutine write_plotfile_flw(mla,flw,ad,phi,sf,istep,dx,time,prob_lo,prob_hi,pfpara)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: flw(:),phi(:),ad(:),sf(:)
    integer        , intent(in   ) :: istep
    real(dp_t)     , intent(in   ) :: dx(:),time
    real(dp_t)     , intent(in   ) :: prob_lo(mla%dim), prob_hi(mla%dim)
    type(pf_para),   intent(in   ) :: pfpara

    ! local variables
    integer           :: n, nlevs, nvarl,i,iter_var,max_iter
    character(len=40)  :: plotfile_name
    character(len=20),allocatable :: variable_names(:)
    character(len=20)  :: f_names(19)
    !character(len=20)  :: k_names(19)
    integer :: un

    type(multifab), allocatable :: plotdata(:)

    ! dimensioned as an array with size dm for fabio_ml_multifab_write_d
    real(dp_t) :: dx_vec(mla%dim)

    if(mla%dim .eq. 2) then
        nvarl=6
        max_iter = 9
    else if(mla%dim .eq. 3) then
        nvarl = 7
        max_iter = 19
    end if

    if(pfpara%plot_mode .eq. 1) then
        if(mla%dim .eq. 2) then
           nvarl = 16
        else if(mla%dim .eq. 3) then
           nvarl = 27
        end if
        do i=1,19
            write(unit=f_names(i),fmt='("f_",i2.2)') i
        end do
    end if

    allocate( variable_names(nvarl) )

    variable_names(1) = "orientation"
    variable_names(2) = "phase_field"
    variable_names(3) = "solute_field"
    variable_names(4) = "temp_field"
    variable_names(5) = "ux"
    variable_names(6) = "uy"

    if(mla%dim .eq. 2) then
        iter_var          = 6
    else if(mla%dim .eq. 3) then
        variable_names(7) = "uz"
        iter_var          = 7
    end if

    if(pfpara%plot_mode .eq. 1) then
        iter_var = iter_var + 1
        variable_names(iter_var) = "rho"
        do i=1,max_iter
            iter_var = iter_var + 1
            variable_names(iter_var) = f_names(i)
        end do
    end if

    nlevs = mla%nlevel
    dx_vec(:) = dx(1)
    allocate(plotdata(nlevs))

    do n=1,nlevs
       ! build plotdata with 1 component and 0 ghost cells
       call multifab_build(plotdata(n),mla%la(n),nvarl,0)

       call multifab_copy_c(plotdata(n),1,ad(n),1,1)
       call multifab_copy_c(plotdata(n),2,phi(n),1,3)
       !call multifab_copy_c(plotdata(n),2,phi(n),1,2)
       !call multifab_copy_c(plotdata(n),4,sf(n),1,1)
       call multifab_copy_c(plotdata(n),5,flw(n),1,2)

       if(mla%dim .eq. 2) then
          iter_var = 6
       else if(mla%dim .eq. 3) then
         call multifab_copy_c(plotdata(n),7,flw(n),3,1)
         iter_var  = 7
       end if

       if(pfpara%plot_mode .eq. 1) then
         iter_var = iter_var + 1
         if(mla%dim .eq. 2) then
            call multifab_copy_c(plotdata(n),iter_var,flw(n),3,max_iter+1)
            !iter_var = iter_var + 1
            !call multifab_copy_c(plotdata(n),iter_var,flw(n),4,max_iter)
         else if(mla%dim .eq. 3) then
            call multifab_copy_c(plotdata(n),iter_var,flw(n),4,max_iter+1)
            !iter_var = iter_var + 1
            !call multifab_copy_c(plotdata(n),iter_var,flw(n),5,max_iter)
         end if
       end if
    end do

    ! define the name of the plotfile that will be written
    write(unit=plotfile_name,fmt='("CAL_DATA/sim",i7.7)') istep

    ! write the plotfile
    call fabio_ml_multifab_write_d(plotdata, mla%mba%rr(:,1), plotfile_name, variable_names, &
                                   mla%mba%pd(1), prob_lo, prob_hi, time, dx_vec)

    ! make sure to destroy the multifab or you'll leak memory
    do n=1,nlevs
       call multifab_destroy(plotdata(n))
    end do

  end subroutine write_plotfile_flw
  !------------------------------------------------------------------------------------------------

  subroutine cal_phi_tag(phi_des, cds, ad,phi, dx,para_uc,para_th,s_k,cal_kind,phi_new,dt)

    type(multifab) , intent(inout) :: phi_des
    type(multifab) , intent(in) :: phi,ad
    real(kind=dp_t), intent(in   ) :: dx,para_uc,para_th,s_k
    integer,         intent(in   ) :: cds,cal_kind
    type(multifab) , intent(in),optional :: phi_new
    real(kind=dp_t), intent(in),optional :: dt

    integer :: lo(phi%dim), hi(phi%dim)
    integer :: dm, ng_p, i

    real(kind=dp_t), pointer ::  pfd(:,:,:,:)
    real(kind=dp_t), pointer ::  pfo(:,:,:,:)
    real(kind=dp_t), pointer ::  pfa(:,:,:,:)
    real(kind=dp_t), pointer ::  pfn(:,:,:,:)


    dm   = phi%dim
    ng_p = phi%ng

     do i=1,nfabs(phi)
       pfd  => dataptr(phi_des,i)
       pfo  => dataptr(phi,i)
       pfa  => dataptr(ad,i)

       if (present(phi_new))  then
          pfn  => dataptr(phi_new,i)
       end if

       lo = lwb(get_box(phi,i))
       hi = upb(get_box(phi,i))

       select case(dm)
       case (2)
          if (present(phi_new))  then
             ! calculate the K1
             call cal_phi_gra_2d_new(pfd(:,:,1,cds),pfo(:,:,1,1),&
                              ng_p,dx,dt,lo,hi,cal_kind,pfn(:,:,1,1))
          else
            ! calculate the K1
             call cal_phi_gra_2d(pfd(:,:,1,cds),pfa(:,:,1,1),pfo(:,:,1,1),pfo(:,:,1,2),pfo(:,:,1,3),&
                              ng_p,dx,lo,hi,para_uc,para_th,s_k,cal_kind)
          end if


       case (3)
          if (present(phi_new))  then
             call cal_phi_gra_3d_new(pfd(:,:,:,cds),pfo(:,:,:,1),&
                              ng_p,dx,dt,lo,hi,cal_kind,pfn(:,:,:,1))
          else
             call cal_phi_gra_3d(pfd(:,:,:,cds),pfa(:,:,:,1),pfo(:,:,:,1),pfo(:,:,:,2),pfo(:,:,:,3),&
                              ng_p,dx,lo,hi,para_uc,para_th,s_k,cal_kind)
          end if
       end select

    end do

  end subroutine cal_phi_tag

  subroutine cal_phi_gra_2d(pfd,ad,pf,uc,th,ng_p,dx,lo,hi,para_uc,para_th,s_k,cal_kind)

    integer           :: lo(2),hi(2),ng_p,cal_kind
    double precision  :: pfd(lo(1):hi(1),lo(2):hi(2))
    double precision  :: ad(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
    double precision  :: pf(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
    double precision  :: uc(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
    double precision  :: th(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
    double precision  :: dx,para_uc,para_th,s_k

    ! local variables
    integer            :: i,j
    real(kind = dp_t)  :: var_e,gra_x,gra_y
    double precision   :: var_target, mod_1
    double precision   :: de1_x, de1_y, de2_xx, de2_yy, de2_xy,thers_gradient=0.05

    select case (cal_kind)
    case (1)
        ! do the loops
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             gra_x = ( pf(i+1,j) - pf(i-1,j) ) / (2.d0 * dx)
             gra_y = ( pf(i,j+1) - pf(i,j-1) ) / (2.d0 * dx)
             var_e = sqrt( gra_x**2 + gra_y**2 )

             gra_x = ( uc(i+1,j) - uc(i-1,j) ) / (2.d0 * dx)
             gra_y = ( uc(i,j+1) - uc(i,j-1) ) / (2.d0 * dx)
             var_e = var_e + para_uc * sqrt( gra_x**2 + gra_y**2 )

             gra_x = ( th(i+1,j) - th(i-1,j) ) / (2.d0 * dx)
             gra_y = ( th(i,j+1) - th(i,j-1) ) / (2.d0 * dx)
             var_e = var_e + para_th * sqrt( gra_x**2 + gra_y**2 )

             !pfd(i,j) = var_e * 2.d0 * dx * dx
             pfd(i,j) = var_e * dx

         end do
      end do

    case (2)
        ! do the loops
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             ! Select those with PF >= 0 and set the value to
             ! the AD, which provides a sharp interface with only solid displaying
             if( pf(i,j) .ge. -0.01d0 ) then
                pfd(i,j) = ad(i,j)
             else
                pfd(i,j) = -1.d0
             end if

         end do
      end do

    case (3)
        ! do the loops
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             pfd(i,j) = 0.d0

             de1_x  = (pf(i+1,j) - pf(i-1,j))/(2.d0 * dx)
             de1_y  = (pf(i,j+1) - pf(i,j-1))/(2.d0 * dx)

             mod_1 = sqrt(de1_x**2 + de1_y**2)
             var_target = mod_1*dx

             if (var_target > thers_gradient ) then

                de2_xx = (pf(i+1,j) - 2.d0*pf(i,j) + pf(i-1,j))/(dx * dx)
                de2_yy = (pf(i,j+1) - 2.d0*pf(i,j) + pf(i,j-1))/(dx * dx)
                de2_xy = (pf(i+1,j+1) + pf(i-1,j-1) - pf(i+1,j-1) - pf(i-1,j+1) )/(dx * dx)

                pfd(i,j) = de2_xx * de1_y**2 + de2_yy * de1_x**2 - 2.d0*de1_x*de1_y*de2_xy
                pfd(i,j) = -0.5*pfd(i,j) / mod_1**3
             end if

         end do
      end do

    case default
       ! do the loops
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             pfd(i,j) = ( (1.d0 - s_k)*uc(i,j) + 1.d0 )*( 1.d0 + s_k - (1.d0 - s_k)*pf(i,j) )
             pfd(i,j) = pfd(i,j) / (2.d0 * s_k)
         end do
      end do

    end select

  end subroutine cal_phi_gra_2d

  subroutine cal_phi_gra_3d(pfd,ad,pf,uc,th,ng_p,dx,lo,hi,para_uc,para_th,s_k,cal_kind)

    integer           :: lo(3),hi(3),ng_p,cal_kind
    double precision  :: pfd(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    double precision  :: pf(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
    double precision  :: ad(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
    double precision  :: uc(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
    double precision  :: th(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
    double precision  :: dx,para_uc,para_th,s_k

    ! local variables
    integer            :: i,j,k
    real(kind = dp_t)  :: var_e,gra_x,gra_y,gra_z
    double precision  :: var_target, mod_1,thers_gradient=0.05
    double precision  :: de1_x,de1_y,de1_z,de2_xx,de2_yy,de2_zz,de2_xy,de2_yz,de2_xz

    select case (cal_kind)
    case (1)
       ! do the loops
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                gra_x = ( pf(i+1,j,k) - pf(i-1,j,k) ) / (2.d0 * dx)
                gra_y = ( pf(i,j+1,k) - pf(i,j-1,k) ) / (2.d0 * dx)
                gra_z = ( pf(i,j,k+1) - pf(i,j,k-1) ) / (2.d0 * dx)
                var_e = sqrt( gra_x**2 + gra_y**2 + gra_z**2 )

                gra_x = ( uc(i+1,j,k) - uc(i-1,j,k) ) / (2.d0 * dx)
                gra_y = ( uc(i,j+1,k) - uc(i,j-1,k) ) / (2.d0 * dx)
                gra_z = ( uc(i,j,k+1) - uc(i,j,k-1) ) / (2.d0 * dx)
                var_e = var_e + para_uc * sqrt( gra_x**2 + gra_y**2 + gra_z**2 )

                gra_x = ( th(i+1,j,k) - th(i-1,j,k) ) / (2.d0 * dx)
                gra_y = ( th(i,j+1,k) - th(i,j-1,k) ) / (2.d0 * dx)
                gra_z = ( th(i,j,k+1) - th(i,j,k-1) ) / (2.d0 * dx)
                var_e = var_e + para_th * sqrt( gra_x**2 + gra_y**2 + gra_z**2 )

                !pfd(i,j,k) = var_e * 3.d0 * dx * dx
                pfd(i,j,k) = var_e * dx

            end do
          end do
       end do

    case (2)
        ! do the loops
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)

             ! Select those with PF >= 0 and set the value to
             ! the AD, which provides a sharp interface with only solid displaying
             if( pf(i,j,k) .ge. -0.01d0 ) then
                pfd(i,j,k) = ad(i,j,k)
             else
                pfd(i,j,k) = -1.d0
             end if

             end do
          end do
       end do

    case (3)
       ! do the loops
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                pfd(i,j,k) = 0.d0

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

                   pfd(i,j,k) = de2_xx * (de1_y**2 + de1_z**2) + de2_yy * (de1_x**2 + de1_z**2) + de2_zz * (de1_x**2 + de1_y**2)
                   pfd(i,j,k) = pfd(i,j,k) - 2.d0*(de1_x*de1_y*de2_xy + de1_x*de1_z*de2_xz + de1_y*de1_z*de2_yz)
                   pfd(i,j,k) = -0.5*pfd(i,j,k) / mod_1**3

               end if

            end do
          end do
       end do

    case default
       ! do the loops
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                pfd(i,j,k) = ( (1.d0 - s_k)*uc(i,j,k) + 1.d0 )*( 1.d0 + s_k - (1.d0 - s_k)*pf(i,j,k) )
                pfd(i,j,k) = pfd(i,j,k) / (2.d0 * s_k)
             end do
          end do
       end do

    end select

  end subroutine cal_phi_gra_3d

  subroutine cal_phi_gra_2d_new(pfd,pf,ng_p,dx,dt,lo,hi,cal_kind,pf_new)

    integer           :: lo(2),hi(2),ng_p,cal_kind
    double precision  :: pfd(lo(1):hi(1),lo(2):hi(2))
    double precision  :: pf(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)
    double precision  :: dx,dt
    double precision  :: pf_new(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p)

    ! local variables
    integer            :: i,j
    real(kind = dp_t)  :: gra_x,gra_y,mod_1,thers=1.0E-5,de1_t,v_mag

    select case (cal_kind)
    case (1)
        ! do the loops
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             pfd(i,j) = 0.d0

             gra_x = ( pf(i+1,j) - pf(i-1,j) ) / (2.d0 * dx)
             gra_y = ( pf(i,j+1) - pf(i,j-1) ) / (2.d0 * dx)

             mod_1 = sqrt(gra_x**2 + gra_y**2)

             if (mod_1 > thers) then

                de1_t = (pf_new(i,j) - pf(i,j)) / dt
                v_mag = abs(de1_t)/mod_1

                pfd(i,j) = v_mag
             end if

         end do
      end do

    end select

  end subroutine cal_phi_gra_2d_new

  subroutine cal_phi_gra_3d_new(pfd,pf,ng_p,dx,dt,lo,hi,cal_kind,pf_new)

    integer           :: lo(3),hi(3),ng_p,cal_kind
    double precision  :: pfd(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    double precision  :: pf(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)
    double precision  :: dx,dt
    double precision  :: pf_new(lo(1)-ng_p:hi(1)+ng_p,lo(2)-ng_p:hi(2)+ng_p,lo(3)-ng_p:hi(3)+ng_p)

    ! local variables
    integer            :: i,j,k
    real(kind = dp_t)  :: gra_x,gra_y,gra_z,mod_1,thers=1.0E-5,de1_t,v_mag

    select case (cal_kind)
    case (1)
       ! do the loops
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                gra_x = ( pf(i+1,j,k) - pf(i-1,j,k) ) / (2.d0 * dx)
                gra_y = ( pf(i,j+1,k) - pf(i,j-1,k) ) / (2.d0 * dx)
                gra_z = ( pf(i,j,k+1) - pf(i,j,k-1) ) / (2.d0 * dx)

                mod_1 = sqrt(gra_x**2 + gra_y**2 + gra_z**2)

                if (mod_1 > thers) then

                   de1_t = (pf_new(i,j,k) - pf(i,j,k)) / dt
                   v_mag = abs(de1_t)/mod_1

                   pfd(i,j,k) = v_mag
                end if

            end do
          end do
       end do

    end select

  end subroutine cal_phi_gra_3d_new

  subroutine cal_phi_tag_test(phi, cds, dx,prob_lo,prob_hi)

    type(multifab) , intent(inout) :: phi
    real(kind=dp_t), intent(in   ) :: dx
    integer,         intent(in   ) :: cds
    real(kind=dp_t), intent(in   ) :: prob_lo(phi%dim)
    real(kind=dp_t), intent(in   ) :: prob_hi(phi%dim)

    integer :: lo(phi%dim), hi(phi%dim)
    integer :: dm, i

    real(kind=dp_t), pointer ::  pf(:,:,:,:)

    dm   = phi%dim

     do i=1,nfabs(phi)
       pf  => dataptr(phi,i)
       lo = lwb(get_box(phi,i))
       hi = upb(get_box(phi,i))

       select case(dm)
       case (2)
          call cal_phi_test_2d(pf(:,:,1,cds),dx,lo,hi,prob_lo,prob_hi)

       case (3)
          call cal_phi_test_3d(pf(:,:,:,cds),dx,lo,hi,prob_lo,prob_hi)

       end select

    end do

  end subroutine cal_phi_tag_test

  subroutine cal_phi_test_2d(phi,dx,lo,hi,prob_lo,prob_hi)

    integer           :: lo(2),hi(2)
    double precision  :: phi(lo(1):hi(1),lo(2):hi(2))
    double precision  :: dx
    double precision  :: prob_lo(2),prob_hi(2)

    ! local variables
    integer            :: i,j
    real(kind = dp_t)  :: x,y,x0,y0,r2,r0 = 50.d0

    ! get the point in the middle
    x0 = prob_lo(1) + 0.5d0*( prob_hi(1) - prob_lo(1) )
    y0 = prob_lo(2) + 0.5d0*( prob_hi(2) - prob_lo(2) )

    ! equiaxed dendrite growth
    do j = lo(2), hi(2)
     y = prob_lo(2) + (dble(j)+0.5d0) * dx
     do i = lo(1), hi(1)
        x = prob_lo(1) + (dble(i)+0.5d0) * dx

        r2 = sqrt( (x - x0)**2 + (y - y0)**2)

        phi(i,j) = -tanh(r2-r0)

     end do
    end do

  end subroutine cal_phi_test_2d

  subroutine cal_phi_test_3d(phi,dx,lo,hi,prob_lo,prob_hi)

    integer           :: lo(3),hi(3)
    double precision  :: phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    double precision  :: dx
    double precision  :: prob_lo(3),prob_hi(3)

    ! local variables
    integer            :: i,j,k
    real(kind = dp_t)  :: x,y,z,x0,y0,z0,r2,r0 = 50.d0

    ! get the point in the middle
    x0 = prob_lo(1) + 0.5d0*( prob_hi(1) - prob_lo(1) )
    y0 = prob_lo(2) + 0.5d0*( prob_hi(2) - prob_lo(2) )
    z0 = prob_lo(3) + 0.5d0*( prob_hi(3) - prob_lo(3) )

    ! equiaxed dendrite growth
    do k = lo(3), hi(3)
      z = prob_lo(3) + (dble(k)+0.5d0) * dx
      do j = lo(2), hi(2)
        y = prob_lo(2) + (dble(j)+0.5d0) * dx
          do i = lo(1), hi(1)
            x = prob_lo(1) + (dble(i)+0.5d0) * dx

            r2 = sqrt( (x - x0)**2 + (y - y0)**2 + (z - z0)**2)

            phi(i,j,k) = -tanh(r2-r0)

         end do
      end do
    end do

  end subroutine cal_phi_test_3d


end module write_plotfile_module
