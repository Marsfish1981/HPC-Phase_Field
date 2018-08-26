module temperature_field_module

  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use multifab_fill_ghost_module
  !use ml_restriction_module
  use mg_module
  use bndry_reg_module
  use cc_stencil_fill_module
  use ml_solve_module
  use pf_utility_module
  use ml_restrict_fill_module
  use ml_cc_restriction_module

  implicit none

  private

  public :: advance_temp_field

contains

  subroutine advance_temp_field(mla,phi_old,phi_new,dx,dt,reset_plotnum,the_bc_tower,pfpara,t_time, &
                                t_l,t_h,prob_lo,prob_hi,istep,cal_mode)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: phi_old(:),phi_new(:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: dt,t_time,t_h,t_l
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(pf_para),   intent(in   ) :: pfpara
    real(kind=dp_t), intent(in   ) :: prob_lo(:),prob_hi(:)
    integer        , intent(in   ) :: istep,cal_mode,reset_plotnum


    ! local variables
    integer :: i,dm,n,nlevs
    logical :: do_implicit_solve

    type(multifab) :: flux(mla%nlevel,mla%dim)

    type(multifab) :: alpha(mla%nlevel)
    type(multifab) :: beta(mla%nlevel,mla%dim)
    type(multifab) :: rhs(mla%nlevel)

    type(multifab), allocatable :: rk_mat(:)

    !type(bndry_reg) :: fine_flx(2:mla%nlevel)

    real(dp_t) ::  xa(mla%dim),  xb(mla%dim)
    real(dp_t) :: pxa(mla%dim), pxb(mla%dim)

    integer :: stencil_order,do_diagnostics

    type(layout) :: la
    type(box) :: pd

    real(dp_t) :: dx_vector(mla%nlevel,mla%dim)
    real(dp_t) :: delta_Cooling
    real(dp_t) :: temp_l, temp_h

    if(pfpara%do_solve_flow == 1 .and. pfpara%solve_flow_only == 1) return;

    dm = mla%dim
    nlevs = mla%nlevel

    if(pfpara%coupled_mode .eq. 0) then! temperature will be skipped

       if(pfpara%cal_tem_mode .eq. 0) then

          ! Copy the data from old phi
          do n=1,nlevs
             call multifab_copy_c(phi_new(n),3,phi_old(n),3,1,1)
          end do

       else if(pfpara%cal_tem_mode .eq. 1) then
          !if( ( pfpara%temp_l - t_time * pfpara%Rate_Cooling ) .ge. pfpara%temp_lowcut ) then
             !delta_Cooling = -dt * pfpara%Rate_Cooling

             if(pfpara%cooling_mode .eq. 0) then
                temp_l = t_l - dt * dble(istep+reset_plotnum) * pfpara%Rate_Cooling

                if(temp_l .le. pfpara%temp_lowcut) then
                    temp_l = pfpara%temp_lowcut
                end if

                temp_h = t_h - t_l + temp_l

             else if(pfpara%cooling_mode .eq. 1) then
                temp_l = t_l - dt * dble(istep+reset_plotnum) * pfpara%Rate_Cooling
                temp_h = t_h

                if(temp_l .le. pfpara%temp_lowcut) then
                    temp_l = pfpara%temp_lowcut
                end if

             else if(pfpara%cooling_mode .eq. 2) then
                temp_l = t_l
                temp_h = t_h - dt * dble(istep+reset_plotnum) * pfpara%Rate_Cooling

                if(temp_h .le. pfpara%temp_lowcut) then
                    temp_h = pfpara%temp_lowcut
                end if

             end if

             call update_th_dev(mla,phi_new,dx,dt,temp_h,temp_l,prob_lo,prob_hi,the_bc_tower,pfpara%temp_lowcut)

        else if(pfpara%cal_tem_mode .eq. 2) then

            if(cal_mode == 0) then
                temp_h = pfpara%T_D(1)
            else
                temp_h = pfpara%T_D(istep)
            end if

            !multifab_setval_c(phi_new, temp_h, 3, 1, all=.true.)
            call update_th_dev(mla,phi_new,dx,dt,temp_h,temp_h,prob_lo,prob_hi,the_bc_tower,pfpara%temp_lowcut)

        end if

       return

    end if

    if(pfpara%cal_tem_mode .eq. 0) then
      do_implicit_solve = .false.
    else
      do_implicit_solve = .true.
    end if

    stencil_order = 2
    do_diagnostics = 0

    if (do_implicit_solve) then

        allocate(rk_mat(nlevs))

        do n=1,nlevs
           call multifab_build(rk_mat(n),mla%la(n),1,1)
        end do

        do n=1,nlevs

          ! Copy all values from old_phi to rk_mat
          ! rk_mat = T0
          call multifab_copy_c(rk_mat(n),1,phi_old(n),3,1,phi_old(1)%ng)

          ! set alpha=1
          call multifab_build(alpha(n),mla%la(n),1,0)
          call setval(alpha(n),1.d0,all=.true.)

          ! set beta=dt
          do i=1,dm
             call multifab_build_edge(beta(n,i),mla%la(n),1,0,i)
             call setval(beta(n,i), dt*(pfpara%t_DM), all=.true.)
          end do

          ! copy phi into rhs
          call multifab_build(rhs(n),mla%la(n),1,0)

          ! rhs = pf_new
          call multifab_copy_c(rhs(n),1,phi_new(n),1,1,0)

          ! rhs = pf_new - pf_old
          call multifab_sub_sub_c(rhs(n),1,phi_old(n),1,1,0)

          ! rhs = (pf_new - pf_old) * r_hyper/2
          call multifab_mult_mult_s_c(rhs(n), 1, (0.5d0*pfpara%r_hyper),1,0)

          ! rhs = (pf_new - pf_old) * r_hyper/2 + T0
          call multifab_plus_plus_c(rhs(n), 1, phi_old(n), 3, 1,0)

       end do

       !do n = 2,nlevs
       !   call bndry_reg_build(fine_flx(n),mla%la(n),ml_layout_get_pd(mla,n))
       !end do

       do n=1,nlevs
          dx_vector(n,1:dm) = dx(n)
       end do

       ! solve (alpha - del dot beta grad) phi = rhs to obtain phi
       ! call ml_cc_solve(mla,rhs,phi,fine_flx,alpha,beta,dx_vector,the_bc_tower,bc_comp)
       call ml_cc_solve(mla, rhs, rk_mat, alpha, beta, dx_vector,the_bc_tower,1)

       ! Copy the solutions
       do n=1,nlevs
          call multifab_copy_c(phi_new(n),3,rk_mat(n),1,1,phi_new(1)%ng)
       end do

       ! deallocate memory
       do n=1,nlevs
          call multifab_destroy(alpha(n))
          call multifab_destroy(  rhs(n))
          do i=1,dm
             call multifab_destroy(beta(n,i))
          end do
       end do
       !do n = 2,nlevs
       !   call bndry_reg_destroy(fine_flx(n))
       !end do
     

       do n=1,nlevs
          call destroy(rk_mat(n))
       end do

    else

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! explicit time advancement
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! build the flux(:,:) multifabs
       do n=1,nlevs
          do i=1,dm
             ! flux(n,i) has one component, zero ghost cells, and is nodal in direction i
             call multifab_build_edge(flux(n,i),mla%la(n),1,0,i)
          end do
       end do

       ! compute the face-centered gradients in each direction
       call compute_flux_th(mla,phi_old,flux,dx,the_bc_tower)

       ! calculate the RK matrices at the first order
       call update_th(mla,phi_new,phi_old,flux,dx,dt,the_bc_tower,pfpara)

       ! make sure to destroy the multifab or you'll leak memory
       do n=1,nlevs
          do i=1,dm
             call multifab_destroy(flux(n,i))
          end do
       end do

    end if

  end subroutine advance_temp_field

  subroutine compute_flux_th(mla,phi_old,flux,dx,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: phi_old(:)
    type(multifab) , intent(inout) :: flux(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: dm, ng_p, ng_f, i, n, nlevs

    real(kind=dp_t), pointer ::  pp(:,:,:,:)
    real(kind=dp_t), pointer :: fxp(:,:,:,:)
    real(kind=dp_t), pointer :: fyp(:,:,:,:)
    real(kind=dp_t), pointer :: fzp(:,:,:,:)

    dm    = mla%dim
    nlevs = mla%nlevel

    ng_p = phi_old(1)%ng
    ng_f = flux(1,1)%ng

    do n=1,nlevs
       do i=1,nfabs(phi_old(n))
          pp  => dataptr(phi_old(n),i)
          fxp => dataptr(flux(n,1),i)
          fyp => dataptr(flux(n,2),i)
          lo = lwb(get_box(phi_old(n),i))
          hi = upb(get_box(phi_old(n),i))
          select case(dm)
          case (2)
             call compute_flux_th_2d(pp(:,:,1,1), ng_p, &
                                     fxp(:,:,1,1),  fyp(:,:,1,1), ng_f, &
                                     lo, hi, dx(n), &
                                     the_bc_tower%bc_tower_array(n)%adv_bc_level_array(i,:,:,1))
          case (3)
             fzp => dataptr(flux(n,3),i)
             call compute_flux_th_3d(pp(:,:,:,1), ng_p, &
                                     fxp(:,:,:,1),  fyp(:,:,:,1), fzp(:,:,:,1), ng_f, &
                                     lo, hi, dx(n), &
                                     the_bc_tower%bc_tower_array(n)%adv_bc_level_array(i,:,:,1))
          end select
       end do
    end do

    ! set level n-1 fluxes to be the average of the level n fluxes covering it
    ! the loop over nlevs must count backwards to make sure the finer grids are done first
    do n=nlevs,2,-1
       do i=1,dm
          call ml_edge_restriction_c(flux(n-1,i),1,flux(n,i),1,mla%mba%rr(n-1,:),i,1)
       end do
    end do

  end subroutine compute_flux_th

  subroutine compute_flux_th_2d(phi, ng_p, fluxx, fluxy, ng_f, lo, hi, dx,adv_bc)

    integer          :: lo(2), hi(2), ng_p, ng_f
    double precision ::   phi(lo(1)-ng_p:,lo(2)-ng_p:)
    double precision :: fluxx(lo(1)-ng_f:,lo(2)-ng_f:)
    double precision :: fluxy(lo(1)-ng_f:,lo(2)-ng_f:)
    double precision :: dx
    integer          :: adv_bc(:,:)

    ! local variables
    integer i,j

    ! x-fluxes

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)+1
          fluxx(i,j) = ( phi(i,j) - phi(i-1,j) ) / dx
       end do
    end do

    ! lo-x boundary conditions
    if (adv_bc(1,1) .eq. EXT_DIR) then
       i=lo(1)
       do j=lo(2),hi(2)
          ! divide by 0.5*dx since the ghost cell value represents
          ! the value at the wall, not the ghost cell-center
          fluxx(i,j) = ( phi(i,j) - phi(i-1,j) ) / (0.5d0*dx)
       end do
    else if (adv_bc(1,1) .eq. FOEXTRAP) then
       ! dphi/dn = 0
       fluxx(lo(1),lo(2):hi(2)) = 0.d0
    end if

    ! hi-x boundary conditions
    if (adv_bc(1,2) .eq. EXT_DIR) then
       i=hi(1)+1
       do j=lo(2),hi(2)
          ! divide by 0.5*dx since the ghost cell value represents
          ! the value at the wall, not the ghost cell-center
          fluxx(i,j) = ( phi(i,j) - phi(i-1,j) ) / (0.5d0*dx)
       end do
    else if (adv_bc(1,2) .eq. FOEXTRAP) then
       ! dphi/dn = 0
       fluxx(hi(1)+1,lo(2):hi(2)) = 0.d0
    end if

    ! y-fluxes
    do j=lo(2),hi(2)+1
       do i=lo(1),hi(1)
          fluxy(i,j) = ( phi(i,j) - phi(i,j-1) ) / dx
       end do
    end do

    ! lo-y boundary conditions
    if (adv_bc(2,1) .eq. EXT_DIR) then
       j=lo(2)
       do i=lo(1),hi(1)
          ! divide by 0.5*dx since the ghost cell value represents
          ! the value at the wall, not the ghost cell-center
          fluxy(i,j) = ( phi(i,j) - phi(i,j-1) ) / (0.5d0*dx)
       end do
    else if (adv_bc(2,1) .eq. FOEXTRAP) then
       ! dphi/dn = 0
       fluxy(lo(1):hi(1),lo(2)) = 0.d0
    end if

    ! hi-y boundary conditions
    if (adv_bc(2,2) .eq. EXT_DIR) then
       j=hi(2)+1
       do i=lo(1),hi(1)
          ! divide by 0.5*dx since the ghost cell value represents
          ! the value at the wall, not the ghost cell-center
          fluxy(i,j) = ( phi(i,j) - phi(i,j-1) ) / (0.5d0*dx)
       end do
    else if (adv_bc(2,2) .eq. FOEXTRAP) then
       ! dphi/dn = 0
       fluxy(lo(1):hi(1),hi(2)+1) = 0.d0
    end if

  end subroutine compute_flux_th_2d

  subroutine compute_flux_th_3d(phi, ng_p, fluxx, fluxy, fluxz, ng_f, &
                             lo, hi, dx, adv_bc)

    integer          :: lo(3), hi(3), ng_p, ng_f
    double precision ::   phi(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    double precision :: fluxx(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: fluxy(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: fluxz(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: dx
    integer          :: adv_bc(:,:)

    ! local variables
    integer i,j,k

    ! x-fluxes
    !$omp parallel do private(i,j,k)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1
             fluxx(i,j,k) = ( phi(i,j,k) - phi(i-1,j,k) ) / dx
          end do
       end do
    end do
    !$omp end parallel do

    ! lo-x boundary conditions
    if (adv_bc(1,1) .eq. EXT_DIR) then
       i=lo(1)
       !$omp parallel do private(j,k)
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             ! divide by 0.5*dx since the ghost cell value represents
             ! the value at the wall, not the ghost cell-center
             fluxx(i,j,k) = ( phi(i,j,k) - phi(i-1,j,k) ) / (0.5d0*dx)
          end do
       end do
       !$omp end parallel do
    else if (adv_bc(1,1) .eq. FOEXTRAP) then
       ! dphi/dn = 0
       fluxx(lo(1),lo(2):hi(2),lo(3):hi(3)) = 0.d0
    end if

    ! hi-x boundary conditions
    if (adv_bc(1,2) .eq. EXT_DIR) then
       i=hi(1)+1
       !$omp parallel do private(j,k)
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             ! divide by 0.5*dx since the ghost cell value represents
             ! the value at the wall, not the ghost cell-center
             fluxx(i,j,k) = ( phi(i,j,k) - phi(i-1,j,k) ) / (0.5d0*dx)
          end do
       end do
       !$omp end parallel do
    else if (adv_bc(1,2) .eq. FOEXTRAP) then
       ! dphi/dn = 0
       fluxx(hi(1)+1,lo(2):hi(2),lo(3):hi(3)) = 0.d0
    end if

    ! y-fluxes
    !$omp parallel do private(i,j,k)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)
             fluxy(i,j,k) = ( phi(i,j,k) - phi(i,j-1,k) ) / dx
          end do
       end do
    end do
    !$omp end parallel do

    ! lo-y boundary conditions
    if (adv_bc(2,1) .eq. EXT_DIR) then
       j=lo(2)
       !$omp parallel do private(i,k)
       do k=lo(3),hi(3)
          do i=lo(1),hi(1)
             ! divide by 0.5*dx since the ghost cell value represents
             ! the value at the wall, not the ghost cell-center
             fluxy(i,j,k) = ( phi(i,j,k) - phi(i,j-1,k) ) / (0.5d0*dx)
          end do
       end do
       !$omp end parallel do
    else if (adv_bc(2,1) .eq. FOEXTRAP) then
       ! dphi/dn = 0
       fluxy(lo(1):hi(1),lo(2),lo(3):hi(3)) = 0.d0
    end if

    ! hi-y boundary conditions
    if (adv_bc(2,2) .eq. EXT_DIR) then
       j=hi(2)+1
       !$omp parallel do private(i,k)
       do k=lo(3),hi(3)
          do i=lo(1),hi(1)
             ! divide by 0.5*dx since the ghost cell value represents
             ! the value at the wall, not the ghost cell-center
             fluxy(i,j,k) = ( phi(i,j,k) - phi(i,j-1,k) ) / (0.5d0*dx)
          end do
       end do
       !$omp end parallel do
    else if (adv_bc(2,2) .eq. FOEXTRAP) then
       ! dphi/dn = 0
       fluxy(lo(1):hi(1),hi(2)+1,lo(3):hi(3)) = 0.d0
    end if

    ! z-fluxes
    !$omp parallel do private(i,j,k)
    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             fluxz(i,j,k) = ( phi(i,j,k) - phi(i,j,k-1) ) / dx
          end do
       end do
    end do
    !$omp end parallel do

    ! lo-z boundary conditions
    if (adv_bc(3,1) .eq. EXT_DIR) then
       k=lo(3)
       !$omp parallel do private(i,j)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             ! divide by 0.5*dx since the ghost cell value represents
             ! the value at the wall, not the ghost cell-center
             fluxz(i,j,k) = ( phi(i,j,k) - phi(i,j,k-1) ) / (0.5d0*dx)
          end do
       end do
       !$omp end parallel do
    else if (adv_bc(3,1) .eq. FOEXTRAP) then
       ! dphi/dn = 0
       fluxz(lo(1):hi(1),lo(2):lo(3),lo(3)) = 0.d0
    end if

    ! hi-z boundary conditions
    if (adv_bc(3,2) .eq. EXT_DIR) then
       k=hi(3)+1
       !$omp parallel do private(i,j)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             ! divide by 0.5*dx since the ghost cell value represents
             ! the value at the wall, not the ghost cell-center
             fluxz(i,j,k) = ( phi(i,j,k) - phi(i,j,k-1) ) / (0.5d0*dx)
          end do
       end do
       !$omp end parallel do
    else if (adv_bc(3,2) .eq. FOEXTRAP) then
       ! dphi/dn = 0
       fluxz(lo(1):hi(1),lo(2):hi(2),hi(3)+1) = 0.d0
    end if

  end subroutine compute_flux_th_3d

  subroutine update_th(mla,phi_new,phi_old,flux,dx,dt,the_bc_tower,pfpara)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: phi_new(:),phi_old(:)
    type(multifab) , intent(in   ) :: flux(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: dt
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(pf_para),   intent(in   ) :: pfpara

    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: dm, ng_p, ng_f, i, n, nlevs

    real(kind=dp_t), pointer ::  pp(:,:,:,:)
    real(kind=dp_t), pointer :: ppn(:,:,:,:)
    real(kind=dp_t), pointer :: fxp(:,:,:,:)
    real(kind=dp_t), pointer :: fyp(:,:,:,:)
    real(kind=dp_t), pointer :: fzp(:,:,:,:)

    dm    = mla%dim
    nlevs = mla%nlevel

    ng_p = phi_old(1)%ng
    ng_f = flux(1,1)%ng

    do n=1,nlevs

       do i=1,nfabs(phi_old(n))
          pp  => dataptr(phi_old(n),i)
          ppn => dataptr(phi_new(n),i)
          fxp => dataptr(flux(n,1),i)
          fyp => dataptr(flux(n,2),i)
          lo = lwb(get_box(phi_old(n),i))
          hi = upb(get_box(phi_old(n),i))
          select case(dm)
          case (2)
             call update_th_2d(pp(:,:,1,3),ppn(:,:,1,3), pp(:,:,1,1),ppn(:,:,1,1), ng_p, &
                                fxp(:,:,1,1), fyp(:,:,1,1), ng_f, &
                                lo, hi, dx(n), dt,pfpara%t_DM, pfpara%r_hyper)
          case (3)
             fzp => dataptr(flux(n,3),i)
             call update_th_3d(pp(:,:,:,3),ppn(:,:,:,3), pp(:,:,:,1),ppn(:,:,:,1), ng_p, &
                                fxp(:,:,:,1), fyp(:,:,:,1), fzp(:,:,:,1), ng_f, &
                                lo, hi, dx(n), dt,pfpara%t_DM, pfpara%r_hyper)
          end select
       end do

    end do

    call ml_restrict_and_fill(nlevs, phi_new, mla%mba%rr, the_bc_tower%bc_tower_array,3,1,1)
    !if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
    !   call multifab_fill_boundary_c(phi_new(nlevs),3,1)

       ! fill non-periodic domain boundary ghost cells
    !   call multifab_physbc(phi_new(nlevs),3,1,1,the_bc_tower%bc_tower_array(nlevs))

    !else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
    !   do n=nlevs,2,-1
          ! set level n-1 data to be the average of the level n data covering it
    !      call ml_cc_restriction_c(phi_new(n-1),3,phi_new(n),3,mla%mba%rr(n-1,:),1)

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
    !      call multifab_fill_ghost_cells(phi_new(n),phi_new(n-1),ng_p,mla%mba%rr(n-1,:), &
    !                                     the_bc_tower%bc_tower_array(n-1), &
    !                                     the_bc_tower%bc_tower_array(n), &
    !                                     3,1,1)
    !   end do

    !end if

  end subroutine update_th

  subroutine update_th_2d(th_old,th_new,pf_old,pf_new, ng_p, fluxx, fluxy, &
                           ng_f, lo, hi, dx, dt,t_DM,r_hyper)

    integer          :: lo(2), hi(2), ng_p, ng_f
    double precision :: th_old(lo(1)-ng_p:,lo(2)-ng_p:),th_new(lo(1)-ng_p:,lo(2)-ng_p:)
    double precision :: pf_old(lo(1)-ng_p:,lo(2)-ng_p:),pf_new(lo(1)-ng_p:,lo(2)-ng_p:)
    double precision :: fluxx(lo(1)-ng_f:,lo(2)-ng_f:)
    double precision :: fluxy(lo(1)-ng_f:,lo(2)-ng_f:)
    double precision :: dx, dt,t_DM,r_hyper

    ! local variables
    integer i,j

    !$omp parallel do private(j,i)
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          th_new(i,j) = th_old(i,j) + dt * t_DM * &
               ( fluxx(i+1,j)-fluxx(i,j) + fluxy(i,j+1)-fluxy(i,j) ) / dx + &
               0.5d0 * r_hyper * ( pf_new(i,j) - pf_old(i,j) )
       end do
    end do
    !$omp end parallel do

  end subroutine update_th_2d

  subroutine update_th_3d(th_old,th_new,pf_old,pf_new, ng_p, fluxx, fluxy, fluxz, &
                           ng_f, lo, hi, dx, dt,t_DM,r_hyper)

    integer          :: lo(3), hi(3), ng_p, ng_f
    double precision :: th_old(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:),th_new(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    double precision :: pf_old(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:),pf_new(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    double precision :: fluxx(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: fluxy(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: fluxz(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: dx, dt,t_DM,r_hyper

    ! local variables
    integer i,j,k

    !$omp parallel do private(i,j,k)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             th_new(i,j,k) = th_old(i,j,k) + dt * t_DM * &
                  ( fluxx(i+1,j,k)-fluxx(i,j,k) &
                   +fluxy(i,j+1,k)-fluxy(i,j,k) &
                   +fluxz(i,j,k+1)-fluxz(i,j,k) ) / dx &
                   +0.5d0*r_hyper*( pf_new(i,j,k) - pf_old(i,j,k) )
          end do
       end do
    end do
    !$omp end parallel do

  end subroutine update_th_3d

  subroutine update_th_dev(mla,phi,dx,dt,t_h,t_l,prob_lo,prob_hi,the_bc_tower,t_cut)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: phi(:)
    real(kind=dp_t), intent(in   ) :: dx(:),prob_lo(:),prob_hi(:)
    real(kind=dp_t), intent(in   ) :: dt,t_h,t_l,t_cut
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: dm, ng_p, i, n, nlevs

    real(kind=dp_t), pointer ::  pp(:,:,:,:)

    dm    = mla%dim
    nlevs = mla%nlevel

    ng_p = phi(1)%ng

    do n=1,nlevs

       do i=1,nfabs(phi(n))
          pp  => dataptr(phi(n),i)

          lo = lwb(get_box(phi(n),i))
          hi = upb(get_box(phi(n),i))

          select case(dm)
          case (2)
             call update_th_dev_2d(pp(:,:,1,3),ng_p,lo,hi,prob_lo,prob_hi,t_l,t_h,dx(n),dt,t_cut)
          case (3)
             call update_th_dev_3d(pp(:,:,:,3),ng_p,lo,hi,prob_lo,prob_hi,t_l,t_h,dx(n),dt,t_cut)
          end select
       end do

    end do

    call ml_restrict_and_fill(nlevs, phi, mla%mba%rr, the_bc_tower%bc_tower_array,3,1,1)

    !if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
    !   call multifab_fill_boundary_c(phi(nlevs),3,1)

       ! fill non-periodic domain boundary ghost cells
    !   call multifab_physbc(phi(nlevs),3,1,1,the_bc_tower%bc_tower_array(nlevs))

    !else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
    !   do n=nlevs,2,-1
          ! set level n-1 data to be the average of the level n data covering it
    !      call ml_cc_restriction_c(phi(n-1),3,phi(n),3,mla%mba%rr(n-1,:),1)

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
    !      call multifab_fill_ghost_cells(phi(n),phi(n-1),ng_p,mla%mba%rr(n-1,:), &
    !                                     the_bc_tower%bc_tower_array(n-1), &
    !                                     the_bc_tower%bc_tower_array(n), &
    !                                     3,1,1)
    !   end do

    !end if

  end subroutine update_th_dev

  subroutine update_th_dev_2d(th,ng_p,lo,hi,prob_lo,prob_hi,t_l,t_h,dx,dt,t_cut)

    integer          :: lo(2), hi(2),ng_p
    double precision :: th(lo(1)-ng_p:,lo(2)-ng_p:)
    double precision :: dx,dt,t_l,t_h,prob_lo(2),prob_hi(2),t_cut

    ! local variables
    integer i,j
    double precision :: y_l,y

    ! The domain top to bottom distance
    y_l = prob_hi(2) - prob_lo(2)

    do j=lo(2),hi(2)
       y = (dble(j)+0.5d0) * dx
       do i=lo(1),hi(1)
          th(i,j) = t_l + (t_h - t_l) * y / y_l  ! setup a temperature gradient

          !if(th(i,j) < t_cut) then
          !  th(i,j) = t_cut
          !end if

       end do
    end do

  end subroutine update_th_dev_2d

  subroutine update_th_dev_3d(th,ng_p,lo,hi,prob_lo,prob_hi,t_l,t_h,dx,dt,t_cut)

    integer          :: lo(3), hi(3),ng_p
    double precision :: th(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    double precision :: dx, dt,t_l,t_h,prob_lo(3),prob_hi(3),t_cut

    ! local variables
    integer i,j,k
    double precision :: y_l,y

    ! The domain top to bottom distance
    y_l = prob_hi(2) - prob_lo(2)

    !$omp parallel do private(i,j,k)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          y = (dble(j)+0.5d0) * dx
          do i=lo(1),hi(1)
             th(i,j,k) = t_l + (t_h - t_l) * y / y_l  ! setup a temperature gradient

             !if(th(i,j,k) < t_cut) then
             !   th(i,j,k) = t_cut
             !end if

          end do
       end do
    end do
    !$omp end parallel do

  end subroutine update_th_dev_3d

end module temperature_field_module

