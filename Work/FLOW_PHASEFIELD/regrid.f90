module regrid_module

  use ml_layout_module
  use multifab_module
  use multifab_physbc_module
  use make_new_grids_module
  !use ml_restriction_module
  use multifab_fill_ghost_module
  use ml_restrict_fill_module

  implicit none

  private

  public :: regrid

contains

  subroutine regrid(mla,flw,flw_0,flw_1,surf,phi,ad,phi_new,noise,nlevs,max_levs,dx, &
                    the_bc_tower,amr_buf_width,max_grid_size,pfpara,nc_fl)

    type(ml_layout), intent(inout) :: mla
    type(multifab) , intent(inout) :: phi(:),phi_new(:),ad(:),noise(:)
    type(multifab) , intent(inout) :: flw(:),flw_0(:),flw_1(:),surf(:)
    integer        , intent(inout) :: nlevs, max_levs,nc_fl
    real(dp_t)     , intent(in   ) :: dx(:)
    type(bc_tower) , intent(inout) :: the_bc_tower
    integer        , intent(in   ) :: amr_buf_width, max_grid_size
    type(pf_para)  , intent(in   ) :: pfpara

    select case(pfpara%do_solve_flow)
    case(0)
      call regrid_den(mla,phi,ad,phi_new,noise,nlevs,max_levs,dx, &
                      the_bc_tower,amr_buf_width,max_grid_size,pfpara)
    case(1)
      call regrid_flw(mla,flw,flw_0,flw_1,surf,phi,ad,phi_new,noise,nlevs,max_levs,dx, &
                    the_bc_tower,amr_buf_width,max_grid_size,pfpara,nc_fl)
    end select

  end subroutine regrid

  subroutine regrid_den(mla,phi,ad,phi_new,noise,nlevs,max_levs,dx, &
                    the_bc_tower,amr_buf_width,max_grid_size,pfpara)

    type(ml_layout), intent(inout) :: mla
    type(multifab) , intent(inout) :: phi(:),phi_new(:),ad(:),noise(:)
    integer        , intent(inout) :: nlevs, max_levs
    real(dp_t)     , intent(in   ) :: dx(:)
    type(bc_tower) , intent(inout) :: the_bc_tower
    integer        , intent(in   ) :: amr_buf_width, max_grid_size
    type(pf_para)  , intent(in   ) :: pfpara

    ! local variables
    type(layout)      :: la_array(max_levs)
    type(ml_layout)   :: mla_old
    type(ml_boxarray) :: mba
    type(multifab)    :: phi_old(max_levs),ad_old(max_levs)

    integer :: dm,n,nl,n_buffer

    logical :: new_grid,properly_nested

    dm = mla%dim

    ! create a copy of the original mla
    call ml_layout_build(mla_old,mla%mba,mla%pmask)

    do n=1,nlevs
       ! create copies of the old data
       ! make sure to use mla_old since we will be destroying mla
       call multifab_build(phi_old(n),mla_old%la(n),3,1)
       call multifab_copy_c(phi_old(n),1,phi(n),1,3)

       call multifab_build(ad_old(n),mla_old%la(n),1,1)
       call multifab_copy_c(ad_old(n),1,ad(n),1,1)

       ! get rid of the original multifab so we can create a new one
       ! with a new grid structure and the same name
       call multifab_destroy(phi(n))
       call multifab_destroy(ad(n))
    end do

    ! destory everything
    do n =1,nlevs
       call multifab_destroy(phi_new(n))
       call multifab_destroy(noise(n))
       !call multifab_destroy(rk_mat(n))
       !call multifab_destroy(fab_t(n))
    end do

    ! get rid of the original mla so we can create a new one
    ! with a new grid structure and the same name
    call destroy(mla)

    ! mba is big enough to hold max_levs levels
    ! even though we know we had nlevs last time, we might
    ! want more or fewer levels after regrid (if nlevs < max_levs)
    call ml_boxarray_build_n(mba,max_levs,dm)

    ! tell mba about the ref_ratio between levels
    ! mba%rr(n-1,i) is the refinement ratio between levels n-1 and n in direction i
    ! we use refinement ratio of 2 in every direction between all levels
    do n=2,max_levs
       mba%rr(n-1,:) = 2
    enddo

    ! copy the level 1 boxarray
    call copy(mba%bas(1),mla_old%mba%bas(1))

    ! set the problem domain at all levels
    mba%pd(1) = mla_old%mba%pd(1)
    do n=2,max_levs
       mba%pd(n) = refine(mba%pd(n-1),mba%rr((n-1),:))
    end do

    ! build the level 1 layout.
    call layout_build_ba(la_array(1),mba%bas(1),mba%pd(1),mla_old%pmask)

    ! This makes sure the boundary conditions are properly defined everywhere
    call bc_tower_level_build(the_bc_tower,1,la_array(1))

    ! build level 1 multifab
    call multifab_build(phi(1),la_array(1),3,1)
    call multifab_build(ad(1),la_array(1),1,1)

    ! copy level 1 data from original multifab
    call multifab_copy_c(phi(1),1,phi_old(1),1,3)
    call multifab_copy_c(ad(1),1,ad_old(1),1,1)

    nl = 1
    new_grid = .true.

    ! this is the number of level n+1 buffer cells we require between levels
    ! n and n+2 for proper nesting
    n_buffer = 4

    do while ( (nl .lt. max_levs) .and. new_grid )

       ! need to fill ghost cells here in case we use them in tagging
       call multifab_fill_boundary(phi(nl))
       call multifab_physbc(phi(nl),1,1,3,the_bc_tower%bc_tower_array(nl))

       call multifab_fill_boundary(ad(nl))
       call multifab_physbc(ad(nl),1,1,1,the_bc_tower%bc_tower_array(nl))

       ! determine whether we need finer grids based on tagging criteria
       ! if so, return new_grid=T and the la_array(nl+1)
       call make_new_grids_den(new_grid,la_array(nl),la_array(nl+1),phi(nl),dx(nl), &
                           amr_buf_width,mba%rr(nl,1),nl,max_grid_size,pfpara)

       if (new_grid) then

          ! set the level nl+1 boxarray
          call copy(mba%bas(nl+1),get_boxarray(la_array(nl+1)))

          ! enforce proper nesting within the grid creation procedure
          if (nl .ge. 2) then

             ! Test on whether grids are already properly nested
             properly_nested = ml_boxarray_properly_nested(mba, n_buffer, mla_old%pmask, &
                                                           max_fine_level=nl+1)

             if (.not. properly_nested) then

                do n = 2,nl
                   ! Delete old multifabs so that we can rebuild them.
                   call destroy(phi(n))
                   call destroy(ad(n))
                end do

                ! Change the layout at levels 2 through nl so the new grid
                ! structure is properly nested
                call enforce_proper_nesting(mba,la_array,max_grid_size)

                ! Loop over all the lower levels which we might have changed
                ! when we enforced proper nesting.
                do n = 2,nl

                   ! This makes sure the boundary conditions are properly defined everywhere
                   call bc_tower_level_build(the_bc_tower,n,la_array(n))

                   ! Rebuild the lower level data again if it changed.
                   call multifab_build(phi(n),la_array(n),3,1)
                   call multifab_build(ad(n),la_array(n),1,1)

                   ! first fill all refined cells by interpolating from coarse
                   ! data underneath...
                   call fillpatch(phi(n),phi(n-1),phi(n)%ng,mba%rr(n-1,:), &
                                  the_bc_tower%bc_tower_array(n-1), &
                                  the_bc_tower%bc_tower_array(n), &
                                  1,1,1,3)
                   call fillpatch(ad(n),ad(n-1),ad(n)%ng,mba%rr(n-1,:), &
                                  the_bc_tower%bc_tower_array(n-1), &
                                  the_bc_tower%bc_tower_array(n), &
                                  1,1,1,1)

                   ! ... then overwrite with the original data at that level, if it existed
                   if (mla_old%nlevel .ge. n) then
                      call multifab_copy_c(phi(n),1,phi_old(n),1,3)
                      call multifab_copy_c(ad(n),1,ad_old(n),1,1)
                   end if

                end do

             end if ! if (.not. properly_nested)

          end if ! if (nl .ge. 2) then

          ! Define bc_tower at level nl+1.
          call bc_tower_level_build(the_bc_tower,nl+1,la_array(nl+1))

          ! Build the level nl+1 data only.
          call multifab_build(phi(nl+1),la_array(nl+1),3,1)
          call multifab_build(ad(nl+1),la_array(nl+1),1,1)

          ! first fill all refined cells by interpolating from coarse
          ! data underneath...
          call fillpatch(phi(nl+1),phi(nl),phi(nl)%ng,mba%rr(nl,:), &
                         the_bc_tower%bc_tower_array(nl), &
                         the_bc_tower%bc_tower_array(nl+1), &
                         1,1,1,3)
          call fillpatch(ad(nl+1),ad(nl),ad(nl)%ng,mba%rr(nl,:), &
                         the_bc_tower%bc_tower_array(nl), &
                         the_bc_tower%bc_tower_array(nl+1), &
                         1,1,1,1)

          ! ... then overwrite with the original data at that level, if it existed
          if (mla_old%nlevel .ge. nl+1) then
             call multifab_copy_c(phi(nl+1),1,phi_old(nl+1),1,3)
             call multifab_copy_c(ad(nl+1),1,ad_old(nl+1),1,1)
          end if

          nl = nl+1
          nlevs = nl

       end if

    end do

    nlevs = nl

    ! Note: This build actually sets mla%la(n) = la_array(n) so we mustn't delete
    !       la_array(n).  Doing the build this way means we don't have to re-create
    !       all the multifabs because we have kept the same layouts.
    call ml_layout_build_la_array(mla,la_array,mba,mla_old%pmask,nlevs)

    ! This makes sure the boundary conditions are properly defined everywhere
    do n=1,nlevs
       call bc_tower_level_build(the_bc_tower,n,la_array(n))
    end do

    
    call ml_restrict_and_fill(nlevs, phi, mla%mba%rr, the_bc_tower%bc_tower_array,1,1,3)
    call ml_restrict_and_fill(nlevs, ad, mla%mba%rr, the_bc_tower%bc_tower_array,1,1,1)

    do n=1,nlevs
      call multifab_build(phi_new(n),mla%la(n),3,1)
      call multifab_build(noise(n),mla%la(n),1,1)
    end do

    do n=1,nlevs
       call multifab_copy(phi_new(n), phi(n), 1)
    end do

    do n=1,mla_old%nlevel
       call destroy(phi_old(n))
       call destroy(ad_old(n))
    end do

    call destroy(mba)
    call destroy(mla_old)

  end subroutine regrid_den

  subroutine regrid_flw(mla,flw,flw_0,flw_1,surf,phi,ad,phi_new,noise,nlevs,max_levs,dx, &
                    the_bc_tower,amr_buf_width,max_grid_size,pfpara,nc_fl)

    type(ml_layout), intent(inout) :: mla
    type(multifab) , intent(inout) :: flw(:),flw_0(:),flw_1(:),surf(:)
    type(multifab) , intent(inout) :: phi(:),phi_new(:),ad(:),noise(:)
    integer        , intent(inout) :: nlevs, max_levs,nc_fl
    real(dp_t)     , intent(in   ) :: dx(:)
    type(bc_tower) , intent(inout) :: the_bc_tower
    integer        , intent(in   ) :: amr_buf_width, max_grid_size
    type(pf_para)  , intent(in   ) :: pfpara

    ! local variables
    type(layout)      :: la_array(max_levs)
    type(ml_layout)   :: mla_old
    type(ml_boxarray) :: mba
    type(multifab)    :: phi_old(max_levs),ad_old(max_levs)
    type(multifab)    :: flw_old(max_levs)


    integer :: dm,n,nl,n_buffer,ncomp_flw
    logical :: new_grid,properly_nested
    double precision :: alpha

    dm = mla%dim

    select case(dm)
      case (2)
        ncomp_flw = 21
      case (3)
        ncomp_flw = 42
    end select

    ! create a copy of the original mla
    call ml_layout_build(mla_old,mla%mba,mla%pmask)

    do n=1,nlevs
       ! create copies of the old data
       ! make sure to use mla_old since we will be destroying mla
       call multifab_build(phi_old(n),mla_old%la(n),3,1)
       call multifab_copy_c(phi_old(n),1,phi(n),1,3)

       call multifab_build(ad_old(n),mla_old%la(n),1,1)
       call multifab_copy_c(ad_old(n),1,ad(n),1,1)

       call multifab_build(flw_old(n),mla_old%la(n),ncomp_flw,1)
       call multifab_copy_c(flw_old(n),1,flw(n),1,ncomp_flw)

       ! get rid of the original multifab so we can create a new one
       ! with a new grid structure and the same name
       call multifab_destroy(phi(n))
       call multifab_destroy(ad(n))
       call multifab_destroy(flw(n))
    end do

    ! destory everything
    do n =1,nlevs
       call multifab_destroy(phi_new(n))
       call multifab_destroy(noise(n))
       call multifab_destroy(flw_0(n))
       call multifab_destroy(flw_1(n))
       call multifab_destroy(surf(n))
    end do

    ! get rid of the original mla so we can create a new one
    ! with a new grid structure and the same name
    call destroy(mla)

    ! mba is big enough to hold max_levs levels
    ! even though we know we had nlevs last time, we might
    ! want more or fewer levels after regrid (if nlevs < max_levs)
    call ml_boxarray_build_n(mba,max_levs,dm)

    ! tell mba about the ref_ratio between levels
    ! mba%rr(n-1,i) is the refinement ratio between levels n-1 and n in direction i
    ! we use refinement ratio of 2 in every direction between all levels
    do n=2,max_levs
       mba%rr(n-1,:) = 2
    enddo

    ! copy the level 1 boxarray
    call copy(mba%bas(1),mla_old%mba%bas(1))

    ! set the problem domain at all levels
    mba%pd(1) = mla_old%mba%pd(1)
    do n=2,max_levs
       mba%pd(n) = refine(mba%pd(n-1),mba%rr((n-1),:))
    end do

    ! build the level 1 layout.
    call layout_build_ba(la_array(1),mba%bas(1),mba%pd(1),mla_old%pmask)

    ! This makes sure the boundary conditions are properly defined everywhere
    call bc_tower_level_build(the_bc_tower,1,la_array(1))

    ! build level 1 multifab
    call multifab_build(phi(1),la_array(1),3,1)
    call multifab_build(ad(1),la_array(1),1,1)
    call multifab_build(flw(1),la_array(1),ncomp_flw,1)

    ! copy level 1 data from original multifab
    call multifab_copy_c(phi(1),1,phi_old(1),1,3)
    call multifab_copy_c(ad(1),1,ad_old(1),1,1)
    call multifab_copy_c(flw(1),1,flw_old(1),1,ncomp_flw)

    nl = 1
    new_grid = .true.

    ! this is the number of level n+1 buffer cells we require between levels
    ! n and n+2 for proper nesting
    n_buffer = 4

    do while ( (nl .lt. max_levs) .and. new_grid )

       ! need to fill ghost cells here in case we use them in tagging
       call multifab_fill_boundary(phi(nl))
       call multifab_physbc(phi(nl),1,1,3,the_bc_tower%bc_tower_array(nl))

       call multifab_fill_boundary(ad(nl))
       call multifab_physbc(ad(nl),1,1,1,the_bc_tower%bc_tower_array(nl))

        ! need to fill ghost cells here in case we use them in tagging
       call multifab_fill_boundary(flw(nl))
       call multifab_physbc(flw(nl),1,1,ncomp_flw,the_bc_tower%bc_tower_array(nl))

       ! determine whether we need finer grids based on tagging criteria
       ! if so, return new_grid=T and the la_array(nl+1)
       call make_new_grids_flw(new_grid,la_array(nl),la_array(nl+1),flw(nl),phi(nl),dx(nl), &
                           amr_buf_width,mba%rr(nl,1),nl,max_grid_size,pfpara)

       if (new_grid) then

          ! set the level nl+1 boxarray
          call copy(mba%bas(nl+1),get_boxarray(la_array(nl+1)))

          ! enforce proper nesting within the grid creation procedure
          if (nl .ge. 2) then

             ! Test on whether grids are already properly nested
             properly_nested = ml_boxarray_properly_nested(mba, n_buffer, mla_old%pmask, &
                                                           max_fine_level=nl+1)

             if (.not. properly_nested) then

                do n = 2,nl
                   ! Delete old multifabs so that we can rebuild them.
                   call destroy(phi(n))
                   call destroy(ad(n))
                   call destroy(flw(n))
                end do

                ! Change the layout at levels 2 through nl so the new grid
                ! structure is properly nested
                call enforce_proper_nesting(mba,la_array,max_grid_size)

                ! Loop over all the lower levels which we might have changed
                ! when we enforced proper nesting.
                do n = 2,nl

                   ! This makes sure the boundary conditions are properly defined everywhere
                   call bc_tower_level_build(the_bc_tower,n,la_array(n))

                   ! Rebuild the lower level data again if it changed.
                   call multifab_build(phi(n),la_array(n),3,1)
                   call multifab_build(ad(n),la_array(n),1,1)
                   call multifab_build(flw(n),la_array(n),ncomp_flw,1)

                   ! first fill all refined cells by interpolating from coarse
                   ! data underneath...
                   call fillpatch(phi(n),phi(n-1),phi(n)%ng,mba%rr(n-1,:), &
                                  the_bc_tower%bc_tower_array(n-1), &
                                  the_bc_tower%bc_tower_array(n), &
                                  1,1,1,3)
                   call fillpatch(ad(n),ad(n-1),ad(n)%ng,mba%rr(n-1,:), &
                                  the_bc_tower%bc_tower_array(n-1), &
                                  the_bc_tower%bc_tower_array(n), &
                                  1,1,1,1)

                   ! calculate the alpha
                   !call cal_tau_alpha(alpha, n-1,pfpara%tau_LBM)
                   call fillpatch_flw(flw(n),flw(n-1),flw(n)%ng,mba%rr(n-1,:), &
                                  the_bc_tower%bc_tower_array(n-1), &
                                  the_bc_tower%bc_tower_array(n), &
                                  1,1,1,ncomp_flw,n-1,pfpara%tau_LBM)

                   ! ... then overwrite with the original data at that level, if it existed
                   if (mla_old%nlevel .ge. n) then
                      call multifab_copy_c(phi(n),1,phi_old(n),1,3)
                      call multifab_copy_c(ad(n),1,ad_old(n),1,1)
                      call multifab_copy_c(flw(n),1,flw_old(n),1,ncomp_flw)
                   end if

                end do

             end if ! if (.not. properly_nested)

          end if ! if (nl .ge. 2) then

          ! Define bc_tower at level nl+1.
          call bc_tower_level_build(the_bc_tower,nl+1,la_array(nl+1))

          ! Build the level nl+1 data only.
          call multifab_build(phi(nl+1),la_array(nl+1),3,1)
          call multifab_build(ad(nl+1),la_array(nl+1),1,1)
          call multifab_build(flw(nl+1),la_array(nl+1),ncomp_flw,1)

          ! first fill all refined cells by interpolating from coarse
          ! data underneath...
          call fillpatch(phi(nl+1),phi(nl),phi(nl)%ng,mba%rr(nl,:), &
                         the_bc_tower%bc_tower_array(nl), &
                         the_bc_tower%bc_tower_array(nl+1), &
                         1,1,1,3)
          call fillpatch(ad(nl+1),ad(nl),ad(nl)%ng,mba%rr(nl,:), &
                         the_bc_tower%bc_tower_array(nl), &
                         the_bc_tower%bc_tower_array(nl+1), &
                         1,1,1,1)

          ! calculate the alpha
          !call cal_tau_alpha(alpha, nl,pfpara%tau_LBM)
          call fillpatch_flw(flw(nl+1),flw(nl),flw(nl)%ng,mba%rr(nl,:), &
                         the_bc_tower%bc_tower_array(nl), &
                         the_bc_tower%bc_tower_array(nl+1), &
                         1,1,1,ncomp_flw,nl,pfpara%tau_LBM)

          ! ... then overwrite with the original data at that level, if it existed
          if (mla_old%nlevel .ge. nl+1) then
             call multifab_copy_c(phi(nl+1),1,phi_old(nl+1),1,3)
             call multifab_copy_c(ad(nl+1),1,ad_old(nl+1),1,1)
             call multifab_copy_c(flw(nl+1),1,flw_old(nl+1),1,ncomp_flw)
          end if

          nl = nl+1
          nlevs = nl

       end if

    end do

    nlevs = nl

    ! Note: This build actually sets mla%la(n) = la_array(n) so we mustn't delete
    !       la_array(n).  Doing the build this way means we don't have to re-create
    !       all the multifabs because we have kept the same layouts.
    call ml_layout_build_la_array(mla,la_array,mba,mla_old%pmask,nlevs)

    ! This makes sure the boundary conditions are properly defined everywhere
    do n=1,nlevs
       call bc_tower_level_build(the_bc_tower,n,la_array(n))
    end do

    call ml_restrict_and_fill(nlevs, phi, mla%mba%rr, the_bc_tower%bc_tower_array,1,1,3)
    call ml_restrict_and_fill(nlevs, ad, mla%mba%rr, the_bc_tower%bc_tower_array,1,1,1)
    call ml_restrict_and_fill_flw(nlevs, flw, mla%mba%rr, the_bc_tower%bc_tower_array,&
                                  pfpara%tau_LBM,1,1,ncomp_flw)

    do n=1,nlevs
      call multifab_build(phi_new(n),mla%la(n),3,1)
      call multifab_build(noise(n),mla%la(n),1,1)
      call multifab_build(flw_0(n),mla%la(n),nc_fl,1)
      call multifab_build(flw_1(n),mla%la(n),nc_fl,1)
      call multifab_build(surf(n),mla%la(n),1,1)
    end do

    do n=1,nlevs
       call multifab_copy(phi_new(n), phi(n), 1)
    end do

    do n=1,mla_old%nlevel
       call destroy(phi_old(n))
       call destroy(ad_old(n))
       call destroy(flw_old(n))
    end do

    call destroy(mba)
    call destroy(mla_old)

  end subroutine regrid_flw

end module regrid_module
