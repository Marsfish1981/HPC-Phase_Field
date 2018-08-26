module ml_cc_restriction_module

  use bl_constants_module
  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bc_module

  implicit none

  private

  public :: ml_cc_restriction, ml_cc_restriction_c,ml_cc_restriction_c_flw
  public :: ml_edge_restriction, ml_edge_restriction_c

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ml_cc_restriction_c(crse, cc, fine,cf, ir, nc)
    use cc_restriction_module
    type(multifab), intent(inout) :: fine
    type(multifab), intent(inout) :: crse
    integer, intent(in)           :: cc, cf, ir(:)
    integer, intent(in), optional :: nc

    integer             :: i, n, lnc, dm
    integer :: lo(get_dim(fine)), hi(get_dim(fine)), loc(get_dim(fine)), lof(get_dim(fine))
    real(dp_t), pointer :: fp(:,:,:,:), cp(:,:,:,:)
    type(layout)        :: lacfine,laf
    type(multifab)      :: cfine
    type(mfiter)        :: mfi
    type(box)           :: bx

    lnc = 1; if ( present(nc) ) lnc = nc

    laf = get_layout(fine)

    call layout_build_coarse(lacfine, laf, ir)

    call multifab_build(cfine, lacfine, nc = lnc, ng = 0)

    dm = get_dim(cfine)

    !$omp parallel private(mfi,n,i,bx,loc,lof,lo,hi,fp,cp)
    call mfiter_build(mfi,cfine,.true.)
    do n = 1, lnc
       do while(next_tile(mfi,i))
          bx = get_tilebox(mfi)

          lo  = lwb(bx)
          hi  = upb(bx)
          loc = lwb(get_pbox(cfine, i))
          lof = lwb(get_pbox( fine, i))

          fp    => dataptr(fine,     i, n+cf-1, 1)
          cp => dataptr(cfine, i, n,      1)
          select case (dm)
          case (1)
             call cc_restriction_1d(cp(:,1,1,1), loc, fp(:,1,1,1), lof, lo, hi, ir)
          case (2)
             call cc_restriction_2d(cp(:,:,1,1), loc, fp(:,:,1,1), lof, lo, hi, ir)
          case (3)
             call cc_restriction_3d(cp(:,:,:,1), loc, fp(:,:,:,1), lof, lo, hi, ir)
          end select
       end do
    end do
    !$omp end parallel

    call copy(crse, cc, cfine, 1, lnc)

    call multifab_destroy(cfine)

  end subroutine ml_cc_restriction_c

  subroutine ml_cc_restriction_c_flw(crse, cc, fine,cf, ir, nc,alfa)
    use cc_restriction_module
    type(multifab), intent(inout) :: fine
    type(multifab), intent(inout) :: crse
    integer, intent(in)           :: cc, cf, ir(:)
    integer, intent(in), optional :: nc
    real(dp_t),intent(in)         :: alfa

    integer             :: i, n, lnc, dm
    integer :: lo(get_dim(fine)), hi(get_dim(fine)), loc(get_dim(fine)), lof(get_dim(fine))
    real(dp_t), pointer :: fp(:,:,:,:), cp(:,:,:,:),pflw(:,:,:,:)
    type(layout)        :: lacfine,laf
    type(multifab)      :: cfine
    type(mfiter)        :: mfi
    type(box)           :: bx

    lnc = 1; if ( present(nc) ) lnc = nc

    laf = get_layout(fine)

    call layout_build_coarse(lacfine, laf, ir)

    call multifab_build(cfine, lacfine, nc = lnc, ng = 0)

    dm = get_dim(cfine)

    !$omp parallel private(mfi,n,i,bx,loc,lof,lo,hi,fp,cp)
    call mfiter_build(mfi,cfine,.true.)
    do n = 1, lnc
       do while(next_tile(mfi,i))
          bx = get_tilebox(mfi)

          lo  = lwb(bx)
          hi  = upb(bx)
          loc = lwb(get_pbox(cfine, i))
          lof = lwb(get_pbox( fine, i))

          fp    => dataptr( fine,   i, n+cf-1, 1)
          cp => dataptr(cfine, i, n,      1)
          select case (dm)
          case (1)
             call cc_restriction_1d(cp(:,1,1,1), loc, fp(:,1,1,1), lof, lo, hi, ir)
          case (2)
             call cc_restriction_2d(cp(:,:,1,1), loc, fp(:,:,1,1), lof, lo, hi, ir)
          case (3)
             call cc_restriction_3d(cp(:,:,:,1), loc, fp(:,:,:,1), lof, lo, hi, ir)
          end select
       end do
    end do
    !$omp end parallel

    do i = 1, nfabs(cfine)
       pflw  => dataptr(cfine, i)

       select case(dm)
       case (2)
           call modify_f_values_flw_2d_c(pflw(:,:,1,1),pflw(:,:,1,2),pflw(:,:,1,3),pflw(:,:,1,4:12),alfa)
       case (3)
           call modify_f_values_flw_3d_c(pflw(:,:,:,1),pflw(:,:,:,2),pflw(:,:,:,3),pflw(:,:,:,4),&
                                        pflw(:,:,:,5:23),alfa)
       end select
    end do

    call copy(crse, cc, cfine, 1, lnc)

    call multifab_destroy(cfine)

  end subroutine ml_cc_restriction_c_flw

  subroutine modify_f_values_flw_2d_c(ux,uy,rho,f_mat,alfa)
  double precision  :: ux(:,:),uy(:,:),rho(:,:)
  double precision  :: f_mat(:,:,:)
  double precision  :: alfa

  integer           :: i, j, iter, nx, ny
  double precision  :: var_stable = 0.d0,var_1,var_2

  integer,parameter :: q=9
  double precision, parameter :: p0_2d = 4.d0/9.d0, p1_2d = 1.d0/9.d0, p2_2d=1.d0/36.d0
  integer           :: ex(q)=(/ 0,1,0,-1,0,1,-1,-1,1 /),ey(q)=(/ 0,0,1,0,-1,1,1,-1,-1 /)
  double precision  :: w(q)=(/ p0_2d,p1_2d,p1_2d,p1_2d,p1_2d,p2_2d,p2_2d,p2_2d,p2_2d /) 
    
  nx = size(ux,1)
  ny = size(ux,2)

  do j = 1, ny
     do i = 1, nx
        do iter=1,q

           var_1 = ex(iter)*ux(i,j)+ey(iter)*uy(i,j)
           var_2 = ux(i,j)**2+uy(i,j)**2

           var_stable = w(iter)*rho(i,j)*(1.0+3.0*var_1+4.5*var_1**2-1.5*var_2) 

           f_mat(i,j,iter) = alfa*f_mat(i,j,iter) + (1.d0-alfa)*var_stable

        end do
     end do
  end do

end subroutine modify_f_values_flw_2d_c

subroutine modify_f_values_flw_3d_c(ux,uy,uz,rho,f_mat,alfa)
  double precision  :: ux(:,:,:),uy(:,:,:),uz(:,:,:),rho(:,:,:)
  double precision  :: f_mat(:,:,:,:)
  double precision  :: alfa
     
  integer           :: i, j, k, iter, nx, ny, nz

  ! local variables
  integer,parameter :: q=19
  double precision, parameter :: p0_3d = 1.d0/3.d0, p1_3d = 1.d0/18.d0, p2_3d=1.d0/36.d0
  integer           :: ex(q)=(/ 0,1,-1,0,0,0,0,1,1,-1,-1,1,-1,1,-1,0,0,0,0 /)
  integer           :: ey(q)=(/ 0,0,0,1,-1,0,0,1,-1,1,-1,0,0,0,0,1,1,-1,-1 /)
  integer           :: ez(q)=(/ 0,0,0,0,0,1,-1,0,0,0,0,1,1,-1,-1,1,-1,1,-1 /)
  double precision  :: w(q)=(/ p0_3d,p1_3d,p1_3d,p1_3d,p1_3d,p1_3d,p1_3d,p2_3d,p2_3d,p2_3d,&
                           p2_3d,p2_3d,p2_3d,p2_3d,p2_3d,p2_3d,p2_3d,p2_3d,p2_3d /)  

  double precision  :: var_stable = 0.d0,var_1 = 0.d0,var_2 = 0.d0
    
  nx = size(ux,1)
  ny = size(uy,2)
  nz = size(uz,3)

  do k = 1, nz
     do j = 1, ny
        do i = 1, nx
           do iter=1,q

              var_1 = ex(iter)*ux(i,j,k)+ey(iter)*uy(i,j,k)+ez(iter)*uz(i,j,k)
              var_2 = ux(i,j,k)**2+uy(i,j,k)**2+uz(i,j,k)**2
              var_stable = w(iter)*rho(i,j,k)*(1.0+3.0*var_1+4.5*var_1**2-1.5*var_2) 

              f_mat(i,j,k,iter) = alfa*f_mat(i,j,k,iter) + (1.d0-alfa)*var_stable

           end do
        end do
     end do
  end do
    
end subroutine modify_f_values_flw_3d_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ml_cc_restriction(crse, fine, ir)
    type(multifab), intent(inout) :: fine
    type(multifab), intent(inout) :: crse
    integer,        intent(in)    :: ir(:)
    if ( ncomp(crse) .ne. ncomp(fine) ) then
       call bl_error('ml_cc_restriction: crse & fine must have same # of components')
    end if
    call ml_cc_restriction_c(crse, 1, fine, 1, ir, ncomp(crse))
  end subroutine ml_cc_restriction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ml_edge_restriction_c(crse, cc, fine, cf, ir, face, nc)
    use edge_restriction_module
    type(multifab), intent(inout) :: fine
    type(multifab), intent(inout) :: crse
    integer,        intent(in)    :: cc, cf, ir(:)
    integer,        intent(in)    :: face
    integer, intent(in), optional :: nc

    integer             :: i, n, lnc, dm, len
    integer             :: lo(get_dim(fine)), hi(get_dim(fine)), loc(get_dim(fine)), lof(get_dim(fine))
    real(dp_t), pointer :: fp(:,:,:,:), cp(:,:,:,:)
    type(box)           :: bx,fine_domain,crse_domain
    type(layout)        :: lacfine, lacfine_lo, lacfine_hi 
    type(layout)        :: la_lo,la_hi,laf
    type(multifab)      :: cfine, fine_lo, fine_hi
    type(list_box)      :: bxs_lo,bxs_hi
    type(boxarray)      :: ba_lo,ba_hi
    logical             :: nodal(get_dim(fine)), pmask(get_dim(fine))
    type(mfiter)        :: mfi

    dm = get_dim(crse)

    lnc = 1; if ( present(nc) ) lnc = nc

    laf = get_layout(fine)

    call layout_build_coarse(lacfine, laf, ir)

    call multifab_build(cfine, lacfine, nc = ncomp(crse), ng = 0, nodal = nodal_flags(crse))

    !$omp parallel private(mfi,n,i,bx,lo,hi,loc,lof,fp,cp)
    call mfiter_build(mfi,cfine,.true.)
    do n = 1, lnc
       do while(next_tile(mfi,i))
          bx = get_tilebox(mfi)

          lo  = lwb(bx)
          hi  = upb(bx)
          loc = lwb(get_pbox(cfine,i))
          lof = lwb(get_pbox(fine, i))

          fp  => dataptr(fine,  i, n+cf-1, 1)
          cp  => dataptr(cfine, i, n,      1)
          select case (dm)
          case (1)
             call edge_restriction_1d(cp(:,1,1,1), loc, fp(:,1,1,1), lof, lo, hi, ir)
          case (2)
             call edge_restriction_2d(cp(:,:,1,1), loc, fp(:,:,1,1), lof, lo, hi, ir, face)
          case (3)
             call edge_restriction_3d(cp(:,:,:,1), loc, fp(:,:,:,1), lof, lo, hi, ir, face)
          end select
       enddo
    end do
    !$omp end parallel

    call copy(crse, cc, cfine, 1, lnc)

    call multifab_destroy(cfine)
    !
    ! Now do periodic fix-up if necessary.
    !
    pmask = get_pmask(get_layout(crse))

    if (pmask(face)) then

       fine_domain = get_pd(get_layout(fine))
       crse_domain = get_pd(get_layout(crse))
       nodal(:)    = .false.
       nodal(face) = .true.
       len         = box_extent_d(fine_domain,face)
       !
       ! First copy from lo edges to hi edges.
       !
       do i = 1, nboxes(fine%la)
          bx = get_box(fine%la,i)
          if (bx%lo(face) == fine_domain%lo(face)) then
             bx = shift(bx, len, face)
             bx = intersection(bx,grow(fine_domain,1,face,+1))
             call push_back(bxs_lo,bx)
          end if
       end do

       if (.not. empty(bxs_lo)) then

          call boxarray_build_l(ba_lo,bxs_lo,sort=.false.)
          call destroy(bxs_lo)
          call layout_build_ba(la_lo,ba_lo,fine_domain,pmask)
          call boxarray_destroy(ba_lo)
          call multifab_build(fine_lo, la_lo, nc = ncomp(fine), ng = 0, nodal = nodal)
   
          call multifab_copy_on_shift(fine_lo, 1, fine, cf, lnc, len, face)

          call layout_build_coarse(lacfine_lo, la_lo, ir)
          call multifab_build(cfine, lacfine_lo, nc = ncomp(crse), ng = 0, nodal = nodal_flags(crse))

          !$OMP PARALLEL DO PRIVATE(i,n,lo,hi,loc,lof,fp,cp)
          do i = 1, nfabs(fine_lo)
             lo  = lwb(get_ibox(cfine,i))
             hi  = upb(get_ibox(cfine,i))
             hi(face) = lo(face)
             loc = lwb(get_pbox(cfine,i))
             lof = lwb(get_pbox(fine_lo, i))
             do n = 1, lnc
                fp  => dataptr(fine_lo, i, n, 1)
                cp  => dataptr(cfine     , i, n, 1)
                select case (dm)
                case (1)
                   call edge_restriction_1d(cp(:,1,1,1), loc, fp(:,1,1,1), lof, lo, hi, ir)
                case (2)
                   call edge_restriction_2d(cp(:,:,1,1), loc, fp(:,:,1,1), lof, lo, hi, ir, face)
                case (3)
                   call edge_restriction_3d(cp(:,:,:,1), loc, fp(:,:,:,1), lof, lo, hi, ir, face)
                end select
             enddo
          end do
          !$OMP END PARALLEL DO
   
          call copy(crse, cc, cfine, 1, lnc)

          call multifab_destroy(cfine)
          call multifab_destroy(fine_lo)
          call layout_destroy(la_lo)
       
       end if
       !
       ! Next copy from hi edges to lo edges.
       !
       do i = 1, nboxes(fine%la)
          bx = get_box(fine%la,i)
          if (bx%hi(face) == fine_domain%hi(face)) then
             bx = shift(bx, -len, face)
             bx = intersection(bx,grow(fine_domain,1,face,-1))
             call push_back(bxs_hi,bx)
          end if
       end do

       if (.not. empty(bxs_hi)) then

          call boxarray_build_l(ba_hi,bxs_hi,sort=.false.)
          call destroy(bxs_hi)
          call layout_build_ba(la_hi,ba_hi,fine_domain,pmask)
          call boxarray_destroy(ba_hi)
          call multifab_build(fine_hi, la_hi, nc = ncomp(fine), ng = 0, nodal = nodal)
   
          call multifab_copy_on_shift(fine_hi, 1, fine, cf, lnc, -len, face)

          call layout_build_coarse(lacfine_hi, la_hi, ir)
          call multifab_build(cfine, lacfine_hi, nc = ncomp(crse), ng = 0, nodal = nodal_flags(crse))

          !$OMP PARALLEL DO PRIVATE(i,n,lo,hi,loc,lof,fp,cp)
          do i = 1, nfabs(fine_hi)
             lo  = lwb(get_ibox(cfine,i))
             hi  = upb(get_ibox(cfine,i))
             lo(face) = hi(face)
             loc = lwb(get_pbox(cfine,i))
             lof = lwb(get_pbox(fine_hi, i))
             do n = 1, lnc
                fp  => dataptr(fine_hi, i, n, 1)
                cp  => dataptr(cfine     , i, n, 1)
                select case (dm)
                case (1)
                   call edge_restriction_1d(cp(:,1,1,1), loc, fp(:,1,1,1), lof, lo, hi, ir)
                case (2)
                   call edge_restriction_2d(cp(:,:,1,1), loc, fp(:,:,1,1), lof, lo, hi, ir, face)
                case (3)
                   call edge_restriction_3d(cp(:,:,:,1), loc, fp(:,:,:,1), lof, lo, hi, ir, face)
                end select
             enddo
          end do
          !$OMP END PARALLEL DO

          call copy(crse, cc, cfine, 1, lnc)

          call multifab_destroy(cfine)
          call multifab_destroy(fine_hi)
          call layout_destroy(la_hi)

       end if ! .not. empty

    end if ! pmask(face)

  end subroutine ml_edge_restriction_c

  subroutine modify_f_values_flw_2d(ux,uy,rho,f_mat,alfa)
  double precision  :: ux(:,:),uy(:,:),rho(:,:)
  double precision  :: f_mat(:,:,:)
  double precision  :: alfa

  integer           :: i, j, iter, nx, ny
  double precision  :: var_stable = 0.d0,var_1,var_2

  integer,parameter :: q=9
  double precision, parameter :: p0_2d = 4.d0/9.d0, p1_2d = 1.d0/9.d0, p2_2d=1.d0/36.d0
  integer           :: ex(q)=(/ 0,1,0,-1,0,1,-1,-1,1 /),ey(q)=(/ 0,0,1,0,-1,1,1,-1,-1 /)
  double precision  :: w(q)=(/ p0_2d,p1_2d,p1_2d,p1_2d,p1_2d,p2_2d,p2_2d,p2_2d,p2_2d /) 
    
  nx = size(ux,1)
  ny = size(ux,2)

  do j = 1, ny
     do i = 1, nx
 
        do iter=1,q

           var_1 = ex(iter)*ux(i,j)+ey(iter)*uy(i,j)
           var_2 = ux(i,j)**2+uy(i,j)**2

           var_stable = w(iter)*rho(i,j)*(1.0+3.0*var_1+4.5*var_1**2-1.5*var_2) 

           f_mat(i,j,iter) = alfa*f_mat(i,j,iter) + (1.d0-alfa)*var_stable

        end do
     end do
  end do

end subroutine modify_f_values_flw_2d

subroutine modify_f_values_flw_3d(ux,uy,uz,rho,f_mat,alfa)
  double precision  :: ux(:,:,:),uy(:,:,:),uz(:,:,:),rho(:,:,:)
  double precision  :: f_mat(:,:,:,:)
  double precision  :: alfa
     
  integer           :: i, j, k, iter, nx, ny, nz

  ! local variables
  integer,parameter :: q=19
  double precision, parameter :: p0_3d = 1.d0/3.d0, p1_3d = 1.d0/18.d0, p2_3d=1.d0/36.d0
  integer           :: ex(q)=(/ 0,1,-1,0,0,0,0,1,1,-1,-1,1,-1,1,-1,0,0,0,0 /)
  integer           :: ey(q)=(/ 0,0,0,1,-1,0,0,1,-1,1,-1,0,0,0,0,1,1,-1,-1 /)
  integer           :: ez(q)=(/ 0,0,0,0,0,1,-1,0,0,0,0,1,1,-1,-1,1,-1,1,-1 /)
  double precision  :: w(q)=(/ p0_3d,p1_3d,p1_3d,p1_3d,p1_3d,p1_3d,p1_3d,p2_3d,p2_3d,p2_3d,&
                           p2_3d,p2_3d,p2_3d,p2_3d,p2_3d,p2_3d,p2_3d,p2_3d,p2_3d /)  

  double precision  :: var_stable = 0.d0,var_1 = 0.d0,var_2 = 0.d0
    
  nx = size(ux,1)
  ny = size(uy,2)
  nz = size(uz,3)

  do k = 1, nz
     do j = 1, ny
        do i = 1, nx

           do iter=1,q

              var_1 = ex(iter)*ux(i,j,k)+ey(iter)*uy(i,j,k)+ez(iter)*uz(i,j,k)
              var_2 = ux(i,j,k)**2+uy(i,j,k)**2+uz(i,j,k)**2
              var_stable = w(iter)*rho(i,j,k)*(1.0+3.0*var_1+4.5*var_1**2-1.5*var_2) 

              f_mat(i,j,k,iter) = alfa*f_mat(i,j,k,iter) + (1.d0-alfa)*var_stable

           end do
        end do
     end do
  end do
    
end subroutine modify_f_values_flw_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ml_edge_restriction(crse, fine, ir, face)
    type(multifab), intent(inout) :: fine
    type(multifab), intent(inout) :: crse
    integer,        intent(in)    :: ir(:)
    integer,        intent(in)    :: face

    if ( ncomp(crse) .ne. ncomp(fine) ) then
       call bl_error('ml_edge_restriction: crse & fine must have same # of components')
    end if
    call ml_edge_restriction_c(crse, 1, fine, 1, ir, face, ncomp(crse))

  end subroutine ml_edge_restriction

end module ml_cc_restriction_module
