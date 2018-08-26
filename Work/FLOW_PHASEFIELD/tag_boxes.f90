module tag_boxes_module

  use multifab_module
  use bl_error_module

  implicit none

contains


  subroutine tag_boxes(mf,tagboxes,dx,lev,para_uc,para_th,amr_thers,aux_tag_mf)

    type( multifab)         , intent(in   ) :: mf
    type(lmultifab)         , intent(inout) :: tagboxes
    real(dp_t)              , intent(in   ) :: dx,para_uc,para_th,amr_thers
    integer                 , intent(in   ) :: lev
    type(multifab), optional, intent(in   ) :: aux_tag_mf
    ! aux_tag_mf allows user to pass in additional multifabs for tagging logic

    ! local variables
    real(kind = dp_t), pointer :: mfp(:,:,:,:)
    logical          , pointer :: tp(:,:,:,:)
    integer           :: i, lo(get_dim(mf)), hi(get_dim(mf)), ng

    if (present(aux_tag_mf)) then
       call bl_error("tag_boxes.f90: aux_tag_mf passed to tag_boxes without implementation")
    end if

    ng = nghost(mf)

    do i = 1, nfabs(mf)
       mfp => dataptr(mf, i)
       tp  => dataptr(tagboxes, i)
       lo =  lwb(get_box(tagboxes, i))
       hi =  upb(get_box(tagboxes, i))
       select case (get_dim(mf))
       case (2)
          call tag_boxes_2d(tp(:,:,1,1),mfp(:,:,1,1),mfp(:,:,1,2),mfp(:,:,1,3),lo,hi,ng,dx,lev,&
                            para_uc,para_th,amr_thers)
       case  (3)
          call tag_boxes_3d(tp(:,:,:,1),mfp(:,:,:,1),mfp(:,:,:,2),mfp(:,:,:,3),lo,hi,ng,dx,lev,&
                            para_uc,para_th,amr_thers)
       end select
    end do

  end subroutine tag_boxes

  subroutine tag_boxes_2d(tagbox,mf_pf,mf_uc,mf_th,lo,hi,ng,dx,lev,para_uc,para_th,amr_thers)

    integer          , intent(in   ) :: lo(2),hi(2),ng
    logical          , intent(  out) :: tagbox(lo(1):hi(1),lo(2):hi(2))
    real(kind = dp_t), intent(in   ) :: mf_pf(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
    real(kind = dp_t), intent(in   ) :: mf_uc(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
    real(kind = dp_t), intent(in   ) :: mf_th(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
    real(dp_t)       , intent(in   ) :: dx,para_uc,para_th,amr_thers
    integer          , intent(in   ) :: lev

    ! local variables
    integer :: i,j
    real(kind = dp_t)  :: var_e,gra_x,gra_y,cri_thers

    ! initially say that we do not want to tag any cells for refinement
    tagbox = .false.
    !cri_thers = amr_thers * ( 4.d0**(lev-1) )
    cri_thers = amr_thers

    ! tag all boxes where the first component of mf >= 1.01
    do j = lo(2),hi(2)
        do i = lo(1),hi(1)
             gra_x = ( mf_pf(i+1,j) - mf_pf(i-1,j) ) / 2.d0
             gra_y = ( mf_pf(i,j+1) - mf_pf(i,j-1) ) / 2.d0
             var_e = sqrt( gra_x**2 + gra_y**2 )

             gra_x = ( mf_uc(i+1,j) - mf_uc(i-1,j) ) / 2.d0
             gra_y = ( mf_uc(i,j+1) - mf_uc(i,j-1) ) / 2.d0
             var_e = var_e + para_uc * sqrt( gra_x**2 + gra_y**2 )

             gra_x = ( mf_th(i+1,j) - mf_th(i-1,j) ) / 2.d0
             gra_y = ( mf_th(i,j+1) - mf_th(i,j-1) ) / 2.d0
             var_e = var_e + para_th * sqrt( gra_x**2 + gra_y**2 )

             !var_e = var_e * 2.d0 * dx * dx
             !var_e = var_e

             if (var_e .ge. cri_thers) then
                tagbox(i,j) = .true.
             end if

        end do
    end do

  end subroutine tag_boxes_2d

  subroutine tag_boxes_3d(tagbox,mf_pf,mf_uc,mf_th,lo,hi,ng,dx,lev,para_uc,para_th,amr_thers)

    integer          , intent(in   ) :: lo(3),hi(3),ng
    logical          , intent(  out) :: tagbox(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    real(kind = dp_t), intent(in   ) :: mf_pf(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
    real(kind = dp_t), intent(in   ) :: mf_uc(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
    real(kind = dp_t), intent(in   ) :: mf_th(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
    real(dp_t)       , intent(in   ) :: dx,para_uc,para_th,amr_thers
    integer          , intent(in   ) :: lev

    ! local variables
    integer :: i,j,k
    real(kind = dp_t)  :: var_e,gra_x,gra_y,gra_z,cri_thers

    ! initially say that we do not want to tag any cells for refinement
    tagbox = .false.
    !cri_thers = amr_thers * ( 4.d0**(lev-1) )
    cri_thers = amr_thers

    ! tag all boxes where the first component of mf >= 1.01
    do k = lo(3),hi(3)
        do j = lo(2),hi(2)
            do i = lo(1),hi(1)
                gra_x = ( mf_pf(i+1,j,k) - mf_pf(i-1,j,k) ) / 2.d0
                gra_y = ( mf_pf(i,j+1,k) - mf_pf(i,j-1,k) ) / 2.d0
                gra_z = ( mf_pf(i,j,k+1) - mf_pf(i,j,k-1) ) / 2.d0
                var_e = sqrt( gra_x**2 + gra_y**2 + gra_z**2 )

                gra_x = ( mf_uc(i+1,j,k) - mf_uc(i-1,j,k) ) / 2.d0
                gra_y = ( mf_uc(i,j+1,k) - mf_uc(i,j-1,k) ) / 2.d0
                gra_z = ( mf_uc(i,j,k+1) - mf_uc(i,j,k-1) ) / 2.d0
                var_e = var_e + para_uc * sqrt( gra_x**2 + gra_y**2 + gra_z**2 )

                gra_x = ( mf_th(i+1,j,k) - mf_th(i-1,j,k) ) / 2.d0
                gra_y = ( mf_th(i,j+1,k) - mf_th(i,j-1,k) ) / 2.d0
                gra_z = ( mf_th(i,j,k+1) - mf_th(i,j,k-1) ) / 2.d0
                var_e = var_e + para_th * sqrt( gra_x**2 + gra_y**2 + gra_z**2 )

                !var_e = var_e * 3.d0 * dx * dx
                !var_e = var_e

                if (var_e .ge. cri_thers) then
                   tagbox(i,j,k) = .true.
                end if

            end do
        end do
    end do

  end subroutine tag_boxes_3d

  !
  ! The *_init routines are for the initialization of the seed related tag
  ! First we label these nearby to the seed center, as specified by safe_n_init
  ! Then do the calculation for some time, normally 1000 steps with dt = 0.1*
  ! aiming at achieving a quasi-steady solute state, then the computation is
  ! proceeded by regriding and advances
  !
  subroutine tag_boxes_init(mf,seed_pos,seed_r,seed_num,tagboxes,dx,prob_lo,prob_hi,safe_n_init,&
                            seed_type, aux_tag_mf)

    type( multifab)         , intent(in   ) :: mf
    integer                 , intent(in   ) :: seed_num,seed_type
    real(kind=dp_t)         , intent(in   ) :: prob_lo(mf%dim),prob_hi(mf%dim)
    real(dp_t)              , intent(in   ) :: seed_pos(seed_num,mf%dim)
    real(dp_t)              , intent(in   ) :: seed_r,safe_n_init
    type(lmultifab)         , intent(inout) :: tagboxes
    real(dp_t)              , intent(in   ) :: dx
    type(multifab), optional, intent(in   ) :: aux_tag_mf
    ! aux_tag_mf allows user to pass in additional multifabs for tagging logic

    ! local variables
    real(kind = dp_t), pointer :: mfp(:,:,:,:)
    logical          , pointer :: tp(:,:,:,:)
    integer           :: i, lo(get_dim(mf)), hi(get_dim(mf)), ng

    if (present(aux_tag_mf)) then
       call bl_error("tag_boxes.f90: aux_tag_mf passed to tag_boxes without implementation")
    end if

    do i = 1, nfabs(mf)
       mfp => dataptr(mf, i)
       tp  => dataptr(tagboxes, i)
       lo =  lwb(get_box(tagboxes, i))
       hi =  upb(get_box(tagboxes, i))
       select case (get_dim(mf))
       case (2)
          call tag_boxes_2d_init(tp(:,:,1,1),seed_pos,seed_r,seed_num,lo,hi,prob_lo,&
                                 prob_hi,dx,safe_n_init,seed_type)
       case  (3)
          call tag_boxes_3d_init(tp(:,:,:,1),seed_pos,seed_r,seed_num,lo,hi,prob_lo,&
                                 prob_hi,dx,safe_n_init,seed_type)
       end select
    end do

  end subroutine tag_boxes_init

  subroutine tag_boxes_2d_init(tagbox,seed_pos,seed_r,seed_num,lo,hi,prob_lo,prob_hi,&
                               dx,safe_n_init,seed_type)

    integer          , intent(in   ) :: lo(2),hi(2)
    real(dp_t)       , intent(in   ) :: prob_lo(2),prob_hi(2)
    integer          , intent(in   ) :: seed_num,seed_type
    logical          , intent(  out) :: tagbox(lo(1):,lo(2):)
    real(dp_t)       , intent(in   ) :: seed_pos(seed_num,2)
    real(dp_t)       , intent(in   ) :: seed_r
    real(dp_t)       , intent(in   ) :: dx,safe_n_init

    ! local variables
    integer :: i,j,s_i
    real(kind = dp_t)  :: dis,safe_sr,y,x,l_m,l_mx,x_l,y_l,l_m2,l_mx2

    y_l = prob_hi(2) - prob_lo(2)
    x_l = prob_hi(1) - prob_lo(1)

    ! initially say that we do not want to tag any cells for refinement
    tagbox = .false.
    !safe_sr = seed_r + safe_n_init * dx
    safe_sr = safe_n_init * dx

    select case(seed_type)
    case(1)  ! columnar
       ! tag all boxes where the first component of mf >= 1.01
       do j = lo(2),hi(2)
          y =  prob_lo(2) +(dble(j)+0.5d0) * dx
          l_m = abs(y - seed_r)
          do i = lo(1),hi(1)
             if( l_m .le. safe_sr ) then
                tagbox(i,j) = .true.
             end if

          end do
       end do

    case(2)
      ! tag all boxes where the first component of mf >= 1.01
       do j = lo(2),hi(2)
          y =  prob_lo(2) +(dble(j)+0.5d0) * dx
          l_m  = abs(y - seed_r)
          l_m2 = abs(y - (y_l-seed_r) )
          do i = lo(1),hi(1)
            x = prob_lo(1) + (dble(i)+0.5d0) * dx
            l_mx  = abs(x - seed_r)
            l_mx2 = abs(x - (x_l-seed_r) )
             if( (l_m .le. safe_sr)  .or. (l_mx .le. safe_sr) .or. &
                 (l_m2 .le. safe_sr) .or. (l_mx2 .le. safe_sr)  ) then
                tagbox(i,j) = .true.
             end if

          end do
       end do

    case default ! Equiaxed
       ! tag all boxes where the first component of mf >= 1.01
       do j = lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+0.5d0) * dx
          do i = lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0) * dx

             do s_i = 1,seed_num
                dis = sqrt( (x - seed_pos(s_i,1))**2 + (y - seed_pos(s_i,2))**2)
                dis = abs(dis - seed_r)

                if( dis .le. safe_sr ) then
                   tagbox(i,j) = .true.
                end if

             end do

           end do
       end do

    end select



  end subroutine tag_boxes_2d_init

  subroutine tag_boxes_3d_init(tagbox,seed_pos,seed_r,seed_num,lo,hi,prob_lo,prob_hi,&
                               dx,safe_n_init,seed_type)

    integer          , intent(in   ) :: lo(3),hi(3)
    real(dp_t)       , intent(in   ) :: prob_lo(3),prob_hi(3)
    integer          , intent(in   ) :: seed_num,seed_type
    logical          , intent(  out) :: tagbox(lo(1):,lo(2):,lo(3):)
    real(dp_t)       , intent(in   ) :: seed_pos(seed_num,3)
    real(dp_t)       , intent(in   ) :: seed_r
    real(dp_t)       , intent(in   ) :: dx,safe_n_init

    ! local variables
    integer :: i,j,k,s_i
    real(kind = dp_t)  :: dis,safe_sr,y,x,z, l_measure
    double precision :: x_left, x_right, y_top, y_bottom

    ! initially say that we do not want to tag any cells for refinement
    tagbox = .false.
    !safe_sr = seed_r + safe_n_init * dx
    safe_sr = safe_n_init * dx

    select case(seed_type)
    case(1)  ! columnar
       ! tag all boxes where the first component of mf >= 1.01
       do k = lo(3),hi(3)
           do j = lo(2),hi(2)
               y = (dble(j)+0.5d0) * dx
               do i = lo(1),hi(1)
                  
                  l_measure = abs(y - seed_r)

                  if( l_measure .le. safe_sr ) then
                     tagbox(i,j,k) = .true.
                   end if

               end do
           end do
       end do

    case (2)
       ! tag all boxes where the first component of mf >= 1.01
       do k = lo(3),hi(3)
           !z = prob_lo(3) + (dble(k)+0.5d0) * dx
           do j = lo(2),hi(2)
               y = prob_lo(2) + (dble(j)+0.5d0) * dx
               do i = lo(1),hi(1)
                  x = prob_lo(1) + (dble(i)+0.5d0) * dx

                  do s_i = 1,seed_num
                     dis = sqrt( (x - seed_pos(s_i,1))**2 + (y - seed_pos(s_i,2))**2 )
                     dis = abs(dis - seed_r)

                     if( dis .le. safe_sr ) then
                        tagbox(i,j,k) = .true.
                     end if

                  end do

               end do
           end do
       end do

    case (3)
     do k = lo(3), hi(3)
       !z = prob_lo(3) + (dble(k)+0.5d0) * dx
       do j = lo(2), hi(2)
         y = prob_lo(2) + (dble(j)+0.5d0) * dx
         do i = lo(1), hi(1)
            x = prob_lo(1) + (dble(i)+0.5d0) * dx

            do s_i = 1,seed_num

               x_left   = seed_pos(s_i,1) - seed_r - safe_sr
               x_right  = seed_pos(s_i,1) + seed_r + safe_sr
               y_bottom = seed_pos(s_i,2) - seed_r - safe_sr 
               y_top    = seed_pos(s_i,2) + seed_r + safe_sr

               if ( (x >= x_left   .and. x <= x_right) .and. &
                    (y >= y_bottom .and. y <= y_top  ) ) then
                  tagbox(i,j,k) = .true.
               end if

               x_left   = seed_pos(s_i,1) - seed_r + safe_sr
               x_right  = seed_pos(s_i,1) + seed_r - safe_sr
               y_bottom = seed_pos(s_i,2) - seed_r + safe_sr 
               y_top    = seed_pos(s_i,2) + seed_r - safe_sr

               if ( (x >= x_left   .and. x <= x_right) .and. &
                    (y >= y_bottom .and. y <= y_top  ) ) then
                  tagbox(i,j,k) = .false.
               end if

            end do

          end do  !i
        end do  !j
     end do  !k

    case default
       ! tag all boxes where the first component of mf >= 1.01
       do k = lo(3),hi(3)
           z = prob_lo(3) + (dble(k)+0.5d0) * dx
           do j = lo(2),hi(2)
               y = prob_lo(2) + (dble(j)+0.5d0) * dx
               do i = lo(1),hi(1)
                  x = prob_lo(1) + (dble(i)+0.5d0) * dx

                  do s_i = 1,seed_num
                     dis = sqrt( (x - seed_pos(s_i,1))**2 + (y - seed_pos(s_i,2))**2 + (z - seed_pos(s_i,3))**2)
                     dis = abs(dis - seed_r)

                     if( dis .le. safe_sr ) then
                        tagbox(i,j,k) = .true.
                     end if

                  end do

               end do
           end do
       end do

    end select



  end subroutine tag_boxes_3d_init

  subroutine tag_boxes_flw(flw,mf,tagboxes,dx,lev,para_uc,para_th,para_u,amr_thers,&
                           multi_add_tau,aux_tag_mf)

    type( multifab)         , intent(in   ) :: flw,mf
    type(lmultifab)         , intent(inout) :: tagboxes
    real(dp_t)              , intent(in   ) :: dx,para_uc,para_th,amr_thers,para_u
    real(dp_t)              , intent(in   ) :: multi_add_tau
    integer                 , intent(in   ) :: lev
    type(multifab), optional, intent(in   ) :: aux_tag_mf
    ! aux_tag_mf allows user to pass in additional multifabs for tagging logic

    ! local variables
    real(kind = dp_t), pointer :: mfp(:,:,:,:)
    real(kind = dp_t), pointer :: flp(:,:,:,:)
    logical          , pointer :: tp(:,:,:,:)
    integer           :: i, lo(get_dim(mf)), hi(get_dim(mf)), ng

    if (present(aux_tag_mf)) then
       call bl_error("tag_boxes.f90: aux_tag_mf passed to tag_boxes without implementation")
    end if

    ng = nghost(mf)

    do i = 1, nfabs(mf)
       mfp => dataptr(mf, i)
       flp => dataptr(flw, i)
       tp  => dataptr(tagboxes, i)
       lo =  lwb(get_box(tagboxes, i))
       hi =  upb(get_box(tagboxes, i))
       select case (get_dim(mf))
       case (2)
          call tag_boxes_2d_flw(tp(:,:,1,1),flp(:,:,1,1:2),mfp(:,:,1,1),mfp(:,:,1,2),mfp(:,:,1,3),lo,hi,ng,dx,lev,&
                            para_uc,para_th,para_u,amr_thers,multi_add_tau)
       case  (3)
          call tag_boxes_3d_flw(tp(:,:,:,1),flp(:,:,:,1:3),mfp(:,:,:,1),mfp(:,:,:,2),mfp(:,:,:,3),lo,hi,ng,dx,lev,&
                            para_uc,para_th,para_u,amr_thers,multi_add_tau)
       end select
    end do

  end subroutine tag_boxes_flw

  subroutine tag_boxes_2d_flw(tagbox,u_vec,mf_pf,mf_uc,mf_th,lo,hi,ng,dx,lev,para_uc,para_th,para_u,amr_thers,&
                              multi_add_tau)

    integer          , intent(in   ) :: lo(2),hi(2),ng
    logical          , intent(  out) :: tagbox(lo(1):hi(1),lo(2):hi(2))
    real(kind = dp_t), intent(in   ) :: u_vec(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,1:2)
    real(kind = dp_t), intent(in   ) :: mf_pf(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
    real(kind = dp_t), intent(in   ) :: mf_uc(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
    real(kind = dp_t), intent(in   ) :: mf_th(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng)
    real(dp_t)       , intent(in   ) :: dx,para_uc,para_th,amr_thers,para_u,multi_add_tau
    integer          , intent(in   ) :: lev

    ! local variables
    integer :: i,j
    real(kind = dp_t)  :: var_e,gra_x,gra_y,cri_thers

    ! initially say that we do not want to tag any cells for refinement
    tagbox = .false.
    !cri_thers = amr_thers * ( 4.d0**(lev-1) )
    cri_thers = amr_thers

    ! tag all boxes where the first component of mf >= 1.01
    do j = lo(2),hi(2)
        do i = lo(1),hi(1)
             gra_x = ( mf_pf(i+1,j) - mf_pf(i-1,j) ) / 2.d0
             gra_y = ( mf_pf(i,j+1) - mf_pf(i,j-1) ) / 2.d0
             var_e = sqrt( gra_x**2 + gra_y**2 )

             gra_x = ( mf_uc(i+1,j) - mf_uc(i-1,j) ) / 2.d0
             gra_y = ( mf_uc(i,j+1) - mf_uc(i,j-1) ) / 2.d0
             var_e = var_e + para_uc * sqrt( gra_x**2 + gra_y**2 )

             gra_x = ( mf_th(i+1,j) - mf_th(i-1,j) ) / 2.d0
             gra_y = ( mf_th(i,j+1) - mf_th(i,j-1) ) / 2.d0
             var_e = var_e + para_th * sqrt( gra_x**2 + gra_y**2 )

             gra_x = ( u_vec(i+1,j,1) - u_vec(i-1,j,1) ) / 2.d0
             gra_y = ( u_vec(i,j+1,2) - u_vec(i,j-1,2) ) / 2.d0
             var_e = var_e + para_u * sqrt( gra_x**2 + gra_y**2 )

             !var_e = var_e * 2.d0 * dx * dx
             !var_e = var_e * dx

             if (var_e .ge. cri_thers) then
                tagbox(i,j) = .true.
             end if

        end do
    end do

  end subroutine tag_boxes_2d_flw

  subroutine tag_boxes_3d_flw(tagbox,u_vec,mf_pf,mf_uc,mf_th,lo,hi,ng,dx,lev,para_uc,para_th,para_u,amr_thers,&
                              multi_add_tau)

    integer          , intent(in   ) :: lo(3),hi(3),ng
    logical          , intent(  out) :: tagbox(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    real(kind = dp_t), intent(in   ) :: u_vec(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,1:3)
    real(kind = dp_t), intent(in   ) :: mf_pf(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
    real(kind = dp_t), intent(in   ) :: mf_uc(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
    real(kind = dp_t), intent(in   ) :: mf_th(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
    real(dp_t)       , intent(in   ) :: dx,para_uc,para_th,amr_thers,para_u,multi_add_tau
    integer          , intent(in   ) :: lev

    ! local variables
    integer :: i,j,k
    real(kind = dp_t)  :: var_e,gra_x,gra_y,gra_z,cri_thers

    ! initially say that we do not want to tag any cells for refinement
    tagbox = .false.
    !cri_thers = amr_thers * ( 4.d0**(lev-1) )
    cri_thers = amr_thers

    ! tag all boxes where the first component of mf >= 1.01
    do k = lo(3),hi(3)
        do j = lo(2),hi(2)
            do i = lo(1),hi(1)
                gra_x = ( mf_pf(i+1,j,k) - mf_pf(i-1,j,k) ) / 2.d0
                gra_y = ( mf_pf(i,j+1,k) - mf_pf(i,j-1,k) ) / 2.d0
                gra_z = ( mf_pf(i,j,k+1) - mf_pf(i,j,k-1) ) / 2.d0
                var_e = sqrt( gra_x**2 + gra_y**2 + gra_z**2 )

                gra_x = ( mf_uc(i+1,j,k) - mf_uc(i-1,j,k) ) / 2.d0
                gra_y = ( mf_uc(i,j+1,k) - mf_uc(i,j-1,k) ) / 2.d0
                gra_z = ( mf_uc(i,j,k+1) - mf_uc(i,j,k-1) ) / 2.d0
                var_e = var_e + para_uc * sqrt( gra_x**2 + gra_y**2 + gra_z**2 )

                gra_x = ( mf_th(i+1,j,k) - mf_th(i-1,j,k) ) / 2.d0
                gra_y = ( mf_th(i,j+1,k) - mf_th(i,j-1,k) ) / 2.d0
                gra_z = ( mf_th(i,j,k+1) - mf_th(i,j,k-1) ) / 2.d0
                var_e = var_e + para_th * sqrt( gra_x**2 + gra_y**2 + gra_z**2 )

                gra_x = ( u_vec(i+1,j,k,1) - u_vec(i-1,j,k,1) ) / 2.d0
                gra_y = ( u_vec(i,j+1,k,2) - u_vec(i,j-1,k,2) ) / 2.d0
                gra_z = ( u_vec(i,j,k+1,3) - u_vec(i,j,k-1,3) ) / 2.d0
                var_e = var_e + para_u * sqrt( gra_x**2 + gra_y**2 + gra_z**2 )

                !var_e = var_e * 3.d0 * dx * dx
                !var_e = var_e * dx

                if (var_e .ge. cri_thers) then
                   tagbox(i,j,k) = .true.
                end if

            end do
        end do
    end do

  end subroutine tag_boxes_3d_flw

  !
  ! The *_init routines are for the initialization of the seed related tag
  ! First we label these nearby to the seed center, as specified by safe_n_init
  ! Then do the calculation for some time, normally 1000 steps with dt = 0.1*
  ! aiming at achieving a quasi-steady solute state, then the computation is
  ! proceeded by regriding and advances
  !
  subroutine tag_boxes_init_flw(flw,mf,seed_pos,seed_r,seed_num,tagboxes,dx,prob_lo,prob_hi,safe_n_init,&
                                seed_type,solve_flow_only,aux_tag_mf)

    type( multifab)         , intent(in   ) :: flw,mf
    integer                 , intent(in   ) :: seed_num,seed_type,solve_flow_only
    real(kind=dp_t)         , intent(in   ) :: prob_lo(mf%dim),prob_hi(mf%dim)
    real(dp_t)              , intent(in   ) :: seed_pos(seed_num,mf%dim)
    real(dp_t)              , intent(in   ) :: seed_r,safe_n_init
    type(lmultifab)         , intent(inout) :: tagboxes
    real(dp_t)              , intent(in   ) :: dx
    type(multifab), optional, intent(in   ) :: aux_tag_mf
    ! aux_tag_mf allows user to pass in additional multifabs for tagging logic

    ! local variables
    real(kind = dp_t), pointer :: mfp(:,:,:,:)
    real(kind = dp_t), pointer :: flp(:,:,:,:)
    logical          , pointer :: tp(:,:,:,:)
    integer           :: i, lo(get_dim(mf)), hi(get_dim(mf)), ng

    ng = nghost(mf)

    if (present(aux_tag_mf)) then
       call bl_error("tag_boxes.f90: aux_tag_mf passed to tag_boxes without implementation")
    end if


    do i = 1, nfabs(mf)
       mfp => dataptr(mf, i)
       tp  => dataptr(tagboxes, i)
       lo =  lwb(get_box(tagboxes, i))
       hi =  upb(get_box(tagboxes, i))

       select case (get_dim(mf))
       case (2)
          flp => dataptr(flw,i,1,2)
          call tag_boxes_2d_init_flw(tp(:,:,1,1),flp(:,:,1,:),seed_pos,seed_r,seed_num,lo,hi,&
                                     prob_lo,prob_hi,ng,dx,safe_n_init,seed_type,solve_flow_only)
       case  (3)
          flp => dataptr(flw,i,1,3)
          call tag_boxes_3d_init_flw(tp(:,:,:,1),flp(:,:,1,:),seed_pos,seed_r,seed_num,lo,hi,&
                                     prob_lo,prob_hi,ng,dx,safe_n_init,seed_type,solve_flow_only)
       end select
    end do

  end subroutine tag_boxes_init_flw

  subroutine tag_boxes_2d_init_flw(tagbox,u_vec,seed_pos,seed_r,seed_num,lo,hi,prob_lo,prob_hi,&
                                   ng,dx,safe_n_init,seed_type,solve_flow_only)

    integer          , intent(in   ) :: lo(2),hi(2),ng
    real(dp_t)       , intent(in   ) :: prob_lo(2),prob_hi(2)
    integer          , intent(in   ) :: seed_num,seed_type,solve_flow_only
    real(kind = dp_t), intent(in   ) :: u_vec(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,1:2)
    logical          , intent(  out) :: tagbox(lo(1):,lo(2):)
    real(dp_t)       , intent(in   ) :: seed_pos(seed_num,2)
    real(dp_t)       , intent(in   ) :: seed_r
    real(dp_t)       , intent(in   ) :: dx,safe_n_init

    ! local variables
    integer :: i,j,s_i
    real(kind = dp_t)  :: dis,safe_sr,x,y,l_m,l_mx,l_m2,l_mx2,x_l,y_l

    y_l = prob_hi(2) - prob_lo(2)
    x_l = prob_hi(1) - prob_lo(1)

    ! initially say that we do not want to tag any cells for refinement
    tagbox = .false.
    !safe_sr = seed_r + safe_n_init * dx
    safe_sr = safe_n_init * dx

    select case(seed_type)
    case(1)  ! columnar
       ! tag all boxes where the first component of mf >= 1.01
       do j = lo(2),hi(2)
          y = (dble(j)+0.5d0) * dx
          l_m = abs(y-seed_r)
          do i = lo(1),hi(1)
             if( l_m .le. safe_sr ) then
                tagbox(i,j) = .true.
             end if

          end do
       end do

    case(2)
      ! tag all boxes where the first component of mf >= 1.01
       do j = lo(2),hi(2)
          y =  prob_lo(2) +(dble(j)+0.5d0) * dx
          l_m  = abs(y - seed_r)
          l_m2 = abs(y - (y_l-seed_r) )
          do i = lo(1),hi(1)
            x = prob_lo(1) + (dble(i)+0.5d0) * dx
            l_mx = abs(x - seed_r)
            l_mx2 = abs(x - (x_l-seed_r) )
             if( (l_m .le. safe_sr)  .or. (l_mx .le. safe_sr) .or. &
                 (l_m2 .le. safe_sr) .or. (l_mx2 .le. safe_sr)  ) then
                tagbox(i,j) = .true.
             end if

          end do
       end do

    case default  ! Equiaxed
       ! tag all boxes where the first component of mf >= 1.01

         do j = lo(2),hi(2)
            y = prob_lo(2) + (dble(j)+0.5d0) * dx
            do i = lo(1),hi(1)
               x = prob_lo(1) + (dble(i)+0.5d0) * dx

               do s_i = 1,seed_num
                  dis = sqrt( (x - seed_pos(s_i,1))**2 + (y - seed_pos(s_i,2))**2)
                  dis = abs(dis - seed_r)

                  if( dis .le. safe_sr ) then
                     tagbox(i,j) = .true.
                  end if

               end do
             end do
         end do
       

    end select


  end subroutine tag_boxes_2d_init_flw

  subroutine tag_boxes_3d_init_flw(tagbox,u_vec,seed_pos,seed_r,seed_num,lo,hi,prob_lo,prob_hi, &
                                   ng,dx,safe_n_init,seed_type,solve_flow_only)

    integer          , intent(in   ) :: lo(3),hi(3),ng
    real(dp_t)       , intent(in   ) :: prob_lo(3),prob_hi(3)
    integer          , intent(in   ) :: seed_num,seed_type,solve_flow_only
    real(kind = dp_t), intent(in   ) :: u_vec(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,1:3)
    logical          , intent(  out) :: tagbox(lo(1):,lo(2):,lo(3):)
    real(dp_t)       , intent(in   ) :: seed_pos(seed_num,3)
    real(dp_t)       , intent(in   ) :: seed_r
    real(dp_t)       , intent(in   ) :: dx,safe_n_init

    ! local variables
    integer :: i,j,k,s_i
    real(kind = dp_t)  :: dis,safe_sr,y,x,z,l_m
    double precision :: x_left, x_right, y_top, y_bottom

    ! initially say that we do not want to tag any cells for refinement
    tagbox = .false.
    !safe_sr = seed_r + safe_n_init * dx
    safe_sr = safe_n_init * dx

    select case(seed_type)
    case(1)  ! columnar
       ! tag all boxes where the first component of mf >= 1.01
       do k = lo(3),hi(3)
           do j = lo(2),hi(2)
               y = (dble(j)+0.5d0) * dx
               l_m = abs(y - seed_r)
               do i = lo(1),hi(1)

                  if( l_m .le. safe_sr ) then
                     tagbox(i,j,k) = .true.
                   end if

               end do
           end do
       end do
 
    case (2)
       ! tag all boxes where the first component of mf >= 1.01
       do k = lo(3),hi(3)
           !z = prob_lo(3) + (dble(k)+0.5d0) * dx
           do j = lo(2),hi(2)
               y = prob_lo(2) + (dble(j)+0.5d0) * dx
               do i = lo(1),hi(1)
                  x = prob_lo(1) + (dble(i)+0.5d0) * dx

                  do s_i = 1,seed_num
                     dis = sqrt( (x - seed_pos(s_i,1))**2 + (y - seed_pos(s_i,2))**2 )
                     dis = abs(dis - seed_r)

                     if( dis .le. safe_sr ) then
                        tagbox(i,j,k) = .true.
                     end if

                  end do

               end do
           end do
       end do

    case (3)
     do k = lo(3), hi(3)
       !z = prob_lo(3) + (dble(k)+0.5d0) * dx
       do j = lo(2), hi(2)
         y = prob_lo(2) + (dble(j)+0.5d0) * dx
         do i = lo(1), hi(1)
            x = prob_lo(1) + (dble(i)+0.5d0) * dx

            do s_i = 1,seed_num

               x_left   = seed_pos(s_i,1) - seed_r - safe_sr
               x_right  = seed_pos(s_i,1) + seed_r + safe_sr
               y_bottom = seed_pos(s_i,2) - seed_r - safe_sr 
               y_top    = seed_pos(s_i,2) + seed_r + safe_sr

               if ( (x >= x_left   .and. x <= x_right) .and. &
                    (y >= y_bottom .and. y <= y_top  ) ) then
                  tagbox(i,j,k) = .true.
               end if

               x_left   = seed_pos(s_i,1) - seed_r + safe_sr
               x_right  = seed_pos(s_i,1) + seed_r - safe_sr
               y_bottom = seed_pos(s_i,2) - seed_r + safe_sr 
               y_top    = seed_pos(s_i,2) + seed_r - safe_sr

               if ( (x >= x_left   .and. x <= x_right) .and. &
                    (y >= y_bottom .and. y <= y_top  ) ) then
                  tagbox(i,j,k) = .false.
               end if

            end do

          end do  !i
        end do  !j
     end do  !k

    case default
       ! tag all boxes where the first component of mf >= 1.01
       do k = lo(3),hi(3)
           z = prob_lo(3) + (dble(k)+0.5d0) * dx
           do j = lo(2),hi(2)
               y = prob_lo(2) + (dble(j)+0.5d0) * dx
               do i = lo(1),hi(1)
                  x = prob_lo(1) + (dble(i)+0.5d0) * dx

                  do s_i = 1,seed_num
                     dis = sqrt( (x - seed_pos(s_i,1))**2 + (y - seed_pos(s_i,2))**2 + (z - seed_pos(s_i,3))**2)
                     dis = abs(dis - seed_r)

                     if( dis .le. safe_sr ) then
                        tagbox(i,j,k) = .true.
                     end if

                  end do

               end do
           end do
       end do

    end select



  end subroutine tag_boxes_3d_init_flw

end module tag_boxes_module
