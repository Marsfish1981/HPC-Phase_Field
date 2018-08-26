! print out the extrema of each variable in the plotfiles 


program fextrema

  use plotfile_module
  use filler_module
  use f2kcli
  use bl_IO_module

  implicit none

  integer f, farg
  integer :: n

  integer :: unit

  type(plotfile) :: pf
  integer narg
  character(len=256) ::  fname, varname

  integer :: ntime, numvars
  integer :: max_level

  character (len = 3) :: minname(50), maxname(50)

  logical :: single

  integer :: ivar

  single = .FALSE.

  data minname / 50*"min"/
  data maxname / 50*"max"/

  narg = command_argument_count()

 
  farg = 1

  varname = ''

  do while (farg <= narg)

     call get_command_argument(farg, value = fname)

     select case (fname)

     case ('-s','--single')
        single = .TRUE.

     case ('-v','--variable')
        farg = farg + 1
        call get_command_argument(farg, value = varname)

     case default
        exit

     end select

     farg = farg + 1
     
  enddo

  if ( farg > narg ) then
     print *, " "
     print *, "Report the extrema (min/max) for each variable in a plotfile"
     print *, "usage: "
     print *, "   fextrema [-s|--single] {[-v|--variable] name} plotfiles"
     print *, " "
     print *, "By default, the variable information is specified in columns, one line"
     print *, "per file.  If -s or --single are specified, then for each plotfile, each"
     print *, "variable's information is printed on a separate line"
     print *
     stop
  endif

  unit = unit_new()

  ! ntime is the number of files to loop over
  ntime = narg - farg  + 1

  do f = 1, ntime

     call get_command_argument(f + farg - 1, value = fname)
     call build(pf, fname, unit)

     ! if we are outputting only a single variable, make sure it exists
     ivar = -1
     if (varname /= '') ivar = plotfile_var_index(pf, trim(varname))

     max_level = pf%flevel

     numvars = pf%nvars

     if (single) then

        write (*,*) 'plotfile = ', trim(fname)
        write (*,*) 'time = ', pf%tm
        write (*,200) "variable", "minimum value", "maximum value"
        do n = 1, numvars
           write (*,201) pf%names(n), plotfile_minval(pf,n,max_level), plotfile_maxval(pf,n,max_level)
        enddo
        write (*,*) " "
     else
        if (ivar == -1) then
           if (f == 1) then
              write (*,100) "time", (pf%names(n), n =1,numvars)
              write (*,101) (minname(n), maxname(n), n=1,numvars)
           endif

           write (*,102) pf%tm, (plotfile_minval(pf,n,max_level), plotfile_maxval(pf,n,max_level), n=1,numvars)
        else
           if (f == 1) then
              write (*,100) "time", pf%names(ivar)
              write (*,101) minname(ivar), maxname(ivar)
           endif
           
           write (*,102) pf%tm, plotfile_minval(pf,ivar,max_level), plotfile_maxval(pf,ivar,max_level)
        endif

     endif

     call destroy(pf)

  enddo

100 format ("#", 1x, a22, 1x, 50("|",a44))
101 format ("#", 1x, 22x, 1x, 50("|",8x,a3,9x,1x,8x,a3,9x,3x))
102 format (1x, g22.11, 1x, 50(g20.10, 1x, g20.10, 3x))

200 format (1x, a22, 1x, a22, 1x, a22)
201 format (1x, a22, 1x, g22.11, 1x, g22.11)

end program fextrema
