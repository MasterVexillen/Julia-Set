module complex_stuffs
  implicit none

  integer, parameter :: dp = selected_real_kind(15,300)
  real(kind=dp), parameter :: pi = 3.141592653589793_dp

contains
  function f_julia(z,c) !P(z)
    implicit none

    complex(kind=dp), intent(in) :: z, c
    complex(kind=dp) :: f_julia

    f_julia = z**2 + c
  end function f_julia

  function f_newraph(z,c) ! f(z) = z - P(z)/P'(z)
    implicit none

    complex(kind=dp), intent(in) :: z, c
    complex(kind=dp) :: f_newraph

    f_newraph = z - f_julia(z,c) / (-2.0_dp*z*sin(z**2))
  end function f_newraph
  
  function cmplog(z)
    implicit none

    complex(kind=dp), intent(in) :: z
    complex(kind=dp) :: cmplog
    real(kind=dp) :: arg

    arg = atan(aimag(z)/real(z,dp))
    cmplog = cmplx(log(abs(z)), arg, dp)
  end function cmplog

  function cmpsinh(z)
    implicit none

    complex(kind=dp), intent(in) :: z
    complex(kind=dp) :: cmpsinh

    cmpsinh = 0.5_dp * (exp(z) - exp(-z))
  end function cmpsinh

  function cmppow(z,pow)
    implicit none

    complex(kind=dp), intent(in) :: z, pow
    complex(kind=dp) :: cmppow, lnpow
    real(kind=dp) :: z_abs, z_arg

    z_abs = abs(z)
    z_arg = atan(aimag(z)/real(z))
    lnpow = cmplx(real(pow)*log(z_abs)-aimag(pow)*z_arg, aimag(pow)*log(z_abs)+real(pow)*z_arg, dp)
    cmppow = exp(lnpow)
  end function cmppow
  
  function cmpsq(z)
    implicit none

    complex(kind=dp), intent(in) :: z
    complex(kind=dp) :: cmpsq

    cmpsq = cmppow(z,cmplx(0.5_dp,0.0_dp,dp))
  end function cmpsq
  
end module complex_stuffs

program julia_all
  use complex_stuffs
  implicit none

  integer, parameter :: maxcol = 256

  integer, dimension(2) :: pixel

  real(kind=dp), dimension(4) :: ranges ! (xmin, xmax, ymin, ymax)
  real(kind=dp), dimension(2) :: cval

  complex(kind=dp) :: z, c
  integer :: i, j, in_unit, out_unit, grad_unit, colour_unit, nr_unit, ierr
  integer :: max_iter, max_rad, currstep, maxstep, nrstep, max_grad, min_grad, max_nr, min_nr
  real(kind=dp) :: re_ctr, im_ctr, dx, dy, hue_offset, re_range, im_range, nrlimit
  integer, dimension(:,:), allocatable :: steps, grad, abs_newraph
  real(kind=dp), dimension(3) :: hsv_curr
  real(kind=dp), dimension(3) :: hsv_freq
  integer, dimension(3) :: rgb_curr
  integer, dimension(2) :: grad_curr
  real(kind=dp), dimension(:,:), allocatable :: hsv
  integer, dimension(:,:), allocatable :: rgb

  namelist /inp/ pixel, re_ctr, im_ctr, re_range, max_iter, max_rad, hsv_freq, hue_offset
  namelist /param/ cval

  ! read in namelist
  in_unit = 15
  open(in_unit, file='nml_julia', status='old', iostat=ierr)
  if (ierr.ne.0) stop 'Error opening nml_julia in main.'

  read(15, nml=inp)
  read(15, nml=param)
  
  ! determine max values of Re(z) and Im(z)
  ranges(1) = re_ctr - 0.5_dp * re_range
  ranges(2) = re_ctr + 0.5_dp * re_range
  im_range = re_range * real(pixel(2),dp) / real(pixel(1),dp)
  ranges(3) = im_ctr - 0.5_dp * im_range
  ranges(4) = im_ctr + 0.5_dp * im_range
  write(*, 10) 'Canvas domain: (', ranges(1), ',', ranges(2), ') x (', ranges(3), ',', ranges(4), ')'
  if (cval(2).lt.0) then
     write(*, 11) 'C = ', cval(1), ' -', abs(cval(2)), 'i'
  else
     write(*, 11) 'C = ', cval(1), ' +', cval(2), 'i'
  end if
10 format(A,F8.5,A,F8.5,A,F8.5,A,F8.5,A)
11 format(A,F8.5,A,F8.5,A)

  ! write ppm file headers
  out_unit = 16
  open(out_unit, file='julia.ppm', status='replace', iostat=ierr)
  if (ierr.ne.0) stop 'Error opening julia.ppm in main'

  write(out_unit, '(A2)') 'P3'
  write(out_unit, '(2(I4.4,1X))') pixel(1), pixel(2)
  write(out_unit, '(I3.3)') maxcol

  grad_unit = 17
  open(grad_unit, file='equi.ppm', status='replace', iostat=ierr)
  if (ierr.ne.0) stop 'Error opening equi.ppm in main'

  write(grad_unit, '(A2)') 'P3'
  write(grad_unit, '(2(I4.4,1X))') pixel(1)-2, pixel(2)-2
  write(grad_unit, '(I3.3)') maxcol

  nr_unit = 18
  open(nr_unit, file='newton.ppm', status='replace', iostat=ierr)
  if (ierr.ne.0) stop 'Error opening newton.ppm in main'

  write(nr_unit, '(A2)') 'P3'
  write(nr_unit, '(2(I4.4,1X))') pixel(1), pixel(2)
  write(nr_unit, '(I3.3)') maxcol

  colour_unit = 19
  open(colour_unit, file='colour', status='replace', iostat=ierr)
  if (ierr.ne.0) stop 'Error opening colour in main'

  ! calculations
  allocate(steps(pixel(2),pixel(1)), stat=ierr)
  if (ierr.ne.0) stop 'Error allocating steps in main.'

  allocate(grad(2:pixel(2)-1,2:pixel(1)-1), stat=ierr)
  if (ierr.ne.0) stop 'Error allocating grad in main.'

  allocate(abs_newraph(pixel(2),pixel(1)), stat=ierr)
  if (ierr.ne.0) stop 'Error allocating abs_newraph in main.'

  c = cmplx(cval(1), cval(2), dp)
  dx = (ranges(2)-ranges(1)) / real(pixel(1),dp)
  dy = (ranges(4)-ranges(3)) / real(pixel(2),dp)

  !$OMP PARALLEL DO &
  !$OMP DEFAULT(NONE) &
  !$OMP PRIVATE(i,j,z,currstep) &
  !$OMP SHARED(c,pixel,ranges,dx,dy,max_iter,max_rad,steps)
   do i = 1, pixel(2)
     do j = 1, pixel(1)
        z = cmplx(ranges(1) + real(j*dx, dp), ranges(4) - real(i*dy, dp), dp)
        currstep = 0
        do while (currstep .lt. max_iter)
           z = f_julia(z,c)
           currstep = currstep + 1
           if (abs(z) .lt. 1.0e-15_dp) then
              currstep = max_iter
              exit
           end if
           if (abs(z) .gt. max_rad) exit
        end do
        steps(i,j) = currstep
     end do
  end do
  !$OMP END PARALLEL DO
  maxstep = maxval(steps)
  
  do i = 2, pixel(2)-1
     do j = 2, pixel(1)-1
        grad_curr = (/ steps(i+1,j)-steps(i-1,j), steps(i,j+1)-steps(i,j-1) /)
        grad(i,j) = nint(norm2(real(grad_curr,dp)))
     end do
  end do
  min_grad = minval(grad)
  max_grad = maxval(grad)

  ! !$OMP PARALLEL DO &
  ! !$OMP DEFAULT(NONE) &
  ! !$OMP PRIVATE(i,j,z,currstep) &
  ! !$OMP SHARED(c,pixel,ranges,dx,dy,nrstep,abs_newraph,nrlimit)
  ! do i = 1, pixel(2)
  !    do j = 1, pixel(1)
  !       z = cmplx(ranges(1) + real(j*dx, dp), ranges(4) - real(i*dy, dp), dp)
  !       do currstep = 1, nrstep
  !          z = f_newraph(z,c)
  !          if (abs(aimag(z)).gt.1.0e10_dp) then
  !             z = cmplx(real(nrlimit,dp),0.0_dp,dp)
  !             exit
  !          end if
  !       end do
  !       abs_newraph(i,j) = nint(abs(z))
  !    end do
  ! end do
  ! !$OMP END PARALLEL DO
  ! min_nr = minval(abs_newraph)
  ! max_nr = maxval(abs_newraph)

  call julia_out()
  call equipot_out()
  ! call nr_out()

  deallocate(steps, stat=ierr)
  if (ierr.ne.0) stop 'Error deallocating steps in main.'

  deallocate(grad, stat=ierr)
  if (ierr.ne.0) stop 'Error deallocating grad in main.'

  deallocate(abs_newraph, stat=ierr)
  if (ierr.ne.0) stop 'Error deallocating abs_newraph in main.'

  close(in_unit, iostat=ierr)
  if (ierr.ne.0) stop 'Error closing julia.in in main.'

  close(out_unit, iostat=ierr)
  if (ierr.ne.0) stop 'Error closing julia.ppm in main.'

  close(grad_unit, iostat=ierr)
  if (ierr.ne.0) stop 'Error closing equi.ppm in main.'

  close(colour_unit, iostat=ierr)
  if (ierr.ne.0) stop 'Error closing colour in main.'

contains
  subroutine julia_out()
    ! HSV colouring scheme for julia (+ conversion to RGB)
    allocate(hsv(maxstep, 3), stat=ierr)
    if (ierr.ne.0) stop 'Error allocating hsv in main.'

    allocate(rgb(maxstep, 3), stat=ierr)
    if (ierr.ne.0) stop 'Error allocating rgb in main.'

    hue_offset = hue_offset * pi / 180.0_dp
    do i = 1, maxstep
       hsv(i,1) = mod(2.0_dp * pi * real(i,dp) * hsv_freq(1) / real(maxstep,dp) + hue_offset, 2.0_dp * pi)
       hsv(i,2) = 0.5_dp * (sin(2.0_dp * pi * i * hsv_freq(2) / real(maxstep,dp) - hue_offset) + 1)
       hsv(i,3) = 0.5_dp * (cos(2.0_dp * pi * i * hsv_freq(3) / real(maxstep,dp) - hue_offset) + 1)
       hsv_curr = hsv(i,:)
       call hsv_to_rgb(hsv_curr, rgb_curr, maxcol)
       rgb(i,:) = rgb_curr(:)
    end do

    ! output julia
    do i = 1, pixel(2)
       do j = 1, pixel(1)
          write(out_unit, '(3(i3.3,1x))') rgb(steps(i,j),:)
       end do
    end do

    deallocate(hsv, stat=ierr)
    if (ierr.ne.0) stop 'Error deallocating hsv in main.'

    deallocate(rgb, stat=ierr)
    if (ierr.ne.0) stop 'Error deallocating rgb in main.'
  end subroutine julia_out

  subroutine equipot_out()
    ! HSV colouring scheme for grad (+ conversion to RGB)
    allocate(hsv(min_grad:max_grad, 3), stat=ierr)
    if (ierr.ne.0) stop 'Error allocating hsv2 in main.'

    allocate(rgb(min_grad:max_grad, 3), stat=ierr)
    if (ierr.ne.0) stop 'Error allocating rgb2 in main.'

    do i = min_grad, max_grad
       hsv(i,1) = mod(2.0_dp * pi * real(i-min_grad,dp) * hsv_freq(1) / real(max_grad-min_grad,dp) &
            & , 2.0_dp * pi)
       hsv(i,2) = 0.5_dp * (sin(2.0_dp * pi * real(i-min_grad,dp) * hsv_freq(2) / real(max_grad-min_grad,dp) &
            & ) + 1)
       hsv(i,3) = 0.5_dp * (cos(2.0_dp * pi * real(i-min_grad,dp) * hsv_freq(3) / real(max_grad-min_grad,dp) &
            & ) + 1)
       hsv_curr = hsv(i,:)
       call hsv_to_rgb(hsv_curr, rgb_curr, maxcol)
       rgb(i,:) = rgb_curr(:)
    end do

    ! output grad
    do i = 2, pixel(2)-1
       do j = 2, pixel(1)-1
          write(grad_unit, '(3(i3.3,1x))') rgb(grad(i,j),:)
       end do
    end do

    deallocate(hsv, stat=ierr)
    if (ierr.ne.0) stop 'Error deallocating hsv in main.'

    deallocate(rgb, stat=ierr)
    if (ierr.ne.0) stop 'Error deallocating rgb in main.'
  end subroutine equipot_out

  subroutine nr_out()
    ! HSV colouring scheme for Newton-Raphson (+ conversion to RGB)
    allocate(hsv(min_grad:max_grad, 3), stat=ierr)
    if (ierr.ne.0) stop 'Error allocating hsv2 in main.'

    allocate(rgb(min_grad:max_grad, 3), stat=ierr)
    if (ierr.ne.0) stop 'Error allocating rgb2 in main.'

    print*, min_nr, max_nr
    do i = min_nr, max_nr
       hsv(i,1) = mod(2.0_dp * pi * real(i-min_nr,dp) * hsv_freq(1) / real(max_nr-min_nr,dp) &
            & , 2.0_dp * pi)
       hsv(i,2) = 0.5_dp * (sin(2.0_dp * pi * real(i-min_nr,dp) * hsv_freq(2) / real(max_nr-min_nr,dp) &
            & ) + 1)
       hsv(i,3) = 0.5_dp * (cos(2.0_dp * pi * real(i-min_nr,dp) * hsv_freq(3) / real(max_nr-min_nr,dp) &
            & ) + 1)
       hsv_curr = hsv(i,:)
       call hsv_to_rgb(hsv_curr, rgb_curr, maxcol)
       rgb(i,:) = rgb_curr(:)
       write(colour_unit,*) i, rgb_curr(:)
    end do

    ! output Newton-Raphson
    do i = 1, pixel(2)
       do j = 1, pixel(1)
          write(nr_unit, '(3(i3.3,1x))') rgb(abs_newraph(i,j),:)
       end do
    end do

    deallocate(hsv, stat=ierr)
    if (ierr.ne.0) stop 'Error deallocating hsv in main.'

    deallocate(rgb, stat=ierr)
    if (ierr.ne.0) stop 'Error deallocating rgb in main.'
  end subroutine nr_out
  
end program julia_all

subroutine hsv_to_rgb (hsv_array, rgb_array, maxcol)
  implicit none

  integer, parameter :: dp = selected_real_kind(15,300)
  real(kind=dp), parameter :: pi = 3.141592653589793_dp

  real(kind=dp), dimension(3), intent(in) :: hsv_array
  real(kind=dp), dimension(3) :: rgb_temp
  integer, dimension(3) :: rgb_array
  integer, intent(in) :: maxcol

  real(kind=dp) :: chroma, modhue, x, match
  integer :: i

  chroma = hsv_array(3) * hsv_array(2)
  modhue = hsv_array(1) * 3.0_dp / pi
  x = chroma * (1.0_dp - abs(mod(modhue,2.0_dp)-1.0_dp))

  if (modhue.gt.0.0 .and. modhue.lt.1.0) then
     rgb_temp = (/chroma, x, 0.0_dp/)
  elseif (modhue.gt.1.0 .and. modhue.lt.2.0) then
     rgb_temp = (/x, chroma, 0.0_dp/)
  elseif (modhue.gt.2.0 .and. modhue.lt.3.0) then
     rgb_temp = (/0.0_dp, chroma, x/)
  elseif (modhue.gt.3.0 .and. modhue.lt.4.0) then
     rgb_temp = (/0.0_dp, x, chroma/)
  elseif (modhue.gt.4.0 .and. modhue.lt.5.0) then
     rgb_temp = (/x, 0.0_dp, chroma/)
  elseif (modhue.gt.5.0 .and. modhue.lt.6.0) then
     rgb_temp = (/chroma, 0.0_dp, x/)
  else
     rgb_temp = (/0.0_dp, 0.0_dp, 0.0_dp/)
  end if

  match = hsv_array(3) - chroma
  do i = 1, 3
     rgb_array(i) = nint(real(maxcol,dp)*(rgb_temp(i) + match))
  end do
end subroutine hsv_to_rgb

