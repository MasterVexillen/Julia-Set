function f(z,c)
  implicit none

  integer, parameter :: dp = selected_real_kind(15,300)
  complex(kind=dp), intent(in) :: z, c
  complex(kind=dp) :: f
  complex(kind=dp), external :: cmpsinh, cmplog

  f = z**2 + c
end function f

program julia
!  use hsv_to_rgb
  implicit none

  integer, parameter :: dp = selected_real_kind(15,300)
  integer, parameter :: maxcol = 256
  real(kind=dp), parameter :: pi = 3.141592653589793_dp

  integer, dimension(2) :: pixel

  real(kind=dp), dimension(4) :: ranges ! (xmin, xmax, ymin, ymax)

  complex(kind=dp) :: z, c
  complex(kind=dp), external :: f
  integer :: i, j, maxiter, maxrad, currstep, maxstep, in_unit, out_unit, colour_unit, ierr
  real(kind=dp) :: dx, dy, creal, cimag, hue_offset, re_range, im_range, hue
  real(kind=dp), dimension(:,:), allocatable :: hsv
  real(kind=dp), dimension(3) :: hsv_curr
  real(kind=dp), dimension(3) :: hsv_freq
  integer, dimension(3) :: rgb_curr
  integer, dimension(:,:), allocatable :: rgb

  ! read in parameters
  in_unit = 15
  open(in_unit, file='julia.in', status='old', iostat=ierr)
  if (ierr.ne.0) stop 'Error opening julia.in in main.'

  read(15, *)
  read(15, *) pixel(1), pixel(2)
  read(15, *)
  read(15, *)
  read(15, *) ranges(1), ranges(3), re_range
  read(15, *)
  read(15, *)
  read(15, *) maxiter, maxrad
  read(15, *)
  read(15, *)
  read(15, *) hsv_freq(1), hsv_freq(2), hsv_freq(3), hue_offset
  read(15, *)
  read(15, *)
  read(15, *) creal, cimag

  ! determine max values of Re(z) and Im(z)
  ranges(2) = ranges(1) + 0.5_dp * re_range
  ranges(1) = ranges(2) - re_range
  im_range = re_range * real(pixel(2),dp) / real(pixel(1),dp)
  ranges(4) = ranges(3) + 0.5_dp * im_range
  ranges(3) = ranges(4) - im_range
  write(*, 10) 'Canvas domain: (', ranges(1), ',', ranges(2), ') x (', ranges(3), ',', ranges(4), ')'
  10 format(A,F8.5,A,F8.5,A,F8.5,A,F8.5,A)
  
  ! write out ppm file
  out_unit = 16
  open(out_unit, file='julia.ppm', status='replace', iostat=ierr)
  if (ierr.ne.0) stop 'Error opening julia.ppm in main'

  write(out_unit, '(A2)') 'P3'
  write(out_unit, '(2(I4.4,1X))') pixel(1), pixel(2)
  write(out_unit, '(I3.3)') maxcol

  ! pre-run to determine max number of steps
  c = cmplx(creal, cimag, dp)
  dx = (ranges(2)-ranges(1)) / real(pixel(1),dp)
  dy = (ranges(4)-ranges(3)) / real(pixel(2),dp)
  
  maxstep = 0
  do i = 1, pixel(2)
     do j = 1, pixel(1)
        z = cmplx(ranges(1) + real(j*dx, dp), ranges(4) - real(i*dy, dp), dp)
        currstep = 0
        do while (currstep .lt. maxiter)
           z = f(z,c)
           currstep = currstep + 1
           if (abs(z) .lt. 1.0e-15_dp) then
              currstep = maxiter
              exit
           end if
           if (abs(z) .gt. maxrad) exit
        end do
        if (currstep.gt.maxstep) maxstep = currstep
     end do
  end do

  ! HSV colouring scheme (+ conversion to RGB)
  allocate(hsv(maxiter, 3), stat=ierr)
  if (ierr.ne.0) stop 'Error allocating hsv in main.'

  allocate(rgb(maxiter, 3), stat=ierr)
  if (ierr.ne.0) stop 'Error allocating rgb in main.'

  colour_unit = 17
  open(colour_unit, file='colour', status='replace', iostat=ierr)
  if (ierr.ne.0) stop 'Error opening colour in main'

  hue_offset = hue_offset * pi / 180.0_dp
  do i = 1, maxstep
     hsv(i,1) = mod(2.0_dp * pi * real(i,dp) * hsv_freq(1) / real(maxstep,dp) + hue_offset, 2.0_dp * pi)
     hsv(i,2) = 0.5_dp * (sin(2.0_dp * pi * i * hsv_freq(2) / real(maxstep,dp) - hue_offset) + 1)
     hsv(i,3) = 0.5_dp * (cos(2.0_dp * pi * i * hsv_freq(3) / real(maxstep,dp) - hue_offset) + 1)
     hsv_curr = hsv(i,:)
     call hsv_to_rgb(hsv_curr, rgb_curr, maxcol)
     rgb(i,:) = rgb_curr(:)
     write(colour_unit, *) i, hsv_curr(:), rgb_curr(:)
  end do

  ! real calculations
  c = cmplx(creal, cimag, dp)
  dx = (ranges(2)-ranges(1)) / real(pixel(1),dp)
  dy = (ranges(4)-ranges(3)) / real(pixel(2),dp)
  
  do i = 1, pixel(2)
     do j = 1, pixel(1)
        z = cmplx(ranges(1) + real(j*dx, dp), ranges(4) - real(i*dy, dp), dp)
        currstep = 0
        do while (currstep .lt. maxiter)
           z = f(z,c)
           currstep = currstep + 1
           if (abs(z) .lt. 1.0e-15_dp) then
              currstep = maxiter
              exit
           end if
           if (abs(z) .gt. maxrad) exit
        end do
        ! hsv_curr = hsv(currstep,:)
        ! call hsv_to_rgb(hsv_curr, rgb_curr, maxcol)
        write(out_unit, '(3(i3.3,1x))') rgb(currstep,:)
     end do
  end do

  deallocate(hsv, stat=ierr)
  if (ierr.ne.0) stop 'Error deallocating hsv in main.'

  deallocate(rgb, stat=ierr)
  if (ierr.ne.0) stop 'Error deallocating rgb in main.'
  
  close(in_unit, iostat=ierr)
  if (ierr.ne.0) stop 'Error closing julia.in in main.'

  close(out_unit, iostat=ierr)
  if (ierr.ne.0) stop 'Error closing julia.ppm in main.'

  close(colour_unit, iostat=ierr)
  if (ierr.ne.0) stop 'Error closing colour in main.'
  
end program julia

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

function cmplog(z)
  implicit none

  integer, parameter :: dp = selected_real_kind(15,300)

  complex(kind=dp), intent(in) :: z
  complex(kind=dp) :: cmplog
  real(kind=dp) :: arg

  arg = atan2(aimag(z), real(z,dp))
  cmplog = cmplx(log(abs(z)), arg, dp)
end function cmplog

function cmpsinh(z)
  implicit none

  integer, parameter :: dp = selected_real_kind(15,300)
  
  complex(kind=dp), intent(in) :: z
  complex(kind=dp) :: cmpsinh

  cmpsinh = 0.5_dp * (exp(z) - exp(-z))
end function cmpsinh
