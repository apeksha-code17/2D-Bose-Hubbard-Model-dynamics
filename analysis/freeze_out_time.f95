program find_freeze_point
  implicit none
  integer :: iounit_in, iounit_out
  integer :: num_points
  real(kind=8) :: TQ, time, Jth, rhoAvg, psi_modAvg, O_t, E_res
  real(kind=8) :: tolerance
  logical :: found

  ! Define tolerance for comparison
  tolerance = 1D-2
  found = .false.

  ! Open the input and output files
  open(unit=10, file='outp_400_1.dat', status='old', action='read')
  open(unit=20, file='plot.dat', status='unknown',position='append', action='write')

  ! Read through the file line by line
  do
     read(10,*,end=100) TQ, time, Jth, rhoAvg, psi_modAvg, O_t, E_res

     ! Check the condition O_t - 1D-4 == O_t
     if (abs(O_t - 1.0D0) > tolerance .and. .not. found) then
        write(20, '(5F15.8)') TQ, O_t,Jth, time, E_res
        found = .true.
     end if
  end do

100 continue
  if (.not. found) then
     print *, 'No matching point found with O_t close to 1.0D0 within tolerance.'
  else
     print *, 'Freeze point data written to freeze_points.dat'
  end if

  ! Close the files
  close(10)
  close(20)

end program find_freeze_point
