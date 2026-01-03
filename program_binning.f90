program test
  use, intrinsic :: iso_fortran_env, only: int64, real64, iostat_end
  implicit none

  ! --------------------
  ! User parameters
  ! --------------------
  character(len=*), parameter :: infile  = "outputs/timeseries_L20_T2.00_MCSTOT100000000.dat"
  character(len=*), parameter :: outfile = "binning_program_L_20_T_2.0_abs_M.dat"
  integer(int64),   parameter :: header_lines = 2_int64
  integer(int64),   parameter :: burnin_valid = 1000_int64   ! ignore first 1000 VALID numeric records

  ! --------------------
  ! Data / variables
  ! --------------------
  real(real64), allocatable :: energy(:), magnet(:)
  real(real64) :: av_E, av_M, sig_E, sig_M
  integer :: i, i_max
  integer(int64) :: bb
  integer(int64) :: ndata

  integer :: uout

  ! Count how many valid numeric records exist after burn-in
  ndata = count_valid_records_after_burnin(infile, header_lines, burnin_valid)

  if (ndata <= 1_int64) then
    error stop "Not enough valid data records after burn-in."
  end if

  allocate(energy(ndata), magnet(ndata))

  call read_timeseries_after_burnin(infile, header_lines, burnin_valid, energy, magnet)

  open(newunit=uout, file=outfile, status="replace", action="write")

  ! Maximum block size exponent: bb = 2^i, with at least 2 points per block-mean estimate
  i_max = int( floor( log(real(ndata, real64)/2.0_real64) / log(2.0_real64) ) )
  if (i_max < 1) i_max = 1

  do i = 1, i_max
    bb = shiftl(1_int64, i)          ! bb = 2^i as int64
    sig_E = 0.0_real64
    sig_M = 0.0_real64

    call block_analysis(ndata, bb, energy, av_E, sig_E)
    call block_analysis(ndata, bb, magnet, av_M, sig_M)

    write(uout, *) bb, av_E, sig_E, av_M, sig_M
  end do

  close(uout)

contains

  integer(int64) function count_valid_records_after_burnin(fname, nheader, nburn) result(n_after)
    use, intrinsic :: iso_fortran_env, only: int64, real64, iostat_end
    implicit none
    character(len=*), intent(in) :: fname
    integer(int64),   intent(in) :: nheader, nburn

    integer :: u, ios
    character(len=1024) :: line
    character(len=256)  :: msg

    integer(int64) :: n_valid
    integer(int64) :: mc
    real(real64)   :: e, m

    n_after = 0_int64
    n_valid = 0_int64

    open(newunit=u, file=fname, status="old", action="read", iostat=ios, iomsg=msg)
    if (ios /= 0) error stop "OPEN failed: " // trim(msg)

    ! Skip header lines (as raw lines)
    call skip_n_lines(u, nheader)

    do
      read(u, '(A)', iostat=ios) line
      if (ios == iostat_end) exit
      if (ios /= 0) cycle

      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#') cycle

      read(line, *, iostat=ios) mc, e, m
      if (ios /= 0) cycle

      n_valid = n_valid + 1_int64
      if (n_valid > nburn) n_after = n_after + 1_int64
    end do

    close(u)
  end function count_valid_records_after_burnin


  subroutine read_timeseries_after_burnin(fname, nheader, nburn, energy, magnet)
    use, intrinsic :: iso_fortran_env, only: int64, real64, iostat_end
    implicit none
    character(len=*), intent(in) :: fname
    integer(int64),   intent(in) :: nheader, nburn
    real(real64),     intent(out) :: energy(:), magnet(:)

    integer :: u, ios
    character(len=1024) :: line
    character(len=256)  :: msg

    integer(int64) :: n_valid, idx
    integer(int64) :: mc
    real(real64)   :: e, m

    n_valid = 0_int64
    idx     = 0_int64

    open(newunit=u, file=fname, status="old", action="read", iostat=ios, iomsg=msg)
    if (ios /= 0) error stop "OPEN failed: " // trim(msg)

    call skip_n_lines(u, nheader)

    do
      read(u, '(A)', iostat=ios) line
      if (ios == iostat_end) exit
      if (ios /= 0) cycle

      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#') cycle

      read(line, *, iostat=ios) mc, e, m
      if (ios /= 0) cycle

      n_valid = n_valid + 1_int64
      if (n_valid <= nburn) cycle

      idx = idx + 1_int64
      if (idx > int(size(energy), int64)) exit   ! safety

      energy(idx) = e
      magnet(idx) = abs(m)
    end do

    close(u)

    if (idx < int(size(energy), int64)) then
      ! If fewer records were read than counted (e.g., file changed between passes),
      ! shrink arrays logically by leaving trailing values unused.
      ! (Kept simple: just stop to avoid silent mismatch.)
      error stop "Read fewer valid records than expected (file may have changed)."
    end if
  end subroutine read_timeseries_after_burnin


  subroutine skip_n_lines(u, n)
    use, intrinsic :: iso_fortran_env, only: int64, iostat_end
    implicit none
    integer,        intent(in) :: u
    integer(int64), intent(in) :: n
    integer(int64) :: k
    integer :: ios
    character(len=1024) :: line

    do k = 1_int64, n
      read(u, '(A)', iostat=ios) line
      if (ios == iostat_end) exit
    end do
  end subroutine skip_n_lines


  subroutine block_analysis(ndata, bb, data, avg, sig)
    use, intrinsic :: iso_fortran_env, only: int64, real64
    implicit none
    integer(int64), intent(in)  :: ndata, bb
    real(real64),   intent(in)  :: data(:)
    real(real64),   intent(out) :: avg, sig

    integer(int64) :: nb, ib
    integer(int64) :: beg, endd
    real(real64)   :: co
    real(real64), allocatable :: blocks(:)

    if (ndata <= 0_int64) then
      avg = 0.0_real64
      sig = 0.0_real64
      return
    end if

    nb = ndata / bb   ! integer division; leftover tail is ignored (same as your original logic)

    if (nb <= 1_int64) then
      avg = sum(data(1:int(ndata))) / real(ndata, real64)
      sig = 0.0_real64
      return
    end if

    allocate(blocks(int(nb)))  ! nb should be safe to fit in default int for your sizes
    blocks = 0.0_real64

    do ib = 1_int64, nb
      beg  = 1_int64 + bb*(ib-1_int64)
      endd = bb*ib
      blocks(int(ib)) = sum(data(int(beg):int(endd))) / real(bb, real64)
    end do

    avg = sum(blocks) / real(nb, real64)

    co = sum( (blocks - avg)**2 ) / ( real(nb, real64) * real(nb-1_int64, real64) )
    sig = sqrt(co)

    deallocate(blocks)
  end subroutine block_analysis

end program test
