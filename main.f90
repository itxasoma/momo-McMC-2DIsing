PROGRAM practice3
  use r1279block
  use GLOBAL
  use lattice
  use ising_energy
  use block
  IMPLICIT NONE
  
  ! ============ VARIABLE DECLARATIONS ============
  integer:: i, iseed, ispin, delta_E, n2DE
  integer:: hh
  integer:: E, M, try_count  
  integer(kind=8):: MCS, MC_count
  integer(kind=8):: nMCS_local, nmeas_local, MCS_discard
  integer(kind=8):: n_measurements, n_after_discard, idx_start
  
  double precision:: avg_E, avg_M, avg_E_data, sigma_E_data
  double precision, allocatable:: table(:), E_data(:), M_data(:), absM_data(:)
  double precision:: time_start, time_end, flips_per_sec
  
  ! File names
  character(len=100):: filename_ts, filename_binning_E, filename_binning_M, inputfile
  character(len=20):: T_str, L_str
  character(len=50):: NMC_str
  

  write(*,*) "============================================="
  write(*,*) "  2D ISING MODEL - MONTE CARLO SIMULATION    "
  write(*,*) "============================================="
  write(*,*)
  
  !----------------------------------------------------------------------------
  ! INPUT PARAMETERS
  !----------------------------------------------------------------------------
  
  ! Manually:
  !L = 20                       ! Lattice size (20 for test, 100 for production)
  !T = 2.d0                     ! Temperature (2.0, 2.27, 2.6)
  !nMCS_local = 10**8           ! Total number of MCS
  !nmeas_local = 10             ! Measurements every nmeas MCS
  !MCS_discard = 10**3          ! Discard first MCS for equilibration
  !iseed = 12345                ! RNG seed
  !z = 4                        ! Coordination number (2D square lattice)
  
  ! Reading input file(s):
  inputfile = "DATA.in"
  call get_command_argument(1, inputfile)
  call readdata(inputfile, L, T, nMCS_local, nmeas_local, MCS_discard, iseed, z)

  write(*,*) "Input Parameters:"
  write(*,*) "  L =", L
  write(*,*) "  T =", T
  write(*,*) "  nMCS =", nMCS_local
  write(*,*) "  nmeas =", nmeas_local
  write(*,*) "  MCS_discard =", MCS_discard
  write(*,*) "  seed =", iseed
  write(*,*)
  
  !----------------------------------------------------------------------------
  ! INITIALIZATION
  !----------------------------------------------------------------------------
  
  ! Initialize RNG
  call setr1279(iseed)
  
  ! Initialize globals (sets N and beta)
  call init_globals()
  
  write(*,*) "System:"
  write(*,*) "  N = L^2 =", N
  write(*,*) "  beta = 1/T =", beta
  write(*,*)
  
  ! Allocate lattice
  call alloc_lattice(N)
  call pbc()
  call square_lattice()
  
  ! Create exponential lookup table
  allocate(table(2*z+1))
  call table_exp(z, table)
  
  ! Allocate data storage for time series
  n_measurements = nMCS_local / nmeas_local
  allocate(E_data(n_measurements))
  allocate(M_data(n_measurements))
  
  ! Random spin initialization
  write(*,*) "Initializing random spins..."
  s(:) = 0
  do i = 1, N
    s(i) = 2*mod(int(2*r1279()),2) - 1
  end do
  
  ! Calculate initial energy and magnetization
  call E_ising(s, nbr, E)
  call M_ising(s, M)
  
  write(*,*) "  Initial E/N =", dble(E)/dble(N)
  write(*,*) "  Initial M/N =", dble(M)/dble(N)
  write(*,*)
  
  !----------------------------------------------------------------------------
  ! PREPARE OUTPUT FILES
  !----------------------------------------------------------------------------
  
  write(T_str, '(F5.2)') T
  write(L_str, '(I3)') L
  write(NMC_str, '(I10)') nMCS_local
  
  filename_ts = 'timeseries_L'//trim(adjustl(L_str))//'_T'//trim(adjustl(T_str))//'_MCSTOT'//trim(adjustl(NMC_str))//'.dat'
  filename_binning_E = 'binning_E_L'//trim(adjustl(L_str))//'_T'//trim(adjustl(T_str))//'_MCTOT'//trim(adjustl(NMC_str))//'.dat'
  filename_binning_M = 'binning_M_L'//trim(adjustl(L_str))//'_T'//trim(adjustl(T_str))//'_MCTOT'//trim(adjustl(NMC_str))//'.dat'
  
  open(71, file=filename_ts)
  write(71,'(A,I4,A,F6.3,A,I12,A,I6)') &
    "# L=", L, " T=", T, " nMCS=", nMCS_local, " nmeas=", nmeas_local
  write(71,*) "# MCS    E/N    M/N"
  
  !----------------------------------------------------------------------------
  ! MONTE CARLO SIMULATION
  !----------------------------------------------------------------------------
  
  write(*,*) "Starting Monte Carlo simulation..."
  write(*,*)
  
  ! Initialize
  avg_E = 0.d0
  avg_M = 0.d0
  MC_count = 0
  
  ! Start timing
  call cpu_time(time_start)
  
  ! Main MC loop
  do MCS = 1, nMCS_local
    
    ! One Monte Carlo Step = N spin flip attempts
    do try_count = 1, N
      
      ! Select random spin
      ispin = mod(int(N*r1279()),N) + 1
      
      ! Calculate local field
      call h_local(z, s, nbr, hh, ispin)
      
      ! Calculate energy change
      delta_E = 2*s(ispin)*hh
      
      ! Index for lookup table
      n2DE = int(dble(delta_E)/2.d0) + z + 1
      
      ! Metropolis algorithm
      if (delta_E .lt. 0) then
        ! Accept: energy decreases
        s(ispin) = -s(ispin)
        E = E + delta_E
        M = M + 2*s(ispin)
        
      else if (r1279() .lt. table(n2DE)) then
        ! Accept with probability exp(-beta*delta_E)
        s(ispin) = -s(ispin)
        E = E + delta_E
        M = M + 2*s(ispin)
      endif
      
    enddo  ! End N attempts = 1 MCS
    
    !--------------------------------------------------------------------------
    ! MEASUREMENTS
    !--------------------------------------------------------------------------
    
    if (mod(MCS, nmeas_local) .eq. 0) then
      
      MC_count = MC_count + 1
      
      ! Store time series (always store, discard later in analysis)
      E_data(MC_count) = E
      M_data(MC_count) = M  ! Store M (not |M|!)
      
      ! Write time series to file
      ! For first 10^3 MCS: write all
      ! For first 10^6 MCS: write every 100th point
      ! if (MCS .le. 10**3 .or. (MCS .le. 10**6 .and. mod(MCS, 100*nmeas_local) .eq. 0)) then
      !   write(71,*) MCS, dble(E)/dble(N), dble(M)/dble(N)
      ! endif 

      ! Plot time series for E and M for the first 10**3 MCS and the first 10**6 MCS.
      if ((MCS.le.nMCS_local).and.(mod(MCS, nmeas_local).eq.0)) write(71,*) MCS, dble(E)/dble(N), dble(M)/dble(N)
      
      ! Accumulate averages (only after discarding initial MCS)
      if (MCS > MCS_discard) then
        avg_E = avg_E + dble(E)
        avg_M = avg_M + dble(abs(M))  ! Use |M| for average
      endif
      
    endif
    
    ! Progress report
    if (mod(MCS, nMCS_local/10) .eq. 0) then
      write(*,'(A,I12,A,F6.2,A)') "  MCS: ", MCS, "  (", &
        100.d0*dble(MCS)/dble(nMCS_local), "% complete)"
    endif
    
  enddo  ! End MCS loop
  
  ! End timing
  call cpu_time(time_end)
  
  close(71)
  
  write(*,*)
  write(*,*) "Simulation complete!"
  write(*,*)
  
  !----------------------------------------------------------------------------
  ! COMPUTE FINAL AVERAGES
  !----------------------------------------------------------------------------
  
  n_after_discard = (nMCS_local - MCS_discard) / nmeas_local
  
  avg_E = avg_E / dble(n_after_discard)
  avg_M = avg_M / dble(n_after_discard)
  
  ! Compute speed (attempted flips per second)
  flips_per_sec = dble(nMCS_local * N) / (time_end - time_start)
  
  write(*,*) "========================================="
  write(*,*) "             FINAL RESULTS               "
  write(*,*) "========================================="
  write(*,'(A,E12.5,A)') "  Speed: ", flips_per_sec, " flips/sec"
  write(*,*)
  write(*,'(A,F14.10)') "  <E>/N   = ", avg_E/dble(N)
  write(*,'(A,F14.10)') "  <|M|>/N = ", avg_M/dble(N)
  write(*,*) "========================================="
  write(*,*)
  
  !----------------------------------------------------------------------------
  ! BINNING ANALYSIS
  !----------------------------------------------------------------------------
  
  write(*,*) "Performing binning analysis..."
  
  ! Calculate index where to start (after discarding)
  idx_start = MCS_discard / nmeas_local + 1
  
  ! Binning for Energy
  open(80, file=filename_binning_E)
  write(80,*) "# m   <E>/N   sigma_E/N"
  
  do i = 0, 25
    if (dble(n_after_discard) / dble(2**i) .lt. 2.d0) exit
    call average_block(E_data(idx_start:MC_count), 2**i, avg_E_data, sigma_E_data)
    write(80,*) 2**i, avg_E_data/dble(N), sigma_E_data/dble(N)
  enddo
  close(80)
  
  ! Binning for |Magnetization|
  allocate(absM_data(n_after_discard))
  absM_data = abs(M_data(idx_start:MC_count))
  
  open(81, file=filename_binning_M)
  write(81,*) "# m   <|M|>/N   sigma_|M|/N"
  
  do i = 0, 25
    if (dble(n_after_discard) / dble(2**i) .lt. 2.d0) exit
    call average_block(absM_data, 2**i, avg_E_data, sigma_E_data)
    write(81,*) 2**i, avg_E_data/dble(N), sigma_E_data/dble(N)
  enddo
  close(81)
  
  deallocate(absM_data)
  
  write(*,*) "Binning analysis complete!"
  write(*,*)
  
  !----------------------------------------------------------------------------
  ! OUTPUT SUMMARY
  !----------------------------------------------------------------------------
  
  write(*,*) "Output files created:"
  write(*,*) "  ", trim(filename_ts)
  write(*,*) "  ", trim(filename_binning_E)
  write(*,*) "  ", trim(filename_binning_M)
  write(*,*)
  
  !----------------------------------------------------------------------------
  ! CLEANUP
  !----------------------------------------------------------------------------
  
  deallocate(E_data)
  deallocate(M_data)
  deallocate(table)
  call dealloc_lattice()
  
  write(*,*) "============================================="
  write(*,*) "       PROGRAM COMPLETED SUCCESSFULLY       "
  write(*,*) "============================================="
  
END PROGRAM practice3


! Allows the main program to read the input parameters from an external *.in file.
SUBROUTINE readdata(input, L, T, nMCS_local, nmeas_local, MCS_discard, iseed, z)
  IMPLICIT NONE
  character(*), intent(in):: input
  integer, intent(out):: L, iseed, z
  integer(kind=8), intent(out) :: nMCS_local, nmeas_local, MCS_discard
  double precision, intent(out):: T

  open(10, file=trim(input), status='old')
    read(10,*) L
    read(10,*) T
    read(10,*) nMCS_local
    read(10,*) nmeas_local
    read(10,*) MCS_discard
    read(10,*) iseed
    read(10,*) z
  close(10)

END SUBROUTINE readdata