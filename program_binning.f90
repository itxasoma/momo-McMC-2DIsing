PROGRAM TEST
IMPLICIT NONE
double precision, allocatable, dimension(:):: energy, magnet
double precision:: av_E, av_M
double precision:: sig_E, sig_M
double precision :: log2
integer:: i, j, k, MCSTEP, io, idx_start,Nb,i_max
integer(kind=8):: bb, nlines
character(len=50):: relleno

! Reads time series files with the MC results and counts lines (MCsteps)
nlines = 0
open(55, file = "outputs/timeseries_L20_T2.00_MCSTOT100000000.dat")
do 
    read(55,*,IOSTAT=IO)
    if (IO/=0) EXIT
    nlines = nlines + 1
enddo
close(55)

! Ignores header (first two lines) and store energy and magnetization
open(55, file = "outputs/timeseries_L20_T2.00_MCSTOT100000000.dat")
read(55,*) relleno
read(55,*) relleno
allocate(energy(nlines-1002))
allocate(magnet(nlines-1002))

! Ignore first 1000 lines (as requested)
do i = 1003, nlines
    read(55,*) MCSTEP, energy(i-1002), magnet(i-1002)
    ! We are intesested in |M|, not M for the binning
    magnet(i-1002)=abs(magnet(i-1002))
enddo
close(55)

open(88, file = "binning_program_L_20_T_2.0_abs_M.dat")

! Compute maximum block size power from data size
i_max = floor(log2((nlines-1002)/2.d0))
do i = 1, i_max
    bb=2**i
    sig_E=0.d0
    sig_M=0.d0
    call block_analysis(nlines-1002, bb, energy, av_E, sig_E)
    call block_analysis(nlines-1002, bb, magnet, av_M, sig_M)
    
    write(88,*) bb, av_E, sig_E, av_M, sig_M
enddo

close(88)


!   37 | call block_analysis(nlines-1002, 2_int8**i, energy, av_E, sig_E)
!      |                                        1
!Error: Missing kind-parameter at (1)

END PROGRAM

double precision function log2(x)
  implicit none
  double precision, intent(in) :: x
  log2=log(x)/log(2.d0)
end function

subroutine block_analysis(ndata, BB, data, avg, sig)
!   ******* BLOCK AVERAGE METHOD ***************
!   Intended workflow: (1) run simulation, (2) then make block analyis FOR 1 OBSERVABLE

!   Total simulation steps: N
!   Size of 1 block : BB
!   Number of blocks : Nb=N/BB
!   Quantity to sample : m
!   Block average of quantity m: m_bar = average(m_i in block of size B)
!   Simulation average of quantity m: <m>=average(m_bar in all blocks)
!   Simulation error bar of quantity m: sigma(m)
!   c_0'=mean square dispalcement of each block from the simulation mean
! **********************************************
    IMPLICIT NONE
    integer:: Nb,ib, beg, end
    integer(kind=8), intent(in):: ndata, BB
    double precision, intent(in):: data(ndata)
    double precision, intent(out):: avg, sig
    double precision:: co
    double precision, allocatable:: blocks(:)

    ! Compute nº of blocks (assumes Nb=10^k and BB up to 2^k or 5^k)
    Nb = ndata / BB
    if (Nb.le.1) then
        print*, "Nb < 2 !!!"
        avg = sum(data) / dble(ndata)
        sig = 0.d0
        return
    end if
    allocate(blocks(Nb))

    ! Initialize
    blocks(:) = 0.d0
    co  = 0.d0
    sig = 0.d0
    avg = 0.d0
    
    ! Compute block averages
    do ib = 1, Nb
        ! Translate nº of block to samples segment (slicing)
        beg = 1+BB*(ib-1)
        end = BB*ib
        blocks(ib) = sum(data(beg:end))/dble(BB)
    end do
    ! Compute simulation average
    avg = sum(blocks(:)) / dble(Nb)

    ! Compute simulation error bar
    do ib=1,Nb
        co = co + (blocks(ib)-avg)**2.d0
    end do
    co = co / (dble(Nb)*dble(Nb-1))
    print*, co
    
    ! Standard error of the mean
    sig = sqrt(co)
end subroutine block_analysis
