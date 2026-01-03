PROGRAM JACKKNIFEe
IMPLICIT NONE
double precision,allocatable,dimension(:)::energy,magnet,abs_magnet
double precision,dimension(2)::mean,mean2,var,mean2_magnet,mean_absmagnet
character(len=50)::relleno
integer::io,MCSTEP,i
integer(kind=8)::bb,nlines

! LOOK FOR SOME ERRORS IN JACKKNIFE, CHI MUST BE REALLY HIGH FOR T LOWER THAN T_C

nlines=0

OPEN(55,file='outputs/timeseries_L100_T2.00_MCSTOT100000000.dat')
do 
    READ(55,*,IOSTAT=IO)
    IF (IO/=0) EXIT
    nlines=nlines+1
enddo
CLOSE(55)

print*,nlines-1002

OPEN(55,file='outputs/timeseries_L100_T2.00_MCSTOT100000000.dat')
read(55,*) relleno
read(55,*) relleno

allocate(energy(nlines-1002))
allocate(magnet(nlines-1002))
allocate(abs_magnet(nlines-1002))


do i=1003,nlines
read(55,*) MCSTEP,energy(i-1002),magnet(i-1002)
abs_magnet(i-1002)=abs(magnet(i-1002))
enddo
CLOSE(55)

bb=2**16

call Jackknife(nlines-1002,energy,bb,mean,mean2,var)

print*,'---------------Jackknife applied to energy----------------------'

print*, 'Mean:', mean(:)
print*, 'Mean2:', mean2(:)
print*, 'Var:', var(:)


print*, 'c', var(1)/dble(4.d0*400.d0)


print*,'---------------Jackknife applied to absolute magnetization----------------'

call Jackknife(nlines-1002,magnet,bb,mean,mean2,var)

mean2_magnet=mean2(1)

call Jackknife(nlines-1002,abs_magnet,bb,mean,mean2,var)

mean_absmagnet=mean(1)


print*, 'Mean:', mean(:)
print*, 'Mean2:', mean2(:)
print*, 'Var:', var(:)

!print*, 'Chi=', (mean2_magnet - mean_absmagnet**2.d0)/dble(2*10**4)
print*, 'Chi=', (mean2_magnet)/dble(2*10**4)

!Mean:  -1.7455891699540105        3.3002498086814551E-005
!Mean2:   3.0470817147246296        1.1522639602815172E-004
!Var:   1.6445668581330084E-007   1.9357948677239489E-008
END PROGRAM

subroutine Jackknife(n_data,data,BB,mean,mean2,var)
    implicit none
    double precision, intent(in), dimension(n_data):: data
    double precision, intent(out), dimension(2) :: mean, mean2, var
    double precision, allocatable :: blocks_data(:),block_mean(:), blocks_mean2(:)
    double precision, allocatable :: jackknife_data(:),jackknife_mean2(:),jackknife_var(:)
    integer(kind=8)::n_series,Nb,ib,beg,end
    integer(kind=8), intent(in) :: n_data, BB
    
    Nb=n_data/BB
    
    allocate(blocks_data(Nb))
    allocate(jackknife_data(Nb))
    allocate(blocks_mean2(Nb))
    allocate(jackknife_mean2(Nb))
    allocate(jackknife_var(Nb))

    blocks_data(:)=0.d0
    jackknife_data(:)=0.d0

! --- Jackknife data ---
    ! Compute block averages for the data
    do ib=1,Nb
        ! Translate nº of block to samples segment (slicing)
        beg=1+BB*(ib-1)
        end=BB*ib
        blocks_data(ib)=sum(data(beg:end))/dble(BB)
    end do

    ! Compute Jackknife for data
    jackknife_data(:)=sum(blocks_data(:))
    do ib=1,Nb
        jackknife_data(ib) = jackknife_data(ib)-blocks_data(ib)
    end do
    jackknife_data(:)=jackknife_data(:)/dble(Nb-1)

! --- Jackknife <X> ---
    mean(1)=sum(jackknife_data(:))/size(jackknife_data)
    mean(2)=0.d0
    do ib=1,Nb
    mean(2)=mean(2)+(jackknife_data(ib)-mean(1))**2
    enddo
    mean(2)=(Nb-1)*mean(2)/dble(Nb)
    mean(2)=sqrt(mean(2))
    
! --- Jackknife <X^2> ---
    ! Change of variables over data-blocks
    blocks_mean2(:)=blocks_data(:)**2
    ! Compute Jackknife array for <X^2>
    jackknife_mean2(:)=sum(blocks_mean2(:))
    do ib=1,Nb
        jackknife_mean2(ib) = jackknife_mean2(ib)-blocks_mean2(ib)
    end do
    jackknife_mean2(:)=jackknife_mean2(:)/dble(Nb-1)
    ! Compute Jackknife estimate and error for <X^2>
    mean2(1)=sum(jackknife_mean2(:))/dble(Nb)   ! Estimate
    mean2(2)=0.d0   ! Error (as the *sqrt* of the Jackknife sample variance)
    do ib=1,Nb
        mean2(2)=mean2(2)+(jackknife_mean2(ib)-mean2(1))**2.d0
    end do
    mean2(2)=(Nb-1)*mean2(2)/dble(Nb)
    mean2(2)=sqrt(mean2(2))

! --- Jackknife Var(X) ---
    do ib=1,Nb
        jackknife_var(ib)=jackknife_mean2(ib)-(jackknife_data(ib))**2
    end do
    var(1)=sum(jackknife_var(:))/dble(Nb)
    var(2)=0.d0   ! Error
    do ib=1,Nb
        var(2)=var(2)+(jackknife_var(ib)-var(1))**2.d0
    end do
    var(2)=(Nb-1)*var(2)/dble(Nb)
    var(2)=sqrt(var(2))

end subroutine Jackknife