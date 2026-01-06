program jackknife_Xi_C 
    ! Compilation:
    ! gfortran -g -Wall -fbounds-check -fbacktrace jackknife.f90 -o jackknife.exe/jackknife.out
    ! Execution:
    ! ./jackknife.exe/jackknife.out L T B filename (where B is some optimal block size)
    implicit none
    integer(kind=8) :: stage, idata, ndata, iostat, i, block_max, block_opt
    double precision, allocatable :: ener(:), magnet(:)
    double precision :: T, L, Xi_norm, C_norm
    double precision :: jackM2_avg, jackM2_err, jackabsM_avg, jackabsM_err, naive_err
    double precision :: jackXi_avg, jackXi_err, jackHcap_avg, jackHcap_err
    double precision, external :: log2, identity, generic_var, Xierr_prop
    character(len=50) :: dummy, inL, inT, inBopt, infile, instage

    ! Read system parameters
    call getarg(1, inL)
    call getarg(2, inT)
    call getarg(3, inBopt)
    call getarg(4, instage)
    call getarg(5, infile)
    read(inL, *) L
    read(inT, *) T
    read(inBopt, *) block_opt
    read(instage, *) stage 
    ! Proportionality constants for susceptibility and heat capacity
    Xi_norm=dble(dble(T)*dble(L)**2.d0)
    C_norm=dble((dble(L)*dble(T))**2.d0)

    ! Count number of data/lines in file
    open(1, file=infile)
        read(1, *, iostat=iostat) dummy
        read(1, *, iostat=iostat) dummy
        ndata=0
        do while (iostat.ge.0)
            read(1, *, iostat=iostat)
            ndata=ndata+1
        end do
    close(1)

    ! Save data in variables
    allocate(ener(ndata))
    allocate(magnet(ndata))

    open(1, file=infile)
        read(1, *, iostat=iostat) dummy
        read(1, *, iostat=iostat) dummy
        ndata=1
        iostat=0
        do while (iostat.ge.0)  
            read(1, *, iostat=iostat) idata, ener(ndata), magnet(ndata)
            ndata=ndata+1
        end do
    close(1)
    ! Multiply E/N and M/N by N to get E and M (for correct Jackknife)
    ener(:)=ener(:)*dble(L**2.d0)
    magnet(:)=magnet(:)*dble(L**2.d0)

    ! Determine maximum block size power from data size
    block_max = floor(log2(ndata/2.d0))

    if (stage.eq.1) then 
        ! --- Xi_err and C_err vs block size---
        ! From this stage we can extract the optimal value of the block size 
        open(3, file="jackerr_Xi-C_vs_block.out", status='unknown', position='append', action='write')
        write(3, '(A4,X,F5.1X,A2,F4.2)') "# L=", L, "T=", T
        do i=1,block_max
            ! Compute "naive" Xi error propagation for fiven block-size
            call jackknife(2**i, ndata, magnet**2.d0, magnet**2.d0, identity, jackM2_avg, jackM2_err)
            call jackknife(2**i, ndata, abs(magnet), abs(magnet), identity, jackabsM_avg, jackabsM_err)
            naive_err=Xierr_prop(Xi_norm, jackabsM_avg, jackM2_err, jackabsM_err)
            ! Compute "good" Jackknife errors for fiven block-size
            call jackknife(2**i, ndata, ener**2.d0, ener, generic_var, jackHcap_avg, jackHcap_err)
            call jackknife(2**i, ndata, magnet**2.d0, abs(magnet), generic_var, jackXi_avg, jackXi_err)
            write(3, *) 2**i, jackXi_err/Xi_norm, jackHcap_err/C_norm, naive_err
        end do
        write(3, '(2/)')
        close(3)
    else if (stage.eq.2) then
        ! ---Xi and C vs T---
        ! Use the otpimal value of the block size from the previous stage
        open(4, file="jack_Xi-C_vs_T.out", status='unknown', position='append', action='write')
        write(4, '(A4,X,F5.1,X,A3,X,I10)') "# L=", L, "BB=", block_opt
        ! -- Heat capacity per spin (Hcap) with optimal block size
        call jackknife(block_opt, ndata, ener**2.d0, ener, generic_var, jackHcap_avg, jackHcap_err)
        ! --Magnetic susceptibility per spin (Xi) with optimal block size
        call jackknife(block_opt, ndata, magnet**2.d0, abs(magnet), generic_var, jackXi_avg, jackXi_err)
        write(4, *) T, jackXi_avg/Xi_norm, jackHcap_avg/C_norm, jackXi_err/Xi_norm, jackHcap_err/C_norm
        write(4, '(2/)')
        close(4)
    end if
    
end program jackknife_Xi_C

double precision function log2(x)
    ! Computes the logarithm of x in base 2.
  implicit none
  double precision, intent(in) :: x
  log2=log(x)/log(2.d0)
end function

double precision function Xierr_prop(Xi_norm, avg_absM, err_M2, err_absM)
    ! Computes the error (std) propagation of the susceptibility Xi,
    ! from error(<M**2>) and error(<|M|>**2) from binning. 
    ! Note that this ignores possible correlations between <M**2> and <|M|>.
    implicit  none
    double precision :: Xi_norm, err_M2, err_absM, avg_absM
    Xierr_prop=sqrt(err_M2**2.d0+4.d0*(avg_absM*err_absM)**2.d0)/dble(Xi_norm)
end function

double precision function identity(x, y)
    ! Maps (X, Y) into X (the first argument)
    implicit none
    double precision :: x, y
    identity=x
end function

double precision function generic_var(x2_avg, y_avg)
    ! Generalized variance function; maps (<X^2>, <Y>) to <X^2> - <Y>^2.
    ! Ignores all constant prefactors.
    implicit none
    double precision, intent(in) :: x2_avg, y_avg
    generic_var=x2_avg-y_avg**2.d0
end function

subroutine jackknife(block_size, ndat, xdat, ydat, funxy, funjack_avg, funjack_err)
    ! Implements the Jackknife resampling method over a general function 'funxy'.
    ! Assumes that the function only depends on 2 variables x, y.
    ! Assumes the same data size 'ndat' for both X and Y.
    ! Assumes that an optimal block size (binning-wise) given 'ndat'.
    implicit none
    integer(kind=8), intent(in) :: block_size, ndat
    integer(kind=8) :: iblock, nblocks, tini, tend
    double precision, intent(in) :: xdat(ndat), ydat(ndat)
    double precision :: funxy
    double precision, intent(out) :: funjack_avg, funjack_err
    double precision, allocatable :: xblocks(:), yblocks(:)
    double precision, allocatable :: xjack(:), yjack(:), funjack(:)

    ! Initialize
    nblocks=ndat/block_size
    allocate(xblocks(nblocks))
    allocate(yblocks(nblocks))
    allocate(xjack(nblocks))
    allocate(yjack(nblocks))
    allocate(funjack(nblocks))

    ! Compute *block averages* over X, Y data
    xblocks(:)=0.d0
    yblocks(:)=0.d0
    do iblock=1,nblocks
        tini=(iblock-1)*block_size+1
        tend=iblock*block_size
        xblocks(iblock)=sum(xdat(tini:tend))
        yblocks(iblock)=sum(ydat(tini:tend))
    end do
    xblocks(:)=xblocks(:)/dble(block_size)
    yblocks(:)=yblocks(:)/dble(block_size)

    ! Compute Jackknife averages over X, Y data
    xjack(:)=sum(xblocks)
    yjack(:)=sum(yblocks)
    do iblock=1,nblocks
        xjack(iblock)=xjack(iblock)-xblocks(iblock)
        yjack(iblock)=yjack(iblock)-yblocks(iblock)
    end do
    xjack(:)=xjack(:)/dble(nblocks)
    yjack(:)=yjack(:)/dble(nblocks)

    ! Apply function to Jackknife averages (construct vector)
    funjack(:)=0.d0
    do iblock=1,nblocks
        funjack(iblock)=funxy(xjack(iblock),yjack(iblock))
    end do

    ! Compute total Jacknife average and error (std)
    funjack_avg=sum(funjack)/dble(nblocks)
    funjack_err=0.d0
    do iblock=1,nblocks
        funjack_err=funjack_err+(funjack(iblock)-funjack_avg)**2.d0
    end do
    funjack_err=sqrt(dble(nblocks-1)*funjack_err/dble(nblocks))
end subroutine jackknife