module block
IMPLICIT NONE
contains
    SUBROUTINE average_block(X,block_size,av_X,sigma_X)
        IMPLICIT NONE 
        double precision, intent(in), dimension(:):: X            
        integer, intent(in):: BLOCK_SIZE
        double precision, intent(out):: AV_X, SIGMA_X
        integer:: n_data, n_blocks, i, start_idx
        double precision, ALLOCATABLE:: block_X(:)
                
        ! Initialize variables
        av_X = 0.d0
        n_data = size(X)

        ! Integer division of how many blocks we have
        ! n_blocks = (n_data) / block_size
        n_blocks = int(dble(n_data) / block_size)
        if (n_blocks .le. 1) then
            av_X = SUM(X( 1 : n_data)) / dble(n_data)
            sigma_X = 0.d0
            return
        end if

        allocate(block_X(n_blocks))
        ! Compute each block average (for every data-block)
        do i = 1, n_blocks
            start_idx = 1 + (i-1)*block_size
            block_x(i) = SUM(x(start_idx : start_idx + block_size - 1)) / block_size
        enddo

        ! Average and error calculation            
        av_x = SUM(block_x) / n_blocks
        sigma_x = SQRT( MAX(0.d0, SUM( (block_x - av_x)**2 ) / dble(n_blocks*(n_blocks-1))) )
                
        deallocate(block_x)
    END SUBROUTINE

END MODULE           
                