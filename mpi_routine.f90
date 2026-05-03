Module mpi_routine
    use math
    use base
    use mpi
    use mod_io
    Implicit None

    contains
        subroutine singal_to_master(mpi_infos, phi)
            Implicit None
            type(mmp), target :: mpi_infos
            Real, dimension(0:nx+1, 0:ny+1) :: phi

            Integer :: mpi_i, mpi_n
            Integer, dimension(:), pointer :: domain
            Real, dimension(:, :), pointer :: trans
            Integer :: ierr, status(mpi_status_size)
            Integer :: iw, ie, js, jn
            Integer :: p_nx, p_ny
            Integer :: i

            mpi_i = mpi_infos%mpi_i
            mpi_n = mpi_infos%mpi_n
            iw = mpi_infos%t_domain(mpi_i+1, 1)
            ie = mpi_infos%t_domain(mpi_i+1, 2)
            js = mpi_infos%t_domain(mpi_i+1, 3)
            jn = mpi_infos%t_domain(mpi_i+1, 4)
            
            p_nx = ie - iw + 1; p_ny = jn - js + 1 
            do i = 1, mpi_n-1
                if(i == mpi_i) then
                    call mpi_send(phi(iw:ie, js:jn), p_nx*p_ny, mpi_real, 0, i, mpi_comm_world, ierr)
                endif
            enddo
            if(mpi_i == 0) then
                do i = 1, mpi_n-1
                    iw = mpi_infos%t_domain(i+1, 1)
                    ie = mpi_infos%t_domain(i+1, 2)
                    js = mpi_infos%t_domain(i+1, 3)
                    jn = mpi_infos%t_domain(i+1, 4)
                    p_nx = ie - iw + 1; p_ny = jn - js + 1
                    allocate(trans(iw:ie, js:jn))
                    call mpi_recv(trans, p_nx*p_ny, mpi_real, i, i, mpi_comm_world, status, ierr)
                    phi(iw:ie, js:jn) = trans
                    deallocate(trans)
                enddo
            endif
        
        end subroutine singal_to_master

        subroutine exchange(mpi_infos, phi)
            Implicit None
            type(mmp), target :: mpi_infos
            Real, dimension(0:nx+1, 0:ny+1) :: phi
            
            Real, dimension(:), pointer :: trans_lat, trans_lon
            Integer, dimension(:), pointer ::mpi_i_lon, mpi_i_lat
            Integer :: p_nx, p_ny
            Integer :: iw, ie, js, jn
            Integer :: i
            Integer :: mpi_i, mpi_n, status(mpi_status_size), ierr
            Logical :: left, right, lower, upper

            mpi_i = mpi_infos%mpi_i
            mpi_n = mpi_infos%mpi_n
            iw = mpi_infos%t_domain(mpi_i+1, 1)
            ie = mpi_infos%t_domain(mpi_i+1, 2)
            js = mpi_infos%t_domain(mpi_i+1, 3)
            jn = mpi_infos%t_domain(mpi_i+1, 4)
            mpi_i_lon => mpi_infos%mpi_loc(:, 1)
            mpi_i_lat => mpi_infos%mpi_loc(:, 2)
            p_nx = ie - iw + 1 
            p_ny = jn - js + 1
            
            left = .true.
            right = .true.
            lower = .true.
            upper = .true.
            if(mpi_i_lon(mpi_i+1) == 0) left = .false.
            if(mpi_i_lon(mpi_i+1) == np_lon-1) right = .false.
            if(mpi_i_lat(mpi_i+1) == 0) lower = .false.
            if(mpi_i_lat(mpi_i+1) == np_lat-1) upper = .false.
            
            if(left .and. lower)  then 
                call mpi_send(phi(iw+1, js+1), 1, mpi_real, mpi_i-np_lon-1, 5*mpi_n+mpi_i*4+1, mpi_comm_world, ierr)
            endif
            if(left .and. upper)  then
                call mpi_send(phi(iw+1, jn-1), 1, mpi_real, mpi_i+np_lon-1, 5*mpi_n+mpi_i*4+2, mpi_comm_world, ierr)
            endif
            if(right .and. upper) then
                call mpi_send(phi(ie-1, jn-1), 1, mpi_real, mpi_i+np_lon+1, 5*mpi_n+mpi_i*4+3, mpi_comm_world, ierr)
            endif
            if(right .and. lower) then
                call mpi_send(phi(ie-1, js+1), 1, mpi_real, mpi_i-np_lon+1, 5*mpi_n+mpi_i*4+4, mpi_comm_world, ierr)
            endif

            if(left .and. lower)  then
                call mpi_recv(phi(iw-1, js-1), 1, mpi_real, mpi_i-np_lon-1, 5*mpi_n+(mpi_i-np_lon)*4-1, mpi_comm_world, status, ierr)
            endif
            if(left .and. upper)  then
                call mpi_recv(phi(iw-1, jn+1), 1, mpi_real, mpi_i+np_lon-1, 5*mpi_n+(mpi_i+np_lon)*4  , mpi_comm_world, status, ierr)
            endif
            if(right .and. upper) then
                call mpi_recv(phi(ie+1, jn+1), 1, mpi_real, mpi_i+np_lon+1, 5*mpi_n+(mpi_i+np_lon)*4+5, mpi_comm_world, status, ierr)
            endif
            if(right .and. lower) then
                call mpi_recv(phi(ie+1, js-1), 1, mpi_real, mpi_i-np_lon+1, 5*mpi_n+(mpi_i-np_lon)*4+6, mpi_comm_world, status, ierr)
            endif
            do i = 0, mpi_n-1
                if(mpi_i == i) then
                    if(mpi_i_lon(i+1) == 0) then
                        call mpi_send(reshape(phi(ie-1, js:jn), (/p_ny/)), p_ny, mpi_real, i+1, mpi_n+i, mpi_comm_world, ierr)
                    elseif(mpi_i_lon(i+1) == np_lon-1) then
                        call mpi_send(reshape(phi(iw+1, js:jn), (/p_ny/)), p_ny, mpi_real, i-1, 2*mpi_n+i, mpi_comm_world, ierr)
                    else                        
                        call mpi_send(reshape(phi(ie-1, js:jn), (/p_ny/)), p_ny, mpi_real, i+1, mpi_n+i, mpi_comm_world, ierr)
                        call mpi_send(reshape(phi(iw+1, js:jn), (/p_ny/)), p_ny, mpi_real, i-1, 2*mpi_n+i, mpi_comm_world, ierr)
                    endif

                    if(mpi_i_lat(i+1) == 0) then
                        call mpi_send(reshape(phi(iw:ie, jn-1), (/p_nx/)), p_nx, mpi_real, i+np_lon, 3*mpi_n+i, mpi_comm_world, ierr)
                    elseif(mpi_i_lat(i+1) == np_lat-1) then
                        call mpi_send(reshape(phi(iw:ie, js+1), (/p_nx/)), p_nx, mpi_real, i-np_lon, 4*mpi_n+i, mpi_comm_world, ierr)
                    else
                        call mpi_send(reshape(phi(iw:ie, jn-1), (/p_nx/)), p_nx, mpi_real, i+np_lon, 3*mpi_n+i, mpi_comm_world, ierr)
                        call mpi_send(reshape(phi(iw:ie, js+1), (/p_nx/)), p_nx, mpi_real, i-np_lon, 4*mpi_n+i, mpi_comm_world, ierr)
                    endif
                endif
            enddo
            do i = 0, mpi_n-1
                if(i == mpi_i) then
                    allocate(trans_lat(p_nx), trans_lon(p_ny))
                    if(mpi_i_lon(i+1) == 0) then
                        call mpi_recv(trans_lon, p_ny, mpi_real, i+1, 2*mpi_n+i+1, mpi_comm_world, status, ierr)
                        phi(ie+1, js:jn) = trans_lon
                    elseif(mpi_i_lon(i+1) == np_lon-1) then
                        call mpi_recv(trans_lon, p_ny, mpi_real, i-1, mpi_n+i-1, mpi_comm_world, status, ierr)
                        phi(iw-1, js:jn) = trans_lon
                    else
                        call mpi_recv(trans_lon, p_ny, mpi_real, i+1, 2*mpi_n+i+1, mpi_comm_world, status, ierr)
                        phi(ie+1, js:jn) = trans_lon
                        call mpi_recv(trans_lon, p_ny, mpi_real, i-1, mpi_n+i-1, mpi_comm_world, status, ierr)
                        phi(iw-1, js:jn) = trans_lon
                    endif

                    if(mpi_i_lat(i+1) == 0) then
                        call mpi_recv(trans_lat, p_nx, mpi_real, i+np_lon, 4*mpi_n+i+np_lon, mpi_comm_world, status, ierr)
                        phi(iw:ie, jn+1) = trans_lat
                    elseif(mpi_i_lat(i+1) == np_lat-1) then
                        call mpi_recv(trans_lat, p_nx, mpi_real, i-np_lon, 3*mpi_n+i-np_lon, mpi_comm_world, status, ierr)
                        phi(iw:ie, js-1) = trans_lat
                    else
                        call mpi_recv(trans_lat, p_nx, mpi_real, i+np_lon, 4*mpi_n+i+np_lon, mpi_comm_world, status, ierr)
                        phi(iw:ie, jn+1) = trans_lat
                        call mpi_recv(trans_lat, p_nx, mpi_real, i-np_lon, 3*mpi_n+i-np_lon, mpi_comm_world, status, ierr)
                        phi(iw:ie, js-1) = trans_lat
                    endif
                    deallocate(trans_lon, trans_lat)
                endif
            enddo

            call mpi_barrier(mpi_comm_world, ierr)
        end subroutine exchange

End Module mpi_routine
