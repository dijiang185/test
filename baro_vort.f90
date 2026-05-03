program main
    use base
    use mod_io
    use mpi
    use solver
    use mpi_routine
    Implicit None
    type(vars), target :: var
    type(mmp), target :: mpi_infos
    Real, dimension(:, :), allocatable :: f
    Integer :: ierr, mpi_i, mpi_n, status(MPI_STATUS_SIZE)
    Real, dimension(:, :), pointer :: phi, vort, u, v, trans, vort0
    Integer :: p_nx, p_ny, mpi_i_lat, mpi_i_lon
    Integer, pointer :: ie, iw, jn, js
    Integer, dimension(:), pointer :: domain
    Integer :: i, j
    Real :: t
    Real, dimension(:, :), pointer :: ten
    double precision :: time_start, time_end

    
    call param_init
    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_world, mpi_i, ierr)
    call mpi_comm_size(mpi_comm_world, mpi_n, ierr)
    if(mpi_i == 0) then
        print*, '********************************************************'
        print*, '*******barotropic vorticity model based on MPI**********'
        print*, '********************************************************'
        print*, 'nx', nx
        print*, 'ny', ny
        print*, 'dx', dx
        time_start = mpi_wtime()
        if(file_type == 0) call grd_init
        if(file_type == 1) call nc_init
        if(np_lon*np_lat /= mpi_n) then
            print*, 'np should be equal to np_lat*np_lon...', mpi_n, np_lat*np_lon
            stop
        endif
    endif
    
    call data_init(mpi_i, mpi_n, var, f, mpi_infos)

    phi => var%phi
    vort => var%vort
    vort0 => var%vort0
    u => var%u
    v => var%v
    if(test_case == 0) call read_nc(vort(1:nx, 1:ny))
    call boundary(mpi_infos, vort)
    vort0 = vort
    call relaxtion(mpi_infos, var)
    call singal_to_master(mpi_infos, phi)
    call phi_to_uv(mpi_infos, var)
    call singal_to_master(mpi_infos, u)
    call singal_to_master(mpi_infos, v)
    if(mpi_i == 0) then
        call write_to_nc(vort(1:nx, ny:1:-1), vortid)
        call write_to_nc(u(1:nx, ny:1:-1), uid)
        call write_to_nc(v(1:nx, ny:1:-1), vid)
        call write_to_nc(phi(1:nx, ny:1:-1), phiid)
    endif
    !Time loop
    t = 0
    call solver1(mpi_infos, var, f)
    t = t + deltat 
    Do while(t <= timax)
        if(mpi_i == 0) write(*, 100) 'integral:', t, '   wtime:', mpi_wtime()-time_start, '  seconds'
        call solver2(mpi_infos, var, f)
        !if(mod(int(t), 8*int(deltat)) == 0) then
            !call smooth(mpi_infos, vort, 0.5)
            !call smooth(mpi_infos, vort, -0.5)
        !endif
        if(mod(int(t), int(write_nt)) == 0) then
            call singal_to_master(mpi_infos, phi)
            call singal_to_master(mpi_infos, vort)
            call phi_to_uv(mpi_infos, var)
            call singal_to_master(mpi_infos, u)
            call singal_to_master(mpi_infos, v)
            if(mpi_i == 0) then
                nt = nt + 1
                print*, 'write to file'
                call write_to_nc(phi(1:nx, ny:1:-1), phiid)
                call write_to_nc(vort(1:nx, ny:1:-1), vortid)
                call write_to_nc(u(1:nx, ny:1:-1), uid)
                call write_to_nc(v(1:nx, ny:1:-1), vid)
            endif
        endif
        t = t + deltat
    Enddo

    if(mpi_i == 0) then
        call close_nc
        time_end = mpi_wtime()
        print*, 'Normal termination'
        print 101, 'Total run time:', time_end-time_start, '  seconds'
    endif
    call mpi_finalize(ierr)
    100 Format('', A, F12.2, A, F12.3, A)
    101 Format('', A, F12.3, A)
end















