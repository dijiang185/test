Module base
    use math
    Implicit None
    Public
    
    type mmp
        Integer :: mpi_i
        Integer :: mpi_n
        Integer, allocatable, dimension(:, :) :: mpi_loc
        Integer, allocatable, dimension(:, :) :: domain, t_domain
    end type mmp
    type vars
        Real, allocatable, dimension(:, :) :: phi, vort, u, v, trans, vort0
    end type vars
    Real, parameter :: omega=7.292e-5, a=6.4e+6, pi=3.1415926
    Real, dimension(:), allocatable :: lat, lon
    Real :: theta, deltat, timax, cen_lon, cen_lat, write_nt, dx
    Real :: SL_its, SL_gdt, eddy_its, V_max
    Integer :: max_it_num

    Real :: r1, r2, r3, r4, zeta1, zeta2, zeta3, vort_amp

    Integer :: nx, ny, test_case, file_type
    Integer :: np_lat, np_lon

    namelist /param/ nx, ny, timax, theta, deltat, write_nt, cen_lon,&
                    cen_lat, file_type, max_it_num, np_lat, np_lon
    namelist /param_case/ r1, r2, r3, r4, zeta1, zeta2, zeta3, vort_amp, test_case
    namelist /sl_nml/ SL_its, SL_gdt, eddy_its, V_max

    contains
        subroutine param_init()
            open(20, file='namelist', form='formatted', status='old', access='sequential')
            read(20, nml=param)
            rewind(20)
            read(20, nml=param_case)
            rewind(20)
            read(20, nml=sl_nml)
            close(20)
            dx = a*pi/180*theta
        end subroutine param_init

        subroutine data_init(mpi_i, mpi_n, var, f, mpi_infos)
            Implicit None
            Integer :: mpi_i, mpi_n
            type(vars) :: var
            type(mmp) :: mpi_infos
            Real, dimension(:, :), allocatable :: f
            
            Integer :: p_nx, p_ny, mpi_i_lat, mpi_i_lon
            Integer :: ie, iw, jn, js
            Integer :: i, j, k
            Integer :: cen_i, cen_j
            Real :: phi, rr, phim
            Real :: x, y
            
            allocate(mpi_infos%mpi_loc(mpi_n, 2))
            allocate(mpi_infos%domain(mpi_n, 4))
            allocate(mpi_infos%t_domain(mpi_n, 4))
            do i = 0, mpi_n-1
                mpi_i_lat = i/np_lon
                mpi_i_lon = mod(i, np_lon)
                p_nx = nx/np_lon ; p_ny = ny/np_lat
                if(mpi_i_lon == 0) then
                    iw = 0
                    ie = p_nx + 1
                elseif(mpi_i_lon == np_lon-1) then
                    iw = mpi_i_lon*p_nx-1; ie = nx+1
                else
                    iw = mpi_i_lon*p_nx-1; ie = (mpi_i_lon+1)*p_nx+1
                endif
                if(mpi_i_lat == 0) then
                    js= 0
                    jn = p_ny + 1
                elseif(mpi_i_lat == np_lat-1) then
                    js = mpi_i_lat*p_ny-1; jn = ny+1
                else
                    js = mpi_i_lat*p_ny-1; jn = (mpi_i_lat+1)*p_ny+1
                endif

                mpi_infos%mpi_loc(i+1, :) = (/mpi_i_lon, mpi_i_lat/)
                mpi_infos%domain(i+1, :) = (/iw, ie, js, jn/)
                mpi_infos%t_domain(i+1, :) = (/iw+1, ie-1, js+1, jn-1/)
                
                if(mpi_i == i) then
                    mpi_infos%mpi_i = mpi_i
                    mpi_infos%mpi_n = mpi_n

                    allocate( var%phi(0:nx+1, 0:ny+1) )
                    allocate( var%vort, mold=var%phi )
                    allocate( var%vort0, mold=var%phi )
                    allocate( var%u, mold=var%phi )
                    allocate( var%v, mold=var%phi )

                    var%phi = 0; var%vort = 0
                    var%u = 0;var%v = 0
                endif
            enddo
            allocate(f(0:nx+1, 0:ny+1))
            do i = 0, nx+1
                if(mod(ny, 2) == 0) then
                    f(i, :) = [(cen_lat+j*theta, j=-ny/2-1, ny/2)]
                else
                    f(i, :) = [(cen_lat+j*theta, j=-ny/2-1, ny/2+1)]
                endif
            enddo
            f = 2*omega*sind(f)
            f = 0.0
            allocate(lat(ny), lon(nx))
            if(mod(nx, 2) == 0) then
                lon = [(cen_lon+i*theta, i=-nx/2, nx/2-1)]
            else
                lon = [(cen_lon+i*theta, i=-nx/2, nx/2)]
            endif
            if(mod(ny, 2) == 0) then
                lat = [(cen_lat+i*theta, i=-ny/2, ny/2-1)]
            else
                lat = [(cen_lat+i*theta, i=-ny/2, ny/2)]
            endif


            if(test_case == 1) then
                if(mod(nx, 2) == 0) then
                    cen_i = nx/2
                else
                    cen_i = nx/2 + 1
                endif
                if(mod(ny, 2) == 0) then
                    cen_j = ny/2
                else
                    cen_j = ny/2 + 1
                endif
                do j = 1, ny
                    do i = 1, nx
                        rr = sqrt((1.0*i-cen_i)**2.0+(1.0*j-cen_j)**2.0)*dx
                        if(rr .eq. 0.0) then
                            var%vort(i, j) = zeta1
                            cycle
                        endif
                        phi = acos((i-cen_i)*dx/rr)
                        if(rr .le. r1) then
                            var%vort(i, j) = zeta1
                        elseif(rr .le. r2) then
                            var%vort(i, j) = zeta1*hermit((rr-r1)/(r2-r1)) + zeta2*hermit((r2-rr)/(r2-r1))
                            do k = 1, 12
                                call random_number(phim)
                                var%vort(i, j) = var%vort(i, j) + vort_amp*cos(phi*k+phim*2.0*pi)*hermit((r2-rr)/(r2-r1))
                            enddo
                        elseif(rr .le. r3) then
                            var%vort(i, j) = zeta2
                            do k = 1, 12
                                call random_number(phim)
                                var%vort(i, j) = var%vort(i, j) + vort_amp*cos(phi*k+phim*2.0*pi)
                            enddo
                        elseif(rr .le. r4) then
                            var%vort(i, j) = zeta2*hermit((rr-r3)/(r4-r3)) + zeta3*hermit((r4-rr)/(r4-r3))
                            do k = 1, 12
                                call random_number(phim)
                                var%vort(i, j) = var%vort(i, j) + vort_amp*cos(phi*k+phim*2.0*pi)*hermit((r4-rr)/(r4-r3))
                            enddo
                        else
                            var%vort(i, j) = zeta3
                        endif
                    enddo
                enddo
                !do j = 1, ny
                !    do i = 1, nx
                !        var%u(i, j) = -SL_its * tanh((lat(j) - cen_lat)*SL_gdt)
                !    enddo
                !enddo
                !var%u(:, 0) = -SL_its * tanh(lat(1)-theta - cen_lat)*SL_gdt
                !var%u(:, ny+1) = -SL_its * tanh(lat(ny)+theta - cen_lat)*SL_gdt
                !var%vort(1:nx, 1:ny) = gradient(var%v(1:nx, 1:ny), axio=1, d=dx) - & 
                !    gradient(var%u(1:nx, 1:ny), axio=2, d=dx)
                !do j = 1, ny
                !    do i = 1, nx
                !        x = (1.0*i-cen_i)*theta
                !        y = (1.0*j-cen_j)*theta
                !        r = sqrt(x**2 + y**2)
                !        if(r <= 2.0) then
                !            var%vort(i, j) = & !var%vort(i, j) + &
                !            eddy_its*exp(-(r**2)/6.0)
                !        endif
                !    enddo
                !enddo
            endif
        end subroutine data_init

        function hermit(s)
            real :: s, hermit
            hermit = 1.0 - 3.0*s**2.0 + 2.0*s**3.0
        end function hermit


End Module base




