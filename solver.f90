Module solver
    use math
    use base
    use mpi
    use mpi_routine
    Implicit None
    Integer, parameter :: nrk=3

    contains
        subroutine solver2(mpi_infos, var, f)
            Implicit None
            type(mmp) :: mpi_infos
            type(vars), target :: var
            Real, dimension(0:nx+1, 0:ny+1) :: f
            
            Real, dimension(:, :), pointer :: phi, vort, u, v, vort0, vort1
            Real, dimension(:, :), pointer :: tdc
            Integer :: mpi_i, mpi_n
            Integer :: iw, ie, js, jn
            Integer :: i, j 
            Real :: gamma=0.1
            
            phi => var%phi
            vort => var%vort
            vort0 => var%vort0
            u => var%u
            v => var%v
            mpi_i = mpi_infos%mpi_i
            mpi_n = mpi_infos%mpi_n
            iw = mpi_infos%domain(mpi_i+1, 1)
            ie = mpi_infos%domain(mpi_i+1, 2)
            js = mpi_infos%domain(mpi_i+1, 3)
            jn = mpi_infos%domain(mpi_i+1, 4)
            allocate(vort1, mold=vort)
            allocate(tdc, mold=vort)
            call jacobina(vort(iw:ie, js:jn)+f(iw:ie, js:jn), phi(iw:ie, js:jn), tdc(iw:ie, js:jn), 7)
            call exchange(mpi_infos, tdc)
            vort1(iw:ie, js:jn) = vort0(iw:ie, js:jn) + 2*deltat*tdc(iw:ie, js:jn)
            vort0(iw:ie, js:jn) = vort(iw:ie, js:jn) + gamma*(vort0(iw:ie, js:jn) - &
                                  2*vort(iw:ie, js:jn) + vort1(iw:ie, js:jn))
            vort(iw:ie, js:jn) = vort1(iw:ie, js:jn)
            call boundary(mpi_infos, vort)
            call relaxtion(mpi_infos, var)

            !print*, mpi_i, maxval(vort(iw:ie, js:jn)), minval(vort(iw:ie, js:jn))
            deallocate(vort1, tdc)
        end subroutine solver2

        subroutine solver1(mpi_infos, var, f)
            Implicit None
            type(mmp), target :: mpi_infos
            type(vars), target :: var
            Real, dimension(0:nx+1, 0:ny+1) :: f
            
            type(vars) :: w
            Real, dimension(:, :), pointer :: phi, vort
            Real, dimension(:, :), pointer :: tdc
            Real, dimension(:), pointer :: wight
            Real, dimension(:, :, :), pointer :: k
            Integer :: mpi_i, mpi_n, iw, ie, js, jn, p_nx, p_ny
            Integer :: i, j
            
            phi => var%phi
            vort => var%vort
            mpi_i = mpi_infos%mpi_i
            mpi_n = mpi_infos%mpi_n
            iw = mpi_infos%domain(mpi_i+1, 1)
            ie = mpi_infos%domain(mpi_i+1, 2)
            js = mpi_infos%domain(mpi_i+1, 3)
            jn = mpi_infos%domain(mpi_i+1, 4)
            
            allocate(wight(nrk))
            allocate(k(nrk, 0:nx+1, 0:ny+1))
            allocate(w%phi, w%u, w%v, mold=phi)
            allocate(w%vort, mold=phi)
            select case(nrk)
            case(3)
                wight = (/1.0, 4.0, 1.0/)
                call jacobina(vort(iw:ie, js:jn)+f(iw:ie, js:jn), phi(iw:ie, js:jn), k(1, iw:ie, js:jn), 7)
                call exchange(mpi_infos, k(1, :, :))

                w%vort(iw:ie, js:jn) = vort(iw:ie, js:jn) + k(1, iw:ie, js:jn)*deltat/2
                call relaxtion(mpi_infos, w)
                call jacobina(w%vort(iw:ie, js:jn)+f(iw:ie, js:jn), w%phi(iw:ie, js:jn), k(2, iw:ie, js:jn), 7)
                call exchange(mpi_infos, k(2, :, :))

                w%vort(iw:ie, js:jn) = w%vort(iw:ie, js:jn) + 2*k(2, iw:ie, js:jn)*deltat - 3/2*k(1, iw:ie, js:jn)*deltat
                call relaxtion(mpi_infos, w)
                call jacobina(w%vort(iw:ie, js:jn)+f(iw:ie, js:jn), w%phi(iw:ie, js:jn), k(3, iw:ie, js:jn), 7)
                call exchange(mpi_infos, k(3, :, :))

                do i = 1, nrk
                    vort(iw:ie, js:jn) = vort(iw:ie, js:jn) + wight(i)*k(i, iw:ie, js:jn)*deltat/6
                enddo
                call boundary(mpi_infos, vort)
                call relaxtion(mpi_infos, var)
            case(4)
                wight = (/1.0, 2.0, 2.0, 1.0/)
                call jacobina(vort(iw:ie, js:jn)+f(iw:ie, js:jn), phi(iw:ie, js:jn), k(1, iw:ie, js:jn), 7)
                call exchange(mpi_infos, k(1, :, :))
                
                w%vort(iw:ie, js:jn) = vort(iw:ie, js:jn) + k(1, iw:ie, js:jn)*deltat/2
                call singal_to_master(mpi_infos, k(1, :, :))
                call relaxtion(mpi_infos, w)
                call jacobina(w%vort(iw:ie, js:jn)+f(iw:ie, js:jn), w%phi(iw:ie, js:jn), k(2, iw:ie, js:jn), 7)
                call exchange(mpi_infos, k(2, :, :))

                w%vort(iw:ie, js:jn) = w%vort(iw:ie, js:jn) + k(2, iw:ie, js:jn)*deltat/2
                call relaxtion(mpi_infos, w)
                call jacobina(w%vort(iw:ie, js:jn)+f(iw:ie, js:jn), w%phi(iw:ie, js:jn), k(3, iw:ie, js:jn), 7)
                call exchange(mpi_infos, k(3, :, :))

                w%vort(iw:ie, js:jn) = w%vort(iw:ie, js:jn) + k(3, iw:ie, js:jn)*deltat
                call relaxtion(mpi_infos, w)
                call jacobina(w%vort(iw:ie, js:jn)+f(iw:ie, js:jn), w%phi(iw:ie, js:jn), k(4, iw:ie, js:jn), 7)
                call exchange(mpi_infos, k(4, :, :))
                
                do i = 1, nrk
                    vort(iw:ie, js:jn) = vort(iw:ie, js:jn) + wight(i)*k(i, iw:ie, js:jn)*deltat/6
                enddo
                call relaxtion(mpi_infos, var)
            end select
            deallocate(wight)
            deallocate(k)
            deallocate(w%phi)
            deallocate(w%vort)
        end subroutine solver1
        
        subroutine phi_to_uv(mpi_infos, var)
            type(mmp) :: mpi_infos
            type(vars), target :: var
            Real, dimension(:, :), pointer :: phi, u, v
            Integer :: mpi_i, mpi_n, iw, ie, js, jn
            Integer :: i, j
            
            mpi_i = mpi_infos%mpi_i
            mpi_n = mpi_infos%mpi_n
            iw = mpi_infos%domain(mpi_i+1, 1)
            ie = mpi_infos%domain(mpi_i+1, 2)
            js = mpi_infos%domain(mpi_i+1, 3)
            jn = mpi_infos%domain(mpi_i+1, 4)
            u => var%u
            v => var%v
            phi => var%phi

            u(iw:ie, js:jn) = -gradient(phi(iw:ie, js:jn), axio=2, d=dx)
            v(iw:ie, js:jn) =  gradient(phi(iw:ie, js:jn), axio=1, d=dx)

        end subroutine phi_to_uv

        subroutine relaxtion(mpi_infos, var)
            type(vars), target :: var
            type(mmp), target :: mpi_infos

            Real, dimension(:, :), pointer :: phi, vort, u
            Real, allocatable, dimension(:, :) :: pre_phi
            Real, dimension(:), pointer :: trans_lat, trans_lon
            Integer :: mpi_i, mpi_n, iw, ie, js, jn
            Integer, dimension(:), pointer ::mpi_i_lon, mpi_i_lat
            Integer :: ierr, status(mpi_status_size)
            Integer :: i, j
            Integer :: p_nx, p_ny, n, total, totals
            Real :: r, epsilon=0.05

            mpi_i = mpi_infos%mpi_i
            mpi_n = mpi_infos%mpi_n
            mpi_i_lon => mpi_infos%mpi_loc(:, 1)
            mpi_i_lat => mpi_infos%mpi_loc(:, 2)
            
            phi => var%phi; vort => var%vort; u => var%u
            phi = 0
            allocate(pre_phi, mold=phi)
            pre_phi = 0
            iw = mpi_infos%t_domain(mpi_i+1, 1)
            ie = mpi_infos%t_domain(mpi_i+1, 2)
            js = mpi_infos%t_domain(mpi_i+1, 3)
            jn = mpi_infos%t_domain(mpi_i+1, 4)
            p_nx = ie - iw + 1 
            p_ny = jn - js + 1
            n = 1; total = 0
            totals = (nx-2)*(ny-2) + (np_lat - 1)*(nx-2) + (np_lon - 1)*(ny-2) + (np_lon-1)*(np_lat - 1)
            do while(total < totals .and. n <= max_it_num)
                do i = iw, ie
                    do j = js, jn
                        pre_phi(i, j) = (phi(i+1, j) + phi(i-1, j) + phi(i, j+1) + phi(i, j-1))/(dx*dx) - 4/(dx*dx)*phi(i, j) - vort(i, j)
                    enddo
                enddo
                call exchange(mpi_infos, pre_phi)
                pre_phi(iw-1:ie+1, js-1:jn+1) = dx*dx/4*pre_phi(iw-1:ie+1, js-1:jn+1)
                phi(iw-1:ie+1, js-1:jn+1) = phi(iw-1:ie+1, js-1:jn+1) + pre_phi(iw-1:ie+1, js-1:jn+1)
                call boundary(mpi_infos, phi)
                !phi(:, 0) = phi(:, 1)+(u(:, 0)+u(:, 1))/2.0*dx
                !phi(:, ny+1) = phi(:, ny)-(u(:, ny)+u(:, ny+1))/2.0*dx
                !phi(:, 0) = 0.0
                !phi(:, ny+1) = 0.0
                call mpi_allreduce(count(abs(pre_phi(iw:ie, js:jn)/phi(iw:ie, js:jn)) < epsilon), total, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)
                n = n + 1
                !if(mpi_i == 0) print*, mpi_i, phi(iw, jn)
                !if(mpi_i == 3) print*, mpi_i, phi(iw, js)
                !if(mpi_i == 2) print*, mpi_i, phi(ie, jn)
                !if(mpi_i == 5) print*, mpi_i, phi(ie, js)
                if(n >= max_it_num) then
                    if(mpi_i == 0) &
                    print*, 'relaxtion too long....', total, totals
                    !call close_nc
                    !stop
                endif
            enddo
            deallocate(pre_phi)
        end subroutine relaxtion

        !subroutine boundary(vort)
        !    Real, dimension(0:nx+1, 0:ny+1) :: vort

        !    vort(:, 0) = 0.0
        !    vort(:, ny+1) = 0.0
        !end subroutine boundary
        subroutine boundary(mpi_infos, phi)
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
            
            !!!!!!! cycle boundary in zonal derection
            if(.not. left .and. upper) then
                call mpi_send(phi(iw+1, jn-1), 1, mpi_real, mpi_i+2*np_lon-1, 6*mpi_n+mpi_i, mpi_comm_world, ierr)
                call mpi_recv(phi(iw-1, jn+1), 1, mpi_real, mpi_i+2*np_lon-1, 6*mpi_n+(mpi_i+2*np_lon-1), mpi_comm_world, status, ierr)
            endif
            if(.not. left .and. lower) then
                call mpi_send(phi(iw+1, js+1), 1, mpi_real, mpi_i-1, 6*mpi_n+mpi_i, mpi_comm_world, ierr)
                call mpi_recv(phi(iw-1, js-1), 1, mpi_real, mpi_i-1, 6*mpi_n+(mpi_i-1), mpi_comm_world, status, ierr)
            endif
            if(.not. right .and. upper) then
                call mpi_send(phi(ie-1, jn-1), 1, mpi_real, mpi_i+1, 6*mpi_n+mpi_i, mpi_comm_world, ierr)
                call mpi_recv(phi(ie+1, jn+1), 1, mpi_real, mpi_i+1, 6*mpi_n+(mpi_i+1), mpi_comm_world, status, ierr)
            endif
            if(.not. right .and. lower) then
                call mpi_send(phi(iw-1, js+1), 1, mpi_real, mpi_i-2*np_lon+1, 6*mpi_n+mpi_i, mpi_comm_world, ierr)
                call mpi_recv(phi(ie+1, js-1), 1, mpi_real, mpi_i-2*np_lon+1, 6*mpi_n+(mpi_i-2*np_lon+1), mpi_comm_world, status, ierr)
            endif
            !call mpi_barrier(mpi_comm_world, ierr)

            if(.not. left) then
                call mpi_send(reshape(phi(iw+1, js:jn), (/p_ny/)), p_ny, mpi_real, mpi_i+np_lon-1, mpi_i, mpi_comm_world, ierr)
                allocate(trans_lon(p_ny))
                call mpi_recv(trans_lon, p_ny, mpi_real, mpi_i+np_lon-1, mpi_i+np_lon-1, mpi_comm_world, status, ierr)
                phi(iw-1, js:jn) = trans_lon
                deallocate(trans_lon)
            endif
            if(.not. right) then
                call mpi_send(reshape(phi(ie-1, js:jn), (/p_ny/)), p_ny, mpi_real, mpi_i-np_lon+1, mpi_i, mpi_comm_world, ierr)
                allocate(trans_lon(p_ny))
                call mpi_recv(trans_lon, p_ny, mpi_real, mpi_i-np_lon+1, mpi_i-np_lon+1, mpi_comm_world, status, ierr)
                phi(ie+1, js:jn) = trans_lon
                deallocate(trans_lon)
            endif
            !!!!!!!!!!!!!!!!!!
        end subroutine boundary 

        
        subroutine smooth(mpi_infos, vort, s)
            Implicit None
            type(mmp) :: mpi_infos
            Real, dimension(0:nx+1, 0:ny+1) :: vort
            Real :: s
            
            Real, dimension(:, :), allocatable :: w
            Integer :: mpi_i
            Integer :: iw, ie, js, jn
            Integer :: i, j

            mpi_i = mpi_infos%mpi_i
            iw = mpi_infos%t_domain(mpi_i+1, 1)
            ie = mpi_infos%t_domain(mpi_i+1, 2)
            js = mpi_infos%t_domain(mpi_i+1, 3)
            jn = mpi_infos%t_domain(mpi_i+1, 4)
            
            allocate(w(iw:ie, js:jn))
            do i = iw, ie
                do j = js, jn
                    w(i, j) = vort(i, j) + s/4*(vort(i-1, j) + vort(i+1, j) + vort(i, j-1) + vort(i, j+1) - 4*vort(i, j))
                enddo
            enddo
            vort(iw:ie, js:jn) = w
            call exchange(mpi_infos, vort)

        end subroutine smooth

        subroutine jacobina(f1, f2, jcb, order)
            Implicit None
            Real, dimension(:, :) :: f1, f2, jcb
            Integer :: order
            
            Real, dimension(:, :), allocatable :: jcb1, jcb2, jcb3
            Integer :: m, n

            if(order < 1 .or. order > 7) then
                print*, 'order should be 1~7...'
                stop
            endif
            select case(order)
            case(1)
                allocate(jcb1, mold=f1)
                jcb1 = jacobina1()
                jcb = jcb1
            case(2)
                allocate(jcb2, mold=f1)
                jcb2 = jacobina2()
                jcb = jcb2
            case(3)
                allocate(jcb3, mold=f1)
                jcb3 = jacobina3()
                jcb = jcb3
            case(4)
                allocate(jcb1, mold=f1)
                jcb1 = jacobina1()
                allocate(jcb2, mold=f1)
                jcb2 = jacobina2()
                jcb = (jcb1 + jcb2)/2
            case(5)
                allocate(jcb2, mold=f1)
                jcb2 = jacobina2()
                allocate(jcb3, mold=f1)
                jcb3 = jacobina3()
                jcb = (jcb2 + jcb3)/2
            case(6)
                allocate(jcb1, mold=f1)
                jcb1 = jacobina1()
                allocate(jcb3, mold=f1)
                jcb3 = jacobina3()
                jcb = (jcb1 + jcb3)/2
            case(7)
                allocate(jcb1, mold=f1)
                jcb1 = jacobina1()
                allocate(jcb2, mold=f1)
                jcb2 = jacobina2()
                allocate(jcb3, mold=f1)
                jcb3 = jacobina3()
                jcb = (jcb1 + jcb2 + jcb3)/3
            end select

            contains
                function jacobina1()
                    Implicit None
                    Real, dimension(:, :), allocatable :: jacobina1
                     
                    allocate(jacobina1, mold=f1)
                    jacobina1 = gradient(f1, axio=1, d=dx)*gradient(f2, axio=2, d=dx) &
                                - gradient(f1, axio=2, d=dx)*gradient(f2, axio=1, d=dx)
                end function jacobina1
                function jacobina2()
                    Implicit None
                    Real, dimension(:, :), allocatable :: jacobina2

                    allocate(jacobina2, mold=f1)
                    jacobina2 = gradient(f2*gradient(f1, axio=1, d=dx), axio=2, d=dx) - &
                                gradient(f2*gradient(f1, axio=2, d=dx), axio=1, d=dx) 
                end function jacobina2
                function jacobina3()
                    Implicit None
                    Real, dimension(:, :), allocatable :: jacobina3
                    
                    allocate(jacobina3, mold=f1)
                    jacobina3 = gradient(f1*gradient(f2, axio=2, d=dx), axio=1, d=dx) - &
                                gradient(f1*gradient(f2, axio=1, d=dx), axio=2, d=dx)
                end function jacobina3
                

        end subroutine jacobina
        

End Module solver
