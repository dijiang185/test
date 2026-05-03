Module Math
        Interface gradient
                module procedure gradient_1d
                module procedure gradient_2d
                module procedure gradient_3d
                module procedure gradient_4d
        End interface gradient
        
        contains
                Function gradient_1d(x, d) result(ans)
                        Implicit None
                        Real, optional :: d
                        Real :: dd
                        Real, dimension(:) :: x
                        Real, allocatable, dimension(:) :: ans
                        Integer :: n
                        Integer :: i

                        n = ubound(x, 1)
                        allocate(ans(n))
                        If(present(d)) Then
                                dd = d
                        Else
                                dd = 1
                        Endif
                        If(dd == 0. .or. dd < 1e-6) Then
                               ans = -9999
                        Else 
                                Do i = 1, n
                                        If (i == 1) Then
                                                ans(i) = (x(i+1) - x(i))/dd
                                        Elseif(i == n) Then
                                                ans(i) = (x(i) - x(i-1))/dd
                                        Else
                                                ans(i) = (x(i+1) - x(i-1))/(2*dd)
                                        Endif
                                Enddo
                        Endif
                End Function gradient_1d
                
                Function gradient_2d(x, axio, d) result(ans)
                        Implicit None
                        Integer  :: d1, d2
                        Integer, optional :: axio
                        Integer :: axios
                        Real, optional :: d
                        Real :: dd
                        Real, dimension(:, :) :: x
                        Real, allocatable, dimension(:, :) :: ans
                        Integer :: i

                        d1 = ubound(x, 1)
                        d2 = ubound(x, 2)
                        allocate(ans(d1, d2))
                        If(present(axio)) Then
                                axios = axio
                        Else
                                axios = 1
                        Endif
                        If(present(d)) Then
                                dd = d
                        Else
                                dd = 1.0
                        Endif
            
                        If(axios /= 1 .and. axios /= 2) Then
                                Print*, 'axios must bs 1 or 2'
                                stop
                        Endif
                        If(axios == 1) Then
                                Do i = 1, d2
                                        ans(:, i) = gradient_1d(x(:, i), d)
                                Enddo
                        Else
                                Do i = 1, d1
                                        ans(i, :) = gradient_1d(x(i, :), d)
                                Enddo
                        Endif


                End Function gradient_2d

                Function gradient_3d(x, axio, d) result(ans)
                        Implicit None
                        Real, dimension(:, :, :) :: x
                        Real, allocatable, dimension(:, :, :) :: ans
                        Integer, optional :: axio
                        Real, optional :: d
                        Integer :: d1, d2, d3, axios
                        Real :: dd
                        Integer :: i, j

                        If(present(axio)) Then
                                axios = axio
                        Else
                                axios = 1
                        Endif

                        If(present(d)) Then
                                dd = d
                        Else
                                dd = 1.0
                        Endif

                        d1 = ubound(x, 1)
                        d2 = ubound(x, 2)
                        d3 = ubound(x, 3)
                        allocate(ans(d1, d2, d3))

                        If(axios == 1) Then
                                Do i = 1, d2
                                        do j = 1, d3
                                                ans(:, i, j) = gradient_1d(x(:, i, j), dd)
                                        Enddo
                                Enddo
                        Elseif(axios == 2) Then
                                Do i = 1, d1
                                        do j = 1, d3
                                                ans(i, :, j) = gradient_1d(x(i, :, j), dd)
                                        Enddo
                                Enddo
                        Else
                                Do i = 1, d1
                                        do j = 1, d2
                                                ans(i, j, :) = gradient_1d(x(i, j, :), dd)
                                        Enddo
                                Enddo
                        Endif
                End Function gradient_3d
                
                Function gradient_4d(x, axio, d) result(ans)
                        Implicit None
                        Real, dimension(:, :, :, :) :: x
                        Real, allocatable, dimension(:, :, :, :) :: ans
                        Integer, optional :: axio
                        real, optional :: d
                        Integer :: i, j, k, axios
                        Integer :: d1, d2, d3, d4
                        Real :: dd

                        If(present(axio)) Then
                                axios = axio
                        Else
                                axios = 1
                        Endif

                        If(present(d)) Then
                                dd = d
                        Else
                                dd = 1.0
                        Endif

                        d1 = ubound(x, 1)
                        d2 = ubound(x, 2)
                        d3 = ubound(x, 3)
                        d4 = ubound(x, 4)
                        allocate(ans(d1, d2, d3, d4))
                        If(axios == 1) Then
                                Do i = 1, d2
                                        do j = 1, d3
                                                Do k = 1, d4
                                                        ans(:, i, j, k) = gradient_1d(x(:, i, j, k), dd)
                                                Enddo
                                        Enddo
                                Enddo
                        Elseif(axio == 2) Then
                                Do i = 1, d1
                                        do j = 1, d3
                                                Do k  =1, d4
                                                        ans(i, :, j, k) = gradient_1d(x(i, :, j, k), dd)
                                                Enddo
                                        Enddo
                                Enddo
                        Elseif(axio == 3) Then
                                Do i = 1, d1
                                        do j = 1, d2
                                                Do k = 1, d4
                                                        ans(i, j, :, k) = gradient_1d(x(i, j, :, k), dd)
                                                Enddo
                                        Enddo
                                Enddo
                        Else
                                Do i = 1, d1
                                        do j = 1, d2
                                                Do k = 1, d3
                                                        ans(i, j, k, :) = gradient_1d(x(i, j, k, :), dd)
                                                Enddo
                                        Enddo
                                Enddo
                        Endif
                End Function gradient_4d 
                
                Function cov(x, y) result(ans)
                        Implicit None
                        Real, dimension(:) :: x, y
                        Real :: ans, Ex, Ey
                        Integer :: n

                        n = ubound(x, 1)
                        If(n /= ubound(y, 1)) Then
                                Print*, 'The length of x and y must be same.'
                                stop
                        Endif
                        Ex = sum(x)/n
                        Ey = sum(y)/n
                        ans = sum((x - Ex)*(x - Ey))/(n - 1)
                End Function cov

                Function corr(x, y) result(ans)
                        Implicit None
                        Real, dimension(:) :: x, y
                        Real :: ans
                        
                        ans = cov(x, y)/sqrt(cov(x, x)*cov(y, y))
                End Function corr

                Function regression(x, y) result(ans)
                        Implicit None
                        Real, dimension(:) :: y, x
                        Real, dimension(2) :: ans
                        Real :: Ex, Ey
                        Integer :: n

                        n = ubound(x, 1)
                        Ex = sum(x)/n
                        Ey = sum(y)/n

                        ans(1) = (sum(x*y)/n - Ex*Ey)/(sum(x**2)/n - Ex**2)
                        ans(2) = Ey - ans(1)*Ex
                End Function regression

                Subroutine write_matrix(x)
                        Implicit None
                        Real, dimension(:, :) :: x
                        Integer :: i, n

                        n = ubound(x, 1)
                        do i = 1, n
                                print*, x(i, :)
                        enddo
                End Subroutine write_matrix
                
                Recursive Function det(x) result(ans)
                        Implicit None
                        Real, dimension(:, :) :: x
                        Real :: ans
                        Real, allocatable, dimension(:, :) :: xx
                        Integer :: n
                        Integer :: i, j

                        n = ubound(x, 1)
                        If(n .ne. ubound(x, 2)) then
                                print*, 'dim1 not equal dim2...'
                                stop
                        Endif

                        If(n .eq. 1) Then
                                ans = x(1, 1)
                        Else
                                allocate(xx(n-1, n-1))
                                ans = 0
                                Do i = 1, n
                                        if(i .eq. 1) then
                                                xx(:, :) = x(2:n, 2:n)
                                        elseif(i .eq. n) then
                                                xx(:, :) = x(2:n, 1:n-1)
                                        else
                                                xx(:, 1:i-1) = x(2:n, 1:i-1)
                                                xx(:, i:n-1) = x(2:n, i+1:n)
                                        endif
                                        ans = ans + (-1)**(1+i)*x(1, i)*det(xx)
                                Enddo
                        Endif
                        End Function det

                                

End Module Math


