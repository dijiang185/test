Module mod_io
    use base
    use netcdf
    Implicit None
    Public 
    Integer :: ncid, phiid, vortid, uid, vid, tid, nt

    contains
        subroutine grd_init()
            open(22, file='vort.ctl', status='unknown')

            write(22, "(A)") 'dset vort.grd'
            write(22, "(A)") 'title vort'
            write(22, "(A)") 'undef -9999'
            write(22, 110  ) 'xdef', nx, 'linear', cen_lon-nx/2*theta, theta
            write(22, 110  ) 'ydef', ny, 'linear', cen_lat-ny/2*theta, theta
            write(22, 111  ) 'xdef', 1, 'levels', 850
            write(22, 112  ) 'tdef', int(timax/write_nt)+1, 'linear', '00Z01JAN1979', 1, 'mn'
            write(22, "(A, I4)") 'vars', 4
            write(22, 113  ) 'phi', 1, 99, 'm**2/s'
            write(22, 113  ) 'vort', 1, 99, '1/s'
            write(22, 113  ) 'u', 1, 99, 'm/s'
            write(22, 113  ) 'v', 1, 99, 'm/s'
            write(22, "(A)")    'endvars'
            close(22)
            110 Format('', A, I6, A10, F10.2, F8.4)
            111 Format('', A, I6, A10, I10)
            112 Format('', A, I6, A10, A15, I4, A)
            113 Format('', A, I4, I4, A8)
        end subroutine grd_init

        subroutine nc_init()
            Implicit None
            Integer :: xid, yid, tvarid, yvarid, xvarid
            Integer :: i, j
            Real :: lat(ny), lon(nx), times(int(timax/write_nt)+1)
            
            nt = 1
            call check(nf90_create('barotropic.nc', NF90_CLOBBER, ncid))
            call check(nf90_def_dim(ncid, 'lon', nx, xid))
            call check(nf90_def_dim(ncid, 'lat', ny, yid))
            call check(nf90_def_dim(ncid, 'time', int(timax/write_nt)+2, tid))

            call check(nf90_def_var(ncid, 'lon', nf90_float, (/xid/), xvarid))
            call check(nf90_def_var(ncid, 'lat', nf90_float, (/yid/), yvarid))
            call check(nf90_def_var(ncid, 'time', nf90_float, (/tid/), tvarid))
            call check(nf90_def_var(ncid, 'phi', nf90_float, (/xid, yid, tid/), phiid))
            call check(nf90_def_var(ncid, 'vort', nf90_float, (/xid, yid, tid/), vortid))
            call check(nf90_def_var(ncid, 'uwnd', nf90_float, (/xid, yid, tid/), uid))
            call check(nf90_def_var(ncid, 'vwnd', nf90_float, (/xid, yid, tid/), vid))

            call check(nf90_put_att(ncid, xid, 'long_name', 'Longitude'))
            call check(nf90_put_att(ncid, xid, 'units', 'degrees_east'))
            call check(nf90_put_att(ncid, yid, 'long_name', 'Latitude'))
            call check(nf90_put_att(ncid, yid, 'units', 'degrees_north'))
            call check(nf90_put_att(ncid, tid, 'long_name', 'Time'))
            call check(nf90_put_att(ncid, tid, 'units', 'minutes since 1979-01-01 00:00'))
            call check(nf90_put_att(ncid, phiid, 'long_name', 'stream unction'))
            call check(nf90_put_att(ncid, phiid, 'units', 'm**2/s'))
            call check(nf90_put_att(ncid, vortid, 'long_name', 'vorticity'))
            call check(nf90_put_att(ncid, vortid, 'units', '1/m'))
            call check(nf90_put_att(ncid, uid, 'long_name', 'uwnd'))
            call check(nf90_put_att(ncid, uid, 'units', 'm/s'))
            call check(nf90_put_att(ncid, vid, 'long_name', 'vwnd'))
            call check(nf90_put_att(ncid, vid, 'units', 'm/s'))
            call check(nf90_enddef(ncid))
            
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
            times = [(write_nt/60*i, i=0, size(times)-1)]
            !print*, lat
            !print*, lon
            !print*, times
            call check(nf90_put_var(ncid, xvarid, lon))
            call check(nf90_put_var(ncid, yvarid, lat(ny:1:-1)))
            call check(nf90_put_var(ncid, tid, times))
             
        end subroutine nc_init

        subroutine write_to_nc(var, varid)
            Implicit None
            Real, dimension(:, :) :: var
            Integer :: varid, i, j
            !do i = 1, nx
            !    call check(nf90_put_var(ncid, varid, var(i, :), start=(/i, 1, nt/), count=(/1, ny, 1/)))
            !enddo
            call check(nf90_put_var(ncid, varid, var, start=(/1, 1, nt/), count=(/nx, ny, 1/)))
        end subroutine write_to_nc

        subroutine close_nc()
            call check(nf90_close(ncid))
        end subroutine close_nc

        subroutine read_nc(vort)
            Implicit None
            Real, dimension(:, :) :: vort
            Integer :: p_ncid, varid
            !call check(nf90_open('var.nc', nf90_nowrite, p_ncid))
            !call check(nf90_inq_varid(p_ncid, 'vort', varid))
            !call check(nf90_get_var(p_ncid, varid, vort))
            call check(nf90_open('restart.nc', nf90_nowrite, p_ncid))
            call check(nf90_inq_varid(p_ncid, 'vort', varid))
            call check(nf90_get_var(p_ncid, varid, vort, start=(/1, 1, 130/), count=(/nx, ny, 1/)))
            vort = vort(:, ny:1:-1)
        end subroutine read_nc
        subroutine check(status)                
            Integer :: status
            If(status /= 0) Then
                Print*, nf90_strerror(status)
                stop
            Endif
        end subroutine check

End Module mod_io


