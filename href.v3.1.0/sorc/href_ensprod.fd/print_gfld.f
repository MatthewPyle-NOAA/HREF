        subroutine print_gfld(gfld)

        use grib_mod

        type(gribfield) :: gfld
        
        write(*,*) 'Print a gfld info ...........'
        write(*,*) 'GRIB eddition#:',gfld%version
        write(*,*) 'GRIB messege discipline',gfld%discipline
        write(*,*) 'GRIB Section 1 ID:',(gfld%idsect(i),i=1,13)
        write(*,*) 'GRIB Section 1 Array Length:',gfld%idsectlen
        write(*,*) 'GRIB locallen:',(gfld%locallen)
        write(*,*) 'GRIB Local:',(gfld%local(i),i=1,gfld%locallen)
        write(*,*) 'GRIB Field#',gfld%ifldnum
        write(*,*) 'GRIB griddef',gfld%griddef
        write(*,*) 'GRIB Total grids',gfld%ngrdpts
        write(*,*) 'GRIB GRID Def Template 3.#',gfld%igdtnum
        write(*,*) 'GRIB GRID ID Array:',
     +                (gfld%igdtmpl(i),i=1,gfld%igdtlen)
        write(*,*) 'GRIB GRID Temp length:',gfld%igdtlen
        write(*,*) 'GRIB numoct_opt=',gfld%numoct_opt
        write(*,*) 'GRIB interp_opt=',gfld%interp_opt
        write(*,*) 'GRIB Product Template 4.#:',gfld%ipdtnum
        write(*,*) 'GRIB Product ID Array:',
     +                (gfld%ipdtmpl(i),i=1,gfld%ipdtlen)
        write(*,*) 'GRIB Product ID Array length:',gfld%ipdtlen
        write(*,*) 'GRIB idrtnum =', gfld%idrtnum
        write(*,*) 'GRIB idrtlen =', gfld%idrtlen
        write(*,*) 'GRIB gfld%idrtmpl=',gfld%idrtmpl
        write(*,*) 'GRIB num_coord=',gfld%num_coord        
        write(*,'(a10,10f9.2)')'fld value:',
     +             (gfld%fld(i),i=738501,738510)
        return
        end
