          program ffg_gen


! Main program that just calls the two subpieces needed to generate
! the final FFG product

! M. Pyle - 20190502 - initial version

        integer :: itag
        character(len=8) PDY

        read(5,*) itag
        read(5,'(A8)') PDY

	if (itag .eq. 1) then
	write(0,*) 'call read_grib1'
        call read_grib1(PDY)
	elseif (itag .eq. 2) then
	write(0,*) 'call stitch'
        call stitch(PDY)
        else
	write(0,*) 'bad itag: ', itag

	endif

	end program ffg_gen
