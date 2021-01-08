cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     subroutine meanomeg: compute average vertical motion between two levels 
c     
c     04-17-2014, B. Zhou
c     09-09-2016, M. Pyle (adapted from meanwind.f)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   	subroutine meanomeg (nv,ifunit,jpdtn,jf,iens,Lm,Lp,Lt,eps,  
     +           derv_mn,derv_sp,wgt,mbrname)


        use grib_mod
        include 'parm.inc'
        parameter (lvl=9) 

c    for derived variables
        Character*4 dvname(maxvar)
        Integer dk5(maxvar), dk6(maxvar), dk4(maxvar)
        Character*1 dMsignal(maxvar), dPsignal(maxvar)
        Integer dMlvl(maxvar), dMeanLevel(maxvar,maxmlvl)
        Integer dPlvl(maxvar), dProbLevel(maxvar,maxplvl)
        Character*1 dop(maxvar)
        Integer dTlvl(maxvar)
        Real    dThrs(maxvar,maxtlvl)
        Integer MPairLevel(maxvar,maxmlvl,2)
        Integer PPairLevel(maxvar,maxplvl,2)
                                                                                                                                                                
        common /dtbl/nderiv,
     +              dvname,dk4,dk5,dk6,dMlvl,dPlvl,dTlvl,
     +              dMeanLevel,dProbLevel,dThrs,
     +              dMsignal,dPsignal,MPairLevel,PPairLevel,dop


        
        INTEGER, intent(IN) :: nv, jf, iens
        REAL,dimension(jf,Lm),intent(INOUT) ::  derv_mn,derv_sp
        real apoint(iens)
                                                                
        INTEGER miss(iens) 
        real wgt(30)

        integer,dimension(iens),intent(IN) :: ifunit
        type(gribfield) :: gfld

         REAL, dimension(jf,iens,lvl)::OMEGpr           
         REAL, dimension(jf,iens)::  meanOMEGpr
         Integer p(50)
         character*5 eps
         character*7 mbrname(50)

         data (p(i),i=1,lvl)
     +  /700,675,650,625,600,575,550,525,500/

CCCCccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Get all paramters from current and previous files 
c   
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        write(*,*) 'In meanomeg (700-500mb omeg average).....'
        write(*,*) 'ifunit=', ifunit, eps

        write(*,*)  nv,jf,iens,Lm,Lp,Lt

        !jpdtn=0
        jp27=-9999

         miss=0

         OMEGpr=0.0

         do 600 irun=1,iens
           do 500 k=1,lvl
            call readGB2(ifunit(irun),jpdtn,2,8,100,p(k),jp27,gfld,
     +              eps,ie)
                OMEGpr(:,irun,k)=gfld%fld
500        continue

            meanOMEGpr(:,irun) = 0.0
            do k=1,lvl
              meanOMEGpr(:,irun)=meanOMEGpr(:,irun)+OMEGpr(:,irun,k)
            end do
              meanOMEGpr(:,irun)= meanOMEGpr(:,irun)/lvl

            write(*,*) 'read/compute meanOMEGpr for irun', irun, ' done'

600      continue 


	write(0,*) 'dMlvl(nv) for meanomeg: ', dMlvl(nv)

           DO lv=1,dMlvl(nv)-1    !one layer between 2 levels 
              do igrid = 1,jf
                 apoint(:)=meanOMEGpr(igrid,:)

	if (igrid .eq. 1) write(0,*) 'wgt(1:5): ', wgt(1:5)

	if (mod(igrid,100000) .eq. 0) then
	write(0,*) 'apoint: ', igrid, apoint
	endif
                 call getmean(apoint,iens,amean,aspread,miss,wgt)
                 derv_mn(igrid,lv)=amean

	if (mod(igrid,100000) .eq. 0) then
	write(0,*) 'derv_mn: ', igrid,lv,derv_mn(igrid,lv)
	endif

                 derv_sp(igrid,lv)=aspread
              end do
           end do

           return
           end



