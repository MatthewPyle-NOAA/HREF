cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Subroutine snow4href: compute 3 hour accumulated snow for HREF
c
c  Author: Binbin Zhou, 2/5/2016
c
c  Input: nv, itime, i00, precip, jf, iens, interval
c  output: derv_mn,  derv_pr
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   	subroutine snow4href(nv,ifunit,ipunit,jf,iens,Lm,Lp,Lt,
     +              derv_mn,derv_sp,derv_pr,eps) 

        use grib_mod
        include 'parm.inc'

        type (gribfield) :: gfld       

c    for derived variables
        Character*4 dvname(maxvar)
        Integer dk5(maxvar), dk6(maxvar),dk4(maxvar)
        Character*1 dMsignal(maxvar), dPsignal(maxvar)
        Integer dMlvl(maxvar), dMeanLevel(maxvar,maxmlvl)
        Integer dPlvl(maxvar), dProbLevel(maxvar,maxplvl)
        Character*1 dop(maxvar)
        Integer dTlvl(maxvar)
        Real    dThrs(maxvar,maxtlvl)
        Integer MPairLevel(maxvar,maxmlvl,2)
        Integer PPairLevel(maxvar,maxplvl,2)
        character*5 eps
        integer miss(iens)
        real    wgt(iens)
                                                                                                                                                                
        common /dtbl/nderiv,
     +              dvname,dk4,dk5,dk6,dMlvl,dPlvl,dTlvl,
     +              dMeanLevel,dProbLevel,dThrs,
     +              dMsignal,dPsignal,MPairLevel,PPairLevel,dop


        integer,dimension(iens),intent(IN) :: ifunit,ipunit
        INTEGER, intent(IN) :: nv,iens
        REAL,dimension(jf) :: apcp
        REAL,dimension(jf) :: ice_type
        REAL,dimension(jf) :: snw_type
        REAL,dimension(jf,iens) :: snow
        REAL,dimension(jf,Lm),intent(INOUT) :: derv_mn,derv_sp
        REAL,dimension(jf,Lp,Lt),intent(INOUT) :: derv_pr

        REAL,dimension(iens)  :: apoint
        REAL                  :: amean,aspread,aprob
        miss=0
        wgt=1.0 

        write(*,*) 'In snow4href:',eps,dvname(nv),dMlvl(nv),dPlvl(nv)
        write(*,*)  nv,jf,iens,Lm,Lp,Lt,dk4(nv),dk5(nv),dk6(nv)
        
        snow=0.
        do 400 irun=1,iens

          !write(*,*) 'ifunit, ipunit=',ifunit(irun),ipunit(irun)
          call readGB2(ifunit(irun),0,1,195,1,0,-9999,gfld,eps,ie) !CSNOW
           if(ie.eq.0) then
            snw_type(:)=gfld%fld(:)
            !write(*,*) 'read CSNOW done for :',irun
           else
            write(*,*) 'read CSNOW has error:',ie
            stop 44
           end if

          call readGB2(ifunit(irun),0,1,194,1,0,-9999,gfld,eps,ie) !CICEP
            if(ie.eq.0) then
            ice_type(:)=gfld%fld(:)
            !write(*,*) 'read CICEP done for :',irun
           else
            write(*,*) 'read CICEP has error:',ie
            stop 45
           end if

          if (dvname(nv).eq.'SN3h') then
             jpd27=3
          else if (dvname(nv).eq.'SN6h') then
             jpd27=6
          else if (dvname(nv).eq.'SN12') then
             jpd27=12
          else if (dvname(nv).eq.'SN24') then
             jpd27=24
          else
             write(*,*) 'Wrong Snowfall request name'
             stop 5555
          end if
         
          call readGB2(ipunit(irun),8,1,8,1,0,jpd27,gfld,eps,ie) !APCP
           if(ie.eq.0) then
            apcp(:)=gfld%fld(:)
            !write(*,*) 'read APCP done for :',irun
           else
            write(*,*) 'read APCP has error:',ie
            stop 46
           end if

           do i=1,jf
            if (snw_type(i).eq.1.0.or.ice_type(i).eq.1.0) then
             snow(i,irun)=apcp(i)*10.
            end if
           end do
 
400        continue
 

             do lv=1,dMlvl(nv)                               
               do igrid = 1,jf
                 apoint(:)=snow(igrid,:) 
                 call getmean(apoint,iens,amean,aspread,miss,wgt)
                 derv_mn(igrid,lv) = amean
                 derv_sp(igrid,lv) = aspread
                end do
              write(*,*) 'mean/spread done'
             end do


             do lv=1,dPlvl(nv)
                do lt =1,dTlvl(nv)
                  do igrid = 1, jf
                   apoint(:)=snow(igrid,:)
                    if(trim(dop(nv)).ne.'-') then
                     thr1 = dThrs(nv,lt)
                     thr2 = 0.
               call getprob(apoint,iens,thr1,thr2,dop(nv),aprob,
     +                  miss,wgt)
                     derv_pr(igrid,lv,lt)=aprob
                    else
                     if(lt.lt.dTlvl(nv)) then
                       thr1 = dThrs(nv,lt)
                       thr2 = dThrs(nv,lt+1)
               call getprob(apoint,iens,thr1,thr2,dop(nv),aprob,
     +                  miss,wgt)
                       derv_pr(igrid,lv,lt)=aprob
                     end if
                    end if
                   
                  end do
                 write(*,*) 'prob done for', lv,lt
                end do
             end do
         

           return
           end
