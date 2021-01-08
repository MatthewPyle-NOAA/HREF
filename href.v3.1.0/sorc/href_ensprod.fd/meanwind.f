cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     subroutine meanwind: compute average wind between two levels 
c     
c     04-17-2014, B. Zhou
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   	subroutine meanwind (nv,ifunit,jpdtn,jf,iens,Lm,Lp,Lt,eps,  
     +           derv_pr,wgt,mbrname)


        use grib_mod
        include 'parm.inc'
        parameter (lvl=23) 

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
        REAL,dimension(jf,Lm,Lt),intent(INOUT) ::  derv_pr
        real apoint(iens)
                                                                
        INTEGER miss(iens) 
        real wgt(30)

        integer,dimension(iens),intent(IN) :: ifunit
        type(gribfield) :: gfld

         REAL, dimension(jf,iens,lvl)::Wpr               
         REAL, dimension(jf,iens):: Upr, Vpr,meanWpr,
     +           meanUpr, meanVpr              
         Integer p(50)
         character*5 eps
         character*7 mbrname(50)

         data (p(i),i=1,lvl)
     +  /850,825,800,775,750,725,700,675,650,625,600,575,550,
     +   525,500,475,450,425,400,375,350,325,300/

CCCCccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Get all paramters from current and previous files 
c   
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        write(*,*) 'In meanwind (850-300mb wind average).....'
        write(*,*) 'ifunit=', ifunit, eps

        write(*,*)  nv,jf,iens,Lm,Lp,Lt

        !jpdtn=0
        jp27=-9999

         miss=0

         Wpr=0.0
         do 600 irun=1,iens
           lvl_loop: do k=1,lvl

          if(mbrname(irun)(1:5).eq.'conus'.and.(p(k).eq.475
     +   .or.p(k).eq.425.or.p(k).eq.375.or.p(k).eq.325)) cycle lvl_loop

          if(mbrname(irun)(1:4).eq.'fv3s'.and.(p(k).eq.475
     +   .or.p(k).eq.425.or.p(k).eq.375.or.p(k).eq.325)) cycle lvl_loop

          if( (mbrname(irun)(1:2).eq.'ak'.or.mbrname(irun)(1:2).eq.'hi'
     +   .or.mbrname(irun)(1:2).eq.'pr') .and.(p(k).eq.475
     +   .or. p(k) .eq. 450 .or.p(k).eq.425.or.p(k).eq.375
     +   .or. p(k) .eq. 350 .or.p(k).eq.325)) cycle lvl_loop

          if( (mbrname(irun)(1:4).eq.'hrrr')
     +    .and.(p(k).eq.475
     +   .or. p(k) .eq. 450 .or.p(k).eq.425.or.p(k).eq.375
     +   .or. p(k) .eq. 350 .or.p(k).eq.325)) cycle lvl_loop

            call readGB2(ifunit(irun),jpdtn,2,2,100,p(k),jp27,gfld,
     +              eps,ie)
                Upr(:,irun)=gfld%fld
            call readGB2(ifunit(irun),jpdtn,2,3,100,p(k),jp27,gfld,
     +              eps,ie)
                Vpr(:,irun)=gfld%fld
            Wpr(:,irun,k)=sqrt(Upr(:,irun)*Upr(:,irun)+
     +                         Vpr(:,irun)*Vpr(:,irun))   
            end do lvl_loop

            meanWpr(:,irun) = 0.0
            do k=1,lvl
              meanWpr(:,irun)=meanWpr(:,irun)+Wpr(:,irun,k)
            end do
            if(mbrname(irun)(1:5).eq.'conus') then
               meanWpr(:,irun)= meanWpr(:,irun)/19
            elseif(mbrname(irun)(1:2).eq.'ak') then
               meanWpr(:,irun)= meanWpr(:,irun)/17
            else
               meanWpr(:,irun)= meanWpr(:,irun)/lvl
            end if

            write(*,*) 'read/compute meanWpr for irun', irun, ' done'

600      continue 


           DO lv=1,dPlvl(nv)-1    !one layer between 2 levels 

             do lh = 1, dTlvl(nv)
                              
              do igrid = 1,jf

                 apoint(:)=meanWpr(igrid,:)

                 if(trim(dop(nv)).ne.'-') then
                  thr1 = dThrs(nv,lh)
                  thr2 = 0.  
                  call getprob(apoint,iens,thr1,thr2,dop(nv),
     +               aprob,miss,wgt)
                  derv_pr(igrid,lv,lh)=aprob

                 else
                  if(lh.lt.dTlvl(nv)) then
                    thr1 = dThrs(nv,lh)
                    thr2 = dThrs(nv,lh+1)
                    call getprob(apoint,iens,thr1,thr2,dop(nv),
     +                   aprob,miss,wgt) 
                    derv_pr(igrid,lv,lh)=aprob
                  end if
                 end if

              end do

             end do
           end do

           return
           end



