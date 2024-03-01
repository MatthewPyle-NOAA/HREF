cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     subroutine getconv: compute convection occurrence probability 
c      using Steve W (GSD) method
c   
c     Author: Binbin Zhou, Sept 1, 2009
c       B.Zhou: Chnaged to GRIB2 I/O, July 5, 2013
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine getconv(nv,ipunit,jpdtn,jf,im,jm,est,iens,Lm,Lp,Lt,
     +      gribid,derv_pr,wgt)

        use grib_mod
        include 'parm.inc'

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
 
        Character*5 eps

        common /dtbl/nderiv,
     +              dvname,dk4,dk5,dk6,dMlvl,dPlvl,dTlvl,
     +              dMeanLevel,dProbLevel,dThrs,
     +              dMsignal,dPsignal,MPairLevel,PPairLevel,dop 

        INTEGER, intent(IN) :: nv, jf, iens
        REAL,dimension(jf,Lp,Lt),intent(INOUT) ::  derv_pr

        real CNVPapnt(iens),CNVP(jf),CNVPSMTH(jf)
        real CNVapnt(iens), cnv(jf,iens)
        real wgt(30),threshold(24),aprob
        integer lon_indx,mid,im,jm,est,wst,gribid

        integer miss(iens)
        integer,dimension(iens),intent(IN) :: ipunit
        type (gribfield) :: gfld


        data (threshold(i),i=1,24)
     +  /0.85,0.80,0.775,0.75,0.725,0.70,0.675,0.65,0.675,0.70,
     +  0.725,0.75,0.775,0.80,0.85, 0.90,0.95,1.0,1.02,1.0,0.975,
     +  0.95,0.925,0.9/

        if (gribid.eq.130 .or. gribid.eq.212) then        !only CONUS does smoothing and using different thresholds
          mid = im/2                                    !for its east and west parts
        else
          mid = 4
        end if 

        wst = est - 2
        if (wst.lt.0) wst = wst + 24                    !use mountain region time as west time

        write(*,*) 'In getconv -->'

        write(*,*) 'wgt=',wgt
        write(*,*)'est=',est,'wst=',wst,'mid=',mid
 
        miss=0
        if (jpdtn.eq.0) then
          jpdtnp=8
        else if(jpdtn.eq.1) then
           jpdtnp=11
        else
           jpdtnp=8
        end if
         

        iensloop: do i = 1, iens

         call readGB2(ipunit(i),jpdtnp,1,10,1,0,1,gfld,eps,ie)

          if (ie.eq.0) then
           CNVP(:)=gfld%fld
          else
           miss(i)=1
           cycle iensloop
          end if

         call smooth_points(CNVP,CNVPSMTH,im,jm,mid,jf)

          write(*,*) i, igrid, CNVP(64393),CNVPSMTH(64393) 

         do igrid = 1,jf

          lon_indx = mod(igrid,im)
          if (lon_indx.eq.0) lon_indx=im

          CNVPapnt(i)=CNVPSMTH(igrid)
        
            if(lon_indx.gt.mid) then                !east
             if(CNVPapnt(i).ge.threshold(est)) then 
                CNVapnt(i) = 1.                                          
             else
                CNVapnt(i) = 0.
             endif   
            else
              if(CNVPapnt(i).ge.0.6*threshold(wst)) then
                CNVapnt(i) = 1.
             else
                CNVapnt(i) = 0.
             endif
            end if

            CNV(igrid,i)=CNVapnt(i)
          end do

         end do iensloop


           do 30 lv=1,dPlvl(nv)

             derv_pr(:,lv,:)=0.0

             write(*,*) 'miss=',miss

             do lh = 1, dTlvl(nv)
 
              thr1 = dThrs(nv,lh)
                             
              do igrid = 1,jf

               CNVapnt(:)=CNV(igrid,:)

               call getprob(CNVapnt,iens,
     +                thr1,thr2,dop(nv),aprob,miss,wgt)
          
               derv_pr(igrid,lv,lh)=aprob 
        
               if(igrid.eq.64393) then
                write(*,*) 'Thre', dop(nv),dThrs(nv,lh),':',
     +              CNVapnt,'prob=',aprob
               end if

              end do
             end do

30          continue  


        return
        end

       subroutine smooth_points(var1,var2,im,jm,mid,jf)

        real var1(jf),var2(jf),var2d(im,jm),a
        integer mid

        var2=var1

        do j=1,jm
         do i=1,im
          var2d(i,j)=var1(i+(j-1)*im)
         end do
        end do

        do j=4,jm-3                    !west reqion: 7x7 points smoothing
         do i=4,mid
           a = 0.
           do j1 = j-3, j+3
             do i1 = i-3, i+3
               a = a + var2d(i1,j1)
             end do
           end do
           var2(i+(j-1)*im)=a/49.0
         end do
        end do

       do j=5,jm-4                    !east reqion: 9x9 points smoothing
         do i=mid+1,im-4
           a = 0.
           do j1 = j-4, j+4
             do i1 = i-4, i+4
               a = a + var2d(i1,j1)
             end do
           end do
           var2(i+(j-1)*im)=a/81.0
         end do
       end do

       return
       end
