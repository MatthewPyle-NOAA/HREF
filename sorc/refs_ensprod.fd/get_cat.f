cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     subroutine get_cat: compute Clear-Air-Turbulence based on Ellrod Algo (1992)
c     
c     Author: Binbin Zhou, Apr, 25, 2009
c      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   	subroutine get_cat (nv,ifunit,jpdtn,jf,iens,Lp,Lt, 
     +              derv_pr,nx,ny,dx,dy,wgt)

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
        character*5 eps                                                                                                                                                               
        common /dtbl/nderiv,
     +              dvname,dk4,dk5,dk6,dMlvl,dPlvl,dTlvl,
     +              dMeanLevel,dProbLevel,dThrs,
     +              dMsignal,dPsignal,MPairLevel,PPairLevel,dop


        INTEGER, intent(IN) :: nv, jf, iens
        REAL,dimension(jf,Lp,Lt),intent(INOUT) :: derv_pr

        real wgt(30)

        integer miss(iens)

        real p(Lp),apoint_cat(iens)
        real,dimension(jf,iens,Lp-1) :: cat3d
        real,dimension(jf,Lp) :: u,v,h
        real,dimension(jf,Lp-1) :: cat

        integer,dimension(iens),intent(IN) :: ifunit
        type(gribfield) :: gfld

       
        jpd10=dk6(nv)
        !jpdtn=0
        jp27=-9999
 
        write(*,*) 'In get_cat .....'
        write(*,*) 'nv,ifunit,jf,iens,Lp,Lt,jpd10,nx,ny,dx,dy',
     +              nv,ifunit,jf,iens,Lp,Lt,jpd10,nx,ny,dx,dy

        do k=1,dPlvl(nv)
          p(k)=dProbLevel(nv,k)
        end do
        write(*,*) 'p(k)',p
        
        do i=1,iens
           do k=1,dPlvl(nv)
              jpd12=dProbLevel(nv,k)
              
             call readGB2(ifunit(i),jpdtn,2,2,jpd10,jpd12,jp27,gfld,
     +            eps,iret)   !u
              if(iret.eq.0) then 
                u(:,k)=gfld%fld
              else
                write(*,*)'U missing at level',jpd12,'in file',ifunit(i)
                stop 222
              end if
             call readGB2(ifunit(i),jpdtn,2,3,jpd10,jpd12,jp27,gfld,
     +            eps,iret)   !v
              if(iret.eq.0) then
                v(:,k)=gfld%fld
              else
                write(*,*)'V missing at level',jpd12,'in file',ifunit(i)
                stop 223
              end if

             call readGB2(ifunit(i),jpdtn,3,5,jpd10,jpd12,jp27,gfld,
     +            eps,iret)   !h
               if(iret.eq.0) then
                 h(:,k)=gfld%fld
               else
                write(*,*)'H missing at level',jpd12,'in file',ifunit(i)
                stop 224
              end if

           end do
                
           ms=0

           call get_turb(i,u,v,h,p,jf,Lp,nx,ny,dx,dy,cat,ms)

           do k=1,Lp-1
             cat3d(:,i,k)=cat(:,k)
           end do
         end do

           miss=0  !No any missing
             do lv=1,dPlvl(nv)-1 
              do lh = 1, dTlvl(nv)      
               do igrid=1,jf
 
                  apoint_cat=cat3d(igrid,:,lv)
                      
                  thr1 = dThrs(nv,lh)
                  thr2 = 0.
                  call getprob(apoint_cat,iens,thr1,thr2,dop(nv),aprob,
     +               miss,wgt)
                  derv_pr(igrid,lv,lh)=aprob

                end do
               end do
              end do
         
         return
         end

      subroutine get_turb(ie,u,v,h,p,jf,Lp,nx,ny,dx,dy,cat,ms)
  
      REAL,dimension(jf,Lp),intent(IN) :: u,v,h
      REAL,dimension(jf,Lp-1),intent(INOUT) :: cat
      INTEGER p(Lp),  jf,nx,ny, ms

      REAL,dimension(nx,ny) :: u0,v0,h0,w0,u1,v1,h1,w1
      REAL dx,dy,tb
      REAL f1,f2       !Note f1 and f2 are two tuned numbers
                       !the higher resolution, the smaller f1, f2

      f1=2500.
      f2=1000.

      cat = 0.
      do 600 k=2,Lp

       do j=1,ny
        do i = 1,nx
          igrid=i+(j-1)*nx
          u0(i,j)=u(igrid,k)
          v0(i,j)=v(igrid,k)
          w0(i,j)=sqrt(u0(i,j)*u0(i,j)+v0(i,j)*v0(i,j))
          h0(i,j)=h(igrid,k)

          u1(i,j)=u(igrid,k-1)
          v1(i,j)=v(igrid,k-1)
          w1(i,j)=sqrt(u1(i,j)*u1(i,j)+v1(i,j)*v1(i,j))
          h1(i,j)=h(igrid,k-1)

        end do
       end do


        do j=1,ny-1
         do i=1,nx-1

          igrid=i+(j-1)*nx

          dsh=(v0(i+1,j)-v0(i,j))*f1/dx   !dsh=dv/dx+du/dy
     +       +(u0(i,j+1)-u0(i,j))*f1/dy

          dst=(u0(i+1,j)-u0(i,j))*f1/dx    !dst=du/dx-dv/dy
     +       -(v0(i,j+1)-v0(i,j))*f1/dy

          def = sqrt(dsh*dsh + dst*dst)

          cvg =-((u0(i+1,j)-u0(i,j))*f1/dx    !cvg=-(du/dx+dv/dy)
     +          +(v0(i,j+1)-v0(i,j))*f1/dy)

          dh=h0(i,j)-h1(i,j)
          dh1=h0(i+1,j)-h1(i+1,j)
          dh2=h0(i-1,j)-h1(i-1,j)
          dh3=h0(i,j+1)-h1(i,j+1)
          dh4=h0(i,j-1)-h1(i,j-1)

          edge=dh*dh1*dh2*dh3*dh4

          if(edge.eq.0.0) then             !edge=0 means it is on domain edge
           vws = 0.0
          else
           vws = (w1(i,j)-w0(i,j))*f2/dh                        !vws=dw/dz
          end if

          turbindx=abs(vws)*(def+abs(cvg))
          if(turbindx.le.4.)  tb=0.
          if(turbindx.le.8.0.and.turbindx.gt.4.0) tb=1.
          if(turbindx.le.12.0.and.turbindx.gt.8.0) tb=2.
          if(turbindx.gt.12) tb=3.
          cat(igrid,k-1)=tb


c         if(igrid.ge.738400.and.igrid.le.738500) then
c          write(*,*) igrid,'level=',k-1,'ens=',ie,
c     +     'miss=',ms,'def=',def,'cvg=',cvg,'vws=',vws,
c     +      'cat=',cat(igrid,k-1)
c         end if

        end do
      end do

600   continue    !end of k

      end

