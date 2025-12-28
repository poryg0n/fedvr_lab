     module math_util
     implicit none 
     contains

!      subroutine composite_simpson_18c(ndim, xx, df, res, ff)
!      implicit none
!      integer :: i, ndim, nn, nm, even
!      real*8 :: hh
!      complex*16 :: res, res1
!      real*8, dimension(*) :: xx
!      complex*16, dimension(*) :: df
!
!      real*8, dimension(4) :: xxl
!      complex*16, dimension(4) :: dfl
!      complex*16, optional, dimension(*) :: ff
!
!      nn=ndim-1
!      if (mod(nn,2).eq.0) then
!         even = 1 
!!        write(*,*) "N is even", nn
!      else
!!        write(*,*) "Number of subintervals is not even"
!         even = 0
!         nn=nn-3
!      end if
!
!      hh = (xx(2)-xx(1))
!!     write(*,*) hh, nn
!      
!      res = 0.d0
!      res1 = 0.d0
!      do i=1,nn/2
!         res = res +                                                  &
!                hh/3*( df(2*i-1) + 4.d0*df(2*i) + df(2*i+1) )
!!        write(*,*) nn/2, i, hh, res
!      enddo
!      if (even.eq.0) then
!!        write(*,*) "Not even even"
!         do i=1,4
!            xxl(i)=xx(ndim-4+i)
!            dfl(i)=df(ndim-4+i)
!         enddo
!         call composite_simpson_38c(4, xxl, dfl, res1) 
!      end if
!      res = res + res1
!!     write(*,*) res
!
!      if(present(ff)) then
!         ff(1)=df(1)
!         do i=1,ndim
!            ff(i+1) = ff(i) + .5d0*(df(i+1)+df(i))*(xx(i+1)-xx(i))
!         enddo
!      end if
!      end


!      subroutine composite_simpson_38c(ndim, xx, df, res, ff)
!      implicit none
!      integer :: i, ndim, mdim, nn, nm, third
!      real*8 :: hh
!      complex*16 :: res, res1
!      real*8, dimension(*) :: xx
!      complex*16, dimension(*) :: df
!
!      real*8, dimension(3) :: xxl
!      complex*16, dimension(3) :: dfl
!      real*8, optional, dimension(*) :: ff
!
!      third=0
!      nn=ndim-1
!      if (mod(nn,3).eq.0) then
!         third = 1
!      else
!!        write(*,*) "Number of subintervals is not multiple by 3"
!         nm=0
!         do while (nn.gt.3)
!            nn=nn-1
!            nm=nm+1
!!           write(*,*) nn, nm
!            if (mod(nn,3).eq.0) then            
!               if (mod(nm,2).eq.0) then
!                  mdim = nm + 1
!                  exit
!               end if
!            end if
!         enddo
!      end if
!
!      hh = (xx(2)-xx(1))
!!     write(*,*) hh, nn
!      
!      res = 0.d0
!      do i=1,nn/3
!         res = res +                                                  &
!                3.d0/8*hh*( df(3*i-2)                                 &
!                    + 3.d0*df(3*i-1) + 3.0d0*df(3*i)                  &
!                          +  df(3*i+1) )
!!        write(*,*) nn/3, i, hh, res
!      enddo
!      if (third.eq.0) then
!         do i=0,nm
!            xxl(i)=xx(ndim-nm+i)
!            dfl(i)=df(ndim-nm+i)
!         enddo
!         call composite_simpson_18c(mdim, xxl, dfl, res1) 
!      end if
!      res = res + res1
!!     write(*,*) res
!
!      end



      subroutine composite_trapeze_c(ndim, xx, df, res, ff)
      implicit none
      integer :: k, ndim, nn, even
      real*8 :: dx                                                    
      complex*16 :: res
      real*8, dimension(*) :: xx
      complex*16, dimension(*) :: df
      complex*16, optional, dimension(*) :: ff


      res = 0.d0
      nn = ndim-1
      do k=2,ndim
         dx   = xx(k)-xx(k-1)
         res = res + 0.5d0* dx * ( (df(k-1) + df(k)))
!        write(101,*) nn+1, k, dx, xx(k), df(k), res
      enddo

      if(present(ff)) then
         ff(1)=df(1)
         do k=1,ndim
            ff(k+1) = ff(k) + .5d0*(df(k+1)+df(k))*(xx(k+1)-xx(k))
         enddo
      end if
      end


!      subroutine composite_simpson_18(ndim, xx, df, res, ff)
!      implicit none
!      integer :: i, ndim, nn, nm, even
!      real*8 :: hh, res, res1
!      real*8, dimension(*) :: xx, df
!      real*8, dimension(4) :: xxl, dfl
!      real*8, optional, dimension(*) :: ff
!
!
!      nn=ndim-1
!      if (mod(nn,2).eq.0) then
!         even = 1 
!!        write(*,*) "N is even", nn
!      else
!!        write(*,*) "Number of subintervals is not even"
!         even = 0
!         nn=nn-3
!      end if
!
!      hh = (xx(2)-xx(1))
!!     write(*,*) hh, nn
!      
!      res = 0.d0
!      res1 = 0.d0
!      do i=1,nn/2
!         res = res +                                                  &
!                hh/3*( df(2*i-1) + 4.d0*df(2*i) + df(2*i+1) )
!!        write(*,*) nn/2, i, hh, res
!      enddo
!      if (even.eq.0) then
!!        write(*,*) "Not even even"
!         do i=1,4
!            xxl(i)=xx(ndim-4+i)
!            dfl(i)=df(ndim-4+i)
!         enddo
!         call composite_simpson_38(4, xxl, dfl, res1) 
!      end if
!      res = res + res1
!!     write(*,*) res
!
!      if(present(ff)) then
!         ff(1)=df(1)
!         do i=1,ndim
!            ff(i+1) = ff(i) + .5d0*(df(i+1)+df(i))*(xx(i+1)-xx(i))
!         enddo
!      end if
!      end



!      subroutine composite_simpson_38(ndim, xx, df, res, ff)
!      implicit none
!      integer :: i, ndim, mdim, nn, nm, third
!      real*8 :: hh, res, res1
!      real*8, dimension(*) :: xx, df
!      real*8, dimension(4) :: xxl, dfl
!      real*8, optional, dimension(*) :: ff
!
!
!      third=0
!      nn=ndim-1
!      if (mod(nn,3).eq.0) then
!         third = 1
!      else
!!        write(*,*) "Number of subintervals is not multiple by 3"
!         nm=0
!         do while (nn.gt.3)
!            nn=nn-1
!            nm=nm+1
!!           write(*,*) nn, nm
!            if (mod(nn,3).eq.0) then            
!               if (mod(nm,2).eq.0) then
!                  mdim = nm + 1
!                  exit
!               end if
!            end if
!         enddo
!      end if
!
!      hh = (xx(2)-xx(1))
!!     write(*,*) hh, nn
!      
!      res = 0.d0
!      do i=1,nn/3
!         res = res +                                                  &
!                3.d0/8*hh*( df(3*i-2)                                 &
!                    + 3.d0*df(3*i-1) + 3.0d0*df(3*i)                  &
!                          +  df(3*i+1) )
!!        write(*,*) nn/3, i, hh, res
!      enddo
!      if (third.eq.0) then
!         do i=0,nm
!            xxl(i)=xx(ndim-nm+i)
!            dfl(i)=df(ndim-nm+i)
!         enddo
!         call composite_simpson_18(mdim, xxl, dfl, res1) 
!      end if
!      res = res + res1
!!     write(*,*) res
!
!      end


!      subroutine composite_simpson(ndim, xx, df, res, ff)
!      implicit none
!      integer :: i, ndim, nn, even
!      real*8 :: h2i, h2ip1, hnm1, hnm2,                               &
!             alpha, beta, eta, res
!      real*8, dimension(*) :: xx, df
!      real*8, optional, dimension(*) :: ff
!
!      even=0
!      if (mod(ndim,2).eq.0) then
!!        write(*,*) "N is even", ndim
!         even = 1
!         nn   = ndim-1
!      else
!!        write(*,*) "N is odd", ndim
!         nn   = ndim-2
!      end if
!      
!      res = 0.d0
!      do i=0,nn/2-1
!         h2i   = xx(2*i+2)-xx(2*i+1) 
!         h2ip1 = xx(2*i+3)-xx(2*i+2) 
!         res = res + (h2i + h2ip1)/6 *       (                          &
!                 (2.d0 - h2ip1/h2i) * df(2*i+1)                         &
!              +  ((h2i + h2ip1)**2/(h2i*h2ip1)) * df(2*i+2)             &
!              +  (2.d0 - h2i/h2ip1) * df(2*i+3) )
!!        write(*,*) nn/2,2*i+1, 2*i, i, xx(2*i+2)-xx(2*i+1), res
!      enddo
!
!      if (even.eq.0) then
!!        write(*,*) i
!         hnm1  = xx(2*i+2)-xx(2*i+1) 
!         hnm2  = xx(2*i+1)-xx(2*i) 
!
!         alpha =  (2.d0*hnm1**2 + 3.0d0*hnm1*hnm2)/(6.d0*(hnm1+hnm2))
!         beta  =  (hnm1**2 + 3.0d0*hnm1*hnm2)/(6.d0*hnm2)
!         eta   =  (hnm1**3)/(6.d0*hnm2*(hnm1+hnm2))
!         res = res +                                                  &
!                 alpha*df(2*i+2) + beta * df(2*i+1) -eta*df(2*i)
!      end if
!
!      if(present(ff)) then
!         ff(1)=df(1)
!         do i=1,ndim
!            ff(i+1) = ff(i) + .5d0*(df(i+1)+df(i))*(xx(i+1)-xx(i))
!         enddo
!      end if
!      end


      subroutine differentiation(ndim, xx, ff, dx, df)
      implicit none
      integer :: i, ndim
      real*8, dimension(*) :: xx, dx
      complex*16, dimension(*) :: ff, df
      complex*16, dimension(ndim) :: dx1, df1, dx2, df2

      dx1 = 0.d0
      dx2 = 0.d0
      df1 = 0.d0
      df2 = 0.d0

!      do i=1,ndim-1
!         dx1(i)=xx(i+1)-xx(i)
!!        df1(i)=(ff(i+1)-ff(i)) / dx(i)
!         df1(i)=(ff(i+1)-ff(i)) 
!      enddo
!
!      do i=2,ndim
!         dx2(i)=xx(i)-xx(i-1)
!!        df2(i)=(ff(i)-ff(i-1)) / dx(i)  
!         df2(i)=(ff(i)-ff(i-1))
!      enddo
!
!      do i=1,ndim
!         dx(i)=.5d0* (dx2(i)+dx1(i))
!         df(i)=.5d0* (df2(i)+df1(i)) / dx(i)
!      enddo

      if (ndim.le.1) then
         dx(1)=xx(1)
         df(1)=ff(1)
         return
      else

         do i=1,ndim-1
            dx1(i)=xx(i+1)-xx(i)
            df1(i)=(ff(i+1)-ff(i))/dx1(i)
         enddo

         do i=2,ndim
            dx2(i)=xx(i)-xx(i-1)
            df2(i)=(ff(i)-ff(i-1))/dx2(i)
         enddo

!        do i=1,ndim
!           dx(i)=.5d0* (dx2(i)+dx1(i))
!           df(i)=.5d0* (df2(i)+df1(i)) 
!        enddo

         do i=1,ndim
            if (i.eq.1) then
               dx(i)=dx1(i)
               df(i)=df1(i)
            else if (i.eq.ndim) then
               dx(i)=dx2(i)
               df(i)=df2(i)
            else
               dx(i)=.5d0*(dx2(i)+dx1(i))
               df(i)=.5d0*(df2(i)+df1(i))
            end if
         enddo
      end if

      end

      subroutine differentiation_real(ndim, xx, ff, dx, df)
      implicit none
      integer :: i, ndim
      real*8, dimension(*) :: xx, dx
      real*8, dimension(*) :: ff, df
      real*8, dimension(ndim) :: dx1, df1, dx2, df2

      dx1 = 0.d0
      dx2 = 0.d0
      df1 = 0.d0
      df2 = 0.d0

      if (ndim.le.1) then
         dx(1)=xx(1)
         df(1)=ff(1)
         return
      else

         do i=1,ndim-1
            dx1(i)=xx(i+1)-xx(i)
            df1(i)=(ff(i+1)-ff(i))/dx1(i)
         enddo

         do i=2,ndim
            dx2(i)=xx(i)-xx(i-1)
            df2(i)=(ff(i)-ff(i-1))/dx2(i)
         enddo

         do i=1,ndim
            if (i.eq.1) then
               dx(i)=dx1(i)
               df(i)=df1(i)
            else if (i.eq.ndim) then
               dx(i)=dx2(i)
               df(i)=df2(i)
            else
               dx(i)=.5d0*(dx2(i)+dx1(i))
               df(i)=.5d0*(df2(i)+df1(i))
            end if
         enddo
      end if

      return
      end


      subroutine inv_fourier_trsfm(ndim, krange, xx, kk, ff, ffk)
      implicit none
      integer :: i, k, ndim, krange
      real*8 :: ppi
      complex*16 :: ci, res
      real*8, dimension(ndim) :: xx
      real*8, dimension(krange) :: kk
      complex*16, dimension(ndim) :: ff, aux
      complex*16, dimension(krange) :: ffk

      ppi = 4.d0*datan(1.d0)
      ci = (0.d0,1.d0)

      res = ppi
      do k=1,krange
         do i=1,ndim
            aux(i) = ff(i) * exp(-ci*xx(i)*kk(k))
         enddo
         call composite_trapeze_c(ndim, xx, aux, res)
         ffk(k) = res
      enddo
      ffk=ffk/(2.d0*ppi)

      return
      end



      subroutine fourier_trsfm(krange, ndim, kk, xx, ffk, ff)
      implicit none
      integer :: i, k, ndim, krange
      real*8 :: ppi
      complex*16 :: ci, res
      real*8, dimension(ndim) :: xx
      real*8, dimension(krange) :: kk
      complex*16, dimension(ndim) :: ff
      complex*16, dimension(krange) :: ffk, aux

      ppi = 4.d0*datan(1.d0)
      ci = (0.d0,1.d0)

      res = ppi
      do i=1,ndim
         do k=1,krange
            aux(k) = ffk(k) * exp(ci*xx(i)*kk(k))
         enddo
         call composite_trapeze_c(krange, kk, aux, res)
         ff(i) = res
      enddo

      return
      end




      subroutine fourier_diff(ndim, krange, kk, xx, ff, dx, df)
      implicit none
      integer :: i, ndim, krange, k
      real*8, dimension(*) :: xx, dx, kk
      complex*16 :: ci
      complex*16, dimension(*) :: ff, df
      complex*16, dimension(ndim) :: dx1, df1, dx2, df2
      complex*16, dimension(krange) :: dkk, ffk, dfk, aux

      ci = (0.d0,1.d0)
      call inv_fourier_trsfm(ndim, krange, xx, kk, ff, ffk)
      do k=1,krange
         aux(k) = ci*kk(k)*ffk(k)
      enddo
      call fourier_trsfm(krange, ndim, kk, xx, aux, df)


      return
      end




      real*8 function sgn(x)
      implicit none
      real*8 :: x

      if(x.gt.0d0) then
         sgn = 1.d0
      else if(x.lt.0d0) then
         sgn = - 1.d0
      else 
         sgn = 0.d0
      end if

      return
      end



      real*8 function pis(np,ns,s,i,w,x,r)
      implicit none
      integer :: np, ns, s, i
      integer :: j, k
      real*8 :: r, eps
      real*8, dimension(*) :: w, x
      real*8, dimension(np) :: rr, ww
      real*8, dimension(ns) :: rs

      call sector_mapper(np,ns,s,w,x,ww,rr)
      eps= 1.d-20
      pis= 1.d0

      if (s.ne.ns) then
         do j=1,np-1
            if (i/=j) then
               pis = pis* (r-rr(j))/(rr(i)-rr(j))
            end if
         enddo
      else
         do j=1,np-2
            if (i/=j) then
               pis = pis* (r-rr(j))/(rr(i)-rr(j))
            end if
         enddo
      end if

!     pis = pis/sqrt(ww(i))

  
      return 
      end


      real*8 function g_pis(np,ns,n,w,x,r)
      implicit none
      integer :: np,ns, n, m, mp
      integer :: i, j, k, s, aux1, sp
      real*8 :: aux0, xmin
      real*8, dimension(*) :: w, x
      real*8, dimension(np) :: rr, ww
      real*8, dimension(ns) :: rs

      real*8 :: r
!     real*8 :: pis


      do sp=1,ns
        if (n.le.((sp-1)*(np-1)+np-1)) then
           i = n - (sp-1)*(np-1)
           exit
        end if
      enddo

      s=sp
      call sector_mapper(np,ns,s,w,x,ww,rr)  

      if ((r.gt.x((s-1)*(np-1)+1)).and.((r.lt.x((s)*(np-1)+np-2)))) then
         if (n.eq.((s-1)*(np-1)+np-1)) then
            if (r.lt.((ww(1)*x(n)+ww(np-1)*x(n+1))/(ww(np-1)+ww(1)))) then
               g_pis = sqrt(ww(np-1)) * pis(np,ns,s,np-1,w,x,r)
            else 
               g_pis = sqrt(ww(1)) * pis(np,ns,s+1,1,w,x,r)
            end if
            g_pis = g_pis / sqrt(ww(1)+ww(np-1))
         else
            g_pis = pis(np,ns,s,i,w,x,r)
         end if
      else
         g_pis = 0.d0
         write(29,*) x((s-1)*(np-1)+1), x((s)*(np-1)+np-2)
      end if


      return 
      end



!     real*8 function dpis(np,ns,s,i,w,x,r)
!     implicit none
!     integer :: np, ns, s, i
!     integer :: j, k, m
!     real*8 :: aux, r
!     real*8, dimension(*) :: w, x
!     real*8, dimension(np) :: rr, ww

!     real*8 :: pis
! 
!     dpis= 0.d0
!     call sector_mapper(np,ns,s,w,x,ww,rr)


!     if (s.ne.ns) then
!        do k=1,np-1
!           aux = 1.d0
!           if ((k/=i)) then
!              do j=1,np-1
!                 if ((j/=i).and.(j/=k)) then
!                    aux = aux* (r-rr(j))/(rr(i)-rr(j))
!                 end if
!              enddo
!              dpis = dpis + aux / (rr(i)-rr(k))
!           end if
!        enddo

!     else
!        do k=1,np-2
!           aux = 1.d0
!           if ((k/=i)) then
!              do j=1,np-2
!                 if ((j/=i).and.(j/=k)) then
!                    aux = aux* (r-rr(j))/(rr(i)-rr(j))
!                 end if
!              enddo
!              dpis = dpis + aux / (rr(i)-rr(k))
!           end if
!        enddo

!     end if
!
!     dpis = dpis/sqrt(ww(i))
!
!     return         
!     end



!     real*8 function dg_pis(np,ns,n,w,x,r)
!     implicit none
!     integer :: np, ns, n, s, i
!     integer :: j, k, m, sp, sm
!     real*8 :: aux, r, aux0, aux1
!     real*8, dimension(*) :: w, x
!     real*8, dimension(np) :: rr, ww
!     real*8, dimension(np) :: rrr, www

!     real*8 :: dpis
!     real*8 :: pis
! 
!     dg_pis = 0.d0

!      do sp=1,ns
!        if (n.le.((sp-1)*(np-1)+np-1)) then
!           i = n - (sp-1)*(np-1)
!           exit
!        end if
!      enddo

!      s=sp

!      if (n.eq.((s-1)*(np-1)+np-1)) then
!         call sector_mapper(np,ns,s,w,x,ww,rr)    ! the ww from this
!
!         if (r.lt.((ww(1)*x(n)+ww(np-1)*x(n+1))/(ww(np-1)+ww(1)))) then
!            dg_pis = sqrt(ww(np-1))* dpis(np,ns,s,np-1,w,x,r)
!         else
!            dg_pis = sqrt(ww(1)) * dpis(np,ns,s+1,1,w,x,r)
!         end if

!         dg_pis = dg_pis /(sqrt(ww(1)+ww(np-1)))
!
!      else
!
!         dg_pis = dpis(np,ns,s,i,w,x,r)
!
!      end if


!     return         
!     end



!     real*8 function dpis2(np,ns,s,i,w,x,j)
!     implicit none
!     integer :: np, ns, s, i
!     integer :: j, k, m
!     real*8 :: aux, r
!     real*8, dimension(*) :: w, x
!     real*8, dimension(np) :: rr, ww

!     real*8 :: pis
! 
!     dpis2= 0.d0
!     call sector_mapper(np,ns,s,w,x,ww,rr)

!     if ((j/=i)) then
!        aux = 1.d0
!        do k=1,np-1
!           if ((k/=i).and.(k/=j)) then
!              aux = aux* (rr(j)-rr(k))/(rr(i)-rr(k))
!           end if
!        enddo
!        dpis2 =  aux / (rr(i)-rr(j))
!     else 
!        dpis2 =  0.5d0*((pis(np,ns,s,i,w,x,rr(np-1)))**2 -           &
!                       (pis(np,ns,s,i,w,x,rr(1)))**2)
!     end if

!     dpis2 = dpis2 / sqrt(ww(i))
!
!     return         
!     end



      subroutine sector_mapper(np,ns,s,w,x,ww,rr)
      implicit none
      integer :: np, m, ns
      integer :: j, s, k
      real*8, dimension(*) :: x, w
      real*8, dimension(np) :: rr, ww


      do k=1,s
         do j=1,np-1
            if ((k.lt.ns)) then
               m = (k-1)*(np-1) + j
               rr(j)=x(m)
               ww(j)=w(m)
            else
               if ((j.lt.(np-1))) then
                  m = (k-1)*(np-1) + j
                  rr(j)=x(m)
                  ww(j)=w(m)
               end if
            end if
         enddo
      enddo

      return
      end



      end



