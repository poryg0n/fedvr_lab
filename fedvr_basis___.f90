      module fedvr_basis
      use fedvr
      use math_util, only :  differentiation_real, fourier_trsfm, inv_fourier_trsfm, fourier_diff
      implicit none
      contains
      



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

       pis = pis/sqrt(ww(i))

  
       return 
       end



       subroutine sector_mapper(np,ns,s,w,x,ww,rr)
       implicit none
       integer :: np, m, ns
       integer :: j, s, k
       real*8, dimension(*) :: x, w
       real*8, dimension(np) :: rr, ww
 
       write(13,*) "results from sector_mapper"
       do k=1,s
          do j=1,np-1
             if ((k.lt.ns)) then
                m = (k-1)*(np-1) + j
                rr(j)=x(m)
                ww(j)=w(m)
                write(13,*) m, ww(j), rr(j), ns*(np-1)+1, ns, np, s, k, j 
             else
                if ((j.lt.(np-1))) then
                   m = (k-1)*(np-1) + j
                   rr(j)=x(m)
                   ww(j)=w(m)
                   write(13,*) m, ww(j), rr(j), ns*(np-1)+1, ns, np, s, k, j 
                end if
             end if
          enddo
          write(13,*) 
       enddo
 
       return
       end

       end 
