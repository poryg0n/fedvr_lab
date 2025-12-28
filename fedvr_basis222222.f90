      module fedvr_basis
      use lobatto
      use math_util, only :  differentiation_real, fourier_trsfm,     &
                            inv_fourier_trsfm, fourier_diff
      implicit none
      contains
      

        subroutine basis_1d(basis,lnbr,nnbr,xa,wa)
        use util
        integer :: i, j, s, k, nnbr, lnbr, n, np, nmax, nn,        &
                         sp
        real*8 :: wa(*), xa(*) 
        real*8 :: basis(lnbr*(nnbr-1)-1,lnbr*(nnbr-1)-1) 

        do s=1,lnbr
           sp=s
           do n=1,nnbr-2
              do np=1,nnbr-2
                 basis((s-1)*(nnbr-1) + n, (sp-1)*(nnbr-1)+np) =      &
                    pis_ij(nnbr,lnbr,s,n+1,wa,xa,np+1)
              enddo 
           enddo

           n=nnbr-1
           do np=1,nnbr-2
             if (((nnbr-1)*(s-1)+ n < (nnbr-1)*lnbr)                 &
                  .and.((nnbr-1)*(sp-1)+np < (nnbr*lnbr))) then
                   basis((s-1)*(nnbr-1) + n, (sp-1)*(nnbr-1)+np) =   &
                  dsqrt(wa(nnbr)/(wa(1)+wa(nnbr))) *                 &
                       pis_ij(nnbr,lnbr,s,n+1,wa,xa,np+1)
              end if
           enddo


           np=nnbr-1
           do n=1,nnbr-2
             if (((nnbr-1)*(s-1)+ n < (nnbr-1)*lnbr)                  &
                .and.((nnbr-1)*(sp-1)+np < ((nnbr-1)*lnbr))) then
                basis((s-1)*(nnbr-1) + n, (sp-1)*(nnbr-1)+np) =       &
                dsqrt(wa(nnbr)/(wa(1)+wa(nnbr)))      *               &
                     pis_ij(nnbr,lnbr,s,n+1,wa,xa,np+1)
              end if
           enddo

           n  = nnbr - 1
           np = nnbr - 1
           if (((nnbr-1)*(s-1) + n < (nnbr-1)*lnbr)                  &
                .and.(((nnbr-1)*(sp-1)+np) < (nnbr*lnbr))) then
                   basis((s-1)*(nnbr-1) + n, (sp-1)*(nnbr-1)+np) =   &
                       (wa(nnbr)/(wa(1)+wa(nnbr))) *                 &
                       pis_ij(nnbr,lnbr,s,nnbr,wa,xa,nnbr)           &
               + (wa(1)/(wa(1)+wa(nnbr))) *                          &
                       pis_ij(nnbr,lnbr,s,1,wa,xa,1)
           end if 

        enddo


!       do s=2,lnbr
!          sp=s-1
!          np=nnbr-1
!          do n=1,nnbr-1
!             if(((nnbr-1)*(s-1)+n < (nnbr-1)*lnbr)                  &
!               .and.((nnbr-1)*(sp-1)+np < nnbr*lnbr)) then
!                  basis((s-1)*(nnbr-1) + n, (sp-1)*(nnbr-1)+np) =   &
!                 dsqrt(wa(1)/(wa(1)+wa(nnbr))) *                    &
!                    pis_ij(nnbr,lnbr,s,n+1,wa,xa,1)
!             end if
!          enddo

!          n=nnbr-1
!          if(((nnbr-1)*(s-1)+n < (nnbr-1)*lnbr)                      &
!            .and.((nnbr-1)*(sp-1)+np < nnbr*lnbr)) then
!               basis((s-1)*(nnbr-1) + n, (sp-1)*(nnbr-1)+np) =       &
!              dsqrt(wa(1)*wa(nnbr))/(wa(1)+wa(nnbr)) *               &
!                 pis_ij(nnbr,lnbr,s,1,wa,xa,nnbr)
!          end if

!       enddo


!       do s=1,lnbr-1
!          sp=s+1
!          n=nnbr-1
!          do np=1,nnbr-1
!             if(((nnbr-1)*(s-1)+n < (nnbr-1)*lnbr)                  &
!               .and.((nnbr-1)*(sp-1)+np < nnbr*lnbr)) then
!                  basis((s-1)*(nnbr-1) + n, (sp-1)*(nnbr-1)+np) =   &
!                 dsqrt(wa(1)/(wa(1)+wa(nnbr))) *                    &
!                    pis_ij(nnbr,lnbr,s,1,wa,xa,np+1)
!             end if
!          enddo

!          np=nnbr-1
!          if(((nnbr-1)*(s-1)+n < (nnbr-1)*lnbr)                      &
!            .and.((nnbr-1)*(sp-1)+np < (nnbr-1)*lnbr)) then
!               basis((s-1)*(nnbr-1) + n, (sp-1)*(nnbr-1)+np) =       &
!              dsqrt(wa(nnbr)*wa(1))/(wa(1)+wa(nnbr)) *               &
!                 pis_ij(nnbr,lnbr,s,1,wa,xa,nnbr)
!          end if

!       enddo


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
 
 
       do j=1,np
          if (i/=j) then
             pis = pis* (r-x(j))/(x(i)-x(j))
          end if
       enddo

       pis = pis / dsqrt(w(i))
 
 
       return 
       end


       real*8 function pis_ij(np,ns,s,i,w,x,j)
       implicit none
       integer :: np, ns, s, i, j
       integer :: k
       real*8 :: eps
       real*8, dimension(*) :: w, x
       real*8, dimension(np) :: rr, ww
       real*8, dimension(ns) :: rs
 
 
 
       if (i.eq.j) then
          pis_ij = 1.d0
       else
          pis_ij = 0.d0
       end if
     
       pis_ij = pis_ij/dsqrt(w(i))
 
       return 
       end


        real*8 function gpis(np,ns,n,w,x,r)
        implicit none
        integer :: np,ns, n, m, mp
        integer :: i, j, k, s, aux1, sp, sm
        real*8 :: aux0, xmin
        real*8, dimension(*) :: w, x
        real*8, dimension(np) :: rr, ww
        real*8, dimension(ns) :: rs
  
        real*8 :: r
  
        do sp=1,ns
          if (n.le.((sp-1)*(np-1)+np)) then
             i = n - (sp-1)*(np-1)
             exit
          end if
        enddo
    
        s=sp
        if (n.eq.((s-1)*(np-1)+np)) then
  
           if (r.lt.x(n)) then
              gpis = dsqrt(w(np)) * pis(np,ns,s,np,w,x,r)
           else
              gpis = dsqrt(w(1))  * pis(np,ns,s+1,1,w,x,r)
           end if

           gpis = gpis / dsqrt(w(1)+w(np))
        else
           gpis = pis(np,ns,s,i,w,x,r)
        end if

        return 
        end


        real*8 function gpis_ij(np,ns,n,w,x,j)
        implicit none
        integer :: np,ns, n, m, mp
        integer :: i, j, k, s, aux1, sp, sm
        real*8 :: aux0, xmin
        real*8, dimension(*) :: w, x
  
  
        do sp=1,ns
          if (n.le.((sp-1)*(np-1)+np)) then
             i = n - (sp-1)*(np-1)
             exit
          end if
        enddo
    
        s=sp
        if (n.eq.((s-1)*(np-1)+np)) then
  
           if (x(j).lt.x(n)) then
              gpis_ij = dsqrt(w(np)) * pis_ij(np,ns,s,np,w,x,j)
           else
              gpis_ij = dsqrt(w(1)) *  pis_ij(np,ns,s,1,w,x,j)
           end if

           gpis_ij = gpis_ij / dsqrt(w(1)+w(np))
        else
           gpis_ij = pis_ij(np,ns,s,i,w,x,j)
        end if

        return 
        end
  
 
 
 
        real*8 function dpis(np,ns,s,i,w,x,r)
        implicit none
        integer :: np, ns, s, i
        integer :: j, k, m
        real*8 :: aux, r
        real*8, dimension(*) :: w, x
        real*8, dimension(np) :: rr, ww
 
 
        do j=1,np
           aux = 1.d0
           if ((j/=i)) then
              do k=1,np
                 if ((k/=i).and.(k/=j)) then
                    aux = aux* (r-x(k))/(x(i)-x(k))
                 end if
              enddo
              dpis = dpis + aux / (x(i)-x(j))
           end if
        enddo
  
        dpis = dpis/sqrt(w(i))
  
        return         
        end
 
 
        real*8 function dg_pis(np,ns,n,w,x,r)
        implicit none
        integer :: np, ns, n, s, i
        integer :: j, k, m, sp, sm
        real*8 :: aux, r, aux0, aux1
        real*8, dimension(*) :: w, x
  
    
         do s=1,ns
           if (n.le.((s-1)*(np-1)+np)) then
              i = n - (s-1)*(np-1)
              exit
           end if
         enddo
  
  
         if (n.eq.((s-1)*(np-1)+np)) then 
            if (r.gt.x(i)) then 
               dg_pis = dsqrt(w(np)) * dpis(np,ns,s,np,w,x,r)
            else  
               dg_pis = dsqrt(w(1)) * dpis(np,ns,s,np,w,x,r)
            end if
               dg_pis = dg_pis / sqrt(w(1)+w(np))
         else 
            dg_pis = dpis(np,ns,s,i,w,x,r)
         end if
  
        return         
        end
  
  
  

        real*8 function dpis_ij(np,ns,s,i,w,x,j)
        implicit none
        integer :: np, ns, s, i
        integer :: j, k, m
        real*8 :: aux, r
        real*8, dimension(*) :: w, x
        real*8, dimension(np) :: rr, ww
 
!       real*8 :: pis

 
        if ((j/=i)) then
           aux = 1.d0
           do k=1,np
              if ((k/=i).and.(k/=j)) then
                 aux = aux* (x(j)-x(k))/(x(i)-x(k))
              end if
           enddo
           dpis_ij =  aux / (x(i)-x(j))
        else
           dpis_ij = 0.5d0 *                                            &
             (pis_ij(np,ns,s,i,w,x,np)**2 - pis_ij(np,ns,s,i,w,x,1)**2)
        end if
 
        dpis_ij = dpis_ij / sqrt(w(i))
  
        return         
        end




!        real*8 function dgpis_ref(np,ns,n,w,x,sp,j)
!        implicit none
!        integer :: np, ns, sp, s, i
!        integer :: j, k, n
!        real*8 :: aux, r
!        real*8, dimension(*) :: w, x
!        real*8, dimension(np) :: rr, ww
! 
!        dgpis_ref = 0.0d0

!        do j=1,np
!           if (j.ne.n)
!              dgpis_ref = dgpis_ref  +                               &
!                 .5d0 * w(j) *                                       &
!                  dpis_ij(np,ns,s,i,w,x,j)
!           end if
!        enddo
!        return         
!        end
 
 



      end module
