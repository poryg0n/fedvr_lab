      module fedvr_basis
      use lobatto
      use math_util, only :  differentiation_real, fourier_trsfm,     &
                            inv_fourier_trsfm, fourier_diff
      implicit none
      contains
      

        subroutine analytical_basis_vectors(workdir, s, i,            &
                                             xa, xs, xt, xx, xy,      &
                                             wa, wx, wy, wz,          &
                                             auxf, auxd, wfc, dwfc)
        use util
        include 'open/files_index'
        integer :: aux0, krange
        integer :: i, j, s, k, np, ns, n, m, mp, nmax1, nmax2, nn,    &
                         sp, npp, ip
        real*8 :: r, step, xx1, xx2, scaling,                         &
                             step_, kstep, y0, kk1, kk2
        real*8, dimension(:), allocatable :: xa, xt, xx, xy, xs,  &
                                     wa, wx, wy, wz

        real*8, dimension(:), allocatable :: aux1, aux2, aux3, aux4
        real*8, dimension(:) :: auxf, auxd

        real*8, dimension(:), allocatable :: ppis, ppis_,             &
                                    gppis, gppis_,                    &
                                    ddppis, ddppis_,                  &
                                    ff, df, dff, dx

        complex*16, dimension(:) :: wfc, dwfc

        complex*16, dimension(:), allocatable :: ffc, ffk, ffc1,      &
                              dfc, wwfc, wwfc_

        character(len=1024) :: filename01='quad_pts_wts'
        character(len=1024) :: filename02='ppi_'
        character(len=1024) :: filename03='ppic_'
        character(len=1024) :: filename04='dpi_'
        character(len=1024) :: filename05='dpic_'

        character(len=1024) :: filename1
        character(len=1024) :: filename2
        character(len=1024) :: filename3
        character(len=1024) :: filename4
        character(len=1024) :: filename5

        character(len=1024) :: filename1_
        character(len=1024) :: filename2_
        character(len=1024) :: filename3_
        character(len=1024) :: filename4_
        character(len=1024) :: filename5_


        character(len=1024) :: format_string
        character(len=1024) :: indx, indx1, indx2 
        character(255) :: workdir


!       real*8 :: pis, dpis, dpis2, g_pis, dg_pis
 
        include 'param'
        include 'common/open_files'
        include 'open/files_specs'

!       include 'common/laser'
        include 'open/files_deriv_analyt'

!       call open_files_to_store(0)

!       nbpts=200
 
!       ut=100
!       ut=200
!       ut=300
!       ut=400
 
        nmax1 = lnbr*(nnbr-1) - 1
        nmax2 = lnbr*(nnbr-1) + 1
        xx1 = -xmax/2.d0
        xx2 =  xmax/2.d0
 
        scaling = 0.5d0*xmax/lnbr
 
 
       allocate(ppis(nbpts), ppis_(nbpts))
       allocate(gppis(nbpts), gppis_(nbpts))
       allocate(wwfc(nbpts), wwfc_(nbpts))

       allocate(aux1(nmax1), aux2(nmax1))

!      allocate(aux3(nbpts))
       step = (xx2 - xx1) / (nbpts-1)

       allocate(ff(nbpts), df(nbpts), dff(nbpts),                   &
                                        ddppis(nbpts), dx(nbpts))

       write(ut80,*) "phase 1 : Find the choke points", lnbr*(nnbr-1)-1
       npp = nnbr-1                           ! for effective version

       m = (s-1)*(nnbr-1) + i

       if (i.ne.nnbr) then
          if ((i.ne.(nnbr-1)).or.(s.ne.lnbr)) then     
             do sp=1,lnbr
                do j=1,nnbr-1
                   if((sp.ne.lnbr).or.(j.ne.(nnbr-1))) then
                      mp = (sp-1)*(nnbr-1)+j
                      write(ut80,*)  mp+1, sp, j+1,                    & 
                          dsqrt(wa(j+1)), xa(j), xx(mp),               &
                          pis(nnbr,lnbr,sp,j,xx,xy,wa,xx(mp)),         &
                          pis_ij(wa,j+1,j+1),                          &
                          g_pis(nnbr,lnbr,mp,xx,xy,wa,xx(mp)),         &
                          dpis_ij(nnbr,lnbr,sp,j+1,xx,xy,wa,j+1),      &
                          dpis(nnbr,lnbr,sp,j,xx,xy,wa,xx(mp)),        &
                          dg_pis(nnbr,lnbr,mp,xx,xy,wa,xx(mp)),        &
                          0.d0


                      aux1(mp) = g_pis(nnbr,lnbr,m,xx,xy,wa,xx(mp))
                      aux2(mp) = dg_pis(nnbr,lnbr,m,xx,xy,wa,xx(mp))
                      write(ut81,*)   m, sp, i+1, j+1, xx(mp),           &
                          pis(nnbr,lnbr,sp,i,xx,xy,wa,xx(mp)),           &
                          pis_ij(wy,i+1,j+1),                            &
                          g_pis(nnbr,lnbr,m,xx,xy,wa,xx(mp)),            &
                          dpis_ij(nnbr,lnbr,s,i+1,xx,xy,wa,j+1),         &
                          dpis(nnbr,lnbr,s,i,xx,xy,wa,xx(mp)),           &
                          dg_pis(nnbr,lnbr,m,xx,xy,wa,xx(mp)),            &
                          0.d0
                   end if
                enddo
                write(ut80,*)
                write(ut81,*)
             enddo


! ***        continuous \pi and \Pi
             do k=1,nbpts
                xt(k) = xx1 + (k-1)*step
                ppis(k)    = pis(nnbr,lnbr,s,i,xx,xy,wa,xt(k))
                gppis(k)   = g_pis(nnbr,lnbr,m,xx,xy,wa,xt(k))
                write(ut82,*) k, lnbr, nnbr, s, i+1, m, xt(k),        &
                   pis(nnbr,lnbr,s,i,xx,xy,wa,xt(k)),                 &
                   g_pis(nnbr,lnbr,m,xx,xy,wa,xt(k)),                 &
                   dpis(nnbr,lnbr,s,i,xx,xy,wa,xt(k)),                &
                   dg_pis(nnbr,lnbr,m,xx,xy,wa,xt(k))
             enddo
          end if

          auxf = auxf + aux1
          auxd = auxd + aux2
          do k=1,nmax1
             write(68,*) s, i,  k, wfc(m), wfc(m) * wx(m) * auxf(k)
             write(69,*) s, i,  k, dwfc(m), dwfc(m) * wx(m) * auxd(k)
          enddo
          write(68,*)
          write(69,*)
       end if


       write(ut90,*) "phase 1 : Find the choke points", lnbr*(nnbr-1)+1
       do sp=1,lnbr
         do j=1,nnbr
            mp = (sp-1)*(nnbr-1)+j
            write(ut90,*)   mp, sp, j,                               &
                dsqrt(wa(j)), xa(j), xy(mp),                         &
                pis_(nnbr,lnbr,sp,j,xx,xy,wa,xy(mp)),                &
                pis_ij(wa,j,j),                                      &
                g_pis_(nnbr,lnbr,mp,xx,xy,wa,xy(mp)),                &
                dpis_ij(nnbr,lnbr,sp,j,xx,xy,wa,j),                  &
                dpis_(nnbr,lnbr,sp,j,xx,xy,wa,xy(mp)),               &
                dg_pis_(nnbr,lnbr,mp,xx,xy,wa,xy(mp)),               &
                0.d0

            write(ut91,*)   m, sp, i, j, xy(mp),                     &
                pis_ij(wa,i,j),                                      &
                pis_(nnbr,lnbr,s,i,xx,xy,wa,xy(mp)),                 &
                g_pis_(nnbr,lnbr,m,xx,xy,wa,xy(mp)),                 &
                dpis_ij(nnbr,lnbr,s,i,xx,xy,wa,j),                   &
                dpis_(nnbr,lnbr,s,i,xx,xy,wa,xy(mp)),                &
                dg_pis_(nnbr,lnbr,m,xx,xy,wa,xy(mp)),                 &
                0.d0
          enddo
          write(ut90,*)
          write(ut91,*)
       enddo

! *** continuous \pi and \Pi
      do k=1,nbpts
         xt(k) = xx1 + (k-1)*step
         ppis_(k)    = pis_(nnbr,lnbr,s,i,xx,xy,wa,xt(k))
         gppis_(k)   = g_pis_(nnbr,lnbr,m,xx,xy,wa,xt(k))
         write(ut92,*) k, lnbr, nnbr, s, i, m,  xt(k),                &
                   pis_(nnbr,lnbr,s,i,xx,xy,wa,xt(k)),                &
                   g_pis_(nnbr,lnbr,m,xx,xy,wa,xt(k)),                &
                   dpis_(nnbr,lnbr,s,i,xx,xy,wa,xt(k)),               &
                   dg_pis_(nnbr,lnbr,m,xx,xy,wa,xt(k))
      enddo


! **  This is the part where we assess the derivative


!!!!! ***  Derivative at the choke points
!!!!       if (m.ne.(lnbr*(nnbr-1)+1)) then
!!!!          do sp=1,lnbr
!!!!             do j=1,nnbr
!!!!                 mp = (sp-1)*(nnbr-1)+j
!!!!                 auxf(mp)= g_pis(nnbr,lnbr,m,wx,xx,xx(mp))
!!!!!                auxd(mp)= dg_pis(nnbr,lnbr,m,wx,xx,xx(mp))
!!!!                 auxd(mp)= dgpis_ref(nnbr,lnbr,m,wx,xx,sp,j)
!!!!                 write(ut43,*)  s, i, j, mp, xx(mp), wx(mp),            &
!!!!!                    pis(nnbr,lnbr,s,i,wx,xx,xx(mp)),                   &
!!!!                     g_pis(nnbr,lnbr,m,wx,xx,xx(mp)),                 &
!!!!                     dpis(nnbr,lnbr,s,i,wx,xx,xx(mp)),                  &
!!!!                     dpis_ij(np,ns,s,i,w,x,j)
!!!!!                    dg_pis(nnbr,lnbr,m,wx,xx,xx(mp))
!!!!                     dgpis_ref(nnbr,lnbr,m,wx,xx,sp,j)
!!!!                 write(422,*) m, sp,j,mp,dg_pis(nnbr,lnbr,m,wx,xx,xx(mp))
!!!!              enddo
!!!!          enddo
!!!!       else
!!!!          auxf(mp)= 0.d0
!!!!          auxd(mp)= 0.d0
!!!!       end if


!!!!         do k=1,nbpts
!!!!            aux4(k) = dg_pis(nnbr,lnbr,m,wx,xx,xt(k))
!!!!            write(ut42,*) k, s, i, xt(k), ddppis(k),                   &
!!!!                     dpis(nnbr,lnbr,s,i,wx,xx,xt(k)),                  &
!!!!                     dg_pis(nnbr,lnbr,m,wx,xx,xt(k)) , real(dfc(k)),   &
!!!!                     (dpis(nnbr,lnbr,s,i,wx,xx,xt(k))-ddppis(k)) /     &
!!!!                      dpis(nnbr,lnbr,s,i,wx,xx,xt(k))*100
!!!!         enddo
!!!!
!!!!        m = (nnbr-1)*(s-1)+i
 
 
        return
        end





        subroutine analytical_basis(workdir, wfc, dwfc)
        use util
        include 'open/files_index'
        integer :: aux0, krange
        integer :: i, j, s, k, np, ns, n, m, mp, nmax1, nn,           &
                         sp, nmax2
        real*8 :: r, step, xx1, xx2, scaling,                         &
                             kstep, y0, kk1, kk2

        real*8, dimension(:), allocatable :: xa, xx, xy, xt, xs,      &
                                   wa, wx, wy, wz

        real*8, dimension(:), allocatable :: aux1, aux2, aux3, aux4,  &
                                             auxf, auxd 

        complex*16, dimension(:), allocatable :: auxc

        complex*16, dimension(:), allocatable :: auxc1, auxc2,        &
                                            auxc3, auxc4
        complex*16, dimension(:) :: wfc, dwfc

        character(255) :: workdir


!       real*8 :: pis, dpis, dpis2, g_pis, dg_pis
 
        include 'param'

        nmax1 = (nnbr-1)*lnbr-1
        nmax2 = (nnbr-1)*lnbr+1
        m = (s-1)*(nnbr-1) + i
        mp = m-1

        allocate(xa(nnbr), wa(nnbr))
        allocate(xx(nmax1), wx(nmax1), wz(nmax1))
        allocate(auxf(nmax1), auxd(nmax1))
!       allocate(auxc(nmax1))
 
        allocate(xy(nmax2), wy(nmax2))
        allocate(xt(nbpts))
        allocate(xs(lnbr+1))
 
        xx1 = -xmax/2.d0
        xx2 =  xmax/2.d0
        step = (xx2 - xx1) / (nbpts-1)
        do k=1,nbpts
           xt(k) =  xx1 + (k-1) * step
        enddo
 
        scaling = 0.5d0*xmax/lnbr

        call lobatto_compute( nnbr, xa, wa )

        do sp=1,lnbr
           do j=1,nnbr
              mp = (nnbr-1)*(sp-1)+j
              xy(mp) = xa(j) + (sp-.5d0-.5d0*lnbr)*2
              xy(mp) = xy(mp)*scaling
              if(j.eq.(nnbr)) then
                 wy(mp) = wa(nnbr)+wa(1)
              else
                 wy(mp) = wa(j)
              end if 
              write(14,*) mp, xy(mp), wy(mp)
              if (j.eq.1) then
                 xs(sp) = xy(mp)
              end if
           enddo
           write(15,*) sp, xs(sp)
        enddo
        xs(sp) = xy(mp)
        write(15,*) sp, xs(sp)
 

! ***   The effective basis
        do sp=1,lnbr
           do j=1,nnbr-1
              if((sp.ne.lnbr).or.(j.ne.(nnbr-1))) then
                 mp = (nnbr-1)*(sp-1)+j
                 xx(mp) = xa(j+1) + (sp-.5d0-.5d0*lnbr)*2
                 xx(mp) = xx(mp)*scaling
                 if(j.eq.(nnbr-1)) then
                    wx(mp) = dsqrt(wa(nnbr)+wa(1))
                    wz(mp) = wa(nnbr)+wa(1)
                 else
                    wx(mp) = dsqrt(wa(j+1))
                    wz(mp) = wa(j+1)
                 end if 
              end if 
           enddo
        enddo


! *** compute all the basis vectors and stock them in files
        do s=1,lnbr
           do i=1,nnbr
              call analytical_basis_vectors(workdir, s, i,            &
                                            xa, xs, xt, xx, xy,       &
                                            wa, wx, wy, wz,           &
                                            auxf, auxd, wfc, dwfc)
            
           enddo

!          if((s.ne.lnbr).or.(i.ne.(nnbr-1))) then
!             m = (nnbr-1)*(s-1)+i
!             auxc = auxc + wfc(m) * aux1
!          end if

        enddo

        do i=1,nmax1
           write(70,*) i, auxf(i), wfc(i),  wfc(i) * auxf(i) * wx(i)
           write(71,*) i, auxd(i), dwfc(i), wfc(i) * auxd(i) * wx(i)
        enddo

!
!
!!       if ((mp.ge.1).and.(mp.le.nmax1)) then
!!          if ((s.ne.(lnbr)).or.(i.ne.(nnbr-1))) then
!!             do sp=1,lnbr
!!                do j=1,nnbr-1
!!                   if((sp.ne.lnbr).or.(j.ne.(nnbr-1))) then
!!                      k = (sp-1)*(nnbr-1) + j
!!                      xx(k)=xy(k+1)
!!                      wx(k)=wy(k+1)
!!                      aux1(k)=auxf(k+1)
!!                      aux2(k)=auxd(k+1)
!!                   end if
!!                enddo
!!             enddo
!!
!!!!             auxc1 = auxc1 +  aux1
!!!!             auxc2 = auxc2 +  aux2
!!!!             auxc3 = auxc3 +  aux3
!!!!             auxc4 = auxc4 +  aux4
!!!
!!!
!!!!             auxc1 = auxc1 + (wfc(mp)/wxx(mp)/dsqrt(scaling)) * aux1
!!!!             auxc2 = auxc2 + (wfc(mp)/wxx(mp)/dsqrt(scaling)) * aux2
!!!!             auxc3 = auxc3 + (wfc(mp)/wxx(mp)/dsqrt(scaling)) * aux3
!!!!             auxc4 = auxc4 + (wfc(mp)/wxx(mp)/dsqrt(scaling)) * aux4
!!!
!!!
!!!!             auxc1 = auxc1 + wfc(mp)*dsqrt(wx(mp)) * aux1
!!!!             auxc2 = auxc2 + wfc(mp)*dsqrt(wx(mp)) * aux2
!!!!             auxc3 = auxc3 + wfc(mp)*dsqrt(wx(mp)) * aux3
!!!!             auxc4 = auxc4 + wfc(mp)*dsqrt(wx(mp)) * aux4
!!!!             write(766,*) mp, wfc(mp), wx(mp), wxx(mp)
!!!!          end if
!!!!       end if
!
!
!!       do k=1,nmax1
!!          auxc1 = auxc1 + wfc(k)* aux1
!!          auxc2 = auxc2 + wfc(k)* aux2
!!          auxc3 = auxc3 + wfc(k)* aux3
!!          auxc4 = auxc4 + wfc(k)* aux4
!!       enddo
 
 
        return
        end



        real*8 function pis_ij(wa,i,j)
        implicit none
        integer :: np, ns, s, i, j
        integer :: k
        real*8 :: eps
        real*8, dimension(*) :: wa
  
  
        if (i.eq.j) then
           pis_ij = 1.d0
        else
           pis_ij = 0.d0
        end if
      
        pis_ij = pis_ij/dsqrt(wa(i))             
  
        return 
        end
 
 
       real*8 function pis_(np,ns,s,i,xx,xy,wa,r)
       implicit none
       integer :: np, ns, s, i, npp
       integer :: j, k
       real*8 :: r, eps
       real*8, dimension(*) :: xx, xy, wa
       real*8, dimension(np) :: rr, ww
       real*8, dimension(ns) :: rs


       include 'param'
 
       call sector_mapper(np,ns,s,xx,xy,wa,rr,2)
       pis_ = 1.d0


       if ((r.ge.xy((s-1)*(nnbr-1)+1)).and.                             &
                  (r.le.xy((s-1)*(nnbr-1)+nnbr))) then
          do j=1,np
             if (i/=j) then
                pis_ = pis_ * (r-rr(j))/(rr(i)-rr(j))
             end if
          enddo
       else
          pis_ = 0.d0
       end if

       pis_ = pis_ / dsqrt(wa(i))
 
       return 
       end


 
 
         real*8 function g_pis_(np,ns,n,xx,xy,wa,r)
         implicit none
         integer :: np,ns, n, m, mp
         integer :: i, j, k, s, aux1, sp, sm
         real*8 :: aux0, xmin
         real*8, dimension(*) :: xx, xy,  wa
         real*8, dimension(np) :: rr, ww
         real*8, dimension(ns) :: rs
   
         real*8 :: r
 !       real*8 :: pis
   
         do s=1,ns
!          if ((n-s*(np-1).le.0)) then
           if (n.le.((s-1)*(np-1)+np)) then
              i = n - (s-1)*(np-1)
              exit
           end if
         enddo
         
         if (n.eq.np) then
            write(46,*) s,n,i, xy((s-1)*(np-1)+np)
         end if
     
         if (n.eq.((s-1)*(np-1)+np)) then 
             call sector_mapper(np,ns,s,xx,xy,wa,rr,2)
             if (r.le.xy((s-1)*(np-1)+np)) then
                g_pis_ = dsqrt(wa(np)) * pis_(np,ns,s,np,xx,xy,wa,r)
             else
                g_pis_ = dsqrt(wa(1))  * pis_(np,ns,s+1,1,xx,xy,wa,r)
             end if
  
             g_pis_ = g_pis_ / dsqrt(wa(1)+wa(np))

         else
            g_pis_ = pis_(np,ns,s,i,xx,xy,wa,r)
         end if
 
 
         return 
         end




         real*8 function g_pis(np,ns,n,xx,xy,wa,r)
         implicit none
         integer :: np,ns, n, m, mp, npp
         integer :: i, j, k, s, aux1, sp, sm
         real*8 :: aux0, xmin
         real*8, dimension(*) :: xx, xy,  wa
         real*8, dimension(np) :: rr, ww
         real*8, dimension(ns) :: rs
   
         real*8 :: r
 !       real*8 :: pis


         do s=1,ns
           if ((s.ne.ns).or.(i.ne.(np-1))) then
              if ((n-s*(np-1).le.0)) then
                 i = n - (s-1)*(np-1)
                 exit
              end if
           end if
         enddo

          

         if (n.eq.(np-1)) then
            write(47,*) s,n,i, xx(s*(np-1))
            do k=1,np-1
               write(48,*) s, i, k, n,  xx(s*(np-1)), xy((s-1)*(np-1)+np)
            enddo
            write(48,*)
         end if



    
!        if (i.lt.(np)) then 
            if (n.eq.(s*(np-1))) then 
                call sector_mapper(np,ns,s,xx,xy,wa,rr,1)
                if (r.le.xy((s-1)*(np-1)+np)) then
                   g_pis = dsqrt(wa(np)) * pis(np,ns,s,np-1,xx,xy,wa,r)
                else
                   g_pis = dsqrt(wa(1))  * pis_(np,ns,s+1,1,xx,xy,wa,r)
                end if
  
                g_pis = g_pis / dsqrt(wa(1)+wa(np))

            else
               g_pis = pis(np,ns,s,i,xx,xy,wa,r)
            end if
!        end if
 
 
         return 
         end



!        real*8 function dpis_ij(np,ns,s,i,xa,wa,j)
!        implicit none
!        integer :: np, ns, s, i
!        integer :: j, k, m
!        real*8 :: aux, r
!        real*8, dimension(*) :: wa, xa
!        real*8, dimension(np) :: rr, ww
! 
!!       real*8 :: pis
!   
!        dpis_ij= 0.d0
! 
!        if ((j/=i)) then
!           aux = 1.d0
!           do k=1,np
!              if ((k/=i).and.(k/=j)) then
!                 aux = aux* (xa(j)-xa(k))/(xa(i)-xa(k))
!              end if
!           enddo
!           dpis_ij =  aux / (xa(i)-xa(j))
!        else
!           dpis_ij = 0.5d0 *                                         &
!             ( pis_ij(wa,i,np)**2 - pis_ij(wa,i,1)**2 )
!        end if
! 
!        dpis_ij  = dpis_ij / dsqrt(wa(i))
!  
!        return         
!        end



         real*8 function dpis_ij(np,ns,s,i,xx,xy,wa,j)
         implicit none
         integer :: np, ns, s, i
         integer :: j, k, m
         real*8 :: aux, r
         real*8, dimension(*) :: xx, xy, wa
         real*8, dimension(np) :: rr, ww
  
 !       real*8 :: pis
    
         call sector_mapper(np,ns,s,xx,xy,wa,rr,2)
         dpis_ij= 0.d0
  
         if ((j/=i)) then
            aux = 1.d0
            do k=1,np
               if ((k/=i).and.(k/=j)) then
                  aux = aux* (rr(j)-rr(k))/(rr(i)-rr(k))
               end if
            enddo
            dpis_ij =  aux / (rr(i)-rr(j))

         else
            dpis_ij = 0.5d0 *                                         &
              ( pis_ij(wa,i,np)**2 - pis_ij(wa,i,1)**2 )
         end if
  
         dpis_ij  = dpis_ij / dsqrt(wa(i))
   
         return         
         end
   
  
  
  
        real*8 function dpis_(np,ns,s,i,xx,xy,wa,r)
        implicit none
        integer :: np, ns, s, i, npp
        integer :: j, k
        real*8 :: r, aux
        real*8, dimension(*) :: xx, xy, wa
        real*8, dimension(np) :: rr, ww
        real*8, dimension(ns) :: rs
 
!       real*8 :: pis
 
        include 'param'
 
        call sector_mapper(np,ns,s,xx,xy,wa,rr,2)
        dpis_ = 0.d0

        if ((r.ge.xy((s-1)*(np-1)+1)).and.                            &
                  (r.le.xy((s-1)*(np-1)+np))) then 
           do j=1,np
              aux = 1.d0
              if ((j/=i)) then
                 do k=1,np
                    if ((k/=i).and.(k/=j)) then
                       aux = aux* (r-rr(k))/(rr(i)-rr(k))
                    end if
                 enddo
                 dpis_ = dpis_ + aux / (rr(i)-rr(j))
              end if
           enddo
  
           dpis_ = dpis_ / dsqrt(wa(i))
        end if
  

        return         
        end



        real*8 function dpis(np,ns,s,i,xx,xy,wa,r)
        implicit none
        integer :: np, ns, s, i, npp
        integer :: j, k
        real*8 :: r, aux
        real*8, dimension(*) :: xx, xy, wa
        real*8, dimension(np) :: rr, ww
        real*8, dimension(ns) :: rs
 
!       real*8 :: pis
 
        include 'param'
 
        call sector_mapper(np,ns,s,xx,xy,wa,rr,1)
        dpis = 0.d0


!        if (s.ne.ns) then
!           npp=np-1
!        else
!           npp=np-2
!        end if

        if ((r.gt.xy((s-1)*(np-1)+1)).and.                            &
                  (r.le.xy((s-1)*(np-1)+np))) then 
           do j=1,np-1
              aux = 1.d0
              if ((j/=i)) then
                 do k=1,np-1
                    if ((k/=i).and.(k/=j)) then
                       aux = aux* (r-rr(k))/(rr(i)-rr(k))
                    end if
                 enddo
                 dpis = dpis + aux / (rr(i)-rr(j))
              end if
           enddo
  
           dpis = dpis / dsqrt(wa(i+1))

        end if
  
        return         
        end



         real*8 function dg_pis_(np,ns,n,xx,xy,wa,r)
         implicit none
         integer :: np,ns, n, m, mp
         integer :: i, j, k, s, aux1, sp, sm
         real*8 :: aux0, xmin
         real*8, dimension(*) :: xx, xy, wa
         real*8, dimension(np) :: rr, ww
         real*8, dimension(ns) :: rs
   
         real*8 :: r
 !       real*8 :: pis
   
         do s=1,ns
           if (n.le.((s-1)*(np-1)+np)) then
              i = n - (s-1)*(np-1)
              exit
           end if
         enddo
         
         if (n.eq.np) then
            write(46,*) s,n,i, xy((s-1)*(np-1)+np)
         end if
     
         if (n.eq.((s-1)*(np-1)+np)) then 
             call sector_mapper(np,ns,s,xx,xy,wa,rr,2)
             if (r.lt.xy((s-1)*(np-1)+np)) then
                dg_pis_ = dsqrt(wa(np)) *                             &
                              dpis_(np,ns,s,np,xx,xy,wa,r)
             else if (r.gt.xy((s-1)*(np-1)+np)) then
                dg_pis_ = dsqrt(wa(1))  *                             &
                              dpis_(np,ns,s+1,1,xx,xy,wa,r)
             end if
  
             dg_pis_ = dg_pis_ / dsqrt(wa(1)+wa(np))

         else
            dg_pis_ = dpis_(np,ns,s,i,xx,xy,wa,r)
         end if
 
 
         return 
         end



         real*8 function dg_pis(np,ns,n,xx,xy,wa,r)
         implicit none
         integer :: np,ns, n, m, mp
         integer :: i, j, k, s, aux1, sp, sm
         real*8 :: aux0, xmin
         real*8, dimension(*) :: xx, xy, wa
         real*8, dimension(np) :: rr, ww
         real*8, dimension(ns) :: rs
   
         real*8 :: r
 !       real*8 :: pis
   
         do s=1,ns
           if ((s.ne.ns).or.(i.ne.(np-1))) then
              if ((n-s*(np-1).le.0)) then
                 i = n - (s-1)*(np-1)
                 exit
              end if
           end if
         enddo
         
         if (n.eq.(np-1)) then
            write(47,*) s,n,i, xx(s*(np-1))
            do k=1,np-1
               write(48,*) s, i, k, n,  xx(s*(np-1)), xy((s-1)*(np-1)+np)
            enddo
            write(48,*)
         end if
     
!        if (i.lt.(np)) then 
            if (n.eq.(s*(np-1))) then 
                call sector_mapper(np,ns,s,xx,xy,wa,rr,1)
                if (r.lt.xy((s-1)*(np-1)+np)) then
                   dg_pis = dsqrt(wa(np)) *                           &
                                 dpis(np,ns,s,np-1,xx,xy,wa,r)
                else if (r.gt.xy((s-1)*(np-1)+np)) then
                   dg_pis = dsqrt(wa(1))  *                           &
                                 dpis_(np,ns,s+1,1,xx,xy,wa,r)
                end if
     
                dg_pis = dg_pis / dsqrt(wa(1)+wa(np))
   
            else
               dg_pis = dpis(np,ns,s,i,xx,xy,wa,r)
            end if
!        end if



         return 
         end









 
 
!        real*8 function dg_pis(np,ns,n,w,x,r)
!        implicit none
!        integer :: np, ns, n, s, i
!        integer :: j, k, m, sp, sm
!        real*8 :: aux, r, aux0, aux1
!        real*8, dimension(*) :: w, x
!        real*8, dimension(np) :: rr, ww
!        real*8, dimension(np) :: rrr, www
!  
!!       real*8 :: dpis
!!       real*8 :: pis
!    
!        dg_pis= 0.d0
!  
!         do s=1,ns
!           if (n.le.((s-1)*(np-1)+np)) then
!              i = n - (s-1)*(np-1)
!              exit
!           end if
!         enddo
!  
!  
!         call sector_mapper(np,ns,s+1,w,x,ww,rr)     
!
!         s=sp
!  
!         if (n.eq.((s-1)*(np-1)+np)) then
!  
!            if (s.ne.ns) then
!!              if (r.ge.x((s-1)*(np-1)+1).and.r.le.x((s+1)*(np-1)+1)) then
!               if (r.gt.x((s-1)*(np-1)+1).and.r.lt.x((s+1)*(np-1)+1)) then
!                  if (r.lt.x(n)) then
!                     dg_pis = dpis(np,ns,s,np,w,x,r)
!                  else
!                     dg_pis = dpis(np,ns,s+1,1,w,x,r)
!                  end if
!               else
!                  dg_pis = 0.d0
!               end if
!            else
!!              if (r.ge.x((s-1)*(np-1)+1).and.r.le.x((s)*(np-1)+1)) then
!               if (r.gt.x((s-1)*(np-1)+1).and.r.lt.x((s)*(np-1)+1)) then
!                  dg_pis = dpis(np,ns,s,np,w,x,r)
!               else
!                  dg_pis = 0.d0
!               end if
!            end if 
!  
!            dg_pis = sqrt(ww(1))  * dg_pis /(sqrt(ww(1)+ww(np)))
!   
!         else
!  
!!          if ((r.gt.x((s-1)*(np-1)+1)).and.(r.lt.x((s)*(np-1)+1))) then
!           if ((r.ge.x((s-1)*(np-1)+1)).and.(r.le.x((s)*(np-1)+1))) then
!              dg_pis = dpis(np,ns,s,i,w,x,r)
!           else
!              dg_pis = 0.d0 
!           end if
!  
!         end if
!  
!        return         
!        end
  
  
  




!         real*8 function dgpis_ref(np,ns,n,w,x,sp,j)
!         implicit none
!         integer :: np, ns, sp, s, i
!         integer :: j, k, n
!         real*8 :: aux, r
!         real*8, dimension(*) :: w, x
!         real*8, dimension(np) :: rr, ww
!  
!
!
!         do s=1,ns
!           if (n.le.((s-1)*(np-1)+np)) then
!              i = n - (s-1)*(np-1)
!              exit
!           end if
!         enddo
!
!
!!        if (n.eq.((s-1)*(np-1)+np)) then         
!!           if ((j.eq.1).or.(j.eq.np)) then         
!!              dgpis_ref = 0.d0
!!           else 
!!              if ((sp-s).eq.1) then                      ! to next sect.
!!                 dgpis_ref = dpis2(np,ns,sp,1,w,x,j) 
!!              else if ((sp-s).eq.1) then                 ! from prev sect.
!!                 dgpis_ref = dpis2(np,ns,s,np,w,x,j) 
!!              else if ((s-sp).eq.0) then
!!                 dgpis_ref = dpis2(np,ns,sp,i,w,x,j) 
!!              else 
!!                 dgpis_ref = 0.d0
!!              end if
!!           end if
!
!!           dgpis_ref = dsqrt(ww(i)) *                             &
!!                          dgpis_ref /dsqrt(ww(1)+ww(np))
!!  
!!        else
!!           if ((j.eq.1)) then         
!!              dgpis_ref = 0.d0
!!           else 
!!              if ((s-sp).eq.0) then
!!                 dgpis_ref = dpis2(np,ns,sp,i,w,x,j) 
!!              else
!!                 dgpis_ref = 0.d0
!!              end if
!!           end if
!!        end if
!
!         dgpis_ref = 1.d0
!         return         
!         end
 
 
        subroutine sector_mapper(np,ns,s,xx,xy,wa,rr,version)
        implicit none
        integer :: np, m, ns, npp, version
        integer :: j, s, sp
        real*8, dimension(*) :: xx, xy, wa
        real*8, dimension(np) :: rr, ww

        include 'param'

        rr=0.0d0
        ww=0.0d0

        if (version.eq.2) then               ! full
          npp=np

           write(72,*) "full"
           do sp=1,s
              do j=1,npp
                 m = (sp-1)*(np-1) + j
                 rr(j)=xy(m)
!                ww(j)=wa(m)
                 if (sp.eq.s) then
                    write(72,*) sp, j, rr(j)
                 end if
              enddo
           enddo
           
        else if (version.eq.1) then     ! effective 

           if (s.ne.ns) then
              npp=np-1
           else
              npp=np-2
           end if

           write(73,*) "effective"
           do sp=1,s
              do j=1,npp
                 if((sp.ne.ns).or.(j.ne.(np-1))) then
                    m = (sp-1)*(np-1) + j
                    rr(j)=xx(m)
!                   ww(j)=wx(m)
                    if (sp.eq.s) then
                       write(73,*) sp, j, rr(j)
                    end if
                 end if
              enddo
           enddo

        end if


       return
       end




        subroutine sector_mapper2(np,ns,s,w,x,ww,rr)
        implicit none
        integer :: np, m, ns, npp
        integer :: j, s, sp
        real*8, dimension(*) :: x, w
        real*8, dimension(np) :: rr, ww
 
        do sp=1,s
           do j=1,np-1
              if((sp.ne.ns).or.(j.ne.(np-1))) then
                 m = (sp-1)*(np-1) + j
                 rr(j)=x(m)
                 ww(j)=w(m)
              end if
           enddo
        enddo



!       if (s.ne.ns) then
!          npp=np-1
!       else
!          npp=np-2
!       end if

!       do sp=1,s
!          do j=1,npp
!             m = (sp-1)*(np-1) + j
!             rr(j)=x(m)
!             ww(j)=w(m)
!          enddo
!       enddo
 


       return
       end



      end module
