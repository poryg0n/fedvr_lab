      module fedvr_basis
      use lobatto
      use math_util, only :  differentiation_real, fourier_trsfm,     &
                            inv_fourier_trsfm, fourier_diff
      implicit none
      contains
      

        subroutine analytical_basis_vectors(workdir, s, i,            &
                                             xa, xs, xt, xx, xy,      &
                                             wa, wx, wy, wz,          &
                                             gppi_d, dgppi_d,         &
                                             gppi, dgppi)
        use util
        include 'open/files_index'
        integer :: aux0, krange
        integer :: i, j, s, k, np, ns, n, m, mp, nmax1, nmax2, nn,    &
                         sp, npp, ip, mk, mq, q
        real*8 :: r, step, xx1, xx2, scaling,                         &
                             step_, kstep, y0, kk1, kk2
        real*8, dimension(:), allocatable :: xa, xt, xx, xy, xs,  &
                                     wa, wx, wy, wz

        real*8, dimension(:) :: gppi_d, dgppi_d, gppi, dgppi

        real*8, dimension(:), allocatable :: aux1, aux2, aux3, aux4,  &
                                             aux
        real*8, dimension(:), allocatable :: ppis, ppis_,             &
                                    gppis, gppis_,                    &
                                    dppis, dppis_,                    &
                                    dgppis, dgppis_,                  &
                                    ff, df, dff, dx


        complex*16, dimension(:), allocatable :: ffc, ffk, ffc1,      &
                              dfc

        character(len=1024) :: filename01='quad_pts_wts'
        character(len=1024) :: filename02='ppi_'
        character(len=1024) :: filename03='ppic_'
        character(len=1024) :: filename04='basis_vect.dat'
        character(len=1024) :: filename05='basis_vect_deriv.dat'


        character(len=1024) :: filename1
        character(len=1024) :: filename2
        character(len=1024) :: filename3
        character(len=1024) :: filename4
        character(len=1024) :: filename5

        character(len=1024) :: filename1_
        character(len=1024) :: filename2_
        character(len=1024) :: filename3_
!       character(len=1024) :: filename4_
!       character(len=1024) :: filename5_


        character(len=1024) :: format_string
        character(len=1024) :: indx, indx1, indx2 
        character(255) :: workdir


!       real*8 :: pis, dpis, dpis2, g_pis, dg_pis
 
        include 'param'
        include 'common/open_files'
        include 'open/files_specs'

!       include 'common/laser'
        include 'open/local_basis_files'

!       call open_files_to_store(0)

!       nbpts=200
 
        nmax1 = lnbr*(nnbr-1) - 1
        nmax2 = lnbr*(nnbr-1) + 1
        xx1 = -xmax/2.d0
        xx2 =  xmax/2.d0
 
        scaling = 0.5d0*xmax/lnbr
 
 
       allocate(ppis(nbpts),  ppis_(nbpts))
       allocate(gppis(nbpts), gppis_(nbpts))
       allocate(dgppis(nbpts),dgppis_(nbpts))

       allocate(aux1(nmax1), aux2(nmax1), aux(nmax1))

       allocate(aux3(nbpts))
       allocate(aux4(nbpts))
       step = (xx2 - xx1) / (nbpts-1)

       allocate(ff(nbpts), df(nbpts), dff(nbpts),                   &
                                        dppis(nbpts), dx(nbpts))

       m = (s-1)*(nnbr-1) + i

       if (i.lt.(nnbr)) then
          if ((i.ne.(nnbr-1)).or.(s.ne.lnbr)) then     
             do sp=1,lnbr
                do j=1,nnbr-1
                   if((sp.ne.lnbr).or.(j.ne.(nnbr-1))) then
                      k = j+1
                      mp = (sp-1)*(nnbr-1)+j
                      mk = mp+1
                      if (m.eq.1) then
                         write(ut80,*)  sp, mp, k,                     & 
                             dsqrt(wa(k)), xa(k), xx(mp),              &
!                            pis(nnbr,lnbr,sp,k,xx,xy,wa,xx(mp)),      &
                             pis_ij(wa,k,k),                           &
!                            gpis_nm(nnbr,lnbr,wa,mk,mk),              &
                             g_pis(nnbr,lnbr,mk,xx,xy,wa,xx(mp)),      &
!                            dpis_ij(nnbr,lnbr,sp,k,xx,xy,wa,k),       &
                             dpis(nnbr,lnbr,sp,k,xx,xy,wa,xx(mp)),     &
                             dg_pis(nnbr,lnbr,mk,xx,xy,wa,xx(mp)),     &
                             dgpis_nm(nnbr,lnbr,mk,xx,xy,wa,mk),       &
                             0.d0
                      end if

                      q = i + 1 
                      mq = m + 1 
                      gppi_d(mp) = g_pis(nnbr,lnbr,mq,xx,xy,wa,xx(mp))
                      dgppi_d(mp) = dg_pis(nnbr,lnbr,mq,xx,xy,wa,xx(mp))

!          ***        not working
!                     gppi_d(mp) = gpis_nm(nnbr,lnbr,wa,mq,mk)
!                     dgppi_d(mp) = dgpis_nm(nnbr,lnbr,mq,xx,xy,wa,mk)

                      write(ut81,*)   mq, sp, q, j, xx(mp),              &
                          pis(nnbr,lnbr,s,q,xx,xy,wa,xx(mp)),            &
                          pis_ij(wa,q,j+1),                              &
!                         gpis_nm(nnbr,lnbr,wa,mq,mk),                   &
                          g_pis(nnbr,lnbr,mq,xx,xy,wa,xx(mp)),           &
                          dpis_ij(nnbr,lnbr,s,q,xx,xy,wa,j+1),           &
                          dpis(nnbr,lnbr,s,q,xx,xy,wa,xx(mp)),           &
                          dg_pis(nnbr,lnbr,mq,xx,xy,wa,xx(mp)),          &
                          dgpis_nm(nnbr,lnbr,mq,xx,xy,wa,mp+1)
                   end if
                enddo
                write(ut81,*)
             enddo


! ***     continuous \pi and \Pi
             q = i + 1 
             mq = m + 1 
             do k=1,nbpts
                xt(k) = xx1 + (k-1)*step
                ppis(k)    = pis(nnbr,lnbr,s,q,xx,xy,wa,xt(k))
                gppis(k)   = g_pis(nnbr,lnbr,mq,xx,xy,wa,xt(k))
                dgppis(k)  = dg_pis(nnbr,lnbr,mq,xx,xy,wa,xt(k))
                write(ut82,*) k, lnbr, nnbr, s, q, mq, xt(k),         &
                   pis(nnbr,lnbr,s,q,xx,xy,wa,xt(k)),                 &
                   g_pis(nnbr,lnbr,mq,xx,xy,wa,xt(k)),                &
                   dpis(nnbr,lnbr,s,q,xx,xy,wa,xt(k)),                &
                   dg_pis(nnbr,lnbr,mq,xx,xy,wa,xt(k))
             enddo
          end if
       end if


!      write(ut90,*) "phase 1 : Find the choke points", lnbr*(nnbr-1)+1

       do sp=1,lnbr
          do j=1,nnbr
             mp = (sp-1)*(nnbr-1)+j
             if (m.eq.1) then
                write(ut90,*)   mp, sp, j,                            &
                   dsqrt(wa(j)), xa(j), xy(mp),                       &
!                  pis(nnbr,lnbr,sp,j,xx,xy,wa,xy(mp)),               &
                   pis_ij(wa,j,j),                                    &
                   gpis_nm(nnbr,lnbr,wa,mp,mp),                       &
!                  g_pis(nnbr,lnbr,mp,xx,xy,wa,xy(mp)),               &
                   dpis_ij(nnbr,lnbr,sp,j,xx,xy,wa,j),                &
!                  dpis(nnbr,lnbr,sp,j,xx,xy,wa,xy(mp)),              &
                   dg_pis(nnbr,lnbr,mp,xx,xy,wa,xy(mp)),              &
!                  dgpis_nm(nnbr,lnbr,mp,xx,xy,wa,mp),                &
                   0.d0
             end if

             write(ut91,*)   m, sp, i, j, xy(mp),                     &
                pis_ij(wa,i,j),                                       &
!               pis(nnbr,lnbr,s,i,xx,xy,wa,xy(mp)),                   &
                gpis_nm(nnbr,lnbr,wa,m,mp),                           &
                g_pis(nnbr,lnbr,m,xx,xy,wa,xy(mp)),                   &
                dpis_ij(nnbr,lnbr,s,i,xx,xy,wa,j),                    &
                dpis(nnbr,lnbr,s,i,xx,xy,wa,xy(mp)),                  &
                dg_pis(nnbr,lnbr,m,xx,xy,wa,xy(mp)),                  &
                dgpis_nm(nnbr,lnbr,m,xx,xy,wa,mp)
          enddo
          write(ut90,*)
          write(ut91,*)
       enddo

! *** continuous \pi and \Pi
      do k=1,nbpts
         xt(k) = xx1 + (k-1)*step
         ppis_(k)    = pis(nnbr,lnbr,s,i,xx,xy,wa,xt(k))
         gppis_(k)   = g_pis(nnbr,lnbr,m,xx,xy,wa,xt(k))
         dgppis_(k)  = dg_pis(nnbr,lnbr,m,xx,xy,wa,xt(k))
         write(ut92,*) k, lnbr, nnbr, s, i, m,  xt(k),                &
                   pis(nnbr,lnbr,s,i,xx,xy,wa,xt(k)),                &
                   g_pis(nnbr,lnbr,m,xx,xy,wa,xt(k)),                &
                   dpis(nnbr,lnbr,s,i,xx,xy,wa,xt(k)),               &
                   dg_pis(nnbr,lnbr,m,xx,xy,wa,xt(k))
      enddo


!     wwfc  =  wwfc  + wfc(m) * wx(m) * gppis_
!     dwwfc  = dwwfc + wfc(m) * wx(m) * dgppis_
 
      return
      end





        subroutine write_analytical_basis(workdir)
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
                                            gppi_d, dgppi_d,          &
                                            gppi, dgppi

        complex*16, dimension(:), allocatable :: auxc

        complex*16, dimension(:), allocatable :: auxc1, auxc2,        &
                                            auxc3, auxc4,             &
                                            wwfc, dwwfc

!       complex*16, dimension(:) :: wfc, dwfc

        character(255) :: workdir
        character(len=1024) :: filename4
        character(len=1024) :: filename5
        character(len=1024) :: filename04='basis_vect.dat'
        character(len=1024) :: filename05='basis_vect_deriv.dat'

        include 'param'

        include 'common/open_files'
        include 'open/files_specs'

!       include 'open/basis_vectors_files'
        open(ut83,file=trim(workdir)//trim(filename04),status='new')
        open(ut84,file=trim(workdir)//trim(filename05),status='new')

        nmax1 = (nnbr-1)*lnbr-1
        nmax2 = (nnbr-1)*lnbr+1
        m = (s-1)*(nnbr-1) + i
        mp = m-1


        allocate(xa(nnbr), wa(nnbr))
        allocate(xx(nmax1), wx(nmax1), wz(nmax1))
!       allocate(wwfc_d(nmax1), dwwfc_d(nmax1))
!       allocate(auxc(nmax1))
 
        allocate(xy(nmax2), wy(nmax2))
        allocate(xt(nbpts))
        allocate(wwfc(nbpts))
        allocate(dwwfc(nbpts))

        allocate(gppi_d(nmax1), dgppi_d(nmax1))
        allocate(gppi(nbpts), dgppi(nbpts))

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
                                            gppi_d, dgppi_d,          &
                                            gppi, dgppi)

              m = (s-1)*(nnbr-1) + i
              if (i.lt.(nnbr)) then
                if ((i.ne.(nnbr-1)).or.(s.ne.lnbr)) then

                    write(ut83,*) s, i, m, gppi_d
                    write(ut84,*) s, i, m, dgppi_d

!                   write(ut85,*) s, i, m, gppi
!                   write(ut86,*) s, i, m, dgppi
                 end if
              end if

           enddo
        enddo

        close(ut83)
        close(ut84)

        return
        end





        subroutine compute_analytic_deriv(workdir, wx, xx, dwwfc,     &
                                                    wwfc, wfc)
        use util
        include 'open/files_index'
        integer :: aux0, krange
        integer :: i, j, s, k, np, ns, n, m, mp, nmax1, nn,           &
                         sp, nmax2, ios
        real*8 :: r, step, xx1, xx2, scaling,                         &
                             kstep, y0, kk1, kk2

        real*8, dimension(:), allocatable :: xa, wa, xx,              &
                                   wx, wy, wz

        real*8, dimension(:), allocatable :: aux1, aux2, aux3, aux4

        complex*16, dimension(:), allocatable :: auxc

        complex*16, dimension(:), allocatable :: auxc1, auxc2,        &
                                            auxc3, auxc4

        complex*16, dimension(:) :: dwwfc, wwfc, wfc

        character(255) :: workdir
        character(len=1024) :: filename04='basis_vect.dat'
        character(len=1024) :: filename05='basis_vect_deriv.dat'
        character(len=1024) :: filename07='analytic_wfc_dwfc_disc.dat'
        character(len=1024) :: filename08='analytic_wfc_dwfc.dat'


        include 'param'

        include 'common/open_files'
        include 'open/files_specs'

        open(ut70,file=trim(workdir)//trim(filename07),status='unknown')
        open(ut71,file=trim(workdir)//trim(filename08),status='unknown')
        open(ut83,file=trim(workdir)//trim(filename04),status='old')
        open(ut84,file=trim(workdir)//trim(filename05),status='old')

        nmax1 = (nnbr-1)*lnbr-1
        nmax2 = (nnbr-1)*lnbr+1

        wwfc = c0
        dwwfc = c0


        allocate(aux1(nmax1), aux2(nmax1))

        do mp=1,nmax1
           read(ut83,*, iostat=ios) s, i, m, aux1
           read(ut84,*, iostat=ios) s, i, m, aux2
           if (ios /= 0) exit
           wwfc(mp) = wx(mp) * wfc(mp) * aux1(mp)
           do i=1,nmax1
              dwwfc(mp) = dwwfc(mp) +  wx(mp) * wfc(mp) * aux2(i)
           enddo
           write(ut70,*) mp, wwfc(mp), dwwfc(mp)
        enddo

!       write(*,*) s, i, m, aux1
      
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


 
 
       real*8 function pis(np,ns,s,i,xx,xy,wa,r)
       implicit none
       integer :: np, ns, s, i, npp
       integer :: j, k
       real*8 :: r, eps
       real*8, dimension(*) :: xx, xy, wa
       real*8, dimension(np) :: rr, ww
       real*8, dimension(ns) :: rs


       include 'param'
 
       call sector_mapper(np,ns,s,xx,xy,wa,rr)
       pis = 1.d0

       if ((r.ge.xy((s-1)*(nnbr-1)+1)).and.                             &
                  (r.le.xy((s-1)*(nnbr-1)+nnbr))) then
          do j=1,np
             if (i/=j) then
                pis = pis * (r-rr(j))/(rr(i)-rr(j))
             end if
          enddo
       else
          pis = 0.d0
       end if

       pis = pis / dsqrt(wa(i))
 
       return 
       end



        real*8 function gpis_nm(np,ns,wa,n,m)
        implicit none
        integer :: np, ns, s, i, j, n, m
        integer :: k, sp
        real*8 :: eps
        real*8, dimension(*) :: wa

        include 'param'


        do s=1,ns
          if (n.le.((s-1)*(np-1)+np)) then
             i = n - (s-1)*(np-1)
             exit
          end if
        enddo

        do sp=1,ns
          if (m.le.((sp-1)*(np-1)+np)) then
             j = m - (sp-1)*(np-1)
             exit
          end if
        enddo


        if ((n-m).eq.0) then
          gpis_nm = 1.d0 / dsqrt(wa(i))
        else 
          gpis_nm = 0.0d0
        end if

        if (n.eq.((s-1)*(np-1)+np)) then 
           if ((sp-s).le.1) then
              gpis_nm = dsqrt(wa(np)) * gpis_nm
           else if ((s-sp).le.1) then
              gpis_nm = dsqrt(wa(1)) * gpis_nm
           end if

           gpis_nm = gpis_nm /dsqrt(wa(1) + wa(np))


        end if
      
        return 
        end


 
 
         real*8 function g_pis(np,ns,n,xx,xy,wa,r)
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
!            call sector_mapper(np,ns,s,xx,xy,wa,rr)
             if (r.le.xy((s-1)*(np-1)+np)) then
                g_pis = dsqrt(wa(np)) * pis(np,ns,s,np,xx,xy,wa,r)
             else
                g_pis = dsqrt(wa(1))  * pis(np,ns,s+1,1,xx,xy,wa,r)
             end if
  
             g_pis = g_pis / dsqrt(wa(1)+wa(np))

         else
            g_pis = pis(np,ns,s,i,xx,xy,wa,r)
         end if
 
 
         return 
         end





         real*8 function dpis_ij(np,ns,s,i,xx,xy,wa,j)
         implicit none
         integer :: np, ns, s, i
         integer :: j, k, m
         real*8 :: aux, r
         real*8, dimension(*) :: xx, xy, wa
         real*8, dimension(np) :: rr, ww
  
 !       real*8 :: pis
    
         call sector_mapper(np,ns,s,xx,xy,wa,rr)
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
 
        call sector_mapper(np,ns,s,xx,xy,wa,rr)
        dpis = 0.d0

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
                 dpis = dpis + aux / (rr(i)-rr(j))
              end if
           enddo
  
           dpis = dpis / dsqrt(wa(i))
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
           if (n.le.((s-1)*(np-1)+np)) then
              i = n - (s-1)*(np-1)
              exit
           end if
         enddo
         
!        if (n.eq.np) then
!           write(46,*) s,n,i, xy((s-1)*(np-1)+np)
!        end if
     
         if (n.eq.((s-1)*(np-1)+np)) then 
!            call sector_mapper(np,ns,s,xx,xy,wa,rr)

             if ((r.lt.xy((s-1)*(np-1)+np)).and.                      &
                                  (r.gt.xy((s-1)*(np-1)+1))) then
                dg_pis = dsqrt(wa(np)) *                              &
                              dpis(np,ns,s,np,xx,xy,wa,r)
             else if ((r.gt.xy((s-1)*(np-1)+np)).and.                 &
                                     (r.lt.xy(s*(np-1)+np))) then
                dg_pis = dsqrt(wa(1))  *                              &
                              dpis(np,ns,s+1,1,xx,xy,wa,r)
             end if
  
             dg_pis = dg_pis / dsqrt(wa(1)+wa(np))

         else
            if (r.ge.xy((s-1)*(np-1)+1).and.                          &
                                (r.le.xy((s-1)*(np-1)+np))) then     
               dg_pis = dpis(np,ns,s,i,xx,xy,wa,r)
            end if
         end if
 
 
         return 
         end



 
 
        real*8 function dgpis_nm(np,ns,n,xx,xy,wa,m)
        implicit none
        integer :: np, ns, sp, s, i
        integer :: j, k, n, m
        real*8 :: aux, r
        real*8, dimension(*) :: wa, xx, xy
        real*8, dimension(np) :: rr, ww
   
 
        do s=1,ns
          if (n.le.((s-1)*(np-1)+np)) then
             i = n - (s-1)*(np-1)
             exit
          end if
        enddo

        do sp=1,ns
          if (m.le.((sp-1)*(np-1)+np)) then
             j = m - (sp-1)*(np-1)
             exit
          end if
        enddo

        write(25,*) s, i, sp, j 

!!       if (abs(s-sp+1).le.1) then
!        if (abs(s-sp).le.1) then
!
!           if (i.eq.np) then
!              if ((sp-s).eq.1) then
!                 dgpis_nm = dsqrt(wa(1)) * dpis_ij(np,ns,sp,1,xx,xy,wa,j) 
!              else if ((s-sp).eq.0) then
!                 dgpis_nm = dsqrt(wa(np)) * dpis_ij(np,ns,s,np,xx,xy,wa,j) 
!              end if
!
!              dgpis_nm = dgpis_nm /dsqrt(wa(1) + wa(np))
!
!           else if (i.eq.1) then
!              if ((s-sp).eq.1) then
!                 dgpis_nm = dsqrt(wa(np)) * dpis_ij(np,ns,sp,np,xx,xy,wa,j) 
!              else if ((s-sp).eq.0) then
!                 dgpis_nm = dsqrt(wa(1)) * dpis_ij(np,ns,s,1,xx,xy,wa,j) 
!              end if
!
!              dgpis_nm = dgpis_nm /dsqrt(wa(1) + wa(np))
!
!           else
!              if ((s-sp).eq.0) then
!                dgpis_nm = dpis_ij(np,ns,s,i,xx,xy,wa,j)
!              end if
!           end if
!        else 
!             dgpis_nm = 0.0d0
!        end if

!       if (i.eq.np) then
!          if ((sp-s).eq.1) then
!             dgpis_nm = dsqrt(wa(1))  * dpis_ij(np,ns,sp,1,xx,xy,wa,j) 
!          else if ((sp-s).eq.0) then
!             dgpis_nm = dsqrt(wa(np)) * dpis_ij(np,ns,s,np,xx,xy,wa,j) 
!          end if

!          dgpis_nm = dgpis_nm /dsqrt(wa(1) + wa(np))

!       else if (i.eq.1) then
!          if ((s-sp).eq.1) then
!             dgpis_nm = dsqrt(wa(np)) * dpis_ij(np,ns,sp,np,xx,xy,wa,j) 
!             dgpis_nm = dgpis_nm /dsqrt(wa(1) + wa(np))
!          else if ((sp-s).eq.0) then
!             dgpis_nm = dpis_ij(np,ns,s,i,xx,xy,wa,j)
!          end if

!       else
!          if ((s-sp).eq.0) then
!             dgpis_nm = dpis_ij(np,ns,s,i,xx,xy,wa,j)
!          else
!             dgpis_nm = 0.d0
!          end if
!       end if


        if ((i.eq.np).and.((sp-s).eq.0)) then
           dgpis_nm = dsqrt(wa(np)) * dpis_ij(np,ns,s,np,xx,xy,wa,j)
           dgpis_nm = dgpis_nm /dsqrt(wa(1) + wa(np))
        else if ((i.eq.np).and.((sp-s).eq.1)) then
           dgpis_nm = dsqrt(wa(1)) * dpis_ij(np,ns,sp,1,xx,xy,wa,j)
           dgpis_nm = dgpis_nm /dsqrt(wa(1) + wa(np))
        else if ((i.eq.1).and.((s-sp).eq.1)) then
           dgpis_nm = dsqrt(wa(np)) * dpis_ij(np,ns,sp,np,xx,xy,wa,j)
           dgpis_nm = dgpis_nm /dsqrt(wa(1) + wa(np))
        else if ((s-sp).eq.0) then
           dgpis_nm = dpis_ij(np,ns,s,i,xx,xy,wa,j)
        else
           dgpis_nm = 0.d0
        end if
 
        return         
        end
 
 
        subroutine sector_mapper(np,ns,s,xx,xy,wa,rr)
        implicit none
        integer :: np, m, ns, npp, version
        integer :: j, s, sp
        real*8, dimension(*) :: xx, xy, wa
        real*8, dimension(np) :: rr, ww

        include 'param'

        rr=0.0d0
        ww=0.0d0

        write(72,*) "full"
        do sp=1,s
           do j=1,np
              m = (sp-1)*(np-1) + j
              rr(j)=xy(m)
              ww(j)=wa(m)
              if (sp.eq.s) then
                 write(72,*) sp, j, rr(j)
              end if
           enddo
        enddo
           

       return
       end



        subroutine correct_the_points_(xx, xt, dwwfc_d, dwwfc)
        implicit none
        integer :: nmax
        integer :: i, j, k, sp
        real*8 :: xtest, aux
        real*8, dimension(*) :: xx, xt

        real*8, dimension(*) :: dwwfc_d
        complex*16, dimension(*) :: dwwfc

        include 'param'

        nmax = (nnbr-1)*lnbr-1

        do i=1,nmax
           write(88,*) i, dwwfc_d(i)
        enddo

        do k=1,nbpts
           write(89,*) k, dwwfc(k)
        enddo

          
 
        do sp=1,lnbr/2-1
           aux = dwwfc_d(sp*(nnbr-1))
           xtest = 1.d0
           do k=1,nbpts/2
              if (abs(xx(sp*(nnbr-1))-xt(k)).lt.xtest) then 
                 dwwfc_d(sp*(nnbr-1)) = dwwfc(k) 
                 xtest = xx(sp*(nnbr-1))-xt(k)
              end if 
           enddo
        enddo


        do sp=lnbr/2+1, lnbr-1
           aux = dwwfc_d(sp*(nnbr-1))
           xtest = 1.d0
           do k=nbpts/2+1,nbpts
              if (abs(xx(sp*(nnbr-1))-xt(k)).lt.xtest) then 
                 dwwfc_d(sp*(nnbr-1)) = dwwfc(k) 
                 xtest = xx(sp*(nnbr-1))-xt(k)
              end if 
           enddo
        enddo

       return
       end


       subroutine correct_the_points(m, xx, xt, dgpis_d, dgpis)
       implicit none
       integer :: nmax
       integer :: i, j, k, sp, m
       real*8 :: xtest, aux
       real*8, dimension(*) :: xx, xt

       real*8, dimension(*) :: dgpis_d
       real*8, dimension(*) :: dgpis

       include 'param'

       nmax = (nnbr-1)*lnbr-1

       if (mod(m,(nnbr-1)).eq.0) then
           aux = dgpis_d(m)
           xtest = 1.d0
           if (xx(m).le.0) then
              do k=1,nbpts/2
                 if ((abs(xx(m))-xt(k)).lt.xtest) then 
                    dgpis_d(m) = dgpis(k) 
                    xtest = xx(m)-xt(k)
                 end if 
              enddo
           else 
              do k=nbpts/2,nbpts
                 if ((abs(xx(m))-xt(k)).lt.xtest) then 
                    dgpis_d(m) = dgpis(k) 
                    xtest = xx(m)-xt(k)
                 end if 
              enddo
           end if 
       end if

       return
       end



      end module
