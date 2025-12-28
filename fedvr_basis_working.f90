      module fedvr_basis
      use lobatto
      use math_util, only :  differentiation_real, fourier_trsfm,     &
                            inv_fourier_trsfm, fourier_diff
      implicit none
      contains
      

        subroutine analytical_basis(workdir,s,i,xx,wx,                &
                            auxf,auxd,aux3,aux4)
        use util
        include 'open/files_index'
        integer :: aux0, krange
        integer :: i, j, s, k, np, ns, n, m, mp, nmax, nn,            &
                         sp
        real*8 :: r, step, xx1, xx2, scaling,                         &
                             step_, kstep, y0, kk1, kk2
        real*8, dimension(:), allocatable :: xa, xs, wa, kk
        real*8, dimension(:) :: xx, wx, auxf, auxd, aux3, aux4
        real*8, dimension(:), allocatable :: ppis, gppis, ddppis,     &
                                    ff, df, dff, dx, xt
        complex*16, dimension(:), allocatable :: ffc, ffk, ffc1,      &
                              dfc, ppisc 

        character(len=1024) :: filename10='quad_pts_wts.dat'
        character(len=1024) :: filename20='ppi_'
        character(len=1024) :: filename30='dpi_'
        character(len=1024) :: filename40='dp2_'

        character(len=1024) :: filename1
        character(len=1024) :: filename2
        character(len=1024) :: filename3
        character(len=1024) :: filename4

        character(len=1024) :: format_string
        character(len=1024) :: indx
        character(255) :: workdir


!       real*8 :: pis, dpis, dpis2, g_pis, dg_pis
 
        include 'param'
        include 'common/open_files'
!       include 'open/files_specs'

!       include 'common/laser'
        include 'open/files_deriv_analyt'

!       call open_files_to_store(0)

        allocate(xa(nnbr),wa(nnbr), xs(lnbr), xt(nbpts))
!       nbpts=200
 
!       ut=100
!       ut=200
!       ut=300
!       ut=400
 
        nmax = lnbr*(nnbr-1) + 1
        xx1 = -xmax/2.d0
        xx2 =  xmax/2.d0
 
        scaling = 0.5d0*xmax/lnbr
 
 
        write(*,*) "np = ", nnbr
        write(*,*) "ns = ", lnbr
        write(*,*) "nmax = ", nmax
        write(*,*)

       call lobatto_compute( nnbr, xa, wa )

       do j=1,nnbr
          write(*,*) nnbr, j, xa(j), wa(j)
       enddo
       write(333,*) 
       write(*,*) "do i get here?", nnbr, lnbr, nmax
 

        do sp=1,lnbr
           do j=1,nnbr
              mp = (nnbr-1)*(sp-1)+j
              xx(mp) = xa(j) + (sp-.5d0-.5d0*lnbr)*2
              xx(mp) = xx(mp)*scaling
              if(j.eq.(nnbr)) then
                 wx(mp) = dsqrt(wa(nnbr)+wa(1))
              else
                 wx(mp) = dsqrt(wa(j))
              end if 
              write(*,*) mp, xx(mp), wx(mp)
           enddo
           write(*,*)
        enddo

       allocate(ppis(nbpts), gppis(nbpts), ppisc(nbpts))
       step = (xx2 - xx1) / (nbpts-1)

       allocate(ff(nbpts), df(nbpts), dff(nbpts),                   &
                                        ddppis(nbpts), dx(nbpts))

      write(ut40,*) "phase 1 : Find the choke points", m
      do sp=1,lnbr
         do j=1,nnbr
            m = (sp-1)*(nnbr-1)+j
            write(ut40,*)   m, sp, j, xx(m), wx(m),                    &
                   pis(nnbr,lnbr,sp,j,wx,xx,xx(m)),                    &
                   g_pis(nnbr,lnbr,m,wx,xx,xx(m)), 0.d0
!           write(*,*)     m, sp, j, xx(m), wx(m), pis(np,ns,sp,j,wx,xx,xx(m)),
!           g_pis(np,ns,m,wx,xx,xx(m))
         enddo
         write(ut40,*)
      enddo


!      write(*,*) "phase 2 : Lets model the sector changes"
!      sp=1
!      do m=1,nmax-1
!         aux0 = m-sp*(np-1)-1
!         aux1 = mod(m,sp*np)
!         j = -(sp-1)*(np-1)+m
!         write(*,*) m, sp, aux0, j
!         if (aux0.ge.-1) then
!            sp=sp+1
!            write(*,*)
!         end if
!      enddo

 
!        write(indx, '(I0,"_",I0)'), s, i
!        indx = trim(indx)//".dat"
!!       print *, trim(indx)
! 
!        filename2 = trim(filename2)//indx
!        filename3 = trim(filename3)//indx
!        filename4 = trim(filename4)//indx
!
!        open(ut41,file=trim(workdir)//trim(filename2),status='unknown')
!        open(ut42,file=trim(workdir)//trim(filename3),status='unknown')
!        open(ut43,file=trim(workdir)//trim(filename4),status='unknown')

        m=(s-1)*(nnbr-1) + i
        do k=1,nbpts
           xt(k) = xx1 + (k-1)*step
           ppis(k) = pis(nnbr,lnbr,s,i,wx,xx,xt(k))                
           aux3(k) = g_pis(nnbr,lnbr,m,wx,xx,xt(k))
           write(ut41,*) k, lnbr, nnbr, s, i, xt(k), ppis(k),         &
                     g_pis(nnbr,lnbr,m,wx,xx,xt(k))
!          write(*,*)   k, lnbr, nnbr, s, i, xt(k), ppis(k),              &
!                    g_pis (nnbr,lnbr,m,wx,xx,xt(k))
        enddo
 
 
! ***  set dvr function
       krange=20000
       allocate(kk(krange), ffc(nbpts), ffk(krange), ffc1(nbpts),   &
                      dfc(nbpts)) 

!!       ppi = 4.d0*datan(1.d0)
!        kk1=-2.d0*ppi
!        kk2=-kk1
!        kstep = (kk2-kk1)/(krange-1)
!        do k=1,krange
!           kk(k) = kk1 + (k-1)*kstep
!        enddo
!        do k=1,nbpts
!           xt(k) = xx1 + (k-1)*step
!           ffc(k) = 1.d0
!        enddo
!        call fourier_trsfm(nbpts, krange, xt, kk, ffc, ffk)
!        call inv_fourier_trsfm(krange, nbpts, kk, xt, ffk, ffc1)
! 
! 
!         do k=1,krange
!            write(23,*) k, kk(k), real(ffk(k))
!         enddo
!         do k=1,nbpts
!            write(24,*) k, xt(k), real(ffc(k)), real(ffc1(k))
!         enddo
  
    
         do k=1,nbpts
            xt(k) = xx1 + (k-1)*step
            ppis(k) = pis(nnbr,lnbr,s,i,wx,xx,xt(k))       ! its  correct
            gppis(k) = g_pis(nnbr,lnbr,m,wx,xx,xt(k))
   
!           ppisc(k) = ppis(k)
            ppisc(k) = gppis(k)
         enddo
   
         call differentiation_real(nbpts, xt, gppis, dx, ddppis)
         call differentiation_real(nbpts, xt, ff, dx, df)
         call fourier_diff(nbpts, krange, kk, xt, ppisc, dx, dfc)
  
! **  This is the part where we assess the derivative
         do k=1,nbpts
            aux4(k) = dg_pis(nnbr,lnbr,m,wx,xx,xt(k))
            write(ut42,*) k, s, i, xt(k), ddppis(k),                   &
                     dpis(nnbr,lnbr,s,i,wx,xx,xt(k)),                  &
                     dg_pis(nnbr,lnbr,m,wx,xx,xt(k)) , real(dfc(k)),   &
                     (dpis(nnbr,lnbr,s,i,wx,xx,xt(k))-ddppis(k)) /     &
                      dpis(nnbr,lnbr,s,i,wx,xx,xt(k))*100
         enddo

        m = (nnbr-1)*(s-1)+i
! ***  Derivative at the choke points
       if (m.ne.(lnbr*(nnbr-1)+1)) then
          do sp=1,lnbr
             do j=1,nnbr
                 mp = (sp-1)*(nnbr-1)+j
!                aux(mp)= dpis2(nnbr,lnbr,s,i,wx,xx,j)
                 auxf(mp)= g_pis(nnbr,lnbr,m,wx,xx,xx(mp))
!                auxd(mp)= dg_pis(nnbr,lnbr,m,wx,xx,xx(mp))
                 auxd(mp)= dgpis_ref(nnbr,lnbr,m,wx,xx,sp,j)
                 write(ut43,*)  s, i, j, mp, xx(mp), wx(mp),            &
!                    pis(nnbr,lnbr,s,i,wx,xx,xx(mp)),                   &
                     g_pis(nnbr,lnbr,m,wx,xx,xx(mp)),                   &
                     dpis(nnbr,lnbr,s,i,wx,xx,xx(mp)),                  &
                     dpis2(nnbr,lnbr,s,i,wx,xx,j),                      &
!                    dg_pis(nnbr,lnbr,m,wx,xx,xx(mp))
                     dgpis_ref(nnbr,lnbr,m,wx,xx,sp,j)
                 write(422,*) m, sp,j,mp,dg_pis(nnbr,lnbr,m,wx,xx,xx(mp))
              enddo
          enddo
       else
          auxf(mp)= 0.d0
          auxd(mp)= 0.d0
       end if
 
 
        return
        end





        subroutine analytical_basis_reduced(workdir, s, i,            &
                           auxc1, auxc2, auxc3, auxc4, wfc)
        use util
        include 'open/files_index'
        integer :: aux0, krange
        integer :: i, j, s, k, np, ns, n, m, mp, nmax1, nn,            &
                         sp, nmax2
        real*8 :: r, step, xx1, xx2, scaling,                         &
                             kstep, y0, kk1, kk2
        real*8, dimension(:), allocatable :: xa, xs, xx, xt, wa, wx,  &
                              kk, aux1, aux2, aux3, aux4, auxf, auxd
        real*8, dimension(:), allocatable :: xy, xtt, wy
        complex*16, dimension(:), allocatable :: auxc

        complex*16, dimension(:) :: auxc1, auxc2, auxc3, auxc4
        complex*16, dimension(:) :: wfc

        character(255) :: workdir


!       real*8 :: pis, dpis, dpis2, g_pis, dg_pis
 
        include 'param'


        nmax1 = (nnbr-1)*lnbr-1
        nmax2 = (nnbr-1)*lnbr+1
        m = (s-1)*(nnbr-1) + i
        mp = m-1

        allocate(xx(nmax1), wx(nmax1))
        allocate(xy(nmax2), xtt(nbpts), wy(nmax2))
        allocate(auxf(nmax2), auxd(nmax2))
        allocate(aux1(nmax1), aux2(nmax1))
        allocate(aux3(nbpts), aux4(nbpts))



        xx1 = -xmax/2.d0
        xx2 =  xmax/2.d0
        step = (xx2 - xx1) / (nbpts-1)
        do k=1,nbpts
           xtt(k) =  xx1 + (k-1) * step
        enddo

        write(*,*) s, i, j
           call analytical_basis(workdir,s,i,xy,wy,                   &
                          auxf,auxd,aux3,aux4)

        if ((m.gt.1).and.(m.lt.nmax2)) then
           do sp=1,lnbr
              do j=1,nnbr-1
                 if((sp.ne.lnbr).or.(j.ne.(nnbr-1))) then
                    k = (sp-1)*(nnbr-1) + j
                    xx(k)=xy(k+1)
                    aux1(k)=auxf(k+1)
                    aux2(k)=auxd(k+1)
                    write(424,*) m,mp,sp,j,k,aux3(k),aux4(k)
!                else
!                   k = (sp-1)*(nnbr-1) + j - 1
!                   xx(k)=xy(k+1)
!                   aux1(k)=0.d0
!                   aux2(k)=0.d0
                 end if
              enddo
           enddo


!          if((i.ne.1)) then 
               auxc1 = auxc1 + wfc(mp)*aux1
               auxc2 = auxc2 + wfc(mp)*aux2
               auxc3 = auxc3 + wfc(mp)*aux3
               auxc4 = auxc4 + wfc(mp)*aux4
!          end if
        end if
 
!         do k=1,nmax1
!            write(888,*) s, i, k, xx(k), aux1(k), aux2(k)
!         enddo
!
!         do k=1,nbpts
!            write(889,*) s, i, k, xtt(k), aux3(k), aux4(k)
!         enddo
!
!        do k=1,nmax2
!           write(1111,*) s, i, k, xy(k), auxf(k), auxd(k)
!        enddo
 
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
 
 
       do k=1,ns-1
          rs(k) = x(k*(np-1))
       enddo
       rs(ns) = x(ns*(np-1))
 
       call sector_mapper(np,ns,s,w,x,ww,rr)
       eps= 1.d-20
       pis= 1.d0
 
       if (r.ge.rr(1).and.(r.le.rr(np))) then
          do j=1,np
             if (i/=j) then
                pis = pis* (r-rr(j))/(rr(i)-rr(j))
             end if
          enddo
       else
             pis = 0.d0
       end if
       pis = pis/sqrt(ww(i))
 
 
       return 
       end


       real*8 function pisij(np,ns,s,i,w,x,j)
       implicit none
       integer :: np, ns, s, i, j
       integer :: k
       real*8 :: eps
       real*8, dimension(*) :: w, x
       real*8, dimension(np) :: rr, ww
       real*8, dimension(ns) :: rs
 
 
       call sector_mapper(np,ns,s,w,x,ww,rr)
 
       if (i.eq.j) then
          pisij = 1.d0
       else
          pisij = 0.d0
       end if
     
       pisij = pisij/sqrt(ww(i))
 
       return 
       end
 
 

        real*8 function g_pis(np,ns,n,w,x,r)
        implicit none
        integer :: np,ns, n, m, mp
        integer :: i, j, k, s, aux1, sp, sm
        real*8 :: aux0, xmin
        real*8, dimension(*) :: w, x
        real*8, dimension(np) :: rr, ww
        real*8, dimension(ns) :: rs
  
        real*8 :: r
!       real*8 :: pis
  
        do sp=1,ns
          if (n.le.((sp-1)*(np-1)+np)) then
             i = n - (sp-1)*(np-1)
             exit
          end if
        enddo
    
        s=sp
        if (n.eq.((s-1)*(np-1)+np)) then
  
           call sector_mapper(np,ns,s,w,x,ww,rr)       
           if (r.lt.x(n)) then
              g_pis = sqrt(ww(np)) *                                   &
                             pis(np,ns,s,np,w,x,r)/(sqrt(ww(1)+ww(np)))
           else
              g_pis = sqrt(ww(1))  *                                   &
                      pis(np,ns,s+1,1,w,x,r)/(sqrt(ww(1)+ww(np)))
           end if
        else
           g_pis = pis(np,ns,s,i,w,x,r)
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
 
!       real*8 :: pis
   
        dpis= 0.d0
        call sector_mapper(np,ns,s,w,x,ww,rr)
 
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
  
        dpis = dpis/sqrt(ww(i))
  
        return         
        end
 
 
        real*8 function dg_pis(np,ns,n,w,x,r)
        implicit none
        integer :: np, ns, n, s, i
        integer :: j, k, m, sp, sm
        real*8 :: aux, r, aux0, aux1
        real*8, dimension(*) :: w, x
        real*8, dimension(np) :: rr, ww
        real*8, dimension(np) :: rrr, www
  
!       real*8 :: dpis
!       real*8 :: pis
    
        dg_pis= 0.d0
  
  
         do sp=1,ns
           if (n.le.((sp-1)*(np-1)+np)) then
              i = n - (sp-1)*(np-1)
              exit
           end if
         enddo
  
  
         s=sp
  
         if (n.eq.((s-1)*(np-1)+np)) then
            call sector_mapper(np,ns,s+1,w,x,ww,rr)     
  
            if (s.ne.ns) then
!              if (r.ge.x((s-1)*(np-1)+1).and.r.le.x((s+1)*(np-1)+1)) then
               if (r.gt.x((s-1)*(np-1)+1).and.r.lt.x((s+1)*(np-1)+1)) then
                  if (r.lt.x(n)) then
                     dg_pis = dpis(np,ns,s,np,w,x,r)
                  else
                     dg_pis = dpis(np,ns,s+1,1,w,x,r)
                  end if
               else
                  dg_pis = 0.d0
               end if
            else
!              if (r.ge.x((s-1)*(np-1)+1).and.r.le.x((s)*(np-1)+1)) then
               if (r.gt.x((s-1)*(np-1)+1).and.r.lt.x((s)*(np-1)+1)) then
                  dg_pis = dpis(np,ns,s,np,w,x,r)
               else
                  dg_pis = 0.d0
               end if
            end if 
  
            dg_pis = sqrt(ww(1))  * dg_pis /(sqrt(ww(1)+ww(np)))
   
         else
  
!          if ((r.gt.x((s-1)*(np-1)+1)).and.(r.lt.x((s)*(np-1)+1))) then
           if ((r.ge.x((s-1)*(np-1)+1)).and.(r.le.x((s)*(np-1)+1))) then
              dg_pis = dpis(np,ns,s,i,w,x,r)
           else
              dg_pis = 0.d0 
           end if
  
         end if
  
        return         
        end
  
  
  
        real*8 function dpis2(np,ns,s,i,w,x,j)
        implicit none
        integer :: np, ns, s, i
        integer :: j, k, m
        real*8 :: aux, r
        real*8, dimension(*) :: w, x
        real*8, dimension(np) :: rr, ww
 
!       real*8 :: pis
   
        dpis2= 0.d0
        call sector_mapper(np,ns,s,w,x,ww,rr)

 
        if ((j/=i)) then
           aux = 1.d0
           do k=1,np
              if ((k/=i).and.(k/=j)) then
                 aux = aux* (rr(j)-rr(k))/(rr(i)-rr(k))
              end if
           enddo
           dpis2 =  aux / (rr(i)-rr(j))
        else
           dpis2 = 0.5d0 *                                         &
             (pisij(np,ns,s,i,w,x,np)**2 -  pisij(np,ns,s,i,w,x,1)**2)
           dpis2 = 0.0d0
        end if
 
        dpis2 = dpis2 / sqrt(ww(i))
  
        return         
        end




         real*8 function dgpis_ref(np,ns,n,w,x,sp,j)
         implicit none
         integer :: np, ns, sp, s, i
         integer :: j, k, n
         real*8 :: aux, r
         real*8, dimension(*) :: w, x
         real*8, dimension(np) :: rr, ww
  
!       real*8 :: pis
    
!         if (s.eq.sp) then
!            dgpis_ref = dpis2(np,ns,s,i,w,x,j)
!         else
!            dgpis_ref = 0.0d0 
!         end if

         do s=1,ns
           if (n.le.((s-1)*(np-1)+np)) then
              i = n - (s-1)*(np-1)
              exit
           end if
         enddo

!        s=sp

!         if (s.eq.sp) then
!            dgpis_ref = dpis2(np,ns,s,i,w,x,j)
!         else
!            dgpis_ref = 0.0d0 
!         end if

         call sector_mapper(np,ns,s,w,x,ww,rr)         ! need this for the ww

         if (n.eq.((s-1)*(np-1)+np)) then         
            dgpis_ref = 1.d0
            if((j.eq.np).or.(j.eq.1)) then
!           if((j.eq.1)) then
               dgpis_ref = 0.d0
            else
               if (((sp-s).eq.1)) then                  !upper segment np to 1
                  dgpis_ref = dpis2(np,ns,sp,1,w,x,j)
               else if ((s-sp).eq.1.and.(i.eq.1)) then  !lower segment 1 to np
                    dgpis_ref = dpis2(np,ns,s,5,w,x,j)
               else if ((s-sp).eq.0) then
                    dgpis_ref = dpis2(np,ns,sp,i,w,x,j)
               else
                 dgpis_ref = 0.d0
               end if
  
  
                dgpis_ref = sqrt(ww(1)) *                             &
                               dgpis_ref /(sqrt(ww(1)+ww(np)))

            end if

            write(656,*) i,j,n,dgpis_ref
   
         else
            if(j.eq.1) then
               dgpis_ref = 0.d0
            else
               if ((s-sp).eq.0) then
                  dgpis_ref = dpis2(np,ns,sp,i,w,x,j) 
               else
                  dgpis_ref = 0.d0
               end if
            end if
         end if

!        if ((n.ne.1).and.(n.ne.(ns*(np-1)-1))) then         
!           if (i.eq.j) then
!               dgpis_ref = 0.d0
!           end if
!        end if



!         do s=1,ns
!           if (n.le.((s-1)*(np-1)+np)) then
!              i = n - (s-1)*(np-1)
!              exit
!           end if
!         enddo
!
!         call sector_mapper(np,ns,s,w,x,ww,rr)         ! need this for the ww
!
!         if (n.eq.((s-1)*(np-1)+np)) then         
!            if ((sp-s).eq.1) then                      ! to next sect.
!               dgpis_ref = dpis2(np,ns,sp,1,w,x,j) 
!            else if ((sp-s).eq.1) then                 ! from prev sect.
!               dgpis_ref = dpis2(np,ns,s,np,w,x,j) 
!            else if ((s-sp).eq.0) then
!               dgpis_ref = dpis2(np,ns,sp,i,w,x,j) 
!            else 
!               dgpis_ref = 0.d0
!            end if
!   
!         else
!            if ((s-sp).eq.0) then
!               dgpis_ref = dpis2(np,ns,sp,i,w,x,j) 
!            else
!               dgpis_ref = 0.d0
!            end if
!         end if

         return         
         end
 
 
        subroutine sector_mapper(np,ns,s,w,x,ww,rr)
        implicit none
        integer :: np, m, ns
        integer :: j, s, k
        real*8, dimension(*) :: x, w
        real*8, dimension(np) :: rr, ww
 
 
        do k=1,s
           do j=1,np
              m = (k-1)*(np-1) + j
              rr(j)=x(m)
              ww(j)=w(m)
           enddo
        enddo
 
       return
       end





        subroutine analytical_basis_save(workdir,s,i,xx,wx,xt)
        use util
        include 'open/files_index'
        integer :: aux0, aux1, krange
        integer :: i, j, s, k, np, ns, n, m, mp, nmax, nn,            &
                         sp
        real*8 :: r, step, xx1, xx2, scaling,                         &
                             step_, kstep, y0, kk1, kk2
        real*8, dimension(:), allocatable :: xa, xs, wa, kk
        real*8, dimension(:) :: xx, xt, wx 
        real*8, dimension(:), allocatable :: ppis, gppis, ddppis,     &
                                    ff, df, dff,  dx
        complex*16, dimension(:), allocatable :: ffc, ffk, ffc1,      &
                              dfc, ppisc 

        character(len=1024) :: filename10='quad_pts_wts.dat'
        character(len=1024) :: filename20='ppi_'
        character(len=1024) :: filename30='dpi_'
        character(len=1024) :: filename40='dp2_'

        character(len=1024) :: filename1
        character(len=1024) :: filename2
        character(len=1024) :: filename3
        character(len=1024) :: filename4

        character(len=1024) :: format_string
        character(len=1024) :: indx
        character(255) :: workdir


!       real*8 :: pis, dpis, dpis2, g_pis, dg_pis
 
        include 'param'
        include 'common/open_files'
!       include 'open/files_specs'

!       include 'common/laser'
        include 'open/files_deriv_analyt'

!       call open_files_to_store(0)

        allocate(xa(nnbr),wa(nnbr), xs(lnbr))
        nbpts=200
 
!       ut=100
!       ut=200
!       ut=300
!       ut=400
 
        nmax = lnbr*(nnbr-1) + 1
        xx1 = -xmax/2.d0
        xx2 =  xmax/2.d0
 
        scaling = 0.5d0*xmax/lnbr
 
 
!       allocate(xx(nmax), xt(nbpts), wx(nmax))
 
        write(*,*) "np = ", nnbr
        write(*,*) "ns = ", lnbr
        write(*,*) "nmax = ", nmax
        write(*,*)
 
        call lobatto_compute( nnbr, xa, wa )
 
!       do j=1,nnbr
!          write(*,*) nnbr, j, xa(j), wa(j)
!       enddo
!       write(333,*) 
!       write(*,*) "do i get here?", nnbr, lnbr, nmax
  
 
        do sp=1,lnbr
           do j=1,nnbr
              m = (nnbr-1)*(sp-1)+j
              xx(m) = xa(j) + (sp-.5d0-.5d0*lnbr)*2
              xx(m) = xx(m)*scaling
              if(j.eq.(nnbr)) then
                 wx(m) = dsqrt(wa(nnbr)+wa(1))
              else
                 wx(m) = dsqrt(wa(j))
              end if 
!             write(*,*) m, xx(m), wx(m)
           enddo
!          write(*,*)
        enddo
 
        allocate(ppis(nbpts), gppis(nbpts), ppisc(nbpts))
        step = (xx2 - xx1) / (nbpts-1)
 
        allocate(ff(nbpts), df(nbpts), dff(nbpts),                   &
                                         ddppis(nbpts),dx(nbpts))

      write(ut40,*) "phase 1 : Find the choke points", m
      do sp=1,lnbr
         do j=1,nnbr
            m = (sp-1)*(nnbr-1)+j
            write(ut40,*)   m, sp, j, xx(m), wx(m),                    &    
                   pis(nnbr,lnbr,sp,j,wx,xx,xx(m)),                    & 
                   g_pis(nnbr,lnbr,m,wx,xx,xx(m)), 0.d0
!           write(*,*)     m, sp, j, xx(m), wx(m), pis(np,ns,sp,j,wx,xx,xx(m)),
!           g_pis(np,ns,m,wx,xx,xx(m))
         enddo
         write(ut40,*)
      enddo


!      write(*,*) "phase 2 : Lets model the sector changes"
!      sp=1
!      do m=1,nmax-1
!         aux0 = m-sp*(np-1)-1
!         aux1 = mod(m,sp*np)
!         j = -(sp-1)*(np-1)+m
!         write(*,*) m, sp, aux0, j
!         if (aux0.ge.-1) then
!            sp=sp+1
!            write(*,*)
!         end if
!      enddo

 
!        write(indx, '(I0,"_",I0)') s, i
!        indx = trim(indx)//".dat"
!!       print *, trim(indx)
! 
!        filename2 = trim(filename2)//indx
!        filename3 = trim(filename3)//indx
!        filename4 = trim(filename4)//indx
!
!        open(ut41,file=trim(workdir)//trim(filename2),status='unknown')
!        open(ut42,file=trim(workdir)//trim(filename3),status='unknown')
!        open(ut43,file=trim(workdir)//trim(filename4),status='unknown')

        m=(s-1)*(nnbr-1) + i
        do k=1,nbpts
           xt(k) = xx1 + (k-1)*step
           ppis(k) = pis(nnbr,lnbr,s,i,wx,xx,xt(k))                
           write(ut41,*) k, lnbr, nnbr, s, i, xt(k), ppis(k),              &
                     g_pis (nnbr,lnbr,m,wx,xx,xt(k))
!          write(*,*)   k, lnbr, nnbr, s, i, xt(k), ppis(k),              &
!                    g_pis (nnbr,lnbr,m,wx,xx,xt(k))
        enddo
 
 
! ***  set dvr function
        krange=20000
        allocate(kk(krange), ffc(nbpts), ffk(krange), ffc1(nbpts),   &
                       dfc(nbpts)) 
 
!       ppi = 4.d0*datan(1.d0)
        kk1=-2.d0*ppi
        kk2=-kk1
        kstep = (kk2-kk1)/(krange-1)
        do k=1,krange
           kk(k) = kk1 + (k-1)*kstep
        enddo
        do k=1,nbpts
           xt(k) = xx1 + (k-1)*step
           ffc(k) = 1.d0
        enddo
        call fourier_trsfm(nbpts, krange, xt, kk, ffc, ffk)
        call inv_fourier_trsfm(krange, nbpts, kk, xt, ffk, ffc1)
 
 
         do k=1,krange
            write(23,*) k, kk(k), real(ffk(k))
         enddo
         do k=1,nbpts
            write(24,*) k, xt(k), real(ffc(k)), real(ffc1(k))
         enddo
 
   
         do k=1,nbpts
            xt(k) = xx1 + (k-1)*step
            ppis(k) = pis(nnbr,lnbr,s,i,wx,xx,xt(k))       ! its  correct
            gppis(k) = g_pis(nnbr,lnbr,m,wx,xx,xt(k))
  
!          ppisc(k) = ppis(k)
            ppisc(k) = gppis(k)
         enddo
  
         call differentiation_real(nbpts, xt, gppis, dx, ddppis)
         call differentiation_real(nbpts, xt, ff, dx, df)
         call fourier_diff(nbpts, krange, kk, xt, ppisc, dx, dfc)
 
! **  This is the part where we assess the derivative
         do k=1,nbpts
            write(ut42,*) k, s, i, xt(k), ddppis(k),                   &
                     dpis(nnbr,lnbr,s,i,wx,xx,xt(k)),                  &
                     dg_pis(nnbr,lnbr,m,wx,xx,xt(k)) , real(dfc(k)),   &
                     (dpis(nnbr,lnbr,s,i,wx,xx,xt(k))-ddppis(k)) /     &
                      dpis(nnbr,lnbr,s,i,wx,xx,xt(k))*100
         enddo
 
! ***  Derivative at the choke points
        do sp=1,lnbr
           do j=1,nnbr
              mp = (sp-1)*(nnbr-1)+j
              write(ut43,*)  s, i, j, mp, xx(mp), wx(mp),              &
                  pis(nnbr,lnbr,s,i,wx,xx,xx(mp)),                     &
                  dpis(nnbr,lnbr,s,i,wx,xx,xx(mp)),                    &
                  dpis2(nnbr,lnbr,s,i,wx,xx,j),                        &
                  dgpis_ref(nnbr,lnbr,m,wx,xx,sp,j)
           enddo
        enddo
 
 
        return
        end


      end module
