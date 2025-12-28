      program fedvr_build
      use util
      use field
      use propagators
      use fedvr_conf_struct
      implicit none
      integer :: nmax, krange, lwork, info, status,                   &
                  ut1, ut2, ut3, ut4, ut5, ut6, ut7, ut8, ut9, ut10,  &
                  ut11, ut12, ut13, ut14, ut15,                       &
                  i, j, k, ij, p, q, r

      real*8 :: t_ini, t_end, emax, tau,                              &
                    step, eft, t,                                     &
                    scaling, abstol, pemd, prob, p_0, p_ion, dk,      &
                    start, finish, lap, aux, aux1, aux2,              &
                    omg_step

      complex*16 :: ci, c0, c1, cnum, a0, cprobk, auxcc

      real*8, allocatable, dimension(:) :: xa, wa, xx, wx,            &
                                   vec_matup, eigval,                 &
                                   kk, init,                          & 
                                   auxr1, auxr2, auxr3,               &
                                   time_, omg,                        &
                                   p_k, dx, dt,                       &
                                   qvc1, qvc2, qvc3, qvc4,            &
                                   matint_v, vextv,                   &
                                   ipiv, work, ifail, dnrm2
  
      complex*16, allocatable, dimension(:) :: wf0, wfc0,             &
                                               wf, wfc, wfc_, dwfc,   &
                                               init2, phi0, phic0,    &
                                               phi, phic,             &
                                               rr_, ak, cprob2, auxc, &
                                               auxc1, auxc2, auxc3,   &
                                               pk0, b0w_1, b0w_2,     &
                                               d_t, dd_t, d_w, dd_w,  &
                                               x_t,                   &
                                               p_t, pp_t, p_w, pp_w,  &
                                               dp_t, ft

      real*8, allocatable, dimension(:,:) :: lu, id, inv,             &
                                       pkk, bkw_1, bkw_2,             &
                                       eigenvec, ss

      complex*16, allocatable, dimension(:,:) :: in_states,           &
                                                 out_states,          &
                                                 cwf

      include 'param'

      common/pulse_par/emax,tau,t_ini,t_end
      common/files/ut1,ut2,ut3,ut4,ut5,ut6,ut7,ut8,ut9,ut10,          &
                   ut11,ut12,ut13,ut14,ut15


      namelist / param / nnbr,lnbr,xmax,nsteps
      nmax = (nnbr-1)*lnbr-1
      lwork=4*nmax
  
      call cpu_time(start)
!     write(6,nml=param) 
      print*, "&PARAMETERS"
      print*, "N_TAU = ", n_tau
      print*, "NNBR = ", nnbr
      print*, "LNBR = ", lnbr
      print*, "XMAX = ", xmax
      print*, "NSTEPS = ", nsteps
!     print*, "INTENSITY  = ", peaki
      print*, "FIELD STRENGTH =", dsqrt(peaki/peaki0)
      print*, "OMEGA  = ", omega
      print*, "KMAX  = ", kmax
      print*, "OMG_MAX  = ", omg_max
      write(6,*) nmax 
      write(6,*) 



      allocate(xa(nmax),wa(nmax),xx(nmax),wx(nmax))
      allocate(eigval(nmax),init(nmax),init2(nmax),                   &
                     wf0(nmax),wfc0(nmax),                            &
                     phic0(nmax),phi0(nmax),                          &
                     wf(nmax),wfc(nmax),wfc_(nmax),                   &
                     auxc(nmax))
      allocate(eigenvec(nmax,nmax))
 
      allocate(matint_v(nmax),vextv(nmax))
      allocate(time_(nsteps+1), ft(nsteps+1))
      allocate(d_t(nsteps+1), x_t(nsteps+1),                          &
                        dd_t(nsteps), d_w(n_omg+1), dd_w(n_omg+1),    &                                 dt(nsteps), p_t(nsteps+1), dp_t(nsteps+1),    &
                        pp_t(nsteps+1), p_w(n_omg+1), pp_w(n_omg+1))

      c0=dcmplx(0.d0,0.d0)
      ci=dcmplx(0.d0,1.d0)
      c1=dcmplx(1.d0,0.d0)
 
      call open_files
      call pulse_chr(emax,tau,t_ini,t_end)
      call plot_field
 
 
      scaling = 0.5d0*xmax/lnbr
!!      call lobatto_compute(nnbr,xa,wa) 
!!
!!      call fedvr_hamilton_conf(scaling,xa,wa,xx,wx,eigval,            &
!!                         eigenvec)
!!      call cpu_time(lap)
!!      write(*,*) "cpu time after diagonalization", lap-start
!!
!! 
!!      do i=1,nmax
!!         j=1
!!         write(ut3,*) i, eigval(i), dsqrt(2.0d0*eigval(i))
!!         write(ut4,*) xx(i),eigenvec(i,j)/wx(i)/dsqrt(scaling),       &
!!             eigenvec(i,j+1)/wx(i)/dsqrt(scaling),                    &
!!             eigenvec(i,j+2)/wx(i)/dsqrt(scaling),                    &
!!             eigenvec(i,j+3)/wx(i)/dsqrt(scaling),                    &
!!             eigenvec(i,j+4)/wx(i)/dsqrt(scaling)                     
!!      enddo
!!
!!        
!!!     allocate(inv(nmax,nmax))
!!!     inv=eigenvec
!!! 
!!!     allocate(ipiv(nmax),work(lwork))
!!!     call dgetrf(nmax,nmax,inv,nmax,ipiv,info) 
!!!     call dgetri(nmax,inv,nmax,ipiv,work,lwork,info) 
!!!     deallocate(ipiv,work)
!!  
!!!     id = matmul(inv,eigenvec) 
!!!     call dgemm('n','n',nmax,nmax,nmax,1.d0,inv,nmax,                &
!!!            eigenvec,nmax,0.d0,id,nmax) 
!!!     stop 
!!   
!!! *** Preparing the initial state
!!      wfc0 = dcmplx(eigenvec(:,1))/wx/dsqrt(scaling)
!!      wf0 = wfc0 * wx * dsqrt(scaling)
!!
!!      auxc=wf0
!!!     wf0 = matmul(inv,auxc)
!!!     deallocate(inv)
!!      wf0 = matmul(transpose(eigenvec),auxc)
!!
!!!     wf0=c0
!!!     wf0(1)=c1
!! 
!!      do i=1,nmax
!!!        last column is the analytical expression
!!         init(i)=dsqrt(vkpa)*exp(-vkpa*abs(xx(i)))
!!         init2(i)=-(ci/omega)*(dsqrt(vkpa))**3 * sgn(xx(i))            &    
!!                  * (exp(-vkpa*abs(xx(i)))                             &        
!!                         -exp(-dsqrt(vkpa**2+2*omega)*abs(xx(i))))
!!         write(ut5,*) xx(i), real(wf0(i)),                            &
!!                      eigenvec(i,1), real(wfc0(i)),                   &
!!                      init(i), imag(init2(i))
!!      enddo
!!
!!
!!
!!      matint_v=0.0d0 
!!      matint_v=xx
!!      eft=0.0d0
!!   
!!      step=(t_end-t_ini)/dble(nsteps)
!!      time_(1)=t_ini 
!!
!!      allocate(omg(n_omg+1))
!!      omg_step = (omg_max - omg_min)/dble(n_omg) 
!!      
!!      do k=1,n_omg+1
!!         omg(k)=omg_min + (k-1)*omg_step
!!      enddo
!!       
!!      allocate(dx(nmax),dwfc(nmax))
!!      allocate(ss(nmax,nsteps+1))
!!
!!!     wfc = wfc0
!!
!!      allocate(auxc1(nsteps+1),auxc2(nsteps+1),auxc3(nsteps+1))
!!
!!!     pp_t = c0
!!      do i=0,nsteps
!!         
!!         ij=i+1
!!         
!!         t=time_(ij)
!!
!!
!!         auxc = - conjg(wfc) * xx * wfc * wx*wx*scaling
!!         do j=1,nmax
!!            d_t(ij) = d_t(ij) + auxc(j)
!!            x_t(ij) = x_t(ij) - auxc(j)
!!         enddo
!!
!!         call differentiation(nmax, xx, wfc, dx, dwfc)
!!         dwfc = -ci*dwfc
!!         auxc =  conjg(wfc) * dwfc * wx*wx*scaling
!!
!!!        do j=1,nmax
!!!           auxc(j) =  conjg(wfc(j)) * dwfc(j) * wx(j)*wx(j)*scaling
!!!        enddo
!!
!!         do j=1,nmax
!!            p_t(ij) = p_t(ij) + auxc(j)
!!         enddo
!!
!!
!!         eft=efield(t)
!!         auxc = conjg(wfc) * eft * wfc * wx*wx*scaling
!!         do j=1,nmax
!!            ft(ij) = ft(ij) + auxc(j)
!!         enddo
!!
!!         do k=1,ij
!!            auxc1(k) = ft(k)
!!         enddo
!!         call composite_simpson_18c(ij, time_, auxc1, auxcc)
!!         pp_t(ij) = auxcc
!!!        write(122,*) ij, time_(ij), p_t(ij), real(p_t(ij)),          &
!!!                      pp_t(ij), real(pp_t(ij))
!!
!!!        write(ut13,*) ij, time_(ij), real(d_t(ij)), real(p_t(ij)),   &
!!!                                     d_t(ij), p_t(ij), ft(ij)
!! 
!!         vextv = eft*matint_v 
!!
!!         call split_op(i,step,nmax,vextv,eigenvec,eigval,             &
!!                        wf0,wf)
!!
!!! ***    To get working for the nontrivial case
!!!        vextv = 0.0d0
!!!        nonlin = eft*matint
!!!        call ifrk4(i,step,nmax,nonlin,vextv,eigenvec,eigval,         &
!!!                           wf0,wf)
!!
!!         wfc_=matmul(eigenvec,wf)
!!         wfc=wfc_/wx/dsqrt(scaling)
!!
!!!        wfc=matmul(eigenvec,wf)
!!!        wfc=wfc/wx/dsqrt(scaling)
!!
!!
!!!!       data storage
!!!        if (mod(i,store_val).eq.0) then
!!!           wfc_=matmul(eigenvec,wf)
!!!           wfc=wfc_/wx/dsqrt(scaling)
!!!           do j=1,nmax
!!!              write(ut6,*) time_(ij),j,eigval(j),                    &
!!!                                real(wf(j)),imag(wf(j))
!!!              write(ut7,*) time_(ij),xx(j),real(wfc(j)),imag(wfc(j))
!!!              aux1 = 0.0d0
!!!              aux2 = 0.0d0
!!
!!!              if (xx(j).eq.0d0) then
!!!                 write(ut8,*) time_(ij),xx(j), &
!!!                      real(wfc(j)), imag(wfc(j))
!!!              end if
!!!           enddo
!!!           do k=1,nmax
!!!              aux1 = aux1 + abs(wf(k))**2
!!!              aux2 = aux2 + abs(wfc(k))**2
!!!           enddo
!!!           aux1 = dsqrt(aux1)
!!!           aux2 = dsqrt(aux2)
!!!           write(ut13,*) time_(ij),aux1, aux2
!!!        end if
!!
!!         p_ion=0.0d0
!!         do j=2,nmax
!!            p_ion = p_ion + abs(wf(j))**2
!!         enddo
!!         write(ut10,*) time_(ij), abs(wf(1))**2, p_ion
!! 
!!
!!         time_(ij+1)=time_(ij)+step 
!!      enddo
!!
!!
!!
!!! *** Data storage at the end of the pulse
!!      wfc_=matmul(eigenvec,wf)
!!      wfc=wfc_/wx/dsqrt(scaling)
!!
!!!     wfc=matmul(eigenvec,wf)
!!!     wfc=wfc/wx/dsqrt(scaling)
!!
!!      do j=1,nmax
!!         write(ut6,*) time_(ij),j,xx(j),dsqrt(2.d0*eigval(j)),       &
!!            real(wf(j)),imag(wf(j))
!!         write(ut7,*) time_(ij), xx(j), real(wfc(j)), imag(wfc(j))
!!         if (xx(j).eq.0d0) then
!!            write(ut8,*) time_(ij), xx(j), real(wfc(j)), imag(wfc(j))
!!         end if
!!      enddo
!!
!!      deallocate(eigenvec)


! *** Post processing
! *** Find the first excited states...
      p=0
      do p=1,nmax
         if (eigval(p).gt.0.0d0) then
            exit
         end if
      enddo
!     print *, p
!     stop

      krange = 4d0*nmax
      dk = kmax/(krange/2-1)

      allocate(rr_(krange), kk(krange))
      allocate(in_states(nmax,krange))
!     allocate(out_states(nmax,krange))
      kk = 0.0d0
      rr_ = c0

      in_states = c0
!     out_states = c0


      write(*,*) krange
      do j=1,krange
         if (j.le.(krange/2-p+1)) then  
            ij = krange/2 + 1 -j
            k = - ij
            kk(j)= -kmax + (j-1)*dk              ! evenly spaced distrib
            
!           write(*,*) j, k, kk((j)),                                 &
!             dsqrt(2.d0*eigval(krange/2+1-j)), eigval(krange/2+1-j)  
         else if (j.ge.(krange/2+p)) then  
            ij = j -krange/2
            k = ij
            kk(j) = (ij-1)*dk                    ! evenly spaced distrib

!           write(*,*) j, k, kk((j)),                                 &
!             dsqrt(2.d0*eigval(j-krange/2)), eigval(j-krange/2)      
         else
            k=0

!           write(*,*) j, k, kk(j), 0.0d0
         end if

         rr_(j) = ci*vkpa/(-abs(kk(j)) - ci*vkpa)
!        write(700,*) j, kk(j), rr_(j), real(rr_(j)), imag(rr_(j))
         do i=1,nmax
            in_states(i,j) = exp(ci*kk(j)*xx(i))                      &
                 +    rr_(j) * exp( -ci*abs(kk(j)*xx(i)))
         enddo

!         rr_(j) = ci*vkpa/(abs(kk(j)) - ci*vkpa)
!!        write(800,*) j, kk(j), rr_(j), real(rr_(j)), imag(rr_(j))
!         do i=1,nmax
!            out_states(i,j) = exp(ci*kk(j)*xx(i))                     &
!                 +    rr_(j) * exp( ci*abs(kk(j)*xx(i)))
!         enddo
!        write(*,*) j, kk(j)
       
      enddo
 
!     do k=1,krange
!        write(*,*) k, rr_(k)
!     enddo

      deallocate(rr_)


! *** write the states From -k2 to k2
      do i=1,nmax
         write(ut11,*) i, xx(i),                                      &
            real(in_states(i,krange/2-p)),                            &
            imag(in_states(i,krange/2-p)),                            &
   
            real(in_states(i,krange/2-p+1)),                          &
            imag(in_states(i,krange/2-p+1)),                          &   

            real(in_states(i,krange/2-p+2)),                          &
            imag(in_states(i,krange/2-p+2)),                          &  !bs 

            real(in_states(i,krange/2-p+3)),                          &
            imag(in_states(i,krange/2-p+3)),                          &  !bs

            real(in_states(i,krange/2-p+4)),                          &
            imag(in_states(i,krange/2-p+4)),                          &         

            real(in_states(i,krange/2-p+5)),                          &
            imag(in_states(i,krange/2-p+5))                           


!        write(ut12,*) i, xx(i),                                      &
!           real(out_states(i,krange/2-p)),                           &
!           imag(out_states(i,krange/2-p)),                           &
!  
!           real(out_states(i,krange/2-p+1)),                         &
!           imag(out_states(i,krange/2-p+1)),                         &   

!           real(out_states(i,krange/2-p+2)),                         &
!           imag(out_states(i,krange/2-p+2)),                         &  !bs 

!           real(out_states(i,krange/2-p+3)),                         &
!           imag(out_states(i,krange/2-p+3)),                         &  !bs

!           real(out_states(i,krange/2-p+4)),                         &
!           imag(out_states(i,krange/2-p+4)),                         &
!        
!           real(out_states(i,krange/2-p+5)),                         &
!           imag(out_states(i,krange/2-p+5))                          
      enddo
 


! *** Compute the PEMD by means of the theoretical scattering states
!  -- Quadrature

      a0=c0
      auxc = init * wfc * wx*wx*scaling

      do i=1,nmax
         a0 = auxc(i) + a0
      enddo
      a0 = exp(ci*eigval(1)*time_(nsteps+1))*a0
      write(*,*) "a0 with init", (abs(a0))**2


      allocate(ak(krange),pk0(krange))
      ak=c0


      call differentiation(nmax, xx, wfc0, dx, dwfc)
      dwfc = -ci*dwfc

      do k=1,krange

!        auxc = conjg(in_states(:,k)) * wfc_ *wx*dsqrt(scaling)  
         auxc = conjg(in_states(:,k)) * wfc *wx*wx*scaling  

!        if (k.lt.nmax) then
!           if(xx(i).lt.0.d0) then
!              auxc = conjg(out_states(:,k)) * wfc_ *wx*dsqrt(scaling)  
!           else
!              auxc = conjg(out_states(:,k)) * wfc_ *wx*dsqrt(scaling)  
!           endif
!        else
!           if(xx(i).le.0.d0) then
!              auxc = conjg(in_states(:,k)) * wfc_ *wx*dsqrt(scaling)  
!           else
!              auxc = conjg(in_states(:,k)) * wfc_ *wx*dsqrt(scaling)  
!           end if
!        endif

         do j=1,nmax
            ak(k) = auxc(j) + ak(k)
         enddo

         ak(k) = exp(ci*(0.5d0*kk(k)**2)*time_(nsteps+1))*ak(k)

         auxc =  conjg(in_states(:,k)) * dwfc * wx*wx*scaling
         do j=1,nmax
            pk0(k) = auxc(j) + pk0(k)
         enddo

         write(ut9,*) k, kk(k), real(ak(k)), imag(ak(k)),             &
                0.5d0*kk(k)**2
!        write(1001,*) k, kk(k), real(pk0(k)),                        &
!                    (2.d0*kk(k)*vkpa**(3.d0/2))/(vkpa**2 + kk(k)**2)
      enddo

      deallocate(in_states)
!     deallocate(out_states)

      allocate(p_k(krange), auxr1(krange/2), auxr2(krange/2))
      p_k = (abs(ak))**2
      p_ion = 0.0d0


      auxr1=0.0d0 
      auxr2=0.0d0 
      aux1=0.0d0 

      do k=1,krange/2
        auxr1(k) = kk(k)
        auxr2(k) = p_k(k)
      enddo
      call composite_simpson_18(krange/2, auxr1, auxr2, aux1)
      
      do k=1,krange/2
        auxr1(k) = kk(k+krange/2)
        auxr2(k) = p_k(k+krange/2)
      enddo
      call composite_simpson_18(krange/2, auxr1, auxr2, p_ion)

      p_ion = p_ion + aux1
      p_ion = p_ion/(2.d0*ppi)
      p_0 = (abs(a0))**2


      deallocate(kk)
      deallocate(auxr1,auxr2)
      allocate(qvc1(n_omg+1),qvc2(n_omg+1),qvc3(n_omg+1))
      allocate(qvc4(n_omg+1))

      call differentiation(nsteps+1, time_, d_t, dt, dd_t)
      call differentiation(nsteps+1, time_, p_t, dt, dp_t)

!          do k=1,nsteps+1
!!            auxr1(k) = time_(k)
!             auxc1(k) = ft(k)
!          enddo
!          call composite_simpson_18c(nsteps+1,time_,auxc1,auxcc,auxc2)
!          do ij=1,nsteps+1 
!             write(125,*) ij, time_(ij), real(p_t(ij)), real(auxc2(ij))
!             write(126,*) ij, time_(ij), real(ft(ij)), real(dp_t(ij))
!             write(127,*) ij, time_(ij), real(d_t(ij)), real(dd_t(ij))
!          enddo


      do ij=1,nsteps+1
         write(ut13,*) ij, time_(ij), real(d_t(ij)), real(dd_t(ij)),   &
                                real(p_t(ij)), real(pp_t(ij)),         &
                                d_t(ij), p_t(ij), dd_t(ij), pp_t(ij),  &
                                ft(ij)
      enddo

      do i=1,n_omg+1
         ! Eq. 80
         auxc1 = exp(ci*omg(i)*time_) * d_t 
         call composite_simpson_18c(nsteps+1, time_, auxc1, auxcc)
         d_w(i) = auxcc
         qvc1(i) =  omg(i)**2 * abs(d_w(i))**2

         ! Eq. 79
         auxc1 = exp(ci*omg(i)*time_) * p_t
         call composite_simpson_18c(nsteps+1, time_, auxc1, auxcc)
         p_w(i) = auxcc
         qvc2(i) =  abs(p_w(i))**2

         auxc1 = -exp(ci*omg(i)*time_) * dd_t
         call composite_simpson_18c(nsteps+1, time_, auxc1, auxcc)
         dd_w(i) = auxcc
         qvc3(i) = abs(dd_w(i))**2

         auxc2 = exp(ci*omg(i)*time_) * x_t 
         call composite_simpson_18c(nsteps+1, time_, auxc2, auxcc)
         pp_w(i) = exp(ci*omg(i)*time_(nsteps+1))*x_t(nsteps+1)       &
                 - exp(ci*omg(i)*time_(1))*x_t(1) - ci*omg(i)*auxcc
         qvc4(i) =  abs(pp_w(i))**2


         write(ut14,*) i, omg(i), d_w(i), dd_w(i), p_w(i), pp_w(i),   &
                                      real(d_w(i)), imag(d_w(i)),     &
                                      real(p_w(i)), imag(p_w(i)),     &
                                      real(dd_w(i)), imag(dd_w(i)),   &
                                      real(pp_w(i)), imag(pp_w(i))

         write(ut15,*) i, omg(i), qvc1(i), qvc2(i), qvc3(i), qvc4(i)
      enddo

     


      write(*,*) "surv. probability", p_0, abs(wf(1))**2
      write(*,*) "ioni. probability", 1.0d0-p_0, 1.d0-abs(wf(1))**2
      write(*,*) "comp. simpson p_ion", p_ion*(2*ppi), p_ion
      write(*,*) "time step", step
      write(*,*) "k step", dk
      write(*,*) "omega step", omg_step
       
      write(*,*) "p_0 + p_ion", p_0 + p_ion


      deallocate(ak)
      deallocate(xx,wx)
      call cpu_time(finish)
      write(*,*) "total cpu time", finish-start

      end 
