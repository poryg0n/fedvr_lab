      program fedvr_build
      use util
      use omp_lib
      use math_util
      use field
      use propagators
      use fedvr_conf_struct
      implicit none
      include 'open/files_index'
      integer :: lwork, info, store_val,                              &
                  i, j, k, ij, p, q, r,                               &
                  jn, ijk,                                            &
                  nmax, krange

      real*8 :: t_ini, t_end, emax, tau,                              &
                    step, eft, t,                                     &
                    scaling, abstol, pemd, prob, p_0, p_ion, dk,      &
                    start, finish, lap, aux, aux1, aux2,              &
                    omg_, omg_step

      complex*16 :: cnum, a0, cprobk, auxcc, auxcc1, auxcc2

      real*8, allocatable, dimension(:) :: xa, wa, xx, wx,            &
                                   vec_matup, eigval,                 &
                                   kk,                                & 
                                   auxr1, auxr2, auxr3,               &
                                   time_, omg,                        &
                                   p_k, dx, dt,                       &
                                   qvc1, qvc2, qvc3, qvc4,            &
                                   matint_v, vextv, vext_,            &
                                   ipiv, work, ifail, dnrm2
  
      complex*16, allocatable, dimension(:) :: wf0, wfc0,             &
                                               wf, wfc, wfc_,         &
                                               wwfc,                  &
                                               dwfc, dwfc2, fftpsi,   &
                                               wfck, dwfck,           &
                                               init1, init2,          &
                                               phi0, phic0,           &
                                               phi, phic,             &
                                               auxc, auxck1, auxck2,  &
                                               rr_, ak, cprob2,       &
                                               auxc1, auxc2,          &
                                               auxc3, auxc4,          &
                                               pk0, b0w_1, b0w_2,     &
                                               d_t, dd_t, d_w, dd_w,  &
                                               x_t,                   &
                                               p_t, pp_t, p_w, pp_w,  &
                                               dp_t, dx_t, ddx_t,     &  
                                               dp_t_, vextvc, dvextv, &
                                               ss, t1, t2, ft, intft

      real*8, allocatable, dimension(:,:) :: lu, id, inv,             &
                                       pkk, bkw_1, bkw_2,             &
                                       eigenvec

      complex*16, allocatable, dimension(:,:) :: in_states,           &
                                                 out_states,          &
                                                 cwf

      include 'param'
      include 'common/open_files'
      include 'common/laser'


!     namelist / param / nnbr,lnbr,xmax,nsteps
      nmax = (nnbr-1)*lnbr-1
      krange = 2.d0*nmax

      lwork=4*nmax
 
      omg_ = 1.d0

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
      print*, "OMG_  = ", omg_
      print*, "KMAX  = ", kmax
      print*, "OMG_MAX  = ", omg_max
      write(6,*) nmax 
      write(6,*) 


      call open_files_to_store(0)


      allocate(xa(nmax),wa(nmax),xx(nmax),wx(nmax))
      allocate(eigval(nmax),init1(nmax),init2(nmax),                  &
                     wf0(nmax),wfc0(nmax),                            &
                     phic0(nmax),phi0(nmax),                          &
                     phic(nmax),phi(nmax),                            &
                     wf(nmax),wfc(nmax),wfc_(nmax),                   &
                     wwfc(nmax),auxc(nmax))
      allocate(eigenvec(nmax,nmax))
 
      allocate(matint_v(nmax),vextv(nmax), vextvc(nmax),ss(nmax))
      allocate(vext_(nmax), dvextv(nmax))
      allocate(time_(nsteps+1), ft(nsteps+1), intft(nsteps+1))
      allocate(kk(krange), auxck1(krange), auxck2(krange))
      allocate(d_t(nsteps+1),x_t(nsteps+1),p_t(nsteps+1),              &
                        dx_t(nsteps+1), dt(nsteps+1),                  &
                        dd_t(nsteps+1),                                &
                        ddx_t(nsteps+1), dp_t(nsteps+1),               &
                        dp_t_(nsteps+1), pp_t(nsteps+1),               &
                        d_w(n_omg+1), dd_w(n_omg+1),                   &
                        p_w(n_omg+1), pp_w(n_omg+1))
      allocate(t1(nsteps+1),t2(nsteps+1))

 
      call pulse_chr(emax,tau,t_ini,t_end)
      call plot_field
 
 
      scaling = 0.5d0*xmax/lnbr
      call lobatto_compute(nnbr,xa,wa) 

      call fedvr_hamilton_conf(scaling,xa,wa,xx,wx,eigval,            &
                         eigenvec)
      call cpu_time(lap)
      write(*,*) "cpu time after diagonalization", lap-start

 
      do i=1,nmax
         j=1
         write(ut11,*) i, eigval(i), dsqrt(2.0d0*eigval(i))
         write(ut12,*) xx(i),eigenvec(i,j)/wx(i)/dsqrt(scaling),      &
             eigenvec(i,j+1)/wx(i)/dsqrt(scaling),                    &
             eigenvec(i,j+2)/wx(i)/dsqrt(scaling),                    &
             eigenvec(i,j+3)/wx(i)/dsqrt(scaling),                    &
             eigenvec(i,j+4)/wx(i)/dsqrt(scaling)                     
      enddo

        
!     allocate(inv(nmax,nmax))
!     inv=eigenvec
! 
!     allocate(ipiv(nmax),work(lwork))
!     call dgetrf(nmax,nmax,inv,nmax,ipiv,info) 
!     call dgetri(nmax,inv,nmax,ipiv,work,lwork,info) 
!     deallocate(ipiv,work)
  
!     id = matmul(inv,eigenvec) 
!     call dgemm('n','n',nmax,nmax,nmax,1.d0,inv,nmax,                &
!            eigenvec,nmax,0.d0,id,nmax) 
!     stop 
   
! *** Preparing the initial state
      do i=1,nmax
!        set up the initial conditions analytically
         init1(i)=dsqrt(vkpa)*exp(-vkpa*abs(xx(i)))
         init2(i)=-(ci/omg_)*(dsqrt(vkpa))**3 * sgn(xx(i))           &    
                  * (exp(-vkpa*abs(xx(i)))                            &        
                         -exp(-dsqrt(vkpa**2+2*omg_)*abs(xx(i))))
      enddo

!     wfc0 = eigenvec(:,1)/wx/dsqrt(scaling)
      wfc0 = init1
      phic0 = init2
      wf0 = wfc0 * wx * dsqrt(scaling)
      phi0 = phic0 * wx * dsqrt(scaling)
      auxc=wf0
      wf0 = matmul(transpose(eigenvec),auxc)
      auxc=phi0
      phi0 = matmul(transpose(eigenvec),auxc)

      allocate(dx(nmax),dwfc(nmax), dwfc2(nmax), wfck(krange))
      call differentiation(nmax, xx, wfc0, dx, dwfc)
      dwfc = -ci * dwfc


      do k=1,krange
         if (k.le.(krange/2-p+1)) then
            ij = krange/2 + 1 -k
            j = - ij
            kk(k)= -kmax + (k-1)*dk              ! evenly spaced distrib
         else if (k.ge.(krange/2+p)) then
            ij = k -krange/2
            j = ij
            kk(k) = (ij-1)*dk                    ! evenly spaced distrib
         end if
      enddo

! ***  fourier transform differentiation
      auxck1 = c0
      do k=1,krange
          auxc = wfc0 * exp(-ci*kk(k)*xx) * wx * wx * scaling
          do i=1,nmax
             auxck1(k) = auxck1(k) + auxc(i)
          enddo
      enddo 

      do i=1,nmax
         auxck2 = ci * kk * auxck1 * exp(ci*kk*xx(i))
         call composite_simpson_18c(krange, kk, auxck2, auxcc1)
         dwfc2(i) = auxcc1/(2.d0*ppi)
      enddo


      do i=1,nmax
         write(ut10,*) xx(i), wx(i), wf0(i), phi0(i),                 &
                      wfc0(i), phic0(i), dwfc(i),                     &
                      real(wfc0(i)), imag(phic0(i)), imag(dwfc(i))              
         write(14,*) xx(i), wx(i), dwfc(i), dwfc2(i)
      enddo


      matint_v=0.0d0 
      matint_v=xx
      eft=0.0d0
    
      step=(t_end-t_ini)/dble(nsteps)
      store_val = nsteps/nb_points
      time_(1)=t_ini 
 
      allocate(omg(n_omg+1))
      omg_step = (omg_max - omg_min)/dble(n_omg) 
      
      do k=1,n_omg+1
         omg(k)=omg_min + (k-1)*omg_step
      enddo
       
      dk = kmax/(krange/2-1)
 
      allocate(auxc1(nmax),auxc2(nmax),auxc3(nmax),auxc4(nmax))
      wfc = wfc0
      phic = phic0

!     pp_t = c0


      vext_ = 0.d0 
      write(*,*) "propagating" 

      do i=0,nsteps
         ij=i+1
         t=time_(ij)

         eft=efield(t)
         vextv = -eft*matint_v
         auxc1 = wfc

         ft(ij) = eft 
         vextvc = vextv
         call differentiation(nmax, xx, vextvc, dx, dvextv)
         auxc = conjg(wfc) * dvextv * wfc * wx*wx*scaling

         dp_t_(ij) = c0
         do j=1,nmax
            dp_t_(ij) = dp_t_(ij) + auxc(j)
         enddo

         auxc = - conjg(wfc) * vextv * wfc * wx*wx*scaling

         do j=1,nmax
            vext_(ij) = vext_(ij) + auxc(j)
         enddo
         
         call composite_trapeze_c(ij, time_, ft, auxcc)

 
         auxc = - conjg(wfc) * xx * wfc * wx*wx*scaling
         d_t(ij) = c0
         x_t(ij) = c0
         do j=1,nmax
            d_t(ij) = d_t(ij) + auxc(j)
            x_t(ij) = x_t(ij) - auxc(j)
         enddo
 
         call differentiation(nmax, xx, wfc, dx, dwfc)
         dwfc = -ci*dwfc
         auxc =  conjg(wfc) * dwfc * wx*wx*scaling
 
         p_t(ij) = c0
         do j=1,nmax
            p_t(ij) = p_t(ij) + auxc(j)
         enddo
 
 
!        call composite_simpson_18c(ij, time_, ft, auxcc)
         call composite_trapeze_c(ij, time_, ft, auxcc)
         pp_t(ij) = auxcc

         auxc =  conjg(wfc) * auxcc* wfc * wx*wx*scaling
         pp_t(ij) = c0
         do j=1,nmax
            pp_t(ij) = pp_t(ij) + auxc(j)
         enddo
  
        call differentiation(ij, time_, x_t, dt, dx_t)
        call differentiation(ij, time_, dx_t, dt, ddx_t)
        call differentiation(ij, time_, p_t, dt, dp_t)

         write(ut22,*) ij, time_(ij), x_t(ij),                        &
                         dx_t(ij), p_t(ij),                           &
                         dp_t(ij), ddx_t(ij), pp_t(ij), vext_(ij)
 
         
!        call composite_trapeze_c(ij, time_, t1, auxcc)


!! ***    test of the ehrenfest theorem here
!         write(222,*) ij, time_(ij), dp_t(ij), dp_t_(ij),             &
!                                      -ft(ij), p_t(ij)
!         write(333,*) ij, time_(ij), dp_t(ij), dp_t_(ij),             &
!                                            dp_t(ij)- dp_t_(ij)
!         write(444,*) ij, time_(ij), dx_t(ij), p_t(ij), pp_t(ij),     &
!                       real(p_t(ij)-pp_t(ij)), imag(p_t(ij)-pp_t(ij))
         
  
         call split_op(i,step,nmax,vextv,eigenvec,eigval,              &
                       wf0,wf)
 
! ***    The inhomogeneous term (trapezoidal rule)
         wfc_=matmul(eigenvec,wf)
         wfc=wfc_/wx/dsqrt(scaling)
         auxc2 = wfc
  
         call differentiation(nmax, xx, auxc1, dx, auxc3)
         call differentiation(nmax, xx, auxc2, dx, auxc4)

         auxc1 = auxc3 * wx * dsqrt(scaling)
         auxc2 = auxc4 * wx * dsqrt(scaling)
         auxc3 = matmul(transpose(eigenvec),auxc1)
         auxc4 = matmul(transpose(eigenvec),auxc2)

         auxc1 = exp(ci*omg_*t) * auxc3
         auxc2 = exp(ci*omg_*(t+step)) * auxc4
         call split_op(i,step,nmax,vextv,eigenvec,eigval,               &
                      wf0,auxc1)


! ***    homogeneous part
!        auxc1 = c0           
!        auxc2 = c0           
         call split_op(i,step,nmax,vextv,eigenvec,eigval,               &
                        phi0,phi)

         phi = phi - 0.5d0 * step* (auxc1 + auxc2)


         time_(ij+1)=time_(ij)+step 
         call flush(ut22) 
      enddo



! *** Data storage at the end of the pulse
      wfc_=matmul(eigenvec,wf)
      wfc=wfc_/wx/dsqrt(scaling)

      wfc_=matmul(eigenvec,phi)
      phic=wfc_/wx/dsqrt(scaling)


      do j=1,nmax
         write(ut30,*) time_(ij), j, xx(j), dsqrt(2.d0*eigval(j)),    &
            wf(j), phi(j), wfc(j), phic(j)
      enddo

      do j=1,nmax
         write(624,*) time_(ij), j, xx(j), dsqrt(2.d0*eigval(j)),     &
            wfc(j), phic(j), real(wfc(j)), imag(wfc(j)),              &
            real(phic(j)), imag(phic(j))
      enddo
      call flush(ut30) 


      deallocate(eigenvec)
      call cpu_time(finish)
      write(*,*) "time step", step
      write(*,*) "total cpu time", finish-start

      end 
