      program projection
      use util
      use math_util
      use field
      use propagators
      use fedvr_conf_struct
      implicit none
      include 'open/files_index'
      integer :: lwork, info, count,                                  &
                  i, j, k, ij, p, q, r,                               &
                  nmax, krange


      character(len=32) :: arg, workdir

      real*8 :: t_ini, t_end, emax, tau,                              &
                    step, eft, t,                                     &
                    scaling, abstol, pemd, prob, p_0, p_ion, dk,      &
                    start, finish, lap, aux, aux1, aux2,              &
                    omg_step, omg0, eps

      complex*16 :: cnum, a0, cprobk, auxcc, auxcc1, auxcc2,          &
                    rr_, pk0, p0k, pkk,ppk_t,  b0w_1, b0w_2
 

      real*8, allocatable, dimension(:) :: xa, wa, xx, wx,            &
                                   vec_matup, eigval, eigval2,        &
                                   auxr1, auxr2, auxr3,               &
                                   time_, omg, kk,                    &
                                   p_k, dx, dt,                       &
                                   qvc1, qvc2, qvc3, qvc4, qvc5,      &
                                   matint_v, vextv,                   &
                                   ipiv, work, ifail, dnrm2
  

      complex*16, allocatable, dimension(:) :: wf0, wfc0,             &
                                               wf, wfc, wfc_, dwfc,   &
                                               init1, init2,          &
                                               phi0, phic0,           &
                                               phi, phic,             &
                                               ak, cprob2,            &
                                               auxc, auxk1, auxk2,    &
                                               auxc1, auxc2, auxc3,   &
                                               in_states, in_states_, &
                                               bkw_1, bkw_2,          &
                                               d_t, dd_t, d_w, dd_w,  &
                                               x_t,                   &
                                               p_t, pp_t, p_w, pp_w,  &
                                               prp_t,                 &
                                               dp_t, ft,              &
                                               v_ext, ss

      real*8, allocatable, dimension(:,:) :: lu, id, inv,             &
                                       eigenvec


      include 'param'
      include 'common/open_files'
      include 'common/laser'

 
      namelist / param / nnbr,lnbr,xmax,nsteps
!     lwork=4*nmax
      abstol = 1.0d-30
      nmax = (nnbr-1)*lnbr-1
      krange = 2.d0*nmax


      call cpu_time(start)
      print*, ""
      print*, "********************  Processing ***********************"
      print*, "&PARAMETERS"
      print*, "N_OC = ", noc
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
      allocate(eigval(nmax),init1(nmax),init2(nmax),                  &
                     wf0(nmax),wfc0(nmax),                            &
                     phic0(nmax),phi0(nmax),                          &
                     phic(nmax),phi(nmax),                            &
                     wf(nmax),wfc(nmax),wfc_(nmax),                   &
                     auxc(nmax),auxc1(nmax),auxc2(nmax))

!     allocate(eigenvec(nmax,nmax))

 
!     allocate(matint_v(nmax),vextv(nmax))
      allocate(time_(nsteps+1), ft(nsteps+1))
      allocate(d_t(nsteps+1), x_t(nsteps+1),                          &
                        dd_t(nsteps), d_w(n_omg+1), dd_w(n_omg+1),    &
                        dt(nsteps), p_t(nsteps+1), dp_t(nsteps+1),    &
                        prp_t(nsteps+1),                              &
                        pp_t(nsteps+1), p_w(n_omg+1), pp_w(n_omg+1))
 

      scaling = 0.5d0*xmax/lnbr

      step=(t_end-t_ini)/dble(nsteps)
      time_(1)=t_ini



!     count = command_argument_count()
!     i = 0
!     print *, count

!     do
!       call get_command_argument(i, arg)
!       if (len_trim(arg) == 0) exit

!       write (*,*) trim(arg)
!       i = i+1
!     end do


       call get_command_argument(1, arg)
       if (len_trim(arg) == 0) stop
       call open_files_to_plot(arg)

!      workdir = "data/run_"
!      workdir = trim(workdir)//trim(arg)//'/'
!      print*, "workdir in post_proc", workdir 
!      print*, "workdir in post_proc", workdir 

       call pulse_chr(emax,tau,t_ini,t_end)
       call plot_field

      step=(t_end-t_ini)/dble(nsteps)
      time_(1)=t_ini


      do i=0,nsteps
         ij=i+1

         read(ut22,*) j, time_(ij), x_t(ij), d_t(ij), p_t(ij)
         prp_t(ij) = real(p_t(ij))
 
         write(ut45,*) time_(ij), real(x_t(ij)), real(d_t(ij)),        &
                            real(p_t(ij)), imag(p_t(ij))
         t=time_(ij)
         time_(ij+1)=time_(ij)+step
       enddo


!      do i=1,2
!        read(ut11,*) j, aux
!        write(*,*) i, aux
!      enddo
 

       omg0 = 1.d0
       eps = 1.d-3
       dk = kmax/(krange/2-1)
  
       print *, krange, dk
       allocate(kk(krange),eigval2(krange))
       allocate(in_states(nmax), in_states_(nmax))

        do i=1,nmax
           read(ut11,*) j, eigval(i)

           read(ut10,*) xx(i), wx(i), wf0(i), phi0(i),                 &
                        wfc0(i), phic0(i)

           read(ut30,*) t_end, j, xx(i), aux,                          &
                  wf(i), phi(i), wfc(i), phic(i)
        enddo


! *** Post processing
! *** Find the first excited states...
      p=0
      do p=1,nmax
         if (eigval(p).gt.0.0d0) then
            exit
         end if
      enddo
! *** show the first excited state
      print *, p, eigval(p), wf(p), wfc(p)
!     stop

! *** probability density
      do i=1,nmax
         aux1 = real(wfc(i))**2 + imag(wfc(i))**2
         aux2 = real(phic(i))**2 + imag(phic(i))**2
         write(ut40,*) i,xx(i),real(wfc(i)),imag(wfc(i)),              &
                       real(phic(i)),imag(phic(i)), aux1, aux2
      enddo

!     uncomment below if you only want prob density!
      write(*,*) "Density prob. only"
      call flush(6)
      stop

      in_states = c0
      in_states_ = c0

      kk = 0.0d0
      rr_ = c0

     allocate(dwfc(nmax),dx(nmax))
     allocate(ak(krange),bkw_1(krange),bkw_2(krange))
     allocate(p_k(krange), auxr1(krange), auxr2(krange))
     allocate(auxk1(krange),auxk2(krange))

! *** Compute the PEMD by means of the theoretical scattering states
!  -- Quadrature

      a0=c0
      auxc = wfc0 * wfc * wx*wx*scaling

      do i=1,nmax
         a0 = auxc(i) + a0
      enddo
      a0 = exp(ci*eigval(1)*time_(nsteps+1))*a0
      write(*,*) "a0 with init", (abs(a0))**2



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

         eigval2(k) = 0.5d0* (kk(k)**2)
      enddo

      do k=1,krange

         rr_ = ci*vkpa/(-abs(kk(k)) - ci*vkpa)
         in_states = exp(ci*kk(k)*xx)                                    &
                    +    rr_ * exp( -ci*abs(kk(k)*xx))
  
         auxc1 = conjg(in_states) * wfc *wx*wx*scaling  
         auxc2 = conjg(in_states) * phic *wx*wx*scaling  

         do j=1,nmax
            ak(k) = auxc1(j) + ak(k)
            bkw_1(k) = auxc2(j) + bkw_1(k)
         enddo

         ak(k) = exp(ci*(0.5d0*kk(k)**2)*time_(nsteps+1))*ak(k)

         bkw_1(k) = exp(ci*(0.5d0*kk(k)**2)*time_(nsteps+1))*bkw_1(k)


! ***    evaluation of pk0, p0k and pkk
         pk0 = c0
         call differentiation(nmax, xx, wfc0, dx, dwfc)
         dwfc = -ci*dwfc
  
         auxc =  conjg(in_states) * dwfc * wx*wx*scaling
         do j=1,nmax
            pk0 = auxc(j) + pk0
         enddo

         p0k = c0
         call differentiation(nmax, xx, in_states, dx, dwfc)
         dwfc = -ci*dwfc
  
         auxc =  conjg(wfc0) * dwfc * wx*wx*scaling
         do j=1,nmax
            p0k = auxc(j) + p0k
         enddo

         auxcc = (2.d0*kk(k)*vkpa**(2.d0/3))/(vkpa**2+ kk(k)**2) 
         write(ut43,*) k, kk(k), pk0, p0k,                             & 
                                 real(pk0), real(p0k), real(auxcc)
 
 
         do j=1,krange
            rr_ = ci*vkpa/(-abs(kk(j)) - ci*vkpa)
            in_states_ = exp(ci*kk(j)*xx)                              &
                       +    rr_ * exp( -ci*abs(kk(j)*xx))
  
            call differentiation(nmax, xx, in_states_, dx, dwfc)
            dwfc = -ci*dwfc

            auxc =  conjg(in_states) * dwfc * wx*wx*scaling
            pkk = c0
            do i=1,nmax
               pkk = auxc(i) + pkk
            enddo

            auxcc1 = c0
            auxc = exp(ci*(kk(k)-kk(j))*xx) * wx*wx*scaling
            do i=1,nmax
               auxcc1 = auxc(i) + auxcc1
            enddo
            auxcc1 = (1.d0)/(2.d0*ppi)* auxcc1

            auxcc1 = (2.d0*ppi) * kk(k) * auxcc1
            auxcc1 = auxcc1                                            &
                - ((2.d0*vkpa)/(kk(k)**2 - kk(j)**2 +ci*eps ))         &
                * ((kk(j)*abs(kk(k)))/(abs(kk(k)) - ci*vkpa)           &
                -  (kk(k)*abs(kk(j)))/(abs(kk(j)) + ci*vkpa))

!           write(44,*) k, j, kk(k), kk(j), pkk, auxcc1

! ***       write some values for testing
            if (k.eq.(krange/2-p+1)) then
               write(ut44,*) k, j, kk(k), kk(j), pkk, auxcc1,          &
                    real(pkk), imag(pkk), real(auxcc1), imag(auxcc1)
            end if

! *** let us ignore the redundant value at kk=0.0d0 for now
!           auxk2(j)=exp(ci*(eigval2(k)+omg0 - eigval2(j))*t_end)      &
!                     * pkk * ak(j) 
            auxk2(j)=exp(ci*(eigval2(k)+omg0 - eigval2(j))*t_end)      &
                      * auxcc1 * ak(j) 

            auxk2(j)=auxk2(j)/(eigval2(k)+omg0 - eigval2(j)+ci*eps) 
            auxr2(j) = kk(j)
                  
         enddo

! *** evaluate the integral term in bkk_2
         call composite_simpson_18c(krange, auxr2, auxk2, auxcc2)
         auxcc2 = auxcc2/(2.d0*ppi)
 
! *** evaluate the nonintegral term in bkk_2
         auxcc1=exp(ci*(eigval2(k)+ omg0 - eigval(1))*t_end)           &
                      * pk0 * a0 

         auxcc1=auxcc1/(eigval2(k)+omg0 - eigval(1)) 
         bkw_2(k) = bkw_1(k)  + auxcc1 +  auxcc2

! *** let us ignore the redundant value at kk=0.0d0 for now
         auxk1(k)=exp(ci*(eigval(1)+omg0 - eigval2(k))*t_end)          &
                   * p0k * ak(k) 
         auxk1(k)=auxk1(k)/(eigval(1)+omg0 - eigval2(k) + ci*eps) 
         auxr1(k) = kk(k)

         p_k(k) = (abs(ak(k)))**2

         write(ut41,*) k, kk(k), eigval2(k),                           &
                            ak(k), bkw_1(k), bkw_2(k), p_k(k)
      enddo

! *** evaluate b0w
      call composite_simpson_18c(krange, auxr1, auxk1, auxcc1)
      auxcc1 = auxcc1/(2.d0*ppi)

      auxc = wfc0 * phic * wx*wx*scaling      
      b0w_1 = c0
      do i=1,nmax
         b0w_1 = auxc(i) + b0w_1
      enddo
      b0w_1 = exp(ci*eigval(1)*time_(nsteps+1))*b0w_1

      b0w_2 = b0w_1  + auxcc1
      write(*,*) "b0w(T), b0w", b0w_1, b0w_2

       p_ion = 0.0d0

!      do k=1,krange
!         write(20000,*) k, kk(k), p_k(k), ak(k), eigval2(k)
!      enddo
 
       auxr1=0.0d0 
       auxr2=0.0d0 
       aux1=0.0d0 
 

! *** let us ignore the redundant value at kk=0.0d0 for now
!      do k=1,krange
!        auxr1(k) = kk(k)
!        auxr2(k) = p_k(k)
!      enddo
!      call composite_simpson_18(krange, auxr1, auxr2, aux1)
       call composite_simpson_18(krange, kk, p_k, p_ion)
 
!      p_ion = p_ion + aux1

       p_ion = p_ion/(2.d0*ppi)
       p_0 = (abs(a0))**2


      write(*,*) "surv. probability", p_0, abs(wf(1))**2
      write(*,*) "ioni. probability", 1.0d0-p_0, 1.d0-abs(wf(1))**2
      write(*,*) "comp. simpson p_ion", p_ion*(2*ppi), p_ion
      write(*,*) "time step", step
      write(*,*) "k step", dk
       
      write(*,*) "p_0 + p_ion", p_0 + p_ion

! *** evaluate the hhg spectra
      deallocate(auxr1,auxr2)
      deallocate(auxc1,auxc2)
      allocate(auxc1(nsteps+1),auxc2(nsteps+1),auxc3(nsteps+1))

      allocate(qvc1(n_omg+1),qvc2(n_omg+1),qvc3(n_omg+1))
      allocate(qvc4(n_omg+1))
      allocate(qvc5(krange))
  
      call differentiation(nsteps+1, time_, d_t, dt, dd_t)
      call differentiation(nsteps+1, time_, p_t, dt, dp_t)

      allocate(omg(n_omg+1))
      omg_step = (omg_max - omg_min)/dble(n_omg)

      do k=1,n_omg+1
          omg(k)=omg_min + (k-1)*omg_step
      enddo


       do i=1,n_omg+1
          ! Eq. 80
!         auxc1 = exp(ci*omg(i)*time_) * d_t 
!         call composite_simpson_18c(nsteps+1, time_, auxc1, auxcc)
!         d_w(i) = auxcc
!         qvc1(i) =  omg(i)**2 * abs(d_w(i))**2
 
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
          pp_w(i) = exp(ci*omg(i)*time_(nsteps+1))*x_t(nsteps+1)      &
                  - exp(ci*omg(i)*time_(1))*x_t(1) - ci*omg(i)*auxcc
          qvc4(i) =  abs(pp_w(i))**2
 
 
!!         write(ut14,*) i, omg(i), d_w(i),                            &
!!                               dd_w(i), p_w(i), pp_w(i),             &
!!                               real(d_w(i)), imag(d_w(i)),           &
!!                               real(p_w(i)), imag(p_w(i)),           &
!!                               real(dd_w(i)), imag(dd_w(i)),         &
!!                               real(pp_w(i)), imag(pp_w(i))
!!
          write(ut42,*) i, omg(i), qvc1(i), qvc2(i), qvc3(i), qvc4(i)
       enddo

      deallocate(auxc1)
      allocate(auxr1(krange),auxc1(krange))
 
       do k=1,krange
          auxr1(k)= kk(k)
          auxc1(k)= (bkw_2(k))**2
       enddo
      
       call composite_simpson_18c(krange, kk, auxc1, auxcc)

       qvc5(1) = b0w_2**2 + auxcc/(2.d0*ppi)
       write(*,*) "omega step", omg_step
       write(*,*) omg0, qvc5(1)
       call cpu_time(finish)
       write(*,*) "total cpu time", finish-start



      end 
