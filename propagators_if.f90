      module propagators_if
!     use constants
      implicit none
      contains

         subroutine split_op(it,dt,nmax,eigval,vext,                  &
                   nonlin,eigenvec,psi_0,psi_n)
         implicit none
         integer :: it, nmax, stage, i
         real*8 :: dt, dt2
         complex*16 :: ci, c1
         real*8, dimension(:) :: eigval, vext
         real*8, dimension(:,:) :: nonlin, eigenvec
         complex*16, dimension(:) :: psi_n, psi_0

         complex*16, dimension(:), allocatable ::  auxc
  
         allocate(auxc(nmax))

         c1 = dcmplx(1.d0,0.d0)
         ci = dcmplx(0.d0,1.d0)
 
         if (it.eq.0) then
            psi_n = psi_0
         else
            dt2 = 0.5d0*dt


   
!           call start_vec(it,dt,nmax,eigval,vext,                    &
!                     nonlin,eigenvec,psi_n,psi_n)
 ! ***   exp(-iEdt/2) * Psi
             psi_n = exp(-ci * eigval * dt2) * psi_n
 
 !     Rotate to coordinate representation
             auxc = matmul(eigenvec,psi_n) 
 !     exp(-idt Vext) * Psi
             auxc = exp(-ci * vext*dt) * auxc
 
 !     Rotate to energy representation
             psi_n = matmul(transpose(eigenvec),auxc)
 !     exp(-iEdt/2) * Psi
             psi_n = exp(-ci * eigval * dt2) * psi_n
             
         end if 

         end subroutine



         subroutine ifrk1(it,dt,nmax,eigval,vext,                     &
                   nonlin,eigenvec,psi_0,psi_n)
         implicit none
         integer :: it, nmax, stage, i
         real*8 :: dt, eta, cc1
         real*8, dimension(:) :: eigval, vext
         real*8, dimension(:,:) :: nonlin, eigenvec
         complex*16, dimension(:) :: psi_n, psi_0

         complex*16, dimension(:), allocatable ::  auxc,              &
                                   ck1
  
         allocate(auxc(nmax),                                         &
                     ck1(nmax))

         cc1=0.0d0
         if (it.eq.0) then
            psi_n = psi_0
         else
            eta = 1.d0 - cc1
            call ifrk_s1(it,dt,nmax,cc1,eigval,vext,                  &
                      nonlin,eigenvec,psi_n,ck1)

            call start_vec(it,dt,nmax,eigval,vext,                    &
                      nonlin,eigenvec,psi_n,auxc)
            
            psi_n = auxc + dt*ck1
         end if 

         end subroutine


         subroutine ifrk2(it,dt,nmax,eigval,vext,                     &
                   nonlin,eigenvec,psi_0,psi_n)
         implicit none
         integer :: it, nmax, stage, i
         real*8 :: dt, eta, cc1, cc2, alpha
         real*8, dimension(:) :: eigval, vext
         real*8, dimension(:,:) :: nonlin, eigenvec
         complex*16, dimension(:) :: psi_n, psi_0

         complex*16, dimension(:), allocatable ::  auxc,              &
                                   ck1, ck2
  
         allocate(auxc(nmax),                                         &
                    ck1(nmax),ck2(nmax))

         cc1=0.0d0
         cc2=1.d0/4
         alpha = cc2
         if (it.eq.0) then
            psi_n = psi_0
         else
            call ifrk_s1(it,dt,nmax,cc1,                              &
                      eigval,vext,nonlin,eigenvec,psi_n,ck1)
            call ifrk_s2(it,dt,nmax,cc1,cc2,                          &
                      eigval,vext,nonlin,eigenvec,psi_n,ck1,ck2)

            call start_vec(it,dt,nmax,                                &
                   eigval,vext,nonlin,eigenvec,psi_n,auxc)
            
            psi_n = auxc + dt*((1.d0-1.d0/(2*alpha))*ck1              &
                                     +    (1.0d0/(2*alpha))*ck2)               
         end if 

         end subroutine


         subroutine ifrk3(it,dt,nmax,eigval,vext,                     &
                   nonlin,eigenvec,psi_0,psi_n)
         implicit none
         integer :: it, nmax, stage, i
         real*8 :: dt, eta, cc1, cc2, cc3, alpha
         real*8, dimension(:) :: eigval, vext
         real*8, dimension(:,:) :: nonlin, eigenvec
         complex*16, dimension(:) :: psi_n, psi_0

         complex*16, dimension(:), allocatable ::  auxc,              &
                                   ck1, ck2, ck3
  
         allocate(auxc(nmax),                                         &
                    ck1(nmax),ck2(nmax),ck3(nmax))

         alpha = 1.0d0/2
         cc1=0.0d0
         cc2=alpha
         cc3=1.d0


         if (it.eq.0) then
            psi_n = psi_0
         else
            call ifrk_s1(it,dt,nmax,cc1,                              &
                      eigval,vext,nonlin,eigenvec,psi_n,ck1)
            call ifrk_s2(it,dt,nmax,cc1,cc2,                          &
                      eigval,vext,nonlin,eigenvec,psi_n,ck1,ck2)
            call ifrk_s3(it,dt,nmax,cc1,cc2,cc3,                      &
                      eigval,vext,nonlin,eigenvec,psi_n,ck1,ck2,ck3)

            call start_vec(it,dt,nmax,                                &
                   eigval,vext,nonlin,eigenvec,psi_n,auxc)


            psi_n = auxc                                              &
                     +           dt*((1.d0/2-1.d0/(6*alpha))*ck1      &
                     +         (1.d0/(6*alpha*(1.d0-alpha)))*ck2      &
                     +     ((2.d0-3*alpha)/(6*(1.d0-alpha)))*ck3 )  
            
         end if 

         end subroutine



         subroutine ifrk4(it,dt,nmax,eigval,vext,                     &
                   nonlin,eigenvec,psi_0,psi_n,rk_38)
         implicit none
         integer :: it, nmax, stage, i, rk4
         integer, optional :: rk_38
         real*8 :: dt, eta, cc1, cc2, cc3, cc4, alpha
         real*8, dimension(:) :: eigval, vext
         real*8, dimension(:,:) :: nonlin, eigenvec
         complex*16, dimension(:) :: psi_n, psi_0

         complex*16, dimension(:), allocatable ::  auxc,              &
                                   ck1, ck2, ck3, ck4
  
         allocate(auxc(nmax),                                         &
                    ck1(nmax),ck2(nmax),ck3(nmax),ck4(nmax))


         if (present(rk_38)) then
            cc1=0.0d0
            cc2=1.d0/3
            cc3=2.d0/3
            cc4=1.0d0
            rk4 = 0
         else 
            cc1=0.0d0
            cc2=1.d0/2
            cc3=1.0d0/2
            cc4=1.0d0
            rk4 = 1
         end if

         if (it.eq.0) then
            psi_n = psi_0
         else
            call ifrk_s1(it,dt,nmax,cc1,                              &
                      eigval,vext,nonlin,eigenvec,psi_n,ck1)
            call ifrk_s2(it,dt,nmax,cc1,cc2,                          &
                      eigval,vext,nonlin,eigenvec,psi_n,ck1,ck2)
            call ifrk_s3(it,dt,nmax,cc1,cc2,cc3,                      &
                     eigval,vext,nonlin,eigenvec,psi_n,ck1,ck2,ck3,rk4)
            call ifrk_s4(it,dt,nmax,cc1,cc2,cc3,cc4,rk4,              &
                     eigval,vext,nonlin,eigenvec,psi_n,ck1,ck2,ck3,ck4)

            call start_vec(it,dt,nmax,                                &
                   eigval,vext,nonlin,eigenvec,psi_n,auxc)

            if (rk4.eq.1) then
               psi_n = auxc                                           &
                     + dt*((1.d0/6)*(ck1+ 2.d0*ck2 + 2.d0*ck3 + ck4)) 
            else
               psi_n = auxc                                           &
                     + dt*((1.d0/8)*(ck1+ 3.d0*ck2 + 3.d0*ck3 + ck4)) 
            end if
         end if 

         end subroutine



! * *****************

 

          subroutine ifrk_s1(it,dt,nmax,cc1,                          &
                    eigval,vext,nonlin,eigenvec,psi_n,ck1)
          implicit none
          integer :: it, nmax, stage, i
          real*8 :: dt, dt2, eta, cc1
          complex*16 :: c0, c1, ci
          real*8, dimension(:) :: eigval, vext
          real*8, dimension(:,:) :: nonlin, eigenvec
          complex*16, dimension(:) :: ck1, psi_n
 
          complex*16, dimension(:), allocatable ::  auxc1, auxc2
 
          c0 = dcmplx(0.d0,0.d0)
          c1 = dcmplx(1.d0,0.d0)
          ci = dcmplx(0.d0,1.d0)
   
          allocate(auxc1(nmax),auxc2(nmax))
 
          dt2 = 0.5d0*dt
          eta = 1.d0 - cc1
 
          auxc1 = psi_n
 !  ***   to coordinate representation
          auxc2 = matmul(eigenvec,auxc1)

 ! *** N(Psi,t) = N * Psi
          auxc1 = auxc2
          auxc2 = -ci*matmul(nonlin,auxc1)
        
 !  ***   to energy representation
          auxc1 = matmul(transpose(eigenvec),auxc2)
 
          call integ_fac(it,dt,nmax,eta,                              &
                            eigval,vext,nonlin,eigenvec,auxc1,ck1)

          end subroutine


          subroutine ifrk_s2(it,dt,nmax,cc1,cc2,                      &
                    eigval,vext,nonlin,eigenvec,psi_n,ck1,ck2)
          implicit none
          integer :: it, nmax, stage, i
          real*8 :: dt, dt2, eta, cc1, cc2
          complex*16 :: c0, c1, ci
          real*8, dimension(:) :: eigval, vext
          real*8, dimension(:,:) :: nonlin, eigenvec
          complex*16, dimension(:) :: ck1, ck2, psi_n
 
          complex*16, dimension(:), allocatable ::  auxc1, auxc2
 
          c0 = dcmplx(0.d0,0.d0)
          c1 = dcmplx(1.d0,0.d0)
          ci = dcmplx(0.d0,1.d0)
   
          allocate(auxc1(nmax),auxc2(nmax))
 
          dt2 = 0.5d0*dt
          eta = 1.d0 - cc2
 
          auxc1 = psi_n + dt*cc2*ck1 

!         call ifrk_s1(it,dt,nmax,cc2,                                &
!                   eigval,vext,nonlin,eigenvec,psi_n,auxc1)

!         auxc1 = psi_n + cc2*dt*auxc1

 !  ***   to coordinate representation
          auxc2 = matmul(eigenvec,auxc1)

 ! *** N(Psi,t) = N * Psi
          auxc1 = auxc2
          auxc2 = -ci*matmul(nonlin,auxc1)
        
 !  ***   to energy representation
          auxc1 = matmul(transpose(eigenvec),auxc2)
 
          call integ_fac(it,dt,nmax,eta,                              &
                            eigval,vext,nonlin,eigenvec,auxc1,ck2)

          end subroutine



          subroutine ifrk_s3(it,dt,nmax,cc1,cc2,cc3,                  &
                    eigval,vext,nonlin,eigenvec,psi_n,ck1,ck2,ck3,rk4)
          implicit none
          integer :: it, nmax, stage, i
          integer, intent(in), optional :: rk4
          real*8 :: dt, dt2, eta, cc1, cc2, cc3, alpha
          complex*16 :: c0, c1, ci
          real*8, dimension(:) :: eigval, vext
          real*8, dimension(:,:) :: nonlin, eigenvec
          complex*16, dimension(:) :: ck1, ck2, ck3, psi_n
 
          complex*16, dimension(:), allocatable ::  auxc1, auxc2
 
          c0 = dcmplx(0.d0,0.d0)
          c1 = dcmplx(1.d0,0.d0)
          ci = dcmplx(0.d0,1.d0)
   
          allocate(auxc1(nmax),auxc2(nmax))
 
          dt2 = 0.5d0*dt
          eta = 1.d0 - cc3
          alpha = cc2

 
          if (present(rk4)) then

             if (rk4.eq.1) then
                cc1 = 0
                cc2 = 1.d0/2
                cc3 = cc2
                auxc1 = psi_n  + dt* cc3 * ck2 
             else
                cc1 = 0
                cc2 = 1.d0/3
                cc3 = 2.0d0/3
                auxc1 = psi_n  + dt* cc2 * ( - ck1 + 3.d0 * ck2 ) 
             end if
          else

             auxc1 = psi_n                                            &
             + dt*((1.d0+(1.d0-alpha)/(alpha*(3.d0*alpha-2.d0)))*ck1  &
                       -     (1.d0-alpha)/(alpha*(3.d0*alpha-2))*ck2 )
          end if

 !  ***   to coordinate representation
          auxc2 = matmul(eigenvec,auxc1)

 ! *** N(Psi,t) = N * Psi
          auxc1 = auxc2
          auxc2 = -ci*matmul(nonlin,auxc1)
        
 !  ***   to energy representation
          auxc1 = matmul(transpose(eigenvec),auxc2)
 
          call integ_fac(it,dt,nmax,eta,                               &
                            eigval,vext,nonlin,eigenvec,auxc1,ck3)

          end subroutine



          subroutine ifrk_s4(it,dt,nmax,cc1,cc2,cc3,cc4,rk4,          &
                    eigval,vext,nonlin,eigenvec,psi_n,ck1,ck2,ck3,ck4)
          implicit none
          integer :: it, nmax, stage, i, rk4
          real*8 :: dt, dt2, eta, cc1, cc2, cc3, cc4, alpha
          complex*16 :: c0, c1, ci
          real*8, dimension(:) :: eigval, vext
          real*8, dimension(:,:) :: nonlin, eigenvec
          complex*16, dimension(:) :: ck1, ck2, ck3, ck4, psi_n
 
          complex*16, dimension(:), allocatable ::  auxc1, auxc2
 
          c0 = dcmplx(0.d0,0.d0)
          c1 = dcmplx(1.d0,0.d0)
          ci = dcmplx(0.d0,1.d0)
   
          allocate(auxc1(nmax),auxc2(nmax))
 
          dt2 = 0.5d0*dt
          eta = 1.d0 - cc4

          if (rk4.eq.1) then 
             auxc1 = psi_n                                            &
                            + dt * cc4 * ck3 
         else
             auxc1 = psi_n                                            &
                            + dt * cc4 * ( ck1 - ck2 + ck3 ) 
         end if

 !  ***   to coordinate representation
          auxc2 = matmul(eigenvec,auxc1)

 ! *** N(Psi,t) = N * Psi
          auxc1 = auxc2
          auxc2 = -ci*matmul(nonlin,auxc1)
        
 !  ***   to energy representation
          auxc1 = matmul(transpose(eigenvec),auxc2)
 
          call integ_fac(it,dt,nmax,eta,                               &
                            eigval,vext,nonlin,eigenvec,auxc1,ck4)

          end subroutine





         subroutine start_vec(it,dt,nmax,eigval,vext,                &
                              nonlin,eigenvec, psi_n,auxc)
         implicit none
         integer :: it, nmax, stage, order, i
         real*8 :: eta, dt, dt2
         complex*16 :: c0, c1, ci
         real*8, dimension(:) :: eigval, vext
         real*8, dimension(:,:) :: nonlin, eigenvec
         complex*16, dimension(:) :: psi_n, auxc

         complex*16, dimension(:), allocatable :: auxc1, auxc2
      
         allocate(auxc1(nmax),auxc2(nmax))  

         c0 = dcmplx(0.d0,0.d0)
         c1 = dcmplx(1.d0,0.d0)
         ci = dcmplx(0.d0,1.d0)

         eta = 1.d0
         dt2 = 0.5d0 * dt

         auxc1 = psi_n
         call integ_fac(it,dt,nmax,eta,                               &
                           eigval,vext,nonlin,eigenvec,auxc1,auxc)

         end subroutine



         subroutine integ_fac(it,dt,nmax,eta,                          &
                           eigval,vext,nonlin,eigenvec,auxc1,vec)
         implicit none
         integer :: it, nmax, i
         real*8 :: eta, dt, dt2
         complex*16 :: c0, c1, ci
         real*8, dimension(:) :: eigval, vext
         real*8, dimension(:,:) :: nonlin, eigenvec
         complex*16, dimension(:) :: vec, auxc1
   
         complex*16, dimension(:), allocatable :: auxc2
   
         c0 = dcmplx(0.d0,0.d0)
         c1 = dcmplx(1.d0,0.d0)
         ci = dcmplx(0.d0,1.d0)
     
         allocate(auxc2(nmax))
   
         dt2 = 0.5d0*dt
   
! ***   exp(-iEdt/2) * Psi
         auxc1 = exp(-ci * eta * eigval * dt2) * auxc1

!     Rotate to coordinate representation
         auxc2 = matmul(eigenvec,auxc1) 
!     exp(-idt Vext) * Psi
          auxc2 = exp(-ci * eta * vext * dt) * auxc2

!     Rotate to energy representation
         auxc1 = matmul(transpose(eigenvec),auxc2)
!     exp(-iEdt/2) * Psi
         vec = exp(-ci * eta * eigval * dt2) * auxc1

         end subroutine


      end module
