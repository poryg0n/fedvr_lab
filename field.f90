      module field
      use util
      implicit none
      contains

      real*8 function efield(t)

!  *********************************************************************
!
!    Computes the electric field E(t) defined as :   E(t)=-d/dt A(t),
!    where A(t) denotes the potential vector given by
!                A(t)=A_0*f(t)*cos(omega*t+phase)  
!    f(t) is the envelop
!
!   Arguments
!   =========
!
!    amax   : input real*8 scalar. On entry, the maximum amplitude of
!             the potential vector. Unchanged on exit.
!
!    omega  : input real*8 scalar. On entry, the angular frequency of 
!             the laser pulse. Unchanged on exit
!
!    phase  : input real*8 scalar. On entry, the phase of the laser
!             pulse. Unchanged on exit.
!
!    tau  :   input real*8 scalar. On entry, the parameter tau of the
!             pulse envelop. Unchanged on exit.
!
!    t    :   input real*8 scalar. On entry, the time. Unchanged on
!             exit.
!
!  *********************************************************************

      implicit none

! ****** INCLUDE FILES ******

      include 'param'

! ******  DEFINES THE TYPE OF VARIABLES  ******

      double precision :: emax, tau, t, tt, &
                          t_ini, t_end
!                         envelop

!  ******   COMMONS    ******

      common/pulse_par/emax,tau,t_ini,t_end


      efield=0.0d0
      if ((t <= t_ini).or.(t >= t_end)) return

!     ppi=4.0d0*datan(1.0d0)       
      tt=ppi*t/tau
!     envelop=(sin(tt))**2
!     efield=emax*envelop*sin(omega*t+phase)

!     emax is absorbed into the envelop
      efield=envelop(t)*cos(omega*t+phase)

      return
      end
      


      real*8 function vec_pot(t)

!  *********************************************************************
!
!    Computes the potential vector A(t) defined as follows :
!              A(t)=A_0*f(t)*sin(omega*t+phase)
!    where f(t) is the envelop
!
!   Arguments
!   =========
!
!    amax   : input real*8 scalar. On entry, the maximum amplitude of
!             the potential vector. Unchanged on exit.
!
!    omega  : input real*8 scalar. On entry, the angular frequency of 
!             the laser pulse. Unchanged on exit
!
!    phase  : input real*8 scalar. On entry, the phase of the laser
!             pulse. Unchanged on exit.
!
!    tau    : input real*8 scalar. On entry, the total duration time
!             of the pulse. Unchanged on exit.
!
!    t      : input real*8 scalar. On entry, the time. Unchanged on
!             exit.
!
!  *********************************************************************

      implicit none

! ****** INCLUDE FILES ******

      include 'param'

! ******  DEFINES THE TYPE OF VARIABLES  ******

      double precision :: emax, tau, t, t_ini, t_end, &
                          np_om, nm_om, aux, auxp, &
                          auxm, ph, php, phm, vpot1, &
                          vpot2
                          

!  ******   COMMONS    ******

      common/pulse_par/emax,tau,t_ini,t_end


      vec_pot=0.0d0
!     if ((t <= t_ini).or.(t >= t_end)) return

!     np_om=dble(noc+1)/dble(noc)
!     np_om=omega*np_om
!     nm_om=dble(noc-1)/dble(noc)
!     nm_om=omega*nm_om
!     aux=0.5d0*emax/omega
!     auxp=0.25d0*emax/np_om
!     auxm=0.25d0*emax/nm_om
!     ph=omega*t+phase
!     php=np_om*t+phase
!     phm=nm_om*t+phase
!     vpot1=aux*cos(ph)-auxp*cos(php)-auxm*cos(phm)
!     vpot2=(auxp+auxm-aux)*cos(phase)

!     vec_pot=vpot1+vpot2

      return
      end
      


      real*8 function envelop(t)

!  *********************************************************************
!
!        Computes the envelop f(t) of the laser pulse
!
!   Arguments
!   =========
!
!    tau    : input real*8 scalar. On entry, the total duration time
!             of the pulse. Unchanged on exit.
!
!    t      : input real*8 scalar. On entry, the time. Unchanged on
!             exit.
!
!  *********************************************************************

      implicit none

! ****** INCLUDE FILES ******

      include 'param'

! ******  DEFINES THE TYPE OF VARIABLES  ******

      double precision :: emax, tau, t, t_ini, t_end, tt

!  ******   COMMONS    ******

      common/pulse_par/emax,tau,t_ini,t_end


      tt=t/tau
      
!     write(*,*) ppi,tau
      envelop=emax*exp(-(2.0d0*tt)**2)

      return
      end

 
      subroutine pulse_chr(emax,tau,t_ini,t_end)

!  *********************************************************************
!
!   Computes few characteristics of the laser pulse whose potential
!   vector writes
!                   A(t)=A_0*f(t)*cos(omega*t+phase)
!
!   Arguments
!   =========
!
!   funit30  : input integer scalar. Unit of the file in which
!              laser pulse informations are printed. Unchanged
!              on exit.
!
!   amax     : output real*8 scalar. On exit, the maximum amplitude
!              of the potential vector (in a.u.).
!
!   tau      : output real*8 scalar. On entry, Total time duration
!              of the laser pulse.
!
!   t_ini    : output real*8 scalar. Specifies the initial time i.e.
!              the time at which starts the interaction with the
!              pulse.
!
!   t_end    : output real*8 scalar. Specifies the final time i.e.
!              the time at which the interaction with the pulse 
!              stops. Unchanged on exit.
!
!   period   : output real*8 scalar. On exit, the period of the laser
!              pulse.
!
!   noc     : output real*8 scalar. On exit, the number of periods
!              of the laser pulse.
!
!  *********************************************************************

      implicit none

! ****** INCLUDE FILES ******

      include 'param'

! ******  DEFINES THE TYPE OF VARIABLES  ******

      integer :: ut1, ut2,ut3, ut4, ut5, ut6, ut7, ut8, ut9, ut10

      double precision :: amax, tau, t_ini, t_end, &
                          emax, i0, ip, up, kpg

! ****** COMMON STATEMENTS ***********

      common/files/ut1,ut2,ut3,ut4,ut5,ut6,ut7,ut8,ut9,ut10

! ***** COMPUTES REQUIRED DATA *******

      i0=peaki/peaki0                     ! laser intensity in a.u.
      ip=0.5d0                            ! ionization potential of H in a.u.
      up=i0/(4.0d0*omega*omega)           ! ponderomotive energy in a.u.
      kpg=sqrt(0.5d0*ip/up)               ! keldysh parameter
      
      emax=dsqrt(i0)
      amax=emax/omega
      tau=dble(noc)*period
      t_ini=-n_tau*0.5d0*tau
      t_end=n_tau*0.5d0*tau

!  *****  WRITES INFORMATIONS ABOUT THE LASER PULSE   *****

      write(ut1,120)
      write(ut1,125) emax,amax,omega,period,phase
      write(ut1,130) noc,n_tau,tau
      write(ut1,135) t_ini
      write(ut1,140) t_end
      write(ut1,145) i0
      write(ut1,150) ip
      write(ut1,155) up
      write(ut1,160) kpg

!  *****   FORMATS   *****

 120  format(///,5x,'OUTPUT INFORMATIONS ABOUT THE LASER PULSE',/, &
             5x,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~',//)
 125  format(2x,'- Laser field characteristics (in a.u.) :',//,6x, &
             'Maximum amplitude of the electric field =',2x,d21.14,//,6x, &
             'Maximum amplitude of the potential vector =',2x,d21.14, &
             //,6x,'Angular frequency =',2x,d21.14,//,6x,'Time period =', &
             2x,d21.14,//,6x,'Angular phase (in rad.) =',2x,d21.14,/)
 130  format(6x,'Number of optical cycles ', &
             '(in units of the laser period) =',2x,i4,/,6x, &
             'within the whole pulse duration',//,6x, &
             'Interaction time scaling parameter =',2x,i4,//,6x, &
             'Total time duration of the laser pulse =',2x,d21.14,/)
 135  format(6x,'initial time of the laser pulse =',2x,d21.14,/)
 140  format(6x,'final time of the laser pulse =',2x,d21.14,/)
 145  format(6x,'laser intensity (in a.u.) =',2x,d21.14,/)
 150  format(6x,'ionization potential (in a.u.) =',2x,d21.14,/)
 155  format(6x,'ponderomotive energy (in a.u.) =',2x,d21.14,/)
 160  format(6x,'keldysh parameter =',2x,d21.14,/)

      return
      end

      subroutine plot_field

!  *********************************************************************
!
!    For many values of the time t, this routine computes and prints :
!         the potential vector  A(t)=A_0*f(t)*sin(omega*t+phase) 
!         the electric field    E(t)=-d/dt A(t)
!         the pulse envelop     f(t)
!
!   Arguments
!   =========
!
!   tau       : input real*8 scalar. On entry, the total duration time
!               of the pulse. Unchanged on exit.
!
!   amax      : input real*8 scalar. On entry, the maximum amplitude
!               of the potential vector. Unchanged on exit.
!
!   nb_points : input integer scalar. On entry, the number of points
!               at which A(t), E(t) and f(t) will be evaluated.
!               Unchanged on exit.
!
!   OTHERS    : see the real*8 functions "efield", "vec_pot" and
!               "envelop" for comments on other arguments.
!
!  *********************************************************************

      implicit none

! ****** INCLUDE FILES ******

      include 'param'

      integer :: i, ut1, ut2, ut3, ut4, ut5, &
                 ut6, ut7, ut8, ut9, ut10

      double precision :: emax, tau, t, step, &
                          t_ini, t_end

!     double precision :: envelop, efield, vec_pot

! ****** COMMON STATEMENTS ***********

      common/pulse_par/emax,tau,t_ini,t_end
      common/files/ut1,ut2,ut3,ut4,ut5,ut6,ut7,ut8,ut9,ut10

    
      step=(t_end-t_ini)/dble(nb_points)

      write(ut2,110) omega,phase,peaki

      write(ut2,120)
      do i=0,nb_points
         t=t_ini+step*i
         write(ut2,*) t,envelop(t)
      enddo

      write(ut2,*) '   '
      write(ut2,*) '   '
      write(ut2,*) '   '

      write(ut2,130)
      do i=0,nb_points
         t=t_ini+step*i
         write(ut2,*) t,efield(t)
      enddo
      write(ut2,*) '   '
      write(ut2,*) '   '
      write(ut2,*) '   '

      write(ut2,140)
      do i=0,nb_points
         t=t_ini+step*i
         write(ut2,*) t,vec_pot(t)
      enddo

 110  format('#',2x,'- Laser pulse characteristics',/, &
             4x,'omega(a.u.) = ',d11.4,/,4x,'phase(rad) = ', &
             d11.4,/,4x,'peaki(W/cm**2) = ',d11.4,/,'#')
 120  format('#',/,'#',9x,'PLOT FILE FOR THE PULSE ENVELOP',/,'#')
 130  format('#',/,'#',9x,'PLOT FILE FOR THE ELECTRIC FIELD',/,'#')
 140  format('#',/,'#',9x,'PLOT FILE FOR THE VECTOR POTENTIAL',/,'#')

      return
      end
      end module
