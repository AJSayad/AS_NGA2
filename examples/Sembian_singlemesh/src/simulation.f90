!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use mpcomp_class,      only: mpcomp
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use cclabel_class,     only: cclabel
   use event_class,       only: event
   use monitor_class,     only: monitor
   use shockgen_class,    only: shockgen
   implicit none
   private; public :: simulation_init,simulation_run,simulation_final,sgen_init,update_shockprofile
   
   !> Multiphase compressible flow solver and corresponding time tracker
   type(mpcomp),      public :: fs
   type(timetracker), public :: time

   !> shock-Gen simulation
   type(shockgen), pointer :: sg=>null()
   real(WP) :: ushock          ! used to store the velocity of the shock (used for finding tmax in shock generator)
   logical, public :: sgenflag ! true for running shockgenerator

   !> CCL for postprocessing
   type(cclabel) :: ccl
   
   !> Ensight postprocessing
   type(surfmesh) :: smesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,consfile,dropfile
   
   !> Private work arrays
   real(WP), dimension(:,:,:,:,:), allocatable :: dQdt
   real(WP), dimension(:,:,:)    , allocatable :: Ui,Vi,Wi,Ma,beta,visc
   
   !> Equations of state
   real(WP) :: PinfL,GammaL,CvL
   real(WP) :: PinfG,GammaG,CvG
   
   !> Flow parameters
   real(WP) :: Ms,Xs
   real(WP) :: Grho0,GP0,u0,M0,muG
   real(WP) :: Grho1,GP1,u1,M1
   real(WP) :: Lrho0,muL
   real(WP) :: relshockvel,ReG,viscG,viscL,visc_ratio

   !> Cylinder parameters 
   real(WP) :: dcyl,xcyl

   !> Various post-processing info
   real(WP) :: Vcore,Mcore,Xcore,Ycore,Zcore !< Drop core data
   real(WP), dimension(3) :: Cmin,Cmax       !< Core extent
   
contains

   !> Function that returns a smooth Heaviside of thickness delta
   real(WP) function Hshock(x,delta)
   real(WP), intent(in) :: x,delta
   ! Goes from 0 to 1 as x goes from negative to positive
   Hshock=1.0_WP/(1.0_WP+exp(-x/delta))
   end function Hshock
   
   !> P=EOS(RHO,I) for liquid
   pure real(WP) function get_PL(RHO,I)
      implicit none
      real(WP), intent(in) :: RHO,I
      get_PL=RHO*I*(GammaL-1.0_WP)-GammaL*PinfL
   end function get_PL
   !> T=f(RHO,P) for liquid
   pure real(WP) function get_TL(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_TL=(P+PinfL)/(CvL*RHO*(GammaL-1.0_WP))
   end function get_TL
   !> C=f(RHO,P) for liquid
   pure real(WP) function get_CL(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_CL=sqrt(GammaL*(P+PinfL)/RHO)
   end function get_CL
   !> S=f(RHO,P) for liquid
   pure real(WP) function get_SL(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_SL=CvL*log((P+PinfL)/RHO**GammaL)
   end function get_SL
   
   
   !> P=EOS(RHO,I) for gas
   pure real(WP) function get_PG(RHO,I)
      implicit none
      real(WP), intent(in) :: RHO,I
      get_PG=RHO*I*(GammaG-1.0_WP)-GammaG*PinfG
   end function get_PG
   !> T=f(RHO,P) for gas
   pure real(WP) function get_TG(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_TG=(P+PinfG)/(CvG*RHO*(GammaG-1.0_WP))
   end function get_TG
   !> C=f(RHO,P) for gas
   pure real(WP) function get_CG(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_CG=sqrt(GammaG*(P+PinfG)/RHO)
   end function get_CG
   !> S=f(RHO,P) for gas
   pure real(WP) function get_SG(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_SG=CvG*log((P+PinfG)/RHO**GammaG)
   end function get_SG

   
   !> Mechanical relaxation model (generalized)
   subroutine P_relaxG(VF,Q,pjump)
      implicit none
      real(WP),                intent(inout) :: VF
      real(WP), dimension(1:), intent(inout) :: Q
      real(WP),                intent(in   ) :: pjump
      real(WP) :: PG,PL,ZG,ZL,Pint,cJ
      real(WP) :: a,b,d,n1,n0,d1,d0,Peq,VFeq
      real(WP), parameter :: RHOGmin=1.0e-3_WP
      real(WP), parameter :: phist=1.0_WP,phi0=0.0_WP   !< Temporal weighting, phist=1 should yield best results
      ! ================ Handle gas flotsams ================
      if (Q(2)/(1.0_WP-VF).lt.RHOGmin) return
      ! ================ First step for mechanical relaxation ================
      ! Get phasic pressures
      PL=get_PL(RHO=Q(1)/(       VF),I=Q(3)/Q(1))
      PG=get_PG(RHO=Q(2)/(1.0_WP-VF),I=Q(4)/Q(2))
      ! Handle limit cases - should mass/energy be tranasfered or lost? - this should probably never happen...
      if (PL.le.-PinfL) then
         print*,"****************** LIQUID CLIPPED!",PL,VF,Q
         VF=0.0_WP; Q(2)=sum(Q(1:2)); Q(1)=0.0_WP; Q(4)=sum(Q(3:4)); Q(3)=0.0_WP; return
      end if
      if (PG.le.-PinfG) then
         print*,"****************** GAS CLIPPED!",PG,VF,Q
         VF=1.0_WP; Q(1)=sum(Q(1:2)); Q(2)=0.0_WP; Q(3)=sum(Q(3:4)); Q(4)=0.0_WP; return
      end if
      ! Get phasic impedances
      ZL=Q(1)/(       VF)*get_CL(RHO=Q(1)/(       VF),P=PL)**2
      ZG=Q(2)/(1.0_WP-VF)*get_CG(RHO=Q(2)/(1.0_WP-VF),P=PG)**2
      cJ=ZL/(ZG+ZL)
      ! Calculate model interface pressure
      Pint=(ZG*PL+ZL*PG)/(ZG+ZL)
      ! Setup quadratic problem
      n1 = VF*phist
      n0 = VF*(phi0*Pint-phist*cJ*pjump)+Q(3) 
      d1 = phist+1.0_WP/(GammaL-1.0_WP)
      d0 = phi0*Pint-phist*cJ*pjump+GammaL/(GammaL-1.0_WP)*PinfL
      a = d1*( 1.0_WP/(GammaG-1.0_WP)+phist*VF) &
        + n1*(-1.0_WP/(GammaG-1.0_WP)-phist)
      b = d1*( (GammaG*PinfG-pjump)/(GammaG-1.0_WP)-Q(4)+VF*(phi0*Pint-phist*cJ*pjump)) &
        + n1*(-(GammaG*PinfG-pjump)/(GammaG-1.0_WP)-phi0*Pint+phist*cJ*pjump) &
        + d0*( 1.0_WP/(GammaG-1.0_WP)+phist*VF) &
        + n0*(-1.0_WP/(GammaG-1.0_WP)-phist   )
      d = d0*( (GammaG*PinfG-pjump)/(GammaG-1.0_WP)-Q(4)+VF*(phi0*Pint-phist*cJ*pjump)) &
        + n0*(-(GammaG*PinfG-pjump)/(GammaG-1.0_WP)-phi0*Pint+phist*cJ*pjump)
      ! Get equilibrium pressure
      Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      ! Check if pressure is sound
      if (Peq.le.max(-PinfG,-PinfL)) return
      ! Get equilibrium volume fraction
      VFeq=(n1*Peq+n0)/(d1*Peq+d0)
      ! Adjust conserved quantities
      Q(3)=Q(3)-(phi0*Pint+phist*Peq)*(VFeq-VF)
      Q(4)=Q(4)+(phi0*Pint+phist*Peq)*(VFeq-VF)
      VF=VFeq
   end subroutine P_relaxG

   !> Mechanical relaxation model
   subroutine P_relax(VF,Q)
      implicit none
      real(WP),                intent(inout) :: VF
      real(WP), dimension(1:), intent(inout) :: Q
      real(WP) :: PG,PL,ZG,ZL,Pint
      real(WP) :: a,b,d,coeffL,coeffG,Peq,VFeq
      real(WP), parameter :: RHOGmin=1.0e-3_WP
      ! ================ Handle gas flotsams ================
      if (Q(2)/(1.0_WP-VF).lt.RHOGmin) return
      ! ================ First step for mechanical relaxation ================
      ! Get phasic pressures
      PL=get_PL(RHO=Q(1)/(       VF),I=Q(3)/Q(1))
      PG=get_PG(RHO=Q(2)/(1.0_WP-VF),I=Q(4)/Q(2))
      ! Handle limit cases - should mass/energy be transfered or lost? - this should probably never happen...
      if (PL.le.-PinfL) then
         print*,"****************** LIQUID CLIPPED!",PL,VF,Q
         VF=0.0_WP; Q(2)=sum(Q(1:2)); Q(1)=0.0_WP; Q(4)=sum(Q(3:4)); Q(3)=0.0_WP; return
      end if
      if (PG.le.-PinfG) then
         print*,"****************** GAS CLIPPED!",PG,VF,Q
         VF=1.0_WP; Q(1)=sum(Q(1:2)); Q(2)=0.0_WP; Q(3)=sum(Q(3:4)); Q(4)=0.0_WP; return
      end if
      ! Get phasic impedances
      ZL=Q(1)/(       VF)*get_CL(RHO=Q(1)/(       VF),P=PL)**2
      ZG=Q(2)/(1.0_WP-VF)*get_CG(RHO=Q(2)/(1.0_WP-VF),P=PG)**2
      ! Calculate model interface pressure
      Pint=(ZG*PL+ZL*PG)/(ZG+ZL)
      ! Setup quadratic problem
      coeffL=(GammaL-1.0_WP)*Pint+2.0_WP*GammaL*PinfL
      coeffG=(GammaG-1.0_WP)*Pint+2.0_WP*GammaG*PinfG
      a=1.0_WP+GammaG*VF+GammaL*(1.0_WP-VF)
      b=coeffL*(1.0_WP-VF)+coeffG*VF-(1.0_WP+GammaG)*VF*PL-(1.0_WP+GammaL)*(1.0_WP-VF)*PG
      d=-(coeffG*VF*PL+coeffL*(1.0_WP-VF)*PG)
      ! Get equilibrium pressure
      Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      ! Check if pressure is sound
      if (Peq.le.max(-PinfG,-PinfL)) return
      ! Get equilibrium volume fraction
      VFeq=VF*((gammaL-1.0_WP)*Peq+2.0_WP*PL+coeffL)/((1.0_WP+gammaL)*Peq+coeffL)
      ! Adjust conserved quantities
      Q(3)=Q(3)-0.5_WP*(Pint+Peq)*(VFeq-VF)
      Q(4)=Q(4)+0.5_WP*(Pint+Peq)*(VFeq-VF)
      VF=VFeq
   end subroutine P_relax
   
   
   !> Mechanical relaxation model (implicit)
   subroutine P_relax_implicit(VF,Q)
      implicit none
      real(WP),                intent(inout) :: VF
      real(WP), dimension(1:), intent(inout) :: Q
      real(WP) :: a,b,d,d1,d0,Peq,VFeq,invG1G,invG1L,facG,facL
      real(WP), parameter :: RHOGmin=1.0e-3_WP
      ! Handle gas flotsams
      if (Q(2)/(1.0_WP-VF).lt.RHOGmin) return
      ! Setup quadratic problem
      invG1G=1.0_WP/(GammaG-1.0_WP); invG1L=1.0_WP/(GammaL-1.0_WP)
      d0=PinfL*GammaL*invG1L; d1=1.0_WP+invG1L
      facG=GammaG*PinfG*invG1G; facL=invG1G+VF
      a=d1*facL-VF*(invG1G+1.0_WP)
      b=d1*(facG-Q(4))-VF*facG+d0*facL-Q(3)*(invG1G+1.0_WP)
      d=d0*(facG-Q(4))-Q(3)*facG
      ! Get equilibrium pressure
      Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      ! Check if pressure is sound
      if (Peq.le.max(-PinfG,-PinfL)) return
      ! Get equilibrium volume fraction
      VFeq=(VF*Peq+Q(3))/(d1*Peq+d0)
      ! Adjust conserved quantities
      Q(3)=Q(3)-Peq*(VFeq-VF)
      Q(4)=Q(4)+Peq*(VFeq-VF)
      VF=VFeq
   end subroutine P_relax_implicit


   !> Various postprocessing
   subroutine postproc()
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE,MPI_MIN,MPI_MAX
      use parallel, only: MPI_REAL_WP
      implicit none
      integer :: i,j,k,n,ierr
      real(WP), parameter :: VFlo_extent=0.1_WP
      ! Update CCL
      call ccl%build(make_label,same_label)
      ! Core is id=1, skip core analysis if no core is present
      if (ccl%nstruct.ge.1) then
         ! Extract core volume, mass, barycenter, and extent
         Vcore=0.0_WP; Mcore=0.0_WP; Xcore=0.0_WP; Ycore=0.0_WP; Zcore=0.0_WP
         Cmin=+[huge(1.0_WP),huge(1.0_WP),huge(1.0_WP)]; Cmax=-[huge(1.0_WP),huge(1.0_WP),huge(1.0_WP)]
         do n=1,ccl%struct(1)%n_
            ! Get cell index
            i=ccl%struct(1)%map(1,n); j=ccl%struct(1)%map(2,n); k=ccl%struct(1)%map(3,n)
            ! Increment volume
            Vcore=Vcore+fs%VF(i,j,k)*fs%cfg%vol(i,j,k)
            ! Increment mass
            Mcore=Mcore+fs%Q(i,j,k,1)*fs%cfg%vol(i,j,k)
            ! Increment barycenter
            Xcore=Xcore+fs%Q(i,j,k,1)*fs%BL(1,i,j,k)*fs%cfg%vol(i,j,k)
            Ycore=Ycore+fs%Q(i,j,k,1)*fs%BL(2,i,j,k)*fs%cfg%vol(i,j,k)
            Zcore=Zcore+fs%Q(i,j,k,1)*fs%BL(3,i,j,k)*fs%cfg%vol(i,j,k)
            ! Increment extent
            Cmin=min(Cmin,[cfg%x(i  ),cfg%y(j  ),cfg%z(k  )])
            Cmax=max(Cmax,[cfg%x(i+1),cfg%y(j+1),cfg%z(k+1)])
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,Vcore,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,Mcore,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,Xcore,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); Xcore=Xcore/Mcore ! Shouldn't ever
         call MPI_ALLREDUCE(MPI_IN_PLACE,Ycore,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); Ycore=Ycore/Mcore ! be dividing by
         call MPI_ALLREDUCE(MPI_IN_PLACE,Zcore,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); Zcore=Zcore/Mcore ! zero here...
         call MPI_ALLREDUCE(MPI_IN_PLACE,Cmin ,3,MPI_REAL_WP,MPI_MIN,fs%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,Cmax ,3,MPI_REAL_WP,MPI_MAX,fs%cfg%comm,ierr)
      end if
      
   contains
      !> Function that identifies cells that need a label
      logical function make_label(i1,j1,k1)
         implicit none
         integer, intent(in) :: i1,j1,k1
         if (fs%VF(i1,j1,k1).gt.0.0_WP) then; make_label=.true.; else; make_label=.false.; end if
      end function make_label
      !> Function that identifies if cell pairs have same label
      logical function same_label(i1,j1,k1,i2,j2,k2)
         use irl_fortran_interface, only: calculateNormal,calculateCentroid
         implicit none
         integer , intent(in) :: i1,j1,k1,i2,j2,k2
         real(WP), dimension(3) :: N1,N2,O1,O2
         ! Big default, assume same label
         same_label=.true.
         ! Look more closely at the local polygon alignment to decide whether to use same labels
         if (fs%VF(i1,j1,k1).gt.0.0_WP.and.fs%VF(i1,j1,k1).lt.1.0_WP.and.fs%VF(i2,j2,k2).gt.0.0_WP.and.fs%VF(i2,j2,k2).lt.1.0_WP) then
            ! Get polygon normals
            N1=calculateNormal(fs%interface_polygon(i1,j1,k1))
            N2=calculateNormal(fs%interface_polygon(i2,j2,k2))
            ! If not at least ~75 degrees, return
            if (dot_product(N1,N2).ge.0.3_WP) return
            ! Get polygon barycenters
            O1=calculateCentroid(fs%interface_polygon(i1,j1,k1))
            O2=calculateCentroid(fs%interface_polygon(i2,j2,k2))
            ! If pointing towards one another, use different labels
            if (dot_product(O1-O2,N1).lt.0.0_WP.and.dot_product(O2-O1,N2).lt.0.0_WP) same_label=.false.
         end if
      end function same_label
   end subroutine postproc
   
   !> constant dynamic viscosity model
   subroutine cst_dyn_visc(mu,visc_cst,T)
      implicit none
      real(WP), intent(inout) :: mu       ! array to be populated 
      real(WP), intent(in)    :: visc_cst ! dynamic viscosity value 
      real(WP), intent(in)    :: T        ! temperature (not used for constant viscosity model)
      mu = visc_cst                       ! set constant viscosity
   end subroutine cst_dyn_visc

   !> sutherland model for viscosity
   subroutine sutherland_air(mu,visc_cst,T)
      implicit none
      real(WP), intent(inout) :: mu            ! viscosity array
      real(WP), intent(in)    :: T             ! temperature array
      real(WP), parameter :: mu0=1.716e-5_WP   ! [Pa*s] https://www.cfd-online.com/Wiki/Sutherland%27s_law 
      real(WP), parameter :: T0=273.15_WP      ! [K] reference temperature
      real(WP), parameter :: S=110.4_WP        ! [K] sutherland constant for air
      real(WP), intent(in) :: visc_cst         ! reference nondim dynamic viscosity (from Re # calc)
      real(WP) :: S_nondim                     ! non-dimensional sutherland constant
      mu = mu0*((T/T0)**1.5)*((T0 + S)/(T + S))! dimensional sutherlands model
   end subroutine sutherland_air

   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      ! Initialize eos and flow parameters
      initialize_parameters: block
         use string,   only: str_long
         use messager, only: log
         character(str_long) :: message
         ! read in case flags
         call param_read('Shock generator' ,sgenflag)
         ! Set PinfG to zero
         PinfG=0.0_WP
         ! Get fluid parameters
         call param_read('Liquid gamma',GammaL)
         call param_read('Gas gamma',GammaG)
         call param_read('Liquid Pinf',PinfL)
         call param_read('Liquid density',Lrho0)
         call param_read('Pre-shock density',Grho0,default=1.204_WP)
         call param_read('Pre-shock pressure',GP0,default=1.01325e5_WP)
         call param_read('Mach number of shock',Ms,default=1.47_WP)
         call param_read('Shock location',XS)
         call param_read('Liquid viscosity',muL)
         call param_read('Gas viscosity',muG)
         ! Use shock relations to get post-shock numbers (1)
         GP1 = GP0 * (2.0_WP*GammaG*Ms**2 - (GammaG-1.0_WP)) / (GammaG+1.0_WP)
         Grho1 = Grho0 * (Ms**2 * (GammaG+1.0_WP) / ((GammaG-1.0_WP)*Ms**2 + 2.0_WP))
         ! Calculate post-shock Mach number
         M1 = sqrt(((GammaG-1.0_WP)*(Ms**2)+2.0_WP)/(2.0_WP*GammaG*(Ms**2)-(GammaG-1.0_WP)))
         ! Calculate post-shock velocity
         u1 = -M1 * sqrt(GammaG*GP1/Grho1) + Ms*sqrt(GammaG*GP0/Grho0)
         ! Velocity at which shock moves
         relshockvel = -Grho1*u1/(Grho0-Grho1)
         if (cfg%amRoot) then
           print*,"===== Problem Setup Description ====="
           print*,'Mach number', Ma
           print*,'Pre-shock:  Density',Grho0,'Pressure',GP0
           print*,'Post-shock: Density',Grho1,'Pressure',GP1,'Velocity',u1
           print*,'Shock velocity', relshockvel
         end if
         ! Initialize shock conditions
         ! Calculate specific heats
         ! CvL=(LP0+PinfL)/(Lrho0*1.0_WP*(GammaL-1.0_WP))
         ! CvG=(GP0)/(Grho0*1.0_WP*(GammaG-1.0_WP))
         CvL=1077.7_WP
         CvG=1956.45_WP
         ! Initialize liquid at left
         call param_read('Cylinder diameter',dcyl)
         call param_read('Cylinder location',xcyl)
      end block initialize_parameters
      
      ! Initialize time tracker
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max time',time%tmax)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         time%dt=time%dtmax
      end block initialize_timetracker
      
      ! Create multipgase compressible flow solver
      create_velocity_solver: block
         ! Initialize solver with required thermodynamic functions
         call fs%initialize(cfg=cfg,name='Compressible NS')
         ! Provide relaxation and thermodynamic models
         fs%relax=>P_relax_implicit
         fs%getPL=>get_PL; fs%getCL=>get_CL; fs%getSL=>get_SL; fs%getTL=>get_TL
         fs%getPG=>get_PG; fs%getCG=>get_CG; fs%getSG=>get_SG; fs%getTG=>get_TG
      end block create_velocity_solver
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(dQdt(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:fs%nQ,1:4))
         allocate(beta(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(visc(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ma(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays

      ! Prepare post-processing
      prep_postprocess: block
         call ccl%initialize(pg=cfg%pgrid,name='ccl')
      end block prep_postprocess
      
      ! Prepare initial conditions
      initial_conditions: block
         use irl_fortran_interface, only: setNumberOfPlanes,setPlane
         use mms_geom,              only: initialize_volume_moments
         use mpcomp_class,          only: VFlo
         integer :: i,j,k
         ! Initialize primary variables
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  fs%BL(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]; fs%BG(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]
                  call setNumberOfPlanes(fs%PLIC(i,j,k),1); call setPlane(fs%PLIC(i,j,k),0,[0.0_WP,0.0_WP,0.0_WP],sign(1.0_WP,fs%VF(i,j,k)-0.5_WP))
                  ! Now set volume moments for a droplet or a slab
                  call initialize_volume_moments(lo=[fs%cfg%x(i),fs%cfg%y(j),fs%cfg%z(k)],hi=[fs%cfg%x(i+1),fs%cfg%y(j+1),fs%cfg%z(k+1)],&
                  levelset=levelset_cyl,time=0.0_WP,level=5,VFlo=VFlo,VF=fs%VF(i,j,k),BL=fs%BL(:,i,j,k),BG=fs%BG(:,i,j,k))
                  fs%U(i,j,k)=u1*Hshock(Xs-fs%cfg%x(i),delta=0.5_WP*fs%dx)
                  fs%V(i,j,k)=0.0_WP; fs%W(i,j,k)=0.0_WP
                  if (fs%VF(i,j,k).lt.1.0_WP) then
                     fs%RHOG(i,j,k)=Grho0+(Grho1-Grho0)*Hshock(Xs-fs%cfg%xm(i),delta=0.5_WP*fs%dx)
                     fs%PG  (i,j,k)=GP0  +(GP1  -GP0  )*Hshock(Xs-fs%cfg%xm(i),delta=0.5_WP*fs%dx)
                     fs%IG  (i,j,k)=(fs%PG(i,j,k)+GammaG*PinfG)/(fs%RHOG(i,j,k)*(GammaG-1.0_WP))
                  end if
                  if (fs%VF(i,j,k).gt.0.0_WP) then
                     fs%PL  (i,j,k)=GP0
                     fs%RHOL(i,j,k)=Lrho0
                     fs%IL  (i,j,k)=(fs%PL(i,j,k)+GammaL*PinfL)/(fs%RHOL(i,j,k)*(GammaL-1.0_WP))
                  end if
               end do
            end do
         end do
         ! Build PLIC interface
         call fs%build_interface()
         ! Initialize conserved variables
         fs%Q(:,:,:,1)=        fs%VF *fs%RHOL
         fs%Q(:,:,:,2)=(1.0_WP-fs%VF)*fs%RHOG
         fs%Q(:,:,:,3)= fs%Q(:,:,:,1)*fs%IL
         fs%Q(:,:,:,4)= fs%Q(:,:,:,2)*fs%IG
         call fs%get_momentum()
         ! Communicate conserved variables (not needed in general, but allows 2D runs without changing loop above...)
         do i=1,fs%nQ; call fs%cfg%sync(fs%Q(:,:,:,i)); end do
         ! Apply user-provided relaxation model
         call fs%apply_relax()
         ! Rebuild primitive variables
         call fs%get_primitive()
         ! Interpolate velocity
         call fs%interp_vel(Ui,Vi,Wi)
         ! Compute local Mach number
         Ma=sqrt(Ui**2+Vi**2+Wi**2)/fs%C
      end block initial_conditions

      shockgen_setup: block
         use param,    only: param_read
         use parallel, only: group,amRoot
         use mpi_f08,  only: MPI_GROUP,MPI_COMM_WORLD,MPI_DOUBLE_PRECISION
         real(WP), dimension(3) :: X0=[0.0_WP,0.0_WP,0.0_WP]
         integer , dimension(3) :: meshsize,sgen_partition
         integer  :: i,j,k,ierr,rank,nshock,shock_index,nx
         integer  :: q=2 ! (mpcomp RHOG is Q(i,j,k,2), spmcomp RHO is Q(i,j,k,1)) used in update_shockprofile
         real(WP) :: dx,Xend,Lx
         type(MPI_GROUP) :: sgen_group
         real(WP), dimension(:), allocatable :: RHOG_profile,IG_profile,PG_profile,U_profile ! shock profile arrays
         ushock = relshockvel ! set shock velocity from main sim
         if (sgenflag.eqv.(.true.))then
            ! allocate shock profile arrays on all procs
            call param_read('nshock',nshock,default=8) ! set number of points for shock profile (left and right of center total pts = 2*nshock+1)
            allocate(RHOG_profile(2*nshock+1),PG_profile(2*nshock+1),IG_profile(2*nshock+1),U_profile(2*nshock+2)) ! add an extra point for U for staggered grid
            RHOG_profile = 0.0_WP; PG_profile = 0.0_WP; IG_profile = 0.0_WP; U_profile = 0.0_WP

            ! create shockgen group
            call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
            call MPI_GROUP_INCL(cfg%group,1,0,sgen_group,ierr)
            
            ! read in mesh size
            call param_read('Lx',Lx); call param_read('nx',nx); dx=Lx/real(nx,WP)
            Xend = dcyl + Xs !> let the shock travel 1 diameter
            ! print*, "DEBUG: sim.f90: shockgen_setup: Xend = ", Xend
            ! print*, "DEBUG: sim.f90: shockgen_setup: dcyl = ", dcyl
            ! print*, "DEBUG: sim.f90: shockgen_setup: Xs = ", Xs
            ! print*, "DEBUG: sim.f90: shockgen_setup: Lx = ", Lx
            if(amRoot)then   !> setup and run sg only on root proc
               call sgen_init(dx,[nx,1,1],X0,sgen_group,viscG,Xs,Xend,ushock)
               !> run shock gen simulation
               do while (.not.sg%time%done())
                  call sg%step()
               end do
               call sg%finalize(nshock,Xend,RHOG_profile,IG_profile,PG_profile,U_profile) 
            end if

            ! broadcast shock profile arrays to all procs
            call MPI_BCAST(RHOG_profile,2*nshock+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
            call MPI_BCAST(IG_profile,  2*nshock+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(PG_profile,  2*nshock+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(U_profile,   2*nshock+2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

            ! find and update shock 
            shock_index = ceiling(abs(Xs - cfg%xm(1))/cfg%dx(1))
            call update_shockprofile(fs%RHOG(:,:,:),fs%PG(:,:,:),fs%IG(:,:,:),fs%U(:,:,:), &
                                    & RHOG_profile,PG_profile,IG_profile,U_profile,cfg%imino_,cfg%imaxo_, &
                                    & cfg%jmino_,cfg%jmaxo_,cfg%kmino_,cfg%kmaxo_,nshock,shock_index) 

            ! rebuild conserved quantites
            call fs%build_interface()
            ! Initialize conserved variables
            fs%Q(:,:,:,1)=        fs%VF *fs%RHOL
            fs%Q(:,:,:,2)=(1.0_WP-fs%VF)*fs%RHOG
            fs%Q(:,:,:,3)= fs%Q(:,:,:,1)*fs%IL
            fs%Q(:,:,:,4)= fs%Q(:,:,:,2)*fs%IG 
            call fs%get_momentum() 
            ! Communicate conserved variables
            do i=1,fs%nQ; call fs%cfg%sync(fs%Q(:,:,:,i)); end do
            ! Rebuild primitive variables
            call fs%get_primitive() 
            ! Compute local Mach number
            Ma=sqrt(Ui**2+Vi**2+Wi**2)/fs%C
            ! deallocate vars
            deallocate(RHOG_profile,IG_profile,PG_profile,U_profile)
         end if
      end block shockgen_setup
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='ShockDrop')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('VOF',fs%VF)
         call ens_out%add_scalar('RHOL',fs%RHOL)
         call ens_out%add_scalar('RHOG',fs%RHOG)
         call ens_out%add_scalar('IL',fs%IL)
         call ens_out%add_scalar('IG',fs%IG)
         call ens_out%add_scalar('PL',fs%PL)
         call ens_out%add_scalar('PG',fs%PG)
         call ens_out%add_scalar('Mach',Ma)
         call ens_out%add_scalar('beta',beta)
         call ens_out%add_scalar('visc',visc)
         call ens_out%add_scalar('visctot',fs%visc)
         ! Create surface mesh for PLIC
         smesh=surfmesh(nvar=0,name='plic')
         call fs%update_surfmesh(smesh)
         call ens_out%add_surface('plic',smesh)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      ! Create monitor files
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         call fs%get_info()
         ! Create simulation monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%RHOLmax,'max(RHOL)')
         call mfile%add_column(fs%RHOLmin,'min(RHOL)')
         call mfile%add_column(fs%ILmax  ,'max(IL)'  )
         call mfile%add_column(fs%ILmin  ,'min(IL)'  )
         call mfile%add_column(fs%PLmax  ,'max(PL)'  )
         call mfile%add_column(fs%PLmin  ,'min(PL)'  )
         call mfile%add_column(fs%TLmax  ,'max(TL)'  )
         call mfile%add_column(fs%TLmin  ,'min(TL)'  )
         call mfile%add_column(fs%RHOGmax,'max(RHOG)')
         call mfile%add_column(fs%RHOGmin,'min(RHOG)')
         call mfile%add_column(fs%IGmax  ,'max(IG)'  )
         call mfile%add_column(fs%IGmin  ,'min(IG)'  )
         call mfile%add_column(fs%PGmax  ,'max(PG)'  )
         call mfile%add_column(fs%PGmin  ,'min(PG)'  )
         call mfile%add_column(fs%TGmax  ,'max(TG)'  )
         call mfile%add_column(fs%TGmin  ,'min(TG)'  )
         call mfile%add_column(fs%VFmax  ,'VFmax'    )
         call mfile%add_column(fs%VFmin  ,'VFmin'    )
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLa_x,'Acoustic xCFL')
         call cflfile%add_column(fs%CFLa_y,'Acoustic yCFL')
         call cflfile%add_column(fs%CFLa_z,'Acoustic zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
         ! Create conservation monitor
         consfile=monitor(fs%cfg%amRoot,'conservation')
         call consfile%add_column(time%n,'Timestep number')
         call consfile%add_column(time%t,'Time')
         call consfile%add_column(fs%VFint  ,'Volume')
         call consfile%add_column(fs%Qint(1),'Liquid mass')
         call consfile%add_column(fs%Qint(2),'Gas mass')
         call consfile%add_column(fs%Qint(3),'Liquid energy')
         call consfile%add_column(fs%Qint(4),'Gas energy')
         call consfile%add_column(fs%Qint(5),'U Momentum')
         call consfile%add_column(fs%Qint(6),'V Momentum')
         call consfile%add_column(fs%Qint(7),'W Momentum')
         call consfile%add_column(fs%RHOKLint,'Liquid KE')
         call consfile%add_column(fs%RHOKGint,'Gas KE')
         call consfile%add_column(fs%RHOSLint,'Liquid entropy')
         call consfile%add_column(fs%RHOSGint,'Gas entropy')
         call consfile%write()
         ! Create drop output
         dropfile=monitor(fs%cfg%amRoot,'Core_info')
         call dropfile%add_column(time%n,'Timestep number')
         call dropfile%add_column(time%t,'Time')
         call dropfile%add_column(fs%VFint  ,'Total volume')
         call dropfile%add_column(fs%Qint(1),'Total mass')
         call dropfile%add_column(ccl%nstruct,'N drops')
         call dropfile%add_column(Vcore,'Core volume')
         call dropfile%add_column(Mcore,'Core mass')
         call dropfile%add_column(Xcore,'Core X')
         call dropfile%add_column(Ycore,'Core Y')
         call dropfile%add_column(Zcore,'Core Z')
         call dropfile%add_column(Cmin(1),'Core Xmin')
         call dropfile%add_column(Cmin(2),'Core Ymin')
         call dropfile%add_column(Cmin(3),'Core Zmin')
         call dropfile%add_column(Cmax(1),'Core Xmax')
         call dropfile%add_column(Cmax(2),'Core Ymax')
         call dropfile%add_column(Cmax(3),'Core Zmax')
         call dropfile%write()
      end block create_monitor
      
   contains
      !> Function that defines a level set function for a cylindrical droplet
      function levelset_cyl(xyz,t) result(G)
         implicit none
         real(WP), dimension(3),intent(in) :: xyz
         real(WP), intent(in) :: t
         real(WP) :: G
         G=1.0_WP-sqrt(((xyz(1)-xcyl)/dcyl*2.0_WP)**2+(xyz(2)/dcyl*2.0_WP)**2)
      end function levelset_cyl
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation (RK4)
   subroutine simulation_run
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Remember conserved variables
         fs%Qold=fs%Q
         
         ! Remember phasic quantities
         fs%RHOLold=fs%RHOL; fs%ILold=fs%IL; fs%PLold=fs%PL
         fs%RHOGold=fs%RHOG; fs%IGold=fs%IG; fs%PGold=fs%PG

         ! Remember volume moments and interface
         fs%VFold=fs%VF
         fs%BLold=fs%BL
         fs%BGold=fs%BG
         copy_plic_to_old: block
            use irl_fortran_interface, only: copy
            integer :: i,j,k
            do k=fs%cfg%kmino_,fs%cfg%kmaxo_; do j=fs%cfg%jmino_,fs%cfg%jmaxo_; do i=fs%cfg%imino_,fs%cfg%imaxo_
               call copy(fs%PLICold(i,j,k),fs%PLIC(i,j,k))
            end do; end do; end do
         end block copy_plic_to_old

         ! Tag cells for semi-Lagrangian transport
         call fs%SLtag()

         ! Prepare SGS viscosity models
         call fs%get_vreman(dt=time%dt,visc=visc)
         call fs%get_viscartif(dt=time%dt,beta=beta)
         mixture_viscosity: block
            integer  :: i,j,k
            real(WP) :: Lvof,Lrho,Gvof,Grho
            real(WP) :: Lvisc,Gvisc,Lbeta,Gbeta
            real(WP), parameter :: eps=1.0e-15_WP
            do k=fs%cfg%kmino_+1,fs%cfg%kmaxo_-1; do j=fs%cfg%jmino_+1,fs%cfg%jmaxo_-1; do i=fs%cfg%imino_+1,fs%cfg%imaxo_-1
               ! Create smooth mass info distribution
               Lvof=sum(       fs%VF(i-1:i+1,j-1:j+1,k-1:k+1)  )
               Gvof=sum(1.0_WP-fs%VF(i-1:i+1,j-1:j+1,k-1:k+1)  )
               Lrho=sum(       fs%Q (i-1:i+1,j-1:j+1,k-1:k+1,1))/(Lvof+eps)
               Grho=sum(       fs%Q (i-1:i+1,j-1:j+1,k-1:k+1,2))/(Gvof+eps)
               ! Harmonic average of VISC
               ! Lvisc=Lrho*(muL+visc(i,j,k)); Gvisc=Grho*(muG+visc(i,j,k)); fs%VISC(i,j,k)=(Lvof+Gvof)/(Lvof/max(Lvisc,eps)+Gvof/max(Gvisc,eps))
               Lvisc=Lrho*visc(i,j,k)+muL; Gvisc=Grho*visc(i,j,k)+muG; fs%VISC(i,j,k)=(Lvof+Gvof)/(Lvof/max(Lvisc,eps)+Gvof/max(Gvisc,eps))
               ! Harmonic average of BETA
               Lbeta=Lrho*beta(i,j,k); Gbeta=Grho*beta(i,j,k); fs%BETA(i,j,k)=(Lvof+Gvof)/(Lvof/max(Lbeta,eps)+Gvof/max(Gbeta,eps))
            end do; end do; end do
         end block mixture_viscosity

         ! Perform first semi-Lagrangian transport step =====================================================
         call fs%SLstep(dt=0.5_WP*time%dt,U=fs%U,V=fs%V,W=fs%W)
         call fs%build_interface()

         ! First RK step ====================================================================================
         ! Get non-SL RHS and increment
         call fs%rhs(dQdt(:,:,:,:,1))
         fs%Q=fs%Qold+0.5_WP*time%dt*dQdt(:,:,:,:,1)
         ! Increment Q with SL terms
         fs%Q=fs%Q+fs%SLdQ
         ! Recompute primitive variables
         call fs%get_primitive()

         ! Second RK step ===================================================================================
         ! Get non-SL RHS and increment
         call fs%rhs(dQdt(:,:,:,:,2))
         fs%Q=fs%Qold+0.5_WP*time%dt*dQdt(:,:,:,:,2)
         ! Increment Q with SL terms
         fs%Q=fs%Q+fs%SLdQ
         ! Recompute primitive variables
         call fs%get_primitive()

         ! Perform second semi-Lagrangian transport step ====================================================
         call fs%SLstep(dt=1.0_WP*time%dt,U=fs%U,V=fs%V,W=fs%W)
         call fs%build_interface()

         ! Third RK step ====================================================================================
         ! Get non-SL RHS and increment
         call fs%rhs(dQdt(:,:,:,:,3))
         fs%Q=fs%Qold+1.0_WP*time%dt*dQdt(:,:,:,:,3)
         ! Increment Q with SL terms
         fs%Q=fs%Q+fs%SLdQ
         ! Recompute primitive variables
         call fs%get_primitive()

         ! Fourth RK step ===================================================================================
         ! Get non-SL RHS and increment
         call fs%rhs(dQdt(:,:,:,:,4))
         fs%Q=fs%Qold+time%dt/6.0_WP*(dQdt(:,:,:,:,1)+2.0_WP*dQdt(:,:,:,:,2)+2.0_WP*dQdt(:,:,:,:,3)+dQdt(:,:,:,:,4))
         ! Increment Q with SL terms
         fs%Q=fs%Q+fs%SLdQ
         ! Apply user-provided relaxation model
         call fs%apply_relax()
         ! Recompute primitive variables
         call fs%get_primitive()
         ! Call BCs
         call apply_bconds()

         ! Interpolate velocity
         call fs%interp_vel(Ui,Vi,Wi)

         ! Compute local Mach number
         Ma=sqrt(Ui**2+Vi**2+Wi**2)/fs%C

         ! Output to ensight
         if (ens_evt%occurs()) then
            call fs%update_surfmesh(smesh)
            call ens_out%write_data(time%t)
         end if

         ! Perform and output monitoring
         call fs%get_info()
         call mfile%write()
         call cflfile%write()
         call consfile%write()
         call postproc()
         call dropfile%write()

      end do
      
   end subroutine simulation_run

   ! AS_SGEN: CLEAN THIS BLOCK UP, ADD VISCOSITY MOFELS FOR SGEN here?
   subroutine sgen_init(dx,meshsize,X0,sgen_group,viscG,Xs,Xend,ushock)
      use param,    only: param_read
      use parallel, only: group,amRoot
      use mpi_f08,  only: MPI_GROUP,MPI_COMM_WORLD,MPI_DOUBLE_PRECISION
      implicit none
      real(WP), intent(in)  :: dx,viscG,ushock,Xs,Xend
      real(WP), dimension(3), intent(in) :: X0
      integer, dimension(3), intent(in) :: meshsize
      type(MPI_GROUP), intent(in) :: sgen_group
      integer, dimension(3) :: sgen_partition=(/1,1,1/) ! run in serial
      integer  :: i,j,k
      ! Allocate and initialize the solver
      allocate(sg); call sg%initialize(dx=dx,meshsize=meshsize,startloc=X0,sgen_group=sgen_group,partition=sgen_partition)
      ! Provide thermodynamic model
      sg%fs%getP=>get_PG; sg%fs%getC=>get_CG; sg%fs%getS=>get_SG; sg%fs%getT=>get_TG
      ! We need to transfer our viscosity explicitly...
      sg%cst_visc=viscG
      ! choose our viscosity model
      sg%fs%visc_model=>cst_dyn_visc
      ! time info
      sg%time%dtmax = time%dtmax; sg%time%dt = time%dtmax; sg%time%cflmax = time%cflmax 
      sg%ens_evt%tper = ens_evt%tper    ! uncomment this line (and ensight lines in shockgen_class to save data files for sgen)
      sg%time%tmax = (Xend - Xs)/ushock ! this is how long it will take to reach Xend
      do k=sg%cfg%kmino_,sg%cfg%kmaxo_
         do j=sg%cfg%jmino_,sg%cfg%jmaxo_
            do i=sg%cfg%imino_,sg%cfg%imaxo_
               sg%fs%U(i,j,k)=u1*Hshock(Xs-sg%fs%cfg%x(i),delta=0.5_WP*sg%fs%dx); sg%fs%V(i,j,k)=0.0_WP; sg%fs%W(i,j,k)=0.0_WP ! velocity
               sg%fs%Q(i,j,k,1)=Grho0+(Grho1-Grho0)*Hshock(Xs-sg%fs%cfg%xm(i),delta=0.5_WP*sg%fs%dx)                              ! density
               sg%fs%P   (i,j,k)=GP0  +(GP1  -GP0 )*Hshock(Xs-sg%fs%cfg%xm(i),delta=0.5_WP*sg%fs%dx)                             ! pressure
               sg%fs%I   (i,j,k)=(sg%fs%P(i,j,k)+GammaG*PinfG)/(sg%fs%Q(i,j,k,1)*(GammaG-1.0_WP))                              ! internal energy
            end do
         end do
      end do
      ! Initialize conserved variables
      sg%fs%Q(:,:,:,2)=sg%fs%Q(:,:,:,1)*sg%fs%I
      call sg%fs%get_momentum()
      ! Rebuild primitive variables
      call sg%fs%get_primitive()
      ! Interpolate velocity
      call sg%fs%interp_vel(sg%Ui,sg%Vi,sg%Wi)
      ! Compute local Mach number
      sg%Ma=sqrt(sg%Ui**2+sg%Vi**2+sg%Wi**2)/sg%fs%C
      call sg%output_ensight(t=sg%time%t)
   end subroutine sgen_init

   subroutine update_shockprofile(rhog,pg,ig,u,rhog_profile,pg_profile,ig_profile,U_profile,IMIN,IMAX,JMIN,JMAX,KMIN,KMAX,nshock,shock_index)
      use parallel, only: amRoot
      implicit none
      integer, intent(in) :: shock_index,nshock,IMIN,IMAX,JMIN,JMAX,KMIN,KMAX
      real(WP), dimension(:), intent(in) :: rhog_profile,pg_profile,ig_profile,u_profile
      real(WP), dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)  , intent(inout) :: rhog,pg,ig,u
      integer :: i
      ! zero everything out we do this to ensure there are no leftover values from Heaviside function (overkill)
      rhog = 0.0_WP; pg = 0.0_WP; ig = 0.0_WP; u = 0.0_WP
      ! set initial values as a discontinuity
      do i = IMIN,IMAX
         if (i.lt.shock_index) then
            rhog(i,:,:) = Grho1
            pg  (i,:,:) = GP1
            ig  (i,:,:) =(GP1 + GammaG*PinfG)/(Grho1*(GammaG-1.0_WP))
            u   (i,:,:) = u1
         else
            rhog(i,:,:) = Grho0
            pg  (i,:,:) = GP0
            ig  (i,:,:) =(GP0 + GammaG*PinfG)/(Grho0*(GammaG-1.0_WP))
            u   (i,:,:) = 0.0_WP
         end if
      end do
      ! overwrite with shock profile where appropriate
      do i = IMIN,IMAX
         if ((i.ge.shock_index - nshock).and.(i.le. shock_index + nshock)) then   ! shock-profile
            rhog(i,:,:) = rhog_profile(i - (shock_index - nshock) + 1)
            pg  (i,:,:) = pg_profile  (i - (shock_index - nshock) + 1)
            ig  (i,:,:) = ig_profile  (i - (shock_index - nshock) + 1)
         end if
      end do
      do i = IMIN,IMAX
         if ((i.ge.shock_index - nshock).and.(i.le. shock_index + nshock+1)) then ! shock-profile for velocity (includes added point for cell faces)
            u(i,:,:) = u_profile(i - (shock_index - nshock) + 1)
         end if
      end do
   end subroutine update_shockprofile

   !> Apply boundary conditions
   subroutine apply_bconds()
      use irl_fortran_interface, only: setPlane
      implicit none
      integer :: i,j,k
      
      ! Apply clipped Neumann on primitive variables in x+
      if (.not.fs%cfg%xper.and.fs%cfg%iproc.eq.fs%cfg%npx) then
         do k=fs%cfg%kmino_,fs%cfg%kmaxo_; do j=fs%cfg%jmino_,fs%cfg%jmaxo_
            ! Copy over from imax to imax+1 and above
            do i=fs%cfg%imax+1,fs%cfg%imaxo
               ! Copy primitive variables
               fs%RHOL(i,j,k)=fs%RHOL(fs%cfg%imax,j,k)
               fs%PL  (i,j,k)=fs%PL  (fs%cfg%imax,j,k)
               fs%IL  (i,j,k)=fs%IL  (fs%cfg%imax,j,k)
               fs%RHOG(i,j,k)=fs%RHOG(fs%cfg%imax,j,k)
               fs%PG  (i,j,k)=fs%PG  (fs%cfg%imax,j,k)
               fs%IG  (i,j,k)=fs%IG  (fs%cfg%imax,j,k)
               fs%U  (i,j,k)=max(fs%U(fs%cfg%imax,j,k),0.0_WP) ! AS apply clipped Neumann for velocity in x+ to prevent inflow 
               ! fs%U   (i,j,k)=fs%U   (fs%cfg%imax,j,k)
               fs%V   (i,j,k)=fs%V   (fs%cfg%imax,j,k)
               fs%W   (i,j,k)=fs%W   (fs%cfg%imax,j,k)
               fs%VF  (i,j,k)=fs%VF  (fs%cfg%imax,j,k)
               ! Also adjust interface data
               call setPlane(fs%PLIC(i,j,k),0,[+1.0_WP,0.0_WP,0.0_WP],fs%cfg%x(i)+fs%dx*fs%VF(i,j,k))
               fs%BL(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]
               fs%BG(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]
            end do
         end do; end do
      end if
      
      ! Apply clipped Neumann on primitive variables in y+
      if (.not.fs%cfg%yper.and.fs%cfg%jproc.eq.fs%cfg%npy) then
         do k=fs%cfg%kmino_,fs%cfg%kmaxo_; do i=fs%cfg%imino_,fs%cfg%imaxo_
            ! Copy over from jmax to jmax+1 and above
            do j=fs%cfg%jmax+1,fs%cfg%jmaxo
               ! Copy primitive variables
               fs%RHOL(i,j,k)=fs%RHOL(i,fs%cfg%jmax,k)
               fs%PL  (i,j,k)=fs%PL  (i,fs%cfg%jmax,k)
               fs%IL  (i,j,k)=fs%IL  (i,fs%cfg%jmax,k)
               fs%RHOG(i,j,k)=fs%RHOG(i,fs%cfg%jmax,k)
               fs%PG  (i,j,k)=fs%PG  (i,fs%cfg%jmax,k)
               fs%IG  (i,j,k)=fs%IG  (i,fs%cfg%jmax,k)
               fs%U   (i,j,k)=fs%U   (i,fs%cfg%jmax,k)
               fs%V  (i,j,k)=max(fs%V(i,fs%cfg%jmax,k),0.0_WP) ! AS apply clipped Neumann for velocity in y+ to prevent inflow 
               !fs%V   (i,j,k)=fs%V   (i,fs%cfg%jmax,k)
               fs%W   (i,j,k)=fs%W   (i,fs%cfg%jmax,k)
               fs%VF  (i,j,k)=fs%VF  (i,fs%cfg%jmax,k)
               ! Also adjust interface data
               call setPlane(fs%PLIC(i,j,k),0,[0.0_WP,+1.0_WP,0.0_WP],fs%cfg%y(j)+fs%dy*fs%VF(i,j,k))
               fs%BL(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]
               fs%BG(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]
            end do
         end do; end do
      end if
      
      ! Apply clipped Neumann on primitive variables in y-
      if (.not.fs%cfg%yper.and.fs%cfg%jproc.eq.1) then
         do k=fs%cfg%kmino_,fs%cfg%kmaxo_; do i=fs%cfg%imino_,fs%cfg%imaxo_
            ! First copy over V from jmin+1 to jmin
            fs%V(i,fs%cfg%jmin,k)=min(fs%V(i,fs%cfg%jmin+1,k),0.0_WP)
            ! Then copy over from jmin to jmin-1 and below
            do j=fs%cfg%jmino,fs%cfg%jmin-1
               ! Copy primitive variables
               fs%RHOL(i,j,k)=fs%RHOL(i,fs%cfg%jmin,k)
               fs%PL  (i,j,k)=fs%PL  (i,fs%cfg%jmin,k)
               fs%IL  (i,j,k)=fs%IL  (i,fs%cfg%jmin,k)
               fs%RHOG(i,j,k)=fs%RHOG(i,fs%cfg%jmin,k)
               fs%PG  (i,j,k)=fs%PG  (i,fs%cfg%jmin,k)
               fs%IG  (i,j,k)=fs%IG  (i,fs%cfg%jmin,k)
               fs%U   (i,j,k)=fs%U   (i,fs%cfg%jmin,k)
               fs%V  (i,j,k)=min(fs%V(i,fs%cfg%jmin,k),0.0_WP) ! AS apply clipped Neumann for velocity in y+ to prevent inflow 
               !fs%V   (i,j,k)=fs%V   (i,fs%cfg%jmin,k)
               fs%W   (i,j,k)=fs%W   (i,fs%cfg%jmin,k)
               fs%VF  (i,j,k)=fs%VF  (i,fs%cfg%jmin,k)
               ! Also adjust interface data
               call setPlane(fs%PLIC(i,j,k),0,[0.0_WP,-1.0_WP,0.0_WP],-fs%cfg%y(j+1)+fs%dy*fs%VF(i,j,k))
               fs%BL(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]
               fs%BG(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]
            end do
         end do; end do
      end if
      
      ! AS: This is intended as a 2D case, we use periodic BCs in z (handeled in geometry.f90)
      ! ! Apply UNclipped Neumann on primitive variables in z+
      ! if (.not.fs%cfg%zper.and.fs%cfg%kproc.eq.fs%cfg%npz) then
      !    do j=fs%cfg%jmino_,fs%cfg%jmaxo_; do i=fs%cfg%imino_,fs%cfg%imaxo_
      !       ! Copy over from kmax to kmax+1 and above
      !       do k=fs%cfg%kmax+1,fs%cfg%kmaxo
      !          ! Copy primitive variables
      !          fs%RHOL(i,j,k)=fs%RHOL(i,j,fs%cfg%kmax)
      !          fs%PL  (i,j,k)=fs%PL  (i,j,fs%cfg%kmax)
      !          fs%IL  (i,j,k)=fs%IL  (i,j,fs%cfg%kmax)
      !          fs%RHOG(i,j,k)=fs%RHOG(i,j,fs%cfg%kmax)
      !          fs%PG  (i,j,k)=fs%PG  (i,j,fs%cfg%kmax)
      !          fs%IG  (i,j,k)=fs%IG  (i,j,fs%cfg%kmax)
      !          fs%U   (i,j,k)=fs%U   (i,j,fs%cfg%kmax)
      !          fs%V   (i,j,k)=fs%V   (i,j,fs%cfg%kmax)
      !          !fs%W  (i,j,k)=max(fs%W(i,j,fs%cfg%kmax),0.0_WP)
      !          fs%W   (i,j,k)=fs%W   (i,j,fs%cfg%kmax)
      !          fs%VF  (i,j,k)=fs%VF  (i,j,fs%cfg%kmax)
      !          ! Also adjust interface data
      !          call setPlane(fs%PLIC(i,j,k),0,[0.0_WP,0.0_WP,+1.0_WP],fs%cfg%z(k)+fs%dz*fs%VF(i,j,k))
      !          fs%BL(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]
      !          fs%BG(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]
      !       end do
      !    end do; end do
      ! end if
      
      ! ! Apply UNclipped Neumann on primitive variables in z-
      ! if (.not.fs%cfg%zper.and.fs%cfg%kproc.eq.1) then
      !    do j=fs%cfg%jmino_,fs%cfg%jmaxo_; do i=fs%cfg%imino_,fs%cfg%imaxo_
      !       ! First copy over W from kmin+1 to kmin
      !       fs%W(i,j,fs%cfg%kmin)=min(fs%W(i,j,fs%cfg%kmin+1),0.0_WP)
      !       ! Then copy over from kmin to kmin-1 and below
      !       do k=fs%cfg%kmino,fs%cfg%kmin-1
      !          ! Copy primitive variables
      !          fs%RHOL(i,j,k)=fs%RHOL(i,j,fs%cfg%kmin)
      !          fs%PL  (i,j,k)=fs%PL  (i,j,fs%cfg%kmin)
      !          fs%IL  (i,j,k)=fs%IL  (i,j,fs%cfg%kmin)
      !          fs%RHOG(i,j,k)=fs%RHOG(i,j,fs%cfg%kmin)
      !          fs%PG  (i,j,k)=fs%PG  (i,j,fs%cfg%kmin)
      !          fs%IG  (i,j,k)=fs%IG  (i,j,fs%cfg%kmin)
      !          fs%U   (i,j,k)=fs%U   (i,j,fs%cfg%kmin)
      !          fs%V   (i,j,k)=fs%V   (i,j,fs%cfg%kmin)
      !          !fs%W  (i,j,k)=min(fs%W(i,j,fs%cfg%kmin),0.0_WP)
      !          fs%W   (i,j,k)=fs%W   (i,j,fs%cfg%kmin)
      !          fs%VF  (i,j,k)=fs%VF  (i,j,fs%cfg%kmin)
      !          ! Also adjust interface data
      !          call setPlane(fs%PLIC(i,j,k),0,[0.0_WP,0.0_WP,-1.0_WP],-fs%cfg%z(k+1)+fs%dz*fs%VF(i,j,k))
      !          fs%BL(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]
      !          fs%BG(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]
      !       end do
      !    end do; end do
      ! end if
      
      ! Rebuild conserved quantities
      fs%Q(:,:,:,1)=        fs%VF *fs%RHOL
      fs%Q(:,:,:,2)=(1.0_WP-fs%VF)*fs%RHOG
      fs%Q(:,:,:,3)= fs%Q(:,:,:,1)*fs%IL
      fs%Q(:,:,:,4)= fs%Q(:,:,:,2)*fs%IG
      call fs%get_momentum()
      
   end subroutine apply_bconds

   ! !> Apply boundary conditions
   ! subroutine apply_bconds()
   !    use irl_fortran_interface, only: setPlane
   !    integer :: i,j,k,incoming
   !    real(WP) :: L_plus,L_minus,Pinf,Uinf,L_minus_target,L_minus_new,sigma
   !    real(WP) :: lambda_1,lambda_2,lambda_3,p_ext,rho_ext,u_ext,I_ext
   !    sigma=0.02_WP

   !    ! Apply characteristic nonreflecting condition at the outflow (attempt 2)
   !    if (.not.fs%cfg%xper.and.fs%cfg%iproc.eq.fs%cfg%npx) then
   !       do k=fs%cfg%kmino_,fs%cfg%kmaxo_; do j=fs%cfg%jmino_,fs%cfg%jmaxo_
   !          ! Copy over from imax to imax+1 and above
   !          do i=fs%cfg%imax+1,fs%cfg%imaxo
   !             lambda_1=fs%U(fs%cfg%imax,j,k)-fs%C(fs%cfg%imax,j,k)
   !             lambda_2=fs%U(fs%cfg%imax,j,k)
   !             lambda_3=fs%U(fs%cfg%imax,j,k)+fs%C(fs%cfg%imax,j,k)
   !             incoming=0
   !             if (lambda_1.le.0.0_WP) incoming=incoming+1
   !             if (lambda_2.le.0.0_WP) incoming=incoming+1
   !             if (lambda_3.le.0.0_WP) incoming=incoming+1
   !             p_ext=238557.89124999996_WP; rho_ext=2.1799412922956609_WP; u_ext=225.89432002347564_WP; I_ext=p_ext/(rho_ext*(GammaG-1.0_WP))
   !             if (fs%U(fs%cfg%imax,j,k).lt.0.5_WP*u_ext) then
   !                ! Copy primitive variables
   !                fs%RHOL(i,j,k)=fs%RHOL(fs%cfg%imax,j,k)
   !                fs%PL  (i,j,k)=fs%PL  (fs%cfg%imax,j,k)
   !                fs%IL  (i,j,k)=fs%IL  (fs%cfg%imax,j,k)
   !                fs%RHOG(i,j,k)=fs%RHOG(fs%cfg%imax,j,k)
   !                ! fs%PG  (i,j,k)=0.5_WP*fs%RHOG(fs%cfg%imax,j,k)*fs%C(fs%cfg%imax,j,k)*(L_plus+L_minus_new)
   !                fs%PG  (i,j,k)=fs%PG  (fs%cfg%imax,j,k)
   !                fs%IG  (i,j,k)=fs%IG  (fs%cfg%imax,j,k)
   !                ! fs%U   (i,j,k)=0.5_WP*(L_plus-L_minus_new)
   !                fs%U   (i,j,k)=fs%U   (fs%cfg%imax,j,k)
   !                fs%V   (i,j,k)=fs%V   (fs%cfg%imax,j,k)
   !                fs%W   (i,j,k)=fs%W   (fs%cfg%imax,j,k)
   !                fs%VF  (i,j,k)=fs%VF  (fs%cfg%imax,j,k)
   !             else 
   !                fs%RHOL(i,j,k)=fs%RHOL(fs%cfg%imax,j,k)
   !                fs%PL  (i,j,k)=fs%PL  (fs%cfg%imax,j,k)
   !                fs%IL  (i,j,k)=fs%IL  (fs%cfg%imax,j,k)
   !                fs%RHOG(i,j,k)=rho_ext
   !                fs%PG  (i,j,k)=p_ext
   !                fs%IG  (i,j,k)=I_ext
   !                fs%U   (i,j,k)=u_ext
   !                fs%V   (i,j,k)=fs%V   (fs%cfg%imax,j,k)
   !                fs%W   (i,j,k)=fs%W   (fs%cfg%imax,j,k)
   !                fs%VF  (i,j,k)=fs%VF  (fs%cfg%imax,j,k)
   !             end if
   !             ! Also adjust interface data
   !             call setPlane(fs%PLIC(i,j,k),0,[+1.0_WP,0.0_WP,0.0_WP],fs%cfg%x(i)+fs%dx*fs%VF(i,j,k))
   !             fs%BL(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]
   !             fs%BG(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]
   !          end do
   !       end do; end do
   !    end if
   !    ! Apply characteristic nonreflecting condition at the outflow (attempt 1)
   !    ! if (.not.fs%cfg%xper.and.fs%cfg%iproc.eq.fs%cfg%npx) then
   !    !    do k=fs%cfg%kmino_,fs%cfg%kmaxo_; do j=fs%cfg%jmino_,fs%cfg%jmaxo_
   !    !       ! Copy over from imax to imax+1 and above
   !    !       do i=fs%cfg%imax,fs%cfg%imaxo
   !    !          ! Get characteristic amplitudes from interior (for now simplify to just gas)
   !    !          L_plus =fs%PG(fs%cfg%imax-1,j,k)/(fs%RHOG(fs%cfg%imax-1,j,k)*fs%C(i,j,k))+fs%U(fs%cfg%imax-1,j,k)
   !    !          L_minus=fs%PG(fs%cfg%imax-1,j,k)/(fs%RHOG(fs%cfg%imax-1,j,k)*fs%C(i,j,k))-fs%U(fs%cfg%imax-1,j,k)
   !    !          ! Get target characteristic from farfield
   !    !          Pinf=238557.89124999996_WP; Uinf=225.89432002347564_WP
   !    !          L_minus_target=Pinf/(fs%RHOG(fs%cfg%imax,j,k)*fs%C(i,j,k))-Uinf
   !    !          ! Relax incomping amplitude to reduce reflections
   !    !          L_minus_new=(1.0_WP-sigma)*L_minus+sigma*L_minus_target
   !    !          ! Copy primitive variables
   !    !          fs%RHOL(i,j,k)=fs%RHOL(fs%cfg%imax-1,j,k)
   !    !          fs%PL  (i,j,k)=fs%PL  (fs%cfg%imax-1,j,k)
   !    !          fs%IL  (i,j,k)=fs%IL  (fs%cfg%imax-1,j,k)
   !    !          fs%RHOG(i,j,k)=fs%RHOG(fs%cfg%imax-1,j,k)
   !    !          fs%PG  (i,j,k)=0.5_WP*fs%RHOG(fs%cfg%imax-1,j,k)*fs%C(fs%cfg%imax-1,j,k)*(L_plus+L_minus_new)
   !    !          fs%IG  (i,j,k)=fs%IG  (fs%cfg%imax-1,j,k)
   !    !          fs%U   (i,j,k)=0.5_WP*(L_plus-L_minus_new)
   !    !          fs%V   (i,j,k)=fs%V   (fs%cfg%imax-1,j,k)
   !    !          fs%W   (i,j,k)=fs%W   (fs%cfg%imax-1,j,k)
   !    !          fs%VF  (i,j,k)=fs%VF  (fs%cfg%imax-1,j,k)
   !    !          ! Also adjust interface data
   !    !          call setPlane(fs%PLIC(i,j,k),0,[+1.0_WP,0.0_WP,0.0_WP],fs%cfg%x(i)+fs%dx*fs%VF(i,j,k))
   !    !          fs%BL(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]
   !    !          fs%BG(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]
   !    !       end do
   !    !    end do; end do
   !    ! end if
   !    ! ! Apply clipped Neumann on primitive variables in x-
   !    ! if (.not.fs%cfg%xper.and.fs%cfg%iproc.eq.1) then
   !    !    do k=fs%cfg%kmino_,fs%cfg%kmaxo_; do j=fs%cfg%jmino_,fs%cfg%jmaxo_
   !    !       ! First copy over U from imin+1 to imin
   !    !       fs%U(fs%cfg%imin,j,k)=min(fs%U(fs%cfg%imin+1,j,k),0.0_WP)
   !    !       ! Then copy over from imin to imin-1 and below
   !    !       do i=fs%cfg%imino,fs%cfg%imin-1
   !    !          ! Copy primitive variables
   !    !          fs%RHOL(i,j,k)=fs%RHOL(fs%cfg%imin,j,k)
   !    !          fs%PL  (i,j,k)=fs%PL  (fs%cfg%imin,j,k)
   !    !          fs%IL  (i,j,k)=fs%IL  (fs%cfg%imin,j,k)
   !    !          fs%RHOG(i,j,k)=fs%RHOG(fs%cfg%imin,j,k)
   !    !          fs%PG  (i,j,k)=fs%PG  (fs%cfg%imin,j,k)
   !    !          fs%IG  (i,j,k)=fs%IG  (fs%cfg%imin,j,k)
   !    !          fs%U   (i,j,k)=min(fs%U(fs%cfg%imin,j,k),0.0_WP)
   !    !          fs%V   (i,j,k)=fs%V   (fs%cfg%imin,j,k)
   !    !          fs%W   (i,j,k)=fs%W   (fs%cfg%imin,j,k)
   !    !          fs%VF  (i,j,k)=fs%VF  (fs%cfg%imin,j,k)
   !    !          ! Also adjust interface data
   !    !          call setPlane(fs%PLIC(i,j,k),0,[-1.0_WP,0.0_WP,0.0_WP],fs%cfg%x(i)+fs%dx*fs%VF(i,j,k))
   !    !          fs%BL(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]
   !    !          fs%BG(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]
   !    !       end do
   !    !    end do; end do
   !    ! end if
   !    ! Apply clipped Neumann on primitive variables in y+
   !    if (.not.fs%cfg%yper.and.fs%cfg%jproc.eq.fs%cfg%npy) then
   !       do k=fs%cfg%kmino_,fs%cfg%kmaxo_; do i=fs%cfg%imino_,fs%cfg%imaxo_ 
   !          ! Copy over from jmax to jmax+1 and above
   !          do j=fs%cfg%jmax+1,fs%cfg%jmaxo
   !             ! Copy primitive variables
   !             fs%RHOL(i,j,k)=fs%RHOL(i,fs%cfg%jmax,k)
   !             fs%PL  (i,j,k)=fs%PL  (i,fs%cfg%jmax,k)
   !             fs%IL  (i,j,k)=fs%IL  (i,fs%cfg%jmax,k)
   !             fs%RHOG(i,j,k)=fs%RHOG(i,fs%cfg%jmax,k)
   !             fs%PG  (i,j,k)=fs%PG  (i,fs%cfg%jmax,k)
   !             fs%IG  (i,j,k)=fs%IG  (i,fs%cfg%jmax,k)
   !             fs%U   (i,j,k)=fs%U   (i,fs%cfg%jmax,k)
   !             fs%V  (i,j,k)=max(fs%V(i,fs%cfg%jmax,k),0.0_WP)
   !             fs%W   (i,j,k)=fs%W   (i,fs%cfg%jmax,k)
   !             fs%VF  (i,j,k)=fs%VF  (i,fs%cfg%jmax,k)
   !             ! Also adjust interface data
   !             call setPlane(fs%PLIC(i,j,k),0,[0.0_WP,+1.0_WP,0.0_WP],fs%cfg%y(j)+fs%dy*fs%VF(i,j,k))
   !             fs%BL(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]
   !             fs%BG(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]
   !          end do
   !       end do; end do
   !    end if
   !    ! Apply clipped Neumann on primitive variables in y-
   !    if (.not.fs%cfg%yper.and.fs%cfg%jproc.eq.1) then
   !       do k=fs%cfg%kmino_,fs%cfg%kmaxo_; do i=fs%cfg%imino_,fs%cfg%imaxo_
   !          ! First copy over V from jmin+1 to jmin
   !          fs%V(i,fs%cfg%jmin,k)=min(fs%V(i,fs%cfg%jmin+1,k),0.0_WP)
   !          ! Then copy over from jmin to jmin-1 and below
   !          do j=fs%cfg%jmino,fs%cfg%jmin-1
   !             ! Copy primitive variables
   !             fs%RHOL(i,j,k)=fs%RHOL(i,fs%cfg%jmin,k)
   !             fs%PL  (i,j,k)=fs%PL  (i,fs%cfg%jmin,k)
   !             fs%IL  (i,j,k)=fs%IL  (i,fs%cfg%jmin,k)
   !             fs%RHOG(i,j,k)=fs%RHOG(i,fs%cfg%jmin,k)
   !             fs%PG  (i,j,k)=fs%PG  (i,fs%cfg%jmin,k)
   !             fs%IG  (i,j,k)=fs%IG  (i,fs%cfg%jmin,k)
   !             fs%U   (i,j,k)=fs%U   (i,fs%cfg%jmin,k)
   !             fs%V   (i,j,k)=min(fs%V(i,fs%cfg%jmin,k),0.0_WP)
   !             fs%W   (i,j,k)=fs%W   (i,fs%cfg%jmin,k)
   !             fs%VF  (i,j,k)=fs%VF  (i,fs%cfg%jmin,k)
   !             ! Also adjust interface data
   !             call setPlane(fs%PLIC(i,j,k),0,[0.0_WP,-1.0_WP,0.0_WP],-fs%cfg%y(j+1)+fs%dy*fs%VF(i,j,k))
   !             fs%BL(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]
   !             fs%BG(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]
   !          end do
   !       end do; end do
   !    end if
   !    ! Rebuild conserved quantities
   !    fs%Q(:,:,:,1)=        fs%VF *fs%RHOL
   !    fs%Q(:,:,:,2)=(1.0_WP-fs%VF)*fs%RHOG
   !    fs%Q(:,:,:,3)= fs%Q(:,:,:,1)*fs%IL
   !    fs%Q(:,:,:,4)= fs%Q(:,:,:,2)*fs%IG
   !    call fs%get_momentum()
      
   !    ! ! Apply clipped Neumann on primitive variables in x+
   !    ! if (.not.fs%cfg%xper.and.fs%cfg%iproc.eq.fs%cfg%npx) then
   !    !    do k=fs%cfg%kmino_,fs%cfg%kmaxo_; do j=fs%cfg%jmino_,fs%cfg%jmaxo_
   !    !       ! Copy over from imax to imax+1 and above
   !    !       do i=fs%cfg%imax+1,fs%cfg%imaxo
   !    !          ! Copy primitive variables
   !    !          fs%RHOL(i,j,k)=fs%RHOL(fs%cfg%imax,j,k)
   !    !          fs%PL  (i,j,k)=fs%PL  (fs%cfg%imax,j,k)
   !    !          fs%IL  (i,j,k)=fs%IL  (fs%cfg%imax,j,k)
   !    !          fs%RHOG(i,j,k)=fs%RHOG(fs%cfg%imax,j,k)
   !    !          fs%PG  (i,j,k)=fs%PG  (fs%cfg%imax,j,k)
   !    !          fs%IG  (i,j,k)=fs%IG  (fs%cfg%imax,j,k)
   !    !          ! fs%U  (i,j,k)=max(fs%U(fs%cfg%imax,j,k),0.0_WP)
   !    !          fs%U   (i,j,k)=fs%U(fs%cfg%imax,j,k)
   !    !          fs%V   (i,j,k)=fs%V   (fs%cfg%imax,j,k)
   !    !          fs%W   (i,j,k)=fs%W   (fs%cfg%imax,j,k)
   !    !          fs%VF  (i,j,k)=fs%VF  (fs%cfg%imax,j,k)
   !    !          ! Also adjust interface data
   !    !          call setPlane(fs%PLIC(i,j,k),0,[+1.0_WP,0.0_WP,0.0_WP],fs%cfg%x(i)+fs%dx*fs%VF(i,j,k))
   !    !          fs%BL(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]
   !    !          fs%BG(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]
   !    !       end do
   !    !       ! do i=fs%cfg%imax,fs%cfg%imaxo
   !    !       !    fs%U   (i,j,k)=fs%U(fs%cfg%imax-1,j,k)
   !    !       ! end do
   !    !    end do; end do
   !    ! end if

   !    ! ! ! Apply clipped Neumann on primitive variables in x-
   !    ! ! if (.not.fs%cfg%xper.and.fs%cfg%iproc.eq.1) then
   !    ! !    do k=fs%cfg%kmino_,fs%cfg%kmaxo_; do j=fs%cfg%jmino_,fs%cfg%jmaxo_
   !    ! !       ! First copy over U from imin+1 to imin
   !    ! !       fs%U(fs%cfg%imin,j,k)=min(fs%U(fs%cfg%imin+1,j,k),0.0_WP)
   !    ! !       ! Then copy over from imin to imin-1 and below
   !    ! !       do i=fs%cfg%imino,fs%cfg%imin-1
   !    ! !          ! Copy primitive variables
   !    ! !          fs%RHOL(i,j,k)=fs%RHOL(fs%cfg%imin,j,k)
   !    ! !          fs%PL  (i,j,k)=fs%PL  (fs%cfg%imin,j,k)
   !    ! !          fs%IL  (i,j,k)=fs%IL  (fs%cfg%imin,j,k)
   !    ! !          fs%RHOG(i,j,k)=fs%RHOG(fs%cfg%imin,j,k)
   !    ! !          fs%PG  (i,j,k)=fs%PG  (fs%cfg%imin,j,k)
   !    ! !          fs%IG  (i,j,k)=fs%IG  (fs%cfg%imin,j,k)
   !    ! !          fs%U   (i,j,k)=min(fs%U(fs%cfg%imin,j,k),0.0_WP)
   !    ! !          fs%V   (i,j,k)=fs%V   (fs%cfg%imin,j,k)
   !    ! !          fs%W   (i,j,k)=fs%W   (fs%cfg%imin,j,k)
   !    ! !          fs%VF  (i,j,k)=fs%VF  (fs%cfg%imin,j,k)
   !    ! !          ! Also adjust interface data
   !    ! !          call setPlane(fs%PLIC(i,j,k),0,[-1.0_WP,0.0_WP,0.0_WP],fs%cfg%x(i)+fs%dx*fs%VF(i,j,k))
   !    ! !          fs%BL(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]
   !    ! !          fs%BG(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]
   !    ! !       end do
   !    ! !    end do; end do
   !    ! ! end if
   !    ! ! Apply clipped Neumann on primitive variables in y+
   !    ! if (.not.fs%cfg%yper.and.fs%cfg%jproc.eq.fs%cfg%npy) then
   !    !    do k=fs%cfg%kmino_,fs%cfg%kmaxo_; do i=fs%cfg%imino_,fs%cfg%imaxo_ 
   !    !       ! Copy over from jmax to jmax+1 and above
   !    !       do j=fs%cfg%jmax+1,fs%cfg%jmaxo
   !    !          ! Copy primitive variables
   !    !          fs%RHOL(i,j,k)=fs%RHOL(i,fs%cfg%jmax,k)
   !    !          fs%PL  (i,j,k)=fs%PL  (i,fs%cfg%jmax,k)
   !    !          fs%IL  (i,j,k)=fs%IL  (i,fs%cfg%jmax,k)
   !    !          fs%RHOG(i,j,k)=fs%RHOG(i,fs%cfg%jmax,k)
   !    !          fs%PG  (i,j,k)=fs%PG  (i,fs%cfg%jmax,k)
   !    !          fs%IG  (i,j,k)=fs%IG  (i,fs%cfg%jmax,k)
   !    !          fs%U   (i,j,k)=fs%U   (i,fs%cfg%jmax,k)
   !    !          fs%V  (i,j,k)=max(fs%V(i,fs%cfg%jmax,k),0.0_WP)
   !    !          fs%W   (i,j,k)=fs%W   (i,fs%cfg%jmax,k)
   !    !          fs%VF  (i,j,k)=fs%VF  (i,fs%cfg%jmax,k)
   !    !          ! Also adjust interface data
   !    !          call setPlane(fs%PLIC(i,j,k),0,[0.0_WP,+1.0_WP,0.0_WP],fs%cfg%y(j)+fs%dy*fs%VF(i,j,k))
   !    !          fs%BL(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]
   !    !          fs%BG(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]
   !    !       end do
   !    !    end do; end do
   !    ! end if

   !    ! ! Apply clipped Neumann on primitive variables in y-
   !    ! if (.not.fs%cfg%yper.and.fs%cfg%jproc.eq.1) then
   !    !    do k=fs%cfg%kmino_,fs%cfg%kmaxo_; do i=fs%cfg%imino_,fs%cfg%imaxo_
   !    !       ! First copy over V from jmin+1 to jmin
   !    !       fs%V(i,fs%cfg%jmin,k)=min(fs%V(i,fs%cfg%jmin+1,k),0.0_WP)
   !    !       ! Then copy over from jmin to jmin-1 and below
   !    !       do j=fs%cfg%jmino,fs%cfg%jmin-1
   !    !          ! Copy primitive variables
   !    !          fs%RHOL(i,j,k)=fs%RHOL(i,fs%cfg%jmin,k)
   !    !          fs%PL  (i,j,k)=fs%PL  (i,fs%cfg%jmin,k)
   !    !          fs%IL  (i,j,k)=fs%IL  (i,fs%cfg%jmin,k)
   !    !          fs%RHOG(i,j,k)=fs%RHOG(i,fs%cfg%jmin,k)
   !    !          fs%PG  (i,j,k)=fs%PG  (i,fs%cfg%jmin,k)
   !    !          fs%IG  (i,j,k)=fs%IG  (i,fs%cfg%jmin,k)
   !    !          fs%U   (i,j,k)=fs%U   (i,fs%cfg%jmin,k)
   !    !          fs%V   (i,j,k)=min(fs%V(i,fs%cfg%jmin,k),0.0_WP)
   !    !          fs%W   (i,j,k)=fs%W   (i,fs%cfg%jmin,k)
   !    !          fs%VF  (i,j,k)=fs%VF  (i,fs%cfg%jmin,k)
   !    !          ! Also adjust interface data
   !    !          call setPlane(fs%PLIC(i,j,k),0,[0.0_WP,-1.0_WP,0.0_WP],-fs%cfg%y(j+1)+fs%dy*fs%VF(i,j,k))
   !    !          fs%BL(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]
   !    !          fs%BG(:,i,j,k)=[fs%cfg%xm(i),fs%cfg%ym(j),fs%cfg%zm(k)]
   !    !       end do
   !    !    end do; end do
   !    ! end if

   !    ! ! Rebuild conserved quantities
   !    ! fs%Q(:,:,:,1)=        fs%VF *fs%RHOL
   !    ! fs%Q(:,:,:,2)=(1.0_WP-fs%VF)*fs%RHOG
   !    ! fs%Q(:,:,:,3)= fs%Q(:,:,:,1)*fs%IL
   !    ! fs%Q(:,:,:,4)= fs%Q(:,:,:,2)*fs%IG
   !    ! call fs%get_momentum()
      
   ! end subroutine apply_bconds
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      ! Deallocate work arrays
      deallocate(dQdt,Ui,Vi,Wi,Ma,beta,visc)
   end subroutine simulation_final
   
   
end module simulation