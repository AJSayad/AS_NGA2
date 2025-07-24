!> Definition for a shockgen class (shock generator)
module shockgen_class
  use precision,         only: WP
  !use geometry,          only: cfg
  use config_class,      only: config
  use mpcomp_class,      only: mpcomp
  use timetracker_class, only: timetracker
  use ensight_class,     only: ensight
  !use surfmesh_class,    only: surfmesh
  !use cclabel_class,     only: cclabel
  use event_class,       only: event
  use monitor_class,     only: monitor
  implicit none
  private
  public :: sgen
  !> Equations of state
  real(WP) :: PinfL,GammaL,CvL
  real(WP) :: PinfG,GammaG,CvG
  
  !> sgen object
  type :: sgen
      !> Multiphase compressible flow solver and corresponding time tracker
      type(mpcomp),      public :: fs
      type(timetracker), public :: time
      !> Single config
      type(config), public :: cfg
      !> Ensight postprocessing
      !type(surfmesh) :: smesh
      type(ensight)  :: ens_out
      type(event)    :: ens_evt
      !> Simulation monitor file
      type(monitor) :: mfile,cflfile,consfile
      !> Private work arrays
      real(WP), dimension(:,:,:,:,:), allocatable :: dQdt
      real(WP), dimension(:,:,:)    , allocatable :: Ui,Vi,Wi,Ma,beta,visc
      !> Flow parameters move these to IC setup
      real(WP) :: Ms,Xstart,Xend
      real(WP) :: rho1,p1,u1,M1
      real(WP) :: rho2,p2,u2,M2
      real(WP) :: rho_ratio,c_ratio
      real(WP) :: rhoL,ML
      real(WP) :: ReG,viscG,viscL,visc_ratio
      !> Drop info (added just for log files, likely to be removed)
      real(WP), public :: ddrop                 ! drop diameter
      real(WP), public, dimension(3) :: dctr    ! drop center (x,y,z)
      !> shock profile
      integer  :: shock_index, nshock                                                      ! i index of shock and number of points to thicken shock
      real(WP), dimension(:), allocatable :: RHOG_global,IG_global,PG_global,Ui_global     ! global arrays to collect centerline values from each proc
      real(WP), dimension(:), allocatable :: RHOG_center,IG_center,PG_center,Ui_center     ! centerline arrays
      real(WP), dimension(:), allocatable :: RHOG_profile,IG_profile,PG_profile,Ui_profile ! shock profile arrays
      !> Flags
      logical :: dim_flag ! Flag for dimensional or nondimensional init
    contains
      procedure :: init  !> initialize sgen simulation
      procedure :: step  !> advance sgen simulation by one timestep
      procedure :: final !> finalize sgen simulation and transfer shock profile
  end type sgen
  
contains

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
   
  !> Thermo-mechanical relaxation model
  subroutine PT_relax(VF,Q)
    implicit none
    real(WP),                intent(inout) :: VF
    real(WP), dimension(1:), intent(inout) :: Q
    real(WP) :: PG,PL,ZG,ZL,Pint
    real(WP) :: a,b,d,coeffL,coeffG,Peq,VFeq
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
    ! Get equilibrium volume fraction
    VFeq=VF*((gammaL-1.0_WP)*Peq+2.0_WP*PL+coeffL)/((1.0_WP+gammaL)*Peq+coeffL)
    ! Adjust conserved quantities
    Q(3)=Q(3)-0.5_WP*(Pint+Peq)*(VFeq-VF)
    Q(4)=Q(4)+0.5_WP*(Pint+Peq)*(VFeq-VF)
    VF=VFeq
    ! ================= Second step for thermal relaxation =================
    ! Setup quadratic problem
    a=Q(1)*CvL+Q(2)*CvG
    b=Q(1)*CvL*(GammaL*PinfL+PinfG)+Q(2)*CvG*(GammaG*PinfG+PinfL)-sum(Q(3:4))*(Q(1)*CvL*(GammaL-1.0_WP)+Q(2)*CvG*(GammaG-1.0_WP))
    d=(Q(1)*CvL*GammaL+Q(2)*CvG*GammaG)*PinfL*PinfG-sum(Q(3:4))*(Q(1)*CvL*(GammaL-1.0_WP)*PinfG+Q(2)*CvG*(GammaG-1.0_WP)*PinfL)
    ! Get equilibrium pressure
    Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
    ! Get equilibrium volume fraction
    VFeq=Q(1)*CvL*(GammaL-1.0_WP)*(Peq+PinfG)/(Q(1)*CvL*(GammaL-1.0_WP)*(Peq+PinfG)+Q(2)*CvG*(GammaG-1.0_WP)*(Peq+PinfL))
    ! Clean up solution
    if (VFeq.lt.0.0_WP) then; VFeq=0.0_WP; Peq=max(Peq,-PinfL); end if
    if (VFeq.gt.1.0_WP) then; VFeq=1.0_WP; Peq=max(Peq,-PinfG); end if
    ! Adjust conserved quantities
    Q(3)=(       VFeq)*(Peq+GammaL*PinfL)/(GammaL-1.0_WP)
    Q(4)=(1.0_WP-VFeq)*(Peq+GammaG*PinfG)/(GammaG-1.0_WP)
    VF=VFeq
    ! Last debugging check... Probably should never happen...
    if (Peq.lt.-PinfG) print*,"****************** NEGATIVE PRESSURE! -","VFeq",VFeq,"Peq",Peq ! AS commented out time%t (duct tape fix)
  end subroutine PT_relax

  !> Initialization of the shock generator (sgen) simulation
  subroutine init(this,shockgen_group)
    use sgrid_class, only: cartesian,sgrid
    use mpi_f08,     only: MPI_Group
    use param,       only: param_read
    implicit none
    class(sgen), intent(inout)  :: this
    type(MPI_Group), intent(in) :: shockgen_group
    type(sgrid) :: grid
    integer  :: i,j,k,nx,ny,nz
    real(WP) :: Lx,Ly,Lz 
    real(WP) :: ddrop,CPD,D0X,D0Y,D0Z
    real(WP) :: dx,dy,dz
    logical  :: xper,yper,zper
    real(WP), dimension(:), allocatable :: x,y,z

    call param_read('Dimensional flag', this%dim_flag)

    !> generate 1D mesh
    create_grid: block
      ! set periodic BCs 
      xper=.true.; yper=.true.; zper=.true.
      ! Read in grid definition
      call param_read('Drop diameter',ddrop)
      call param_read('CPD',CPD)
      call param_read('D0X',D0X)

      ! compute physical domain length
      Lx = D0X*ddrop                   ! compute domain lengths
      nx = ceiling((CPD*Lx)/ddrop)     ! compute number of uniform cells
      dx = Lx/nx                       ! uniform mesh spacing
      dy = dx; ny = 1; dz = dx; nz = 1 ! 1D case

      ! allocate arrays
      allocate(x(nx+1));allocate(y(ny+1));allocate(z(nz+1))

      ! create simple uniform rectilinear mesh
      do i=1,nx+1
        x(i) = real(i-1,WP)*dx
      end do
      y(1) = -0.5_WP*dy; y(2) = 0.5_WP*dy ! 1D mesh
      z(1) = -0.5_WP*dz; z(2) = 0.5_WP*dz ! 1D mesh

      ! General serial grid object
      grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=xper,yper=yper,zper=zper,name='ShockGen')
    end block create_grid
    
    ! Create a config from that grid on our entire group
    create_cfg: block
      use parallel, only: group
      integer, dimension(3) :: partition
      ! Read in partition
      call param_read('Partition',partition,short='p'); partition(2) = 1; partition(3) = 1; ! manually overwrite partition in y and z directions
      ! Create partitioned grid
      this%cfg=config(grp=shockgen_group,decomp=partition,grid=grid)
    end block create_cfg

        ! Create multipgase compressible flow solver
    create_velocity_solver: block
      ! Initialize solver with required thermodynamic functions
      call this%fs%initialize(cfg=this%cfg,getPL=get_PL,getCL=get_CL,getPG=get_PG,getCG=get_CG,name='Compressible NS')
      ! Provide relaxation model
      this%fs%relax=>P_relax
      ! Provide entropy calculation functions
      this%fs%getSL=>get_SL; this%fs%getSG=>get_SG
      ! Provide temperature calculation functions
      this%fs%getTL=>get_TL; this%fs%getTG=>get_TG
    end block create_velocity_solver

    !> create timetracker
    create_timetracker:block
      use param, only: param_read
      real(WP) :: vshock ! velocity at which the shock moves through unshocked air
      ! initialize timetracker
      this%time=timetracker(amRoot=this%cfg%amRoot)
      call param_read('Max timestep size',this%time%dtmax)
      call param_read('Max cfl number',this%time%cflmax)
      this%time%dt=this%time%dtmax
      ! Set PinfG to zero
      PinfG=0.0_WP
      ! Read in Gammas
      call param_read('Liquid gamma',GammaL)
      call param_read('Gas gamma'   ,GammaG)
      ! Read in shock Mach number and location
      call param_read('Shock Mach number',this%Ms)
      call param_read('Shockgen start',this%Xstart)
      call param_read('Shockgen end',this%Xend)

      if (this%dim_flag.eqv.(.false.))then ! non dimensional
         ! First generate static shock with normalized pre-shock conditions
         this%M1=this%Ms
         this%rho1=1.0_WP
         this%rho2=this%rho1*(GammaG+1.0_WP)*this%M1**2/((GammaG-1.0_WP)*this%M1**2+2.0_WP)
         this%p1=0.25_WP*this%rho1/GammaG*((GammaG+1.0_WP)*this%M1/(this%M1**2-1.0_WP))**2 ! Ensures that |u2-u1|=1
         this%p2=this%p1*(2.0_WP*GammaG/(GammaG+1.0_WP)*(this%M1**2-1.0_WP)+1.0_WP)
         this%u1=this%M1*sqrt(GammaG*this%p1/this%rho1); ! shock speed through unshocked air
         ! calculate final shock generator time based on how far the shock moves
         this%time%tmax = (this%Xend - this%Xstart)/this%u1
      else ! dimensional
         ! Read in Gas variables
         call param_read('Pre-shock density',this%rho1)
         call param_read('Pre-shock pressure',this%p1)
         ! Use shock relations (shock fixed frame) to calculate post-shock conditions 
         this%M1 = this%Ms ! set M1 equal to shock Mach number for now 
         this%rho2=this%rho1*(GammaG+1.0_WP)*this%M1**2/((GammaG-1.0_WP)*this%M1**2+2.0_WP) ! post shock density  (Anderson 3.53)
         this%p2=this%p1*(2.0_WP*GammaG/(GammaG+1.0_WP)*(this%M1**2-1.0_WP)+1.0_WP)         ! post shock pressure (Anderson 3.57)
         this%u1=this%M1*sqrt(GammaG*this%p1/this%rho1)                                     ! shock speed through unshocked air
         ! calculate final shock generator time based on how far the shock moves
         this%time%tmax = (this%Xend - this%Xstart)/this%u1         
      end if               
    end block create_timetracker

    !> initialize EOS and flow parameters
    initialize_parameters: block
      use string,   only: str_long
      use messager, only: log
      character(str_long) :: message
      ! Set PinfG to zero
      PinfG=0.0_WP
      ! Read in Gammas
      call param_read('Liquid gamma',GammaL)
      call param_read('Gas gamma'   ,GammaG)
      ! Read in shock Mach number and location
      call param_read('Shock Mach number',this%Ms)
      call param_read('Shockgen start',this%Xstart)
      if (this%dim_flag.eqv.(.false.))then ! non dimensional
        ! First generate static shock with normalized pre-shock conditions
        this%M1=this%Ms
        this%rho1=1.0_WP
        this%rho2=this%rho1*(GammaG+1.0_WP)*this%M1**2/((GammaG-1.0_WP)*this%M1**2+2.0_WP)
        this%p1=0.25_WP*this%rho1/GammaG*((GammaG+1.0_WP)*this%M1/(this%M1**2-1.0_WP))**2 ! Ensures that |u2-u1|=1
        this%p2=this%p1*(2.0_WP*GammaG/(GammaG+1.0_WP)*(this%M1**2-1.0_WP)+1.0_WP)
        this%u1=this%M1*sqrt(GammaG*this%p1/this%rho1)
        this%u2=this%u1*this%rho1/this%rho2
        ! Now shift frame of reference to obtain moving shock
        this%u2=abs(this%u2-this%u1); this%M2=this%u2/sqrt(GammaG*this%p2/this%rho2); this%u1=0.0_WP; this%M1=this%u1/sqrt(GammaG*this%p1/this%rho1)
        ! Read in density ratio and use it to set liquid density
        call param_read('Density ratio',this%rho_ratio); this%rhoL=this%rho_ratio*this%rho1
        ! Read in sound speed ratio and use it to set PinfL
        call param_read('Sound speed ratio',this%c_ratio)
        PinfL=this%p1*(this%rho_ratio*this%c_ratio**2*GammaG/GammaL-1.0_WP)
        this%ML=this%u2/sqrt(GammaL*(this%p1+PinfL)/this%rhoL)
        ! Set heat capacities corresponding to a normalized pre-shock and liquid temperature
        CvL=(this%p1+PinfL)/(this%rhoL*(GammaL-1.0_WP))
        CvG=(this%p1+PinfG)/(this%rho1*(GammaG-1.0_WP))
        ! Viscous parameters
        call param_read('Gas Reynolds number',this%ReG); this%viscG=this%rho1*1.0_WP*this%u2/this%ReG 
        call param_read('Viscosity ratio',this%visc_ratio); this%viscL=this%visc_ratio*this%viscG
      else ! dimensional
        ! Read in Liquid variables
        call param_read('Drop diameter',this%ddrop) ! AS this might be moved in the future
        call param_read('Drop center',this%dctr)    ! AS this might be moved in the future
        call param_read('Liquid Pinf',PinfL)
        call param_read('Liquid density',this%rhoL)
        call param_read('Liquid dynamic viscosity',this%viscL)
        call param_read('Liquid specific heat (constant vol)',CvL)

        ! Read in Gas variables
        call param_read('Pre-shock density',this%rho1)
        call param_read('Pre-shock pressure',this%p1)
        call param_read('Gas dynamic viscosity',this%viscG)
        call param_read('Gas specific heat (constant vol)',CvG)

        ! Use shock relations (shock fixed frame) to calculate post-shock conditions 
        this%M1 = this%Ms ! set M1 equal to shock Mach number for now 
        this%rho2=this%rho1*(GammaG+1.0_WP)*this%M1**2/((GammaG-1.0_WP)*this%M1**2+2.0_WP) ! post shock density  (Anderson 3.53)
        this%p2=this%p1*(2.0_WP*GammaG/(GammaG+1.0_WP)*(this%M1**2-1.0_WP)+1.0_WP)         ! post shock pressure (Anderson 3.57)
        this%u1=this%M1*sqrt(GammaG*this%p1/this%rho1)                                     ! velocity in state 1 (left side of shock in fixed frame)
        this%u2=this%u1*this%rho1/this%rho2                                                ! velocity in state 2 (Anderson 3.53, right side of shock in fixed frame)
        ! we now shift to accomodate a moving shock (converting from shock fixed frame to lab frame)
        this%u2=abs(this%u2-this%u1); this%M2=this%u2/sqrt(GammaG*this%p2/this%rho2) ! post-shock gas velocity and post shock Mach number
        this%u1=0.0_WP; this%M1=this%u1/sqrt(GammaG*this%p1/this%rho1)               ! set pre-shock gas velocity to zero and update pre-shock Mach number

        ! compute some non-dimensional parameters for log files
        this%rho_ratio=this%rho2/this%rho1       ! density ratio
        this%visc_ratio=this%viscL/this%viscG    ! viscosity ratio
        this%ReG = this%rho2*this%u2*this%ddrop/this%viscG ! Reynolds number based on post-shock conditions 
        this%c_ratio=sqrt(GammaL*(this%p1+PinfL)/this%rhoL)/sqrt(GammaG*this%p1/this%rho1) ! sound speed ratio
        this%ML=this%u2/sqrt(GammaL*(this%p1+PinfL)/this%rhoL)
      end if
      ! Output case info
      if (this%cfg%amRoot) then
        write(message,'("[Liquid EOS] => Gamma=",es12.5)') GammaL; call log(message)
        write(message,'("[Liquid EOS] =>  Pinf=",es12.5)')  PinfL; call log(message)
        write(message,'("[Liquid EOS] =>    Cv=",es12.5)')    CvL; call log(message)
        write(message,'("[Gas EOS]    => Gamma=",es12.5)') GammaG; call log(message)
        write(message,'("[Gas EOS]    =>    Cv=",es12.5)')    CvG; call log(message)
        write(message,'("[Shock Mach number]     =>     Ms=",es12.5)')     this%Ms; call log(message)
        write(message,'("[Pre -shock conditions] =>   rho1=",es12.5)')   this%rho1; call log(message)
        write(message,'("[Pre -shock conditions] =>     p1=",es12.5)')     this%p1; call log(message)
        write(message,'("[Pre -shock conditions] =>     u1=",es12.5)')     this%u1; call log(message)
        write(message,'("[Pre -shock conditions] =>     M1=",es12.5)')     this%M1; call log(message)
        write(message,'("[Post-shock conditions] =>   rho2=",es12.5)')   this%rho2; call log(message)
        write(message,'("[Post-shock conditions] =>     p2=",es12.5)')     this%p2; call log(message)
        write(message,'("[Post-shock conditions] =>     u2=",es12.5)')     this%u2; call log(message)
        write(message,'("[Post-shock conditions] =>     M2=",es12.5)')     this%M2; call log(message)
        write(message,'("[Liquid Mach number] =>        ML=",es12.5)')     this%Ml; call log(message)
        write(message,'("[Density ratio]      => rhoL/rho1=",es12.5)') this%rho_ratio; call log(message)
        write(message,'("[Sound speed ratio]  =>     cl/c1=",es12.5)')   this%c_ratio; call log(message)
        write(message,'("[Gas Reynolds]     =>     ReG=",es12.5)')        this%ReG; call log(message)
        write(message,'("[Viscosity ratio]  => muL/muG=",es12.5)') this%visc_ratio; call log(message)
        write(message,'("[Gas    viscosity] =>     muG=",es12.5)')      this%viscG; call log(message)
        write(message,'("[Liquid viscosity] =>     muL=",es12.5)')      this%viscL; call log(message)
      end if
    end block initialize_parameters

    ! Allocate work arrays
    allocate_work_arrays: block
      allocate(this%dQdt(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:this%fs%nQ,1:4))
      allocate(this%beta(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(this%visc(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(this%Ui(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(this%Vi(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(this%Wi(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(this%Ma(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
    end block allocate_work_arrays
    
    ! Prepare initial conditions
    initial_conditions: block
      use irl_fortran_interface, only: setNumberOfPlanes,setPlane
      use mms_geom,              only: initialize_volume_moments
      use mpcomp_class,          only: VFlo
      integer :: i,j,k
      ! Initialize primary variables
      this%fs%VF=0.0_WP;this%fs%RHOL=1.0_WP; this%fs%PL=1.0_WP; this%fs%IL=1.0_WP! singlephase, no need for levelset, set VOF=0 everywhere, liquid props=1.0
      do k=this%cfg%kmino_,this%cfg%kmaxo_
        do j=this%cfg%jmino_,this%cfg%jmaxo_
          do i=this%cfg%imino_,this%cfg%imaxo_
            ! even in singlephase, init barycenters and set # of planes (as instructed by Chase, related to BCs)
            this%fs%BL(:,i,j,k)=[this%fs%cfg%xm(i),this%fs%cfg%ym(j),this%fs%cfg%zm(k)]; this%fs%BG(:,i,j,k)=[this%fs%cfg%xm(i),this%fs%cfg%ym(j),this%fs%cfg%zm(k)]
            call setNumberOfPlanes(this%fs%PLIC(i,j,k),1); call setPlane(this%fs%PLIC(i,j,k),0,[0.0_WP,0.0_WP,0.0_WP],sign(1.0_WP,this%fs%VF(i,j,k)-0.5_WP))

            this%fs%U(i,j,k) = 0.0_WP; this%fs%V(i,j,k)=0.0_WP; this%fs%W(i,j,k)=0.0_WP ! set post shock gas velocity
            if(this%cfg%xm(i).lt.this%Xstart)then ! post-shock region
              this%fs%U   (i,j,k) = this%u2   ! velocity
              this%fs%RHOG(i,j,k) = this%rho2 ! density
              this%fs%PG  (i,j,k) = this%p2   ! pressure 
              this%fs%IG  (i,j,k) = (this%fs%PG(i,j,k)+GammaG*PinfG)/(this%fs%RHOG(i,j,k)*(GammaG-1.0_WP)) ! internal energy 
            else ! pre-shock region
              this%fs%U   (i,j,k) = this%u1   ! velocity (should be zero for most cases)
              this%fs%RHOG(i,j,k) = this%rho1 ! density
              this%fs%PG  (i,j,k) = this%p1   ! pressure 
              this%fs%IG  (i,j,k) = (this%fs%PG(i,j,k)+GammaG*PinfG)/(this%fs%RHOG(i,j,k)*(GammaG-1.0_WP)) ! internal energy
            end if
          end do
        end do
      end do

      ! AS I'm pretty sure we don't need to call build interface
      ! Build PLIC interface
      !call fs%build_interface()
      ! Initialize conserved variables
      this%fs%Q(:,:,:,1)=         this%fs%VF*this%fs%RHOL
      this%fs%Q(:,:,:,2)=(1.0_WP-this%fs%VF)*this%fs%RHOG
      this%fs%Q(:,:,:,3)= this%fs%Q(:,:,:,1)*this%fs%IL
      this%fs%Q(:,:,:,4)= this%fs%Q(:,:,:,2)*this%fs%IG
      call this%fs%get_momentum()
      ! Communicate conserved variables (not needed in general, but allows 2D runs without changing loop above...)
      do i=1,this%fs%nQ; call this%fs%cfg%sync(this%fs%Q(:,:,:,i)); end do
      ! Rebuild primitive variables
      call this%fs%get_primitive()
      ! Interpolate velocity
      call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
      ! Compute local Mach number
      this%Ma=sqrt(this%Ui**2+this%Vi**2+this%Wi**2)/this%fs%C
    end block initial_conditions

    ! Add Ensight output
    create_ensight: block
      ! Create Ensight output from cfg
      this%ens_out=ensight(cfg=this%cfg,name='ShockGen')
      ! Create event for Ensight output
      this%ens_evt=event(time=this%time,name='Ensight output')
      call param_read('Ensight output period',this%ens_evt%tper)
      ! Add variables to output
      call this%ens_out%add_vector('velocity',this%Ui,this%Vi,this%Wi)
      call this%ens_out%add_scalar('VOF',this%fs%VF)
      call this%ens_out%add_scalar('RHOL',this%fs%RHOL)
      call this%ens_out%add_scalar('RHOG',this%fs%RHOG)
      call this%ens_out%add_scalar('IL',this%fs%IL)
      call this%ens_out%add_scalar('IG',this%fs%IG)
      call this%ens_out%add_scalar('PL',this%fs%PL)
      call this%ens_out%add_scalar('PG',this%fs%PG)
      call this%ens_out%add_scalar('Mach',this%Ma)
      call this%ens_out%add_scalar('beta',this%beta)
      call this%ens_out%add_scalar('visc',this%visc)
      call this%ens_out%add_scalar('TL',this%fs%TL)
      call this%ens_out%add_scalar('TG',this%fs%TG)
      ! AS debugging
      call this%ens_out%add_scalar('sound_speed',this%fs%C) 
      ! ! Create surface mesh for PLIC
      ! this%smesh=surfmesh(nvar=0,name='plic')
      ! call this%fs%update_surfmesh(smesh)
      ! call this%ens_out%add_surface('plic',this%smesh)
      ! Output to ensight
      if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
    end block create_ensight

    ! Create monitor files
    create_monitor: block
      ! Prepare some info about fields
      call this%fs%get_cfl(dt=this%time%dt,cfl=this%time%cfl)
      call this%fs%get_info()
      ! Create simulation monitor
      this%mfile=monitor(this%fs%cfg%amRoot,'simulation')
      call this%mfile%add_column(this%time%n,'Timestep number')
      call this%mfile%add_column(this%time%t,'Time')
      call this%mfile%add_column(this%time%dt,'Timestep size')
      call this%mfile%add_column(this%time%cfl,'Maximum CFL')
      call this%mfile%add_column(this%fs%Umax,'Umax')
      call this%mfile%add_column(this%fs%Vmax,'Vmax')
      call this%mfile%add_column(this%fs%Wmax,'Wmax')
      call this%mfile%add_column(this%fs%RHOLmax,'max(RHOL)')
      call this%mfile%add_column(this%fs%RHOLmin,'min(RHOL)')
      call this%mfile%add_column(this%fs%ILmax  ,'max(IL)'  )
      call this%mfile%add_column(this%fs%ILmin  ,'min(IL)'  )
      call this%mfile%add_column(this%fs%PLmax  ,'max(PL)'  )
      call this%mfile%add_column(this%fs%PLmin  ,'min(PL)'  )
      call this%mfile%add_column(this%fs%TLmax  ,'max(TL)'  )
      call this%mfile%add_column(this%fs%TLmin  ,'min(TL)'  )
      call this%mfile%add_column(this%fs%RHOGmax,'max(RHOG)')
      call this%mfile%add_column(this%fs%RHOGmin,'min(RHOG)')
      call this%mfile%add_column(this%fs%IGmax  ,'max(IG)'  )
      call this%mfile%add_column(this%fs%IGmin  ,'min(IG)'  )
      call this%mfile%add_column(this%fs%PGmax  ,'max(PG)'  )
      call this%mfile%add_column(this%fs%PGmin  ,'min(PG)'  )
      call this%mfile%add_column(this%fs%TGmax  ,'max(TG)'  )
      call this%mfile%add_column(this%fs%TGmin  ,'min(TG)'  )
      call this%mfile%add_column(this%fs%VFmax  ,'VFmax'    )
      call this%mfile%add_column(this%fs%VFmin  ,'VFmin'    )
      call this%mfile%write()
      ! Create CFL monitor
      this%cflfile=monitor(this%fs%cfg%amRoot,'cfl')
      call this%cflfile%add_column(this%time%n,'Timestep number')
      call this%cflfile%add_column(this%time%t,'Time')
      call this%cflfile%add_column(this%fs%CFLc_x,'Convective xCFL')
      call this%cflfile%add_column(this%fs%CFLc_y,'Convective yCFL')
      call this%cflfile%add_column(this%fs%CFLc_z,'Convective zCFL')
      call this%cflfile%add_column(this%fs%CFLa_x,'Acoustic xCFL')
      call this%cflfile%add_column(this%fs%CFLa_y,'Acoustic yCFL')
      call this%cflfile%add_column(this%fs%CFLa_z,'Acoustic zCFL')
      call this%cflfile%add_column(this%fs%CFLv_x,'Viscous xCFL')
      call this%cflfile%add_column(this%fs%CFLv_y,'Viscous yCFL')
      call this%cflfile%add_column(this%fs%CFLv_z,'Viscous zCFL')
      call this%cflfile%write()
      ! Create conservation monitor
      this%consfile=monitor(this%fs%cfg%amRoot,'conservation')
      call this%consfile%add_column(this%time%n,'Timestep number')
      call this%consfile%add_column(this%time%t,'Time')
      call this%consfile%add_column(this%fs%VFint  ,'Volume')
      call this%consfile%add_column(this%fs%Qint(1),'Liquid mass')
      call this%consfile%add_column(this%fs%Qint(2),'Gas mass')
      call this%consfile%add_column(this%fs%Qint(3),'Liquid energy')
      call this%consfile%add_column(this%fs%Qint(4),'Gas energy')
      call this%consfile%add_column(this%fs%Qint(5),'U Momentum')
      call this%consfile%add_column(this%fs%Qint(6),'V Momentum')
      call this%consfile%add_column(this%fs%Qint(7),'W Momentum')
      call this%consfile%add_column(this%fs%RHOKLint,'Liquid KE')
      call this%consfile%add_column(this%fs%RHOKGint,'Gas KE')
      call this%consfile%add_column(this%fs%RHOSLint,'Liquid entropy')
      call this%consfile%add_column(this%fs%RHOSGint,'Gas entropy')
      call this%consfile%write()
    end block create_monitor
  end subroutine init

  !> take one step with specified dt
  subroutine step(this)
    implicit none
    class(sgen), intent(inout) :: this
    ! do while (.not. this%time%done()) ! AS this loop needs to be moved to sim.f90 eventually
    ! Increment time
    call this%fs%get_cfl(dt=this%time%dt,cfl=this%time%cfl)
    call this%time%adjust_dt()
    call this%time%increment()
    
    ! Remember conserved variables
    this%fs%Qold=this%fs%Q
    
    ! Remember phasic quantities
    this%fs%RHOLold=this%fs%RHOL; this%fs%ILold=this%fs%IL; this%fs%PLold=this%fs%PL
    this%fs%RHOGold=this%fs%RHOG; this%fs%IGold=this%fs%IG; this%fs%PGold=this%fs%PG
    
    ! Remember volume moments and interface
    this%fs%VFold = this%fs%VF
    this%fs%BLold = this%fs%BL
    this%fs%BGold = this%fs%BG
    
    copy_plic_to_old: block
      use irl_fortran_interface, only: copy
      integer :: i,j,k
      do k=this%fs%cfg%kmino_,this%fs%cfg%kmaxo_; do j=this%fs%cfg%jmino_,this%fs%cfg%jmaxo_; do i=this%fs%cfg%imino_,this%fs%cfg%imaxo_
         call copy(this%fs%PLICold(i,j,k),this%fs%PLIC(i,j,k))
      end do; end do; end do
    end block copy_plic_to_old
    
    ! Tag cells for semi-Lagrangian transport
    call this%fs%SLtag()
    
    ! Prepare SGS viscosity models
    call this%fs%get_viscartif(dt=this%time%dt,beta=this%beta)
    call this%fs%get_vreman   (dt=this%time%dt,visc=this%visc)
    mixture_viscosity: block
      integer  :: i,j,k
      real(WP) :: Lvof,Lrho,Gvof,Grho
      real(WP) :: Lvisc,Gvisc,Lbeta,Gbeta
      real(WP), parameter :: eps=1.0e-15_WP
      do k=this%fs%cfg%kmino_+1,this%fs%cfg%kmaxo_-1; do j=this%fs%cfg%jmino_+1,this%fs%cfg%jmaxo_-1; do i=this%fs%cfg%imino_+1,this%fs%cfg%imaxo_-1
         ! Create smooth mass info distribution
         Lvof=sum(       this%fs%VF(i-1:i+1,j-1:j+1,k-1:k+1)  )
         Gvof=sum(1.0_WP-this%fs%VF(i-1:i+1,j-1:j+1,k-1:k+1)  )
         Lrho=sum(       this%fs%Q (i-1:i+1,j-1:j+1,k-1:k+1,1))/(Lvof+eps)
         Grho=sum(       this%fs%Q (i-1:i+1,j-1:j+1,k-1:k+1,2))/(Gvof+eps)
         ! Harmonic average of VISC
         Lvisc=Lrho*(this%viscL+this%visc(i,j,k)); Gvisc=Grho*(this%viscG+this%visc(i,j,k)); this%fs%VISC(i,j,k)=(Lvof+Gvof)/(Lvof/max(Lvisc,eps)+Gvof/max(Gvisc,eps))
         ! Harmonic average of BETA
         Lbeta=Lrho*this%beta(i,j,k); Gbeta=Grho*this%beta(i,j,k); this%fs%BETA(i,j,k)=(Lvof+Gvof)/(Lvof/max(Lbeta,eps)+Gvof/max(Gbeta,eps))
      end do; end do; end do
    end block mixture_viscosity
    
    ! Perform first semi-Lagrangian transport step =====================================================
    call this%fs%SLstep(dt=0.5_WP*this%time%dt,U=this%fs%U,V=this%fs%V,W=this%fs%W)
    call this%fs%build_interface()
    
    ! First RK step ====================================================================================
    ! Get non-SL RHS and increment
    call this%fs%rhs(this%dQdt(:,:,:,:,1))
    this%fs%Q=this%fs%Qold+0.5_WP*this%time%dt*this%dQdt(:,:,:,:,1)
    ! Increment Q with SL terms
    this%fs%Q=this%fs%Q+this%fs%SLdQ
    ! Recompute primitive variables
    call this%fs%get_primitive()
    
    ! Second RK step ===================================================================================
    ! Get non-SL RHS and increment
    call this%fs%rhs(this%dQdt(:,:,:,:,2))
    this%fs%Q=this%fs%Qold+0.5_WP*this%time%dt*this%dQdt(:,:,:,:,2)
    ! Increment Q with SL terms
    this%fs%Q=this%fs%Q+this%fs%SLdQ
    ! Recompute primitive variables
    call this%fs%get_primitive()
    
    ! Perform second semi-Lagrangian transport step ====================================================
    call this%fs%SLstep(dt=1.0_WP*this%time%dt,U=this%fs%U,V=this%fs%V,W=this%fs%W)
    call this%fs%build_interface()
    
    ! Third RK step ====================================================================================
    ! Get non-SL RHS and increment
    call this%fs%rhs(this%dQdt(:,:,:,:,3))
    this%fs%Q=this%fs%Qold+1.0_WP*this%time%dt*this%dQdt(:,:,:,:,3)
    ! Increment Q with SL terms
    this%fs%Q=this%fs%Q+this%fs%SLdQ
    ! Recompute primitive variables
    call this%fs%get_primitive()
    
    ! Fourth RK step ===================================================================================
    ! Get non-SL RHS and increment
    call this%fs%rhs(this%dQdt(:,:,:,:,4))
    this%fs%Q=this%fs%Qold+this%time%dt/6.0_WP*(this%dQdt(:,:,:,:,1)+2.0_WP*this%dQdt(:,:,:,:,2)+2.0_WP*this%dQdt(:,:,:,:,3)+this%dQdt(:,:,:,:,4))
    ! Increment Q with SL terms
    this%fs%Q=this%fs%Q+this%fs%SLdQ
    ! Apply user-provided relaxation model
    call this%fs%apply_relax()
    ! Recompute primitive variables
    call this%fs%get_primitive()
    ! Apply Neumann condition at the outflow
    neumann_outflow: block
      use irl_fortran_interface, only: setPlane
      integer :: i,j,k
      ! Apply clipped Neumann on primitive variables in x+
      if (.not.this%cfg%xper.and.this%cfg%iproc.eq.this%cfg%npx) then
         do k=this%cfg%kmino_,this%cfg%kmaxo_; do j=this%cfg%jmino_,this%cfg%jmaxo_
            ! Copy over from imax to imax+1 and above
            do i=this%cfg%imax+1,this%cfg%imaxo
               ! Copy primitive variables
               this%fs%RHOL(i,j,k)=this%fs%RHOL(this%cfg%imax,j,k)
               this%fs%PL  (i,j,k)=this%fs%PL  (this%cfg%imax,j,k)
               this%fs%IL  (i,j,k)=this%fs%IL  (this%cfg%imax,j,k)
               this%fs%RHOG(i,j,k)=this%fs%RHOG(this%cfg%imax,j,k)
               this%fs%PG  (i,j,k)=this%fs%PG  (this%cfg%imax,j,k)
               this%fs%IG  (i,j,k)=this%fs%IG  (this%cfg%imax,j,k)
               this%fs%U  (i,j,k)=max(this%fs%U(this%cfg%imax,j,k),0.0_WP)
               this%fs%V   (i,j,k)=this%fs%V   (this%cfg%imax,j,k)
               this%fs%W   (i,j,k)=this%fs%W   (this%cfg%imax,j,k)
               this%fs%VF  (i,j,k)=this%fs%VF  (this%cfg%imax,j,k)
               ! Also adjust interface data
               call setPlane(this%fs%PLIC(i,j,k),0,[+1.0_WP,0.0_WP,0.0_WP],this%cfg%x(i)+this%fs%dx*this%fs%VF(i,j,k))
               this%fs%BL(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
               this%fs%BG(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
            end do
         end do; end do
      end if
      ! Apply clipped Neumann on primitive variables in y+
      if (.not.this%cfg%yper.and.this%cfg%jproc.eq.this%cfg%npy) then
         do k=this%cfg%kmino_,this%cfg%kmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
            ! Copy over from jmax to jmax+1 and above
            do j=this%cfg%jmax+1,this%cfg%jmaxo
               ! Copy primitive variables
               this%fs%RHOL(i,j,k)=this%fs%RHOL(i,this%cfg%jmax,k)
               this%fs%PL  (i,j,k)=this%fs%PL  (i,this%cfg%jmax,k)
               this%fs%IL  (i,j,k)=this%fs%IL  (i,this%cfg%jmax,k)
               this%fs%RHOG(i,j,k)=this%fs%RHOG(i,this%cfg%jmax,k)
               this%fs%PG  (i,j,k)=this%fs%PG  (i,this%cfg%jmax,k)
               this%fs%IG  (i,j,k)=this%fs%IG  (i,this%cfg%jmax,k)
               this%fs%U   (i,j,k)=this%fs%U   (i,this%cfg%jmax,k)
               this%fs%V  (i,j,k)=max(this%fs%V(i,this%cfg%jmax,k),0.0_WP)
               this%fs%W   (i,j,k)=this%fs%W   (i,this%cfg%jmax,k)
               this%fs%VF  (i,j,k)=this%fs%VF  (i,this%cfg%jmax,k)
               ! Also adjust interface data
               call setPlane(this%fs%PLIC(i,j,k),0,[0.0_WP,+1.0_WP,0.0_WP],this%cfg%y(j)+this%fs%dy*this%fs%VF(i,j,k))
               this%fs%BL(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
               this%fs%BG(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
            end do
         end do; end do
      end if
      ! Apply clipped Neumann on primitive variables in y-
      if (.not.this%cfg%yper.and.this%cfg%jproc.eq.1) then
         do k=this%cfg%kmino_,this%cfg%kmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
            ! First copy over V from jmin+1 to jmin
            this%fs%V(i,this%cfg%jmin,k)=min(this%fs%V(i,this%cfg%jmin+1,k),0.0_WP)
            ! Then copy over from jmin to jmin-1 and below
            do j=this%cfg%jmino,this%cfg%jmin-1
               ! Copy primitive variables
               this%fs%RHOL(i,j,k)=this%fs%RHOL(i,this%cfg%jmin,k)
               this%fs%PL  (i,j,k)=this%fs%PL  (i,this%cfg%jmin,k)
               this%fs%IL  (i,j,k)=this%fs%IL  (i,this%cfg%jmin,k)
               this%fs%RHOG(i,j,k)=this%fs%RHOG(i,this%cfg%jmin,k)
               this%fs%PG  (i,j,k)=this%fs%PG  (i,this%cfg%jmin,k)
               this%fs%IG  (i,j,k)=this%fs%IG  (i,this%cfg%jmin,k)
               this%fs%U   (i,j,k)=this%fs%U   (i,this%cfg%jmin,k)
               this%fs%V  (i,j,k)=min(this%fs%V(i,this%cfg%jmin,k),0.0_WP)
               this%fs%W   (i,j,k)=this%fs%W   (i,this%cfg%jmin,k)
               this%fs%VF  (i,j,k)=this%fs%VF  (i,this%cfg%jmin,k)
               ! Also adjust interface data
               call setPlane(this%fs%PLIC(i,j,k),0,[0.0_WP,-1.0_WP,0.0_WP],this%cfg%y(j)+this%fs%dy*this%fs%VF(i,j,k))
               this%fs%BL(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
               this%fs%BG(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
            end do
         end do; end do
      end if
      ! Apply clipped Neumann on primitive variables in z+
      if (.not.this%cfg%zper.and.this%cfg%kproc.eq.this%cfg%npz) then
         do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
            ! Copy over from kmax to kmax+1 and above
            do k=this%cfg%kmax+1,this%cfg%kmaxo
               ! Copy primitive variables
               this%fs%RHOL(i,j,k)=this%fs%RHOL(i,j,this%cfg%kmax)
               this%fs%PL  (i,j,k)=this%fs%PL  (i,j,this%cfg%kmax)
               this%fs%IL  (i,j,k)=this%fs%IL  (i,j,this%cfg%kmax)
               this%fs%RHOG(i,j,k)=this%fs%RHOG(i,j,this%cfg%kmax)
               this%fs%PG  (i,j,k)=this%fs%PG  (i,j,this%cfg%kmax)
               this%fs%IG  (i,j,k)=this%fs%IG  (i,j,this%cfg%kmax)
               this%fs%U   (i,j,k)=this%fs%U   (i,j,this%cfg%kmax)
               this%fs%V   (i,j,k)=this%fs%V   (i,j,this%cfg%kmax)
               this%fs%W  (i,j,k)=max(this%fs%W(i,j,this%cfg%kmax),0.0_WP)
               this%fs%VF  (i,j,k)=this%fs%VF  (i,j,this%cfg%kmax)
               ! Also adjust interface data
               call setPlane(this%fs%PLIC(i,j,k),0,[0.0_WP,0.0_WP,+1.0_WP],this%cfg%z(k)+this%fs%dz*this%fs%VF(i,j,k))
               this%fs%BL(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
               this%fs%BG(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
            end do
         end do; end do
      end if
      ! Apply clipped Neumann on primitive variables in z-
      if (.not.this%cfg%zper.and.this%cfg%kproc.eq.1) then
         do j=this%cfg%jmino_,this%cfg%jmaxo_; do i=this%cfg%imino_,this%cfg%imaxo_
            ! First copy over W from kmin+1 to kmin
            this%fs%W(i,j,this%cfg%kmin)=min(this%fs%W(i,j,this%cfg%kmin+1),0.0_WP)
            ! Then copy over from kmin to kmin-1 and below
            do k=this%cfg%kmino,this%cfg%kmin-1
               ! Copy primitive variables
               this%fs%RHOL(i,j,k)=this%fs%RHOL(i,j,this%cfg%kmin)
               this%fs%PL  (i,j,k)=this%fs%PL  (i,j,this%cfg%kmin)
               this%fs%IL  (i,j,k)=this%fs%IL  (i,j,this%cfg%kmin)
               this%fs%RHOG(i,j,k)=this%fs%RHOG(i,j,this%cfg%kmin)
               this%fs%PG  (i,j,k)=this%fs%PG  (i,j,this%cfg%kmin)
               this%fs%IG  (i,j,k)=this%fs%IG  (i,j,this%cfg%kmin)
               this%fs%U   (i,j,k)=this%fs%U   (i,j,this%cfg%kmin)
               this%fs%V   (i,j,k)=this%fs%V   (i,j,this%cfg%kmin)
               this%fs%W  (i,j,k)=min(this%fs%W(i,j,this%cfg%kmin),0.0_WP)
               this%fs%VF  (i,j,k)=this%fs%VF  (i,j,this%cfg%kmin)
               ! Also adjust interface data
               call setPlane(this%fs%PLIC(i,j,k),0,[0.0_WP,0.0_WP,-1.0_WP],this%cfg%z(k)+this%fs%dz*this%fs%VF(i,j,k))
               this%fs%BL(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
               this%fs%BG(:,i,j,k)=[this%cfg%xm(i),this%cfg%ym(j),this%cfg%zm(k)]
            end do
         end do; end do
      end if
      ! Rebuild conserved quantities
      this%fs%Q(:,:,:,1)=        this%fs%VF *this%fs%RHOL
      this%fs%Q(:,:,:,2)=(1.0_WP-this%fs%VF)*this%fs%RHOG
      this%fs%Q(:,:,:,3)= this%fs%Q(:,:,:,1)*this%fs%IL
      this%fs%Q(:,:,:,4)= this%fs%Q(:,:,:,2)*this%fs%IG
      call this%fs%get_momentum()
    end block neumann_outflow

    ! Interpolate velocity
    call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
         
    ! Compute local Mach number
    this%Ma=sqrt(this%Ui**2+this%Vi**2+this%Wi**2)/this%fs%C

    ! Output to ensight
    if (this%ens_evt%occurs()) then
       !call this%fs%update_surfmesh(this%fs%smesh)
       call this%ens_out%write_data(this%time%t)
    end if

    call this%fs%get_info()
    call this%mfile%write()
    call this%cflfile%write()
    call this%consfile%write()
    !end do
  end subroutine step
  
  !> Finalize shock generator (sgen) simulation and transfer shock profile
  subroutine final(this,shockgen_group)
    use param,    only: param_read
    use parallel, only: amRoot
    use mpi_f08,  only: MPI_Group,MPI_ALLREDUCE,MPI_SUM,MPI_DOUBLE_PRECISION
    use parallel
    implicit none
    class (sgen), intent(inout) :: this
    type(MPI_Group), intent(in) :: shockgen_group
    integer  :: i,ierr
    real(WP) :: tol ! tolerance for finding shock

    ! Allocate and initialize centerline arrays (these are local to each proc)
    allocate(this%RHOG_center(this%cfg%imin:this%cfg%imax)); this%RHOG_center = 0.0_WP ! gas density
    allocate(  this%PG_center(this%cfg%imin:this%cfg%imax)); this%PG_center   = 0.0_WP ! gas pressure 
    allocate(  this%IG_center(this%cfg%imin:this%cfg%imax)); this%IG_center   = 0.0_WP ! gas internal energy
    allocate(  this%Ui_center(this%cfg%imin:this%cfg%imax)); this%Ui_center   = 0.0_WP ! velocity
    ! Allocate and initialize global arrays (this is where each proc stores their section of the centerline)
    allocate(this%RHOG_global(this%cfg%imin:this%cfg%imax)); this%RHOG_global = 0.0_WP ! gas density
    allocate(  this%PG_global(this%cfg%imin:this%cfg%imax)); this%PG_global   = 0.0_WP ! gas pressure 
    allocate(  this%IG_global(this%cfg%imin:this%cfg%imax)); this%IG_global   = 0.0_WP ! gas internal energy
    allocate(  this%Ui_global(this%cfg%imin:this%cfg%imax)); this%Ui_global   = 0.0_WP ! velocity

    ! read in nshock and allocate profile arrays (these are to be transfered to shockdrop)
    call param_read('nshock',this%nshock)        ! number of points left and right of shock center for thickening
    allocate(this%RHOG_profile(2*this%nshock+1)) ! density
    allocate( this%PG_profile(2*this%nshock+1))  ! pressure
    allocate( this%IG_profile(2*this%nshock+1))  ! internal energy
    allocate( this%Ui_profile(2*this%nshock+1))  ! velocity 
    
    ! if our proc is on the centerline, store each procs centerline variables
    ! AS: we're running in 1D now so I don't think we need to worry about the if statement for finding the centerline
    this%RHOG_center(this%cfg%imin_:this%cfg%imax_) = this%fs%RHOG(this%cfg%imin_:this%cfg%imax_,1,1)
    this%PG_center  (this%cfg%imin_:this%cfg%imax_) = this%fs%PG(  this%cfg%imin_:this%cfg%imax_,1,1)
    this%IG_center  (this%cfg%imin_:this%cfg%imax_) = this%fs%IG(  this%cfg%imin_:this%cfg%imax_,1,1)
    this%Ui_center  (this%cfg%imin_:this%cfg%imax_) = this%   Ui(  this%cfg%imin_:this%cfg%imax_,1,1) ! does Ui belong to sgen or mpcomp?

    ! use mpi_allreduce sum to 'concatenate' arrays (this works as a sum since the arrays have a value of 0 outside of their local proc domain)
    call MPI_ALLREDUCE(this%RHOG_center,this%RHOG_global,this%cfg%imax,MPI_DOUBLE_PRECISION,MPI_SUM,this%cfg%comm,ierr)
    call MPI_ALLREDUCE(this%PG_center  ,this%PG_global  ,this%cfg%imax,MPI_DOUBLE_PRECISION,MPI_SUM,this%cfg%comm,ierr)
    call MPI_ALLREDUCE(this%IG_center  ,this%IG_global  ,this%cfg%imax,MPI_DOUBLE_PRECISION,MPI_SUM,this%cfg%comm,ierr)
    call MPI_ALLREDUCE(this%Ui_center  ,this%Ui_global  ,this%cfg%imax,MPI_DOUBLE_PRECISION,MPI_SUM,this%cfg%comm,ierr)

    ! loop over domain and find the shock
    do i=this%cfg%imin,this%cfg%imax
       !print*, "sgen: final: index =  ", i
       if ((this%cfg%xm(i).lt.(this%Xend+tol)).and.(this%cfg%xm(i).gt.(this%Xend-tol)))then
          !print*, "sgen: found the shock index at i = ", i
          this%shock_index = i
       end if
    end do

    ! save the shock profile
    do i=this%shock_index-this%nshock,this%shock_index+this%nshock
       this%RHOG_profile(i-this%shock_index+1) = this%RHOG_global(i)
       this%PG_profile  (i-this%shock_index+1) = this%PG_global  (i)
       this%IG_profile  (i-this%shock_index+1) = this%IG_global  (i)
       this%Ui_profile  (i-this%shock_index+1) = this%Ui_global  (i)
    end do
  end subroutine final

end module shockgen_class
