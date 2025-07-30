!> Definition of shock generator class
! AS intended to run in serial for now 
module shockgen_class  
  use precision,         only: WP
  use string,            only: str_medium
  use inputfile_class,   only: inputfile
  use config_class,      only: config
  use spcomp_class,      only: spcomp
  use timetracker_class, only: timetracker
  use ensight_class,     only: ensight
  !use event_class,       only: event
  !use monitor_class,     only: monitor
  implicit none
  private
  public :: sgen
  !> Equations of state
  real(WP) :: PinfL,GammaL,CvL
  real(WP) :: PinfG,GammaG,CvG

  !>sgen object
  type :: sgen
     !> config
     type(config) :: cfg
     !> Flow solver
     type(spcomp) :: fs        !< Single-phase compressible solver
     type(timetracker) :: time !< Time info
     !> Ensight postprocessing
     !type(ensight)  :: ens_out
     !> Simulation monitor file
     !type(monitor) :: mfile,cflfile,consfile
     !> Work arrays
     real(WP), dimension(:,:,:,:,:), allocatable :: dQdt
     real(WP), dimension(:,:,:)    , allocatable :: Ui,Vi,Wi,Ma,beta,visc
     !> Constant phasic kinematic viscosity
     real(WP) :: cst_visc
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
     integer  :: shock_index,nshock                                             ! i index of shock and number of points to thicken shock
     real(WP), dimension(:), allocatable :: RHOG_profile,IG_profile,PG_profile,Ui_profile ! shock profile arrays
     !> domain partiion
     integer, dimension(3) :: partition
     !> Flags
     logical :: dim_flag ! Flag for dimensional or nondimensional init
     logical :: inviscid_flag  ! Flag for running inviscid
   contains
     procedure :: init                            !< Initialize shock-drop simulation
     procedure :: step                            !< Advance shock-drop simulation by one time step
     procedure :: finalize                        !< Finalize shock-drop simulation
     !procedure :: output_monitor                  !< Monitoring for shock-drop case
     !procedure :: output_ensight                  !< Ensight output for shock-drop case
     !procedure, private :: prepare_viscosities    !< Prepare viscosities
     procedure, private :: apply_bconds           !< Apply boundary conditions
  end type sgen

contains

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

  !> initialization of shock generator simulation
  subroutine init(this,sgen_group)
    use sgrid_class, only: cartesian,sgrid
    use param,       only: param_read
    use mpi_f08,     only: MPI_Group
    implicit none
    class(sgen), intent(inout)  :: this
    type(MPI_Group), intent(in) :: sgen_group
    type(sgrid) :: grid
    integer  :: i,nx,ny,nz
    real(WP) :: Lx 
    real(WP) :: CPD,D0X
    real(WP) :: dx,dy,dz
    logical  :: xper,yper,zper
    real(WP), dimension(:), allocatable :: x,y,z

    call param_read('Dimensional flag', this%dim_flag)
    call param_read('Inviscid flag',this%inviscid_flag)
    call param_read('Drop diameter', this%ddrop) ! usd for mesh generation
    call param_read('CPD', CPD)
    call param_read('D0X',D0X)

    create_grid: block
      !> generate 1D mesh
      ! set periodic BCs
      ! AS note: when we set periodic in x, we see a strange interpolated velocity value at the right boundary
      xper=.false.; yper=.true.; zper=.true.
      ! compute physical domain length
      Lx = D0X*this%ddrop                   ! compute domain lengths
      nx = ceiling((CPD*Lx)/this%ddrop)     ! compute number of uniform cells
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
      !General serial grid object
      grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=xper,yper=yper,zper=zper,name='ShockGen')
    end block create_grid
    create_cfg: block
      use parallel, only: group,comm
      ! Read in partition
      this%partition = (/1,1,1/) ! run in serial
      ! Create partitioned grid
      this%cfg=config(grp=sgen_group,decomp=this%partition,grid=grid) 
    end block create_cfg

    !> create timetracker
    create_timetracker:block
      use parallel, only : amRoot
      use param,    only : param_read
      ! initialize timetracker
      this%time=timetracker(amRoot=this%cfg%amRoot,name='shockgen_timer')
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
      call param_read('Shockgen start',this%Xstart) ! drop diameter lengths
      call param_read('Shockgen end',this%Xend)     ! drop diameter 
      if (this%dim_flag.eqv.(.false.))then ! non dimensional
         ! First generate static shock with normalized pre-shock conditions
         this%M1=this%Ms
         this%rho1=1.0_WP
         this%rho2=this%rho1*(GammaG+1.0_WP)*this%M1**2/((GammaG-1.0_WP)*this%M1**2+2.0_WP)
         this%p1=0.25_WP*this%rho1/GammaG*((GammaG+1.0_WP)*this%M1/(this%M1**2-1.0_WP))**2 ! Ensures that |u2-u1|=1
         this%p2=this%p1*(2.0_WP*GammaG/(GammaG+1.0_WP)*(this%M1**2-1.0_WP)+1.0_WP)
         this%u1=this%M1*sqrt(GammaG*this%p1/this%rho1); ! shock speed through unshocked air
         ! calculate final shock generator time based on how far the shock moves
         this%time%tmax = (this%Xend - this%Xstart)/this%u1 ! assumes ddrop = 1
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
         this%time%tmax = (this%Xend*this%ddrop - this%Xstart*this%ddrop)/this%u1
      end if
    end block create_timetracker

    !> Create singelphase compressible flow solver
    create_velocity_solver: block
      call this%fs%initialize(cfg=this%cfg,name='Compressible NS')
      this%fs%getP=>get_PG; this%fs%getC=>get_CG; this%fs%getS=>get_SG; this%fs%getT=>get_TG
    end block create_velocity_solver

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

        ! Viscous parameters
        if (this%inviscid_flag.eqv.(.true.))then
           this%viscG = 0.0_WP
        else
           call param_read('Gas dynamic viscosity',this%viscG)        
        end if

        ! Use shock relations (shock fixed frame) to calculate post-shock conditions 
        this%M1 = this%Ms ! set M1 equal to shock Mach number for now 
        this%rho2=this%rho1*(GammaG+1.0_WP)*this%M1**2/((GammaG-1.0_WP)*this%M1**2+2.0_WP) ! post shock density  (Anderson 3.53)
        this%p2=this%p1*(2.0_WP*GammaG/(GammaG+1.0_WP)*(this%M1**2-1.0_WP)+1.0_WP)         ! post shock pressure (Anderson 3.57)
        this%u1=this%M1*sqrt(GammaG*this%p1/this%rho1)                                     ! velocity in state 1 (left side of shock in fixed frame)
        this%u2=this%u1*this%rho1/this%rho2                                ! velocity in state 2 (Anderson 3.53, right side of shock in fixed frame)
        ! we now shift to accomodate a moving shock (converting from shock fixed frame to lab frame)
        this%u2=abs(this%u2-this%u1); this%M2=this%u2/sqrt(GammaG*this%p2/this%rho2) ! post-shock gas velocity and post shock Mach number
        this%u1=0.0_WP; this%M1=this%u1/sqrt(GammaG*this%p1/this%rho1)         ! set pre-shock gas velocity to zero and update pre-shock Mach number

        ! compute some non-dimensional parameters for log files
        this%rho_ratio=this%rho2/this%rho1       ! density ratio
        this%visc_ratio=this%viscL/this%viscG    ! viscosity ratio
        this%ReG = this%rho2*this%u2*this%ddrop/this%viscG ! Reynolds number based on post-shock conditions 
        this%c_ratio=sqrt(GammaL*(this%p1+PinfL)/this%rhoL)/sqrt(GammaG*this%p1/this%rho1) ! sound speed ratio
        this%ML=this%u2/sqrt(GammaL*(this%p1+PinfL)/this%rhoL)
     end if
    end block initialize_parameters

    ! Allocate work arrays
    allocate_work_arrays: block
      allocate(this%dQdt(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:this%fs%nQ,1:4))
      allocate(this%beta(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(this%visc(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(this%Ui  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(this%Vi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(this%Wi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(this%Ma  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
    end block allocate_work_arrays

    ! Prepare initial conditions
    initial_conditions: block
      integer :: i,j,k
      ! Initialize primary variables
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if ((this%cfg%xm(i).lt.this%Xstart*this%ddrop)) then    ! post shock region
                  this%fs%U(i,j,k)   = this%u2; this%fs%V(i,j,k) = 0.0_WP; this%fs%W(i,j,k) = 0.0_WP ! velocity
                  this%fs%Q(i,j,k,1) = this%rho2 ! density
                  this%fs%P(i,j,k)   = this%p2   ! pressure
                  this%fs%I(i,j,k)   = (this%fs%P(i,j,k)+GammaG*PinfG)/(this%fs%Q(i,j,k,1)*(GammaG-1.0_WP)) ! internal energy
               else
                  this%fs%U(i,j,k)   = this%u1; this%fs%V(i,j,k) = 0.0_WP; this%fs%W(i,j,k) = 0.0_WP ! velocity
                  this%fs%Q(i,j,k,1) = this%rho1 ! density
                  this%fs%P(i,j,k)   = this%p1   ! pressure
                  this%fs%I(i,j,k)   = (this%fs%P(i,j,k)+GammaG*PinfG)/(this%fs%Q(i,j,k,1)*(GammaG-1.0_WP)) ! internal energy
               end if
            end do
         end do
      end do

      ! Initialize conserved variables
      this%fs%Q(:,:,:,2)=this%fs%Q(:,:,:,1)*this%fs%I
      call this%fs%get_momentum()
      ! Rebuild primitive variables
      call this%fs%get_primitive()
      ! Interpolate velocity
      call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
      ! Compute local Mach number
      this%Ma=sqrt(this%Ui**2+this%Vi**2+this%Wi**2)/this%fs%C
    end block initial_conditions

    ! ! Add Ensight output
    ! create_ensight: block
    !   ! Create Ensight output from cfg
    !   this%ens_out=ensight(cfg=this%cfg,name='Shockgen')
    !   ! Add variables to output
    !   call this%ens_out%add_vector('velocity',this%Ui,this%Vi,this%Wi)
    !   call this%ens_out%add_scalar('RHO',this%fs%Q(:,:,:,1))
    !   call this%ens_out%add_scalar('I',this%fs%I)
    !   call this%ens_out%add_scalar('P',this%fs%P)
    !   call this%ens_out%add_scalar('Mach',this%Ma)
    !   call this%ens_out%add_scalar('beta',this%beta)
    !   call this%ens_out%add_scalar('visc',this%visc)
    !   call this%ens_out%add_vector('cell_vel',this%fs%U,this%fs%V,this%fs%W)
    !   call this%output_ensight()
    ! end block create_ensight

    ! ! Create monitor files
    ! create_monitor: block
    !   ! Create simulation monitor
    !   this%mfile=monitor(this%fs%cfg%amRoot,'shockgen_sim')
    !   call this%mfile%add_column(this%time%n,'Timestep number')
    !   call this%mfile%add_column(this%time%t,'Time')
    !   call this%mfile%add_column(this%time%dt,'Timestep size')
    !   call this%mfile%add_column(this%time%cfl,'Maximum CFL')
    !   call this%mfile%add_column(this%fs%Umax,'Umax')
    !   call this%mfile%add_column(this%fs%Vmax,'Vmax')
    !   call this%mfile%add_column(this%fs%Wmax,'Wmax')
    !   call this%mfile%add_column(this%fs%RHOmax,'max(RHO)')
    !   call this%mfile%add_column(this%fs%RHOmin,'min(RHO)')
    !   call this%mfile%add_column(this%fs%Imax  ,'max(I)'  )
    !   call this%mfile%add_column(this%fs%Imin  ,'min(I)'  )
    !   call this%mfile%add_column(this%fs%Pmax  ,'max(P)'  )
    !   call this%mfile%add_column(this%fs%Pmin  ,'min(P)'  )
    !   call this%mfile%add_column(this%fs%Tmax  ,'max(T)'  )
    !   call this%mfile%add_column(this%fs%Tmin  ,'min(T)'  )
    !   ! Create CFL monitor
    !   this%cflfile=monitor(this%fs%cfg%amRoot,'shockgen_cfl')
    !   call this%cflfile%add_column(this%time%n,'Timestep number')
    !   call this%cflfile%add_column(this%time%t,'Time')
    !   call this%cflfile%add_column(this%fs%CFLc_x,'Convective xCFL')
    !   call this%cflfile%add_column(this%fs%CFLc_y,'Convective yCFL')
    !   call this%cflfile%add_column(this%fs%CFLc_z,'Convective zCFL')
    !   call this%cflfile%add_column(this%fs%CFLa_x,'Acoustic xCFL')
    !   call this%cflfile%add_column(this%fs%CFLa_y,'Acoustic yCFL')
    !   call this%cflfile%add_column(this%fs%CFLa_z,'Acoustic zCFL')
    !   call this%cflfile%add_column(this%fs%CFLv_x,'Viscous xCFL')
    !   call this%cflfile%add_column(this%fs%CFLv_y,'Viscous yCFL')
    !   call this%cflfile%add_column(this%fs%CFLv_z,'Viscous zCFL')
    !   ! Create conservation monitor
    !   this%consfile=monitor(this%fs%cfg%amRoot,'shockgen_cons')
    !   call this%consfile%add_column(this%time%n,'Timestep number')
    !   call this%consfile%add_column(this%time%t,'Time')
    !   call this%consfile%add_column(this%fs%Qint(1),'Mass')
    !   call this%consfile%add_column(this%fs%Qint(2),'Energy')
    !   call this%consfile%add_column(this%fs%Qint(3),'U Momentum')
    !   call this%consfile%add_column(this%fs%Qint(4),'V Momentum')
    !   call this%consfile%add_column(this%fs%Qint(5),'W Momentum')
    !   call this%consfile%add_column(this%fs%RHOKint,'Kinetic Energy')
    !   call this%consfile%add_column(this%fs%RHOSint,'Entropy')
    ! end block create_monitor
    ! call this%output_monitor()
  end subroutine init

  ! we need to create the step subroutine, then we need to add final subroutine and then we can start adding these into simulaiton.f90

  subroutine step(this)
    implicit none
    class(sgen), intent(inout)  :: this
    
    ! Increment time
    call this%fs%get_cfl(dt=this%time%dt,cfl=this%time%cfl)
    call this%time%increment()

    ! Remember conserved variables
    this%fs%Qold=this%fs%Q

    ! Prepare SGS viscosity models
    !call this%prepare_viscosities()
    ! Get LAD
    call this%fs%get_viscartif(dt=this%time%dt,beta=this%beta); this%fs%BETA=this%fs%Q(:,:,:,1)*(this%beta              )
    if(this%inviscid_flag.eqv.(.false.))then ! only get vreman model if we're running viscous
       ! Get eddy viscosity
       call this%fs%get_vreman   (dt=this%time%dt,visc=this%visc); this%fs%VISC=this%fs%Q(:,:,:,1)*(this%visc+this%cst_visc)
    end if

    ! --> we're getting the error from the trhs timer already being started, I'm not sure how to fix it besides commenting out
    ! First RK step ====================================================================================
    ! Get non-SL RHS and increment
    call this%fs%rhs(this%dQdt(:,:,:,:,1))
    this%fs%Q=this%fs%Qold+0.5_WP*this%time%dt*this%dQdt(:,:,:,:,1)
    ! Recompute primitive variables
    call this%fs%get_primitive()
    ! Second RK step ===================================================================================
    ! Get non-SL RHS and increment
    call this%fs%rhs(this%dQdt(:,:,:,:,2))
    this%fs%Q=this%fs%Qold+0.5_WP*this%time%dt*this%dQdt(:,:,:,:,2)
    ! Recompute primitive variables
    call this%fs%get_primitive()
    ! Third RK step ====================================================================================
    ! Get non-SL RHS and increment
    call this%fs%rhs(this%dQdt(:,:,:,:,3))
    this%fs%Q=this%fs%Qold+1.0_WP*this%time%dt*this%dQdt(:,:,:,:,3)
    ! Recompute primitive variables
    call this%fs%get_primitive()
    ! Fourth RK step ===================================================================================
    ! Get non-SL RHS and increment
    call this%fs%rhs(this%dQdt(:,:,:,:,4))
    this%fs%Q=this%fs%Qold+this%time%dt/6.0_WP*(this%dQdt(:,:,:,:,1)+2.0_WP*this%dQdt(:,:,:,:,2)+2.0_WP*this%dQdt(:,:,:,:,3)+this%dQdt(:,:,:,:,4))
    ! Recompute primitive variables
    call this%fs%get_primitive()

    ! Apply boundary conditions
    call this%apply_bconds()

    ! Interpolate velocity
    call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)

    ! Compute local Mach number
    this%Ma=sqrt(this%Ui**2+this%Vi**2+this%Wi**2)/this%fs%C
    !call this%output_monitor()
    !call this%output_ensight()
  end subroutine step

  !> Finalize shockgenerator simulation (get shock profile for simulation.f90)
  subroutine finalize(this)
    use param,    only : param_read
    use mpi_f08,  only : MPI_DOUBLE_PRECISION,MPI_BCAST,MPI_COMM_WORLD
    implicit none
    class(sgen), intent(inout) :: this
    integer  :: i,ierr,initial_index
    real(WP) :: min_distance

    call param_read('nshock',this%nshock)
    ! allocate shock profile arrays
    ! these need to be allocated on every proc --> allocated in sim.f90
    allocate(this%RHOG_profile(2*this%nshock+1),this%PG_profile(2*this%nshock+1),this%IG_profile(2*this%nshock+1),this%Ui_profile(2*this%nshock+1))
    this%RHOG_profile = 0.0_WP; this%PG_profile = 0.0_WP; this%IG_profile = 0.0_WP; this%Ui_profile = 0.0_WP
    ! we're setup to run this simulaiton in serial for now, so no need to worry about proc subdomains (unless this becomes too expensive)
    this%shock_index = -10; initial_index = -10; min_distance=HUGE(1.0_WP) ! initialize shock_index and min_distance
    do i = this%cfg%imin, this%cfg%imax
       if (abs(this%cfg%xm(i) - this%Xend*this%ddrop).lt.min_distance) then
          min_distance = abs(this%cfg%xm(i) - this%Xend*this%ddrop)
          this%shock_index = i          
       end if
    end do
    
    ! now that we've found the shock, lets store the profile variables
    do i=this%shock_index-this%nshock,this%shock_index+this%nshock              ! saving 2*n_shock+1 points (accounts for center point)
       this%RHOG_profile(i-this%shock_index+this%nshock+1) = this%fs%Q(i,1,1,1) ! density 
       this%IG_profile  (i-this%shock_index+this%nshock+1) = this%fs%I(i,1,1)   ! internal energy
       this%PG_profile  (i-this%shock_index+this%nshock+1) = this%fs%P(i,1,1)   ! pressure
       this%Ui_profile  (i-this%shock_index+this%nshock+1) = this%Ui  (i,1,1)   ! velocity
    end do
    if (this%shock_index.eq.-10)then
       print*, "WARNING: shock generator did not find the shock."
    else
       print*, "sgen: stored shock index", this%shock_index
       print*, "RHOG profile: ", this%RHOG_profile
       print*, "Internal energy profile: ", this%IG_profile
       print*, "Pressure profile: ", this%PG_profile
       print*, "Velocity profile: ", this%Ui_profile
    end if

    ! communicate profile arrays to all other processors (0 --> root proc) ! this was moved to simulation.f90
    ! call MPI_BCAST(this%RHOG_profile,this%nshock+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    ! call MPI_BCAST(this%IG_profile  ,this%nshock+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    ! call MPI_BCAST(this%PG_profile  ,this%nshock+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    ! call MPI_BCAST(this%Ui_profile  ,this%nshock+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    
    ! Deallocate work arrays
    deallocate(this%dQdt,this%Ui,this%Vi,this%Wi,this%Ma,this%beta,this%visc)
  end subroutine finalize

  !> Calculate viscosities
  !subroutine prepare_viscosities(this)
  !  implicit none
  !  class(sgen), intent(inout) :: this
  !  ! Get LAD
  !  call this%fs%get_viscartif(dt=this%time%dt,beta=this%beta); this%fs%BETA=this%fs%Q(:,:,:,1)*(this%beta              )
  !  ! Get eddy viscosity
  !  call this%fs%get_vreman   (dt=this%time%dt,visc=this%visc); this%fs%VISC=this%fs%Q(:,:,:,1)*(this%visc+this%cst_visc)
  !end subroutine prepare_viscosities

  !> Apply boundary conditions
  subroutine apply_bconds(this)
    implicit none
    class(sgen), intent(inout) :: this
    integer :: i,j,k

    ! Apply clipped Neumann on primitive variables in x+
    if (.not.this%fs%cfg%xper.and.this%fs%cfg%iproc.eq.this%fs%cfg%npx) then
       do k=this%fs%cfg%kmino_,this%fs%cfg%kmaxo_; do j=this%fs%cfg%jmino_,this%fs%cfg%jmaxo_
          ! Copy over from imax to imax+1 and above
          do i=this%fs%cfg%imax+1,this%fs%cfg%imaxo
             ! Copy primitive variables
             this%fs%Q(i,j,k,1)=this%fs%Q(this%fs%cfg%imax,j,k,1)
             this%fs%P(i,j,k)=this%fs%P(this%fs%cfg%imax,j,k)
             this%fs%I(i,j,k)=this%fs%I(this%fs%cfg%imax,j,k)
             this%fs%U(i,j,k)=max(this%fs%U(this%fs%cfg%imax,j,k),0.0_WP)
             this%fs%V(i,j,k)=this%fs%V(this%fs%cfg%imax,j,k)
             this%fs%W(i,j,k)=this%fs%W(this%fs%cfg%imax,j,k)
          end do
       end do; end do
    end if

    ! Apply clipped Neumann on primitive variables in y+
    if (.not.this%fs%cfg%yper.and.this%fs%cfg%jproc.eq.this%fs%cfg%npy) then
       do k=this%fs%cfg%kmino_,this%fs%cfg%kmaxo_; do i=this%fs%cfg%imino_,this%fs%cfg%imaxo_
          ! Copy over from jmax to jmax+1 and above
          do j=this%fs%cfg%jmax+1,this%fs%cfg%jmaxo
             ! Copy primitive variables
             this%fs%Q(i,j,k,1)=this%fs%Q(i,this%fs%cfg%jmax,k,1)
             this%fs%P(i,j,k)=this%fs%P(i,this%fs%cfg%jmax,k)
             this%fs%I(i,j,k)=this%fs%I(i,this%fs%cfg%jmax,k)
             this%fs%U(i,j,k)=this%fs%U(i,this%fs%cfg%jmax,k)
             this%fs%V(i,j,k)=max(this%fs%V(i,this%fs%cfg%jmax,k),0.0_WP)
             this%fs%W(i,j,k)=this%fs%W(i,this%fs%cfg%jmax,k)
          end do
       end do; end do
    end if

    ! Apply clipped Neumann on primitive variables in y-
    if (.not.this%fs%cfg%yper.and.this%fs%cfg%jproc.eq.1) then
       do k=this%fs%cfg%kmino_,this%fs%cfg%kmaxo_; do i=this%fs%cfg%imino_,this%fs%cfg%imaxo_
          ! First copy over V from jmin+1 to jmin
          this%fs%V(i,this%fs%cfg%jmin,k)=min(this%fs%V(i,this%fs%cfg%jmin+1,k),0.0_WP)
          ! Then copy over from jmin to jmin-1 and below
          do j=this%fs%cfg%jmino,this%fs%cfg%jmin-1
             ! Copy primitive variables
             this%fs%Q(i,j,k,1)=this%fs%Q(i,this%fs%cfg%jmin,k,1)
             this%fs%P(i,j,k)=this%fs%P(i,this%fs%cfg%jmin,k)
             this%fs%I(i,j,k)=this%fs%I(i,this%fs%cfg%jmin,k)
             this%fs%U(i,j,k)=this%fs%U(i,this%fs%cfg%jmin,k)
             this%fs%V(i,j,k)=min(this%fs%V(i,this%fs%cfg%jmin,k),0.0_WP)
             this%fs%W(i,j,k)=this%fs%W(i,this%fs%cfg%jmin,k)
          end do
       end do; end do
    end if

    ! Apply clipped Neumann on primitive variables in z+
    if (.not.this%fs%cfg%zper.and.this%fs%cfg%kproc.eq.this%fs%cfg%npz) then
       do j=this%fs%cfg%jmino_,this%fs%cfg%jmaxo_; do i=this%fs%cfg%imino_,this%fs%cfg%imaxo_
          ! Copy over from kmax to kmax+1 and above
          do k=this%fs%cfg%kmax+1,this%fs%cfg%kmaxo
             ! Copy primitive variables
             this%fs%Q(i,j,k,1)=this%fs%Q(i,j,this%fs%cfg%kmax,1)
             this%fs%P(i,j,k)=this%fs%P(i,j,this%fs%cfg%kmax)
             this%fs%I(i,j,k)=this%fs%I(i,j,this%fs%cfg%kmax)
             this%fs%U(i,j,k)=this%fs%U(i,j,this%fs%cfg%kmax)
             this%fs%V(i,j,k)=this%fs%V(i,j,this%fs%cfg%kmax)
             this%fs%W(i,j,k)=max(this%fs%W(i,j,this%fs%cfg%kmax),0.0_WP)
          end do
       end do; end do
    end if

    ! Apply clipped Neumann on primitive variables in z-
    if (.not.this%fs%cfg%zper.and.this%fs%cfg%kproc.eq.1) then
       do j=this%fs%cfg%jmino_,this%fs%cfg%jmaxo_; do i=this%fs%cfg%imino_,this%fs%cfg%imaxo_
          ! First copy over W from kmin+1 to kmin
          this%fs%W(i,j,this%fs%cfg%kmin)=min(this%fs%W(i,j,this%fs%cfg%kmin+1),0.0_WP)
          ! Then copy over from kmin to kmin-1 and below
          do k=this%fs%cfg%kmino,this%fs%cfg%kmin-1
             ! Copy primitive variables
             this%fs%Q(i,j,k,1)=this%fs%Q(i,j,this%fs%cfg%kmin,1)
             this%fs%P(i,j,k)=this%fs%P(i,j,this%fs%cfg%kmin)
             this%fs%I(i,j,k)=this%fs%I(i,j,this%fs%cfg%kmin)
             this%fs%U(i,j,k)=this%fs%U(i,j,this%fs%cfg%kmin)
             this%fs%V(i,j,k)=this%fs%V(i,j,this%fs%cfg%kmin)
             this%fs%W(i,j,k)=min(this%fs%W(i,j,this%fs%cfg%kmin),0.0_WP)
          end do
       end do; end do
    end if

    ! Rebuild conserved quantities
    this%fs%Q(:,:,:,2)=this%fs%Q(:,:,:,1)*this%fs%I
    call this%fs%get_momentum()

  end subroutine apply_bconds

  ! !> Perform and output monitoring for the far-field shock problem
  ! subroutine output_monitor(this)
  !   implicit none
  !   class(sgen), intent(inout) :: this
  !   call this%fs%get_info()
  !   call this%mfile%write()
  !   call this%cflfile%write()
  !   call this%consfile%write()
  ! end subroutine output_monitor

  ! !> Output ensight files for the far-field shock problem
  ! subroutine output_ensight(this,t)
  !   implicit none
  !   class(sgen), intent(inout) :: this
  !   real(WP), intent(in), optional :: t
  !   if (present(t)) then
  !      call this%ens_out%write_data(t)
  !   else
  !      call this%ens_out%write_data(this%time%t)
  !   end if
  ! end subroutine output_ensight
  
end module shockgen_class
