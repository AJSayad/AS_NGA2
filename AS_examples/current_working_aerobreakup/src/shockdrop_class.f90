!> Definition for a shockdrop_class (shock generator)
module shockdrop_class
  use precision,         only: WP
  use config_class,      only: config
  use mast_class,        only: mast
  use vfs_class,         only: vfs
  use matm_class,        only: matm
  use timetracker_class, only: timetracker
  use ensight_class,     only: ensight
  use surfmesh_class,    only: surfmesh
  use event_class,       only: event
  use monitor_class,     only: monitor
  use hypre_str_class,   only: hypre_str
  use pardata_class,     only: pardata
  use param,             only: param_read
  use shockgen_class,    only: sgen
  implicit none
  private

  real(WP), public :: ddrop
  real(WP), public, dimension(3) :: dctr

  public :: sdrop

  !> sdrop object
  type :: sdrop
     !> Config
     type(config)       :: cfg !> Mesh for solver
     !> Flow solver
     type(mast)         :: fs
     type(vfs)          :: vf
     type(matm)         :: matmod
     type(timetracker)  :: time
     type(hypre_str)    :: ps
     type(hypre_str)    :: vs
     !> surface mesh
     type(surfmesh)     :: smesh
     !> Ensight postprocessing
     type(ensight)      :: ens_out, ens_out_smesh 
     type(event)        :: ens_evt, ens_evt_smesh
     !> Simulation monitor file
     type(monitor)      :: mfile,cflfile,cvgfile
     !> Fluid parameters
     integer            :: relax_model
     integer            :: shock_index,n_shock
     real(WP)           :: xshock
     !> Provide a pardata and an event tracker for saving restarts
     type(event)   :: save_evt
     type(pardata) :: df
     logical       :: restarted     
   contains
     procedure, public :: init_grid                !> initialize grid for shockdrop simulation
     procedure, public :: init                     !> initialize sdrop simulation
     procedure, public :: update_mixture_variables !> update mixture variables after shock profile
     procedure :: writeIC                          !> write initial conditions to ensight
     procedure :: step                             !> advance sgen simulation by one timestep
     procedure :: final                            !> finalize sgen simulation
     procedure :: restart                          !> restart simulation from timestamp
  end type sdrop

contains

  !> Function that defines a level set function for a cylindrical droplet (2D)
  function levelset_cyl(xyz,t) result(G)
    implicit none
    real(WP), dimension(3),intent(in) :: xyz
    real(WP), intent(in) :: t
    real(WP) :: G
    G=1.0_WP-sqrt((xyz(1)-dctr(1))**2+(xyz(2)-dctr(2))**2)/(ddrop/2.0)
  end function levelset_cyl
  
  !> Function that defines a level set function for a spherical droplet (3D)
  function levelset_sphere(xyz,t) result(G)
    implicit none
    real(WP), dimension(3),intent(in) :: xyz
    real(WP), intent(in) :: t
    real(WP) :: G
    G=1.0_WP-sqrt((xyz(1)-dctr(1))**2+(xyz(2)-dctr(2))**2+(xyz(3)-dctr(3))**2)/(ddrop/2.0)
  end function levelset_sphere
  
  subroutine init_grid(this)
    implicit none
    class(sdrop), intent(inout) :: this
    
    ! Create mesh for sdrop
    create_sdrop_config: block
      use sgrid_class, only: cartesian,sgrid
      use param,       only: param_read, param_exists
      use parallel,    only: amRoot,group
      use messager,    only: die
      type(sgrid) :: grid
      integer, dimension(3) :: partition
      integer :: i,j,k,nx,ny,nz
      real(WP) :: ddrop, dx
      real(WP) :: Lx,Ly,Lz
      real(WP), dimension(:), allocatable :: x,y,z
      
      ! variables for stretching in x
      integer ::  nx_stretchL,nx_stretchR
      real(WP) :: dx_old,alpha,dx_ref,start_ref
      
      ! variables for stretching in y
      real(WP) :: dy,dy_old,dy_stretch
      integer ::  ny_stretch
      
      ! variables for stretching in z
      real(WP) :: dz,dz_old,dz_stretch
      integer :: nz_stretch

      alpha=1.03_WP ! mesh stretching ratio
      ! Read in grid definition
      call param_read('Lx',Lx); call param_read('Lx ref', start_ref, default=0.0_WP);
      call param_read('nx',nx); call param_read('nx stretch left',nx_stretchL); call param_read('nx stretch right',nx_stretchR);
      dx = Lx/nx; dx_ref = (Lx - start_ref)/real(nx,WP)
      call param_read('Ly',Ly); call param_read('ny',ny); call param_read('ny stretch',ny_stretch);
      call param_read('nz',nz,default=1); call param_read('nz stretch',nz_stretch)
      if (nz.eq.1) then
         Lz = dx
      else
         call param_read('Lz',Lz)
      end if      
      allocate(x(nx+nx_stretchL+nx_stretchR+1));allocate(y(ny+2*ny_stretch+1));allocate(z(nz+2*nz_stretch+1));

      ! Read in droplet information
      call param_read('Droplet diameter',ddrop)
      
      !uniform mesh x
      do i=nx_stretchL+1,nx+nx_stretchL+1
         x(i) = start_ref + real(i-1-nx_stretchL,WP)*dx_ref
      end do
      
      ! stretch left of domain
      do i=nx_stretchL,1,-1
         dx_old = abs(x(i+2) - x(i+1))
         x(i) = x(i+1) - dx_old*alpha
      end do
      
      ! stretch right of domain
      do i=nx+nx_stretchL+2,nx+nx_stretchL+nx_stretchR+1
         dx_old = x(i-1)-x(i-2)
         x(i) = x(i-1)+dx_old*alpha
      end do      

      ! y mesh
      do j = 1,ny+2*ny_stretch+1 !initialize y mesh array
         y(j) = 0.0_WP
      end do
      
      dy = Ly/ny ! define uniform grid spacing
      y(ny/2+ny_stretch+1) = 0.0_WP !define the centerline of the domain
      
      !y array uniform region
      do j = ny/2+ny_stretch+2,ny+ny_stretch+1
         y(j) = y(j-1) + dy
      end do
      
      !stretching in y
      do j = ny+ny_stretch+2,ny+(2*ny_stretch)+1
         dy_old = y(j-2) - y(j-3)
         dy_stretch = alpha*dy_old
         y(j) = y(j-1) + dy_stretch
      end do
      
      ! mirror y across y=0 line
      do j = 1,ny/2+ny_stretch
         y(j) = -y(ny-j+(2*ny_stretch)+2)
      end do
      
      ! z mesh
      dz = Lz/nz ! define uniform grid spacing
      z(nz/2+nz_stretch+1) = 0.0_WP !define centerline
      
      do k = nz/2+nz_stretch+2,nz+nz_stretch+1
         z(k) = z(k-1) + dz
      end do
      
      !stretching in z
      do k = nz+nz_stretch+2,nz+(2*nz_stretch)+1
         dz_old = z(k-2) - z(k-3)
         dz_stretch = alpha*dz_old
         z(k) = z(k-1) + dz_stretch
      end do
      
      ! mirror z across centerline
      do k = 1,nz/2+nz_stretch
         z(k) = -z(nz-k+(2*nz_stretch)+2)
      end do
      
      if(amRoot)then
         print*, "======== MESH DESCRIPTION IN x ========"
         print*, "Uniform region length: ", Lx
         print*, "Number of cells in the uniform region: ", nx
         print*, "Number of cells added to left of domain: ", nx_stretchL
         print*, "Number of cells added to right of domain: ", nx_stretchR
         print*, "Total number of cells in domain: ", nx+nx_stretchL+nx_stretchR
         print*, "Stretching ratio in x: ", alpha
         print*, "Leftmost point (stretched region left): ", x(1)
         print*, "Start of uniform region: ", start_ref
         print*, "End of uniform region: ",x(nx+nx_stretchL+1)
         print*, "Rightmost point (stretched region right): ", x(nx_stretchL+nx+nx_stretchR+1)
         print*, "Aspect ratio in uniform region: ", dx/dy
         print*, "Number of cells per droplet diameter: ", nx*(ddrop/Lx)
         print*, "======================================="
         print*, "======== MESH DESCRIPTION IN y ========"
         print*, 'Uniform region height: ', Ly
         print*, 'Stretching in y starts at: +- ', Ly/2
         print*, 'Number of cells added to the top and to the bottom: ', ny_stretch/2
         print*, 'Stretching ratio in y: ', alpha
         print*, "Aspect ratio in uniform region: ", dx/dy
         print*, "Number of cells per droplet diameter: ", ny*(ddrop/Ly)
         print*, "========================================"
         print*, "======== MESH DESCRIPTION IN z ========"
         print*, 'Uniform region height: ', Lz
         print*, 'Stretching in z starts at: +- ', Lz/2
         print*, 'Number of cells added to the front and to the back: ', nz_stretch/2
         print*, 'Stretching ratio in z: ', alpha
         print*, "Aspect ratio in uniform region: ", dx/dz
         print*, "Number of cells per droplet diameter: ", nz*(ddrop/Lz)
         print*, "========================================"
      end if
      
      ! General serial grid object
      if (nz.gt.1)then
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='ShockDrop')
      else
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.true.,name='ShockDrop')
      end if
      ! Read in partition
      call param_read('Partition',partition,short='p')
      ! Create partitioned grid
      this%cfg=config(grp=group,decomp=partition,grid=grid)
    end block  create_sdrop_config
  end subroutine init_grid

  subroutine init(this,dt_in,dtmax_in)
    implicit none
    class(sdrop), intent(inout) :: this
    real(WP), intent(in) :: dt_in, dtmax_in
    !> set up timetracker
    initialize_timetracker: block
      use param,             only: param_read
      use shockgen_class,    only: sgen
      this%time=timetracker(amRoot=this%cfg%amRoot)
      call param_read('Max time',this%time%tmax)
      call param_read('Max cfl number',this%time%cflmax)
      call param_read('Max steps',this%time%nmax)
      this%time%itmax=2
      this%time%dt = dt_in       ! shockgen%time%dt passed as arguments for init subroutine
      this%time%dtmax = dtmax_in ! shockgen%time%dtmax passed as arguments for init subroutine
    end block initialize_timetracker
    
    !> handle saves here
    save: block
      use param,   only: param_read
      use string,  only: str_medium
      use filesys, only: makedir,isdir
      
      character(len=str_medium) :: timestamp
      integer, dimension(3)     :: iopartition
      
      ! create event for saving restart files
      this%save_evt=event(this%time,'Restart output')
      call param_read('Restart output period', this%save_evt%tper)
      ! check if we are restarting
      call param_read('Restart from', timestamp, default='')
      this%restarted=.false.; if(len_trim(timestamp).gt.0) this%restarted=.true.
      ! read in I/O partition
      call param_read('I/O partition',iopartition)

      if (this%cfg%amRoot) then
         if(.not.isdir('restart')) call makedir('restart') 
      end if
      ! prepare pardata object for saving restart files
      call this%df%initialize(pg=this%cfg,iopartition=iopartition,filename=trim(this%cfg%name),nval=2,nvar=33)
      this%df%valname=['t ','dt']
      
      this%df%varname=['Grho   ','Lrho   ','RHO    ','Ui     ','Vi     ','Wi     ','U      ','V      ','W      ','rhoUi  ','rhoVi  ','rhoWi  ','GrhoE  ','LrhoE  ','GP     ','LP     ','P      ','PA     ','Tmptr  ','GrhoSS2','LrhoSS2','RHOSS2 ','P11    ','P12    ','P13    ','P14    ','P21    ','P22    ','P23    ','P24    ','SL_x   ','SL_y   ','SL_z   ']

    end block save
    
    ! Initialize our VOF solver and field
    create_VOF_solver: block
      use mms_geom, only: cube_refine_vol
      use vfs_class, only: VFhi,VFlo,plicnet,flux,neumann
      use irl_fortran_interface 
      
      integer   :: i,j,k,n,si,sj,sk
      real(WP), dimension(3,8) :: cube_vertex
      real(WP), dimension(3) :: v_cent,a_cent
      real(WP)  :: vol,area
      integer,  parameter :: amr_ref_lvl=4

      ! Create a VOF solver with PLICnet reconstruction
      call this%vf%initialize(cfg=this%cfg,reconstruction_method=plicnet,transport_method=flux,name='VOF')

      ! Initialize liquid
      call param_read('Droplet diameter',ddrop)
      call param_read('Droplet location',dctr)
      do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
         do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
            do i=this%vf%cfg%imino_,this%vf%cfg%imaxo_
               ! Set cube vertices
               n=0
               do sk=0,1
                  do sj=0,1
                     do si=0,1
                        n=n+1; cube_vertex(:,n)=[this%vf%cfg%x(i+si),this%vf%cfg%y(j+sj),this%vf%cfg%z(k+sk)]
                     end do
                  end do
               end do
               ! Call adaptive refinement code to get volume and barycenters recursively
               vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
               if (this%vf%cfg%nz.eq.1) then
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_cyl,0.0_WP,amr_ref_lvl)
               else
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_sphere,0.0_WP,amr_ref_lvl)
               end if
               this%vf%VF(i,j,k)=vol/this%vf%cfg%vol(i,j,k)
               if (this%vf%VF(i,j,k).ge.VFlo.and.this%vf%VF(i,j,k).le.VFhi) then
                  this%vf%Lbary(:,i,j,k)=v_cent
                  this%vf%Gbary(:,i,j,k)=([this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]-this%vf%VF(i,j,k)*this%vf%Lbary(:,i,j,k))/(1.0_WP-this%vf%VF(i,j,k))
                  if (this%vf%cfg%nz.eq.1) this%vf%Gbary(3,i,j,k)=v_cent(3);
               else
                  this%vf%Lbary(:,i,j,k)=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
                  this%vf%Gbary(:,i,j,k)=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
               end if
            end do
         end do
      end do
      
      ! Update the band
      call this%vf%update_band()
      ! Perform interface reconstruction from VOF field
      call this%vf%build_interface()
      ! Set initial interface at the boundaries
      call this%vf%set_full_bcond()
      ! Create discontinuous polygon mesh from IRL interface
      call this%vf%polygonalize_interface()
      ! Calculate distance from polygons
      call this%vf%distance_from_polygon()
      ! Calculate subcell phasic volumes
      call this%vf%subcell_vol()
      ! Calculate curvature
      call this%vf%get_curvature()
      ! Reset moments to guarantee compatibility with interface reconstruction
      call this%vf%reset_volume_moments()
      
      ! add boundary conditions on VOF
      call this%vf%add_bcond(name='xright',type=neumann,locator=right_of_domain,dir='xp')
      call this%vf%add_bcond(name='ytop',type=neumann,locator=top_of_domain,dir='yp')
      call this%vf%add_bcond(name='ybottom',type=neumann,locator=bot_of_domain,dir='ym')
      if (this%cfg%nz.gt.1)then
         call this%vf%add_bcond(name='zfront',type=neumann,locator=fnt_of_domain,dir='zp')
         call this%vf%add_bcond(name='zback',type=neumann,locator=bck_of_domain,dir='zm')
      end if
      ! apply boundary conditions on VOF
      call this%vf%apply_bcond(this%time%t,this%time%dt)
    end block create_VOF_solver
    
    !> create two-phase compressible flow solver
    create_flow_solver: block
      use hypre_str_class, only: pcg_pfmg   ! preconditioned conjugate gradient method for pressure and velocity
      use param,           only: param_read
      ! Create flow solver
      this%fs=mast(cfg=this%cfg,name='Two-phase All-Mach',vf=this%vf)
      ! Configure pressure solver
      this%ps=hypre_str(cfg=this%cfg,name='Pressure',method=pcg_pfmg,nst=7)
      this%ps%maxlevel=10
      call param_read('Pressure iteration',this%ps%maxit)
      call param_read('Pressure tolerance',this%ps%rcvg)
      ! Configure implicit velocity solver
      this%vs=hypre_str(cfg=this%cfg,name='Velocity',method=pcg_pfmg,nst=7)
      call param_read('Implicit iteration',this%vs%maxit)
      call param_read('Implicit tolerance',this%vs%rcvg)
      ! Setup the solver
      call this%fs%setup(pressure_solver=this%ps,implicit_solver=this%vs)
    end block create_flow_solver
    
    !> set initial and boundary  conditions
    set_IC_BC: block
      use mast_class,      only: bc_scope,bcond,mech_egy_mech_hhz,dirichlet,neumann,clipped_neumann ! boundary conditions
      use param,           only: param_read
      use mathtools,       only: Pi
      use parallel,        only: amRoot
      use messager,        only: die
      integer :: i,j,k,n,nx
      real(WP), dimension(3) :: xyz
      real(WP) :: gamm_l,Pref_l,gamm_g,visc_l,visc_g,Pref,cv_l0,cv_g0,kappa_l,kappa_g
      real(WP) :: vshock,relshockvel,Lx,dx
      real(WP) :: Grho0, GP0, Grho1, GP1, ST, Ma1, Ma, Lrho0, LP0, Mas
      type(bcond), pointer :: mybc
      
      ! Create material model class
      this%matmod=matm(cfg=this%cfg,name='Liquid-gas models')
      
      ! Get EOS parameters from input
      call param_read('Liquid Pref', Pref_l)
      call param_read('Liquid gamma',gamm_l)
      call param_read('Gas gamma',gamm_g)
      
      ! Register equations of state
      call this%matmod%register_stiffenedgas('liquid',gamm_l,Pref_l)
      call this%matmod%register_idealgas('gas',gamm_g)
      
      ! get viscosity, thermal conductivity, and specific heat for phases
      call param_read('Liquid dynamic viscosity',visc_l)
      call param_read('Gas dynamic viscosity',visc_g)
      call param_read('Liquid thermal conductivity',kappa_l)
      call param_read('Gas thermal conductivity',kappa_g)
      call param_read('Liquid specific heat (constant vol)',cv_l0)
      call param_read('Gas specific heat (constant vol)',cv_g0)
      
      ! Register flow solver variables with material models
      call this%matmod%register_thermoflow_variables('liquid',this%fs%Lrho,this%fs%Ui,this%fs%Vi,this%fs%Wi,this%fs%LrhoE,this%fs%LP)
      call this%matmod%register_thermoflow_variables('gas'   ,this%fs%Grho,this%fs%Ui,this%fs%Vi,this%fs%Wi,this%fs%GrhoE,this%fs%GP)
      call this%matmod%register_diffusion_thermo_models(viscconst_gas=visc_g, viscconst_liquid=visc_l,hdffconst_gas=kappa_g, hdffconst_liquid=kappa_l,sphtconst_gas=cv_g0,sphtconst_liquid=cv_l0)
      
      ! Read in surface tension coefficient
      call param_read('Surface tension coefficient',this%fs%sigma)
      
      ! Liquid and gas density
      call param_read('Liquid density',Lrho0);
      call param_read('Pre-shock density',Grho0,default=1.204_WP)
      call param_read('Pre-shock pressure',GP0,default=1.01325e5_WP)
      call param_read('Mach number of shock',Ma,default=1.47_WP)
      
      !use shock relations to get post shock numbers
      GP1 = GP0 * (2.0_WP*gamm_g*Ma**2 - (gamm_g-1.0_WP)) / (gamm_g+1.0_WP)
      Grho1 = Grho0 * (Ma**2 * (gamm_g+1.0_WP) / ((gamm_g-1.0_WP)*Ma**2 + 2.0_WP))
      !calculate post shock Mach number (mach number of gas behind shock)
      Ma1 = sqrt(((gamm_g-1.0_WP)*(Ma**2)+2.0_WP)/(2.0_WP*gamm_g*(Ma**2)-(gamm_g-1.0_WP)))
      !calculate post shock velocity (velocity of the gas behind the shock)
      vshock = -Ma1 * sqrt(gamm_g*GP1/Grho1) + Ma*sqrt(gamm_g*GP0/Grho0)
      !velocity at which the shock moves
      relshockvel = -Grho1*vshock/(Grho0-Grho1)
      dx = Lx/nx ! mesh spacing in uniform region
      
      if (amRoot) then
         print*, "===== Problem Setup Description ====="
         print*, 'Mach number', Ma
         print*, 'Pre-shock:  Density',Grho0,'Pressure',GP0
         print*, 'Post-shock: Density',Grho1,'Pressure',GP1,'Gas Velocity',vshock
         print*, 'Shock velocity', relshockvel
         print*, "===================================="
      end if
      
      ! if no restart, initialize cell centered velocities in y and z to zero 
      this%fs%Vi = 0.0_WP; this%fs%Wi = 0.0_WP  
      ! Zero face velocities as well for the sake of dirichlet boundaries
      this%fs%V = 0.0_WP; this%fs%W = 0.0_WP
      
      call param_read('Shock location',this%xshock)
      call param_read('Droplet diameter',ddrop)
      
      ! shock discontinuity initialization
      do i=this%fs%cfg%imino_,this%fs%cfg%imaxo_
         if (this%cfg%xm(i).le.this%xshock) then
            this%fs%Grho(i,:,:) = Grho1
            this%fs%Ui(i,:,:) = vshock
            this%fs%GP(i,:,:) = GP1
            this%fs%GrhoE(i,:,:) = this%matmod%EOS_energy(GP1,Grho1,vshock,0.0_WP,0.0_WP,'gas')
         elseif (this%cfg%xm(i).ge.this%xshock) then
            this%fs%Grho(i,:,:) = Grho0
            this%fs%Ui(i,:,:) = 0.0_WP
            this%fs%GP(i,:,:) = GP0
            this%fs%GrhoE(i,:,:) = this%matmod%EOS_energy(GP0,Grho0,0.0_WP,0.0_WP,0.0_WP,'gas')
         end if
      end do

      !!! shock profile is handled in simulation.f90
      
      if (this%fs%cfg%nz.eq.1) then
         ! Cylinder configuration, curv = 1/r
         LP0 = GP0 + 2.0/ddrop*this%fs%sigma
      else
         ! Sphere configuration, curv = 1/r + 1/r
         LP0 = GP0 + 4.0/ddrop*this%fs%sigma
      end if
      
      !initialize liquid quantities
      this%fs%Lrho = Lrho0
      this%fs%LP = LP0
      
      ! Initialize liquid energy, with surface tension
      this%fs%LrhoE = this%matmod%EOS_energy(LP0,Lrho0,0.0_WP,0.0_WP,0.0_WP,'liquid')
      
      ! Define boundary conditions - initialized values are intended dirichlet values too, for the cell centers
      call this%fs%add_bcond(name= 'inflow',type=dirichlet      ,locator=left_of_domain ,face='x',dir=-1)
      call this%fs%add_bcond(name='outflow',type=clipped_neumann,locator=right_of_domain,face='x',dir=+1)
      call this%fs%add_bcond(name='outflow',type=clipped_neumann,locator=bot_of_domain,face='y',dir=-1)
      call this%fs%add_bcond(name='outflow',type=clipped_neumann,locator=top_of_domain,face='y',dir=+1)
      if (this%cfg%nz.gt.1)then
         call this%fs%add_bcond(name='outflow',type=clipped_neumann,locator=bck_of_domain,face='z',dir=-1)
         call this%fs%add_bcond(name='outflow',type=clipped_neumann,locator=fnt_of_domain,face='z',dir=+1)
      end if

      ! Calculate face velocities
      call this%fs%interp_vel_basic(this%vf,this%fs%Ui,this%fs%Vi,this%fs%Wi,this%fs%U,this%fs%V,this%fs%W)
      
      ! Apply face BC - inflow
      call this%fs%get_bcond('inflow',mybc)
      do n=1,mybc%itr%n_
         i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
         this%fs%U(i,j,k)=vshock
      end do

      ! Apply face BC - outflow
      bc_scope = 'velocity'
      call this%fs%apply_bcond(this%time%dt,bc_scope)
      
      ! Calculate mixture density and momenta
      this%fs%RHO   = (1.0_WP-this%vf%VF)*this%fs%Grho  + this%vf%VF*this%fs%Lrho
      this%fs%rhoUi = this%fs%RHO*this%fs%Ui; this%fs%rhoVi = this%fs%RHO*this%fs%Vi; this%fs%rhoWi = this%fs%RHO*this%fs%Wi
      
      ! set pressure relax model
      this%relax_model = mech_egy_mech_hhz
      ! Perform initial pressure relax
      call this%fs%pressure_relax(this%vf,this%matmod,this%relax_model)
      
      ! Calculate initial phase and bulk moduli
      call this%fs%init_phase_bulkmod(this%vf,this%matmod)
      call this%fs%reinit_phase_pressure(this%vf,this%matmod)
      call this%fs%harmonize_advpressure_bulkmod(this%vf,this%matmod)
      
      ! Set initial pressure to harmonized field based on internal energy
      this%fs%P = this%fs%PA
      
    end block set_IC_BC
    
    !> create surfmesh object for interface polygon output
    create_smesh: block
      use irl_fortran_interface
      integer :: i,j,k,nplane,np
      this%smesh=surfmesh(nvar=1,name='plic')
      this%smesh%varname(1)='curv'
      call this%vf%update_surfmesh(this%smesh)
      this%smesh%var(1,:)=0.0_WP
      np=0;
      do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
         do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
            do i=this%vf%cfg%imin_,this%vf%cfg%imax_
               do nplane=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                  if (getNumberOfVertices(this%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                     np=np+1;
                     this%smesh%var(1,np)=this%vf%curv(i,j,k)
                  end if
               end do
            end do
         end do
      end do
    end block create_smesh
    
    !> Add Ensight output
    create_ensight: block
      use param,           only: param_read
      ! Create Ensight output from cfg
      this%ens_out=ensight(cfg=this%cfg,name='ShockDroplet')
      ! Create event for Ensight output
      this%ens_evt=event(time=this%time,name='Ensight output')
      call param_read('Ensight output period',this%ens_evt%tper)
      ! Add variables to output
      call this%ens_out%add_vector('velocity',this%fs%Ui,this%fs%Vi,this%fs%Wi)
      call this%ens_out%add_scalar('P',this%fs%P)
      call this%ens_out%add_scalar('PA',this%fs%PA)
      call this%ens_out%add_scalar('Grho',this%fs%Grho)
      call this%ens_out%add_scalar('Lrho',this%fs%Lrho)
      call this%ens_out%add_scalar('Density',this%fs%RHO)
      call this%ens_out%add_scalar('Bulkmod',this%fs%RHOSS2)
      call this%ens_out%add_scalar('VOF',this%vf%VF)
      call this%ens_out%add_scalar('curvature',this%vf%curv)
      call this%ens_out%add_scalar('Mach',this%fs%Mach)
      call this%ens_out%add_scalar('fvf',this%cfg%VF)
      call this%ens_out%add_scalar('Tmptr',this%fs%Tmptr)
      call this%ens_out%add_scalar('SL_x',this%fs%sl_x) 
      call this%ens_out%add_scalar('SL_y',this%fs%sl_y) 
      call this%ens_out%add_scalar('SL_z',this%fs%sl_z) 
      call this%ens_out%add_scalar('LP',this%fs%LP) 
      call this%ens_out%add_scalar('GP',this%fs%GP) 
      call this%ens_out%add_scalar('LrhoE',this%fs%LrhoE) 
      call this%ens_out%add_scalar('GrhoE',this%fs%GrhoE)         
      ! Output to ensight
      !if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
    end block create_ensight
    
    !> block for writing smesh data
    create_ensight_smesh: block
      use param,           only: param_read
      real(WP) :: smesh_tper ! declare variable for smesh output frequency
      call param_read('Ensight smesh output period', smesh_tper)      
      ! create ensight output from cfg for smesh surface reconstruction
      this%ens_out_smesh=ensight(cfg=this%cfg,name='droplet_smesh')      
      ! create event for ensight output
      this%ens_evt_smesh=event(time=this%time,name='Ensight output smesh')
      this%ens_evt_smesh%tper = smesh_tper
      ! add variables to output
      call this%ens_out_smesh%add_surface('smesh',this%smesh)
      ! output to ensight
      !if (this%ens_evt_smesh%occurs()) call this%ens_out_smesh%write_data(this%time%t)
    end block create_ensight_smesh
    
    !> Create a monitor file
    create_monitor: block
      ! Prepare some info about fields
      call this%fs%get_cfl(this%time%dt,this%time%cfl)
      call this%fs%get_max()
      call this%vf%get_max()
      ! Create simulation monitor
      this%mfile=monitor(this%fs%cfg%amRoot,'simulation')
      call this%mfile%add_column(this%time%n,'Timestep number')
      call this%mfile%add_column(this%time%t,'Time')
      call this%mfile%add_column(this%time%dt,'Timestep size')
      call this%mfile%add_column(this%time%cfl,'Maximum CFL')
      call this%mfile%add_column(this%fs%RHOmin,'RHOmin')
      call this%mfile%add_column(this%fs%RHOmax,'RHOmax')
      call this%mfile%add_column(this%fs%Umax,'Umax')
      call this%mfile%add_column(this%fs%Vmax,'Vmax')
      call this%mfile%add_column(this%fs%Wmax,'Wmax')
      call this%mfile%add_column(this%fs%Pmax,'Pmax')
      call this%mfile%add_column(this%fs%Tmax,'Tmax')
      call this%mfile%write()
      ! Create CFL monitor
      this%cflfile=monitor(this%fs%cfg%amRoot,'cfl')
      call this%cflfile%add_column(this%time%n,'Timestep number')
      call this%cflfile%add_column(this%time%t,'Time')
      call this%cflfile%add_column(this%fs%CFLst,'STension CFL')
      call this%cflfile%add_column(this%fs%CFLc_x,'Convective xCFL')
      call this%cflfile%add_column(this%fs%CFLc_y,'Convective yCFL')
      call this%cflfile%add_column(this%fs%CFLc_z,'Convective zCFL')
      call this%cflfile%add_column(this%fs%CFLv_x,'Viscous xCFL')
      call this%cflfile%add_column(this%fs%CFLv_y,'Viscous yCFL')
      call this%cflfile%add_column(this%fs%CFLv_z,'Viscous zCFL')
      call this%cflfile%add_column(this%fs%CFLa_x,'Acoustic xCFL')
      call this%cflfile%add_column(this%fs%CFLa_y,'Acoustic yCFL')
      call this%cflfile%add_column(this%fs%CFLa_z,'Acoustic zCFL')
      call this%cflfile%write()
      ! Create convergence monitor
      this%cvgfile=monitor(this%fs%cfg%amRoot,'cvg')
      call this%cvgfile%add_column(this%time%n,'Timestep number')
      call this%cvgfile%add_column(this%time%it,'Iteration')
      call this%cvgfile%add_column(this%time%t,'Time')
      call this%cvgfile%add_column(this%fs%impl_it_x,'Impl_x iteration')
      call this%cvgfile%add_column(this%fs%impl_rerr_x,'Impl_x error')
      call this%cvgfile%add_column(this%fs%impl_it_y,'Impl_y iteration')
      call this%cvgfile%add_column(this%fs%impl_rerr_y,'Impl_y error')
      call this%cvgfile%add_column(this%fs%implicit%it,'Impl_z iteration')
      call this%cvgfile%add_column(this%fs%implicit%rerr,'Impl_z error')
      call this%cvgfile%add_column(this%fs%psolv%it,'Pressure iteration')
      call this%cvgfile%add_column(this%fs%psolv%rerr,'Pressure error')
    end block create_monitor
  end subroutine init
  
  !> write initial condition to ensight
  subroutine writeIC(this)
    implicit none
    class(sdrop), intent(inout) :: this
    if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
    if (this%ens_evt_smesh%occurs()) call this%ens_out_smesh%write_data(this%time%t)
  end subroutine writeIC
  
  subroutine update_mixture_variables(this)
    use mast_class, only: mech_egy_mech_hhz
    implicit none
    class(sdrop), intent(inout) :: this
    
    ! Calculate mixture density and momenta
    this%fs%RHO   = (1.0_WP-this%vf%VF)*this%fs%Grho  + this%vf%VF*this%fs%Lrho
    this%fs%rhoUi = this%fs%RHO*this%fs%Ui; this%fs%rhoVi = this%fs%RHO*this%fs%Vi; this%fs%rhoWi = this%fs%RHO*this%fs%Wi
    
    ! set pressure relax model
    this%relax_model = mech_egy_mech_hhz
    ! Perform initial pressure relax
    call this%fs%pressure_relax(this%vf,this%matmod,this%relax_model)
    
    ! Calculate initial phase and bulk moduli
    call this%fs%init_phase_bulkmod(this%vf,this%matmod)
    call this%fs%reinit_phase_pressure(this%vf,this%matmod)
    call this%fs%harmonize_advpressure_bulkmod(this%vf,this%matmod)
  end subroutine update_mixture_variables
  
  !> advance shock drop simualation by dt
  subroutine step(this)
    use messager, only: die
    implicit none
    class(sdrop), intent(inout) :: this
    
    ! Perform time integration in simulation.f90 file 
    ! Increment time
    call this%fs%get_cfl(this%time%dt,this%time%cfl)
    call this%time%adjust_dt()
    call this%time%increment()
    
    ! Reinitialize phase pressure by syncing it with conserved phase energy
    call this%fs%reinit_phase_pressure(this%vf,this%matmod)
    ! remember old velocity and density
    this%fs%Uiold=this%fs%Ui; this%fs%Viold=this%fs%Vi; this%fs%Wiold=this%fs%Wi;this%fs%RHOold = this%fs%RHO
    
    ! Remember old flow variables (phase)
    this%fs%Grhoold = this%fs%Grho; this%fs%Lrhoold = this%fs%Lrho
    this%fs%GrhoEold=this%fs%GrhoE; this%fs%LrhoEold=this%fs%LrhoE
    this%fs%GPold   =   this%fs%GP; this%fs%LPold   =   this%fs%LP
    
    ! Remember old interface, including VF and barycenters
    call this%vf%copy_interface_to_old()
    
    ! Create in-cell reconstruction
    call this%fs%flow_reconstruct(this%vf)
    
    ! Zero variables that will change during subiterations
    this%fs%P = 0.0_WP;this%fs%Pjx = 0.0_WP;this%fs%Pjy = 0.0_WP;this%fs%Pjz = 0.0_WP;this%fs%Hpjump = 0.0_WP
    
    ! Determine semi-Lagrangian advection flag
    call this%fs%flag_sl(this%time%dt,this%vf)
    
    ! Perform sub-iterations
    do while (this%time%it.le.this%time%itmax)
       
       ! Predictor step, involving advection and pressure terms
       call this%fs%advection_step(this%time%dt,this%vf,this%matmod)
       
       ! Viscous step
       call this%fs%diffusion_src_explicit_step(this%time%dt,this%vf,this%matmod)
       
       ! Prepare pressure projection
       call this%fs%pressureproj_prepare(this%time%dt,this%vf,this%matmod)
       
          ! Initialize and solve Helmholtz equation
       call this%fs%psolv%setup()
       this%fs%psolv%sol=this%fs%PA-this%fs%P
       call this%fs%psolv%solve()
       call this%fs%cfg%sync(this%fs%psolv%sol)
       
       ! Perform corrector step using solution
       this%fs%P=this%fs%P+this%fs%psolv%sol
       
       call this%fs%pressureproj_correct(this%time%dt,this%vf,this%fs%psolv%sol)
       
       ! Record convergence monitor
       call this%cvgfile%write()
       ! Increment sub-iteration counter
       this%time%it=this%time%it+1
       
    end do

    ! Pressure relaxation
    call this%fs%pressure_relax(this%vf,this%matmod,this%relax_model)
    
    ! Output to ensight
    if (this%ens_evt%occurs()) then
       call this%ens_out%write_data(this%time%t)            
    end if
    
    if (this%ens_evt_smesh%occurs()) then
       !update surfmesh object
       update_smesh: block
         use irl_fortran_interface
         integer :: i,j,k,nplane,np
         ! Transfer polygons to smesh
         call this%vf%update_surfmesh(this%smesh)
         ! Also populate nplane variable
         this%smesh%var(1,:)=0.0_WP
         np=0
         do k= this%vf%cfg%kmin_, this%vf%cfg%kmax_
            do j= this%vf%cfg%jmin_, this%vf%cfg%jmax_
               do i= this%vf%cfg%imin_, this%vf%cfg%imax_
                  do nplane=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                     if (getNumberOfVertices(this%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                        np=np+1;
                        this%smesh%var(1,np)=this%vf%curv(i,j,k)
                     end if
                  end do
               end do
            end do
         end do
       end block update_smesh
       call this%ens_out_smesh%write_data(this%time%t)
    end if
    
    ! Perform and output monitoring
    call this%fs%get_max()
    call this%vf%get_max()
    call this%fs%get_viz()
    call this%mfile%write()
    call this%cflfile%write()
    
    if (this%save_evt%occurs())then
       save_restart: block
         use irl_fortran_interface
         use string, only: str_medium
         character(len=str_medium) :: timestamp
         real(WP), dimension(:,:,:), allocatable :: P11,P12,P13,P14
         real(WP), dimension(:,:,:), allocatable :: P21,P22,P23,P24
         integer :: i,j,k
         real(WP), dimension(4) :: plane
         
         ! Handle IRL data
         allocate(P11(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(P12(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(P13(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(P14(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(P21(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(P22(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(P23(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(P24(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         ! Store IRL data
         do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
            do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
               do i=this%vf%cfg%imino_,this%vf%cfg%imaxo_
                  ! First plane
                  plane=getPlane(this%vf%liquid_gas_interface(i,j,k),0)
                  P11(i,j,k)=plane(1); P12(i,j,k)=plane(2); P13(i,j,k)=plane(3); P14(i,j,k)=plane(4)
                  ! Second plane
                  plane=0.0_WP
                  if (getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k)).eq.2) plane=getPlane(this%vf%liquid_gas_interface(i,j,k),1)
                  P21(i,j,k)=plane(1); P22(i,j,k)=plane(2); P23(i,j,k)=plane(3); P24(i,j,k)=plane(4)
               end do
            end do
         end do
         ! Prefix for files
         write(timestamp,'(es12.5)') this%time%t
         ! Populate df and write it
         call this%df%push(name='t'      ,val=this%time%t    )
         call this%df%push(name='dt'     ,val=this%time%dt   )
         call this%df%push(name='Grho'   ,var=this%fs%Grho   )
         call this%df%push(name='Lrho'   ,var=this%fs%Lrho   )
         call this%df%push(name='RHO'    ,var=this%fs%RHO    )
         call this%df%push(name='Ui'     ,var=this%fs%Ui     )
         call this%df%push(name='Vi'     ,var=this%fs%Vi     )
         call this%df%push(name='Wi'     ,var=this%fs%Wi     )
         call this%df%push(name='U'      ,var=this%fs%U      )
         call this%df%push(name='V'      ,var=this%fs%V      )
         call this%df%push(name='W'      ,var=this%fs%W      )
         call this%df%push(name='rhoUi'  ,var=this%fs%rhoUi  )
         call this%df%push(name='rhoVi'  ,var=this%fs%rhoVi  )
         call this%df%push(name='rhoWi'  ,var=this%fs%rhoWi  )
         call this%df%push(name='GrhoE'  ,var=this%fs%GrhoE  )
         call this%df%push(name='LrhoE'  ,var=this%fs%LrhoE  )
         call this%df%push(name='GP'     ,var=this%fs%GP     )
         call this%df%push(name='LP'     ,var=this%fs%LP     )
         call this%df%push(name='P'      ,var=this%fs%P      ) 
         call this%df%push(name='PA'     ,var=this%fs%PA     )
         call this%df%push(name='Tmptr'  ,var=this%fs%Tmptr  )
         call this%df%push(name='GrhoSS2',var=this%fs%GrhoSS2) 
         call this%df%push(name='LrhoSS2',var=this%fs%LrhoSS2) 
         call this%df%push(name='RHOSS2' ,var=this%fs%RHOSS2 )
         call this%df%push(name='P11'    ,var=P11       )
         call this%df%push(name='P12'    ,var=P12       )
         call this%df%push(name='P13'    ,var=P13       )
         call this%df%push(name='P14'    ,var=P14       )
         call this%df%push(name='P21'    ,var=P21       )
         call this%df%push(name='P22'    ,var=P22       )
         call this%df%push(name='P23'    ,var=P23       )
         call this%df%push(name='P24'    ,var=P24       )
         call this%df%push(name='SL_x'   ,var=this%fs%sl_x   )
         call this%df%push(name='SL_y'   ,var=this%fs%sl_y   )
         call this%df%push(name='SL_z'   ,var=this%fs%sl_z   )
         call this%df%write(fdata='restart/data_'//trim(adjustl(timestamp)))
         ! Deallocate
         deallocate(P11,P12,P13,P14,P21,P22,P23,P24)
       end block save_restart
    end if
  end subroutine step
  
  !> restart the simulation
  ! (this will probably need to be cleaned up a bit, but this is a self contained restart capability)
  ! data file saving logic will remain inside of the main shockdrop subroutines since we always want to save data
  subroutine restart(this)
    use param,   only: param_read
    use string,  only: str_medium
    use filesys, only: makedir,isdir
    implicit none
    class(sdrop), intent(inout) :: this
    
    character(len=str_medium) :: timestamp
    integer, dimension(3)     :: iopartition

    call param_read('Restart from', timestamp, default=''); this%restarted=.true.
    ! read in I/O partition
    call param_read('I/O partition',iopartition)
    ! we are restarting, read file
    call this%df%initialize(pg=this%cfg,iopartition=iopartition,fdata='restart/data_'//trim(adjustl(timestamp)))
    !> setup saved data files (we always want to save data files)
    ! create event for saving restart files
    this%save_evt=event(this%time,'Restart output')
    call param_read('Restart output period', this%save_evt%tper)

    !> set up timetracker
    restart_timetracker: block
      use param,             only: param_read
      use shockgen_class,    only: sgen
      this%time=timetracker(amRoot=this%cfg%amRoot)
      call param_read('Max time',this%time%tmax)
      call param_read('Max timestep size',this%time%dtmax)
      call param_read('Max cfl number',this%time%cflmax)
      call param_read('Max steps',this%time%nmax)
      this%time%itmax=2
      ! update timetracker
      call this%df%pull(name='t' ,val=this%time%t )
      call this%df%pull(name='dt',val=this%time%dt)
      this%time%told=this%time%t-this%time%dt
    end block restart_timetracker

    !> Initialize our VOF solver and field
    restart_VOF_solver: block
      use mms_geom, only: cube_refine_vol
      use vfs_class, only: VFhi,VFlo,plicnet,flux,neumann
      use irl_fortran_interface 
      
      integer   :: i,j,k,n,si,sj,sk
      real(WP), dimension(3,8) :: cube_vertex
      real(WP), dimension(3) :: v_cent,a_cent
      real(WP)  :: vol,area
      integer,  parameter :: amr_ref_lvl=4
      real(WP), dimension(:,:,:), allocatable :: P11,P12,P13,P14 
      real(WP), dimension(:,:,:), allocatable :: P21,P22,P23,P24
      
      ! Create a VOF solver with PLICnet reconstruction
      call this%vf%initialize(cfg=this%cfg,reconstruction_method=plicnet,transport_method=flux,name='VOF')
      
      ! initialize the interface including restarts         
      ! Read in the planes directly and set the IRL interface
      allocate(P11(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P11',var=P11)
      allocate(P12(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P12',var=P12)
      allocate(P13(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P13',var=P13)
      allocate(P14(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P14',var=P14)
      allocate(P21(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P21',var=P21)
      allocate(P22(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P22',var=P22)
      allocate(P23(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P23',var=P23)
      allocate(P24(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P24',var=P24)
      do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
         do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
            do i=this%vf%cfg%imin_,this%vf%cfg%imax_
               ! Check if the second plane is meaningful
               if (this%vf%two_planes.and.P21(i,j,k)**2+P22(i,j,k)**2+P23(i,j,k)**2.gt.0.0_WP) then
                  call setNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k),2)
                  call setPlane(this%vf%liquid_gas_interface(i,j,k),0,[P11(i,j,k),P12(i,j,k),P13(i,j,k)],P14(i,j,k))
                  call setPlane(this%vf%liquid_gas_interface(i,j,k),1,[P21(i,j,k),P22(i,j,k),P23(i,j,k)],P24(i,j,k))
               else
                  call setNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k),1)
                  call setPlane(this%vf%liquid_gas_interface(i,j,k),0,[P11(i,j,k),P12(i,j,k),P13(i,j,k)],P14(i,j,k))
               end if
            end do
         end do
      end do
      
      call this%vf%sync_interface() ! this syncs the interface across processors
      deallocate(P11,P12,P13,P14,P21,P22,P23,P24)
      ! Reset moments to guarantee compatibility with interface reconstruction
      call this%vf%reset_volume_moments()! this resets volumetric moments based on the interface
      
      ! ensure ghost cells are correct
      if (this%cfg%iproc.eq.1) this%vf%VF(this%cfg%imino:this%cfg%imin-1,:,:)=0.0_WP 
      if (this%cfg%jproc.eq.this%cfg%npy) this%vf%VF(:,this%cfg%jmaxo:this%cfg%jmax+1,:)=0.0_WP
      if (this%cfg%jproc.eq.1) this%vf%VF(:,this%cfg%jmino:this%cfg%jmin-1,:)=0.0_WP
      if (this%cfg%kproc.eq.1) this%vf%VF(:,:,this%cfg%kmino:this%cfg%kmin-1)=0.0_WP
      if (this%cfg%kproc.eq.this%cfg%npz) this%vf%VF(:,:,this%cfg%kmaxo:this%cfg%kmax+1)=0.0_WP
      
      ! add boundary conditions on VOF
      call this%vf%add_bcond(name='xright',type=neumann,locator=right_of_domain,dir='xp')
      call this%vf%add_bcond(name='ytop',type=neumann,locator=top_of_domain,dir='yp')
      call this%vf%add_bcond(name='ybottom',type=neumann,locator=bot_of_domain,dir='ym')
      if (this%cfg%nz.gt.1)then
         call this%vf%add_bcond(name='zfront',type=neumann,locator=fnt_of_domain,dir='zp')
         call this%vf%add_bcond(name='zback',type=neumann,locator=bck_of_domain,dir='zm')
      end if
      ! apply boundary conditions on VOF
      call this%vf%apply_bcond(this%time%t,this%time%dt)
      
      ! Update the band
      call this%vf%update_band() ! searches for a fluid interface and updates the band 
      ! set interface at the boundaries
      call this%vf%set_full_bcond() ! sets the liquid/gas plane boundaries based on the VOF field
      ! Create discontinuous polygon mesh from IRL interface
      call this%vf%polygonalize_interface() ! creates a polygonal representation of the interface
      ! Calculate distance from polygons
      call this%vf%distance_from_polygon() ! calculates the distance from the interface within the band
      ! Calculate subcell phasic volumes
      call this%vf%subcell_vol() ! calcualtes the phase volumes within each cell based on the interface
      ! Calculate curvature
      call this%vf%get_curvature() ! calculates the curvature of the interface using a least squares fit
    end block restart_VOF_solver

    !> create two-phase compressible flow solver
    restart_flow_solver: block
      use hypre_str_class, only: pcg_pfmg   ! preconditioned conjugate gradient method for pressure and velocity
      use param,           only: param_read
      ! Create flow solver
      this%fs=mast(cfg=this%cfg,name='Two-phase All-Mach',vf=this%vf)
      ! Configure pressure solver
      this%ps=hypre_str(cfg=this%cfg,name='Pressure',method=pcg_pfmg,nst=7)
      this%ps%maxlevel=10
      call param_read('Pressure iteration',this%ps%maxit)
      call param_read('Pressure tolerance',this%ps%rcvg)
      ! Configure implicit velocity solver
      this%vs=hypre_str(cfg=this%cfg,name='Velocity',method=pcg_pfmg,nst=7)
      call param_read('Implicit iteration',this%vs%maxit)
      call param_read('Implicit tolerance',this%vs%rcvg)
      ! Setup the solver
      call this%fs%setup(pressure_solver=this%ps,implicit_solver=this%vs)
    end block restart_flow_solver

    !> apply our restarted boundary conditions and field variables
    restart_IC_BC: block
      use mast_class,      only: bc_scope,bcond,mech_egy_mech_hhz,dirichlet,neumann,clipped_neumann ! boundary conditions
      use param,           only: param_read
      use mathtools,       only: Pi
      use parallel,        only: amRoot
      use messager,        only: die
      integer :: i,j,k,n,nx
      real(WP), dimension(3) :: xyz
      real(WP) :: gamm_l,Pref_l,gamm_g,visc_l,visc_g,Pref,cv_l0,cv_g0,kappa_l,kappa_g
      real(WP) :: Grho0, GP0, Grho1, GP1, ST, Ma1, Ma, Lrho0, LP0, Mas, vshock, relshockvel
      type(bcond), pointer :: mybc

      if (amRoot) then
         print*, "===================================="
         print*, "============ RESTARTING ============"
         print*, "===================================="
      end if
      
      ! Create material model class
      this%matmod=matm(cfg=this%cfg,name='Liquid-gas models')
      
      ! Get EOS parameters from input
      call param_read('Liquid Pref', Pref_l)
      call param_read('Liquid gamma',gamm_l)
      call param_read('Gas gamma',gamm_g)
      
      ! Register equations of state
      call this%matmod%register_stiffenedgas('liquid',gamm_l,Pref_l)
      call this%matmod%register_idealgas('gas',gamm_g)
      
      ! get viscosity, thermal conductivity, and specific heat for phases
      call param_read('Liquid dynamic viscosity',visc_l)
      call param_read('Gas dynamic viscosity',visc_g)
      call param_read('Liquid thermal conductivity',kappa_l)
      call param_read('Gas thermal conductivity',kappa_g)
      call param_read('Liquid specific heat (constant vol)',cv_l0)
      call param_read('Gas specific heat (constant vol)',cv_g0)
      
      ! Register flow solver variables with material models
      call this%matmod%register_thermoflow_variables('liquid',this%fs%Lrho,this%fs%Ui,this%fs%Vi,this%fs%Wi,this%fs%LrhoE,this%fs%LP)
      call this%matmod%register_thermoflow_variables('gas'   ,this%fs%Grho,this%fs%Ui,this%fs%Vi,this%fs%Wi,this%fs%GrhoE,this%fs%GP)
      call this%matmod%register_diffusion_thermo_models(viscconst_gas=visc_g, viscconst_liquid=visc_l,hdffconst_gas=kappa_g, hdffconst_liquid=kappa_l,sphtconst_gas=cv_g0,sphtconst_liquid=cv_l0)
      
      ! Read in surface tension coefficient
      call param_read('Surface tension coefficient',this%fs%sigma)
      ! Liquid and gas density
      call param_read('Liquid density',Lrho0);
      call param_read('Pre-shock density',Grho0,default=1.204_WP)
      call param_read('Pre-shock pressure',GP0,default=1.01325e5_WP)
      call param_read('Mach number of shock',Ma,default=1.47_WP)

      ! set pressure relax model
      this%relax_model = mech_egy_mech_hhz

      !use shock relations to get post shock numbers
      GP1 = GP0 * (2.0_WP*gamm_g*Ma**2 - (gamm_g-1.0_WP)) / (gamm_g+1.0_WP)
      Grho1 = Grho0 * (Ma**2 * (gamm_g+1.0_WP) / ((gamm_g-1.0_WP)*Ma**2 + 2.0_WP))
      !calculate post shock Mach number (mach number of gas behind shock)
      Ma1 = sqrt(((gamm_g-1.0_WP)*(Ma**2)+2.0_WP)/(2.0_WP*gamm_g*(Ma**2)-(gamm_g-1.0_WP)))
      !calculate post shock velocity (velocity of the gas behind the shock)
      vshock = -Ma1 * sqrt(gamm_g*GP1/Grho1) + Ma*sqrt(gamm_g*GP0/Grho0)
      !velocity at which the shock moves
      relshockvel = -Grho1*vshock/(Grho0-Grho1)

      ! Read data
      ! when df%pull is called, it syncs the overlapping cells that are in the physical domain
      ! the overlapping ghost cells on the boundaries are not saved with the restart files and must be set 
      call this%df%pull(name='Grho'   ,var=this%fs%Grho   );
      call this%df%pull(name='Lrho'   ,var=this%fs%Lrho   ); 
      call this%df%pull(name='RHO'    ,var=this%fs%RHO    ); 
      call this%df%pull(name='Ui'     ,var=this%fs%Ui     ); ! cell centered value (see mast_class.f90)
      call this%df%pull(name='Vi'     ,var=this%fs%Vi     ); 
      call this%df%pull(name='Wi'     ,var=this%fs%Wi     ); 
      call this%df%pull(name='U'      ,var=this%fs%U      ); ! face value (see mast_class.f90)
      call this%df%pull(name='V'      ,var=this%fs%V      ); 
      call this%df%pull(name='W'      ,var=this%fs%W      ); 
      call this%df%pull(name='rhoUi'  ,var=this%fs%rhoUi  ); 
      call this%df%pull(name='rhoVi'  ,var=this%fs%rhoVi  ); 
      call this%df%pull(name='rhoWi'  ,var=this%fs%rhoWi  ); 
      call this%df%pull(name='GrhoE'  ,var=this%fs%GrhoE  ); 
      call this%df%pull(name='LrhoE'  ,var=this%fs%LrhoE  ); 
      call this%df%pull(name='GP'     ,var=this%fs%GP     ); 
      call this%df%pull(name='LP'     ,var=this%fs%LP     ); 
      call this%df%pull(name='P'      ,var=this%fs%P      ); 
      call this%df%pull(name='PA'     ,var=this%fs%PA     );
      call this%df%pull(name='Tmptr'  ,var=this%fs%Tmptr  );
      call this%df%pull(name='GrhoSS2',var=this%fs%GrhoSS2); 
      call this%df%pull(name='LrhoSS2',var=this%fs%LrhoSS2); 
      call this%df%pull(name='RHOSS2' ,var=this%fs%RHOSS2 );
      call this%df%pull(name='SL_x'   ,var=this%fs%sl_x   ); ! AS pull sensor data
      call this%df%pull(name='SL_y'   ,var=this%fs%sl_y   );
      call this%df%pull(name='SL_z'   ,var=this%fs%sl_z   );
      
      ! define boundary conditions
      call this%fs%add_bcond(name= 'inflow',type=dirichlet      ,locator=left_of_domain ,face='x',dir=-1)
      call this%fs%add_bcond(name='outflow',type=clipped_neumann,locator=right_of_domain,face='x',dir=+1)
      call this%fs%add_bcond(name='outflow',type=clipped_neumann,locator=bot_of_domain,face='y',dir=-1)
      call this%fs%add_bcond(name='outflow',type=clipped_neumann,locator=top_of_domain,face='y',dir=+1)
      if (this%cfg%nz.gt.1)then
         call this%fs%add_bcond(name='outflow',type=clipped_neumann,locator=bck_of_domain,face='z',dir=-1)
         call this%fs%add_bcond(name='outflow',type=clipped_neumann,locator=fnt_of_domain,face='z',dir=+1)
      end if
      
      ! Ensure that we are only on the boundaries
      ! inlet at left of domain
      if (this%cfg%iproc.eq.1) then ! apply inlet (dirichlet) condition
         
         !try using post shock values for simplicity
         this%fs%U(this%cfg%imino:this%cfg%imin-1,:,:) = vshock; this%fs%V(this%cfg%imino:this%cfg%imin-1,:,:) = 0.0_WP; this%fs%W(this%cfg%imino:this%cfg%imin-1,:,:) = 0.0_WP
         this%fs%Ui(this%cfg%imino:this%cfg%imin-1,:,:) = vshock; this%fs%Vi(this%cfg%imino:this%cfg%imin-1,:,:) = 0.0_WP; this%fs%Wi(this%cfg%imino:this%cfg%imin-1,:,:) = 0.0_WP
         this%fs%rhoUi(this%cfg%imino:this%cfg%imin-1,:,:) = Grho1*vshock; this%fs%rhoVi(this%cfg%imino:this%cfg%imin-1,:,:) = 0.0_WP; this%fs%rhoWi(this%cfg%imino:this%cfg%imin-1,:,:) = 0.0_WP
         this%fs%Grho(this%cfg%imino:this%cfg%imin-1,:,:) = Grho1; this%fs%GP(this%cfg%imino:this%cfg%imin-1,:,:) = GP1
         this%fs%GrhoE(this%cfg%imino:this%cfg%imin-1,:,:) = this%matmod%EOS_energy(GP1,Grho1,vshock,0.0_WP,0.0_WP,'gas')
         
         this%fs%P(this%cfg%imino:this%cfg%imin-1,:,:) = GP1
         this%fs%PA(this%cfg%imino:this%cfg%imin-1,:,:) = GP1
         
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino,this%cfg%imin-1
                  this%fs%GrhoSS2(i,j,k) = this%matmod%EOS_gas(i,j,k,'M')
               end do
            end do
         end do
         
         ! set VOF to zero at inlet 
         this%vf%VF(this%cfg%imino:this%cfg%imin-1,:,:) = 0.0_WP            
      end if
      
      ! outlet at right of domain
      if (this%cfg%iproc.eq.this%cfg%npx) then ! apply (neumann) outlet condition
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imax+1,this%cfg%imaxo
                  this%fs%U(i,j,k) = this%fs%U(this%cfg%imax_,j,k); this%fs%V(i,j,k) = this%fs%V(this%cfg%imax_,j,k); this%fs%W(i,j,k) = this%fs%W(this%cfg%imax_,j,k);
                  this%fs%Ui(i,j,k) = this%fs%Ui(this%cfg%imax_,j,k); this%fs%Vi(i,j,k) = this%fs%Vi(this%cfg%imax_,j,k); this%fs%Wi(i,j,k) = this%fs%Wi(this%cfg%imax_,j,k);
                  this%fs%rhoUi(i,j,k) = this%fs%rhoUi(this%cfg%imax_,j,k); this%fs%rhoVi(i,j,k) = this%fs%rhoVi(this%cfg%imax_,j,k); this%fs%rhoWi(i,j,k) = this%fs%rhoWi(this%cfg%imax_,j,k);
                  this%fs%Grho(i,j,k) = this%fs%Grho(this%cfg%imax_,j,k); this%fs%GP(i,j,k) = this%fs%GP(this%cfg%imax_,j,k)                        
                  this%fs%GrhoE(i,j,k) = this%matmod%EOS_energy(this%fs%GP(this%cfg%imax_,j,k),this%fs%Grho(this%cfg%imax_,j,k),this%fs%U(this%cfg%imax_,j,k),this%fs%V(this%cfg%imax_,j,k),this%fs%W(this%cfg%imax_,j,k),'gas')
                  this%fs%GrhoSS2(i,j,k) = this%matmod%EOS_gas(i,j,k,'M')
                  
                  this%fs%P(i,j,k) = this%fs%P(this%cfg%imax_,j,k)
                  this%fs%PA(i,j,k) = this%fs%PA(this%cfg%imax_,j,k)
                  
                  ! set nuemann condition on VOF
                  this%vf%VF(i,j,k) = this%vf%VF(this%cfg%imax_,j,k)
               end do
            end do
         end do
      end if
      
      ! outlet at bottom of domain
      if (this%cfg%jproc.eq.1) then ! apply (neumann) outlet condition
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino,this%cfg%jmin-1
               do i=this%cfg%imino_,this%cfg%imaxo_
                  this%fs%U(i,j,k) = this%fs%U(i,this%cfg%jmin_,k); this%fs%V(i,j,k) = this%fs%V(i,this%cfg%jmin_,k); this%fs%W(i,j,k) = this%fs%W(i,this%cfg%jmin_,k);
                  this%fs%Ui(i,j,k) = this%fs%Ui(i,this%cfg%jmin_,k); this%fs%Vi(i,j,k) = this%fs%Vi(i,this%cfg%jmin_,k); this%fs%Wi(i,j,k) = this%fs%Wi(i,this%cfg%jmin_,k);
                  this%fs%rhoUi(i,j,k) = this%fs%rhoUi(i,this%cfg%jmin_,k); this%fs%rhoVi(i,j,k) = this%fs%rhoVi(i,this%cfg%jmin_,k); this%fs%rhoWi(i,j,k) = this%fs%rhoWi(i,this%cfg%jmin_,k);
                  this%fs%Grho(i,j,k) = this%fs%Grho(i,this%cfg%jmin_,k); this%fs%GP(i,j,k) = this%fs%GP(i,this%cfg%jmin_,k)
                  this%fs%GrhoE(i,j,k) = this%matmod%EOS_energy(this%fs%GP(i,this%cfg%jmin_,k),this%fs%Grho(i,this%cfg%jmin_,k),this%fs%U(i,this%cfg%jmin_,k),this%fs%V(i,this%cfg%jmin_,k),this%fs%W(i,this%cfg%jmin_,k),'gas')                        
                  this%fs%GrhoSS2(i,j,k) = this%matmod%EOS_gas(i,j,k,'M')
                  
                  this%fs%P(i,j,k) = this%fs%P(i,this%cfg%jmin_,k)
                  this%fs%PA(i,j,k) = this%fs%PA(i,this%cfg%jmin_,k)
                  
                  ! set nuemann condition on VOF
                  this%vf%VF(i,j,k) = this%vf%VF(i,this%cfg%jmin_,k)
                  
               end do
            end do
         end do
      end if
      
      ! outlet at top of domain
      if (this%cfg%jproc.eq.this%cfg%npy) then ! apply (neumann) outlet condition
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmax+1,this%cfg%jmaxo
               do i=this%cfg%imino_,this%cfg%imaxo_
                  this%fs%U(i,j,k) = this%fs%U(i,this%cfg%jmax_,k); this%fs%V(i,j,k) = this%fs%V(i,this%cfg%jmax_,k); this%fs%W(i,j,k) = this%fs%W(i,this%cfg%jmax_,k);
                  this%fs%Ui(i,j,k) = this%fs%Ui(i,this%cfg%jmax_,k); this%fs%Vi(i,j,k) = this%fs%Vi(i,this%cfg%jmax_,k); this%fs%Wi(i,j,k) = this%fs%Wi(i,this%cfg%jmax_,k);
                  this%fs%rhoUi(i,j,k) =this%fs%rhoUi(i,this%cfg%jmax_,k); this%fs%rhoVi(i,j,k) = this%fs%rhoVi(i,this%cfg%jmax_,k); this%fs%rhoWi(i,j,k) = this%fs%rhoWi(i,this%cfg%jmax_,k);
                  this%fs%Grho(i,j,k) = this%fs%Grho(i,this%cfg%jmax_,k); this%fs%GP(i,j,k) = this%fs%GP(i,this%cfg%jmax_,k)
                  this%fs%GrhoE(i,j,k) = this%matmod%EOS_energy(this%fs%GP(i,this%cfg%jmax_,k),this%fs%Grho(i,this%cfg%jmax_,k),this%fs%U(i,this%cfg%jmax_,k),this%fs%V(i,this%cfg%jmax_,k),this%fs%W(i,this%cfg%jmax_,k),'gas')
                  this%fs%GrhoSS2(i,j,k) = this%matmod%EOS_gas(i,j,k,'M')
                  
                  this%fs%P(i,j,k) = this%fs%P(i,this%cfg%jmax_,k)
                  this%fs%PA(i,j,k) = this%fs%PA(i,this%cfg%jmax_,k)
                  
                  ! set nuemann condition on VOF
                  this%vf%VF(i,j,k) = this%vf%VF(i,this%cfg%jmax_,k)
                  
               end do
            end do
         end do
      end if
      
      if (this%cfg%nz.gt.1)then
         ! outlet at back of domain
         if (this%cfg%kproc.eq.1) then ! apply (neumann) outlet condition
            do k=this%cfg%kmino,this%cfg%kmin-1
               do j=this%cfg%jmino_,this%cfg%jmaxo_
                  do i=this%cfg%imino_,this%cfg%imaxo_
                     this%fs%U(i,j,k) = this%fs%U(i,j,this%cfg%kmin_); this%fs%V(i,j,k) = this%fs%V(i,j,this%cfg%kmin_); this%fs%W(i,j,k) = this%fs%W(i,j,this%cfg%kmin_);
                     this%fs%Ui(i,j,k) = this%fs%Ui(i,j,this%cfg%kmin_); this%fs%Vi(i,j,k) = this%fs%Vi(i,j,this%cfg%kmin_); this%fs%Wi(i,j,k) = this%fs%Wi(i,j,this%cfg%kmin_);
                     this%fs%rhoUi(i,j,k) = this%fs%rhoUi(i,j,this%cfg%kmin_); this%fs%rhoVi(i,j,k) = this%fs%rhoVi(i,j,this%cfg%kmin_); this%fs%rhoWi(i,j,k) = this%fs%rhoWi(i,j,this%cfg%kmin_);
                     this%fs%Grho(i,j,k) = this%fs%Grho(i,j,this%cfg%kmin_); this%fs%GP(i,j,k) = this%fs%GP(i,j,this%cfg%kmin_)
                     this%fs%GrhoE(i,j,k) = this%matmod%EOS_energy(this%fs%GP(i,j,this%cfg%kmin_),this%fs%Grho(i,j,this%cfg%kmin_),this%fs%U(i,j,this%cfg%kmin_),this%fs%V(i,j,this%cfg%kmin_),this%fs%W(i,j,this%cfg%kmin_),'gas')
                     this%fs%GrhoSS2(i,j,k) = this%matmod%EOS_gas(i,j,k,'M')
                     
                     this%fs%P(i,j,k) = this%fs%P(i,j,this%cfg%kmin_)
                     this%fs%PA(i,j,k) = this%fs%PA(i,j,this%cfg%kmin_)
                     
                     ! set nuemann condition on VOF
                     this%vf%VF(i,j,k) = this%vf%VF(i,j,this%cfg%kmin_)
                  end do
               end do
            end do
         end if
         
         ! outlet at front of domain
         if (this%cfg%kproc.eq.this%cfg%npz) then ! apply (neumann) outlet condition
            do k=this%cfg%kmax+1,this%cfg%kmaxo
               do j=this%cfg%jmino_,this%cfg%jmaxo_
                  do i=this%cfg%imino_,this%cfg%imaxo_
                     this%fs%U(i,j,k) = this%fs%U(i,j,this%cfg%kmax_); this%fs%V(i,j,k) = this%fs%V(i,j,this%cfg%kmax_); this%fs%W(i,j,k) = this%fs%W(i,j,this%cfg%kmax_);
                     this%fs%Ui(i,j,k) = this%fs%Ui(i,j,this%cfg%kmax_); this%fs%Vi(i,j,k) = this%fs%Vi(i,j,this%cfg%kmax_); this%fs%Wi(i,j,k) = this%fs%Wi(i,j,this%cfg%kmax_);
                     this%fs%rhoUi(i,j,k) = this%fs%rhoUi(i,j,this%cfg%kmax_); this%fs%rhoVi(i,j,k) = this%fs%rhoVi(i,j,this%cfg%kmax_); this%fs%rhoWi(i,j,k) = this%fs%rhoWi(i,j,this%cfg%kmax_);
                     this%fs%Grho(i,j,k) = this%fs%Grho(i,j,this%cfg%kmax_); this%fs%GP(i,j,k) = this%fs%GP(i,j,this%cfg%kmax_)
                     this%fs%GrhoE(i,j,k) = this%matmod%EOS_energy(this%fs%GP(i,j,this%cfg%kmax_),this%fs%Grho(i,j,this%cfg%kmax_),this%fs%U(i,j,this%cfg%kmax_),this%fs%V(i,j,this%cfg%kmax_),this%fs%W(i,j,this%cfg%kmax_),'gas')
                     this%fs%GrhoSS2(i,j,k) = this%matmod%EOS_gas(i,j,k,'M')
                     
                     this%fs%P(i,j,k) = this%fs%P(i,j,this%cfg%kmax_)
                     this%fs%PA(i,j,k) = this%fs%PA(i,j,this%cfg%kmax_)
                     
                     ! set nuemann condition on VOF
                     this%vf%VF(i,j,k) = this%vf%VF(i,j,this%cfg%kmax_)
                  end do
               end do
            end do
         end if
      end if
      
      call this%matmod%update_temperature(this%VF,this%fs%Tmptr)
      
      ! Apply face BC - outflow
      bc_scope = 'velocity'
      call this%fs%apply_bcond(this%time%dt,bc_scope)
    end block restart_IC_BC

    !> create surfmesh object for interface polygon output
    restart_smesh: block
      use irl_fortran_interface
      integer :: i,j,k,nplane,np
      this%smesh=surfmesh(nvar=1,name='plic')
      this%smesh%varname(1)='curv'
      call this%vf%update_surfmesh(this%smesh)
      this%smesh%var(1,:)=0.0_WP
      np=0;
      do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
         do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
            do i=this%vf%cfg%imin_,this%vf%cfg%imax_
               do nplane=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                  if (getNumberOfVertices(this%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                     np=np+1;
                     this%smesh%var(1,np)=this%vf%curv(i,j,k)
                  end if
               end do
            end do
         end do
      end do
    end block restart_smesh
    
    !> Add Ensight output
    restart_ensight: block
      use param,           only: param_read
      ! Create Ensight output from cfg
      this%ens_out=ensight(cfg=this%cfg,name='ShockDroplet')
      ! Create event for Ensight output
      this%ens_evt=event(time=this%time,name='Ensight output')
      call param_read('Ensight output period',this%ens_evt%tper)
      ! Add variables to output
      call this%ens_out%add_vector('velocity',this%fs%Ui,this%fs%Vi,this%fs%Wi)
      call this%ens_out%add_scalar('P',this%fs%P)
      call this%ens_out%add_scalar('PA',this%fs%PA)
      call this%ens_out%add_scalar('Grho',this%fs%Grho)
      call this%ens_out%add_scalar('Lrho',this%fs%Lrho)
      call this%ens_out%add_scalar('Density',this%fs%RHO)
      call this%ens_out%add_scalar('Bulkmod',this%fs%RHOSS2)
      call this%ens_out%add_scalar('VOF',this%vf%VF)
      call this%ens_out%add_scalar('curvature',this%vf%curv)
      call this%ens_out%add_scalar('Mach',this%fs%Mach)
      call this%ens_out%add_scalar('fvf',this%cfg%VF)
      call this%ens_out%add_scalar('Tmptr',this%fs%Tmptr)
      call this%ens_out%add_scalar('SL_x',this%fs%sl_x) 
      call this%ens_out%add_scalar('SL_y',this%fs%sl_y) 
      call this%ens_out%add_scalar('SL_z',this%fs%sl_z) 
      call this%ens_out%add_scalar('LP',this%fs%LP) 
      call this%ens_out%add_scalar('GP',this%fs%GP) 
      call this%ens_out%add_scalar('LrhoE',this%fs%LrhoE) 
      call this%ens_out%add_scalar('GrhoE',this%fs%GrhoE)         
      ! Output to ensight
      !if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
    end block restart_ensight
    
    !> block for writing smesh data
    restart_ensight_smesh: block
      use param,           only: param_read
      real(WP) :: smesh_tper ! declare variable for smesh output frequency
      call param_read('Ensight smesh output period', smesh_tper)
      
      ! create ensight output from cfg for smesh surface reconstruction
      this%ens_out_smesh=ensight(cfg=this%cfg,name='droplet_smesh')
      
      ! create event for ensight output
      this%ens_evt_smesh=event(time=this%time,name='Ensight output smesh')
      this%ens_evt_smesh%tper = smesh_tper
         
      ! add variables to output
      call this%ens_out_smesh%add_surface('smesh',this%smesh)
      
      ! output to ensight
      !if (this%ens_evt_smesh%occurs()) call this%ens_out_smesh%write_data(this%time%t)
    end block restart_ensight_smesh
    
    !> Create a monitor file
    restart_monitor: block
      ! Prepare some info about fields
      call this%fs%get_cfl(this%time%dt,this%time%cfl)
      call this%fs%get_max()
      call this%vf%get_max()
      ! Create simulation monitor
      this%mfile=monitor(this%fs%cfg%amRoot,'simulation')
      call this%mfile%add_column(this%time%n,'Timestep number')
      call this%mfile%add_column(this%time%t,'Time')
      call this%mfile%add_column(this%time%dt,'Timestep size')
      call this%mfile%add_column(this%time%cfl,'Maximum CFL')
      call this%mfile%add_column(this%fs%RHOmin,'RHOmin')
      call this%mfile%add_column(this%fs%RHOmax,'RHOmax')
      call this%mfile%add_column(this%fs%Umax,'Umax')
      call this%mfile%add_column(this%fs%Vmax,'Vmax')
      call this%mfile%add_column(this%fs%Wmax,'Wmax')
      call this%mfile%add_column(this%fs%Pmax,'Pmax')
      call this%mfile%add_column(this%fs%Tmax,'Tmax')
      call this%mfile%write()
      ! Create CFL monitor
      this%cflfile=monitor(this%fs%cfg%amRoot,'cfl')
      call this%cflfile%add_column(this%time%n,'Timestep number')
      call this%cflfile%add_column(this%time%t,'Time')
      call this%cflfile%add_column(this%fs%CFLst,'STension CFL')
      call this%cflfile%add_column(this%fs%CFLc_x,'Convective xCFL')
      call this%cflfile%add_column(this%fs%CFLc_y,'Convective yCFL')
      call this%cflfile%add_column(this%fs%CFLc_z,'Convective zCFL')
      call this%cflfile%add_column(this%fs%CFLv_x,'Viscous xCFL')
      call this%cflfile%add_column(this%fs%CFLv_y,'Viscous yCFL')
      call this%cflfile%add_column(this%fs%CFLv_z,'Viscous zCFL')
      call this%cflfile%add_column(this%fs%CFLa_x,'Acoustic xCFL')
      call this%cflfile%add_column(this%fs%CFLa_y,'Acoustic yCFL')
      call this%cflfile%add_column(this%fs%CFLa_z,'Acoustic zCFL')
      call this%cflfile%write()
      ! Create convergence monitor
      this%cvgfile=monitor(this%fs%cfg%amRoot,'cvg')
      call this%cvgfile%add_column(this%time%n,'Timestep number')
      call this%cvgfile%add_column(this%time%it,'Iteration')
      call this%cvgfile%add_column(this%time%t,'Time')
      call this%cvgfile%add_column(this%fs%impl_it_x,'Impl_x iteration')
      call this%cvgfile%add_column(this%fs%impl_rerr_x,'Impl_x error')
      call this%cvgfile%add_column(this%fs%impl_it_y,'Impl_y iteration')
      call this%cvgfile%add_column(this%fs%impl_rerr_y,'Impl_y error')
      call this%cvgfile%add_column(this%fs%implicit%it,'Impl_z iteration')
      call this%cvgfile%add_column(this%fs%implicit%rerr,'Impl_z error')
      call this%cvgfile%add_column(this%fs%psolv%it,'Pressure iteration')
      call this%cvgfile%add_column(this%fs%psolv%rerr,'Pressure error')
    end block restart_monitor
    
  end subroutine restart
  
  !> Finalize shock droplet simulation
  subroutine final(this)
    implicit none
    class(sdrop), intent(inout) :: this
  end subroutine final
  
  !> functions for localizing domain and setting up levelsets
  !> Function that localizes the left (x-) of the domain
  function left_of_domain(pg,i,j,k) result(isIn)
    use pgrid_class, only: pgrid
    implicit none
    class(pgrid), intent(in) :: pg
    integer, intent(in) :: i,j,k
    logical :: isIn
    isIn=.false.
    if (i.eq.pg%imin) isIn=.true.
  end function left_of_domain
  
  !> Function that localizes the right (x+) of the domain
  function right_of_domain(pg,i,j,k) result(isIn)
    use pgrid_class, only: pgrid
    implicit none
    class(pgrid), intent(in) :: pg
    integer, intent(in) :: i,j,k
    logical :: isIn
    isIn=.false.
    if (i.eq.pg%imax+1) isIn=.true.
  end function right_of_domain
     
  !> Function that localizes the bottom (y-) of the domain
  function bot_of_domain(pg,i,j,k) result(isIn)
    use pgrid_class, only: pgrid
    implicit none
    class(pgrid), intent(in) :: pg
    integer, intent(in) :: i,j,k
    logical :: isIn
    isIn=.false.
    if (j.eq.pg%jmin) isIn=.true.
  end function bot_of_domain
  
  !> Function that localizes the top (y+) of the domain
  function top_of_domain(pg,i,j,k) result(isIn)
    use pgrid_class, only: pgrid
    implicit none
    class(pgrid), intent(in) :: pg
    integer, intent(in) :: i,j,k
    logical :: isIn
    isIn=.false.
    if (j.eq.pg%jmax+1) isIn=.true.
  end function top_of_domain
  
  !> Function that localizes the back (z-) of the domain
  function bck_of_domain(pg,i,j,k) result(isIn)
    use pgrid_class, only: pgrid
    implicit none
    class(pgrid), intent(in) :: pg
    integer, intent(in) :: i,j,k
    logical :: isIn
    isIn=.false.
    if (k.eq.pg%kmin) isIn=.true.
  end function bck_of_domain
  
  !> Function that localizes the front (z+) of the domain
  function fnt_of_domain(pg,i,j,k) result(isIn)
    use pgrid_class, only: pgrid
    implicit none
    class(pgrid), intent(in) :: pg
    integer, intent(in) :: i,j,k
    logical :: isIn
    isIn=.false.
    if (k.eq.pg%kmax+1) isIn=.true.
  end function fnt_of_domain
  
end module shockdrop_class
   
