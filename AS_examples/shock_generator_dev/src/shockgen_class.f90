!> Definition for a shockgen class (shock generator)
module shockgen_class
  use precision,         only: WP
  use config_class,      only: config
  use mast_class,        only: mast
  use vfs_class,         only: vfs
  use matm_class,        only: matm
  use timetracker_class, only: timetracker
  use ensight_class,     only: ensight
  use event_class,       only: event
  use monitor_class,     only: monitor
  use hypre_str_class,   only: hypre_str
  use pardata_class,     only: pardata
  implicit none
  private

  public :: sgen

  !> sgen object
  type :: sgen
     !> Config
     type(config)       :: cfg !> Mesh for solver
     !> Flow solver
     type(mast)         :: fs
     type(vfs)          :: vf
     type(matm)         :: matmod
     type(timetracker)  :: time
     type(hypre_str)    :: ps
     type(hypre_str)    :: vs
     !> Ensight postprocessing
     type(ensight)      :: ens_out
     type(event)        :: ens_evt
     !> Simulation monitor file
     type(monitor)      :: mfile,cflfile,cvgfile
     !> Fluid parameters
     real(WP)           :: visc !AS is this needed? I borrowed this from Chase's shear layer pre-sim sml_class.f90
     real(WP), dimension(:),  allocatable :: Grho_profile, GrhoE_profile, Ui_profile, GP_profile
     integer            :: relax_model
   contains
     procedure :: init  !> initialize sgen simulation
     procedure :: step  !> advance sgen simulation by one timestep
     procedure :: final !> finalize sgen simulation
  end type sgen

contains

  !> Initialization of the shock generator (sgen) simulation
  subroutine init(this)
    implicit none
    class(sgen), intent(inout) :: this

    ! Create mesh for sgen (we do this twice, once here and once in shockdrop_class.f90
    ! becuase the sgen mesh is a reduced size since we only need the shock profile data in x direction)
    create_config: block
      use sgrid_class, only: cartesian,sgrid
      use param,       only: param_read, param_exists
      use parallel,    only: amRoot,group
      use messager,    only: die
      type(sgrid) :: grid
      integer, dimension(3) :: partition
      integer  :: i,j,k,nx,ny,nz
      real(WP) :: Lx,dx,Ly,dy,Lz,dz,alpha
      real(WP), dimension(:), allocatable :: x,y,z

      ! variables for stretching
      integer  :: nx_stretchL,nx_stretchR,ny_stretch,nz_stretch
      real(WP) :: dx_old,dx_ref,start_ref,dy_old,dy_stretch,dz_old,dz_stretch

      ! stretching ratio
      alpha = 1.03_WP

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

    end block create_config
    
    initialize_timetracker: block
      use param, only: param_read

      real(WP) :: start_xshock,final_xshock,vshock,relshockvel ! how far the shock will travel
      real(WP) :: gamm_g,Ma,Grho0,GP0,Grho1,GP1,Ma1            ! gas properties
      
      this%time=timetracker(amRoot=this%cfg%amRoot)
      call param_read('Shockgen Max timestep size',this%time%dtmax)
      call param_read('Max cfl number',this%time%cflmax)

      ! use shock values to calculate final simulation time
      call param_read('Single phase shock location',start_xshock)
      call param_read('Gas gamma',gamm_g)
      call param_read('Final shock location',final_xshock) !final singlephase shock location
      call param_read('Pre-shock density',Grho0,default=1.204_WP)
      call param_read('Pre-shock pressure',GP0,default=1.01325e5_WP)
      call param_read('Mach number of shock',Ma,default=1.47_WP)
      
      !use shock relations to get post shock numbers
      GP1 = GP0 * (2.0_WP*gamm_g*Ma**2 - (gamm_g-1.0_WP)) / (gamm_g+1.0_WP)
      Grho1 = Grho0 * (Ma**2 * (gamm_g+1.0_WP) / ((gamm_g-1.0_WP)*Ma**2 + 2.0_WP))
      Ma1 = sqrt(((gamm_g-1.0_WP)*(Ma**2)+2.0_WP)/(2.0_WP*gamm_g*(Ma**2)-(gamm_g-1.0_WP)))
      vshock = -Ma1 * sqrt(gamm_g*GP1/Grho1) + Ma*sqrt(gamm_g*GP0/Grho0)
      relshockvel = -Grho1*vshock/(Grho0-Grho1)

      ! calculate final shock time based on shock speed
      this%time%tmax = (final_xshock - start_xshock) / relshockvel

      call param_read('Max steps',this%time%nmax)
      this%time%dt=this%time%dtmax
      this%time%itmax=2
    end block initialize_timetracker

    !> create VOF solver
    create_VOF_solver: block
      ! even though this is singlephase, MAST depends on vfs_class, so we still need to initailze it
      use mms_geom, only: cube_refine_vol
      use vfs_class, only: VFhi,VFlo,plicnet,flux,neumann
      use irl_fortran_interface

      integer :: i,j,k,n,si,sj,sk
      real(WP), dimension(3,8) :: cube_vertex
      real(WP), dimension(3) :: v_cent,a_cent
      real(WP) :: vol,area
      integer, parameter :: amr_ref_lvl=4
  
      ! set VOF to zero everywhere
      this%vf%VF = 0.0_WP
      
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

      ! AS do we need this for a singlephase sim?
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
      use hypre_str_class, only: pcg_pfmg ! preconditioned conjugate gradient method for pressure and velocity
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
      real(WP) :: xshock,vshock,relshockvel,Lx,start_xshock
      real(WP) :: Grho0, GP0, Grho1, GP1, ST, Ma1, Ma, Lrho0, LP0, Mas
      type(bcond), pointer :: mybc
      
      ! variables for shock generation
      integer  :: n_shock,shock_index
      real(WP) :: final_xshock,delta,dx,tol,shock_loc
      
      ! set up for shock profile
      call param_read('n_shock',n_shock) ! number of points to capture shock profile
      call param_read('Lx',Lx); call param_read('nx',nx)
      dx = Lx/nx ! mesh spacing in uniform region
      tol = dx/2 ! set tolerance for reading in shock profile
      delta = 2*dx*n_shock ! shock thickness
      
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

      ! singlephase, set liquid properties to unity, VOF has been set to zero in VOF solver setup already
      this%fs%Lrho = 1.0_WP; this%fs%LP = 1.0_WP; this%fs%LrhoE = 1.0_WP

      ! get gas properties from input file
      call param_read('Pre-shock density',Grho0,default=1.204_WP)
      call param_read('Pre-shock pressure',GP0,default=1.01325e5_WP)
      call param_read('Mach number of shock',Ma,default=1.47_WP)
      call param_read('Single phase shock location',start_xshock) 
      call param_read('Final shock location',final_xshock) 

      !use shock relations to get post shock numbers
      GP1 = GP0 * (2.0_WP*gamm_g*Ma**2 - (gamm_g-1.0_WP)) / (gamm_g+1.0_WP)
      Grho1 = Grho0 * (Ma**2 * (gamm_g+1.0_WP) / ((gamm_g-1.0_WP)*Ma**2 + 2.0_WP))
      !calculate post shock Mach number (mach number of gas behind shock)
      Ma1 = sqrt(((gamm_g-1.0_WP)*(Ma**2)+2.0_WP)/(2.0_WP*gamm_g*(Ma**2)-(gamm_g-1.0_WP)))
      !calculate post shock velocity (velocity of the gas behind the shock)
      vshock = -Ma1 * sqrt(gamm_g*GP1/Grho1) + Ma*sqrt(gamm_g*GP0/Grho0)
      !velocity at which the shock moves
      relshockvel = -Grho1*vshock/(Grho0-Grho1)

      if (amRoot) then
         print*, "===== Shockgen Problem Setup Description ====="
         print*, 'Mach number', Ma
         print*, 'Pre-shock:  Density',Grho0,'Pressure',GP0
         print*, 'Post-shock: Density',Grho1,'Pressure',GP1,'Gas Velocity',vshock
         print*, 'Shock velocity', relshockvel
         print*, "Total shock profile points: ", 2*n_shock
         print*, "Shock thickness: ", delta
         print*, "Tolerance for finding shock center: ", tol
         print*, "=============================================="
      end if

      ! zero out cell centered velocities in y and z 
      this%fs%Vi = 0.0_WP; this%fs%Wi = 0.0_WP  
      ! Zero out face velocities as well for the sake of dirichlet boundaries
      this%fs%V = 0.0_WP; this%fs%W = 0.0_WP

      ! set our initial conditions (shock discontinuity)
      do i=this%fs%cfg%imino_,this%fs%cfg%imaxo_
         if (this%cfg%xm(i).lt.start_xshock) then ! post shock values
            this%fs%Grho(i,:,:) = Grho1 
            this%fs%Ui(i,:,:) = vshock
            this%fs%GP(i,:,:) = GP1
            this%fs%GrhoE(i,:,:) = this%matmod%EOS_energy(GP1,Grho1,vshock,0.0_WP,0.0_WP,'gas')
         else ! pre shock values
            this%fs%Grho(i,:,:) = Grho0
            this%fs%Ui(i,:,:) = 0.0_WP
            this%fs%GP(i,:,:) = GP0
            this%fs%GrhoE(i,:,:) = this%matmod%EOS_energy(GP0,Grho0,0.0_WP,0.0_WP,0.0_WP,'gas')
         end if
      end do

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
      
      ! choose pressure relax model
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
      
    !> singlephase, no need for smesh
    !> create ensight output
    create_ensight: block
      use param,           only: param_read
      ! Create Ensight output from cfg
      this%ens_out=ensight(cfg=this%cfg,name='Shockgen')
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
      if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
    end block create_ensight

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

  !> Take one time step with specified dt
  subroutine step(this)
    use messager, only: die
    implicit none
    class(sgen), intent(inout) :: this

    ! Perform time integration in simulation.f90 file 
    ! Increment time
    call this%fs%get_cfl(this%time%dt,this%time%cfl)
    call this%time%adjust_dt()
    call this%time%increment()
    
    ! Reinitialize phase pressure by syncing it with conserved phase energy
    call this%fs%reinit_phase_pressure(this%vf,this%matmod)
    ! remember old velocity and density
    this%fs%Uiold=this%fs%Ui; this%fs%Viold=this%fs%Vi; this%fs%Wiold=this%fs%Wi;this%fs%RHOold = this%fs%RHO
    
    ! AS do we need these steps for singlephase?
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
    
    ! Perform and output monitoring
    call this%fs%get_max()
    call this%vf%get_max()
    call this%fs%get_viz()
    call this%mfile%write()
    call this%cflfile%write()
    
  end subroutine step
    
  !> Finalize shock generator (sgen) simulation
  subroutine final(this)
    use param, only: param_read
    use parallel,   only: amRoot
    implicit none
    class (sgen), intent(inout) :: this    
    integer :: i,j,shock_index,n_shock, nx
    real(WP) :: final_xshock, delta, start_ref,Lx
    real(WP) :: tol ! tolerance for finding final shock location in singlephase
    real(WP), dimension(:),  allocatable :: Grho_profile, GrhoE_profile, Ui_profile, GP_profile
    ! allocate shock profile arrays
    allocate(Grho_profile(2*n_shock+1));allocate(GrhoE_profile(2*n_shock+1));allocate(Ui_profile(2*n_shock+1));allocate(GP_profile(2*n_shock+1))
    
    call param_read('nx',nx)
    call param_read('Lx',Lx);  call param_read('Lx ref',start_ref);
    call param_read('n_shock',n_shock)
    call param_read('Final shock location',final_xshock)
    call param_read('Lx ref', start_ref, default=0.0_WP);
    
    delta = 2*this%cfg%dx(1)*n_shock !shock thickness
    tol = (Lx - start_ref)/nx ! set the tolerance to the mesh spacing in the uniform region
    
    do i=this%cfg%imino_,this%cfg%imaxo_
       if ((this%cfg%xm(i).lt.(final_xshock+tol)).and.(this%cfg%xm(i).gt.(final_xshock-tol))) then
          print*, "The shock has been found at index: ", i
          shock_index=i
          print*, "stored shock index", shock_index
       end if
    end do
    
    do i=shock_index-n_shock,shock_index+n_shock ! save 2*n_shock points for shock profile
       Grho_profile(i+n_shock-shock_index+1)=this%fs%Grho(i,1,1)
       GrhoE_profile(i+n_shock-shock_index+1)=this%fs%GrhoE(i,1,1)
       GP_profile(i+n_shock-shock_index+1)=this%fs%GP(i,1,1)
       Ui_profile(i+n_shock-shock_index+1)=this%fs%Ui(i,1,1)
    end do

  end subroutine final
  
  !> functions for localizing domain
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
  
end module shockgen_class
