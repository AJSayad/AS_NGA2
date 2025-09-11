!> Definition for a counterflow class
module counterflow_class
   use precision,         only: WP
   use config_class,      only: config
   use iterator_class,    only: iterator
   use surfmesh_class,    only: surfmesh
   use ensight_class,     only: ensight
   use hypre_str_class,   only: hypre_str
   use ddadi_class,       only: ddadi
   use vfs_class,         only: vfs
   use tpns_class,        only: tpns
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use monitor_class,     only: monitor
   use pardata_class,     only: pardata
   use stracker_class,    only: stracker
   implicit none
   private
   
   public :: counterflow
   
   !> counterflow object
   type :: counterflow
      
      !> Config
      type(config) :: cfg
      
      !> Flow solver
      type(vfs)         :: vf     !< Volume fraction solver
      type(tpns)        :: fs     !< Two-phase flow solver
      type(hypre_str)   :: ps     !< Structured Hypre linear solver for pressure
      type(timetracker) :: time   !< Time info
      
      !> Implicit solver
      logical     :: use_implicit !< Is an implicit solver used?
      type(ddadi) :: vs           !< DDADI solver for velocity
      
      !> Ensight postprocessing
      type(surfmesh) :: smesh     !< Surface mesh for interface
      type(ensight)  :: ens_out   !< Ensight output for flow variables
      type(event)    :: ens_evt   !< Event trigger for Ensight output
      
      !> Simulation monitoring files
      type(monitor) :: mfile      !< General simulation monitoring
      type(monitor) :: cflfile    !< CFL monitoring
      
      !> Work arrays
      real(WP), dimension(:,:,:), allocatable :: resU,resV,resW      !< Residuals
      real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi            !< Cell-centered velocities
      
      !> Iterator for VOF removal
      type(iterator) :: vof_removal_layer  !< Edge of domain where we actively remove VOF
      real(WP)       :: vof_removed        !< Integral of VOF removed
      integer        :: nlayer=4           !< Size of buffer layer for VOF removal
      
      !> Provide a pardata object for restarts
      logical       :: restarted  !< Is the simulation restarted?
      type(pardata) :: df         !< Pardata object for restart I/O
      type(event)   :: save_evt   !< Event to trigger restart I/O
      
   contains
      procedure :: init                            !< Initialize round jet simulation
      procedure :: step                            !< Advance round jet simulation by one time step
      procedure :: final                           !< Finalize round jet simulation
   end type counterflow
   
   
contains
   
   
   !> Initialization of roundjet simulation
   subroutine init(this)
      use param, only: param_read
      implicit none
      class(counterflow), intent(inout) :: this
      
      
      ! Initialize the config
      initialize_config: block
         use sgrid_class, only: sgrid,cartesian
         use parallel,    only: group
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz
         real(WP), dimension(:), allocatable :: x,y,z
         type(sgrid) :: grid
         integer, dimension(3) :: partition
         ! Read in grid definition
         call param_read('Lx',Lx); call param_read('nx',nx); allocate(x(nx+1))
         call param_read('Ly',Ly); call param_read('ny',ny); allocate(y(ny+1))
         call param_read('Lz',Lz); call param_read('nz',nz); allocate(z(nz+1))
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx-0.5_WP*Lx
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='cflow')
         ! Read in partition
         call param_read('Partition',partition)
         ! Create partitioned grid
         this%cfg=config(grp=group,decomp=partition,grid=grid)
         ! No walls in the domain
         this%cfg%VF=1.0_WP
      end block initialize_config
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         this%time=timetracker(amRoot=this%cfg%amRoot)
         call param_read('Max timestep size',this%time%dtmax)
         call param_read('Max cfl number',this%time%cflmax)
         call param_read('Max time',this%time%tmax)
         this%time%dt=this%time%dtmax
         this%time%itmax=2
      end block initialize_timetracker
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(this%resU(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resV(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resW(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Ui  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Vi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Wi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use vfs_class, only: VFlo,VFhi,remap,plicnet,r2pnet
         use mms_geom,  only: cube_refine_vol
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver with plicnet reconstruction
         call this%vf%initialize(cfg=this%cfg,reconstruction_method=r2pnet,transport_method=remap,name='VOF')
         this%vf%thin_thld_min=0.0_WP
         this%vf%flotsam_thld=0.0_WP
         this%vf%maxcurv_times_mesh=1.0_WP
         this%vf%smoothing_maxite=3
         ! Initialize to cylindrical interface
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
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_cylinder,0.0_WP,amr_ref_lvl)
                  this%vf%VF(i,j,k)=vol/this%vf%cfg%vol(i,j,k)
                  if (this%vf%VF(i,j,k).ge.VFlo.and.this%vf%VF(i,j,k).le.VFhi) then
                     this%vf%Lbary(:,i,j,k)=v_cent
                     this%vf%Gbary(:,i,j,k)=([this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]-this%vf%VF(i,j,k)*this%vf%Lbary(:,i,j,k))/(1.0_WP-this%vf%VF(i,j,k))
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
         ! Set simple full-liquid/full-gas interface planes in geometric overlap cells
         call this%vf%set_full_bcond()
         ! Now apply Neumann condition on interface at inlets to have proper round injection
         neumann_irl: block
            use irl_fortran_interface, only: getPlane,new,construct_2pt,RectCub_type,&
            &                                setNumberOfPlanes,setPlane,matchVolumeFraction
            real(WP), dimension(1:4) :: plane
            type(RectCub_type) :: cell
            call new(cell)
            if (this%vf%cfg%iproc.eq.1) then
               do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
                  do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
                     do i=this%vf%cfg%imino,this%vf%cfg%imin-1
                        ! Extract plane data and copy in overlap
                        plane=getPlane(this%vf%liquid_gas_interface(this%vf%cfg%imin,j,k),0)
                        call construct_2pt(cell,[this%vf%cfg%x(i  ),this%vf%cfg%y(j  ),this%vf%cfg%z(k  )],&
                        &                       [this%vf%cfg%x(i+1),this%vf%cfg%y(j+1),this%vf%cfg%z(k+1)])
                        plane(4)=dot_product(plane(1:3),[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)])
                        call setNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k),1)
                        call setPlane(this%vf%liquid_gas_interface(i,j,k),0,plane(1:3),plane(4))
                        call matchVolumeFraction(cell,this%vf%VF(i,j,k),this%vf%liquid_gas_interface(i,j,k))
                     end do
                  end do
               end do
            end if
            if (this%vf%cfg%iproc.eq.this%vf%cfg%npx) then
               do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
                  do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
                     do i=this%vf%cfg%imax+1,this%vf%cfg%imaxo
                        ! Extract plane data and copy in overlap
                        plane=getPlane(this%vf%liquid_gas_interface(this%vf%cfg%imax,j,k),0)
                        call construct_2pt(cell,[this%vf%cfg%x(i  ),this%vf%cfg%y(j  ),this%vf%cfg%z(k  )],&
                        &                       [this%vf%cfg%x(i+1),this%vf%cfg%y(j+1),this%vf%cfg%z(k+1)])
                        plane(4)=dot_product(plane(1:3),[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)])
                        call setNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k),1)
                        call setPlane(this%vf%liquid_gas_interface(i,j,k),0,plane(1:3),plane(4))
                        call matchVolumeFraction(cell,this%vf%VF(i,j,k),this%vf%liquid_gas_interface(i,j,k))
                     end do
                  end do
               end do
            end if
         end block neumann_irl
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
      end block create_and_initialize_vof
      
      
      ! Create an iterator for removing VOF at edges
      create_iterator: block
         this%vof_removal_layer=iterator(this%cfg,'VOF removal',vof_removal_layer_locator)
         this%vof_removed=0.0_WP
      end block create_iterator
      
      
      ! Create an incompressible flow solver with bconds
      create_flow_solver: block
         use hypre_str_class, only: pcg_pfmg2
         use tpns_class,      only: dirichlet,clipped_neumann
         real(WP) :: Froude
         ! Create flow solver
         this%fs=tpns(cfg=this%cfg,name='Two-phase NS')
         ! Set fluid properties
         this%fs%rho_g=1.0_WP; call param_read('Density ratio',this%fs%rho_l)
         call param_read('Reynolds number',this%fs%visc_l); this%fs%visc_l=this%fs%rho_l/this%fs%visc_l
         call param_read('Viscosity ratio',this%fs%visc_g); this%fs%visc_g=this%fs%visc_l/this%fs%visc_g
         call param_read('Weber number',this%fs%sigma); this%fs%sigma=this%fs%rho_l/this%fs%sigma
         call param_read('Froude number',Froude,default=0.0_WP); if (Froude.gt.0.0_WP) this%fs%gravity=[Froude**(-2.0_WP),0.0_WP,0.0_WP]
         ! Inflows in xm/xp and outflows in ym/yp/zm/zp
         call this%fs%add_bcond(name='inflow_xm',type=dirichlet,face='x',dir=-1,canCorrect=.false.,locator=xm_locator)
         call this%fs%add_bcond(name='inflow_xp',type=dirichlet,face='x',dir=+1,canCorrect=.false.,locator=xp_locator)
         call this%fs%add_bcond(name='outlet_ym',type=clipped_neumann,face='y',dir=-1,canCorrect=.true.,locator=ym_locator)
         call this%fs%add_bcond(name='outlet_yp',type=clipped_neumann,face='y',dir=+1,canCorrect=.true.,locator=yp_locator)
         call this%fs%add_bcond(name='outlet_zm',type=clipped_neumann,face='z',dir=-1,canCorrect=.true.,locator=zm_locator)
         call this%fs%add_bcond(name='outlet_zp',type=clipped_neumann,face='z',dir=+1,canCorrect=.true.,locator=zp_locator)
         ! Configure pressure solver
         this%ps=hypre_str(cfg=this%cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         this%ps%maxlevel=16
         call param_read('Pressure iteration',this%ps%maxit)
         call param_read('Pressure tolerance',this%ps%rcvg)
         ! Check if we want to use an implicit solver
         call param_read('Use implicit solver',this%use_implicit)
         if (this%use_implicit) then
            ! Configure implicit velocity solver
            this%vs=ddadi(cfg=this%cfg,name='Velocity',nst=7)
            ! Setup the solver
            call this%fs%setup(pressure_solver=this%ps,implicit_solver=this%vs)
         else
            ! Setup the solver
            call this%fs%setup(pressure_solver=this%ps)
         end if
      end block create_flow_solver
      
      
      ! Initialize our velocity field
      initialize_velocity: block
         use tpns_class, only: bcond
         type(bcond), pointer :: mybc
         integer :: n,i,j,k
         ! Zero initial field
         this%fs%U=0.0_WP; this%fs%V=0.0_WP; this%fs%W=0.0_WP
         ! Apply inflow condition on velocity
         call this%fs%get_bcond('inflow_xm',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            if (sqrt(this%cfg%ym(j)**2+this%cfg%zm(k)**2).le.0.5_WP) this%fs%U(i,j,k)=+1.0_WP
         end do
         call this%fs%get_bcond('inflow_xp',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            if (sqrt(this%cfg%ym(j)**2+this%cfg%zm(k)**2).le.0.5_WP) this%fs%U(i,j,k)=-1.0_WP
         end do
         ! Apply all other boundary conditions
         call this%fs%apply_bcond(this%time%t,this%time%dt)
         ! Compute MFR through all boundary conditions
         call this%fs%get_mfr()
         ! Adjust MFR for global mass balance
         call this%fs%correct_mfr()
         ! Compute cell-centered velocity
         call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
         ! Compute divergence
         call this%fs%get_div()
      end block initialize_velocity
      
      
      ! Handle restart/saves here
      handle_restart: block
         use string,                only: str_medium
         use filesys,               only: makedir,isdir
         use irl_fortran_interface, only: setNumberOfPlanes,setPlane
         character(len=str_medium) :: timestamp
         integer, dimension(3) :: iopartition
         real(WP), dimension(:,:,:), allocatable :: P11,P12,P13,P14
         real(WP), dimension(:,:,:), allocatable :: P21,P22,P23,P24
         integer :: i,j,k
         ! Create event for saving restart files
         this%save_evt=event(this%time,'Jet restart output')
         call param_read('Restart output period',this%save_evt%tper)
         ! Read in the I/O partition
         call param_read('I/O partition',iopartition)
         ! Check if we are restarting
         call param_read('Restart from',timestamp,default='')
         this%restarted=.false.; if (len_trim(timestamp).gt.0) this%restarted=.true.
         ! Perform pardata initialization
         if (this%restarted) then
            ! Read in the file
            call this%df%initialize(pg=this%cfg,iopartition=iopartition,fdata='restart/data_'//trim(timestamp))
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
            call this%vf%sync_interface()
            deallocate(P11,P12,P13,P14,P21,P22,P23,P24)
            ! Reset moments
            call this%vf%reset_volume_moments()
            ! Update the band
            call this%vf%update_band()
            ! Create discontinuous polygon mesh from IRL interface
            call this%vf%polygonalize_interface()
            ! Calculate distance from polygons
            call this%vf%distance_from_polygon()
            ! Calculate subcell phasic volumes
            call this%vf%subcell_vol()
            ! Calculate curvature
            call this%vf%get_curvature()
            ! Now read in the velocity solver data
            call this%df%pull(name='U',var=this%fs%U)
            call this%df%pull(name='V',var=this%fs%V)
            call this%df%pull(name='W',var=this%fs%W)
            call this%df%pull(name='P',var=this%fs%P)
            call this%df%pull(name='Pjx',var=this%fs%Pjx)
            call this%df%pull(name='Pjy',var=this%fs%Pjy)
            call this%df%pull(name='Pjz',var=this%fs%Pjz)
            ! Apply all other boundary conditions
            call this%fs%apply_bcond(this%time%t,this%time%dt)
            ! Compute MFR through all boundary conditions
            call this%fs%get_mfr()
            ! Adjust MFR for global mass balance
            call this%fs%correct_mfr()
            ! Compute cell-centered velocity
            call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
            ! Compute divergence
            call this%fs%get_div()
            ! Also update time
            call this%df%pull(name='t' ,val=this%time%t )
            call this%df%pull(name='dt',val=this%time%dt)
            this%time%told=this%time%t-this%time%dt
            !this%time%dt=this%time%dtmax !< Force max timestep size anyway
         else
            ! We are not restarting, prepare a new directory for storing restart files
            if (this%cfg%amRoot) then
               if (.not.isdir('restart')) call makedir('restart')
            end if
            ! Prepare pardata object for saving restart files
            call this%df%initialize(pg=this%cfg,iopartition=iopartition,filename=trim(this%cfg%name),nval=2,nvar=15)
            this%df%valname=['t ','dt']
            this%df%varname=['U  ','V  ','W  ','P  ','Pjx','Pjy','Pjz','P11','P12','P13','P14','P21','P22','P23','P24']
         end if
      end block handle_restart
      
      
      ! Create surfmesh object for interface polygon output
      create_smesh: block
         use irl_fortran_interface, only: getNumberOfPlanes,getNumberOfVertices
         integer :: i,j,k,np,nplane
         this%smesh=surfmesh(nvar=3,name='plic')
         this%smesh%varname(1)='nplane'
         this%smesh%varname(2)='thickness'
         this%smesh%varname(3)='curvature'
         ! Transfer polygons to smesh
         call this%vf%update_surfmesh(this%smesh)
         ! Calculate curv2p and thickness even for plic
         if (.not.this%vf%two_planes) then
            allocate(this%vf%curv2p(1:2,this%vf%cfg%imino_:this%vf%cfg%imaxo_,this%vf%cfg%jmino_:this%vf%cfg%jmaxo_,this%vf%cfg%kmino_:this%vf%cfg%kmaxo_)); this%vf%curv2p   =0.0_WP
            allocate(this%vf%thickness (this%vf%cfg%imino_:this%vf%cfg%imaxo_,this%vf%cfg%jmino_:this%vf%cfg%jmaxo_,this%vf%cfg%kmino_:this%vf%cfg%kmaxo_)); this%vf%thickness=0.0_WP
            this%vf%curv2p(1,:,:,:)=this%vf%curv(:,:,:)
            call this%vf%get_thickness()
         end if
         ! Populate surface variables
         np=0
         do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
            do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
               do i=this%vf%cfg%imin_,this%vf%cfg%imax_
                  do nplane=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                     if (getNumberOfVertices(this%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                        np=np+1; this%smesh%var(1,np)=real(getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k)),WP)
                        this%smesh%var(2,np)=this%vf%thickness(i,j,k)
                        this%smesh%var(3,np)=this%vf%curv2p(nplane,i,j,k)
                     end if
                  end do
               end do
            end do
         end do
      end block create_smesh
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         this%ens_out=ensight(cfg=this%cfg,name='cflow')
         ! Create event for Ensight output
         this%ens_evt=event(time=this%time,name='Ensight output')
         call param_read('Ensight output period',this%ens_evt%tper)
         ! Add variables to output
         call this%ens_out%add_vector('velocity',this%Ui,this%Vi,this%Wi)
         call this%ens_out%add_scalar('pressure',this%fs%P)
         call this%ens_out%add_scalar('VOF',this%vf%VF)
         call this%ens_out%add_scalar('curvature',this%vf%curv)
         call this%ens_out%add_surface('plic',this%smesh)
         ! Output to ensight
         if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
      end block create_ensight
      
      
      ! Create a monitor file
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
         call this%mfile%add_column(this%fs%Umax,'Umax')
         call this%mfile%add_column(this%fs%Vmax,'Vmax')
         call this%mfile%add_column(this%fs%Wmax,'Wmax')
         call this%mfile%add_column(this%fs%Pmax,'Pmax')
         call this%mfile%add_column(this%vf%VFmax,'VOF maximum')
         call this%mfile%add_column(this%vf%VFmin,'VOF minimum')
         call this%mfile%add_column(this%vf%VFint,'VOF integral')
         call this%mfile%add_column(this%vof_removed,'VOF removed')
         call this%mfile%add_column(this%vf%SDint,'SD integral')
         call this%mfile%add_column(this%fs%divmax,'Maximum divergence')
         call this%mfile%add_column(this%fs%psolv%it,'Pressure iteration')
         call this%mfile%add_column(this%fs%psolv%rerr,'Pressure error')
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
         call this%cflfile%write()
      end block create_monitor
      
      
   contains
      
      
      !> Function that defines a level set function for a cylinder
      function levelset_cylinder(xyz,t) result(G)
         implicit none
         real(WP), dimension(3),intent(in) :: xyz
         real(WP), intent(in) :: t
         real(WP) :: G
         G=0.5_WP-sqrt(xyz(2)**2+xyz(3)**2)
      end function levelset_cylinder
      
      
      !> Function that localizes region of VOF removal
      function vof_removal_layer_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (j.le.pg%jmin+this%nlayer.or.&
         &   j.ge.pg%jmax-this%nlayer.or.&
         &   k.le.pg%kmin+this%nlayer.or.&
         &   k.ge.pg%kmax-this%nlayer) isIn=.true.
      end function vof_removal_layer_locator
      
      
      !> Function that localizes the right (x+) of the domain
      function xp_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (i.eq.pg%imax+1) isIn=.true.
      end function xp_locator
      
      
      !> Function that localizes the left (x-) of the domain
      function xm_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         implicit none
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (i.eq.pg%imin) isIn=.true.
      end function xm_locator
      
      
      !> Function that localizes the top (y+) of the domain
      function yp_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         implicit none
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (j.eq.pg%jmax+1) isIn=.true.
      end function yp_locator
      
      
      !> Function that localizes the bottom (y-) of the domain
      function ym_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         implicit none
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (j.eq.pg%jmin) isIn=.true.
      end function ym_locator
      
      
      !> Function that localizes the front (z+) of the domain
      function zp_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         implicit none
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (k.eq.pg%kmax+1) isIn=.true.
      end function zp_locator
      
      
      !> Function that localizes the back (z-) of the domain
      function zm_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         implicit none
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (k.eq.pg%kmin) isIn=.true.
      end function zm_locator
      
      
   end subroutine init
   
   
   !> Take one time step
   subroutine step(this)
      use tpns_class, only: arithmetic_visc
      implicit none
      class(counterflow), intent(inout) :: this
      
      ! Increment time
      call this%fs%get_cfl(this%time%dt,this%time%cfl)
      call this%time%adjust_dt()
      call this%time%increment()
      
      ! Remember old VOF
      this%vf%VFold=this%vf%VF
      
      ! Remember old velocity
      this%fs%Uold=this%fs%U
      this%fs%Vold=this%fs%V
      this%fs%Wold=this%fs%W
      
      ! Prepare old staggered density (at n)
      call this%fs%get_olddensity(vf=this%vf)
      
      ! VOF solver step
      call this%vf%advance(dt=this%time%dt,U=this%fs%U,V=this%fs%V,W=this%fs%W)
      
      ! Prepare new staggered viscosity (at n+1)
      call this%fs%get_viscosity(vf=this%vf,strat=arithmetic_visc)
      
      ! Perform sub-iterations
      do while (this%time%it.le.this%time%itmax)
         
         ! Build mid-time velocity
         this%fs%U=0.5_WP*(this%fs%U+this%fs%Uold)
         this%fs%V=0.5_WP*(this%fs%V+this%fs%Vold)
         this%fs%W=0.5_WP*(this%fs%W+this%fs%Wold)
         
         ! Preliminary mass and momentum transport step at the interface
         call this%fs%prepare_advection_upwind(dt=this%time%dt)
         
         ! Explicit calculation of drho*u/dt from NS
         call this%fs%get_dmomdt(this%resU,this%resV,this%resW)
         
         ! Add momentum source terms
         call this%fs%addsrc_gravity(this%resU,this%resV,this%resW)
         
         ! Assemble explicit residual
         this%resU=-2.0_WP*this%fs%rho_U*this%fs%U+(this%fs%rho_Uold+this%fs%rho_U)*this%fs%Uold+this%time%dt*this%resU
         this%resV=-2.0_WP*this%fs%rho_V*this%fs%V+(this%fs%rho_Vold+this%fs%rho_V)*this%fs%Vold+this%time%dt*this%resV
         this%resW=-2.0_WP*this%fs%rho_W*this%fs%W+(this%fs%rho_Wold+this%fs%rho_W)*this%fs%Wold+this%time%dt*this%resW   
         
         ! Form implicit residuals
         call this%fs%solve_implicit(this%time%dt,this%resU,this%resV,this%resW)
         
         ! Apply these residuals
         this%fs%U=2.0_WP*this%fs%U-this%fs%Uold+this%resU
         this%fs%V=2.0_WP*this%fs%V-this%fs%Vold+this%resV
         this%fs%W=2.0_WP*this%fs%W-this%fs%Wold+this%resW
         
         ! Apply other boundary conditions on the resulting fields
         call this%fs%apply_bcond(this%time%t,this%time%dt)
         
         ! Solve Poisson equation
         call this%fs%update_laplacian()
         call this%fs%correct_mfr()
         call this%fs%get_div()
         if (this%vf%two_planes) then
            call this%fs%add_surface_tension_jump_twoVF(dt=this%time%dt,div=this%fs%div,vf=this%vf)
         else
            call this%fs%add_surface_tension_jump(dt=this%time%dt,div=this%fs%div,vf=this%vf)
         end if
         this%fs%psolv%rhs=-this%fs%cfg%vol*this%fs%div/this%time%dt
         this%fs%psolv%sol=0.0_WP
         call this%fs%psolv%solve()
         call this%fs%shift_p(this%fs%psolv%sol)
         
         ! Correct velocity
         call this%fs%get_pgrad(this%fs%psolv%sol,this%resU,this%resV,this%resW)
         this%fs%P=this%fs%P+this%fs%psolv%sol
         this%fs%U=this%fs%U-this%time%dt*this%resU/this%fs%rho_U
         this%fs%V=this%fs%V-this%time%dt*this%resV/this%fs%rho_V
         this%fs%W=this%fs%W-this%time%dt*this%resW/this%fs%rho_W
         
         ! Increment sub-iteration counter
         this%time%it=this%time%it+1
         
      end do
      
      ! Recompute interpolated velocity and divergence
      call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
      call this%fs%get_div()
      
      ! Remove VOF at edge of domain
      remove_vof: block
         use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
         use parallel, only: MPI_REAL_WP
         integer :: n,i,j,k,ierr
         real(WP) :: my_vof_removed
         this%vof_removed=0.0_WP
         do n=1,this%vof_removal_layer%no_
            i=this%vof_removal_layer%map(1,n)
            j=this%vof_removal_layer%map(2,n)
            k=this%vof_removal_layer%map(3,n)
            if (n.le.this%vof_removal_layer%n_) this%vof_removed=this%vof_removed+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
            this%vf%VF(i,j,k)=0.0_WP
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%vof_removed,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
         call this%vf%clean_irl_and_band()
      end block remove_vof
      
      ! Output to ensight
      if (this%ens_evt%occurs()) then
         update_smesh: block
            use irl_fortran_interface, only: getNumberOfPlanes,getNumberOfVertices
            integer :: i,j,k,np,nplane
            ! Transfer polygons to smesh
            call this%vf%update_surfmesh(this%smesh)
            ! Calculate thickness and curv2p even for plic
            if (.not.this%vf%two_planes) then
               this%vf%curv2p(1,:,:,:)=this%vf%curv(:,:,:)
               call this%vf%get_thickness()
            end if
            ! Populate surface variables
            np=0
            do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
               do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
                  do i=this%vf%cfg%imin_,this%vf%cfg%imax_
                     do nplane=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                        if (getNumberOfVertices(this%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                           np=np+1; this%smesh%var(1,np)=real(getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k)),WP)
                           this%smesh%var(2,np)=this%vf%thickness(i,j,k)
                           this%smesh%var(3,np)=this%vf%curv2p(nplane,i,j,k)
                        end if
                     end do
                  end do
               end do
            end do
         end block update_smesh
         call this%ens_out%write_data(this%time%t)
      end if
      
      ! Perform and output monitoring
      call this%fs%get_max()
      call this%vf%get_max()
      call this%mfile%write()
      call this%cflfile%write()
      
      ! Finally, see if it's time to save restart files
      if (this%save_evt%occurs()) then
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
            call this%df%push(name='t'  ,val=this%time%t )
            call this%df%push(name='dt' ,val=this%time%dt)
            call this%df%push(name='U'  ,var=this%fs%U   )
            call this%df%push(name='V'  ,var=this%fs%V   )
            call this%df%push(name='W'  ,var=this%fs%W   )
            call this%df%push(name='P'  ,var=this%fs%P   )
            call this%df%push(name='Pjx',var=this%fs%Pjx )
            call this%df%push(name='Pjy',var=this%fs%Pjy )
            call this%df%push(name='Pjz',var=this%fs%Pjz )
            call this%df%push(name='P11',var=P11         )
            call this%df%push(name='P12',var=P12         )
            call this%df%push(name='P13',var=P13         )
            call this%df%push(name='P14',var=P14         )
            call this%df%push(name='P21',var=P21         )
            call this%df%push(name='P22',var=P22         )
            call this%df%push(name='P23',var=P23         )
            call this%df%push(name='P24',var=P24         )
            call this%df%write(fdata='restart/data_'//trim(adjustl(timestamp)))
            ! Deallocate
            deallocate(P11,P12,P13,P14,P21,P22,P23,P24)
         end block save_restart
      end if
      
   end subroutine step
   

   !> Finalize nozzle simulation
   subroutine final(this)
      implicit none
      class(counterflow), intent(inout) :: this
      
      ! Deallocate work arrays
      deallocate(this%resU,this%resV,this%resW,this%Ui,this%Vi,this%Wi)
      
   end subroutine final
   
   
end module counterflow_class