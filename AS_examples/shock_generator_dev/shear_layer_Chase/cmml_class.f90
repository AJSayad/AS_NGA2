!> Definition for an ml class
module ml_class
   use precision,         only: WP
   use config_class,      only: config
   use mast_class,        only: mast
   use matm_class,        only: matm
   use vfs_class,         only: vfs
   use timetracker_class, only: timetracker
   use monitor_class,     only: monitor
   use ensight_class,     only: ensight
   use event_class,       only: event
   use string,            only: str_medium
   use hypre_str_class,   only: hypre_str
   use hypre_uns_class,   only: hypre_uns
   use mathtools,         only: Pi
   use iterator_class,    only: iterator
   use surfmesh_class,    only: surfmesh
   use pardata_class,     only: pardata
   use inputfile_class,   only: inputfile
   implicit none
   private
   
   public :: ml
    
   !> SML object
   type :: ml
      !> Config
      type(config)      :: cfg   !< Mesh for solver
      !> Flow solver
      type(mast)        :: fs    !< Incompressible flow solver
      type(vfs)         :: vf
      type(hypre_str)   :: ps    
      type(hypre_str)   :: vs
      ! type(hypre_uns)   :: ps    
      ! type(hypre_uns)   :: vs    
      type(timetracker) :: time  !< Time info
      type(matm)        :: matmod
      type(iterator) :: top_layer, btm_layer 
      !> Ensight postprocessing
      type(ensight)     :: ens_out
      type(event)       :: ens_evt
      type(event)       :: ppevt
      type(surfmesh)    :: smesh
      !> Simulation monitor file
      type(monitor)     :: mfile,cflfile,turbfile
      !> Work arrays
      real(WP), dimension(:,:,:), allocatable :: recU,recV,recW
      !> Fluid parameters
      real(WP) :: visc
      !> Choice of relaxation model
      integer :: relax_model
      !> Postproc output
      integer :: junit
      !> Provide a pardata and an event tracker for saving restarts
      type(event)     :: save_evt
      type(pardata)   :: df
      logical         :: restarted
      !> Input file for the simulation
      type(inputfile) :: input
   contains
      procedure :: init          !< Initialize SML simulation
      procedure :: step          !< Advance SML simulation by one time step
      procedure :: final         !< Finalize SML simulation
   end type ml
    
   !> Problem Definition
   real(WP) :: Ls,Lyg,Lyl,Ma_c,Re_d,m_thick,TKE,Re_m,r_rho,r_vel,r_visc,turb_lengthV,turb_lengthP,r_kappa
   real(WP) :: GP,Grho,Ma_g,gamm_g,Reg,Weg,visc_g,g_thick,Pr_g,cv_g0,kappa_g
   real(WP) :: LP,Lrho,Ma_l,gamm_l,Pref_l,visc_l,l_thick,Pr_l,cv_l0,kappa_l
   !> Turbulence stuff
   real(WP) :: myTKE,EPS_min,myEPS,eta,EPS_check
   real(WP), dimension(:), allocatable :: KE_old
   !> Postproc output
   integer :: junit
   character(len=str_medium) :: filename,timestamp

   contains

   !> Function that defines a level set function for a initial wavy interface
   function levelset_wavy(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      ! G=-xyz(2)
      G=Lyg-xyz(2)
   end function levelset_wavy
   
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
 
   !> Function that localizes the bottom (y-) of the domain
   function btm_of_domain(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmin-1) isIn=.true.
   end function btm_of_domain
 
   !> Function that localizes top sponge
   function top_sponge(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (pg%y(pg%jmax+1)-pg%ym(j).le.Ls) isIn=.true.
   end function top_sponge
 
   !> Function that localizes bottom sponge
   function btm_sponge(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (pg%ym(j)-pg%y(pg%jmin).le.Ls) isIn=.true.
   end function btm_sponge
    
   !> Initialization of SML simulation
   subroutine init(this)
      implicit none
      class(ml), intent(inout) :: this

      ! Setup an input file
      read_input: block
         use parallel, only: amRoot
         this%input=inputfile(amRoot=amRoot,filename='input')
      end block read_input
      
      ! Create the SML mesh
      create_config: block
         use sgrid_class, only: cartesian,sgrid
         use param,       only: param_read,param_exists
         use parallel,    only: group
         real(WP), dimension(:), allocatable :: x,y,z
         integer, dimension(3) :: partition
         type(sgrid) :: grid
         integer :: i,j,k,nx,ny,nz,ny_uni
         real(WP) :: Lx,Ly,Lz,dy,Lyg,Lyl,dy_tmp,y_uni,ry,y_tmp

         ! Read in grid definition
         call param_read('Lx',Lx); call param_read('nx',nx); allocate(x(nx+1))
         call param_read('Liquid Ly',Lyl); call param_read('Gas Ly',Lyg); call param_read('ny',ny);allocate(y(ny+1))
         call param_read('Lz',Lz); call param_read('nz',nz); allocate(z(nz+1))
         ! Create simple rectilinear grid for x and z
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx-0.5_WP*Lx
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         do j=1,ny+1
           ! y(j)=real(j-1,WP)/real(ny,WP)*(Lyl+Lyg)-0.2_WP*(Lyl+Lyg)
           y(j)=real(j-1,WP)/real(ny,WP)*(Lyl+Lyg)-(Lyl/(Lyl+Lyg))*(Lyl+Lyg)
         end do
        !  ! Create grid in y, with stretching if specified
        !  dy = (Lyl + Lyg)/ny
        !  if (param_exists('Stretching ratio in y')) then
        !     call param_read('Stretching ratio in y',ry)
        !     call param_read('Stretching begins at y',y_uni)
         !     ! Determine actual number of points to be used
        !     ny = 0
        !     y_tmp = 0.0_WP
        !     dy_tmp = dy
        !     do while (y_tmp.lt.0.5_WP*(Lyl+Lyg))
        !        ny = ny+2 ! Symmetric
        !        if (y_tmp.gt.y_uni) dy_tmp = dy_tmp*ry
        !        y_tmp = y_tmp + dy_tmp
        !     end do
        !     ny_uni = ceiling(y_uni/dy)
        !     if (this%cfg%amRoot) print*,'ny with stretched mesh',ny
        !     allocate(y(ny+1))
        !     y(ny/2+1) = 0.0_WP
        !     do j=ny/2+2,ny/2+ny_uni+1
        !        y(j) = y(j-1)+dy
        !     end do
        !     do j=ny/2+ny_uni+2,ny+1
        !        dy = dy*ry
        !        y(j) = y(j-1)+dy
        !     end do
        !     do j=1,ny/2
        !        y(j) = -y(ny+2-j)
        !     end do
        !  else
        !     ! Uniform mesh
        !     allocate(y(ny+1))
        !     do j=1,ny+1
        !        y(j)=real(j-1,WP)/real(ny,WP)*(Lyl+Lyg)-Lyl
        !     end do
        !  end if
         ! Read in sponge length
         call param_read('Sponge Zone',Ls)
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.true.,yper=.false.,zper=.true.,name='ml')
         ! Read in partition
         call param_read('Partition',partition,short='p')
         ! Create partitioned grid without walls
         this%cfg=config(grp=group,decomp=partition,grid=grid)
      end block create_config
      
      ! Initialize the work arrays
      allocate_work_arrays: block
         allocate(this%recU(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%recV(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%recW(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(KE_old(this%cfg%jmin:this%cfg%jmax)); KE_old=0.0_WP
      end block allocate_work_arrays
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         use param, only: param_read
         this%time=timetracker(amRoot=this%cfg%amRoot)
         call param_read('Max timestep size',this%time%dtmax)
         call param_read('Max cfl number',this%time%cflmax)
         this%time%dt=this%time%dtmax
         this%time%itmax=2
      end block initialize_timetracker
      
      ! Handle restart/saves here
      restart_and_save: block
      use string,   only: str_medium
      use filesys,  only: makedir,isdir
      character(len=str_medium)  :: timestamp
      integer, dimension(3)  :: iopartition
      ! Create event for saving restart files
      this%save_evt=event(this%time,'Restart output')
      call this%input%read('Restart output period',this%save_evt%tper)
      ! Check if we are restarting
      call this%input%read('Restart from',timestamp,default='')
      this%restarted=.false.; if (len_trim(timestamp).gt.0) this%restarted=.true.
      ! Read in the I/O partition
      call this%input%read('I/O partition',iopartition)
      ! Perform pardata initialization
      if (this%restarted) then
         ! We are restarting, read the file
         call this%df%initialize(pg=this%cfg,iopartition=iopartition,fdata='restart/data_'//trim(adjustl(timestamp)))
      else
         ! We are not restarting, prepare a new directory for storing restart files
         if (this%cfg%amRoot) then
            if (.not.isdir('restart')) call makedir('restart')
         end if
         ! Prepare pardata object for saving restart files
         call this%df%initialize(pg=this%cfg,iopartition=iopartition,filename=trim(this%cfg%name),nval=2,nvar=29)
         this%df%valname=['t ','dt']
         this%df%varname=['Grho   ','Lrho   ','RHO    ','Ui     ','Vi     ','Wi     ','U      ','V      ','W      ','rhoUi  ','rhoVi  ','rhoWi  ','GrhoE  ','LrhoE  ','GP     ','LP     ','P      ','PA     ','GrhoSS2','LrhoSS2','RHOSS2 ','P11    ','P12    ','P13    ','P14    ','VOF    '] ! 'GP     ','LP     ','P      ','PA     ','Pjx    ','Pjy    ','Pjz    ','Tmptr  ','GrhoSS2','LrhoSS2'
      end if
      end block restart_and_save

      ! Revisit timetracker to adjust time and time step values if this is a restart
      update_timetracker: block
      if (this%restarted) then
         call this%df%pull(name='t' ,val=this%time%t )
         call this%df%pull(name='dt',val=this%time%dt)
         this%time%told=this%time%t-this%time%dt
      end if
      end block update_timetracker

      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
      use mms_geom,  only: cube_refine_vol
      use vfs_class, only: lvira,VFhi,VFlo,plicnet,flux,neumann
      use mathtools, only: Pi,twoPi
      use random,    only: random_uniform
      use parallel,  only: MPI_REAL_WP
      use mpi_f08
      use irl_fortran_interface
      integer :: i,j,k,n,si,sj,sk,ierr
      real(WP), dimension(3,8) :: cube_vertex
      real(WP), dimension(3) :: v_cent,a_cent
      real(WP) :: vol,area
      integer, parameter :: amr_ref_lvl=4
      real(WP), dimension(:,:,:), allocatable :: P11,P12,P13,P14
      ! Create a VOF solver
      call this%vf%initialize(cfg=this%cfg,reconstruction_method=plicnet,transport_method=flux,name='VOF')
      ! ! Initialize the interface including restarts
      if (this%restarted) then
        ! Read in the planes directly and set the IRL interface
        allocate(P11(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P11',var=P11)
        allocate(P12(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P12',var=P12)
        allocate(P13(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P13',var=P13)
        allocate(P14(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P14',var=P14)
        do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
           do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
              do i=this%vf%cfg%imino_,this%vf%cfg%imaxo_
                    call setNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k),1)
                    call setPlane(this%vf%liquid_gas_interface(i,j,k),0,[P11(i,j,k),P12(i,j,k),P13(i,j,k)],P14(i,j,k))
              end do
           end do
        end do
        deallocate(P11,P12,P13,P14)
      !   call this%cfg%sync(this%vf%VF)
        call this%df%pull(name='VOF',var=this%vf%VF)
        ! Ensure that boundaries are correct
        this%vf%VF(:,this%vf%cfg%jmino:this%vf%cfg%jmin+1,:)=1.0_WP
        this%vf%VF(:,this%vf%cfg%jmax:this%vf%cfg%jmaxo,:)=0.0_WP
        call this%vf%add_bcond(name='ybottom',type=neumann,locator=btm_of_domain,dir='ym')
        call this%vf%apply_bcond(this%time%t,this%time%dt)
        call this%vf%sync_interface()
        call this%vf%set_full_bcond()
        ! Reset moments
        call this%vf%reset_volume_moments()
        ! Ensure that boundaries are correct
      !   if (this%vf%cfg%jproc.eq.1)               this%vf%VF(:,this%vf%cfg%jmino:this%vf%cfg%jmin+1,:)=1.0_WP
      !   if (this%vf%cfg%jproc.eq.this%vf%cfg%npy) this%vf%VF(:,this%vf%cfg%jmax:this%vf%cfg%jmaxo,:)=0.0_WP
      !   call this%cfg%sync(this%vf%VF)
      !   call this%df%pull(name='VOF',var=this%vf%VF)
        ! Update the band
        call this%vf%update_band()
        ! Set interface planes at the boundaries
        call this%vf%set_full_bcond()
        ! Create discontinuous polygon mesh from IRL interface
        call this%vf%polygonalize_interface()
        ! Calculate distance from polygons
        call this%vf%distance_from_polygon()
        ! Calculate subcell phasic volumes
        call this%vf%subcell_vol()
        ! Calculate curvature
        call this%vf%get_curvature()
      else
      ! Initialize unpeturbed interface
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
               call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_wavy,0.0_WP,amr_ref_lvl)
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
      ! Set interface at the boundaries
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
      end if
      end block create_and_initialize_vof
       
       
      ! Create a compressible, multiphase flow solver with clipped Neumann conditions top and bottom
      create_flow_solver: block
         use param,           only: param_read
         use hypre_str_class, only: pcg_pfmg
         use hypre_uns_class, only: pcg_amg
         use mast_class,      only: mech_egy_mech_hhz,clipped_neumann,bc_scope,bcond,thermmech_egy_mech_hhz
         use matm_class,      only: none
         integer :: i,j,k
         real(WP) :: c
         real(WP) :: liq_vol,gas_vol,tot_vol
         ! Create material model class
         this%matmod=matm(cfg=this%cfg,name='Liquid-gas models')
         ! Get nondimensional parameters from input (compressible multiphase case)
         call param_read('Gas gamma',gamm_g)
         call param_read('Liquid gamma',gamm_l);
         call param_read('Gas Reynolds number',Reg); visc_g=(1.0_WP)/(Reg+epsilon(Reg))
         call param_read('Viscosity ratio',r_visc); ! visc_l=r_visc*visc_g
         call param_read('Density ratio',r_rho); Grho=1.0_WP;Lrho = r_rho*Grho
         call param_read('Velocity ratio',r_vel); visc_l = visc_g*(1.0_WP-r_vel)/r_vel
         call param_read('Gas Mach number',Ma_g); GP = (1.0_WP**2.0_WP)/(gamm_g*Ma_g**2.0_WP); LP = GP
         call param_read('Liquid Mach number',Ma_l); Pref_l = ((r_rho*1.0_WP**2.0_WP)/(gamm_l*Ma_l**2.0_WP)) - LP
         cv_g0 = GP/(Grho*1.0_WP*(gamm_g-1.0_WP)); cv_l0 = (LP+Pref_l)/(Lrho*1.0_WP*(gamm_l-1.0_WP))
         call param_read('Gas Prandtl number',Pr_g); kappa_g = gamm_g*cv_g0*visc_g/Pr_g
         call param_read('Liquid Prandtl number',Pr_l); kappa_l = gamm_l*cv_l0*visc_l/Pr_l
         ! kappa_g = 0.0_WP; kappa_l = 0.0_WP
         ! print*, GP, Pref_l, cv_g0, cv_l0, kappa_g, kappa_l
         ! Register equations of state
         call this%matmod%register_stiffenedgas('liquid',gamm_l,Pref_l)
         call this%matmod%register_idealgas('gas',gamm_g)
         ! Create flow solver
         this%fs=mast(cfg=this%cfg,name='Two-phase All-Mach',vf=this%vf)
         ! Register flow solver variables with material models
         call this%matmod%register_thermoflow_variables('liquid',this%fs%Lrho,this%fs%Ui,this%fs%Vi,this%fs%Wi,this%fs%LrhoE,this%fs%LP)
         call this%matmod%register_thermoflow_variables('gas'   ,this%fs%Grho,this%fs%Ui,this%fs%Vi,this%fs%Wi,this%fs%GrhoE,this%fs%GP)
         call this%matmod%register_diffusion_thermo_models(viscconst_gas=visc_g,viscconst_liquid=visc_l,hdffconst_gas=kappa_g,hdffconst_liquid=kappa_l,sphtconst_gas=cv_g0,sphtconst_liquid=cv_l0)
         call param_read('Gas Weber number',Weg); this%fs%sigma=(1.0_WP**3.0_WP)/(Weg+epsilon(Weg)); 
         ! Configure pressure solver
         this%ps=hypre_str(cfg=this%cfg,name='Pressure',method=pcg_pfmg,nst=7)
         this%ps%maxlevel=12
         call param_read('Pressure iteration',this%ps%maxit)
         call param_read('Pressure tolerance',this%ps%rcvg)
         ! Configure implicit velocity solver
         this%vs=hypre_str(cfg=this%cfg,name='Velocity',method=pcg_pfmg,nst=7)
         call param_read('Implicit iteration',this%vs%maxit)
         call param_read('Implicit tolerance',this%vs%rcvg)
         ! Setup the solver
         call this%fs%setup(pressure_solver=this%ps,implicit_solver=this%vs)
         this%fs%Lrho = Lrho; this%fs%LP = LP
         this%fs%Grho = Grho; this%fs%GP = GP
         ! Solve for initial constant to enforce d_g = 1.0
         !  c = -((1.0_WP + r_vel)**2.0_WP)/(r_vel*LOG(2.0_WP) + LOG(2.0_WP) - 1.0_WP)
         c = 1.0_WP / ((1.0_WP - LOG(2.0_WP))*(r_vel-1.0_WP)**2.0_WP)
         if (this%restarted) then
            ! Read data
            call this%df%pull(name='Grho'   ,var=this%fs%Grho   ); 
            call this%df%pull(name='Lrho'   ,var=this%fs%Lrho   ); 
            call this%df%pull(name='RHO'    ,var=this%fs%RHO   ); 
            call this%df%pull(name='Ui'     ,var=this%fs%Ui     ); 
            call this%df%pull(name='Vi'     ,var=this%fs%Vi     ); 
            call this%df%pull(name='Wi'     ,var=this%fs%Wi     ); 
            call this%df%pull(name='U'      ,var=this%fs%U      ); 
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
            ! call this%df%pull(name='Pjx'    ,var=this%fs%Pjx    ); 
            ! call this%df%pull(name='Pjy'    ,var=this%fs%Pjy    ); 
            ! call this%df%pull(name='Pjz'    ,var=this%fs%Pjz    ); 
            ! call this%df%pull(name='Tmptr'  ,var=this%fs%Tmptr  ); 
            call this%df%pull(name='GrhoSS2',var=this%fs%GrhoSS2); 
            call this%df%pull(name='LrhoSS2',var=this%fs%LrhoSS2); 
            call this%df%pull(name='RHOSS2',var=this%fs%RHOSS2); 
            ! Apply boundary conditions
            ! Define BCs at top and bottom - though sponges will determine bdy behavior
            call this%fs%add_bcond(name='top_y'   ,type=clipped_neumann,locator=top_of_domain  ,celldir='yp')
            call this%fs%add_bcond(name='btm_y'   ,type=clipped_neumann,locator=btm_of_domain  ,celldir='ym')
            ! call this%fs%apply_bcond(this%time%dt,'density')
            ! call this%fs%apply_bcond(this%time%dt,'momentum')
            ! call this%fs%apply_bcond(this%time%dt,'energy')
            ! call this%fs%apply_bcond(this%time%dt,'velocity')
            ! call apply_sponges(this)
            ! ! Treat top boundaries
            do k=this%fs%cfg%kmino_,this%fs%cfg%kmaxo_
               do j=this%fs%cfg%jmax_-1,this%fs%cfg%jmaxo_
                  do i=this%fs%cfg%imino_,this%fs%cfg%imaxo_
                     this%fs%U(i,j,k) = 1.0_WP; this%fs%V(i,j,k) = 0.0_WP; this%fs%W(i,j,k) = 0.0_WP
                     this%fs%Ui(i,j,k) = 1.0_WP; this%fs%Vi(i,j,k) = 0.0_WP; this%fs%Wi(i,j,k) = 0.0_WP
                     this%fs%rhoUi(i,j,k) = 1.0_WP; this%fs%rhoVi(i,j,k) = 0.0_WP; this%fs%rhoWi(i,j,k) = 0.0_WP
                     this%fs%Grho(i,j,k) = Grho; this%fs%GP(i,j,k) = GP; this%fs%Tmptr(i,j,k) = 1.0_WP;
                     this%fs%GrhoE(i,j,k) = this%matmod%EOS_energy(GP,Grho,1.0_WP,0.0_WP,0.0_WP,'gas')
                     this%fs%GrhoSS2(i,j,k) = this%matmod%EOS_gas(i,j,k,'M')
                  end do
               end do
            end do
            ! Treat bottom boundaries
            do k=this%fs%cfg%kmino_,this%fs%cfg%kmaxo_
               do j=this%fs%cfg%jmino_,this%fs%cfg%jmin_+1
                  do i=this%fs%cfg%imino_,this%fs%cfg%imaxo_
                     this%fs%U(i,j,k) = 0.0_WP; this%fs%V(i,j,k) = 0.0_WP; this%fs%W(i,j,k) = 0.0_WP
                     this%fs%Ui(i,j,k) = 0.0_WP; this%fs%Vi(i,j,k) = 0.0_WP; this%fs%Wi(i,j,k) = 0.0_WP
                     this%fs%rhoUi(i,j,k) = 0.0_WP; this%fs%rhoVi(i,j,k) = 0.0_WP; this%fs%rhoWi(i,j,k) = 0.0_WP
                     this%fs%Lrho(i,j,k) = Lrho; this%fs%LP(i,j,k) = LP; this%fs%Tmptr(i,j,k) = 1.0_WP;
                     this%fs%LrhoE(i,j,k) = this%matmod%EOS_energy(LP,Lrho,0.0_WP,0.0_WP,0.0_WP,'liquid')
                     this%fs%LrhoSS2(i,j,k) = this%matmod%EOS_liquid(i,j,k,'M')
                  end do
               end do
            end do
            ! Calculate mixture density and momenta
            this%fs%RHO   = (1.0_WP-this%vf%VF)*this%fs%Grho  + this%vf%VF*this%fs%Lrho
            ! Choose relaxation procedure
            this%relax_model = thermmech_egy_mech_hhz  ! thermmech_egy_mech_hhz
            ! Calculate cell center and face viscosity
            do k=this%fs%cfg%kmin_-1,this%fs%cfg%kmax_+2
               do j=this%fs%cfg%jmin_-1,this%fs%cfg%jmax_+2
                  do i=this%fs%cfg%imin_-1,this%fs%cfg%imax_+2
                     !! -- CELL CENTER -- !!
                     liq_vol=sum(this%vf%Lvol(:,:,:,i,j,k))
                     gas_vol=sum(this%vf%Gvol(:,:,:,i,j,k))
                     tot_vol=gas_vol+liq_vol
                     if (tot_vol.gt.0.0_WP) then
                        this%fs%therm_cond(i,j,k)=kappa_g*kappa_l/(kappa_l*gas_vol/tot_vol+kappa_g*liq_vol/tot_vol+epsilon(1.0_WP))
                        this%fs%visc(i,j,k)=visc_g*visc_l/(visc_l*gas_vol/tot_vol+visc_g*liq_vol/tot_vol+epsilon(1.0_WP))
                     end if
                  end do
               end do
            end do
            this%fs%Tmptr = 1.0_WP;
            call this%matmod%update_temperature(this%VF,this%fs%Tmptr)
            ! Perform initial pressure relax
            ! call this%fs%pressure_relax(this%vf,this%matmod,this%relax_model)
            ! Calculate initial phase and bulk moduli
            ! call this%fs%init_phase_bulkmod(this%vf,this%matmod)
            ! call this%fs%reinit_phase_pressure(this%vf,this%matmod)
            ! call this%fs%harmonize_advpressure_bulkmod(this%vf,this%matmod)
            ! Set initial pressure to harmonized field based on internal energy
            ! this%fs%P = this%fs%PA
            ! Initialize first guess for pressure (0 works best)
            ! this%fs%psolv%sol=0.0_WP

         else
            ! Set initial velocity field
            this%fs%Vi=0.0_WP;this%fs%Wi=0.0_WP
            do k=this%fs%cfg%kmino_,this%fs%cfg%kmaxo_
               do j=this%fs%cfg%jmino_,this%fs%cfg%jmaxo_
                  do i=this%fs%cfg%imino_,this%fs%cfg%imaxo_
                     if (this%fs%cfg%y(j).le.0.0_WP) then
                        this%fs%Ui(i,j,k)    = r_vel*tanh(this%fs%cfg%ym(j)/c) + r_vel
                        ! this%fs%LrhoE(i,j,k) = this%matmod%EOS_energy(LP,Lrho,this%fs%Ui(i,j,k),this%fs%Vi(i,j,k),this%fs%Wi(i,j,k),'liquid')
                        ! this%fs%LP(i,j,k) = LP
                     else
                        ! this%fs%Ui(i,j,k)    = tanh(this%fs%cfg%ym(j)/c)
                        this%fs%Ui(i,j,k)    = (1.0_WP - r_vel)*tanh(this%fs%cfg%ym(j)/c) + r_vel
                        ! this%fs%GrhoE(i,j,k) = this%matmod%EOS_energy(GP,Grho,this%fs%Ui(i,j,k),this%fs%Vi(i,j,k),this%fs%Wi(i,j,k),'gas')
                        ! this%fs%GP(i,j,k) = GP
                     end if
                     this%fs%GrhoE(i,j,k) = this%matmod%EOS_energy(GP,Grho,this%fs%Ui(i,j,k),this%fs%Vi(i,j,k),this%fs%Wi(i,j,k),'gas')
                     this%fs%LrhoE(i,j,k) = this%matmod%EOS_energy(LP,Lrho,this%fs%Ui(i,j,k),this%fs%Vi(i,j,k),this%fs%Wi(i,j,k),'liquid')
                  end do
               end do
            end do
            ! this%fs%Lrho = Lrho; this%fs%LP = LP
            ! this%fs%Grho = Grho; this%fs%GP = GP
            ! Define BCs at top and bottom - though sponges will determine bdy behavior
            call this%fs%add_bcond(name='top_y'   ,type=clipped_neumann,locator=top_of_domain  ,celldir='yp')
            call this%fs%add_bcond(name='btm_y'   ,type=clipped_neumann,locator=btm_of_domain  ,celldir='ym')
            ! Calculate face velocities
            call this%fs%interp_vel_basic(this%vf,this%fs%Ui,this%fs%Vi,this%fs%Wi,this%fs%U,this%fs%V,this%fs%W)
            call this%matmod%update_temperature(this%vf,this%fs%Tmptr)
            ! Calculate mixture density and momenta
            this%fs%RHO   = (1.0_WP-this%vf%VF)*this%fs%Grho  + this%vf%VF*this%fs%Lrho
            this%fs%rhoUi = this%fs%RHO*this%fs%Ui; this%fs%rhoVi = this%fs%RHO*this%fs%Vi; this%fs%rhoWi = this%fs%RHO*this%fs%Wi
            ! Choose relaxation procedure
            this%relax_model = thermmech_egy_mech_hhz  ! thermmech_egy_mech_hhz
            ! Perform initial pressure relax
            call this%fs%pressure_relax(this%vf,this%matmod,this%relax_model)
            ! Calculate initial phase and bulk moduli
            call this%fs%init_phase_bulkmod(this%vf,this%matmod)
            call this%fs%reinit_phase_pressure(this%vf,this%matmod)
            call this%fs%harmonize_advpressure_bulkmod(this%vf,this%matmod)
            ! Set initial pressure to harmonized field based on internal energy
            this%fs%P = this%fs%PA
            ! Initialize first guess for pressure (0 works best)
            this%fs%psolv%sol=0.0_WP
            ! call param_read('Liquid Ly',Lyl); call param_read('Gas Ly',Lyg)
            ! call apply_sponges(this)
            ! call turb_stats(this,m_thick,g_thick,l_thick,TKE,Re_m)
            ! Get Mach number
            ! call this%fs%get_viz()
            ! Make directory for postprocessing data
            if (this%fs%cfg%amRoot) then
               call execute_command_line('mkdir -p PostProc')
            end if
         end if
      end block create_flow_solver

      ! Create surfmesh object for interface polygon output
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

      ! Add Ensight output
      create_ensight: block
         use param, only: param_read
         ! Create Ensight output from cfg
         this%ens_out=ensight(cfg=this%cfg,name='cmml')
         ! Create event for Ensight output
         this%ens_evt=event(time=this%time,name='Ensight output')
         call param_read('Ensight output period',this%ens_evt%tper)
         ! Add variables to output
          call this%ens_out%add_scalar('Density',this%fs%RHO)
         call this%ens_out%add_scalar('Grho',this%fs%Grho)
         call this%ens_out%add_scalar('Lrho',this%fs%Lrho)
         call this%ens_out%add_vector('velocity',this%fs%Ui,this%fs%Vi,this%fs%Wi)
         call this%ens_out%add_scalar('pressure',this%fs%P)
         call this%ens_out%add_scalar('GP',this%fs%GP)
         call this%ens_out%add_scalar('LP',this%fs%LP)
         call this%ens_out%add_scalar('VOF',this%vf%VF)
         call this%ens_out%add_scalar('curvature',this%vf%curv)
         call this%ens_out%add_scalar('LrhoE',this%fs%LrhoE)
         call this%ens_out%add_scalar('GrhoE',this%fs%GrhoE)
         ! call this%ens_out%add_scalar('visc',this%fs%visc)
         !  call this%ens_out%add_scalar('SL_x',this%fs%sl_x)
         !  call this%ens_out%add_scalar('SL_y',this%fs%sl_y)
         !  call this%ens_out%add_scalar('SL_z',this%fs%sl_z)
         call this%ens_out%add_surface('plic',this%smesh)
         !  call this%ens_out%add_vector('SL_trans',this%fs%sl_x,this%fs%sl_y,this%fs%sl_z)
         call this%ens_out%add_scalar('Mach',this%fs%Mach)
         call this%ens_out%add_scalar('Temp',this%fs%Tmptr)
         !  call this%ens_out%add_scalar('GBulkMod',this%fs%GrhoSS2)
         !  call this%ens_out%add_scalar('LBulkMod',this%fs%LrhoSS2)
         ! Output to ensight
         if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
      end block create_ensight
 
      ! Create monitoring file
      create_monitor: block
         ! Prepare some info about fields
         call this%fs%get_cfl(this%time%dt,this%time%cfl)
         call this%fs%get_max()
         call this%vf%get_max()
         call this%fs%get_viz()
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
         !  call this%mfile%add_column(this%fs%divmax,'Maximum divergence')
         call this%mfile%add_column(this%fs%psolv%it,'Pressure iteration')
         call this%mfile%add_column(this%fs%psolv%rerr,'Pressure error')
         call this%mfile%write()
         ! Create turbulent statistics monitor
         this%turbfile=monitor(this%fs%cfg%amRoot,'Turb_Stats')
         call this%turbfile%add_column(this%time%n,'Timestep number')
         call this%turbfile%add_column(this%time%t,'Time')
         call this%turbfile%add_column(m_thick,'M_thick')
         call this%turbfile%add_column(TKE,'TKE')
         call this%turbfile%add_column(EPS_min,'EPS')
         call this%turbfile%add_column(EPS_check,'EPS_check')
         call this%turbfile%add_column(eta,'eta')
         call this%turbfile%add_column(Re_m,'Re_m')
         call this%turbfile%add_column(g_thick,'g_thick')
         call this%turbfile%add_column(l_thick,'l_thick')
         call this%turbfile%add_column(turb_lengthV,'V Correlation')
         call this%turbfile%add_column(turb_lengthP,'P Correlation')
         call this%turbfile%write()
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
      end block create_monitor

      ! Create a specialized post-processing file
      create_postproc: block
         use param,      only: param_read
         ! Create event for data postprocessing
         this%ppevt=event(time=this%time,name='Postproc output')
         call param_read('Postproc output period',this%ppevt%tper)
      end block create_postproc 

      ! Create an iterator for imposing the bottom and top sponges
      create_iterator_top: block
         this%top_layer = iterator(pg=this%fs%cfg,name='Top Sponge',locator=top_sponge)
      end block create_iterator_top
  
      create_iterator_btm: block
         this%btm_layer = iterator(pg=this%fs%cfg,name='Bottom Sponge',locator=btm_sponge)
      end block create_iterator_btm
       
   end subroutine init
    
 
   !> Take one time step with specified dt
   subroutine step(this)
      implicit none
      class(ml), intent(inout) :: this
      
      ! Increment time
      call this%fs%get_cfl(this%time%dt,this%time%cfl)
      call this%time%adjust_dt()
      call this%time%increment()
      
      ! Reinitialize phase pressure by syncing it with conserved phase energy
      call this%fs%reinit_phase_pressure(this%vf,this%matmod)
      this%fs%Uiold=this%fs%Ui; this%fs%Viold=this%fs%Vi; this%fs%Wiold=this%fs%Wi
      this%fs%RHOold = this%fs%RHO
      ! Remember old flow variables (phase)
      this%fs%Grhoold = this%fs%Grho;  this%fs%Lrhoold = this%fs%Lrho
      this%fs%GrhoEold= this%fs%GrhoE; this%fs%LrhoEold= this%fs%LrhoE
      this%fs%GPold   = this%fs%GP;    this%fs%LPold   = this%fs%LP
      
      ! Remember old interface, including VF and barycenters
      call this%vf%copy_interface_to_old()

      ! Create in-cell reconstruction
      call this%fs%flow_reconstruct(this%vf)

      ! Zero variables that will change during subiterations
      this%fs%P = 0.0_WP
      this%fs%Pjx = 0.0_WP; this%fs%Pjy = 0.0_WP; this%fs%Pjz = 0.0_WP
      this%fs%Hpjump = 0.0_WP

      ! Determine semi-Lagrangian advection flag
      call this%fs%flag_sl(this%time%dt,this%vf)
       
      ! Perform sub-iterations
      do while (this%time%it.le.this%time%itmax)
         
         ! Predictor step, involving advection and pressure terms
         call this%fs%advection_step(this%time%dt,this%vf,this%matmod)

         ! Diffusion and built-in source term (gravity) step
         call this%fs%diffusion_src_explicit_step(this%time%dt,this%vf,this%matmod)
       
         ! Perform sponge forcing
         call apply_sponges(this)
       
         ! Prepare pressure projection
         call this%fs%pressureproj_prepare(this%time%dt,this%vf,this%matmod) 
       
         ! Initialize and solve Helmholtz equation
         call this%fs%psolv%setup()
         call this%fs%psolv%solve()
         call this%fs%cfg%sync(this%fs%psolv%sol)
       
         ! Perform corrector step using solution
         this%fs%P=this%fs%P+this%fs%psolv%sol
         call this%fs%pressureproj_correct(this%time%dt,this%vf,this%fs%psolv%sol)
       
         ! Record convergence monitor
         !   call this%cvgfile%write()
       
         ! Increment sub-iteration counter
         this%time%it=this%time%it+1
           
      end do
       
      ! Pressure relaxation
      call this%fs%pressure_relax(this%vf,this%matmod,this%relax_model)
 
      ! Output to ensight
      this%fs%PA = this%matmod%EOS_all(this%vf);
      if (this%ens_evt%occurs()) then
         ! update surfmesh object
         update_smesh: block
            use irl_fortran_interface
				integer :: i,j,k,nplane,np
				! Transfer polygons to smesh
				call this%vf%update_surfmesh(this%smesh)
				! Also populate nplane variable
				this%smesh%var(1,:)=0.0_WP
				np=0
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
				end block update_smesh
         call this%ens_out%write_data(this%time%t)
      end if
       
      ! Perform and output monitoring
      call this%fs%get_max()
      call this%vf%get_max()
      call this%fs%get_viz()
      call this%mfile%write()
      call this%cflfile%write()
      ! Perform and output monitoring
      if (this%ppevt%occurs()) then
         call turb_stats(this,m_thick,g_thick,l_thick,Re_m,turb_lengthV,turb_lengthP,TKE,EPS_min,eta,EPS_check,KE_old)
         call this%turbfile%write()
      end if

      ! Finally, see if it's time to save restart files
      if (this%save_evt%occurs()) then
         save_restart: block
         use irl_fortran_interface
         use string, only: str_medium
         character(len=str_medium) :: timestamp
         real(WP), dimension(:,:,:), allocatable :: P11,P12,P13,P14
         integer :: i,j,k
         real(WP), dimension(4) :: plane
         ! Handle IRL data
         allocate(P11(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(P12(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(P13(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(P14(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         ! Store IRL data
         do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
            do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
               do i=this%vf%cfg%imino_,this%vf%cfg%imaxo_
                  ! First plane
                  plane=getPlane(this%vf%liquid_gas_interface(i,j,k),0)
                  P11(i,j,k)=plane(1); P12(i,j,k)=plane(2); P13(i,j,k)=plane(3); P14(i,j,k)=plane(4)
                  plane=0.0_WP
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
         ! call this%df%push(name='Pjx'    ,var=this%fs%Pjx    ) 
         ! call this%df%push(name='Pjy'    ,var=this%fs%Pjy    ) 
         ! call this%df%push(name='Pjz'    ,var=this%fs%Pjz    )
         ! call this%df%push(name='Tmptr'  ,var=this%fs%Tmptr  ) 
         call this%df%push(name='GrhoSS2',var=this%fs%GrhoSS2) 
         call this%df%push(name='LrhoSS2',var=this%fs%LrhoSS2) 
         call this%df%push(name='RHOSS2',var=this%fs%RHOSS2) 
         call this%df%push(name='P11'    ,var=P11            )
         call this%df%push(name='P12'    ,var=P12            )
         call this%df%push(name='P13'    ,var=P13            )
         call this%df%push(name='P14'    ,var=P14            )
         call this%df%push(name='VOF'    ,var=this%vf%VF     )
         call this%df%write(fdata='restart/data_'//trim(adjustl(timestamp)))
         end block save_restart
      end if
       
   end subroutine step
 
   !> Finalize nozzle simulation
   subroutine final(this)
      implicit none
      class(ml), intent(inout) :: this
      
      ! Deallocate work arrays
      deallocate(this%recU,this%recV,this%recW)
      
   end subroutine final
    
    
   ! !> Function that defines a level set function for a initial wavy interface
   ! function levelset_wavy(xyz,t) result(G)
   !    implicit none
   !    real(WP), dimension(3),intent(in) :: xyz
   !    real(WP), intent(in) :: t
   !    real(WP) :: G
   !    ! G=-xyz(2)
   !    G=max(G,Lyg-xyz(2))
   ! end function levelset_wavy
   
   ! !> Function that localizes the top (y+) of the domain
   ! function top_of_domain(pg,i,j,k) result(isIn)
   !    use pgrid_class, only: pgrid
   !    implicit none
   !    class(pgrid), intent(in) :: pg
   !    integer, intent(in) :: i,j,k
   !    logical :: isIn
   !    isIn=.false.
   !    if (j.eq.pg%jmax+1) isIn=.true.
   ! end function top_of_domain
 
   ! !> Function that localizes the bottom (y-) of the domain
   ! function btm_of_domain(pg,i,j,k) result(isIn)
   !    use pgrid_class, only: pgrid
   !    implicit none
   !    class(pgrid), intent(in) :: pg
   !    integer, intent(in) :: i,j,k
   !    logical :: isIn
   !    isIn=.false.
   !    if (j.eq.pg%jmin-1) isIn=.true.
   ! end function btm_of_domain
 
   ! !> Function that localizes top sponge
   ! function top_sponge(pg,i,j,k) result(isIn)
   !    use pgrid_class, only: pgrid
   !    implicit none
   !    class(pgrid), intent(in) :: pg
   !    integer, intent(in) :: i,j,k
   !    logical :: isIn
   !    isIn=.false.
   !    if (pg%y(pg%jmax+1)-pg%ym(j).le.Ls) isIn=.true.
   ! end function top_sponge
 
   ! !> Function that localizes bottom sponge
   ! function btm_sponge(pg,i,j,k) result(isIn)
   !    use pgrid_class, only: pgrid
   !    implicit none
   !    class(pgrid), intent(in) :: pg
   !    integer, intent(in) :: i,j,k
   !    logical :: isIn
   !    isIn=.false.
   !    if (pg%ym(j)-pg%y(pg%jmin).le.Ls) isIn=.true.
   ! end function btm_sponge
 
   subroutine apply_sponges(this)
      use mathtools,  only: Pi
      implicit none
      class(ml), intent(inout) :: this
      integer :: i,j,k,n,m
      logical :: in_sponge_top, in_sponge_btm
      real(WP) :: swt_top, swt_btm, psponge, rhos, usponge, vsponge, wsponge

      ! Apply sponges using iterators
      do n=1,this%top_layer%n_
         swt_top = min(1.0_WP,max(0.0_WP,(this%fs%cfg%ym(this%top_layer%map(2,n)) - Lyg + Ls)/Ls))**2
         ! Get sponge solution if within a sponge
         psponge = GP
         rhos = Grho
         usponge = 1.0_WP; wsponge = 0.0_WP; vsponge = 0.0_WP
         ! Apply changes to variables
         this%fs%Grho (this%top_layer%map(1,n),this%top_layer%map(2,n),this%top_layer%map(3,n)) = this%fs%Grho (this%top_layer%map(1,n),this%top_layer%map(2,n),this%top_layer%map(3,n)) & 
            - swt_top*(this%fs%Grho(this%top_layer%map(1,n),this%top_layer%map(2,n),this%top_layer%map(3,n))-rhos)
         this%fs%Ui   (this%top_layer%map(1,n),this%top_layer%map(2,n),this%top_layer%map(3,n)) = this%fs%Ui   (this%top_layer%map(1,n),this%top_layer%map(2,n),this%top_layer%map(3,n)) & 
            - swt_top*(this%fs%Ui(this%top_layer%map(1,n),this%top_layer%map(2,n),this%top_layer%map(3,n))-usponge)
         this%fs%Vi   (this%top_layer%map(1,n),this%top_layer%map(2,n),this%top_layer%map(3,n)) = this%fs%Vi   (this%top_layer%map(1,n),this%top_layer%map(2,n),this%top_layer%map(3,n)) & 
            - swt_top*(this%fs%Vi(this%top_layer%map(1,n),this%top_layer%map(2,n),this%top_layer%map(3,n))-vsponge)
         this%fs%Wi   (this%top_layer%map(1,n),this%top_layer%map(2,n),this%top_layer%map(3,n)) = this%fs%Wi   (this%top_layer%map(1,n),this%top_layer%map(2,n),this%top_layer%map(3,n)) & 
            - swt_top*(this%fs%Wi(this%top_layer%map(1,n),this%top_layer%map(2,n),this%top_layer%map(3,n))-wsponge)
         this%fs%GP   (this%top_layer%map(1,n),this%top_layer%map(2,n),this%top_layer%map(3,n)) = this%fs%GP   (this%top_layer%map(1,n),this%top_layer%map(2,n),this%top_layer%map(3,n)) & 
            - swt_top*(this%fs%GP(this%top_layer%map(1,n),this%top_layer%map(2,n),this%top_layer%map(3,n))-psponge)
         this%fs%GrhoE(this%top_layer%map(1,n),this%top_layer%map(2,n),this%top_layer%map(3,n)) = this%fs%GrhoE(this%top_layer%map(1,n),this%top_layer%map(2,n),this%top_layer%map(3,n)) & 
            - swt_top*(this%fs%GrhoE(this%top_layer%map(1,n),this%top_layer%map(2,n),this%top_layer%map(3,n))-this%matmod%EOS_energy(psponge,rhos,usponge,vsponge,wsponge,'gas'))
         ! Update related quantities
         this%fs%RHO  (this%top_layer%map(1,n),this%top_layer%map(2,n),this%top_layer%map(3,n)) = this%fs%Grho(this%top_layer%map(1,n),this%top_layer%map(2,n),this%top_layer%map(3,n))
         this%fs%rhoUi(this%top_layer%map(1,n),this%top_layer%map(2,n),this%top_layer%map(3,n)) = this%fs%RHO(this%top_layer%map(1,n),this%top_layer%map(2,n),this%top_layer%map(3,n)) &
            *this%fs%Ui(this%top_layer%map(1,n),this%top_layer%map(2,n),this%top_layer%map(3,n))
         this%fs%rhoVi(this%top_layer%map(1,n),this%top_layer%map(2,n),this%top_layer%map(3,n)) = this%fs%RHO(this%top_layer%map(1,n),this%top_layer%map(2,n),this%top_layer%map(3,n)) &
            *this%fs%Vi(this%top_layer%map(1,n),this%top_layer%map(2,n),this%top_layer%map(3,n))
         this%fs%rhoWi(this%top_layer%map(1,n),this%top_layer%map(2,n),this%top_layer%map(3,n)) = this%fs%RHO(this%top_layer%map(1,n),this%top_layer%map(2,n),this%top_layer%map(3,n)) &
            *this%fs%Wi(this%top_layer%map(1,n),this%top_layer%map(2,n),this%top_layer%map(3,n))
         ! Remove liquid
         this%vf%VF   (this%top_layer%map(1,n),this%top_layer%map(2,n),this%top_layer%map(3,n)) = 0.0_WP
         this%fs%Lrho (this%top_layer%map(1,n),this%top_layer%map(2,n),this%top_layer%map(3,n)) = 0.0_WP
         this%fs%LP   (this%top_layer%map(1,n),this%top_layer%map(2,n),this%top_layer%map(3,n)) = 0.0_WP
         this%fs%LrhoE(this%top_layer%map(1,n),this%top_layer%map(2,n),this%top_layer%map(3,n)) = 0.0_WP
      end do

      do n=1,this%btm_layer%n_
         swt_btm = min(1.0_WP,max(0.0_WP,(-this%fs%cfg%ym(this%btm_layer%map(2,n)) - Lyl + Ls)/Ls))**2
         ! Get sponge solution if within a sponge
         psponge = LP
         rhos = Lrho
         usponge = 0.0_WP; wsponge = 0.0_WP; vsponge = 0.0_WP
         ! Apply changes to variables
         this%fs%Lrho (this%btm_layer%map(1,n),this%btm_layer%map(2,n),this%btm_layer%map(3,n)) = this%fs%Lrho (this%btm_layer%map(1,n),this%btm_layer%map(2,n),this%btm_layer%map(3,n)) &
            - swt_btm*(this%fs%Lrho(this%btm_layer%map(1,n),this%btm_layer%map(2,n),this%btm_layer%map(3,n))-rhos)
         this%fs%Ui   (this%btm_layer%map(1,n),this%btm_layer%map(2,n),this%btm_layer%map(3,n)) = this%fs%Ui   (this%btm_layer%map(1,n),this%btm_layer%map(2,n),this%btm_layer%map(3,n)) &
            - swt_btm*(this%fs%Ui(this%btm_layer%map(1,n),this%btm_layer%map(2,n),this%btm_layer%map(3,n))-usponge)
         this%fs%Vi   (this%btm_layer%map(1,n),this%btm_layer%map(2,n),this%btm_layer%map(3,n)) = this%fs%Vi   (this%btm_layer%map(1,n),this%btm_layer%map(2,n),this%btm_layer%map(3,n)) &
            - swt_btm*(this%fs%Vi(this%btm_layer%map(1,n),this%btm_layer%map(2,n),this%btm_layer%map(3,n))-vsponge)
         this%fs%Wi   (this%btm_layer%map(1,n),this%btm_layer%map(2,n),this%btm_layer%map(3,n)) = this%fs%Wi   (this%btm_layer%map(1,n),this%btm_layer%map(2,n),this%btm_layer%map(3,n)) &
            - swt_btm*(this%fs%Wi(this%btm_layer%map(1,n),this%btm_layer%map(2,n),this%btm_layer%map(3,n))-wsponge)
         this%fs%LP   (this%btm_layer%map(1,n),this%btm_layer%map(2,n),this%btm_layer%map(3,n)) = this%fs%LP   (this%btm_layer%map(1,n),this%btm_layer%map(2,n),this%btm_layer%map(3,n)) &
            - swt_btm*(this%fs%LP(this%btm_layer%map(1,n),this%btm_layer%map(2,n),this%btm_layer%map(3,n))-psponge)
         this%fs%LrhoE(this%btm_layer%map(1,n),this%btm_layer%map(2,n),this%btm_layer%map(3,n)) = this%fs%LrhoE(this%btm_layer%map(1,n),this%btm_layer%map(2,n),this%btm_layer%map(3,n)) &
            - swt_btm*(this%fs%LrhoE(this%btm_layer%map(1,n),this%btm_layer%map(2,n),this%btm_layer%map(3,n))-this%matmod%EOS_energy(psponge,rhos,usponge,vsponge,wsponge,'liquid'))
         ! Update related quantities
         this%fs%RHO  (this%btm_layer%map(1,n),this%btm_layer%map(2,n),this%btm_layer%map(3,n)) = this%fs%Lrho(this%btm_layer%map(1,n),this%btm_layer%map(2,n),this%btm_layer%map(3,n))
         this%fs%rhoUi(this%btm_layer%map(1,n),this%btm_layer%map(2,n),this%btm_layer%map(3,n)) = this%fs%RHO(this%btm_layer%map(1,n),this%btm_layer%map(2,n),this%btm_layer%map(3,n)) &
            *this%fs%Ui(this%btm_layer%map(1,n),this%btm_layer%map(2,n),this%btm_layer%map(3,n))
         this%fs%rhoVi(this%btm_layer%map(1,n),this%btm_layer%map(2,n),this%btm_layer%map(3,n)) = this%fs%RHO(this%btm_layer%map(1,n),this%btm_layer%map(2,n),this%btm_layer%map(3,n)) &
            *this%fs%Vi(this%btm_layer%map(1,n),this%btm_layer%map(2,n),this%btm_layer%map(3,n))
         this%fs%rhoWi(this%btm_layer%map(1,n),this%btm_layer%map(2,n),this%btm_layer%map(3,n)) = this%fs%RHO(this%btm_layer%map(1,n),this%btm_layer%map(2,n),this%btm_layer%map(3,n)) &
            *this%fs%Wi(this%btm_layer%map(1,n),this%btm_layer%map(2,n),this%btm_layer%map(3,n))
         ! Remove gas
         this%vf%VF   (this%btm_layer%map(1,n),this%btm_layer%map(2,n),this%btm_layer%map(3,n)) = 1.0_WP
         this%fs%Grho (this%btm_layer%map(1,n),this%btm_layer%map(2,n),this%btm_layer%map(3,n)) = 0.0_WP
         this%fs%GP   (this%btm_layer%map(1,n),this%btm_layer%map(2,n),this%btm_layer%map(3,n)) = 0.0_WP
         this%fs%GrhoE(this%btm_layer%map(1,n),this%btm_layer%map(2,n),this%btm_layer%map(3,n)) = 0.0_WP
      end do
 
      ! Reset interface-based quantities due to VF changes
      call this%vf%advect_interface(0.0_WP,this%fs%U,this%fs%V,this%fs%W)
      call this%vf%remove_flotsams()
      call this%vf%remove_thinstruct()
      call this%vf%sync_and_clean_barycenters()
      call this%vf%update_band()
      call this%vf%build_interface()
      call this%vf%polygonalize_interface()
      call this%vf%distance_from_polygon()
      call this%vf%subcell_vol()
      call this%vf%get_curvature()
      call this%vf%reset_moments()
 
   end subroutine    

   subroutine turb_stats(this,m_thick,g_thick,l_thick,Re_m,turb_lengthV,turb_lengthP,TKE,EPS_min,eta,EPS_check,KE_old)
      use mathtools,  only: Pi
      use mpi_f08,    only: MPI_ALLREDUCE,MPI_SUM 
      use parallel,   only: MPI_REAL_WP 
      use param,      only: param_read
      use string,     only: str_medium
      implicit none
      class(ml), intent(inout) :: this
      integer :: i,j,k,ierr
      real(WP), dimension(:), allocatable :: rho_pavg,vol_pavg,tmp_quick
      real(WP), dimension(:), allocatable :: rhoU_pavg,rhoV_pavg,rhoW_pavg,P_pavg,Tmpr,rhoE_pavg
      real(WP), dimension(:), allocatable :: vf_pavg,curv_pavg
      real(WP), dimension(:), allocatable :: Grho_pavg, GrhoU_pavg, GrhoE_pavg
      real(WP), dimension(:), allocatable :: Lrho_pavg, LrhoU_pavg, LrhoE_pavg
      real(WP), dimension(:), allocatable :: P_work,P_dil,Dissp,Prod,Transp,STen,Flux,KE,dKdt
      real(WP), dimension(:), intent(inout) :: KE_old
      real(WP) :: nu_avg,turb_lengthV,turb_lengthP
      real(WP) :: m_thick,g_thick,l_thick,Re_m
      real(WP) :: TKE,EPS_min,eta,EPS_check
      ! Allocate 
      allocate(rho_pavg(this%fs%cfg%jmin:this%fs%cfg%jmax));   rho_pavg=0.0_WP
      allocate(vol_pavg(this%fs%cfg%jmin:this%fs%cfg%jmax));   vol_pavg=0.0_WP
      allocate(rhoU_pavg(this%fs%cfg%jmin:this%fs%cfg%jmax));  rhoU_pavg=0.0_WP
      allocate(rhoV_pavg(this%fs%cfg%jmin:this%fs%cfg%jmax));  rhoV_pavg=0.0_WP
      allocate(rhoW_pavg(this%fs%cfg%jmin:this%fs%cfg%jmax));  rhoW_pavg=0.0_WP
      allocate(rhoE_pavg(this%fs%cfg%jmin:this%fs%cfg%jmax));  rhoE_pavg=0.0_WP
      allocate(P_pavg(this%fs%cfg%jmin:this%fs%cfg%jmax));     P_pavg=0.0_WP
      allocate(vf_pavg(this%fs%cfg%jmin:this%fs%cfg%jmax));    vf_pavg=0.0_WP
      allocate(curv_pavg(this%fs%cfg%jmin:this%fs%cfg%jmax));  curv_pavg=0.0_WP
      allocate(Grho_pavg(this%fs%cfg%jmin:this%fs%cfg%jmax));  Grho_pavg=0.0_WP
      allocate(Lrho_pavg(this%fs%cfg%jmin:this%fs%cfg%jmax));  Lrho_pavg=0.0_WP
      allocate(GrhoU_pavg(this%fs%cfg%jmin:this%fs%cfg%jmax)); GrhoU_pavg=0.0_WP
      allocate(LrhoU_pavg(this%fs%cfg%jmin:this%fs%cfg%jmax)); LrhoU_pavg=0.0_WP
      allocate(GrhoE_pavg(this%fs%cfg%jmin:this%fs%cfg%jmax)); GrhoE_pavg=0.0_WP
      allocate(LrhoE_pavg(this%fs%cfg%jmin:this%fs%cfg%jmax)); LrhoE_pavg=0.0_WP
      allocate(Tmpr(this%fs%cfg%jmin:this%fs%cfg%jmax));       Tmpr=0.0_WP
      allocate(dKdt(this%fs%cfg%jmin:this%fs%cfg%jmax));       dKdt=0.0_WP
      allocate(KE(this%fs%cfg%jmin:this%fs%cfg%jmax));         KE=0.0_WP
      allocate(P_work(this%fs%cfg%jmin:this%fs%cfg%jmax));     P_work=0.0_WP
      allocate(P_dil(this%fs%cfg%jmin:this%fs%cfg%jmax));      P_dil=0.0_WP
      allocate(Dissp(this%fs%cfg%jmin:this%fs%cfg%jmax));      Dissp=0.0_WP
      allocate(Prod(this%fs%cfg%jmin:this%fs%cfg%jmax));       Prod=0.0_WP
      allocate(Transp(this%fs%cfg%jmin:this%fs%cfg%jmax));     Transp=0.0_WP
      allocate(Flux(this%fs%cfg%jmin:this%fs%cfg%jmax));       Flux=0.0_WP
      allocate(STen(this%fs%cfg%jmin:this%fs%cfg%jmax));       STen=0.0_WP
      allocate(tmp_quick(this%fs%cfg%jmin:this%fs%cfg%jmax));  tmp_quick=0.0_WP
      nu_avg=0.0_WP;m_thick=0.0_WP;g_thick=0.0_WP;l_thick=0.0_WP
      ! Calculate plane averages
      call plane_avg_init(this,rho_pavg,vol_pavg,rhoU_pavg,rhoV_pavg,rhoW_pavg,P_pavg,vf_pavg,Grho_pavg,Lrho_pavg,GrhoU_pavg,LrhoU_pavg,nu_avg,curv_pavg)
      ! Calculate momentum thickness
      call momentum_thickness(this,m_thick,g_thick,l_thick,TKE,Re_m,rho_pavg,rhoU_pavg,rhoV_pavg,rhoW_pavg,nu_avg,Grho_pavg,Lrho_pavg,GrhoU_pavg,LrhoU_pavg,vf_pavg)
      ! Calculate TKE budget
      call TKE_budget(this,KE,P_work,P_dil,Dissp,Prod,Transp,Flux,STen,TKE,EPS_min,eta,EPS_check,rho_pavg,rhoU_pavg,rhoV_pavg,rhoW_pavg,P_pavg,vol_pavg,curv_pavg)
      dKdt = (KE - KE_old)/this%time%dt
      eta = (nu_avg**3.0_WP/EPS_min)**0.25_WP
      ! Calculate a few more averages
       ! Calculate spatial averages
      do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
         do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
            do i=this%fs%cfg%imin_,this%fs%cfg%imax_
               ! Total averaged calculations
               Tmpr(j)      = Tmpr(j) + this%fs%Tmptr(i,j,k)*this%fs%cfg%vol(i,j,k)
               rhoE_pavg(j) = rhoE_pavg(j) + ((1.0_WP-this%vf%VF(i,j,k))*this%fs%GrhoE(i,j,k) + this%vf%VF(i,j,k)*this%fs%LrhoE(i,j,k))*this%fs%cfg%vol(i,j,k)
               ! Volume fraction specific calculations
               GrhoE_pavg(j) = GrhoE_pavg(j) + this%fs%GrhoE(i,j,k)*this%fs%cfg%vol(i,j,k)
               LrhoE_pavg(j) = LrhoE_pavg(j) + this%fs%LrhoE(i,j,k)*this%fs%cfg%vol(i,j,k)
            end do
         end do
      end do
      ! All-reduce the data
      call MPI_ALLREDUCE(Tmpr,tmp_quick,this%fs%cfg%ny,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
      Tmpr = tmp_quick/vol_pavg
      call MPI_ALLREDUCE(rhoE_pavg,tmp_quick,this%fs%cfg%ny,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
      rhoE_pavg = tmp_quick/vol_pavg
      call MPI_ALLREDUCE(GrhoE_pavg,tmp_quick,this%fs%cfg%ny,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
      GrhoE_pavg = tmp_quick/vol_pavg
      call MPI_ALLREDUCE(LrhoE_pavg,tmp_quick,this%fs%cfg%ny,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
      LrhoE_pavg = tmp_quick/vol_pavg
      ! Post Process Data
      if (this%fs%cfg%amRoot) then
         ! Reynolds Stresses
         filename='PostProc_'
         write(timestamp,'(f6.0)') this%time%t
         open(newunit=junit,file='PostProc/'//trim(adjustl(filename))//trim(adjustl(timestamp)),form='formatted',status='replace',access='stream',iostat=ierr)
         write(junit,'(a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x)') &
            'Height','dKdt','Transp','Prod','Dissp','P Dil','P Work','Flux','STension','rho','U_fav','V_fav','W_fav','KE','Grho','Lrho','VF','GrhoU','LrhoU','GrhoE','LrhoE','rhoE','Tmptr'
         do j=this%fs%cfg%jmin,this%fs%cfg%jmax
            write(junit,'(es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x)') & 
               this%fs%cfg%ym(j),dKdt(j),Transp(j),Prod(j),Dissp(j),P_dil(j),P_work(j),Flux(j),STen(j),rho_pavg(j),rhoU_pavg(j)/rho_pavg(j),rhoV_pavg(j)/rho_pavg(j),rhoW_pavg(j)/rho_pavg(j),KE(j),Grho_pavg(j),Lrho_pavg(j),vf_pavg(j),GrhoU_pavg(j),LrhoU_pavg(j),GrhoE_pavg(j),LrhoE_pavg(j),rhoE_pavg(j),Tmpr(j)
         end do
         close(junit)
      end if
      KE_old = KE
      deallocate(rho_pavg,vol_pavg,rhoU_pavg,rhoV_pavg,rhoW_pavg,P_pavg,Tmpr,rhoE_pavg,tmp_quick)
      deallocate(vf_pavg,curv_pavg,Grho_pavg,GrhoU_pavg,GrhoE_pavg,Lrho_pavg,LrhoU_pavg,LrhoE_pavg)
      deallocate(P_work,P_dil,Dissp,Prod,Transp,Flux,STen,dKdt,KE)

      contains

         subroutine plane_avg_init(this,rho_pavg,vol_pavg,rhoU_pavg,rhoV_pavg,rhoW_pavg,P_pavg,vf_pavg,Grho_pavg,Lrho_pavg,GrhoU_pavg,LrhoU_pavg,nu_avg,curv_pavg)
            implicit none
            class(ml), intent(inout) :: this
            integer :: i,j,k,ierr,iunit
            real(WP), dimension(:), allocatable :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10,tmp11,tmp12
            real(WP) :: tmp13
            real(WP), dimension(:), intent(inout) :: rho_pavg,vol_pavg
            real(WP), dimension(:), intent(inout) :: rhoU_pavg,rhoV_pavg,rhoW_pavg,P_pavg
            real(WP), dimension(:), intent(inout) :: vf_pavg,curv_pavg
            real(WP), dimension(:), intent(inout) :: Grho_pavg, GrhoU_pavg
            real(WP), dimension(:), intent(inout) :: Lrho_pavg, LrhoU_pavg
            real(WP), intent(inout) :: nu_avg
            
            ! Allocate temporary arrays for MPI
            allocate(tmp1(this%fs%cfg%jmin:this%fs%cfg%jmax));   tmp1=0.0_WP
            allocate(tmp2(this%fs%cfg%jmin:this%fs%cfg%jmax));   tmp2=0.0_WP
            allocate(tmp3(this%fs%cfg%jmin:this%fs%cfg%jmax));   tmp3=0.0_WP
            allocate(tmp4(this%fs%cfg%jmin:this%fs%cfg%jmax));   tmp4=0.0_WP
            allocate(tmp5(this%fs%cfg%jmin:this%fs%cfg%jmax));   tmp5=0.0_WP
            allocate(tmp6(this%fs%cfg%jmin:this%fs%cfg%jmax));   tmp6=0.0_WP
            allocate(tmp7(this%fs%cfg%jmin:this%fs%cfg%jmax));   tmp7=0.0_WP
            allocate(tmp8(this%fs%cfg%jmin:this%fs%cfg%jmax));   tmp8=0.0_WP
            allocate(tmp9(this%fs%cfg%jmin:this%fs%cfg%jmax));   tmp9=0.0_WP
            allocate(tmp10(this%fs%cfg%jmin:this%fs%cfg%jmax));  tmp10=0.0_WP
            allocate(tmp11(this%fs%cfg%jmin:this%fs%cfg%jmax));  tmp11=0.0_WP
            allocate(tmp12(this%fs%cfg%jmin:this%fs%cfg%jmax));  tmp12=0.0_WP
            ! Calculate spatial averages
            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                     ! Total averaged calculations
                     rho_pavg(j)  = rho_pavg(j) + this%fs%RHO(i,j,k)*this%fs%cfg%vol(i,j,k)
                     vol_pavg(j)  = vol_pavg(j) + this%fs%cfg%vol(i,j,k)
                     rhoU_pavg(j) = rhoU_pavg(j) + this%fs%RHO(i,j,k)*this%fs%Ui(i,j,k)*this%fs%cfg%vol(i,j,k)
                     rhoV_pavg(j) = rhoV_pavg(j) + this%fs%RHO(i,j,k)*this%fs%Vi(i,j,k)*this%fs%cfg%vol(i,j,k)
                     rhoW_pavg(j) = rhoW_pavg(j) + this%fs%RHO(i,j,k)*this%fs%Wi(i,j,k)*this%fs%cfg%vol(i,j,k)
                     P_pavg(j)    = P_pavg(j) + this%fs%P(i,j,k)*this%fs%cfg%vol(i,j,k)
                     nu_avg       = nu_avg + this%fs%visc(i,j,k)*this%fs%cfg%vol(i,j,k)/this%fs%RHO(i,j,k)
                     ! Volume fraction specific calculations
                     vf_pavg(j)    = vf_pavg(j) + this%vf%VF(i,j,k)*this%fs%cfg%vol(i,j,k)
                     curv_pavg(j)  = curv_pavg(j) + this%vf%curv(i,j,k)*this%fs%cfg%vol(i,j,k)
                     Grho_pavg(j)  = Grho_pavg(j) + (1.0_WP - this%vf%VF(i,j,k))*this%fs%Grho(i,j,k)*this%fs%cfg%vol(i,j,k)
                     Lrho_pavg(j)  = Lrho_pavg(j) + this%vf%VF(i,j,k)*this%fs%Lrho(i,j,k)*this%fs%cfg%vol(i,j,k)
                     GrhoU_pavg(j) = GrhoU_pavg(j) + (1.0_WP - this%vf%VF(i,j,k))*this%fs%Grho(i,j,k)*this%fs%Ui(i,j,k)*this%fs%cfg%vol(i,j,k)
                     LrhoU_pavg(j) = LrhoU_pavg(j) + this%vf%VF(i,j,k)*this%fs%Lrho(i,j,k)*this%fs%Ui(i,j,k)*this%fs%cfg%vol(i,j,k)
                  end do
               end do
            end do
            ! All-reduce the data
            call MPI_ALLREDUCE(rho_pavg,tmp1,this%fs%cfg%ny,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
            call MPI_ALLREDUCE(vol_pavg,tmp2,this%fs%cfg%ny,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
            call MPI_ALLREDUCE(rhoU_pavg,tmp3,this%fs%cfg%ny,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
            call MPI_ALLREDUCE(rhoV_pavg,tmp4,this%fs%cfg%ny,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
            call MPI_ALLREDUCE(rhoW_pavg,tmp5,this%fs%cfg%ny,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
            call MPI_ALLREDUCE(P_pavg,tmp6,this%fs%cfg%ny,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
            call MPI_ALLREDUCE(nu_avg,tmp13,1,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
            call MPI_ALLREDUCE(vf_pavg,tmp7,this%fs%cfg%ny,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
            call MPI_ALLREDUCE(curv_pavg,tmp8,this%fs%cfg%ny,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
            call MPI_ALLREDUCE(Grho_pavg,tmp9,this%fs%cfg%ny,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
            call MPI_ALLREDUCE(Lrho_pavg,tmp10,this%fs%cfg%ny,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
            call MPI_ALLREDUCE(GrhoU_pavg,tmp11,this%fs%cfg%ny,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
            call MPI_ALLREDUCE(LrhoU_pavg,tmp12,this%fs%cfg%ny,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
            ! Volume-average for each plane
            vol_pavg   = tmp2
            rho_pavg   = tmp1/vol_pavg
            rhoU_pavg  = tmp3/vol_pavg
            rhoV_pavg  = tmp4/vol_pavg
            rhoW_pavg  = tmp5/vol_pavg
            P_pavg     = tmp6/vol_pavg
            nu_avg     = tmp13/this%fs%cfg%vol_total
            vf_pavg    = tmp7/vol_pavg
            curv_pavg  = tmp8/vol_pavg
            Grho_pavg  = tmp9/vol_pavg
            Lrho_pavg  = tmp10/vol_pavg
            GrhoU_pavg = tmp11/vol_pavg
            LrhoU_pavg = tmp12/vol_pavg
            ! ! Phase-specific Favre averaging (to prevent division by zero)
            ! do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
            !    if ((1-vf_pavg(j)) .ne. 0.0_WP) then
            !       Grho_pavg(j) = Grho_pavg(j)/(1.0_WP-vf_pavg(j))
            !    end if
            !    if (vf_pavg(j) .ne. 0.0_WP) then
            !       Lrho_pavg(j) = Lrho_pavg(j)/vf_pavg(j)
            !    end if
            !    if (Grho_pavg(j) .ne. 0.0_WP) then
            !       GrhoU_pavg(j) = GrhoU_pavg(j)/(1.0_WP - vf_pavg(j))
            !    end if
            !    if (Lrho_pavg(j) .ne. 0.0_WP) then
            !       LrhoU_pavg(j) = LrhoU_pavg(j)/(vf_pavg(j))
            !    end if
            ! end do
            ! Deallocate temporary arrays
            deallocate(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10,tmp11,tmp12)
         end subroutine plane_avg_init

         subroutine momentum_thickness(this,m_thick,g_thick,l_thick,TKE,Re_m,rho_pavg,rhoU_pavg,rhoV_pavg,rhoW_pavg,nu_avg,Grho_pavg,Lrho_pavg,GrhoU_pavg,LrhoU_pavg,vf_pavg)
            implicit none
            class(ml), intent(inout) :: this
            integer :: j
            real(WP) :: m_thick,TKE,myTKE,Re_m,nu_avg
            real(WP), intent(out) :: g_thick, l_thick
            real(WP), dimension(:), intent(in) :: rho_pavg,rhoU_pavg,rhoV_pavg,rhoW_pavg,Grho_pavg,Lrho_pavg,GrhoU_pavg,LrhoU_pavg,vf_pavg
            
            ! ! Calculate Momentum thicknesses (KEEP FOR NOW)
            ! do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
            !    if (Grho_pavg(j) .ne. 0.0_WP) then
            !       g_thick = g_thick + (1.0_WP/(Grho*(1.0_WP+r_vel)**2.0_WP))*((1.0_WP-vf_pavg(j))*Grho_pavg(j)*(1.0_WP - GrhoU_pavg(j)/Grho_pavg(j))*(GrhoU_pavg(j)/Grho_pavg(j) + r_vel))*this%vf%cfg%dy(j)
            !    end if
            !    if (Lrho_pavg(j) .ne. 0.0_WP) then
            !       l_thick = l_thick + (1.0_WP/(Lrho*(1.0_WP+r_vel)**2.0_WP))*(vf_pavg(j)*Lrho_pavg(j)*(1.0_WP - LrhoU_pavg(j)/Lrho_pavg(j))*(LrhoU_pavg(j)/Lrho_pavg(j) + r_vel))*this%vf%cfg%dy(j)
            !    end if
            !    m_thick = m_thick + (1.0_WP/((1.0_WP)**2.0_WP))*(rho_pavg(j)*(0.5_WP - (rhoU_pavg(j)/rho_pavg(j)))*((rhoU_pavg(j)/rho_pavg(j)) + 0.5_WP))*this%fs%cfg%dy(j)
            ! end do
            ! Revised momentum thicknesses (from strictly splitting up the momentum thickness)
            do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
               if (Grho_pavg(j) .ne. 0.0_WP) then
                  g_thick = g_thick + (1.0_WP/(1.0_WP*(1.0_WP+r_vel)**2.0_WP))*(Grho_pavg(j)*(1.0_WP - rhoU_pavg(j)/rho_pavg(j))*(rhoU_pavg(j)/rho_pavg(j)))*this%vf%cfg%dy(j)
               end if
               if (Lrho_pavg(j) .ne. 0.0_WP) then
                  l_thick = l_thick + (1.0_WP/(1.0_WP*(1.0_WP+r_vel)**2.0_WP))*(Lrho_pavg(j)*(1.0_WP - rhoU_pavg(j)/rho_pavg(j))*(rhoU_pavg(j)/rho_pavg(j)))*this%vf%cfg%dy(j)
               end if
               m_thick = m_thick + (1.0_WP/(1.0_WP*(1.0_WP)**2.0_WP))*(rho_pavg(j)*(1.0_WP - (rhoU_pavg(j)/rho_pavg(j)))*((rhoU_pavg(j)/rho_pavg(j))))*this%fs%cfg%dy(j)
            end do

            ! Calculate current TKE
            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                     myTKE=myTKE+0.5_WP*((this%fs%Ui(i,j,k)-(rhoU_pavg(j)/rho_pavg(j)))**2+(this%fs%Vi(i,j,k)-(rhoV_pavg(j)/rho_pavg(j)))**2+(this%fs%Wi(i,j,k)-(rhoW_pavg(j)/rho_pavg(j)))**2)*this%fs%cfg%vol(i,j,k)
                  end do
               end do
            end do
            call MPI_ALLREDUCE(myTKE,TKE,1,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr); TKE=TKE/this%fs%cfg%vol_total

            ! Calculate momentum thickness Reynolds Number
            Re_m = m_thick*1.0_WP/(nu_avg+epsilon(1.0_WP))

         end subroutine momentum_thickness

         subroutine TKE_budget(this,KE,P_work,P_dil,Dissp,Prod,Transp,Flux,STen,TKE,EPS_min,eta,EPS_check,rho_pavg,rhoU_pavg,rhoV_pavg,rhoW_pavg,P_pavg,vol_pavg,curv_pavg)
            implicit none
            class(ml), intent(inout) :: this
            integer :: i,j,k
            real(WP) :: TKE,EPS_min,EPS_check,eta
            real(WP), dimension(:), intent(inout)   :: KE,P_work,P_dil,Dissp,Prod,Transp,Flux,STen
            real(WP), dimension(:), intent(in)      :: rho_pavg,rhoU_pavg,rhoV_pavg,rhoW_pavg,P_pavg,vol_pavg,curv_pavg
            real(WP), dimension(:,:,:), allocatable :: Transp_local,Prod_local,Dissp_local,Pdil_local,Pwork_local,UU_local,STen_local
            real(WP), dimension(:,:,:), allocatable :: Ui_fluc,Vi_fluc,Wi_fluc,P_fluc
            real(WP), dimension(:,:,:,:), allocatable :: Pgrad,Pgrad_fluc
            ! real(WP), dimension(:,:,:), allocatable :: Pgradx,Pgrady,Pgradz
            ! real(WP), dimension(:,:,:), allocatable :: dPdx_fluc,dPdy_fluc,dPdz_fluc
            real(WP), dimension(:,:,:,:,:), allocatable :: gradU,gradU_fluc
            real(WP), dimension(:,:,:,:), allocatable :: SR,SS,Re_stress,SS_fluc,dUU_local
            real(WP), dimension(:,:), allocatable :: Re_stress_pavg,Pgrad_pavg,Ufluc_pavg,dUU_pavg
            real(WP), dimension(:,:,:), allocatable :: gradU_fav
            real(WP), dimension(:), allocatable :: eps_int,UU_pavg

            ! Allocate Arrays
            allocate(Ui_fluc(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_));  Ui_fluc=0.0_WP
            allocate(Vi_fluc(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_));  Vi_fluc=0.0_WP
            allocate(Wi_fluc(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_));  Wi_fluc=0.0_WP
            allocate(P_fluc (this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_));  P_fluc=0.0_WP
            allocate(gradU(1:3,1:3,this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_));  gradU=0.0_WP
            allocate(gradU_fluc(1:3,1:3,this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_));  gradU_fluc=0.0_WP
            allocate(SR(1:6,this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_));   SR=0.0_WP
            allocate(SS(1:6,this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_));   SS=0.0_WP
            allocate(SS_fluc(1:6,this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_));   SS_fluc=0.0_WP
            allocate(Re_stress(1:6,this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_));   Re_stress=0.0_WP
            allocate(Pgrad(1:3,this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_)); Pgrad=0.0_WP
            allocate(Pgrad_fluc(1:3,this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_)); Pgrad_fluc=0.0_WP
            allocate(Transp_local(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_)); Transp_local=0.0_WP
            allocate(Prod_local(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_)); Prod_local=0.0_WP
            allocate(Dissp_local(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_)); Dissp_local=0.0_WP
            allocate(Pdil_local(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_)); Pdil_local=0.0_WP
            allocate(Pwork_local(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_)); Pwork_local=0.0_WP
            allocate(UU_local(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_)); UU_local=0.0_WP
            allocate(dUU_local(1:3,this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_)); dUU_local=0.0_WP
            allocate(STen_local(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_)); STen_local=0.0_WP
            allocate(gradU_fav (1:3,1:3,this%fs%cfg%jmin:this%fs%cfg%jmax));   gradU_fav=0.0_WP
            allocate(Ufluc_pavg(1:3,this%fs%cfg%jmin:this%fs%cfg%jmax));   Ufluc_pavg=0.0_WP
            allocate(Pgrad_pavg(1:3,this%fs%cfg%jmin:this%fs%cfg%jmax));   Pgrad_pavg=0.0_WP
            allocate(Re_stress_pavg(1:6,this%fs%cfg%jmin:this%fs%cfg%jmax));   Re_stress_pavg=0.0_WP
            allocate(eps_int(this%fs%cfg%jmin:this%fs%cfg%jmax));   eps_int=0.0_WP
            allocate(UU_pavg(this%fs%cfg%jmin:this%fs%cfg%jmax));   UU_pavg=0.0_WP
            allocate(dUU_pavg(1:3,this%fs%cfg%jmin:this%fs%cfg%jmax));   dUU_pavg=0.0_WP
            EPS_check=0.0_WP
            ! Calculate fluctuating quantities from plane favre-averaged velocities
            call fluctuating_quant(this,Ui_fluc,Vi_fluc,Wi_fluc,P_fluc,rho_pavg,rhoU_pavg,rhoV_pavg,rhoW_pavg,P_pavg,vol_pavg,curv_pavg)
            ! Calculate intermediate local quantities
            call this%fs%get_gradu(gradU)
            call this%fs%get_strainrate(SR)
            call get_shearstress(this,SS,SR)
            call get_gradUfluc(this,gradU_fluc,Ui_fluc,Vi_fluc,Wi_fluc)
            call get_gradPfluc(this,Pgrad_fluc(1,:,:,:),Pgrad_fluc(2,:,:,:),Pgrad_fluc(1,:,:,:),P_fluc)
            call this%fs%get_pgradm(this%fs%P,Pgrad(1,:,:,:),Pgrad(2,:,:,:),Pgrad(3,:,:,:))
            call get_Re_stress(this,Re_stress,rho_pavg,rhoU_pavg,rhoV_pavg,rhoW_pavg)
            call get_shearstressfluc(this,SS_fluc,SS,vol_pavg)
            ! Calculate local Transport
            ! call get_transport_local(this,Transp_local,Ui_fluc,Vi_fluc,Wi_fluc,P_fluc,SS_fluc,vol_pavg)
            call get_transport_local_test(this,Transp_local,Ui_fluc,Vi_fluc,Wi_fluc,P_fluc,SS,vol_pavg)
            ! Calculate local Dissipation
            call get_dissipation_local(this,Dissp_local,SS,gradU_fluc)
            ! Calculate local Pressure Dilatation
            call get_p_dilatation_local(this,Pdil_local,gradU_fluc,P_fluc)
            ! Calculate local momentum flux
            call get_momentum_local(this,UU_local,dUU_local,Ui_fluc,Vi_fluc,Wi_fluc)
            ! Calculate local surface tension
            call get_stension_local(this,STen_local,Ui_fluc,Vi_fluc,Wi_fluc)
            ! ! Get total TKE, eps, and eta
            ! do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
            !    do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
            !       do i=this%fs%cfg%imin_,this%fs%cfg%imax_
            !          myTKE=myTKE+0.5_WP*((this%fs%Ui(i,j,k)-(rhoU_pavg(j)/rho_pavg(j)))**2+(this%fs%Vi(i,j,k)-(rhoV_pavg(j)/rho_pavg(j)))**2+(this%fs%Wi(i,j,k)-(rhoW_pavg(j)/rho_pavg(j)))**2)*this%fs%cfg%vol(i,j,k)
            !       end do
            !    end do
            ! end do
            ! Plane average quantities for terms that make up production, pressure work, and momentum flux
            call get_gradU_fav(this,gradU_fav,gradU,vol_pavg,rho_pavg)
            call get_plane_avg(this,Ufluc_pavg(1,:),Ui_fluc(:,:,:),vol_pavg)
            call get_plane_avg(this,Ufluc_pavg(2,:),Vi_fluc(:,:,:),vol_pavg)
            call get_plane_avg(this,Ufluc_pavg(3,:),Wi_fluc(:,:,:),vol_pavg)
            call get_plane_avg(this,UU_pavg(:),UU_local(:,:,:),vol_pavg)
            do i=1,3
               call get_plane_avg(this,Pgrad_pavg(i,:),Pgrad(i,:,:,:),vol_pavg)
               call get_plane_avg(this,dUU_pavg(i,:),dUU_local(i,:,:,:),vol_pavg)
            end do
            do i=1,6
               call get_plane_avg(this,Re_stress_pavg(i,:),Re_stress(i,:,:,:),vol_pavg)
            end do
            ! Calculate final production, pressure work, and momentum flux
            P_work = Ufluc_pavg(1,:)*Pgrad_pavg(1,:) + Ufluc_pavg(2,:)*Pgrad_pavg(2,:) + Ufluc_pavg(3,:)*Pgrad_pavg(3,:)
            Prod = Re_stress_pavg(1,:)*gradU_fav(1,1,:) + Re_stress_pavg(2,:)*gradU_fav(2,2,:) + Re_stress_pavg(3,:)*gradU_fav(3,3,:) + &
            & Re_stress_pavg(4,:)*(gradU_fav(1,2,:)+gradU_fav(2,1,:))+Re_stress_pavg(5,:)*(gradU_fav(2,3,:)+gradU_fav(3,2,:))+Re_stress_pavg(6,:)*(gradU_fav(1,3,:)+gradU_fav(3,1,:))
            Flux = 0.5_WP*(dUU_pavg(1,:)*(rhoU_pavg(:)/rho_pavg(:)) + UU_pavg*gradU_fav(1,1,:)) + 0.5_WP*(dUU_pavg(2,:)*(rhoV_pavg(:)/rho_pavg(:)) + UU_pavg*gradU_fav(2,2,:)) + &
                 & 0.5_WP*(dUU_pavg(3,:)*(rhoW_pavg(:)/rho_pavg(:)) + UU_pavg*gradU_fav(3,3,:))
            ! Plane average rest of TKE budget terms using local values
            call get_plane_avg(this,Dissp,Dissp_local,vol_pavg)
            call get_plane_avg(this,P_dil,Pdil_local,vol_pavg)
            call get_plane_avg(this,Transp,Transp_local,vol_pavg)
            call get_plane_avg(this,STen,STen_local,vol_pavg)
            ! Quick method to get integrated dissipation rate for checking self similarity
            call get_plane_avg(this,eps_int,2.0_WP*this%fs%visc(:,:,:)*(SR(1,:,:,:)**2+SR(2,:,:,:)**2+SR(3,:,:,:)**2+2.0_WP*(SR(4,:,:,:)**2+SR(5,:,:,:)**2+SR(6,:,:,:)**2)),vol_pavg)
            do j=this%fs%cfg%jmin,this%fs%cfg%jmax
               EPS_check = EPS_check + eps_int(j)*this%fs%cfg%dy(j)
            end do
            ! Get Kolmogorov length scale
            call get_kolmogorov(this,eta,Dissp,rho_pavg,vol_pavg)
            call get_KE(this,KE,rho_pavg,rhoU_pavg,rhoV_pavg,rhoW_pavg,vol_pavg)
            ! call MPI_ALLREDUCE(myTKE,TKE,1,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr); TKE=TKE/this%fs%cfg%vol_total
            ! Deallocate Arrays
            deallocate(Ui_fluc,Vi_fluc,Wi_fluc,P_fluc,gradU_fluc,Pgrad_fluc,SS_fluc)
            deallocate(gradU_fav,Ufluc_pavg,Re_stress_pavg,Pgrad_pavg,eps_int,UU_pavg,dUU_pavg)
            deallocate(SR,SS,Re_stress,Pgrad,gradU)
            deallocate(Transp_local,Prod_local,Dissp_local,Pdil_local,Pwork_local,UU_local,dUU_local,STen_local)

         end subroutine

         subroutine fluctuating_quant(this,Ui_fluc,Vi_fluc,Wi_fluc,P_fluc,rho_pavg,rhoU_pavg,rhoV_pavg,rhoW_pavg,P_pavg,vol_pavg,curv_pavg)
            implicit none
            class(ml), intent(inout) :: this
            integer :: i,j,K
            real(WP), dimension(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_), intent(inout) :: Ui_fluc,Vi_fluc,Wi_fluc,P_fluc
            real(WP), dimension(this%fs%cfg%jmin:this%fs%cfg%jmax), intent(in)  :: rho_pavg,rhoU_pavg,rhoV_pavg,rhoW_pavg,P_pavg,vol_pavg,curv_pavg

            ! Calculate Fluctuating quantities from plane favre-averaged velocities
            do k=this%fs%cfg%kmino_,this%fs%cfg%kmaxo_
               do j=this%fs%cfg%jmino_,this%fs%cfg%jmaxo_
                  do i=this%fs%cfg%imino_,this%fs%cfg%imaxo_
                     Ui_fluc(i,j,k) = this%fs%Ui(i,j,k) - rhoU_pavg(j)/rho_pavg(j)
                     Vi_fluc(i,j,k) = this%fs%Vi(i,j,k) - rhoV_pavg(j)/rho_pavg(j)
                     Wi_fluc(i,j,k) = this%fs%Wi(i,j,k) - rhoW_pavg(j)/rho_pavg(j)
                     P_fluc(i,j,k)  = this%fs%P(i,j,k)  - P_pavg(j)
                  end do
               end do
            end do
            ! Communicate
            call this%cfg%sync(Ui_fluc)
            call this%cfg%sync(Vi_fluc)
            call this%cfg%sync(Wi_fluc)
            call this%cfg%sync(P_fluc)

         end subroutine fluctuating_quant

         subroutine get_gradUfluc(this,dUdx,Ui,Vi,Wi)
            implicit none
            class(ml), intent(inout) :: this
            integer :: i,j,k
            real(WP), dimension(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_), intent(in) :: Ui,Vi,Wi
            real(WP), dimension(:,:,:), allocatable :: U,V,W,Vf_x,Wf_y,Uf_z,Wf_x,Uf_y,Vf_z
            real(WP), dimension(1:3,1:3,this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_), intent(out) :: dUdx

            ! Allocate
            allocate(U(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_));     U=0.0_WP
            allocate(V(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_));     V=0.0_WP
            allocate(W(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_));     W=0.0_WP
            allocate(Vf_x(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_));  Vf_x=0.0_WP
            allocate(Wf_y(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_));  Wf_y=0.0_WP
            allocate(Uf_z(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_));  Uf_z=0.0_WP
            allocate(Wf_x(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_));  Wf_x=0.0_WP
            allocate(Uf_y(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_));  Uf_y=0.0_WP
            allocate(Vf_z(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_));  Vf_z=0.0_WP            
            ! Try interpolation with no density weighting
            call interp_no_rho(this,Ui,Vi,Wi,U,V,W)
            call interp_no_rho(this,Vi,Wi,Ui,Vf_x,Wf_y,Uf_z)
            call interp_no_rho(this,Wi,Ui,Vi,Wf_x,Uf_y,Vf_z)
            ! Calculate derivatives
            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                     dUdx(1,1,i,j,k) = sum(this%fs%divp_x(:,i,j,k) * U   (i:i+1,j,k))
                     dUdx(2,1,i,j,k) = sum(this%fs%divp_x(:,i,j,k) * Vf_x(i:i+1,j,k))
                     dUdx(3,1,i,j,k) = sum(this%fs%divp_x(:,i,j,k) * Wf_x(i:i+1,j,k))
                     dUdx(1,2,i,j,k) = sum(this%fs%divp_y(:,i,j,k) * Uf_y(i,j:j+1,k))
                     dUdx(2,2,i,j,k) = sum(this%fs%divp_y(:,i,j,k) * V   (i,j:j+1,k))
                     dUdx(3,2,i,j,k) = sum(this%fs%divp_y(:,i,j,k) * Wf_y(i,j:j+1,k))
                     dUdx(1,3,i,j,k) = sum(this%fs%divp_z(:,i,j,k) * Uf_z(i,j,k:k+1))
                     dUdx(2,3,i,j,k) = sum(this%fs%divp_z(:,i,j,k) * Vf_z(i,j,k:k+1))
                     dUdx(3,3,i,j,k) = sum(this%fs%divp_z(:,i,j,k) * W   (i,j,k:k+1))
                  end do
               end do
            end do
            ! Communicate
            call this%fs%cfg%sync(dUdx)
            ! Deallocate
            deallocate(U,V,W,Vf_x,Wf_y,Uf_z,Wf_x,Uf_y,Vf_z)

         end subroutine get_gradUfluc

         subroutine get_gradPfluc(this,dQ1dx,dQ1dy,dQ1dz,Q1)
            implicit none
            class(ml), intent(inout) :: this
            integer :: i,j,K
            real(WP), dimension(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_), intent(in) :: Q1
            real(WP), dimension(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_), intent(out) :: dQ1dx,dQ1dy,dQ1dz

            ! Calculate derivatives
            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                     dQ1dx(i,j,k) = sum(this%fs%divp_x(:,i,j,k) * Q1(i:i+1,j,k))
                     dQ1dy(i,j,k) = sum(this%fs%divp_y(:,i,j,k) * Q1(i,j:j+1,k))
                     dQ1dz(i,j,k) = sum(this%fs%divp_z(:,i,j,k) * Q1(i,j,k:k+1))
                  end do
               end do
            end do
            ! Communicate
            call this%fs%cfg%sync(dQ1dx)
            call this%fs%cfg%sync(dQ1dy)
            call this%fs%cfg%sync(dQ1dz)

         end subroutine get_gradPfluc

         subroutine get_shearstress(this,SS,SR)
            implicit none
            class(ml), intent(inout) :: this
            real(WP), dimension(1:6,this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_), intent(in) :: SR
            real(WP), dimension(1:6,this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_), intent(out) :: SS

            ! Calculate shear stress using the strain rate
            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                     SS(1,i,j,k) = 2.0_WP*this%fs%visc(i,j,k)*SR(1,i,j,k)
                     SS(2,i,j,k) = 2.0_WP*this%fs%visc(i,j,k)*SR(2,i,j,k)
                     SS(3,i,j,k) = 2.0_WP*this%fs%visc(i,j,k)*SR(3,i,j,k)
                     SS(4,i,j,k) = 2.0_WP*this%fs%visc(i,j,k)*SR(4,i,j,k)
                     SS(5,i,j,k) = 2.0_WP*this%fs%visc(i,j,k)*SR(5,i,j,k)
                     SS(6,i,j,k) = 2.0_WP*this%fs%visc(i,j,k)*SR(6,i,j,k)
                  end do
               end do
            end do
            call this%fs%cfg%sync(SS)
         end subroutine get_shearstress

         subroutine get_shearstressfluc(this,SS_fluc,SS,vol_pavg)
            implicit none
            class(ml), intent(inout) :: this
            real(WP), dimension(1:6,this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_), intent(out) :: SS_fluc
            real(WP), dimension(1:6,this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_), intent(in) :: SS
            real(WP), dimension(this%fs%cfg%jmin:this%fs%cfg%jmax), intent(in) :: vol_pavg
            real(WP), dimension(:,:), allocatable :: SS_pavg

            ! Allocate
            allocate(SS_pavg(1:6,this%fs%cfg%jmino:this%fs%cfg%jmax))
            ! Plane average SS
            call get_plane_avg(this,SS_pavg(1,:),SS(1,:,:,:),vol_pavg)
            call get_plane_avg(this,SS_pavg(2,:),SS(2,:,:,:),vol_pavg)
            call get_plane_avg(this,SS_pavg(3,:),SS(3,:,:,:),vol_pavg)
            call get_plane_avg(this,SS_pavg(4,:),SS(4,:,:,:),vol_pavg)
            call get_plane_avg(this,SS_pavg(5,:),SS(5,:,:,:),vol_pavg)
            call get_plane_avg(this,SS_pavg(6,:),SS(6,:,:,:),vol_pavg)
            ! Do SS' = SS - \bar{SS}
            ! Calculate shear stress using the strain rate
            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                     SS_fluc(1,i,j,k) = SS(1,i,j,k) - SS_pavg(1,j)
                     SS_fluc(2,i,j,k) = SS(2,i,j,k) - SS_pavg(2,j)
                     SS_fluc(3,i,j,k) = SS(3,i,j,k) - SS_pavg(3,j)
                     SS_fluc(4,i,j,k) = SS(4,i,j,k) - SS_pavg(4,j)
                     SS_fluc(5,i,j,k) = SS(5,i,j,k) - SS_pavg(5,j)
                     SS_fluc(6,i,j,k) = SS(6,i,j,k) - SS_pavg(6,j)
                  end do
               end do
            end do
            ! Communicate
            call this%cfg%sync(SS_fluc)
            ! Deallocate
            deallocate(SS_pavg)

         end subroutine get_shearstressfluc

         subroutine get_KE(this,KE,rho_pavg,rhoU_pavg,rhoV_pavg,rhoW_pavg,vol_pavg)
            implicit none
            class(ml), intent(inout) :: this
            integer :: i,j,k
            real(WP), dimension(this%fs%cfg%jmin:this%fs%cfg%jmax), intent(in)  :: rho_pavg,rhoU_pavg,rhoV_pavg,rhoW_pavg,vol_pavg
            real(WP), dimension(this%fs%cfg%jmin:this%fs%cfg%jmax), intent(out) :: KE
            real(WP), dimension(:), allocatable :: my_KE

            allocate(my_KE(this%fs%cfg%jmin:this%fs%cfg%jmax)); my_KE = 0.0_WP

            ! Get total TKE, eps, and eta
            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                    my_KE(j) = my_KE(j) + 0.5_WP*this%fs%RHO(i,j,k)*((this%fs%Ui(i,j,k)-(rhoU_pavg(j)/rho_pavg(j)))**2.0_WP+(this%fs%Vi(i,j,k)-(rhoV_pavg(j)/rho_pavg(j)))**2.0_WP+(this%fs%Wi(i,j,k)-(rhoW_pavg(j)/rho_pavg(j)))**2.0_WP)*this%fs%cfg%vol(i,j,k)
                  end do
               end do
            end do
            call MPI_ALLREDUCE(my_KE,KE,this%fs%cfg%ny,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr); 
            KE = KE/(vol_pavg)

         end subroutine get_KE

         subroutine get_Re_stress(this,Re_stress,rho_pavg,rhoU_pavg,rhoV_pavg,rhoW_pavg)
            implicit none
            class(ml), intent(inout) :: this
            integer :: i,j,k
            real(WP), dimension(this%fs%cfg%jmin:this%fs%cfg%jmax), intent(in) :: rho_pavg,rhoU_pavg,rhoV_pavg,rhoW_pavg
            real(WP), dimension(1:6,this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_), intent(out) :: Re_stress
            
            ! Calculate plane-averages
            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                     ! Local Reynolds Stresses
                     Re_stress(1,i,j,k) = this%fs%RHO(i,j,k)*((this%fs%Ui(i,j,k) - rhoU_pavg(j)/rho_pavg(j))**2.0_WP)
                     Re_stress(2,i,j,k) = this%fs%RHO(i,j,k)*((this%fs%Vi(i,j,k) - rhoV_pavg(j)/rho_pavg(j))**2.0_WP)
                     Re_stress(3,i,j,k) = this%fs%RHO(i,j,k)*((this%fs%Wi(i,j,k) - rhoW_pavg(j)/rho_pavg(j))**2.0_WP)
                     Re_stress(4,i,j,k) = this%fs%RHO(i,j,k)*(this%fs%Ui(i,j,k) - rhoU_pavg(j)/rho_pavg(j))*(this%fs%Vi(i,j,k) - rhoV_pavg(j)/rho_pavg(j))
                     Re_stress(5,i,j,k) = this%fs%RHO(i,j,k)*(this%fs%Ui(i,j,k) - rhoU_pavg(j)/rho_pavg(j))*(this%fs%Wi(i,j,k) - rhoW_pavg(j)/rho_pavg(j))
                     Re_stress(6,i,j,k) = this%fs%RHO(i,j,k)*(this%fs%Vi(i,j,k) - rhoV_pavg(j)/rho_pavg(j))*(this%fs%Wi(i,j,k) - rhoW_pavg(j)/rho_pavg(j))
                  end do
               end do
            end do
            ! Communicate
            call this%cfg%sync(Re_stress)

         end subroutine get_Re_stress

         subroutine get_transport_local(this,Transp_local,Ui_fluc,Vi_fluc,Wi_fluc,P_fluc,SS_fluc,vol_pavg)
            implicit none
            class(ml), intent(inout) :: this
            integer :: i,j,k
            real(WP), dimension(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_), intent(out) :: Transp_local
            real(WP), dimension(this%fs%cfg%jmin:this%fs%cfg%jmax), intent(in) :: vol_pavg
            real(WP), dimension(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_), intent(in) :: Ui_fluc,Vi_fluc,Wi_fluc,P_fluc
            real(WP), dimension(1:6,this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_), intent(in) :: SS_fluc
            real(WP), dimension(:,:,:), allocatable :: dTrdx,dTrdy,dTrdz

            ! Allocation
            allocate(dTrdx(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_));  dTrdx=0.0_WP
            allocate(dTrdy(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_));  dTrdy=0.0_WP
            allocate(dTrdz(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_));  dTrdz=0.0_WP
            ! Intermediate calculations at cell centers
            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                     dTrdx(i,j,k) = SS_fluc(1,i,j,k)*Ui_fluc(i,j,k) + SS_fluc(4,i,j,k)*Vi_fluc(i,j,k) + SS_fluc(6,i,j,k)*Wi_fluc(i,j,k) &
                        - 0.5_WP*(this%fs%RHO(i,j,k)*Ui_fluc(i,j,k)**3.0_WP + this%fs%RHO(i,j,k)*Ui_fluc(i,j,k)*Vi_fluc(i,j,k)**2.0_WP &
                        + this%fs%RHO(i,j,k)*Ui_fluc(i,j,k)*Wi_fluc(i,j,k)**2.0_WP) - P_fluc(i,j,k)*Ui_fluc(i,j,k)
                     dTrdy(i,j,k) = SS_fluc(4,i,j,k)*Ui_fluc(i,j,k) + SS_fluc(2,i,j,k)*Vi_fluc(i,j,k) + SS_fluc(5,i,j,k)*Wi_fluc(i,j,k) &
                        - 0.5_WP*(this%fs%RHO(i,j,k)*Vi_fluc(i,j,k)*Ui_fluc(i,j,k)**2.0_WP + this%fs%RHO(i,j,k)*Vi_fluc(i,j,k)**3.0_WP &
                        + this%fs%RHO(i,j,k)*Vi_fluc(i,j,k)*Wi_fluc(i,j,k)**2.0_WP) - P_fluc(i,j,k)*Vi_fluc(i,j,k)
                     dTrdz(i,j,k) = SS_fluc(6,i,j,k)*Ui_fluc(i,j,k) + SS_fluc(5,i,j,k)*Vi_fluc(i,j,k) + SS_fluc(3,i,j,k)*Wi_fluc(i,j,k) & 
                        - 0.5_WP*(this%fs%RHO(i,j,k)*Wi_fluc(i,j,k)*Ui_fluc(i,j,k)**2.0_WP + this%fs%RHO(i,j,k)*Wi_fluc(i,j,k)*Vi_fluc(i,j,k)**2.0_WP &
                        + this%fs%RHO(i,j,k)*Wi_fluc(i,j,k)**3.0_WP) - P_fluc(i,j,k)*Wi_fluc(i,j,k)
                  end do
               end do
            end do
            ! Communicate
            call this%fs%cfg%sync(dTrdx)
            call this%fs%cfg%sync(dTrdy)
            call this%fs%cfg%sync(dTrdz)
            ! Interpolate to cell centers (this was added/changed!!!)
            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                     dTrdx(i,j,k) = sum(this%fs%itpi_x(:,i,j,k) * dTrdx(i-1:i,j,k))
                     dTrdy(i,j,k) = sum(this%fs%itpi_y(:,i,j,k) * dTrdy(i,j-1:j,k))
                     dTrdz(i,j,k) = sum(this%fs%itpi_z(:,i,j,k) * dTrdz(i,j,k-1:k))
                  end do
               end do
            end do
            ! Communicate
            call this%fs%cfg%sync(dTrdx)
            call this%fs%cfg%sync(dTrdy)
            call this%fs%cfg%sync(dTrdz)
            ! Calculate derivatives
            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                     ! dTrdx(i,j,k) = sum(this%fs%divp_x(:,i,j,k) * dTrdx(i:i+1,j,k))
                     ! dTrdy(i,j,k) = sum(this%fs%divp_y(:,i,j,k) * dTrdy(i,j:j+1,k))
                     ! dTrdz(i,j,k) = sum(this%fs%divp_z(:,i,j,k) * dTrdz(i,j,k:k+1))
                     Transp_local(i,j,k) = sum(this%fs%divp_x(:,i,j,k) * dTrdx(i:i+1,j,k)) + sum(this%fs%divp_y(:,i,j,k) * dTrdy(i,j:j+1,k)) + sum(this%fs%divp_z(:,i,j,k) * dTrdz(i,j,k:k+1))
                  end do
               end do
            end do
            ! Communicate
            call this%cfg%sync(Transp_local)
            ! Deallocate
            deallocate(dTrdx,dTrdy,dTrdz)

         end subroutine get_transport_local

         subroutine get_transport_local_test(this,Transp_local,Ui_fluc,Vi_fluc,Wi_fluc,P_fluc,SS,vol_pavg)
            implicit none
            class(ml), intent(inout) :: this
            integer :: i,j,k
            real(WP), dimension(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_), intent(out) :: Transp_local
            real(WP), dimension(this%fs%cfg%jmin:this%fs%cfg%jmax), intent(in) :: vol_pavg
            real(WP), dimension(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_), intent(in) :: Ui_fluc,Vi_fluc,Wi_fluc,P_fluc
            real(WP), dimension(1:6,this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_), intent(in) :: SS
            real(WP), dimension(:,:,:), allocatable :: dTrdx,dTrdy,dTrdz

            ! Allocation
            allocate(dTrdx(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_));  dTrdx=0.0_WP
            allocate(dTrdy(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_));  dTrdy=0.0_WP
            allocate(dTrdz(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_));  dTrdz=0.0_WP
            ! Intermediate calculations at cell centers
            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                     dTrdx(i,j,k) = SS(1,i,j,k)*Ui_fluc(i,j,k) + SS(4,i,j,k)*Vi_fluc(i,j,k) + SS(6,i,j,k)*Wi_fluc(i,j,k) &
                        - 0.5_WP*(this%fs%RHO(i,j,k)*Ui_fluc(i,j,k)**3.0_WP + this%fs%RHO(i,j,k)*Ui_fluc(i,j,k)*Vi_fluc(i,j,k)**2.0_WP &
                        + this%fs%RHO(i,j,k)*Ui_fluc(i,j,k)*Wi_fluc(i,j,k)**2.0_WP) - P_fluc(i,j,k)*Ui_fluc(i,j,k)
                     dTrdy(i,j,k) = SS(4,i,j,k)*Ui_fluc(i,j,k) + SS(2,i,j,k)*Vi_fluc(i,j,k) + SS(5,i,j,k)*Wi_fluc(i,j,k) &
                        - 0.5_WP*(this%fs%RHO(i,j,k)*Vi_fluc(i,j,k)*Ui_fluc(i,j,k)**2.0_WP + this%fs%RHO(i,j,k)*Vi_fluc(i,j,k)**3.0_WP &
                        + this%fs%RHO(i,j,k)*Vi_fluc(i,j,k)*Wi_fluc(i,j,k)**2.0_WP) - P_fluc(i,j,k)*Vi_fluc(i,j,k)
                     dTrdz(i,j,k) = SS(6,i,j,k)*Ui_fluc(i,j,k) + SS(5,i,j,k)*Vi_fluc(i,j,k) + SS(3,i,j,k)*Wi_fluc(i,j,k) & 
                        - 0.5_WP*(this%fs%RHO(i,j,k)*Wi_fluc(i,j,k)*Ui_fluc(i,j,k)**2.0_WP + this%fs%RHO(i,j,k)*Wi_fluc(i,j,k)*Vi_fluc(i,j,k)**2.0_WP &
                        + this%fs%RHO(i,j,k)*Wi_fluc(i,j,k)**3.0_WP) - P_fluc(i,j,k)*Wi_fluc(i,j,k)
                  end do
               end do
            end do
            ! Communicate
            call this%fs%cfg%sync(dTrdx)
            call this%fs%cfg%sync(dTrdy)
            call this%fs%cfg%sync(dTrdz)
            ! Interpolate to cell centers (this was added/changed!!!)
            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                     dTrdx(i,j,k) = sum(this%fs%itpi_x(:,i,j,k) * dTrdx(i-1:i,j,k))
                     dTrdy(i,j,k) = sum(this%fs%itpi_y(:,i,j,k) * dTrdy(i,j-1:j,k))
                     dTrdz(i,j,k) = sum(this%fs%itpi_z(:,i,j,k) * dTrdz(i,j,k-1:k))
                  end do
               end do
            end do
            ! Communicate
            call this%fs%cfg%sync(dTrdx)
            call this%fs%cfg%sync(dTrdy)
            call this%fs%cfg%sync(dTrdz)
            ! Calculate derivatives
            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                     ! dTrdx(i,j,k) = sum(this%fs%divp_x(:,i,j,k) * dTrdx(i:i+1,j,k))
                     ! dTrdy(i,j,k) = sum(this%fs%divp_y(:,i,j,k) * dTrdy(i,j:j+1,k))
                     ! dTrdz(i,j,k) = sum(this%fs%divp_z(:,i,j,k) * dTrdz(i,j,k:k+1))
                     Transp_local(i,j,k) = sum(this%fs%divp_x(:,i,j,k) * dTrdx(i:i+1,j,k)) + sum(this%fs%divp_y(:,i,j,k) * dTrdy(i,j:j+1,k)) + sum(this%fs%divp_z(:,i,j,k) * dTrdz(i,j,k:k+1))
                  end do
               end do
            end do
            ! Communicate
            call this%cfg%sync(Transp_local)
            ! Deallocate
            deallocate(dTrdx,dTrdy,dTrdz)

         end subroutine get_transport_local_test

         subroutine get_dissipation_local(this,Dissp_local,SS,gradU_fluc)
            implicit none
            class(ml), intent(inout) :: this
            integer :: i,j,K
            real(WP), dimension(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_), intent(out) :: Dissp_local
            real(WP), dimension(1:6,this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_), intent(in) :: SS
            real(WP), dimension(1:3,1:3,this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_), intent(in) :: gradU_fluc

            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                     Dissp_local(i,j,k) = SS(1,i,j,k)*gradU_fluc(1,1,i,j,k) + SS(2,i,j,k)*gradU_fluc(2,2,i,j,k) + SS(3,i,j,k)*gradU_fluc(3,3,i,j,k) &
                                      & + SS(4,i,j,k)*(gradU_fluc(2,1,i,j,k)+gradU_fluc(1,2,i,j,k)) + SS(5,i,j,k)*(gradU_fluc(3,2,i,j,k)+gradU_fluc(2,3,i,j,k)) &
                                      & + SS(6,i,j,k)*(gradU_fluc(3,1,i,j,k)+gradU_fluc(1,3,i,j,k))
                  end do
               end do
            end do
            ! Communicate
            call this%cfg%sync(Dissp_local)

         end subroutine get_dissipation_local

         subroutine get_p_dilatation_local(this,Pdil_local,gradU_fluc,P_fluc)
            implicit none
            class(ml), intent(inout) :: this
            integer :: i,j,K
            real(WP), dimension(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_), intent(out) :: Pdil_local
            real(WP), dimension(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_), intent(in) :: P_fluc
            real(WP), dimension(1:3,1:3,this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_), intent(in) :: gradU_fluc

            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                     Pdil_local(i,j,k) = P_fluc(i,j,k)*(gradU_fluc(1,1,i,j,k)+gradU_fluc(2,2,i,j,k)+gradU_fluc(3,3,i,j,k))
                  end do
               end do
            end do
            ! Communicate
            call this%cfg%sync(Pdil_local)

         end subroutine get_p_dilatation_local

         subroutine get_momentum_local(this,UU_local,dUU_local,Ui_fluc,Vi_fluc,Wi_fluc)
            implicit none
            class(ml), intent(inout) :: this
            integer :: i,j,k
            real(WP),dimension(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_),intent(out) :: UU_local
            real(WP),dimension(1:3,this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_),intent(out) :: dUU_local
            real(WP),dimension(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_),intent(in) :: Ui_fluc,Vi_fluc,Wi_fluc
            real(WP),dimension(:,:,:),allocatable :: UU_x,UU_y,UU_z
            ! Allocate
            allocate(UU_x(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_)); UU_x=0.0_WP
            allocate(UU_y(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_)); UU_y=0.0_WP
            allocate(UU_z(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_)); UU_z=0.0_WP
            ! Calculate \rho*u"u" at cell centers
            do k=this%fs%cfg%kmino_,this%fs%cfg%kmaxo_
               do j=this%fs%cfg%jmino_,this%fs%cfg%jmaxo_
                  do i=this%fs%cfg%imino_,this%fs%cfg%imaxo_
                     UU_local(i,j,k) = this%fs%RHO(i,j,k)*(Ui_fluc(i,j,k)**2.0_WP + Vi_fluc(i,j,k)**2.0_WP + Wi_fluc(i,j,k)**2.0_WP)
                  end do
               end do
            end do
            ! Communicate
            call this%fs%cfg%sync(UU_local)
            ! Interpolate \rho*u"u" to cell faces (this was changed!!!)
            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                     UU_x(i,j,k) = sum(this%fs%itpi_x(:,i,j,k) * UU_local(i-1:i,j,k))
                     UU_y(i,j,k) = sum(this%fs%itpi_y(:,i,j,k) * UU_local(i,j-1:j,k))
                     UU_z(i,j,k) = sum(this%fs%itpi_z(:,i,j,k) * UU_local(i,j,k-1:k))
                  end do
               end do
            end do
            ! Communicate
            call this%fs%cfg%sync(UU_x)
            call this%fs%cfg%sync(UU_y)
            call this%fs%cfg%sync(UU_z)
            ! Calculate d\rho*u"u"/dx (this was changed from UU_local to UU_x,UU_y,UU_z)
            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                     dUU_local(1,i,j,k) = sum(this%fs%divp_x(:,i,j,k) * UU_x(i:i+1,j,k))
                     dUU_local(2,i,j,k) = sum(this%fs%divp_y(:,i,j,k) * UU_y(i,j:j+1,k))
                     dUU_local(3,i,j,k) = sum(this%fs%divp_x(:,i,j,k) * UU_z(i,j,k:k+1))
                  end do
               end do
            end do
            ! Communicate
            call this%fs%cfg%sync(dUU_local)
            ! Deallocate
            deallocate(UU_x,UU_y,UU_z)

         end subroutine get_momentum_local

         subroutine get_stension_local(this,ST_local,Ui_fluc,Vi_fluc,Wi_fluc)
            implicit none
            class(ml), intent(inout) :: this
            integer :: i,j,k
            real(WP),dimension(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_),intent(out) :: ST_local
            real(WP),dimension(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_),intent(in) :: Ui_fluc,Vi_fluc,Wi_fluc
            real(WP),dimension(:,:,:),allocatable :: dadx,dady,dadz
            ! Allocate
            allocate(dadx(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_)); dadx=0.0_WP
            allocate(dady(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_)); dady=0.0_WP
            allocate(dadz(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_)); dadz=0.0_WP
            ! Calculate d\alpha/dx_i
            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                     dadx(i,j,k) = (this%vf%VF(i+1,j,k) - this%vf%VF(i-1,j,k)) / (2.0_WP*this%fs%cfg%dx(i))
                     dady(i,j,k) = (this%vf%VF(i,j+1,k) - this%vf%VF(i,j-1,k)) / (2.0_WP*this%fs%cfg%dy(j))
                     dadz(i,j,k) = (this%vf%VF(i,j,k+1) - this%vf%VF(i,j,k-1)) / (2.0_WP*this%fs%cfg%dz(k))
                  end do
               end do
            end do
            ! Communicate
            call this%fs%cfg%sync(dadx)
            call this%fs%cfg%sync(dady)
            call this%fs%cfg%sync(dadz)
            ! Calculate \sigma*\kappa*u_fluc_i*d\alpha/dx_i
            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                     ST_local(i,j,k) = this%fs%sigma*this%vf%curv(i,j,k)*(Ui_fluc(i,j,k)*dadx(i,j,k) + Vi_fluc(i,j,k)*dady(i,j,k) + Wi_fluc(i,j,k)*dadz(i,j,k))
                  end do
               end do
            end do
            ! Communicate
            call this%fs%cfg%sync(ST_local)
            ! Deallocate
            deallocate(dadx,dady,dadz)

         end subroutine get_stension_local

         subroutine interp_no_rho(this,Ui,Vi,Wi,U,V,W)
            implicit none
            class(ml), intent(inout) :: this
            integer :: i,j,k
            real(WP), dimension(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_), intent(in) :: Ui,Vi,Wi
            real(WP), dimension(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_), intent(out) :: U,V,W

            ! Interpolate without density weighting (MAST does density weighting by default)
            do k=this%cfg%kmin_,this%cfg%kmax_+1
               do j=this%cfg%jmin_,this%cfg%jmax_+1
                  do i=this%cfg%imin_,this%cfg%imax_+1
                     ! X face
                     U(i,j,k) = sum(this%fs%itpi_x(:,i,j,k) * Ui(i-1:i,j,k))
                     ! Y face
                     V(i,j,k) = sum(this%fs%itpi_y(:,i,j,k) * Vi(i,j-1:j,k))
                     ! Z face
                     W(i,j,k) = sum(this%fs%itpi_z(:,i,j,k) * Wi(i,j,k-1:k))
                  end do
               end do
            end do

         end subroutine interp_no_rho

         subroutine get_gradU_fav(this,gradU_fav,gradU,vol_pavg,rho_pavg)
            implicit none
            class(ml), intent(inout) :: this
            integer :: i
            real(WP),dimension(1:3,1:3,this%fs%cfg%jmin:this%fs%cfg%jmax),intent(out) :: gradU_fav
            real(WP),dimension(1:3,1:3,this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_),intent(in) :: gradU
            real(WP),dimension(this%fs%cfg%jmin:this%fs%cfg%jmax),intent(in) :: vol_pavg,rho_pavg
            real(WP),dimension(:,:,:),allocatable :: gradU_fav_tmp

            ! Allocate temporary arrays for MPI
            allocate(gradU_fav_tmp(1:3,1:3,this%fs%cfg%jmin:this%fs%cfg%jmax));   gradU_fav_tmp=0.0_WP
            ! Calculate components for TKE terms
            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                     ! Sum Favre-averaged velocity derivatives
                     gradU_fav_tmp(1,1,j) = gradU_fav_tmp(1,1,j) + (this%fs%RHO(i,j,k) * gradU(1,1,i,j,k) * this%fs%cfg%vol(i,j,k)) 
                     gradU_fav_tmp(2,1,j) = gradU_fav_tmp(2,1,j) + (this%fs%RHO(i,j,k) * gradU(2,1,i,j,k) * this%fs%cfg%vol(i,j,k))
                     gradU_fav_tmp(3,1,j) = gradU_fav_tmp(3,1,j) + (this%fs%RHO(i,j,k) * gradU(3,1,i,j,k) * this%fs%cfg%vol(i,j,k))
                     gradU_fav_tmp(1,2,j) = gradU_fav_tmp(1,2,j) + (this%fs%RHO(i,j,k) * gradU(1,2,i,j,k) * this%fs%cfg%vol(i,j,k))
                     gradU_fav_tmp(2,2,j) = gradU_fav_tmp(2,2,j) + (this%fs%RHO(i,j,k) * gradU(2,2,i,j,k) * this%fs%cfg%vol(i,j,k))
                     gradU_fav_tmp(3,2,j) = gradU_fav_tmp(3,2,j) + (this%fs%RHO(i,j,k) * gradU(3,2,i,j,k) * this%fs%cfg%vol(i,j,k))
                     gradU_fav_tmp(1,3,j) = gradU_fav_tmp(1,3,j) + (this%fs%RHO(i,j,k) * gradU(1,3,i,j,k) * this%fs%cfg%vol(i,j,k))
                     gradU_fav_tmp(2,3,j) = gradU_fav_tmp(2,3,j) + (this%fs%RHO(i,j,k) * gradU(2,3,i,j,k) * this%fs%cfg%vol(i,j,k))
                     gradU_fav_tmp(3,3,j) = gradU_fav_tmp(3,3,j) + (this%fs%RHO(i,j,k) * gradU(3,3,i,j,k) * this%fs%cfg%vol(i,j,k))
                  end do
               end do
            end do
            ! All-reduce the data
            call MPI_ALLREDUCE(gradU_fav_tmp(1,1,:),gradU_fav(1,1,:),this%fs%cfg%ny,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
            call MPI_ALLREDUCE(gradU_fav_tmp(2,1,:),gradU_fav(2,1,:),this%fs%cfg%ny,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
            call MPI_ALLREDUCE(gradU_fav_tmp(3,1,:),gradU_fav(3,1,:),this%fs%cfg%ny,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
            call MPI_ALLREDUCE(gradU_fav_tmp(1,2,:),gradU_fav(1,2,:),this%fs%cfg%ny,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
            call MPI_ALLREDUCE(gradU_fav_tmp(2,2,:),gradU_fav(2,2,:),this%fs%cfg%ny,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
            call MPI_ALLREDUCE(gradU_fav_tmp(3,2,:),gradU_fav(3,2,:),this%fs%cfg%ny,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
            call MPI_ALLREDUCE(gradU_fav_tmp(1,3,:),gradU_fav(1,3,:),this%fs%cfg%ny,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
            call MPI_ALLREDUCE(gradU_fav_tmp(2,3,:),gradU_fav(2,3,:),this%fs%cfg%ny,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
            call MPI_ALLREDUCE(gradU_fav_tmp(3,3,:),gradU_fav(3,3,:),this%fs%cfg%ny,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
            ! Plane-average the data
            gradU_fav(1,1,:) = gradU_fav(1,1,:)/(vol_pavg*rho_pavg)
            gradU_fav(2,1,:) = gradU_fav(2,1,:)/(vol_pavg*rho_pavg)
            gradU_fav(3,1,:) = gradU_fav(3,1,:)/(vol_pavg*rho_pavg)
            gradU_fav(1,2,:) = gradU_fav(1,2,:)/(vol_pavg*rho_pavg)
            gradU_fav(2,2,:) = gradU_fav(2,2,:)/(vol_pavg*rho_pavg)
            gradU_fav(3,2,:) = gradU_fav(3,2,:)/(vol_pavg*rho_pavg)
            gradU_fav(1,3,:) = gradU_fav(1,3,:)/(vol_pavg*rho_pavg)
            gradU_fav(2,3,:) = gradU_fav(2,3,:)/(vol_pavg*rho_pavg)
            gradU_fav(3,3,:) = gradU_fav(3,3,:)/(vol_pavg*rho_pavg)
            ! Deallocate temporary arrays
            deallocate(gradU_fav_tmp)

         end subroutine get_gradU_fav

         subroutine get_plane_avg(this,Qa1,Q1,vol_pavg)
            implicit none
            class(ml), intent(inout) :: this
            integer :: i,j,k
            real(WP),dimension(this%fs%cfg%imino_:this%fs%cfg%imaxo_,this%fs%cfg%jmino_:this%fs%cfg%jmaxo_,this%fs%cfg%kmino_:this%fs%cfg%kmaxo_),intent(in) :: Q1
            real(WP),dimension(this%fs%cfg%jmin:this%fs%cfg%jmax),intent(in) :: vol_pavg
            real(WP),dimension(this%fs%cfg%jmin:this%fs%cfg%jmax),intent(out) :: Qa1
            real(WP),dimension(:), allocatable :: tmp

            ! Allocate temporary arrays for MPI
            allocate(tmp(this%fs%cfg%jmin:this%fs%cfg%jmax));   tmp=0.0_WP
            Qa1 = 0.0_WP;
            ! Calculate plane averages
            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                     ! Sum all entries
                     Qa1(j) = Qa1(j) + Q1(i,j,k)*this%fs%cfg%vol(i,j,k)
                  end do
               end do
            end do
            ! All-reduce the data
            call MPI_ALLREDUCE(Qa1,tmp,this%fs%cfg%ny,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
            ! Plane-average the data
            Qa1 = tmp/vol_pavg
            ! Deallocate temporary arrays
            deallocate(tmp)

         end subroutine get_plane_avg

         subroutine get_kolmogorov(this,eta,Dissp,rho_pavg,vol_pavg)
            implicit none
            class(ml), intent(inout) :: this
            integer :: i,j,k
            real(WP),intent(out) :: eta
            real(WP),dimension(this%fs%cfg%jmin:this%fs%cfg%jmax),intent(in) :: Dissp,rho_pavg,vol_pavg
            real(WP),dimension(:),allocatable :: eta_pavg,visc_pavg,visc_tmp

            ! Allocate
            allocate(eta_pavg(this%fs%cfg%jmin:this%fs%cfg%jmax)); eta_pavg=0.0_WP
            allocate(visc_pavg(this%fs%cfg%jmin:this%fs%cfg%jmax)); visc_pavg=0.0_WP
            allocate(visc_tmp(this%fs%cfg%jmin:this%fs%cfg%jmax)); visc_tmp=0.0_WP
            ! Calculate favre averaged viscosity
            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                     ! Sum all entries
                     visc_pavg(j) = visc_pavg(j) + this%fs%visc(i,j,k)*this%fs%cfg%vol(i,j,k)*this%fs%RHO(i,j,k)
                  end do
               end do
            end do
            ! All-reduce the data
            call MPI_ALLREDUCE(visc_pavg,visc_tmp,this%fs%cfg%ny,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
            ! Plane average
            visc_pavg = visc_tmp/(vol_pavg*rho_pavg)
            ! Calculate dissipation rate for each plane
            do j=this%fs%cfg%jmin+1,this%fs%cfg%jmax-1
               eta_pavg(j) = (((visc_pavg(j)/rho_pavg(j))**3.0_WP)/Dissp(j))**0.25_WP
               if (eta_pavg(j).eq.0.0_WP) eta_pavg(j) = 1.0_WP
            end do

            do j=this%fs%cfg%jmin,this%fs%cfg%jmax
               if (eta_pavg(j).eq.0.0_WP) eta_pavg(j) = 1.0_WP
            end do

            eta = MINVAL(eta_pavg)

            ! Deallocate
            deallocate(eta_pavg,visc_pavg,visc_tmp)
         
         end subroutine get_kolmogorov
      
   end subroutine turb_stats
 
 
end module ml_class
