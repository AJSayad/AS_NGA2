!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use mast_class,        only: mast
   use vfs_class,         only: vfs
   use matm_class,        only: matm
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   use hypre_str_class,   only: hypre_str
   use pardata_class,     only: pardata !AS restart
   use inputfile_class,   only: inputfile !AS restart 
   implicit none
   private

   !> Single two-phase flow solver, volume fraction solver, and material model set
   !> With corresponding time tracker
   type(mast),        public :: fs
   type(vfs),         public :: vf
   type(matm),        public :: matmod
   type(timetracker), public :: time
   type(hypre_str),   public :: ps
   type(hypre_str),   public :: vs

   !> Ensight postprocessing
   type(surfmesh) :: smesh
   type(ensight) :: ens_out, ens_out_smesh 
   type(event)   :: ens_evt, ens_evt_smesh 

   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,cvgfile

   public :: simulation_init,simulation_run,simulation_final

   !> Problem definition
   real(WP) :: ddrop
   real(WP), dimension(3) :: dctr
   integer :: relax_model
   real(WP) :: xshock,vshock,relshockvel
   real(WP) :: Grho0, GP0, Grho1, GP1, ST, Ma1, Ma, Lrho0, LP0, Mas, gamm_g

   !> variables for shock extraction
   logical :: extract_flag
   integer :: n_shock
   real(WP) :: tshock,final_xshock,delta,start_xshock

   !AS restart
   !>Provide a pardata and an event tracker for saving restarts
   type(event)   :: save_evt
   type(pardata) :: df
   logical       :: restarted
   
   !AS this input file stuff might not be needed for restarting our case?
   !> input file for the simulation
   type(inputfile) :: input
   
contains

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

   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none

      !AS restart
      ! set up an input file
      read_input: block
        use parallel, only: amRoot
        input=inputfile(amRoot=amRoot,filename='input')
      end block read_input

      ! if extraction flag is true --> run singlephase and then multiphase
      call param_read('Profile extraction flag',extract_flag)

      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)

         if (extract_flag.eqv.(.true.)) then ! singlephase simulation
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

            time%tmax = (final_xshock - start_xshock) / relshockvel !calculate final singlephase time based on final shock position
            print*, "Singlephase ending time: ", time%tmax

         else if (extract_flag.eqv.(.false.)) then ! multiphase simulation with extracted shock profile
            call param_read('Max time',time%tmax)
            call param_read('Multiphase max timestep size',time%dtmax)
            if (cfg%amRoot)then
               print*, "Multiphase ending time: ", time%tmax
            end if
         end if

         call param_read('Max steps',time%nmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker

      !AS restart
      ! Handle restart and saves here
      restart_and_save: block
        use string,  only: str_medium
        use filesys, only: makedir,isdir
        
        character(len=str_medium) :: timestamp
        integer, dimension(3)     :: iopartition

        ! create event for saving restart files
        save_evt=event(time,'Restart output')
        call param_read('Restart output period', save_evt%tper)

        ! check if we are restarting
        call param_read('Restart from', timestamp, default='')
        restarted=.false.; if(len_trim(timestamp).gt.0) restarted=.true.

        ! read in I/O partition
        call param_read('I/O partition',iopartition)

        print*, 'restart output period from input: ',save_evt%tper
        print*, 'I/O partition from input: ', iopartition !AS working 

        ! initialize pardata
        if (restarted) then
           print*, 'Defensive prog check. sim.f90 Line 176: init pardata running, restarting.' !AS working
           ! we are restarting, read file
           call df%initialize(pg=cfg,iopartition=iopartition,fdata='restart/data_'//trim(adjustl(timestamp))) !AS THIS NEEDS TO BE VERIFIED 1/17
        else
           ! we are not restarting, create a directory for the restart files
           print*, 'Defensive prog check. sim.f90 Line 181: creating directory for restart files, not restarting.' !AS working
           if (cfg%amRoot) then
              if(.not.isdir('restart')) call makedir('restart') !AS working
           end if

           ! prepare pardata object for saving restart files
           call df%initialize(pg=cfg,iopartition=iopartition,filename=trim(cfg%name),nval=2,nvar=29)
           df%valname=['t ','dt']
           df%varname=['Grho   ','Lrho   ','RHO    ','Ui     ','Vi     ','Wi     ','U      ','V      ','W      ','rhoUi  ','rhoVi  ','rhoWi  ','GrhoE  ','LrhoE  ','GP     ','LP     ','P      ','PA     ','GrhoSS2','LrhoSS2','RHOSS2 ','P11    ','P12    ','P13    ','P14    ','VOF    ']           
        end if
      end block restart_and_save

      !AS restart
      ! revisit timetracker to adjust the time and time step values if restarting
      update_timetracker: block
        if (restarted) then
           call df%pull(name='t' ,val=time%t )
           call df%pull(name='dt',val=time%dt)
           time%told=time%t-time%dt
        end if
      end block update_timetracker
        
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use mms_geom, only: cube_refine_vol
         use vfs_class, only: r2p,lvira,elvira,VFhi,VFlo,plicnet,flux
         use mast_class, only: dirichlet,clipped_neumann,bc_scope !AS restart
         use irl_fortran_interface !AS restart
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         real(WP), dimension(:,:,:), allocatable :: P11,P12,P13,P14 !AS restart, planes for interface reconstruction
         
         ! Create a VOF solver with lvira reconstruction
         call vf%initialize(cfg=cfg,reconstruction_method=lvira,name='VOF')

         ! Create a VOF solver with plicnet reconstruction with flux transport method
         !call vf%initialize(cfg=cfg,reconstruction_method=plicnet,transport_method=flux,name='VOF')

         !AS restart
         ! initialize the interface including restarts
         if (restarted)then
            !restart
            !read in the planes directly and set the IRL interface
            allocate(P11(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); call df%pull(name='P11',var=P11)
            allocate(P12(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); call df%pull(name='P12',var=P12)
            allocate(P13(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); call df%pull(name='P13',var=P13)
            allocate(P14(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); call df%pull(name='P14',var=P14)

            do k=vf%cfg%kmino_,vf%cfg%kmaxo_
               do j=vf%cfg%jmino_,vf%cfg%jmaxo_
                  do i=vf%cfg%imino_,vf%cfg%imaxo_
                     call setNumberOfPlanes(vf%liquid_gas_interface(i,j,k),1)
                     call setPlane(vf%liquid_gas_interface(i,j,k),0,[P11(i,j,k),P12(i,j,k),P13(i,j,k)],P14(i,j,k))
                  end do
               end do
            end do
            deallocate(P11,P12,P13,P14)
            
            call df%pull(name='VOF',var=vf%VF) !read in VOF
            ! ensure that boundaries are correct !!!BE CAREFUL HERE!!!
            ! AS assuming that the droplet hasn't reached the boundaries yet, there should only be gas on the edges of the domain
            vf%VF(vf%cfg%imino:vf%cfg%imin+1,:,:) = 0.0_WP  ! cells on the left 'slice' of the domain
            vf%VF(vf%cfg%imax:vf%cfg%imaxo,:,:) = 0.0_WP    ! cells on the right 'slice' of the domain
            vf%VF(:,vf%cfg%jmino:vf%cfg%jmin+1,:) = 0.0_WP  ! cells on the bottom 'slice' of the domain
            vf%VF(:,vf%cfg%jmax:vf%cfg%jmaxo,:) = 0.0_WP    ! cells on the top 'slice' of the domain
            vf%VF(:,:,vf%cfg%kmino:vf%cfg%kmin+1) = 0.0_WP  ! cells on the back 'slice' of the domain
            vf%VF(:,:,vf%cfg%kmax:vf%cfg%kmaxo) = 0.0_WP    ! cells on the front 'slice' of domain
            call vf%sync_interface()
            
            !AS set BC's on vof using MAST? Chase mentioned it might not be necessary to do it
            !boundary conditions using MAST apply_bcond_vf subroutine
            !don't forget to include mast in this block
            
            ! Boundary conditions on VF are built into the mast solver
            ! Update the band
            call vf%update_band()
            ! Perform interface reconstruction from VOF field
            ! removed build interface ---> missing from Chase's example
            !call vf%build_interface()
            ! Set initial interface at the boundaries
            call vf%set_full_bcond()
            ! Create discontinuous polygon mesh from IRL interface
            call vf%polygonalize_interface()
            ! Calculate distance from polygons
            call vf%distance_from_polygon()
            ! Calculate subcell phasic volumes
            call vf%subcell_vol()
            ! Calculate curvature
            call vf%get_curvature()
            ! Reset moments to guarantee compatibility with interface reconstruction
            call vf%reset_volume_moments()

         else
            !no restart
            ! Initialize liquid
            if (extract_flag.eqv.(.true.)) then ! if singlephase, set drop diameter to 0.0
               ddrop = 0.0_WP
            else if (extract_flag.eqv.(.false.)) then ! otherwise, initialize droplet based on diameter
               call param_read('Droplet diameter',ddrop)
               call param_read('Droplet location',dctr)
               do k=vf%cfg%kmino_,vf%cfg%kmaxo_
                  do j=vf%cfg%jmino_,vf%cfg%jmaxo_
                     do i=vf%cfg%imino_,vf%cfg%imaxo_
                        ! Set cube vertices
                        n=0
                        do sk=0,1
                           do sj=0,1
                              do si=0,1
                                 n=n+1; cube_vertex(:,n)=[vf%cfg%x(i+si),vf%cfg%y(j+sj),vf%cfg%z(k+sk)]
                              end do
                           end do
                        end do
                        ! Call adaptive refinement code to get volume and barycenters recursively
                        vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
                        if (vf%cfg%nz.eq.1) then
                           call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_cyl,0.0_WP,amr_ref_lvl)
                        else
                           call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_sphere,0.0_WP,amr_ref_lvl)
                        end if
                        vf%VF(i,j,k)=vol/vf%cfg%vol(i,j,k)
                        if (vf%VF(i,j,k).ge.VFlo.and.vf%VF(i,j,k).le.VFhi) then
                           vf%Lbary(:,i,j,k)=v_cent
                           vf%Gbary(:,i,j,k)=([vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]-vf%VF(i,j,k)*vf%Lbary(:,i,j,k))/(1.0_WP-vf%VF(i,j,k))
                           if (vf%cfg%nz.eq.1) vf%Gbary(3,i,j,k)=v_cent(3);
                        else
                           vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                           vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                        end if
                     end do
                  end do
               end do
            end if
            
            ! Boundary conditions on VF are built into the mast solver
            ! Update the band
            call vf%update_band()
            ! Perform interface reconstruction from VOF field
            call vf%build_interface()
            ! Set initial interface at the boundaries
            call vf%set_full_bcond()
            ! Create discontinuous polygon mesh from IRL interface
            call vf%polygonalize_interface()
            ! Calculate distance from polygons
            call vf%distance_from_polygon()
            ! Calculate subcell phasic volumes
            call vf%subcell_vol()
            ! Calculate curvature
            call vf%get_curvature()
            ! Reset moments to guarantee compatibility with interface reconstruction
            call vf%reset_volume_moments()
         end if
      end block create_and_initialize_vof

      ! Create a compressible two-phase flow solver
      create_and_initialize_flow_solver: block
         use mast_class,      only: clipped_neumann,dirichlet,bc_scope,bcond,mech_egy_mech_hhz
         use hypre_str_class, only: pcg_pfmg
         use mathtools,       only: Pi
         use parallel,        only: amRoot
         use messager,        only: die
         use param,           only: param_read, param_exists

         integer :: i,j,k,n,nx
         real(WP), dimension(3) :: xyz
         real(WP) :: gamm_l,Pref_l,gamm_g,visc_l,visc_g,Pref,cv_l0,cv_g0,kappa_l,kappa_g
         real(WP) :: xshock,vshock,relshockvel,Lx
         real(WP) :: Grho0, GP0, Grho1, GP1, ST, Ma1, Ma, Lrho0, LP0, Mas
         type(bcond), pointer :: mybc

         ! variables for shock extraction
         integer :: n_shock,q,shock_index
         real(WP) :: final_xshock,delta, dx,tol,shock_loc

         ! variables for reading in shock profile
         real(WP), dimension(:),  allocatable :: Grho_profile, GrhoE_profile, Ui_profile, GP_profile

         ! set up for shock profile
         call param_read('n_shock',n_shock) ! number of points to the left and right of shock for profile
         call param_read('Lx',Lx); call param_read('nx',nx)
         
         dx = Lx/nx ! mesh spacing in uniform region
         tol = dx/2 ! set tolerance for reading in shock profile
         delta = 2*dx*n_shock ! shock thickness
         if (amRoot) then
            print*, "Total shock profile points: ", 2*n_shock
            print*, "Shock thickness: ", delta
         end if
         
         ! Initialize conditions
         if (.not.restarted)then ! if we are not restarting then read in the numerical shock 
            if (extract_flag.eqv.(.true.)) then ! singlephase simulation
               call param_read('Single phase shock location',start_xshock); print*,"start_xshock", start_xshock !singlephase shock starting location
               call param_read('Final shock location',final_xshock) !final shock location
               ddrop = 0.0_WP !set to singelphase
            else if (extract_flag.eqv.(.false.)) then ! multiphase simulation
               call param_read('Shock location',xshock)
               call param_read('Droplet diameter',ddrop) !set to multiphase
               
               !allocate profile arrays
               allocate(Grho_profile(2*n_shock+1))
               allocate(GrhoE_profile(2*n_shock+1))
               allocate(Ui_profile(2*n_shock+1))
               allocate(GP_profile(2*n_shock+1))
               
               ! read in singlephase profile data and store in variables
               open(unit=1, file='Grho_profile.dat')
               read(1,*) Grho_profile
               close(1)
               
               open(unit=2, file='GrhoE_profile.dat')
               read(2,*) GrhoE_profile
               close(2)
               
               open(unit=3, file='Ui_profile.dat')
               read(3,*) Ui_profile
               close(3)
               
               open(unit=4, file='GP_profile.dat')
               read(4,*) GP_profile
               close(4)
            end if
         end if

         ! Create material model class
         matmod=matm(cfg=cfg,name='Liquid-gas models')

         ! Get EOS parameters from input
         call param_read('Liquid Pref', Pref_l)
         call param_read('Liquid gamma',gamm_l)
         call param_read('Gas gamma',gamm_g)

         ! Register equations of state
         call matmod%register_stiffenedgas('liquid',gamm_l,Pref_l)
         call matmod%register_idealgas('gas',gamm_g)

         ! Create flow solver
         fs=mast(cfg=cfg,name='Two-phase All-Mach',vf=vf)

         call param_read('Liquid dynamic viscosity',visc_l)
         call param_read('Gas dynamic viscosity',visc_g)
         call param_read('Liquid thermal conductivity',kappa_l)
         call param_read('Gas thermal conductivity',kappa_g)
         call param_read('Liquid specific heat (constant vol)',cv_l0)
         call param_read('Gas specific heat (constant vol)',cv_g0)

         ! Register flow solver variables with material models
         call matmod%register_thermoflow_variables('liquid',fs%Lrho,fs%Ui,fs%Vi,fs%Wi,fs%LrhoE,fs%LP)
         call matmod%register_thermoflow_variables('gas'   ,fs%Grho,fs%Ui,fs%Vi,fs%Wi,fs%GrhoE,fs%GP)
         call matmod%register_diffusion_thermo_models(viscconst_gas=visc_g, viscconst_liquid=visc_l,hdffconst_gas=kappa_g, hdffconst_liquid=kappa_l,sphtconst_gas=cv_g0,sphtconst_liquid=cv_l0)

         ! Read in surface tension coefficient
         call param_read('Surface tension coefficient',fs%sigma)

         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg,nst=7)
         ps%maxlevel=10
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)

         ! Configure implicit velocity solver
         vs=hypre_str(cfg=cfg,name='Velocity',method=pcg_pfmg,nst=7)
         call param_read('Implicit iteration',vs%maxit)
         call param_read('Implicit tolerance',vs%rcvg)

         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)

         ! Liquid and gas density
         call param_read('Liquid density',Lrho0); fs%Lrho = Lrho0
         call param_read('Pre-shock density',Grho0,default=1.204_WP)
         call param_read('Pre-shock pressure',GP0,default=1.01325e5_WP)
         call param_read('Mach number of shock',Ma,default=1.47_WP)

         ! Initially 0 velocity in y and z
         fs%Vi = 0.0_WP; fs%Wi = 0.0_WP

         ! Zero face velocities as well for the sake of dirichlet boundaries
         fs%V = 0.0_WP; fs%W = 0.0_WP

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
            print*,"===== Problem Setup Description ====="
            print*,'Mach number', Ma
            print*,'Pre-shock:  Density',Grho0,'Pressure',GP0
            print*,'Post-shock: Density',Grho1,'Pressure',GP1,'Gas Velocity',vshock
            print*,'Shock velocity', relshockvel
         end if
         
         if (.not.restarted)then !AS restart   
            ! Initialize gas phase quantities
            if (extract_flag.eqv.(.true.))then
               !singlephase initialization
               do  i=fs%cfg%imino_,fs%cfg%imaxo_
                  if (cfg%xm(i).lt.start_xshock) then ! post shock values
                     fs%Grho(i,:,:) = Grho1
                     fs%Ui(i,:,:) = vshock
                     fs%GP(i,:,:) = GP1
                     fs%GrhoE(i,:,:) = matmod%EOS_energy(GP1,Grho1,vshock,0.0_WP,0.0_WP,'gas')
                  else ! pre shock values
                     fs%Grho(i,:,:) = Grho0
                     fs%Ui(i,:,:) = 0.0_WP
                     fs%GP(i,:,:) = GP0
                  fs%GrhoE(i,:,:) = matmod%EOS_energy(GP0,Grho0,0.0_WP,0.0_WP,0.0_WP,'gas')
               end if
            end do
            
         else
            !multiphase numerical initialization
            do i=fs%cfg%imino_,fs%cfg%imaxo_
               if (cfg%xm(i).le.xshock) then
                  fs%Grho(i,:,:) = Grho1
                  fs%Ui(i,:,:) = vshock
                  fs%GP(i,:,:) = GP1
                  fs%GrhoE(i,:,:) = matmod%EOS_energy(GP1,Grho1,vshock,0.0_WP,0.0_WP,'gas')
               elseif (cfg%xm(i).ge.xshock) then
                  fs%Grho(i,:,:) = Grho0
                  fs%Ui(i,:,:) = 0.0_WP
                  fs%GP(i,:,:) = GP0
                  fs%GrhoE(i,:,:) = matmod%EOS_energy(GP0,Grho0,0.0_WP,0.0_WP,0.0_WP,'gas')
               end if
            end do
            
            ! find shock index
            do i=cfg%imino_,cfg%imaxo_
               if ((cfg%xm(i).lt.(xshock+tol)).and.(cfg%xm(i).gt.(xshock-tol))) then
                  print*, "The shock has been found at index: ", i
                  shock_index=i
                  print*, "stored shock index", shock_index
                  shock_loc = cfg%xm(shock_index) ! store location of shock corresponding to index
                  print*, "The found shock location (cell center) is: ", shock_loc
               end if
            end do
            
            do i=cfg%imino_,cfg%imaxo_ ! read in shock profile
               if ((cfg%xm(i).le.cfg%xm(shock_index+n_shock)).and.(cfg%xm(i).ge.cfg%xm(shock_index-n_shock)))then
                  fs%Grho(i,:,:) = Grho_profile(i+n_shock-shock_index+1)
                  fs%GP(i,:,:) = GP_profile(i+n_shock-shock_index+1)
                  fs%GrhoE(i,:,:) = GrhoE_profile(i+n_shock-shock_index+1)
                  fs%Ui(i,:,:) = Ui_profile(i+n_shock-shock_index+1)
               end if
            end do
         end if
         
         ! Calculate liquid pressure
         if (fs%cfg%nz.eq.1) then
            ! Cylinder configuration, curv = 1/r
            ! LP0 = Pref + 2.0/ddrop*fs%sigma
            LP0 = GP0 + 2.0/ddrop*fs%sigma
         else
            ! Sphere configuration, curv = 1/r + 1/r
            ! LP0 = Pref + 4.0/ddrop*fs%sigma
            LP0 = GP0 + 4.0/ddrop*fs%sigma
         end if
         
         !initialize liquid quantities
         fs%Lrho = Lrho0
         fs%LP = LP0
         ! Initialize liquid energy, with surface tension
         fs%LrhoE = matmod%EOS_energy(LP0,Lrho0,0.0_WP,0.0_WP,0.0_WP,'liquid')
         
         ! Define boundary conditions - initialized values are intended dirichlet values too, for the cell centers
         call fs%add_bcond(name= 'inflow',type=dirichlet      ,locator=left_of_domain ,face='x',dir=-1)
         call fs%add_bcond(name='outflow',type=clipped_neumann,locator=right_of_domain,face='x',dir=+1)
         
         ! Calculate face velocities
         call fs%interp_vel_basic(vf,fs%Ui,fs%Vi,fs%Wi,fs%U,fs%V,fs%W)
         
         ! Apply face BC - inflow
         call fs%get_bcond('inflow',mybc)
         do n=1,mybc%itr%n_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%U(i,j,k)=vshock
         end do
         
         ! Apply face BC - outflow
         bc_scope = 'velocity'
         call fs%apply_bcond(time%dt,bc_scope)
         
         ! Calculate mixture density and momenta
         fs%RHO   = (1.0_WP-vf%VF)*fs%Grho  + vf%VF*fs%Lrho
         fs%rhoUi = fs%RHO*fs%Ui; fs%rhoVi = fs%RHO*fs%Vi; fs%rhoWi = fs%RHO*fs%Wi
         
         ! Perform initial pressure relax
         relax_model = mech_egy_mech_hhz
         call fs%pressure_relax(vf,matmod,relax_model)
         
         ! Calculate initial phase and bulk moduli
         call fs%init_phase_bulkmod(vf,matmod)
         call fs%reinit_phase_pressure(vf,matmod)
         call fs%harmonize_advpressure_bulkmod(vf,matmod)
         
         ! Set initial pressure to harmonized field based on internal energy
         fs%P = fs%PA
      else !AS restart
         ! Read data
         call df%pull(name='Grho'   ,var=fs%Grho   ); 
         call df%pull(name='Lrho'   ,var=fs%Lrho   ); 
         call df%pull(name='RHO'    ,var=fs%RHO   ); 
         call df%pull(name='Ui'     ,var=fs%Ui     ); 
         call df%pull(name='Vi'     ,var=fs%Vi     ); 
         call df%pull(name='Wi'     ,var=fs%Wi     ); 
         call df%pull(name='U'      ,var=fs%U      ); 
         call df%pull(name='V'      ,var=fs%V      ); 
         call df%pull(name='W'      ,var=fs%W      ); 
         call df%pull(name='rhoUi'  ,var=fs%rhoUi  ); 
         call df%pull(name='rhoVi'  ,var=fs%rhoVi  ); 
         call df%pull(name='rhoWi'  ,var=fs%rhoWi  ); 
         call df%pull(name='GrhoE'  ,var=fs%GrhoE  ); 
         call df%pull(name='LrhoE'  ,var=fs%LrhoE  ); 
         call df%pull(name='GP'     ,var=fs%GP     ); 
         call df%pull(name='LP'     ,var=fs%LP     ); 
         call df%pull(name='P'      ,var=fs%P      ); 
         call df%pull(name='PA'     ,var=fs%PA     ); 
         call df%pull(name='GrhoSS2',var=fs%GrhoSS2); 
         call df%pull(name='LrhoSS2',var=fs%LrhoSS2); 
         call df%pull(name='RHOSS2',var=fs%RHOSS2); 
         
         ! define boundary conditions
         call fs%add_bcond(name= 'inflow',type=dirichlet      ,locator=left_of_domain ,face='x',dir=-1)
         call fs%add_bcond(name='outflow',type=clipped_neumann,locator=right_of_domain,face='x',dir=+1)

         ! Calculate face velocities
         call fs%interp_vel_basic(vf,fs%Ui,fs%Vi,fs%Wi,fs%U,fs%V,fs%W)
         
         ! Apply face BC - inflow
         call fs%get_bcond('inflow',mybc)
         do n=1,mybc%itr%n_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%U(i,j,k)=vshock
         end do
         
         ! Apply face BC - outflow
         bc_scope = 'velocity'
         call fs%apply_bcond(time%dt,bc_scope)
         
         !calculate mixture density and momenta
         fs%RHO   = (1.0_WP-vf%VF)*fs%Grho  + vf%VF*fs%Lrho
         fs%rhoUi = fs%RHO*fs%Ui; fs%rhoVi = fs%RHO*fs%Vi; fs%rhoWi = fs%RHO*fs%Wi
         
         ! choose relaxation model and Perform initial pressure relax
         relax_model = mech_egy_mech_hhz
         call fs%pressure_relax(vf,matmod,relax_model)
         
      end if
      
    end block create_and_initialize_flow_solver

       ! create surfmesh object for interface polygon output
       create_smesh: block
         use irl_fortran_interface
         integer :: i,j,k,nplane,np
         smesh=surfmesh(nvar=1,name='plic')
         smesh%varname(1)='curv'
         call vf%update_surfmesh(smesh)
         smesh%var(1,:)=0.0_WP
         np=0;
         do k=vf%cfg%kmin_,vf%cfg%kmax_
            do j=vf%cfg%jmin_,vf%cfg%jmax_
               do i=vf%cfg%imin_,vf%cfg%imax_
                  do nplane=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j,k))
                     if (getNumberOfVertices(vf%interface_polygon(nplane,i,j,k)).gt.0) then
                        np=np+1;
                        smesh%var(1,np)=vf%curv(i,j,k)
                     end if
                  end do
               end do
            end do
         end do
       end block create_smesh

      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='ShockDroplet')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',fs%Ui,fs%Vi,fs%Wi)
         call ens_out%add_scalar('P',fs%P)
         call ens_out%add_scalar('PA',fs%PA)
         call ens_out%add_scalar('Grho',fs%Grho)
         call ens_out%add_scalar('Lrho',fs%Lrho)
         call ens_out%add_scalar('Density',fs%RHO)
         call ens_out%add_scalar('Bulkmod',fs%RHOSS2)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('curvature',vf%curv)
         call ens_out%add_scalar('Mach',fs%Mach)
         call ens_out%add_scalar('fvf',cfg%VF)
         call ens_out%add_scalar('Tmptr',fs%Tmptr)
         call ens_out%add_scalar('SL_x',fs%sl_x) 
         call ens_out%add_scalar('SL_y',fs%sl_y) 
         call ens_out%add_scalar('SL_z',fs%sl_z) 
         call ens_out%add_scalar('LP',fs%LP) 
         call ens_out%add_scalar('GP',fs%GP) 
         call ens_out%add_scalar('LrhoE',fs%LrhoE) 
         call ens_out%add_scalar('GrhoE',fs%GrhoE) 
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
       end block create_ensight

       ! block for writing smesh data more frequently than field variables
       create_ensight_smesh: block
         real(WP) :: smesh_tper ! declare variable for smesh output frequency
         call param_read('Ensight smesh output period', smesh_tper)

         ! create ensight output from cfg for smesh surface reconstruction
         ens_out_smesh=ensight(cfg=cfg,name='droplet_smesh')

         ! create event for ensight output
         ens_evt_smesh=event(time=time,name='Ensight output smesh')
         ens_evt_smesh%tper = smesh_tper

         ! add variables to output
         call ens_out_smesh%add_surface('smesh',smesh)

         ! output to ensight
         if (ens_evt_smesh%occurs()) call ens_out_smesh%write_data(time%t)
       end block create_ensight_smesh

      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call vf%get_max()
         ! Create simulation monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%RHOmin,'RHOmin')
         call mfile%add_column(fs%RHOmax,'RHOmax')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(fs%Tmax,'Tmax')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLst,'STension CFL')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%add_column(fs%CFLa_x,'Acoustic xCFL')
         call cflfile%add_column(fs%CFLa_y,'Acoustic yCFL')
         call cflfile%add_column(fs%CFLa_z,'Acoustic zCFL')
         call cflfile%write()
         ! Create convergence monitor
         cvgfile=monitor(fs%cfg%amRoot,'cvg')
         call cvgfile%add_column(time%n,'Timestep number')
         call cvgfile%add_column(time%it,'Iteration')
         call cvgfile%add_column(time%t,'Time')
         call cvgfile%add_column(fs%impl_it_x,'Impl_x iteration')
         call cvgfile%add_column(fs%impl_rerr_x,'Impl_x error')
         call cvgfile%add_column(fs%impl_it_y,'Impl_y iteration')
         call cvgfile%add_column(fs%impl_rerr_y,'Impl_y error')
         call cvgfile%add_column(fs%implicit%it,'Impl_z iteration')
         call cvgfile%add_column(fs%implicit%rerr,'Impl_z error')
         call cvgfile%add_column(fs%psolv%it,'Pressure iteration')
         call cvgfile%add_column(fs%psolv%rerr,'Pressure error')
      end block create_monitor

   end subroutine simulation_init

   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
      use messager, only: die
      implicit none

      ! Perform time integration
      do while (.not.time%done())

         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Reinitialize phase pressure by syncing it with conserved phase energy
         call fs%reinit_phase_pressure(vf,matmod)
         fs%Uiold=fs%Ui; fs%Viold=fs%Vi; fs%Wiold=fs%Wi
         fs%RHOold = fs%RHO

         ! Remember old flow variables (phase)
         fs%Grhoold = fs%Grho; fs%Lrhoold = fs%Lrho
         fs%GrhoEold=fs%GrhoE; fs%LrhoEold=fs%LrhoE
         fs%GPold   =   fs%GP; fs%LPold   =   fs%LP

         ! Remember old interface, including VF and barycenters
         call vf%copy_interface_to_old()

         ! Create in-cell reconstruction
         call fs%flow_reconstruct(vf)

         ! Zero variables that will change during subiterations
         fs%P = 0.0_WP
         fs%Pjx = 0.0_WP; fs%Pjy = 0.0_WP; fs%Pjz = 0.0_WP
         fs%Hpjump = 0.0_WP

         ! Determine semi-Lagrangian advection flag
         call fs%flag_sl(time%dt,vf)

         ! Perform sub-iterations
         do while (time%it.le.time%itmax)

            ! Predictor step, involving advection and pressure terms
            call fs%advection_step(time%dt,vf,matmod)

            ! Viscous step
            call fs%diffusion_src_explicit_step(time%dt,vf,matmod)

            ! Prepare pressure projection
            call fs%pressureproj_prepare(time%dt,vf,matmod)

            ! Initialize and solve Helmholtz equation
            call fs%psolv%setup()
            fs%psolv%sol=fs%PA-fs%P
            call fs%psolv%solve()
            call fs%cfg%sync(fs%psolv%sol)

            ! Perform corrector step using solution
            fs%P=fs%P+fs%psolv%sol
            call fs%pressureproj_correct(time%dt,vf,fs%psolv%sol)

            ! Record convergence monitor
            call cvgfile%write()
            ! Increment sub-iteration counter
            time%it=time%it+1

         end do

         ! Pressure relaxation
         call fs%pressure_relax(vf,matmod,relax_model)

         ! Output to ensight
         if (ens_evt%occurs()) then

            !update surfmesh object
            update_smesh: block
              use irl_fortran_interface
              integer :: i,j,k,nplane,np
              ! Transfer polygons to smesh
              call vf%update_surfmesh(smesh)
              ! Also populate nplane variable
              smesh%var(1,:)=0.0_WP
              np=0
              do k=vf%cfg%kmin_,vf%cfg%kmax_
                 do j=vf%cfg%jmin_,vf%cfg%jmax_
                    do i=vf%cfg%imin_,vf%cfg%imax_
                       do nplane=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j,k))
                          if (getNumberOfVertices(vf%interface_polygon(nplane,i,j,k)).gt.0) then
                             np=np+1;
                             smesh%var(1,np)=vf%curv(i,j,k)
                          end if
                       end do
                    end do
                 end do
              end do
            end block update_smesh

            !AS write data
            call ens_out%write_data(time%t)

         end if

         if (ens_evt_smesh%occurs()) then
            call ens_out_smesh%write_data(time%t)
         end if

         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         call fs%get_viz()
         call mfile%write()
         call cflfile%write()

         if (save_evt%occurs())then
            save_restart: block
              use irl_fortran_interface
              use string, only: str_medium
              character(len=str_medium) :: timestamp
              real(WP), dimension(:,:,:), allocatable :: P11,P12,P13,P14
              integer :: i,j,k
              real(WP), dimension(4) :: plane

              ! Handle IRL data
              allocate(P11(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
              allocate(P12(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
              allocate(P13(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
              allocate(P14(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
              ! Store IRL data
              do k=vf%cfg%kmino_,vf%cfg%kmaxo_
                 do j=vf%cfg%jmino_,vf%cfg%jmaxo_
                    do i=vf%cfg%imino_,vf%cfg%imaxo_
                       ! First plane
                       plane=getPlane(vf%liquid_gas_interface(i,j,k),0)
                       P11(i,j,k)=plane(1); P12(i,j,k)=plane(2); P13(i,j,k)=plane(3); P14(i,j,k)=plane(4)
                       plane=0.0_WP
                    end do
                 end do
              end do

              ! Prefix for files
              write(timestamp,'(es12.5)') time%t
              ! Populate df and write it
              call df%push(name='t'      ,val=time%t    )
              call df%push(name='dt'     ,val=time%dt   )
              call df%push(name='Grho'   ,var=fs%Grho   )
              call df%push(name='Lrho'   ,var=fs%Lrho   )
              call df%push(name='RHO'    ,var=fs%RHO    )
              call df%push(name='Ui'     ,var=fs%Ui     )
              call df%push(name='Vi'     ,var=fs%Vi     )
              call df%push(name='Wi'     ,var=fs%Wi     )
              call df%push(name='U'      ,var=fs%U      )
              call df%push(name='V'      ,var=fs%V      )
              call df%push(name='W'      ,var=fs%W      )
              call df%push(name='rhoUi'  ,var=fs%rhoUi  )
              call df%push(name='rhoVi'  ,var=fs%rhoVi  )
              call df%push(name='rhoWi'  ,var=fs%rhoWi  )
              call df%push(name='GrhoE'  ,var=fs%GrhoE  )
              call df%push(name='LrhoE'  ,var=fs%LrhoE  )
              call df%push(name='GP'     ,var=fs%GP     )
              call df%push(name='LP'     ,var=fs%LP     )
              call df%push(name='P'      ,var=fs%P      ) 
              call df%push(name='PA'     ,var=fs%PA     ) 
              call df%push(name='GrhoSS2',var=fs%GrhoSS2) 
              call df%push(name='LrhoSS2',var=fs%LrhoSS2) 
              call df%push(name='RHOSS2',var=fs%RHOSS2) 
              call df%push(name='P11'    ,var=P11            )
              call df%push(name='P12'    ,var=P12            )
              call df%push(name='P13'    ,var=P13            )
              call df%push(name='P14'    ,var=P14            )
              call df%push(name='VOF'    ,var=vf%VF     )
              call df%write(fdata='restart/data_'//trim(adjustl(timestamp)))
              
            end block save_restart

         end if

      end do
   end subroutine simulation_run

   !> Finalize the NGA2 simulation
   subroutine simulation_final
     use param, only: param_read
     use parallel,   only: amRoot
     implicit none
     integer :: i,j,shock_index,n_shock, nx
     real(WP) :: final_xshock, delta, start_ref,Lx
     real(WP) :: tol ! tolerance for finding final shock location in singlephase

     call param_read('nx',nx)
     call param_read('Lx',Lx);  call param_read('Lx ref',start_ref);
     call param_read('n_shock',n_shock)
     call param_read('Final shock location',final_xshock)
     call param_read('Lx ref', start_ref, default=0.0_WP);
     
     delta = 2*cfg%dx(1)*n_shock !shock thickness
     tol = (Lx - start_ref)/nx ! set the tolerance to the mesh spacing in the uniform region
     
     if (extract_flag.eqv.(.true.)) then
        !set up shock profile data files
        open(1, file='Grho_profile.dat')
        open(2, file='GrhoE_profile.dat')
        open(3, file='Ui_profile.dat')
        open(4, file='GP_profile.dat')

        do i=cfg%imino_,cfg%imaxo_
           if ((cfg%xm(i).lt.(final_xshock+tol)).and.(cfg%xm(i).gt.(final_xshock-tol))) then
              print*, "The shock has been found at index: ", i
              shock_index=i
              print*, "stored shock index", shock_index
           end if
        end do

        do i=shock_index-n_shock,shock_index+n_shock ! this is now writing 2*n_shock points
           write(1,*) fs%Grho(i,1,1)
           write(2,*) fs%GrhoE(i,1,1)
           write(3,*) fs%Ui(i,1,1)
           write(4,*) fs%GP(i,1,1)
        end do

        close(1)
        close(2)
        close(3)
        close(4)
     end if

   end subroutine simulation_final

end module simulation
