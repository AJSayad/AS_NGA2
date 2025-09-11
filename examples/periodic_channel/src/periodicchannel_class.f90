module periodicchannel_class
   use precision,         only: WP
   use config_class,      only: config
   use fft2d_class,       only: fft2d
   use ddadi_class,       only: ddadi
   use incomp_class,      only: incomp
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   use pardata_class,     only: pardata
   implicit none
   private
   
   public :: periodicchannel
   
   !> periodicchannel object
   type :: periodicchannel
      
      !> Config
      type(config) :: cfg
      
      !> Flow solver
      type(incomp)      :: fs     !< Incompressible flow solver
      type(fft2d)       :: ps     !< Fourier-accelerated pressure solver
      type(timetracker) :: time   !< Time info
      
      !> Implicit solver
      logical     :: use_implicit !< Is an implicit solver used?
      type(ddadi) :: vs           !< DDADI solver for velocity
      
      !> SGS model
      logical        :: use_sgs   !< Is an LES model used?
      type(sgsmodel) :: sgs       !< SGS model for eddy viscosity
      real(WP), dimension(:,:,:,:,:), allocatable :: gradU !< Velocity gradient
      
      !> Ensight postprocessing
      type(ensight) :: ens_out    !< Event trigger for Ensight output
      type(event)   :: ens_evt    !< Ensight output for flow variables
      
      !> Simulation monitoring files
      type(monitor) :: mfile      !< General simulation monitoring
      type(monitor) :: cflfile    !< CFL monitoring

      !> Work arrays
      real(WP), dimension(:,:,:), allocatable :: resU,resV,resW      !< Residuals
      real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi            !< Cell-centered velocities
      
      !> Flow conditions
      real(WP) :: visc,meanU,bforce
      
      !> Provide a pardata object for restarts
      logical       :: restarted
      type(pardata) :: df
      type(event)   :: save_evt
      
   contains
      procedure :: init                            !< Initialize pipe simulation
      procedure :: step                            !< Advance pipe simulation by one time step
      procedure :: final                           !< Finalize pipe simulation
   end type periodicchannel
   
contains
   
   
   !> Initialization of periodicchannel simulation
   subroutine init(this)
      use parallel, only: group
      use param,    only: param_read
      implicit none
      class(periodicchannel), intent(inout) :: this
      
      
      ! Initialize the config
      initialize_config: block
         use sgrid_class, only: sgrid,cartesian
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz,stretch
         real(WP), dimension(:), allocatable :: x,y,z
         type(sgrid) :: grid
         integer, dimension(3) :: partition
         ! Grid specification
         call param_read('Lx',Lx); call param_read('nx',nx); allocate(x(nx+1))
         Ly=1.0_WP;                call param_read('ny',ny); allocate(y(ny+1)); call param_read('Stretching',stretch)
         call param_read('Lz',Lz); call param_read('nz',nz); allocate(z(nz+1))
         ! Create simple rectilinear grid in x and z, tanh-stretched grid in y
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx
         end do
         do j=1,ny+1
            y(j)=0.5_WP*Ly*tanh(stretch*(2.0_WP*real(j-1,WP)/real(ny,WP)-1.0_WP))/tanh(stretch)
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=1,x=x,y=y,z=z,xper=.true.,yper=.false.,zper=.true.,name='channel')
         ! Read in partition
         call param_read('Partition',partition,short='p')
         ! Create partitioned grid
         this%cfg=config(grp=group,decomp=partition,grid=grid)
         ! Create walls top and bottom
         this%cfg%VF=0.0_WP
         this%cfg%VF(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)=1.0_WP
         call this%cfg%sync(this%cfg%VF)
         ! Output the grid
         if (this%cfg%amRoot) call this%cfg%pgrid%write('channel.grid')
      end block initialize_config
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         this%time=timetracker(amRoot=this%cfg%amRoot)
         call param_read('Max time',this%time%tmax)
         call param_read('Max timestep size',this%time%dtmax)
         call param_read('Max cfl number',this%time%cflmax)
         this%time%dt=this%time%dtmax
         this%time%itmax=2
      end block initialize_timetracker
      
      
      ! Create an incompressible flow solver without bconds
      create_flow_solver: block
         ! Create flow solver
         this%fs=incomp(cfg=this%cfg,name='Incompressible NS')
         ! Set the flow properties
         this%fs%rho=1.0_WP
         call param_read('Reynolds number',this%visc)
         this%visc=1.0_WP/this%visc
         this%fs%visc=this%visc
         ! Configure pressure solver
         this%ps=fft2d(cfg=this%cfg,name='Pressure',nst=7)
         ! Check if implicit velocity solver is used
         call param_read('Use implicit solver',this%use_implicit)
         if (this%use_implicit) then
            ! Configure implicit solver
            this%vs=ddadi(cfg=this%cfg,name='Velocity',nst=7)
            ! Finish flow solver setup
            call this%fs%setup(pressure_solver=this%ps,implicit_solver=this%vs)
         else
            ! Finish flow solver setup
            call this%fs%setup(pressure_solver=this%ps)
         end if
      end block create_flow_solver
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(this%resU(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resV(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resW(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Ui(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Vi(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Wi(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize our velocity field
      initialize_velocity: block
         use mathtools, only: twoPi
         use random,    only: random_uniform
         integer :: i,j,k
         real(WP) :: amp
         ! Initial fields
         where (this%fs%umask.eq.0) this%fs%U=1.0_WP
         this%fs%V=0.0_WP; this%fs%W=0.0_WP; this%fs%P=0.0_WP
         this%meanU=1.0_WP; this%bforce=0.0_WP
         ! Add fluctuations for faster transition
         call param_read('Fluctuation amp',amp,default=0.2_WP)
         do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
            do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
               do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                  this%fs%U(i,j,k)=this%fs%U(i,j,k)+random_uniform(lo=-0.5_WP*amp,hi=0.5_WP*amp)+amp*cos(8.0_WP*twoPi*this%fs%cfg%xm(i)/this%fs%cfg%xL)*cos(8.0_WP*twoPi*this%fs%cfg%zm(k)/this%fs%cfg%zL)
                  this%fs%W(i,j,k)=this%fs%W(i,j,k)+random_uniform(lo=-0.5_WP*amp,hi=0.5_WP*amp)+amp*cos(8.0_WP*twoPi*this%fs%cfg%xm(i)/this%fs%cfg%xL)*cos(8.0_WP*twoPi*this%fs%cfg%zm(k)/this%fs%cfg%zL)
               end do
            end do
         end do
         call this%fs%cfg%sync(this%fs%U)
         call this%fs%cfg%sync(this%fs%W)
         ! Compute cell-centered velocity
         call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
         ! Compute divergence
         call this%fs%get_div()
      end block initialize_velocity
      
      
      ! Create an LES model
      create_sgs: block
         call param_read('Use SGS model',this%use_sgs)
         if (this%use_sgs) then
            allocate(this%gradU(1:3,1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            this%sgs=sgsmodel(cfg=this%fs%cfg,umask=this%fs%umask,vmask=this%fs%vmask,wmask=this%fs%wmask)
         end if
      end block create_sgs
      
      
      ! Handle restart here
      perform_restart: block
         use string,  only: str_medium
         use filesys, only: makedir,isdir
         character(len=str_medium) :: filename
         integer, dimension(3) :: iopartition
         ! Create event for saving restart files
         this%save_evt=event(this%time,'Channel restart output')
         call param_read('Restart output period',this%save_evt%tper)
         ! Read in the partition for I/O
         call param_read('I/O partition',iopartition)
         ! Prepare a restart directory for storing files for restart
         if (this%cfg%amRoot.and..not.isdir('restart')) call makedir('restart')
         ! Check if a restart file was provided
         call param_read('Restart from',filename,default='')
         this%restarted=.false.; if (len_trim(filename).gt.0) this%restarted=.true.
         ! Perform pardata initialization
         if (this%restarted) then
            ! Read in the file
            call this%df%initialize(pg=this%cfg,iopartition=iopartition,fdata=trim(filename))
            ! Put the data at the right place
            call this%df%pull(name='U',var=this%fs%U)
            call this%df%pull(name='V',var=this%fs%V)
            call this%df%pull(name='W',var=this%fs%W)
            call this%df%pull(name='P',var=this%fs%P)
            ! Update cell-centered velocity
            call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
            ! Update divergence
            call this%fs%get_div()
            ! Also update time
            call this%df%pull(name='t' ,val=this%time%t )
            call this%df%pull(name='dt',val=this%time%dt)
            this%time%told=this%time%t-this%time%dt
            !this%time%dt=this%time%dtmax !< Force max timestep size anyway
            ! Get bodyforce
            call this%df%pull(name='bf',val=this%bforce)
         else
            ! If we are not restarting, we will still need a datafile for saving restart files
            call this%df%initialize(pg=this%cfg,iopartition=iopartition,filename=trim(this%cfg%name),nval=3,nvar=4)
            this%df%valname=['dt','t ','bf']; this%df%varname=['U','V','W','P']
         end if
      end block perform_restart
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         this%ens_out=ensight(cfg=this%cfg,name='channel')
         ! Create event for Ensight output
         this%ens_evt=event(time=this%time,name='Ensight output')
         call param_read('Ensight output period',this%ens_evt%tper)
         ! Add variables to output
         call this%ens_out%add_vector('velocity',this%Ui,this%Vi,this%Wi)
         call this%ens_out%add_scalar('pressure',this%fs%P)
         if (this%use_sgs) call this%ens_out%add_scalar('visc_sgs',this%sgs%visc)
         ! Output to ensight
         if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call this%fs%get_cfl(this%time%dt,this%time%cfl)
         call this%fs%get_max()
         ! Create simulation monitor
         this%mfile=monitor(this%fs%cfg%amRoot,'simulation')
         call this%mfile%add_column(this%time%n,'Timestep number')
         call this%mfile%add_column(this%time%t,'Time')
         call this%mfile%add_column(this%time%dt,'Timestep size')
         call this%mfile%add_column(this%time%cfl,'Maximum CFL')
         call this%mfile%add_column(this%meanU,'Mean U')
         call this%mfile%add_column(this%bforce,'Body force')
         call this%mfile%add_column(this%fs%Umax,'Umax')
         call this%mfile%add_column(this%fs%Vmax,'Vmax')
         call this%mfile%add_column(this%fs%Wmax,'Wmax')
         call this%mfile%add_column(this%fs%Pmax,'Pmax')
         call this%mfile%add_column(this%fs%divmax,'Maximum divergence')
         call this%mfile%add_column(this%fs%psolv%it,'Pressure iteration')
         call this%mfile%add_column(this%fs%psolv%rerr,'Pressure error')
         call this%mfile%write()
         ! Create CFL monitor
         this%cflfile=monitor(this%fs%cfg%amRoot,'cfl')
         call this%cflfile%add_column(this%time%n,'Timestep number')
         call this%cflfile%add_column(this%time%t,'Time')
         call this%cflfile%add_column(this%fs%CFLc_x,'Convective xCFL')
         call this%cflfile%add_column(this%fs%CFLc_y,'Convective yCFL')
         call this%cflfile%add_column(this%fs%CFLc_z,'Convective zCFL')
         call this%cflfile%add_column(this%fs%CFLv_x,'Viscous xCFL')
         call this%cflfile%add_column(this%fs%CFLv_y,'Viscous yCFL')
         call this%cflfile%add_column(this%fs%CFLv_z,'Viscous zCFL')
         call this%cflfile%write()
      end block create_monitor
      
      
   end subroutine init
   
   
   !> Take one time step
   subroutine step(this)
      implicit none
      class(periodicchannel), intent(inout) :: this
      
      ! Increment time
      call this%fs%get_cfl(this%time%dt,this%time%cfl)
      call this%time%adjust_dt()
      call this%time%increment()
      
      ! Remember old velocity
      this%fs%Uold=this%fs%U
      this%fs%Vold=this%fs%V
      this%fs%Wold=this%fs%W
      
      ! Turbulence modeling
      if (this%use_sgs) then
         sgs_modeling: block
            use sgsmodel_class, only: vreman
            this%resU=this%fs%rho
            call this%fs%get_gradu(this%gradU)
            call this%sgs%get_visc(type=vreman,dt=this%time%dtold,rho=this%resU,gradu=this%gradU)
            this%fs%visc=this%visc+this%sgs%visc
         end block sgs_modeling
      end if
      
      ! Perform sub-iterations
      do while (this%time%it.le.this%time%itmax)
         
         ! Build mid-time velocity
         this%fs%U=0.5_WP*(this%fs%U+this%fs%Uold)
         this%fs%V=0.5_WP*(this%fs%V+this%fs%Vold)
         this%fs%W=0.5_WP*(this%fs%W+this%fs%Wold)
         
         ! Explicit calculation of drho*u/dt from NS
         call this%fs%get_dmomdt(this%resU,this%resV,this%resW)
         
         ! Assemble explicit residual
         this%resU=-2.0_WP*(this%fs%rho*this%fs%U-this%fs%rho*this%fs%Uold)+this%time%dt*this%resU
         this%resV=-2.0_WP*(this%fs%rho*this%fs%V-this%fs%rho*this%fs%Vold)+this%time%dt*this%resV
         this%resW=-2.0_WP*(this%fs%rho*this%fs%W-this%fs%rho*this%fs%Wold)+this%time%dt*this%resW
         
         ! Add body forcing
         forcing: block
            use mpi_f08,  only: MPI_SUM,MPI_ALLREDUCE,MPI_IN_PLACE
            use parallel, only: MPI_REAL_WP
            integer :: i,j,k,ierr
            real(WP) :: Uvol,meanUold
            real(WP), parameter :: coeff=1.0e-2_WP
            ! Compute MFRs at t^n and t^(n+1)
            Uvol=0.0_WP; this%meanU=0.0_WP; meanUold=0.0_WP
            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                     this%meanU=this%meanU+this%fs%cfg%dxm(i)*this%fs%cfg%dy(j)*this%fs%cfg%dz(k)*(2.0_WP*this%fs%U(i,j,k)-this%fs%Uold(i,j,k))
                     meanUold  =meanUold  +this%fs%cfg%dxm(i)*this%fs%cfg%dy(j)*this%fs%cfg%dz(k)*this%fs%Uold(i,j,k)
                     Uvol      =Uvol      +this%fs%cfg%dxm(i)*this%fs%cfg%dy(j)*this%fs%cfg%dz(k)
                  end do
               end do
            end do
            call MPI_ALLREDUCE(MPI_IN_PLACE,Uvol      ,1,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,this%meanU,1,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr); this%meanU=this%meanU/Uvol
            call MPI_ALLREDUCE(MPI_IN_PLACE,meanUold  ,1,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr); meanUold  =meanUold/Uvol
            ! Adjust bodyforce
            this%bforce=this%bforce+coeff*(this%fs%rho*(1.0_WP-this%meanU)/this%time%dt-this%fs%rho*(this%meanU-meanUold)/this%time%dtold)
            where (this%fs%umask.eq.0) this%resU=this%resU+this%time%dt*this%bforce
         end block forcing
         
         ! Form implicit residuals
         call this%fs%solve_implicit(this%time%dt,this%resU,this%resV,this%resW)
         
         ! Apply these residuals
         this%fs%U=2.0_WP*this%fs%U-this%fs%Uold+this%resU
         this%fs%V=2.0_WP*this%fs%V-this%fs%Vold+this%resV
         this%fs%W=2.0_WP*this%fs%W-this%fs%Wold+this%resW
         
         ! Apply other boundary conditions on the resulting fields
         call this%fs%apply_bcond(this%time%t,this%time%dt)
         
         ! Solve Poisson equation
         call this%fs%correct_mfr()
         call this%fs%get_div()
         this%fs%psolv%rhs=-this%fs%cfg%vol*this%fs%div*this%fs%rho/this%time%dt
         this%fs%psolv%sol=0.0_WP
         call this%fs%psolv%solve()
         call this%fs%shift_p(this%fs%psolv%sol)
         
         ! Correct velocity
         call this%fs%get_pgrad(this%fs%psolv%sol,this%resU,this%resV,this%resW)
         this%fs%P=this%fs%P+this%fs%psolv%sol
         this%fs%U=this%fs%U-this%time%dt*this%resU/this%fs%rho
         this%fs%V=this%fs%V-this%time%dt*this%resV/this%fs%rho
         this%fs%W=this%fs%W-this%time%dt*this%resW/this%fs%rho
         
         ! Increment sub-iteration counter
         this%time%it=this%time%it+1
         
      end do
      
      ! Recompute interpolated velocity and divergence
      call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
      call this%fs%get_div()
      
      ! Output to ensight
      if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
      
      ! Perform and output monitoring
      call this%fs%get_max()
      call this%mfile%write()
      call this%cflfile%write()
      
      ! Finally, see if it's time to save restart files
      if (this%save_evt%occurs()) then
         save_restart: block
            use string, only: str_medium
            character(len=str_medium) :: timestamp
            ! Prefix for files
            write(timestamp,'(es12.5)') this%time%t
            ! Populate df and write it
            call this%df%push(name='t' ,val=this%time%t )
            call this%df%push(name='dt',val=this%time%dt)
            call this%df%push(name='bf',val=this%bforce )
            call this%df%push(name='U' ,var=this%fs%U   )
            call this%df%push(name='V' ,var=this%fs%V   )
            call this%df%push(name='W' ,var=this%fs%W   )
            call this%df%push(name='P' ,var=this%fs%P   )
            call this%df%write(fdata='restart/channel_'//trim(adjustl(timestamp)))
         end block save_restart
      end if
      
   end subroutine step
   
   
   !> Finalize simulation
   subroutine final(this)
      implicit none
      class(periodicchannel), intent(inout) :: this
      
      ! Deallocate work arrays
      deallocate(this%resU,this%resV,this%resW,this%Ui,this%Vi,this%Wi)
      if (this%use_sgs) deallocate(this%gradU)
      
   end subroutine final
   
   
end module periodicchannel_class