!> Various definitions and tools for initializing NGA2 config
module geometry
   use config_class, only: config
   use precision,    only: WP
   implicit none
   private

   !> Single config
   type(config), public :: cfg
   public :: geometry_init

contains

  !> Initialization of problem geometry
   subroutine geometry_init
      use sgrid_class, only: sgrid
      use param,       only: param_read, param_exists
      use parallel,    only: amRoot
      use messager,    only: die
      implicit none
      type(sgrid) :: grid

      ! Create a grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
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
         
         ! variable for shock profile extraction
         logical :: extract_flag

         call param_read('Profile extraction flag',extract_flag)
         
         ! Read in grid definition
         call param_read('Lx',Lx); call param_read('Lx ref', start_ref, default=0.0_WP);
         
         call param_read('nx',nx); call param_read('nx stretch left',nx_stretchL); call param_read('nx stretch right',nx_stretchR);
         allocate(x(nx+nx_stretchL+nx_stretchR+1));

         dx = Lx/nx
         dx_ref = (Lx - start_ref)/real(nx,WP)

         ! if profile is being extracted, run the simulation in 1D
         if (extract_flag.eqv.(.true.)) then
            call param_read('Ly',Ly);
            ny = 1 !if singlephase, run in 1D
            nz = 1
            Lz = dx
            allocate(y(ny+1))
            allocate(z(nz+1))
         else if (extract_flag.eqv.(.false.)) then !if multiphase, run sim in 2D or 3D
            call param_read('Ly',Ly); call param_read('ny',ny); call param_read('ny stretch',ny_stretch);
            allocate(y(ny+2*ny_stretch+1))

            call param_read('nz',nz,default=1)
            call param_read('nz stretch',nz_stretch)
            allocate(z(nz+2*nz_stretch+1))

            if (nz.eq.1) then
               Lz = dx
            else
               call param_read('Lz',Lz)
            end if
         end if
        
         ! Read in droplet information
         call param_read('Droplet diameter',ddrop)
         
            !uniform mesh x
            alpha=1.05_WP
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

            ! mesh in y
            if (extract_flag.eqv.(.true.)) then
               ! uniform 1 cell mesh in y for 1D sim
               do j=1,ny+1
                  ny_stretch = 0 ! 1D
                  dy = Ly/ny
                  y(j) = real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
               end do

               ! uniform 1 cell mesh in z for 1D sim ***
               do k=1,nz+1
                  nz_stretch = 0 ! 1D
                  dz = Lz/nz
                  z(k) = real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
               end do
            else
               do j=1,ny+2*ny_stretch+1 !initialize y mesh array
                  y(j) = 0.0_WP
               end do

               dy = Ly/ny ! define uniform grid spacing
               y(ny/2+ny_stretch+1) = 0.0_WP !define the centerline of the domain

               !y array uniform region
               do j=ny/2+ny_stretch+2,ny+ny_stretch+1
                  y(j) = y(j-1) + dy
               end do

               !stretching in y
               do j=ny+ny_stretch+2,ny+(2*ny_stretch)+1
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
            end if

            if(amRoot)then
               print*, "======== MESH DESCRIPTION IN x ========"
               print*, "Uniform region length: ", Lx
               print*, "Number of cells in the uniform region: ", nx
               print*, "Number of cells added to left of domain: ", nx_stretchL
               print*, "Number of cells added to left of domain: ", nx_stretchR
               print*, "Total number of cells in domain: ", nx+nx_stretchL+nx_stretchR
               print*, "Stretching ratio in x: ", alpha
               print*, "Leftmost point (stretched region left): ", x(1)
               print*, "Start of uniform region: ", start_ref
               print*, "End of uniform region: ",x(nx+nx_stretchL+1)
               print*, "Rightmost point (stretched region right): ", x(nx_stretchL+nx+nx_stretchR+1)
               print*, "Aspect ratio in uniform region: ", dx/dy
               print*, "Number of cells per droplet diameter: ", nx*(ddrop/Lx)
               print*, "======================================="             
            end if

            if (amRoot) then
               print*, "======== MESH DESCRIPTION IN y ========"
               print*, 'Uniform region height: ', Ly
               print*, 'Stretching in y starts at: +- ', Ly/2
               print*, 'Number of cells added to the top and to the bottom: ', ny_stretch/2
               print*, 'Stretching ratio in y: ', alpha
               print*, "Aspect ratio in uniform region: ", dx/dy
               print*, "Number of cells per droplet diameter: ", ny*(ddrop/Ly)
               print*, "========================================"                 
            end if

            if (amRoot) then
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
            grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.true.,zper=.true.,name='ShockDrop')

      end block create_grid

      ! Create a config from that grid on our entire group
      create_cfg: block
         use parallel, only: group
         integer, dimension(3) :: partition

         ! Read in partition
         call param_read('Partition',partition,short='p')

         ! Create partitioned grid
         cfg=config(grp=group,decomp=partition,grid=grid)

      end block create_cfg


   end subroutine geometry_init


end module geometry
