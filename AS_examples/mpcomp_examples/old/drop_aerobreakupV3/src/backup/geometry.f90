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
      use param,       only: param_read
      implicit none
      type(sgrid) :: grid
      
      ! Create a grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
         integer  :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz 
         real(WP) :: ddrop,CPD,D0X,D0Y,D0Z
         real(WP) :: dx,dy,dz
         logical  :: xper,yper,zper
         real(WP), dimension(:), allocatable :: x,y,z

         ! AS note: the origin (0,0,0) is the 'back-left corner' (from default paraview view x+ right, y+ up, z+ pointing out of screen)     
         ! the case is fully non-periodic but we assume a 2D sim (hence zper) until we read input file
         xper=.false.; yper=.false.; zper=.true.

         ! Read in grid definition
         call param_read('Drop diameter',ddrop)
         call param_read('CPD',CPD)
         call param_read('D0X',D0X)
         call param_read('D0Y',D0Y)
         call param_read('D0Z',D0Z)

         ! compute physical domain length
         Lx = D0X*ddrop; Ly = D0Y*ddrop                                 ! compute domain lengths
         nx = ceiling((CPD*Lx)/ddrop); ny = ceiling((CPD*Ly)/ddrop)     ! compute number of uniform cells
         dx = Lx/nx; dy = Ly/ny                                         ! uniform mesh spacing
         if (D0Z.gt.0) then ! handle 3D case
            Lz = D0Z*ddrop; nz = ceiling((CPD*Lz)/ddrop); dz = Lz/nz; zper = .false.;
         else
            nz = 1; Lz = dx; dz = dx
         end if 

         ! allocate arrays
         allocate(x(nx+1));allocate(y(ny+1));allocate(z(nz+1))

         ! create simple uniform rectilinear mesh
         do i=1,nx+1; x(i) = real(i-1,WP)*dx;end do
         do j=1,ny+1; y(j) = real(j-1,WP)*dy;end do
         if (D0Z.gt.0)then
               do k=1,nz+1; z(k) = real(k-1,WP)*dz;end do
            else
               z(1) = -dz; z(2) = dz ! 2D mesh
         end if
         
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=xper,yper=yper,zper=zper,name='ShockDrop')
         
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
