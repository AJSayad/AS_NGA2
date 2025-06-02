!> Various definitions and tools for running an NGA2 simulation
module simulation
  use precision,         only: WP
  use shockgen_class,    only: sgen
  use shockdrop_class,   only: sdrop
  use param,             only: param_read
  implicit none
  private

  ! modul level storage for profile arrays
  real(WP),public,dimension(:), allocatable :: savedGrho_profile, savedGP_profile, savedGrhoE_profile, savedUi_profile
  real(WP) :: saved_dt,saved_dtmax
  !> shock droplet simulation
  type(sdrop) :: shockdrop

  public :: run_shock_generator,simulation_init,simulation_run,simulation_final

contains

  subroutine run_shock_generator
    shock_generator: block
      !> shock generator simulation
      type(sgen) ::  shockgen
      integer :: n_shock

      call param_read('n_shock',n_shock)
      allocate(savedGrho_profile(2*n_shock+1));savedGrho_profile   = 0.0_WP
      allocate(savedGP_profile(2*n_shock+1));savedGP_profile       = 0.0_WP
      allocate(savedGrhoE_profile(2*n_shock+1));savedGrhoE_profile = 0.0_WP
      allocate(savedUi_profile(2*n_shock+1));savedUi_profile       = 0.0_WP
      
      ! initialize the shock generator sim
      call shockgen%init()
      ! run shock generator
      do while (.not.shockgen%time%done())
         ! advance shock generator sim by one step
         call shockgen%step()
      end do
      call shockgen%final()

      ! explicitly copy variables here
      savedGrho_profile  = shockgen%Grho_profile
      savedGP_profile    = shockgen%GP_profile
      savedGrhoE_profile = shockgen%GrhoE_profile
      savedUi_profile    = shockgen%Ui_profile
      saved_dt = shockgen%time%dt
      saved_dtmax = shockgen%time%dtmax
    end block shock_generator

  end subroutine run_shock_generator

  !> initialize full simulation
  subroutine simulation_init

    if(.not.shockdrop%restarted)then
       call run_shock_generator
    end if
    
    ! initialize shock droplet sim
    ! note: restart logic is still built into shockdrop%init subroutine
    call shockdrop%init(saved_dt,saved_dtmax)
    ! add coupling block if needed
    shock_profile: block
      integer :: i,n_shock,shock_index
      real(WP) :: tol,Lx,dx,shock_loc
      call param_read('Lx',Lx)
      call param_read('n_shock',n_shock)
      dx = Lx/shockdrop%cfg%nx
      tol = dx/2

      ! 1. first loop through each proc subdomain and find physical shock locations
      ! 2. Loop again through each proc subdomain and determine which points need to be updated (based on physical values)
      shock_loc = -10.0_WP ! initialize to non-physical number 
      ! find shock index
      if (.not.shockdrop%restarted)then
         ! every processor looks for the shock
         do i=shockdrop%cfg%imin,shockdrop%cfg%imax
            if ((shockdrop%cfg%xm(i).lt.(shockdrop%xshock+tol)).and.(shockdrop%cfg%xm(i).gt.(shockdrop%xshock-tol))) then
               shock_index=i
               shock_loc = shockdrop%cfg%xm(shock_index) ! store location of shock corresponding to index
               !print*, "sim.f90: Rank: ", shockdrop%cfg%rank
               !print*, "sim.f90: The shock has been found at index: ", i
               !print*, "sim.f90: The found shock location (cell center) is: ", shock_loc
            end if
         end do

         do i=shockdrop%cfg%imin_,shockdrop%cfg%imax_
            if ((shockdrop%cfg%xm(i).ge.(shock_loc-n_shock*dx).and.(shockdrop%cfg%xm(i).le.(shock_loc+n_shock*dx))))then
               shockdrop%fs%Grho(i,:,:)  = savedGrho_profile(i+n_shock-shock_index+1)
               shockdrop%fs%GP(i,:,:)    = savedGP_profile(i+n_shock-shock_index+1)
               shockdrop%fs%GrhoE(i,:,:) = savedGrhoE_profile(i+n_shock-shock_index+1)
               shockdrop%fs%Ui(i,:,:)    = savedUi_profile(i+n_shock-shock_index+1) 
            end if
         end do
         call shockdrop%update_mixture_variables() ! update mixture density, bulkmod, and momenta
         call shockdrop%writeIC() ! write IC
      else
         call shockdrop%writeIC()
      end if
      
    end block shock_profile
  end subroutine simulation_init
  
  !> run full simulation
  subroutine simulation_run
    do while (.not.shockdrop%time%done()) 
       ! advance shock generator sim by one step
       call shockdrop%step()
    end do    
  end subroutine simulation_run
  
  !> finalize simulation
  subroutine simulation_final
    ! deallocate work arrays
    deallocate(savedGrho_profile);deallocate(savedGP_profile);deallocate(savedGrhoE_profile);deallocate(savedUi_profile)
    deallocate(saved_dt);deallocated(saved_dtmax)
    
    call shockdrop%final()
  end subroutine simulation_final
  
end module simulation
