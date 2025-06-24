!> Various definitions and tools for running an NGA2 simulation
module simulation
  use precision,         only: WP
  use shockgen_class,    only: sgen
  use shockdrop_class,   only: sdrop
  use ensight_class,     only: ensight
  use event_class,       only: event
  use mast_class,        only: mast
  use matm_class,        only: matm
  use param,             only: param_read
  implicit none
  private

  !> shock generator simulation
  type(sgen) ::  shockgen

  !> shock droplet simulation
  type(sdrop) :: shockdrop

  public :: simulation_init,simulation_run,simulation_final

contains

  !> initialize full simulation
  subroutine simulation_init
    
    ! initialize the shock generator sim
    if (.not.shockdrop%restarted)then
       call shockgen%init()
    end if

    ! initialize coupler if needed
    if (.not.shockdrop%restarted)then
       do while (.not.shockgen%time%done())
          ! advance shock generator sim by one step
          call shockgen%step()
       end do
       call shockgen%final()
    end if

    ! initialize shock droplet sim
    ! note: restart logic is still built into shockdrop%init subroutine
    call shockdrop%init(shockgen%time%dt,shockgen%time%dtmax)

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
               shockdrop%fs%Grho(i,:,:)  = shockgen%Grho_profile(i+n_shock-shock_index+1)
               shockdrop%fs%GP(i,:,:)    = shockgen%GP_profile(i+n_shock-shock_index+1)
               shockdrop%fs%GrhoE(i,:,:) = shockgen%GrhoE_profile(i+n_shock-shock_index+1)
               shockdrop%fs%Ui(i,:,:)    = shockgen%Ui_profile(i+n_shock-shock_index+1)
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
    call shockdrop%final()
  end subroutine simulation_final
  
end module simulation
