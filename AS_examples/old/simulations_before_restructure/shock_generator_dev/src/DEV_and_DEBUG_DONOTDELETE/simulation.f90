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
  !use coupler_class,     only: coupler ! AS used for coupling sims with two different domains, see Chase's notes
  implicit none
  private

  !> shock generator simulation
  type(sgen) ::  shockgen

  !> shock droplet simulation
  type(sdrop) :: shockdrop

  !> couplers from shockgen to shockdrop
  !type(coupler) ::

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
    if (.not.shockdrop%restarted)then
       call shockdrop%init(shockgen%time%dt,shockgen%time%dtmax)
    end if

    ! add coupling block if needed
    shock_profile: block
      integer :: i,n_shock,shock_index
      real(WP) :: tol,Lx,dx,shock_loc
      call param_read('Lx',Lx)
      call param_read('n_shock',n_shock)
      dx = Lx/shockdrop%cfg%nx
      tol = dx/2
      
      ! find shock index
      if (.not.shockdrop%restarted)then
         do i=shockdrop%cfg%imin,shockdrop%cfg%imax
            if ((shockdrop%cfg%xm(i).lt.(shockdrop%xshock+tol)).and.(shockdrop%cfg%xm(i).gt.(shockdrop%xshock-tol))) then
               shock_index=i
               shock_loc = shockdrop%cfg%xm(shock_index) ! store location of shock corresponding to index
               print*, "sim.f90: The shock has been found at index: ", i
               print*, "sim.f90: The found shock location (cell center) is: ", shock_loc
            end if
         end do
         
         !print*, "shock index    = ", shock_index
         !print*, "shock location = ", shock_loc
         !print*, "n_shock        = ", n_shock
         !print*, "tol            = ", tol
         !print*, "dx             = ", dx
         !print*, "shock_index+n_shock = ", shock_index+n_shock
         !print*, "shock_index-n_shock = ", shock_index-n_shock
         !print*, "xm(shock_index+n_shock) = ", shockdrop%cfg%xm(shock_index+n_shock)
         !print*, "xm(shock_index-n_shock) = ", shockdrop%cfg%xm(shock_index-n_shock)
         
         ! read in numerical shock profile 
         do i=shockdrop%cfg%imin_,shockdrop%cfg%imax_
            if ((shockdrop%cfg%xm(i).le.shockdrop%cfg%xm(shock_index+n_shock)).and.(shockdrop%cfg%xm(i).ge.shockdrop%cfg%xm(shock_index-n_shock)))then
               shockdrop%fs%Grho(i,:,:)  = shockgen%Grho_profile(i+n_shock-shock_index+1)
               shockdrop%fs%GP(i,:,:)    = shockgen%GP_profile(i+n_shock-shock_index+1)
               shockdrop%fs%GrhoE(i,:,:) = shockgen%GrhoE_profile(i+n_shock-shock_index+1)
               shockdrop%fs%Ui(i,:,:)    = shockgen%Ui_profile(i+n_shock-shock_index+1)
               
               !print*, "profile index    : ", i+n_shock-shock_index+1
               !print*, "field array index: ", i
               !print*, "Profile Grho : ", shockgen%Grho_profile(i+n_shock-shock_index+1)
               !print*, "Field Grho   : ", shockdrop%fs%Grho(i,50,1)
               !print*, "Profile GP   : ", shockgen%GP_profile(i+n_shock-shock_index+1)
               !print*, "Field GP     : ", shockdrop%fs%GP(i,50,1)
               !print*, "Profile GrhoE: ", shockgen%GrhoE_profile(i+n_shock-shock_index+1)
               !print*, "Field GrhoE  : ", shockdrop%fs%GrhoE(i,50,1)
               !print*, "Profile Ui   : ", shockgen%Ui_profile(i+n_shock-shock_index+1)
               !print*, "Field Ui     : ", shockdrop%fs%Ui(i,50,1)
            end if
         end do
         
         !print*, "CENTERLINE Field Grho     : ", shockdrop%fs%Grho(:,50,1)
         !print*, "CENTERLINE Field GrhoE     : ", shockdrop%fs%GrhoE(:,50,1)
         !print*, "CENTERLINE Field GP        : ", shockdrop%fs%GP(:,50,1)
         !print*, "CENTERLINE Field Ui        : ", shockdrop%fs%Ui(:,50,1)
         
         call shockdrop%update_mixture_variables()
         call shockdrop%writeIC()
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
