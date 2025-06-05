!> Various definitions and tools for running an NGA2 simulation
module simulation
  use precision,         only: WP
  use shockgen_class,    only: sgen
  use shockdrop_class,   only: sdrop
  use ensight_class,     only: ensight
  use event_class,       only: event
  use mast_class,        only: mast
  use matm_class,        only: matm
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
    !if (.not.shockdrop%restarted)then
    print*, "simulation.f90: calling shockgen init."
    call shockgen%init()
    !end if

    ! initialize shock droplet sim
    !call shockdrop%init()

    ! initialize coupler

    !if (.not.shockdrop%restarted)then
    do while (.not.shockgen%time%done())
       ! advance shock generator sim by one step
       print*, "simulation.f90: calling shockgen step."
       print*, "simulation.f90 shockgen%time%done(): ", shockgen%time%done()
       call shockgen%step()
    end do
    call shockgen%final()
    print*, "SIMULATION.F90"
    print*, "shockgen%Grho_profile: ", shockgen%Grho_profile
    print*, "shockgen%GrhoE_profile: ", shockgen%GrhoE_profile
    print*, "shockgen%GP_profile: ", shockgen%GP_profile
    print*, "shockgen%Ui_profile: ", shockgen%Ui_profile
    !end if

       ! add coupling block
       shock_profile: block
         integer :: i
         
         ! find shock index
         if (.not.shockdrop%restarted)then
            do i=shockdrop%cfg%imin,shockdrop%cfg%imax
               if ((shockdrop%cfg%xm(i).lt.(shockdrop%xshock+shockdrop%tol)).and.(shockdrop%cfg%xm(i).gt.(shockdrop%xshock-shockdrop%tol))) then
                  shockdrop%shock_index=i
                  shockdrop%shock_loc = shockdrop%cfg%xm(shockdrop%shock_index) ! store location of shock corresponding to index
                  if(shockdrop%cfg%amRoot)then
                     print*, "The shock has been found at index: ", i
                     print*, "The found shock location (cell center) is: ", shockdrop%shock_loc
                  end if
               end if
            end do
            
            ! read in numerical shock profile 
            do i=shockdrop%cfg%imino_,shockdrop%cfg%imaxo_ 
               if ((shockdrop%cfg%xm(i).le.shockdrop%cfg%xm(shockdrop%shock_index+shockdrop%n_shock)).and.(shockdrop%cfg%xm(i).ge.shockdrop%cfg%xm(shockdrop%shock_index-shockdrop%n_shock)))then
                  shockdrop%fs%Grho(i,:,:) = shockgen%Grho_profile(i+shockdrop%n_shock-shockdrop%shock_index+1)
                  shockdrop%fs%GP(i,:,:) = shockgen%GP_profile(i+shockdrop%n_shock-shockdrop%shock_index+1)
                  shockdrop%fs%GrhoE(i,:,:) = shockgen%GrhoE_profile(i+shockdrop%n_shock-shockdrop%shock_index+1)
                  shockdrop%fs%Ui(i,:,:) = shockgen%Ui_profile(i+shockdrop%n_shock-shockdrop%shock_index+1)
               end if
            end do
         end if

         !update timestep from shockgen to shockdrop
         shockdrop%time%dtmax = shockgen%time%dtmax
       end block shock_profile
         
     end subroutine simulation_init

     !> run full simulation
     subroutine simulation_run

     end subroutine simulation_run
     
     !> finalize simulation
     subroutine simulation_final

     end subroutine simulation_final

end module simulation
