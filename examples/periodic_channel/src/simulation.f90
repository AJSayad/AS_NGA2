!> Various definitions and tools for running an NGA2 simulation
module simulation
   use periodicchannel_class, only: periodicchannel
   implicit none
   private
   
   !> Periodic channel simulation
   type(periodicchannel) :: channel
   
   public :: simulation_init,simulation_run,simulation_final
   
contains
   
   
   !> Initialization of our simulation
   subroutine simulation_init
      implicit none
      ! Initialize periodic channel simulation
      call channel%init()
   end subroutine simulation_init
   
   
   !> Run the simulation
   subroutine simulation_run
      implicit none
      ! Time-integrate channel simulation
      do while (.not.channel%time%done())
         call channel%step()
      end do
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      ! Finalize channel simulation
      call channel%final()
   end subroutine simulation_final
   
   
end module simulation
