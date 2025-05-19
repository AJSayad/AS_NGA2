!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use sml_class,         only: sml
   use ml_class,          only: ml
   use coupler_class,     only: coupler
   use ensight_class,     only: ensight
   use event_class,       only: event
   use mast_class,        only: mast
   use matm_class,        only: matm
   implicit none
   private
   
   !> SML simulation
   type(sml) :: mixing

   !> Mixing Layer simulation
   type(ml)  :: shear
   
   !> Couplers from injector to atomization
   type(coupler) :: xcpl_i2a,ycpl_i2a,zcpl_i2a

   public :: simulation_init,simulation_run,simulation_final
   
contains
   
   
   !> Initialization of our simulation
   subroutine simulation_init
      implicit none
      
      ! Initialize mixing layer that generates the turbulent field
      if (.not.shear%restarted) then
         call mixing%init()
      end if

      ! Initialize the actual shear layer simulation I want to run
      call shear%init()

      ! ! Initialize coupler: this allows for a proper and generalized coupling between "mixing" and "shear" cases
      ! create_coupler_i2a: block
      !    use parallel, only: group
      !    xcpl_i2a=coupler(src_grp=group,dst_grp=group,name='sml_to_ml'); call xcpl_i2a%set_src(mixing%cfg,'x'); call xcpl_i2a%set_dst(shear%cfg,'x'); call xcpl_i2a%initialize()
      !    ycpl_i2a=coupler(src_grp=group,dst_grp=group,name='sml_to_ml'); call ycpl_i2a%set_src(mixing%cfg,'y'); call ycpl_i2a%set_dst(shear%cfg,'y'); call ycpl_i2a%initialize()
      !    zcpl_i2a=coupler(src_grp=group,dst_grp=group,name='sml_to_ml'); call zcpl_i2a%set_src(mixing%cfg,'z'); call zcpl_i2a%set_dst(shear%cfg,'z'); call zcpl_i2a%initialize()
      ! end block create_coupler_i2a

      if (.not.shear%restarted) then
         ! Perform time integration for fluctuations first
         do while (mixing%time%t < 125) 
   
            ! Advance mixing layer simulation - this is where I run my first simulation by calling this "mixing" class
            call mixing%step()
   
         end do

         ! Once I finish running "mixing", I want to add the resulting velocity and pressure fluctuations into my initial conditions 
         ! of the actual shear layer case. Because both domains happen to be the exact same, I can get away with a simple transfer
         ! based on the indices. However, if the domains are different or the setup is more complicated, then I would have to use 
         ! the more general coupler class. 
         couple_simple: block
         use param,           only: param_read
         integer :: n,i,j,k,ierr
         real(WP) :: GP,LP,gamm_g,gamm_l,Grho,Lrho,Ma_g,r_rho
         ! Interpolate velocities to cell centers
         call mixing%fs%interp_vel(mixing%Ui,mixing%Vi,mixing%Wi)
         ! Add fluctuation field onto the actual mixing layer (compressible multiphase)
         do k=shear%fs%cfg%kmino_,shear%fs%cfg%kmaxo_
            do j=shear%fs%cfg%jmino_,shear%fs%cfg%jmaxo_
               do i=shear%fs%cfg%imino_,shear%fs%cfg%imaxo_
                  if (shear%vf%VF(i,j,k).eq.0.0_WP) then 
                     shear%fs%Ui(i,j,k)=shear%fs%Ui(i,j,k)+mixing%Ui(i,j,k)
                     shear%fs%Vi(i,j,k)=shear%fs%Vi(i,j,k)+mixing%Vi(i,j,k)
                     shear%fs%Wi(i,j,k)=shear%fs%Wi(i,j,k)+mixing%Wi(i,j,k)
                     shear%fs%P (i,j,k)=shear%fs%P (i,j,k)+mixing%fs%P(i,j,k) 
                  end if
               end do
            end do
         end do
   
         end block couple_simple
   
         ! ! Handle coupling between "mixing" and "shear" here. This would be the general structure for using the coupler class
         ! ! However, the ligament class is likely a better example if you need to figure out more info on the coupler class 
         ! ! for more complicated situations.
         ! coupling_i2a: block
         ! integer :: n,i,j,k,ierr
         ! ! Interpolate velocities to cell centers
         ! call mixing%fs%interp_vel(mixing%Ui,mixing%Vi,mixing%Wi)
         ! ! Exchange data using cpl12x/y/z couplers
         ! call xcpl_i2a%push(mixing%Ui); call xcpl_i2a%transfer(); call xcpl_i2a%pull(shear%recU)
         ! call ycpl_i2a%push(mixing%Vi); call ycpl_i2a%transfer(); call ycpl_i2a%pull(shear%recV)
         ! call zcpl_i2a%push(mixing%Wi); call zcpl_i2a%transfer(); call zcpl_i2a%pull(shear%recW)
         ! ! Add fluctuation field onto the actual mixing layer
         ! do k=shear%fs%cfg%kmino_,shear%fs%cfg%kmaxo_
         !    do j=shear%fs%cfg%jmino_,shear%fs%cfg%jmaxo_
         !       do i=shear%fs%cfg%imino_,shear%fs%cfg%imaxo_
         !          shear%fs%Ui(i,j,k)=shear%fs%Ui(i,j,k)+shear%recU(i,j,k)
         !          shear%fs%Vi(i,j,k)=shear%fs%Vi(i,j,k)+shear%recV(i,j,k)
         !          shear%fs%Wi(i,j,k)=shear%fs%Wi(i,j,k)+shear%recW(i,j,k)
         !       end do
         !    end do
         ! end do
         ! ! Need something like MPI_Barrier
         ! call MPI_Barrier(shear%fs%cfg%comm,ierr)
         ! call shear%cfg%sync(shear%fs%Ui)
         ! call shear%cfg%sync(shear%fs%Vi)
         ! call shear%cfg%sync(shear%fs%Wi)
         ! end block coupling_i2a
   
         end if
      
   end subroutine simulation_init
   
   
   !> Run the simulation - here, we now call the "shear" class that houses the main simulation I want to run
   subroutine simulation_run
      use param,           only: param_read
      real(WP) :: max_time
      ! implicit none
      
      call param_read('Max time',max_time)
      ! Perform time integration for actual mixing layer simulation
      do while (shear%time%t < max_time)

         ! Advance mixing layer simulation
         call shear%step()

      end do
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none

      ! Finalize fluctuation simulation
      call mixing%final()

      ! Finalize mixing layer simulation
      call shear%final()
   end subroutine simulation_final
   
   
end module simulation
