!> Definition for an sgen class (shock generator)
module sgen_class
  use precision,         only: WP
  use geometry,          only: cfg
  use mast_class,        only: mast
  use vfs_class,         only: vfs
  use matm_class,        only: matm
  use timetracker_class, only: timetracker
  use ensight_class,     only: ensight
  use event_class,       only: event
  use monitor_class,     only: monitor
  use hypre_str_class,   only: hypre_str
  use pardata_class,     only: pardata
  implicit none
  private

  public :: sgen

  !> sgen object
  type :: sgen
     !> Config
     type(config) :: cfg !> Mesh for solver
     !> Flow solver
     type(mast),        :: fs
     type(vfs),         :: vf
     type(matm),        :: matmod
     type(timetracker), :: time
     type(hypre_str),   :: ps
     type(hypre_str),   :: vs
     !> Ensight postprocessing
     type(ensight)      :: ens_out
     type(event)        :: ens_evt
     !> Simulation monitor file
     type(monitor)      :: mfile,cflfile,cvgfile
     !> Fluid parameters
     real(WP)           :: visc !AS is this needed? I borrowed this from Chase's shear layer pre-sim sml_class.f90
   contains
     procedure :: init  !> initialize sgen simulation
     procedure :: step  !> advance sgen simulation by one timestep
     procedure :: final !> finalize sgen simulation
  end type sgen


contains

  !> Initialization of the shock generator (sgen) simulation



  !> Take one time step with specified dt


  !> Finalize shock generator (sgen) simulation
