!module 
!
!contains


subroutine NEWTONLINODESOLVE(spec0,nspec,Temp,press,dtime)

!==================================================================================================================================!
!
!
!       this is a dummy routine
!       needed to build dummy library libNewtonLinODE.so
!
!==================================================================================================================================!

  implicit none

  integer,intent(in) :: nspec

  double precision,intent(in) :: dtime
  double precision :: spec0(1:nspec),Temp,press
! ---

!----------------------------------------------------------------------------------------------------------------------------------!
  write(*,*) 'This is a dummy routine for the chemistry integration by Newton Linearization - NEWTONLINODESOLVE'
  write(*,*) 'source code not provided - use different option for integration of chemical mechanism'
  write(*,*) 'modify ./control/cloudproperites'
  stop

  return
!----------------------------------------------------------------------------------------------------------------------------------!


end subroutine

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

subroutine NEWTONLINODEINIT
!==================================================================================================================================!
!
!
!       this is a dummy routine
!       needed to build dummy library libNewtonLinODE.so
!
!==================================================================================================================================!
  implicit none

!----------------------------------------------------------------------------------------------------------------------------------!
  write(*,*) 'This is a dummy routine for the chemistry integration by Newton Linearization - NEWTONLINODESOLVE'
  write(*,*) 'source code not provided - use different option for integration of chemical mechanism'
  write(*,*) 'modify ./control/cloudproperites'
  stop

  return
!----------------------------------------------------------------------------------------------------------------------------------!


end subroutine
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

