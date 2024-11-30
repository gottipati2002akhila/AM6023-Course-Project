!--------------------------------------------------------
!--------------------------------------------------------
!   etica: Solution to the elastica with boundary force
!   By - 
!--------------------------------------------------------
!--------------------------------------------------------

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      REAL(8), PARAMETER :: PI = 4 * atan(1.0)
      DOUBLE PRECISION theta, thetap, nx, ny, x, y, S
      DOUBLE PRECISION eta
      
      	nx =U(1)
      	ny =U(2)
      	theta = U(3)
      	thetap =U(4)
      	x =U(5)
      	y =U(6)
      	eta =7.00
      	F(1) =-eta**4*y*SIN(theta)
      	F(2) = eta**4*y*COS(theta)
      	F(3) = thetap
      	F(4) = nx*SIN(theta) - ny*COS(theta)
      	F(5) = COS(theta)
      	F(6) = SIN(theta) 
  
      	
      END SUBROUTINE FUNC
!----------------------------------------------------------------------
! setting initial guess to be zero 
! ----------------------------------------------------
      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- -----

      IMPLICIT NONE
      REAL(8), PARAMETER :: PI = 4 * atan(1.0)
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T
	
       PAR(1) =1
       U(1)=0
       U(2)=0
       U(3)=0
       U(4)=0
       U(5)=0
       U(6)=0

       
      END SUBROUTINE STPNT
!----------------------------------------------------------------------

      SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC)
!     ---------- ----
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
      DOUBLE PRECISION, INTENT(IN) :: PAR(*), U0(NDIM), U1(NDIM)
      DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
      DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)
       
       FB(1)=U0(3)	    ! theta(0) =0
       FB(2)=U1(3)         !theta(1)=0
       FB(3)=U0(5)          ! x(0) =0
       FB(4)=U1(5)-PAR(1)+1! x(1)-1+delta =0
       FB(5)=U0(6)          ! y(0) =0
       FB(6)=U1(6)          ! y(1) =0
      
      END SUBROUTINE BCND
      
!----------------------------------------------------------------------
      SUBROUTINE ICND(NDIM,PAR,ICP,NINT,U,UOLD,UDOT,UPOLD,FI,IJAC,DINT)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), NINT, IJAC
      DOUBLE PRECISION, INTENT(IN) :: PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), UOLD(NDIM), UDOT(NDIM), UPOLD(NDIM)
      DOUBLE PRECISION, INTENT(OUT) :: FI(NINT)
      DOUBLE PRECISION, INTENT(INOUT) :: DINT(NINT,*)

     

      END SUBROUTINE ICND

      SUBROUTINE FOPT
      END SUBROUTINE FOPT
      
      ! ----------------------------------------------------------------------

      SUBROUTINE PVLS(NDIM,U,PAR)
      !     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)

      DOUBLE PRECISION, EXTERNAL :: GETP,GETU2
      INTEGER NDX,NCOL,NTST
      DOUBLE PRECISION :: s
      
	PAR(2) =GETP('BV1',1,U)
	PAR(3) =GETP('BV0',1,U)

      END SUBROUTINE PVLS

!----------------------------------------------------------------------

      
      
      

