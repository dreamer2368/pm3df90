module modPM3D

	use modPlasma
	use modMesh
	use modAssign
	use modRecord

	implicit none

	type PM3D
		integer :: nt, ni, n, ng(3)
		real(mp) :: L(3), eps0, wp
		real(mp) :: dt, A0, B0

		type(plasma) :: p
		type(mesh) :: m
		type(recordData) :: r
		type(pmAssign) :: a
	end type

	! parameters setup
!	integer, parameter :: Nd(3) = (/ 2**0, 2**7, 2**7 /)

!	real(mp) :: Ti = 40.0_mp
!	real(mp) :: Tp = 2.0_mp*pi/wp
!	real(mp) :: T = 40.0_mp
!	real(mp) :: T = 40.0_mp + 10.0_mp*(2.0_mp*pi/wp)

	!initial spatial distribution
	integer, parameter :: mode = 1

	!initial velocity distribution
	real(mp), parameter :: vT = 0.0_mp
	real(mp), parameter :: v0 = 0.0_mp
!	real(mp) :: vp0(N,3)
!	integer :: pm(N)

	!grid and operators setup
!	integer, parameter :: Ng(3) = (/ 64, 64, 64 /)

	!error convergence variable
	real(mp) :: fDA(30)
	real(mp) :: ek(30)
	real(mp) :: e(9,30)

contains

	subroutine buildPM3D(this,Tf,Ti,Ng,N,dt,L,A,B)
		type(PM3D), intent(out) :: this
		real(mp), intent(in) :: Tf,Ti
		integer, intent(in) :: Ng(3), N
		real(mp), intent(in), optional :: dt, A, B, L(3)
		if( present(dt) ) then
			this%dt = dt
		else
			this%dt = 0.2_mp
		end if
		if( present(A) ) then
			this%A0 = A
		else
			this%A0 = 1.0_mp
		end if
		if( present(B) ) then
			this%B0 = B
		else
			this%B0 = 0.0_mp
		end if
		if( present(L) ) then
			this%L = L
		else
			this%L = 2*pi/( sqrt(3.0_mp)/2.0_mp/sqrt(2.0_mp)/0.2_mp )
		end if
		this%ng = Ng
		this%n = N
		this%nt = CEILING(Tf/this%dt)
		this%dt = Tf/this%nt
		this%ni = FLOOR(Ti/this%dt) + 1
		print *, 'Ni = ',this%ni,', Nt = ',this%nt,', dt = ',this%dt

		this%eps0 = 1.0_mp
		this%wp = 1.0_mp
		call buildPlasma(this%p,N)
		call buildMesh(this%m,this%L,Ng)
		call buildAssign(this%a,N,Ng)
		call buildRecord(this%r,this%nt,this%n,this%L,this%ng)
	end subroutine

	subroutine destroyPM3D(this)
		type(PM3D), intent(inout) :: this

		call destroyPlasma(this%p)
		call destroyMesh(this%m)
		call destroyAssign(this%a)
		call destroyRecord(this%r)
	end subroutine

end module