module init

	use modPM3D
	use random
	implicit none

contains

	function randn(N) result(x)
		integer, intent(in) :: N
		real(mp) :: x(N)
		integer :: i

		do i = 1,N
!			x(i) = SQRT(-2.0_mp*LOG(RAND()))*COS(2.0_mp*pi*RAND())
			x(i) = random_normal()
		end do
	end function

	subroutine twostream_initialize(this,Nd,v0)			!generate initial distribution
		type(PM3D), intent(inout) :: this
		integer, intent(in) :: Nd(3)
		real(mp), intent(in) :: v0
		real(mp) :: xp0(PRODUCT(Nd),3), vp0(PRODUCT(Nd),3), qs,ms,spwt
		real(mp), dimension(this%m%ng(1),this%m%ng(2),this%m%ng(3)) :: rho_back
		real(mp) :: L(3),qe,me,mode
		integer :: N,i,i1,i2,i3,pm(PRODUCT(Nd))
		integer :: nseed,clock
		integer, allocatable :: seed(:)

		L = this%L
		N = PRODUCT(Nd)
		mode=1.0_mp

		qe = -this%wp*this%wp/(N/PRODUCT(L))
		qs = qe
		me = -qe
		ms = me
		spwt = 1.0_mp
		rho_back = -qe*N/PRODUCT(L)

		call buildSpecies(this%p(1),N,qs,ms,spwt)
		call setMesh(this%m,rho_back)

		!spatial distribution initialize
!		do i3 = 1,Nd(3)
!			do i2 = 1,Nd(2)
!				do i1 = 1,Nd(1)
!					xp0(i1+Nd(1)*(i2-1)+Nd(1)*Nd(2)*(i3-1),:) =	&
!								(/ (i1-1)*L(1)/Nd(1),(i2-1)*L(2)/Nd(2),(i3-1)*L(3)/Nd(3) /)
!				end do
!			end do
!		end do

		call RANDOM_SEED(size=nseed)
		allocate(seed(nseed))
!		call SYSTEM_CLOCK(COUNT=clock)
!		seed = clock + 127*(/ ( i, i=1,nseed ) /)
		seed = 127*(/ ( i, i=1,nseed ) /)
		call RANDOM_SEED(put=seed)
		deallocate(seed)

		call RANDOM_NUMBER(xp0)
		do i1=1,3
			xp0(:,i1) = xp0(:,i1)*L(i1)
		end do
		xp0(:,1) = xp0(:,1) - this%A0(1)*L(1)/N*SIN( 2.0_mp*pi*xp0(:,1)/L(1)*mode ) + 0.5_mp*L(1)

		!velocity distribution initialize
!		vp0 = vT*randn(N)
		vp0 = 0.0_mp
		pm = (/ ( i, i=1,N ) /)
		pm = 1 - 2*MOD(pm,2)
		vp0(:,1) = vp0(:,1) + pm*v0

		call setSpecies(this%p(1),N,xp0,vp0)
	end subroutine

end module
