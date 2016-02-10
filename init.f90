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
			x(i) = SQRT(-2.0_mp*LOG(RAND()))*COS(2.0_mp*pi*RAND())
		end do
	end function

	subroutine particle_initialize(this,Nd,v0,xp0,vp0,qs,ms,rho_back)			!generate initial distribution
		type(PM3D), intent(inout) :: this
		integer, intent(in) :: Nd(3)
		real(mp), intent(in) :: v0
		real(mp), intent(out) :: xp0(this%n,3), vp0(this%n,3), qs(this%n),ms(this%n),rho_back
		real(mp) :: L(3),qe,me
		integer :: i,i1,i2,i3,pm(this%n)
		L = this%L

		if( this%n .ne. PRODUCT(Nd) ) then
			print *, '!!FAULT :: Nd array size is not equal to the number of particles.'
			stop
		end if

		qe = -this%wp*this%wp/(this%n/PRODUCT(L))
		qs = qe
		me = -qe
		ms = me
		rho_back = -qe*this%n/PRODUCT(L)

		!spatial distribution initialize
		do i3 = 1,Nd(3)
			do i2 = 1,Nd(2)
				do i1 = 1,Nd(1)
					xp0(i1+Nd(1)*(i2-1)+Nd(1)*Nd(2)*(i3-1),:) =	&
								(/ (i1-1)*L(1)/Nd(1),(i2-1)*L(2)/Nd(2),(i3-1)*L(3)/Nd(3) /)
				end do
			end do
		end do
		xp0(:,1) = xp0(:,1) - this%A0*L(1)/this%n*SIN( 2.0_mp*pi*xp0(:,1)/L(1)*mode ) + 0.5_mp*L(1)

		!velocity distribution initialize
!		vp0 = vT*randn(N)
		vp0 = 0.0_mp
		pm = (/ ( i, i=1,this%n ) /)
		pm = 1 - 2*MOD(pm,2)
		vp0(:,1) = vp0(:,1) + pm*v0
	end subroutine

end module
