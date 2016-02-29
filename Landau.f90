module Landau

	use init
	use timeStep

	implicit none

contains

	subroutine LandauTraj(Ng, Nd)
		integer, intent(in) :: Ng(3), Nd(3)
		type(PM3D) :: this
		type(adjoint) :: adj
		real(mp) :: wp,eps0,Ti=0.0_mp,Tf,Tp,rho_back
		integer :: N
		real(mp) :: xp0(Nd(1)*Nd(2)*Nd(3),3), vp0(Nd(1)*Nd(2)*Nd(3),3), qs(Nd(1)*Nd(2)*Nd(3)), ms(Nd(1)*Nd(2)*Nd(3))
		real(mp) :: B0=0.0_mp, dB
		integer :: i,j
		N = Nd(1)*Nd(2)*Nd(3)

		wp = 1.0_mp
		Tp = 2.0_mp*pi/wp
		eps0 = 1.0_mp

		Tf = Ti + Tp*10.0_mp

		call buildPM3D(this,Tf,Ti,Ng,N,dt=0.1_mp,L=4.0_mp*pi*(/1,1,1/),A=0.1_mp)
		call buildAdjoint(adj,this)

		call LandauInit(this,Nd,0.0_mp,1.0_mp,xp0,vp0,qs,ms,rho_back)

		call forwardsweep(this,xp0,vp0,qs,ms,rho_back,Orbit_radius)
		call printPlasma(this%r)

		call destroyAdjoint(adj)
		call destroyPM3D(this)
	end subroutine

	subroutine LandauInit(this,Nd,v0,vT,xp0,vp0,qs,ms,rho_back)			!generate initial distribution
		type(PM3D), intent(inout) :: this
		integer, intent(in) :: Nd(3)
		real(mp), intent(in) :: v0, vT
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
		xp0(:,1) = xp0(:,1) - this%A0*SIN( 2.0_mp*pi*xp0(:,1)/L(1)*mode ) + 0.5_mp*L(1)

		!velocity distribution initialize
		vp0(:,1) = vT*randn(this%n)
		vp0(:,2) = vT*randn(this%n)
		vp0(:,3) = vT*randn(this%n)
		pm = (/ ( i, i=1,this%n ) /)
		pm = 1 - 2*MOD(pm,2)
		vp0(:,1) = vp0(:,1) + pm*v0
	end subroutine

end module
