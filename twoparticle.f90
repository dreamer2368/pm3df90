module twoparticle

	use init
	use timeStep

	implicit none

contains

	subroutine twoParticleTest(v0, Ng)
		real(mp), intent(in) :: v0
		integer, intent(in) :: Ng(3)
		type(PM3D) :: this
		type(adjoint) :: adj
		real(mp) :: wp,eps0,Tf,Ti=20.0_mp,Tp,rho_back
		integer, parameter :: N = 2
		real(mp) :: xp0(N,3), vp0(N,3), qs(N), ms(N)
		real(mp) :: B0=0.0_mp, dB, dJdA, fdB(20)
		integer :: i

		wp = 1.0_mp
		Tp = 2.0_mp*pi/wp
		eps0 = 1.0_mp
		Tf = Ti + 10.0_mp*Tp

		call buildPM3D(this,Tf,Ti,Ng,N,B=B0)
		call buildAdjoint(adj,this)

		call twoParticleInit(this,v0,xp0,vp0,qs,ms,rho_back)
		call forwardsweep(this,xp0,vp0,qs,ms,rho_back)
		call printPlasma(this%r)

		call destroyAdjoint(adj)
		call destroyPM3D(this)
	end subroutine

	subroutine twoParticleInit(this,v0,xp0,vp0,qs,ms,rho_back)
		type(PM3D) :: this
		real(mp), intent(in) :: v0
		real(mp), intent(out), dimension(this%n,3) :: xp0, vp0
		real(mp), intent(out), dimension(this%n) :: qs, ms
		real(mp), intent(out) :: rho_back
		real(mp) :: qe, me, L(3)
		real(mp) :: Q										!Virial ratio = Potential/Kinetic
		L = this%L

		if( this%n .ne. 2 ) then
			print *, '!!FAULT :: Not two particles.'
			stop
		end if

		qe = -this%wp*this%wp/(this%n/PRODUCT(L))
		qs = qe
		me = -qe
		ms = me
		rho_back = -qe*this%n/PRODUCT(L)

		xp0(1,:) = 0.5_mp*this%L
		xp0(2,:) = 0.5_mp*this%L
		xp0(1,1) = xp0(1,1) - 0.25_mp*this%L(1)
		xp0(2,1) = xp0(2,1) + 0.25_mp*this%L(1)
		vp0(1,:) = v0*(/ 1.0_mp, 0.0_mp, 0.0_mp /)
		vp0(2,:) = -vp0(1,:)

		print *, '======Charges======'
		print *, 'particle 1 : ', qs(1)
		print *, 'particle 2 : ', qs(2)
		print *, '======Masses======'
		print *, 'particle 1 : ', ms(1)
		print *, 'particle 2 : ', ms(2)
		print *, '======Positions======'
		print *, 'particle 1 : ', xp0(1,:)
		print *, 'particle 2 : ', xp0(2,:)
		print *, '======Velocities======'
		print *, 'particle 1 : ', vp0(1,:)
		print *, 'particle 2 : ', vp0(2,:)
	end subroutine

end module
