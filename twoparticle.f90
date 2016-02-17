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
		real(mp) :: wp,eps0,Ti=30.0_mp,Tp,rho_back
		real(mp) :: Tf(6)
		integer, parameter :: N = 2
		real(mp) :: xp0(N,3), vp0(N,3), qs(N), ms(N)
		real(mp) :: B0=1.0_mp, dB, dJdA, fdB(20), ek(20)
		integer :: i,j

		wp = 1.0_mp
		Tp = 2.0_mp*pi/wp
		eps0 = 1.0_mp

		Tf = Ti + Tp*(/ 0.5_mp/pi, 1.0_mp, 2.0_mp, 4.0_mp, 10.0_mp, 20.0_mp /)
		fdB = (/ ( EXP(-i*1.0_mp), i=1,20 ) /)
		open(unit=301,file='data/dA.bin',status='replace',form='unformatted',access='stream')
		write(301) B0*fdB
		close(301)

!		Tf = Ti + Tp*20.0_mp
!		i = 6

		open(unit=401,file='data/ek.bin',status='replace',form='unformatted',access='stream')
		do i=1,size(Tf)
			call buildPM3D(this,Tf(i),Ti,Ng,N,B=B0)
			call buildAdjoint(adj,this)

			call twoParticleInit(this,v0,xp0,vp0,qs,ms,rho_back)
			call forwardsweep(this,xp0,vp0,qs,ms,rho_back)
			call printPlasma(this%r)
			call QoI(adj,this,0)

			call backward_sweep(adj,this,dJdvp,dJdA)
			print *, dJdA

			do j=1,size(fdB)
				dB = B0*fdB(j)
				this%B0 = B0 + dB
				print *, 'B = ',this%B0
				call twoParticleInit(this,v0,xp0,vp0,qs,ms,rho_back)
				call forwardsweep(this,xp0,vp0,qs,ms,rho_back)
				call QoI(adj,this,1)
				print *, (adj%J1 - adj%J0)/dB

				ek(j) = abs( (adj%J1 - adj%J0)/dB - dJdA )
				print *, 'error = ', ek(j)
			end do

			write(401) ek

			call destroyAdjoint(adj)
			call destroyPM3D(this)
		end do
		close(401)
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
!		qs(1) = qe
!		qs(2) = -qe
		me = -qe
		ms = me
		rho_back = -qe*this%n/PRODUCT(L)
!		rho_back = 0.0_mp

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
