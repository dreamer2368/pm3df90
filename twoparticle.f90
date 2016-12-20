module twoparticle

	use init
	use timeStep

	implicit none

contains

!	subroutine twoParticle_rep_AdjTest(v0, Ng,QoI,dQoI,source,Dsource,dQoI_dsource)
!		real(mp), intent(in) :: v0
!		integer, intent(in) :: Ng(3)
!		type(PM3D) :: this
!		type(adjoint) :: adj
!		real(mp) :: wp,eps0,Ti=30.0_mp,Tp,rho_back,dt
!		real(mp) :: Tf(6)
!		integer, parameter :: N = 2
!		real(mp) :: xp0(N,3), vp0(N,3), qs(N), ms(N)
!		real(mp) :: B0=1.0_mp, dB, dJdA, fdB(20), ek(20)
!		integer :: i,j
!		interface
!			subroutine QoI(adj,pm,i)
!				use modPM3D
!				use modAdj
!				type(adjoint), intent(inout) :: adj
!				type(PM3D), intent(in) :: pm
!				integer, intent(in), optional :: i
!			end subroutine
!		end interface
!		interface
!			subroutine dQoI(adj,pm,k)
!				use modPM3D
!				use modAdj
!				type(adjoint), intent(inout) :: adj
!				type(PM3D), intent(in) :: pm
!				integer, intent(in) :: k
!			end subroutine
!		end interface
!		interface
!			subroutine source(pm,k,str)
!				use modPM3D
!				type(PM3D), intent(inout) :: pm
!				integer, intent(in) :: k
!				character(len=*), intent(in) :: str
!			end subroutine
!		end interface
!		interface
!			subroutine Dsource(adj,pm,k,str)
!				use modPM3D
!				use modAdj
!				type(adjoint), intent(inout) :: adj
!				type(PM3D), intent(in) :: pm
!				integer, intent(in) :: k
!				character(len=*), intent(in) :: str
!			end subroutine
!		end interface
!		interface
!			subroutine dQoI_dsource(adj,pm,dJdA)
!				use modPM3D
!				use modAdj
!				use constants
!				type(adjoint), intent(in) :: adj
!				type(PM3D), intent(in) :: pm
!				real(mp), intent(out) :: dJdA
!			end subroutine
!		end interface
!		wp = 1.0_mp
!		Tp = 2.0_mp*pi/wp
!		eps0 = 1.0_mp
!		dt = 0.2_mp
!
!		Tf = Ti + Tp*(/ 0.5_mp/pi, 1.0_mp, 2.0_mp, 4.0_mp, 10.0_mp, 20.0_mp /)
!		fdB = (/ ( EXP(-i*1.0_mp), i=1,20 ) /)
!		open(unit=301,file='data/dA.bin',status='replace',form='unformatted',access='stream')
!		write(301) B0*fdB
!		close(301)
!
!		open(unit=401,file='data/ek.bin',status='replace',form='unformatted',access='stream')
!		open(unit=402,file='data/dJdA.bin',status='replace',form='unformatted',access='stream')
!		do i=1,size(Tf)
!			call buildPM3D(this,Tf(i),Ti,Ng,N,dt=dt,B=B0)
!			call buildAdjoint(adj,this)
!
!			B0 = 0.1_mp*this%L(1)
!			this%B0 = B0
!
!			call twoParticleInit(this,v0,xp0,vp0,qs,ms,rho_back)
!			call forwardsweep(this,xp0,vp0,qs,ms,rho_back,source)
!!			call printPlasma(this%r)
!			call QoI(adj,this,0)
!
!			call backward_sweep(adj,this,dQoI,Dsource)
!			call dQoI_dsource(adj,this,dJdA)
!			print *, dJdA
!			write(402) dJdA
!
!			do j=1,size(fdB)
!				dB = B0*fdB(j)
!				this%B0 = B0 + dB
!				print *, 'B = ',this%B0
!				call twoParticleInit(this,v0,xp0,vp0,qs,ms,rho_back)
!				call forwardsweep(this,xp0,vp0,qs,ms,rho_back,source)
!				call QoI(adj,this,1)
!				print *, (adj%J1 - adj%J0)/dB
!
!				ek(j) = abs( (adj%J1 - adj%J0)/dB - dJdA )/abs(dJdA)
!				print *, 'error = ', ek(j)
!			end do
!
!			write(401) ek
!
!			call destroyAdjoint(adj)
!			call destroyPM3D(this)
!		end do
!		close(401)
!		close(402)
!	end subroutine
!
!	subroutine twoParticle_att_AdjTest(v0, Ng,QoI,dQoI,source,Dsource,dQoI_dsource)
!		real(mp), intent(in) :: v0
!		integer, intent(in) :: Ng(3)
!		type(PM3D) :: this
!		type(adjoint) :: adj
!		real(mp) :: wp,eps0,Ti=0.0_mp,Tp,rho_back,dt
!		real(mp) :: Tf(6)
!		integer, parameter :: N = 2
!		real(mp) :: xp0(N,3), vp0(N,3), qs(N), ms(N)
!		real(mp) :: A0, B0=0.0_mp, dJdA, fdB(26), ek(26)
!		integer :: i,j
!		interface
!			subroutine QoI(adj,pm,i)
!				use modPM3D
!				use modAdj
!				type(adjoint), intent(inout) :: adj
!				type(PM3D), intent(in) :: pm
!				integer, intent(in), optional :: i
!			end subroutine
!		end interface
!		interface
!			subroutine dQoI(adj,pm,k)
!				use modPM3D
!				use modAdj
!				type(adjoint), intent(inout) :: adj
!				type(PM3D), intent(in) :: pm
!				integer, intent(in) :: k
!			end subroutine
!		end interface
!		interface
!			subroutine source(pm,k,str)
!				use modPM3D
!				type(PM3D), intent(inout) :: pm
!				integer, intent(in) :: k
!				character(len=*), intent(in) :: str
!			end subroutine
!		end interface
!		interface
!			subroutine Dsource(adj,pm,k,str)
!				use modPM3D
!				use modAdj
!				type(adjoint), intent(inout) :: adj
!				type(PM3D), intent(in) :: pm
!				integer, intent(in) :: k
!				character(len=*), intent(in) :: str
!			end subroutine
!		end interface
!		interface
!			subroutine dQoI_dsource(adj,pm,dJdA)
!				use modPM3D
!				use modAdj
!				use constants
!				type(adjoint), intent(in) :: adj
!				type(PM3D), intent(in) :: pm
!				real(mp), intent(out) :: dJdA
!			end subroutine
!		end interface
!		wp = 1.0_mp
!		Tp = 2.0_mp*pi/wp
!		eps0 = 1.0_mp
!		dt = 0.01_mp
!
!		Tf = Ti + Tp*(/ 0.5_mp/pi, 1.0_mp, 2.0_mp, 4.0_mp, 10.0_mp, 20.0_mp /)
!		fdB = (/ ( EXP(-i*1.0_mp), i=1,26 ) /)
!		open(unit=301,file='data/dA.bin',status='replace',form='unformatted',access='stream')
!		write(301) fdB
!		close(301)
!
!		open(unit=401,file='data/ek.bin',status='replace',form='unformatted',access='stream')
!		open(unit=402,file='data/dJdA.bin',status='replace',form='unformatted',access='stream')
!		do i=1,size(Tf)
!			call buildPM3D(this,Tf(i),Ti,Ng,N,dt=dt,B=B0)
!			call buildAdjoint(adj,this)
!
!!			A0 = 0.1_mp*this%L(1)
!			A0 = 0.9_mp*this%m%dx(1)
!			this%A0 = A0
!
!			call twoParticleInit(this,v0,xp0,vp0,qs,ms,rho_back)
!			call forwardsweep(this,xp0,vp0,qs,ms,rho_back,source)
!!			call printPlasma(this%r)
!			call QoI(adj,this,0)
!
!			call backward_sweep(adj,this,dQoI,Dsource)
!			call dQoI_dsource(adj,this,dJdA)
!			print *, dJdA
!			write(402) dJdA
!
!			do j=1,size(fdB)
!				this%B0 = B0 + fdB(j)
!				print *, 'B = ',this%B0
!				call twoParticleInit(this,v0,xp0,vp0,qs,ms,rho_back)
!				call forwardsweep(this,xp0,vp0,qs,ms,rho_back,source)
!				call QoI(adj,this,1)
!				print *, (adj%J1 - adj%J0)/fdB(j)
!
!				ek(j) = abs( (adj%J1 - adj%J0)/fdB(j) - dJdA )/abs(dJdA)
!				print *, 'error = ', ek(j)
!			end do
!
!			write(401) ek
!
!			call destroyAdjoint(adj)
!			call destroyPM3D(this)
!		end do
!		close(401)
!		close(402)
!	end subroutine
!
!	subroutine twoParticleTraj(v0, Ng)
!		real(mp), intent(in) :: v0
!		integer, intent(in) :: Ng(3)
!		type(PM3D) :: this
!		type(adjoint) :: adj
!		real(mp) :: wp,eps0,Ti=0.0_mp,Tf,Tp,rho_back
!		integer, parameter :: N = 2
!		real(mp) :: xp0(N,3), vp0(N,3), qs(N), ms(N)
!		real(mp) :: B0=0.0_mp, dB
!		integer :: i,j
!
!		wp = 1.0_mp
!		Tp = 2.0_mp*pi/wp
!		eps0 = 1.0_mp
!
!		Tf = Ti + Tp*10.0_mp
!
!		call buildPM3D(this,Tf,Ti,Ng,N,dt=0.01_mp,B=B0)
!		call buildAdjoint(adj,this)
!
!		this%A0 = 0.9_mp*this%m%dx(1)
!
!		call twoParticleInit(this,v0,xp0,vp0,qs,ms,rho_back)
!		print *, '======Charges======'
!		print *, 'particle 1 : ', qs(1)
!		print *, 'particle 2 : ', qs(2)
!		print *, '======Masses======'
!		print *, 'particle 1 : ', ms(1)
!		print *, 'particle 2 : ', ms(2)
!		print *, '======Positions======'
!		print *, 'particle 1 : ', xp0(1,:)
!		print *, 'particle 2 : ', xp0(2,:)
!		print *, '======Velocities======'
!		print *, 'particle 1 : ', vp0(1,:)
!		print *, 'particle 2 : ', vp0(2,:)
!
!		call forwardsweep(this,xp0,vp0,qs,ms,rho_back,Orbit_radius)
!		call printPlasma(this%r)
!
!		call destroyAdjoint(adj)
!		call destroyPM3D(this)
!	end subroutine
!
!	subroutine twoParticleInit(this,v0,xp0,vp0,qs,ms,rho_back)
!		type(PM3D) :: this
!		real(mp), intent(in) :: v0
!		real(mp), intent(out), dimension(this%n,3) :: xp0, vp0
!		real(mp), intent(out), dimension(this%n) :: qs, ms
!		real(mp), intent(out) :: rho_back
!		real(mp) :: qe, me, L(3), r0
!		real(mp) :: Q										!Virial ratio = Potential/Kinetic
!		L = this%L
!
!		if( this%n .ne. 2 ) then
!			print *, '!!FAULT :: Not two particles.'
!			stop
!		end if
!
!		qe = -this%wp*this%wp/(this%n/PRODUCT(L))
!		me = -qe
!
!		!Repulsive particles
!		qs = qe
!		ms = me
!		rho_back = -qe*this%n/PRODUCT(L)
!
!!		xp0(1,:) = this%m%dx*( 0.5_mp + 32.0_mp )
!!		xp0(2,:) = this%m%dx*( 0.5_mp + 32.0_mp )
!		xp0(1,:) = 0.5_mp*this%L
!		xp0(2,:) = 0.5_mp*this%L
!		xp0(1,1) = xp0(1,1) - 0.25_mp*this%L(1)
!		xp0(2,1) = xp0(2,1) + 0.25_mp*this%L(1)
!		vp0(1,:) = v0*(/ 1.0_mp, 0.0_mp, 0.0_mp /)
!		vp0(2,:) = -vp0(1,:)
!
!		!Attractive particles
!!		qs(1) = qe
!!		qs(2) = -qe
!!		rho_back = 0.0_mp
!!		ms = me
!!
!!		r0 = this%A0
!!
!!		xp0(1,:) = this%m%dx*( 0.5_mp + 16.0_mp )
!!		xp0(2,:) = this%m%dx*( 0.5_mp + 16.0_mp )
!!		xp0(1,1) = xp0(1,1) - 0.5_mp*r0
!!		xp0(2,1) = xp0(2,1) + 0.5_mp*r0
!!		vp0(1,:) = v0*(/ 0.0_mp, 1.0_mp, 0.0_mp /)
!!		vp0(2,:) = -vp0(1,:)
!	end subroutine
!
!	subroutine EfieldKernel()
!		type(PM3D) :: pm
!		real(mp) :: Ti=0.0_mp, Tf=10.0_mp, qe, me, qs(2), ms(2), rho_back, L(3), rhs(32,32,32)
!		real(mp) :: xp0(2,3), vp0(2,3), xd(1000)
!		integer :: Ng(3), N=2, Nd=1000
!		integer :: i
!
!		Ng = (/32,32,32/)
!		call buildPM3D(pm,Tf,Ti,Ng,N)
!
!		L = pm%L
!		qe = -pm%wp*pm%wp/(pm%n/PRODUCT(L))
!		me = -qe
!		xd = (/ ( L(1)*i*1.0_mp/Nd, i=1,1000 ) /)
!		vp0 = 0.0_mp
!		xp0(1,:) = 0.5_mp*L
!		xp0(2,:) = 0.5_mp*L
!		!Attractive particles
!!		qs(1) = qe
!!		qs(2) = -qe
!!		rho_back = 0.0_mp
!!		ms = me
!		!Repulsive particles
!		qs = qe
!		rho_back = -qe*pm%n/PRODUCT(L)
!		ms = me
!
!		open(unit=301,file='data/Fpx.bin',status='replace',form='unformatted',access='stream')
!		do i=1,Nd
!			xp0(2,1) = xd(i)
!			call setPlasma(pm%p,xp0,vp0,qs,ms)
!			call setMesh(pm%m, rho_back)
!
!			call assignMatrix(pm%a,pm%p,pm%m,pm%p%xp)
!			call chargeAssign(pm%a,pm%p,pm%m)
!			rhs = -pm%m%rho/pm%eps0
!			call FFTPoisson(pm%m%phi,rhs,pm%m%W)
!			pm%m%E = - Gradient(pm%m%phi,pm%m%dx,pm%ng)
!			call forceAssign(pm%a,pm%p,pm%m)
!			write(301) pm%p%qs(2)/pm%p%ms(2)*pm%p%Ep(2,1)
!		end do
!		close(301)
!
!		open(unit=302,file='data/xd.bin',status='replace',form='unformatted',access='stream')
!		write(302) xd
!		close(302)
!
!		call destroyPM3D(pm)
!	end subroutine

end module
