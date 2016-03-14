module Landau

	use init
	use timeStep

	implicit none

contains

	subroutine Landau_AdjTest(vT,Ng,Nd,QoI,dQoI,source,Dsource,dQoI_dsource)
		real(mp), intent(in) :: vT
		integer, intent(in) :: Ng(3),Nd(3)
		type(PM3D) :: this
		type(adjoint) :: adj
		real(mp) :: wp,eps0,Ti,Tp,rho_back,L(3),dt
		real(mp) :: Tf(6)
		integer:: N
		real(mp) :: xp0(Nd(1)*Nd(2)*Nd(3),3), vp0(Nd(1)*Nd(2)*Nd(3),3), qs(Nd(1)*Nd(2)*Nd(3)), ms(Nd(1)*Nd(2)*Nd(3))
		real(mp) :: B0, dB, dJdA, fdB(20), ek(20)
		integer :: i,j
		interface
			subroutine QoI(adj,pm,i)
				use modPM3D
				use modAdj
				type(adjoint), intent(inout) :: adj
				type(PM3D), intent(in) :: pm
				integer, intent(in), optional :: i
			end subroutine
		end interface
		interface
			subroutine dQoI(adj,pm,k)
				use modPM3D
				use modAdj
				type(adjoint), intent(inout) :: adj
				type(PM3D), intent(in) :: pm
				integer, intent(in) :: k
			end subroutine
		end interface
		interface
			subroutine source(pm,k,str)
				use modPM3D
				type(PM3D), intent(inout) :: pm
				integer, intent(in) :: k
				character(len=*), intent(in) :: str
			end subroutine
		end interface
		interface
			subroutine Dsource(adj,pm,k,str)
				use modPM3D
				use modAdj
				type(adjoint), intent(inout) :: adj
				type(PM3D), intent(in) :: pm
				integer, intent(in) :: k
				character(len=*), intent(in) :: str
			end subroutine
		end interface
		interface
			subroutine dQoI_dsource(adj,pm,dJdA)
				use modPM3D
				use modAdj
				use constants
				type(adjoint), intent(in) :: adj
				type(PM3D), intent(in) :: pm
				real(mp), intent(out) :: dJdA
			end subroutine
		end interface
		N = PRODUCT(Nd)
		wp = 1.0_mp
		Tp = 2.0_mp*pi/wp
		eps0 = 1.0_mp
		dt = 0.1_mp
		Ti = dt
		L = 4.0_mp*pi

!		Tf = Ti + Tp*(/ 0.5_mp/pi, 1.0_mp, 2.0_mp, 4.0_mp, 10.0_mp, 20.0_mp /)
		Tf = Ti + Tp*(/ 0.25_mp/pi, 0.5_mp/pi, 1.0_mp, 4.0_mp, 10.0_mp, 20.0_mp /)
		B0 = 1.0_mp
		fdB = (/ ( EXP(-i*1.0_mp), i=1,20 ) /)
		open(unit=301,file='data/dA.bin',status='replace',form='unformatted',access='stream')
		write(301) B0*fdB
		close(301)

		open(unit=401,file='data/ek.bin',status='replace',form='unformatted',access='stream')
		open(unit=402,file='data/dJdA.bin',status='replace',form='unformatted',access='stream')
		do i=1,2!size(Tf)
			call buildPM3D(this,Tf(i),Ti,Ng,N,dt=dt,L=L,A=0.1_mp,B=B0)
			call buildAdjoint(adj,this)

			this%B0 = B0

			call LandauInit(this,Nd,0.0_mp,vT,xp0,vp0,qs,ms,rho_back)
			call forwardsweep(this,xp0,vp0,qs,ms,rho_back,source)
			call QoI(adj,this,0)

			call backward_sweep(adj,this,dQoI,Dsource)
			call dQoI_dsource(adj,this,dJdA)
			print *, dJdA
			write(402) dJdA

			do j=1,size(fdB)
				dB = B0*fdB(j)
				this%B0 = B0 + dB
				print *, 'B = ',this%B0
				call forwardsweep(this,xp0,vp0,qs,ms,rho_back,source)
				call QoI(adj,this,1)
				print *, (adj%J1 - adj%J0)/dB

				ek(j) = abs( (adj%J1 - adj%J0)/dB - dJdA )/abs(dJdA)
				print *, 'error = ', ek(j)
			end do

			write(401) ek

			call destroyAdjoint(adj)
			call destroyPM3D(this)
		end do
		close(401)
		close(402)
	end subroutine

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
