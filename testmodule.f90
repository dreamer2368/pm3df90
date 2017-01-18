module testmodule

	use init
	use timeStep

	implicit none

contains

!	subroutine Compare_Traj_twostream(v0,Ng,Nd,source,B1,B2)
!		real(mp), intent(in) :: v0, B1, B2
!		integer, intent(in) :: Ng(3), Nd(3)
!		type(PM3D) :: this
!		type(adjoint) :: adj
!		real(mp) :: Tf=1.6_mp,Ti=0.2_mp,rho_back
!		integer :: N
!		real(mp) :: xp0(PRODUCT(Nd),3), vp0(PRODUCT(Nd),3), qs(PRODUCT(Nd)), ms(PRODUCT(Nd))
!		interface
!			subroutine source(pm,k,str)
!				use modPM3D
!				type(PM3D), intent(inout) :: pm
!				integer, intent(in) :: k
!				character(len=*), intent(in) :: str
!			end subroutine
!		end interface
!		N = PRODUCT(Nd)
!
!		call buildPM3D(this,Tf,Ti,Ng,N,dt=0.2_mp,B=0.0_mp)
!		call buildAdjoint(adj,this)
!
!		call particle_initialize(this,Nd,v0,xp0,vp0,qs,ms,rho_back)
!		call forwardsweep(this,xp0,vp0,qs,ms,rho_back,source)
!		call printPlasma(this%r,'0')
!
!		this%B0 = 0.0_mp+B1
!		call particle_initialize(this,Nd,v0,xp0,vp0,qs,ms,rho_back)
!		call forwardsweep(this,xp0,vp0,qs,ms,rho_back,source)
!		call printPlasma(this%r,'1')
!
!		this%B0 = 0.0_mp+B2
!		call particle_initialize(this,Nd,v0,xp0,vp0,qs,ms,rho_back)
!		call forwardsweep(this,xp0,vp0,qs,ms,rho_back,source)
!		call printPlasma(this%r,'2')
!
!		call destroyAdjoint(adj)
!		call destroyPM3D(this)
!	end subroutine
!
!	subroutine Twostream_AdjTest(v0,Ng,Nd,QoI,dQoI,source,Dsource,dQoI_dsource)
!		real(mp), intent(in) :: v0
!		integer, intent(in) :: Ng(3),Nd(3)
!		type(PM3D) :: this
!		type(adjoint) :: adj
!		real(mp) :: wp,eps0,Ti,Tp,rho_back,L(3),dt
!		real(mp) :: Tf(6)
!		integer:: N
!		real(mp) :: xp0(Nd(1)*Nd(2)*Nd(3),3), vp0(Nd(1)*Nd(2)*Nd(3),3), qs(Nd(1)*Nd(2)*Nd(3)), ms(Nd(1)*Nd(2)*Nd(3))
!??????		real(mp) :: B0, dB, dJdA, fdB(20), ek(20)
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
!		N = PRODUCT(Nd)
!		wp = 1.0_mp
!		Tp = 2.0_mp*pi/wp
!		eps0 = 1.0_mp
!		dt = 0.1_mp
!		Ti = dt
!
!		Tf = Ti + Tp*(/ 0.5_mp/pi, 1.0_mp, 2.0_mp, 4.0_mp, 10.0_mp, 20.0_mp /)
!		B0 = 0.0_mp
!		fdB = (/ ( EXP(-i*1.0_mp), i=1,20 ) /)
!		open(unit=301,file='data/dA.bin',status='replace',form='unformatted',access='stream')
!		write(301) fdB
!		close(301)
!
!		open(unit=401,file='data/ek.bin',status='replace',form='unformatted',access='stream')
!		open(unit=402,file='data/dJdA.bin',status='replace',form='unformatted',access='stream')
!		do i=1,size(Tf)
!			call buildPM3D(this,Tf(i),Ti,Ng,N,dt=dt,A=1.0_mp,B=B0)
!			call buildAdjoint(adj,this)
!
!			this%B0 = B0
!
!			call particle_initialize(this,Nd,v0,xp0,vp0,qs,ms,rho_back)
!			call forwardsweep(this,xp0,vp0,qs,ms,rho_back,source)
!			call QoI(adj,this,0)
!
!			call backward_sweep(adj,this,dQoI,Dsource)
!			call dQoI_dsource(adj,this,dJdA)
!			print *, dJdA
!			write(402) dJdA
!
!			do j=1,size(fdB)
!				dB = fdB(j)
!				this%B0 = B0 + dB
!				print *, 'B = ',this%B0
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

	subroutine test_fullAdjoint(v0, Ng, Nd,QoI,dQoI,control,Dcontrol,dQoI_dcontrol)
		real(mp), intent(in) :: v0
		integer, intent(in) :: Ng(3), Nd(3)
		type(PM3D) :: this
		type(adjoint) :: adj
		real(mp) :: Tf=0.8_mp,Ti=0.4_mp,rho_back
		integer :: N
		real(mp) :: J0, J1
		real(mp) :: A0(1), grad(1), dA , fdA(20)
		integer :: i
		interface
			subroutine QoI(pm,k,J)
				use modPM3D
				type(PM3D), intent(in) :: pm
				integer, intent(in) :: k
				real(mp), intent(inout) :: J
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
			subroutine control(pm,k,str)
				use modPM3D
				type(PM3D), intent(inout) :: pm
				integer, intent(in) :: k
				character(len=*), intent(in) :: str
			end subroutine
		end interface
		interface
			subroutine Dcontrol(adj,pm,k,str)
				use modPM3D
				use modAdj
				type(adjoint), intent(inout) :: adj
				type(PM3D), intent(in) :: pm
				integer, intent(in) :: k
				character(len=*), intent(in) :: str
			end subroutine
		end interface
		interface
			subroutine dQoI_dcontrol(adj,pm,k,str,grad)
				use modPM3D
				use modAdj
				use constants
				type(adjoint), intent(in) :: adj
				type(PM3D), intent(in) :: pm
				integer, intent(in) :: k
				character(len=*), intent(in) :: str
				real(mp), intent(inout) :: grad(:)
			end subroutine
		end interface

		N = PRODUCT(Nd)
		A0 = 0.0_mp
		call buildPM3D(this,Tf,Ti,Ng,1,A=A0,dir='test')

		call twostream_initialize(this,Nd,v0)
		call forwardsweep(this,control,QoI,J0)
		call printPlasma(this%r)
		print *, 'J0: ',J0

		call buildAdjoint(adj,this)
		call backward_sweep(adj,this,grad,dQoI,Dcontrol,dQoI_dcontrol,control)
		print *, 'grad: ', grad
		open(unit=301,file='data/dJdA.bin',status='replace',form='unformatted',access='stream')
		write(301) grad
		close(301)

		fdA = (/ ( 10.0_mp**(-i*1.0_mp), i=1,20 ) /)
		open(unit=301,file='data/dA.bin',status='replace',form='unformatted',access='stream')
		write(301) fdA
		close(301)

		open(unit=401,file='data/dJdAFD.bin',status='replace',form='unformatted',access='stream')
		do i=1,size(fdA)
			dA = fdA(i)
			this%A0 = A0 + dA
			print *, 'A = ',this%A0
			call twostream_initialize(this,Nd,v0)
			call forwardsweep(this,control,QoI,J1)
			print *, 'J1: ',J1
			print *, 'FD: ',(J1 - J0)/dA
			print *, 'error = ', ABS( (grad - (J1 - J0)/dA)/grad )
			write(401) (J1 - J0)/dA
		end do
		close(401)

		call destroyAdjoint(adj)
		call destroyPM3D(this)
	end subroutine

	subroutine test_checkpoint(N,Np,QoI)
		integer, intent(in) :: N(3)										!!grid number
		integer, intent(in) :: Np										!!number of particles
		type(PM3D) :: pm, pm1
		type(adjoint) :: adj

		real(mp) :: Tf = 1.0_mp, Ti = 1.0_mp
		real(mp) :: dx(3)												!!grid size
		real(mp) :: qe, qs, ms
		real(mp) :: L(3) = (/ 2.0_mp, 2.0_mp, 2.0_mp /)					!cubic size
		real(mp), dimension(N(1),N(2),N(3)) :: rho_back, rhs, phi
		real(mp), dimension(N(1),N(2),N(3),3) :: Es
		real(mp), dimension(Np,3) :: xp, vp, xp1, vp1, dxp, dxps1, dxps2

		real(mp) :: J0,J1,dJdxp(Np,3)
		real(mp) :: fxp(20)
		integer :: i,j,k(2)

		interface
			subroutine QoI(pm,k,J)
				use modPM3D
				type(PM3D), intent(in) :: pm
				integer, intent(in) :: k
                real(mp), intent(inout) :: J
			end subroutine
		end interface

		call buildPM3D(pm,Tf,Ti,N,1,L=L,dir='adj_test',mod_input=3)
		call buildPM3D(pm1,Tf,Ti,N,1,L=L,dir='adj_test',mod_input=3)

		!particle, mesh setup
!		print *, 'Original xp'
		do j=1,Np
			xp(j,:) = 0.1_mp*L*(/ j, 7*j, 4*j /)
!			print *, xp(j,:)
		end do
		vp = 0.1_mp
		qe = -0.1_mp
		qs = qe
		ms = -qs
        call buildSpecies(pm%p(1),Np,qs,ms,1.0_mp)
		call setSpecies(pm%p(1),Np,xp,vp)
		rho_back = -qe*pm%p(1)%np/PRODUCT(L)
		call setMesh(pm%m,rho_back)

        call buildSpecies(pm1%p(1),Np,qs,ms,1.0_mp)
		rho_back = -qe*pm1%p(1)%np/PRODUCT(L)
		call setMesh(pm1%m,rho_back)

        call forwardsweep(pm,Null_input,QoI,J0)
        call printPlasma(pm%r)
        print *, pm%p(1)%xp
        print *, pm%p(1)%vp
        print *, pm%p(1)%Ep

        call checkpoint(pm1,pm%r,pm%nt,Null_input)
        print *, pm%p(1)%xp
        print *, pm%p(1)%vp
        print *, pm%p(1)%Ep

		call destroyPM3D(pm)
		call destroyPM3D(pm1)
	end subroutine

	subroutine test_particle_adj3(N,Np,QoI,dQoI)
		integer, intent(in) :: N(3)										!!grid number
		integer, intent(in) :: Np										!!number of particles
		type(PM3D) :: pm
		type(adjoint) :: adj

		real(mp) :: Tf = 0.2_mp, Ti = 0.2_mp
		real(mp) :: dx(3)												!!grid size
		real(mp) :: qe, qs, ms
		real(mp) :: L(3) = (/ 2.0_mp, 2.0_mp, 2.0_mp /)					!cubic size
		real(mp), dimension(N(1),N(2),N(3)) :: rho_back, rhs, phi
		real(mp), dimension(N(1),N(2),N(3),3) :: Es
		real(mp), dimension(Np,3) :: xp, vp, xp1, vp1, dxp, dxps1, dxps2

		real(mp) :: A(3), J0,J1,grad(12),dJdxp(Np,3)
		real(mp) :: fxp(20)
		integer :: i,j,k(2)

		interface
			subroutine QoI(pm,k,J)
				use modPM3D
				type(PM3D), intent(in) :: pm
				integer, intent(in) :: k
                real(mp), intent(inout) :: J
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

    A=(/1.0_mp,1.0_mp,0.0_mp/)

		call buildPM3D(pm,Tf,Ti,N,1,L=L,A=A,dir='adj_test')

		!particle, mesh setup
!		print *, 'Original xp'
		do j=1,Np
			xp(j,:) = 0.1_mp*L*(/ j, 7*j, 4*j /)
!			print *, xp(j,:)
		end do
		vp = 0.1_mp
		qe = -0.1_mp
		qs = qe
		ms = -qs
        call buildSpecies(pm%p(1),Np,qs,ms,1.0_mp)
		call setSpecies(pm%p(1),Np,xp,vp)
		rho_back = -qe*pm%p(1)%np/PRODUCT(L)
		call setMesh(pm%m,rho_back)

		call forwardsweep(pm,Null_input,QoI,J0)

		call buildAdjoint(adj,pm)
    call backward_sweep(adj,pm,grad,dQoI,Null_Dinput,testdJdA,testControl)

		print *, 'dJdxp'
    print *, grad(1:3)
    print *, grad(4:6)
		print *, 'dJdvp'
    print *, grad(7:9)
    print *, grad(10:12)

		!FD approximation - choose the component that you want to measure
		k = (/2,2/)
		open(unit=301,file='data/dJdA.bin',status='replace',form='unformatted',access='stream')
		write(301) -adj%p(1)%xp(k(1),k(2))/pm%dt
		close(301)
		fxp = (/ ( 10.0_mp**(-i*1.0_mp), i=1,20 ) /)
		open(unit=301,file='data/dA.bin',status='replace',form='unformatted',access='stream')
		write(301) ABS( xp(k(1),k(2))*fxp )
		close(301)

		open(unit=301,file='data/dJdAFD.bin',status='replace',form='unformatted',access='stream')
		do i=1,20
        A=(/1.0_mp*k(1),1.0_mp*k(2),fxp(i)/)
        deallocate(pm%A0)
        allocate(pm%A0(3))
        pm%A0 = A
        print *, pm%A0
        call setSpecies(pm%p(1),Np,xp,vp)

        call forwardsweep(pm,testControl,QoI,J1)

        dxp = 0.0_mp
        dxp(k(1),k(2)) = xp(k(1),k(2))*fxp(i)

			print *, 'dJdxp(',k(1),',',k(2),')=', (J1-J0)/dxp(k(1),k(2))
			print *, 'error = ', ABS( ( -adj%p(1)%xp(k(1),k(2))/pm%dt - (J1-J0)/dxp(k(1),k(2)) ) )
			write(301) (J1-J0)/dxp(k(1),k(2))
		end do
		close(301)

		call destroyPM3D(pm)
		call destroyAdjoint(adj)
	end subroutine

	subroutine test_particle_adj2(N,Np,QoI)
		integer, intent(in) :: N(3)										!!grid number
		integer, intent(in) :: Np										!!number of particles
		type(PM3D) :: pm
		type(adjoint) :: adj

		real(mp) :: Tf = 1.0_mp, Ti = 1.0_mp
		real(mp) :: dx(3)												!!grid size
		real(mp) :: qe, qs, ms
		real(mp) :: L(3) = (/ 2.0_mp, 2.0_mp, 2.0_mp /)					!cubic size
		real(mp), dimension(N(1),N(2),N(3)) :: rho_back, rhs, phi
		real(mp), dimension(N(1),N(2),N(3),3) :: Es
		real(mp), dimension(Np,3) :: xp, vp, xp1, vp1, dxp, dxps1, dxps2

		real(mp) :: J0,J1,dJdxp(Np,3)
		real(mp) :: fxp(20)
		integer :: i,j,k(2)

		interface
			subroutine QoI(pm,k,J)
				use modPM3D
				type(PM3D), intent(in) :: pm
				integer, intent(in) :: k
                real(mp), intent(inout) :: J
			end subroutine
		end interface

		call buildPM3D(pm,Tf,Ti,N,1,L=L,dir='adj_test')

		!particle, mesh setup
!		print *, 'Original xp'
		do j=1,Np
			xp(j,:) = 0.1_mp*L*(/ j, 7*j, 4*j /)
!			print *, xp(j,:)
		end do
		vp = 0.1_mp
		qe = -0.1_mp
		qs = qe
		ms = -qs
        call buildSpecies(pm%p(1),Np,qs,ms,1.0_mp)
		call setSpecies(pm%p(1),Np,xp,vp)
		rho_back = -qe*pm%p(1)%np/PRODUCT(L)
		call setMesh(pm%m,rho_back)

		!particle move
        call move(pm%p(1),pm%dt)
		!assignment
        call assignMatrix(pm%a(1),pm%p(1),pm%m,pm%p(1)%xp)
		call chargeAssign(pm%a,pm%p,pm%m)
		!field solver
		rhs = -pm%m%rho/pm%eps0
		call FFTPoisson(pm%m%phi,rhs,pm%m%W)
		pm%m%E = - Gradient(pm%m%phi,pm%m%dx,pm%m%ng)
		!force assignment
		call forceAssign(pm%a(1),pm%p(1),pm%m)
		!particle accel
        call accel(pm%p(1),pm%dt)
		!QoI evaluation
		call QoI(pm,pm%nt,J0)

		call buildAdjoint(adj,pm)
        call reset_Dadj(adj)

		!Adjoint sensitivity solver
		adj%p(1)%vp = - 2.0_mp*pm%p(1)%vp*pm%dt
        call Adj_accel(adj)

		adj%p(1)%Ep = pm%p(1)%qs/pm%p(1)%ms*adj%p(1)%vp

        adj%m%E = 0.0_mp
		call Adj_forceAssign_E(pm%a(1),adj%p(1)%Ep,adj%m%E)
        adj%m%E = adj%m%E + adj%dm%E

		call FFTAdj(adj%m%E,adj%m%rho,pm%m%W,pm%m%dx)
		adj%m%rho = -adj%m%rho/pm%eps0

		call Adj_chargeAssign(pm%a(1),pm%p(1),pm%m,adj%m%rho,adj%dp(1)%xp)
		call Adj_forceAssign_xp(pm%a(1),pm%m,pm%m%E,adj%p(1)%Ep,adj%dp(1)%xp)
      call Adj_move(adj)
!		adj%xps = - pm%dt*( dxps1 + dxps2 )

		print *, 'dJdxp'
		do i=1,Np
			print *, -adj%p(1)%xp(i,:)/pm%dt
		end do
		print *, 'dJdvp'
		do i=1,Np
			print *, -adj%p(1)%vp(i,:)/pm%dt - adj%p(1)%xp(i,:)
		end do

		!FD approximation - choose the component that you want to measure
		k = (/1,3/)
		open(unit=301,file='data/dJdA.bin',status='replace',form='unformatted',access='stream')
		write(301) -adj%p(1)%xp(k(1),k(2))/pm%dt
		close(301)
		fxp = (/ ( 10.0_mp**(-i*1.0_mp), i=1,20 ) /)
		open(unit=301,file='data/dA.bin',status='replace',form='unformatted',access='stream')
		write(301) ABS( xp(k(1),k(2))*fxp )
		close(301)

		open(unit=301,file='data/dJdAFD.bin',status='replace',form='unformatted',access='stream')
		do i=1,20
			dxp = 0.0_mp
			dxp(k(1),k(2)) = xp(k(1),k(2))*fxp(i)
			xp1 = xp + dxp
			call setSpecies(pm%p(1),Np,xp1,vp)

			!particle move
            call move(pm%p(1),pm%dt)
			!assignment
			call assignMatrix(pm%a(1),pm%p(1),pm%m,pm%p(1)%xp)
			call chargeAssign(pm%a,pm%p,pm%m)
			!field solver
			rhs = -pm%m%rho/pm%eps0
			call FFTPoisson(pm%m%phi,rhs,pm%m%W)
			pm%m%E = - Gradient(pm%m%phi,pm%m%dx,pm%m%ng)
			!force assignment
			call forceAssign(pm%a(1),pm%p(1),pm%m)
			!particle accel
            call accel(pm%p(1),pm%dt)
			!QoI evaluation
            call QoI(pm,pm%nt,J1)

			print *, 'dJdxp(',k(1),',',k(2),')=', (J1-J0)/dxp(k(1),k(2))
			print *, 'error = ', ABS( ( -adj%p(1)%xp(k(1),k(2))/pm%dt - (J1-J0)/dxp(k(1),k(2)) ) )
			write(301) (J1-J0)/dxp(k(1),k(2))
		end do
		close(301)

		call destroyPM3D(pm)
		call destroyAdjoint(adj)
	end subroutine

	subroutine twostream(v0,Ng,Nd)
		real(mp), intent(in) :: v0
		integer, intent(in) :: Ng(3), Nd(3)
		type(PM3D) :: this
		real(mp) :: Tf=100.0_mp,Ti=5.0_mp,rho_back
		integer :: N
		real(mp) :: xp0(PRODUCT(Nd),3), vp0(PRODUCT(Nd),3), qs, ms

		N = PRODUCT(Nd)
		call buildPM3D(this,Tf,Ti,Ng,1,dir='test')

		call twostream_initialize(this,Nd,v0)

		call forwardsweep(this,Null_input)
		call printPlasma(this%r)

		call destroyPM3D(this)
	end subroutine
!
!	subroutine test_FFTPoisson_adj(N,Nf)
!		integer, intent(in) :: N(3)										!!grid number
!		integer, intent(in) :: Nf										!!number of fourier modes
!
!		real(mp) :: dx(3)												!!grid size
!		real(mp) :: xg(N(1)), yg(N(2)), zg(N(3))
!		real(mp) :: eps0 = 1.0_mp
!		real(mp) :: L(3) = (/ 1.0_mp, 1.0_mp, 1.0_mp /)					!cubic size
!		real(mp), dimension(N(1),N(2),N(3)) :: rho, rhs, phi, weight, rhos
!		real(mp), dimension(N(1),N(2),N(3),3) :: E,Es
!		real(mp) :: rhok(N(1),N(2),N(3),Nf)
!
!		real(mp) :: A(Nf), dA(Nf)										!coefficients of fourier modes
!		real(mp) :: J0,J1,dJdA(Nf)
!
!		complex(mp) :: W(N(1),N(2),N(3))
!		integer :: i,j,k,m
!
!		dx = L/N
!		xg = (/ (i*L(1)/N(1), i=0,N(1)-1) /)
!		yg = (/ (i*L(2)/N(2), i=0,N(2)-1) /)
!		zg = (/ (i*L(3)/N(3), i=0,N(3)-1) /)
!
!		A = 1.0_mp
!		rho = 0.0_mp
!		rhok = 0.0_mp
!		do i=1,Nf
!			do k=1,N(3)
!				do j=1,N(2)
!					do m=1,N(1)
!						rhok(m,j,k,i) = SIN(2.0_mp*pi*xg(m)/L(1)*i)	&
!									*SIN(2.0_mp*pi*yg(j)/L(2)*i)	&
!									*SIN(2.0_mp*pi*zg(k)/L(3)*i)
!					end do
!				end do
!			end do
!			rho = rho + A(i)*rhok(:,:,:,i)
!		end do
!
!		weight = 1.0_mp
!
!		call FFTPoisson_setup(N,W,L)
!
!		rhs = -rho/eps0
!
!		call FFTPoisson(phi,rhs,W)
!
!		E = - Gradient(phi,dx,N)
!
!		J0 = SUM( PRODUCT(dx)*weight*(E(:,:,:,1)**2 + E(:,:,:,2)**2 + E(:,:,:,3)**2) )
!		print *, 'J0 = ', J0
!
!		do i=1,3
!			Es(:,:,:,i) = -2.0_mp*PRODUCT(dx)*weight*E(:,:,:,i)
!		end do
!
!		call FFTAdj(Es,rhos,W,dx)
!		rhos = rhos/eps0
!
!		do i=1,Nf
!			dJdA(i) = SUM( rhos*rhok(:,:,:,i) )
!		end do
!		print *, 'dJdA = ', dJdA
!
!		dA = 0.0_mp
!		j = 1
!		dA(j) = (0.1_mp)**9
!		A = A+dA
!		rho = 0.0_mp
!		do i=1,Nf
!			rho = rho + A(i)*rhok(:,:,:,i)
!		end do
!
!		rhs = -rho/eps0
!
!		call FFTPoisson(phi,rhs,W)
!
!		E = - Gradient(phi,dx,N)
!
!		J1 = SUM( PRODUCT(dx)*weight*(E(:,:,:,1)**2 + E(:,:,:,2)**2 + E(:,:,:,3)**2) )
!		print *, 'J1 = ', J1
!		print *, 'actual dJdA = ', (J1-J0)/dA(j)
!	end subroutine
!
!	subroutine verify_assignment()
!		type(PM3D) :: this
!		real(mp) :: Tf=40.0_mp,Ti=20.0_mp,rho_back
!		integer :: Ng(3), N
!
!		real(mp) :: testqs(2), testms(2), testxp(2,3), testvp(2,3)
!		integer :: testN
!		integer :: i,j, g(2,3)
!
!		Ng = (/ 64, 64, 64 /)
!		testN = 2
!
!		call buildPM3D(this,Tf,Ti,Ng,testN)
!		testqs(1) = -this%wp*this%wp/(2/PRODUCT(this%L))
!		testqs(2) = -testqs(1)
!		rho_back = 0.0_mp
!		testms = ABS(testqs(1))
!		testxp(1,:) = (/ 0.4_mp, 0.6_mp, 0.7_mp /)*this%L
!		testxp(2,:) = (/ 1.0_mp-1.0_mp/Ng(1), 0.0_mp, 0.0_mp /)*this%L
!		testvp = 0.0_mp
!
!		call setPlasma(this%p,testxp,testvp,testqs,testms)
!		call setMesh(this%m,rho_back)
!		call assignMatrix(this%a,this%p,this%m,this%p%xp)
!
!		do i=1,testN
!			print *, 'particle ',i,' (column-major order)'
!			do j=1,8
!!				print *, 2-MOD(j,2),MOD((j-1)/2,2)+1,MOD((j-1)/4,2)+1
!				print *, 'g(',this%a%g(i,2-MOD(j,2),1),',',this%a%g(i,MOD((j-1)/2,2)+1,2),',',this%a%g(i,MOD((j-1)/4,2)+1,3),')=',	&
!				this%a%frac(i,2-MOD(j,2),MOD((j-1)/2,2)+1,MOD((j-1)/4,2)+1)
!			end do
!			print *, 'sum of fraction = ', SUM(this%a%frac(i,:,:,:))
!		end do
!
!		call chargeAssign(this%a,this%p,this%m)
!
!		print *, 'Check on the grid'
!		do i=1, testN
!			print *, 'Particle ',i
!			g = this%a%g(i,:,:)
!			print *, 'Xg=',g(:,1)
!			print *, 'Yg=',g(:,2)
!			print *, 'Zg=',g(:,3)
!			print *, 'Particle fraction(column-major order)'
!			print *, this%m%rho( g(:,1), g(:,2), g(:,3) )*PRODUCT(this%m%dx)/this%p%qs(i)
!		end do
!
!		call destroyPM3D(this)
!	end subroutine

end module
