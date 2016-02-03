module testmodule

	use init
	use timeStep

	implicit none

contains

	subroutine test_particle_adj2(N,Np)
		integer, intent(in) :: N(3)										!!grid number
		integer, intent(in) :: Np										!!number of particles
		type(plasma) :: p
		type(mesh) :: m
		type(pmassign) ::a

		real(mp) :: dx(3)												!!grid size
		real(mp) :: eps0 = 1.0_mp, rho_back, qe, qs(Np), ms(Np)
		real(mp) :: L(3) = (/ 2.0_mp, 2.0_mp, 2.0_mp /)					!cubic size
		real(mp), dimension(N(1),N(2),N(3)) :: rhs, phi, weight, rhos
		real(mp), dimension(N(1),N(2),N(3),3) :: Es
		real(mp), dimension(Np,3) :: xp, vp, xps, xps1, xps2, Eps
		real(mp) :: dxp

		real(mp) :: J0,J1,dJdxp(Np,3)
		real(mp) :: fxp = (0.1_mp)**9
		integer :: i,j,k(2)

		call buildPlasma(p,Np)
		call buildMesh(m,L,N)
		call buildAssign(a,Np,N)

		!particle, mesh setup
		print *, 'Original xp'
		do j=1,Np
		xp(j,:) = 0.1_mp*L*(/ j, 7*j, 4*j /)
		print *, xp(j,:)
		end do
		qe = -0.1_mp
		qs = qe
		ms = -qs
		rho_back = -qe*p%n/PRODUCT(L)
		vp = 0.0_mp
		call setPlasma(p,xp,vp,qs,ms)
		call setMesh(m,rho_back)
		!assignment
		call assignMatrix(a,p,m,p%xp)
		call chargeAssign(a,p,m)
		!field solver
		rhs = -m%rho/eps0
		call FFTPoisson(m%phi,rhs,m%W)
		m%E = - Gradient(m%phi,m%dx,m%ng)
		!force assignment
		call forceAssign(a,p,m)
		!QoI evaluation
		J0 = SUM( (p%Ep(:,1)**2 + p%Ep(:,2)**2 + p%Ep(:,3)**2) )
		print *, 'J0 = ', J0

		!Adjoint sensitivity solver
		Eps = - 2.0_mp*p%Ep

		call Adj_chargeAssign(a,Eps,Es)

		call FFTAdj(Es,rhos,m%W,m%dx)
		rhos = -rhos/eps0

		call Adj_rhosAssign(a,p,m,rhos,xps1)
		call Adj_EpsAssign(a,m,Eps,xps2)
		xps = xps1 + xps2

		print *, 'dJdxp'
		do i=1,Np
			print *, xps(i,:)
		end do

		!FD approximation - choose the component that you want to measure
		k = (/2,3/)
		dxp = xp(k(1),k(2))*fxp
		xp(k(1),k(2)) = xp(k(1),k(2)) + dxp
		print *, 'Perturbed xp'
		do j=1,Np
		print *, xp(j,:)
		end do
		call setPlasma(p,xp,vp,qs,ms)
		!assignment
		call assignMatrix(a,p,m,p%xp)
		call chargeAssign(a,p,m)
		!field solver
		rhs = -m%rho/eps0
		call FFTPoisson(m%phi,rhs,m%W)
		m%E = - Gradient(m%phi,m%dx,m%ng)
		!force assignment
		call forceAssign(a,p,m)
		!QoI evaluation
		J1 = SUM( (p%Ep(:,1)**2 + p%Ep(:,2)**2 + p%Ep(:,3)**2) )
		print *, 'J1 = ', J1
		print *, 'dJdxp(',k(1),',',k(2),')=', (J1-J0)/dxp

		call destroyPlasma(p)
		call destroyMesh(m)
		call destroyAssign(a)
	end subroutine

	subroutine test_particle_adj(N,Np)
		integer, intent(in) :: N(3)										!!grid number
		integer, intent(in) :: Np										!!number of particles
		type(plasma) :: p
		type(mesh) :: m
		type(pmassign) ::a

		real(mp) :: dx(3)												!!grid size
		real(mp) :: eps0 = 1.0_mp, rho_back, qe, qs(Np), ms(Np)
		real(mp) :: L(3) = (/ 2.0_mp, 2.0_mp, 2.0_mp /)					!cubic size
		real(mp), dimension(N(1),N(2),N(3)) :: rhs, phi, weight, rhos
		real(mp), dimension(N(1),N(2),N(3),3) :: Es
		real(mp) :: xp(Np,3), vp(Np,3), dxp, xps(Np,3)

		real(mp) :: J0,J1,dJdxp(Np,3)
		real(mp) :: fxp = (0.1_mp)**9
		integer :: i,j,k(2)

		call buildPlasma(p,Np)
		call buildMesh(m,L,N)
		call buildAssign(a,Np,N)

	!particle, mesh setup
		print *, 'Original xp'
		do j=1,Np
			xp(j,:) = 0.1_mp*L*(/ j, 7*j, 4*j /)
			print *, xp(j,:)
		end do
		qe = -0.1_mp
		qs = qe
		ms = -qs
		rho_back = -qe*p%n/PRODUCT(L)
		vp = 0.0_mp
		call setPlasma(p,xp,vp,qs,ms)
		call setMesh(m,rho_back)
	!assignment
		call assignMatrix(a,p,m,p%xp)
		call chargeAssign(a,p,m)
	!field solver
		rhs = -m%rho/eps0
		call FFTPoisson(m%phi,rhs,m%W)
		m%E = - Gradient(m%phi,m%dx,m%ng)
	!QoI evaluation
		weight = 0.0_mp
		weight(2*N(1)/5:3*N(1)/5,2*N(2)/5:3*N(2)/5,2*N(3)/5:3*N(3)/5) = 1.0_mp
		J0 = SUM( PRODUCT(m%dx)*weight*(m%E(:,:,:,1)**2 + m%E(:,:,:,2)**2 + m%E(:,:,:,3)**2) )
		print *, 'J0 = ', J0

	!Adjoint sensitivity solver
		do i=1,3
			Es(:,:,:,i) = -2.0_mp*PRODUCT(m%dx)*weight*m%E(:,:,:,i)
		end do

		call FFTAdj(Es,rhos,m%W,m%dx)
		rhos = -rhos/eps0

		call Adj_rhosAssign(a,p,m,rhos,xps)

		print *, 'dJdxp'
		do i=1,Np
			print *, xps(i,:)
		end do

	!FD approximation - choose the component that you want to measure
		k = (/1,2/)
		dxp = xp(k(1),k(2))*fxp
		xp(k(1),k(2)) = xp(k(1),k(2)) + dxp
		print *, 'Perturbed xp'
		do j=1,Np
			print *, xp(j,:)
		end do
		call setPlasma(p,xp,vp,qs,ms)
	!assignment
		call assignMatrix(a,p,m,p%xp)
		call chargeAssign(a,p,m)
	!field solver
		rhs = -m%rho/eps0
		call FFTPoisson(m%phi,rhs,m%W)
		m%E = - Gradient(m%phi,m%dx,m%ng)
	!QoI re-evaluation, Sensitivity approximation
		J1 = SUM( PRODUCT(m%dx)*weight*(m%E(:,:,:,1)**2 + m%E(:,:,:,2)**2 + m%E(:,:,:,3)**2) )
		print *, 'J1 = ', J1
		print *, 'dJdxp(',k(1),',',k(2),')=', (J1-J0)/dxp

		call destroyPlasma(p)
		call destroyMesh(m)
		call destroyAssign(a)
	end subroutine

	subroutine twostream(v0,Ng,Nd)
		real(mp), intent(in) :: v0
		integer, intent(in) :: Ng(3), Nd(3)
		type(PM3D) :: this
		real(mp) :: Tf=60.0_mp,Ti=20.0_mp,rho_back
		integer :: N
		real(mp) :: xp0(PRODUCT(Nd),3), vp0(PRODUCT(Nd),3), qs(PRODUCT(Nd)), ms(PRODUCT(Nd))

		N = PRODUCT(Nd)
		call buildPM3D(this,Tf,Ti,Ng,N)

		call particle_initialize(this,Nd,v0,xp0,vp0,qs,ms,rho_back)
		call forwardsweep(this,xp0,vp0,qs,ms,rho_back)

		call destroyPM3D(this)
	end subroutine

	subroutine test_FFTPoisson_adj(N,Nf)
		integer, intent(in) :: N(3)										!!grid number
		integer, intent(in) :: Nf										!!number of fourier modes

		real(mp) :: dx(3)												!!grid size
		real(mp) :: xg(N(1)), yg(N(2)), zg(N(3))
		real(mp) :: eps0 = 1.0_mp
		real(mp) :: L(3) = (/ 1.0_mp, 1.0_mp, 1.0_mp /)					!cubic size
		real(mp), dimension(N(1),N(2),N(3)) :: rho, rhs, phi, weight, rhos
		real(mp), dimension(N(1),N(2),N(3),3) :: E,Es
		real(mp) :: rhok(N(1),N(2),N(3),Nf)

		real(mp) :: A(Nf), dA(Nf)										!coefficients of fourier modes
		real(mp) :: J0,J1,dJdA(Nf)

		complex(mp) :: W(N(1),N(2),N(3))
		integer :: i,j,k,m

		dx = L/N
		xg = (/ (i*L(1)/N(1), i=0,N(1)-1) /)
		yg = (/ (i*L(2)/N(2), i=0,N(2)-1) /)
		zg = (/ (i*L(3)/N(3), i=0,N(3)-1) /)

		A = 1.0_mp
		rho = 0.0_mp
		rhok = 0.0_mp
		do i=1,Nf
			do k=1,N(3)
				do j=1,N(2)
					do m=1,N(1)
						rhok(m,j,k,i) = SIN(2.0_mp*pi*xg(m)/L(1)*i)	&
									*SIN(2.0_mp*pi*yg(j)/L(2)*i)	&
									*SIN(2.0_mp*pi*zg(k)/L(3)*i)
					end do
				end do
			end do
			rho = rho + A(i)*rhok(:,:,:,i)
		end do

		weight = 1.0_mp

		call FFTPoisson_setup(N,W,L)

		rhs = -rho/eps0

		call FFTPoisson(phi,rhs,W)

		E = - Gradient(phi,dx,N)

		J0 = SUM( PRODUCT(dx)*weight*(E(:,:,:,1)**2 + E(:,:,:,2)**2 + E(:,:,:,3)**2) )
		print *, 'J0 = ', J0

		do i=1,3
			Es(:,:,:,i) = -2.0_mp*PRODUCT(dx)*weight*E(:,:,:,i)
		end do

		call FFTAdj(Es,rhos,W,dx)
		rhos = rhos/eps0

		do i=1,Nf
			dJdA(i) = SUM( rhos*rhok(:,:,:,i) )
		end do
		print *, 'dJdA = ', dJdA

		dA = 0.0_mp
		j = 1
		dA(j) = (0.1_mp)**9
		A = A+dA
		rho = 0.0_mp
		do i=1,Nf
			rho = rho + A(i)*rhok(:,:,:,i)
		end do

		rhs = -rho/eps0

		call FFTPoisson(phi,rhs,W)

		E = - Gradient(phi,dx,N)

		J1 = SUM( PRODUCT(dx)*weight*(E(:,:,:,1)**2 + E(:,:,:,2)**2 + E(:,:,:,3)**2) )
		print *, 'J1 = ', J1
		print *, 'actual dJdA = ', (J1-J0)/dA(j)
	end subroutine

	subroutine verify_assignment()
		type(PM3D) :: this
		real(mp) :: Tf=40.0_mp,Ti=20.0_mp,rho_back
		integer :: Ng(3), N

		real(mp) :: testqs(2), testms(2), testxp(2,3), testvp(2,3)
		integer :: testN
		integer :: i,j, g(2,3)

		Ng = (/ 64, 64, 64 /)
		testN = 2

		call buildPM3D(this,Tf,Ti,Ng,testN)
		testqs(1) = -this%wp*this%wp/(2/PRODUCT(this%L))
		testqs(2) = -testqs(1)
		rho_back = 0.0_mp
		testms = ABS(testqs(1))
		testxp(1,:) = (/ 0.4_mp, 0.6_mp, 0.7_mp /)*this%L
		testxp(2,:) = (/ 1.0_mp-1.0_mp/Ng(1), 0.0_mp, 0.0_mp /)*this%L
		testvp = 0.0_mp

		call setPlasma(this%p,testxp,testvp,testqs,testms)
		call setMesh(this%m,rho_back)
		call assignMatrix(this%a,this%p,this%m,this%p%xp)

		do i=1,testN
			print *, 'particle ',i,' (column-major order)'
			do j=1,8
!				print *, 2-MOD(j,2),MOD((j-1)/2,2)+1,MOD((j-1)/4,2)+1
				print *, 'g(',this%a%g(i,2-MOD(j,2),1),',',this%a%g(i,MOD((j-1)/2,2)+1,2),',',this%a%g(i,MOD((j-1)/4,2)+1,3),')=',	&
				this%a%frac(i,2-MOD(j,2),MOD((j-1)/2,2)+1,MOD((j-1)/4,2)+1)
			end do
			print *, 'sum of fraction = ', SUM(this%a%frac(i,:,:,:))
		end do

		call chargeAssign(this%a,this%p,this%m)

		print *, 'Check on the grid'
		do i=1, testN
			print *, 'Particle ',i
			g = this%a%g(i,:,:)
			print *, 'Xg=',g(:,1)
			print *, 'Yg=',g(:,2)
			print *, 'Zg=',g(:,3)
			print *, 'Particle fraction(column-major order)'
			print *, this%m%rho( g(:,1), g(:,2), g(:,3) )*PRODUCT(this%m%dx)/this%p%qs(i)
		end do

		call destroyPM3D(this)
	end subroutine

end module
