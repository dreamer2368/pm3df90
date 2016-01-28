module timeStep

	use modPM3D

	implicit none

contains

	subroutine halfStep(this)
		type(PM3D), intent(inout) :: this
		integer :: i, j
		real(mp) :: rhs(this%ng(1),this%ng(2),this%ng(3))
		real(mp) :: dt,L(3),N
		dt = this%dt
		L = this%L
		N = this%n

		call assignMatrix(this%a,this%p,this%m,this%p%xp)

		!charge assignment
		call chargeAssign(this%a,this%p,this%m)

		!field equation rhs
		rhs = -this%m%rho/this%eps0
		!solve field equation
		call FFTPoisson(this%m%phi,rhs,this%m%W)				!FFT-based Poisson solver
!
!		!Electric field : -D*phi
		this%m%E = - Gradient(this%m%phi,this%m%dx, this%ng)

		!Force assignment : mat'*E
		call forceAssign(this%a,this%p,this%m)

		!Half time step advancement in velocity
		this%p%vp(:,1) = this%p%vp(:,1) + dt/2.0_mp*this%p%qs/this%p%ms*this%p%Ep(:,1)
		this%p%vp(:,2) = this%p%vp(:,2) + dt/2.0_mp*this%p%qs/this%p%ms*this%p%Ep(:,2)
		this%p%vp(:,3) = this%p%vp(:,3) + dt/2.0_mp*this%p%qs/this%p%ms*this%p%Ep(:,3)
	end subroutine

	subroutine updatePlasma(this)
		type(PM3D), intent(inout) :: this
		integer :: i, j, k
		real(mp) :: rhs(this%ng(1),this%ng(2),this%ng(3))
		real(mp) :: dt,L(3),N,B
		real(mp) :: start, finish
		dt = this%dt
		L = this%L
		N = this%n
		B = this%B0

		do k = 1,this%nt
			call recordPlasma(this%r,this%p,this%m,k)									!record for n=0~(Nt-1)

			if(k==this%ni) then
				this%p%xp(:,1) = this%p%xp(:,1) + dt*B*L(1)/N*SIN(4.0_mp*pi*this%p%xp(:,1)/L(1)*mode)		!xp_(Ni-1) : k=Ni, xp_Ni : k=(Ni+1)
			end if
			this%p%xp = this%p%xp + dt*this%p%vp

			call assignMatrix(this%a,this%p,this%m,this%p%xp)

			!Sometimes time stepping continues even though particles are assigned way outside the boundary. Even though it continues without error message, it's wrong.
			if( MAXVAL(this%a%g)>MAXVAL(this%ng) .or. MINVAL(this%a%g)<1 ) then
				print *, 'n= ', k, MAXVAL(this%a%g), MINVAL(this%a%g)
				print *, 'Boundary handling is failed. particle is way outside BC. stopped time stepping.'
				stop
			end if

			!charge assignment
			call chargeAssign(this%a,this%p,this%m)

			!field equation rhs
			rhs = -this%m%rho/this%eps0
			!solve field equation
			call FFTPoisson(this%m%phi,rhs,this%m%W)				!FFT-based Poisson solver

			!Electric field : -D*phi
			this%m%E = - Gradient(this%m%phi,this%m%dx, this%ng)

			!Force assignment : mat'*E
			call forceAssign(this%a,this%p,this%m)

			!Single time step advancement in velocity
			this%p%vp(:,1) = this%p%vp(:,1) + dt*this%p%qs/this%p%ms*this%p%Ep(:,1)
			this%p%vp(:,2) = this%p%vp(:,2) + dt*this%p%qs/this%p%ms*this%p%Ep(:,2)
			this%p%vp(:,3) = this%p%vp(:,3) + dt*this%p%qs/this%p%ms*this%p%Ep(:,3)
		end do
	end subroutine

	subroutine QOI(this,J)
		type(PM3D), intent(in) :: this
		real(mp), intent(out) :: J
		integer :: N,Nt,Ni
		N = this%n
		Nt = this%nt
		Ni = this%ni

		J = 1.0_mp/N/(Nt-Ni)*SUM(this%r%vpdata(:,:,Ni+1:Nt)**2)		!omitted 1/N/T for the sake of machine precision
	end subroutine

	subroutine QOItwo(this,A,B,J)
		type(PM3D), intent(in) :: this
		integer, intent(in) :: A, B
		real(mp), intent(out) :: J
		integer :: N
		N = this%n

		J = 1.0_mp/N/(B-A)*SUM(this%r%vpdata(:,:,A+1:B)**2)		!omitted 1/N/T for the sake of machine precision
	end subroutine

	subroutine forwardsweep(this,xp0,vp0,qs,ms,rho_back)
		type(PM3D), intent(inout) :: this
		real(mp), dimension(this%n,3), intent(in) :: xp0, vp0
		real(mp), intent(in) :: qs(this%n), ms(this%n), rho_back

		!initial condition
		call setPlasma(this%p,xp0,vp0,qs,ms)
		call setMesh(this%m,rho_back)

		!Time stepping
		call halfStep(this)
		call updatePlasma(this)
		!export the data - when needed
		call printPlasma(this%r)
	end subroutine
!
!	subroutine adjoint(this,delJdelA)
!		type(plasma), intent(inout) :: this
!		real(mp), intent(out) :: delJdelA
!		real(mp) :: dts
!		real(mp) :: xps(size(this%xp)), vps(size(this%vp)), Es(size(this%E)), phis(size(this%phi))
!		real(mp) :: dvps(size(this%vp))
!		real(mp) :: rhs(size(this%E)), phis1(size(this%phi)-1)
!		integer :: i,k, nk
!		real(mp) :: coeff(this%n)
!		dts = -dt
!		xps = 0.0_mp
!		vps = 0.0_mp
!
!		do k=1,this%nt
!			!vps update
!			nk = this%nt+1-k
!			dvps = merge(2.0_mp/N/(Nt-Ni)/dt*this%vpdata(:,this%nt+1-k),0.0_mp,nk>=Ni+1)
!			if (k==this%nt) then
!				vps = 1.0_mp/2.0_mp*( vps + dts*( -xps + dvps ) )	!omitted 1/N/T(see QOI function)
!			else
!				vps = vps + dts*( -xps + dvps )
!			end if
!			!assignment
!			call assignMatrix(this,this%xpdata(:,nk))
!
!			!Es update : (qe/me/dx)*mat*vps
!			call vpsAssign(this,Es,vps)
!
!			!phis update
!			!rhs = D*Es
!			rhs = multiplyD(Es,dx)
!
!			!K*phis = D*Es
!			call CG_K(phis1,rhs(1:Ng-1),dx)				!Conjugate-Gradient Iteration
!			phis(1:Ng-1) = phis1
!			phis(Ng) = 0.0_mp
!
!			!xps update
!			if(nk==this%ni) then
!				xps = xps - dts*xps*(B*L/N*4.0_mp*pi*mode/L)*COS(4.0_mp*pi*mode*this%xpdata(:,nk)/L)
!			end if
!			xps = xps + xpsUpdate(this,vps,phis,dts,k)
!
!			call recordAdjoint(this,nk,xps)					!recorde xps from n=(Nt-1) to n=0
!			if( nk<this%ni ) then
!				exit
!			end if
!		end do
!
!		coeff = - L/N*SIN( 4.0_mp*pi*mode*this%xpdata(:,this%ni)/L )
!		delJdelA = 0.0_mp
!		do i = 1,N
!			delJdelA = delJdelA + dt*coeff(i)*this%xpsdata(i,this%ni+1)
!		end do
!	end subroutine

end module
