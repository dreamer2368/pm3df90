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
		dt = this%dt
		L = this%L
		N = this%n
		B = this%B0

		do k = 1,this%nt
!			call recordPlasma(this%r,this%p,this%m,k)									!record for n=0~(Nt-1) : xp_(k-1 and vp_(k-1/2)

			if(k==this%ni) then
				this%p%xp(:,1) = this%p%xp(:,1) + dt*B*L(1)/N*SIN(4.0_mp*pi*this%p%xp(:,1)/L(1)*mode)		!xp_(Ni-1) : k=Ni, xp_Ni : k=(Ni+1)
			end if
			this%p%xp = this%p%xp + dt*this%p%vp

			call assignMatrix(this%a,this%p,this%m,this%p%xp)

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

			call recordPlasma(this%r,this%p,this%m,k)									!record for n=1~Nt : xp_k and vp_(k+1/2)
		end do
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
!		call printPlasma(this%r)
	end subroutine

end module
