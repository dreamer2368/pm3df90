module timeStep

	use modQoI

	implicit none

contains

	subroutine forwardsweep(this,xp0,vp0,qs,ms,rho_back,source)
		type(PM3D), intent(inout) :: this
		real(mp), dimension(this%n,3), intent(in) :: xp0, vp0
		real(mp), intent(in) :: qs(this%n), ms(this%n), rho_back
		integer :: k
		interface
			subroutine source(pm,k,str)
				use modPM3D
				type(PM3D), intent(inout) :: pm
				integer, intent(in) :: k
				character(len=*), intent(in) :: str
			end subroutine
		end interface

		!initial condition
		call setPlasma(this%p,xp0,vp0,qs,ms)
		call setMesh(this%m,rho_back)

		!Time stepping
		call halfStep(this)
		do k=1,this%nt
			call updatePlasma(this,source,k)
		end do
	end subroutine

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

!		!Electric field : -D*phi
		this%m%E = - Gradient(this%m%phi,this%m%dx, this%ng)

		!Force assignment : mat'*E
		call forceAssign(this%a,this%p,this%m)

		!Half time step advancement in velocity
		this%p%vp(:,1) = this%p%vp(:,1) + dt/2.0_mp*this%p%qs/this%p%ms*this%p%Ep(:,1)
		this%p%vp(:,2) = this%p%vp(:,2) + dt/2.0_mp*this%p%qs/this%p%ms*this%p%Ep(:,2)
		this%p%vp(:,3) = this%p%vp(:,3) + dt/2.0_mp*this%p%qs/this%p%ms*this%p%Ep(:,3)
	end subroutine

	subroutine updatePlasma(this,source,k)
		type(PM3D), intent(inout) :: this
		integer, intent(in) :: k
		integer :: i, j
		real(mp) :: rhs(this%ng(1),this%ng(2),this%ng(3))
		real(mp) :: dt,L(3),N,B
		interface
			subroutine source(pm,k,str)
				use modPM3D
				type(PM3D), intent(inout) :: pm
				integer, intent(in) :: k
				character(len=*), intent(in) :: str
			end subroutine
		end interface
		dt = this%dt
		L = this%L
		N = this%n
		B = this%B0

!		call recordPlasma(this%r,this%p,this%m,k)									!record for n=0~(Nt-1) : xp_(k-1 and vp_(k-1/2)

		call source(this,k,'xp')
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
		call source(this,k,'E')

		!Force assignment : mat'*E
		call forceAssign(this%a,this%p,this%m)

		!Single time step advancement in velocity
		call source(this,k,'vp')
		this%p%vp(:,1) = this%p%vp(:,1) + dt*this%p%qs/this%p%ms*this%p%Ep(:,1)
		this%p%vp(:,2) = this%p%vp(:,2) + dt*this%p%qs/this%p%ms*this%p%Ep(:,2)
		this%p%vp(:,3) = this%p%vp(:,3) + dt*this%p%qs/this%p%ms*this%p%Ep(:,3)

		call recordPlasma(this%r,this%p,this%m,k)									!record for n=1~Nt : xp_k and vp_(k+1/2)
	end subroutine

	subroutine Forcing(this,k,str)
		type(PM3D), intent(inout) :: this
		integer, intent(in) :: k
		character(len=*), intent(in) :: str

		SELECT CASE (str)
			CASE('xp')
				if(k==this%ni) then
					this%p%xp(:,1) = this%p%xp(:,1) + this%dt*this%B0*this%L(1)/this%n*SIN(4.0_mp*pi*this%p%xp(:,1)/this%L(1)*mode)		!xp_(Ni-1) : k=Ni, xp_Ni : k=(Ni+1)
				end if
		END SELECT
	end subroutine

	subroutine QOItwo(this,A,B,J)
		type(PM3D), intent(in) :: this
		integer, intent(in) :: A, B
		real(mp), intent(out) :: J
		integer :: N
		N = this%n

		J = 1.0_mp/N/(B-A)*SUM(this%r%vpdata(:,:,A+1:B)**2)		!omitted 1/N/T for the sake of machine precision
	end subroutine

!===================Adjoint time stepping==========================

	subroutine backward_sweep(adj,pm, dJ, Dsource, dJdA)
		type(adjoint), intent(inout) :: adj
		type(PM3D), intent(inout) :: pm
		real(mp), intent(out) :: dJdA
		integer :: k, nk
		interface
			subroutine dJ(adj,pm,k)
				use modPM3D
				use modAdj
				type(adjoint), intent(inout) :: adj
				type(PM3D), intent(in) :: pm
				integer, intent(in) :: k
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

		adj%xps = 0.0_mp
		adj%vps = 0.0_mp
		do k=1,pm%nt
			nk = pm%nt+1-k

			call reset_Dadj(adj)

			!======= dJ : source term ==========
			call dJ(adj,pm,nk)

			!======= dv_p =============
			call Dsource(adj,pm,nk,'vp')
			call Adj_accel(adj,pm,adj%dvps)

!			!Check when adjoint reach to the initial step
!			if( k .eq. pm%nt ) then
!				adj%vps = 2.0_mp*adj%vps
!			end if

			!======= dE_p =============
			adj%Eps(:,1) = pm%p%qs/pm%p%ms*adj%vps(:,1)
			adj%Eps(:,2) = pm%p%qs/pm%p%ms*adj%vps(:,2)
			adj%Eps(:,3) = pm%p%qs/pm%p%ms*adj%vps(:,3)

			call assignMatrix(pm%a,pm%p,pm%m,pm%r%xpdata(:,:,nk))

			!======= dE_g =============
			call Adj_forceAssign_E(pm%a,adj%Eps,adj%Es)

			!======= dPhi_g, dRho_g =============
			call FFTAdj(adj%Es,adj%rhos,pm%m%W,pm%m%dx)
			adj%rhos = -adj%rhos/pm%eps0

			!======= dx_p =============
			call Adj_chargeAssign(pm%a,pm%p,pm%m,adj%rhos,adj%dxps)
			call Adj_forceAssign_xp(pm%a,pm%m,pm%r%Edata(:,:,:,:,nk),adj%Eps,adj%dxps)
			call Dsource(adj,pm,nk,'xp')
			call Adj_move(adj,pm,adj%dxps)

			call recordAdjoint(pm%r,nk,adj%xps)									!record for nk=1~Nt : xps_nk and vp_(nk+1/2)

			if( nk<pm%ni ) then
				exit
			end if
		end do

		dJdA = SUM( - pm%L(1)/pm%n*SIN( 4.0_mp*pi*mode*pm%r%xpdata(:,1,pm%ni-1)/pm%L(1) )*pm%r%xpsdata(:,1,pm%ni) )
	end subroutine

	subroutine Dforcing(adj,pm,k,str)
		type(adjoint), intent(inout) :: adj
		type(PM3D), intent(in) :: pm
		integer, intent(in) :: k
		character(len=*), intent(in) :: str

		select case (str)
			case ('xp')
				if( k .eq. pm%ni-1 ) then
					adj%dxps = adj%dxps - adj%xps*(pm%B0*pm%L(1)/pm%n*4.0_mp*pi*mode/pm%L(1))*COS(4.0_mp*pi*mode*pm%r%xpdata(:,:,k)/pm%L(1))
				end if
		end select
	end subroutine

end module
