module timeStep

	use modQoI
	use modControl
	use modRecord

	implicit none

contains

	subroutine forwardsweep(this,control,QoI,J)
		type(PM3D), intent(inout) :: this
		real(mp), intent(out), optional :: J
!		real(mp), dimension(this%n,3), intent(in) :: xp0, vp0
!		real(mp), intent(in) :: qs(this%n), ms(this%n), rho_back
		integer :: k
		interface
			subroutine control(pm,k,str)
				use modPM3D
				type(PM3D), intent(inout) :: pm
				integer, intent(in) :: k
				character(len=*), intent(in) :: str
			end subroutine
		end interface
		interface
			subroutine QoI(pm,k,J)
				use modQoI
				type(PM3D), intent(in) :: pm
				real(mp), intent(inout) :: J
				integer, intent(in) :: k
			end subroutine
		end interface
		optional :: QoI
		if( present(J) ) then
			J = 0.0_mp
		end if
		k=0

		!initial condition
!		call setPlasma(this%p,xp0,vp0,qs,ms)
!		call setMesh(this%m,rho_back)

		!Time stepping
!		call halfStep(this)
		if( present(J) ) then
			call QoI(this,k,J)
		end if
		call recordPlasma(this%r,this,k)							!x_0, v_1/2							
		do k=1,this%nt
			call updatePlasma(this,control,k)
			if( present(J) ) then
				call QoI(this,k,J)
			end if
			call recordPlasma(this%r,this,k)						!x_k, v_{k+1/2}
		end do
	end subroutine

	subroutine halfStep(this)
		type(PM3D), intent(inout) :: this
		integer :: i, j, Ns
		real(mp) :: rhs(this%ng(1),this%ng(2),this%ng(3))
		real(mp) :: dt,L(3)
		dt = this%dt
		L = this%L
		Ns = this%ns

		do i=1,Ns
			call assignMatrix(this%a(i),this%p(i),this%m,this%p(i)%xp)
		end do

		!charge assignment
		call chargeAssign(this%a,this%p,this%m)

		!field equation rhs
		rhs = -this%m%rho/this%eps0

		!solve field equation
		call FFTPoisson(this%m%phi,rhs,this%m%W)				!FFT-based Poisson solver

!		!Electric field : -D*phi
		this%m%E = - Gradient(this%m%phi,this%m%dx, this%ng)

		!Force assignment : mat'*E
		do i=1,Ns
			call forceAssign(this%a(i),this%p(i),this%m)
		end do
		!Half time step advancement in velocity
		do i=1,Ns
			call accel(this%p(i),0.5_mp*dt)
		end do
	end subroutine

	subroutine updatePlasma(this,control,k)
		type(PM3D), intent(inout) :: this
		integer, intent(in) :: k
		integer :: i, j,Ns
		real(mp) :: rhs(this%ng(1),this%ng(2),this%ng(3))
		real(mp) :: dt,L(3)
		interface
			subroutine control(pm,k,str)
				use modPM3D
				type(PM3D), intent(inout) :: pm
				integer, intent(in) :: k
				character(len=*), intent(in) :: str
			end subroutine
		end interface
		dt = this%dt
		L = this%L
		Ns = this%ns

!		call recordPlasma(this%r,this%p,this%m,k)									!record for n=0~(Nt-1) : xp_(k-1 and vp_(k-1/2)

		call control(this,k,'xp')
		do i=1,Ns
			call move(this%p(i),dt)
		end do

		do i=1,Ns
			call assignMatrix(this%a(i),this%p(i),this%m,this%p(i)%xp)
		end do

		!charge assignment
		call chargeAssign(this%a,this%p,this%m)

		!field equation rhs
		rhs = -this%m%rho/this%eps0
		!solve field equation
		call FFTPoisson(this%m%phi,rhs,this%m%W)				!FFT-based Poisson solver

		!Electric field : -D*phi
		this%m%E = - Gradient(this%m%phi,this%m%dx, this%ng)
		call control(this,k,'E')

		!Force assignment : mat'*E
		do i=1,Ns
			call forceAssign(this%a(i),this%p(i),this%m)
		end do

		!Single time step advancement in velocity
		call control(this,k,'vp')
		do i=1,Ns
			call accel(this%p(i),dt)
		end do
	end subroutine

!===================Adjoint time stepping==========================

	subroutine backward_sweep(adj,pm, grad, dJ, Dcontrol, dJdA, control)
		type(adjoint), intent(inout) :: adj
		type(PM3D), intent(inout) :: pm
		real(mp), intent(out) :: grad(:)
		integer :: k, nk, i
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
			subroutine dJdA(adj,pm,k,str,grad)
				use modPM3D
				use modAdj
				type(adjoint), intent(in) :: adj
				type(PM3D), intent(in) :: pm
				integer, intent(in) :: k
				character(len=*), intent(in) :: str
				real(mp), intent(inout) :: grad(:)
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
		grad = 0.0_mp

!		adj%xps = 0.0_mp
!		adj%vps = 0.0_mp
		do k=1,pm%nt
			nk = pm%nt+1-k

			call reset_Dadj(adj)

			!=====  Checkpointing  =====
			call checkpoint(pm,pm%r,nk,control)

			!===== dJdA : 1st sensitivity calculation ======
			call dJdA(adj,pm,nk,'before',grad)

			!======= dJ : source term ==========
			call dJ(adj,pm,nk)

			!======= dv_p =============
			call Dcontrol(adj,pm,nk,'vp')
			call Adj_accel(adj)

!			!Check when adjoint reach to the initial step
!			if( k .eq. pm%nt ) then
!				adj%vps = 2.0_mp*adj%vps
!			end if

			!======= dE_p =============
			do i=1,adj%ns
				adj%p(i)%Ep = pm%p(i)%qs/pm%p(i)%ms*adj%p(i)%vp
			end do

!			call assignMatrix(pm%a,pm%p,pm%m,pm%r%xpdata(:,:,nk))

			!======= dE_g =============
			adj%m%E = 0.0_mp
			do i=1,adj%ns
				call Adj_forceAssign_E(pm%a(i),adj%p(i)%Ep,adj%m%E)
			end do
			adj%m%E = adj%m%E + adj%dm%E

			!======= dPhi_g, dRho_g =============
			call FFTAdj(adj%m%E,adj%m%rho,pm%m%W,pm%m%dx)
			adj%m%rho = -adj%m%rho/pm%eps0

			!======= dx_p =============
			do i=1,adj%ns
				call Adj_chargeAssign(pm%a(i),pm%p(i),pm%m,adj%m%rho,adj%dp(i)%xp)
				call Adj_forceAssign_xp(pm%a(i),pm%m,pm%m%E,adj%p(i)%Ep,adj%dp(i)%xp)
			end do
			call Dcontrol(adj,pm,nk,'xp')
			call Adj_move(adj)
!			call recordAdjoint(pm%r,adj,nk)									!record for nk=1~Nt : xps_nk and vp_(nk+1/2)

			!===== dJdA : 2nd sensitivity calculation ======
			call dJdA(adj,pm,nk,'after',grad)
		end do
		nk = 0
		call reset_Dadj(adj)
		!=====  Checkpointing  =====
		call checkpoint(pm,pm%r,nk,control)

		!===== dJdA : 1st sensitivity calculation ======
		call dJdA(adj,pm,nk,'before',grad)

		!======= dJ : source term ==========
		call dJ(adj,pm,nk)

		!======= dv_p =============
		call Dcontrol(adj,pm,nk,'vp')
		call Adj_accel(adj)

		!======= dx_p =============
		call Dcontrol(adj,pm,nk,'xp')
		call Adj_move(adj)

		!===== dJdA : 2nd sensitivity calculation ======
		call dJdA(adj,pm,nk,'after',grad)
	end subroutine

	subroutine checkpoint(pm,r,nk,control)
		type(PM3D), intent(inout) :: pm
		type(recordData), intent(inout) :: r
		integer, intent(in) :: nk
		real(mp) :: rhs(pm%ng(1),pm%ng(2),pm%ng(3))
		integer :: kr,i
		character(len=1000) :: istr, kstr, dir_temp
		real(mp), allocatable :: xp0(:,:), vp0(:,:)
		interface
			subroutine control(pm,k,str)
				use modPM3D
				type(PM3D), intent(inout) :: pm
				integer, intent(in) :: k
				character(len=*), intent(in) :: str
			end subroutine
		end interface

		kr = merge(nk,nk/r%mod_r,r%mod_r.eq.1)
		write(kstr,*) kr
		do i=1,pm%ns
			allocate(xp0(r%np(i,kr+1),3))
			allocate(vp0(r%np(i,kr+1),3))
			write(istr,*) i
			dir_temp='data/'//r%dir//'/xp/'//trim(adjustl(kstr))//'_'//trim(adjustl(istr))//'.bin'
			open(unit=305,file=trim(dir_temp),form='unformatted',access='stream')
			dir_temp='data/'//r%dir//'/vp/'//trim(adjustl(kstr))//'_'//trim(adjustl(istr))//'.bin'
			open(unit=306,file=trim(dir_temp),form='unformatted',access='stream')
			read(305) xp0
			read(306) vp0
			close(305)
			close(306)
			call destroySpecies(pm%p(i))
			call setSpecies(pm%p(i),r%np(i,kr+1),xp0,vp0)
			deallocate(xp0)
			deallocate(vp0)
		end do
		if( nk-kr*r%mod_r.eq.0 ) then
			do i=1, pm%ns
				call assignMatrix(pm%a(i),pm%p(i),pm%m,pm%p(i)%xp)
			end do
!			call adjustGrid(pm)

			!charge assignment
			call chargeAssign(pm%a,pm%p,pm%m)

			!field equation rhs
			rhs = -pm%m%rho/pm%eps0
			!solve field equation
			call FFTPoisson(pm%m%phi,rhs,pm%m%W)				!FFT-based Poisson solver

			!Electric field : -D*phi
			pm%m%E = - Gradient(pm%m%phi,pm%m%dx, pm%ng)
			call control(pm,nk,'E')

			!Force assignment : mat'*E
			do i=1, pm%ns
				call forceAssign(pm%a(i), pm%p(i), pm%m)
			end do
		else
			do i=1,nk-kr*r%mod_r
				call updatePlasma(pm,control,i)
			end do
		end if
	end subroutine

end module
