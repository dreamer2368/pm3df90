module assignFunctions

	use declaration

	implicit none

contains

	subroutine assignMatrix(this,xp)	!apply BC and create assignment matrix
		type(plasma), intent(inout) :: this
		real(mp), intent(inout) :: xp(this%n,3)
		integer :: i, k
		integer :: g1
		integer :: gl		!left grid point
		integer :: gr		!right grid point
		real(mp) :: fracl(3)		!fraction for left grid point
		real(mp) :: fracr(3)		!fraction for right grid point

		!apply BC
		do i=1,this%n
			!dimension
			do k=1,3
				if( xp(i,k)<0 ) then
					xp(i,k) = xp(i,k) + L
				elseif( xp(i,k)>=L ) then
					xp(i,k) = xp(i,k) - L
				end if
			end do
		end do

		this%frac = 1.0_mp
		!assignment matrix
		do i=1,this%n
			!dimension
			do k=1,3
				g1 = FLOOR(xp(i,k)/dx(k) - 0.5_mp)+1
	!			g1 = FLOOR(xp(i)/dx)+1
				gl = g1
				gr = gl+1
				if (gl<1) then
					gl = gl + this%ng(k)
				elseif (gl>this%ng(k)) then
					gl = gl - this%ng(k)
				end if
				if (gr<1) then
					gr = gr + this%ng(k)
				elseif (gr>this%ng(k)) then
					gr = gr - this%ng(k)
				end if

				fracl(k) = 1.0_mp - ABS(xp(i,k)/dx(k) - g1 + 0.5_mp)
	!			fracl = 1.0_mp - ABS(xp(i)/dx - g1 + 1.0_mp)
				fracr(k) = 1.0_mp - fracl(k)
				this%g(i,:,k) = (/ gl, gr /)
			end do
			this%frac(i,1,:,:) = this%frac(i,1,:,:)*fracl(1)
			this%frac(i,2,:,:) = this%frac(i,2,:,:)*fracr(1)
			this%frac(i,:,1,:) = this%frac(i,:,1,:)*fracl(2)
			this%frac(i,:,2,:) = this%frac(i,:,2,:)*fracr(2)
			this%frac(i,:,:,1) = this%frac(i,:,:,1)*fracl(3)
			this%frac(i,:,:,2) = this%frac(i,:,:,2)*fracr(3)
		end do
	end subroutine

	subroutine chargeAssign(this)
		type(plasma), intent(inout) :: this
		integer :: i, g(2,3)
		real(mp) :: dV

		dV = dx(1)*dx(2)*dx(3)

		this%rho = 0.0_mp
		do i=1,this%n
			g = this%g(i,:,:)
			this%rho( g(:,1), g(:,2), g(:,3) ) = this%rho( g(:,1), g(:,2), g(:,3) ) + this%qs(i)/dV*this%frac(i,:,:,:)
		end do
		this%rho = this%rho + rho_back
	end subroutine

	subroutine forceAssign(this)
		type(plasma), intent(inout) :: this
		integer :: i

		this%Ep = 0.0_mp
		do i=1,this%n
			this%Ep(i,1) = SUM( this%frac(i,:,:,:)*this%E( this%g(i,:,1), this%g(i,:,2), this%g(i,:,3), 1 ) )
			this%Ep(i,2) = SUM( this%frac(i,:,:,:)*this%E( this%g(i,:,1), this%g(i,:,2), this%g(i,:,3), 2 ) )
			this%Ep(i,3) = SUM( this%frac(i,:,:,:)*this%E( this%g(i,:,1), this%g(i,:,2), this%g(i,:,3), 3 ) )
		end do
	end subroutine
!
!	subroutine vpsAssign(this,Es,vps)
!		type(plasma), intent(in) :: this
!		real(mp), intent(in) :: vps(:)
!		real(mp), intent(out) :: Es( size(this%E) )
!		integer :: i
!
!		Es = 0.0_mp
!		do i=1,this%n
!			Es( this%gl(i) ) = Es( this%gl(i) ) + this%fracl(i)*qe/me/dx*vps(i)
!			Es( this%gr(i) ) = Es( this%gr(i) ) + this%fracr(i)*qe/me/dx*vps(i)
!		end do
!	end subroutine
!
!	function xpsUpdate(this,vps,phis,dts,k) result(dxps)
!		type(plasma), intent(in) :: this
!		real(mp), intent(in) :: vps(:)
!		real(mp), intent(in) :: phis(:)
!		real(mp), intent(in) :: dts
!		integer, intent(in) :: k					!current time step
!		real(mp) :: dxps( size(vps) )
!		integer :: i
!
!		dxps = 0.0_mp
!		do i=1,this%n
!			dxps(i) = dts*( -qe/me/dx*vps(i)*( this%Edata(this%gr(i),this%nt+1-k) - this%Edata(this%gl(i),this%nt+1-k) ) )
!			dxps(i) = dxps(i) + dts*qe/dx/eps0*( phis(this%gr(i)) - phis(this%gl(i)) )
!		end do
!	end function
!
end module