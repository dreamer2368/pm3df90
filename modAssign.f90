module modAssign

	use modSpecies
	use modMesh

	implicit none

	type pmAssign
		integer :: np
		integer :: ng(3)

		integer, allocatable :: g(:,:,:)
		real(mp), allocatable :: frac(:,:,:,:)
	end type

contains

	subroutine buildAssign(this,ng)
		type(pmAssign), intent(out) :: this
		integer, intent(in) :: ng(3)

		this%ng = ng

		allocate(this%g(1,2,3))
		allocate(this%frac(1,2,2,2))
	end subroutine

	subroutine destroyAssign(this)
		type(pmAssign), intent(inout) :: this

		deallocate(this%g)
		deallocate(this%frac)
	end subroutine

	subroutine assignMatrix(this,p,m,xp)	!apply BC and create assignment matrix
		type(pmAssign), intent(inout) :: this
		type(species), intent(inout) :: p
		type(mesh), intent(in) :: m
		real(mp), intent(inout) :: xp(:,:)
		integer :: i, k, np
		integer :: g1
		integer :: gl		!left grid point
		integer :: gr		!right grid point
		real(mp) :: fracl(3)		!fraction for left grid point
		real(mp) :: fracr(3)		!fraction for right grid point

		np = size(xp,1)
		this%np = np
		deallocate(this%g)
		deallocate(this%frac)
		allocate(this%g(np,2,3))
		allocate(this%frac(np,2,2,2))
		!apply BC
		do i=1,np
			!dimension
			do k=1,3
				if( xp(i,k)<0 ) then
					xp(i,k) = xp(i,k) + m%L(k)
				elseif( xp(i,k)>=m%L(k) ) then
					xp(i,k) = xp(i,k) - m%L(k)
				end if
				if( xp(i,k)< -m%dx(k) .or. xp(i,k)>=m%L(k) + m%dx(k) ) then
					print *, 'FAULT!! particles outside the boundary : check if particles are too fast'
					print *, 'xp(',i,', ',k,') =',xp(i,k)
					stop
				end if
			end do
		end do

		this%frac = 1.0_mp
		!assignment matrix
		do i=1,this%np
			!dimension
			do k=1,3
				g1 = FLOOR(xp(i,k)/m%dx(k) - 0.5_mp)+1
	!			g1 = FLOOR(xp(i)/m%dx)+1
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

				fracl(k) = 1.0_mp - ABS(xp(i,k)/m%dx(k) - g1 + 0.5_mp)
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

	subroutine chargeAssign(this,p,m)
		type(pmAssign), intent(inout) :: this(:)
		type(species), intent(inout) :: p(:)
		type(mesh), intent(inout) :: m
		integer :: i,ip, g(2,3)
		real(mp) :: dV

		dV = PRODUCT(m%dx)

		m%rho = 0.0_mp
		do ip=1,size(p)
			do i=1,p(ip)%np
				g = this(ip)%g(i,:,:)
				m%rho( g(:,1), g(:,2), g(:,3) ) = m%rho( g(:,1), g(:,2), g(:,3) ) + p(ip)%qs*p(ip)%spwt/dV*this(ip)%frac(i,:,:,:)
			end do
		end do
		m%rho = m%rho + m%rho_back
	end subroutine

	subroutine forceAssign(this,p,m)
		type(pmAssign), intent(inout) :: this
		type(species), intent(inout) :: p
		type(mesh), intent(in) :: m
		integer :: i,j

		p%Ep = 0.0_mp
		do i=1,this%np
			p%Ep(i,1) = SUM( this%frac(i,:,:,:)*m%E( this%g(i,:,1), this%g(i,:,2), this%g(i,:,3), 1 ) )
			p%Ep(i,2) = SUM( this%frac(i,:,:,:)*m%E( this%g(i,:,1), this%g(i,:,2), this%g(i,:,3), 2 ) )
			p%Ep(i,3) = SUM( this%frac(i,:,:,:)*m%E( this%g(i,:,1), this%g(i,:,2), this%g(i,:,3), 3 ) )
		end do
	end subroutine

end module
