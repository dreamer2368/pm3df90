module modSpecies

	use constants

	implicit none

	type species
		integer :: np

		real(mp) :: qs, ms
		real(mp), allocatable :: spwt(:)

		real(mp), allocatable :: xp(:,:)
		real(mp), allocatable :: vp(:,:)

		real(mp), allocatable :: Ep(:,:)
	end type

contains

	subroutine buildSpecies(this,qs,ms)
		type(species), intent(out) :: this
		real(mp), intent(in) :: qs, ms

		this%qs = qs
		this%ms = ms

		print *, 'Species built up: ms=',ms,', qs=',qs
	end subroutine

	subroutine setSpecies(this,np0,xp0,vp0,spwt0)
		type(species), intent(inout) :: this
		integer, intent(in) :: np0
		real(mp), intent(in), dimension(np0,3) :: xp0, vp0
		real(mp), intent(in), dimension(np0) :: spwt0
		real(mp) :: temp
		if( allocated(this%xp) ) then
			deallocate(this%xp)
		end if
		if( allocated(this%vp) ) then
			deallocate(this%vp)
		end if
		if( allocated(this%Ep) ) then
			deallocate(this%Ep)
		end if
		if( allocated(this%spwt) ) deallocate(this%spwt)

		this%np = np0
		allocate(this%xp(np0,3))
		allocate(this%vp(np0,3))
		allocate(this%Ep(np0,3))
		this%xp = xp0
		this%vp = vp0
		this%Ep = 0.0_mp

		allocate(this%spwt(np0))
		this%spwt = spwt0
	end subroutine

	subroutine destroySpecies(this)
		type(species), intent(inout) :: this

		if( allocated(this%xp) ) then
			deallocate(this%xp)
		end if
		if( allocated(this%vp) ) then
			deallocate(this%vp)
		end if
		if( allocated(this%Ep) ) then
			deallocate(this%Ep)
		end if
		if( allocated(this%spwt) ) deallocate(this%spwt)
	end subroutine

	subroutine move(this,dt)
		type(species), intent(inout) :: this
		real(mp), intent(in) :: dt

		this%xp = this%xp + dt*this%vp
	end subroutine

	subroutine accel(this,dt)
		type(species), intent(inout) :: this
		real(mp), intent(in) :: dt

		this%vp = this%vp + dt*this%qs/this%ms*this%Ep
	end subroutine

end module
