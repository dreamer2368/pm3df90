module modSpecies

	use constants

	implicit none

	type species
		integer :: np

		real(mp) :: qs, ms, spwt

		real(mp), allocatable :: xp(:,:)
		real(mp), allocatable :: vp(:,:)

		real(mp), allocatable :: Ep(:,:)
	end type

contains

	subroutine buildSpecies(this,np,qs,ms,spwt)
		type(species), intent(out) :: this
		integer, intent(in) :: np
		real(mp), intent(in) :: qs, ms, spwt

		this%np = np
		this%qs = qs
		this%ms = ms
		this%spwt = spwt

		print *, 'Species built up: ms=',ms,', qs=',qs,', specific weight=',spwt
	end subroutine

	subroutine setSpecies(this,np0,xp0,vp0)
		type(species), intent(inout) :: this
		integer, intent(in) :: np0
		real(mp), intent(in), dimension(np0,3) :: xp0, vp0
		if( allocated(this%xp) ) then
			deallocate(this%xp)
		end if
		if( allocated(this%vp) ) then
			deallocate(this%vp)
		end if
		if( allocated(this%Ep) ) then
			deallocate(this%Ep)
		end if

		this%np = np0
		allocate(this%xp(np0,3))
		allocate(this%vp(np0,3))
		allocate(this%Ep(np0,3))
		this%xp = xp0
		this%vp = vp0
		this%Ep = 0.0_mp
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
