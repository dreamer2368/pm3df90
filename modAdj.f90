module modAdj

	use modPM3D

	implicit none

	type adjoint
		integer :: nt, ni, n, ng(3)
		real(mp) :: J0, J1
		real(mp), allocatable :: xps(:,:), vps(:,:), Eps(:,:), Es(:,:,:,:), rhos(:,:,:)
	end type

contains

	subroutine buildAdjoint(this,pm)
		type(adjoint), intent(out) :: this
		type(PM3D), intent(in) :: pm

		this%ng = pm%ng
		this%n = pm%n
		this%nt = pm%nt
		this%ni = pm%ni

		allocate(this%xps(this%n,3))
		allocate(this%vps(this%n,3))
		allocate(this%Eps(this%n,3))
		allocate(this%Es(this%ng(1),this%ng(2),this%ng(3),3))
		allocate(this%rhos(this%ng(1),this%ng(2),this%ng(3)))
	end subroutine

	subroutine destroyAdjoint(this)
		type(adjoint), intent(inout) :: this

		deallocate(this%xps)
		deallocate(this%vps)
		deallocate(this%Eps)
		deallocate(this%Es)
		deallocate(this%rhos)
	end subroutine

	subroutine QoI(this,pm,i)
		type(adjoint), intent(inout) :: this
		type(PM3D), intent(in) :: pm
		integer, intent(in), optional :: i
		integer :: input
		if( present(i) ) then
			input = i
		else
			input = 0
		end if

		if( input.eq.0 ) then
			this%J0 = SUM( (pm%p%vp(:,1)**2 + pm%p%vp(:,2)**2 + pm%p%vp(:,3)**2) )
			print *, 'J0 = ', this%J0
		elseif( input.eq.1 ) then
			this%J1 = SUM( (pm%p%vp(:,1)**2 + pm%p%vp(:,2)**2 + pm%p%vp(:,3)**2) )
			print *, 'J1 = ', this%J1
		end if
	end subroutine

!subroutine QOI(this,J)
!type(PM3D), intent(in) :: this
!real(mp), intent(out) :: J
!integer :: N,Nt,Ni
!N = this%n
!Nt = this%nt
!Ni = this%ni
!
!J = 1.0_mp/N/(Nt-Ni)*SUM(this%r%vpdata(:,:,Ni+1:Nt)**2)		!omitted 1/N/T for the sake of machine precision
!end subroutine

end module
