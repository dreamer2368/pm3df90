module modQoI

	use modAdj

	implicit none

contains

!=============Mean Kinetic Energy====================

	subroutine MKE(this,pm,i)
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
			this%J0 = 1.0_mp/pm%n/(pm%nt-pm%ni)*SUM( pm%r%vpdata(:,:,pm%ni+1:pm%nt)**2 )
			print *, 'J0 = ', this%J0
		elseif( input.eq.1 ) then
			this%J1 = 1.0_mp/pm%n/(pm%nt-pm%ni)*SUM( pm%r%vpdata(:,:,pm%ni+1:pm%nt)**2 )
			print *, 'J1 = ', this%J1
		end if
	end subroutine

	subroutine dMKE(adj,pm,k)
		type(adjoint), intent(inout) :: adj
		type(PM3D), intent(in) :: pm
		integer, intent(in) :: k

		if( k >= pm%ni+1 ) then
			adj%dvps = adj%dvps + 2.0_mp/pm%n/(pm%nt-pm%ni)*pm%r%vpdata(:,:,k)
		end if
	end subroutine

!============Testmodule : Single step Total Kinetic Energy

	subroutine TestQoI_singleStep(this,pm,i)
		type(adjoint), intent(inout) :: this
		type(PM3D), intent(in) :: pm
		integer, intent(in), optional :: i
		integer :: input
		if( present(i) ) then
			input = i
		else
			input = 0
		end if

		!For testmodule
		if( input.eq.0 ) then
			this%J0 = SUM( (pm%p%vp(:,1)**2 + pm%p%vp(:,2)**2 + pm%p%vp(:,3)**2) )
			print *, 'J0 = ', this%J0
		elseif( input.eq.1 ) then
			this%J1 = SUM( (pm%p%vp(:,1)**2 + pm%p%vp(:,2)**2 + pm%p%vp(:,3)**2) )
			print *, 'J1 = ', this%J1
		end if
	end subroutine

end module