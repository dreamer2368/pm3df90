module modQoI

	use modAdj

	implicit none

contains

!=============Mean Kinetic Energy====================

	subroutine MKE(pm,k,J)
		type(PM3D), intent(in) :: pm
		integer, intent(in) :: k
		real(mp), intent(inout) :: J
		integer :: input

		if( k.ge.pm%ni ) then
			J = J + 1.0_mp/pm%p(1)%np/(pm%nt-pm%ni)*SUM( pm%p(1)%vp**2 )
		end if
	end subroutine

	subroutine dMKE(adj,pm,k)
		type(adjoint), intent(inout) :: adj
		type(PM3D), intent(in) :: pm
		integer, intent(in) :: k

		if( k.ge.pm%ni ) then
			adj%dp(1)%vp = adj%dp(1)%vp + 2.0_mp/pm%p(1)%np/(pm%nt-pm%ni)*pm%p(1)%vp
		end if
	end subroutine

!============Mean Efield Energy======================

	subroutine MPE(pm,k,J)
		type(PM3D), intent(in) :: pm
		integer, intent(in) :: k
		real(mp), intent(inout) :: J

		if( k.ge.pm%ni ) then
			J = J + 1.0_mp/PRODUCT(pm%ng)/(pm%nt-pm%ni)*SUM( pm%m%E**2 )
		end if
	end subroutine

	subroutine dMPE(adj,pm,k)
		type(adjoint), intent(inout) :: adj
		type(PM3D), intent(in) :: pm
		integer, intent(in) :: k

		if( k .ge. pm%ni ) then
			adj%dm%E = - 2.0_mp/PRODUCT(pm%ng)/(pm%nt-pm%ni)*pm%m%E
		end if
	end subroutine

!============Testmodule : Single step Total Kinetic Energy

	subroutine TestQoI_singleStep(pm,k,J)
		type(PM3D), intent(in) :: pm
        integer, intent(in) :: k
        real(mp), intent(inout) :: J

		!For testmodule
		if( k.eq.pm%nt ) then
			J = SUM( (pm%p(1)%vp(:,1)**2 + pm%p(1)%vp(:,2)**2 + pm%p(1)%vp(:,3)**2) )
		end if
	end subroutine

	subroutine dTestQoI_singleStep(adj,pm,k)
		type(adjoint), intent(inout) :: adj
		type(PM3D), intent(in) :: pm
		integer, intent(in) :: k

		if( k.eq.pm%nt ) then
			adj%p(1)%vp = - 2.0_mp*pm%p(1)%vp*pm%dt
		end if
	end subroutine

end module
