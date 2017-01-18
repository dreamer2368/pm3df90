module modControl

	use modAdj

	implicit none

contains

!==============Default=====================================================

	subroutine Null_input(this,k,str)
		type(PM3D), intent(inout) :: this
		integer, intent(in) :: k
		character(len=*), intent(in) :: str
	end subroutine

	subroutine Null_Dinput(adj,pm,k,str)
		type(adjoint), intent(inout) :: adj
		type(PM3D), intent(in) :: pm
		integer, intent(in) :: k
		character(len=*), intent(in) :: str
	end subroutine

	subroutine Null_dJdA(adj,pm,k,str,grad)
		type(adjoint), intent(in) :: adj
		type(PM3D), intent(in) :: pm
		integer, intent(in) :: k
		character(len=*), intent(in) :: str
		real(mp), intent(inout) :: grad(:)
	end subroutine

!==============Increase thermal velocity===================================

	subroutine V_Th(this,k,str)
		type(PM3D), intent(inout) :: this
		integer, intent(in) :: k
		character(len=*), intent(in) :: str

		SELECT CASE (str)
			CASE('vp')
				if(k==this%ni) then
					this%p(1)%vp = this%p(1)%vp + this%A0(1)*this%p(1)%vp
				end if
		END SELECT
	end subroutine

	subroutine dV_Th(adj,pm,k,str)
		type(adjoint), intent(inout) :: adj
		type(PM3D), intent(in) :: pm
		integer, intent(in) :: k
		character(len=*), intent(in) :: str

		select case (str)
			case ('vp')
				if( k .eq. pm%ni-1 ) then
					adj%dp(1)%vp = adj%p(1)%vp*pm%A0(1)/(-pm%dt)
				end if
		end select
	end subroutine

	subroutine dV_Th_dB(adj,pm,k,str,grad)
		type(adjoint), intent(in) :: adj
		type(PM3D), intent(in) :: pm
		integer, intent(in) :: k
		character(len=*), intent(in) :: str
		real(mp), intent(inout) :: grad(:)

		if( str.eq.'before' .and. k.eq.pm%ni-1 ) then
			grad(1) = -sum( adj%p(1)%vp*pm%p(1)%vp )/pm%dt
		end if

!		dJdA = SUM( - pm%r%vpdata(:,:,pm%ni-1)*pm%r%vpsdata(:,:,pm%ni)/pm%dt )
	end subroutine

!==============Wave perturbation on position===============================

	subroutine IC_wave(this,k,str)
		type(PM3D), intent(inout) :: this
		integer, intent(in) :: k
		character(len=*), intent(in) :: str
		integer :: mode = 1

		SELECT CASE (str)
			CASE('xp')
				if(k==this%ni) then
					this%p(1)%xp(:,1) = this%p(1)%xp(:,1)	&
											+ this%dt*this%A0(1)*this%L(1)/this%p(1)%np	&
												*SIN(4.0_mp*pi*this%p(1)%xp(:,1)/this%L(1)*mode)		!xp_(Ni-1) : k=Ni, xp_Ni : k=(Ni+1)
				end if
		END SELECT
	end subroutine

	subroutine dIC_wave(adj,pm,k,str)
		type(adjoint), intent(inout) :: adj
		type(PM3D), intent(in) :: pm
		integer, intent(in) :: k
		character(len=*), intent(in) :: str

		select case (str)
			case ('xp')
				if( k .eq. pm%ni-1 ) then
!					adj%dp(1)%xp = adj%dp(1)%xp - adj%p(1)%xp*(pm%B0*pm%L(1)/pm%p(1)%np*4.0_mp*pi*mode/pm%L(1))*COS(4.0_mp*pi*mode*pm%p(1)%xpdata(:,:,k)/pm%L(1))
				end if
		end select
	end subroutine

	subroutine dIC_wave_dB(adj,pm,dJdA)
		type(adjoint), intent(in) :: adj
		type(PM3D), intent(in) :: pm
		real(mp), intent(out) :: dJdA

!		dJdA = SUM( - pm%L(1)/pm%n*SIN( 4.0_mp*pi*mode*pm%r%xpdata(:,1,pm%ni-1)/pm%L(1) )*pm%r%xpsdata(:,1,pm%ni) )
	end subroutine

!================Two Particle Orbit radius===========================

	subroutine Orbit_radius(pm,k,str)
		type(PM3D), intent(inout) :: pm
		integer, intent(in) :: k
		character(len=*), intent(in) :: str

		select case (str)
			case('xp')
				if( k.eq.2 ) then
!					pm%p%xp(1,1) = pm%p%xp(1,1) - pm%B0*0.5_mp*(pm%p%xp(2,1)-pm%p%xp(1,1))
!					pm%p%xp(2,1) = pm%p%xp(2,1) + pm%B0*0.5_mp*(pm%p%xp(2,1)-pm%p%xp(1,1))
				end if
		end select
	end subroutine

	subroutine dOrbit(adj,pm,k,str)
		type(adjoint), intent(inout) :: adj
		type(PM3D), intent(in) :: pm
		integer, intent(in) :: k
		character(len=*), intent(in) :: str

	end subroutine

	subroutine dOrbit_dB(adj,pm,dJdA)
		type(adjoint), intent(in) :: adj
		type(PM3D), intent(in) :: pm
		real(mp), intent(out) :: dJdA
		real(mp) :: r0

!		r0 = pm%r%xpdata(2,1,1) - pm%r%xpdata(1,1,1)
!		dJdA = - ( -0.5_mp*r0*pm%r%xpsdata(1,1,2) + 0.5_mp*r0*pm%r%xpsdata(2,1,2) )/pm%dt
	end subroutine

!=================Testmodule =======================================

	subroutine testControl(pm,k,str)
		type(PM3D), intent(inout) :: pm
		integer, intent(in) :: k
		character(len=*), intent(in) :: str

		select case (str)
			case('xp')
				if( k.eq.1 ) then
					pm%p(1)%xp(int(pm%A0(1)),int(pm%A0(2))) = pm%p(1)%xp(int(pm%A0(1)),int(pm%A0(2)))*(1.0_mp+pm%A0(3))
				end if
		end select
	end subroutine

	subroutine testdJdA(adj,pm,k,str,grad)
		type(PM3D), intent(in) :: pm
		type(adjoint), intent(in) :: adj
		integer, intent(in) :: k
		character(len=*), intent(in) :: str
		real(mp), intent(inout) :: grad(:)

		if( k.eq.1 ) then
			grad(1:3) = -adj%p(1)%xp(1,:)/pm%dt
			grad(4:6) = -adj%p(1)%xp(2,:)/pm%dt
			grad(7:9) = -adj%p(1)%vp(1,:)/pm%dt - adj%p(1)%xp(1,:)
			grad(10:12) = -adj%p(1)%vp(2,:)/pm%dt - adj%p(1)%xp(2,:)
		end if
	end subroutine

end module
