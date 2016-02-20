module modSource

	use modAdj

	implicit none

contains

!==============Wave perturbation on position===============================

	subroutine IC_wave(this,k,str)
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

	subroutine dIC_wave(adj,pm,k,str)
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

	subroutine dIC_wave_dB(adj,pm,dJdA)
		type(adjoint), intent(in) :: adj
		type(PM3D), intent(in) :: pm
		real(mp), intent(out) :: dJdA

		dJdA = SUM( - pm%L(1)/pm%n*SIN( 4.0_mp*pi*mode*pm%r%xpdata(:,1,pm%ni-1)/pm%L(1) )*pm%r%xpsdata(:,1,pm%ni) )
	end subroutine

end module