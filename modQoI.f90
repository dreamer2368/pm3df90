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
		type(PM3D), intent(inout) :: pm
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

!=============E-field Probe : measure E-field at positions

	subroutine E_probe(pm,i,J)
		type(PM3D), intent(in) :: pm
		integer, intent(in) :: i
		real(mp), intent(inout) :: J
		type(species) :: probe
		type(pmAssign) :: a
		real(mp), dimension(3,3) :: xp0, vp0
		integer :: k
		!first probe position
		xp0(1,:) = (/pm%L(1)/8.0_mp,0.0_mp,0.0_mp/)
		!second probe position
		xp0(2,:) = (/0.0_mp,pm%L(1)/8.0_mp,0.0_mp/)
		!third probe position
		xp0(3,:) = (/0.0_mp,0.0_mp,pm%L(1)/8.0_mp/)
		!set zero velocity
		vp0 = 0.0_mp

		call buildSpecies(probe,3,1.0_mp,1.0_mp,1.0_mp)
		call setSpecies(probe,3,xp0,vp0)

		!locate probes
		call buildAssign(a,pm%m%ng)
		call assignMatrix(a,probe,pm%m,probe%xp)

		!measure E-field
		call forceAssign(a,probe,pm%m)

		do k=1,3
			print *, k,'-th probe at (',probe%xp(k,1),',',probe%xp(k,2),',',probe%xp(k,3),')'
			print *, '            = E(',probe%Ep(k,1),',',probe%Ep(k,2),',',probe%Ep(k,3),')'
		end do

		call destroySpecies(probe)
		call destroyAssign(a)
	end subroutine

!===================Debye_length=============================

   subroutine Debye_length(pm,k,J)
      type(PM3D), intent(in) :: pm
      integer, intent(in) :: k
      real(mp), intent(inout) :: J
      real(mp), dimension(pm%ng(1),pm%ng(2),pm%ng(3)) :: Xg2      !squared distance from center
      integer :: i1,i2,i3
      
      do i3 = 1,pm%ng(3)
         do i2 = 1,pm%ng(2)
            do i1 = 1,pm%ng(1)
               Xg2(i1,i2,i3) =  SUM( ((/ (i1-0.5_mp),(i2-0.5_mp),(i3-0.5_mp) /)*pm%L/pm%ng - 0.5_mp*pm%L)**2 )
            end do
         end do
      end do

		J = J + PRODUCT(pm%m%dx)/pm%nt*SUM( Xg2*pm%m%rho )
   end subroutine

end module
