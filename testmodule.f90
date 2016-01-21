module testmodule

	use init
	use timeStep

	implicit none

contains

	subroutine verify_assignment(this)
		type(plasma), intent(inout) :: this
		real(mp) :: testqs(2), testms(2), testxp(2,3), testvp(2,3)
		integer :: testN
		real(mp), allocatable :: rhs(:,:,:)
		integer :: i,j, g(2,3)

		testN = 2

		call setup
		testqs(1) = -wp*wp/(2/L/L/L)
		testqs(2) = -testqs(1)
		rho_back = 0.0_mp
		testms = wp*wp/(2/L/L/L)
		testxp(1,:) = (/ 0.4_mp*L, 0.6_mp*L, 0.7_mp*L /)
		testxp(2,:) = (/ L-dx(1), 0.0_mp*L, 0.0_mp*L /)
		testvp = 0.0_mp

		call buildPlasma(this,Nt,Ni,2,Ng,testxp,testvp,testqs,testms)
		call assignMatrix(this,this%xp)

		do i=1,testN
			print *, 'particle ',i,' (column-major order)'
			do j=1,8
!				print *, 2-MOD(j,2),MOD((j-1)/2,2)+1,MOD((j-1)/4,2)+1
				print *, 'g(',this%g(i,2-MOD(j,2),1),',',this%g(i,MOD((j-1)/2,2)+1,2),',',this%g(i,MOD((j-1)/4,2)+1,3),')=',	&
				this%frac(i,2-MOD(j,2),MOD((j-1)/2,2)+1,MOD((j-1)/4,2)+1)
			end do
			print *, 'sum of fraction = ', SUM(this%frac(i,:,:,:))
		end do

		call chargeAssign(this)

		print *, 'Check on the grid'

		do i=1, testN
			print *, 'Particle ',i
			g = this%g(i,:,:)
			print *, 'Xg=',g(:,1)
			print *, 'Yg=',g(:,2)
			print *, 'Zg=',g(:,3)
			print *, 'Particle fraction(column-major order)'
			print *, this%rho( g(:,1), g(:,2), g(:,3) )*dx(1)*dx(2)*dx(3)/this%qs(i)
		end do

		allocate(rhs(this%ng(1),this%ng(2),this%ng(3)))
		rhs = -this%rho/eps0
		call FFTPoisson(this%phi,rhs,W)
		this%E = - Gradient(this%phi,dx,Ng)

		open(unit=301,file='data/testphi.bin',status='replace',form='unformatted',access='stream')
		write(301) this%phi
		close(301)

		open(unit=301,file='data/testE.bin',status='replace',form='unformatted',access='stream')
		write(301) this%E
		close(301)

		call destroyPlasma(this)
	end subroutine

end module