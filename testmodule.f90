module testmodule

	use init
	use timeStep

	implicit none

contains

	subroutine test_FFTPoisson_adj(N,Nf)
		integer, intent(in) :: N(3)										!!grid number
		integer, intent(in) :: Nf										!!number of fourier modes

		real(mp) :: dx(3)												!!grid size
		real(mp) :: xg(N(1)), yg(N(2)), zg(N(3))
		real(mp) :: eps0 = 1.0_mp
		real(mp) :: L(3) = (/ 1.0_mp, 1.0_mp, 1.0_mp /)					!cubic size
		real(mp), dimension(N(1),N(2),N(3)) :: rho, rhs, weight, rhos
		real(mp), dimension(N(1),N(2),N(3),3) :: E,Es
		real(mp) :: rhok(N(1),N(2),N(3),Nf)

		real(mp) :: A(Nf), dA(Nf)										!coefficients of fourier modes
		real(mp) :: J0,J1,dJdA(Nf)

		complex(mp) :: W(N(1),N(2),N(3))
		integer :: i,j,k,m

		dx = L/N
		xg = (/ (i*L(1)/N(1), i=0,N(1)-1) /)
		yg = (/ (i*L(2)/N(2), i=0,N(2)-1) /)
		zg = (/ (i*L(3)/N(3), i=0,N(3)-1) /)

		A = 1.0_mp
		rho = 0.0_mp
		rhok = 0.0_mp
		do i=1,Nf
			do k=1,N(3)
				do j=1,N(2)
					do m=1,N(1)
						rhok(m,j,k,i) = SIN(2.0_mp*pi*xg(m)/L(1)*i)	&
									*SIN(2.0_mp*pi*yg(j)/L(2)*i)	&
									*SIN(2.0_mp*pi*zg(k)/L(3)*i)
					end do
				end do
			end do
			rho = rho + A(i)*rhok(:,:,:,i)
		end do

		weight = 1.0_mp

		call FFTPoisson_setup(N,dx,W)

		rhs = -rho/eps0

		call FFTEfield(E,rhs,W)

		J0 = SUM( PRODUCT(dx)*weight*(E(:,:,:,1)**2 + E(:,:,:,2)**2 + E(:,:,:,3)**2) )
		print *, 'J0 = ', J0

		do i=1,3
			Es(:,:,:,i) = -2.0_mp*PRODUCT(dx)*weight*E(:,:,:,i)
		end do

		call FFTAdj(Es,rhos,W)
		rhos = rhos/eps0

		do i=1,Nf
			dJdA(i) = SUM( rhos*rhok(:,:,:,i) )
		end do
		print *, 'dJdA = ', dJdA

		dA = 0.0_mp
		dA(2) = (0.1_mp)**9
		A = A+dA
		rho = 0.0_mp
		do i=1,Nf
			rho = rho + A(i)*rhok(:,:,:,i)
		end do

		rhs = -rho/eps0

		call FFTEfield(E,rhs,W)

		J1 = SUM( PRODUCT(dx)*weight*(E(:,:,:,1)**2 + E(:,:,:,2)**2 + E(:,:,:,3)**2) )
		print *, 'J1 = ', J1
		print *, 'actual dJdA = ', (J1-J0)/dA(2)
	end subroutine

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