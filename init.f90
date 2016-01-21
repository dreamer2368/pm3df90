module init

	use declaration
	use MatrixSolver
	use random
	implicit none

contains

	subroutine setup()
		!reference perturbation amplitude
		A = A0
		B = B0
		!time step parameter
		dt = 0.2_mp
		Nt = CEILING(T/dt)
		dt = T/Nt
		Ni = FLOOR(Ti/dt) + 1
		print *, 'Ni=',Ni,', Nt=',Nt,', dt=',dt

		!Grid setup
		dx = L/Ng
		!FFTSolver setup
		call FFTPoisson_setup(Ng,dx,W)

		xg = (/ ( i*L/Ng(1), i=0,Ng(1)-1 ) /) + 0.5_mp*dx(1)

		open(unit=301,file='data/L.out',status='replace')
		write(301,*) L
		close(301)
	end subroutine

	function randn(N) result(x)
		integer, intent(in) :: N
		real(mp) :: x(N)
		integer :: i

		do i = 1,N
			x(i) = SQRT(-2.0_mp*LOG(RAND()))*COS(2.0_mp*pi*RAND())
		end do
	end function

	subroutine particle_initialize()
		integer :: i1,i2,i3
		real(mp) :: start, finish

		!species charge and mass assign - single species with electrons right now
		qs = qe
		ms = me
		rho_back = -qe*N/L/L/L						!For uniform background charge

		!spatial distribution initialize
		do i3 = 1,Nd(3)
			do i2 = 1,Nd(2)
				do i1 = 1,Nd(1)
					xp0(i1+Nd(1)*(i2-1)+Nd(1)*Nd(2)*(i3-1),:) = (/ (i1-1)*L/Nd(1),(i2-1)*L/Nd(2),(i3-1)*L/Nd(3) /)
				end do
			end do
		end do
		xp0(:,1) = xp0(:,1) - A*L/N*SIN( 2.0_mp*pi*xp0(:,1)/L*mode ) + 0.5_mp*L

		!velocity distribution initialize
!		vp0 = vT*randn(N)
		vp0 = 0.0_mp
		pm = (/ ( i, i=1,N ) /)
		pm = 1 - 2*MOD(pm,2)
		vp0(:,1) = vp0(:,1) + pm*v0

		open(unit=301,file='data/xp0.bin',status='replace',form='unformatted',access='stream')
		write(301) xp0
		close(301)
	end subroutine

end module