module declaration

	use modPlasma

	implicit none

	! parameters setup
	integer, parameter :: Nd(3) = (/ 2**0, 2**7, 2**7 /)
	integer, parameter :: N = Nd(1)*Nd(2)*Nd(3)
	real(mp), parameter :: L = 2*pi/( sqrt(3.0_mp)/2.0_mp/sqrt(2.0_mp)/0.2_mp )
	real(mp), parameter :: eps0 = 1.0_mp
	real(mp), parameter :: wp = 1.0_mp
	real(mp), parameter :: qe = -wp*wp/(N/L/L/L)			!electron charge
	real(mp), parameter :: me = -qe							!electron mass
	real(mp) :: rho_back = -qe*N/L/L/L						!For uniform background charge

	real(mp) :: Ti = 40.0_mp
	real(mp) :: Tp = 2.0_mp*pi/wp
!	real(mp) :: T = 40.0_mp
	real(mp) :: T = 40.0_mp + 10.0_mp*(2.0_mp*pi/wp)

	!for loop index
	integer :: i

	! time step parameter
	real(mp) :: dt
	integer :: Nt, Ni

	!initial spatial distribution
	real(mp) :: xp0(N,3)
	real(mp), parameter :: A0 = 1.0_mp, B0 = 0.0_mp
	real(mp) :: A = A0
	real(mp) :: B = B0
	integer, parameter :: mode = 1

	!initial velocity distribution
	real(mp), parameter :: vT = 0.0_mp
	real(mp), parameter :: v0 = 0.0_mp
	real(mp) :: vp0(N,3)
	integer :: pm(N)

	!grid and operators setup
	integer, parameter :: Ng(3) = (/ 64, 64, 64 /)
	real(mp) :: dx(3)
	real(mp) :: xg(Ng(1))
	!FFT_Poisson
	complex(mp) :: W(Ng(1)/2+1,Ng(2),Ng(3))
!	real(mp) :: rhs(Ng-1)

	!time-stepping particle & field variable are composed of type fluid
	type(plasma) :: langmuir
	real(mp) :: qs(N), ms(N)			!charge and mass of species

	!error convergence variable
	real(mp) :: fDA(30)
	real(mp) :: ek(30)
	real(mp) :: e(9,30)

	!check grid passing
	integer :: pass
	integer :: coincide

!	character(32) :: filename ! You can make this longer than 32 characters, if needed

end module
