module init

	use modPM3D
	use random
	implicit none

contains

	subroutine buildDebye(this,Tf,Ld,Nd,Ng,vT,dir_input,mod_input)
		type(PM3D), intent(inout) :: this
		real(mp), intent(in) :: Tf, Ld, vT
		integer, intent(in) :: Nd,Ng(3)
      character(len=*), intent(in) :: dir_input
		integer, intent(in) :: mod_input
		real(mp) :: Q = 1.0_mp
		real(mp), dimension(Nd*Nd*Nd,3) :: xp0, vp0
		real(mp), dimension(Nd*Nd*Nd) :: spwt
		real(mp), dimension(Ng(1),Ng(2),Ng(3)) :: rho_back
		real(mp) :: xg(Ng(1)), yg(Ng(2)), zg(Ng(3))
		real(mp) :: qs,ms
		real(mp) :: w
		integer :: i,i1,i2,i3

		!Parameter setup
		call buildPM3D(this,Tf,1.0_mp,Ng,1,dt=0.1_mp,L=(/Ld,Ld,Ld/),dir=dir_input,mod_input=mod_input)
		qs = -1.0_mp
		ms = 1.0_mp
		spwt = Ld*Ld*Ld/Nd/Nd/Nd														!rho_0 = 1
		call buildSpecies(this%p(1),qs,ms)

		!Random seed
		call init_random_seed

		!spatial distribution initialize
		do i3 = 1,Nd
			do i2 = 1,Nd
				do i1 = 1,Nd
					xp0(i1+Nd*(i2-1)+Nd*Nd*(i3-1),:) =	&
								(/ (i1-1)*Ld/Nd,(i2-1)*Ld/Nd,(i3-1)*Ld/Nd /)
				end do
			end do
		end do

		!Initial spatial distribution: uniformly-random
!		call RANDOM_NUMBER(xp0)
!		xp0 = xp0*Ld

		!Initial velocity distribution: Maxwell distribution ~ N(0,0.1^2)
		vp0 = randn(Nd*Nd*Nd,3)
		vp0 = vp0*vT

		!Set up species with initial distribution
		call setSpecies(this%p(1),Nd*Nd*Nd,xp0,vp0,spwt)

		!Grid space
		xg =	(/ ( (i-0.5_mp)*Ld/Ng(1),i=1,Ng(1) ) /)
		yg =	(/ ( (i-0.5_mp)*Ld/Ng(2),i=1,Ng(2) ) /)
		zg =	(/ ( (i-0.5_mp)*Ld/Ng(3),i=1,Ng(3) ) /)

		!Set up fixed background charge: Gaussian bump
		w = Ld/20.0_mp
		rho_back = 1.0_mp - Q/Ld/Ld/Ld
		do i3=1,Ng(3)
			do i2=1,Ng(2)
				do i1=1,Ng(1)
					rho_back(i1,i2,i3) = rho_back(i1,i2,i3) + Q/SQRT(8.0_mp*pi*pi*pi)/w/w/w*	&
														EXP( -((xg(i1)-Ld/2.0_mp)**2+(yg(i2)-Ld/2.0_mp)**2+(zg(i3)-Ld/2.0_mp)**2)/2.0_mp/w/w )
				end do
			end do
		end do
		call setMesh(this%m,rho_back)
	end subroutine

	subroutine twostream_initialize(this,Nd,v0)			!generate initial distribution
		type(PM3D), intent(inout) :: this
		integer, intent(in) :: Nd(3)
		real(mp), intent(in) :: v0
		real(mp) :: xp0(PRODUCT(Nd),3), vp0(PRODUCT(Nd),3), qs,ms,spwt(PRODUCT(Nd))
		real(mp), dimension(this%m%ng(1),this%m%ng(2),this%m%ng(3)) :: rho_back
		real(mp) :: L(3),qe,me,mode
		integer :: N,i,i1,i2,i3,pm(PRODUCT(Nd))
		integer :: nseed,clock
		integer, allocatable :: seed(:)

		L = this%L
		N = PRODUCT(Nd)
		mode=1.0_mp

		qe = -this%wp*this%wp/(N/PRODUCT(L))
		qs = qe
		me = -qe
		ms = me
		spwt = 1.0_mp
		rho_back = -qe*N/PRODUCT(L)

		call buildSpecies(this%p(1),qs,ms)
		call setMesh(this%m,rho_back)

		!spatial distribution initialize
!		do i3 = 1,Nd(3)
!			do i2 = 1,Nd(2)
!				do i1 = 1,Nd(1)
!					xp0(i1+Nd(1)*(i2-1)+Nd(1)*Nd(2)*(i3-1),:) =	&
!								(/ (i1-1)*L(1)/Nd(1),(i2-1)*L(2)/Nd(2),(i3-1)*L(3)/Nd(3) /)
!				end do
!			end do
!		end do

		call RANDOM_SEED(size=nseed)
		allocate(seed(nseed))
!		call SYSTEM_CLOCK(COUNT=clock)
!		seed = clock + 127*(/ ( i, i=1,nseed ) /)
		seed = 127*(/ ( i, i=1,nseed ) /)
		call RANDOM_SEED(put=seed)
		deallocate(seed)

		call RANDOM_NUMBER(xp0)
		do i1=1,3
			xp0(:,i1) = xp0(:,i1)*L(i1)
		end do
		xp0(:,1) = xp0(:,1) - this%A0(1)*L(1)/N*SIN( 2.0_mp*pi*xp0(:,1)/L(1)*mode ) + 0.5_mp*L(1)

		!velocity distribution initialize
!		vp0 = vT*randn(N)
		vp0 = 0.0_mp
		pm = (/ ( i, i=1,N ) /)
		pm = 1 - 2*MOD(pm,2)
		vp0(:,1) = vp0(:,1) + pm*v0

		call setSpecies(this%p(1),N,xp0,vp0,spwt)
	end subroutine

end module
