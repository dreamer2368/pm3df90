module modAdj

	use modPM3D

	implicit none

	type adjoint
		integer :: nt, ni, n, ng(3)
		real(mp) :: J0, J1
		real(mp), allocatable :: xps(:,:), vps(:,:), Eps(:,:), Es(:,:,:,:), rhos(:,:,:)
	end type

contains

	subroutine buildAdjoint(this,pm)
		type(adjoint), intent(out) :: this
		type(PM3D), intent(in) :: pm

		this%ng = pm%ng
		this%n = pm%n
		this%nt = pm%nt
		this%ni = pm%ni

		allocate(this%xps(this%n,3))
		allocate(this%vps(this%n,3))
		allocate(this%Eps(this%n,3))
		allocate(this%Es(this%ng(1),this%ng(2),this%ng(3),3))
		allocate(this%rhos(this%ng(1),this%ng(2),this%ng(3)))
	end subroutine

	subroutine destroyAdjoint(this)
		type(adjoint), intent(inout) :: this

		deallocate(this%xps)
		deallocate(this%vps)
		deallocate(this%Eps)
		deallocate(this%Es)
		deallocate(this%rhos)
	end subroutine

	subroutine QoI(this,pm,i)
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
		!For testmodule
!		if( input.eq.0 ) then
!			this%J0 = SUM( (pm%p%vp(:,1)**2 + pm%p%vp(:,2)**2 + pm%p%vp(:,3)**2) )
!			print *, 'J0 = ', this%J0
!		elseif( input.eq.1 ) then
!			this%J1 = SUM( (pm%p%vp(:,1)**2 + pm%p%vp(:,2)**2 + pm%p%vp(:,3)**2) )
!			print *, 'J1 = ', this%J1
!		end if
	end subroutine

	subroutine backward_sweep(adj,pm, dJdA)
		type(adjoint), intent(inout) :: adj
		type(PM3D), intent(inout) :: pm
		real(mp), intent(out) :: dJdA
		real(mp), dimension(pm%n,3) :: dvps, dxps1, dxps2, dxps3
		integer :: k, nk

		adj%xps = 0.0_mp
		adj%vps = 0.0_mp
		do k=1,pm%nt
			nk = pm%nt+1-k

			!======= dv_p =============
			dvps = merge(	2.0_mp/pm%n/(pm%nt-pm%ni)*pm%r%vpdata(:,:,nk),	&
							0.0_mp,	&
							nk >= pm%ni+1	)
			call Adj_accel(adj,pm,dvps)

			!Check when adjoint reach to the initial step
!			if( k .eq. pm%nt ) then
!				adj%vps = 0.5_mp*adj%vps
!			end if

			!======= dE_p =============
			adj%Eps(:,1) = pm%p%qs/pm%p%ms*adj%vps(:,1)
			adj%Eps(:,2) = pm%p%qs/pm%p%ms*adj%vps(:,2)
			adj%Eps(:,3) = pm%p%qs/pm%p%ms*adj%vps(:,3)

			call assignMatrix(pm%a,pm%p,pm%m,pm%r%xpdata(:,:,nk))

			!======= dE_g =============
			call Adj_forceAssign_E(pm%a,adj%Eps,adj%Es)

			!======= dPhi_g, dRho_g =============
			call FFTAdj(adj%Es,adj%rhos,pm%m%W,pm%m%dx)
			adj%rhos = -adj%rhos/pm%eps0

			!======= dx_p =============
			call Adj_chargeAssign(pm%a,pm%p,pm%m,adj%rhos,dxps1)
			call Adj_forceAssign_xp(pm%a,pm%m,pm%r%Edata(:,:,:,:,nk),adj%Eps,dxps2)
			if( nk .eq. pm%ni-1 ) then
				dxps3 = -adj%xps*(pm%B0*pm%L(1)/pm%n*4.0_mp*pi*mode/pm%L(1))*COS(4.0_mp*pi*mode*pm%r%xpdata(:,:,nk)/pm%L(1))
			else
				dxps3 = 0.0_mp
			end if
			call Adj_move(adj,pm,dxps1+dxps2+dxps3)

			call recordAdjoint(pm%r,nk,adj%xps)									!record for nk=1~Nt : xps_nk and vp_(nk+1/2)

			if( nk<pm%ni ) then
				exit
			end if
		end do

		dJdA = SUM( - pm%L(1)/pm%n*SIN( 4.0_mp*pi*mode*pm%r%xpdata(:,1,pm%ni-1)/pm%L(1) )*pm%r%xpsdata(:,1,pm%ni) )
	end subroutine

	subroutine Adj_chargeAssign(this,p,m,rhos,xps)
		type(pmAssign), intent(inout) :: this
		type(plasma), intent(in) :: p
		type(mesh), intent(in) :: m
		real(mp), intent(in) :: rhos(this%ng(1),this%ng(2),this%ng(3))
		real(mp), intent(out) :: xps(this%n,3)
		real(mp) :: dV
		integer :: i

		dV = PRODUCT(m%dx)
		xps = 0.0_mp
		do i=1,this%n
			xps(i,1) = 1.0_mp/m%dx(1)*SUM( SUM( this%frac(i,:,:,:), 1)*rhos( this%g(i,2,1), this%g(i,:,2), this%g(i,:,3) )	&
											- SUM( this%frac(i,:,:,:), 1)*rhos( this%g(i,1,1), this%g(i,:,2), this%g(i,:,3) ) )
			xps(i,2) = 1.0_mp/m%dx(2)*SUM( SUM( this%frac(i,:,:,:), 2)*rhos( this%g(i,:,1), this%g(i,2,2), this%g(i,:,3) )	&
											- SUM( this%frac(i,:,:,:), 2)*rhos( this%g(i,:,1), this%g(i,1,2), this%g(i,:,3) ) )
			xps(i,3) = 1.0_mp/m%dx(3)*SUM( SUM( this%frac(i,:,:,:), 3)*rhos( this%g(i,:,1), this%g(i,:,2), this%g(i,2,3) )	&
											- SUM( this%frac(i,:,:,:), 3)*rhos( this%g(i,:,1), this%g(i,:,2), this%g(i,1,3) ) )
			xps(i,:) = - p%qs(i)/dV*xps(i,:)
		end do
	end subroutine

	subroutine Adj_forceAssign_E(this,Eps,Es)
		type(pmAssign), intent(inout) :: this
		real(mp), intent(in) :: Eps(this%n,3)
		real(mp), intent(out) :: Es(this%ng(1),this%ng(2),this%ng(3),3)
		integer :: i, g(2,3)

		Es = 0.0_mp
		do i=1,this%n
			g = this%g(i,:,:)
			Es( g(:,1), g(:,2), g(:,3), 1 ) = Es( g(:,1), g(:,2), g(:,3), 1 ) + Eps(i,1)*this%frac(i,:,:,:)
			Es( g(:,1), g(:,2), g(:,3), 2 ) = Es( g(:,1), g(:,2), g(:,3), 2 ) + Eps(i,2)*this%frac(i,:,:,:)
			Es( g(:,1), g(:,2), g(:,3), 3 ) = Es( g(:,1), g(:,2), g(:,3), 3 ) + Eps(i,3)*this%frac(i,:,:,:)
		end do
	end subroutine

	subroutine Adj_forceAssign_xp(this,m,E,Eps,xps)
		type(pmAssign), intent(inout) :: this
		type(mesh), intent(in) :: m
		real(mp), intent(in) :: E(this%ng(1),this%ng(2),this%ng(3),3)
		real(mp), intent(in) :: Eps(this%n,3)
		real(mp), intent(out) :: xps(this%n,3)
		real(mp) :: frac(2,2), Epsum(this%n), Esum(this%ng(1),this%ng(2),this%ng(3))
		integer :: i,j

		xps = 0.0_mp
		frac = 0.0_mp
		!sum : sum in each direction --- this rank will be added by the gradient of assignment
		do i=1,this%n
			frac = 1.0_mp/m%dx(1)*SUM( this%frac(i,:,:,:), 1 )
			do j=1,3
				xps(i,1) = xps(i,1) + Eps(i,j)*SUM( frac*( E(this%g(i,2,1), this%g(i,:,2), this%g(i,:,3), j)		&
														- E(this%g(i,1,1), this%g(i,:,2), this%g(i,:,3), j) )		)
			end do
			frac = 1.0_mp/m%dx(2)*SUM( this%frac(i,:,:,:), 2 )
			do j=1,3
				xps(i,2) = xps(i,2) + Eps(i,j)*SUM( frac*( E(this%g(i,:,1), this%g(i,2,2), this%g(i,:,3), j)		&
														- E(this%g(i,:,1), this%g(i,1,2), this%g(i,:,3), j) )		)
			end do
			frac = 1.0_mp/m%dx(3)*SUM( this%frac(i,:,:,:), 3 )
			do j=1,3
				xps(i,3) = xps(i,3) + Eps(i,j)*SUM( frac*( E(this%g(i,:,1), this%g(i,:,2), this%g(i,2,3), j)		&
														- E(this%g(i,:,1), this%g(i,:,2), this%g(i,1,3), j) )		)
			end do
		end do
		xps = - xps
	end subroutine

	subroutine FFTAdj(Es,rhos,W,dx)
		real(mp), intent(in) :: Es(:,:,:,:), dx(3)
		complex(mp), intent(in) :: W(size(Es,1),size(Es,2),size(Es,3))
		real(mp), intent(out) :: rhos(size(Es,1),size(Es,2),size(Es,3))
		integer*8 :: plan
		complex(mp), dimension(size(Es,1),size(Es,2),size(Es,3)) :: Esb, EsFFT, phisFFT, rhsFFT, phisb, rhosb
		integer :: L,M,N, i,wi

		if( size(Es,4).ne.3 ) then
			print *, 'FAULT!! :: Adjoint Electric field is not 3 dimension'
			stop
		end if

		L = size(Es,1)
		M = size(Es,2)
		N = size(Es,3)

		phisb =  Divergence( Es, dx, (/L,M,N/) )

		phisb = phisb*1.0_mp/L/M/N
		call dfftw_plan_dft_3d(plan,L,M,N,phisb,phisFFT,FFTW_BACKWARD,FFTW_ESTIMATE)
		call dfftw_execute_dft(plan,phisb,phisFFT)

		rhsFFT = phisFFT/W
		call dfftw_plan_dft_3d(plan,L,M,N,rhsFFT,rhosb,FFTW_FORWARD,FFTW_ESTIMATE)
		call dfftw_execute_dft(plan,rhsFFT,rhosb)
		call dfftw_destroy_plan(plan)
		rhos = REALPART(rhosb)
	end subroutine

	subroutine Adj_accel(adj,pm,dvps)
		type(adjoint), intent(inout) :: adj
		type(PM3D), intent(in) :: pm
		real(mp), intent(in) :: dvps(adj%n,3)

		adj%vps = adj%vps + (-pm%dt)*( -adj%xps +dvps )
	end subroutine

	subroutine Adj_move(adj,pm,dxps)
		type(adjoint), intent(inout) :: adj
		type(PM3D), intent(in) :: pm
		real(mp), intent(in) :: dxps(adj%n,3)

		adj%xps = adj%xps + (-pm%dt)*( dxps )
	end subroutine

end module
