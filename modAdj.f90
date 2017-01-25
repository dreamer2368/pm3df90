module modAdj

	use modPM3D

	implicit none

	type adjoint
		integer :: nt, ni, ns, ng(3)
		real(mp) :: J0, J1
		real(mp) :: dts

		type(species), allocatable :: p(:), dp(:)
		type(mesh) :: m, dm

!		real(mp), allocatable :: xps(:,:), vps(:,:), Eps(:,:), Es(:,:,:,:), rhos(:,:,:)
!		real(mp), allocatable :: dvps(:,:), dxps(:,:), dEs(:,:,:,:)
	end type

contains

	subroutine buildAdjoint(this,pm)
		type(adjoint), intent(out) :: this
		type(PM3D), intent(in) :: pm
		integer :: i
		real(mp), allocatable :: vec(:,:)

		this%ng = pm%ng
		this%ns = pm%ns
		this%nt = pm%nt
		this%ni = pm%ni
		this%dts = -pm%dt

		allocate(this%p(pm%ns))
		allocate(this%dp(pm%ns))
		do i=1,pm%ns
			call buildSpecies(this%p(i),pm%p(i)%qs,pm%p(i)%ms)
			call buildSpecies(this%dp(i),0.0_mp,0.0_mp)
			allocate(vec(pm%p(i)%np,3))
			vec = 0.0_mp
			call setSpecies(this%p(i),pm%p(i)%np,vec,vec,pm%p(i)%spwt)
			call setSpecies(this%dp(i),pm%p(i)%np,vec,vec,pm%p(i)%spwt)
			deallocate(vec)
		end do
		call buildMesh(this%m,pm%L,pm%ng)
		call buildMesh(this%dm,pm%L,pm%ng)

!		allocate(this%xps(this%n,3))
!		allocate(this%vps(this%n,3))
!		allocate(this%Eps(this%n,3))
!		allocate(this%Es(this%ng(1),this%ng(2),this%ng(3),3))
!		allocate(this%rhos(this%ng(1),this%ng(2),this%ng(3)))

!		allocate(this%dvps(this%n,3))
!		allocate(this%dxps(this%n,3))
!		allocate(this%dEs(this%ng(1),this%ng(2),this%ng(3),3))
	end subroutine

	subroutine destroyAdjoint(this)
		type(adjoint), intent(inout) :: this
		integer :: i

		call destroyMesh(this%m)
		call destroyMesh(this%dm)
		do i=1,this%ns
			call destroySpecies(this%p(i))
			call destroySpecies(this%dp(i))
		end do
		deallocate(this%p)
		deallocate(this%dp)

!		deallocate(this%xps)
!		deallocate(this%vps)
!		deallocate(this%Eps)
!		deallocate(this%Es)
!		deallocate(this%rhos)

!		deallocate(this%dvps)
!		deallocate(this%dxps)
!		deallocate(this%dEs)
	end subroutine

!	subroutine recordAdjoint(this,adj,k)
!		type(recordData), intent(inout) :: this
!		type(adjoint), intent(in) :: adj
!		integer, intent(in) :: k

!		this%xpsdata(:,:,k) = adj%xps
!		this%vpsdata(:,:,k) = adj%vps
!	end subroutine

	subroutine reset_Dadj(adj)
		type(adjoint), intent(inout) :: adj
		integer :: i

		do i=1,adj%ns
			adj%dp(i)%xp = 0.0_mp
			adj%dp(i)%vp = 0.0_mp
			adj%dp(i)%Ep = 0.0_mp
		end do
		adj%dm%E = 0.0_mp
		adj%dm%rho = 0.0_mp
		adj%dm%phi = 0.0_mp

!		adj%dvps = 0.0_mp
!		adj%dxps = 0.0_mp
!		adj%dEs = 0.0_mp
	end subroutine

	subroutine Adj_chargeAssign(this,p,m,rhos,xps)
		type(pmAssign), intent(inout) :: this
		type(species), intent(in) :: p
		type(mesh), intent(in) :: m
		real(mp), intent(in) :: rhos(this%ng(1),this%ng(2),this%ng(3))
		real(mp), intent(out) :: xps(this%np,3)
		real(mp) :: dxps(this%np,3)
		real(mp) :: dV
		integer :: i

		dV = PRODUCT(m%dx)
		dxps = 0.0_mp
		do i=1,this%np
			dxps(i,1) = p%spwt(i)/m%dx(1)*SUM( SUM( this%frac(i,:,:,:), 1)*rhos( this%g(i,2,1), this%g(i,:,2), this%g(i,:,3) )	&
											- SUM( this%frac(i,:,:,:), 1)*rhos( this%g(i,1,1), this%g(i,:,2), this%g(i,:,3) ) )
			dxps(i,2) = p%spwt(i)/m%dx(2)*SUM( SUM( this%frac(i,:,:,:), 2)*rhos( this%g(i,:,1), this%g(i,2,2), this%g(i,:,3) )	&
											- SUM( this%frac(i,:,:,:), 2)*rhos( this%g(i,:,1), this%g(i,1,2), this%g(i,:,3) ) )
			dxps(i,3) = p%spwt(i)/m%dx(3)*SUM( SUM( this%frac(i,:,:,:), 3)*rhos( this%g(i,:,1), this%g(i,:,2), this%g(i,2,3) )	&
											- SUM( this%frac(i,:,:,:), 3)*rhos( this%g(i,:,1), this%g(i,:,2), this%g(i,1,3) ) )
		end do
		dxps = - p%qs/dV*dxps
		xps = xps + dxps
	end subroutine

	subroutine Adj_forceAssign_E(this,Eps,Es)
		type(pmAssign), intent(inout) :: this
		real(mp), intent(in) :: Eps(this%np,3)
		real(mp), intent(inout) :: Es(this%ng(1),this%ng(2),this%ng(3),3)
		integer :: i, g(2,3)

!		Es = 0.0_mp
		do i=1,this%np
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
		real(mp), intent(in) :: Eps(this%np,3)
		real(mp), intent(out) :: xps(this%np,3)
		real(mp) :: dxps(this%np,3)
		real(mp) :: frac(2,2), Epsum(this%np), Esum(this%ng(1),this%ng(2),this%ng(3))
		integer :: i,j

		dxps = 0.0_mp
		frac = 0.0_mp
		!sum : sum in each direction --- this rank will be added by the gradient of assignment
		do i=1,this%np
			frac = 1.0_mp/m%dx(1)*SUM( this%frac(i,:,:,:), 1 )
			do j=1,3
				dxps(i,1) = dxps(i,1) + Eps(i,j)*SUM( frac*( E(this%g(i,2,1), this%g(i,:,2), this%g(i,:,3), j)		&
														- E(this%g(i,1,1), this%g(i,:,2), this%g(i,:,3), j) )		)
			end do
			frac = 1.0_mp/m%dx(2)*SUM( this%frac(i,:,:,:), 2 )
			do j=1,3
				dxps(i,2) = dxps(i,2) + Eps(i,j)*SUM( frac*( E(this%g(i,:,1), this%g(i,2,2), this%g(i,:,3), j)		&
														- E(this%g(i,:,1), this%g(i,1,2), this%g(i,:,3), j) )		)
			end do
			frac = 1.0_mp/m%dx(3)*SUM( this%frac(i,:,:,:), 3 )
			do j=1,3
				dxps(i,3) = dxps(i,3) + Eps(i,j)*SUM( frac*( E(this%g(i,:,1), this%g(i,:,2), this%g(i,2,3), j)		&
														- E(this%g(i,:,1), this%g(i,:,2), this%g(i,1,3), j) )		)
			end do
		end do
		dxps = - dxps

		xps = xps + dxps
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

	subroutine Adj_accel(adj)
		type(adjoint), intent(inout) :: adj
!		type(PM3D), intent(in) :: pm
!		real(mp), intent(in) :: dvps(adj%n,3)
		integer :: i

		do i=1,adj%ns
			adj%p(i)%vp = adj%p(i)%vp + adj%dts*( -adj%p(i)%xp + adj%dp(i)%vp )
		end do
	end subroutine

	subroutine Adj_move(adj)
		type(adjoint), intent(inout) :: adj
!		type(PM3D), intent(in) :: pm
!		real(mp), intent(in) :: dxps(adj%n,3)
		integer :: i

		do i=1,adj%ns
			adj%p(i)%xp = adj%p(i)%xp + adj%dts*( adj%dp(i)%xp )
		end do
	end subroutine

end module
