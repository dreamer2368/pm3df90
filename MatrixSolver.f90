module MatrixSolver

!	use modPlasma										!machine precision mp should be assigned

	implicit none

	include 'fftw3.f'

contains

	subroutine FFTPoissonGrad(y,rhs,W)
		real(mp), intent(in) :: rhs(:,:,:)
		complex(mp), intent(in) :: W(size(rhs,1)/2+1,size(rhs,2),size(rhs,3))
		real(mp), intent(out) :: y(size(rhs,1),size(rhs,2),size(rhs,3),3)			!!Gradient of the solution
		integer*8 :: plan
		complex(mp) :: rhsFFT(size(rhs,1)/2+1,size(rhs,2),size(rhs,3))
		complex(mp) :: xFFT(size(rhs,1)/2+1,size(rhs,2),size(rhs,3))				!!Spectral solution of Poisson equation
		complex(mp) :: yFFT(size(rhs,1)/2+1,size(rhs,2),size(rhs,3),3)
		integer :: L,M,N, i,wi

		L = size(rhs,1)
		M = size(rhs,2)
		N = size(rhs,3)

		call dfftw_plan_dft_r2c_3d(plan,L,M,N,rhs,rhsFFT,FFTW_ESTIMATE)
		call dfftw_execute_dft_r2c(plan,rhs,rhsFFT)
		call dfftw_destroy_plan(plan)

		xFFT = rhsFFT/W
		xFFT(1,1,1) = (0.0_mp, 0.0_mp)

		!gradient in z direction
		do i=0,N-1
			if( i.le.N/2 ) then
				wi = i
			else
				wi = - ( N-i )
			end if
			yFFT(:,:,i+1,3) = xFFT(:,:,i+1)*2.0_mp*pi*eye*wi
		end do

		!gradient in y direction
		do i=0,M-1
			if( i.le.M/2 ) then
				wi = i
			else
				wi = - ( M-i )
			end if
			yFFT(:,i+1,:,2) = xFFT(:,i+1,:)*2.0_mp*pi*eye*wi
		end do

		!gradient in x direction
		do i=0,L/2
			wi = i
			yFFT(i+1,:,:,1) = xFFT(i+1,:,:)*2.0_mp*pi*eye*wi
		end do

		call dfftw_plan_dft_c2r_3d(plan,L,M,N,yFFT(:,:,:,1),y(:,:,:,1),FFTW_ESTIMATE)
		call dfftw_execute_dft_c2r(plan,yFFT(:,:,:,1),y(:,:,:,1))
		call dfftw_destroy_plan(plan)

		call dfftw_plan_dft_c2r_3d(plan,L,M,N,yFFT(:,:,:,2),y(:,:,:,2),FFTW_ESTIMATE)
		call dfftw_execute_dft_c2r(plan,yFFT(:,:,:,2),y(:,:,:,2))
		call dfftw_destroy_plan(plan)

		call dfftw_plan_dft_c2r_3d(plan,L,M,N,yFFT(:,:,:,3),y(:,:,:,3),FFTW_ESTIMATE)
		call dfftw_execute_dft_c2r(plan,yFFT(:,:,:,3),y(:,:,:,3))
		call dfftw_destroy_plan(plan)

		y = y*1.0_mp/L/M/N
	end subroutine

	function Gradient(x,dx,Ng) result(y)						!Derivative with periodic BC
		real(mp), intent(in) :: x(:,:,:)						!(nx,ny,nz)
		real(mp), intent(in) :: dx(3)							!dx, dy, dz
		integer, intent(in) :: Ng(3)							!nx, ny, nz
		real(mp) :: y(size(x,1),size(x,2),size(x,3),3)
		integer :: i

		y=0.0_mp
		do i=2,Ng(1)-1
			y(i,:,:,1) = 0.5_mp/dx(1)*( x(i+1,:,:) - x(i-1,:,:) )
		end do
		y(1,:,:,1) = 0.5_mp/dx(1)*( x(2,:,:) - x(Ng(1),:,:) )
		y(Ng(1),:,:,1) = 0.5_mp/dx(1)*( x(1,:,:) - x(Ng(1)-1,:,:) )

		do i=2,Ng(2)-1
			y(:,i,:,2) = 0.5_mp/dx(2)*( x(:,i+1,:) - x(:,i-1,:) )
		end do
		y(:,1,:,2) = 0.5_mp/dx(2)*( x(:,2,:) - x(:,Ng(2),:) )
		y(:,Ng(2),:,2) = 0.5_mp/dx(2)*( x(:,1,:) - x(:,Ng(2)-1,:) )

		do i=2,Ng(3)-1
			y(:,:,i,3) = 0.5_mp/dx(3)*( x(:,:,i+1) - x(:,:,i-1) )
		end do
		y(:,:,1,3) = 0.5_mp/dx(3)*( x(:,:,2) - x(:,:,Ng(3)) )
		y(:,:,Ng(3),3) = 0.5_mp/dx(3)*( x(:,:,1) - x(:,:,Ng(3)-1) )
	end function

	subroutine FFTPoisson_setup(N,dx,W)
		integer, intent(in) :: N(3)
		real(mp), intent(in) :: dx(3)
		complex(mp), intent(out) :: W(N(1)/2+1,N(2),N(3))
		integer :: i,j,k, wi,wj,wk
		complex(mp) :: wx,wy,wz

!		wx = exp(2.0_mp*pi*eye/N(1))
!		wy = exp(2.0_mp*pi*eye/N(2))
!		wz = exp(2.0_mp*pi*eye/N(3))
		wx = 2.0_mp*pi*eye
		wy = 2.0_mp*pi*eye
		wz = 2.0_mp*pi*eye
		do k=0,N(3)-1
			if( k.le.N(3)/2 ) then
				wk = k
			else
				wk = - ( N(3)-k )
			end if
			do j=0,N(2)-1
				if( j.le.N(2)/2 ) then
					wj = j
				else
					wj = - ( N(2)-j )
				end if
				do i=0,N(1)/2
!					W(i+1,j+1,k+1) = ( wx**i - 2.0_mp + wx**(-i) )/dx(1)/dx(1)	&
!								+ ( wy**j - 2.0_mp + wy**(-j) )/dx(2)/dx(2)	&
!								+ ( wz**k - 2.0_mp + wz**(-k) )/dx(3)/dx(3)
					if( i.le.N(1)/2 ) then
						wi = i
					else
						wi = - ( N(1)-i )
					end if
					W(i+1,j+1,k+1) = + (wx*wi)**2.0_mp + (wy*wj)**2.0_mp + (wz*wk)**2.0_mp
				end do
			end do
		end do
		W(1,1,1) = 1.0_mp
	end subroutine

	subroutine FFTPoisson(x,rhs,W)
		real(mp), intent(in) :: rhs(:,:,:)
		complex(mp), intent(in) :: W(size(rhs,1)/2+1,size(rhs,2),size(rhs,3))
		real(mp), intent(out) :: x(size(rhs,1),size(rhs,2),size(rhs,3))
		integer*8 :: plan
		complex(mp) :: rhsFFT(size(rhs,1)/2+1,size(rhs,2),size(rhs,3))
		complex(mp) :: xFFT(size(rhs,1)/2+1,size(rhs,2),size(rhs,3))
		integer :: L,M,N

		L = size(rhs,1)
		M = size(rhs,2)
		N = size(rhs,3)

		call dfftw_plan_dft_r2c_3d(plan,L,M,N,rhs,rhsFFT,FFTW_ESTIMATE)
		call dfftw_execute_dft_r2c(plan,rhs,rhsFFT)
		call dfftw_destroy_plan(plan)

		xFFT = rhsFFT/W
		xFFT(1,1,1) = (0.0_mp, 0.0_mp)

		call dfftw_plan_dft_c2r_3d(plan,L,M,N,xFFT,x,FFTW_ESTIMATE)
		call dfftw_execute_dft_c2r(plan,xFFT,x)
		call dfftw_destroy_plan(plan)

		x = x*1.0_mp/L/M/N
	end subroutine

end module