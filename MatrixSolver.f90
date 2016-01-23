module MatrixSolver

	use constants										!machine precision mp should be assigned

	implicit none

	include 'fftw3.f'

contains

	subroutine FFTAdj(Es,rhos,W)
		real(mp), intent(in) :: Es(:,:,:,:)
		complex(mp), intent(in) :: W(size(Es,1),size(Es,2),size(Es,3))
		real(mp), intent(out) :: rhos(size(Es,1),size(Es,2),size(Es,3))
		integer*8 :: plan
		complex(mp), dimension(size(Es,1),size(Es,2),size(Es,3)) :: Esb, EsFFT, phisFFT, rhsFFT, rhosb
		integer :: L,M,N, i,wi

		if( size(Es,4).ne.3 ) then
			print *, 'FAULT!! :: Adjoint Electric field is not 3 dimension'
			stop
		end if

		L = size(Es,1)
		M = size(Es,2)
		N = size(Es,3)

		phisFFT = 0.0_mp
		Esb = Es(:,:,:,1)
		call dfftw_plan_dft_3d(plan,L,M,N,Esb,EsFFT,FFTW_BACKWARD,FFTW_ESTIMATE)
		call dfftw_execute_dft(plan,Esb,EsFFT)
		call dfftw_destroy_plan(plan)
		do i=0,L-1
			if( i.le.L/2 ) then
				wi = i
			else
				wi = - ( L-i )
			end if
			phisFFT(i+1,:,:) = phisFFT(i+1,:,:) + EsFFT(i+1,:,:)*2.0_mp*pi*eye*wi
		end do

		Esb = Es(:,:,:,2)
		call dfftw_plan_dft_3d(plan,L,M,N,Esb,EsFFT,FFTW_BACKWARD,FFTW_ESTIMATE)
		call dfftw_execute_dft(plan,Esb,EsFFT)
		call dfftw_destroy_plan(plan)
		do i=0,M-1
			if( i.le.M/2 ) then
				wi = i
			else
				wi = - ( M-i )
			end if
			phisFFT(:,i+1,:) = phisFFT(:,i+1,:) + EsFFT(:,i+1,:)*2.0_mp*pi*eye*wi
		end do

		Esb = Es(:,:,:,3)
		call dfftw_plan_dft_3d(plan,L,M,N,Esb,EsFFT,FFTW_BACKWARD,FFTW_ESTIMATE)
		call dfftw_execute_dft(plan,Esb,EsFFT)
		call dfftw_destroy_plan(plan)
		do i=0,N-1
			if( i.le.N/2 ) then
				wi = i
			else
				wi = - ( N-i )
			end if
			phisFFT(:,:,i+1) = phisFFT(:,:,i+1) + EsFFT(:,:,i+1)*2.0_mp*pi*eye*wi
		end do

		phisFFT = phisFFT*1.0_mp/L/M/N

		rhsFFT = phisFFT/W
		call dfftw_plan_dft_3d(plan,L,M,N,rhsFFT,rhosb,FFTW_FORWARD,FFTW_ESTIMATE)
		call dfftw_execute_dft(plan,rhsFFT,rhosb)
		call dfftw_destroy_plan(plan)
		rhos = REALPART(rhosb)
	end subroutine

	subroutine FFTEfield(y,rhs,W)
		real(mp), intent(in) :: rhs(:,:,:)
		complex(mp), intent(in) :: W(size(rhs,1),size(rhs,2),size(rhs,3))
		real(mp), intent(out) :: y(size(rhs,1),size(rhs,2),size(rhs,3),3)			!!Gradient of the solution
		integer*8 :: plan
		complex(mp), dimension(size(rhs,1),size(rhs,2),size(rhs,3)) :: rhsb, rhsFFT, xFFT, yb, yFFT
		integer :: L,M,N, i,wi

		L = size(rhs,1)
		M = size(rhs,2)
		N = size(rhs,3)

		rhsb = rhs
		call dfftw_plan_dft_3d(plan,L,M,N,rhsb,rhsFFT,FFTW_FORWARD,FFTW_ESTIMATE)
		call dfftw_execute_dft(plan,rhsb,rhsFFT)
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
			yFFT(:,:,i+1) = xFFT(:,:,i+1)*2.0_mp*pi*eye*wi
		end do

		call dfftw_plan_dft_3d(plan,L,M,N,yFFT,yb,FFTW_BACKWARD,FFTW_ESTIMATE)
		call dfftw_execute_dft(plan,yFFT,yb)
		call dfftw_destroy_plan(plan)
		y(:,:,:,3) = REALPART(yb)*1.0_mp/L/M/N

		!gradient in y direction
		do i=0,M-1
			if( i.le.M/2 ) then
				wi = i
			else
				wi = - ( M-i )
			end if
			yFFT(:,i+1,:) = xFFT(:,i+1,:)*2.0_mp*pi*eye*wi
		end do

		call dfftw_plan_dft_3d(plan,L,M,N,yFFT,yb,FFTW_BACKWARD,FFTW_ESTIMATE)
		call dfftw_execute_dft(plan,yFFT,yb)
		call dfftw_destroy_plan(plan)
		y(:,:,:,2) = REALPART(yb)*1.0_mp/L/M/N

		!gradient in x direction
		do i=0,L-1
			if( i.le.L/2 ) then
				wi = i
			else
				wi = - ( L-i )
			end if
			yFFT(i+1,:,:) = xFFT(i+1,:,:)*2.0_mp*pi*eye*wi
		end do

		call dfftw_plan_dft_3d(plan,L,M,N,yFFT,yb,FFTW_BACKWARD,FFTW_ESTIMATE)
		call dfftw_execute_dft(plan,yFFT,yb)
		call dfftw_destroy_plan(plan)
		y(:,:,:,1) = REALPART(yb)*1.0_mp/L/M/N
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
		complex(mp), intent(out) :: W(N(1),N(2),N(3))
		integer :: i,j,k, wi,wj,wk
		complex(mp) :: wx,wy,wz

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
				do i=0,N(1)
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
		complex(mp), intent(in) :: W(size(rhs,1),size(rhs,2),size(rhs,3))
		real(mp), intent(out) :: x(size(rhs,1),size(rhs,2),size(rhs,3))
		integer*8 :: plan
		complex(mp), dimension(size(rhs,1),size(rhs,2),size(rhs,3)) :: rhsFFT, xFFT, rhsb, xb
		integer :: L,M,N

		L = size(rhs,1)
		M = size(rhs,2)
		N = size(rhs,3)

		rhsb = rhs
		call dfftw_plan_dft_3d(plan,L,M,N,rhsb,rhsFFT,FFTW_FORWARD,FFTW_ESTIMATE)
		call dfftw_execute_dft(plan,rhsb,rhsFFT)
		call dfftw_destroy_plan(plan)

		xFFT = rhsFFT/W
		xFFT(1,1,1) = (0.0_mp, 0.0_mp)

		call dfftw_plan_dft_3d(plan,L,M,N,xFFT,xb,FFTW_BACKWARD,FFTW_ESTIMATE)
		call dfftw_execute_dft(plan,xFFT,xb)
		call dfftw_destroy_plan(plan)

		x = REALPART(xb)*1.0_mp/L/M/N
	end subroutine

end module