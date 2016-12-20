module modMesh

	use MatrixSolver

	implicit none

	type mesh
		integer :: ng(3)
		real(mp) :: L(3), dx(3)
		complex(mp), allocatable :: W(:,:,:)

		real(mp), allocatable :: E(:,:,:,:)
		real(mp), allocatable :: rho(:,:,:), rho_back(:,:,:)
		real(mp), allocatable :: phi(:,:,:)
	end type

contains

	subroutine buildMesh(this,L,ng)
		type(mesh), intent(out) :: this
		integer, intent(in) :: ng(3)
		real(mp), intent(in) :: L(3)

		this%L = L
		this%ng = ng
		this%dx = L/ng

		allocate(this%E(ng(1),ng(2),ng(3),3))
		allocate(this%phi(ng(1),ng(2),ng(3)))
		allocate(this%rho(ng(1),ng(2),ng(3)))
		allocate(this%rho_back(ng(1),ng(2),ng(3)))
		this%rho_back = 0.0_mp

		allocate(this%W(ng(1),ng(2),ng(3)))
!		allocate(this%W(ng(1)/2+1,ng(2),ng(3)))

		call FFTPoisson_setup(ng,this%W,L)
!		call FFTPoisson_setup2(ng,this%dx,this%W)
	end subroutine

	subroutine setMesh(this,rho_back)
		type(mesh), intent(inout) :: this
		real(mp), intent(in) :: rho_back(this%ng(1),this%ng(2),this%ng(3))

		this%rho_back = rho_back
	end subroutine

	subroutine destroyMesh(this)
		type(mesh), intent(inout) :: this

		deallocate(this%E)
		deallocate(this%rho)
		deallocate(this%phi)
		deallocate(this%rho_back)

		deallocate(this%W)
	end subroutine

end module
