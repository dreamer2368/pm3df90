module modMesh

	use MatrixSolver

	implicit none

	type mesh
		integer :: ng(3)
		real(mp) :: L(3), dx(3), rho_back
		complex(mp), allocatable :: W(:,:,:)

		real(mp), allocatable :: E(:,:,:,:)
		real(mp), allocatable :: rho(:,:,:)
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
		this%rho_back = 0.0_mp

		allocate(this%E(ng(1),ng(2),ng(3),3))
		allocate(this%phi(ng(1),ng(2),ng(3)))
		allocate(this%rho(ng(1),ng(2),ng(3)))

		allocate(this%W(ng(1),ng(2),ng(3)))

		call FFTPoisson_setup(ng,this%W)
	end subroutine

	subroutine setMesh(this,rho_back)
		type(mesh), intent(out) :: this
		real(mp), intent(in) :: rho_back

		this%rho_back = rho_back
	end subroutine

	subroutine destroyMesh(this)
		type(mesh), intent(inout) :: this

		deallocate(this%E)
		deallocate(this%rho)
		deallocate(this%phi)

		deallocate(this%W)
	end subroutine

end module