module modRecord

	use modPlasma
	use modMesh

	implicit none

	type recordData
		integer :: nt, n, ng(3)
		real(mp) :: L(3), dx(3)

		real(mp), allocatable :: xpdata(:,:,:)
		real(mp), allocatable :: vpdata(:,:,:)
		real(mp), allocatable :: xpsdata(:,:,:)
		real(mp), allocatable :: Exdata(:,:,:,:)
	end type

contains

	subroutine buildRecord(this,nt,n,L,ng)

		type(recordData), intent(out) :: this
		integer, intent(in) :: nt, n, ng(3)
		real(mp), intent(in) :: L(3)

		this%nt = nt
		this%n = n
		this%L = L
		this%ng = ng

		allocate(this%xpdata(n,3,nt))
		allocate(this%vpdata(n,3,nt))
		allocate(this%xpsdata(n,3,nt))
		allocate(this%Exdata(ng(1),ng(2),ng(3),nt))
	end subroutine

	subroutine destroyRecord(this)
		type(recordData), intent(inout) :: this

		deallocate(this%xpdata)
		deallocate(this%vpdata)
		deallocate(this%xpsdata)
		deallocate(this%Exdata)
	end subroutine

	subroutine recordPlasma(this,p,m,k)
		type(recordData), intent(inout) :: this
		type(plasma), intent(in) :: p
		type(mesh), intent(in) :: m
		integer, intent(in) :: k

		this%xpdata(:,:,k) = p%xp
		this%vpdata(:,:,k) = p%vp
		this%Exdata(:,:,:,k) = m%E(:,:,:,1)
	end subroutine

	subroutine printPlasma(this)
		type(recordData), intent(in) :: this
		integer :: i

		open(unit=301,file='data/record.out',status='replace')
		write(301,*) this%n, this%ng, this%nt, this%L
		close(301)

		open(unit=301,file='data/xp.bin',status='replace',form='unformatted',access='stream')
		open(unit=302,file='data/vp.bin',status='replace',form='unformatted',access='stream')
		open(unit=303,file='data/Ex.bin',status='replace',form='unformatted',access='stream')
		do i = 1,this%nt
			write(301) this%xpdata(:,:,i)
			write(302) this%vpdata(:,:,i)
			write(303) this%Exdata(:,:,:,i)
		end do
		close(301)
		close(302)
		close(303)
	end subroutine

	subroutine recordAdjoint(this,k,xps)
		type(recordData), intent(inout) :: this
		integer, intent(in) :: k
		real(mp), intent(in) :: xps(this%n,3)

		this%xpsdata(:,:,k) = xps
	end subroutine

end module