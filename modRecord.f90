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
		real(mp), allocatable :: vpsdata(:,:,:)

		real(mp), allocatable :: Epdata(:,:,:)
		real(mp), allocatable :: Edata(:,:,:,:,:)
		real(mp), allocatable :: PE(:), KE(:)
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
		allocate(this%vpsdata(n,3,nt))

		allocate(this%Epdata(n,3,nt))
		allocate(this%Edata(ng(1),ng(2),ng(3),3,nt))

		allocate(this%PE(nt))
		allocate(this%KE(nt))
	end subroutine

	subroutine destroyRecord(this)
		type(recordData), intent(inout) :: this

		deallocate(this%xpdata)
		deallocate(this%vpdata)
		deallocate(this%xpsdata)
		deallocate(this%vpsdata)

		deallocate(this%Epdata)
		deallocate(this%Edata)

		deallocate(this%PE)
		deallocate(this%KE)
	end subroutine

	subroutine recordPlasma(this,p,m,k)
		type(recordData), intent(inout) :: this
		type(plasma), intent(in) :: p
		type(mesh), intent(in) :: m
		integer, intent(in) :: k

		this%xpdata(:,:,k) = p%xp
		this%vpdata(:,:,k) = p%vp
		this%Epdata(:,:,k) = p%Ep
		this%Edata(:,:,:,:,k) = m%E
		this%PE(k) = 0.5_mp*SUM(m%E**2)*PRODUCT(m%dx)
		this%KE(k) = 0.5_mp*SUM(p%ms*SUM(p%vp**2,2))
	end subroutine

	subroutine printPlasma(this,str)
		type(recordData), intent(in) :: this
		character(len=*), intent(in), optional :: str
		integer :: i

		if( present(str) ) then
			open(unit=300,file='data/record'//str,status='replace')
			open(unit=301,file='data/xp'//str//'.bin',status='replace',form='unformatted',access='stream')
			open(unit=302,file='data/vp'//str//'.bin',status='replace',form='unformatted',access='stream')
			open(unit=303,file='data/E'//str//'.bin',status='replace',form='unformatted',access='stream')
			open(unit=304,file='data/PE'//str//'.bin',status='replace',form='unformatted',access='stream')
			open(unit=305,file='data/KE'//str//'.bin',status='replace',form='unformatted',access='stream')
		else
			open(unit=300,file='data/record',status='replace')
			open(unit=301,file='data/xp.bin',status='replace',form='unformatted',access='stream')
			open(unit=302,file='data/vp.bin',status='replace',form='unformatted',access='stream')
			open(unit=303,file='data/E.bin',status='replace',form='unformatted',access='stream')
			open(unit=304,file='data/PE.bin',status='replace',form='unformatted',access='stream')
			open(unit=305,file='data/KE.bin',status='replace',form='unformatted',access='stream')
		end if

		write(300,*) this%n, this%ng, this%nt, this%L
		close(300)

		do i = 1,this%nt
			write(301) this%xpdata(:,:,i)
			write(302) this%vpdata(:,:,i)

!			write(303) this%Epdata(:,:,i)
			write(303) this%Edata(:,:,:,1,i)
			write(303) this%Edata(:,:,:,2,i)
			write(303) this%Edata(:,:,:,3,i)
			write(304) this%PE(i)
			write(305) this%KE(i)
		end do
		close(301)
		close(302)
		close(303)
		close(304)
		close(305)
	end subroutine

end module
