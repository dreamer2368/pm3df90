module modPlasma

	use constants

	implicit none

	type plasma
		integer :: nt
		integer :: ni
		integer :: n
		integer :: ng(3)

		real(mp), allocatable :: qs(:), ms(:)

		real(mp), allocatable :: time(:)
		real(mp), allocatable :: xp(:,:)
		real(mp), allocatable :: vp(:,:)
		real(mp), allocatable :: E(:,:,:,:)
		real(mp), allocatable :: rho(:,:,:)
		real(mp), allocatable :: phi(:,:,:)

		integer, allocatable :: g(:,:,:)
		real(mp), allocatable :: frac(:,:,:,:)

		real(mp), allocatable :: Ep(:,:)
		real(mp), allocatable :: xpdata(:,:,:)
		real(mp), allocatable :: vpdata(:,:,:)
		real(mp), allocatable :: Exdata(:,:,:,:)

		real(mp), allocatable :: xpsdata(:,:,:)
	end type

contains

	subroutine buildPlasma(this,nt,ni,n,ng,xp0,vp0,qs,ms)

		type(plasma), intent(out) :: this
		integer, intent(in) :: nt
		integer, intent(in) :: ni
		integer, intent(in) :: n
		integer, intent(in) :: ng(3)
		real(mp), intent(in) :: xp0(n,3)
		real(mp), intent(in) :: vp0(n,3)
		real(mp), intent(in) :: qs(n), ms(n)

		this%nt = nt
		this%ni = ni
		this%n = n
		this%ng = ng

		allocate(this%qs(n))
		allocate(this%ms(n))

		this%qs = qs
		this%ms = ms

		allocate(this%time(nt))
		allocate(this%xp(n,3))
		allocate(this%vp(n,3))
		allocate(this%E(ng(1),ng(2),ng(3),3))
		allocate(this%phi(ng(1),ng(2),ng(3)))
		allocate(this%rho(ng(1),ng(2),ng(3)))

		allocate(this%g(n,2,3))
		allocate(this%frac(n,2,2,2))

		allocate(this%Ep(n,3))

		this%xp = xp0
		this%vp = vp0

		allocate(this%xpdata(n,3,nt))
		allocate(this%vpdata(n,3,nt))
		allocate(this%Exdata(ng(1),ng(2),ng(3),nt))

		allocate(this%xpsdata(n,3,nt))

	end subroutine

	subroutine destroyPlasma(this)
		type(plasma), intent(inout) :: this

		deallocate(this%qs)
		deallocate(this%ms)
		deallocate(this%time)
		deallocate(this%xp)
		deallocate(this%vp)
		deallocate(this%E)
		deallocate(this%rho)
		deallocate(this%phi)

		deallocate(this%g)
		deallocate(this%frac)

		deallocate(this%Ep)
		deallocate(this%xpdata)
		deallocate(this%vpdata)
		deallocate(this%Exdata)

		deallocate(this%xpsdata)
	end subroutine

	subroutine recordPlasma(this,k)

		type(plasma), intent(inout) :: this
		integer, intent(in) :: k

		this%xpdata(:,:,k) = this%xp
		this%vpdata(:,:,k) = this%vp
		this%Exdata(:,:,:,k) = this%E(:,:,:,1)

	end subroutine

	subroutine printPlasma(this)
		type(plasma), intent(in) :: this
		integer :: i

		open(unit=301,file='data/record.out',status='replace')
		write(301,*) this%n, this%ng, this%nt
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
		type(plasma), intent(inout) :: this
		integer, intent(in) :: k
		real(mp), intent(in) :: xps(this%n,3)

		this%xpsdata(:,:,k) = xps
	end subroutine

end module