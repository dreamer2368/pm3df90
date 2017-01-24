module modPM3D

	use modSpecies
	use modMesh
	use modAssign

	implicit none

	!Defined here to resolve circular dependency
	type recordData
		integer :: nt, ns, ng(3), mod_r		!mod_r: number of timesteps to save data
		real(mp) :: L(3), dx(3)
		character(len=:), allocatable :: dir

		integer, allocatable :: np(:,:)

		real(mp), allocatable :: Edata(:,:,:,:,:)
		real(mp), allocatable :: rhodata(:,:,:,:)
		real(mp), allocatable :: phidata(:,:,:,:)
		real(mp), allocatable :: PE(:), KE(:,:)
	end type

	type PM3D
		integer :: nt, ni, ns, ng(3)
		real(mp) :: L(3), eps0, wp
		real(mp) :: dt
		real(mp), allocatable :: A0(:)

		type(species), allocatable :: p(:)
		type(mesh) :: m
		type(recordData) :: r
		type(pmAssign), allocatable :: a(:)
	end type

contains

	subroutine buildPM3D(this,Tf,Ti,Ng,Ns,dt,L,A,dir,mod_input)
		type(PM3D), intent(out) :: this
		real(mp), intent(in) :: Tf,Ti
		integer, intent(in) :: Ng(3), Ns
		real(mp), intent(in), optional :: dt, A(:), L(3)
		character(len=*), intent(in), optional :: dir
		integer, intent(in), optional :: mod_input
		integer :: i, mod_r
		if( present(dt) ) then
			this%dt = dt
		else
			this%dt = 0.2_mp
		end if
		if( present(A) ) then
			allocate(this%A0(size(A)))
			this%A0 = A
		else
			allocate(this%A0(1))
			this%A0 = 1.0_mp
		end if
		if( present(L) ) then
			this%L = L
		else
			this%L = 2*pi/( sqrt(3.0_mp)/2.0_mp/sqrt(2.0_mp)/0.2_mp )
		end if
		if( present(mod_input) ) then
			mod_r = mod_input
		else
			mod_r = 1
		end if
		this%ng = Ng
		this%ns = Ns
		this%nt = CEILING(Tf/this%dt)
		this%dt = Tf/this%nt
		this%ni = FLOOR(Ti/this%dt)+1
      print *, 'Plasma is created'
      print *, 'L = (',this%L,')'
      print *, 'Ng = (',Ng,')'
      print *, 'Ns = ',Ns,', A = ',this%A0
		print *, 'Ni = ',this%ni,', Nt = ',this%nt,', dt = ',this%dt

		this%eps0 = 1.0_mp
		this%wp = 1.0_mp
		allocate(this%p(Ns))
		allocate(this%a(Ns))
!		call buildPlasma(this%p,N)
		call buildMesh(this%m,this%L,Ng)
		do i=1,Ns
			call buildAssign(this%a(i),Ng)
		end do
		if( present(dir) ) then
			call buildRecord(this%r,this%nt,this%ns,this%L,this%ng,mod_r,dir)
		else
			call buildRecord(this%r,this%nt,this%ns,this%L,this%ng,mod_r)
		end if
	end subroutine

	subroutine destroyPM3D(this)
		type(PM3D), intent(inout) :: this
		integer :: i

		deallocate(this%A0)

		do i=1,this%ns
			call destroySpecies(this%p(i))
			call destroyAssign(this%a(i))
		end do
		deallocate(this%p)
		deallocate(this%a)
		call destroyMesh(this%m)
		call destroyRecord(this%r)
	end subroutine

	subroutine buildRecord(this,nt,ns,L,ng,mod_r,input_dir)
		type(recordData), intent(out) :: this
		integer, intent(in) :: nt, ns, ng(3), mod_r
		real(mp), intent(in) :: L(3)
		character(len=*), intent(in), optional :: input_dir
		integer :: nr
		nr = nt/mod_r+1

		this%nt = nt
		this%ns = ns
		this%L = L
		this%ng = ng
		this%mod_r = mod_r

		allocate(this%np(ns,nr))
		allocate(this%Edata(ng(1),ng(2),ng(3),3,nr))
		allocate(this%rhodata(ng(1),ng(2),ng(3),nr))
		allocate(this%phidata(ng(1),ng(2),ng(3),nr))
		allocate(this%PE(nr))
		allocate(this%KE(ns,nr))

		this%np = 0
		this%Edata = 0.0_mp
		this%rhodata = 0.0_mp
		this%phidata = 0.0_mp
		this%PE = 0.0_mp
		this%KE = 0.0_mp

		if( present(input_dir) ) then
			allocate(character(len=len(input_dir)) :: this%dir)
			this%dir = input_dir
		else
			allocate(character(len=0) :: this%dir)
			this%dir = ''
		end if

		call system('mkdir -p data/'//this%dir//'/xp')
		call system('mkdir -p data/'//this%dir//'/vp')
!		call system('mkdir -p data/'//this%dir//'/Ep')

		call system('rm data/'//this%dir//'/xp/*.*')
		call system('rm data/'//this%dir//'/vp/*.*')
!		call system('rm data/'//this%dir//'/Ep/*.*')
	end subroutine

	subroutine destroyRecord(this)
		type(recordData), intent(inout) :: this

		deallocate(this%np)
		deallocate(this%Edata)
		deallocate(this%rhodata)
		deallocate(this%phidata)
		deallocate(this%PE)
		deallocate(this%KE)
		deallocate(this%dir)
	end subroutine

end module
