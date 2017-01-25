module modRecord

	use modPM3D	

	implicit none

contains

	subroutine recordPlasma(this,pm,k)
		type(recordData), intent(inout) :: this
		type(PM3D), intent(in) :: pm
		integer, intent(in) :: k					!k : time step
		integer :: ns,j, kr								!ns : species
		character(len=100) :: nstr, kstr

		if( (this%mod_r.eq.1) .or. (mod(k,this%mod_r).eq.0) ) then
			kr = merge(k,k/this%mod_r,this%mod_r.eq.1)
			do ns=1,pm%ns
				write(nstr,*) ns
				write(kstr,*) kr
				open(unit=305,file='data/'//this%dir//'/xp/'//trim(adjustl(kstr))//'_'	&
					//trim(adjustl(nstr))//'.bin',status='replace',form='unformatted',access='stream')
				open(unit=306,file='data/'//this%dir//'/vp/'//trim(adjustl(kstr))//'_'	&
					//trim(adjustl(nstr))//'.bin',status='replace',form='unformatted',access='stream')
				open(unit=307,file='data/'//this%dir//'/spwt/'//trim(adjustl(kstr))//'_'	&
					//trim(adjustl(nstr))//'.bin',status='replace',form='unformatted',access='stream')
!				open(unit=308,file='data/'//this%dir//'/Ep/'//trim(adjustl(kstr))//'_'	&
!					//trim(adjustl(nstr))//'.bin',status='replace',form='unformatted',access='stream')
				write(305) pm%p(ns)%xp
				write(306) pm%p(ns)%vp
				write(307) pm%p(ns)%spwt
!				write(308) pm%p(ns)%Ep
				close(305)
				close(306)
				close(307)
!				close(308)
				this%np(ns,kr+1) = pm%p(ns)%np
				this%KE(ns,kr+1) = 0.5_mp*pm%p(ns)%ms*SUM(pm%p(ns)%spwt*SUM(pm%p(ns)%vp**2,2))
			end do

			this%phidata(:,:,:,kr+1) = pm%m%phi
			this%Edata(:,:,:,:,kr+1) = pm%m%E
			this%rhodata(:,:,:,kr+1) = pm%m%rho
			this%PE(kr+1) = 0.5_mp*SUM(pm%m%E**2)*PRODUCT(pm%m%dx)

			print *, '============= ',k,'-th Time Step ================='
			do ns=1,pm%ns
				print *, 'Species(',ns,'): ',pm%p(ns)%np, ', KE: ', this%KE(ns,kr+1),'J'
			end do
		end if
	end subroutine

	subroutine printPlasma(this)
		type(recordData), intent(in) :: this
		character(len=100) :: s
		integer :: i,j

		open(unit=300,file='data/'//this%dir//'/record',status='replace')
		open(unit=301,file='data/'//this%dir//'/E.bin',status='replace',form='unformatted',access='stream')
		open(unit=302,file='data/'//this%dir//'/rho.bin',status='replace',form='unformatted',access='stream')
		open(unit=303,file='data/'//this%dir//'/PE.bin',status='replace',form='unformatted',access='stream')
		open(unit=304,file='data/'//this%dir//'/Np.bin',status='replace',form='unformatted',access='stream')
		open(unit=305,file='data/'//this%dir//'/phi.bin',status='replace',form='unformatted',access='stream')
		do i=1,this%ns
			write(s,*) i
			open(unit=307+i,file='data/'//this%dir//'/KE_'//trim(adjustl(s))//'.bin',status='replace',form='unformatted',access='stream')
		end do

		write(300,*) this%ns, this%ng, this%nt, this%L, this%mod_r
		close(300)

		do i = 1,this%nt/this%mod_r+1
			write(301) this%Edata(:,:,:,1,i)
			write(301) this%Edata(:,:,:,2,i)
			write(301) this%Edata(:,:,:,3,i)
			write(302) this%rhodata(:,:,:,i)
			write(303) this%PE(i)
			write(304) this%np(:,i)
			write(305) this%phidata(:,:,:,i)
			do j=1,this%ns
				write(307+j) this%KE(j,i)
			end do
		end do

		close(301)
		close(302)
		close(303)
		close(304)
		close(305)
		do i=1,this%ns
			close(307+i)
		end do
	end subroutine

end module
