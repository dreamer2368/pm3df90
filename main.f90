program main

	use testmodule
	use twoparticle
	use Landau

	implicit none

	real(mp) :: start, finish

	! print to screen
	print *, 'calling program main'

	call cpu_time(start)

!	call twoParticle_rep_AdjTest(0.2_mp,(/32,32,32/),MKE,dMKE,IC_wave,dIC_wave,dIC_wave_dB)
!	call twoParticleTraj(1.0_mp,(/32,32,32/))
!	call twostream(0.2_mp,(/64,64,64/),(/2**8,2**5,2**5/))
	call test_fullAdjoint(0.2_mp,(/32,32,32/),(/2**8,2**5,2**5/),MKE,dMKE,V_Th,dV_Th,dV_Th_dB)
!	call test_particle_adj3((/64,64,64/),2,TestQoI_singleStep,dTestQoI_singleStep)
!	call test_particle_adj2((/64,64,64/),2,TestQoI_singleStep)
!	call test_checkpoint((/64,64,64/),2,TestQoI_singleStep)
!	call test_checkpoint((/64,64,64/),2,TestQoI_singleStep)
!	call twoParticle_att_AdjTest(1.0_mp,(/32,32,32/),MKE,dMKE,Orbit_radius,dOrbit,dOrbit_dB)
!	call EfieldKernel
!	call LandauTraj((/64,64,64/),(/2**7,2**7,2**7/))
!	call Landau_AdjTest(1.0_mp,(/32,32,32/),(/2**6,2**6,2**6/),MPE,dMPE,IC_wave,dIC_wave,dIC_wave_dB)
!	call Compare_Traj_twostream(0.2_mp,(/32,32,32/),(/2**8,2**5,2**5/),V_Th ,EXP(-5.0_mp),EXP(-6.0_mp))
!	call Twostream_AdjTest(0.2_mp,(/32,32,32/),(/2**8,2**5,2**5/),MPE,dMPE,V_Th,dV_Th,dV_Th_dB)
!	call Compare_Traj_Landau(1.0_mp,(/32,32,32/),(/2**7,2**7,2**7/),V_Th ,EXP(-5.0_mp),EXP(-10.0_mp))
!	call Debye_shielding
!	call test_FFTPoisson
!   call Debye_characterization(Debye_length)

	call cpu_time(finish)

	print *, 'Time = ',finish-start

	! print to screen
	print *, 'program main...done.'

contains

	! You can add custom subroutines/functions here later, if you want

	subroutine Debye_shielding
		type(PM3D) :: pm
		integer, parameter :: Nd = 50, Ng(3) = (/32,32,32/)
		real(mp), parameter :: Tf = 250.0_mp, Ld = 5.0_mp, vT = 0.1_mp

		call buildDebye(pm,Tf,Ld,Nd,Ng,vT,dir_input='Debye',mod_input=5)

		open(unit=405,file='data/'//pm%r%dir//'/rho_back.bin',status='replace',form='unformatted',access='stream')
		write(405) pm%m%rho_back
		close(405)

		call forwardsweep(pm,Null_input)
		call printPlasma(pm%r)

		call destroyPM3D(pm)
	end subroutine

	subroutine Debye_characterization(QoI)
		interface
			subroutine QoI(pm,k,J)
				use modQoI
				type(PM3D), intent(in) :: pm
				real(mp), intent(inout) :: J
				integer, intent(in) :: k
			end subroutine
		end interface
		optional :: QoI
		integer :: ierr, my_rank, s
		type(PM3D) :: pm
		integer, parameter :: Nsample = 1001
		integer :: sample_per_core, sendcnt
		integer, allocatable :: recvcnt(:), displc(:)
		real(mp), allocatable :: sendbuf(:,:)
		real(mp) :: recvbuf(Nsample,2)                     !(/ A, J(A) /)
      real(mp) :: A0(Nsample), Ak
		real(mp), parameter :: Tf=200.0_mp, Ld=5.0_mp
		integer, parameter :: Ng(3)=(/32,32,32/), Np=75
		character(len=100)::dir,rank_str
		real(mp) :: J0,J1,grad(1)
		integer :: i
      A0 = (/ ( 0.5_mp*(i-1)/(Nsample-1) + 0.05_mp , i=1,Nsample ) /)

		call MPI_INIT(ierr)
		call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
		call MPI_COMM_SIZE(MPI_COMM_WORLD,s,ierr)
		sample_per_core = Nsample/s
		write(rank_str,*) my_rank

		allocate(recvcnt(0:s-1))
      allocate(displc(0:s-1))
		recvcnt(0:MOD(Nsample,s)-1) = sample_per_core+1
		recvcnt(MOD(Nsample,s):s-1) = sample_per_core
	   displc = 0
	   do i=0,s-1
		  displc(i) = SUM(recvcnt(0:i-1))
	   end do

      if( my_rank.eq.s-1 ) then
		   print *, 'size: ',s
		   print *, 'sample/core: ',sample_per_core
		   print *, 'remainder: ',MOD(Nsample,s)
		   print *, 'recvcnt: ', recvcnt
		   print *, 'displc: ', displc
      end if

		if(my_rank<MOD(Nsample,s) ) then
			sendcnt = sample_per_core+1
			allocate(sendbuf(sendcnt,2))
		else
			sendcnt = sample_per_core
			allocate(sendbuf(sendcnt,2))
		end if
		sendbuf = 0.0_mp

!		call init_random_seed(my_rank)

		do i=1,sendcnt
         Ak = A0( displc(my_rank)+i )
			dir = 'Debye_observable'//trim(adjustl(rank_str))
		   call buildDebye(pm,Tf,Ld,Np,Ng,Ak,dir_input=trim(dir),mod_input=50)

			call forwardsweep(pm,Null_input,QoI,J0)
			print *, 'vT=',Ak,', J0=',J0

			call destroyPM3D(pm)

			sendbuf(i,:) = (/ Ak, J0 /)
		end do

		do i=1,2
			call MPI_GATHERV(sendbuf(:,i),sendcnt,MPI_DOUBLE,recvbuf(:,i),recvcnt,displc,MPI_DOUBLE,s-1,MPI_COMM_WORLD,ierr)
		end do

		call MPI_FINALIZE(ierr)
		if( my_rank.eq.s-1 ) then
			open(unit=301,file='data/Debye_observable/Ak.bin',status='replace',form='unformatted',access='stream')
			open(unit=302,file='data/Debye_observable/Jk.bin',status='replace',form='unformatted',access='stream')
		   write(301) recvbuf(:,1)
		   write(302) recvbuf(:,2)
		   close(301)
		   close(302)
		end if

		deallocate(sendbuf)
		if( my_rank.eq.s-1) then
		   deallocate(recvcnt)
		   deallocate(displc)
		end if
	end subroutine


end program
