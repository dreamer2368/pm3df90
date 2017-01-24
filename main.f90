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
!	call test_fullAdjoint(0.2_mp,(/32,32,32/),(/2**8,2**5,2**5/),MKE,dMKE,V_Th,dV_Th,dV_Th_dB)
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
	call Debye_shielding
!	call test_FFTPoisson

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

		call buildDebye(pm,Tf,Ld,Nd,Ng,vT,mod_input=5)

		open(unit=405,file='data/'//pm%r%dir//'/rho_back.bin',status='replace',form='unformatted',access='stream')
		write(405) pm%m%rho_back
		close(405)

		call forwardsweep(pm,Null_input)
		call printPlasma(pm%r)

		call destroyPM3D(pm)
	end subroutine

end program
