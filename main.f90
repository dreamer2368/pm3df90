program main

	use testmodule
	use twoparticle
	use Landau

	implicit none

	real(mp) :: start, finish

	! print to screen
	print *, 'calling program main'

	call cpu_time(start)

!	call twoParticle_rep_AdjTest(0.4_mp,(/64,64,64/),MKE,dMKE,IC_wave,dIC_wave,dIC_wave_dB)
!	call twoParticleTraj(1.0_mp,(/32,32,32/))
!	call twostream(0.2_mp,(/64,64,64/),(/2**8,2**5,2**5/))
!	call test_fullAdjoint(0.2_mp,(/32,32,32/),(/2**8,2**5,2**5/),MPE,dMPE,IC_wave,dIC_wave,dIC_wave_dB)
!	call test_particle_adj2((/64,64,64/),2,TestQoI_singleStep)
!	call twoParticle_att_AdjTest(1.0_mp,(/32,32,32/),MKE,dMKE,Orbit_radius,dOrbit,dOrbit_dB)
!	call EfieldKernel
!	call LandauTraj((/64,64,64/),(/2**7,2**7,2**7/))
!	call Landau_AdjTest(1.0_mp,(/32,32,32/),(/2**6,2**6,2**6/),MPE,dMPE,IC_wave,dIC_wave,dIC_wave_dB)
!	call Compare_Traj_twostream(0.2_mp,(/32,32,32/),(/2**8,2**5,2**5/),V_Th ,EXP(-5.0_mp),EXP(-6.0_mp))
	call Twostream_AdjTest(0.2_mp,(/32,32,32/),(/2**8,2**5,2**5/),MPE,dMPE,V_Th,dV_Th,dV_Th_dB)

	call cpu_time(finish)

	print *, 'Time = ',finish-start

	! print to screen
	print *, 'program main...done.'

contains

	! You can add custom subroutines/functions here later, if you want

end program
