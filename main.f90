program main

	use testmodule
	use twoparticle

	implicit none

	real(mp) :: start, finish

	! print to screen
	print *, 'calling program main'

	call cpu_time(start)

	call twoParticleTest(0.2_mp,(/64,64,64/))
!	call twostream(0.2_mp,(/64,64,64/),(/2**8,2**5,2**5/))

	call cpu_time(finish)

	print *, 'Time = ',finish-start

	! print to screen
	print *, 'program main...done.'

contains

	! You can add custom subroutines/functions here later, if you want

end program
