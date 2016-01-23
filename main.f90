program main

	use testmodule

	implicit none

	real(mp) :: start, finish

	! print to screen
	print *, 'calling program main'

	call cpu_time(start)

!	call setup					!A,B and QOI are initialized here
!	call particle_initialize
!	call forwardsweep(langmuir)
!	call destroyPlasma(langmuir)

	call test_FFTPoisson_adj((/ 32,32,32 /),2)

	call cpu_time(finish)

	print *, 'Time = ',finish-start

	! print to screen
	print *, 'program main...done.'

contains

	! You can add custom subroutines/functions here later, if you want

end program
