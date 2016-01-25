program main

	use testmodule

	implicit none

	real(mp) :: start, finish

	! print to screen
	print *, 'calling program main'

	call cpu_time(start)

	call twostream(0.2_mp,(/ 64,64,64 /),(/ 2**7,2**4,2**4 /))
!	call verify_assignment
!	call test_FFTPoisson_adj((/ 32,32,32 /),2)

	call cpu_time(finish)

	print *, 'Time = ',finish-start

	! print to screen
	print *, 'program main...done.'

contains

	! You can add custom subroutines/functions here later, if you want

end program
