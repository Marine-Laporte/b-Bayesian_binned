SUBROUTINE sort_shell(arr)

IMPLICIT NONE
integer, DIMENSION(:), INTENT(INOUT) :: arr
!Sorts an array arr into ascending numerical order by Shell’s method (diminishing increment
!sort). arr is replaced on output by its sorted rearrangement.
INTEGER :: i,j,inc,n
REAL*8 :: v
n=size(arr)
inc=1
do !Determine the starting increment.
	inc=3*inc+1
	if (inc > n) exit
end do
do !Loop over the partial sorts.
	inc=inc/3
	do i=inc+1,n !Outer loop of straight insertion.
		v=arr(i)
		j=i
		do !Inner loop of straight insertion.
			if (arr(j-inc) <= v) exit
			arr(j)=arr(j-inc)
			j=j-inc
			if (j <= inc) exit
		end do
		arr(j)=v
	end do
	if (inc <= 1) exit
end do
END SUBROUTINE sort_shell

SUBROUTINE sort_shell_real(arr)
IMPLICIT NONE
real*8, DIMENSION(:), INTENT(INOUT) :: arr
!Sorts an array arr into ascending numerical order by Shell’s method (diminishing increment
!sort). arr is replaced on output by its sorted rearrangement.
INTEGER :: i,j,inc,n
REAL*8 :: v
n=size(arr)
inc=1
do !Determine the starting increment.
	inc=3*inc+1
	if (inc > n) exit
end do
do !Loop over the partial sorts.
	inc=inc/3
	do i=inc+1,n !Outer loop of straight insertion.
		v=arr(i)
		j=i
		do !Inner loop of straight insertion.
			if (arr(j-inc) <= v) exit
			arr(j)=arr(j-inc)
			j=j-inc
			if (j <= inc) exit
		end do
		arr(j)=v
	end do
	if (inc <= 1) exit
end do
END SUBROUTINE sort_shell_real


