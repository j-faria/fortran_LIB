! Sorting related routines
! João Faria (c) 2013
!   Based on code written by John E. Pask, LLNL.

module lib_sort
    
    implicit none
    save
    
    private
    
    integer,parameter :: dp = selected_real_kind(p=15,r=307)

    public sort, sortpairs, argsort

    ! overload argsort
    interface argsort
        module procedure iargsort, &
                         rargsort
    end interface

    ! overload sort
    interface sort
        module procedure sortNums, &
                         sortINums, &
                         sortVecs
    end interface

    ! overload sortpairs
    interface sortpairs
        module procedure sortNumNumPairs, &
                         sortINumCNumPairs, &
                         sortNumVecPairs, &
                         sortINumVecPairs, &
                         sortNumCVecPairs, &
                         sortNumMatPairs, &
                         sortINumMatPairs
    end interface

contains

    subroutine sortNums(array)
    ! sorts array of reals, from smallest to largest
        real(dp), intent(inout):: array(:) ! array of reals
        array = array(argsort(array))
    end subroutine sortNums

    subroutine sortINums(array)
    ! sorts array of integers, from smallest to largest.
        integer, intent(inout):: array(:) ! array of integers
        array = array(argsort(array))
    end subroutine sortINums

subroutine sortVecs(vecs)
! sorts array of vectors, vecs, by length, from smallest to largest
real(dp), intent(inout):: vecs(:,:) ! array of vectors: vecs(i,j) = ith comp. of jth vec.
real(dp) len2(size(vecs,2)) ! array of squares of vector lengths
integer i
do i=1,size(len2)
   len2(i)=dot_product(vecs(:,i),vecs(:,i))
end do
call sortpairs(len2,vecs)
end subroutine

    subroutine sortNumNumPairs(array1, array2)
    ! sorts arrays of reals, array1 and array2, according to increasing array1
        real(dp), intent(inout):: array1(:), array2(:) ! arrays of reals
        integer :: a(size(array1))
        if (size(array1) /= size(array2)) then
        call error('SORTPAIRS', 'arrays must be of same length.')
        end if
        a = argsort(array1)
        array1 = array1(a)
        array2 = array2(a)
    end subroutine sortNumNumPairs

    subroutine sortINumCNumPairs(array1, array2)
    ! sorts arrays of integers, array1, and complex, array2, according to increasing array1
        integer, intent(inout):: array1(:) ! array of integers
        complex(dp), intent(inout):: array2(:) ! array of complex numbers
        integer :: a(size(array1))
        if (size(array1) /= size(array2)) then
        call error('SORTPAIRS', 'arrays must be of same length.')
        end if
        a = argsort(array1)
        array1 = array1(a)
        array2 = array2(a)
    end subroutine sortINumCNumPairs

subroutine sortNumVecPairs(nums, vecs)
! sorts arrays of numbers, nums, and vectors, vecs, according to increasing nums
real(dp), intent(inout):: nums(:) ! array of numbers
real(dp), intent(inout):: vecs(:,:) ! array of vectors: vecs(i,j) = ith comp. of jth vec.
integer :: a(size(nums))
if (size(nums) /= size(vecs,2)) then
call error('SORTPAIRS', 'arrays must be of same length.')
end if
a = argsort(nums)
nums = nums(a)
vecs = vecs(:, a)
end subroutine

subroutine sortINumVecPairs(nums, vecs)
! sorts arrays of integers, nums, and vectors, vecs, according to increasing nums
integer, intent(inout):: nums(:) ! array of numbers
real(dp), intent(inout):: vecs(:,:) ! array of vectors: vecs(i,j) = ith comp. of jth vec.
integer :: a(size(nums))
if (size(nums) /= size(vecs,2)) then
call stop_error("SORTPAIRS ERROR: arrays must be of same length.")
end if
a = argsort(nums)
nums = nums(a)
vecs = vecs(:, a)
end subroutine

subroutine sortNumCVecPairs(nums, vecs)
! sorts arrays of numbers, nums, and complex vectors, vecs, according to increasing nums
real(dp), intent(inout):: nums(:) ! array of numbers
complex(dp), intent(inout):: vecs(:,:) ! array of vectors: vecs(i,j) = ith comp. of jth vec.
integer :: a(size(nums))
if (size(nums) /= size(vecs,2)) then
call stop_error("SORTPAIRS ERROR: arrays must be of same length.")
end if
a = argsort(nums)
nums = nums(a)
vecs = vecs(:, a)
end subroutine

subroutine sortNumMatPairs(nums, mats)
! sorts arrays of numbers, nums, and matrices, mats, according to increasing nums
real(dp), intent(inout):: nums(:) ! array of numbers
real(dp), intent(inout):: mats(:,:,:) ! array of matrices: mats(i,j,n) = (i,j) comp. of nth mat.
integer :: a(size(nums))
if (size(nums) /= size(mats,3)) then
call stop_error("SORTPAIRS ERROR: arrays must be of same length.")
end if
a = argsort(nums)
nums = nums(a)
mats = mats(:, :, a)
end subroutine

subroutine sortINumMatPairs(nums,mats)
! sorts arrays of integers, nums, and matrices, mats, according to increasing nums
integer, intent(inout):: nums(:) ! array of numbers
real(dp), intent(inout):: mats(:,:,:) ! array of matrices: mats(i,j,n) = (i,j) comp. of nth mat.
integer :: a(size(nums))
if (size(nums) /= size(mats,3)) then
call stop_error("SORTPAIRS ERROR: arrays must be of same length.")
end if
a = argsort(nums)
nums = nums(a)
mats = mats(:, :, a)
end subroutine

    function iargsort(a) result(b)
    ! Returns the indices that would sort an array of integers a.
    !
    ! Arguments
    ! ---------
    !
        integer, intent(in):: a(:) ! array of integers
        integer :: b(size(a)) ! indices into the array 'a' that sort it
    !
    ! Example
    ! -------
    !
    ! iargsort([10, 9, 8, 7, 6]) ! Returns [5, 4, 3, 2, 1]
    integer :: N ! number of numbers/vectors
        integer :: i,imin ! indices: i, i of smallest
        integer :: temp ! temporary
        integer :: a2(size(a))
        a2 = a
        N=size(a)
        do i = 1, N
            b(i) = i
        end do
        do i = 1, N-1
            ! find ith smallest in 'a'
            imin = minloc(a2(i:),1) + i - 1

            ! swap to position i in 'a' and 'b', if not already there
            if (imin /= i) then
                temp = a2(i); 
                a2(i) = a2(imin); 
                a2(imin) = temp
                temp = b(i); 
                b(i) = b(imin); 
                b(imin) = temp
            end if
        end do
    end function iargsort

    function rargsort(a) result(b)
    ! Returns the indices that would sort an array of reals.
    !
    ! Arguments
    ! ---------
    !
        real(dp), intent(in):: a(:) ! array of reals
        integer :: b(size(a)) ! indices into the array 'a' that sort it
    !
    ! Example
    ! -------
    !
    ! rargsort([4.1_dp, 2.1_dp, 2.05_dp, -1.5_dp, 4.2_dp]) ! Returns [4, 3, 2, 1, 5]
        integer :: N ! number of numbers/vectors
        integer :: i,imin ! indices: i, i of smallest
        integer :: temp1 ! temporary
        real(dp) :: temp2
        real(dp) :: a2(size(a))
        a2 = a
        N=size(a)
        do i = 1, N
            b(i) = i
        end do
        do i = 1, N-1
            ! find ith smallest in 'a'
            imin = minloc(a2(i:),1) + i - 1
            ! swap to position i in 'a' and 'b', if not already there
            if (imin /= i) then
                temp2 = a2(i); 
                a2(i) = a2(imin); 
                a2(imin) = temp2
                temp1 = b(i); 
                b(i) = b(imin); 
                b(imin) = temp1
            end if
        end do
    end function rargsort

end module
