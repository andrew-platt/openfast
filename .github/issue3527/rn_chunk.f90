! Does RANDOM_NUMBER give the same stream when an array is filled in pieces?
program rn_chunk
   implicit none
   integer, parameter :: n = 5000000
   real(8), allocatable :: a(:), b(:)
   real(4), allocatable :: a4(:), b4(:)
   integer, allocatable :: seed(:)
   integer :: ns, i, i0, m, chunks(5) = (/ 1, 7, 1009, 36000, 262144 /), k
   call random_seed(size=ns); allocate(seed(ns)); seed(1::2) = 2318573; seed(2::2) = 239492
   allocate(a(n), b(n), a4(n), b4(n))
   call random_seed(put=seed); call random_number(a)
   call random_seed(put=seed); call random_number(a4)
   do k = 1,size(chunks)
      call random_seed(put=seed)
      do i0 = 0, n-1, chunks(k)
         m = min(chunks(k), n-i0); call random_number(b(i0+1:i0+m))
      end do
      print '(A,I7,A,I0)', ' real(8) chunk ', chunks(k), ': mismatches = ', count(a /= b)
      call random_seed(put=seed)
      do i0 = 0, n-1, chunks(k)
         m = min(chunks(k), n-i0); call random_number(b4(i0+1:i0+m))
      end do
      print '(A,I7,A,I0)', ' real(4) chunk ', chunks(k), ': mismatches = ', count(a4 /= b4)
   end do
end program
