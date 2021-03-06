module interfaces
  implicit none

  interface
    function dU(U, ng, dx)
      implicit none
      double precision, dimension(:,:,:,:)              :: U
      double precision, allocatable, dimension(:,:,:,:) :: dU
      double precision                                  :: dx
      integer                                           :: ng
    end function

    function partial(U, idx, i)
      implicit none
      double precision, dimension(:,:,:,:) :: U
      double precision                     :: partial
      integer                              :: idx, i
    end function

    function partial2(U, idx, i, j)
      implicit none
      double precision, dimension(:,:,:,:) :: U
      double precision                     :: partial2
      integer                              :: idx, i, j
    end function

    function DiDj(U, idx, i, j, ginv, Chris)
      implicit none
      double precision, dimension(:,:,:,:) :: U
      double precision, dimension(6)       :: ginv
      double precision, dimension(18)      :: Chris
      double precision                     :: DiDj
      integer                              :: idx, i, j
    end function

    function D2(U, idx, ginv, Chris)
      implicit none
      double precision, dimension(:,:,:,:) :: U
      double precision, dimension(6)       :: ginv
      double precision, dimension(18)      :: Chris
      double precision                     :: D2
      integer                              :: idx
    end function
    
  end interface
end module
