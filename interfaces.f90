! use interfaces, only: dU, partial, partial2, DiDj, D2
module interfaces
  implicit none
  double precision :: dx
  integer, dimension(4) :: us !the shape of U

  interface
    function dU(U, ng, dx)
      implicit none
      double precision, dimension(:,:,:,:)              :: U
      double precision, allocatable, dimension(:,:,:,:) :: dU
      double precision                                  :: dx
      integer                                           :: ng
    end function

    function partial(U, idx, i, dx)
      implicit none
      double precision, dimension(:,:,:,:) :: U
      double precision                     :: partial, dx
      integer                              :: idx, i
    end function

    function partial2(U, idx, i, j, dx)
      implicit none
      double precision, dimension(:,:,:,:) :: U
      double precision                     :: partial2, dx
      integer                              :: idx, i, j
    end function

    function DiDj(U, idx, i, j, ginv, Chris, dx)
      implicit none
      double precision, dimension(:,:,:,:) :: U
      double precision, dimension(6)       :: ginv
      double precision, dimension(18)      :: Chris
      double precision                     :: DiDj, dx
      integer                              :: idx, i, j
    end function

    function D2(U, idx, ginv, Chris, dx)
      implicit none
      double precision, dimension(:,:,:,:) :: U
      double precision, dimension(6)       :: ginv
      double precision, dimension(18)      :: Chris
      double precision                     :: D2, dx
      integer                              :: idx
    end function
    
  end interface
end module
