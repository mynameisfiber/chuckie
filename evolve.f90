module evolve
  USE interfaces
  implicit none

  include "omp_lib.h"
  
  CONTAINS

  subroutine advance (U, Unew, n, m, o, ncomp, ng, dx, dt)
    implicit none
    integer, dimension(3)                    :: lo, hi
    integer                                  :: ncomp,ng, n,m,o, i,j,k
    double precision                         :: dx,dt
    double precision, dimension(Ncomp,n,m,o) :: U, Unew

    n = hi(1) - lo(1) + 2 * ng
    m = hi(2) - lo(2) + 2 * ng
    o = hi(3) - lo(3) + 2 * ng

    ! Shu/Osher 3rd order method:
    ! |C.-W. Shu and S. Osher, Efficient implementation of essentially
    ! | non-oscillatory shock-capturing schemes, Journal of Computational
    ! | Physics, 77 (1988), pp. 439–471.
    ! More information:
    ! http://www.amath.washington.edu/~ketch/David_Ketcheson/Publications_files/lowstorage.pdf
    Unew =         U                +         dt * dU(U   , ng, dx)
    Unew = 0.75  * U + 0.25  * Unew + 0.25  * dt * dU(Unew, ng, dx)
    Unew = 1./3. * U + 2./3. * Unew + 2./3. * dt * dU(Unew, ng, dx)
  end subroutine
end module


  function dU(U, ng, dx)
    implicit none
    double precision, dimension(:,:,:,:)              :: U
    double precision, allocatable, dimension(:,:,:,:) :: dU
    integer                                           :: n, m, o, ncomp, ng
    double precision                                  :: dx, g, trA
    double precision, dimension(6)                    :: invg, Aup, Ricci
    double precision, dimension(18)                   :: Chris, Chrislow
    integer                                           :: i,j,k,l,ijk,oi,oj,ok

    ncomp = size(U,1)
    n     = size(U,2)
    m     = size(U,3)
    o     = size(U,4)
    allocate(dU(ncomp, n, m, o))


    !$omp parallel do default(none) shared(dU,U,n,m,o,ncomp,ng,dx) &
    !$omp private(g,trA,invg,Aup,Ricci,Chris,Chrislow,i,j,k,l,ijk,oi,oj,ok)
    do oi=ng,n-ng
      do oj=ng,m-ng
        do ok=ng,o-ng
          dU(:,oi,oj,ok) = 1.

        end do
      end do
    end do

    deallocate(dU)

  end function


  function inverse_metric(g) result(ginv)
    implicit none
    double precision, dimension(6) :: g, ginv

    return
  end function


  function christoffel(g, ginv) result(Chris)
    implicit none
    double precision, dimension(6)  :: g, ginv
    double precision, dimension(18) :: Chris

    return
  end function


  function lower_symbol(g, Chris) result(Chrislow)
    implicit none
    double precision, dimension(6)  :: g
    double precision, dimension(18) :: Chris, Chrislow

    return
  end function


  function raise_matrix(ginv, Alow) result(Aup)
    implicit none
    double precision, dimension(6)  :: ginv
    double precision, dimension(18) :: Alow, Aup
    integer :: i, j, k, l

    do i=1,3
      do j=i,3
        Aup(iMat(i,j)) = 0.0
        do k=1,3
          do l=1,3
            Aup(iMat(i,j)) = Aup(iMat(i,j)) + &
                ginv(iMat(i,k)) * ginv(iMat(j,l)) * Alow(iMat(k,l))
          end do
        end do
      end do
    end do

    return
  end function


  function ricci(g, phi, Gam, Chris, Chrislow, ginv) result(R) 
    implicit none
    double precision                :: phi
    double precision, dimension(6)  :: g, ginv, Gam, R
    double precision, dimension(18) :: Chris, Chrislow

    return
  end function


  function partial(U, idx, i)
    implicit none
    double precision, dimension(:,:,:,:) :: U
    double precision                     :: partial
    integer                              :: idx, i

    return
  end function


  function partial2(U, idx, i, j)
    implicit none
    double precision, dimension(:,:,:,:) :: U
    double precision                     :: partial2
    integer                              :: idx, i, j

    return
  end function


  function DiDj(U, idx, i, j, ginv, Chris)
    implicit none
    double precision, dimension(:,:,:,:) :: U
    double precision, dimension(6)       :: ginv
    double precision, dimension(18)      :: Chris
    double precision                     :: DiDj
    integer                              :: idx, i, j

    if (i .eq. j) then
      DiDj = D2(U, idx, ginv, Chris)
      return
    end if

    return
  end function


  function D2(U, idx, ginv, Chris)
    implicit none
    double precision, dimension(:,:,:,:) :: U
    double precision, dimension(6)       :: ginv
    double precision, dimension(18)      :: Chris
    double precision                     :: D2
    integer                              :: idx

    return
  end function


  function iMat(i, j)
    implicit none
    integer :: i, j, k, iMat

    return
  end function


  function iSym(i, j)
    implicit none
    integer :: i, j, iSym

    return
  end function


  function get_ijk(i, j, k) result(ijk)
    implicit none
    integer :: i, j, k, ijk

    return
  end function
end module
