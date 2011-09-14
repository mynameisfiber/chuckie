
#define iMat(i, j) ( (i) + (j) + MIN(1,(i),(j)) ) 
#define iSym(i, j) ( (i) + (j)*6 )

module evolve
  USE interfaces
  implicit none

  include "omp_lib.h"
  
  CONTAINS

  subroutine advancechombo()
    ! create this function to take in the chombo style data definition (ie:
    ! double precision    U(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,Ncomp)), calc
    ! n, m, o and send to the normal advance function.
    !
    ! this function should _simply_ translate the indicies
  end subroutine

  subroutine advance (U, Unew, n, m, o, ncomp, ng, dx, dt)
    !f2py double precision, dimension(Ncomp, n, m, o), intent(in)  :: U
    !f2py double precision, dimension(Ncomp, n, m, o), intent(out) :: Unew
    !f2py integer, intent(in)                                      :: n, m, o, ncomp, ng
    !f2py double precision, intent(in)                             :: dx, dt
    implicit none
    integer                                  :: ncomp,ng, n,m,o
    double precision                         :: dx,dt
    double precision, dimension(Ncomp,n,m,o) :: U, Unew

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
    use interfaces
    implicit none
    double precision, dimension(:,:,:,:),intent(in)               :: U
    double precision, allocatable, dimension(:,:,:,:),intent(out) :: dU
    integer,intent(in)                                            :: n, m, o, ncomp, ng
    double precision, intent(in)                                  :: dx, g, trA
    double precision, dimension(6),intent(in)                     :: invg, Aup, Ricci
    double precision, dimension(18),intent(in)                    :: Chris, Chrislow
    integer,intent(in)                                            :: i,j,k,l,ijk,oi,oj,ok

    ncomp = size(U,1)
    n     = size(U,2)
    m     = size(U,3)
    o     = size(U,4)

    if (allocated(dU)) then
      deallocate(dU)
    endif
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
    use interfaces, only: D2
    implicit none
    double precision, dimension(:,:,:,:),intent(in) :: U
    double precision, dimension(6),intent(in)       :: ginv
    double precision, dimension(18),intent(in)      :: Chris
    double precision,intent(out)                    :: DiDj
    integer, intent(in)                             :: idx, i, j

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


  function get_ijk(i, j, k) result(ijk)
    implicit none
    integer :: i, j, k, ijk

    return
  end function
