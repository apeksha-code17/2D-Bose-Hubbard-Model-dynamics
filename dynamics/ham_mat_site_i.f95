SUBROUTINE ham_i(u,mu,t,psi_m,&
     psi_l,psi_r,psi_u,psi_d,&
     m)
  IMPLICIT NONE
  INCLUDE 'common.h'
  
  integer :: n,i,j,IErr,lda,lwork
  
  real(kind=8) :: t,u,mu
  real(kind=8):: eigenvalues(dim),groundstateE
  
  complex(kind=8):: m(dim,dim),ground_state(dim),newpsi
  complex(kind=8):: psi_m,psi_u, psi_l, psi_d, psi_r
  
  !print*,psi_m,psi_u,psi_l, psi_d, psi_r
  !STOP
  m=(0.0D0,0.0D0)

  do i=1,nmax+1
     do j=1,nmax+1
        if (i==j)then
           m(i,j) = t*(&
                (psi_r*conjg(psi_m) + conjg(psi_r)*psi_m) / 2.0 + &
                (psi_l*conjg(psi_m) + conjg(psi_l)*psi_m) / 2.0  + &
                (psi_u*conjg(psi_m) + conjg(psi_u)*psi_m) / 2.0  + &
                (psi_d*conjg(psi_m) + conjg(psi_d)*psi_m) / 2.0 ) + &
                ((u / 2.0) * dble((j-1)*(j-2))) - mu * dble(j-1)
           
        else if (i==(j+1))then
           m(i,j)=-(t*(psi_l + psi_r + psi_u + psi_d)*dsqrt(dble(j)))
           
        else if (i==(j-1))then
           m(i,j)=-(t*(conjg(psi_r) + conjg(psi_l)+ conjg(psi_u) + conjg(psi_d))*dsqrt(dble(j-1)))
           
        else
           m(i,j)=(0.0,0.0)
           
        end if
     end do
  end do
  RETURN
END SUBROUTINE ham_i
