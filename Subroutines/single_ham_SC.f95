SUBROUTINE onesiteham(u,mu,t,psi_m,&
     psi_l,psi_r,psi_u,psi_d,&
     newpsi,GS_Eval,GS)
  IMPLICIT NONE
  INCLUDE 'common.h'
  
  integer :: n,i,j,IErr,lda,lwork
  
  real(kind=8) :: t,u,mu
  real(kind=8):: eigenvalues(dim),GS_Eval
  
  complex(kind=8):: m(dim,dim),GS(dim),newpsi
  complex(kind=8):: psi_m,psi_u, psi_l, psi_d, psi_r
  
  !print*,psi_m,psi_u,psi_l, psi_d, psi_r
  !STOP

  CALL ham_i(u,mu,t,psi_m,&
     psi_l,psi_r,psi_u,psi_d,&
     m)
  
  call ceigensolver(m,eigenvalues)
  
  GS = m(:, 1)
  GS_Eval=eigenvalues(1)
  
  newpsi = (0.0,0.0)
  do i = 1, nmax+1
     do j = 1, nmax+1
        if (i == j - 1) then
           newpsi = newpsi + conjg(GS(i)) * GS(j) * dsqrt(dble(j - 1))
        end if
     end do
  end do
  
  
  
  RETURN
END SUBROUTINE onesiteham
