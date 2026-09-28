program comp
  !10/1/2025
  implicit none
  include 'common.h'
  
  integer :: n,i,j,IErr,lda,lwork,k,il,ir,id,iu,l
  
  real(kind=8) :: t,u,mu,z,psiold,r,rn,energy_error,V,muold
  real(kind=8):: psi_error,rho_error,psi_j_error
  real(kind=8) ::new_energy(d,d),old_energy(d,d)
  real(kind=8):: eigenvalues(dim)
  real(kind=8):: rho_real(d,d),psi_mod(d,d)
  
  complex(kind=8) ::psi_i(d,d)
  complex(kind=8) ::psi_new(d,d)
  complex(kind=8) ::rho_new(d,d),rho_i(d,d)
  complex(kind=8)::oneham(dim,dim),ground_state(dim),psi_u, psi_l, psi_d, psi_r, rho_u, rho_l, rho_d, rho_r
  
  
  
 
  V=0.0*u
  t=0.0D0
  
  
  
  do while(t<=0.05)
  
  muold = 0.0D0
  u = 1.0;  mu = 0.0D0
  
  
  
  do while (mu.lt.5.0)
     rho_i=0.2
     psi_i= (0.9,0.3)
     old_energy = 0.0
     mu = mu + 0.01
     
1235 CONTINUE
        
     do l = 1,d
	do k = 1,d
    
    !print*,l,k
    if(l==1)then
       psi_u = psi_i(d,k)
       rho_u = rho_i(d,k)
    end if
    
    if(l==d) then
       psi_d = psi_i(1,k)
       rho_d = rho_i(1,k)
    end if
    
    if(k==1) then
       psi_l = psi_i(l,d)
       rho_l = rho_i(l,d)
    end if
    
    if(k==d) then
       psi_r = psi_i(l,1)
       rho_r = rho_i(l,1)   
    endif
   
    CALL onesiteham(u,mu,t,v,psi_i(l,k),rho_i(l,k),&
         psi_l,psi_r,psi_u,psi_d,&
         rho_l,rho_r,rho_u,rho_d,&
         psi_new(l,k),rho_new(l,k),new_energy(l,k))
    
    psi_mod(l,k) = sqrt(real(psi_i(l,k))**2 + aimag(psi_i(l,k))**2)
    rho_real(l,k) = real(rho_new(l,k))
    
 end do
end do


do l =1,d
   DO k=1,D
      
      
      IF(abs(new_energy(l,k)-old_energy(l,k)) > 1.0E-10 )    THEN
         psi_i = (1.0-0.618)*psi_i + 0.618*psi_new
         rho_i = (1.0-0.618)*rho_i + 0.618*rho_new
         old_energy=new_energy
         GOTO 1235
      ENDIF
      
   ENDDO
end do

!WRITE(*,100)U,mu, psi_mod, rho_real


if (psi_mod(1,1) > 0.003)then
if ((mu-muold) > 0.1)then
print*,t,mu
print*,t,muold
end if
muold=mu
end if

end do
t = t + 0.001

end do

end program comp


!_____********************************_____________________________________________***************************************_________________________**********************************************************_____________________________________________________________________

SUBROUTINE onesiteham(u,mu,t,v,psi_m,rho_m,&
     psi_l,psi_r,psi_u,psi_d,&
     rho_l,rho_r,rho_u,rho_d,&
     newpsi,newrho,groundstateE)
  IMPLICIT NONE
  INCLUDE 'common.h'
  
  integer :: n,i,j,IErr,lda,lwork
  
  real(kind=8) :: t,u,mu,V
  real(kind=8):: eigenvalues(dim),groundstateE
  
  complex(kind=8):: m(dim,dim),ground_state(dim),newpsi,newrho
  complex(kind=8):: psi_m,rho_m,psi_u, psi_l, psi_d, psi_r, rho_u, rho_l, rho_d, rho_r
  
  !print*,psi_m,psi_u,psi_l, psi_d, psi_r
  !STOP
  m=(0.0D0,0.0D0)
  n=dim
  do i=1,n
     do j=1,n
        if (i==j)then
           m(i,j) = t*(&
                (psi_r*conjg(psi_m) + conjg(psi_r)*psi_m) / 2.0 + &
                (psi_l*conjg(psi_m) + conjg(psi_l)*psi_m) / 2.0  + &
                (psi_u*conjg(psi_m) + conjg(psi_u)*psi_m) / 2.0  + &
                (psi_d*conjg(psi_m) + conjg(psi_d)*psi_m) / 2.0 ) + &

                V*(dble(j-1)*rho_r - (rho_m*rho_r)/2.0 + &
                dble(j-1)*rho_l - (rho_m*rho_l)/2.0 + & 
                dble(j-1)*rho_u - (rho_m*rho_u)/2.0 + &
                dble(j-1)*rho_d - (rho_m*rho_d)/2.0) + &
                
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
  
  call ceigensolver(m,eigenvalues)
  
  ground_state = m(:, 1)
  groundstateE=eigenvalues(1)
      
  newpsi = (0.0,0.0)
  do i = 1, n
     do j = 1, n
        if (i == j - 1) then
           newpsi = newpsi + conjg(ground_state(i)) * ground_state(j) * dsqrt(dble(j - 1))
        end if
     end do
  end do
  
  newrho = (0.0, 0.0)
  do i = 1, n
     newrho = newrho + conjg(ground_state(i)) * ground_state(i) * dble(i - 1)
  end do
  
  RETURN
END SUBROUTINE onesiteham
    
     
     
 !*************************COMPLEX***EIGENSOLVER***************************************     
  SUBROUTINE ceigensolver(A,W)
  IMPLICIT NONE

  INCLUDE 'common.h'

  INTEGER :: lwork=3*dim,INFO,lda=dim,N=dim,i,ind,IL,IU,M,LDZ,ISUPPZ(2),LRWORK,LIWORK
  COMPLEX(kind=8) :: A(dim,dim),Work(3*dim),Z(dim,dim)
  REAL(KIND=8) :: W(dim),test,VL,VU,RWORK(24*dim),IWORK(10*dim)
  CHARACTER :: jobz,uplo
  REAL(KIND=8) :: ABSTOL

  
  jobz ='V'       ! Compute eigenvalues and eigenvectors
  uplo ='L'       ! Lower triangle of A is stored
  IL=1
  IU=dim
  M=IU-IL+1
  LDZ=dim
  LRWORK=24*dim
  LIWORK=10*dim
  
  CALL zheev(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,RWORK,INFO)
  IF(INFO.NE.0)PRINT*,"Error in ZHEEV ",Info
  
 RETURN
END SUBROUTINE ceigensolver

