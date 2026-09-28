SUBROUTINE RK5(TQ,time,t_spacing,t_i,&
     u,mu,Jth,J_i,Jc,&
     psi_i,GS_mat,GS_matN)
  
  IMPLICIT NONE
  include 'common.h' ! This includes following
  
  integer :: l,k ! To loop over (l,k)  site positions
  integer :: i,j ! To loop over bra,kets of single site states 
  
  real(kind=8) :: u,mu,Jth ! On-site Repulsion,Chemical Potential, Hopping strength
  
  complex(kind=8) ::psi_i(d,d),psi_new(d,d)  ! All initial and new psi values of sites (k,l)
  complex(kind=8)::psi_u, psi_l, psi_d, psi_r ! psi of : top, left, right, down, and right

  real(kind=8) :: GS_energy(d,d) ! Ground state energy of sites (k,l)
  complex(kind=8)::GS_mat(d,d,dim)! Ground state vector of sites (k,l)
  
  real(kind=8):: rho(d,d),psi_mod(d,d)   ! densities and |psi| of sites (k,l)

  ! Declaration for dynamic calculations
  integer :: t_panels,t_ind  ! No.of Panels for time and time index
  real(kind=8) :: J_i,Jc,J_f ! Initial J,Critical J, and final J (hopping ampl.)
  real(kind=8) :: t_i,t_f,time ! Initial and Final and  time
  real(kind=8) :: t_spacing,TQ ! Time increment (delat t) and Quench time
  real(kind=8) :: time1,time2,time3,time4,time5,time6

  real(kind=8) :: normsum ! To do normalization

  real(kind=8) ::JK1,JK2,JK3,JK4,JK5,JK6 !The hop.amp. corresponding to time for 4th order RK 
  complex(kind=8) :: K1_dGS(dim),K2_dGS(dim),K3_dGS(dim),K4_dGS(dim) ! Corresponding GS coeffients
  complex(kind=8) :: K5_dGS(dim),K6_dGS(dim)
  
  complex(kind=8) :: dham(dim,dim) ! one site hamiltonian mat. to do dynamcics
  complex(kind=8) :: dGS(dim),ndGS(dim),mul(dim) ! The old and new dynamical GS
  complex(kind=8)::GS_matN(d,d,dim)!
  !PRINT*,J_i,Jc
  !PRINT*,"time",Jth,time1,time2,time3,time4,time5
  !PRINT*,t_spacing,Jth,time1,time2,time3,time4,time5
  
  time1=time
  JK1=J_i+(Jc-J_i)*(time1+TQ)/TQ

  time2=time+t_spacing/5.0D0
  JK2=J_i+(Jc-J_i)*(time2+TQ)/TQ
  
  time3=time+3.0D0*t_spacing/10.0D0
  JK3=J_i+(Jc-J_i)*(time3+TQ)/TQ

  time4=time+3.0D0*t_spacing/5.0D0
  JK4=J_i+(Jc-J_i)*(time4+TQ)/TQ
  
  time5=time+t_spacing
  JK5=J_i+(Jc-J_i)*(time5+TQ)/TQ
  
  time6=time+7.0D0*t_spacing/8.0D0
  JK6=J_i+(Jc-J_i)*(time6+TQ)/TQ


  

  
  do l=1,d ! Do RK calculation for each site
     do k=1,d
        
        dGS=GS_mat(l,k,:) !The Previous step Ground State
        
        ! Periodic Boundary conditions 
        if(l==1)then
           psi_u = psi_i(d,k)
        else
           psi_u = psi_i(l-1,k)
        end if
        
        if(l==d) then
           psi_d = psi_i(1,k)
        else
           psi_d = psi_i(l+1,k)
        end if
        
        if(k==1) then
           psi_l = psi_i(l,d)
        else
           psi_l = psi_i(l,k-1)
        end if
        
        if(k==d) then
           psi_r = psi_i(l,1)
        else
           psi_r = psi_i(l,k+1)
        endif
        
      
        CALL ham_i(u,mu,JK1,psi_i(l,k),&
             psi_l,psi_r,psi_u,psi_d,&
             dham)
        K1_dGS=-COMPLEX(0.0D0,1.0D0)*matmul(dham,dGS)
        
        CALL ham_i(u,mu,JK2,psi_i(l,k),&
             psi_l,psi_r,psi_u,psi_d,&
             dham)
        K2_dGS=-COMPLEX(0.0D0,1.0D0)*matmul(dham,dGS&
             +K1_dGS*t_spacing/5.0D0)
           
        CALL ham_i(u,mu,JK3,psi_i(l,k),&
             psi_l,psi_r,psi_u,psi_d,&
             dham)
        K3_dGS=-COMPLEX(0.0D0,1.0D0)*matmul(dham,dGS&
             +3.0D0*K1_dGS*t_spacing/40.0D0&
             +9.0D0*K2_dGS*t_spacing/40.0D0)
           
        CALL ham_i(u,mu,JK4,psi_i(l,k),&
             psi_l,psi_r,psi_u,psi_d,&
             dham)
        K4_dGS=-COMPLEX(0.0D0,1.0D0)*matmul(dham,dGS&
             +3.0D0*K1_dGS*t_spacing/10.0D0&
             -9.0D0*K2_dGS*t_spacing/10.0D0&
             +6.0D0*K3_dGS*t_spacing/5.0D0)

        CALL ham_i(u,mu,JK5,psi_i(l,k),&
             psi_l,psi_r,psi_u,psi_d,&
             dham)
        K5_dGS=-COMPLEX(0.0D0,1.0D0)*matmul(dham,dGS&
             -11.0D0*K1_dGS*t_spacing/54.0D0&
             +5.0D0*K2_dGS*t_spacing/2.0D0&
             -70.0D0*K3_dGS*t_spacing/27.0D0&
             +35.0D0*K4_dGS*t_spacing/27.0D0)

        CALL ham_i(u,mu,JK6,psi_i(l,k),&
             psi_l,psi_r,psi_u,psi_d,&
             dham)
        K6_dGS=-COMPLEX(0.0D0,1.0D0)*matmul(dham,dGS&
             +16131.0D0*K1_dGS*t_spacing/55296.0D0&
             +175.0D0*K2_dGS*t_spacing/512.0D0&
             +575.0D0*K3_dGS*t_spacing/13824.0D0&
             +44275.0D0*K4_dGS*t_spacing/110592.0D0&
             +253.0D0*K5_dGS*t_spacing/4096.0D0)
        

           ! Get the new GS using the 4 terms
        ndGS=dGS+t_spacing*(&
             2825.0D0*K1_dGS/27648.0D0&
             +18575.0D0*K3_dGS/48384.0D0&
             +13525.0D0*K4_dGS/55296.0D0&
             +277.0D0*K5_dGS/14336.0D0&
              +1.0D0*K6_dGS/4.0D0)
        
         !PRINT*,ndGS
           !  Normalized the new GS
           normsum=0.0D0
           do i=1,dim
              normsum=normsum+ndGS(i)*CONJG(ndGS(i))
           end do
           ndGS=ndGS/DSQRT(normsum)
           !PRINT*,normsum,ndGS*CONJG(ndGS)
           GS_matN(l,k,:)=ndGS ! update the previous GS by new GS
        enddo
     end do
  RETURN
END SUBROUTINE RK5
