program comp
  !6/2/2025
  ! Fourth order RK
  implicit none
  
  include 'common.h' ! This includes following
  ! d = number of sites
  ! nmax= max. no. of bosons at a site
  ! dim = nmax+1  dimension of MF Hamiltonian matrix at a site
  
  integer :: l,k ! To loop over (l,k)  site positions
  integer :: i,j ! To loop over bra,kets of single site states 
  
  real(kind=8) :: u,mu,Jth ! On-site Repulsion,Chemical Potential, Hopping strength
  
  real(kind=8):: rnum,inum,alpha
  
  complex(kind=8) ::psi_i(d,d),psi_new(d,d)  ! All initial and new psi values of sites (k,l)
  complex(kind=8)::psi_u, psi_l, psi_d, psi_r ! psi of : top, left, right, down, and right

  real(kind=8) :: GS_energy(d,d) ! Ground state energy of sites (k,l)
  complex(kind=8)::GS_mat(d,d,dim)! Ground state vector of sites (k,l)
  
  
  real(kind=8):: TOT_Eng_old,TOT_Eng_new ! GS Energies per sites for minimization
  real(kind=8):: rho(d,d),psi_mod(d,d)   ! densities and |psi| of sites (k,l)


  ! Declaration for dynamic calculations
  integer :: t_panels,t_ind  ! No.of Panels for time and time index
  real(kind=8) :: J_i,Jc,J_f ! Initial J,Critical J, and final J (hopping ampl.)
  real(kind=8) :: t_i,t_f,time ! Initial and Final and  time
  real(kind=8) :: t_spacing,TQ ! Time increment (delat t) and Quench time

  real(kind=8) :: normsum,delta,psi_modAvg,rhoAvg ! To do normalization

  real(kind=8) ::JK1,JK2,JK3,JK4,JK5,JK6 !The hop.amp. corresponding to time for 4th order RK 
  complex(kind=8) :: K1_dGS(dim),K2_dGS(dim),K3_dGS(dim),K4_dGS(dim) ! Corresponding GS coeffients
  complex(kind=8) :: K5_dGS(dim),K6_dGS(dim)
  
  complex(kind=8) :: dham(dim,dim) ! one site hamiltonian mat. to do dynamcics
  complex(kind=8) :: dGS(dim),ndGS(dim),mul(dim) ! The old and new dynamical GS

  complex(kind=8)::GS_matN(d,d,dim)

  integer :: noofwantedfiles,wi
  real(kind=8) ::Jw(100),temp_value
  real(kind=8) ::tw(100)
  character(len=1024) :: filename

   
  
  u=1.0
  mu=0.5
  
  alpha=0.6 ! Factor to minimize
  
  J_i=0.030D0
  Jth=J_i

  J_f=0.10D0

  TQ=50
  Jc=0.042

  t_panels=500000

  ! Get the minimized GS for initial J vlaue
  !___________________________________________________________________-
  ! Initialize all psis to random values
  do l = 1,d
     do k = 1,d
        CALL random_number(rnum)
        CALL random_number(inum)
        psi_i(l,k)= COMPLEX(rnum,inum)
     end do
  end do
  
  TOT_Eng_old=100.0D0 ! Just for first iteration
  
1235 CONTINUE ! for minimizing energy
  
  do l = 1,d ! index of sites in 2D configuration
     do k = 1,d

        ! Periodic boundary conditions
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

        ! Calculate the one site Hamiltonian and diagonalize to get
        ! GS Energy and GS vectors 
        ! Psi value for the site (l,k)
        CALL onesiteham(u,mu,Jth,psi_i(l,k),&
             psi_l,psi_r,psi_u,psi_d,&
             psi_new(l,k),GS_energy(l,k),GS_mat(l,k,:))

        ! Calculate Rho
        rho(l,k) =0.0D0
        do i = 1, nmax+1
           rho(l,k) = rho(l,k) + conjg(GS_mat(l,k,i)) * GS_mat(l,k,i) * dble(i - 1)
        end do
        ! |psi|
        psi_mod(l,k) = sqrt(real(psi_new(l,k))**2 + aimag(psi_new(l,k))**2)  
     end do
  end do


  ! Calculate GS Energy per site
  TOT_Eng_new=0.0D0
  DO l =1,d
     DO k=1,d
        TOT_Eng_new=TOT_Eng_new+GS_energy(l,k)
     END DO
  END do
  TOT_Eng_new=1000*TOT_Eng_new/DBLE(d**2)

  ! Compare with older and do the Self consitent minimization interms of psi
  IF(abs(TOT_Eng_new-TOT_Eng_old) > 1.0E-12)THEN
     psi_i = (1.0-alpha)*psi_i+alpha*psi_new
     TOT_Eng_old=TOT_Eng_new
     GOTO 1235
  ENDIF

  !OPEN(UNIT=193,FILE='psi_real',status="unknown", position="append", action="write")
  !DO l=1,d
  !   WRITE(193,*)REAL(psi_new(l,:))
  !END DO
  !CLOSE(193)

  !OPEN(UNIT=194,FILE='psi_imag',status="unknown", position="append", action="write")
  !DO l=1,d
  !   WRITE(194,*)AIMAG(psi_new(l,:))
  !END DO
  !CLOSE(194)

  !OPEN(UNIT=195,FILE='rho',status="unknown", position="append", action="write")
  !DO l=1,d
  !   WRITE(195,*)AIMAG(psi_new(l,:))
  !END DO
  !CLOSE(195)
    

  !STOP

  
 
  !______________________________________________________________________
  !_____________ Start the time dependent part_____________      
  t_i=-TQ ! Initial scaled time 
  t_f=TQ*(-1.0D0+(J_f-J_i)/(Jc-J_i)) ! Final scaled time 
 
  t_spacing=(t_f-t_i)/DBLE(t_panels) ! time spacing 
  time=t_i

  temp_value=J_i
  wi=1
  DO WHILE (temp_value.LE.J_f+0.001)
     Jw(wi)=temp_value
     tw(wi)=TQ*(-1.0D0+(Jw(wi)-J_i)/(Jc-J_i))
     temp_value=temp_value+0.0025
     !PRINT*,wi,Jw(wi),tw(wi)
     wi=wi+1
  END DO
  noofwantedfiles=wi-1
  !PRINT*,noofwantedfiles
  !STOP
  

  
  delta=0.000001
  do l = 1,d
     do k = 1,d
        CALL random_number(rnum)
        CALL random_number(inum)
        do i = 1, nmax+1
           GS_mat(l,k,i)=GS_mat(l,k,i)+delta*COMPLEX(rnum,inum)
        end do
     end do
  end do
  
  
  !DO t_ind=1,t_panels ! start the dynamcis  
  DO WHILE(time.LE.t_f)       
     ! First RK term has this time and corresponding J
     Jth=J_i+(Jc-J_i)*(time+TQ)/TQ
     psi_modAvg=0.0D0
     rhoAvg=0.0D0
9876 CONTINUE
     CALL RK5(TQ,time,t_spacing,t_i,&
          u,mu,Jth,J_i,Jc,&
          psi_i,GS_mat,GS_matN)

     ! Calculate the new psi and rho at this time for all sites
     do l=1,d
        do k=1,d         
           
           dGS=GS_matN(l,k,:)
           psi_new(l,k)=COMPLEX(0.0D0,0.0D0)
           rho(l,k) =0.0D0
           do i=1,dim
              rho(l,k) = rho(l,k) + conjg(dGS(i))*dGS(i)*dble(i - 1)
              do j=1,dim
                 if (i==j-1) then
                    psi_new(l,k)=psi_new(l,k)+DSQRT(DBLE(j-1))*conjg(dGS(i))*dGS(j)
                 end if
              enddo
           enddo
           psi_mod(l,k)=DSQRT(REAL(psi_new(l,k)*CONJG(psi_new(l,k))))
           psi_modAvg= psi_modAvg+psi_mod(l,k)
           rhoAvg=rhoAvg+rho(l,k)
        end do
     end do
     psi_modAvg= psi_modAvg/DBLE(d**2)
     rhoAvg=rhoAvg/DBLE(d**2)
     !IF(DABS(rho(d/2,d/2)-1.0D0).GT.0.0001)THEN
     !   t_spacing=t_spacing*0.5D0
     !    PRINT*,t_spacing,rho(d/2,d/2)
     !   GOTO 9876
     !END IF
     
     psi_i=psi_new
     GS_mat=GS_matN

     !OPEN(UNIT=193,FILE='psi_real',status="unknown", position="append", action="write")
     !DO l=1,d
     !   WRITE(193,*)REAL(psi_new(l,:))
     !END DO
     !CLOSE(193)
     
  !OPEN(UNIT=194,FILE='psi_imag',status="unknown", position="append", action="write")
  !DO l=1,d
  !   WRITE(194,*)AIMAG(psi_new(l,:))
  !END DO
  !CLOSE(194)

  !OPEN(UNIT=195,FILE='rho',status="unknown", position="append", action="write")
  !DO l=1,d
  !   WRITE(195,*)AIMAG(psi_new(l,:))
  !END DO
  !CLOSE(195)

     ! Update the old psis as new psis
     !IF(rho(1,1)>1.5)STOP
     ! OPEN(UNIT=193,FILE='psi_real',status="unknown", position="append", action="write")
     !DO l=1,d
     ! WRITE(193,*)time,Jth,psi_mod
     !END DO
     !CLOSE(193)
     !open(unit=10, file="psi_mod.dat", form="unformatted", access="stream",status="unknown", position="append", action="write")
     !DO l=1,d
     !   WRITE(10)time,Jth,psi_mod
     !END DO
     !CLOSE(10)
     !OPEN(UNIT=193,FILE='psi1',status="unknown", position="append", action="write")
     !WRITE(193,'(5F15.8)')time,Jth,rho(1,1),SQRT(REAL(psi_new(1,1)*CONJG(psi_new(1,1))))
     !CLOSE(193)
     !OPEN(UNIT=194,FILE='psih',status="unknown", position="append", action="write")
     !WRITE(194,'(5F15.8)')time,Jth,rho(d/2,d/2),SQRT(REAL(psi_new(d/2,d/2)*CONJG(psi_new(d/2,d/2))))
     !CLOSE(194)
     OPEN(UNIT=195,FILE='avg',status="unknown", position="append", action="write")
     WRITE(195,'(5F15.8)')time,Jth,rhoAvg,psi_modAvg
     CLOSE(195)
     !IF(t_ind>50)STOP
     
     DO wi=1,noofwantedfiles
        !PRint*,"INSIDE",time,wi
        IF(DABS(time-tw(wi)).LE.0.5*t_spacing)THEN
           !PRint*,"INSIDE"
           write (filename,"(I2)")wi
           filename=trim(filename)
           OPEN(UNIT=291,FILE=filename,status="unknown", position="append", action="write")
           DO l=1,d
              WRITE(291,*)psi_mod(l,:)
           END DO
           CLOSE(291)
        END IF
        !STOP
     ENDDO
     
     time=time+t_spacing
     Jth=J_i+(Jc-J_i)*(time+TQ)/TQ
     !STOP
  end DO
 
100 FORMAT(2F9.2,34F9.4)
STOP
end program comp
