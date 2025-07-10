!Program is working in the OpenMP system. This contains upto the step of calculation of eigen-values and eigen_vector
!This program have the additional features to calculate the eigenfuction from the desired independent variable
!060719

program real
implicit none
double precision,dimension(:),allocatable::b,initial_array,M,SD
double precision,dimension(:,:),allocatable::final_array,EF,U_N
integer,dimension(:),allocatable::cont,k,ek
integer::i,j,i1,l10,j1,I12,L66,l100,nr,num_param,INFO, LWORK,I2,LDA,NEW13
double precision::r1,dp,CHI,R,ERROR2,SUM1993,SUM2,HP,H0,H0_V,HT,OMEGAM,OMEGAX,HS,CHI2,R01,S,RS
INTEGER,parameter ::lwmax = 5000
double precision,DIMENSION(:,:),ALLOCATABLE::C,DC,MATRIX,MATRIX_OUT,NEW_MATRIX
double precision,DIMENSION(:),ALLOCATABLE::W,DW,WORK,COFF
integer,parameter::data_points=32,actual_data_points=32,total_rows=1339                !total_rows=1339
double precision::dat(3,actual_data_points),U(2,total_rows)
character(len=90)::filename,coffname,file_id


l100=10
nr=1339                         !nr=1339
num_param=7

lda = NUM_PARAM


ALLOCATE(W(NUM_PARAM))
ALLOCATE(WORK(LWMAX))
ALLOCATE(DW(NUM_PARAM))

ALLOCATE(MATRIX(NUM_PARAM,NUM_PARAM+1))
ALLOCATE(MATRIX_OUT(LDA,NUM_PARAM))
ALLOCATE(NEW_MATRIX(NUM_PARAM,NUM_PARAM+1))

allocate(ek(num_param))
allocate(k(num_param))
allocate(cont(num_param))
allocate(b(l100))
allocate(initial_array(num_param+1))   
allocate(final_array(NUM_PARAM+1,nr))

ALLOCATE(SD(NUM_PARAM))
ALLOCATE(C(NUM_PARAM,NUM_PARAM))
ALLOCATE(DC(LDA,NUM_PARAM))
ALLOCATE(M(NUM_PARAM))

open(0,file='./beta_data.dat',status='old')
open(1,file='./input.dat',status='old')
read(0,*)dat
read(1,*)u
 close(0)
 close(1)
 
open(69,file='./DATA_p07.txt',status='replace')
open(6969,file='./COV_MATRIX_p07.txt',status='replace')
OPEN(56,FILE='./EIGENVALUES_p07.txt',status='replace')
OPEN(65,FILE='./EIGENVECTORS_p07.txt',status='replace')
OPEN(96,FILE='./EIGENFUNCTIONS_p07.txt',STATUS='replace')

OPEN(9696,FILE='./NEIGENFNS_p07.txt',STATUS='replace')
OPEN(13,FILE='./DATA_FINAL_N_p07.txt',STATUS='replace')
OPEN(1993,FILE='./FINAL_FILE_N_p07.txt',STATUS='replace')

final_array = 0
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(NUM_PARAM,OMEGAM,OMEGAX,H0,H0_V,L100,U,dat,final_array,nr) NUM_THREADS(12)
!$OMP DO

DO L66 = 1,nr                                           !nr is the number of patches
r1 = U(1,L66)                                           !Initial Prior, Upper limit
dp = U(2,L66)                                           !step size	
I12 = 0
k=0
cont = 0
ek = 0
initial_array = 0
        do i=1,l100                                  !creation of net on a particular patch
        b(i) = r1 - i*dp
        end do

        do l10 = 1 ,l100**num_param
                do j=1,num_param
                if(j .ne. 1)then  

                        if(ek(j) .eq. l100**(j-1))then
                        cont(j) = cont(j) + 1 
                        ek(j) = 0
                        end if

        10               if((cont(j-1) - l100*cont(j)) .eq. (l100 - 1))then
                         ek(j) = ek(j) + 1
                         end if

                  k(j) = cont(j-1) - cont(j)*l100 + 1  
                  else
                           if(mod((l10 - 1),l100) .eq. 0 .and. l10 .ge. l100)then
                            cont(1) = cont(1) + 1
                           end if
                 k(j) = l10 - cont(1)*l100
                 end if
                end do

!* Chi calculation	
chi=0
                do i=1,data_points        !38              
                r=dat(1,i)
                error2=dat(3,i)
               sum2=0
                                do j=1,NUM_PARAM
                                sum2=sum2+b(K(J))*(funct(r)**(j-1))
                                end do

                HP=sum2
                hs = dat(2,i)
                chi=chi+((HS-HP)**2)/(error2**2)
                end do

!* Finding the Minimum chi and the corresponding point in parameter space

          I12=I12+1
 
          if(I12 == 1)then
          chi2=chi
          end if
 
                         if (chi .le. chi2)then 
                         chi2=chi
                         initial_array(1)=chi

                          DO J=1,NUM_PARAM
                          initial_array(J+1) = b(K(J))
                          END DO 

                         end if

          END DO

        final_array(1,l66) = chi2
        DO J = 1,NUM_PARAM
        final_array(J+1,l66) = initial_array(J+1)
        END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

9999 FORMAT( 11(:,5X,F12.7) )

DO L66=1,nr
WRITE(69,*)(final_array(J,l66),J=1,NUM_PARAM+1)
END DO

! print*,'complete the first DR'
!*===================================Creating Covariance Matrix===============================================
 do j=1,NUM_PARAM
 R01=0
 S=0
        do i=1,nr
        R01=final_array(j+1,i)+R01                        !finding mean
        end do
 m(j)=R01/nr                                        !value of mean


                do i=1,nr
                S=S+((final_array(j+1,i)-m(j))**2)/(nr-1)   !finding variance
                end do
 
 SD(j)=S                                                       !value of variance
 C(j,j)=SD(j)
 end do
  
 do I1=1,NUM_PARAM-1
        do j=1,NUM_PARAM-I1
        RS=0
                do i=1,nr
                RS=RS+final_array(j+1,i)*final_array(j+1+I1,i)          !finding covariance matrix
                end do
!          C(j,j+I1)=(RS-m(j)*m(j+I1))/(nr-1)
                 C(j,j+I1)=(RS/(nr-1)) - (m(j)*m(j+I1))
         C(j+I1,j)=C(j,j+I1)
         end do
 end do

 do i=1,NUM_PARAM
!write(6969,9999)(c(i,j),j=1,NUM_PARAM)
write(6969,*)(c(i,j),j=1,NUM_PARAM)
 end do
 !print*,'complete calculating the covariance matrix'
!*=============================Creating The Eigenvectors and Eigenvalues======================================
!open(1893,file='cov_matrix_rd_1-a.txt',status='old')
!read(1893,*)C
! close(1893)

!9901 continue

      LWORK = -1
      CALL DSYEV( 'Vectors', 'U', NUM_PARAM, C, LDA, W, WORK, LWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      CALL DSYEV( 'Vectors', 'U', NUM_PARAM, C, LDA, W, WORK, LWORK, INFO )

!*     Check for convergence.

      IF( INFO.GT.0 ) THEN
         ! WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         STOP
      END IF

!*     Print eigenvalues.

     ! CALL PRINT_MATRIX( 'Eigenvalues', 1, NUM_PARAM, W, 1, 0)

!*     Print eigenvectors.

DO I2=1,NUM_PARAM
         DO I1=1,NUM_PARAM
         MATRIX_OUT(I1,I2) = C(I2,I1)
         END DO
END DO

     ! CALL PRINT_MATRIX( 'Eigenvectors (stored columnwise)', NUM_PARAM, NUM_PARAM, MATRIX_OUT,LDA,1)

!*================================ARRANGING EIGENVECTORS AS REQUIRED==========================================
DO I2 = 1,NUM_PARAM
MATRIX(I2,1) = W(I2)
END DO

DO I1=1,NUM_PARAM
  DO I2=1,NUM_PARAM
  MATRIX(I2,I1+1) = C(I1,I2)
  END DO
END DO

CALL ASCENDING(NUM_PARAM,MATRIX,NEW_MATRIX)

DO I1=1,NUM_PARAM
DW(I1)=NEW_MATRIX(I1,1)
END DO

DO I1=1,NUM_PARAM
        DO I2=1,NUM_PARAM
        DC(I2,I1) = NEW_MATRIX(I2,I1+1)
        END DO
END DO

CALL PRINT_MATRIX( 'Eigenvalues', 1, NUM_PARAM, DW, 1,0)
CALL PRINT_MATRIX( 'Eigenvectors (stored column-wise)', NUM_PARAM, NUM_PARAM, DC,LDA,1)

DEALLOCATE(WORK)
DEALLOCATE(final_array)
! print*,'complete calculating the eigenvalues'
!*======================================EIGENFUNCTONS=========================================================

ALLOCATE(EF(NUM_PARAM, actual_data_points))
ALLOCATE(U_N(NUM_PARAM, actual_data_points))

DO i=1,data_points
r=dat(1,i)
error2=dat(3,i)

  DO i1=1,NUM_PARAM
  sum2=0
    DO j=1,NUM_PARAM
    sum2=sum2 + DC(i1,j)*((funct(r))**(j-1))
    END DO
  EF(i1,i)=sum2
  ENDDO

WRITE(96,*)(EF(j,i),j=1,num_param)

SUM2 = 0
DO i1=1,num_param
SUM2 = SUM2 + EF(i1,i)*EF(i1,i)/(error2*error2) !from_Clarkson_and_Zenkle
ENDDO

DO i1=1,num_param
u_n(i1,i) = EF(i1,i)/SQRT(SUM2)
ENDDO
WRITE(9696,*)(u_n(j,i),j=1,num_param)
ENDDO
! PRINT*,'complete calculating the eigenfunction'
!*=====================================FINAL RUN FOR THE CHI SQUARE==========================================
DEALLOCATE(C)
ALLOCATE(final_array(NUM_PARAM+1,nr))
!$OMP DO

DO L66 = 1,nr
r1 = U(1,L66)                                            !Initial Prior, Upper limit
dp = U(2,L66)                                            !step size	
I12 = 0
k=0
cont = 0
ek = 0
initial_array = 0

        DO i=1,l100
        b(i) = r1 - i*dp
        ENDDO

        DO l10 = 1 ,l100**num_param
                do j=1,num_param
                if(j .ne. 1)then 

                        if(ek(j) .eq. l100**(j-1))then
                        cont(j) = cont(j) + 1 
                        ek(j) = 0
                        end if

        100             if((cont(j-1) - l100*cont(j)) .eq. (l100 - 1))then
                        ek(j) = ek(j) + 1
                        end if

                   k(j) = cont(j-1) - cont(j)*l100 + 1
                    else
                           if(mod((l10 - 1),l100) .eq. 0 .and. l10 .ge. l100)then
                            cont(1) = cont(1) + 1
                              end if
                   k(j) = l10 - cont(1)*l100
                   end if
                   end do                !ending loop of j=1,num_param

       !* Chi calculation	
 chi=0
                do i=1,data_points              
                r=dat(1,i)
                error2=dat(3,i)
                 sum2=0
                                do j=1,NUM_PARAM
                                sum2=sum2+b(K(J))*EF(J,I)
                                end do
                HP=sum2
                hs = dat(2,i)
                chi=chi+((HS-HP)**2)/(error2**2)
                end do                                       !i= 1,data_points

                !* Finding the Minimum chi and the corresponding point in parameter space
                !print*,chi	
         I12=I12+1
 
          if(I12 == 1)then
          chi2=chi
          end if

                          if (chi .le. chi2)then 
                          chi2=chi
                           initial_array(1)=chi

                          DO J=1,NUM_PARAM
                          initial_array(J+1) = b(K(J))
                          END DO 
                          
                          end if
                          END DO                             !ending loop of l10

         final_array(1,l66) = chi2
         DO J = 1,NUM_PARAM
         final_array(J+1,l66) = initial_array(J+1)
         END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

DO L66=1,nr
WRITE(13,*)(final_array(J,l66),J=1,NUM_PARAM+1)
END DO

! print*,'complete calculating the second DR'
!*==============================MINIMISING THE FINAL DATA-SET=================================================

 sum2 = final_array(1,1)

DO L66 = 1,nr
        IF (sum2 .GE. final_array(1,L66)) THEN
        J=L66
        SUM2=final_array(1,L66)
        END IF 
END DO
WRITE(1993,*)(final_array(I,J),I=1,NUM_PARAM + 1)

DO I = 1,NUM_PARAM + 1
initial_array(I)=final_array(I,J)
END DO

! print*,'complete finding the minimum'
!*========================================FINDING THE COEFFICIENTS AND THE FUNCTOINAL FORM==================
NEW13 = 2
DEALLOCATE(final_array)
ALLOCATE(final_array(NEW13,data_points))
DO I1 = NUM_PARAM - 5, NUM_PARAM

write(file_id,'(i10)')I1
filename = './resultant_p07_fn1_' // trim(adjustl(file_id)) // '.txt'
coffname = './resultant_coff_p07_fn1_' // trim(adjustl(file_id)) // '.txt'
open(file = trim(filename), unit = 1932, status = 'replace')
open(file = trim(coffname), unit = 1947, status = 'replace')

ALLOCATE(COFF(I1))
        
        DO J = 1,I1
        COFF(J) = 0
                DO I = 1, NUM_PARAM
                COFF(J) = initial_array(I+1)*DC(I,J) + COFF(J)
                END DO
        END DO

         write(1947, *)(coff(i),i=1,i1)

        DO I = 1,data_points
        SUM2 = 0
        r = dat(1,I)
                DO J = 1,I1
                SUM2 = SUM2 + COFF(J)*((funct(r))**(j-1))
                END DO
        final_array(2,I) = SUM2
        final_array(1,I) = DAT(1,I)       
        END DO
! print*,'one step forward'
        DO I=1,data_points
        WRITE(1932,*)(final_array(J,I),J=1,NEW13)
        END DO

DEALLOCATE(COFF)
 close(1932)
 close(1947)
END DO
! print*,'complete the program_calculation of functional form'
!*==============================SUBROUTINES AND FUNCTIONS=====================================================
1001 continue
CONTAINS

         FUNCTION funct(t)
	      double precision, intent(in) :: t
	      double precision::funct
              funct = t / (1 + t )
        end function funct

!*============================================================================================================
      SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA,new )
      CHARACTER*(*) ::   DESC
      INTEGER       ::   M,N,LDA
      double precision          :: A(LDA,*)

      INTEGER     ::     I93, J93, new

     ! WRITE(*,*)
     ! WRITE(*,*) DESC
      
      IF (NEW == 0) THEN
       DO I93 = 1, M
       WRITE(56, *) ( A( I93, J93 ), J93 = 1, N )     
       END DO
9998  FORMAT(11(:,1X,E13.3))
      ELSE
        DO I93 = 1, M
        WRITE(65, *) ( A( J93, I93 ), J93 = 1, N )                 !MAKING THE MATRIX COLUMN-WISE, BUT ACTUALLY IT IS STORED IN ROW-WISE
        END DO
9999 FORMAT( 11(:,5X,F12.7) )
      END IF
       
      END subroutine print_matrix

!=============================================================================================================
      SUBROUTINE ASCENDING(N1,X,X_OUT)
                INTEGER,INTENT(IN)::N1
                INTEGER:: flag,i93,j93,count
    		double precision:: temp,X(N1,N1+1)
    		double precision,INTENT(OUT)::X_OUT(N1,N1+1)

                   count=N1
        
         do
         flag=0
        
            do i93=1,count-1
                if (ABS(x(i93,1))>ABS(x(i93+1,1)))then
            
                        do j93 = 1,NUM_PARAM+1
                        temp=x(i93,j93)
                        x(i93,j93)=x(i93+1,j93)
                        x(i93+1,j93)=temp 
                        end do
                    
                    flag=1    
                    end if
                end do
                if (flag==0)  then 

                X_OUT = X

           exit
          end if
        end do
    
    END SUBROUTINE ASCENDING

end program real
