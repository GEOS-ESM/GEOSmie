      implicit none
      include 'params.h' 
      real*8 LINTERPOL, ERREVAL, ONE_CALC

C****************************************************************
      INTEGER NANG
      real*8 angl(NANG_MAX)
      real*8 F11ANG(NANG_MAX), F22ANG(NANG_MAX), F33ANG(NANG_MAX)
      real*8 F44ANG(NANG_MAX), F12ANG(NANG_MAX), F34ANG(NANG_MAX)
      real*8 alt_angl(NANG_MAX)
      real*8 F11ALT(NANG_MAX), F22ALT(NANG_MAX), F33ALT(NANG_MAX)
      real*8 F44ALT(NANG_MAX), F12ALT(NANG_MAX), F34ALT(NANG_MAX)
      COMMON  /FANG/ angl, F11ANG, F22ANG, F33ANG, F44ANG,
     &  F12ANG, F34ANG, alt_angl, F11ALT, F22ALT, F33ALT,
     &  F44ALT, F12ALT, F34ALT, NANG
C****************************************************************
      real*8 F11OUT(NANG_MAX), F22OUT(NANG_MAX), F33OUT(NANG_MAX)
      real*8 F44OUT(NANG_MAX), F12OUT(NANG_MAX), F34OUT(NANG_MAX)
      COMMON /MATR_OUT/ F11OUT,F22OUT,F33OUT,F44OUT,F12OUT,F34OUT
C****************************************************************
      INTEGER NG
      real*8 X(NG_MAX), W(NG_MAX)
      real*8 F11(NG_MAX), F22(NG_MAX), F33(NG_MAX)
      real*8 F44(NG_MAX), F12(NG_MAX), F34(NG_MAX)
      COMMON /GAUSS_COMM/  F11,F22,F33,F44,F12,F34,X,W,NG
C****************************************************************
      INTEGER LMAX, L1MAX
      REAL*8 AL1(NG_MAX),AL2(NG_MAX),AL3(NG_MAX),AL4(NG_MAX)
      REAL*8 BET1(NG_MAX),BET2(NG_MAX)
      COMMON /EXPCOEFS/ AL1,AL2,AL3,AL4,BET1,BET2,LMAX,L1MAX
C****************************************************************

      real*8 fiterr, prev_fiterr

      integer NG_BRA, NG_KET, NG_GUESS, iter
      real*8 ERR_BRA, ERR_KET
      real*8 DDELT,ang
      real*8 CNORM
      integer narg, iargc, i, IERR, OUT_UNIT, L1
      character*1000 arg, ifn, outfname

      integer ng_step, ng_num, ing
      ng_step = 50
      ng_num = NSPHER/ng_step

      OUT_UNIT = 12

C Check the number of command line arguments. 
C If not one then print usage and stop.
      narg = iargc()
      if (narg .ne. 1) then
         write(*,*) 'Wrong number of arguments:', narg
         call getarg(0,arg)
         write(*,*) 'Usage: ', trim(arg), ' <filename>'
         stop
      endif 

C Read in the scattering matrix to be expanded using the
C first comand line argument.
      call getarg(1,ifn)
      call READMATRIX(ifn)

      if (NSPHER .gt. 0) then
         NG = NSPHER
         fiterr = one_calc()
         write(*,*) , NG, fiterr
      else
         NG = MIN_NSPHER
         fiterr = one_calc()
         write(*,*) , NG, fiterr
         if (fiterr .gt. DESIRED_ERR) then
            do 
               prev_fiterr = fiterr
               NG = NG + DELTA_NG
               fiterr = one_calc()
               write(*,*) , NG, fiterr
               if (fiterr .le. DESIRED_ERR) then
C               if (NG .gt. 5000) then
                  exit
               endif
            enddo
         endif
      endif
      print *, 'FINAL NG', NG, 'Error', fiterr
      outfname = trim(ifn)//'.expan_matr'
      OPEN (UNIT=OUT_UNIT, FILE=outfname)
      do i=1,NANG
         write (OUT_UNIT, '(F6.2,X,4E15.5,2F11.5)'), angl(i)*R2D,
     &    F11OUT(i),F22OUT(i),F33OUT(i),F44OUT(i),F12OUT(i),F34OUT(i)
      end do
      close(UNIT=OUT_UNIT)
      
      outfname = trim(ifn)//'.expan_coeff'
      OPEN (UNIT=OUT_UNIT, FILE=outfname)
      CNORM=1D0/AL1(1)
      write (OUT_UNIT,'(X,I5,6F17.10)'), L1MAX-1, CNORM 
      DO L1=1,L1MAX
         AL1(L1)=AL1(L1)*CNORM
         AL2(L1)=AL2(L1)*CNORM
         AL3(L1)=AL3(L1)*CNORM
         AL4(L1)=AL4(L1)*CNORM
         BET1(L1)=BET1(L1)*CNORM
         BET2(L1)=BET2(L1)*CNORM
         write (OUT_UNIT,'(X,I5,6F17.10)'),
     &    L1-1,AL1(L1),AL2(L1),AL3(L1),AL4(L1),BET1(L1),BET2(L1)
      end do
      close(UNIT=OUT_UNIT)
      
      stop
      end
 
C***************************************************************
C***************************************************************

C one_clac calls SPHER_EXPAN and calculates the expansion coefficients
C for the specified number of Gaussian nodes, then calls MATR 
C to obtain the expanded matrix and returns the discrepancy 
C between the original and expanded matrices.

      DOUBLE PRECISION FUNCTION one_calc() result(fiterr)
      IMPLICIT NONE
      include 'params.h' ! NANG_MAX, NG_MAX
      real*8 LINTERPOL, ERREVAL
 
C****************************************************************
      integer NANG
      real*8 angl(NANG_MAX)
      real*8 F11ANG(NANG_MAX), F22ANG(NANG_MAX), F33ANG(NANG_MAX)
      real*8 F44ANG(NANG_MAX), F12ANG(NANG_MAX), F34ANG(NANG_MAX)
      real*8 alt_angl(NANG_MAX)
      real*8 F11ALT(NANG_MAX), F22ALT(NANG_MAX), F33ALT(NANG_MAX)
      real*8 F44ALT(NANG_MAX), F12ALT(NANG_MAX), F34ALT(NANG_MAX)
      COMMON  /FANG/ angl, F11ANG, F22ANG, F33ANG, F44ANG,
     &  F12ANG, F34ANG, alt_angl, F11ALT, F22ALT, F33ALT,
     &  F44ALT, F12ALT, F34ALT, NANG
C****************************************************************
      real*8 F11OUT(NANG_MAX), F22OUT(NANG_MAX), F33OUT(NANG_MAX)
      real*8 F44OUT(NANG_MAX), F12OUT(NANG_MAX), F34OUT(NANG_MAX)
      COMMON /MATR_OUT/ F11OUT,F22OUT,F33OUT,F44OUT,F12OUT,F34OUT
C****************************************************************
      INTEGER NG
      real*8 X(NG_MAX), W(NG_MAX)
      real*8 F11(NG_MAX), F22(NG_MAX), F33(NG_MAX)
      real*8 F44(NG_MAX), F12(NG_MAX), F34(NG_MAX)
      COMMON /GAUSS_COMM/  F11, F22, F33, F44, F12, F34, X, W, NG
C****************************************************************
      INTEGER LMAX, L1MAX
      REAL*8 AL1(NG_MAX),AL2(NG_MAX),AL3(NG_MAX),AL4(NG_MAX)
      REAL*8 BET1(NG_MAX),BET2(NG_MAX)
      COMMON /EXPCOEFS/AL1,AL2,AL3,AL4,BET1,BET2, LMAX, L1MAX
C****************************************************************

      real*8 fiterr_alt,  ang
      integer i,ierr
      fiterr_alt = 0.0D0

      L1MAX = NG
      LMAX=L1MAX-1
      call gauss (NG,0,0,X,W)
      do i=1,NG
         ang=dacos(X(i))
         F11(i) =  LINTERPOL(NANG, angl, F11ANG, ang, 1, IERR)
         F22(i) =  LINTERPOL(NANG, angl, F22ANG, ang, 1, IERR)
         F33(i) =  LINTERPOL(NANG, angl, F33ANG, ang, 1, IERR)
         F44(i) =  LINTERPOL(NANG, angl, F44ANG, ang, 1, IERR)
         F12(i) =  LINTERPOL(NANG, angl, F12ANG, ang, 1, IERR)
         F34(i) =  LINTERPOL(NANG, angl, F34ANG, ang, 1, IERR)
      end do 
      call SPHER_EXPAN (AL1,AL2,AL3,AL4,BET1,BET2,L1MAX)
      if (USE_ALT_ANG .eq. 1) then
         call MATR (AL1,AL2,AL3,AL4,BET1,BET2,LMAX,alt_angl, NANG)
         fiterr_alt = ERREVAL(NANG,angl,F11ALT,F11OUT,
     &        ang_min,ang_max,ERRTYP)
      endif
      call MATR (AL1,AL2,AL3,AL4,BET1,BET2,LMAX,angl, NANG)
      fiterr = ERREVAL(NANG,angl,F11ANG,F11OUT,
     &     ang_min,ang_max,ERRTYP)
      fiterr = max(fiterr, fiterr_alt)
C      print *, fiterr, fiterr_alt
      end

C***************************************************************
C***************************************************************
C READMATRIX: Given the input file name, ifn reads in the 
C scattering matrix elements and scattering angles. 
C The data should be arranged in seven columns as follows:
C
C angle, F11, F22, F33, F44, F12, F34
C
C If USE_ALT_ANG .eq. 1 then the matrix elements are lineary 
C interpolated to the alternative angle grid.

      SUBROUTINE READMATRIX(ifn)
      IMPLICIT NONE
      include 'params.h' ! NANG_MAX, NG_MAX
      real*8 LINTERPOL
      integer IN_UNIT, io, imtx_elem, i, IERR, L1
      character*1000 ifn, infname
C****************************************************************
      integer NANG
      real*8 angl(NANG_MAX)
      real*8 F11ANG(NANG_MAX), F22ANG(NANG_MAX), F33ANG(NANG_MAX)
      real*8 F44ANG(NANG_MAX), F12ANG(NANG_MAX), F34ANG(NANG_MAX)
      real*8 alt_angl(NANG_MAX)
      real*8 F11ALT(NANG_MAX), F22ALT(NANG_MAX), F33ALT(NANG_MAX)
      real*8 F44ALT(NANG_MAX), F12ALT(NANG_MAX), F34ALT(NANG_MAX)
      COMMON  /FANG/ angl, F11ANG, F22ANG, F33ANG, F44ANG,
     &  F12ANG, F34ANG, alt_angl, F11ALT, F22ALT, F33ALT,
     &  F44ALT, F12ALT, F34ALT, NANG
C****************************************************************
      real*8 ang
      character(1000) :: line
      IN_UNIT = 11
      
      infname = trim(ifn)
      print *, trim(ifn)
      OPEN (UNIT=IN_UNIT, FILE=trim(ifn), STATUS='OLD')
      DO i=1,NANG_MAX
         READ(IN_UNIT,*,IOSTAT=io) angl(i),F11ANG(i),F22ANG(i),
     &        F33ANG(i),F44ANG(i),F12ANG(i),F34ANG(i)
         IF (io > 0) THEN
            WRITE(*,*) 'READMATRIX: Problem reading input file', io
            STOP
         ELSE IF (io < 0) THEN
            NANG = i-1
            EXIT
         END IF
      END DO
      close(UNIT=IN_UNIT)
      if (i .eq. NANG_MAX) then 
         WRITE(*,*) 'READMATRIX: Too many angles in input', infname
         STOP
      endif
C      do i=1,NANG
C         print '(F6.2,X,6F11.5)', angl(i),
C     &    F11ANG(i),F22ANG(i),F33ANG(i),F44ANG(i),F12ANG(i),F34ANG(i)
C      end do

      if (USE_ALT_ANG .eq. 1) then

         alt_angl(1) = angl(1)
         alt_angl(NANG) = angl(NANG)
         do i=2,NANG/2
            alt_angl(i)=0.5*(angl(i-1) + angl(i))
         end do
         do i=NANG/2+1, NANG-1
            alt_angl(i)=0.5*(angl(i+1) + angl(i))
         end do
         
         do i=1,NANG
            ang = alt_angl(i)
            F11ALT(i) =  LINTERPOL(NANG, angl, F11ANG, ang, 1, IERR)
            F22ALT(i) =  LINTERPOL(NANG, angl, F22ANG, ang, 1, IERR)
            F33ALT(i) =  LINTERPOL(NANG, angl, F33ANG, ang, 1, IERR)
            F44ALT(i) =  LINTERPOL(NANG, angl, F44ANG, ang, 1, IERR)
            F12ALT(i) =  LINTERPOL(NANG, angl, F12ANG, ang, 1, IERR)
            F34ALT(i) =  LINTERPOL(NANG, angl, F34ANG, ang, 1, IERR)
         end do
      endif

      do i=1,NANG
         angl(i)=angl(i)*D2R
         alt_angl(i)=alt_angl(i)*D2R
      end do
      end 
C***************************************************************
C  CALCULATION OF THE EXPANSION COEFFICIENTS  *******************               
                                                                                
      SUBROUTINE SPHER_EXPAN (AL1,AL2,AL3,AL4,BET1,BET2,L1MAX)

C     IMPLICIT REAL*8 (A-H,O-Z)                                                 
      IMPLICIT NONE
      include 'params.h' ! NANG_MAX, NG_MAX
      REAL*8 AL1(NG_MAX),AL2(NG_MAX),AL3(NG_MAX),AL4(NG_MAX)
      REAL*8 BET1(NG_MAX),BET2(NG_MAX)
      REAL*8 P1(NG_MAX),P2(NG_MAX),P3(NG_MAX),P4(NG_MAX),
     *       COEF1(NG_MAX),COEF2(NG_MAX),COEF3(NG_MAX),COEF4(NG_MAX),               
     *       COEF5(NG_MAX),COEF6(NG_MAX),COEF7(NG_MAX),COEF8(NG_MAX)                
      REAL*8 D6
      COMMON /COEF/ COEF1,COEF2,COEF3,COEF4,COEF5,COEF6,COEF7,COEF8,D6          
      COMMON /P/ P1,P2,P3,P4                                                    

      INTEGER NG
      real*8 X(NG_MAX), W(NG_MAX)
      real*8 F11(NG_MAX), F22(NG_MAX), F33(NG_MAX)
      real*8 F44(NG_MAX), F12(NG_MAX), F34(NG_MAX)
      COMMON /GAUSS_COMM/  F11, F22, F33, F44, F12, F34, X, W, NG

      INTEGER*4 L1, L1MAX, L, I
      real*8 A2,A3,CL
      real*8 FF11,FF22,FF33,FF44,FF12,FF34,FM,FP,P1L1,P4L1,WI

      DO 150 L1=3,L1MAX                                                         
          L=L1-1                                                                
          COEF1(L1)=1D0/DFLOAT(L+1)                                             
          COEF2(L1)=DFLOAT(2*L+1)                                               
          COEF3(L1)=1D0/DSQRT(DFLOAT((L+1)*(L+1)-4))                            
          COEF4(L1)=DSQRT(DFLOAT(L*L-4))                                        
          COEF5(L1)=1D0/(DFLOAT(L)*DFLOAT((L+1)*(L+1)-4))                       
          COEF6(L1)=DFLOAT(2*L+1)*DFLOAT(L*(L+1))                               
          COEF7(L1)=DFLOAT((2*L+1)*4)                                           
          COEF8(L1)=DFLOAT(L+1)*DFLOAT(L*L-4)                                   
  150 CONTINUE
      DO 160 L1=1,L1MAX
          AL1(L1)=0D0                                                           
          AL2(L1)=0D0                                                           
          AL3(L1)=0D0                                                           
          AL4(L1)=0D0                                                           
          BET1(L1)=0D0                                                          
          BET2(L1)=0D0                                                          
  160 CONTINUE                                                                  
      D6=0.25D0*DSQRT(6D0)                                                      
      DO 300 I=1,NG                                                             
          CALL GENER (X(I),L1MAX)                                               
          WI=W(I)                                                               
          FF11=F11(I)*WI                                                        
          FF22=F22(I)*WI                                                        
          FF33=F33(I)*WI                                                        
          FF44=F44(I)*WI                                                        
          FF12=F12(I)*WI                                                        
          FF34=F34(I)*WI                                                        
          FP=FF22+FF33                                                          
          FM=FF22-FF33                                                          
          DO 260 L1=1,L1MAX                                                     
              P1L1=P1(L1)                                                       
              P4L1=P4(L1)                                                       
              AL1(L1)=AL1(L1)+FF11*P1L1                                         
              AL4(L1)=AL4(L1)+FF44*P1L1                                         
              AL2(L1)=AL2(L1)+FP*P2(L1)                                         
              AL3(L1)=AL3(L1)+FM*P3(L1)                                         
              BET1(L1)=BET1(L1)+FF12*P4L1                                       
              BET2(L1)=BET2(L1)+FF34*P4L1                                       
  260     CONTINUE                                                              
  300 CONTINUE                                                                  
      DO 350 L1=1,L1MAX                                                         
          CL=DFLOAT(L1-1)+0.5D0                                                 
          L=L1                                                                  
          AL1(L1)=AL1(L1)*CL                                                    
          A2=AL2(L1)*CL*0.5D0                                                   
          A3=AL3(L1)*CL*0.5D0                                                   
          AL2(L1)=A2+A3                                                         
          AL3(L1)=A2-A3                                                         
          AL4(L1)=AL4(L1)*CL                                                    
          BET1(L1)=BET1(L1)*CL                                                  
          BET2(L1)=BET2(L1)*CL                                                  
C          IF (DABS(AL1(L1)).LE.DDELT) GO TO 400                     
  350 CONTINUE                                                                  
C  400 L1MAX=L                                                                   
C      PRINT *, '*********   EXPANSION COEFFICIENTS   *********'
C      PRINT *, '  S     ALPHA 1    ALPHA 2    ALPHA 3',
C     *      '    ALPHA 4     BETA 1     BETA 2'                  
C      DO L=1,L1MAX                                                          
C         PRINT '(I4,1X,6F11.5)', L-1,AL1(L),AL2(L),AL3(L),AL4(L),
C     *       BET1(L),BET2(L)               
C      END DO

      return
      end
C***********************************************************************        
                                                                                
C  CALCULATION OF THE REQUISITE GENERALIZED SPHERICAL FUNCTIONS                           
                                                                                
      SUBROUTINE GENER (U,L1MAX)                                                
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      include 'params.h' ! NANG_MAX, NG_MAX
      REAL*8 P1(NG_MAX),P2(NG_MAX),P3(NG_MAX),P4(NG_MAX),                           
     *       COEF1(NG_MAX),COEF2(NG_MAX),COEF3(NG_MAX),COEF4(NG_MAX),               
     *       COEF5(NG_MAX),COEF6(NG_MAX),COEF7(NG_MAX),COEF8(NG_MAX)                
      COMMON /COEF/ COEF1,COEF2,COEF3,COEF4,COEF5,COEF6,COEF7,COEF8,D6          
      COMMON /P/ P1,P2,P3,P4                                                    
      DUP=1D0+U                                                                 
      DUM=1D0-U                                                                 
      DU=U*U                                                                    
      P1(1)=1D0                                                                 
      P1(2)=U                                                                   
      P1(3)=0.5D0*(3D0*DU-1D0)                                                  
      P2(1)=0D0                                                                 
      P2(2)=0D0                                                                 
      P2(3)=0.25D0*DUP*DUP                                                      
      P3(1)=0D0                                                                 
      P3(2)=0D0                                                                 
      P3(3)=0.25D0*DUM*DUM                                                      
      P4(1)=0D0                                                                 
      P4(2)=0D0                                                                 
      P4(3)=D6*(DU-1D0)                                                         
      LMAX=L1MAX-1                                                              
      DO 100 L1=3,LMAX                                                          
         C1=COEF1(L1)                                                           
         C2=COEF2(L1)                                                           
         C3=COEF3(L1)                                                           
         C4=COEF4(L1)                                                           
         C5=COEF5(L1)                                                           
         C6=COEF6(L1)                                                           
         C7=COEF7(L1)                                                           
         C8=COEF8(L1)                                                           
         CU1=C2*U                                                               
         CU2=C6*U                                                               
         L2=L1+1                                                                
         L3=L1-1                                                                
         DL=DFLOAT(L3)                                                          
         P1(L2)=C1*(CU1*P1(L1)-DL*P1(L3))                                       
         P2(L2)=C5*((CU2-C7)*P2(L1)-C8*P2(L3))                                  
         P3(L2)=C5*((CU2+C7)*P3(L1)-C8*P3(L3))                                  
         P4(L2)=C3*(CU1*P4(L1)-C4*P4(L3))                                       
  100 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
C****************************************************************
 
C    CALCULATION OF THE SCATTERING MATRIX FOR GIVEN EXPANSION
C    COEFFICIENTS
 
C    A1,...,B2 - EXPANSION COEFFICIENTS
C    LMAX - NUMBER OF COEFFICIENTS MINUS 1
C    N - NUMBER OF SCATTERING ANGLES
C        THE CORRESPONDING SCATTERING ANGLES ARE GIVEN BY
C        180*(I-1)/(N-1) (DEGREES), WHERE I NUMBERS THE ANGLES
 
      SUBROUTINE MATR(A1,A2,A3,A4,B1,B2,LMAX,angl,NANG)

      IMPLICIT REAL*8 (A-H,O-Z)
      include 'params.h'
      REAL*8 A1(NG_MAX),A2(NG_MAX),A3(NG_MAX),A4(NG_MAX)
      REAL*8 B1(NG_MAX),B2(NG_MAX)
      real*8 angl(NANG_MAX)

      real*8 F11OUT(NANG_MAX), F22OUT(NANG_MAX), F33OUT(NANG_MAX)
      real*8 F44OUT(NANG_MAX), F12OUT(NANG_MAX), F34OUT(NANG_MAX)
      COMMON /MATR_OUT/ F11OUT,F22OUT,F33OUT,F44OUT,F12OUT,F34OUT

      L1MAX=LMAX+1

C      PRINT 1000
C 1000 FORMAT(' ')
C      PRINT 1001
C 1001 FORMAT(' ',2X,'S',6X,'ALPHA1',6X,'ALPHA2',6X,'ALPHA3',
C     &       6X,'ALPHA4',7X,'BETA1',7X,'BETA2')
C      DO 10 L1=1,L1MAX
C         L=L1-1
C         PRINT 1002,L,A1(L1),A2(L1),A3(L1),A4(L1),B1(L1),B2(L1)
C   10 CONTINUE
C 1002 FORMAT(' ',I5,6F12.5)
C      PRINT 1000
C      PRINT 1003
C 1003 FORMAT(' ',5X,'<',8X,'F11',8X,'F22',8X,'F33',
C     & 8X,'F44',8X,'F12',8X,'F34')

      D6=DSQRT(6D0)*0.25D0
C      PI=dacos(-1d0)
C      D2R=PI/180d0
C      R2D=180d0/PI
      DO 500 I1=1,NANG
         TAA=angl(I1)
         TB=R2D*angl(I1)
         U=DCOS(TAA)
         F11=0D0
         F2=0D0
         F3=0D0
         F44=0D0
         F12=0D0
         F34=0D0
         P1=0D0
         P2=0D0
         P3=0D0
         P4=0D0
         PP1=1D0
         PP2=0.25D0*(1D0+U)*(1D0+U)
         PP3=0.25D0*(1D0-U)*(1D0-U)
         PP4=D6*(U*U-1D0)
         DO 400 L1=1,L1MAX
            L=L1-1
            DL=DFLOAT(L)
            DL1=DFLOAT(L1)
            F11=F11+A1(L1)*PP1
            F44=F44+A4(L1)*PP1
            IF(L.EQ.LMAX) GO TO 350
            PL1=DFLOAT(2*L+1)
            P=(PL1*U*PP1-DL*P1)/DL1
            P1=PP1
            PP1=P
  350       IF(L.LT.2) GO TO 400
            F2=F2+(A2(L1)+A3(L1))*PP2
            F3=F3+(A2(L1)-A3(L1))*PP3
            F12=F12+B1(L1)*PP4
            F34=F34+B2(L1)*PP4
            IF(L.EQ.LMAX) GO TO 400
            PL2=DL*DL1*U
            PL3=DL1*(DL*DL-4D0)
            PL4=1D0/(DL*(DL1*DL1-4D0))
            P=(PL1*(PL2-4D0)*PP2-PL3*P2)*PL4
            P2=PP2
            PP2=P
            P=(PL1*(PL2+4D0)*PP3-PL3*P3)*PL4
            P3=PP3
            PP3=P
            P=(PL1*U*PP4-DSQRT(DL*DL-4D0)*P4)/DSQRT(DL1*DL1-4D0)
            P4=PP4
            PP4=P
  400    CONTINUE
         F22=(F2+F3)*0.5D0
         F33=(F2-F3)*0.5D0
C        F22=F22/F11
C        F33=F33/F11
C        F44=F44/F11
C        F12=-F12/F11
C        F34=F34/F11
         F11OUT(I1) = F11
         F22OUT(I1) = F22
         F33OUT(I1) = F33
         F44OUT(I1) = F44
         F12OUT(I1) = F12
         F34OUT(I1) = F34
C         PRINT '(' ',F6.2,6F11.4)',TB,F11,F22,F33,F44,F12,F34
  500 CONTINUE
C      PRINT *, ' ' 
      RETURN
      END                                                                                

C********************************************************
      SUBROUTINE GAUSS ( N,IND1,IND2,Z,W )
      IMPLICIT REAL*8 (A-H,P-Z)
      REAL*8 Z(N),W(N)
      A=1D0
      B=2D0
      C=3D0
      IND=MOD(N,2)
      K=N/2+IND
      F=DFLOAT(N)
      DO 100 I=1,K
          M=N+1-I
          IF(I.EQ.1) X=A-B/((F+A)*F)
          IF(I.EQ.2) X=(Z(N)-A)*4D0+Z(N)
          IF(I.EQ.3) X=(Z(N-1)-Z(N))*1.6D0+Z(N-1)
          IF(I.GT.3) X=(Z(M+1)-Z(M+2))*C+Z(M+3)
          IF(I.EQ.K.AND.IND.EQ.1) X=0D0
          NITER=0
          CHECK=1D-16
   10     PB=1D0
          NITER=NITER+1
          IF (NITER.LE.100) GO TO 15
          CHECK=CHECK*10D0
   15     PC=X
          DJ=A
          DO 20 J=2,N
              DJ=DJ+A
              PA=PB
              PB=PC
   20         PC=X*PB+(X*PB-PA)*(DJ-A)/DJ
          PA=A/((PB-X*PC)*F)
          PB=PA*PC*(A-X*X)
          X=X-PB
          IF(DABS(PB).GT.check*DABS(X)) GO TO 10
          Z(M)=X
          W(M)=PA*PA*(A-X*X)
          IF(IND1.EQ.0) W(M)=B*W(M)
          IF(I.EQ.K.AND.IND.EQ.1) GO TO 100
          Z(I)=-Z(M)
          W(I)=W(M)
  100 CONTINUE
      IF(IND2.NE.1) GO TO 110
      PRINT 1100,N
 1100 FORMAT(' ***  POINTS AND WEIGHTS OF GAUSSIAN QUADRATURE FORMULA',
     * ' OF ',I4,'-TH ORDER')
      DO 105 I=1,K
          ZZ=-Z(I)
  105     PRINT 1200,I,ZZ,I,W(I)
 1200 FORMAT(' ',4X,'X(',I4,') = ',F17.14,5X,'W(',I4,') = ',F17.14)
      GO TO 115
  110 CONTINUE
C     PRINT 1300,N
 1300 FORMAT(' GAUSSIAN QUADRATURE FORMULA OF ',I4,'-TH ORDER IS USED')
  115 CONTINUE
      IF(IND1.EQ.0) GO TO 140
      DO 120 I=1,N
  120     Z(I)=(A+Z(I))/B
  140 CONTINUE
      !write(*,*) z
      RETURN
      END
 
                 
C**************************************************************
C LINTERPOL returns a linearly interpolated value of array YY of size NN,
C defined on the grid XX for a new value X. IEX determines whether
C extrapolation outside of the grid XX is allowed. If IEX is set to 1 and 
C the value of X is outside the range of XX the returned value will be 
C extrapolated from the two closest YY values. If IEX is not 1
C and X is below low end of XX then the first YY value is returned and
C IERR error code is set to -1. If IEX is not 1
C and X is above the high end of XX then the last YY value is returned 
C and the IERR error code is set to 1.

      DOUBLE PRECISION FUNCTION LINTERPOL(NN,XX,YY,X,IEX,IERR) result(Y)

      IMPLICIT NONE
      INTEGER NN , IERR, I, IEX
      REAL*8 XX(NN), YY(NN), X     

      IF (X .LT. XX(1)) THEN                                                       ! if below low end of XX (error)..
         if (IEX .eq. 1) THEN
            Y = (YY(1)-YY(2))/(XX(1)-XX(2))*(X-XX(1))+YY(1)                        !extrapolate
         ELSE                   
            Y = YY(1)                                                              !  set Y = first YY value
         ENDIF
         IERR = -1                                                                  !  return error code -1
      ELSE IF (X .GT. XX(NN)) THEN                                                  ! if above high end of XX (error)..
         if (IEX .eq. 1) THEN
            Y = (YY(NN)-YY(NN-1))/(XX(NN)-XX(NN-1))*(X-XX(NN))+YY(NN)               ! extrapolate
         ELSE
            Y = YY(NN)                                                              !  set Y = last YY value
         ENDIF
         IERR = +1                                                                  !  return error code +1
      ELSE                                                                          ! if OK
         DO I = 2, NN                                                               ! loop to find first XX > X
            IF (XX(I) .GT. X) EXIT
         END DO
         Y = (YY(I)-YY(I-1))/(XX(I)-XX(I-1))*(X-XX(I-1))+YY(I-1)                    ! interpolate
         IERR = 0                                                                   ! set error code to 0 (no error)
      END IF

      RETURN 

      END FUNCTION LINTERPOL

C**************************************************************
C ERRVAL returns the discrepancy between two arrays YY1 and YY2
C of the dimension NN specified on grid XX. Only the array values 
C corresponding to XX grid values between XMIN and XMAX
C are considered. The type of discrepancy is determined 
C by the ITYP argument. The supported values are 
C MAXABS, MAXPERCENT and MSRE as defined in the file params.h.
C They correspond to maximum absolute difference, maximum absolute  
C percent difference, and mean square root error, respectively.

      DOUBLE PRECISION FUNCTION ERREVAL
     &                    (NN,XX,YY1,YY2,XMIN,XMAX,ITYP) result(ERR)

      IMPLICIT NONE
      include 'params.h'
      INTEGER NN, I, ITYP, xcnt
      REAL*8 XMIN, XMAX, tmp,sq,mean
      REAL*8 XX(NN), YY1(NN), YY2(NN)

      if (xmin .gt. xmax) then
         tmp = xmax
         xmax = xmin
         xmin = tmp
      end if
      if (XX(1) .gt. xmax .or. XX(NN) .lt. xmin) then
         print *, 'ERRVAL: no points in range ', xmin, xmax
      end if
      err = 0.0D0
      xcnt  = 0
      mean = 0.0
      sq = 0.0
      do i=1,NN 
         if (xx(i) .lt. XMIN .or. xx(i) .gt. XMAX) then 
            cycle
         end if
         xcnt = xcnt + 1
         select case (ITYP)
            case (MAXABS)
               tmp = abs(yy1(i) - yy2(i))
            case (MAXPERCENT)
               tmp = 200D0*abs(yy1(i) - yy2(i))/abs(yy1(i) + yy2(i))
            case (MSRE)
               mean = mean + (yy1(i) - yy2(i))
               sq = sq + (yy1(i) - yy2(i))**2 
            case default
               WRITE(*,*) 'ERRVAL: Wrong error type specified', ITYP
               stop
         end select
         if (tmp .gt. err) then 
            err = tmp
         end if
      end do
      if (ITYP .eq. MSRE) then
        if (xcnt .lt. 2) then
          WRITE(*,*) 'ERRVAL: less than 2 points available for MSRE'
          stop
        end if 
        err  = sqrt((sq - mean*mean/xcnt)/(xcnt-1.0D0))
      end if

      RETURN 

      END FUNCTION ERREVAL
