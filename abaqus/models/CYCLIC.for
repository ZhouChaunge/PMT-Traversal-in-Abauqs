       subroutine UEXTERNALDB(LOP, LRESTART, TIME, DTIME, KSTEP, KINC)
c     !! open files at the beginning of the analysis and close them at the end of the analysis

      include 'aba_param.inc'
      DIMENSION TIME(2)
		
      integer :: lenoutdir
      character(len=256) :: outdir, flname

      if (LOP == 0) then
	    call GETOUTDIR(outdir, lenoutdir)
        flname = trim(outdir) // '\info.txt'
        open(unit=100, file=flname, status='UNKNOWN')
        close(unit=100, status='DELETE')
        open(unit=100, file=flname, status='NEW')
       
        flname = trim(outdir) // '\stat.txt'
        open(unit=200, file=flname, status='UNKNOWN')
        close(unit=200, status='DELETE')
        open(unit=200, file=flname, status='NEW')
      else if (LOP ==3) then
	    close (unit=100, status='KEEP')
	    close (unit=200, status='KEEP')
      end if
      end subroutine UEXTERNALDB

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,DRPLDE
     1 ,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     2 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,
     3 DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)
      REAL*8   YYY, ZZZ,ZZETAI
      Real*8, dimension(NTENS) ::  DES,SIGT,DSIGT,DSIGT1,DSIGT2,ZEE,ZGG
      Real*8, dimension(NTENS)::ZHARBETA, SIG, DE, DEPSON, DSIG,HARBETA,
     1 DHARBETA,DSIG1,DSIG2,ZHARBETAI ,DHARBETA1,DHARBETA2
  
      Real(8), dimension(NTENS,NTENS) :: DEPT
      Real(8) :: MM, EPS0, POI, ZRAMDA, ZKAPA, ZM, ZMS, ZBR, 
     1 TNONE, ZR, ZRS, ZZETA, VP, VP0
      Real(8) :: SIGMM, E0, DLAMDA, TAVAL, VPT, PMSTAR, PMBAR, 
     1 DRS1,DVP1,DRS2,DVP2,DVP,DRS,EET,EET1,PCC,DV,EET2
      Integer :: I, J,NSS,NNN,NTENS,NDI,NPROPS !iSubStep, NSubStep, ICLAY   
C      !propoetry constants
      EPS0 = PROPS(1)         
C      !e0 initial void ratio
      POI = PROPS(2)          
C      !possoin ratio
      MM = PROPS(3)    !3.0D0*(RF-1.0D0)/(RF+2.0D0) MM=(M=q/p')critical      
C      !Rf= (sigma1'/sigma3')|failure   Mf=3(Rf-1)/(Rf+2)
      ZRAMDA = PROPS(4)       
C      !lambad
      ZKAPA = PROPS(5)        
C      !kappa
C      ! 3 special constant parameters:   over consolidation ratio; structure; unisotropic
      ZM = PROPS(6)
      ZMS = PROPS(7)
      ZBR = PROPS(8)
       !ICLAY = PROPS(9)        
c      !flag, 1-clay; 0-sand
      TNONE = PROPS(9)
c      !stae variables
c      !initial epsilonvp  
C                *** *** *** *** *** *** *** *** ** 给定参数
      PCC = STATEV(1)        !初始前期固结压力but 输入时代表OCR
c      !R
      ZRSI = STATEV(2)         
c      !R*   
      VPS = STATEV(3)         
c      !epsilonvp   plastic volume strain
      !VP0 = STATEV(4)
      EET = STATEV(4)
      NNN = STATEV(5)  !nnn=0时初始化状态参数
      !zeta
      ZZETAI = STATEV(6)
c      !beta ij
      ZHARBETAI(1) = STATEV(7)
      ZHARBETAI(2) = STATEV(8)
      ZHARBETAI(3) = STATEV(9)
      ZHARBETAI(4) = STATEV(10)
      IF(NTENS > 4)THEN
        ZHARBETAI(5) = STATEV(11)
        ZHARBETAI(6) = STATEV(12)
      END IF
c      !Pass on 
      SIG = 0.0D0
      DE = 0.0D0
       DO I=1,NTENS
          SIG(I)=-STRESS(I)
       END DO
C                *** *** *** *** *** *** *** *** ** 获取初始应力
       IF(NNN.EQ.0)THEN
          SIGMA0=0.0D0
          DO I=1,NDI
             SIGMA0=SIGMA0+SIG(I)
          END DO
          
          SIGMA0 = SIGMA0/NDI
          PM0 = SIGMA0*(MM*MM)/(MM*MM-ZZETAI*ZZETAI)
          PMBAR  =PM0*PCC
          PCC =ZRSI*PMBAR
          EET = EPS0-ZRAMDA*LOG(PCC/TNONE)+ZKAPA*LOG(PCC/PM0)
          STATEV(4) = EET
          ZHARBETAI(1)= 2.D0/3.D0*ZZETAI
          ZHARBETAI(2)= -1.D0/3.D0*ZZETAI
          ZHARBETAI(3)= -1.D0/3.D0*ZZETAI
          ZHARBETAI(4)= 0.D0
          
          
         IF(NTENS > 4)THEN
            ZHARBETAI(5)= 0.D0    
            ZHARBETAI(6)= 0.D0
         END IF
          NNN = 1
       END IF 
       ZDT=1.0D0
       ZT=0.0D0
C                *** *** *** *** *** *** *** *** ** T 以及 DT
       DO I=1,NTENS
          DES(I)=0.0D0
       END DO
C                *** *** *** *** *** *** *** *** ** 获取子步应变增量
       
          
117    CONTINUE
         ZDVP=0.0D0
         VP=VPS
         DRS=0.0
         ZRS=ZRSI
         ZZETA=ZZETAI
         ZHARBETA(:) = ZHARBETAI(:)
         EET1 = EET
       DO I=1,NTENS
         DES(I)=ZDT*(-DSTRAN(I))
      END DO 
C
      IF(Kstep > 3)THEN
      write(unit=200,fmt='(A,F12.6,A,I,A,I,A,I)') 'DSTRAN2:', DSTRAN(2),
     1'   Kstep:', Kstep,'   Kinc:', Kinc,'   NPT:', NPT
      END IF
      
          NSS=0
       DO 10 I=1,NTENS
10    SIGT(I)=SIG(I)
C
       
115    CONTINUE
       NSS=NSS+1 
       DHARBETA=0.0
C
       CALL GET_DMatrix(PROPS,DEPT,SIGT,DES,ZRS,ZZETA,ZHARBETA,
     1   DHARBETA,VP,DVP,DRS,NPROPS,NDI,NTENS,PMSTAR,PMBAR,EET1,
     2   PCC)
          DO 50 J=1,NTENS
50        DSIGT(J)=0.0D0
C
       DO 260 I=1,NTENS
       DO 261 J=1,NTENS
          DSIGT(I)=DSIGT(I)+DEPT(I,J)*DES(J)
261    CONTINUE
          SIGT(I)=SIGT(I)+DSIGT(I)
260   CONTINUE
C 
c      IF(Kstep > 3)THEN
c      write(unit=200,fmt='(A,F12.6,A,I,A,I,A,I)') 'DSIGT(2):', DSIGT(2),
c    1'   Kstep:', Kstep,'   Kinc:', Kinc,'   NPT:', NPT
c      write(unit=100,fmt='(A,F12.6,A,I,A,I,A,I)') 'DEPT(2,2):',
c    1 DEPT(2,2),'    Kstep:', Kstep,'    Kinc:', Kinc,'    NPT:', NPT
c      END IF

       DV=0.0D0
       DO I=1,NDI 
          DV=DV+DES(I)
       END DO
       EET1=EET1-(1.0D0+EET1)*DV
C
      IF(NSS.EQ.1)THEN
          DSIGT1(:)=DSIGT(:)
          DVP1=DVP
          DRS1=DRS
		  DHARBETA1 = 0.0D0
          DHARBETA1(:) = DHARBETA(:) 
          PMSTAR1 = PMSTAR
          PMBAR1 = PMBAR
          EET2 = EET1
       GOTO 115
      ELSE IF(NSS.EQ.2)THEN
          DSIGT2(:)=DSIGT(:)
          DVP2=DVP
          DRS2=DRS
		  DHARBETA2 = 0.0D0
          DHARBETA2(:) = DHARBETA(:)
          PMSTAR2 = PMSTAR
          PMBAR2 = PMBAR
      END IF
C
      DO 7 I=1,NTENS
7        ZEE(I)=(DSIGT2(I)+DSIGT1(I))/2.0D0
      DO 8 I=1,NTENS
8     ZGG(I)=(DSIGT2(I)-DSIGT1(I))/2.0D0
C

      ZDVP=(DVP1+DVP2)/2.0D0
      ZDRS=(DRS2+DRS1)/2.0D0
      DO I =1, NTENS
         DHARBETA(I) = (DHARBETA1(I) + DHARBETA2(I))/2.0D0
      END DO
C
      YYY=0.0D0
      ZZZ=0.0D0
      DO I=1,NTENS
         YYY=YYY+ZGG(I)**2
         ZZZ=ZZZ+(SIG(I)+ZEE(I))**2
      END DO
C     
      YYY=DSQRT(YYY)
      ZZZ=DSQRT(ZZZ)
C
      ZR1=YYY/ZZZ
      ZR2=DABS((PMBAR2-PMBAR1))/DABS(PMBAR)/2.0D0
      ZR3=DABS((PMSTAR2-PMSTAR2))/DABS(PMSTAR)/2.0D0
      ZRR=MAX(ZR1,ZR2,ZR3)
C      
      SSTOL=0.0001
C       
      ZB=0.8*DSQRT(SSTOL/ZRR)
      IF(ZRR.GT.SSTOL) THEN
         IF(ZB.LT.0.1)THEN
            ZB=0.1
         END IF
         ZDT=ZDT*ZB
         GOTO 117  
      ELSE
         IF(ZB.GT.2.0)THEN
            ZB=2.0
         END IF
         ZT=ZT+ZDT
         ZDT=ZDT*ZB
         GOTO 118   
      END IF 
C                *** *** *** *** *** *** *** *** ** 判断试算的子步应变增量是否符合
118   CONTINUE
      DO 11 I=1,NTENS
11       SIG(I)=SIG(I)+ZEE(I)
C
         VPS=VPS+ZDVP
         ZRSI=ZRSI+ZDRS
         IF(ZRSI > 1.0) ZRSI = 1.0
        DO I=1,NTENS
            ZHARBETAI(I)=ZHARBETAI(I)+DHARBETA(I)
        END DO
        
        ZZETAI=0.D0
      DO I=1,NDI
         ZZETAI=ZZETAI+ZHARBETAI(I)**2
      END DO
      DO I=NDI+1,NTENS
           ZZETAI=ZZETAI +2.D0*ZHARBETAI(I)**2
      END DO
        ZZETAI=DSQRT(3.D0/2.D0*ZZETAI)        
        EET=EET2
C
C      
      IF(ZT.LT.1)THEN
          IF(ZDT.GT.(1-ZT))THEN
          ZDT=1.0-ZT
          END IF
          GOTO 117
      END IF
C                *** *** *** *** *** *** *** *** **  继续下一个应变增量计算
      DO I=1,NTENS
         STRESS(I)=-SIG(I)
      END DO
C
      
      DO I=1,NTENS
         DE(I)=-DSTRAN(I)
      END DO

      STATEV(1)=  PCC          
c      !R
      STATEV(2) = ZRSI        
c      !R*   
      STATEV(3)= VPS     
c     !epsilonvp   plastic volume strain
      STATEV(4) = EET
      STATEV(5) = 1  !nnn=0时初始化状态参数
      STATEV(6)  = ZZETAI
c      
      STATEV(7)=ZHARBETAI(1)
      STATEV(8)=ZHARBETAI(2)
      STATEV(9)=ZHARBETAI(3) 
      STATEV(10)=ZHARBETAI(4) 
      IF(NTENS > 4)THEN
        STATEV(11)=ZHARBETAI(5) 
        STATEV(12)=ZHARBETAI(6) 
      END IF            
      CALL GET_DMatrix(PROPS,DEPT,SIG,DE,ZRSI,ZZETAI,ZHARBETAI,
     1   DHARBETA,VPS,DVP,DRS,NPROPS,NDI,NTENS,PMSTAR,PMBAR,
     2   EET,PCC)
C
      DO I=1,NTENS
      DO J=1,NTENS
         DDSDDE(I,J)=DEPT(I,J)
      END DO
      END DO
C     
      RETURN
      END
      
      
      
      SUBROUTINE GET_DMatrix(PROPS,DEPT,SIG,DDSTRAN,ZRS,ZZETA,
     1ZHARBETA,DHARBETA,VP,DVP,DRS,NPROPS,NDI,NTENS,
     2 PMSTAR,PMBAR,EET,PCC) 
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION  PROPS(NPROPS)
      Real*8, dimension(NTENS) :: ZHARBETA, SIG,DEPSON,
     1 HARBETA, ETA, ETASTAR, DHARBETA, ETABB, DDSTRAN,FF,WW,ZZ
      Real(8), dimension(NTENS,NTENS) :: DEPT
      Real(8) :: MM, EPS0, POI, RF, ZRAMDA, ZKAPA, ZM, ZMS, ZBR, 
     1TNONE, ZR, ZRS, ZZETA, VP, VP0, RT3, RT2, ZMB, HEE, SRT1, SRT2, 
     2SJ1, SJ2, SJ3,ZBO,EET,PCC
      Real(8) :: SIGMM, E0, DLAMDA, DMU, DMU2, ZKT0, SJJ, SJ, 
     1ROOT1, QC, SX, SY, SZ, DSX, DSY, DSZ, DSXY, DSYZ, DSZX, ETAVAL, 
     2ETASTARVAL, EBS, PM, VPT, PMSTAR, PMBAR, TT1, TT2, TT3, 
     3ETABBVALUE,SETAB,SBETAB
      Real(8) :: HDD,TMS, SS1, STRSTD, ZMUS, ZMU, HHH, FIJ, D, DF, 
     2DG, RAMDA, GAMA, DRS, DVP
      Integer :: I, J,  ICLAY
      DIMENSION DXX(NDI)
!      
      EPS0 = PROPS(1)         
c      !e0 initial void ratio
      POI = PROPS(2)          
c      !possoin ratio
      MM = PROPS(3)   !3.0D0*(RF-1.0D0)/(RF+2.0D0)        
c      !Rf= (sigma1'/sigma3')|failure   Mf=3(Rf-1)/(Rf+2)
      ZRAMDA = PROPS(4)   ! compression index      
c      !lambad
      ZKAPA = PROPS(5)        ! swelling index
c      !kappa
c      ! 3 special constant parameters:   over consolidation ratio; structure; unisotropic
      ZM = PROPS(6)! parameter controlling development of OCR
      ZMS = PROPS(7)! parameter controlling development of structure
      ZBR = PROPS(8)! parameter controlling development of anisotropy 
      TNONE = PROPS(9)! reference stress (usually taken as atmospheric stress 98kPa)
!
      RT3=1.7320508100D0! ROOT 3
      RT2=1.4142135624D0! ROOT 2
      ZMB=0.95
      HEE=(ZRAMDA-ZKAPA)/(1.0D0+EET) 
C
      SIGMM=0.D0
      DO I=1, NDI
          SIGMM = SIGMM + SIG(I)
      END DO
      SIGMM = SIGMM/NDI
C      
      E0=3.0D0*(1.0D0-2.0D0*POI)*(1.0D0+EET)*SIGMM/ZKAPA
      DLAMDA=E0*POI/(1.0D0+POI)/(1.0D0-2.0D0*POI)
      DMU=E0/2.0D0/(1.0D0+POI)
      DMU2=DMU*2.0D0
      ZKT0=E0/3.0D0/(1.0D0-2.0D0*POI)
C
       DO 126 I=1,NTENS
       DO 126 J=1,NTENS
126       DEPT(I,J)=0.0D0
C
      DO K1=1,NDI
         DO K2=1,NDI
            DEPT(K1,K2)=DLAMDA
            IF(K1.EQ.K2) DEPT(K1,K2)=DLAMDA+DMU2
         END DO
      END DO
c
      DO K1=NDI+1,NTENS
         DEPT(K1,K1)=DMU
      END DO
C
      DXX =0.0D0
      DO I=1,NDI
          DXX(I)=SIG(I)-SIGMM
      END DO
      ETA = 0.0D0
      DO I=1,NDI
         ETA(I)=DXX(I)/SIGMM
      END DO
      DO I=NDI+1,NTENS
         ETA(I)=SIG(I)/SIGMM
      ENDDO
      ETAVAL=0.D0
      DO I=1,NDI
          ETAVAL=ETAVAL+ETA(I)*ETA(I)
      END DO
      DO I= NDI+1, NTENS
          ETAVAL = ETAVAL + 2.D0 * ETA(I)*ETA(I)
      END DO
      ETAVAL=DSQRT(3.0D0/2.0*ETAVAL)
      ETASTAR=0.0D0
      DO I=1,NTENS
          HARBETA(I)=ZHARBETA(I)  
          ETASTAR(I)=ETA(I)-HARBETA(I)  
      END DO
      ETASTARVAL=0.D0
      DO I=1,NDI
          ETASTARVAL=ETASTARVAL+ETASTAR(I)**2
      END DO
      DO I=NDI+1,NTENS
        ETASTARVAL=ETASTARVAL +2.D0*ETASTAR(I)**2
      END DO
      ETASTARVAL=DSQRT(3.D0/2.D0*ETASTARVAL)
C
C     Calculate etabb, modified by yebin on 2008.4.24
      ZBO=ZMB*RT2/RT3*MM
      ETABB = 0.0D0
      DO I=1,NTENS
        IF(ETASTARVAL.LE.1.0E-6) THEN
           ETABB(I)=0.D0
        ELSE
           ETABB(I)=ZBO*RT3/RT2*ETASTAR(I)/ETASTARVAL-HARBETA(I)
        END IF
      END DO
C
      ETABBVALUE=0.D0
      DO I=1,NDI
        ETABBVALUE=ETABBVALUE+ETABB(I)*ETABB(I)
      END DO
      DO I =NDI+1,NTENS
          ETABBVALUE=ETABBVALUE +2.D0*ETABB(I)*ETABB(I)
      END DO     
	ETABBVALUE=DSQRT(ETABBVALUE)
C
	IF(ETABBVALUE.NE.0.0D0)THEN
        DO I=1,NTENS
	    ETABB(I)=ETABB(I)/ETABBVALUE
	  END DO
	END IF
C     calculate etabb end
      SETAB=0.D0
      DO I=1,NDI
        SETAB=SETAB+ETASTAR(I)*ETABB(I)
      END DO
      DO I=NDI+1,NTENS
        SETAB=SETAB +2.D0*ETASTAR(I)*ETABB(I)
      END DO
C
      SBETAB=0.D0
      DO I=1,NDI
        SBETAB=SBETAB+HARBETA(I)*ETABB(I)
      END DO        
      DO I=NDI+1,NTENS
         SBETAB=SBETAB +2.D0*HARBETA(I)*ETABB(I)
      END DO
      PM=SIGMM*(MM**2-ZZETA**2+ETASTARVAL**2)/(MM**2-ZZETA**2)
      PMSTAR=PCC*EXP(VP/HEE)
      PMBAR=PMSTAR/ZRS
      ZR=PM/PMBAR
      IF(ZR > 1.0)ZR = 1.0
      TT1=1.0D0/(MM*MM-ZZETA*ZZETA+ETASTARVAL*ETASTARVAL)/SIGMM
      TT2=MM*MM-ETAVAL*ETAVAL
      FF =0.0D0
      DO I=1,NDI
         FF(I)=TT1*(3.D0*ETASTAR(I)+TT2/3.D0)
      END DO
      DO I=NDI+1,NTENS
          FF(I)=6.D0*TT1*ETASTAR(I)
      END DO   
      HDD=TT1*TT2 
      SS1=MM*MM*SETAB-ETASTARVAL**2*SBETAB-ZZETA**2*SETAB
      STRSTD=(SIGMM/98.D0)**2+0.01
      TMS=TT2+6.D0*MM*ZBR*(1.D0-(ETAVAL/MM))*(ZMB*MM-ZZETA)
     1    *SS1*ETASTARVAL/(MM*MM-ZZETA**2)
     2  /(MM*MM-ZZETA**2+ETASTARVAL**2)
     3  -2.D0*ZMS*MM*(1.D0-ZRS)*ETASTARVAL
     4  -ZM*MM*DLOG(ZR)/ZR
     5  *DSQRT(6.D0*ETASTARVAL**2+TT2*TT2/3.D0)*(STRSTD/(STRSTD+1.D0))
c     
      HHH=TT1*TMS/HEE
      FSS=0D0
      DO I=1,NDI
         FSS=FSS+FF(I)*FF(I)
      END DO
       FDD=0D0
      DO I=NDI+1,NTENS
        FDD=FDD+FF(I)*FF(I)
      END DO 
C      
      FIJ= 2.0D0*FSS + FDD
      HHH=DLAMDA*HDD*HDD+DMU*FIJ+HHH
      D=1.0D0/HHH
      DF=DLAMDA*HDD
      DG=DLAMDA*HDD
      FPP=0D0
      DO I=1,NDI
         FPP=FPP+(DMU2*FF(I)+DLAMDA*HDD)*DDSTRAN(I)
      END DO
      FQQ=0D0
      DO I=NDI+1,NTENS
         FQQ=FQQ+DMU*FF(I)*DDSTRAN(I)
      END DO       
C     
      RAMDA=FPP+FQQ
C
      GAMA=RAMDA*D
      IF(GAMA.LT.0.0D0) GOTO 720
         
      DRS=2.D0*MM*ZMS*ZRS*(1.D0-ZRS)*ETASTARVAL*GAMA*TT1/HEE  
      ZRS = ZRS + DRS
      TT3=2.D0*MM*ZBR*TT1*GAMA/HEE*ETASTARVAL*(ZMB*MM-ZZETA)
      DO I=1,NTENS
          DHARBETA(I)=TT3*ETABB(I)
          ZHARBETA(I)=ZHARBETA(I)+DHARBETA(I)
      END DO
C      
      ZZETA=0.D0
      DO I=1,NDI
          ZZETA=ZZETA+ZHARBETA(I)**2 
      END DO
      DO I=NDI+1,NTENS
        ZZETA=ZZETA+2.D0*ZHARBETA(I)**2
      END DO
C
      ZZETA= DSQRT(3.D0/2.D0*ZZETA)
      DEPSON = 0.0D0
      DO I=1,NTENS
        DEPSON(I)=GAMA*FF(I)
      END DO
      DVP=0.0D0
      DO I=1,NDI
         DVP=DVP+DEPSON(I)
      END DO
      VP=VP+DVP
	  WW = 0.0D0
      DO I=1,NDI
         WW(I)=DG+DMU2*FF(I)
      END DO
      DO I=NDI+1,NTENS
         WW(I)=DMU*FF(I)
      END DO
C
      ZZ = 0.0D0
      DO I=1,NDI
         ZZ(I)=DF+DMU2*FF(I)
      END DO
      DO I=NDI+1,NTENS
         ZZ(I)=DMU*FF(I)
      END DO
      DO I=1,NTENS
         DO J=1,NTENS
             DEPT(I,J)=DEPT(I,J)-D*WW(I)*ZZ(J)
         END DO
      ENDDO
C
720   RETURN 
      END SUBROUTINE 
      
      

