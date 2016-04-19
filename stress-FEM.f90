!本程序计算弹性力学平面应力问题（带厚度），在F90上调试通过
!包含主程序STRESS和四个子程序ELSTMX(KK)、MODIFY、DCMPBD和SLVBD
!主程序主要进行数据输入、结果输出、最大半带宽的计算、组装总体刚度矩阵以及计算单元内
!的应力、应变值
!子程序ELSTMX(KK)计算刚度矩阵系数，生成第KK单元的刚度矩阵K
!子程序MODIFY输入载荷节点的载荷值，位移边界节点处的位移值，对总体刚度矩阵、位移数组
!和节点力数组进行修改
!子程序DEMPBD用高斯消元法将对称等带宽总刚度矩阵的上半部分化为上三角阵
!子程序SLVBD对子程序DEMPBD得到的上三角矩阵进行回代计算，得到节点位移向量DD
!
	PROGRAM STRESS
	COMMON /ELMATX/ESM(6,6),X(3),Y(3),D(3,3)
	COMMON / GRAD/ B(3,6),AR2
	COMMON /MTL/ EM,PR,TH
	COMMON /AV/A(8500),JGF,JGSM,NP,NBW,JEND
	DIMENSION NS(6),U(6),STRA(3),STRE(3)
	DATA IN/60/,IO/61/
!
!将单元以及组成单元的节点数组按照动态数组分配方法
	REAL,ALLOCATABLE::XC(:),YC(:)
	INTEGER,ALLOCATABLE::NEL(:,:)
!
!ESM(6,6)——单元刚度矩阵，X(3),Y(3)——单元节点坐标，D(3,3)材料性质矩阵
!B(3,6)——几何矩阵，AR2——三角形面积的二倍，NP——自由度总数，NBW——最大半带宽
!A(8500)——存储节点位移向量，节点力向量和总体刚度矩阵的数组A
!JGF，JEND为计数单元
!JGF=NP-节点位移向量在数组A中的位置，JGSM=JGF+NP-节点力向量在数组A中的位置
!JEND=JGSM+NP*NBW——刚度矩阵在数组A中的位置，数组A总长度
!NS(6)——一个单元节点的自由度编号数组，U(6)——一个单元的节点自由度
!IN/60/,IO/61/一个输入输出文件设备号，STRA(3),STRE(3)——存储单元的应变、应力
!
!---------------------------程序输入段-------------------------------------
! TITLE——存储计算内容标题
! NN——节点总数，NE——单元总数，EM——杨氏模量，PR——泊松比，TH——板的厚度，N——单元编号
! NEL(I)——节点的X轴坐标，YC(I)——节点的Y轴坐标
! NEL(N,I)——组成第N个三角形的第I个节点的编号(I=1,2,3)
!
	OPEN(60,FILE='INPUT.DAT',STATUS='UNKNOWN')
	OPEN(61,FILE='OUTPUT.DAT',STATUS='UNKNOWN')
	READ(IN,1) TITLE
1	FORMAT(900A)
	WRITE(IO,1) TITLE
	WRITE(IO,*)
! 
! 输入计算模型的节点总数NN和单元总数NE
	READ(IN,*)NN
	READ(IN,*)NE
	WRITE(IO,2)NN,NE
2	FORMAT(/,'输入数据:',/,'节点数=',I3,5X,'单元数=',I4,/)
	NP=2*NN
	ALLOCATE(NEL(1:NE,1:3),XC(NN),YC(NN))
! 
! 输入材料的杨氏模量EM，泊松比PR，平板厚度TH，节点坐标XC(I),YC(I)和组成单元的
! 节点NEL(N,I)
! 组成单元的节点的编号都按逆时针顺序输入
! 
	READ(IN,*)EM
	READ(IN,*)PR
	READ(IN,*)TH
	READ(IN,*)(XC(I),I=1,NN)
	READ(IN,*)(YC(I),I=1,NN)
! 
! 输出材料的性质和计算模型拓扑数据便于检查时对照
! 
	WRITE(IO,3)
3	FORMAT('材料常数为:')
	WRITE(IO,4)EM,PR,TH
4	FORMAT('弹性模量为:',E12.5,5X,'泊松比为:',E12.5,5X,'厚度为:',E12.5)
	WRITE(IO,*)
	DO 6 I=1,NN
	WRITE(IO,5)I,XC(I),YC(I)
5	FORMAT('节点',I3,'坐标为:X=',F8.3,3X,'Y=',F8.3)
6	CONTINUE
	WRITE(IO,*)
	DO 7 KK=1,NE
	READ(IN,*)N,(NEL(N,I),I=1,3)
	WRITE(IO,8)N,(NEL(N,I),I=1,3)
8	FORMAT('单元号码为:',I3,'组成单元的节点号码为:',I3,I3,I3)	
7	CONTINUE
	WRITE(IO,*)
! 
! ---------------------------计算开始-------------------------------
! 计算最大半带宽，B=MAX(D+1)*F
! D是一个单元各节点遍点号之差的最大值，F是一个节点的自由度数
! 
	INBW=0
	NBW=0
	DO 20 KK=1,NE
	DO 25 I=1,3
25	NS(I)=NEL(KK,I)
	DO 21 I=1,2
		IJ=I+1
	DO 21 J=IJ,3
		NB=IABS(NS(I)-NS(J))    !节点号之差的绝对值
		IF(NB.LE.NBW) GOTO 21
		INBW=KK
		NBW=NB
21 CONTINUE
20 CONTINUE
	NBW=(NBW+1)*2  !平面问题节点自由度F=2,NBW此时为最大半带宽
!
!数组A中数据的安排:A(1,2,3....NP|NP+1.....2NP|2NP+1.....JEND)
! 					节点位移向量  节点力向量   总体刚度矩阵
!初始化数组A
	JGF=NP     !JGF=2*NN
	JGSM=JGF+NP  !JGSM=4*NN
	JEND=JGSM+NP*NBW     !JEND=4*NN+2*NN*12
! NP为K中KNN的N，NBW为最大半带宽
!  	|XXXX|
! 	 |XXXX|
! 	  |XXXX|
! 	   |XXX*|
! 	    |XX**|
! 	     |X***|
!用等带宽二维数组方法存储的刚度矩阵中包含*表示的位置，但在实际程序中未使用
! 
	JL=JEND-JGF
	DO 24 I=1,JEND
24	A(I)=0.0
	GOTO 30
! 
! 生成材料性质矩阵
! E为杨氏模量
! V为泊松比
! 
30	R=EM/(1.-PR**2)
	D(1,1)=R
	D(2,2)=D(1,1)
	D(3,3)=R*(1.-PR)/2.
	D(1,2)=PR*R
	D(2,1)=D(1,2)
	D(1,3)=0.0
	D(3,1)=0.0
	D(2,3)=0.0
	D(3,2)=0.0
! 
! 单元矩阵循环的开始
! 
	KK=1
! KK为单元号
! 节点自由度的生成，节点坐标的局部变化
! 
32	DO 31 I=1,3
	J=NEL(KK,I)   !元件KK的三个节点编号
	NS(2*I-1)=J*2-1
	NS(2*I)=J*2   !元件KK各个自由度在总刚度矩阵中的方程号
	X(I)=XC(J)
31	Y(I)=YC(J)  !元件KK的节点坐标值
! 
! 调用子程序ELSTMX计算单元刚度矩阵ESM(6,6)并输出
! 
	CALL ELSTMX(KK)
! 
! 单元刚度矩阵组装成总体刚度矩阵
	DO 33 I=1,6
		II=NS(I)
		DO 34 J=1,6
			JJ=NS(J)+1-II
			IF(JJ.LE.0)GOTO 34
			J1=JGSM+(JJ-1)*NP+II-(JJ-1)*(JJ-2)/2
			A(J1)=A(J1)+ESM(I,J)
34		CONTINUE
33	CONTINUE
	KK=KK+1
	IF(KK.LE.NE)GOTO 32
! 
! 调用子程序MODIFY输入载荷节点处的载荷值，位移边界节点处的位移值，对总体刚度矩阵
! 位移数组和节点力数组进行相应的修改
	CALL MODIFY
		WRITE(IO,35)
35 FORMAT(/,'计算结果为:',/)
		WRITE(IO,36)NBW
36 FORMAT('最大半带宽为',I3)
		WRITE(IO,37)JEND
37 FORMAT('总数组大小为:',I3)
!
!调用子程序DEMPBD，用高斯消元法将对称等带宽总刚度矩阵化为上三角矩阵
	CALL DCMPBD
!
! 调用子程序SLVBD对子程序DEMPBD得到的上三角矩阵进行回代计算，得到节点位移向量DD
	CALL SLVBD
! 
! ---------------------------输出各个节点的位移向量----------------------------
! 
	WRITE(IO,*)
	DO 45 I=1,NP/2
		WRITE(IO,43)I,A(2*I-1),A(2*I)
43 	FORMAT('节点号',I3,5X,'X方向的位移UX=',E12.5,5X,'Y方向的位移UY=',E12.5)
45 	CONTINUE
! 
! ---------------------------附加计算-----------------------------------------
! 计算节点处的应变
! 
	DO 96 KK=1,NE
! 生成节点的自由度
! 节点坐标局部变化
	DO 51 I=1,3
		J=NEL(KK,I)
		NS(2*I-1)=2*J-1
		NS(2*I)=2*J
		X(I)=XC(J)
51		Y(I)=YC(J)
! 
! 单元节点位移局部变化
! 
65 DO 73 I=1,6,2
	NS1=NS(I)
	NS2=NS(I+1)
	U(I)=A(NS1)   !节点X方向的位移量
73	U(I+1)=A(NS2)	!节点Y方向的位移量
!
!计算单元应变
! 
	CALL ELSTMX(KK)
	DO 52 I=1,3
		STRA(I)=0.0
		DO 52 K=1,6
52			STRA(I)=STRA(I)+B(I,K)*U(K)/AR2
!
! 计算单元应力值
! 
	DO 58 I=1,3
		STRE(I)=0
		DO 58 K=1,3
58	STRE(I)=STRE(I)+D(I,K)*(STRA(K))
! 
! 计算主应力S1，S2，TM
! 
	AA=(STRE(1)+STRE(2))/2.
	AB=SQRT((ABS(STRE(1)-STRE(2))/2.)**2+STRE(3)**2)
	S1=AA+AB
	S2=AA-AB
	TM=AB
! 
! 计算主应力方向与X轴的夹角
! 
	IF(ABS(STRE(1)-STRE(2)).LT.0.0001) GOTO 93
	AC=ATAN2(2*STRE(3),(STRE(1)-STRE(2)))
	THM=(180/3.1415926*AC)/2
	GOTO 94
93	THM=90
! 
! ------------------------输出单元应变和应力的计算结果--------------------------
! 
94	WRITE(IO,57)KK
57	FORMAT(/,'单元',I4)
	WRITE(IO,95)STRA(1),STRA(2),STRA(3)
95	FORMAT('EPTOX=',E12.5,2X,'EPTOY=',E12.5,2X,'EPTOXY=',E12.5)
	WRITE(IO,97)STRE(1),STRE(2),STRE(3)
97	FORMAT('SX=',E12.5,5X,'SY=',E12.5,5X,'SXY=',E12.5)
	WRITE(IO,98)S1,S2,TM
98	FORMAT('S1=',E12.5,5X,'S2=',E12.5,5X,'TMAX=',E12.5)
	WRITE(IO,99)THM
99	FORMAT(44X,'ANGEL=',F8.2,'°')
96	CONTINUE
	CLOSE(60)
	CLOSE(61)
	STOP
	END

! 
! 计算刚度系数矩阵，计算第KK单元的刚度矩阵K
! 
	SUBROUTINE ELSTMX(KK)
		COMMON/MTL/EM,PR,TH
		COMMON/GRAD/B(3,6),AR2
		COMMON/ELMATX/ESM(6,6),X(3),Y(3),D(3,3)
		DIMENSION C(6,3)
		DATA IO/61/ IN/60/
! 
! 计算矩阵B
!
	DO 20 I=1,3
	DO 20 J=1,6
20	B(I,J)=0.0
	B(1,1)=Y(2)-Y(3)
	B(1,3)=Y(3)-Y(1)
	B(1,5)=Y(1)-Y(2)
	B(2,2)=X(3)-X(2)
	B(2,4)=X(1)-X(3)
	B(2,6)=X(2)-X(1)
	B(3,1)=B(2,2)
	B(3,2)=B(1,1)
	B(3,3)=B(2,4)
	B(3,4)=B(1,3)
	B(3,5)=B(2,6)
	B(3,6)=B(1,5)
	AR2=X(2)*Y(3)+X(3)*Y(1)+X(1)*Y(2)-X(2)*Y(1)-X(3)*Y(2)-X(1)*Y(3)
! 
! AR2=2*DLT
! 
! 计算矩阵C=(B3*6)*(D3*3)
	DO 22 I=1,6
	DO 22 J=1,3
	C(I,J)=0.0
	DO 22 K=1,3
22	C(I,J)=C(I,J)+B(K,I)*D(K,J)
! 
! 计算矩阵ESM=DTDB=CB
! 
	DO 27 I=1,6
	DO 27 J=1,6
	SUM=0.0
	DO 28 K=1,3
28	SUM=SUM+C(I,K)*B(K,J)
		ESM(I,J)=SUM*TH/(2.*AR2)
27	CONTINUE
	RETURN
	END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!输入载荷节点处的载荷值，位移边界节点处的位移值，对总体刚度矩阵、位移数组和节点力
! 数组进行修改
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE MODIFY
	 	COMMON/AV/A(8500),JGF,JGSM,NP,NBW,JEND
	 	DATA IN/60/,IO/61/
!
!输入节点的集中载荷，放到A数组中节点力向量的相应位置
!IB为施加外载荷节点的自由度(X,Y),BV为该自由度上的载荷值
202	READ(IN,*)IB	
	IF(IB.LE.0) GOTO 208
	READ(IN,*)BV
	IF(MOD(IB,2).EQ.1) GOTO 204
	WRITE(IO,203) IB/2,BV
203	FORMAT('节点',I3,'载荷为:PY=',F11.3)
	GOTO 206
204	WRITE(IO,205) IB/2+1,BV
205	FORMAT('节点',I3,'载荷为:PX=',F8.3)
206		A(JGF+IB)=A(JGF+IB)+BV
		GOTO 202

!
!输入位移边界节点处的位移值，放到A数组中节点位移向量的相对位置
!IB为节点位移的自由度，BV为位移值
!
208	READ(IN,*)IB
		IF(IB.LE.0) RETURN
		READ(IN,*)BV
		IF(MOD(IB,2).EQ.1) GOTO 214
		WRITE(IO,213) (IB+1)/2,BV
213	FORMAT('节点',I3,'位移约束为:V=',F8.3)
	GOTO 209
214	WRITE(IO,215)(IB+1)/2,BV
215	FORMAT('节点',I3,'位移约束为:U=',F8.3)
209		K=IB-1
		DO 211 J=2,NBW
			M=IB+J-1
			IF(M.GT.NP) GOTO 210
			IJ=JGSM+(J-1)*NP+IB-(J-1)*(J-2)/2
			A(JGF+M)=A(JGF+M)-A(IJ)*BV
			A(IJ)=0.0
210	IF(K.LE.0) GOTO 211
		KJ=JGSM+(J-1)*NP+K-(J-1)*(J-2)/2
		A(JGF+k)=A(JGF+K)-A(KJ)*BV
		A(KJ)=0.0
		K=K-1
211	CONTINUE
		A(JGF+IB)=A(JGSM+IB)*BV
221	CONTINUE
		GOTO 208
		END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     GUASS  消元子程序
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE DCMPBD
		COMMON /AV/A(8500),JGF,JGSM,NP,NBW,JEND
!
!用高斯消元法将对称等带宽刚度矩阵化为上三角矩阵
!
		NP1=NP-1
		DO 226 I=1,NP1
		MJ=I+NBW-1
		IF(MJ.GT.NP) MJ=NP
		NJ=I+1
		MK=NBW
		IF((NP-I+1).LT.NBW) MK=NP-I+1
		ND=0
		DO 225 J=NJ,MJ
		MK=MK-1
		ND=ND+1
		NL=ND+1
		DO 225 K=1,MK
		NK=ND+k
		JK=JGSM+(K-1)*NP+J-(K-1)*(K-2)/2
		INL=JGSM+(NL-1)*NP+I-(NL-1)*(NL-2)/2
		INK=JGSM+(NK-1)*NP+I-(NK-1)*(NK-2)/2
		II=JGSM+I
225		A(JK)=A(JK)-A(INL)*A(INK)/A(II)
226		CONTINUE
			RETURN
			END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		回代求解子程序
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE SLVBD
		COMMON /AV/A(8500),JGF,JGSM,NP,NBW,JEND
		DATA IN/60/,IO/61/
		NP1=NP-1
!
! 对子程序DEMPBD得到的上三角矩阵进行回代计算，得到节点位移向量DD
!
		DO 250 I=1,NP1
		MJ=I+NBW-1
		IF(MJ.GT.NP)MJ=NP
		NJ=I+1
		L=1
		DO 250 J=NJ,MJ
		L=L+1
		IL=JGSM+(L-1)*NP+I-(L-1)*(L-2)/2
250		A(JGF+J)=A(JGF+J)-A(IL)*A(JGF+I)/A(JGSM+I)
!
!求解节点自由度的回代计算，从下往上迭代
!
		A(NP)=A(JGF+NP)/A(JGSM+NP)
		DO 252 K=1,NP1
		I=NP-K
		MJ=NBW
		IF((I+NBW-1).GT.NP) MJ=NP-I+1
		SUM=0.0
		DO 251 J=2,MJ
		N=I+J-1
		IJ=JGSM+(J-1)*NP+I-(J-1)*(J-2)/2
251	SUM=SUM+A(IJ)*A(N)
252	A(I)=(A(JGF+I)-SUM)/A(JGSM+I)
		RETURN
	END



