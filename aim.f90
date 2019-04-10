PROGRAM IQA
!
INTEGER :: natom, inte, i
!
REAL :: Energy, Error
!
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: TA, VneAA, VeeCAA, VeeXAA
!
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: VneAB, VenAB, VnnAB, VeeCAB, VeeXAB
!
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: T_self, Vc_self, Vx_self
!
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Vc_int, Vx_int
!
DOUBLE PRECISION :: AA_T, AA_Vc, AA_Vx, AB_Vc, AB_Vx, Tot_Vc, Tot_Vx, Tot_IQA
!
OPEN(UNIT=1, FILE='AA',STATUS='OLD')
OPEN(UNIT=2, FILE='AB',STATUS='OLD')
OPEN(UNIT=3, FILE='out',STATUS='NEW')
OPEN(UNIT=5, FILE='numer', STATUS='OLD')
OPEN(UNIT=6, FILE='wfn', STATUS='OLD')
!
READ(5,*) natom
!
READ(6,*) Energy
!
inte=((natom*(natom+1))/2)-natom
!
ALLOCATE(TA(natom))
ALLOCATE(VneAA(natom))
ALLOCATE(VeeCAA(natom))
ALLOCATE(VeeXAA(natom))
!
ALLOCATE(VneAB(inte))
ALLOCATE(VenAB(inte))
ALLOCATE(VnnAB(inte))
ALLOCATE(VeeCAB(inte))
ALLOCATE(VeeXAB(inte))
!
ALLOCATE(T_self(natom))
ALLOCATE(Vc_self(natom))
ALLOCATE(Vx_self(natom))
!
ALLOCATE(Vc_int(inte))
ALLOCATE(Vx_int(inte))
!
	DO i=1, natom
	READ(1,*) TA(i), VneAA(i), VeeCAA(i), VeeXAA(i)
	END DO
!
	DO i=1, inte
	READ(2,*) VneAB(i), VenAB(i), VnnAB(i), VeeCAB(i), VeeXAB(i)
	END DO
!
	DO i=1, natom
	T_self(i)=TA(i)
	Vc_self(i)=VneAA(i)+VeeCAA(i)
	Vx_self(i)=VeeXAA(i)
	END DO
!
	DO i=1, inte
	Vc_int(i)=(VneAB(i)+VenAB(i)+VnnAB(i)+VeeCAB(i))*2
	Vx_int(i)=VeeXAB(i)*2
	END DO
!
AA_T=sum(T_self)
AA_Vc=sum(Vc_self)
AA_Vx=sum(Vx_self)
!
AB_Vc=sum(Vc_int)
AB_Vx=sum(Vx_int)
!
Tot_Vc=AA_Vc+AB_Vc
Tot_Vx=AA_Vx+AB_Vx
!
Tot_IQA=AA_T+Tot_Vc+Tot_Vx
!
Error=(Energy-Tot_IQA)*2625.5
!
	DO i=1, natom
	 WRITE(3,*) "Atom ", i 
	 WRITE(3,*) "Kinetic= ", T_self(i)
	 WRITE(3,*) "Coulomb= ", Vc_self(i)
	 WRITE(3,*) "Exchange= ", Vx_self(i)
	END DO
!
	DO i=1, inte
	 WRITE(3,*) "interaction ", i
	 WRITE(3,*) "Coulomb= ", Vc_int(i)
	 WRITE(3,*) "Exchange= ", Vx_int(i)
	END DO
!
WRITE(3,*) "Total Kinetic (self) = ", AA_T 
WRITE(3,*) "Total self Coulomb = ", AA_Vc 
WRITE(3,*) "Total self exchange = ", AA_Vx
!
WRITE(3,*) "Total BOND/INT Coulomb = ", AB_Vc 
WRITE(3,*) "Total BOND/INT exchange = ", AB_Vx
!
WRITE(3,*) "Total Kinetic = ",  AA_T
WRITE(3,*) "Total Coulomb = ",  Tot_Vc 
WRITE(3,*) "Total exchange = ", Tot_vx  
WRITE(3,*) "Total IQA (except correlation) = ", Tot_IQA
WRITE(3,*) "The Gaussian energy is = ", Energy
WRITE(3,*) "The integration error is", Error
!
END PROGRAM IQA
