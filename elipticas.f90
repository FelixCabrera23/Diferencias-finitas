! 2020 - 25 - 5
! penduloD2.f90
! Félix Cabrera (walberto.cabrera@gmail.com)

! SOLUCION A LA ECUACIÓN DE POISSON
! POR DIFERENCIAS FINITAS

! Codificación del texto: UTF-8
! Compiladores probados: GNU Fortran (Ubuntu 9.2.1-9ubuntu2) 9.2.1 2019008

! Algoritomo tomado del libro de Analisis Numerico de Burden
! Algoritmo 12.1 pgs 555-556

! Rquiere: Fx.f90, funciones.f90
! Instrucciones de compilacion:
! gfortran -Wall -pedantic -std=f95 -c FX.f90
! gfortran -Wall -pedantic -std=f95 -c GX.f90
! gfortran -Wall -pedantic -std=f95 -c funciones.f90
! gfortran -Wall -pedantic -std=f95 -c elipticas.f90
! gfortran -Wall -pedantic -std=f95 -o poisson elipticas.o funciones.o FX.o GX.o
! ./poisson

PROGRAM poisson
  USE funciones
  IMPLICIT NONE
  
  ! Definimos variables principales
  INTEGER(8):: m, n ! Puntos de frontera y pasos
  REAL(8):: a,b,c,d,tol
  REAL(8), ALLOCATABLE:: x(:),y(:),w(:,:)
  ! Definimos variables auxiliares
  INTEGER(8):: i,j,l, pasos
  REAL(8):: lambda,mu,z,norm,h,k
  ! Variables de control
  INTEGER(4) :: err
  
  ! Abrimos el archivo de configuración (12)
  OPEN (12, FILE='parametros.config', STATUS='old', IOSTAT=err)
  IF (err .ne. 0) STOP 'parametros.config is missing'

  READ(12,*) a
  READ(12,*) b
  READ(12,*) c
  READ(12,*) d
  READ(12,*) n
  READ(12,*) m
  READ(12,*) pasos
  READ(12,*) tol
  CLOSE(12)
  
  ! Guardamos espacio en memoria
  ALLOCATE(x(n),y(m),w(n,m))
  
  ! Abrimos el acribo para los resultados
  OPEN(14, FILE='Resultados.dat', STATUS='new', IOSTAT=err)
  IF (err .ne. 0) STOP 'Resultados.dat already exists'

  PRINT *, '***********************************************************************'
  PRINT *, 'Resolviendo la ecuación de Poisson'
  PRINT *, 'Método= differencias finitas'
  PRINT *, 'Tolerancia =',tol,', Pasos =',pasos
  PRINT *, 'Procesando...'



  ! paso 1
  h = (b-a)/REAL(n,8)
  k = (d-c)/REAL(m,8)
  

  ! paso 2
  DO i=1, n-1
    x(i)=a+i*h
  END DO
  ! paso 3
  DO j=1, m-1
    y(j)=c+j*k
  END DO
  ! Preparamos las variables
  ! paso 4
  w=0
  ! paso 5
  lambda = (h**2)/(k**2)
  mu = 2*(1+lambda)
  l = 1
  
  ! paso 6
  DO WHILE (l<=pasos) ! para los pasos 7-20
    ! paso 7
    ! Iniciamos las variables
    z=(-(h**2)*Fx(x(1),y(m-1))+Gx(a,y(m-1))+lambda*Gx(x(1),d)+lambda*w(1,m-2)+w(2,m-1))/mu
    
    norm=ABS(z-w(1,m-1))
   
    w(1,m-1)=z
    
    ! paso 8
    DO i=2,n-2
      z=(-(h**2)*Fx(x(i),y(m-1))+lambda*Gx(x(i),d)+w(i-1,m-1)+w(i+1,m-1)+lambda*w(i,m-2))/mu
      IF (ABS(w(i,m-1)-z)>norm) norm=ABS(w(i,m-1)-z)
      w(i,m-1)=z
    END DO
    
    ! paso 9
    z=(-(h**2)*Fx(x(n-1),y(m-1))+Gx(b,y(m-1))+lambda*Gx(x(n-1),d)+w(n-2,m-1)+lambda*w(n-1,m-2))/mu
    IF(ABS(w(n-1,m-1)-z)>norm) norm=ABS(w(n-1,m-1)-z)
    w(n-1,m-1) = z
    
    ! paso 10
    j=m-2
    DO WHILE (j>=2) ! pasos 11,12 y 13
      ! paso 11
      z=(-(h**2)*Fx(x(1),y(j))+Gx(a,y(j))+lambda*w(1,j+1)+lambda*w(1,j-1)+w(2,j))/mu
      IF (ABS(w(1,j)-z)>norm) norm=ABS(w(1,j)-z)
      w(1,j)=z
      
      ! paso 12
      DO i=2, n-2
        z=(-(h**2)*Fx(x(i),y(j))+w(i-1,j)+lambda*w(i,j+1)+w(i+1,j)+lambda*w(i,j-1))/mu
        IF (ABS(w(i,j)-z)>norm) norm = ABS(w(i,j)-z)
        w(i,j)=z
      END DO
      
      ! paso 13
      z = (-(h**2)*Fx(x(n-1),y(j))+Gx(b,y(j))+w(n-2,j)+lambda*w(n-1,j+1)+lambda*w(n-1,j-1))/mu
      IF (ABS(w(n-1,j)-z)>norm) norm = ABS(w(n-1,j)-z)
      w(n-1,j)=z
      
      j=j-1
    END DO
    ! paso 14
    z = (-(h**2)*Fx(x(1),y(1))+Gx(a,y(1))+lambda*Gx(x(1),c)+lambda*w(1,2)+w(2,1))/mu
    IF (ABS(w(1,1)-z)>norm) norm = ABS(w(1,1)-z)
    w(1,1)=z
    
    ! paso 15
    DO i=2, n-2
      z=(-(h**2)*Fx(x(i),y(1))+lambda*Gx(x(i),c)+w(i-1,1)+lambda*w(i,2)+w(i+1,1))/mu
      IF (ABS(w(i,1)-z)>norm) norm=ABS(w(i,1)-z)
      w(i,1)=z
    END DO
    
    ! paso 16
    z=(-(h**2)*Fx(x(n-1),y(1))+Gx(b,y(1))+lambda*Gx(x(n-1),c)+w(n-2,1)+lambda*w(n-1,2))/mu
    IF (ABS(w(n-1,1)-z)>norm) norm=ABS(w(n-1,1)-z)
    w(n-1,1)=z
    
    ! paso 17
    IF (norm<=tol) THEN
      ! paso 18
      DO i=1, n-1
        DO j=1, m-1
          WRITE (12,*) i,j,x(i),y(j),w(i,j)
        END DO
      END DO
      !paso 19
      PRINT *, 'Calculado con exito despues de',l,'iteraciones'
      DEALLOCATE (x,y,w)
      PRINT *, '***********************************************************************'
      STOP
    END IF
    ! paso 20
    IF (MOD(l,(pasos/20))==0) THEN
      WRITE(*, fmt='(1x,a,i0)', advance='no') '%'
    END IF
    
    l=l+1
  END DO
  ! paso 21
  PRINT *, 'El proceso fallo despues de',l-1,'iteraciones'
  PRINT *, '***********************************************************************'
  DEALLOCATE (x,y,w)


END PROGRAM poisson













