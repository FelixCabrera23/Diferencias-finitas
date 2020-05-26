! 2020 - 25 - 5
! FX.f90
! Félix Cabrera (walberto.cabrera@gmail.com)

! FORMA PARTE ELIPTICAS

! Codificación del texto: UTF-8
! Compiladores probados: GNU Fortran (Ubuntu 9.2.1-9ubuntu2) 9.2.1 2019008
! Instrucciones de compilación: 
! gfortran -Wall -pedantic -std=f95 -c FX.f90

FUNCTION Fx(xi,yi)
  IMPLICIT NONE
  
  ! Definimos variables de entrada
  REAL(8), INTENT(IN):: xi,yi
  REAL(8):: Fx
  
  Fx= xi*EXP(yi)

END FUNCTION
