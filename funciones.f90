! 2020 - 25 - 5
! funciones.f90
! Félix Cabrera (walberto.cabrera@gmail.com)

! FORMA PARTE DE ELIPTICAS
! MODULO DE INTERFASE DE FUNCIÓNES

! Codificación del texto: UTF-8
! Compiladores probados: GNU Fortran (Ubuntu 9.2.1-9ubuntu2) 9.2.1 2019008

! REQUIERE:
! FX.f90

! Instrucciones de compilación:
! gfortran -Wall -pedantic -std=f95 -c FX.f90
! gfortran -Wall -pedantic -std=f95 -c funciones.f90

MODULE funciones

  ! Funcion F(x,y)
  INTERFACE
    FUNCTION Fx(xi,yi)
      IMPLICIT NONE
        REAL(8), INTENT(IN):: xi,yi
        REAL(8):: Fx
    END FUNCTION Fx
  END INTERFACE

END MODULE
