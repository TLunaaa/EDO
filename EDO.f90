Module primas

Implicit none

INTEGER, PARAMETER :: CantEc = 4

Contains

Function v_prima(v)
	Real(8), Dimension (0:CantEc) :: v_prima, v
	v_prima(0)=1.0 !Le asignamos la derivada de x, la cual es siempre 1 (REAL)
	!Recordar que x=v(0), y=v(1)
	v_prima(1)= v(2)
	v_prima(2)= (-6.10*v(1)+8.80*(v(3)-v(1))-41*v(2))/174.
	v_prima(3)= v(4)
	v_prima(4)= (-5.80*v(2)-8.80*(v(3)-v(1))-21*v(4))/132.
	!v_prima(5)= V(3)-2*V(2)-0.1*V(5)+V(1)+0.1*V(4)
	!v_prima(6)= -2.*V(3) - 0.1*V(6)+V(2)

End function

End module


PROGRAM EDOs

Use primas

 IMPLICIT NONE
    
Real(8), Dimension (0:CantEc) :: v,E=0.
Real(8):: h, Ls,tol
Logical lastIteration
Integer:: iter=0
!>>>Paremetros<<<!
	
   !Salto/Paso
	h=0.5
	
   !Limite Superior
	Ls=30.
	
   !Tolerancia
	tol=0.0001
	
   !Valores Iniciales!
	v(0)=0. !Limite inferior
	v(1)=-5.0
	v(2)= 0.
	v(3)=7.0
	v(4)=0.
	!v(5)=0.
	!v(6)=0.
!>>>>>>><<<<<<<<!
lastIteration = .false.
Open(2, file='Datos.dat', status='replace') 
Do while ((v(0) < Ls) .or. (lastIteration .eqv. .false.))
	iter=iter+1
	write(*,'(a11 i6)') "Iteracion: ",iter

	if(v(0) >= Ls) then
		lastIteration = .true.
	end if

	Call Escribe(1,v) !Escribe en la pantalla
	Call Escribe(2,v) !Escribe en el archivo

	!Call EulerSimple (v,h)
	!Call EulerModificado (v,h)
	!Call EulerMejorado (v,h)
		!Call Nuevo_h(v,h,tol)
	!Call RK4 (v,h)
		Call Nuevo_h_RKF(h,tol,v)
	Call RKF(v,h,E)
	Write(*,*) "{----------------------<<>>-----------------------}"
END DO

Close (2, status='keep')
Call system("gnuplot -persist  'Script.p'")


CONTAINS

Subroutine EulerSimple(v,h)
Real(8), Dimension (0:CantEc) :: v
real(8) h

v=v+h*v_prima(v)

End Subroutine

Subroutine EulerModificado(v,h)
Real(8), Dimension (0:CantEc) :: v, vp
real(8) h

vp=h*v_prima(v)

v= v + h*(v_prima(v) + v_prima(v+vp) )/2.0

End Subroutine

Subroutine EulerMejorado(v,h)
Real(8), Dimension (0:CantEc) :: v, vp
real(8) h

vp=h/2.*v_prima(v)

v= v + h*(v_prima(v+vp))

End Subroutine

Subroutine RK4(v,h)
Real(8), Dimension (0:CantEc) :: v, k1, k2, k3, k4
Real(8) h

k1 = h*v_prima(v)

k2 = h*v_prima(v+k1/2.0)

k3 = h*v_prima(v+k2/2.0)

k4 = h*v_prima(v+k3)

v = v + (k1+2.0*k2+2.0*k3+k4)/6.0

End Subroutine RK4

Subroutine RKF(v, h, E)
Real(8), Dimension (0:CantEc) :: v, k1, k2, k3, k4, k5, k6, E
Real(8) h

k1 = h*v_prima(v)

k2 = h*v_prima(v+k1/4.0)

k3 = h*v_prima(v+(3*k1/32.0)+(9*k2/32.0))

k4 = h*v_prima(v+(1932*k1/2197.0)-(7200*k2/2197.0)+(7296*k3/2197.0))

k5 = h*v_prima(v+(439*k1/216.0)-(8.0*k2)+(3680*k3/513.0)-(845*k4/4104.0))

k6 = h*v_prima(v-(8*k1/27.0)+(2.0*k2)-(3544*k3/2565.0)+(1859*k4/4104.0)-(11*k5/40.0))

v = v + ((25*k1/216.0)+(1408*k3/2565.0)+(2197*k4/4104.0)-(0.2*k5))

E=K1/(360.)-K3*(128./4275.)-K4*(2197./75240.)+K5/(50.)+K6*(2./55.)

End Subroutine RKF

Subroutine Escribe(Disp,v)
Real(8), Dimension (0:CantEc) :: v
Integer Disp,i
	i = 1
	if (Disp == 1) then
		Write(*,'(a2 f25.12)') "x=",v(0)
		Do while(i <= CantEc)
			Write(*,'(a1 i1 a1 f14.5)') "v",i,"=",v(i)
			i = i+1
		End Do
		Write(*,'(a2 f15.5)') "h=",h
	else
		Write(2,'(f25.12)', Advance='no') v(0)
		Do while(i <= CantEc)
			Write(2,'(f25.12)', Advance='no') v(i)
			i = i+1
		End do
		Write(2,*)
	end if
End Subroutine

SUBROUTINE Nuevo_h (v, h, tol)

!Variables
REAL(8), INTENT(INOUT) :: h
REAL(8), INTENT(IN) :: tol
REAL(8), DIMENSION (0:CantEc) :: v, VAux1, VAux2
REAL(8) Error
LOGICAL h_adecuado

	!Cuerpo
	VAux1 = v
	VAux2 = v
	h_adecuado = .false.
	DO WHILE (h_adecuado .eqv. .false.)
		CALL RK4 (VAux1, h)
		CALL RK4 (VAux2, h / 2.0) !Llamo nuevamente para obtener el mismo Y*n+1 pero ahora con 2 pasos h/2
		CALL RK4 (VAux2, h / 2.0) 
		Error = (abs(VAux1(1) )- abs(VAux2(1)))
		IF (abs(Error) > tol) THEN
				h = h / 2.0
		!ELSE IF (abs(Error) < tol / 2) THEN
			 !h = h * 1.99
			ELSE          ! Optamos por nunca agrandar el h, en el caso de contemplar el caso de agrandar el h, la tolerancia debería ser muy pequeña
				h_adecuado = .true.
		END IF
		
		VAux1 = v !Vuelvo a dejar los vectores como estaban, para una posible nueva iteración del ciclo.
		VAux2 = v
	END DO
	write(*,'(A7 F10.4)') "Error :",Error
END SUBROUTINE

SUBROUTINE Nuevo_h_RKF (h, tol, v)

!Variables
REAL(8), INTENT(INOUT) :: h
REAL(8), INTENT(IN) :: tol
REAL(8) a, aux, Error
REAL(8), DIMENSION (0:CantEc) :: E,v,vaux

	vaux = v
	Call RKF(vaux, h, E)
	Error = maxval(abs(E))  !Tomamos el mayor valor del vector
	write(*,*) E
	IF (tol >= Error) THEN
		a = 0.2
	ELSE
		a = 0.22
	END IF
	IF(Error .NE. 0.) THEN
		aux = tol / Error
		aux = abs(aux)
		aux = aux**a
		write(*,'(a6 f20.12)') "Error: ",Error
		write(*,'(a6 f20.12)') "Paso: ",h
		h = h * aux
	END IF
	
END SUBROUTINE

END PROGRAM
