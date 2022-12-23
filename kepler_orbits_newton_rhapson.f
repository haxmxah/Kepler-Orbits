c		Marta Xiulan Aribó Herrera
c 		Mètode de Newton-Raphson i órbitas de Kepler
C 		2019
c-----------------------------------------------------------------------
		
		PROGRAM orbites_kepler
		implicit none
		integer ndat,i
		parameter (ndat = 120)
		double precision x,y,a0,e0,E,D(ndat),pi,E_vector(ndat),dD(ndat)
		parameter (pi = 4.d0* atan(1.d0))
		double precision F,A,B,eps,emax
		integer niter,ndat1
		double precision t,xarrel
		double precision E_0(5)
		external fun
		external fun1
c---------------------------- APARTAT 1 --------------------------------
c		Calcul de la distancia a l'origen de coordenadas i la 
C		seva derivada numèrica amb la subrutina derfun.
c-----------------------------------------------------------------------
		E = 0.d0
		a0 = 189.857d0
		e0 = 0.995086d0

		do i = 1, ndat
			x = a0 * (cos(E)-e0)
			y = a0*sqrt(1.d0-e0**2.d0)*sin(E)
			D(i) = sqrt(x**2.d0+y**2.d0)
			E_vector(i) = E
			E = E + (2.d0*pi/ndat)
		enddo

		call derfun(ndat,E_vector,D,dD)

		open (1,file= "P3-20-21-res.dat")
		do i = 1, ndat
			write(1,200) E_vector(i), D(i), dD(i)
		enddo
		write (1,"(a1)")
		write (1,"(a1)")
200 	format(f20.12,2x,f20.12,2x,f20.12)

		call system ("gnuplot -p plot1.gnu")

C--------------------------- APARTAT 2 ---------------------------------
c		Calcul de les arrels de la funció amb l'algorisme 
C		bisecció que donen els valors extrems.
c-----------------------------------------------------------------------

		A = 0.2d0
		B = 6.1d0
		eps = 1.d-9
		call Bisection(A,B,eps,fun,niter,emax)
		x = a0 * (cos(emax)-e0)
		y = a0*sqrt(1.d0-e0**2.d0)*sin(emax)
		
		write (1, "(a14,5x,a14)") "#Emax", "distància"
		write (1,"(f20.12,2x,f20.12)") emax, sqrt(x**2.d0+y**2.d0)
		write (1,"(a1)")
		write (1,"(a1)")

C--------------------------- APARTAT 3 ---------------------------------
C		Calcul de l'excentricité anomal ami el métope de Newton
C 		Rhapson.
c-----------------------------------------------------------------------
		eps = 1.d-11
		ndat1 = nint(2526.5d0/100)+1
		t = 0.d0
		write (1,"(a20)") "#Per NewtonRap tenim"
		do i = 1, 100
			E = (2.d0*pi/2526.5d0)*t - (e0*sin(pi/6.d0))
!			call NewtonRap(E,eps,fun,niter,xarrel)
			x = a0*(cos(E)-e0)
			y = a0*sqrt(1.d0-e0**2.d0)*sin(E)
!			write (1, "(a20)") "per t=", t
			write (1,400) t,E,x,y
			t = t + (2526.5d0/100.d0)
		enddo

400		format(f20.12,2x,f20.12,2x,f20.12,2x,f20.12)
		write (1,"(a1)")
		write (1,"(a1)")
		call system ("gnuplot -p plot2.gnu")

C--------------------------- EXTRA APARTAT 3 ---------------------------
C		Estudi de la convergència del metode de Newton-Raphson.
c-----------------------------------------------------------------------

		Write (1, "(a25)") "# Convergència Newtonrap"
		E_0(1) = 0.4d0
		E_0(2) = 1.7d0
		E_0(3) = 2.3d0
		E_0(4) = 3.6d0
		E_0(5) = 6.1d0
		call NewtonRap(E,eps,fun,niter,xarrel)
	
		do i = 1,5
			t = 3.d0*2526.5d0/5.d0
			E = (2.d0*pi/2526.5d0)*t - (e0*sin(E_0(i)))
			call NewtonRap(E,eps,fun,niter,xarrel)
			write (1,500) E_0(i), niter
		enddo
500		format(f20.12,2x,i4)
		END PROGRAM


C--------------------------- SUBROUTINES ---------------------------------

		subroutine fun (x,fu,dfu) 
		implicit none
		double precision x, fu, dfu, dfu1, dfu2, dfu3, fu1,fu2
		double precision e0
		e0 = 0.995086d0
		fu1 =sin(2.d0*x)*(1.d0-e0**2.d0)
		fu2 =(cos(x)*(2.d0-e0**2.d0)-e0)*sin(x)
		fu=fu1-fu2
		dfu1 =2.d0*(1.d0-e0**2)*cos(2.d0*x)
		dfu2=cos(x)*((cos(x)*(2.d0-e0**2.d0))-e0)
		dfu3=sin(x)*(sin(x)*(2.d0-e0**2.d0))
		dfu = dfu1 - dfu2-dfu3
		return
		end

		subroutine fun1 (x,fu,dfu) 
		implicit none
		double precision x, fu, dfu, dfu1, dfu2, dfu3, fu1,fu2
		double precision e0,th,e00,pi,e,de
		parameter (pi =4.0*atan(1.d0))
		e0 = 0.995086d0
		e00 = pi/6.d0
		th = 2526.5d0
		e = (2.d0*pi/th)*x - (e0*sin(e00*x))
		de = (2.d0*pi/th) -(e0*e00*cos(e00*x))

		fu1 =sin(2.d0*e)*(1.d0-e0**2.d0)
		fu2 =(cos(e)*(2.d0-e0**2.d0)-e0)*sin(e)
		fu=fu1-fu2
		dfu1 =2.d0*de*(1.d0-e0**2)*cos(2.d0*x)
		dfu2=cos(e)*de*((cos(x)*(2.d0-e0**2.d0))-e0)
		dfu3=sin(e)*(sin(e)*de*(2.d0-e0**2.d0))
		dfu = dfu1 - dfu2-dfu3
		return
		end

				subroutine NewtonRap(x0,eps,fun,niter,xarrel)
		implicit none
		external fun
		double precision x0, eps, xarrel, x
		integer niter, i, maxiter
		double precision fu, dfu, diff
		maxiter = 2000
		x = x0
!		open (1, file = "P3-20-21-res.dat")
!		write (1, "(a13,2x,f20.12)") "# Pel punt t", x0
		do i = 1, maxiter
			niter = i
			call fun1(x,fu,dfu)
!			print*,"hehe", x,fu,dfu
			xarrel = x - (fu/dfu)
!			print*,"xarrel", xarrel
			diff = abs(fu/dfu)
!			x = a0 * (cos(fu)-e0)
!			y = a0*sqrt(1.d0-e0**2.d0)*sin(fu)
!			write (2,"(i4,f20.12)") x, fu, x,y 
			if (diff.LE.eps) then
				print *, "S'ha arribat a la precisio desitjada"
				print *, "El valor de c =", xarrel
				EXIT
			endif
			x = xarrel
		enddo
!		write (1, "(a1)") !espais perquè es reconeixin al gnuplot
!		write (1, "(a1)") 
		return
		end subroutine

		subroutine Bisection(A,B,eps,fun,niter,xarrel)
		implicit none
		double precision A,B, C,eps,xarrel
		double precision fa,dfa, fb,dfb,fc,dfc, diff !funció fun
		integer niter, i

		niter = nint (log((abs(B)-abs(A))/eps)/log(2.d0)) + 1 !iteracions

		Do i = 1, niter
			C = (A+B)/2.d0

			call fun (A,fa,dfa)
			call fun (B,fb,dfb)
			call fun (C,fc,dfc)

			if (fc.eq.0.d0) then
				xarrel = c
				print*, "la solució exacta es x =", xarrel
				exit
			endif
			if (fa*fb.GE.0.d0) then
				print*, "La funció no té canvis de signes en l'interval"
				print*, A,B,fa,fb
				exit
			endif
			if (fa*fc.LT.0.d0) then
				B = C
			else 
				A = C
			endif
			diff = B-A
			if (diff.LE.eps) then
				xarrel = c
				print*, "solució aproximada x =", c
				print*, "error <", diff
				exit
			endif
		enddo
		return
		end

		subroutine derfun(ndat,x,fu,dfu)
		implicit none 
		integer ndat, i
		double precision x(ndat), fu(ndat), dfu(ndat), h
		h = x(ndat)-(x(ndat-1)) !intervals constants

		do i = 1, ndat
			if ((i.LT. ndat).and.(i.GT.1)) then
				dfu(i)= (fu(i+1)-fu(i-1))/(2*h)
			endif
			if (i.EQ.1) then
				dfu(i) = (fu(i+1)- fu(i))/h
			endif
			if (i.EQ.ndat) then
				dfu(i) = (fu(i)-fu(i-1))/h
			endif
		enddo
		return
		end

