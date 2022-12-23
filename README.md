# Kepler-Orbits
Study of Kepler orbits using Newton Rhapson method (numerical methods). The script is separetd in differents sections.

Johannes Kepler could determine the orbits of bodies revolving around the Sun. In this program we are gonna be able to compute different and interesting parameters to put on practice the application of this method.

## 1 Computation of the distance to the origin


$$x(E) = a(\cos(E) -\epsilon)$$  $$y(E)=a\sqrt{1-\epsilon^2} \sin (E)$$


With $a=189.857$ U.A and $\epsilon= 0.995086$.

Where we can compute the distance with $D(E)=\sqrt{x^2+y^2}$ and its derivative.

## 2 Computation of the 0's of the distance function via bisection method.

$$F(E)=\sin(2E)(1-\epsilon^2 -(cos(E)(e-\epsilon^2)-\epsilon)\sin(E)$$

Where we are able to find the maximum energy and maximum distance.

## 3 Abnormal eccentricity

$$E= \frac{2\pi}{T_H}-\epsilon \sin(E)$$

With $E_0=\pi/6$ as the initial value, we are able to compute 100 values of the abnormal eccentricity via Newton-Rhapson method in one full orbit of the comet arround the sun.

The period of the comet is $T_H= 2526.5$ years.

## 4 Study of the convergence of Newton-Rhapson method.

Starting with different values of $E_0= 0.4,1.7,2.3,3.6,6.1$

# Compilation and execution of the program
This program was written in _Fortran_ 77 and the graphics were plotted with _Gnuplot_.
## Linux and Mac
### Compilation

```
gfortran -name_of_the_file.f -o name_of_the_output_file.out
```
### Execution
```
./name_of_the_output_file.out
```

## Windows
### Compilation
```
gfortran -name_of_the_file.f -o name_of_the_output_file.exe
```
### Execution
```
./name_of_the_output_file.exe
```
