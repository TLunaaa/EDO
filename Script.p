# Ejemplo de script Gnuplot, para graficar los datos 
# que se encuentran en un archivo llamado "ejercicio_7.dat"
# Este archivo se llama   ejemploScript.p
set   autoscale                        # escala los ejes automaticamente
unset log                              # quita la escala logaritmica (si la hubiera)
unset label                            # quita los titulos anteriores
set xtic auto                          # establece automaticamente las divisiones del eje x
set ytic auto                          # establece automaticamente las divisiones del eje y
set grid
set title "EDO - PVI"
set xlabel "x"
set ylabel "y"

#linespoints
#lines
#points

plot  "Datos.dat" using 1:2 title 'Y1 aprox' with lines,\
      "Datos.dat" using 1:4 title 'Y2 aprox' with lines  

