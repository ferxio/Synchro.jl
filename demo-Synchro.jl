#========================================

 Modelo para construccion de problemas de
sincronizacion de SD's usando Synchro.jl

J. F. Rojas 2018
========================================#

nOsc = 100
nDim = 3

# variables de acoplamiento 1:2:3->[1 0 1]
thevars = 1:1                 # indices
vars = [0 for i in 1:nDim]'   # row vector
for i in thevars              # active variables
   vars[i] = 2
end
varnames = ["x" "y" "z"]      # names of active variables

nettype = "watts_strogatz"

F = Oscillators.Rossler

#filename = "LogisticMap"
filename = "/home/fernando/MEGAsync/Proyectos/Programas/DATA/Rossler"

rangeeps = 0.0:0.1:0.5

miumbral = 0.0    # umbral de ruid
