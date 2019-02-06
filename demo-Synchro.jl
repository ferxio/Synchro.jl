
#include("Synchro.jl")
#include("Array4.jl")

using Synchro
using Oscillators
using LightGraphs
using Printf
using DelimitedFiles
#========================================

 Modelo para construccion de problemas de
sincronizacion de SD's 

Usando DSys3.jl con RungeKuttaMultiD() en
la solucion numerica de los SD acoplados.

========================================#

include("demo-Synchro-param.jl")

initSynchro(nOsc,nDim,vars)

### construction of the adjacency matrix
### makeMadj(nettype, directed?)
#   "one"    "random"    "all"    "watts_strogatz"
# [0 1 1 1]  [0 1 0 0]  [0 1 1 1]
# [1 0 0 0]  [1 0 1 1]  [1 0 1 1]
# [1 0 0 0]  [0 1 0 1]  [1 1 0 1]
# [1 0 0 0]  [0 1 1 0]  [1 1 1 0]
makeMadj(nettype) #,true)

### or manual construction
#M = zeros(nOsc,nOsc)
#M[2,1] = 1
#M[3,1] = 1
#makeMadj(M) or 
#makeMadj(M,true) for directed network

println("M---",Synchro.Madj,"===") #,Madj)
writedlm("$filename-$nettype-Madj-$nSD",Synchro.Madj)

# the oscillators...
F = Rossler

# construccion de vectores: el de SD y el de funciones
MU = [[rand(-10.:10.) rand(-10.:10.) rand(-10.:10.)] for i in 1:nOsc]   # condiciones iniciales
MF = [F for i in 1:nOsc]                           # funciones que definen al oscilador F
#MF = [ODef[rand(1:4)] for i in 1:nOsc]  # lista osciladores predefinidos (aleatorio)

println("definidos ",MU," ",MF)

# la solucion!!!
for mieps in rangeeps
   # define Eps=peso uniforme de terminos de acoplo, y umbral de ruido
   setEpsWeight(mieps)
   setNoiseTop(miumbral)
   setPhaseDist(midist)
   # resuelve todos los osciladores con cond inicial MU definidos por MF guardando cada 1 datos
   SolveRK(MU,MF,npts,0.01,"$filename-$nettype-$nOsc-$mieps-$miumbral","dat",each) 
   # distribucion de distancias
   println("---------------------------------------------------------")
   println(distances_list())
   println(distances_dist())
   println("---------------------------------------------------------")
   writedlm("$filename-$nettype-distancesList-$nOsc-$mieps-$miumbral.dat",distances_list())
   writedlm("$filename-$nettype-distancesHist-$nOsc-$mieps-$miumbral.dat",distances_dist())
   # distribucion de amplitudes
   writedlm("$filename-$nettype-amplitudesList-$nOsc-$mieps-$miumbral.dat",amplitudes_list())
   writedlm("$filename-$nettype-amplitudesHist-$nOsc-$mieps-$miumbral.dat",amplitudes_dist())
end  # for
mean_shortest_path()
println("it has been solved..!")
