module Synchro

using LightGraphs
using Printf
using DelimitedFiles
using Array4
using Oscillators

#===================================================
  Modulo para resolver sistemas dinamicos en Julia
  para estudios de sincronizacion

  mayo 2017

  Cada array dinamico i tiene vinculo con varios SD
  j. El acoplamiento tendra terminos agregados de la
  forma (x_i-x_j)
  
  necesito agregar una funcion Coupling() que use la
  matriz de adyacencia
====================================================#

   export initSynchro
   export distances
   export distances_list
   export distances_dist
   export amplitudes_list
   export amplitudes_dist
   export SolveRK
   export SolveMap
   export makeMadj
   export mean_shortest_path
   export setEpsWeight
   export setNoiseTop
   export setPhaseDist
   export ODef, MDef
   export nSD, sSD, Madj, Mlap, Eps, Meps, Rnd, dimSD, U, gammaU
   export signal
   nSD = 1             # numer of DS's
   #export nSD 
   dimSD = 1         # dimension by default
   Eps = [0.0]        # weight of coupling
   Rnd = 0.0          # umbral of noise
   gammaU = 0.0      # distancia en espacio fase
   sSD = [1]
   Madj = zeros(nSD,nSD)
   MEps = ones(nSD,nSD).*Eps
   U = [zeros(dimSD) for i in 1:nSD]
   ## abajo ODef = [Lorenz Rossler Sprott1 Sprott2 Sprott3 Sprott4 Sprott5]
   
   ODef = [Lorenz Rossler Sprott1 Sprott2 Sprott3 Sprott4 Sprott5]
   MDef = [Logistic]
   
   function initSynchro(minSD::Int,midim::Int,misSD::Any)
      global nSD,sSD,dimSD,Madj,Eps,Meps
      nSD = minSD
      dimSD = midim
      sSD = misSD
      Madj = zeros(nSD,nSD)
      #Meps = ones(nSD,nSD)*Eps
      println("init nsd ",nSD," --- ",Synchro.nSD,"--",dimSD)
      #println("synchro ",Synchro.nSD)
   end
   
   function initSynchro(miU::Any,miF::Any,misSD::Any)
      global U,nSD,sSD,dimSD,Madj,Eps,Meps
      U = deepcopy(miU)
      F = miF
      nSD = size(U)[1]     # numero de SD's acoplados
      dimSD = size(U[1])[2]
      sSD = misSD
      Madj = zeros(nSD,nSD)
      #Meps = ones(nSD,nSD)*Eps
      println("init nsd ",nSD," --- ",Synchro.nSD,"--",dimSD)
      #println("synchro ",Synchro.nSD)
   end
   
   function dij(u,v) # distancia euclidiana entre dos SD's
      return sqrt(sum((u.-v).*(u.-v)))
   end

   function distances()          # matriz de distancias
      global U, nSD
      D = zeros(nSD,nSD)
      for i in 1:nSD
         for j in i+1:nSD  # derecha de la diagonal
            D[i,j] = dij(U[i],U[j])
         end
      end
      return D
   end
   
   function distances_list()     # distances list
      global U, nSD
      damp = Float32[]
      for i in 1:nSD-1
         for j in i+1:nSD  # derecha de la diagonal
            push!(damp,dij(U[i],U[j]))
         end
      end
      return damp
   end
   
   function distances_dist()     # distances distribution
      dl = distances_list()
      x,y = Histogram(dl)
      return hcat(x,y)
   end
   
   function amplitudes_list()    # distances list
      global U, nSD
      amp = Float32[]
      for i in 1:nSD
         push!(amp,norm(U[i]))
      end
      return amp
   end
   
   function amplitudes_dist()    # amplitudes distribution
      da = amplitudes_list()
      x,y = Histogram(da)
      return hcat(x,y)
   end
   
   function Theta(x)
      if x < 0.0 
         return 0.0 
      end
      return 1.0
   end
      
   function r(gamma) # fracción de pares más cercanos que gamma
      global nSD,U
      suma = 0.0
      N1 = nSD*(nSD-1)
      for i in 1:nSD-1
         for j in i+1:nSD
            suma += Theta(gamma-dij(U[i],U[j]))
         end
      end
      return 2*suma/N1
   end

   function s_(gamma) # fracción de elementos alejados: d_ij>gamma
      global nSD,U
      dijs = []
      suma = 0.0
      for i in 1:nSD
         prod = 1.0
         for j in 1:nSD
            if j != i
               d = dij(U[i],U[j])
               prod = prod*Theta(d - gamma) # cuenta los i que tienen todos alejados
               append!(dijs,d)
               #print("prod ",i," ",j," ",prod,"  ")
            end  #  if
         end  #for
         suma += prod
         #println("suma ",suma)
      end  # for
      return suma/nSD
   end

   function s(gamma) # fracción con al menos un elemento cercano: d_ij<gamma
      #print("s_",s_(gamma),"  ")
      return 1.0-s_(gamma)
   end
   
   #===== SOLUCION NUMERICA =====#
   function makeMadj(nettype::String,directed=false)
      global Madj,Mlap,nSD,Meps
      #   "one"    "random"    "all"    "watts_strogatz"
      # [0 1 1 1]  [0 1 0 0]  [0 1 1 1]
      # [1 0 0 0]  [1 0 1 1]  [1 0 1 1]
      # [1 0 0 0]  [0 1 0 1]  [1 1 0 1]
      # [1 0 0 0]  [0 1 1 0]  [1 1 1 0]
      if nettype == "random"
         Madj = rand(0:1,nSD,nSD)
      elseif nettype == "all"
         Madj = ones(nSD,nSD).-eye(nSD)
      elseif nettype == "one"
         Madj = zeros(nSD,nSD)
         Madj[2:end,1] = 1
      elseif nettype == "watts_strogatz"
         g = watts_strogatz(nSD,div(nSD,4),0.5)
         Madj = zeros(nSD,nSD) + adjacency_matrix(g)
      end
      if !directed   # se hace simetrica
         for i in 1:nSD
            for j in i:nSD
               Madj[i,j] = 0
               Madj[i,j] = Madj[j,i]
            end  # for
         end  # for
      end  # if
      Mlap = deepcopy(Madj)
      for i in 1:nSD
         Mlap[i,i] = -sum(Madj[i,:])
      end # for
      #return M
   end

   function makeMadj(miM::Array,directed=false)
      global Madj,Mlap,nSD,Meps
      Madj = deepcopy(miM)
      if !directed
         for i in 1:nSD
            for j in i:nSD
               Madj[i,j] = 0
               Madj[i,j] = Madj[j,i]
            end  # for
         end  # for
      end  # if
      Mlap = deepcopy(Madj)
      for i in 1:nSD
         Mlap[i,i] = -sum(Madj[i,:])
      end # for
      println(" local Madj ",Madj)
      #return M
   end
   
   function mean_shortest_path()
      global Madj,nSD
      g = Graph(Madj)
      CL = [ sum(dijkstra_shortest_paths(g,i).dists) for i in vertices(g) ]
      CL = sum(CL)/nSD/(nSD-1)
      println("mean shortest path: ",CL)
      return CL
   end
   
   function signal()
      #=========================
      la señal: cómo la elegimos
      =========================#
      global U,nSD,sSD,nSD
      l = div(nSD,5)
      rosc = l:l:nSD
      #println("kepex...",U[rosc],sSD,sum(U[rosc]))
      #return dot(sum(U[rosc])/length(U[rosc]),float(sSD))
      return dot(sum(U),float(sSD'))
   end
   
   function setEpsWeight(miEps)
      global Eps
      Eps = miEps
   end
   
   function setNoiseTop(miNoise)
      global Rnd
      Rnd = miNoise
   end
   
   function setPhaseDist(miPhaseDist)
      global gammaU
      gammaU = miPhaseDist
   end
      
   function makeEpsRnd(valE,valR)
      global Eps,Rnd
      Eps,Rnd = valE,valR
   end
   
   function coupling(U,i)
      #====================================
      evaluacion de termino de acoplamiento
      para el oscilador U[i]
      ====================================#
      global Madj,nSD,sSD,dimSD,Meps,Eps
      tadd = [0.0 for j in 1:dimSD]
      tadd = tadd'
      for j in 1:nSD
         if i!=j
            #println("tadd eps i j M_ij Uj-Ui sSD: ",tadd," / ",Eps," / ",i," ",j," / ",Madj[i,j]," / ",U[j].-U[i]," / ",sSD)
            tadd = tadd .+ Eps*Madj[i,j]*(U[j].-U[i]).*sSD  # Complete Synchronisation?
         end  # if
      end  # for
      #println("tadd final: ",tadd .+ [rand(-Rnd:Rnd) for i in 1:dimSD]')
      return tadd .+ [rand(-Rnd:Rnd) for i in 1:dimSD]'
   end
 
   function SolveMap(X, F, n=50000, name="datosMAP", ext="dat", every=1)
   global U,nSD,dimSD,sSD,gammaU
   #====================================================================
   Solucion numerica de un sistema de SD's discretos (mapeos)
   (al estilo de los Sistemas Dinamicos; pensando en SINCRONIZACION)
   Requiere:
      - un array de SD's (con sus array de parámetros)
      - un array de funciones tipo F(U) que se asocian a cada SD (1 a 1)
       
      Parámetros:  X, el arreglo de valores iniciales de los SD
      F, el array de nombres de funciones de los SD
      n, número de puntos a evaluar
      arch='datos.dat', el archivo de datos
   ====================================================================#
       
      U = deepcopy(X)
      Ulocal = deepcopy(X)
      # se guarda la condición inicial
      ff=open("$name.$ext","w")
      println(ff,"# Synchro.jl Mapeos acoplados")
      println(ff,"# $nSD sistemas (mapeos) de dimension $dimSD:")
      print(ff,"# t ")
      gg=open("$name-tRS","w")
      println(gg,"# Synchro.jl Parametros de Orden t, r(t) y s(t)")
      hh = open("$name-SG.dat","w")
      for sd in 1:nSD
         for i in 1:dimSD
            print(ff,"X$sd[$i] ") 
         end
      end
      println(ff)
      # aplicacion directa del mapeo
      for t in 0:n
         if t % every == 0 
            writeMultiM(ff, U, t) 
            println(gg,t," ",r(gammaU)," ",s(gammaU))
            println(hh,signal())
            #println("U-RK ",U)
         end
         for sd in 1:nSD
            #println("Mapping... ", U[sd], F[sd](Ulocal[sd]), coupling(Ulocal,sd))
            U[sd] = F[sd](Ulocal[sd]) .+ coupling(Ulocal,sd)
            U[sd] = (U[sd] .> 1.0).*1.0 .+ (U[sd] .<= 1.0).*U[sd]
            U[sd] = U[sd] .- (U[sd] .< 0.0).*U[sd]
         end # for sd
         Ulocal = deepcopy(U)
      end # for t
      close(ff)
      close(gg)
      close(hh)
      ###
   end
   
   function SolveRK(X, F, n=50000, h=0.01, name="datosRK", ext="dat", every=10)
      global U,nSD,dimSD,sSD,gammaU
   #====================================================================
   Solucion numerica de un sistema de SD's
   (al estilo de los Sistemas Dinamicos; pensando en SINCRONIZACION)
   Requiere:
      - un array de SD's (con sus array de parámetros)
      - un array de funciones tipo F(U) que se asocian a cada SD (1 a 1)
   Usando el método de Runge-Kutta de Cuarto Orden:
       
      Parámetros:  X, el arreglo de valores iniciales de los SD
      F, el array de nombres de funciones de los SD
      n, número de puntos a evaluar
      h=0.0001, precisión (tamaño de paso)
      arch='datos.dat', el archivo de datos
   =====================================================================#
       
      U = deepcopy(X)
      Ulocal = deepcopy(X)
      #nSD = size(U)[1]     # numero de SD's acoplados
      #dimSD = size(U[1])[2]  # dimension (variables) de un SD
      #ndimSD = ndims(U)
      
      # se guarda la condicin inicial
      println("creando archivo $name.$ext")
      ff=open("$name.$ext","w")
      gg=open("$name-tRS","w")
      hh = open("$name-SG.dat","w")
      println(ff,"# Synchro.jl Solucion por metodo de Runge-Kutta de orden 4")
      println(ff,"# $nSD sistemas, de dimension $dimSD:")
      println(gg,"# Synchro.jl Parametros de Orden r(t) y s(t)")
      print(ff,"# t ")
      for sd in 1:nSD
         for i in 1:dimSD
            print(ff,"X$sd[$i] ") 
         end
      end
      println(ff)
      # método de Runge-Kutta a cada SD
      for t in 0:n
         #print(t," ")
         if t % every == 0 
            writeMultiF(ff, U, t, h) 
            println(gg,t*h," ",r(gammaU)," ",s(gammaU))
            println(hh,signal())
            #println(t*h," ",r(1.)," ",s(1.))
            #println(t)
            #println("U-RK ",U)
         end
         for sd in 1:nSD
            #println(F[sd](U[sd],sd),coupling(sd))
            # aplica la función [sd] al ascilador [sd] + acoplamiento con todos los demas
            k1 = h*( F[sd](Ulocal[sd]) .+ coupling(Ulocal,sd) )
            #println("k1 $sd",k1)
            k1x = Ulocal[sd] .+ 0.5*k1
            k2 = h*( F[sd](k1x) .+ coupling(Ulocal,sd) )
            #println("k2 ",k2)
            k2x = Ulocal[sd] .+ 0.5*k2
            k3 = h*( F[sd](k2x) .+ coupling(Ulocal,sd) )
            #println("k3 ",k3)
            k3x = Ulocal[sd] .+ k3
            k4 = h*( F[sd](k3x) .+ coupling(Ulocal,sd) )
            #println("k4 ",k4)
            #println("el RK added ",(k1 .+ 2*k2 .+ 2*k3 .+ k4))/6
            # nuevos valores del oscilador real calculados con todos los estados de todos
            U[sd] = U[sd] .+ (k1 .+ 2*k2 .+ 2*k3 .+ k4)/6
         end # for sd
         # actualiza el conjunto "local"
         Ulocal = deepcopy(U)
      end # for t
      close(ff)
      close(gg)
      close(hh)
      ###
   end
   
      #===== ESCRITURA =====#
   function writeMultiF(ff, X, t, h)
      global nSD,dim
      #nSD = size(X)[1]     # numero de SD's acoplados
      #dimSD = size(X[1])[2]  # dimension (variables) de un SD
      @printf(ff,"%.10f ",t*h)
      for sd in 1:nSD
         for i in 1:dimSD 
            @printf(ff,"%.10f ",X[sd][i])
         end
      end
      @printf(ff, "\n")
   end

   function writeMultiM(ff, X, t)
      global nSD,dim
      
      @printf(ff,"%.10f ",t)
      for sd in 1:nSD
         for i in 1:dimSD 
            @printf(ff,"%.10f ",X[sd][i])
         end
      end
      @printf(ff, "\n")
   end

end #module
    
#U = Any[[.1 .2], [0.3 0.4], [2.5 2.6], [.7 .8]]
