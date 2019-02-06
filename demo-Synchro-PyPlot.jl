#include("Synchro.jl")
#include("Array4.jl")

using Synchro
using Array4
using Statistics
using Printf
using DelimitedFiles
#========================================
MODULO DE GRAFICOS

 Modelo para construccion de problemas de
sincronizacion de SD's

Usando DSys3.jl con RungeKuttaMultiD() en
la solucion numerica de los SD acoplados.
========================================#

include("demo-Synchro-param.jl")

#M = loadarray("$filename-$nettype-Madj-$nOsc.dat")
#println("M---",size(M)) 

#===== GRAFICAS =====#
using PyPlot 

# lectura de archivos en un array
println("==========================================================================================")
println("  leyendo archivos en z, r, S...")
println()
println("figure(k) muestra la ventana k")
println("saveim(k,\"opcional\",X,Y) guarda imagen k, donde")
println("     \"opcional\" es un agregado al nombre de la imagen")
println("     X,Y el tamaño de la imagen en pulgadas.")
println("==========================================================================================")

zf = [loadarray("$filename-$nettype-$nOsc-$s-$miumbral.dat") for s in rangeeps]    # oscilaciones
rf = [loadarray("$filename-$nettype-$nOsc-$s-$miumbral-tRS") for s in rangeeps]    # r(t), s(t)
Sf = [loadarray("$filename-$nettype-$nOsc-$s-$miumbral-SG.dat") for s in rangeeps] # la señal
dL = [loadarray("$filename-$nettype-distancesList-$nOsc-$s-$miumbral.dat") for s in rangeeps] # distancias list
dH = [loadarray("$filename-$nettype-distancesHist-$nOsc-$s-$miumbral.dat") for s in rangeeps] # distancias dist
aL = [loadarray("$filename-$nettype-amplitudesList-$nOsc-$s-$miumbral.dat") for s in rangeeps] # amplitudes list
aH = [loadarray("$filename-$nettype-amplitudesHist-$nOsc-$s-$miumbral.dat") for s in rangeeps] # amplitudes dist

#hf = loadarray("$filename-$nettype-Hist-$nOsc-0.0-$miumbral.dat")

println("  zf[i]->osciladores $(size(zf[1])); rf[i]->r(t),s(t) $(size(rf[1])); Sf[i]->señal $(size(Sf[1]))")
println("  dL[i]->distancias $(size(zf[1])); rf[i]->r(t),s(t) $(size(rf[1])); Sf[i]->señal $(size(Sf[1]))")
println("  para cada valor i de peso de acoplamiento.")
println("==========================================================================================")

# funcion para guardar GRAFICAS
function saveim(i,addname="", x::Int=6, y::Int=4)
   figure(i)
   savefig("$filename-$nettype-$nOsc-w$i-$miumbral-$addname.eps",
            format="eps",
            dpi=600,
            orientation="landscape",
            papertype="A6",
            transparent=true,
            frameon=false)
   savefig("$filename-$nettype-$nOsc-w$i-$miumbral-$addname.svg",
            format="svg",
            dpi=600,
            orientation="landscape",
            papertype="A6",
            transparent=true,
            frameon=false)
   savefig("$filename-$nettype-$nOsc-w$i-$miumbral-$addname.png",
            format="png",
            dpi=600,
            orientation="landscape",
            papertype="A6",
            transparent=true,
            frameon=false)
   savefig("$filename-$nettype-$nOsc-w$i-$miumbral-$addname.jpg",
            format="jpg",
            dpi=600,
            orientation="landscape",
            papertype="A6",
            transparent=true,
            frameon=false)
end;

function saveimG(i,addname="", x::Int=6, y::Int=4)
   figure(i)
   set_filename("$filename-$nettype-$nOsc-w$i-$miumbral-$addname.png")
   set_print_size("$(x)in,$(y)in")
   printfigure("png")
end

function saveall()
   for i in 1:nOsc
      saveim(i)
   end
end;

function closeall()
   PyPlot.close("all")
   #Gaston.closeall()
end

function showOsc(k::Int, mivar=1)
   figure(k)
   
end

function showRS()

end

function showSignal()

end

closeall()

function plotwithpyplot()
   # creacion de imagenes
   figcounter = 1    # contador de figura
   acounter = 1      # array counter
   thevar = 1        # por ahora
   for s in rangeeps
      # osciladores en la variable X
      fig1 = PyPlot.figure(figcounter,figsize=(12,5))
      #fig = figure()
      PyPlot.subplot(211)
      PyPlot.xlim(0.0,npts)
      p = PyPlot.plot(zf[acounter][:,thevar+1:nOsc+1])
      PyPlot.axis("tight")
      PyPlot.grid("on")
      PyPlot.title(string("OSCILACIONES $filename-$nettype nOsc=$nOsc ",L"\epsilon","=$s Umbral=$miumbral"))
      PyPlot.xlabel(L"tiempo")
      label_y = varnames[thevar]
      PyPlot.ylabel(string("$label_y",L"(t)"))
   
      # la señal
      PyPlot.subplot(212)
      PyPlot.xlim(0.0,npts)
      H = Hurst(Sf[acounter][:])
      PyPlot.plot(Sf[acounter][:],"r-",label="H suma = $H")
      xeeg,yeeg = FindLocalExtrema(Sf[acounter][:])
      H = Hurst(yeeg)
      PyPlot.plot(xeeg,yeeg,"k-",label="H señal = $H")
      PyPlot.legend()
      #PyPlot.title("SEÑAL                                                             ")
      #PyPlot.xlabel(L"tiempo")
      PyPlot.ylabel(L"señal \ EEG")

      figcounter += 1
   
      # parametros r(t) y s(t)
      fig2 = PyPlot.figure(figcounter,figsize=(10,12))
      PyPlot.subplot(311)
      PyPlot.xlim(0.0,npts)
      PyPlot.ylim(0.0,1.1)
      PyPlot.plot(rf[acounter][:,2],"g-",markersize=1.5,label=L"r(t)")
      PyPlot.plot(rf[acounter][:,3],"b-",markersize=1.5,label=L"s(t)")
      legend()
      #PyPlot.title("$filename-$nettype nOsc=$nOsc Eps=$s Umbral=$miumbral")
      PyPlot.xlabel(L"tiempo")
      PyPlot.ylabel(L"r(t),s(t)")
   
      # las distribuciones de distancias
      PyPlot.subplot(312)
      PyPlot.xlim(0.0,maximum(dH[2])+midist/2)
      PyPlot.ylim(0.0,1.0)
      PyPlot.bar(dH[acounter][:,1],dH[acounter][:,2],width=midist,alpha=0.4,color="red",
                 label=string("distribución de distancias y"));
      PyPlot.bar(aH[acounter][:,1],aH[acounter][:,2],width=midist,alpha=0.4,color="blue",
                 label=string(L"amplitudes $\epsilon$=","$s  ",L"$d_{umbral}$=","$midist"));
      PyPlot.plot([midist midist],[0.0 1.0],"ro-",linewidth=2.0)
      PyPlot.legend()
      #PyPlot.title("Distribución de distancias $filename nOsc=$nOsc Eps=$s Umbral=$miumbral")
      #PyPlot.xlabel(L"$distancia \ (d_{ij})$")
      #PyPlot.ylabel(L"$f(d_{ij}) // d_{ij}$")
      
      # las distribuciones de distancias
      PyPlot.subplot(313)
      #PyPlot.xlim(0.0,maximum(dH[2])+midist/2)
      #PyPlot.ylim(0.0,1.0)
      PyPlot.plot(dL[acounter][:],alpha=0.4,label=string("distancias y"));
      PyPlot.plot(aL[acounter][:],alpha=0.4,label=string(L"amplitudes $\epsilon$=","$s  ",L"$d_{umbral}$=","$midist"));
      PyPlot.plot([midist midist],[0.0 1.0],"ro-",linewidth=2.0)
      PyPlot.legend()
      PyPlot.title("Distancias y Amplitudes $filename nOsc=$nOsc Eps=$s Umbral=$miumbral")
      #PyPlot.xlabel(L"$distancia \ (d_{ij})$")
      #PyPlot.ylabel(L"$f(d_{ij}) // d_{ij}$")
   
      figcounter += 1
   
      suptitle = "Eps=$s"
   
      acounter += 1
   end  # for
end   # plotwithpyplot   
   

#plotwithpyplot()

sizez1 = size(zf[1])
sizer1 = size(rf[1])
sizeS1 = size(Sf[1])
println("==========================================================================================")
println("figure(k) muestra la ventana k")
println("saveim(k,\"opcional\",X,Y) guarda imagen k, donde")
println("     \"opcional\" es un agregado al nombre de la imagen")
println("     X,Y el tamaño de la imagen en pulgadas. Se guarda en png, eps, svg y eps.")
println("  zf[i]->osciladores $(sizez1); rf[i]->r(t),s(t) $(sizer1); Sf[i]->señal $(sizeS1)")
println("  para cada valor i de peso de acoplamiento.")
println("==========================================================================================")

