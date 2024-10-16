include("GRASP.jl")
import DataFrames
using DataFrames


n=9

c=geraCusto(n)
#println("C=",c)
inicio = 1
ganacioso =greedy(c,inicio) 
#println(ganacioso)

# k = interações no GRASP
k=1000

#alpha a =0.85
alfa=0.85
graspPure= GRASP(c,inicio,alfa,k)
#println(graspPure)

#Alfa mini aMin e Alfa Maximo aMax
aMin=0.70
aMax=0.97
graspEvo= evoGRASP(c,inicio,aMin,aMax,k)
#println(graspEvo)


algoritmo = ["Ganacioso", "GRASP α fixo", "Grasp α variado"]
df_a = DataFrame(algoritmo=algoritmo,α=[string("|   0   "),string("|  ",alfa,"  "),string("|",aMin," - ",aMax)],Complexidade=["|  n²","|  kn²","|  kn³"])
println(df_a)
println("")
df=DataFrame(algoritmo=algoritmo, Custo=[ganacioso[2],graspPure[2],graspEvo[2]],Caminho_Fechado=[ganacioso[1],graspPure[1],graspEvo[1]])
print(df)

