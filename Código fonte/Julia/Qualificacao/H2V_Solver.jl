

import DataFrames, XLSX
using DataFrames, XLSX
#include("FunctionsH2V.jl")
include("FunctionsH2V_02.jl")
include("H2V_MILP_Function.jl")
include("H2V_Heuristica.jl")
include("Arvore_Minima.jl")
include("GRASP.jl")
nome_planilha="H2V_Results_V0002.xlsx"
nnx=50
try
    XLSX.writetable(nome_planilha)
catch
    println("Já existe planilha com mesmo nome")
end


XLSX.openxlsx(nome_planilha, mode="w") do xf
    sheet = xf[1]
    sheet[1,1]="numero"
    sheet[1,2]="Cidades Portuarias(parametro)"
    sheet[1,3]="Cenectividade Terrestre(parametro)"
    sheet[1,4]="Instancia"
    sheet[1,5]="(n) Size of the sample"
    sheet[1,6]="(R) Modais de Transporte"
    sheet[1,7]=" Tempo (s))"
    sheet[1,8]="GAP"
    sheet[1,9]="Solver Objetive"#"Tempo Total"
    sheet[1,10]="CG lowBound Objetive"
    l=1
    f = 180
    
    for i in l:f
        cc, p, d, nmr, n, R, g = readInstance(i)
        #saida = read_data(i) 
        #i, n, R, MIP_value,MIP_time,MIP_gap,GV_value,GV_time,GV_gap,ML_value,ML_time,ML_gap
        if n<=nnx
            
            #obj_val_Solver,x_val,y_val,h_val,w_val,modelo =SolverH2V(cc,p,d,R,n,g,3600)
            obj_val_Solver,x_val,y_val,h_val,w_val,modelo =SolverH2V_03(cc,p,R,n,g,3600)
            #obj_val_Heuristic,x_val,y_val,h_val,w_val =H2V_Heuristic_Construtive(i)
            aMin,aMax=0.55,0.85
            #obj_val_Heuristic,gg,hh = H2V_Grasp(cc,g,p,aMin,aMax)
            mmx,interacoes = GeraColumLowBound(i,100)
            sheet[1+i,1]=i
            sheet[1+i,2]=string("p",i)
            sheet[1+i,3]=string("g",i)
            sheet[1+i,4]=string("H2V_Results_",i)
            sheet[1+i,5]=n
            sheet[1+i,6]=R
            sheet[1+i,7]=string(round(100*solve_time(modelo))/100)
            sheet[1+i,8]=interacoes#relative_gap(modelo)
            sheet[1+i,9]=obj_val_Solver#"Fail"#string(round(100*solve_time(modelo))/100," s")
            sheet[1+i,10]=mmx #obj_val_Heuristic#"Fail"#relative_gap(modelo)
        
            #df = DataFrame(XLSX.readtable("H2V_Results_V000.xlsx","Sheet1"))
        end
    end
end
df = DataFrame(XLSX.readtable(nome_planilha,"Sheet1"))