

import DataFrames, XLSX
using DataFrames, XLSX
include("FunctionsH2V.jl")
include("H2V_MILP_Function.jl")
include("H2V_Heuristica.jl")
try
    XLSX.writetable("H2V_Results.xlsx")
catch
    println("Já existe planilha com mesmo nome")
end


XLSX.openxlsx("H2V_Results.xlsx", mode="w") do xf
    sheet = xf[1]
    sheet[1,1]="numero"
    sheet[1,2]="Cidades Portuarias(parametro)"
    sheet[1,3]="Cidades Produtoras(parametro)"
    sheet[1,4]="Cenectividade Terrestre(parametro)"
    sheet[1,5]="Instancia"
    sheet[1,6]="(n) Size of the sample"
    sheet[1,7]="(R) Modais de Transporte"
    sheet[1,8]=" Tempo Total"
    sheet[1,9]="GAP"
    sheet[1,10]="Solver Objetive"#"Tempo Total"
    sheet[1,11]="Heuristic Objetive"
    l=1
    f = 180
    for i in l:f
        cc, p, d, nmr, n, R, g = readInstance(i) 
        if n<=f#100
            
            obj_val_Solver,x_val,y_val,h_val,w_val,modelo =SolverH2V(cc,p,d,R,n,g,300)
            obj_val_Heuristic,x_val,y_val,h_val,w_val =H2V_Heuristic_Construtive(i)
            sheet[1+i,1]=i
            sheet[1+i,2]=string("p",i)
            sheet[1+i,3]=string("d",i)
            sheet[1+i,4]=string("g",i)
            sheet[1+i,5]=string("H2V_n",n,"_r",R,"_nmr",i)
            sheet[1+i,6]=n
            sheet[1+i,7]=R
            sheet[1+i,8]=string(round(100*solve_time(modelo))/100," s")
            sheet[1+i,9]=relative_gap(modelo)
            sheet[1+i,10]=obj_val_Solver#"Fail"#string(round(100*solve_time(modelo))/100," s")
            sheet[1+i,11]=obj_val_Heuristic#"Fail"#relative_gap(modelo)
            #df = DataFrame(XLSX.readtable("H2V_Results_Solver.xlsx","Sheet1"))
        end
    end
end
df = DataFrame(XLSX.readtable("H2V_Results.xlsx","Sheet1"))