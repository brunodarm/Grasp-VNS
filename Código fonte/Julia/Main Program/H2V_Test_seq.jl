include("FunctionsH2V.jl")
include("H2V_MILP_Function.jl")

t=1


for f in 1:t
    cc, p, d, nmr, n, R, g = readInstance(3)  


    obj_val,x_val,y_val,h_val,w_val,modelo=SolverH2V(cc,p,d,R,n,g)
    

    ShowResults(obj_val,x_val,y_val,h_val,w_val,n,R)

    println(g)
end

H2V_Heuristic_Construtive(41)