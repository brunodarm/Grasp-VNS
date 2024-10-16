include("FunctionsH2V.jl")
include("H2V_MILP_Function.jl")

t=1

for f in 1:t
    cc, p, d, nmr, n, R, g = readInstance(11)  


    obj_val,x_val,y_val,h_val,w_val=SolverH2V_02(cc,p,d,R,n,g)
    





    println(" Valor Obejetivo z = ",obj_val)
    println("Decision variable x:")
    for i in 1:n
        for j in 1:n
            for r in 1:R
                if Int(round(x_val[i,j,r]))==1
                    println(i,"   ",j,"    ",r,"  x[",i,",",j,",",r,"] = ",Int(round(x_val[i,j,r])))
                end
            end    
        end
    end
    
    
    println()
    println("Decision variable y:")
    for i in 1:n
        for j in 1:n
            for r in 1:R
                if Int(round(y_val[i,j,r]))==1
                    println(i,"   ",j,"    ",r,"  y[",i,",",j,",",r,"] = ",Int(round(y_val[i,j,r])))
                end
            end
        end
    end
    
    println()
    println("Decision variable h:")
    for i in 1:n
        if Int(round(h_val[i]))==1
            println(i,"   h[",i,"] = ",Int(round(h_val[i])))
        end
    end
    
    println()
    println("Decision variable W:")
    for i in 1:n
        for j in 1:n
            
            if Int(round(w_val[i,j]))==1
                println(i,"   ",j,"    ","  w[",i,",",j,",","] = ",Int(round(w_val[i,j])))
            end
             
        end
    end
    
    println("The model was run!")
    
end
