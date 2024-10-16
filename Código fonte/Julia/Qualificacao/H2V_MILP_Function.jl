using JuMP#, CPLEX
import HiGHS

function SolverH2V(cc,p,d,R,n,g,time_Limit=Inf)
    
    assignment = Model(HiGHS.Optimizer)
    #assignment = Model(optimizer_with_attributes(CPLEX.Optimizer,"CPX_PARAM_TILIM" => time_Limit, "CPX_PARAM_THREADS" => 1))
    set_time_limit_sec(assignment,time_Limit)
    
    set_silent(assignment)
    @variable(assignment, x[1:n, 1:n, 1:R], Bin)
    @variable(assignment, y[1:n, 1:n, 1:R], Bin)
    @variable(assignment, h[1:n], Bin)
    @variable(assignment, w[1:n, 1:n, 1:R], Bin)
    #@variable(assignment, s[1:2, 1:n], Bin)
    #--------------------------------------
    @constraint(assignment,[i = 1:n,j = 1:n,r = 1:R], (x[i,j,r])==x[j,i,r])
    @constraint(assignment,[i = 1:n,j = 1:n,r = 1:R], (y[i,j,r])==y[j,i,r])
    @constraint(assignment,[i = 1:n, j = 1:n; i!=j],sum(x[i,j,:])<= g[i,j])
    @constraint(assignment,[i = 1:n, j = 1:n; i!=j],sum(y[i,j,:])+sum(x[i,j,:])<=1)
    @constraint(assignment, [i = 1:n, j = 1:n; i!=j], h[i] + h[j] + sum(x[i,j,:]) <= 2)
    @constraint(assignment,[i = 1:n, j = 1:n, r = 1:R; i!=j], 2*y[i,j,r]<= h[i]+h[j])
    @constraint(assignment,[i = 1:n, j = 1:n, k = 1:n; k!=i && k!=j && i<j ], sum(x[k,i,:])+sum(x[k,j,:])<=sum(x[i,j,:])+1)#+s[1,k])
    #@constraint(assignment,[i = 1:n, j = 1:n, r = 1:R, k = 1:n; k!=i && k!=j && i<j ], (x[k,i,r])+(x[k,j,r])<=(x[i,j,r])+1)#+s[2,k])
    #@constraint(assignment,[i = 1:n], s[1,i]+s[2,i]==1)

    #--Hub

    @constraint(assignment, [i = 1:n, j = 1:n; i!=j], p[i]*h[i] + p[j]*h[j] <=y[i,j,1] + 1)
    @constraint(assignment, [i = 1:n, j = 1:n; i!=j], h[i] + h[j] + g[i,j] <=sum(y[i,j,:]) + 2)
    #@constraint(assignment, [i = 1:n, j = 1:n, r = 2:R; i!=j], y[i,j,r] <=g[i,j])
    #@constraint(assignment, [i = 1:n, j = 1:n, r = 2:R; i!=j], x[i,j,r] <=g[i,j])
    @constraint(assignment,[i = 1:n, j = 1:n, r = 1:R; i!=j], w[i,j,r]<= x[i,j,r])
    @constraint(assignment,[i = 1:n, j = 1:n, r = 1:R; i!=j], w[i,j,r]<= h[i])
    @constraint(assignment,[i = 1:n, j = 1:n, r = 1:R; i!=j], h[i] + x[i,j,r] <=w[i,j,r]+ 1)
    @constraint(assignment, [j = 1:n], h[j] + sum(w[:, j,:]) == 1 )
    @constraint(assignment, [ j = 1:n, r = 1:R], x[j, j,r]+y[j, j,r]+w[j, j,r]  == 0)
    @constraint(assignment, [ i = 1:n, j = 1:n; i!=j], x[i,j,1]<=p[i]*p[j])
    @constraint(assignment, [ i = 1:n, j = 1:n; i!=j], y[i,j,1]<=p[i]*p[j])
    @constraint(assignment, sum(w[:,:,:])== (n-sum(h[:])))
    #@constraint(assignment, sum(x[:,:,:])>= (n+(R-1)*sum(h[:])))
    #@constraint(assignment,[ j = 1:n], sum(x[j,:,:])>= 1)
    
    @constraint(assignment, sum(y[:,:,:])>= 2*(sum(h)-1))
    #@constraint(assignment, sum(h[:])>= sum(d[:]))
    @constraint(assignment,[ i = 1:n],d[i]<=h[i])

    #@constraint(assignment, sum(h[:]) - sum(d[:]) <= (30/100)*(n - sum(d[:]) ))
    # ---- OBJEcTIVE 
    @objective(assignment, Min, sum(cc .* y)+sum(cc .* x))

    optimize!(assignment) 
   # if termination_status(assignment)==1
   #     obj_val = objective_value(assignment)
   # else
   #     obj_val = "No_Optimal"
   # end
    obj_val = objective_value(assignment)
    x_val = value.(x)
    y_val = value.(y)
    h_val = value.(h)
    w_val = value.(w)
    println()
    modelo=assignment
    println("The model was run!")

    return(obj_val,x_val,y_val,h_val,w_val,modelo)
end

function SolverH2V_02(cc,p,R,n,g,UperBound,time_Limit=Inf)
    include("FunctionsH2V.jl")
    m = Cluster_Ground(g)
    assignment = Model(HiGHS.Optimizer)
    #assignment = Model(optimizer_with_attributes(CPLEX.Optimizer,"CPX_PARAM_TILIM" => time_Limit, "CPX_PARAM_THREADS" => 1))
    set_time_limit_sec(assignment,time_Limit)
    
    set_silent(assignment)
    @variable(assignment, x[1:n, 1:n, 1:R], Bin)
    @variable(assignment, y[1:n, 1:n, 1:R], Bin)
    @variable(assignment, h[1:n], Bin)
    @variable(assignment, H[1:n, 1:n], Bin)
    #@variable(assignment, s[1:2, 1:n], Bin)
   

    #-----------------------------------------------
    @constraint(assignment,[i = 1:n, j = 1:n; i!=j],sum(x[i,j,:])<= g[i,j])
    @constraint(assignment,[i = 1:n, j = 1:n; i!=j],sum(y[i,j,2:R])<= g[i,j])
    @constraint(assignment,[i = 1:n, j = 1:n; i!=j],sum(y[i,j,:])+sum(x[i,j,:])<=1)
    @constraint(assignment, [i = 1:n, j = 1:n; i!=j], h[i] + h[j] + sum(x[i,j,:]) <= 2)
    @constraint(assignment,[i = 1:n, j = 1:n, r = 1:R; i!=j], 2*y[i,j,r]<= h[i]+h[j])

    @constraint(assignment,[i = 1:n, j = 1:n, k = 1:n; k!=i && k!=j && i!=j],(H[i,j]+H[i,k])<=sum(x[j,k,:])+1)#+s[1,k])

    #@constraint(assignment,[i = 1:n, j = 1:n, k = 1:n; k!=i && k!=j && i<j ], sum(x[k,i,:])+sum(x[k,j,:])<=sum(x[i,j,:])+1)#+s[1,k])
    #@constraint(assignment,[i = 1:n, j = 1:n, r = 1:R, k = 1:n; k!=i && k!=j && i<j ], (x[k,i,r])+(x[k,j,r])<=(x[i,j,r])+1)#+s[2,k])
    #@constraint(assignment,[i = 1:n], s[1,i]+s[2,i]==1)

    #--Hub

    
    @constraint(assignment, [i = 1:n, j = 1:n; i!=j], h[i] + h[j] + g[i,j] <=sum(y[i,j,1:R]) + 2)
    #@constraint(assignment, [i = 1:n, j = 1:n; i!=j], h[i] + h[j] <=sum(y[i,j,:]) + 1)
    #@constraint(assignment, [i = 1:n, j = 1:n; i!=j], h[i] + h[j] + g[i,j] >=sum(y[i,j,2:R]) + 2)

    

    @constraint(assignment,[i = 1:n, j = 1:n; i!=j], H[i,j]<= sum(x[i,j,:]))
    
    @constraint(assignment,[i = 1:n, j = 1:n; i!=j], H[i,j]<= (h[i]))
    @constraint(assignment,[i = 1:n, j = 1:n; i!=j], H[i,j]<= (1-h[j]))
    @constraint(assignment,[j = 1:n], sum(H[:,j])<=1)
    

    @constraint(assignment,[i = 1:n,j = 1:n,r = 1:R], (x[i,j,r])==x[j,i,r])
    @constraint(assignment,[i = 1:n,j = 1:n,r = 1:R], (y[i,j,r])==y[j,i,r])
   
   
    @constraint(assignment, [ j = 1:n, r = 1:R], x[j, j,r]+y[j, j,r]+H[j, j]  == 0)
    @constraint(assignment, [ i = 1:n, j = 1:n; i!=j], x[i,j,1]<=p[i]*p[j])
    @constraint(assignment, [ i = 1:n, j = 1:n; i!=j], y[i,j,1]<=p[i]*p[j])
    @constraint(assignment, sum(H[:,:])== (n-sum(h[:])))
   

   
    #@constraint(assignment, sum(h)>= 5)
    @constraint(assignment, sum(y[:,:,1])>= m*(m-1))
    
    @constraint(assignment, [i = 1:n], sum((transpose(p).*h).*g[i,:]) >=1)#sum(p))

    @constraint(assignment, [i = 1:n, j = 1:n; i!=j], p[i]*h[i] + p[j]*h[j] <=y[i,j,1] + 1)
    #--------------------------------------
    #@constraint(assignment,[i = 1:n, j = 1:n, r=1:R],0<=(x[i,j,r]) <=1)
   # @constraint(assignment,[i = 1:n, j = 1:n, r=1:R],0<=(y[i,j,r]) <=1)
    #@constraint(assignment,[i = 1:n, j = 1:n, r=1:R],0<=(H[i,j]) <=1)
   # @constraint(assignment,[i = 1:n],(h[i]) <=1)
   
   
   
    # ---- OBJEcTIVE 
    #@objective(assignment, Min, sum(cc .* y) +sum(cc .* x))
    @objective(assignment, Min, sum(cc .* (y+ x)))
    optimize!(assignment) 
   # if termination_status(assignment)==1
   #     obj_val = objective_value(assignment)
   # else
   #     obj_val = "No_Optimal"
   # end
    obj_val = objective_value(assignment)
    x_val = value.(x)
    y_val = value.(y)
    h_val = value.(h)
    w_val = value.(H)
    println()
    modelo=assignment
    println("The model was run!")
    #c,r,cz =H2V_graph(cc,p,g,n)
    #obj_val =custo_Grafo(c,x_val,y_val,R,n)

    return(obj_val,x_val,y_val,h_val,w_val,modelo)
end

function SolverH2V_03(cc,p,R,n,g,UperBound,time_Limit=Inf)
    include("FunctionsH2V.jl")
    m = Cluster_Ground(g)
    assignment = Model(HiGHS.Optimizer)
    #assignment = Model(optimizer_with_attributes(CPLEX.Optimizer,"CPX_PARAM_TILIM" => time_Limit, "CPX_PARAM_THREADS" => 1))
    set_time_limit_sec(assignment,time_Limit)
    
    set_silent(assignment)
    @variable(assignment, x[1:n, 1:n, 1:R] >=0)
    @variable(assignment, y[1:n, 1:n, 1:R] >=0)
    @variable(assignment, h[1:n] >=0)
    @variable(assignment, H[1:n, 1:n] >=0)
    #@variable(assignment, s[1:2, 1:n], Bin)
   

    #-----------------------------------------------
    @constraint(assignment,[i = 1:n, j = 1:n; i!=j],sum(x[i,j,:])<= g[i,j])
    @constraint(assignment,[i = 1:n, j = 1:n; i!=j],sum(y[i,j,2:R])<= g[i,j])
    @constraint(assignment,[i = 1:n, j = 1:n; i!=j],sum(y[i,j,:])+sum(x[i,j,:])<=1)
    @constraint(assignment, [i = 1:n, j = 1:n; i!=j], h[i] + h[j] + sum(x[i,j,:]) <= 2)
    @constraint(assignment,[i = 1:n, j = 1:n, r = 1:R; i!=j], 2*y[i,j,r]<= h[i]+h[j])

    @constraint(assignment,[i = 1:n, j = 1:n, k = 1:n; k!=i && k!=j && i!=j],(H[i,j]+H[i,k])<=sum(x[j,k,:])+1)#+s[1,k])

    #@constraint(assignment,[i = 1:n, j = 1:n, k = 1:n; k!=i && k!=j && i<j ], sum(x[k,i,:])+sum(x[k,j,:])<=sum(x[i,j,:])+1)#+s[1,k])
    #@constraint(assignment,[i = 1:n, j = 1:n, r = 1:R, k = 1:n; k!=i && k!=j && i<j ], (x[k,i,r])+(x[k,j,r])<=(x[i,j,r])+1)#+s[2,k])
    #@constraint(assignment,[i = 1:n], s[1,i]+s[2,i]==1)

    #--Hub

    
    @constraint(assignment, [i = 1:n, j = 1:n; i!=j], h[i] + h[j] + g[i,j] <=sum(y[i,j,1:R]) + 2)
    #@constraint(assignment, [i = 1:n, j = 1:n; i!=j], h[i] + h[j] <=sum(y[i,j,:]) + 1)
    #@constraint(assignment, [i = 1:n, j = 1:n; i!=j], h[i] + h[j] + g[i,j] >=sum(y[i,j,2:R]) + 2)

    

    @constraint(assignment,[i = 1:n, j = 1:n; i!=j], H[i,j]<= sum(x[i,j,:]))
    
    @constraint(assignment,[i = 1:n, j = 1:n; i!=j], H[i,j]<= (h[i]))
    @constraint(assignment,[i = 1:n, j = 1:n; i!=j], H[i,j]<= (1-h[j]))
    @constraint(assignment,[j = 1:n], sum(H[:,j])<=1)
    

    @constraint(assignment,[i = 1:n,j = 1:n,r = 1:R], (x[i,j,r])==x[j,i,r])
    @constraint(assignment,[i = 1:n,j = 1:n,r = 1:R], (y[i,j,r])==y[j,i,r])
   
   
    @constraint(assignment, [ j = 1:n, r = 1:R], x[j, j,r]+y[j, j,r]+H[j, j]  == 0)
    @constraint(assignment, [ i = 1:n, j = 1:n; i!=j], x[i,j,1]<=p[i]*p[j])
    @constraint(assignment, [ i = 1:n, j = 1:n; i!=j], y[i,j,1]<=p[i]*p[j])
    @constraint(assignment, sum(H[:,:])== (n-sum(h[:])))
   

   
    #@constraint(assignment, sum(h)>= 5)
    @constraint(assignment, sum(y[:,:,1])>= m*(m-1))
    
    @constraint(assignment, [i = 1:n], sum((transpose(p).*h).*g[i,:]) >=1)#sum(p))

    @constraint(assignment, [i = 1:n, j = 1:n; i!=j], p[i]*h[i] + p[j]*h[j] <=y[i,j,1] + 1)
    #--------------------------------------
    #@constraint(assignment,[i = 1:n, j = 1:n, r=1:R],0<=(x[i,j,r]) <=1)
   # @constraint(assignment,[i = 1:n, j = 1:n, r=1:R],0<=(y[i,j,r]) <=1)
    #@constraint(assignment,[i = 1:n, j = 1:n, r=1:R],0<=(H[i,j]) <=1)
   # @constraint(assignment,[i = 1:n],(h[i]) <=1)
   
   
   
    # ---- OBJEcTIVE 
    #@objective(assignment, Min, sum(cc .* y) +sum(cc .* x))
    @objective(assignment, Min, sum(cc .* (y+ x)))
    optimize!(assignment) 
   # if termination_status(assignment)==1
   #     obj_val = objective_value(assignment)
   # else
   #     obj_val = "No_Optimal"
   # end
    obj_val = objective_value(assignment)
    x_val = value.(x)
    y_val = value.(y)
    h_val = value.(h)
    w_val = value.(H)
    println()
    modelo=assignment
    println("The model was run!")
    #c,r,cz =H2V_graph(cc,p,g,n)
    #obj_val =custo_Grafo(c,x_val,y_val,R,n)

    return(obj_val,x_val,y_val,h_val,w_val,modelo)
end