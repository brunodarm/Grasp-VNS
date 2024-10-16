##########################################  Teste 01   ##############################################
    i=1

    cc, p, d, nmr, n, R, g = readInstance(i)
    n

    obj_val_Solver,x_val,y_val,h_val,w_val,modelo =SolverH2V_03(cc,p,R,n,g,3600)
    #sum(x_val)+sum(y_val)+sum(w_val)+sum(h_val)
    #obtj=[800000]
    #obj_val,x_val,y_val,h_val,w_val,modeloas= SolverH2V_02(cc,p,R,n,g,Inf)
    #modeloMestre,x_cc=Relaxado(A,gx,l,ux,cx)

    #subProblema,x_c=Relaxado(As,gxs,ls,uxs,cxs)


    #obj_val,x_val,y_val,h_val,w_val,gr = Solution_Construtive(i)
    #x_c=[]
    #x_c=hcat(reshape(x_val,(1,n*n*R)),reshape(y_val,(1,n*n*R)),reshape(w_val,(1,n*n)),h_val)
    #x_c

    #col=A*x_c'
    modelos=ProbRestrito(cc,p,R,n,g)
    A,gx,cx,l,ux =Empacote(modelos)
    s_modelos=Subproblem(cc,p,R,n,g)
    As,gxs,cxs,ls,uxs =Empacote(s_modelos)


    flag,m,objetivo_valor = geraPrimeiraColuna(A,l,ux)
    objetivo_valor
    obtj=[4070.00]
    #obtj=[objetivo_valor]
    col=m
    x, objt= geradorColunas(modelos,s_modelos,obtj,col)
    #obtj=[1000.0]
    col
    obtj
    objt
    objt[end]
    sum(x)




####################################################################################################################


    i=1

    cc, p, d, nmr, n, R, g = readInstance(i)
    n

    obj_val,x_val,y_val,h_val,w_val,modelo = Subproblem(cc,p,R,n,g,2)


    # if termination_status(assignment)==1
    #     obj_val = objective_value(assignment)
    # else
    #     obj_val = "No_Optimal"
    # end


    sum(x_val)
    sum(y_val)
    sum(w_val)
    sum(h_val)
    println()
    


    println("The model was run!")
    #include("Graphs.jl")
    Graph_Plot(g,x_val,y_val,h_val,n)
    

###########################Geração de Colunas de Thomadsen and Larsen ###########################################
    function MasterProblem(a_c,s_b,n,K,B)
        c_c=reshape(zeros(n*n),(n,n))
        c_b=reshape(zeros(n*n),(n,n))
        for i=1:n
            c_c[i,j]=sum(a_c[i]*a_c[j]*c[i,j])
            c_b[i,j]=sum(s_b[i]*s_b[j]*c[i,j])
        end

        modelo = Model(HiGHS.Optimizer)
        @variable(modelo, u_c[1:K], Bin)
        @variable(modelo, v_b[1:B], Bin)

        @objective(modelo, Min,sum(c_c'*u_c) +sum(c_b'*v_b))

        @constraint(modelo,[i=1:n, k=1:K], sum(a_c[i].*u_c[k])==1)
        @constraint(modelo,[i=1:n, k=1:K, b=1:B], -sum(s_c[i].*u_c[k]) + sum(s_b[i].*v_b[b])==0)

        return u_c,v_b
    end
    function ClusterGeneration(c,n,u,v,v_min,v_max)
        B=u
        alfa=v
        modelo = Model(HiGHS.Optimizer)
        #assignment = Model(optimizer_with_attributes(CPLEX.Optimizer,"CPX_PARAM_TILIM" => time_Limit, "CPX_PARAM_THREADS" => 1))
        
        
        set_silent(modelo)
        @variable(modelo, a[1:n], Bin)
        @variable(modelo, x[1:n, 1:n], Bin)
        @variable(modelo, S[1:n], Bin)
        
        @objective(modelo, Max,sum(alfa*a) -sum(B*S)+v- sum(c*y))

        @constraint(modelo,sum(S)==1)
        @constraint(modelo,[i=1:n, j=1:n; i!=j], S[i]<=a[i])
        @constraint(modelo,[i=1:n, j=1:n; i!=j], x[i,j]<=a[i])
        @constraint(modelo,[i=1:n, j=1:n; i!=j], x[i,j]<=a[j])
        @constraint(modelo,[i=1:n, j=1:n; i!=j], a[i]+a[j]<=x[i,j]+1)
        @constraint(modelo, v_min<=sum(a)<=v_max)

        return x,S

    end
    function BackboneGeneration(c,n,u,v,b_min,b_max)
        B=u
        modelo = Model(HiGHS.Optimizer)
        #assignment = Model(optimizer_with_attributes(CPLEX.Optimizer,"CPX_PARAM_TILIM" => time_Limit, "CPX_PARAM_THREADS" => 1))
        set_time_limit_sec(modelo,time_Limit)
        
        set_silent(modelo)
        @variable(modelo, y[1:n, 1:n], Bin)
        @variable(modelo, S[1:n], Bin)
        
        @objective(modelo, Max, sum(B*S)+v- sum(c*y))

        @constraint(modelo,[i=1:n, j=1:n; i!=j], y[i,j]<=S[i])
        @constraint(modelo,[i=1:n, j=1:n; i!=j], y[i,j]<=S[j])
        @constraint(modelo,[i=1:n, j=1:n; i!=j], S[i]+S[j]<=y[i,j]+1)
        @constraint(modelo, b_min<=sum(S)<=b_max)

        return y,S
    end
    function QKP_model(P)
        modelo = Model(HiGHS.Optimizer)
        #assignment = Model(optimizer_with_attributes(CPLEX.Optimizer,"CPX_PARAM_TILIM" => time_Limit, "CPX_PARAM_THREADS" => 1))
        set_silent(modelo)
        @variable(modelo, q[1:n, 1:n], Bin)
        @variable(modelo, qi[1:n], Bin)
        @variable(modelo, w[1:n])
        @objective(modelo, Max,sum( [i=1:n, j=1:n; i<j], P[i,i]*q[i,i])+sum([i=1:n, j=1:n; i<j],P[i,j]*q[i,j]))

    end
    function MasterProblem_QKP()

    end

