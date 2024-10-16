include("FunctionsH2V_02.jl")
include("H2V_MILP_Function.jl")
include("H2V_Heuristica.jl")
include("Arvore_Minima.jl")
include("GRASP.jl")
using JuMP, NLPModelsJuMP, NLPModels
function Empacote(modelo)
    modelo=modelo
    nlp = MathOptNLPModel(modelo) 
    x = zeros(nlp.meta.nvar)
    cx=grad(nlp, x) 
    A=jac(nlp, x) 
    gx = cons(nlp, x)

    nlp.meta.lcon, nlp.meta.ucon # l <= Ax + g <= u
    l=nlp.meta.lcon
    u=nlp.meta.ucon
    for i in 1:length(gx)
        if u[i]==Inf
            u[i]=100000
        end
        if l[i]==-Inf
            l[i]=-100000
        end
    end
    

    return A,gx,cx,l,u
end

function Relaxado(A,gx,l,ux,cx)
    nx=Int(length(A)/length(gx))
    mdl = Model(HiGHS.Optimizer)
    @variable(mdl,x_c[1:nx] .>=0)
    Ax=A*x_c
    eq=Ax+gx
    @constraint(mdl,Lx[i=1:size(A)[1]],l[i]<=eq[i]<=ux[i])
    #@objective(mdl,Min, cx'*x_c)
    #optimize!(mdl)
    #value.(x_c)
    #dual.(Lx)
    
    return mdl  
end

function Relaxado_02(A,gx,l,ux,cx)
    nx=Int(length(A)/length(gx))
    mdl = Model(HiGHS.Optimizer)
    @variable(mdl,x_c[1:nx]>=0)
    Ax=A*x_c
    eq=Ax+gx
    @constraint(mdl,Lx[i=1:size(A)[1]],l[i]<=eq[i]<=ux[i])
    @objective(mdl,Min, (cx'*x_c))
    #optimize!(mdl)
    #value.(x_c)
    #dual.(Lx)
    
    return mdl        
end

function ProblemaMestreRestrito(col,l,ux,obtj,conta)
    colu=col
    conta=conta



    obtj=obtj
    tm=size(colu)[1]
    problemaMestre= Model(HiGHS.Optimizer)
    @variable(problemaMestre,m[1:conta] .>=0)
    
    @constraint(problemaMestre,L[i=1:tm],l[i]<=colu[i,:]'*m<=ux[i])
    @constraint(problemaMestre,v,sum(m)==1)
    @constraint(problemaMestre,m .>=0)
    #println("obtj=",obtj)
    #println("m=",m)
    @objective(problemaMestre,Min,sum(obtj.*m))
    set_silent(problemaMestre)
    optimize!(problemaMestre)
    
    solution_summary(problemaMestre)
    #println("m=",value.(m))
    u=dual.(L)
    #println("u= ",u)
   # println("v= ",dual(v))
    #println("objetivo Problema mestre = ",sum(obtj.*m))
    
    
    return u,dual(v),value.(m), sum(obtj.*value.(m))
    
end


function SubProb(c,u,v,A,subProblem)
    v=v
    u=u
    Aq,gxq,cxq,lq,uxq =Empacote(subProblem)
    Ab=Aq
    gb=gxq
    nx=Int(length(Ab)/length(gb))
    mdl = Model(HiGHS.Optimizer)
    set_silent(mdl)
    @variable(mdl,x_c[1:nx] .>=0)
    Ax=Ab*x_c
    eq=Ax+gb
    @constraint(mdl,Lx[i=1:size(Ab)[1]],lq[i]<=eq[i]<=uxq[i])
    
    #println("Lx= ",x_c)
    #set_silent(modelo_s)
    #modelo_s= Model(HiGHS.Optimizer)
    #@variable(modelo_s,x[1:nq].>=0)
    #println(" c =",c," tamanho c = ", size(c)[1])
    #println(" x =",x," tamanho x = ", size(x)[1])
    #println("size of x_c = ",size(x_c))
    #println("size of c =",size(c))
    #println(" c'*y= ",(c'*x_c))
    #x_c=reshape(x_c,(size(x_c)[1],1))
    part_1=c'*x_c
    #println(" -u'*(A*(x)) = ",(-u'*(A*(x_c))))
    part_2=(-u'*(A*(x_c)))
    @objective(mdl,Min, part_1 +part_2)
    #@objective(subProblema,Min, c'*x -sum(u'*(Aq*(x))))
    optimize!(mdl)
    value.(x_c)
    objective_value(mdl)

    ob=c'*value.(x_c)
    ob=Int(ob-mod(ob,1))
    #println(" Cy=",ob)
    #println("y= ",value.(y))
    #println(" Cy -v =", -c'*value.(y) +v*1.0)
    #println(" teste",)
    #colu[:,conta]= Aq*value.(y)
    #colu
    f=objective_value(mdl)
    #println("ob=",ob)
    #println("A*value.(x') = ",(Aq*value.(x_c)))
    #println("value.(x')=",(value.(x'))[1,:])
    #println("f=",f)
    return ob, A*value.(x_c),value.(x_c),f
end

function geradorColunas(modeloMestre,subProblema,ob_inicial,colu_inicial,tempo)
   
    obtj=ob_inicial
    colu=colu_inicial

    subProblem=subProblema
    permi="."
    yv=[]
    mx=[]
    conta=1
    vx=[]
    global rp=[]
    custoRed=0.0
    u=[]
    #somas=[0,0]
    
    global resu=0.0
    A,gx,cx,l,ux =Empacote(modeloMestre)
    Aq,gxq,cxq,lq,uxq =Empacote(subProblem)
    nq=size(colu)[1]
    modelo_s=Relaxado(Aq,gxq,lq,uxq,cxq)
    myy=[]
    yvv=[]
    last=-10000000000000000
    ob=10000
    while (permi==".")
        
         
        u,v,my,resu = ProblemaMestreRestrito(colu,l,ux,obtj,conta)
        push!(myy,my)
        push!(mx,[my])
        push!(rp,resu)
        push!(vx,v)
        println("LS =",rp[end])

        LS=vx[end]
        if conta>1
            LI=rp[end-1] + custoRed
            println("LI =",LI)
            v=ob
        else
        
            println("LI =",-Inf)
            LI=-Inf
        end        
        println("conta = ",conta)
        ob, Aa,yi,f =SubProb(cx,u,v,A,subProblem)
        #println("Aa=",Aa)
        conta=conta+1
        println("ob=",ob)
        push!(obtj,ob)
        LS=ob
        println("###########################################")
        println(" f =",f," and v = ",v)
    
        tempos=0
        
        println(" CUSTO RELATIVO ",conta-1," = ",f-v)
        if (f-v)<0 && conta<100 && tempos<=tempo
            permi="."
            colu=hcat(colu,Aa)
            push!(yv,yi)
            #println("yv =",yv)
        else
            permi="para"
        end
        #last=f-v
        #push!(custoRed,f-v)
        custoRed=f-v
    end
    #println("yv=",yv) 
    println("######################################################################################")
    push!(yvv,yv)
    #println(" tamanho de ",size(yv[end-i+1]))
    #print("mx[end][end][end-i+1]=",mx[end][end][end-i+1]) 
    #println(" tamanho de ",mx[end][end][end-i+1])
    #somas=[]
    somas=sum(yv[end-i+1].*mx[end][end][end-i+1] for i=1:length(yv))
    println("x=",somas)
    #println("objetivo =",rp[end])
    #somas=sum(yvv[end-i+1]*mx[end][end][end-i+1] for i=1:length(yv))
    #somas=sum(yvv[i]*mx[end][end][i] for i=1:length(yv))
    #println("yv=",yv*mx[end][end][1]) 
    #println("my =",mx[end][end])
    println("somas =",somas)
    return somas,rp,conta
    

end

################################ Problema Restrito ###########################################
function ProbRestrito(cc,p,R,n,g)
    include("FunctionsH2V.jl")
    m = Cluster_Ground(g)
    assignment = Model(HiGHS.Optimizer)
    #assignment = Model(optimizer_with_attributes(CPLEX.Optimizer,"CPX_PARAM_TILIM" => time_Limit, "CPX_PARAM_THREADS" => 1))
    #set_time_limit_sec(assignment,time_Limit)
    
    set_silent(assignment)
    @variable(assignment, x[1:n, 1:n, 1:R] >=0)
    @variable(assignment, y[1:n, 1:n, 1:R] >=0)
    @variable(assignment, H[1:n, 1:n] >=0)
    @variable(assignment, h[1:n] >=0)

    
   
   #-----------------------------------------------
   @constraint(assignment,[i = 1:n, j = 1:n; i!=j],sum(x[i,j,:])<= g[i,j])
   @constraint(assignment,[i = 1:n, j = 1:n; i!=j],sum(y[i,j,2:R])<= g[i,j])
   @constraint(assignment,[i = 1:n, j = 1:n; i!=j],sum(y[i,j,:])+sum(x[i,j,:])<=1)
   @constraint(assignment, [i = 1:n, j = 1:n; i!=j], h[i] + h[j] + sum(x[i,j,:]) <= 2)
   @constraint(assignment,[i = 1:n, j = 1:n, r = 1:R; i!=j], 2*y[i,j,r]<= h[i]+h[j])

   @constraint(assignment,[i = 1:n, j = 1:n, k = 1:n; k!=i && k!=j && i!=j],(H[i,j]+H[i,k])<=sum(x[j,k,:])+1)#+s[1,k])

  

   #--Hub

   
   @constraint(assignment, [i = 1:n, j = 1:n; i!=j], h[i] + h[j] + g[i,j] <=sum(y[i,j,1:R]) + 2)
   

   

   #@constraint(assignment,[i = 1:n, j = 1:n; i!=j], H[i,j]<= sum(x[i,j,:]))
   
   @constraint(assignment,[i = 1:n, j = 1:n; i!=j], H[i,j]<= (h[i]))
   @constraint(assignment,[i = 1:n, j = 1:n; i!=j], H[i,j]<= (1-h[j]))
   @constraint(assignment,[j = 1:n], sum(H[:,j])<=1)
   

   #@constraint(assignment,[i = 1:n,j = 1:n,r = 1:R], (x[i,j,r])==x[j,i,r])
  # @constraint(assignment,[i = 1:n,j = 1:n,r = 1:R], (y[i,j,r])==y[j,i,r])
  
  
   #@constraint(assignment, [ j = 1:n, r = 1:R], x[j, j,r]+y[j, j,r]+H[j, j]  == 0)
   @constraint(assignment, [ i = 1:n, j = 1:n; i!=j], x[i,j,1]<=p[i]*p[j])
   @constraint(assignment, [ i = 1:n, j = 1:n; i!=j], y[i,j,1]<=p[i]*p[j])
   #@constraint(assignment, sum(H[:,:])== (n-sum(h[:])))
  

  
  
   @constraint(assignment, sum(y[:,:,1])>= m*(m-1))
   
   #@constraint(assignment, [i = 1:n], sum((transpose(p).*h).*g[i,:]) >=1)#sum(p))

   @constraint(assignment, [i = 1:n, j = 1:n; i!=j], p[i]*h[i] + p[j]*h[j] <=y[i,j,1] + 1)
 
   # ---- OBJEcTIVE 
   
   @objective(assignment, Min, sum(cc .* (y+ x)))
   # if termination_status(assignment)==1
   #     obj_val = objective_value(assignment)
   # else
   #     obj_val = "No_Optimal"
   # end
    #obj_val = objective_value(assignment)
    #x_val = value.(x)
    #y_val = value.(y)
    #h_val = value.(h)
    #w_val = value.(H)
    #println()
    #modelo=assignment
    #println("The model was run!")
    #c,r,cz =H2V_graph(cc,p,g,n)
    #obj_val =custo_Grafo(c,x_val,y_val,R,n)

    #return(obj_val,x_val,y_val,h_val,w_val,modelo)
    return assignment
end

############################### Subproblema   ###############################################
function Subproblem(cc,p,R,n,g,lula=0)
    include("FunctionsH2V.jl")
    m = Cluster_Ground(g)
    assignment = Model(HiGHS.Optimizer)
    #assignment = Model(optimizer_with_attributes(CPLEX.Optimizer,"CPX_PARAM_TILIM" => time_Limit, "CPX_PARAM_THREADS" => 1))

    
    #set_silent(assignment)
    
    @variable(assignment, x[1:n, 1:n, 1:R], Bin)
    @variable(assignment, y[1:n, 1:n, 1:R], Bin)
    @variable(assignment, h[1:n], Bin)
    @variable(assignment, H[1:n, 1:n] ,Bin)
    
   @constraint(assignment, sum(y[:,:,1])>= m*(m-1))
  
    #____________________________________________________________________________________________
    
    #@constraint(assignment, [i = 1:n, j = 1:n; i!=j], p[i]*h[i] + p[j]*h[j] <=y[i,j,1] + 1)
    @constraint(assignment, [i = 1:n, j = 1:n; i!=j], h[i] + h[j] <=sum(y[i,j,:]) + 1)
    @constraint(assignment, [ j = 1:n, r = 1:R], x[j, j,r]+y[j, j,r]+H[j, j]  == 0)
    @constraint(assignment,[i = 1:n, j = 1:n; i!=j], H[i,j]<= (1-h[j]))
    @constraint(assignment, [i = 1:n, j = 1:n; i!=j], h[i] + h[j] + g[i,j] <=sum(y[i,j,1:R]) + 2)
    @constraint(assignment, [i = 1:n], sum((transpose(p).*h).*g[i,:]) >=1)#sum(p))
    @constraint(assignment,[i = 1:n, j = 1:n; i!=j], H[i,j]<= (h[i]))
    @constraint(assignment,[i = 1:n, j = 1:n, r = 1:R; i!=j], 2*y[i,j,r]<= h[i]+h[j])
    @constraint(assignment, sum(H[:,:])== (n-sum(h[:])))
   
    @constraint(assignment, sum(h)>= m)
    #____________________________________________________________________________________________
    
   
    # ---- OBJEcTIVE 
    #@objective(assignment, Min, sum(cc .* y) +sum(cc .* x))
    @objective(assignment, Min, sum(cc .* (y+ x)))
    #optimize!(assignment) 
    if lula==0
        return assignment
    end
    if lula==1
        optimize!(assignment) 
        obj_val = objective_value(assignment)
        x_val = value.(x)
        y_val = value.(y)
        h_val = value.(h)
        w_val = value.(H)
        return(obj_val,x_val,y_val,h_val,w_val,modelo)

    end
    if lula==2
        optimize!(assignment) 
        obj_val = objective_value(assignment)
        x_val = value.(x)
        y_val = value.(y)
        h_val = value.(h)
        w_val = value.(H)
        return(obj_val,x_val,y_val,h_val,w_val,modelo)
    end
end



########################################################################################

function Solution_Construtive(nmr)
    include("FunctionsH2V.jl")
    cc, p,d, nmr, n, R, g = readInstance(nmr)
    h =Matrix{Int}(reshape(zeros(Int(n)),(1,n)))

    gr=[]
    va=[]
    v=[]

    for i in 1:n
        v=[]
        for j in 1:n
            if i<=j && !(j in va) && g[i,:]==g[j,:]
                push!(va,j)
                push!(v,j)
            end
        end
        if v!=[]
            push!(gr,v)
        end
    end
    gr
    flg=0
    for grx in gr
        flg=0
        for j in grx               
            if p[j]==1
                flg=1
            end
        end
        for j in grx
            if (flg==0 && p[j]==1 && h[j]==0)
                h[j]=1
                flg+=1
            end
        end
        
    end
    for i in 1:n
        if h[i]==0 && sum(h)<1
            h[i]=1
        end
    end
    
    #sum(p.*h)

    #h=[0 1 0 0 0 1 1 0 1 1]
    x,y,w,n,R = Cennectionxy02(nmr,h,cc, p, n, R, g)


    obj_val,x_val,y_val,h_val,w_val=sum(cc.*x+cc.*y),x,y,h,w

    #ShowResults(obj_val,x_val,y_val,h_val,w_val,n,R)
    return(obj_val,x_val,y_val,h_val,w_val,gr)
end

function Cennectionxy02(nmr,h,cc, p, n, R, g)
    #cc, p, d, nmr, n, R, g = readInstance(nmr)
    x = Array{Int}(reshape(zeros(n*n*R),(n,n,R)))
    y = Array{Int}(reshape(zeros(n*n*R),(n,n,R)))
    w = Array{Int}(reshape(zeros(n*n),(n,n)))
    gr=[]
    va=[]
    v=[]
    h=h
    for i in 1:n
        v=[]
        for j in 1:n
            if i<=j && !(j in va) && g[i,:]==g[j,:]
                push!(va,j)
                push!(v,j)
            end
        end
        if v!=[]
            push!(gr,v)
        end
    end
    gr
    #flg=0
    #for grx in gr
    #    flg=0
    #    for j in grx               
    #
    #        if (flg==0 && p[j]==1 ) || (d[j]==1)
    #            h[j]=1
    #            flg+=1
    #        end
    #    end
    #
    #end

    for i in 1:n
        for j in 1:n
            if i!=j
                y[i,j,1]=h[i]*h[j]*p[i]*p[j]
            end
        end
    end

    gh=Vector{Int}[]
    gi=Vector{CartesianIndex}[]
    ord=Vector{Int}[]
    ch=cc
    for i in 1:n
        for j in 1:n
            ch[i,j,1]=cc[i,j,1] +200*(1-p[i]*p[j])
        end
    end


    for i in 1:n
        push!(gh,vec(findmin(cc[i,:,:],dims=2)[1]))
        push!(gi,vec(findmin(cc[i,:,:],dims=2)[2]))
    end
    for grx in gr
        for i in grx
            ord=sortperm(gh[i])
            for j in ord
                r = gi[i][j][2]
                if i!=j && g[i,j]==1 && sum(y[i,j,:])==0 && h[i]*h[j]==1
                    if r==1
                        y[i,j,r]=p[i]*p[j]
                    end
                    if r>1 || y[i,j,r]==0
                        y[i,j,2]=1
                    end
                    y[j,i,r]=y[i,j,r]
                end

            end 
        end
    end




    for i in 1:n
        ord=sortperm(gh[i])
        for j in ord
            r = gi[i][j][2]
            #println("h[",j,"]=",h[j]," g[",i,",",j,"]= ",g[i,j])
            if h[j]==1 && i!=j && h[i]==0 && g[i,j]==1 &&  sum(w[i, :])+sum(w[:, i]) ==0 && sum(x[j,i,:])<1
                #x[i,j,r]=1
                if r==1
                    x[j,i,r]=p[i]*p[j]
                else    
                    x[j,i,r]=1
                end
                #w[i,j,r]=1
                w[j,i]=maximum([x[j,i,r],w[j,i]])
            end 
        end 
    end
    for i in 1:n
        for j in 1:n
            if i!=j && g[i,j]==0
                y[i,j,1]=h[i]*h[j]*p[i]*p[j]
            end
        end
    end


    for i in 1:n
        for j in 1:n
            for k in 1:n
                
                if i<j && k!=j && k!=i && g[k,i]*g[k,j]*g[i,j]==1
                    for r in 1:R
                        if sum(x[i,j,:])<1
                            if r==1
                                x[i,j,r]=sum(x[k,i,:])*sum(x[k,j,:])*p[i]*p[j]
                            else
                                x[i,j,r]=sum(x[k,i,:])*sum(x[k,j,:])
                            end
                        end
                    end
                end
            end
        end
    end
    return (x,y,w,n,R,)
end


##########################################################################################


function geraPrimeiraColuna(col,l,ux)
    colu=col
    conta=1
    
    tm=size(colu)[1]
    problemaMestre= Model(HiGHS.Optimizer)
    @variable(problemaMestre,m[1:tm] .>=0)
    
    @constraint(problemaMestre,L[i in 1:tm],l[i]<=m[i]<=ux[i])
    solution_summary(problemaMestre)
    #@constraint(problemaMestre,v,sum(m)==1)
    #@constraint(problemaMestre,m .>=0)
    @objective(problemaMestre,Max,sum(m))
    set_silent(problemaMestre)
    
    try
        optimize!(problemaMestre)    
        objetivo=objective_value(problemaMestre)
        return 1,value.(m),objetivo
    catch
        return 0,[]
    end
    #optimize!(problemaMestre)
    
    #solution_summary(problemaMestre)
    #println("m=",value.(m))
    #u=dual.(L)
    #println("u= ",u)
    #println("v= ",dual(v))
    #println("objetivo Problema mestre = ",sum(obtj.*m))  
    
end


function GeraColumLowBound(i,tempo)

    i=i

    cc, p, d, nmr, n, R, g = readInstance(i)
    

    #obj_val_Solver,x_val,y_val,h_val,w_val,modelo =SolverH2V_03(cc,p,R,n,g,3600)
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
    #obj_val_Solver

    flag,m,objetivo_valor = geraPrimeiraColuna(A,l,ux)
    #objetivo_valor
    #obtj=[309.00]
    obtj=[objetivo_valor]
    col=m
    x, objt, intera= geradorColunas(modelos,s_modelos,obtj,col,tempo)
    #obtj=[1000.0]
    #col
    #obtj
    #objt
    #objt[end]
    #sum(x)
    #minimum(obtj)
    return minimum(obtj),intera
end