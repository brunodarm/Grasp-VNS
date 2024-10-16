include("FunctionsH2V.jl")
include("H2V_MILP_Function.jl")
include("H2V_Heuristica.jl")
include("Arvore_Minima.jl")
include("GRASP.jl")
using JuMP, NLPModelsJuMP, NLPModels
i=1
cc, p, d, nmr, n, R, g = readInstance(i)
n
cc[1,1,1]

c=reshape(cc,(n*n*R))

modelo=Model(HiGHS.Optimizer)
@variable(modelo, x[1:n*n*R], Bin)
@constraint(modelo, x[1]+2x[2]<=1)

using JuMP, NLPModelsJuMP, NLPModels

nlp = MathOptNLPModel(modelo) # NLPModelsJuMP enters here
x = zeros(nlp.meta.nvar)
cx=grad(nlp, x) # = c = [1, 2, 3, -2, -4, -8]
A=jac(nlp, x) # = A = 5x6 12 entries
gx = cons(nlp, x) # = g = zeros(5)
cx[1]
A[1,200]


#Modelo de extração de coeficientes de um modelo e seus respectivas restrições no modelo A*x <= b
###################################################################################################
    using JuMP, NLPModelsJuMP, NLPModels
    model = Model()
    @variable(model, x[1:3])
    @variable(model, y[1:3])
    @objective(model, Min, sum(x[i] * i - y[i] * 2^i for i = 1:3))
    @constraint(model, [i=1:2], x[i+1] == x[i] + y[i])
    @constraint(model, [i=1:3], x[i] + y[i] <= 1)

    nlp = MathOptNLPModel(model) # NLPModelsJuMP enters here
    x = zeros(nlp.meta.nvar)
    grad(nlp, x) # = c = [1, 2, 3, -2, -4, -8]
    jac(nlp, x) # = A = 5x6 12 entries
    cons(nlp, x) # = g = zeros(5)
    nlp.meta.lcon, nlp.meta.ucon # l <= Ax + g <= u = ([0,0,-Inf,-Inf,-Inf], [0, 0, 1, 1, 1])
    # constraint indexes for each situation
    # [1, 2], [], [3, 4, 5], []
    nlp.meta.jfix, nlp.meta.jlow, nlp.meta.jupp, nlp.meta.jrng
#################################################################################################

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


modelo = Model(HiGHS.Optimizer)#dual_optimizer(HiGHS.Optimizer))#
n=2
u=[]
@variable(modelo,x[1:n])
@constraint(modelo,u1, 2*x[1]+2x[2]>=3)
push!(u,u1)
@constraint(modelo,u2, -2*x[1]+2x[2]<=3)
push!(u,u2)
@constraint(modelo,u3, 2*x[1]+x[2]<=10)
push!(u,u3)
@constraint(modelo,u4, 1<=x[1]<=3)
push!(u,u4)
@constraint(modelo,u5, 1<=x[2]<=3)
push!(u,u5)

@objective(modelo,Min, x[1]-2x[2])
optimize!(modelo)
value.(x)
dual.(u)
A,gx,cx,l,ux =Empacote(modelo)
nx=Int(length(A)/length(gx))
mdl = Model(HiGHS.Optimizer)
@variable(mdl,x_c[1:nx])
Ax=A*x_c

nulo=zeros(length(ux))
eq=Ax+gx

k=0

l
ux
@constraint(mdl,Lx[i=1:length(gx)],l[i]<=eq[i]<=ux[i])

@objective(mdl,Min, cx'*x_c)

optimize!(mdl)
value.(x_c)
dual.(Lx)



######################################################
i=1
cc, p, d, nmr, n, R, g = readInstance(i)
n

obj_val,x_val,y_val,h_val,w_val,modeloas= SolverH2V_02(cc,p,R,n,g,Inf)

A,gx,cx,l,ux =Empacote(modeloas)

Relaxado(A,gx,l,ux,cx)

function Relaxado(A,gx,l,ux,cx)
    nx=Int(length(A)/length(gx))
    mdl = Model(HiGHS.Optimizer)
    @variable(mdl,x_c[1:nx])
    Ax=A*x_c
    eq=Ax+gx
    @constraint(mdl,Lx[i=1:length(gx)],l[i]<=eq[i]<=ux[i])
    #@objective(mdl,Min, cx'*x_c)
    #optimize!(mdl)
    #value.(x_c)
    #dual.(Lx)
    
    return mdl        
end

function Relaxado_03(A,gx,l,ux)
    nx=Int(length(A)/length(gx))
    mdl = Model(HiGHS.Optimizer)
    @variable(mdl,x_c[1:nx])
    Ax=A*x_c
    eq=Ax+gx
    @constraint(mdl,Lx[i=1:length(gx)],l[i]<=eq[i]<=ux[i])
    @objective(mdl,Min,x_c[1])
    #optimize!(mdl)
    #value.(x_c)
    #dual.(Lx)
    
    return mdl        
end



###############################################################################################
A[2331,408]
A[2332,:]
#######################################
A[2331,408]
A[2331,410]
A[2331,280]
ux[2331]
gx[2331]

#h_10 +h_8 - y_10,8,1<=1
####################################

A[2332,410]
A[2332,290]
ux[2332]
gx[2332]

#h_10  - y_10,9,1<=1

####################################

p
# Revisado, faltando colocar Geração de Colunas no formato de Empacote.




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
    @objective(problemaMestre,Min,sum(obtj.*m))
    set_silent(problemaMestre)
    optimize!(problemaMestre)
    
    solution_summary(problemaMestre)
    println("m=",value.(m))
    u=dual.(L)
    println("u= ",u)
    println("v= ",dual(v))
    println("objetivo Problema mestre = ",sum(obtj.*m))
    
    
    return u,dual(v),value.(m), sum(obtj.*value.(m))
    
end


function SubProb(nq,c,u,v,A,modelo)
    v=v
    u=u
    
    modelo_s=modelo
    x=x_c
    set_silent(modelo_s)
    #modelo_s= Model(HiGHS.Optimizer)
    #@variable(modelo_s,x[1:nq].>=0)
    #println(" c =",c," tamanho c = ", size(c)[1])
    #println(" x =",x," tamanho x = ", size(x)[1])
    println(" x = ",x)
    println(" c =",c)
    println(" c'*y= ",c'*x)
    println(" -u'*(A*(x)) = ",-u'*(A*(x)))
    @objective(modelo_s,Min, c'*x -u'*(A*(x)))
    #@objective(subProblema,Min, c'*x -sum(u'*(Aq*(x))))
    optimize!(modelo_s)
    value.(x)
    objective_value(modelo_s)

    ob=c'*value.(x)
    
    #println(" Cy=",ob)
    #println("y= ",value.(y))
    #println(" Cy -v =", -c'*value.(y) +v*1.0)
    #println(" teste",)
    #colu[:,conta]= Aq*value.(y)
    #colu
    f=objective_value(subProblema)
    
    return ob, A*value.(x),value.(x),f
end

#############################################################################################

########################################################################
function geradorColunas(modeloMestre,subProblema,ob_inicial,colu_inicial)
   
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
    modelo_s =Relaxado(Aq,gxq,lq,uxq,cxq)
    while (permi==".")
        
         
        u,v,my,resu = ProblemaMestreRestrito(colu,l,ux,obtj,conta)
        
        push!(mx,[my])
        push!(rp,resu)
        push!(vx,v)
        println("LS =",rp[end])

        LS=vx[end]
        if conta>1
            LI=rp[end-1] + custoRed
            println("LI =",LI)
            
        else
        
            println("LI =",-Inf)
            LI=-Inf
        end
        
        
        println("conta = ",conta)
        ob, Aa,yi,f =SubProb(nq,cx,u,v,A,modelo_s)
        

        push!(obtj,ob)
        println("###########################################")
        println(" f =",f," and v = ",v)
    
        conta=conta+1
        
        println(" CUSTO RELATIVO ",conta-1," = ",f-v)
        if (f-v)<0 && conta<100
            permi="."
            colu=hcat(colu,Aa)
            push!(yv,yi)
            println("yv =",yv)
        else
            permi="para"

        end
        #push!(custoRed,f-v)
        custoRed=f-v
    end
   
    somas=sum(yv[end-i+1].*mx[end][end][end-i+1] for i=1:length(yv))
    println("x=",somas)
    println("objetivo =",rp[end])

    return somas,rp[end]

end




###################### EXEMPLO   #########################

modeloMestre= Model(HiGHS.Optimizer)
n=2
@variable(modeloMestre,x[1:n])
@constraint(modeloMestre,u1, 2*x[1]+2x[2]>=3)
@constraint(modeloMestre,u2, -2*x[1]+2x[2]<=3)
@constraint(modeloMestre,u3, 2*x[1]+x[2]<=10)

@objective(modeloMestre,Min, x[1]-2x[2])
A,gx,cx,l,ux =Empacote(modeloMestre)


subProblema = Model(HiGHS.Optimizer)
n=2
@variable(subProblema,x[1:n])
@constraint(subProblema,u4, 1<=x[1]<=3)
@constraint(subProblema,u5, 1<=x[2]<=3)
#@constraint(subProblema,u6, 1*x[1]+3x[2]<=10)
As,gxs,cxs,ls,uxs =Empacote(subProblema)

#obtj=[200.0]
#colu=[12;0;0]



A
flag,m,objetivo_valor = geraPrimeiraColuna(A,l,ux)
objetivo_valor
obtj=[objetivo_valor]
col=m
x, objt= geradorColunas(modeloMestre,subProblema,obtj,col)

obtj
objt
objt[end]







################################ Problema Restrito ###########################################
function ProbRestrito(cc,p,R,n,g)
    include("FunctionsH2V.jl")
    m = Cluster_Ground(g)
    assignment = Model(HiGHS.Optimizer)
    #assignment = Model(optimizer_with_attributes(CPLEX.Optimizer,"CPX_PARAM_TILIM" => time_Limit, "CPX_PARAM_THREADS" => 1))
    #set_time_limit_sec(assignment,time_Limit)
    
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

    
   
   
   
    # ---- OBJEcTIVE 
    #@objective(assignment, Min, sum(cc .* y) +sum(cc .* x))
    @objective(assignment, Min, sum(cc .* (y+ x)))
    #optimize!(assignment) 
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
    println("The model was run!")
    #c,r,cz =H2V_graph(cc,p,g,n)
    #obj_val =custo_Grafo(c,x_val,y_val,R,n)

    #return(obj_val,x_val,y_val,h_val,w_val,modelo)
    return assignment
end

############################### Subproblema   ###############################################
function Subproblem(cc,p,R,n,g)
    include("FunctionsH2V.jl")
    m = Cluster_Ground(g)
    assignment = Model(HiGHS.Optimizer)
    #assignment = Model(optimizer_with_attributes(CPLEX.Optimizer,"CPX_PARAM_TILIM" => time_Limit, "CPX_PARAM_THREADS" => 1))
    #set_time_limit_sec(assignment,time_Limit)
    
    set_silent(assignment)
    @variable(assignment, x[1:n, 1:n, 1:R], Bin)
    @variable(assignment, y[1:n, 1:n, 1:R], Bin)
    @variable(assignment, h[1:n], Bin)
    @variable(assignment, H[1:n, 1:n], Bin)


    @constraint(assignment, [i = 1:n, j = 1:n; i!=j], p[i]*h[i] + p[j]*h[j] <=y[i,j,1] + 1)
    #--------------------------------------
    @constraint(assignment,[i = 1:n, j = 1:n, r=1:R],0<=(x[i,j,r]) <=1)
    @constraint(assignment,[i = 1:n, j = 1:n, r=1:R],0<=(y[i,j,r]) <=1)
    @constraint(assignment,[i = 1:n, j = 1:n, r=1:R],0<=(H[i,j]) <=1)
    @constraint(assignment,[i = 1:n],(h[i]) <=1)

    # ---- OBJEcTIVE 
    #@objective(assignment, Min, sum(cc .* y) +sum(cc .* x))
    @objective(assignment, Min, sum(cc .* (y+ x)))

    return assignment
end



########################################################################################

for i in 1:2
    baixo=0
    alto=0

    if round(x[i])-x[i]<0
        baixo=Int(round(x[1]))
        alto=Int(baixo+1)
    else
        alto=Int(round(x[i]))
        baixo=Int(alto-1)
    end

    alto
    baixo
end


###############################################################################

################## Gerar solução qualquer #####################################


i=1
cc, p, d, nmr, n, R, g = readInstance(i)
n

#obj_val,x_val,y_val,h_val,w_val,modeloas= SolverH2V_02(cc,p,R,n,g,Inf)
modelos=ProbRestrito(cc,p,R,n,g)
A,gx,cx,l,ux =Empacote(modelos)
modeloMestre=Relaxado(A,gx,l,ux,cx)

mdl=Relaxado_02(A,gx,l,ux)
optimize!(mdl)

value.(x_c)
sum(value.(x_c))

A*value.(x_c)





###################################
    

colu=colu
    conta=size(colu)[2]



    obtj=obtj
    tm=size(colu)[1]
    problemaMestre= Model(HiGHS.Optimizer)
    @variable(problemaMestre,m[1:conta] .>=0)
    eq=colu*m
    @constraint(problemaMestre,L[i=1:tm],l[i]<=colu[i,:]'*m<=ux[i])
    @constraint(problemaMestre,v,sum(m)==1)
    @constraint(problemaMestre,m .>=0)
    @objective(problemaMestre,Min,sum(obtj.*m))
    set_silent(problemaMestre)
    optimize!(problemaMestre)
    
    solution_summary(problemaMestre)
    println("m=",value.(m))
    u=dual.(L)
    println("u= ",u)
    println("v= ",dual(v))
    println("objetivo Problema mestre = ",sum(obtj.*m))

#####################################

modelos=ProbRestrito(cc,p,R,n,g)
A,gx,cx,l,ux =Empacote(modelos)
modeloMestre=Relaxado(A,gx,l,ux,cx)

xo = zeros(size(A)[2])
wo=5
xo[wo]=0
xo[wo]=1

Ao=A*xo
zo=0
for i in 1:size(A)[2]
    
    
end
println(zo)
Ao[i] <=ux[i]
Ao<=ux
randn(4)


xo = zeros(size(A)[2])
qo=[]
zo=0
zi=[10000]
#xo[1]=1
Ao=A*xo
for i in 1:size(A)[2]
    zo=0
    xo[i]=1
    Ao=A*xo
    for j in 1:size(A)[2]
        if !(Ao[j] <=ux[j] && Ao[j] >=l[j])
            zo=zo+1
        end
    end
    if(zi[end]<zo)
        xo[i]=0
    end
    push!(zi,zo)
end
zi
argmin(zi)
sum(xo)
argmax(ux)


xo=rand(size(A)[2])
Ao=A*xo
Ao <=ux && Ao >=l



##################################################################


col=Matrix(reshape(zeros(size(A)[2]),(size(A)[2],1)))
col[2]=1
flag,m = geraPrimeiraColuna(A,l,ux)
m

sum(m)
A*m

function geraPrimeiraColuna(col,l,ux)
    colu=col
    conta=1
    
    tm=size(colu)[1]
    problemaMestre= Model(HiGHS.Optimizer)
    @variable(problemaMestre,m[1:tm] )
    
    @constraint(problemaMestre,L[i in 1:tm],l[i]<=m[i]<=ux[i])
    solution_summary(problemaMestre)
    #@constraint(problemaMestre,v,sum(m)==1)
    #@constraint(problemaMestre,m .>=0)
    @objective(problemaMestre,Max,m[1])
    set_silent(problemaMestre)
    try
        optimize!(problemaMestre)    
        return 1,value.(m)
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


colu=col
    conta=1
    
    tm=size(colu)[1]
    problemaMestre= Model(HiGHS.Optimizer)
    @variable(problemaMestre,m[1:tm] )
    
    @constraint(problemaMestre,L[i in 1:tm],l[i]<=colu[i]'.*m[i]<=ux[i])
    solution_summary(problemaMestre)
    #@constraint(problemaMestre,v,sum(m)==1)
    #@constraint(problemaMestre,m .>=0)
    @objective(problemaMestre,Max,sum(m))
    set_silent(problemaMestre)
    try
        optimize!(problemaMestre)    
        return 1,value.(m)
    catch
        return 0,[]
    end



    import Pkg; Pkg.add("Karnak")


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
#obj_val_Solver

flag,m,objetivo_valor = geraPrimeiraColuna(A,l,ux)
#objetivo_valor
#obtj=[309.00]
obtj=[objetivo_valor]
col=m
x, objt, interacoes= geradorColunas(modelos,s_modelos,obtj,col,100)
#obtj=[1000.0]
col
obtj
objt
objt[end]
sum(x)
minimum(obtj)

interacoes
i=1
tempo=100

mmx,interacoes = GeraColumLowBound(i,tempo)










modelos=ProbRestrito(cc,p,R,n,g)
A,gx,cx,l,ux =Empacote(modelos)
modeloMestre= Relaxado_02(A,gx,l,ux,cx)
modelos


#@objective(modelos,Min, sum(cx'.*x_c))
optimize!(modelos)

objective_value(modelos)


optimize!(modeloMestre)
#@objective(modeloMestre,Min, sum(cx'.*x_c))
objective_value(modeloMestre)
argmax(value.(x_c))
#value.(x_c[220])
x_c.modeloMestre
objective_value(modeloMestre)

A[end,:]
l[end]
ux[end]

    nx=Int(length(A)/length(gx))
    mdl = Model(HiGHS.Optimizer)
    @variable(mdl,x_c[1:nx]>=0)
    Ax=A*x_c
    eq=Ax+gx
    @constraint(mdl,Lx[i=1:size(A)[1]],l[i]<=eq[i]<=ux[i])
    #cx'*x_c
    @objective(mdl,Min, (cx'*x_c))
    optimize!(mdl)
    objective_value(mdl)






flag,m,objetivo_valor = geraPrimeiraColuna(A,l,ux)
objetivo_valor
#obtj=[objetivo_valor]
col=m
x, objt= geradorColunas(modelos,s_modelos,obtj,col)