using JuMP, Dualization
import DataFrames
import HiGHS
import Plots
import SparseArrays

modelo = Model(HiGHS.Optimizer)#dual_optimizer(HiGHS.Optimizer))#
n=2

@variable(modelo,x[1:n])
@constraint(modelo,u1, 2*x[1]+2x[2]>=3)

@constraint(modelo,u2, -2*x[1]+2x[2]<=3)

@constraint(modelo,u3, 2*x[1]+x[2]<=10)

@constraint(modelo,u4, 1<=x[1]<=3)

@constraint(modelo,u5, 1<=x[2]<=3)


@objective(modelo,Min, x[1]-2x[2])




#model=dualize(modelo)




optimize!(modelo)
value.(x)
value.(u)
objective_value(modelo)


modelo = Model(HiGHS.Optimizer)#dual_optimizer(HiGHS.Optimizer))#
n=3
u=[]
@variable(modelo,x[1:n])
@constraint(modelo,u1, 2*x[1]-2x[2]+2x[3]<=1)
push!(u,u1)
@constraint(modelo,u2, 2*x[1]+2x[2]+x[3]<=-2)
push!(u,u2)

@objective(modelo,Min, 3*x[1] -3x[2] -10x[3])




#model=dualize(modelo)




optimize!(modelo)
value.(x)
value.(u)
objective_value(modelo)

b=[3 ;-3; -10;-3;1;-3;1]
A=[[2 2];[2 -2];[-2 -1];[-1 0];[1 0];[0 -1];[0 1]]
c=[1 ;-2]

modelo = Model(HiGHS.Optimizer)#dual_optimizer(HiGHS.Optimizer))#
n=2
L=[1 ;2 ;3;4;5]
@variable(modelo,x[1:n])
@constraint(modelo,L[i in 1:length(b)], A[i,:]'*x .>= b[i])


@objective(modelo,Min, c'*x)
optimize!(modelo)
solution_summary(modelo)
u=[]
for i in 1:length(b)
    push!(u,dual(L[i]))
end
value.(x)
u



subProblema= Model(HiGHS.Optimizer)
@variable(subProblema,y[1:n])
@constraint(subProblema,y.>=0)
@objective(subProblema,Min, c'*y -sum(u.*(A*y)))

c'*y
u'
(A*y)
c'*y +sum(u.*(A*y))

#############################################################################################

using JuMP, Dualization
import DataFrames
import HiGHS
import Plots
import SparseArrays

#obtj
b=[3 ;3; 10;3;1;3;1]
A=[[2 2];[-2 +2];[2 +1];[1 0];[1 0];[0 1];[0 1]]
c=[1 ;-2]
fator=[1;-1;-1;-1;1;-1;1]
mp=1:3
sp=4:7
#(colu.*fator[1:tm])[1,1:conta]
bq=b[mp]
Aq=A[mp,:]
conta=1
bs=b[sp]
As=A[sp,:]
obtj=[1000.0]
n_col=(length(mp))
colu=reshape(zeros(n_col^2),(n_col,n_col))
n_max=1
colu[:,1]=[3;0;0]
L=[1:(length(mp))]
# Problema Metre Restrito #####################################################################
 
 bq[3].*fator[3]
 conta
 colu[:,end]
 colu
 length(colu)
 tm=size(colu)[1]
 problemaMestre= Model(HiGHS.Optimizer)
 @variable(problemaMestre,m[1:conta] .>=0)
 @constraint(problemaMestre,L[i in 1:tm],(colu.*fator[1:tm])[i,1:conta]'*m[1:conta] .>=bq[i].*fator[i])
 #@constraint(problemaMestre,L[1],colu[1,1:conta]'*m[1:conta] .>=bq[1])
 

 @constraint(problemaMestre,v,sum(m)==1)
 @constraint(problemaMestre,m .>=0)
 @objective(problemaMestre,Min,sum(obtj.*m))
 optimize!(problemaMestre)
 solution_summary(problemaMestre)
 π=0
 value.(m)
 u=dual.(L)
 dual(v)
 conta=conta+1
 objective_value(problemaMestre)
 colu[1,1:conta]
 #colu[1,1:conta]'*m[1:conta]
 sum(obtj.*value.(m))
 println("LS = ",objective_value(problemaMestre))
 println("LI =",objective_value(problemaMestre) - (obtj[end-1]-dual(v)))
 (obtj[end-1]-dual(v))
# Subproblema ###############################################################################
 n=2
 subProblema= Model(HiGHS.Optimizer)
 @variable(subProblema,y[1:n].>=0)
 @constraint(subProblema,y.>=0)
 @constraint(subProblema, 1<=y[1]<=3)
 @constraint(subProblema, 1<=y[2]<=3)
 @objective(subProblema,Min, c'*y -sum(u'*(Aq*y)))
 optimize!(subProblema)
 value.(y)
 objective_value(subProblema)

 ob=c'*value.(y)
 
 push!(obtj,ob)
 colu
 [ Aq*value.(y)]
 [ Aq*y]
 value.(y)
 colu[:,conta]= Aq*value.(y)
 colu

 u'*(Aq*y)
 "x=" 
 
 
 #######################################################

 modelo = Model(HiGHS.Optimizer)
 @variable(modelo,x >=0)
 @objective(modelo,Min, 1000*x)
 @constraint(modelo,ll,3x>=3)
 @constraint(modelo,x=1)
 optimize!(modelo)
 value.(x)
 dual(ll)


 Matr = Matrix{Float64}
 Matr=[ 1 2; 3 4]
 vet= [5 ;6]
 hcat(Matr,vet)

 colu=hcat(colu,Aq*value.(y))
colu


function ProblemaMestreRestrito(colu,conta,fator,bq,obtj)
    colu=colu
    conta=conta
    obtj=obtj
    tm=size(colu)[1]
    problemaMestre= Model(HiGHS.Optimizer)
    @variable(problemaMestre,m[1:conta] .>=0)
    @constraint(problemaMestre,L[i in 1:tm],(colu.*fator[1:tm])[i,1:conta]'*m[1:conta] .>=bq[i].*fator[i])
    @constraint(problemaMestre,v,sum(m)==1)
    @constraint(problemaMestre,m .>=0)
    @objective(problemaMestre,Min,sum(obtj.*m))
    set_silent(problemaMestre)
    optimize!(problemaMestre)
    
    solution_summary(problemaMestre)
    println("m=",value.(m))
    u=dual.(L)
    println("u= ",u.*fator[1:tm])
    println("v= ",dual(v))
    println("objetivo Problema mestre = ",sum(obtj.*m))
    
    
    return u.*fator[1:tm],dual(v),value.(m), sum(obtj.*value.(m))
    
end


function SubProb(modelo,c,Aq,u,v)
    v=v
    u=u
    subProblema=modelo
    set_silent(subProblema)
    @objective(subProblema,Min, c'*y -sum(u'*(Aq*(y))))
    optimize!(subProblema)
    value.(y)
    objective_value(subProblema)

    ob=c'*value.(y)
    
    #println(" Cy=",ob)
    #println("y= ",value.(y))
    #println(" Cy -v =", -c'*value.(y) +v*1.0)
    #println(" teste",)
    #colu[:,conta]= Aq*value.(y)
    #colu
    f=objective_value(subProblema)
    
    return ob, Aq*value.(y),value.(y),f
end

#############################################################################################
1+1


using JuMP, Dualization
import DataFrames
import HiGHS
import Plots
import SparseArrays

#obtj

b=[3 ;3; 10;3;1;3;1]
A=[[2 2];[-2 +2];[2 +1];[1 0];[1 0];[0 1];[0 1]]
c=[1 ;-2]
fator=[1;-1;-1;-1;1;-1;1]

mp=1:3
sp=4:7
fatorq=fator[mp]
bq=b[mp]
Aq=A[mp,:]
conta=1
bs=b[sp]
As=A[sp,:]
obtj=[200.0]
colu=[12;0;9]
permi="."

n=2
subProblema= Model(HiGHS.Optimizer)
@variable(subProblema,y[1:n])
@constraint(subProblema,y.>=0)
@constraint(subProblema, 1<=y[1]<=3)
@constraint(subProblema, 1<=y[2]<=3)
#@objective(subProblema,Min, c'*y -sum(u'*(Aq*y)))
yv=[]
mx=[]
conta=1
vx=[]
LSS=[]
rp=[]
custoRed=0.0
u=[]
while (permi==".")
    

    u,v,my,resu = ProblemaMestreRestrito(colu,conta,fator,bq,obtj)
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
    ob, Aa,yi,f =SubProb(subProblema,c,Aq,u,v)
    

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
println("Fim!!!!!!!!!!!!!!!!!!!!!!!")
colu
mx[end][end]
colu
mx[end][end]
(colu*mx[end][end])[1,:]
bq
conta
yv
value.(y)
obtj
mx
yv
vx
rp
soma=[0,0]
for i in (1:length(yv))
    soma=yv[end-i+1].*mx[end][end][end-i+1]+soma
end
soma
rp[end]

colu
size(colu)[2]

yv.*mx[end]
mx[end][1][1]

for i in 1:3
    zz=yv[i].*mx[end-1][1][i]
    println(zz)
end
zz
########################################################################
function geradorColunas(subProblema,c,Aq,bq,fator,ob_inicial,colu_inicial)
    conta=1
    obtj=ob_inicial
    bq=bq
    #colu=zeros(length(Aq))
    fator=fator
    colu=colu_inicial
    subProblema=subProblema
    permi=="."
    yv=[]
    mx=[]
    conta=1
    vx=[]
    rp=[]
    custoRed=0.0
    u=[]
    somas=[0,0]
    while (permi==".")
        
         
        u,v,my,resu = ProblemaMestreRestrito(colu,conta,fator,bq,obtj)
        
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
        ob, Aa,yi,f =SubProb(subProblema,c,Aq,u,v)
        

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
    
    for i in (1:length(yv))
        somas=yv[end-i+1].*mx[end][end][end-i+1]+somas
    end
    println("x=",soma)
    println("objetivo =",rp[end])

    return soma,rp[end]

end

b=[3 ;3; 10;3;1;3;1]
A=[[2 2];[-2 +2];[2 +1];[1 0];[1 0];[0 1];[0 1]]
c=[1 ;-2]
fator=[1;-1;-1;-1;1;-1;1]

mp=1:3
sp=4:7
fatorq=fator[mp]
bq=b[mp]
Aq=A[mp,:]
conta=1
bs=b[sp]
As=A[sp,:]
obtj=[200.0]
colu=[12;0;9]


using JuMP, Dualization
import DataFrames
import HiGHS
import Plots
import SparseArrays

geradorColunas(subProblema,c,Aq,bq,fator,obtj,colu)
