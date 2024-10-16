include("Arvore_Minima.jl")
include("FunctionsH2V.jl")
include("ShowResults.jl")
include("H2V_MILP_Function.jl")
include("GRASP.jl")
using Plots, GraphRecipes
c =[00 30 30 02 06 30 30 30 30 30 06;
    30 00 30 30 30 30 06 03 06 30 30;
    30 30 00 30 30 02 30 30 30 06 30;
    02 30 30 00 02 08 30 06 30 30 03;
    06 30 30 02 00 30 30 30 30 30 06;
    30 30 02 08 30 00 30 08 30 02 30;
    30 06 30 30 30 30 00 02 30 30 30;
    30 03 30 06 30 08 02 00 03 30 30;
    30 06 30 30 30 30 30 03 00 30 30;
    30 30 06 30 30 02 30 30 30 00 30;
    06 30 30 03 06 30 30 30 30 30 00
]

nmr=1
cx,px, dx ,nmr, nx, Rx,gx = readInstance(nmr)
n=nx
c=cx
#ct,r,cz =H2V_graph(cx,px,gx,nx)
c
grd=gx
aMin,aMax=0.55,0.85
menor,g,h = H2V_Grasp(cx,grd,px,0,0)
tempo=@elapsed menor,g,h = H2V_Grasp(cx,grd,px,aMin,aMax)
menor
tempo
Graph_Plot(g,[],[],h,n)
gt=Unado(nx)
gt=reshape(zeros(n*n),(n,n))
obj_val,x_val,y_val,h_val,w_val,modelo=SolverH2V_02(cx,px,Rx,nx,gt,nx)
obj_val

Graph_Plot(g,x_val,y_val,h,n)

m = Cluster(gx)
sum(gx[1,:])

c,r,cz =H2V_graph(cx,px,grd,n)
#p=nh
time= @elapsed st=arvoreGeradora(c,2)

Unado(10)

#LST,st,grafo =arvoreMin(c)

g=Array{Int}(reshape(zeros(n*n),(n,n)))

h = hub_selection(4,st,grd,c,px,85,1)
sum(h)
hg,g3,h_gr = X_acess(h,c,gx,px)
g1=conectaGraph(h_gr,nx)
g2=hubAcess(h,g,grd,px)
sum(c.*g3)
a=(1,4)
#g2[a[1],a[2]]=1
#g2[a[2],a[1]]=1
g3=g1+g2
Graph_Plot(g3,[],[],h,n)

################################################################################
c,r,cz =H2V_graph(cx,px,grd,n)
#h=[]
#h=[1 1 0 0 0 1 1 0 0 1]
n=size(h)[1]
g=Matrix{Int}(reshape(zeros(n*n),(n,n)))
cx=cx
hg=[]
ht=[]
vt=[]
lista=Array{Int}(1:n)
for i in 1:n
    if h[i]==1
        push!(ht,i)
    else
        push!(vt,i)
    end
end
ht
vt
filter!(e->!(e in ht),lista)

size(vt)[1]
while size(vt)[1]>0
    j=rand(vt)
    ca =testeConexao(cx,j,ht,gx,px,h)
    i=argmin(ca)[1]
    r=argmin(ca)[2]
    println(j,", ",ca," argmin=",ht[i]," -- ",ca[argmin(ca)])
    if ca[argmin(ca)]<9999
        push!(hg,[ht[i],j,r])
        g[ht[i],j]=1
        g[j,ht[i]]=1
    end
    vt=selecao(j,vt)
   
end


g
hg
Graph_Plot(g,[],[],h,n)
#################################################################################



######################################################################################
c,r,cz =H2V_graph(cx,px,grd,n)
g=Matrix{Int}(reshape(zeros(n*n),(n,n)))
c=c
stx=[]
n=Int(n)

for i in 1:n
    stt,st,g0=arvoreMin(c,i)
    for j in st
        push!(stx,j)
    end
    g=g+g0
end
px
h=[1 1 0 0 0 1 1 0 0 1]
FreqNode(stx,n)
g
g[10,1]=0
g[1,10]=0
g[10,8]=1
g[8,10]=1

g=grd.*g
Graph_Plot(g,[],[],h,n)

g[1,4]
gr =[1,3,4,6,8,9]
stt,st,g0=arvoreMin(c[gr,gr],1)
Graph_Plot(g0,[],[],h[gr],size(g0)[1])

#########################################################################################