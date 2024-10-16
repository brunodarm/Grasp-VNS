
function greedy(c,i,l)
    l=l
    c=c
    local lista=[]
    for j in 1:size(c)[1]
        push!(lista,j)
    end
    x=i
    resulta=[]
    push!(resulta,x)
    custo=0
    for i in 2:size(lista)[1]
        lista=selecao(x,lista)
        xi=copy(x)
        x=lista[argmin(c[x,lista])]
        custo=custo+c[xi,x]
        push!(resulta,x)
    end
    push!(resulta,l)
    return resulta,custo+c[x,l]
    
end

function greedy_2(c,i)

    c=c
    local lista=[]
    for j in 1:size(c)[1]
        push!(lista,j)
    end
    x=i
    resulta=[]
    push!(resulta,x)
    custo=0
    for i in 2:size(lista)[1]
        lista=selecao(x,lista)
        xi=copy(x)
        x=lista[argmin(c[x,lista])]
        custo=custo+c[xi,x]
        push!(resulta,x)
    end
    push!(resulta,i)
    return resulta,custo+c[x,i]
    
end





function GRASP_02(c,i,a,n)
    n=n
    c=copy(c)
    lista=[]
    for j in 1:size(c)[1]
        push!(lista,j)
    end
    
    x=i
    resulta=[]
    push!(resulta,x)
    custo=0
    theBest=[],99999999
    j=1;
    for j in 1:n
        #x=i
        #resulta=[]
        #push!(resulta,x)
        #custo=0
        #theBest=[],99999999
        for kl in 2:size(lista)[1]
            lista=selecao(x,lista)
            xi=copy(x)
            d=sortperm(c[x,lista])
            xc = c[x,lista][d[1:(size(d)[1]-1)]]
            k=maximum(xc)-minimum(xc)
            axc=minimum(xc)
            x=axc
            if rand()<=a
                axc=axc+k*a
                xcc=xc[findall(e->e<=axc,xc)]
                x=rand(d[1:size(xcc)[1]])
            end  
            custo=custo+c[xi,x]
            push!(resulta,x)
            
        end
        if theBest[2]> custo+c[x,i] 
            push!(resulta,i)
            theBest= resulta,custo+c[x,i] 
        end
    end
   
    return Array{Int}(theBest[1]),theBest[2]
    
    
end

function GRASP(c,i,a,n)
    lista=[]
    for i in 1:size(c)[1]
        push!(lista,i)
    end
    x=i
    resulta=[]
    push!(resulta,x)
    custo=0
    theBest=[],99999999
    for j in 1:n
        for kl in 2:size(lista)[1]
            lista=selecao(x,lista)
            xi=copy(x)
            #x=lista[argmin(c[x,lista])]
            v = copy(c[x,lista])
            filter!(e->e<99,v)
            k=maximum(v)-minimum(v)
            maxV=minimum(v)+k*a
            vv = copy(v)
            filter!(e->e<=maxV,vv)
            #rand(vv)
            x = lista[indexin(rand(vv),c[x,lista])[1]]
            custo=custo+c[xi,x]
            push!(resulta,x)
        end
        
        if theBest[2]> custo+c[x,i] 
            push!(resulta,i)
            theBest= resulta,custo+c[x,i] 
        end
    end
    
    return Array{Int}(theBest[1]),theBest[2]
end



function evoGRASP(c,i,aMin,aMax,n)
    list_alpha = geraAlpha(n,aMin,aMax)
    b=copy(list_alpha)
    theBest=[],99999999
    
    for k in 1:n
        a=rand(b)
        selecao(a,b)
        x = GRASP(c,i,a,n)
        if theBest[2]> x[2] 
            theBest= x 
        end
    end
    return Array{Int}(theBest[1]),theBest[2]
end






function selecao(x,lista)
    if (x in lista)
        filter!(e->e!=x,lista)
    end
    return lista
end








function geraCusto(n)

    c=reshape(zeros(n*n),(n,n))
    for i in 1:n
        for j in 1:n
            if i<j 
                c[i,j]=rand(10:50) 
            end
            if i==j
                c[i,i]=99
            end   
            if i>j
                c[i,j]=c[j,i]
            end
                   
        end
    end
    return c
end

function geraAlpha(n,nMin,nMax)
    num_inter=n
    numMin=nMin
    numMax=nMax
    listaAlpha = []
    x=0
    for i in 1:num_inter
        flg=false
        while flg ==false
            x =rand()
            flg=false
            if x>=numMin && x<=numMax
                flg=true
            end
        end
        push!(listaAlpha,x)
    
    end

    return listaAlpha  
end

function H2V_Grasp(cx,grd,px,aMin,aMax)
    n=size(grd)[1]
    cx=cx
    px=px
    c,r,cz =H2V_graph(cx,px,grd,n)
    st=arvoreGeradora(c,round(log10(n))+1)
    #LST,st,grafo =arvoreMin(c)
    n_a=round(log10(n))
    k_a=4
    k_h=10
    list_alpha =[]
    if (aMin<=0 && aMax==0) || aMax<=0
        list_alpha =[]
    else    
        list_alpha = geraAlpha(n_a+k_a,aMin,aMax)
        
    end
    b=copy(list_alpha)
    
    gm=Matrix{Int}(reshape(zeros(Int(n*n)),(n,n)))
    hm=Array{Int}(zeros(n))

    
    t1=@elapsed ns, gm, menor,hm=Heuristica_GRASP(cx,c,grd,px,st,0)
    a=0
    
    if size(list_alpha)[1]>0
        for k in 1:(n_a+k_a)
            if size(b)[1]>0
                a=rand(b)
                selecao(a,b)
            end
            
            
            for i in 1:(n_a+k_h)
                t2=@elapsed ns, g, custo,h=Heuristica_GRASP(cx,c,grd,px,st,a*100)
                t1=t2+t1
                if custo<menor 
                    menor=copy(custo)
                    gm=copy(g)
                    hm=copy(h)
                end
                if t1>=3600
                    break
                end
            end
            if t1>=3600
                break
            end
            #menor,x,y,hm=custoso(gm,cx,px,grd,hm)
        end
        
    end
    menor,x,y,hmm=custoso(gm,cx,px,grd,hm)
    return menor,gm,hm


    
end



function GraspSequence(stx,n,alfa)
    stx=stx
    n=n
    sty=copy(stx)
    fq = FreqNode_02(stx)
    xf=copy(fq)
    #fqq= FreqNode_02(fq[2,:])
    sequeny=[]
    alfa=alfa
    for h in 1:n
        iz=reverse(sortperm(xf[2,:]))
        xz=copy(xf[2,iz])
        xx=copy(xz )

        xx
        k=maximum(xx) -minimum(xx)
        k*alfa
        yy=maximum(xx)-k*alfa
        filter!(e->e>=yy,xx)
        xx
        aa=rand(xx)
        ab=findfirst(e->e==aa,xf)[2]
        push!(sequeny,xf[1,ab])
        selecao(xf[1,ab],sty)
        xf=FreqNode_02(sty)
    end

    return sequeny

    
end







