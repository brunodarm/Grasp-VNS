
function greedy(c,i)
    lista=[]
    for i in 1:size(c)[1]
        push!(lista,i)
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
        for k in 2:size(lista)[1]
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












