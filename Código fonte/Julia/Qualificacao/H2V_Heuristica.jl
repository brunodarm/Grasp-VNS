

function H2V_Heuristic_Construtive(nmr)
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
    x,y,w,n,R = Cennectionxy(nmr,h,cc, p, n, R, g)


    obj_val,x_val,y_val,h_val,w_val=sum(cc.*x+cc.*y),x,y,h,w

    #ShowResults(obj_val,x_val,y_val,h_val,w_val,n,R)
    return(obj_val,x_val,y_val,h_val,w_val,gr)
end


function H2V_Heuristic_Best()
    for i in 0:(Int(exp2(10))-1)
        vet = Transpose(digits(Int(exp2(10))-i-1,base=2,pad=10))
    
        if vet == h
            println(vet)
        end
    end
    
end


function Heuristic_H2V_V01(nmr)

    include("FunctionsH2V.jl")
    cc, p,d, nmr, n, R, g = readInstance(nmr)
    h =Matrix{Int}(reshape(zeros(Int(n)),(1,n)))
    gr=[]
    va=[]
    v=[]
    for i in 1:n
        v=[]
        for j in 1:n
            if (i<=j) && !(j in va) && (g[i,:]==g[j,:])
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

    #x,y,w,n,R = Cennectionxy(nmr,h,cc, p, n, R, g)
    #obj_val,x_val,y_val,h_val,w_val=sum(cc.*x+cc.*y),x,y,h,w
    ##ShowResults(obj_val,x_val,y_val,h_val,w_val,n,R)
    #return(obj_val,x_val,y_val,h_val,w_val,gr)
    return gr, h
end