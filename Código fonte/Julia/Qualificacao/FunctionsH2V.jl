function processingGenerator(R,n,cont,nmr)
    #output_path =("/home/bruno/Dropbox/Pesquisa/Problemas/Green parallel machines/Source/Instances generator/Testbed_tardiness/")
    output_path =("C:/Users/User/Instancias/H2V")
    instance = string(output_path,"_r",R,"_n",n,"_inst_", cont,"_num",nmr, ".txt")
    #instance = string("p", cont, ".txt")
    #instance = "p" * string(cont) * ".txt"
    process = open(instance, "w")
    global m = reshape(Vector(1:(R*n*n)),(n,n,R))
    for r in 1:R
        for i in 1:n
            for j in 1:n
                m[i,j,r]  =0
            end
        end
    end
    
    for r in 1:R
        for i in 1:n
            for j in 1:n
                if i<j 
                    m[i,j,r] = rand(1:100)
                    print(process, m[i,j,r], "   ")
                end
                if i==j
                    m[i,j,r]=0
                    print(process, m[i,j,r], "   ")
                end
                if j<i
                    m[i,j,r] = m[j,i,r]
                    print(process, m[i,j,r], "   ")
                end

            end
        end
    end
    print(process,"\n")
    close(process)
  # Gerando as instancias Portuarias 
    g,p =  groundGenerator(n)
    instance1 = string(output_path,"_p", nmr, ".txt")
    process1 = open(instance1, "w")
    for i in 1:n
       # pp = rand(0:1)
        print(process1, p[i], "   ")
    end
    print(process1,"\n")
    close(process1)

 # Gerando as instancias de Demanda 
    instance2 = string(output_path,"_d", nmr, ".txt")
    process2 = open(instance2, "w")
    for i in 1:n
        dd = rand(0:1)
        print(process2, dd, "   ")
    end
    print(process2,"\n")
    close(process2)
 # Gerando as instancias de ilhas 
    instance3 = string(output_path,"_g", nmr, ".txt")
    process3 = open(instance3, "w")
    
    for i in 1:n
        for j in 1:n
            print(process3, g[i,j], "   ")
        end
    end
    print(process3,"\n")
    close(process3)
end
    

using DelimitedFiles
function readProcessing(R,n,cont,nmr)
    #input_path =("/home/bruno/Dropbox/Pesquisa/Problemas/Green parallel machines/Source/Instances generator/Testbed_tardiness/")
    input_path =("C:/Users/User/Instancias/H2V")
    process_instance_Custo = string(input_path,"_r",R,"_n",n,"_inst_", cont,"_num",nmr, ".txt")
    process_instance_Porto = string(input_path,"_p", nmr, ".txt")
    process_instance_Demanda = string(input_path,"_d", nmr, ".txt")
    process_instance_Ground = string(input_path,"_g", nmr, ".txt")
    cx,px,dx,gx = readdlm(process_instance_Custo, Int),readdlm(process_instance_Porto, Int),readdlm(process_instance_Demanda, Int),readdlm(process_instance_Ground, Int)
    return(cx,px,dx,gx)
end


using DelimitedFiles
function readInstance(nmr)    
    px =0
    dx=0
    cx=0
    gx=0
    nx=0
    Rx=0
    #instance = string("C:/Users/Pial/Dropbox/Mestrado Bruno Mendonça/Dados/Instancias/Instances.txt")
    instance = string("C:/Users/User/Instancias/Instances.txt")
    file_inst = readdlm(instance, Int)
    f_inst = reshape(file_inst,(:,4))
    #print(f_inst)
    #println(" ")
    for i in 1:size(f_inst)[1]
        #print(i)
        if (f_inst[i,1] == nmr)# || (f_inst[i,2] == R && f_inst[i,3] == n && f_inst[i,4] == cont)
            print(" ",f_inst[i,:])
            println(" .")
            cx,px,dx,gx = readProcessing(f_inst[i,2], f_inst[i,3], f_inst[i,4], f_inst[i,1])
            #println(px)
            nx = f_inst[i,3]
            Rx = f_inst[i,2]
        end
    end
    
    println(size(cx))
    return(reshape(cx,(nx,nx,Rx)), px, dx ,nmr, nx, Rx, reshape(gx,(nx,nx)))  
end

using LinearAlgebra


function groundGenerator(n)

    g =Matrix{Int}(reshape(zeros(n*n),(n,n)))
    for i in 1:n g[i,i]=1 end
    
    for j in 2:n
        g[1,j]=rand(0:1)
    end
    g[1,:]
    for i in 2:n
        for j in 1:n
            if i<j
                t=1
                k=0
                l=0
                for r in 1:i
                    if r< i
                        #println(" r=", r,"  g[",r,",",i,"]= ",g[r,i])
                        if g[r,i]==1
                            k=1
                            t=r
                        end
                        if g[r,j]==1
                            l=1
                            #t=r
                        end
                    end
                end

            
            
                g[i,j] =k*(g[t,j]) + (1-k)*(1-l)*rand(0:1)
                #println(" g[",i,",",j,"]= k*g[",t,",",j,"] + (1-k)*rand(0:1)")
                #rintln(". ",g[i,j]," = ",k,"*",g[t,j]," + (",1-k,")*rand(0:1)")
            end
            if i>j
                g[i,j]=g[j,i]
            end
        end
    end
    gr=Vector(Vector())
    sl=Vector{Int}()
    sle =Vector{Int}()
    #push!(sl,1)
    
    for i in 1:n
        for j in 1:n
            if (g[i,:]== g[j,:]) && !(j in sle) && (i<=j) 
                push!(sl,j)
                push!(sle,j)
                
            end
        end
        if !(sl in gr) && (sl!=[])
            push!(gr,sl)
        end
        sl=[]
    end

    return g
    #gg.append(g)
    #if R >1 
    #    gg.append(groundGenerator(n,R-1,g))
    #else
    #    return (gg)
    #end
    
end


function g_generator(n,R)
    gg=[]
    g = groundGenerator(n)
    push!(gg,g)
    if R>1
        for i in 1:R-1
            gx = groundGenerator(n)
            while sum(g.*gx)<=n
                gx = groundGenerator(n)
            end
            push!(gg,gx)
        end
    end
    gx =[]
    for j in gg
        for i in j
            push!(gx,i)
        end        
    end
    return gx
end


function Cennectionxy(nmr,h,cc, p, n, R, g)
    #cc, p, d, nmr, n, R, g = readInstance(nmr)
    x = Array{Int}(reshape(zeros(n*n*R),(n,n,R)))
    y = Array{Int}(reshape(zeros(n*n*R),(n,n,R)))
    w = Array{Int}(reshape(zeros(n*n*R),(n,n,R)))
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
            if h[j]==1 && i!=j && h[i]==0 && g[i,j]==1 &&  sum(w[i, :,:])+sum(w[:, i,:]) ==0 && sum(x[j,i,:])<1
                #x[i,j,r]=1
                if r==1
                    x[j,i,r]=p[i]*p[j]
                else    
                    x[j,i,r]=1
                end
                #w[i,j,r]=1
                w[j,i,r]=x[j,i,r]
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

function ShowResults(obj_val,x_val,y_val,h_val,w_val,n,R)
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
    println("Decision variable H:")
    for i in 1:n
        for j in 1:n

                if Int(round(w_val[i,j]))==1
                    println(i,"   ",j,"     H[",i,",",j,"] = ",Int(round(w_val[i,j])))
                end
            
        end
    end
    
    println("The model was run!")
    
        
end


function Cluster_Ground(g)
    gr=Cluster(g)

    return size(gr)[1]
    
end
function Cluster(g)
    n=size(g)[1]
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

    return gr
    
end

using Graphs
using Colors
function Graph_Plot(g,x,y,h,n)
    fonte1=""
    fonte2=""
    g=copy(g)
    if size(x)[1]>0 && size(y)[1]>0
        g=Array{Int}(reshape(zeros(n*n),(n,n)))
        for i in 1:R
            g =x[:,:,i]+y[:,:,i]+g
        end
    else
        g=g;
    end
    for i in 1:n
        for j in 1:n

            g[i,j]=round(g[i,j])

        end 

    end

    markercolor =[]
    node_weights =[]
    nodeshape=[]
    for i in 1:n
        push!(node_weights,2)
        if round(h[i])==0
            push!(markercolor,colorant"white")
            push!(nodeshape,:circle)
        else
            push!(markercolor,colorant"lightgreen")
            push!(nodeshape,:rect)
        end
    
    end
    
    node_weights =[]
    for i in 1:n
        push!(node_weights,10)
    end
    
    names = []
    if log(n)/log(10)>=2 && log(n)/log(10)<3
        fonte1="00"
        fonte2="0"
    end
    if log(n)/log(10)<2 && log(n)/log(10)>=1
        fonte1="0"
        fonte2=""
    end
    names
    for i in 1:n 
        if i<10
            push!(names,string(fonte1,i))
        end
        if i<100 && i>=10
            push!(names,string(fonte2,i))
        end
    
    end
    
    graphplot(g,markersize = 0.2, markercolor = markercolor,
    node_weights = node_weights, names = names, nodeshape=nodeshape,
     fontsize = 8,curves=false)
    
end

