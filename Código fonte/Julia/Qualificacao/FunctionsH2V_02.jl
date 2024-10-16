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

function read_data(i)
    #input_path =("/home/bruno/Dropbox/Pesquisa/Problemas/Green parallel machines/Source/Instances generator/Testbed_tardiness/")
    endereco = "C:/Users/User/Instancias/H2V"
    input_path =(endereco)
    process = string(input_path,"_Results_",i, ".txt")
    #MIP_value,MIP_time,MIP_gap,GV_value,GV_time,GV_gap,ML_value,ML_time,ML_gap
    saida = readdlm(process, Int)
    return(saida)
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
    g
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
    p = Vector{Int}()
    for i in 1:n
        push!(p,rand(0:1))
        for j in 1:size(gr)[1]
            if i == gr[j][1]
                p[i]=1       
            end
        end
    end

    return (g,p)
end

