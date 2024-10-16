function processing_generator(R,n,cont,nmr)
    #output_path =("/home/bruno/Dropbox/Pesquisa/Problemas/Green parallel machines/Source/Instances generator/Testbed_tardiness/")
    output_path =("C:/Users/Pial/Dropbox/Mestrado Bruno Mendonça/Dados/Instancias/H2V")
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
    instance1 = string(output_path,"_p", nmr, ".txt")
    process1 = open(instance1, "w")
    for i in 1:n
        pp = rand(0:1)
        print(process1, pp, "   ")
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

end
    
