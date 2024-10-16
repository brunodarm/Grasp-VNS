using DelimitedFiles
function read_instance(R,n,cont,nmr)    
    px =0
    include("H2V_read_processing.jl")
    instance = string("C:/Users/Pial/Dropbox/Mestrado Bruno Mendonça/Dados/Instancias/Instances.txt")
    file_inst = readdlm(instance, Int)
    f_inst = reshape(file_inst,(:,4))
    #print(f_inst)
    #println(" ")
    for i in 1:size(f_inst)[1]
        #print(i)
        if f_inst[i,2] == R && f_inst[i,3] == n && f_inst[i,4] == cont
            print(" ",f_inst[i,:])
            println(" .")
            px = read_processing(R,n,cont)
            #println(px)
            
        else
            println(" ")
        end
    end

    println(size(px))
    return(reshape(px,(n,n,R)) )  
end