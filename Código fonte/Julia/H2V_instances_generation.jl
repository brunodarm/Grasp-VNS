include("H2V_processing_generator.jl")
include("H2V_read_processing.jl")


R =[2,3,4]
n =[10, 25, 50, 100, 250, 500]
qtd = 10 
println("The generator begins...")

instance = string("C:/Users/Pial/Dropbox/Mestrado Bruno Mendonça/Dados/Instancias/Instances.txt")
file = open(instance, "w")
cont = 1
for r in R
    for i in n
        for j in 1:qtd
            println(file, cont, "   ",r, "   ",i, "   ",j)
            processing_generator(r,i,j,cont)
            #p = read_processing(r,i,j,cont)   
            global cont+=1    
            if cont<=2
                #println(p)
            end
        end
    end
end



close(file)
println(cont-1, " generated test instances!")

