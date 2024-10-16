include("FunctionsH2V.jl")



R =[2,3,4]
n =[10, 25, 50, 100, 250, 500]
qtd = 10 
println("The generator begins...")

instance = string("C:/Users/User/Instancias/instances.txt")
file = open(instance, "w")
cont = 1
for r in R
    for i in n
        for j in 1:qtd
            println(file, cont, "   ",r, "   ",i, "   ",j)
            processingGenerator(r,i,j,cont)
            #p = readProcessing(r,i,j,cont)   
            global cont+=1    
            if cont<=2
                #println(p)
            end
        end
    end
end



close(file)
println(cont-1, " generated test instances!")

