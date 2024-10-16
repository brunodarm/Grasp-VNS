using DelimitedFiles
function read_processing(R,n,cont,nmr)
    #input_path =("/home/bruno/Dropbox/Pesquisa/Problemas/Green parallel machines/Source/Instances generator/Testbed_tardiness/")
    input_path =("C:/Users/Pial/Dropbox/Mestrado Bruno Mendonça/Dados/Instancias/H2V")
    process_instance = string(input_path,"r",R,"_n",n,"_inst_", cont, ".txt")
    p = readdlm(process_instance, Int)
    return(p)
end
