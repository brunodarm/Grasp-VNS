using LinearAlgebra


function groundGenerator(n)
    

    g =reshape(zeros(n*n),(n,n))
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
    return (reshape(g,(n,n)))
end

#groundGenerator(25)