function arvoreGeradora(c,n)
    c=c
    stx=[]
    n=Int(n)
    for i in 1:n
        stt,st,g=arvoreMin(c,i)
        for j in st
            push!(stx,j)
        end
    end
    return stx    
end 

function arvoreMin(c=[],a=1)
    n = size(c)[1]
    N =Array(1:n)

    grafo=Array{Int}(reshape(zeros(n*n),(n,n)))
    #seq = seqc(N)
    C =[]
    append!(C,a)
    Cx = N
    filter!(e->e≠C[1],Cx)
    LST=0
    ST=[]
    j=0
    i=1
    while (size(Cx)[1]>0 )
        l = argmin(c[C,Cx])
        a=getindex(c[C,Cx],l)
        b1=findall(e->e==a,c[:,:])
        [b1]
        for i in 1:size(b1)[1]
            if (b1[i][1] in C) && (b1[i][2] in Cx) 
                k = b1[i][1]
                #println(" k = ",k)
                    
                j = b1[i][2]
                #println(" j = ",j)
                append!(C,j)
                filter!(e->e≠j,Cx)
                LST=LST +c[k,j]
                grafo[k,j]=1.0
                grafo[j,k]=1.0
                ST=append!(ST,(k,j))
                i=i+1
            end

        end
    end
    
    return LST,ST,grafo#,seq
end


function seqc(c=[],a=1)
    kx=[]
    LST,ST,grafo = arvoreMin(c,a)
    n = size(c)[1]
    N =Array(1:n)
    for i in 1:n
        #append!(kx,size(findall(e-> e==i,st))[1])
        if size(findall(e-> e==i,ST))[1]> n/10 && sum(kx)<(n)/2
            append!(kx,1)
        else
            append!(kx,0)
        end
        
    end
    return LST,ST,grafo, Array{Int}(kx)
end

function hub_selection_01(p,ST,grd,c,px,Alfa=0,a=1)
    px=px
    n = size(c)[1]
    #print(p," ,")
    #LST,ST,g = arvoreMin(c,a)
    st =copy(ST)
    h=zeros(n)
    gr=Cluster(grd)
    fq = FreqNode(st,n)#frequencia em que os nós aparecem na arvore
    hub_cand=reverse(sortperm(fq[2,:]))#deixa em ordem decrescente de frequencia como retorno o numoro do nó
    
    for i in gr
        hq=[]
        if size(i)[1]==1
            h[i[1]]=1
        end
        h_c=copy(hub_cand)
        filter!(e->e in i,h_c)
        for k in h_c#i#hub_cand[i]
            if px[k]==1
                push!(hq,k)
            end
        end
        # quando alfa=0 entao será ecolha gulosa
        if rand(0:100)<=Alfa 
            h[rand(hq)]=1
        else
            h[hq[1]]=1
        end

    end
    m = Cluster_Ground(grd)  

    for k in 1:size(hub_cand)[1]
        if (fq[2,hub_cand[k]]>0 && sum(h)<= p) #|| px*h<=m
            if rand(0:100)<=Alfa 
                h[hub_cand[rand(k:size(hub_cand)[1])]]=1
            else
                h[hub_cand[k]]=1
            end
            

        end
    end


    return Array{Int}(h)
end


function FreqNode(ST,n)#Conta a frequencia de nós apos MSP
    n=n
    y=copy(sort(ST))
    #n=size(ST)[1]
    freq=Array{Int}(reshape(zeros(n*2),(n,2)))
    for i in 1:n
        freq[i,1]=i
        soma=0
        for j in y 
            if i==j
                soma=soma+1
            end
        end
        freq[i,2]=soma
    end
    return transpose(freq)
end


function ConectAcess(h,c,g,grd,px)
    n=size(c)[1]
    C=copy(c)
    gx=copy(g)
    #Conectando hubs
    gx=hubAcess(h,gx,grd,px)

    #sendo g grafo de elemento binario
    for j in 1:n
        for k in 1:n
            if gx[j,k]>=1
                gx[j,k]=1
            end
            if j==k
                gx[j,k]=0
            end
        end
    end
    #gx = X_conect(gx,h,n)
    
    
    #Custo=sum(C.*(gx))
    return gx#,Custo
    
end

function X_conect(g,h,n)
    gx=g
    for i in 1:n
        flg_1=0
        for j in 1:n
            if h[j]*gx[i,j]==1 && h[i]==0
                flg_1=1
                for k in (i+1):n 
                    if gx[k,j]==1 && h[k]==0
                        gx[i,k]=1
                        gx[k,i]=1
                    end
                end
            end
        end
    end
    return gx
end

function hubAcess(h,g,grd,px)
    gx=copy(g)
    n=size(h)[1]
    if size(px)[1]==0
        px=Array{Int}(zeros(n))
    
        if size(grd)[1]==0
            grd=Array{Int}(reshape(zeros(n*n),(n,n)))
            for i in 1:n
                px[i]=1
                for j in 1:n 
                    grd[i,j]=1
                end 
                
            end
        end
    end
    #Conectando hubs
    for i in 1:n
        if h[i]==1
            gx[i,:]=gx[i,:]+h
        end
    end
     #sendo g grafo de elemento binario
    for j in 1:n
        for k in 1:n
            if round(gx[j,k])>=1 && (grd[j,k]==1 || px[j]+px[k]==2)
                gx[j,k]=1
            else
                gx[j,k]=0
            end
            if j==k
                gx[j,k]=0
            end
        end
    end
    return gx
end

function nodeAcess(h,c,grd,Alfa)
    gr=[]
    dentro=[]
    n=size(c)[1]
    for i in 1:n
         if h[i]==1
              push!(gr,[i])
              push!(dentro,i)
         end
    end
    fora=[]
    for i in 1:n
         if !(i in dentro)
              #print(i,",")
              push!(fora,i)
         end
    end
    g,custo=nodeAcess_1(c,grd,gr,fora)
    #Graph_Plot(g,[],[],h,n)
    fora=[]
    select =0
    for k in 1:size(gr)[1]
         for i in 1:size(gr[k])[1]
              if i>1 && size(gr[k])[1]>=i
                   if rand(0:100)<=Alfa  
                        select=gr[k][rand(2:size(gr[k])[1])]
                   else
                        select=gr[k][i]
                   end
                   if !(select in fora)
                        filter!(e->e≠select,gr[k])
                        push!(fora,select)
                        g,custo=nodeAcess_1(c,grd,gr,select)
                   end
              end
         end
    end
    
    gx=conectaGraph(gr,n)
    gx=gx
    return gx
end

function nodeAcess_1(c,grd,gr=[],fora=[])
    fora=fora
    n=size(c)[1]
    g=Matrix{Int}(reshape(zeros(n*n),(n,n)))
    #sum(c.*g)
    menor=9999
    jx=0;ix=0
    for i in fora
         menor=9999
         for j in 1:size(gr)[1]
              push!(gr[j],i)
              g=conectaGraph(gr,n)
              custo=sum(c.*g)
              if custo<menor && grd[i,gr[j][1]]==1
                   menor=custo
                   jx=j
                   ix=i
              end
              filter!(e->e≠i,gr[j])
              #println(" Custo=",custo," gr=",gr,", i=",i)
         end
         push!(gr[jx],ix)
         #filter!(e->e≠ix,fora)
         #push!(dentro,ix)
         #println("gr=",gr)
         g=conectaGraph(gr,n)
    end
    return g,menor     
end




function H2V_graph(cx,px,gx,nx)
    n=nx
    grd=gx
    gx = hubAcess(transpose(px),gx,[],[])
    cy=cx
    Rx=size(cx)[3]
    
    for r in 1:Rx
         if r==1
              cy[:,:,r]=cx[:,:,r].*px
         else
              cy[:,:,r]=cx[:,:,r].*gx
         end
    end

    for i in 1:nx
         for j in 1:nx
              for r in 1:Rx
                   if cy[i,j,r]==0 #|| grd[i,j]==0
                        cy[i,j,r]=999
                        
                   end
              end
         end
    end

    cx=cy
    cz=Array{Int}(reshape(zeros(n*n),(n,n)))
    c=Array{Int}(reshape(zeros(n*n),(n,n)))
    r=Array{Int}(reshape(zeros(n*n),(n,n)))
    for i in 1:n
         for j in  1:n
              if i!=j
                if grd[i,j]==1
                   c[i,j]=minimum(cx[i,j,1:Rx])
                   r[i,j]=argmin(cx[i,j,1:Rx])
                else
                    c[i,j]=cx[i,j,1]
                    r[i,j]=1
                end
                if round(c[i,j])==0 || c[i,j]>=999
                    cz[i,j]=0
                else
                    cz[i,j]=1
                end
                
              end
         end

    end
    #a=1
    #stt,st,g=arvoreMin(c,a)
    #p, g, custo, h = Varredura(c)
    return c,r,cz
end

function somaCol_Lin(g,ori)
    gx=copy(g)
    n=size(g)[1]
    col=[]
    for i in 1:n
        if ori=="linha"
            push!(col,sum(gx[i,:]))
        end
        if ori=="coluna"
            push!(col,sum(gx[i,:]))
        end

    end
    return col
end

function Varredura(c,grd,px,Alfa)
    stt,st,g=arvoreMin(c,1)
    n=size(c)[1]
    if size(px)[1]==0
        px=Array{Int}(zeros(n))
        
        if size(grd)[1]==0
            grd=Array{Int}(reshape(zeros(n*n),(n,n)))
            for i in 1:n
                px[i]=1
                for j in 1:n 
                    grd[i,j]=1
                end 
                
            end
        end
    end

    
    p=0
    g_o=g
    h=[]
    menor=99999999999999
    for j in 1:n
         g,custo,hx=Revisado(c,j,grd,px,Alfa)
         if custo<menor
              menor=custo
              g_o =g
              p= (j*100-mod(j*100,n))/n
              h=hx
         end
         println(" p=",p,"% custo=",menor," | ")

    end
    return p, g_o, menor,h

end


function nodeAcess_2(h,c,grd)
    n=size(c)[1]
    gx=Array{Int}(reshape(zeros(n*n),(n,n)))
    h=h; ix=0;jx=0;
    menor=99999
    #clus=Matrix{Int}(reshape(zeros(n*n),(n,n)))
    for i in 1:n
        ix=0;jx=0;
        menor=99999
        for j in 1:n
            if h[j]==1 && j!=i && h[i]==0 && c[i,j]>0 
                if c[i,j]<=menor && grd[i,j]==1
                    menor=c[i,j]    
                    ix =i
                    jx=j
                    #gx[ix,jx]=1
                    #gx[jx,ix]=1               
                end       
            end
        end
        if ix>0 && jx>0
            gx[ix,jx]=1
            gx[jx,ix]=1
        end
    end  

    return gx
end

function Revisado(c,nh,grd,px,Alfa)
    p=nh
    stt,st,g=arvoreMin(c,1)
    h = hub_selection(p,st,grd,c,px,Alfa,1)
    


    g =nodeAcess(h,c,grd,Alfa)

    g,custo =ConectAcess(h,c,g,grd,px)

    return g, custo,h
end

function C_cluster(g)
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
    saida = Matrix{Int}(reshape(zeros(n*n),(n,n)))
    k=1
    for i in gr
        for j in i  
            saida[k,j]=1
        end 
        k=k+1
    end

    return saida,gr,size(gr)[1]
    
end
function Heuristic_03(cx,c,nh,grd,px,st,Alfa)
    cx=cx
    px=px
    n=size(cx)[1]
    nh=nh
    #c,r,cz =H2V_graph(cx,px,grd,n)
    #p=nh
    #st=arvoreGeradora(c)
    h = hub_selection(nh,st,grd,c,px,Alfa,1)
    #g =nodeAcess(h,c,grd,Alfa)
    #g =ConectAcess(h,c,g,grd,px)
    #custo,x,y,H =custoso(g,cx,px,grd,h)
    g=Array{Int}(reshape(zeros(n*n),(n,n)))

    hg,g3,h_gr = X_acess(h,c,grd,px)
    g1=conectaGraph(h_gr,nx)
    g2=hubAcess(h,g,grd,px)
    g=g1+g2
    custo = sum(c.*g)
    return g, custo,h
    
end


function Heuristica_GRASP(cx,c,grd,px,st,Alfa)
    px=px
    grd=grd
    cx=cx
    st=st
    c=c
    n=size(cx)[1]
    g=Array{Int}(reshape(zeros(n*n),(n,n)))
    #stt,st,g=arvoreMin(c,1)
    n=size(c)[1]
    if size(px)[1]==0
        px=Array{Int}(zeros(n))
        
        if size(grd)[1]==0
            grd=Array{Int}(reshape(zeros(n*n),(n,n)))
            for i in 1:n
                px[i]=1
                for j in 1:n 
                    grd[i,j]=1
                end 
                
            end
        end
    end

    m = Cluster_Ground(grd)
    ns=0
    upBound=round(n/2)+1
    g_o=Array{Int}(reshape(zeros(n*n),(n,n)))
    h=Array{Int}(reshape(zeros(n),(1,n)))
    menor=99999999999999
    t1=0
    z=0
    if m<upBound/2
        m=Int(round(upBound/2)-1)
    end
    for j in m:upBound
        if t1>=600
            break
        end
        if  z>=(upBound-m)/3 +1
            break
        end
        
        t2=@elapsed g,custo,hx=Heuristic_03(cx,c,j,grd,px,st,Alfa)
        if custo<menor
                menor=custo
                g_o =g
                ns= (j*100-mod(j*100,n))/n
                h=copy(hx)
                z=z+1
        else
            z=0
        end
        t1=t2+t1
        #println(" p=",p,"% custo=",menor," | ")
        
        
        
        
        
        
    end
    return ns, g_o, menor,h
    
end

function custoso(g,cx,px,grd,h)
    n=size(cx)[1]
    Rx=size(cx)[3]
    cx=cx
    c,r,cz =H2V_graph(cx,px,grd,n)
    H=Array(reshape(zeros(n*n),(n,n)))
    x=Array(reshape(zeros(n*n*Rx),(n,n,Rx)))
    y=Array(reshape(zeros(n*n*Rx),(n,n,Rx)))

    for i in 1:n
        for j in 1:n
            if j!=i && r[i,j]!=0
                if h[i]==0 && g[i,j]==1
                    x[i,j,r[i,j]]=1
                    x[j,i,r[i,j]]=1  
                    if h[j]==1 
                       H[j,i]=sum(x[j,i,:]) 
                    end

                end
                if h[i]==1 && h[j]==1 && g[i,j]==1
                    y[i,j,r[i,j]]=1
                    y[j,i,r[i,j]]=1          
                end
            end
        end
    end
    cx
    custo=custo_Grafo(c,x,y,Rx,n)#sum(cx.*(x+y))

    return custo,x,y,H
end


function conectaGraph(gr,n)
    n=n    
    g=Matrix{Int}(reshape(zeros(n*n),(n,n)))
    for k in gr
         g[k,k]=Unado(size(k)[1])
    end
    return g
end




function Unado(n)
    g=Matrix{Int}(reshape(zeros(n*n),(n,n)))
    for i in 1:n 
         for j in 1:n
              if i!=j 
                   g[i,j]=1
              else
                   g[i,j]=0
              end
         end
    end
    return g
end


function X_acess(h,cx,gx,px)
    n=size(h)[1]
    g=Matrix{Int}(reshape(zeros(n*n),(n,n)))
    cx=cx
    gr=Cluster(gx)
    hg=[]
    ht=[]
    h_gr=[]
    lista=Array{Int}(1:n)
    for i in 1:n
        if h[i]==1
            push!(ht,i)
            for k in gr
                if i in k
                    filter!(e->e!=i,k)
                end
            end
        end
    end
    for i in ht 
        push!(h_gr,[i])
    end
    filter!(e->!(e in ht),lista)
    ww=0
    while size(lista)[1]>0 && ww<150
        
        ww=ww+1
        for i in ht
            iix=[]
            flag=1
            #print("i =",i," lista = " ,lista," tamanho de ",size(lista)[1])
            for k in gr
                #print(", gr =",k)
                if size(lista)[1]>0  && size(k)[1]>=1 && sum(gx[i,k])>=1
                    #ca =#testeConexao(cx,i,lis,gx,px)
                    #filter!(e->e!=i,k)
                    delta =argmin(cx[i,k])#ca)
                    ii =k[delta]#[1]]
                    #push!(k,i)
                    #r=delta[2]

                
                    g[i,ii]=1
                    g[ii,i]=1
                    #push!(hg,[i,ii,r])
                    filter!(e->e!=ii,lista)
                    filter!(e->e!=ii,k)
                    for j in 1:size(ht)[1] 
                        if i in h_gr[j]
                            push!(h_gr[j],ii) 
                        end
                    end
                
                end
                if size(lista)[1]==0
                    break
                end
            end
            #println("")
            
            

        end
    end
    #println(h_gr)
    return hg,g,h_gr
end
function Y_acess()
    
end

function testeConexao(cx,a,bb,gx,px,h)
    R=size(cx)[3]
    #bx=bb
    cx=cx
    ca=cx
    saida=true
    for b in bb
        if a!=b
            if px[a]*px[b]==0
                if gx[a,b]==1     
                    #filter!(e->e>1,r)
                    ca[a,b,1]=9999
                else
                    ca[a,b,:]=(ca[a,b,:]./ca[a,b,:]).*9999
                    saida=false
                end
            else
                if gx[a,b]==0
                    #filter!(e->e<=1,r)

                    ca[a,b,:]=(ca[a,b,:]./ca[a,b,:]).*9999
                    saida=false
                end
                if h[a]*h[b]==0
                    ca[a,b,1]=9999
                end
            end
        else
            ca[a,b,:]=(ca[a,b,:]./ca[a,b,:]).*9999
            saida=false
        end
    end
    return ca[a,bb,:]
    
end

function hub_selection(p,ST,grd,c,px,Alfa=0,a=1)
    px=px
    n = size(c)[1]

    #print(p," ,")
    #LST,ST,g = arvoreMin(c,a)
    st =copy(ST)
    h=zeros(n)
    gr=Cluster(grd)
    fq = FreqNode(st,n)#frequencia em que os nós aparecem na arvore
    hub_cand=GraspSequence(st,n,Alfa)#deixa em ordem decrescente de frequencia como retorno o numoro do nó
    
    for i in gr
        hq=[]
        if size(i)[1]==1
            h[i[1]]=1
        end
        h_c=copy(hub_cand)
        filter!(e->e in i,h_c)
        for k in h_c#i#hub_cand[i]
            if px[k]==1
                push!(hq,k)
            end
        end
        # quando alfa=0 entao será ecolha gulosa
        h[hq[1]]=1
        #if rand(0:100)<=Alfa 
        #    h[rand(hq)]=1
        #else
        #    h[hq[1]]=1
        #end

    end

    for k in 1:size(hub_cand)[1]
        if  sum(h)<=n/2 -1
            if (fq[2,hub_cand[k]]>0 && sum(h)<= p ) #|| px*h<=m
                h[hub_cand[k]]=1
                #if rand(0:100)<=Alfa 
                #    h[hub_cand[rand(k:size(hub_cand)[1])]]=1
                #else
                #    h[hub_cand[k]]=1
                #end


            end
        end
    end


    return Array{Int}(h)

end


function FreqNode_02(ST)#Conta a frequencia de nós apos MSP
    y=copy(sort(ST))
    zeta=[]
    while size(y)[1]>0
        a_i=y[1]
        selecao(a_i,y)
        push!(zeta,a_i)
    end
    y=copy(sort(ST))
    n=size(zeta)[1]
    freq=Array{Int}(reshape(zeros(n*2),(n,2)))
    for i in 1:n
        freq[i,1]=zeta[i]
        soma=0
        for j in y 
            if zeta[i]==j
                soma=soma+1
            end
        end
        freq[i,2]=soma
    end
    return transpose(freq)
end


function custo_Grafo(c,x_val,y_val,R,n)
    c=c
    R=R
    n=n
    g=Array{Int}(reshape(zeros(n*n),(n,n)))
    for i in 1:R
        g =x_val[:,:,i]+y_val[:,:,i]+g
    end
    for i in 1:n
        for j in 1:n
            g[i,j]=round(g[i,j])
        end 
    end

    return sum(c.*g)
end

