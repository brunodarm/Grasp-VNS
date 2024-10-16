
using Plots, GraphRecipes
#import Pkg; Pkg.add("LightGraphs")
#import Pkg; Pkg.add("NetworkLayout")
g = wheel_graph(10)
graphplot(g, curves=false)

using Karnak
using Graphs
using NetworkLayout
g = barabasi_albert(10, 1)
@drawsvg begin
    background("black")
    sethue("white")
    drawgraph(g, layout=stress, vertexlabels = 1:nv(g))
end


using Graphs
using Karnak
using NetworkLayout
using Colors
g = barabasi_albert(60, 1)

nv(g)

@drawsvg begin
    background("black")
    sethue("grey40")
    fontsize(8)
    drawgraph(g, 
        layout=stress, 
        vertexlabels = 1:nv(g),
        vertexfillcolors = 
            [RGB(rand(3)/2...) 
               for i in 1:nv(g)]
    )
end 600 400


using LightGraphs
nodelabel = 1:size(g)[1]
gplot(g, nodelabel=nodelabel)