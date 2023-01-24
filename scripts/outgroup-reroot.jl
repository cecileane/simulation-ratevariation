#=
to use the functions in this file, read it within julia with:
include("path/to/outgroup-reroot.jl");
but replace 'path/to/' by the actual path to this file on your machine.
main function defined in this file:

rootgenetrees
rootgenetree!
outgrouplabels

and helper functions:

clustersize
clustersize!
=#

try # install package PhyloNetworks, if not already installed, and load it
    using PhyloNetworks
catch
    @info "need PhyloNetworks package: will add & use it"
    import Pkg; Pkg.add("PhyloNetworks"); using PhyloNetworks
end

"""
    rootgenetrees(genetreefile, speciesnet::HybridNetwork; outgroup)

Read in the input phylogenies in `genetreefile`, root them all with the
species network (and possible outgroups), then write the rooted
phylogenies to an output file. See help on `rootgenetrees!` for how this
re-rooting is done, e.g. when the input phylogenies have a subset of the
taxa in the species network (possibly excluding outgroups).

The output file name is the same as the gene tree file name,
except for an extra "_rooted" before the extension.
For example: if the input phylogenies are in "path/to/genetrees.tre",
then the output file name will be "path/to/genetrees_rooted.tre".
If the input phylogenies are in "genetrees", then the output file name
will be "genetrees_rooted".

Assumption: input phylogenies are trees. In fact, the function would work
sensibly with input phylogenies that are reticulate (non-tree) networks.
But it hasn't been tested.

# example

In the example below, we read a species network (from a file "astralfile")
that is unrooted. We first root it with a known outgroup clade.
Then we read the *rooted* species network and use it to root all phylogenies
in another file (variable "genetreefile").

```julia
julia> astralfile = "path/to/astral.tre"; # contains 102 trees, say

julia> astraltree_unrooted = readTopology(astralfile);

julia> rootgenetrees(astralfile, astraltree_unrooted; outgroup=["ornAna1","hg19"]);
[ Info: rooted input phylogenies written to file path/to/astral_rooted.tre

julia> astraltree_rooted = readMultiTopology("path/to/astral_rooted.tre")[end];

julia> genetreefile = "path/to/besttrees.tre";

julia> rootgenetrees(genetreefile, astraltree_rooted);
[ Info: rooted input phylogenies written to file path/to/besttrees_rooted.tre

```
"""
function rootgenetrees(genetreefile, snet::HybridNetwork;
        outgroup=String[]::AbstractVector)
    gtrees = readMultiTopology(genetreefile)
    match(r"(?<ext>\.[^.]+)$", genetreefile)
    if occursin(".", basename(genetreefile))
        outputfile = replace(genetreefile, r"(\.[^.]+)$" => s"_rooted\1")
    else
        outputfile = genetreefile * "_rooted"
    end
    directEdges!(snet)
    rootedgtrees = map(g -> rootgenetree!(g,snet;outgroup=outgroup,verbose=false),
                       gtrees)
    writeMultiTopology(rootedgtrees, outputfile)
    @info "rooted input phylogenies written to file $outputfile"
    return gtrees
end

"""
    rootgenetree!(genetree::HybridNetwork, speciesnet::HybridNetwork;
        outgroup=String[], directedges=true::Bool)

Re-root one phylogeny `genetree` using the outgroup in `speciesnet`.
The species network is assumed to be rooted, and with edges correctly directed
if `directedges` is set to false.

If `outgroup` is specified, it should form a clade stemming from the root
in the species network, but this is not checked.
If `outgroup` is not specified, it is taken to be the smallest daughter
"clade" of the root in the species network. In a reticulate network, a "clade"
is defined here as a hardwired cluster.

If the gene tree has taxa not present in the species network, an error is thrown.

There are 2 sources of complications:
the gene tree may be missing some taxa (like the outgroup),
and the outgroup may consist of clade of multiple taxa in the species network,
but the gene tree may fail to have a split separating the ingroup from the outgroup.

1. If the gene tree is missing taxa, then the function prunes the species to
   the taxa present in the gene tree. The outgroup is recalculated as the set
   of outgroup taxa present in the gene tree. If this is empty, then the
   outgroup is recalculated as smallest daughter clade of the root in the
   pruned species network.
2. If the outgroup has more than 1 taxon and if the gene tree fails to have
   a spit outgroup | ingroup, then the outgroup is re-set, to be a first
   taxon in the outgroup list.

# note

An alternative criterion would be to choose the rooting that minimizes
the MDC criterion: read Yu Warnow & Nakhleh (2011)

# example

```julia
julia> speciesnet = readTopology("(anoCar2,pytMol0,((chrPic0,(allMis0,((galGal3,taeGut1):2.02)#H9:10.0::0.6):0.61):0.6,((hg19,ornAna1):1.54,#H9:0.0::0.4):0.22):2.32);");

julia> rootonedge!(speciesnet, 13) # to root along the edge separating mammals (hg19,ornAna1) from reptiles

julia> speciesnet = readTopology("((hg19,ornAna1):0.77,(#H9:0.0::0.4,((chrPic0,(allMis0,((galGal3,taeGut1):2.02)#H9:10.0::0.6):0.61):0.6,(anoCar2,pytMol0):2.32):0.22):0.77);");

julia> genetree = readTopology("(anoCar2,(pytMol0,((chrPic0,(allMis0,(galGal3,taeGut1):12.02):0.61):0.6,(hg19,ornAna1):1.76):2.32));");

julia> rootgenetree!(genetree, speciesnet; verbose=false);

julia> writeTopology(genetree)
"((hg19,ornAna1):0.88,((chrPic0,(allMis0,(galGal3,taeGut1):12.02):0.61):0.6,(pytMol0,anoCar2):2.32):0.88);"

julia> genetree = readTopology("((hg19,chrPic0),((ornAna1,(allMis0,galGal3)),(anoCar2,pytMol0)));");

julia> rootgenetree!(genetree, speciesnet);
[ Info: pruned 1 taxa
[ Info: outgroups: hg19, ornAna1
[ Info: outgroups do not form a clade in the gene tree. will use hg19 only

julia> writeTopology(genetree)
"(hg19,(chrPic0,((ornAna1,(allMis0,galGal3)),(anoCar2,pytMol0))));"
````
"""
function rootgenetree!(gtree::HybridNetwork, snet::HybridNetwork;
        outgroup=String[]::AbstractVector, directedges=true::Bool,
        verbose=true::Bool)
    gtree.numHybrids == 0 || @warn "not tested on reticulate (non-tree) input networks"
    directedges && directEdges!(snet)

    # prune the reference network (snet) to have the
    # same of taxa as the input gene tree (gtree)
    gtaxa = tipLabels(gtree)
    staxa = tipLabels(snet)
    issubset(gtaxa, staxa) || error("there are extra taxa in the input tree")
    extratips = setdiff(staxa, gtaxa)
    if !isempty(extratips)
        snet = deepcopy(snet)
        for taxon in extratips
            deleteleaf!(snet, taxon)
        end
        outgroup = setdiff(outgroup, extratips)
        verbose && @info "pruned $(length(extratips)) taxa"
    end

    # build the outgroup set, if empty
    if isempty(outgroup)
        outgroup = outgrouplabels(snet, false)
    end
    verbose && @info "outgroups: $(join(outgroup, ", "))"

    # find the edge, in the input tree, that corresponds to the outgroup
    if length(outgroup)>1 # multiple outgroup taxa: need to find the internal edge
        taxa = setdiff(gtaxa, outgroup) # ingroup only, vector
        ningroup = length(taxa)
        append!(taxa,outgroup) # tacking outgroups at the end
        bipartition1 = vcat([1 for i in 1:ningroup],[0 for o in outgroup])
        bipartition2 = vcat([0 for i in 1:ningroup],[1 for o in outgroup])
        mat = hardwiredClusters(gtree, taxa)
        size(mat,2) == length(bipartition1) + 2 || error("bipartition of wrong size...")
        i = findfirst(
            i -> mat[i,2:end-1] == bipartition1 || mat[i,2:end-1] == bipartition2,
            1:size(mat,1))
        if !isnothing(i) # the outgroup forms a clade in the input tree
            edgenum = mat[i,1]
            rootonedge!(gtree, edgenum)
            return gtree
        end
        verbose && @info "outgroups do not form a clade in the input tree. will use $(outgroup[1]) only"
        # the outgroup does NOT form a clade: the first one will be used.
    end
    outgroup = outgroup[1] # node label, string
    rootatnode!(gtree, outgroup)
    return gtree
end

"""
    outgrouplabels(net::HybridNetwork, directedges=true::Bool;
                    subtreeroot = net.node[net.root],
                    ingroup_also = false::Bool)

Set of labels for tips in the outgroup, which is chosen to be the
descendants of the (first) edges of the root with the minimum number
of descendants.

# example
```julia
julia> net = readTopology("((hg19,ornAna1),(#H9:::0.4,((chrPic0,(allMis0,((galGal3,taeGut1))#H9:::0.6)),(anoCar2,pytMol0))));");

julia> outgrouplabels(net)
Set{String} with 2 elements:
  "ornAna1"
  "hg19"

julia> outgrouplabels(net; ingroup_also=true)
(["hg19", "ornAna1"], ["chrPic0", "allMis0", "galGal3", "taeGut1", "anoCar2", "pytMol0"])

julia> outgrouplabels(net; subtreeroot=net.node[16], ingroup_also=true) # outgroup subset of ingroup!
(["galGal3", "taeGut1"], ["chrPic0", "allMis0", "galGal3", "taeGut1", "anoCar2", "pytMol0"])

julia> outgrouplabels(net; subtreeroot=net.node[15], ingroup_also=true)
(["anoCar2", "pytMol0"], ["chrPic0", "allMis0", "galGal3", "taeGut1"])

julia> outgrouplabels(net; subtreeroot=net.node[9]) # this subtreeroot is a hybrid node
ERROR: the (subnetwork) root has a single daughter

julia> rootatnode!(net, "allMis0");

julia> outgrouplabels(net)
Set{String} with 1 element:
  "allMis0"

```
"""
function outgrouplabels(net::HybridNetwork, directedges=true::Bool;
                    subtreeroot::PhyloNetworks.Node = net.node[net.root],
                    ingroup_also::Bool = false)
    directedges && directEdges!(net)
    daughteredges = [e for e in subtreeroot.edge if PhyloNetworks.getParent(e)===subtreeroot]
    length(daughteredges) > 1 || error("the (subnetwork) root has a single daughter")
    if ingroup_also && length(daughteredges)>2
        @error("the (subnetwork) root has more than 2 daughters: no single ingroup. Turning off 'ingroup_only'.")
        ingroup_also = false
    end
    daughtersizes = [clustersize(ce) for ce in daughteredges]
    i = findmin(daughtersizes)[2] # first cluster of smallest size
    nodenum = PhyloNetworks.descendants(daughteredges[i])
    out_tiplabels = [n.name for n in net.node if n.number in nodenum]
    if ingroup_also
        i = 3-i # changes i: 2->1 or 1->2
        nodenum = PhyloNetworks.descendants(daughteredges[i])
        in_tiplabels = [n.name for n in net.node if n.number in nodenum]
        return (out_tiplabels, in_tiplabels)
    end
    return out_tiplabels
end

"""
    clustersize(edge::Edge)

Size of the edge's hardwired cluster (set of descendants).

`edge` should belong in a rooted network for which `isChild1` is up-to-date.
Run `directEdges!` beforehand on this network.
This is very important, otherwise one might enter an infinite loop,
and the function does not test for this.

# example

```julia
julia> net = readTopology("((hg19,ornAna1),(#H9:::0.4,((chrPic0,(allMis0,((galGal3,taeGut1))#H9:::0.6)),(anoCar2,pytMol0))));");

julia> clustersize(net.edge[3])
2

julia> clustersize(net.edge[17])
6

julia> # printEdges(net) # to see the list of edges with their parent & child nodes
```
"""
function clustersize(edge::PhyloNetworks.Edge)
    visited = Int[]
    ntips = clustersize!(visited, edge)
    return ntips
end
function clustersize!(visited::Vector{Int}, edge::PhyloNetworks.Edge)
    ntips = 0
    n = PhyloNetworks.getChild(edge)
    if n.hybrid # only need to check previous visits for hybrid nodes
        n.number in visited && return ntips
        push!(visited, n.number)
    end
    if n.leaf
        ntips += 1
    end
    for ce in n.edge
        if n == PhyloNetworks.getParent(ce)
            ntips += clustersize!(visited, ce)
        end
    end
    return ntips
end
