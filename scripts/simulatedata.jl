#= code to simulate:
1. networks using SiPhyNetwork (via RCall)
2. gene trees within a network, using PhyloCoalSimulations
3. rate variation: scale branch lengths in gene trees
4. sequences along gene trees, using seqgen

import the file to use its functions:

include("simulatedata.jl")
=#

using Random
using PhyloNetworks
using PhyloCoalSimulations
using RCall
using Distributions

# below: add with full path to executable if "seq-gen" is not on the path
host = chomp(read(`hostname`, String))
if occursin("Lauren", host)
    seqgen = "/Users/laurenfrankel/Desktop/simulation-ratevariation/bin/Seq-Gen-1.3.4/source/seq-gen"
elseif occursin("darwin", host)
    seqgen = "/u/l/e/lefrankel/bin/Seq-Gen-1.3.4/source/seq-gen"
elseif occursin("franklin", host)
    seqgen = "/u/l/e/lefrankel/bin/Seq-Gen-1.3.4/source/seq-gen"
elseif occursin("cecile", host)
    seqgen = "/Users/ane/bin/seq-gen"
else
    @warn "Need to specify Seq-Gen location"
end
ispath(seqgen) || error("can't find seqgen at $seqgen")

# R"install.packages('devtools')"
# R"require(devtools)";
# R"devtools::install_github('jjustison/SiPhyNetwork')";
# use 10 as the precision for the disribution of γ values
# setting genetic distance function using example values
R"require(ape, quietly=T)";
R"require(SiPhyNetwork, quietly=T)";
R"inheritance_fxn <- make.beta.draw(2, 2)";
R"gen_dist_fxn = make.stepwise(probs = c(1,0),distances = c(0.75,Inf))";
hybridproportions = [0.5,0.25,0.25]  # generative, degenerative, neutral
@rput hybridproportions

using Pkg
Pkg.status(["PhyloNetworks","RCall","PhyloCoalSimulations"]);
R"rversion <- R.Version()$version.string";
@rget(rversion);
R"spnversion <- sessionInfo(package='SiPhyNetwork')$otherPkgs$SiPhyNetwork$Version";
R"spnsha     <- sessionInfo(package='SiPhyNetwork')$otherPkgs$SiPhyNetwork$RemoteSha";
@rget(spnversion);
@rget(spnsha);
python_version = ""
try
    global python_version = readchomp(`python --version`)
catch _
    @warn "'python' not in your path"
end
print("""julia version: $VERSION
      $rversion
      SiPhyNetwork version: v$spnversion
      SiPhyNetwork commit sha: $spnsha
      Python version: $python_version
      """)

"""
    simulatenetwork(ntaxa::Int, nsims::Int)

calls `simulatenetwork(ntaxa)` many times, to simulate `nsims` networks
and return them as a tuple with vector of `HybridNetworks` and a vector of 
levels of the networks. The keyword arguments are passed, except for the seed. 
"""
function simulatenetwork(ntaxa::Int, nsims::Int; seed=nothing, kwargs...)
    nets = HybridNetwork[]
    levels = Int[]
    if isnothing(seed)
        seed = rand(1:10000)
        @info "seed for the generation of random seeds: $seed"
    end
    Random.seed!(seed)
    s = rand(1:100_000, nsims) # each s value will be used as a seed within R for 1 network
    for i in 1:nsims
        n, l = simulatenetwork(ntaxa; seed = s[i], kwargs...)
        push!(nets, n)
        push!(levels, l)
    end
    return nets, levels
end

"""
    simulatenetwork(ntaxa::Int, file::AbstractString; seed=nothing,
        lambda=1, mu=0.2, nu=0.2,
        shrink=true, hmin=1, hmax=ntaxa/2, levelmax=hmax)

Simulate 1 network using SiPhyNetwork repeatedly until the network
satisfies the filters, append it to `file` and return it as a tuple of 
`HybridNetwork` and level of network.

The filters are:
- exclude networks with fewer than `hmin` reticulations (trees by default)
- exclude networks with more than `hmax` reticulations
- exclude networks of level higher than `levelmax`

Before checking the level, there is the option to `shrink` 2-cycles and 3-cycles
(which are not detectable by the ABBA-BABA test or by SNaQ).

By default, the seed is not controlled. Otherwise, it is passed to R to seed
the random number generator for the SiPhyNetwork simulation.

external variables (which could be turned into options):
- proportion of the 3 kinds of reticulations: `hybridproportions`
- distribution of γs: R function `inheritance_fxn`
- probability that a proposed hybridization is successful, given the
  average genetic distance between the 2 taxa: R function `gen_dist_fxn`

# example

```repl
julia> network, level = simulatenetwork(7, file="test.phy", seed=121, verbose=false)

julia> network
HybridNetwork, Rooted Network
15 edges
15 nodes: 7 tips, 1 hybrid nodes, 7 internal tree nodes.
tip labels: t7, t10, t9, t12, ...
((((t7:0.018,t10:0.018)I1:0.358,#H15:0.0::0.155)I2:0.195,(t9:0.376,(t12:0.376)#H15:0.0::0.845)I3:0.195)I4:0.984,((t5:0.011,t6:0.011)I5:0.374,t1:0.385)I6:1.17)I7;

julia> level
1
```

R example with dangling nodes, errors with write.net etc:
```r
require(SiPhyNetwork)
inheritance_fxn = make.beta.draw(10,10)
hybridproportions = c(0.5,0.25,0.25)
set.seed(5282)
net <- sim.bdh.taxa.ssa(n=17, numbsim=1,
    lambda=1,mu=0.2, nu=0.8, hybprops=hybridproportions,
    hyb.inher.fxn = inheritance_fxn)[[1]]
write.net(net, digits=1)

library(SiPhyNetwork)
gen_dist_fxn = make.stepwise(probs = c(1,0),distances = c(0.75,Inf))
set.seed(16)
net <- sim.bdh.taxa.ssa(n=4, numbsim=1,
    lambda=1, mu=0.2, nu=0.3, hybprops=c(0.5,0.25,0.25),
    hyb.rate.fxn = gen_dist_fxn, hyb.inher.fxn = make.beta.draw(2,2), complete=FALSE)[[1]]
write.net(net, tol=1e-12)

# how to find a seed for a minimal working example:
for (i in 1:1000){
  cat("i =",i,"\n")
  set.seed(i)
  net <- sim.bdh.taxa.ssa(n=4, numbsim=1,
    lambda=1, mu=0.2, nu=0.3, hybprops=c(0.5,0.25,0.25),
    hyb.rate.fxn = gen_dist_fxn, hyb.inher.fxn = make.beta.draw(2,2), complete=FALSE)[[1]]
  if (length(net) < 2){
    cat("went extinct\n")
    next
  }
  write.net(net, tol=1e-12)
}
```
"""
function simulatenetwork(ntaxa::Int; file=nothing, seed=nothing,
        lambda=1, mu=0.2, nu=0.3,
        shrink=true, hmin=1, hmax=ntaxa/2, levelmax=hmax, verbose=true)
    levelmax>0 || error("a non-tree network must have level 1 or more")
    hmin <= hmax || error("hmin may not be higher than hmax")
    isnothing(seed) || R"set.seed"(seed)
    net = nothing
    level = nothing
    for i in 1:1000
        i += 1
        R"""
        net <- sim.bdh.taxa.ssa(n=$ntaxa, numbsim=1,
            lambda=$lambda, mu=$mu, nu=$nu, hybprops=hybridproportions,
            hyb.rate.fxn = gen_dist_fxn, hyb.inher.fxn = inheritance_fxn, complete = FALSE)[[1]]
        """
        # continue to next iteration if the lineage went completely extinct
        rcopy(R"length(net) > 1") || continue
        # convert R network to Julia network
        net_string = rcopy(R"write.net(net, tol=1e-12)")
        net = readTopology(net_string)
        shrink && PhyloNetworks.shrink3cycles!(net)
        h = net.numHybrids
        (h>=hmin && h<=hmax) || continue # to next iteration if h outside range
        bi = PhyloNetworks.blobInfo(net)[3]
        level = maximum(length.(bi))
        level <= levelmax || continue
        verbose && @info "$i iterations, level=$level"
        break
    end
    PhyloNetworks.nameinternalnodes!(net, "I") # to get be able to map internal node names in gene trees later
    isnothing(file) || writeTopology(net, file; append=true)
    return net, level
end

"""
    scale_network(net::HybridNetwork, file, med_coal_unit::AbstractFloat)

Given a network `net`, median coalescent unit to scale all branches to `med_coal_unit`,
and file path to write network to `file`, `scale_network` scales all branch lengths
in `net` by float `scaler` and writes to `file`.
"""
function scale_network(net::HybridNetwork, file, med_coal_unit::AbstractFloat)
    open(file, "w") do fout
        scaler = med_coal_unit / median([e.length for e in net.edge])
        for e in net.edge
            e.length *= scaler
        end
        writeTopology(net, fout)
    end # close file cleanly
    return net
end

"""
    simulategenetrees(net, nloci, genefile;
        lineage_distribution, lineagebygene_distribution,
        nodemapping=false)

Simulate `nloci` gene trees (1 individual/species) along `net`
using PhyloCoalSimulations.
The gene trees with branch lengths in coalescent units are written to file
"`genefile`_coal_unit". These gene trees have extra degree-2 nodes to map
them into the species network.

The model for rate variation is as follows:
* for each lineage `l` in the network, the lineage's rate `r_l` is drawn from
  `lineage_distribution`, independently across lineages
* for each lineage `l` and each gene `g`, a lineage-by-gene specific rate
  `r_lg` is drawn from `lineagebygene_distribution` (again, independently)
Then, the length of each edge in each gene tree is multiplied by `r_l × r_lg`.

The lineage rates `r_l` are saved as edge lengths in the network, in a new file
named `scalednetwork_lineagerates.phy` in the same directory as `genefile`
(see `save_lineagerates` below). `net` is *not* modified.

Output: gene trees, with edge lengths reflecting rate variation,
without any degree-2 nodes unless `nodemapping=true`.

# example

```repl
julia> net, level = simulatenetwork(4; seed=123, levelmax=1, verbose=false);

julia> writeTopology(net, round=true)
"(t2:0.448,((t1:0.202,#H10:0.0::0.06)I1:0.105,(t4:0.202,(t3:0.202)#H10:0.0::0.94)I2:0.105)I3:0.141)I4;"

julia> Random.seed!(520)

julia> gtrees = simulategenetrees(net, 3, "trial_genetrees", lineage_distribution = lognormal_meanone(0.1), lineagebygene_distribution = lognormal_meanone(0.0));

julia> writeTopology(gtrees[3], round=true)
"((t1:0.88,t2:0.928):0.52,t3:1.399,t4:1.445);"

julia> Random.seed!(520)

julia> gtrees = simulategenetrees(net, 3, "trial_genetrees"; lineage_distribution = lognormal_meanone(0.1), lineagebygene_distribution = lognormal_meanone(0.0), nodemapping=true);

julia> writeTopology(gtrees[3], round=true)
"(((t2:0.391)I4:0.182,((((t1:0.181)I1:0.11)I3:0.019,(((t3:0.176)H10:0.0)I2:0.105)I3:0.019):0.117)I4:0.182):0.33,(((t4:0.215)I2:0.105)I3:0.136)I4:0.513);"

```
"""
function simulategenetrees(net, nloci, genefile;
      lineage_distribution, lineagebygene_distribution,
      nodemapping=false)
    gtcu = genefile * "_coal_unit" # before rate variation
    treelist = simulatecoalescent(net, nloci, 1; nodemapping=true)
    writeMultiTopology(treelist, gtcu) # gene trees with lengths in coalescent units
    length(treelist) == nloci || @warn "unexpected # of gene trees" # sanity check
    lineage_rate = Dict(e.number => rand(lineage_distribution) for e in net.edge)
    # add rate for population above the network's root. find its number first.
    rootedgenumber = max(0, maximum(e.number for e in net.edge)) + 1
    push!(lineage_rate, rootedgenumber => rand(lineage_distribution))
    # initialize memory for lineage-by-gene rates
    lineagebygene_rate = Dict(enum => 0.0 for enum in keys(lineage_rate))
    for tree in treelist
        # draw lineage rates specific to this gene
        for enum in keys(lineagebygene_rate)
            lineagebygene_rate[enum] = rand(lineagebygene_distribution)
        end
        # add rate variation across lineages, and also lineage-by-gene
        for e in tree.edge
            e.length *= lineage_rate[e.inCycle] * lineagebygene_rate[e.inCycle]
        end
        # cleanup gene trees
        # with PhyloNetworks v0.15: the degree-2 root is removed too, but that's fine for seq-gen
        nodemapping || PhyloNetworks.removedegree2nodes!(tree)
    end
    netfile = joinpath(dirname(genefile), "scalednetwork_lineagerates.phy")
    save_lineagerates(net, lineage_rate, netfile) # saved inside edge lengths
    return treelist
end

"""
    save_lineagerates(net, lineage_rate, netfile)

Save the rates in `lineage_rate` in file `netfile`, as edge lengths in `net`
in two ways, resulting in writing 2 networks:
1. edge length = lineage rate * original length in `net` (in coal units)
2. edge length = lineage rate itself
"""
function save_lineagerates(net, lineage_rate, netfile)
    netcp = deepcopy(net)
    for e in netcp.edge e.length *= lineage_rate[e.number]; end # rate * length
    writeTopology(netcp, netfile)
    for e in netcp.edge e.length  = lineage_rate[e.number]; end # rate itself only
    writeTopology(netcp, netfile; append=true)
    return nothing
end

"""
    simulategenetrees_hybridlambda(net, nloci, genefile;
        distribution=lognormal_meanone(0.1), seed=nothing, verbose=false)


**old function, from before PhyloCoalSimulations**

Simulate `nloci` gene trees along `net` using HybridLambda which writes several
output files whose names will start with `genefile`;
then read the gene trees from file and modifies the length of it external,
then return them as a list of HybridNetwork's.

**Warning**: the trees in the genefile are *not* equal to the trees returned
by the function. The file has the trees as generated by HybridLambda, without
rate variation. The trees output by the function have rate variation.

The model for rate variation is as follows: each external edge in the network
is given a lineage rate according to `distribution`. Then, each external edge
in each tree is mapped to the species network: the part that maps to the
external edge in the species tree is multiplied by the lineage rate. The part
that maps to the parent edge is left as is.

**Assumption**: the calculation above assumes that there is a single individual
per species, such that each external edge in a gene tree must map to the entire
corresponding external edge in the species network (plus at least some part of
the parent edge in the network).

# example

```repl
julia> net, level = simulatenetwork(4; seed=123, levelmax=1, verbose=false);

julia> writeTopology(net, round=true)
"((t4:0.504,((t3:0.268,t1:0.268):0.057,#H8:0.0::0.424):0.179):0.16,(t2:0.325)#H8:0.34::0.576);"

julia> gtrees = simulategenetrees_hybridlambda(net, 3, "trial_genetrees"; seed=456, verbose=true);
Default Kingman coalescent on all branches.
Default population size of 10000 on all branches. 
Random seed: 456
Produced gene tree files: 
trial_genetrees_coal_unit

julia> length(gtrees)
3

julia> writeTopology(gtrees[3], round=true)
"(((t4:0.525,t2:0.544):0.27,t3:0.754):0.156,t1:0.948);"
```
"""
function simulategenetrees_hybridlambda(net, nloci, genefile; distribution=lognormal_meanone(0.1), seed=nothing, verbose=false)
    netHL = hybridlambdaformat(net) # string. weird format required by hybrid-lambda
    netHLquoted = "'$netHL'"
    hlout = ( verbose ? stdout : devnull )
    isnothing(seed) || Random.seed!(seed)
    gtcu = genefile * "_coal_unit" # name of output file created by hybrid-Lambda
    # simulate gene trees with hybrid-lambda
    hl = `$hybridlambda -spcu $netHLquoted -num $nloci -seed $seed -o $genefile`
    run(pipeline(hl; stderr = hlout))
    run(pipeline(`sed -e 's/_1//g' $gtcu`, genefile)) # get rid of _1 in names: 1 individual / species
    treelist = readMultiTopology(genefile)
    rm(gtcu)
    length(treelist) == nloci || @warn "unexpected # of trees in $gtcu" # sanity check
    taxon_edge = Dict()
    for e in net.edge
        c = PhyloNetworks.getChild(e)
        if c.leaf
            edge_rate = rand(distribution)
            push!(taxon_edge, c.name => (edge_rate, e.length, edge_rate*e.length))
        end
    end
    for tree in treelist
        for e in tree.edge
            c = PhyloNetworks.getChild(e)
            if c.leaf
                edge_info = taxon_edge[c.name]
                #second value is edge length is species tree and third value is edge rate * length of the network
                e.length = (e.length - edge_info[2]) + edge_info[3]
            end
        end
    end
    return treelist
end


"""
    varyrates!(genetrees, file, scaler::AbstractFloat, nsites=1;
                gene_distribution, edge_distribution)

Multiply the branch lengths in each gene tree from the `genetrees` list
by 2 random rates and a scaling factor:
1. a gene-specific rate, common to all the edges in the same gene tree,
   drawn from `gene_distribution` 
2. an external edge rate drawn from `edge_distribution` 
3. scaling rate `scaler`, coming from function `calc_scaling_factor`
   based on ingroup major tree length and given substitution rate.
The scaled gene trees are then written to `file`, each one preceded by the
number of sites `nsites` to simulate along this gene tree (for a format that
`seqgen` can use).
The `genetrees` list is returned after modification.

The default `nsites` of 1 to simulate unlinked SNPs, each with an independent gene tree.

# example

```repl
julia> net = readTopology("((t4:0.504,((t3:0.268,t1:0.268):0.057,#H8:0.0::0.424):0.179):0.16,(t2:0.325)#H8:0.34::0.576);");

julia> Random.seed!(456);

julia> gtrees = simulategenetrees(net, 3, "trial_genetrees", lineage_distribution = lognormal_meanone(0.7), lineagebygene_distribution = lognormal_meanone(0.7));

julia> writeTopology(gtrees[3], round=true)
"((t2:1.273,t3:1.158):0.144,t1:1.071,t4:8.827);"

julia> varyrates!(gtrees, "trial_genetrees_subst", 1.0, 500, gene_distribution = lognormal_meanone(0.2), edge_distribution = lognormal_meanone(0.0));

julia> writeTopology(gtrees[3], round=true)
"((t2:1.003,t3:0.913):0.114,t1:0.844,t4:6.956);"

```
"""
function varyrates!(genetrees, file, scaler::AbstractFloat, nsites=1;
        gene_distribution, edge_distribution)
    open(file, "w") do fout
        for tree in genetrees
            gene_rate = rand(gene_distribution)
            for e in tree.edge
                edgerate = rand(edge_distribution)
                e.length *= scaler * gene_rate * edgerate
            end
            write(fout, "[$nsites]")
            writeTopology(tree, fout)
            write(fout, "\n")
        end
    end # close file cleanly
    # writeMultiTopology(genetrees, treefile2)
    return genetrees
end

"""
    lognormal_meanone(σ)

Lognormal distribution with mean one on the positive scale and
standard deviation σ on the log scale. The mean, on the log scale,
must be μ = -σ^2/2 for the mean to be 1 on the positive scale.
"""
function lognormal_meanone(sigma)
    mu = -sigma^2 / 2
    return LogNormal(mu, sigma)
end

"""
    function calc_scaling_factor(net::HybridNetwork; sub_rate = 0.03)

Given a network `net`, this will extract the major tree, prune any outgroups, 
then calculate the ingroup major tree length. Dividing a substitution rate 
`sub_rate` by the ingroup major tree length, this returns a scaling factor 
which can be used within the `varyrates!` function.
"""
function calc_scaling_factor(net::HybridNetwork; sub_rate = 0.03)
    #extract major tree from network
    tree = majorTree(net)
    #get outgroups of the network
    outgroups = outgrouplabels(net)
    
    #prune outgroups
    for taxon in outgroups
        deleteleaf!(tree,taxon)
    end
    
    #get ingroup major tree length (imtl)
    imtl = sum(e.length for e in tree.edge)
    
    #calculate scaling factor
    scaler = sub_rate / imtl
    
    return scaler
end

"""
    simulatesequences(ngenes, nsites, treefile, sequencefile, sgseed=1234)

Simulate sequences for `ngenes` along the gene trees in `treefile`
using SeqGen, and save the sequences in `sequencefile`.

Set parameters for seq-gen:
`-m HKY`: HKY model of nucleotide substitution- different rate of 
    transitions and transversions as well as unequal frequencies 
    of the four nucleotides (base frequencies)
`-f`: base frequences of A, C, G and T respectively
`-t`: transition/transversion ratio
`-g`: the number of categories for the discrete gamma rate heterogeneity model
`-a`: a real number >0 that specifies the shape of the gamma distribution to use with gamma rate heterogeneity
`-of`: fasta output
`-q`: quiet
Parameters for Seq-Gen (-m, -a, -t, -g, -f) all adapted to reflect empirical data from the reptiles project.

Parameters altered by arguments:
`-l`: the sequence length in nucleotides of each locus, not the entire alignment ("nsites")
`-z`: seed ("sgseed")
"""
function simulatesequences(nsites, treefile, sequencefile, sgseed=1234)
    sgcmd = `$seqgen -m HKY -f 0.2 0.3 0.3 0.2 -t 3 -g 10 -a 0.35 -l $nsites -of -z $sgseed -q $treefile`;
    run(pipeline(sgcmd, stdout=sequencefile, stderr=devnull));
end

println("done reading this awesome file")
