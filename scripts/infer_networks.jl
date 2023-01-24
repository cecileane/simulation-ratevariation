#= helper code to infer networks:
1. HyDe: https://hybridization-detection.readthedocs.io/
    (suggested order)
    - Function to prune extra outgroups from FASTA `prune_fasta`
    - Function to create map file given FASTA `hyde_input`
    - Function to create HyDe command and run it, given phylip and mapfile `hyde_cmd`
2. ABBA-BABA:
    - Function to determine which site counts are ABBA & BABA, & BBAA when unordered `detect_pattern_order`
3. General:
    - Function to determine outgroups & ingroups are distinct `outgroup_check`
    - Function to check whether there is hybridization in a subnetwork given a network and list of 3 taxa `return_hybrid_number`
    - Function to run through simulation pipeline `simulation_pipeline`
    - Function to summarize individual rep results into paramset set results for D3, Dstat, & Hyde `summarize_hyde_dstat_d3`

import the file to use its functions:
include("infer_networks.jl")
=#

using PhyloNetworks
using CSV
using DataFrames
using Dates

"""
    hyde_input(phylipfile, mapfile)

Write the mapping file required by HyDe associated with the alignment in
`phylipfile`: with taxa listed in the same order. The mapping assumes a single
allele per population, so it maps each allele name to a population with the same name.

The output mapping file can then be used by function `hyde_cmd` to call the
python HyDe command `run_hyde.py` with all necessary arguments
(alignment, map, outgroup, number of sites etc.)

# example

```julia
julia> hyde_input("scripts/examples/sequential.phylip", "scripts/examples/hydemap.txt")
shell> diff scripts/examples/hydemap.txt scripts/examples/map.txt # no difference
```
"""
function hyde_input(phylipfile::AbstractString, mapfile::AbstractString)
    open(mapfile, "w") do map_fh
        open(phylipfile) do algn_fh
            readline(algn_fh) # read header (#sequences, #sites) and do nothing
            for line in eachline(algn_fh)
                taxonname = split(line, "\t")[1]
                println(map_fh, taxonname * " " * taxonname) # print twice: allele_name population_name
            end
        end
    end
end

"""
    detect_pattern_order(pattern_counts::NTuple)

Given results of `calcPatternCounts` (from calcD.jl):
7 vectors for the following number of sites:
    total, ABBA sites, BABA sites, BBAA sites,
    d13, d23, d12 number of site differences
ex. if a single locus: ([6], [1], [2], [0], [2], [3], [3])
ex. if multiple loci:  ([250, 250, 250], [7, 10, 8], [6, 8, 9], [4, 7, 4], ...)
returns ordered indices from smallest to largest values in pattern_counts[1:4],
which should be ABBA and BABA, followed by BBAA, which should be the largest of the site counts.
This resolves the need to enter taxa names in "correct" order of species tree sister pairings
and get correct order of ABBA and BABA.

# example
```julia
julia> file = "scripts/examples/tinyexample.fasta"

julia> x = calcPatternCounts(file, ["P2","P5","P3","P1"])
([6], [0], [2], [1], [2], [4], [2]))

julia> detect_pattern_order(x)
3-element view(::Vector{Int64}, 1:3) with eltype Int64:
 2
 4
 3

julia> x[detect_pattern_order(x)] # the first 2 will be used as ABBA and BABA
([0], [1], [2])

julia> x[5:7] # pairwise distances: use largest 2 for D3
([2], [4], [2])

julia> detect_pattern_order_D3(x)
3-element Vector{Int64}:
 6
 7
 5

julia> x[detect_pattern_order_D3(x)] # the first 2 will be used for D3, as d13 and d23
([4], [2], [2])
```
"""
function detect_pattern_order(pattern_counts::NTuple)
    summed_counts = [sum(pattern_counts[i]) for i in 1:4]
    #sort vector_calc from smallest to largest
    #should indices in ABBA, BABA, and BBAA order for pattern_counts
    return partialsortperm(summed_counts, 1:3)
end
function detect_pattern_order_D3(pattern_counts::NTuple)
    sc = [sum(pattern_counts[i]) for i in 5:7] # d_13, d_23, d_12
    i = argmin(sc) # sc[i] contains the smallest distance: defines the sister pair
    # use the other 2 distances (to the non-sister) for the D3 test
    idx = (i==1 ? [2,3,1] : (i==2 ? [1,3,2] : [1,2,3]) )
    return idx .+ 4 # +4 to get indices in pattern_counts
end

"""
    prune_fasta(net, input_fasta::AbstractString, output_fasta::AbstractString)

Given network `net` (either a string or a HybridNetwork object),
`prune_fasta` will determine outgroup taxa, and if more
than 1 then will prune those taxa from the FASTA alignment `input_fasta` and write 
to new fasta `output_fasta`.

Warning: assumes that ingroups and outgroups are distinct
(can check this with `outgroup_check`).
"""
function prune_fasta(net::AbstractString, input_fasta::AbstractString, output_fasta::AbstractString)
    net = readTopology(net)
    prune_fasta(net, input_fasta, output_fasta)
end
function prune_fasta(net::HybridNetwork, input_fasta::AbstractString, output_fasta::AbstractString)
    #get outgroups of the network
    outgroups = outgrouplabels(net)
    
    #check that there is more than 1 outgroup
    if length(outgroups) > 1
        #get list of outgroups to prune
        to_prune = outgroups[2:end]
        
        #open input fasta
        open(input_fasta) do f
            open(output_fasta, "w") do io
                for line in eachline(f)
                    #if taxa name matches that in to_prune list, skip that line and next (sequence)
                    if line[2:end] in to_prune
                        readline(f) # read next line with the sequence and do nothing
                    #if taxa name not in to_prune list, write to output_fasta
                    else
                        write(io, line) # >taxonname
                        write(io, "\n")
                        write(io, readline(f, keep=true)) # sequence on the next line
                    end
                end
            end
        end
    else
        print("only one outgroup, no need to prune fasta")
    end

end

"""
    prune_phylip(net::HybridNetwork, input_phylip::AbstractString, output_phylip::AbstractString)

Given network `net` (either a string or a HybridNetwork object),
`prune_phylip` will determine outgroup taxa, and if more
than 1 then will prune those taxa from the FASTA alignment `input_fasta` and write 
to new fasta `output_phylip`.

Warning: assumes that ingroups and outgroups are distinct
(can check this with `outgroup_check`).
"""
function prune_phylip(net::AbstractString, input_phylip::AbstractString, output_phylip::AbstractString)
    net = readTopology(net)
    prune_phylip(net, input_phylip, output_phylip)
end
function prune_phylip(net::HybridNetwork, input_phylip::AbstractString, output_phylip::AbstractString)
    outgroups = outgrouplabels(net)
    if length(outgroups) == 1
        @info "only one outgroup, no need to create a reduced phylip file"
        return nothing
    end
    # more than 1 outgroup: build list of outgroups to prune: all but the first
    to_prune = outgroups[2:end]
    open(input_phylip) do f
        open(output_phylip, "w") do io
            nline = 0
            for line in eachline(f)
                nline += 1
                if nline == 1
                    ntaxa = parse(Int, split(line, " ")[1]) - length(to_prune)
                    seqlen = split(line, " ")[2]
                    println(io, "$ntaxa $seqlen")
                    continue
                end
                # if the taxon is not in the list to prune: copy the line to output file
                if !(split(line, "\t")[1] in to_prune)
                    println(io, line)
                end
            end
        end
    end
    return nothing
end

"""
    outgroup_check(network::HybridNetwork)

Given network `net`, determine whether ingroups and outgroups overlap.
Returns Boolean `true` if there is no overlap, `false` if there is overlap.

# example
```julia
julia> seventaxa_network = readTopology("output/seventaxa_network")
julia> outgroup_check(seventaxa_network)
true
```
"""
function outgroup_check(network::HybridNetwork)
    groups = outgrouplabels(network, ingroup_also = true)
    outgroup_check(groups)
end

"""
    outgroup_check(outgroup_tuple::Tuple{Vector{String}, Vector{String}})

Given results of PhyloNetworks outgrouplabels(network, ingroup_also = true) stored to
`outgroup_tuple`, determine whether ingroups and outgroups overlap.
Returns Boolean `true` if there is no overlap, `false` if there is overlap.

# example
```julia
julia> seventaxa_network = readTopology("output/seventaxa_network")
julia> outgroups_ingroups = outgrouplabels(seventaxa_network, ingroup_also = true)
julia> outgroup_check(outgroups_ingroups)
true
```
"""
function outgroup_check(outgroup_tuple::Tuple{Vector{String}, Vector{String}})
    #intersects the arrays for outgroups and ingroups
    #if the length of the overlap is zero, there is no overlap
    return isempty(intersect(outgroup_tuple[1], outgroup_tuple[2]))
end

"""
    hyde_cmd(input_phylip, map_file, outgroup, prefix)

Given a phylip alignment file `input_phylip` (can be generated with `piranha_convert`),
a mapping in file `map_file` (can be generated with `hyde_input`) where the taxa
are listed in the same order as in the alignment file (HyDe requirement) and
a given outgroup, `hyde_cmd` parses out the number of taxa and sites from the
phylip file and then run the HyDe command `run_hyde.py`. Assumes that
`run_hyde.py` is in PATH.
`prefix` is used to name the output files (and location).

```julia
julia> hyde_cmd("scripts/examples/sequential.phylip", "scripts/examples/map.txt", "P1", "scripts/examples/example")

Running run_hyde.py

Reading input file.....Done.
Reading map file  .....Done.

Analyzing 12 triple(s).
Process(`run_hyde.py -i scripts/examples/sequential.phylip -m scripts/examples/map.txt -o P1 -n 5 -t 5 -s 12 --prefix scripts/examples/example`, ProcessExited(0))
```
"""
function hyde_cmd(input_phylip::AbstractString, map_file::AbstractString, outgroup::AbstractString, prefix::AbstractString)
    phylip_line = readline(input_phylip) # read first line only to get number of sites & taxa
    header = split(phylip_line, " ")
    if length(header) != 2
        error("incorrectly formatted phylip header. check that there is one space between taxa and sites")
    end
    taxa = parse(Int64, header[1])
    sites = parse(Int64, header[2])
    # don't get the outgroup taxon name from map file, bc it needs to list taxa
    # in the same order in which they appear in the alignment, for run_hyde
    run(`run_hyde.py -i $input_phylip -m $map_file -o $outgroup -n $taxa -t $taxa -s $sites --prefix $prefix`)
end

"""
    return_hybrid_number(net::HybridNetwork, taxonset::Vector)

Number of reticulations in the subnetwork of `net`, pruned to the taxa in `taxonset`.
While pruning, 2-cycles are *not* simplified.
Output:
1. H1 (integer): the number of hybridization events in the subnetwork
   *without* shrinking 2-cycles
2. H2 (integer): the number of hybridization events in the subnetwork
   after shrinking all 3-cycles
3. the gammas (vector): the gamma value(s) after shrinking 3-cycles,
   sorted in descending order. This will be an empty vector
   if there are no hybridization events (H2=0).

# example
```julia
julia> net = readTopology("((t9:0.087,t7:0.087):1.783,((t4:0.074,(t11:0.074)#H22:0.0::0.821):0.109,(t5:0.074,#H22:0.0::0.179):0.109):1.687);")
julia> return_hybrid_number(net, ["t4", "t11", "t5"])
(1, 1, [0.179])
```
"""
function return_hybrid_number(net::HybridNetwork, threetaxa::Vector)
    net = deepcopy(net)
    #first check whether taxa in threetaxa are in net
    for tip in threetaxa
        if tip ∉ tipLabels(net)
            error("Tip(s) given in threetaxa are not in network")
        end
    end
    for tip in tipLabels(net)
        if tip ∉ threetaxa # then prune it from the network
            # by default, simplify=true shrinks some 2 cycles, but could possibly miss some
            deleteleaf!(net, tip, simplify=false)
        end
    end
    # do *not* shrink 2-cycles: commented out below, bc shrinking 2/3-cycles
    # is even stronger. our results show none of these cycles are identifiable.
    # PhyloNetworks.shrink2cycles!(net)
    h1 = net.numHybrids
    if h1 == 0 # then h2=0 also
        return (0, 0, Float64[])
    end
    # next: calculate h2 and gammas that are hopefully identifiable
    PhyloNetworks.shrink3cycles!(net) # this also shrinks 2-cycles
    gammas = sort([e.gamma for e in net.edge if !e.isMajor], rev= true)
    h2 = net.numHybrids
    return (h1, h2, gammas)
end

"""
    summarize_hyde_dstat_d3(paramset::AbstractString)

Given a string `paramset` that contains the path to the directory for a parameter set (parent directory
for all the reps), summarize_hyde_dstat_d3 will construct DataFrames concatenating and analyzing the results
from HyDe, ABBA-BABA (D-statistic), and D3, if they don't already exist, for each replicate. 

For ABBA-BABA, it simply concatenates all of the dstat.csv files created in each replicate. 

For HyDe, it inserts the outgroup and rep number into the DataFrame, calculates whether there was truly hybridization 
before and after shrinking 3-cycles (H1 and H2) for each triplet plus outgroup, and inserts that and the true gamma value(s)
if applicable into the DataFrame as well. 

For D3, it does the same as above for HyDe. 

It then writes the relevant DataFrames (D-statistic, HyDe, D3) to CSVs in the `paramset` directory if they did not exist prior.
"""

function summarize_hyde_dstat_d3(paramset::AbstractString)
    #DSTAT:
    if !isfile("$(paramset)/concat_dstat.csv")
        #initialize concat_dstat file
        concat_dstat = DataFrame(taxa1= String[], taxa2= String[], taxa3= String[], taxa4= String[],
        Rep= Int[], TotalSites= Int[], ABBAcount= Int[], BABAcount= Int[], BBAAcount= Int[], 
        D= Float64[], D_sd= Float64[], Z= Float64[], P= Float64[],  
        H1=Int[], H2=Int[], Gamma1=Float64[], Gamma2=Float64[], Gamma3=Float64[])

        for (root, dirs, files) in walkdir(paramset)
            for dir in dirs
                if isfile("$(paramset)/$dir/dstat.csv")        
                #append dstat df for each rep to concat_dstat df
                    concat_dstat = vcat(concat_dstat, DataFrame(CSV.File("$(paramset)/$dir/dstat.csv")))
                else
                    error("Missing dstat.csv for rep $dir")
                end
            end
        end

        #write concatenated dstat results to CSV in paramset dir
        CSV.write("$(paramset)/concat_dstat.csv", concat_dstat)

    end
    
    #HYDE:
    if !isfile("$(paramset)/concat_hyde.csv")
        #initialize concat_hyde file
        concat_hyde = DataFrame(Outgroup= String[], P1= String[], Hybrid= String[], P2= String[], Rep=Int[],
                            Zscore= Float64[], Pvalue= Float64[], Gamma= Float64[], 
                            H1=Int[], H2=Int[], Gamma1= Float64[], Gamma2= Float64[], Gamma3= Float64[])

        #loop over all the reps in one paramset
        for (root, dirs, files) in walkdir(paramset)
            for dir in dirs
                if isfile("$(paramset)/$dir/hyde-out.txt")
                    #get outgroup and network
                    #outgroup from the map file is outgroups[1] from outgroups, ingroups = outgrouplabels(net, ingroup_also = true)
                    net = readTopology("$(paramset)/$dir/scalednetwork")
                    outgroups, ingroups = outgrouplabels(net, ingroup_also = true)
                    outgroup = outgroups[1]
                    
                    #convert hyde output csv to df (not taking all the pattern counts columns)
                    #takes 1st six columns- p1, hybrid, p2, zscore, pvalue, gamma
                    hyde_df = DataFrame(CSV.File("$(paramset)/$dir/hyde-out.txt"))[:, 1:6]
                    
                    #insert outgroup column for all rows as 1st column
                    insertcols!(hyde_df, 1, :Outgroup => outgroup)
                    #insert column for rep number after the taxa combos
                    insertcols!(hyde_df, 5, :Rep => dir)
                    #insert column for TrueHybrid with Int type and missing values 
                    hyde_df[!,:H1] = missings(Int, nrow(hyde_df))
                    hyde_df[!,:H2] = missings(Int, nrow(hyde_df))
                    #insert column for true gamma values based on network- float type and missing vals
                    hyde_df[!,:Gamma1] = missings(Float64, nrow(hyde_df))
                    hyde_df[!,:Gamma2] = missings(Float64, nrow(hyde_df))
                    hyde_df[!,:Gamma3] = missings(Float64, nrow(hyde_df))
                    
                    row_number = 1
                    for row in eachrow(hyde_df)[1:end] 
                        #initialize taxa vector- needs to reset for each line
                        taxa = []

                        #push p1, hybrid, p2 to list of taxa
                        push!(taxa, row[2], row[3], row[4])

                        #check whether there's truly hybridization for that line
                        true_hybridization = return_hybrid_number(net, taxa)

                        #add h1 and h2 booleans to their respective columns
                        hyde_df[row_number,:H1] = true_hybridization[1]
                        hyde_df[row_number,:H2] = true_hybridization[2]

                        #assign true gamma values to columns- either missings or float vals
                        gammas = true_hybridization[3]
                        gamma1 = missing
                        gamma2 = missing
                        gamma3 = missing
                        #if gamma values exist (are assigned), replace missings with their vals
                        isassigned(gammas, 1) && (gamma1 = gammas[1])
                        isassigned(gammas, 2) && (gamma2 = gammas[2])
                        isassigned(gammas, 3) && (gamma3 = gammas[3])
                        #add gamma values (or missings) to appropriate cols in df
                        hyde_df[row_number,:Gamma1] = gamma1
                        hyde_df[row_number,:Gamma2] = gamma2
                        hyde_df[row_number,:Gamma3] = gamma3

                        row_number += 1
                    end
                    
                    #append hyde df for each rep to concat_dstat df- now that there's true hybrid and outgroup in there
                    concat_hyde = vcat(concat_hyde, hyde_df)
                else
                    error("Missing hyde-out.txt for rep $dir")
                end        
            end
        end

        #rename Gamma col to HyDeGamma since it has to match HyDe col names when inserting cols from a CSV (line 364)
        rename!(concat_hyde, Dict(:Gamma => "HyDeGamma"))
        #write concatenated hyde results to CSV in paramset dir
        CSV.write("$(paramset)/concat_hyde.csv", concat_hyde)
    end
    
    #D3:
    if !isfile("$(paramset)/concat_d3.csv")
        concat_d3 = DataFrame(taxa1= String[], taxa2= String[], taxa3= String[], taxa4= String[], Rep= Int[],
                    TotalSites= Int[], D13 = Int[], D23= Int[], 
                    D3= Float64[], D_sd= Float64[], Z= Float64[], P= Float64[], 
                    H1=Int[], H2=Int[], Gamma1=Float64[], Gamma2=Float64[], Gamma3=Float64[], D3_Nonsister =AbstractString[])

        for (root, dirs, files) in walkdir(paramset)
            for dir in dirs        
                if isfile("$(paramset)/$dir/d3.csv")        
                #append dstat df for each rep to concat_d3 df
                    concat_d3 = vcat(concat_d3, DataFrame(CSV.File("$(paramset)/$dir/d3.csv")))
                else
                    error("Missing d3.csv for rep $dir")
                end
            end
        end

        #write concatenated d3 results to CSV in paramset dir
        CSV.write("$(paramset)/concat_d3.csv", concat_d3)
    end
end

"""
    parameterset_directory(ntaxa, ngenes, nsites,
        lineage_dist, lineagebygene_dist, gene_dist, edge_dist,
        sub_rate, coal_median)

Directory for output files for all replicates given a set of parameters.
Used by several functions, defined once here for consistency across the pipeline.
"""
function parameterset_directory(ntaxa, ngenes, nsites,
        lineage_dist, lineagebygene_dist, gene_dist, edge_dist,
        sub_rate, coal_median)
    joinpath("output",
      join(["$ntaxa-taxa", "$ngenes-genes", "$nsites-sites",
            "$lineage_dist-lindist", "$lineagebygene_dist-linbygenedist",
            "$gene_dist-genedist", "$edge_dist-edgedist",
            "$sub_rate-subrate", "$coal_median-coalmedian"],
            "_")
    )
end

"""
    simulation_pipeline(ntaxa::Int, ngenes::Int, nsites::Int,
            lineage_dist::AbstractFloat, lineagebygene_dist::AbstractFloat,
            gene_dist::AbstractFloat, edge_dist::AbstractFloat, 
            sub_rate::AbstractFloat, coal_median::AbstractFloat,
            nrep::Int, rep_id_start::Int, 
            nbootstrap_Dsd::Int=5000, seed::Int=nothing)

Given the number of taxa `ntaxa`, number of genes `ngenes`, number of sites `nsites`, 
the substitution rate `sub_rate` in units substitutions/site/edge, the sigma values for 
the lognormal distribution with mean one for 
    1. lineage rate shared across all genes `lineage_dist`
    2. lineage rate specific to each gene `lineagebygene_dist`
    3. individual gene rate `gene_dist`
    4. edge-within-individual-gene rate `edge_dist`
this function:

1. simulates a network with `ntaxa`, 3+ ingroups, and distinct ingroups/outgroups. Network written to `network` output file
2. scales the network to given median coalescent unit `coal_median`
3. simulates `ngenes` number of gene trees within that network, implements species specific variation across all gene 
   trees depending on distribution `lineage_dist`. Gene trees written to `gts` output file
4. calculates scaling rate for gene trees based on ingroup major tree length of network and `sub_rate`
5. varies rates of gene trees based on `lineagebygene_dist`, `edge_dist`, `gene_dist`,
   and scaling rate calculated in #3. Gene trees written to `scaledgts` output file.
6. simulates sequences with `nsites`. Alignment written to `seq.fasta`
7. converts alignment to interleaved fasta if 2+ sites `seq_linear.fasta`, converts fasta to phylip `seq.phylip`
8. runs the D-statistic on each possible 4-taxa grouping with the 1st outgroup, 
   calculates the stdev with `nbootstrap_Dsd`, writes CSV dataframe with statistics `dstat.csv`
9. runs the D3-statistic on each possible 4-taxa grouping with the 1st outgroup, 
   calculates the stdev with `nbootstrap_Dsd`, writes CSV dataframe with statistics `d3.csv`
10. runs HyDe with 1st outgroup, writes output files `hyde-out.txt` & `hyde-out-filtered.txt`

See example use in file `simpipeline.jl`, where replicates are parallelized.
"""
function simulation_pipeline(ntaxa::Int, ngenes::Int, nsites::Int,
        lineage_dist::AbstractFloat, lineagebygene_dist::AbstractFloat,
        gene_dist::AbstractFloat, edge_dist::AbstractFloat,
        sub_rate::AbstractFloat, coal_median::AbstractFloat,
        nrep::Int, rep_id_start::Int, 
        nbootstrap_Dsd::Int=5000, seed::Int=nothing)
    
    #get starting datetime
    datetime = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")

    #create folder for output if it doesn't exist already
    paramset_dir = parameterset_directory(ntaxa, ngenes, nsites,
        lineage_dist, lineagebygene_dist, gene_dist, edge_dist,
        sub_rate, coal_median)
    isdir(paramset_dir) || mkdir(paramset_dir)
    isnothing(seed) || Random.seed!(seed)

    rep_id = rep_id_start
    rep_id_end = rep_id_start + nrep
    while rep_id < rep_id_end
        network_seed = rand(1:10000)
        @show(network_seed)

        # create subfolder for this rep if it doesn't exist already
        rep_dir = joinpath(paramset_dir, string(rep_id, pad = 3))
        isdir(rep_dir) || mkdir(rep_dir)
        
        #check whether network exists already
        #if so, read in network
        if isfile("$rep_dir/network")
            net = readTopology("$rep_dir/network")
            outgroups, ingroups = outgrouplabels(net, ingroup_also = true)
            
            #start logging in append mode
            io = open("$rep_dir/log.txt", "a+")
            logger = SimpleLogger(io)
            defaultlogger = global_logger(logger)
            println(io, "\n $datetime")
            println(io, "network file already exists")
        else
            #simulate network
            #write levels to log file?
            net, levels = simulatenetwork(ntaxa, seed = network_seed)
        
            #save outgroups & ingroups in vector
            #skip if less than three ingroups
            #outgroups[1] is going to be used by ABBA-BABA & HyDe
            outgroups, ingroups = outgrouplabels(net, ingroup_also = true)
            if length(ingroups) < 3
                @info "Error: less than three ingroups"
                continue
            end
            
            #check whether ingroups vs. outgroups distinct
            #if they are not distinct, generate new network seed and simulate network again
            if !outgroup_check((outgroups, ingroups))
                @info "Error: ingroups and outgroups are not distinct for rep"
                continue
            end
            
            #write network to rep_dir folder if it doesn't exist already
            writeTopology(net, "$rep_dir/network")
            
            #now with successful rep, start logging
            io = open("$rep_dir/log.txt", "w+")
            logger = SimpleLogger(io)
            defaultlogger = global_logger(logger)
            write(io, """
            $datetime
            ntaxa = $ntaxa
            ngenes = $ngenes
            nsites = $nsites
            lineage_dist = $lineage_dist
            lineagebygene_dist = $lineagebygene_dist
            gene_dist = $gene_dist
            edge_dist = $edge_dist
            sub_rate = $sub_rate
            network_seed = $network_seed
            """)
        end
        
        #scale network
        if isfile("$rep_dir/scalednetwork") 
            scaled_net = readTopology("$rep_dir/scalednetwork")
            println(io, "scaled network already exists")
        else
            scaled_net = scale_network(net, "$rep_dir/scalednetwork", coal_median)
        end

        #simulate gts, write to rep_dir folder
        genetrees_seed = rand(1:10000)
        @show(genetrees_seed)
        if isfile("$rep_dir/gts_coal_unit")
            gts = readMultiTopology("$rep_dir/gts_coal_unit")
            println(io, "gene trees already exist")
        else
            gts = nothing
            try # was important when using HybridLambda, which often segfaulted
                gts = simulategenetrees(scaled_net, ngenes, "$rep_dir/gts",
                            lineage_distribution = lognormal_meanone(lineage_dist),
                            lineagebygene_distribution = lognormal_meanone(lineagebygene_dist))
            catch e 
                @error "Simulating gene trees didn't work with this network, retrying new topology" writeTopology(scaled_net)
                continue
            end
        end

        if isfile("$rep_dir/scaledgts")
            println(io, "scaled gene trees already exist")
        else
            #calculate scaling factor based on network, for use in varyrates!
            scaler = calc_scaling_factor(scaled_net, sub_rate = sub_rate)
            # vary rates on gene trees: write 'scaledgts' with lengths in subst/site
            varyrates!(gts, "$rep_dir/scaledgts", scaler, nsites, gene_distribution = lognormal_meanone(gene_dist), edge_distribution = lognormal_meanone(edge_dist))
        end

        # simulate sequences with seq-gen
        if isfile("$rep_dir/seq.fasta")
            println(io, "sequence already exists")
        else
            sequence_seed = rand(1:10000)
            @show(sequence_seed)
            println(io, "sequence_seed = $sequence_seed")
            simulatesequences(nsites, "$rep_dir/scaledgts", "$rep_dir/seq.fasta", sequence_seed)
        end
        
        #convert interleaved fasta to linear if 2+ loci
        #rename seq.fasta to seq_linear.fasta if only 1 locus (automatically in linear format)
        if ngenes == 1 && !isfile("$rep_dir/seq_linear.fasta")
            mv("$rep_dir/seq.fasta", "$rep_dir/seq_linear.fasta")
        elseif ngenes > 1 && !isfile("$rep_dir/seq_linear.fasta")
            reformat_fasta_linear("$rep_dir/seq.fasta", "$rep_dir/seq_linear.fasta")
        end
        
        #convert fasta alignment to phylip
        #should output in $rep_dir
        #rename phylip to seq.phylip instead of seq_linear.phylip if seq_linear.fasta is input
        if !isfile("$rep_dir/seq.phylip")
            run(`python3 scripts/fasta2phylip.py -f $rep_dir/seq_linear.fasta`)
            mv("$rep_dir/seq_linear.phylip", "$rep_dir/seq.phylip")
        end
        
        # ABBA-BABA & D3
        if !isfile("$rep_dir/d3.csv") || !isfile("$rep_dir/dstat.csv")
            ningroups = length(ingroups)
            dstat_df = DataFrame(taxa1= String[], taxa2= String[], taxa3= String[], taxa4= String[], Rep= Int[],
                        TotalSites= Int[], ABBAcount= Int[], BABAcount= Int[], BBAAcount= Int[], 
                        D= Float64[], D_sd= Float64[], Z= Float64[], P= Float64[], 
                        H1=Int[], H2=Int[], Gamma1=Float64[], Gamma2=Float64[], Gamma3=Float64[])
            allowmissing!(dstat_df)

            d3_df = DataFrame(taxa1= String[], taxa2= String[], taxa3= String[], taxa4= String[], Rep= Int[],
                        TotalSites= Int[], D13 = Int[], D23= Int[], 
                        D3= Float64[], D_sd= Float64[], Z= Float64[], P= Float64[], 
                        H1=Int[], H2=Int[], Gamma1=Float64[], Gamma2=Float64[], Gamma3=Float64[], D3_Nonsister=AbstractString[])
            allowmissing!(d3_df)

            for in1 in 1:(ningroups-2) for in2 in (in1+1):(ningroups-1) for in3 in (in2+1):ningroups 
                taxa_names = [outgroups[1], ingroups[in1], ingroups[in2], ingroups[in3]]
                taxa_outgrouplast = taxa_names[[2,3,4,1]]
                
                #check whether there is true hybridization
                h1 = nothing
                h2 = nothing
                #gamma values will be missing unless there are hybridization events
                gamma1 = missing
                gamma2 = missing
                gamma3 = missing
                try
                    h1, h2, gammas = return_hybrid_number(scaled_net, taxa_names[2:4])
                    #if there are gamma values, assign to gamma 1-3
                    isassigned(gammas, 1) && (gamma1 = gammas[1])
                    isassigned(gammas, 2) && (gamma2 = gammas[2])
                    isassigned(gammas, 3) && (gamma3 = gammas[3])
                catch e
                    rep_id -= 1
                    @error "return_hybrid_number failed, perhaps due to segfault with Hybrid-Lambda & hybrid ladder network" writeTopology(scaled_net)
                    continue        
                end
                        
                pattern_counts = calcPatternCounts("$rep_dir/seq_linear.fasta", taxa_outgrouplast, repeat([nsites],ngenes))
                #sum the arrays within the pattern_counts tuple
                summed_pattern_counts = sum.(pattern_counts)
                total_counts = summed_pattern_counts[1]
                
                if !isfile("$rep_dir/dstat.csv")
                    #get index of ABBA, BABA, BBBA in pattern_counts
                    #use that to save count #s
                    abba_baba_indices = detect_pattern_order(pattern_counts)
                    abba_counts = summed_pattern_counts[abba_baba_indices[1]]
                    baba_counts = summed_pattern_counts[abba_baba_indices[2]]
                    bbaa_counts = summed_pattern_counts[abba_baba_indices[3]]
                    
                    if abba_counts + baba_counts == 0
                        d_value = missing
                        d_sd = missing
                        d_z = missing
                        d_p = missing
                    else
                        # calculate d, d_sd, z, p_val
                        d_value = calcD(abba_counts, baba_counts)
                        d_sd = (calcD_sd_locusbootstrap(pattern_counts[abba_baba_indices[1:2]]..., nbootstrap_Dsd))
                        d_z = d_value/d_sd
                        d_p = ccdf(Normal(), abs(d_z))*2
                        #= same for Dp, if desired (todo): (would need to stick in a diff elseif block where abba + baba + bbaa != 0)
                        dp = calcDp(pattern_counts[abba_baba_indices[1:3]]...)
                        dp_sd = calcDp_sd_locusbootstrap(pattern_counts[abba_baba_indices[1:3]]..., nbootstrap_Dsd)
                        z, pvalue, save to file
                        =#
                    end

                    #write D-stat to dataframe
                    push!(dstat_df, (taxa_names..., rep_id, total_counts, abba_counts, baba_counts, bbaa_counts,
                            d_value, d_sd, d_z, d_p, h1, h2, gamma1, gamma2, gamma3))
                end

                if !isfile("$rep_dir/d3.csv")
                    #get indices for d23 and d13 in pattern_counts
                    d3_indices = detect_pattern_order_D3(pattern_counts) # permutation of 5:7
                    d13 = summed_pattern_counts[d3_indices[1]]
                    d23 = summed_pattern_counts[d3_indices[2]]
                    d3_nonsister = taxa_outgrouplast[d3_indices[3]-4]

                    # calculate D3 and its bootstrap-based p-value
                    if d13 + d23 == 0
                        d3_value = missing
                        d3_sd = missing
                        d3_z = missing
                        d3_p = missing
                    else 
                        d3_value = calcD3(pattern_counts[d3_indices[1:2]]...)
                        d3_sd = calcD3_sd_locusbootstrap(pattern_counts[d3_indices[1:2]]..., nbootstrap_Dsd)
                        d3_z = d3_value/d3_sd
                        d3_p = ccdf(Normal(), abs(d3_z))*2
                    end

                    #write D3 to dataframe
                    push!(d3_df, (taxa_names..., rep_id, total_counts, d13, d23, d3_value, d3_sd, d3_z, d3_p, h1, h2, gamma1, gamma2, gamma3, d3_nonsister))
                end
                
            end end end

            #write Dstat dataframe to CSV
            isfile("$rep_dir/dstat.csv") || CSV.write("$rep_dir/dstat.csv", dstat_df)

            #write D3 dataframe to CSV
            isfile("$rep_dir/d3.csv") || CSV.write("$rep_dir/d3.csv", d3_df)

            # now that all the criteria passed: successful rep
            # make sure to increment rep *after* rep number is included in csvs for dstat, d3
            @show(rep_id)
            rep_id += 1
        end

        # run HyDe on all ingroup triplets + first outgroup
        if !isfile("$rep_dir/hyde-out.txt")
            phylip = "seq.phylip"
            if length(outgroups) > 1 # then remove sequences for all but the first outgroup
                prune_phylip(scaled_net, "$rep_dir/seq.phylip", "$rep_dir/seq_pruned.phylip")
                phylip = "seq_pruned.phylip"
            end
            hyde_input("$rep_dir/$phylip", "$rep_dir/map.txt") # generate mapping file: taxa in matching order
            hyde_cmd("$rep_dir/$phylip", "$rep_dir/map.txt", outgroups[1], "$rep_dir/hyde")
        end
        
        close(io)
        global_logger(defaultlogger)

    end
end
