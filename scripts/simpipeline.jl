@everywhere begin
	include("outgroup-reroot.jl")
	include("calcD.jl")
	include("infer_networks.jl")
	include("simulatedata.jl")
	using PhyloNetworks
	using Logging
	using Statistics
	using Plots
	using StatsPlots
	using ProgressMeter
end


#ten taxa
ntaxa_list = [10]
#1000 loci of 1000 sites each
nloci_nsites = [(1000,1000)]
# parameters for rate variation
lineage_dist_list       = (0.0, 0.7)
lineagebygene_dist_list = (0.0, 0.7)
gene_dist_list          = (0.0, 0.7)
edge_dist_list          = (0.0, 0.7)
substitution_rate_list  = (0.03,)
# scaling of network to make their total length more homogeneous
coalescent_median_list  = (1.0,)
# collection of all parameter sets:
rates_list = collect( (s1,s2,s3,s4,s5,cm)
    for s1 in lineage_dist_list
    for s2 in lineagebygene_dist_list
    for s3 in gene_dist_list
    for s4 in edge_dist_list
    for s5 in substitution_rate_list
    for cm in coalescent_median_list)
@show rates_list
nreplicates = 100

for ntaxa in ntaxa_list
    @show(ntaxa)

    for (ngenes, nsites) in nloci_nsites
        @show(ngenes)
        @show(nsites)

        for (lineage_dist, lineagebygene_dist, gene_dist, edge_dist, sub_rate, coal_median) in rates_list
            @show(lineage_dist)
            @show(lineagebygene_dist)
            @show(gene_dist)
            @show(edge_dist)
            @show(sub_rate)
            @show(coal_median)

            # paramset: match paramset_dir defined within simulation_pipeline
            paramset = parameterset_directory(ntaxa, ngenes, nsites,
                lineage_dist, lineagebygene_dist, gene_dist, edge_dist,
                sub_rate, coal_median)
            if isfile("$paramset/concat_dstat.csv") &&
                  isfile("$paramset/concat_d3.csv") &&
                  isfile("$paramset/concat_hyde.csv")
		        @info "$paramset is already complete- d3, dstat, hyde concatened"
                continue
            end

            replicate_seed = round.(Int, rand(nreplicates) .* 100_000_000)
            p = Progress(nreplicates, desc = "Rep #: ", barglyphs = BarGlyphs("[=> ]"), showspeed = true)
            progress_pmap(1:nreplicates, progress=p) do irep
              simulation_pipeline(ntaxa, ngenes, nsites,
                lineage_dist, lineagebygene_dist, gene_dist, edge_dist, sub_rate, coal_median,
                1, irep, 5000, replicate_seed[irep])
            end
            
            summarize_hyde_dstat_d3(paramset)
        end
    end
end

