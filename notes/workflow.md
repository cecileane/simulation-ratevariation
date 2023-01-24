# Workflow outline

Steps in simulation pipeline (`simulation_pipeline` function in `scripts/infer_networks.jl`) 
and parameters (either varied or fixed) at each step:

1. Simulate networks: `simulatenetwork`
    - Number of taxa: `ntaxa` (`10`)
    - Lambda (`1`), mu (`0.2`), and nu (`0.2`) values
    - Proportion of 3 kinds of reticulations `hybridproportions` (external, `[0.5,0.25,0.25]` = generative, degenerative, neutral)
    - Distribution of γs: R function `inheritance_fxn` (external, `make.beta.draw(2, 2)`)
    - Probability of hybridization success given the average genetic distance:
      R function `gen_dist_fxn` (external, `make.stepwise(probs = c(1,0),distances = c(0.75,Inf))`)
    - Minimum number of reticulations: `hmin` (`1`)
    - Maximum number of reticulations: `hmax` (ntaxa/2 = `5`)
    - Maximum level: `levelmax` (5, same as `hmax`)

2. Scale network `scale_network` to have a median coalescent unit `med_coalunit` (`1.0`)

3. Simulate gene trees: `simulategenetrees`
    - Number of loci gene trees along given network `nloci` (`1000`)
    - Lineage rates shared across all gene trees drawn from `lineage_distribution` (`0.0` or `0.7`)
    - Lineage rates specific to each gene drawn from `lineagebygene_distribution` (`0.0` or `0.7`)
    - Remove degree-2 nodes (`nodemapping = false`)

4. Vary rates on gene trees: `varyrates!`
    - Scaler (calculated by `calc_scaling_factor(scaled_net, sub_rate = sub_rate)` where sub_rate = `0.03`
      and using ingroup major tree length from `scaled_net`)
    - Number of sites (`1000`)
    - Gene rate distribution `gene_distribution` (`0.0` or `0.7`)
    - Edge rate distribution (specific to a single gene) `edge_distribution` (`0.0` or `0.7`)

5. Simulate sequences: `simulatesequences`
    (parameters below are all fixed)
    - Model of nucleotide substitution `-m HKY` (HKY model = different rate of 
      transitions and transversions as well as unequal frequencies 
      of the four nucleotides (base frequencies))
    - Base frequencies `-f 0.2 0.3 0.3 0.2` (A, C, G and T respectively)
    - Transition/transversion ratio `-t 3`
    - Number of categories for discrete gamma rate heterogeneity model `-g 10`
    - Shape of gamma distribution `-a 0.35`

6. Hybrid detection / network inference
    - ABBA-BABA
    - HyDe
    - D3

Compilation of reps per parameter set happens with function `summarize_hyde_dstat_d3`
and compilation of various parameter sets in R afterward.
