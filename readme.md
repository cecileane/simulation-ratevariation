# pipeline to simulate data with rate variation

Code to reproduce the simulation study in: Frankel & Ané (2023) "Summary tests of introgression are highly sensitive to rate variation across lineages".

Preprint available here: <https://www.biorxiv.org/content/10.1101/2023.01.26.525396v1>

Author information:

- [Lauren Frankel](https://orcid.org/0000-0002-3002-6784): <lefrankel@wisc.edu>
- [Cécile Ané](https://orcid.org/0000-0002-4702-8217)

## goal and model

Goal: simulate sequence data under a network model with ILS
and rate variation, and quantify the performance of various methods
and assess their robustness to various types of rate variation.

Model of rate variation: there are different types of branch lengths.

- branch lengths in **coalescent units** in the species network:
  those were *not* varied.
  If they varied across lineages, then the network could be
  time inconsistent (2 paths from the root to a given hybrid
  nodes could have different lengths) and non-ultrametric
  (2 paths from the root to 2 different tips could have different lengths).
- branch lengths in **substitutions per site** on the species network
  and in gene trees: this type of rate-variation was assessed.
  * they could vary across lineages, making the network (in subst/site)
    either time-inconsistent and/or non-ultrametric
  * they could also vary across genes: some genes overall slow,
    others overall fast
  * they could also vary beyond variation across lineages and across genes.
    for example, the length in substitution/site of a given branch e
    in a given gene g, from some lineage l in the network, could be this:

    t_e × r_l × r_g × r_lg × r_e

    where:  
    t_e is the original length of edge e,  
    r_l = rate for lineage l, shared across all genes for edges in this lineage,  
    r_g = rate for gene g, shared across all edges for this gene tree  
    r_lg = rate specific for edge e in gene l (shared across all individual
    copies of that gene in the same lineage, until and after they
    coalesce within this lineage),  
    r_e = residual rate variation for edge e in gene g, independent of lineage
    if evolves in.

## code

in `scripts`:

- `simulatedata.jl` to simulate:
  * networks using (via RCall) the R package
    [SiPhyNetwork](https://github.com/jjustison/SiPhyNetwork)  
    warning: this is *not* simulating any rate variation in coalescent units.
  * gene trees within the network under the coalescent using
    [PhyloCoalSimulations](https://github.com/cecileane/PhyloCoalSimulations.jl)
  * rate variation in substitutions / site: scale branch lengths in gene trees
  * sequences along these gene trees using `seqgen`.
- `calcD.jl` calculates the ABBA-BABA test and the D3 test.
  earlier version on this [gist](https://gist.github.com/cecileane/0894d70161a2b490dbb330b96a29550a)
  to share it publically, after LF shared it with Claudia.
- `outgroup_reroot.jl` to find the outgroup of a simulated network.
  also on this [gist](https://gist.github.com/cecileane/b427b085b44ab1b968d1f5a5f1af346c)
- `infer_networks.jl` with helper functions for inferring networks/hybridization.
  currently contains functions to set up & convert files for HyDe input.
- `fasta2phylip.py` file conversion script, 
  called by `simulation_pipeline` function in `infer_networks.jl`.
- `simpipeline.jl` used for parallelization `simulation_pipeline` function
  with multiple parameter sets, and then concatenating/summarizing results with
  `summarize_hyde_dstat_d3` from `infer_networks.jl`.

## dependencies

- julia packages: see `Project.toml`
- seqgen for the simulation of sequences: Seq-Gen v1.3.4 downloaded from
  [github](git@github.com:rambaut/Seq-Gen.git) at commit from 2019-08-29.
- [HyDe](https://hybridization-detection.readthedocs.io/installation.html):
  python package, using Python 3, see
  [notes/software_install.md](notes/software_install.md#hyde)

## preliminary runs

See `choosing_params.ipynb` for preliminary runs to look at networks,
and to choose parameters to simulate networks using `SiPhyNetworks`:
`lambda`, `mu`, `nu`, the proportions of the 3 hybrid types, and the success
hybridization probability curve given the average genetic distance between the 2 taxa.

## data processing & visualization

Relevant files are in [data_processing_viz](https://github.com/cecileane/simulation-ratevariation/tree/main/data_processing_viz):

1. `data_processing.Rmd`: combines individual parameter sets' D, D3, HyDe results into `alldatasets_dstat.csv`, `alldatasets_d3.csv`, and `alldatasets_hyde_filtered.csv`, all found in [output](https://github.com/cecileane/simulation-ratevariation/tree/main/output)
2. `resultviz.Rmd`: creates the type-1 error and power plots in the paper, gets max, min, mean # subnetworks per category for plots
3. `paper_plots.jl`: creates example h1 (h) and h2 (h') plot for the paper

## output

### folder structure

the `output` folder has structure:

```
output
├── alldatasets_d3.csv
├── alldatasets_dstat.csv
├── alldatasets_hyde_filtered.csv
├── 10-taxa_1000-genes_1000-sites_0.0-lindist_0.0-linbygenedist_0.0-genedist_0.0-edgedist_0.03-subrate_1.0-coalmedian
│    └──  scalednets.txt
├── ...
└── 10-taxa_1000-genes_1000-sites_0.7-lindist_0.7-linbygenedist_0.7-genedist_0.7-edgedist_0.03-subrate_1.0-coalmedian
    └── scalednets.txt
```

``...`` collapses the 14 other folders from the 16 parameter sets ran. There were 16 parameter sets, derived from all of the possible combinations of no rate variation / yes rate variation for r_l (lindist), r_g (genedist), r_lg (linbygenedist), r_e (edgedist). Directories are named accordingly. Rates were pulled from a log-normal distribution with mean 1 and standard deviation `0.0` (no lineage rate variation) or `0.7` (yes lineage rate variation) on the log scale.

In the three csv files, they all share some column structure:

- `H1`: number of reticulations in that taxon triplet before iteratively shrinking any 2/3-cycle
- `H2`: number of reticulations in that taxon triplet after iteratively shrinking any 2/3-cycle
- `Gamma1`: true gamma (proportion of genomic contribution from a parent species to the hybrid species) value of the first hybridization event, if applicable (`NA` value if not)
- `Gamma2`: true gamma value of the second hybridization event, if applicable (`NA` value if not)
- `Gamma3`: true gamma value of the third hybridization event, if applicable (`NA` value if not)
- `coalmedian`: median edge length of branches species network was scaled to (`1.0` for all parameter sets)
- `subrate`: substitution rate in units substitutions/site/edge (`0.03` for all parameter sets)
- `lindist`: value for standard deviation on the log scale of log-normal distribution with mean 1, used to pull r_l values. `0.0` = no lineage rate variation, `0.7` = yes lineage rate variation
- `linbygenedist`: value for standard deviation on the log scale of log-normal distribution with mean 1, , used to pull r_lg values. `0.0` = no lineage-by-gene rate variation, `0.7` = yes lineage-by-gene rate variation
- `genedist`: value for standard deviation on the log scale of log-normal distribution with mean 1, , used to pull r_g values. `0.0` = no gene rate variation, `0.7` = yes gene rate variation
- `edgedist`: value for standard deviation on the log scale of log-normal distribution with mean 1, , used to pull r_e values. `0.0` = edge lineage rate variation, `0.7` = yes edge rate variation
- `simid`: unique value for each parameter set combination (values range from `1-16`)


### alldatasets_d3.csv

This is a comma-delimited file containing the results for  D3 for every possible ingroup triplet across 100 10-taxon replicate networks for each parameter set. Each row contains a triplet and its values (detailed below). Its size is 15.8MB with 79903 rows and 27 columns.

**D3 specific column structure**:

- `taxa1`: outgroup taxon, not used in these calculations
- `taxa2`, `taxa3` `taxa4`: taxa used in calculation
- `Rep`: replicate number the triplet originated from
- `TotalSites`: number of sites used
- `D13`: pairwise genetic distance between t1 (one of the taxa in `taxa2`, `taxa3`, `taxa4`) and t3 (`D3_Nonsister` below)
- `D23`: pairwise genetic distance between t2 (one of the taxa in `taxa2`, `taxa3`, `taxa4`) and t3 (`D3_Nonsister` below)
- `D3`: D3 value calculated from `D13` and `D23`
- `D3_sd`: standard deviation of D3 across 5000 bootstrap samples
- `Z`: Z-value
- `P`: P-value
- `D3_Nonsister`: taxon from `taxa2`, `taxa3`, `taxa4` that is treated as nonsister`t3` in calculations
- `taxonset`: combination of `taxa2`, `taxa3`, and `taxa4` values

### alldatasets_dstat.csv

This is a comma-delimited file containing the results for the D-statistic for every possible ingroup triplet across 100 10-taxon replicate networks for each parameter set. Each row is a triplet and its values (detailed below). Its size is 12.3MB with 79903 rows and 26 columns.

**D specific column structure**:

- `taxa1`: outgroup taxon
- `taxa2`, `taxa3`, `taxa4`: ingroup taxa
- `Rep`: replicate number the triplet originated from
- `TotalSites`: number of sites used
- `ABBAcount`: number of ABBA patterns
- `BABAcount`: number of BABA patterns
- `BBAAcount`: number of BBAA patterns
- `D`: D-value
- `Dsd`: standard deviation of D across 5000 bootstrap samples
- `Z`: Z-value
- `P`: P-value

### alldatasets_hyde_filtered.csv

This is a comma-delimited file containing the results for  HyDe for every possible ingroup triplet across 100 10-taxon replicate networks for each parameter set. We calculated the HyDe statistic for each of the 3 possible choices of the putative hybrid species for a given ingroup triplet, and we only retained the choice with the smallest p-value and applied a Bonferroni correction for multiple testing. Each row is a triplet and its values (detailed below). The file's size is 13.7MB with 79903 rows and 24 columns.

**HyDe specific column structure**:

- `outgroup`: outgroup used for HyDe statistic
- `P1`: parent 1 taxon
- `Hybrid`: hybrid taxon
- `P2`: parent 2 taxon
- `Rep`: see above
- `Zscore`: Z-score
- `origPvalue`: P-value before Bonferroni correction
- `HyDeGamma`: gamma value inferred from HyDe
- `GammaDiff`: difference between `Gamma1` and `HyDeGamma`, if applicable
- `taxonset`: combination of `P1`, `Hybrid`, and `P2` values, used for Bonferroni correction/filtering
- `Pvalue`:  P-value after Bonferroni correction

### scaled nets

`scalednets.txt`, found in each individual parameter set's directory, are files containing 100 (1 per line, for each replicate) 10-taxa networks in `HybridNetwork` format for each parameter set. These networks were simulated using `SiPhyNetwork`, then scaled all branches to have a median edge length of 1.0, which we used as measuring coalescent units (to later simulate a moderate level of ILS).
