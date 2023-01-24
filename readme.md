# pipeline to simulate data with rate variation

Code to reproduce the simulation study in:
Frankel & Ané (2023)
"Summary tests of introgression are highly sensitive to rate variation across lineages".

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
