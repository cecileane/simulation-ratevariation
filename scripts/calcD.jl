#=
inspired from evobiR::CalcD but requires package "seqinr" only
used in [Blair & Ané 2020](https://doi.org/10.1093/sysbio/syz056)
=#

using BioSymbols
using BioSequences
using FASTX
using PhyloNetworks
using Random
using Statistics
using Distributions
using DataFrames
using CSV

"""
    calcPatternCounts(filename, taxa, nsites_perlocus=nothing; ambig=:D)

Take an alignment in a file in fasta format, a vector `taxa` of 4 taxon names,
the number of sites at each locus (optionally), and calculate
7 vectors for the following number of sites:

1. total
2. ABBA sites
3. BABA sites
4. BBAA sites
5. d_13
6. d_23
7. d_12
   where d_ij is the number of sites at which taxa i and j differ.

These numbers of sites are calculated by locus, so are stored in a vector
for each type of site. If `nsites_perlocus` is not provided, it is assumed
that the alignment has a single locus and each vector has a single value.
Note that the BBAA sites are *not* used by the D statistic, but are used by
the Dp statistic. The d_ij counts are used by the D3 statistic.

Warning: the outgroup should be the last of the 4 taxa for later calculation of D3.

By default, sites with any ambiguity are dropped (`ambig=:D`) before calculating
these 4 numbers. Alternatively, they could be resolved (`ambig=:R`) randomly.
Otherwise, ambiguous sites are just left alone (e.g. ignored: `ambig=:I`).
In that case for instance, `N--N` would count as an ABBA pattern.

**Warnings**:
1. this is assuming that the fasta file is in *sequential* format,
*not interleaved*. More specifically, each taxon name should appear only once.
For this one time, the corresponding sequence could be broken across several lines though.
Use `reformat_fasta_linear` to reformat an interleaved fasta file.

2. use PhyloNetworks v0.15.0 (or above) to avoid a bug with some ambiguous sites.

# example

```julia
julia> file = "examples/tinyexample.fasta";

julia> calcPatternCounts(file, ["P2","P1","P3","typo"])
(missing, missing, missing, missing, missing, missing, missing)

julia> calcPatternCounts(file, ["P2","P1","P3","P5"])
([6], [1], [2], [0], [2], [3], [3])

julia> calcPatternCounts(file, ["P2","P1","P3","P5"], [4,3]) # 2 loci: with 4 sites and 3 sites
([3, 3], [1, 0], [0, 2], [0, 0], [1, 1], [0, 3], [1, 2])

julia> calcPatternCounts(file, ["P2","P1","P3","P5"], repeat([4],2)) # 2 loci each with 4 sites
ERROR: the sum of locus lengths is not the total alignment length

julia> calcPatternCounts(file, ["P1","P2","P3","P5"]) # switching taxa 1 & 2 changes abba into baba
([6], [2], [1], [0], [3], [2], [3])

julia> calcPatternCounts(file, ["P1","P2","P5","P3"]) # switching taxa 3 & 4 changes abba into baba
([6], [1], [2], [0], [1], [2], [3])

julia> calcPatternCounts(file, ["P5","P3","P1","P2"]) # switching 1,2 with 3,4 has no effect
([6], [1], [2], [0], [1], [3], [4])

julia> calcPatternCounts(file, ["P1","P2","P3","P4"]) # one more site excluded, gap for P4
([5], [0], [0], [1], [2], [2], [2])

julia> calcPatternCounts(file, ["P1","P3","P2","P4"]) # switching 2 & 3 switches bbaa into baba
([5], [0], [1], [0], [2], [2], [2])

julia> calcPatternCounts(file, ["P1","P2","P3","P4"]; ambig=:R) # ambiguous Y resolved
([6], [0], [0], [1], [2], [3], [3])
```
"""
function calcPatternCounts(file::String, taxa::AbstractVector,
            nsites_perlocus = nothing; ambig=:D)
  length(Set(taxa)) == 4 || error("we need 4 distinct taxa for ABBA and BABA patterns")

  # read alignment and make it a matrix
  species, alg = PhyloNetworks.readFastaToArray(file)

  # check that the 4 taxa are present, and all have same alignment length
  # sometimes, 1 or more taxa might be missing a gene of interest
  ind = indexin(taxa, species) # species[ind] in correct order to look at ABBA and BABA
  # if any of the 4 taxa is missing: return missing values
  any(isnothing, ind) && return (missing,missing,missing,missing,missing,missing,missing)

  # alignment for the (sub)set of 4 taxa in the appropriate order
  alg_4 = view(alg, ind)
  # check that all 4 sequences have the same # of sites, otherwise return missing values
  nsites = length(alg_4[1])
  if !isnothing(nsites_perlocus)
    nsites == sum(nsites_perlocus) ||
      error("the sum of locus lengths is not the total alignment length")
    nloci = length(nsites_perlocus)
  else
    nloci = 1
    nsites_perlocus = [nsites]
  end
  all(x -> length(x)==nsites, alg_4) || return (missing,missing,missing,missing,missing,missing,missing)

  ## deal with ambiguity in sequence data
  # R: A or G
  # Y: C or T
  # S: G or C
  # W: A or T
  # K: G or T
  # M: A or C
  if ambig in [:D, :R]
    # D: drop all ambiguous sites
    # R: resolve, but first drop any site having a 3- or 4-ambiguity
    target = (BioSymbols.ACGT..., DNA_R, DNA_Y, DNA_S, DNA_W, DNA_K, DNA_M)
    j = nsites # start from the last site, to keep indices valid after deletion
    # for j in nsites:-1:1
    for locus in nloci:-1:1
      for _ in nsites_perlocus[locus]:-1:1
        keepj = (ambig==:D ?
               all(s -> BioSymbols.iscertain(s[j]), alg_4) :
               all(s -> s[j] in target, alg_4) )
        if !keepj
          for s in alg_4 # delete the site
            deleteat!(s,j)
          end
          nsites_perlocus[locus] -= 1 # decrement # sites in this locus
        end
        j -= 1
      end
    end
  end
  nsites = length(alg_4[1]) # possibly less than before
  sum(nsites_perlocus) == nsites ||
    error("bug: after removing uncertain sites, the loci and total length don't match...")
  ## resolve remaining ambigous sites randomly
  if ambig == :R
    function resolver(x)
      if iscertain(x) return x; end
      if x==DNA_R return rand([DNA_A, DNA_G]); end
      if x==DNA_Y return rand([DNA_C, DNA_T]); end
      if x==DNA_S return rand([DNA_G, DNA_C]); end
      if x==DNA_W return rand([DNA_A, DNA_T]); end
      if x==DNA_K return rand([DNA_G, DNA_T]); end
      if x==DNA_M return rand([DNA_A, DNA_C]);
      else error("weird base: $x"); end
    end
    for a in alg_4
      for j in 1:nsites
        isambiguous(a[j]) || continue
        a[j] = resolver(a[j])
      end
    end
  end
  ## calculate the number of sites with pattern ABBA or BABA
  abba_all = Int[] # all loci
  baba_all = Int[]
  bbaa_all = Int[]
  d13_all = Int[]
  d23_all = Int[]
  d12_all = Int[]
  pattern = Vector{DNA}(undef, 4) # to avoid garbage collection
  j = 0 # index in the concatenated alignment
  for locus in 1:nloci
    abba = 0
    baba = 0
    bbaa = 0
    d13 = 0
    d23 = 0
    d12 = 0
    for _ in 1:nsites_perlocus[locus]
      j += 1
      for i in 1:4 pattern[i] = alg_4[i][j]; end
      if pattern[1] != pattern[3] d13 += 1; end
      if pattern[2] != pattern[3] d23 += 1; end
      if pattern[1] != pattern[2] d12 += 1; end
      length(unique(pattern)) == 2 || continue   # focus on biallelic sites
      if pattern[1] == pattern[2]
        if pattern[3] == pattern[4] # else bbba or bbab: do not count
          bbaa += 1
        end
        continue
      end
      # at this point: pattern[1] != pattern[2] so we focus on xy.. patterns
      pattern[3] != pattern[4] || continue   # focus on ..ba patterns
      if pattern[3] == pattern[1]
        baba += 1
      else # pattern[3] == pattern[2]
        abba += 1
      end
    end
    push!(abba_all, abba)
    push!(baba_all, baba)
    push!(bbaa_all, bbaa)
    push!(d13_all, d13)
    push!(d23_all, d23)
    push!(d12_all, d12)
  end
  return nsites_perlocus, abba_all,baba_all,bbaa_all, d13_all,d23_all,d12_all
end

"""
    calcD(abba, baba)
    calcDp(abba, baba, bbaa)
    calcD3(d13, d23)
    calcD3v2(d13, d23, d12)

Take a vector `x` of ABBA counts and a vector `y` of BABA counts, assuming that
in each list: 1 count corresponds to 1 gene. Then calculate the D statistics for
the concatenated data, that is: [sum(ABBA) - sum(BABA)] / [sum(ABBA) + sum(BABA)]

For other statistics:
Dp = [sum(ABBA) - sum(BABA)] / [sum(ABBA) + sum(BABA) + sum(BBAA)]
D3   = [sum(d23) - sum(d13)] / [sum(d23) + sum(d13)]
D3v2 = [sum(d23) - sum(d13)] /  sum(d12)

# examples

```julia
julia> calcD(20, 5)
0.6

julia> calcD([5,4,6,5], [3,0,0,2])
0.6

julia> file = "examples/tinyexample.fasta";

julia> pattern_counts = calcPatternCounts(file, ["P2","P1","P3","P5"], [4,3]) # 2 loci
([3, 3], [1, 0], [0, 2], [0, 0], [1, 1], [0, 3], [1, 2])

julia> d_value = calcD(pattern_counts[2:3]...)
-0.3333333333333333

julia> dp_value = calcDp(pattern_counts[2:4]...)
-0.3333333333333333

julia> d3_value = calcD3(pattern_counts[5:6]...)
0.2

julia> d3v2_value = calcD3v2(pattern_counts[5:7]...)
0.3333333333333333
```
"""
function calcD(x, y)
  sx = sum(x); sy = sum(y)
  return (sx - sy)/(sx + sy)
end
function calcDp(x, y, z)
  sx = sum(x); sy = sum(y)
  return (sx - sy)/(sx + sy + sum(z))
end
function calcD3(x, y)
  sx = sum(x); sy = sum(y)
  return (sy - sx)/(sx + sy)
end
function calcD3v2(x, y, z)
  sx = sum(x); sy = sum(y)
  return (sy - sx)/sum(z)
end

"""
     calcDztest_concatenated(x,y)
     calcDsd_concatenated(x, y, nbootstrap=5000)

Take a number of ABBA sites `x` and number of BABA sites `y` (possibly
calculated by `calcPatternCounts` from a concatenated alignment)
and test whether the data are consistent with mean D=0.

Output: D statistic, its standard deviation, and test p-value (2-tailed).

Assumptions:

1. independent sites in the single (concatenated) alignment
   from which `x,y` come from.
2. the number of "informative" sites for the ABBA-BABA test, that is, the
   number of ABBA + BABA sites is considered fixed. In other words, the test
   conditions on the number of informative sites.

Note that D = (x-y)/(x+y) becomes phat - (1-phat) = 2*(phat -0.5)
where phat = estimated proportion of ABBA sites among informative sites,
that is, among ABBA & BABA sites.
Testing "mean D=0" is equivalent to testing p=0.5, which amounts to doing
is simple binomial test, or chi-square goodness-of-fit test if all sites
are assumed independent.

The first version does a z test (equivalent to chi-square goodness-of-fit test).

The second version calculates the standard deviation of the D statistics
using bootstrap under the null hypothesis that p_abba = p_baba, re-sampling the
ABBA-BABA sites only (*not* resampling all sites: see assumption 2 above).
Then calculates z = D/sd and pvalue = two-tailed probability beyond z under N(0,1).

This second version should be modified to allow for linkage between sites
(in which case the z-test becomes invalid).
For now, this version is equivalent to the z-test but slower for large counts,
and more reliable for small count provided that `nbootstrap` is large.

# examples

```julia
julia> calcDztest_concatenated(20,5) # D=0.6, z=3, p-value=0.0027
(0.6, 3.0000000000000004, 0.002699796063260186)

julia> calcDztest_concatenated(5,20) # D=-0.6, same p-value as above
(-0.6, -3.0, 0.002699796063260186)

julia> calcDsd_concatenated(20, 5, 10000) # same D, similar z & p.
(0.6, 3.0081772816666144, 0.0026281977308263323)
```
"""
function calcDztest_concatenated(x,y)
  n = x+y
  D = (x-y)/n
  phat = x/n
  z = (phat-0.5)*sqrt(n)*2 #  = (phat - p0) / sqrt(p0 * (1-p0) * 1/n) when p0=0.5
  pvalue = ccdf(Normal(), abs(z))*2
  return D, z, pvalue
end

function calcDsd_concatenated(x, y, nbootstrap=5000)
  n = x+y
  D = (x-y)/n
  d = Distributions.Binomial(n, 0.5) # null hypothesis for D=0: x=y=n/2
  bx = rand(d, nbootstrap)
  bD = (2*bx .- n) ./ n # y=n-x so x-y = x-(n-x) = 2x-n
  z = D / Statistics.std(bD)
  pvalue = ccdf(Normal(), abs(z))*2
  return D, z, pvalue
end

"""
    calcD_sd_locusbootstrap( abba, baba,       nbootstrap=5000)
    calcDp_sd_locusbootstrap(abba, baba, bbaa, nbootstrap=5000)
    calcD3_sd_locusbootstrap(  d13, d23,      nbootstrap=5000)
    calcD3v2_sd_locusbootstrap(d13, d23, d12, nbootstrap=5000)

Take a list of ABBA counts and a list of BABA counts, like `calcD`.
For Dp, also take a list of BBAA counts.
For D3 (original and version 2): take a list of pairwise distances.
Calculate the standard deviation of the named D statistics using bootstrap,
re-sampling loci. The sites are *not* resampled: only the loci are.
This is assuming independent loci, but possibly linked sites within loci.

There is one case of concern: if all resampled loci have no informative
(ABBA & BABA; or d13 & d23) sites.
In this case, the D statistic is undefined: 0-0/(0+0)
If that happens, the corresponding bootstrap replicate is dropped.
`nbootstrap` is the total # of valid bootstrap replicates (not counting
dropped reps).

Warning: the bootstrap procedure requires a large sample size to be valid.
Here, this means a large number of loci. In the example below, we only have
4 loci, so the results from this bootstrap cannot be trusted. But 4 loci
makes the example run fast.

# examples with Patterson's D

```julia
julia> abba_4loci = [5,4,6,5]; baba_4loci = [3,0,0,2];

julia> D_value = calcD(abba_4loci, baba_4loci)
0.6

julia> D_sd = calcD_sd_locusbootstrap(abba_4loci, baba_4loci, 10_000)
0.17212358417715387

julia> z = D_value/D_sd
3.485867453134517

julia> pval = ccdf(Normal(), abs(z))*2
0.0004905439979869396

julia> abba_4loci = [11,0,0,9]; baba_4loci = [3,0,0,2];

julia> D_value = calcD(abba_4loci, baba_4loci)
0.6

julia> D_sd = calcD_sd_locusbootstrap(abba_4loci, baba_4loci, 10_000)
0.02462521664116594

julia> calcD_sd_locusbootstrap([1,2], [1,2,3], 10)
ERROR: different loci for ABBA and BABA counts

julia> calcD_sd_locusbootstrap([0,-1,0], [0,0,0], 10)
ERROR: some ABBA counts are negative

julia> calcD_sd_locusbootstrap([0,0,0], [0,0,0], 10)
ERROR: all ABBA and BABA are 0
```

# examples using a file, with Hahn & Hibbins's D3

```julia
julia> file = "examples/tinyexample.fasta";

julia> pattern_counts = calcPatternCounts(file, ["P2","P1","P3","P5"], [4,3]) # 2 loci

julia> d3_value = calcD3(pattern_counts[5:6]...)
0.2

julia> d3_sd = calcD3_sd_locusbootstrap(pattern_counts[5:6]..., 10_000)
0.5740339703361726

julia> z = d3_value/d3_sd
0.34841143614353276

julia> p_value = ccdf(Normal(), abs(z))*2
0.7275312151400344
```
"""
function calcD_sd_locusbootstrap(x, y, nbootstrap=5000)
  nloci = length(x)
  length(y) == nloci || error("different loci for ABBA and BABA counts")
  bD = Vector{Float64}(undef, nbootstrap)
  b = 1
  all(x .>= 0) || error("some ABBA counts are negative")
  all(y .>= 0) || error("some BABA counts are negative")
  any(x .> 0) || any(y .> 0) || error("all ABBA and BABA are 0")
  while b <= nbootstrap
    bsx=0; bsy=0
    bloci = rand(1:nloci, nloci) # sample loci with replacement
    for locus in bloci
      bsx += x[locus]
      bsy += y[locus]
    end
    sum_abbababa = bsx + bsy # none of the loci have any of the ABBA or BABA sites
    if sum_abbababa ≈ 0      # D would be NaN, then SD would be NaN
      continue # try again: don't increment b (bootstrap rep)
    end
    bD[b] = (bsx - bsy)/sum_abbababa
    b += 1
  end
  return Statistics.std(bD)
end
function calcDp_sd_locusbootstrap(x, y, z, nbootstrap=5000)
  nloci = length(x)
  length(y) == nloci || error("different loci for ABBA and BABA counts")
  length(z) == nloci || error("different loci for ABBA and BBAA counts")
  bD = Vector{Float64}(undef, nbootstrap)
  b = 1
  all(x .>= 0) || error("some ABBA counts are negative")
  all(y .>= 0) || error("some BABA counts are negative")
  all(z .>= 0) || error("some BBAA counts are negative")
  any(x .> 0) || any(y .> 0) || any(z .> 0) || error("all ABBA BABA & BBAA are 0")
  while b <= nbootstrap
    bsx=0; bsy=0; bsz=0;
    bloci = rand(1:nloci, nloci) # sample loci with replacement
    for locus in bloci
      bsx += x[locus]
      bsy += y[locus]
      bsz += z[locus]
    end
    sum_denom = bsx + bsy + bsz
    sum_denom ≈ 0  && continue # Dp would be NaN, then SD would be NaN: try again
    bD[b] = (bsx - bsy)/sum_denom
    b += 1
  end
  return Statistics.std(bD)
end
function calcD3_sd_locusbootstrap(x, y, nbootstrap=5000)
  nloci = length(x)
  length(y) == nloci || error("different loci for d23 and d13 counts")
  bD = Vector{Float64}(undef, nbootstrap)
  b = 1
  all(x .>= 0) || error("some d13 counts are negative")
  all(y .>= 0) || error("some d23 counts are negative")
  any(x .> 0) || any(y .> 0) || error("all d13 & d23 are 0")
  while b <= nbootstrap
    bsx=0; bsy=0
    bloci = rand(1:nloci, nloci) # sample loci with replacement
    for locus in bloci
      bsx += x[locus]
      bsy += y[locus]
    end
    sum_denom = bsx + bsy
    sum_denom ≈ 0 && continue # D3 would be NaN, then SD would be NaN
    bD[b] = (bsy - bsx)/sum_denom
    b += 1
  end
  return Statistics.std(bD)
end
function calcD3v2_sd_locusbootstrap(x, y, z, nbootstrap=5000)
  nloci = length(x)
  length(y) == nloci || error("different loci for d23 and d13 counts")
  length(z) == nloci || error("different loci for d12 and d13 counts")
  bD = Vector{Float64}(undef, nbootstrap)
  b = 1
  all(x .>= 0) || error("some d13 counts are negative")
  all(y .>= 0) || error("some d23 counts are negative")
  all(z .>= 0) || error("some d12 counts are negative")
  any(z .> 0)  || error("all d12 are 0")
  while b <= nbootstrap
    bsx=0; bsy=0; bsz=0
    bloci = rand(1:nloci, nloci) # sample loci with replacement
    for locus in bloci
      bsx += x[locus]
      bsy += y[locus]
      bsz += z[locus]
    end
    bsz ≈ 0 && continue # Dp3v2 would be NaN, then SD would be NaN: try again
    bD[b] = (bsy - bsx)/bsz
    b += 1
  end
  return Statistics.std(bD)
end

"""
    reformat_fasta_linear(infile, outfile; sequencetype=LongDNA{4})

Reformat from interleaved (multiple lines and/or multiple records) to sequential
format and write to `outputfile`. Error if `outputfile` already exists.
Return nothing.

# example

```julia
julia> file = "examples/interleaved.fasta";

julia> reformat_fasta_linear(file, "output.fasta")

shell> diff output.fasta examples/sequential.fasta

```
"""
function reformat_fasta_linear(infile::AbstractString, outfile::AbstractString;
                               sequencetype=LongDNA{4})
  ispath(outfile) && error("I don't want to overwrite $outfile")
  seqs = Array{BioSequences.BioSequence}(undef, 0)
  taxa = String[]
  reader = FASTA.Reader(open(infile, "r"))
  for record in reader
    seq1 = sequence(sequencetype, record)
    push!(seqs, seq1)
    push!(taxa, FASTA.identifier(record))
  end
  close(reader)
  taxaunique = unique(taxa)
  writer = FASTA.Writer(open(outfile, "w"), width=-1)
  for taxon in taxaunique
    seq1 = sequencetype() # initialize with empty sequence
    idx = findall(isequal(taxon), taxa) # identify all sequences for this taxon
    for i in idx seq1 *= seqs[i]; end   # concatenate sequences
    record = FASTA.Record(taxon, seq1)
    write(writer, record)
  end
  close(writer)
  return nothing
end

"""
    calcD_locuscount(x, y)

Take a vector `x` of ABBA counts and a vector `y` of BABA counts, assuming that
in each list: 1 count corresponds to 1 gene. Then calculate the a statistic
similar to D, but that counts a number of loci instead of a number of sites.
Namely: each locus is labelled as an
- ABBA-biased locus if it has more ABBA sites than BABA sites
- BABA-biased locus if it has more BABA sites than ABBA sites.
Then we define:
nABBA = number of ABBA-biased loci,
nBABA = number of BABA-biased loci, and
locus-based D = [nABBA - nBABA] / [nABBA + nBABA]

Note that a locus with an equal number of ABBA sites and BABA sites
(or even 0 ABBA and BABA sites) is not biased in either direction,
so it does not enter the counts above.

If we assume that the loci are independent (unlinked with independent gene trees)
then we can carry out a binomial test to see if the locus-based D statistic
is significantly different from 0, that is, if the proportions of ABBA-biased
and BABA-biased loci is different from 0.50.

Output: tuple of (nABBA_loci, nBABA_loci, locus-based D, p-value)

# example

```julia
julia> file = "examples/tinyexample.fasta";

julia> abba_sites = calcPatternCounts(file, ["P2","P1","P3","P5"], [4,3]) # 2 loci
([3, 3], [1, 0], [0, 2], [0, 0])

julia> calcD_locuscount(abba_sites[2], abba_sites[3])
(1, 1, 0.0, 1.0)

julia> calcD_locuscount([5,4,6,5,3, 1], [3,0,0,2,1, 1]) # 6 loci, only 5 contribute
(5, 0, 1.0, 0.0625)


```
"""
function calcD_locuscount(x,y)
  nloci = length(x)
  length(y) == nloci || error("different loci for ABBA and BABA counts")
  nABBAbiased = sum(x .> y)
  nBABAbiased = sum(y .> x)
  nbiased = nABBAbiased + nBABAbiased
  d = (nABBAbiased - nBABAbiased)/nbiased
  pval = 2 * cdf(Binomial(nbiased, 0.5), min(nABBAbiased, nBABAbiased))
  # pval could be >1 if nABBA = nBABA: probability counted twice for that value
  return nABBAbiased, nBABAbiased, d, min(pval,1.0)
end
