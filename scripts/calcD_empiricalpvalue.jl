#=
Alternate functions to calculate D, Dp, D3 and D3v2, and the corresponding tests.
There were NOT used in our pipeline, but could be of interest to us or others later.

Compared to functions in `calcD.jl`:
- The code is better here (applies DRY principle: don't repeat yourself),
  but has a different interface: functions take different arguments
- Option to calculate an empirical p-value, in addition to the normal
  approximation used in the functions in `calcD.jl`.

The normal approximation consists of estimating the SE of D (or D3 etc.)
by bootstrapping loci, then getting a z-score z=D/SE(D) then getting a
p-value by comparing the z-score to the standard normal distribution.

The extra option here calculates an empirical one-sided p-value as the
number of bootstrap samples with a D statistic as extreme as or more extreme
than the D statistic from the original data.
The output p-value is twice the one-sided p-value.
The downside of this empirical p-value is that it requires a very large
bootstrap sample to accurately estimate small p-values. And since bootstrap
is locus-based, it requires a large number of loci also.
=#

"""
    calcD(data, D_statistic_type::Symbol)

Take count data and the type of desired D statistic: either `:D`, `:Dp`, `D3`
or `D3v2`. The count data needs to be:
(ABBA,BABA) counts for the D statistic,
(ABBA,BABA,BBAA) counts for the Dp statistic,
(d23,d13) counts for D3, or
(d23,d13,d12) counts for the D3v2 statistic.
Each one should be a list: 1 count corresponds to 1 gene.
Then output is the D statistics for
the concatenated data, that is:

    D    = [sum(ABBA) - sum(BABA)] / [sum(ABBA) + sum(BABA)]
    Dp   = [sum(ABBA) - sum(BABA)] / [sum(ABBA) + sum(BABA) + sum(BBAA)]
    D3   = [sum(d23) - sum(d13)] / [sum(d23) + sum(d13)]
    D3v2 = [sum(d23) - sum(d13)] /  sum(d12)

# examples

```julia
julia> calcD([20, 5], :D)
0.6

julia> calcD(([5,4,6,5], [3,0,0,2]), :D)
0.6

julia> file = "examples/tinyexample.fasta";

julia> pattern_counts = calcPatternCounts(file, ["P2","P1","P3","P5"], [4,3]) # 2 loci
([3, 3], [1, 0], [0, 2], [0, 0], [1, 1], [0, 3], [1, 2])

julia> d_value = calcD(pattern_counts[2:3], :D)
-0.3333333333333333

julia> dp_value = calcD(pattern_counts[2:4], :Dp)
-0.3333333333333333

julia> d3_value = calcD(pattern_counts[5:6], :D3)
-0.2

julia> d3v2_value = calcD(pattern_counts[5:7], :D3v2)
-0.3333333333333333
```
"""
function calcD(data, stattype::Symbol)
  x = data[1]
  y = data[2]
  z = (stattype in [:Dp,:D3v2] ? data[3] : nothing)
  sx = sum(x); sy = sum(y)
  denom = ( stattype in [:D, :D3] ? sx + sy :
           (stattype == :Dp ? sx + sy + sum(z) : sum(z)))
  return (sx - sy)/denom
end

"""
    calcD_locusbootstrap(data, stattype, nbootstrap=5000)
    bootstrap_sd_pval(original_D, bootstrap_D)
    bootstrap_empirical_pval(original_D, bootstrap_D)

`calcD_locusbootstrap` takes a list of counts of type ABBA-BABA or d23-d13,
like `calcD` does, then returns a bootstrap sample of the desired statistic type
based on resampling loci. The sites are *not* resampled: only the loci are.
This is assuming independent loci, but possibly linked sites within loci.

There is one case of concern: if all resampled loci have no informative
(ABBA & BABA; or d13 & d23) sites.
In this case, the D statistic is undefined: 0-0/(0+0)
If that happens, the corresponding bootstrap replicate is dropped.
`nbootstrap` is the total # of valid bootstrap replicates (not counting
dropped reps).

`bootstrap_sd_pval` and `bootstrap_empirical_pval` take the statistic (D or any
other statistic, really) from the original data,
a bootstrap sample of the same test statistic,
and perform a test of the null hypothesis that the population D is 0.
`bootstrap_sd_pval` assumes that the bootstrap sample has a normal distribution.
`bootstrap_empirical_pval` doesn't make this assumption, but requires a very
large bootstrap sample to accurately estimate small p-values.

Warning: the bootstrap procedure requires a large sample size to be valid.
Here, this means a large number of loci. In the example below, we only have
4 loci, so the results from this bootstrap cannot be trusted. But 4 loci
makes the example run fast.

# examples with Patterson's D

```julia
julia> abba_4loci = [5,4,6,5]; baba_4loci = [3,0,0,2];

julia> d_value = calcD((abba_4loci, baba_4loci), :D)
0.6

julia> d_bootstrap = calcD_locusbootstrap((abba_4loci, baba_4loci), :D, 10_000);

julia> bootstrap_sd_pval(d_value, d_bootstrap)
0.00047227005488971745

julia> bootstrap_empirical_pval(dvalue, dboot)
0.002

julia> abba_4loci = [11,0,0,9]; baba_4loci = [3,0,0,2];

julia> d_value = calcD((abba_4loci, baba_4loci), :D)
0.6

julia> d_bootstrap = calcD_locusbootstrap((abba_4loci, baba_4loci), :D, 10_000);

julia> bootstrap_sd_pval(d_value, d_bootstrap)
2.399890083011516e-130

julia> calcD_locusbootstrap(([1,2], [1,2,3], 10), :D)
ERROR: different numbers of loci for ABBA and BABA counts

julia> calcD_locusbootstrap(([0,-1,0], [0,0,0], 10), :D)
ERROR: some ABBA counts are negative

julia> calcD_locusbootstrap(([0,0,0], [0,0,0], 10), :D)
ERROR: all ABBA and BABA are 0
```

# examples using a file, with Hahn & Hibbins's D3

```julia
julia> file = "examples/tinyexample.fasta";

julia> pattern_counts = calcPatternCounts(file, ["P2","P1","P3","P5"], [4,3]) # 2 loci

julia> d3_value = calcD(pattern_counts[5:6], :D3)
-0.2

julia> d3_boot = calcD_locusbootstrap(pattern_counts[5:6], :D3, 10_000);

julia> bootstrap_sd_pval(d3_value, d3_boot)
0.744300905776986

julia> bootstrap_empirical_pval(d3_value, d3_boot)
0.4904

julia> # should be centered at -0.2.
       # only 3 bars: because 2 loci only --the bootstrap is unreliable then!
       using StatsPlots; histogram(d3_boot)
```
"""
function calcD_locusbootstrap(data, stattype::Symbol, nbootstrap=5000)
  x = data[1]
  y = data[2]
  z = (stattype in [:Dp,:D3v2] ? data[3] : nothing)
  nloci = length(x)
  msg = (stattype in [:D,:Dp] ? "ABBA and BABA" : "d23 and d13")
  length(y) == nloci || error("different numbers of loci for $msg counts")
  msg = (stattype in [:D,:Dp] ? "ABBA" : "d13")
  all(x .>= 0) || error("some $msg counts are negative")
  msg = (stattype in [:D,:Dp] ? "BABA" : "d23")
  all(y .>= 0) || error("some $msg counts are negative")
  if stattype == :Dp
    length(z) == nloci || error("different numbers of loci for ABBA and BBAA counts")
    all(z .>= 0) || error("some BBAA counts are negative")
    any(x .> 0) || any(y .> 0) || any(z .> 0) || error("all ABBA BABA & BBAA are 0")
  elseif stattype == :D3v2
    length(z) == nloci || error("different loci for d12 and d13 counts")
    all(z .>= 0) || error("some d12 counts are negative")
    any(z .> 0)  || error("all d12 are 0")
  else
    msg = (stattype == :D ? "ABBA and BABA" : "d13 & d23")
    any(x .> 0) || any(y .> 0) || error("all $msg are 0")
  end
  bD = Vector{Float64}(undef, nbootstrap)
  b = 1
  while b <= nbootstrap
    bsx=0; bsy=0; bsz=0;
    bloci = rand(1:nloci, nloci) # sample loci with replacement
    for locus in bloci
      bsx += x[locus]
      bsy += y[locus]
      if stattype in [:Dp, :D3v2]
        bsz += z[locus]
      end
    end
    sum_denom = bsx + bsy
    if stattype == :Dp
      sum_denom += bsz
    end
    if stattype == :D3v2
      sum_denom = bsz
    end
    sum_denom ≈ 0 && continue # D would be NaN, then SD would be NaN. Try again, don't increment b
    sum_numerator = bsx - bsy
    bD[b] = sum_numerator/sum_denom
    b += 1
  end
  return bD
end


function bootstrap_sd_pval(D_original, D_bootstrap)
  D_std = std(D_bootstrap, corrected=false, mean=D_original)
  zvalue = D_original / D_std
  pvalue = ccdf(Normal(), abs(zvalue))*2
  return pvalue
end
function bootstrap_empirical_pval(D_original, D_bootstrap)
  p_onesided = (D_original > 0.0 ? mean(D_bootstrap .<= 0) : mean(D_bootstrap .>= 0))
  return 2 * p_onesided
end
