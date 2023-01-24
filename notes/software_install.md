# Software installation notes

Organized by software package.

## HyDe - done on Darwin

```bash
pip install cython numpy matplotlib seaborn multiprocess
pip install phyde
```

Have to edit my local bashrc ( ~/.bashrc.local) to include this line:
`export PATH=$PATH:/u/l/e/lefrankel/.local/bin`
(or `export PATH=$PATH:~/.local/bin` for any user).
Then can call HyDe with `run_hyde.py` in any directory.
got `phyde-0.4.3` on 2022-05-16.

## PhyloNetworks - done on Franklin

Julia 1.8 available on Darwin as of 2022-9-20.

1. Set up Franklin Julia environment by transferring local Manifest.toml to `~/.julia/environments/v1.8/`
2. Activate julia from `simulation-ratevariation` directory with `julia --project` to download packages from Project.toml like so:

```julia package mode
activate .
instantiate
```

## SiPhyNetwork - done on Franklin

### downloading 

R version 4.2.0 (2022-04-22) -- "Vigorous Calisthenics"
installed on the cluster, default in `/usr/bin/R` as of 2022-05-16.

Install packages within R, locally to `~/R/x86_64-pc-linux-gnu-library/4.2`:

```r
install.packages("devtools") # takes forever!!
install.packages('ape')
devtools::install_github("jjustison/SiPhyNetwork")
```

### updating 

```r
library(remotes)
remotes::update_packages()
```
