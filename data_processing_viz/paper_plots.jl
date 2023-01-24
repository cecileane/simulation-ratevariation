using Distributions
using CairoMakie # best for publication-quality vector graphics
using PhyloNetworks
using PhyloPlots
using RCall

const PP = PhyloPlots
const PN = PhyloNetworks

#=
Log-normal distribution used to simulate rate variation
=#

## copy-pasted from def in scripts/simulatedata.jl
function lognormal_meanone(sigma)
    mu = -sigma^2 / 2
    return LogNormal(mu, sigma)
end
# using sigma=0.7
ratedist = lognormal_meanone(0.7);
# find quantiles, to show on the plot: 0.32, 1.92
percentile = round.(quantile(ratedist, [0.10,0.90]), digits=2)

size_inches = (4.5,3) # width,height, in inches. 72 pts per inch
fig = Figure(resolution = 72 .* size_inches, fontsize = 12)
ax = Axis(fig[1, 1],
  xlabel = "rate R\nmean 1, percentiles: $(percentile[1]) (10%) - $(percentile[2]) (90%)",
  ylabel = "density",
  xticks = [0,0.5,1,2,3,4,5])
xlims!(ax, 0,4)
ylims!(ax, 0,1)
hideydecorations!(ax, label=false)
plot!(ratedist)

fig
save("fig_lognormal.pdf", fig)

##

#=
samples of simulated networks to illustrate:

1. h=1 or higher originally but h=0 after shrinking 2 and 3 cycles (H1>0, H2=0),
   without any rate variation, and for which all tests are non significant
   yet gamma was large and there were enough SNPs:
   - rep 44: t16 t3 t7 t15 (case H1=1, H2=0)
   - rep 46: t8  t4 t1 t10 (same type of case)

2. reticulations may cancel each other, such that increasing H2 may decrease
   the power to detect the presence of reticulations.
   example simulations without rate variation, H1=H2, gamma's are ~ 0.4, enough SNPs
   and where D and D3 are significant when h=1 and not significant when h=2.
   - rep 96: t3 t7 t4  t10 (h=1, detected by both D and D3)
   - rep 46: t8 t4 t11 t10 (h=2, not detected by either D and D3, but similar gamma ~ 0.4)

see `resultviz.Rmd` for finding the choice of these replicates & taxon sets.
=#
function quarnet(net::HybridNetwork, fourtaxa::Vector)
    net = deepcopy(net)
    for tip in tipLabels(net)
        if tip ∉ fourtaxa # then prune it from the copied network
            # options simplify and nofuse: to show 2-cycles and nodes from the larger network
            deleteleaf!(net, tip, simplify=false, nofuse=true)
        end
    end
    return net
end

# networks copy-pasted from
# franklin00:/nobackup2/lefrankel/simulation-ratevariation/output_nov21
# from simid=1 (no rate variation) which corresponds to folder
# 10-taxa_1000-genes_1000-sites_0.0-lindist_0.0-linbygenedist_0.0-genedist_0.0-edgedist_0.03-subrate_1.0-coalmedian/
# from file xxx/scalednetwork (xxx = 044 for rep 44) to have
# branch lengths in coalescent units: after SiPhyNetwork + rescaling
net44 = readTopology("((t5:5.103767650391385,(t14:3.435796700756024,((((t2:1.073283269955339,(t3:1.073283269955339)#H25:0.0::0.5959050689)I1:1.2880300656445158,(t11:2.361313335037435)#H23:0.0::0.9108372462)I2:0.08302791918033138,(t7:1.073283269955339,#H25:0.0::0.4040949311)I3:1.3710579846561213)I4:0.03526865005861225,(t15:2.361313335037435,#H23:0.0::0.08916275384)I5:0.11829656923894363)I6:0.926716730044661)I7:1.6532471530280382)I8:4.701077427607387,((t16:0.32346583717508987,t10:0.32346583717508987)I9:0.4679149107541766,t9:0.7913807477605406)I10:9.01346432855097)I11;") 
net46 = readTopology("((((((t4:0.52091928919845,t7:0.52091928919845)I1:1.4569035794981082)#H22:0.0::0.5587859276,#H25:0.08609266787695974::0.4115165525)I2:1.519570161497752,(t5:2.0377313279863976,(((t11:1.8917302007232017)#H25:0.05496130929393801::0.5884834475,(t1:0.7234156477535107,t13:0.7234156477535107)I3:1.2232758619744384)I4:0.031131358602301117,#H22:0.0::0.4412140724)I5:0.059908459820022346)I6:1.4596617017317117)I7:0.7767241380255616,t10:4.274117167737887)I8:3.8367775105672504,((t8:1.8006292489827749,t12:1.8006292489827749)I9:6.306975070287122,t6:8.10760431734196)I10:0.003290360715438531)I11;")
net96 = readTopology("((((t7:0.987699501577265,#H20:0.0::0.4359759431)I1:1.0,(((t4:0.000985080224584983,t1:0.000985080224584983)I2:0.3439639250514807,t2:0.344949005273544)I3:0.64275049605155,(t10:0.987699501577265)#H20:0.0::0.5640240569)I4:1.0)I5:1.4428663726577315,(t11:1.7948360111289163,t9:1.7948360111289163)I6:1.6357298631060806)I7:4.719566120909923,((t3:3.1497211447134297,t8:3.1497211447134297)I8:0.07121744721442134,t6:3.2209385915916227)I9:4.929193400611302)I10;")
net90 = readTopology("(((((t11:2.919937788086447,((t10:1.8200901613505132,#H26:0.0::0.1745316785)I1:0.9483818907159166,#H23:0.09329367604890394::0.3025135972)I2:0.1514657356017557)I3:0.15525863154026753,#H18:0.0::0.4039336612)I4:2.4259296277491984,t8:5.50112604768961)I5:1.0,(((t7:0.6560808885353053,t9:0.6560808885353053)I6:2.419115530986844)#H18:1.745465911994403::0.5960663388,(((t1:1.8200901613505132,(t5:1.8200901613505132)#H26:0.0::0.8254683215)I7:0.855088214656556)#H23:0.0::0.6974864028,t4:2.6751783765298973)I8:2.1454839554049174)I9:1.6804637162776233)I10:2.558973813536396,(t2:0.5408730073129575,t3:0.5408730073129575)I11:8.519226853703918)I12;")
# same networks, but after scaling to get substitutions / site (before seq-gen)
for i in [12,14,6] PN.rotate!(net44, i); end
for i in [13,14]   PN.rotate!(net46, i); end
for i in [11,3]    PN.rotate!(net96, i); end
for i in [21,9,7,6,4, 20,19] PN.rotate!(net90, i); end

net44_q1 = quarnet(net44, ["t16","t3","t7","t15"]);
net46_q1 = quarnet(net46, ["t8","t4","t1","t10"]);
net46_q2 = quarnet(net46, ["t8","t4","t11","t10"]);
net96_q1 = quarnet(net96, ["t3","t7","t4","t10"]);
net90_q1 = quarnet(net90, ["t2","t7","t1","t5"]);
net90_q2 = quarnet(net90, ["t2","t8","t7","t5"])
net90_q2shrunk = PP.removedegree2nodes!(deepcopy(net90_q2), true)
shrink3cycles!(net90_q2shrunk)
PN.rotate!(net90_q2shrunk, 20);
net90_q3 = quarnet(net90, ["t2","t10","t7","t1"]);

#= --- trying a few plots first, to rotate edges for nicer final plots ---#
PP.plot(net44, shownodenumber=true, tipoffset=0.1);
PP.plot(net44_q1; showgamma=true, useedgelength=true)

res = PP.plot(net46, tipoffset=0.1, shownodenumber=true, showedgenumber=true);
res[14][24,:] # to find x and y for edge 24 to outgroups
PP.plot(net46_q1; showgamma=true, useedgelength=true)

res = PP.plot(net96, tipoffset=0.1, shownodenumber=true, showedgenumber=true);
res[14][21,:] # to find x and y for edge 21 to outgroups

res = PP.plot(net90, tipoffset=0.1, shownodenumber=true, showedgenumber=true);
res[14][27,:] # to find x and y for edge 24 to outgroups: 5, 12.5
shrink3cycles!(PP.removedegree2nodes!(deepcopy(net90_q1), true)) # false: h=h'
PP.plot(net90_q1; useedgelength=true, tipoffset=0.1);

PP.plot(net90_q2; useedgelength=true, tipoffset=0.1);
PP.plot(net90_q2shrunk; useedgelength=true, tipoffset=0.1);
shrink3cycles!(PP.removedegree2nodes!(deepcopy(net90_q3), true)) # false: h=h'
PP.plot(net90_q3; useedgelength=true, tipoffset=0.1);
=#

## ------ preliminary plots -----------------------#
# but not retained: reticulations between sister species
R"pdf"("fig_illustrativenets.pdf", height=6, width=6)
R"par"(mar=[.1,.1,.1,.1]); R"layout"([1 2; 1 3; 1 4]);
PP.plot(net46; showgamma=true, tipoffset=0.1);
R"text"(x=3.5, y=11.25+0.1, "outgroups");
R"mtext"("a)", side=3, line=-3, adj=0.1);
PP.plot(net46_q1; useedgelength=true, tipoffset=0.1);
R"mtext"("b)", side=3, line=-3, adj=0.01);
R"mtext"("h=1, h'=0", side=1, line=-5, adj=0.2);
PP.plot(net46_q2; useedgelength=true, tipoffset=0.1);
R"mtext"("c)", side=3, line=-3, adj=0.01);
R"mtext"("h=h'=2", side=1, line=-5, adj=0.2);
# PP.plot(net96; showgamma=true, tipoffset=0.1);
# R"text"(x=3, y=10.25+0.2, "outgroups");
PP.plot(net96_q1; useedgelength=true, tipoffset=0.1);
R"mtext"("d)", side=3, line=-3, adj=0.01)
R"mtext"("h=h'=1", side=1, line=-5, adj=0.2);
R"dev.off"()

## ------ plots for final figures -----------------------#
#         reticulations that hide each other
hcol = "darkorange3" # "deepskyblue3" # color for hybrid edge, major
mcol = "darkorange1" # color for minor hybrid edge
R"pdf"("fig_illustrativenets.pdf", height=5, width=8.5)
R"par"(mar=[.1,.1,.1,.1]);
# 7th plot: empty plot in the middle, for space between plots 1-3 and 4-6
R"layout"([1 2 7 4 5; 1 3 7 4 6], widths=[1,1,0.25,1,1]);
PP.plot(net46; showgamma=true, tipoffset=0.1,
        majorhybridedgecolor=hcol, minorhybridedgecolor=mcol,);
R"text"(x=3.5, y=11.25+0.15, "outgroups");
R"mtext"("a)", side=3, line=-3, adj=0.1);
PP.plot(net46_q1; useedgelength=true, tipoffset=0.1,
        majorhybridedgecolor=hcol, minorhybridedgecolor=mcol,);
R"mtext"("b)", side=3, line=-3, adj=0.01);
R"mtext"("h =1", side=1, line=-6.5, adj=0.15, col=hcol);
R"mtext"("h'=0", side=1, line=-5  , adj=0.15);
PP.plot(net46_q2; useedgelength=true, tipoffset=0.1,
        majorhybridedgecolor=hcol, minorhybridedgecolor=mcol,);
R"mtext"("c)", side=3, line=-3, adj=0.01);
R"mtext"("h=h'=2", side=1, line=-6, adj=0.2, col=hcol);
PP.plot(net90; showgamma=true, tipoffset=0.1,
        majorhybridedgecolor=hcol, minorhybridedgecolor=mcol,);
R"text"(x=5, y=12.5+0.2, "outgroups");
R"mtext"("d)", side=3, line=-3, adj=0.1);
PP.plot(net90_q3; useedgelength=true, tipoffset=0.1,
        majorhybridedgecolor=hcol, minorhybridedgecolor=mcol,);
R"mtext"("e)", side=3, line=-3, adj=0.01)
R"mtext"("h =3", side=1, line=-5  , adj=0.15);
R"mtext"("h'=2", side=1, line=-3.5, adj=0.15, col=hcol);
PP.plot(net90_q2shrunk; useedgelength=true, tipoffset=0.1,
        majorhybridedgecolor=hcol, minorhybridedgecolor=mcol,);
R"mtext"("f)", side=3, line=-3, adj=0.01)
R"mtext"("h =3", side=1, line=-6.5, adj=0.15);
R"mtext"("h'=2", side=1, line=-5  , adj=0.15, col=hcol);
R"dev.off"()
