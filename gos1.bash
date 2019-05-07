#!/usr/bin/env bash

# download and unpack mother v1.41.1
bash code/get_mothur_linux_64_1_41_1.bash

# get references
# generate a customized version of the SILVA reference database that targets the V4 region
# using Silva.seed_v128 and Trainset16_022016
bash code/silva.seed.align.bash

# run mothur through quality control steps
code/mothur/mothur code/get_good_seqs.batch

# calculate pcr/sequencing error rate
code/mothur/mothur code/get_error.batch

# get shared otus
code/mothur/mothur code/get_shared_otus.batch

# draw figure 1 indicating:
# growth rate (performance) 
# slaughter weight (35 d)
draw_fig1_performance.R

# get alpha diversity
code/mothur/mothur "#summary.single(shared=data/mothur/gos1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.shared, calc=invsimpson-chao-shannon, subsample=T)"

# get beta diversity
code/mothur/mothur code/get_betadiversity.batch

# analyze beta diversity
code/mothur/mothur code/get_differential_OTUs.batch

# draw figure S1 
# inverse Simpson plot
# Bray-Curtis ordination
 code/draw_figS1_diversity.R

# draw figure 2 indicating:
# stacked top 10 OTUs
# LEfSE plots
# absolute abundance plot
draw_fig2_differential_OTUs.R

# get rep fasta
code/mothur/mothur code/get_repfasta.batch

# draw figure 3 indicating:
# cecal IL10, IL17A, IL17F
# ileal IL10, IL17A, IL17F
draw_fig3_cytokine_exp_log.R

# draw figure S2 indicating:
# relationship of ileal and cecal cytokines
# IL17A vs IL10
# IL17F vs IL10
# IL17F vs IL17A
draw_figS4_cytokine_corr.R

# draw figure 3:
# histology
draw_figS4_histo.R

# draw figure 4
# system
draw_fig4_mass_cor.R

# draw figure 5
draw_fig5_ci.R

# draw rarefaction curves
draw_figS5_rarefaction.R



