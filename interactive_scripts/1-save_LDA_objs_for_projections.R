# Jake Yeung
# Date of Creation: 2021-04-15
# File: ~/projects/scChIX/interactive_scripts/1-save_LDA_objs_for_projections.R
#


hubprefix <- "/home/jyeung/hub_oudenaarden"

inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/from_cluster/LDA_outputs/ldaAnalysisBins_TopBins_autosomesOnly/lda_outputs.count_mat.K27m3.KeepTopBins_500.KeepAllPlates.K-30.binarize.FALSE/ldaOut.count_mat.K27m3.KeepTopBins_500.KeepAllPlates.K-30.Robj"
load(inf, v=T)
out.lda.k27me3 <- out.lda
tm.k27me3 <- posterior(out.lda.k27me3)

inf2 <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/from_cluster/LDA_outputs/ldaAnalysisBins_TopBins_autosomesOnly/lda_outputs.count_mat.K9m3.KeepTopBins_500.KeepAllPlates.K-30.binarize.FALSE/ldaOut.count_mat.K9m3.KeepTopBins_500.KeepAllPlates.K-30.Robj"
load(inf2, v=T)
out.lda.k9me3 <- out.lda
tm.k9me3 <- posterior(out.lda.k9me3)


outf1 <- "/home/jyeung/projects/scChIX/data/H3K27me3_LdaOutputs.RData"
outf2 <- "/home/jyeung/projects/scChIX/data/H3K9me3_LdaOutputs.RData"
# outf2 <- "/home/jyeung/projects/scChIX/data/SingleIncubLdaOutputsPosteriors.RData"

save(out.lda.k27me3, file = outf1)
save(out.lda.k9me3, file = outf2)
# save(tm.k27me3, out.lda.k9me3, file = outf2)

