# Jake Yeung
# Date of Creation: 2021-04-15
# File: ~/projects/scChIX/interactive_scripts/1-save_projection_objs_all_cells.R
#

rm(list=ls())

jsuffix <- ".FewerTopics2"
jprefix.base <- "TopBins_autosomesOnly"
jprefix <- paste0(jprefix.base, jsuffix, "_KeepAllPlates_clstrbytopics")
jprefix.input <- jprefix
jprefix.output <- paste0(jprefix, "_0_1")
jprefix.ldaoutput <- paste0(jprefix, "_0_1_RemoveBadMixings")

jmarks <- c("K27m3", "K9m3")
names(jmarks) <- jmarks
indir.lda <-  paste0("/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/from_cluster/projections_LDA_outputs/afterUnmixing_", jprefix.ldaoutput, "_K27m-K9m3")

out.proj.lst <- lapply(jmarks, function(jmark){
  inf.lda <- file.path(indir.lda, paste0("project_unmixed_", jmark, ".RData"))
  load(inf.lda)
  return(list(out.lda = out.objs$out.lda, out.lda.predict = out.lda.predict, count.mat.proj = count.mat.proj))
})

outf1 <- "/home/jyeung/projects/scChIX/data/H3K27me3_ProjectionOutput.RData"
outf2 <- "/home/jyeung/projects/scChIX/data/H3K9me3_ProjectionOutput.RData"

out.proj.K27m3 <- out.proj.lst$K27m3
out.proj.K9m3 <- out.proj.lst$K9m3
save(out.proj.K27m3, file = outf1)
save(out.proj.K9m3, file = outf2)

