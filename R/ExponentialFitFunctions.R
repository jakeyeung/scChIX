
IsUpOrDown <- function(dat.impute.lst, refmark, jgene, dat.meta.long){
  
  dat.refmark <- data.frame(signal = dat.impute.lst[[refmark]][jgene, ], 
                            cell = colnames(dat.impute.lst[[refmark]]), 
                            mark = rep(refmark, ncol(dat.impute.lst[[refmark]])), 
                            stringsAsFactors = FALSE)
  
  dat.refmark <- left_join(dat.refmark, dat.meta.long, by = c("cell", "mark")) %>%
    rowwise() %>%
    mutate(type = ifelse(grepl(pattern = jmark.dbl, x = cell), "dbl", "single"), 
           gene = jgene) %>%
    arrange(ptime.exp.norm)
  
  fit.lm <- lm(formula = signal ~ ptime.exp.norm, data = dat.refmark)
  slope <- summary(fit.lm)$coefficients["ptime.exp.norm", "Estimate"]
  # MaxOrMin <- ifelse(slope > 0, min, max)
  UpOrDown <- ifelse(slope > 0, "up", "down")
  return(UpOrDown)
}

expo.relax <- function(ptime.exp.norm, mugamma, gamma, s0){
  s0 + mugamma * (1 - exp(-gamma * ptime.exp.norm))
}


SetupDatForGene.UpOrDown <- function(dat.impute.lst, dat.meta.long, jgene, jmark1 = "H3K4me1", jmark2 = "H3K36me3", refmark = "H3K36me3", returnlst = TRUE, offset.zero = TRUE, frac.cells.filter = 0.01, MaxOrMin = max){
  jmark.dbl <- paste(jmark1, jmark2, collapse = "-")
  dat.exprs.long <- data.frame(signal = c(dat.impute.lst[[jmark1]][jgene, ], dat.impute.lst[[jmark2]][jgene, ]), 
                               cell = c(colnames(dat.impute.lst[[jmark1]]), colnames(dat.impute.lst[[jmark2]])), 
                               mark = c(rep(jmark1, ncol(dat.impute.lst[[jmark1]])), rep(jmark2, ncol(dat.impute.lst[[jmark2]]))), 
                               stringsAsFactors = FALSE)
  dat.exprs.long.annot <- left_join(dat.exprs.long, dat.meta.long, by = c("cell", "mark")) %>%
    rowwise() %>%
    mutate(type = ifelse(grepl(pattern = jmark.dbl, x = cell), "dbl", "single"), 
           gene = jgene) %>%
    arrange(ptime.exp.norm)
  
  # ggplot(dat.exprs.long.annot, aes(x = ptime.exp.norm, y = signal, color = mark)) + 
  #   geom_point()   + 
  #   theme_bw() + 
  #   theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  if (offset.zero){
    # check up or down
    # take average of top 1% of pseudotime as offset, this handles Up or Down
    # joffset <- min(dat.exprs.long.annot$signal)  # only works if up
    dat.refmark <- subset(dat.exprs.long.annot, mark == refmark) %>%
      arrange(ptime.exp.norm)
    
    
    
    # print(MaxOrMin)
    
    # ncells <- nrow(dat.refmark)
    # ncells.filter <- round(frac.cells.filter * ncells) # number of cells to keep
    # cells.select <- (dat.refmark %>% arrange(ptime.exp.norm))$cell[1:ncells.filter]
    # joffset <- mean(subset(dat.refmark, cell %in% cells.select)$signal)
    
    joffset <- MaxOrMin(dat.refmark$signal)
    
    dat.exprs.long.annot <- dat.exprs.long.annot %>%
      ungroup() %>%
      mutate(signal = signal - joffset)
    
    # jtest <- subset(dat.exprs.long.annot, mark == jmark2)
    # 
    # ggplot(jtest, aes(x = ptime.exp.norm, y = signal, color = mark)) + 
    #   geom_point()   + 
    #   theme_bw() + 
    #   theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    # 
    # asymp.init <- min( jtest )$signal
    # loggamma.init <- 0
    # jfit <- nls(formula = signal ~ SSasympOrig(ptime.exp.norm, asymp, loggamma), 
    #             data = jtest, 
    #             start=list(asymp = asymp.init, loggamma = loggamma.init), 
    #             control = nls.control(maxiter = jiter))
    
  }
  if (returnlst){
    dat.exprs.long.annot <- split(x = dat.exprs.long.annot, f = dat.exprs.long.annot$mark)
    # dat.exprs.long.annot$log <- list("UpOrDown" = UpOrDown, "MaxOrMin" = MaxOrMin)
  }
  return(dat.exprs.long.annot)
}

FitAsympTryCatch <- function(dat.exprs.long.gene, asymp.init, loggamma.init, jiter = 1000){
  jfit <- tryCatch({
    jfit <- nls(formula = signal ~ SSasympOrig(ptime.exp.norm, asymp, loggamma), 
                data = dat.exprs.long.gene, 
                start=list(asymp = asymp.init, loggamma = loggamma.init), 
                control = nls.control(maxiter = jiter))
  }, error = function(e) {
    jfit <- e
  })
  return(jfit)
}


FitAsympMarks.UpOrDown <- function(dat.exprs.long.gene.lst, jmark1 = "H3K4me1", jmark2 = "H3K36me3", asymp.init = "auto", loggamma.init = 0, jiter = 1000, MaxOrMin = max){
  jgene <- unique(dat.exprs.long.gene.lst[[1]]$gene)
  if (asymp.init == "auto"){
    asymp.init <- MaxOrMin( (dat.exprs.long.gene.lst %>% bind_rows() )$signal)
  }
  
  # ggplot(dat.exprs.long.gene.lst[[jmark2]], aes(x = ptime.exp.norm, y = signal - 1)) + 
  #   geom_point() + 
  #   geom_hline(yintercept = 0) + 
  #   theme_bw() + 
  #   theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  jfit1 <- FitAsympTryCatch(dat.exprs.long.gene.lst[[jmark1]], asymp.init = asymp.init, loggamma.init = loggamma.init, jiter = jiter)
  jfit2 <- FitAsympTryCatch(dat.exprs.long.gene.lst[[jmark2]], asymp.init = asymp.init, loggamma.init = loggamma.init, jiter = jiter)
  # assmemble output
  dat.out <- data.table(gene = jgene, fitout = list(jfit1, jfit2), mark = c(jmark1, jmark2), stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(is.error = any(class(fitout) %in% c("simpleError")))
  return(dat.out)
}

FitExpoRelaxMarks <- function(dat.exprs.long.gene.lst, jmark1 = "H3K4me1", jmark2 = "H3K36me3", s0.init = "auto", mugamma.init = "auto", gamma.init = 0){
  jgene <- unique(dat.exprs.long.gene.lst[[1]]$gene)
  if (s0.init == "auto"){
    s0.init <- min( (dat.exprs.long.gene.lst %>% bind_rows() )$signal)
  }
  if (mugamma.init == "auto"){
    mugamma.init <- max( (dat.exprs.long.gene.lst %>% bind_rows() )$signal)
  }
  jfit2 <- nls(formula = signal ~ expo.relax(ptime.exp.norm, mugamma, gamma, s0), 
               data = dat.exprs.long.gene.lst[[jmark2]], 
               start=list(s0 = s0.init, mugamma = mugamma.init, gamma = gamma.init))
  s0.from.fit2 <- coef(jfit2)["s0"]
  jfit1 <- nls(formula = signal ~ expo.relax(ptime.exp.norm, mugamma, gamma, s0 = s0.from.fit2),
               data = dat.exprs.long.gene.lst[[jmark1]], 
               start=list(mugamma = mugamma.init, gamma = gamma.init))
  # assmemble output
  dat.out <- data.table(gene = jgene, fitout = list(jfit1, jfit2), mark = c(jmark1, jmark2), stringsAsFactors = FALSE)
  return(dat.out)
}