#Utility functions
library(latex2exp)
library(ggbeeswarm)
#set the global ggplot theme
theme_full <- theme_bw() + theme(axis.text = element_text(size=14),
                             axis.title = element_text(size=16),
                             axis.line = element_blank(),
                             panel.border = element_rect(size=1.5),
                             axis.ticks = element_line(size=1.5),
                             plot.title = element_text(size = 16, hjust =0.5, face="bold"),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank())

theme_half <- theme_bw() + theme(axis.text = element_text(size=14),
                                 axis.title = element_text(size=16),
                                 axis.line =element_line(size=0.8),
                                 panel.border = element_blank(),
                                 axis.ticks = element_line(size=1.5),
                                 plot.title = element_text(size = 16, hjust =0.5, face="bold"),
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank())

# Defien a color scheme, based on ggsci_NEJM panel, for the paper
colList <- c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF","#FFDC91FF","#EE4C97FF")

# format pathway names
setMap <- read_tsv("../data/setToPathway.txt", col_types = "cc")

# Change to mutational status to more readable ones
formatStatus <- function(patBack) {
  patTab <- mutate_all(patBack, as.character) %>%
    pivot_longer(-Patient.ID) %>%
    mutate(value = case_when(
      str_detect(name, "IGHV|Methylation") ~ value,
      str_detect(name, "gain|del|trisomy") & value %in% "1" ~ "yes",
      str_detect(name, "gain|del|trisomy") & value %in% "0" ~ "no",
      value %in% "1" ~ "Mut",
      value %in% "0" ~ "WT"
    )) %>%
    pivot_wider(names_from = "name", values_from = "value")
}

# Prepare column annotations for heatmap
genAnnoCol <- function(varList) {
  annoCol <- lapply(varList, function(name) {
    if (name == "onChr12") {
      c("yes" = colList[1],"no" = "white")
    } else if (str_detect(name, "IGHV")) {
      c(M = colList[4], U = colList[3])
    } else if (str_detect(name, "gain|del|trisomy")) {
      c(yes = "black",no = "white")
    } else {
      c(Mut = "black",WT = "white")
    }
  })
  names(annoCol) = varList
  return(annoCol)
}


# Function for scale a data matrix
mscale <- function(x, center = TRUE, scale = TRUE, censor = NULL, useMad = FALSE){
  if (scale & center) {
    if (useMad) {
      x.scaled <- apply(x, 1, function(y) (y-median(y))/(1.4826*mad(y)))
    } else {
      x.scaled <- apply(x, 1, function(y) (y-mean(y))/(sd(y)))
    }
  } else if (center & !scale) {
    if (useMad) {
      x.scaled <- apply(x, 1, function(y) (y-median(y)))
    } else {
      x.scaled <- apply(x, 1, function(y) (y-mean(y)))
    }
  } else if (!center & scale) {
    if (useMad) {
      x.scaled <- apply(x, 1, function(y) y/(1.4826*mad(y)))
    } else {
      x.scaled <- apply(x, 1, function(y) y/(sd(y)))
    }
  } else {
    x.scaled <- x
  }

  if (!is.null(censor)) {
    x.scaled[x.scaled > censor] <- censor
    x.scaled[x.scaled < -censor] <- -censor
  }
  return(t(as.matrix(x.scaled)))
}


plotCorScatter <- function(inputTab, x, y, x_lab = "X", y_lab = "Y", title = "",
                           col = NULL, shape = NULL, showR2 = TRUE, annoPos = "right", legendPos = "left",
                           dotCol = colList, dotShape = c(16,1,17,2), textCol="darkred") {

  #prepare table for plotting
  plotTab <- tibble(x = inputTab[[x]],y=inputTab[[y]])
  if (!is.null(col)) plotTab <- mutate(plotTab, status = inputTab[[col]])
  if (!is.null(shape)) plotTab <- mutate(plotTab, statusShape = inputTab[[shape]])

  plotTab <- filter(plotTab, !is.na(x), !is.na(y))

  #prepare annotation values
  corRes <- cor.test(plotTab$x, plotTab$y)
  pval <- formatNum(corRes$p.value, digits = 1, format = "e")
  Rval <- formatNum(corRes$estimate, digits = 2)
  R2val <- formatNum(corRes$estimate^2, digits = 2)
  Nval <- nrow(plotTab)
  annoP <- bquote(italic("P")~"="~.(pval))

  if (showR2) {
    annoCoef <-  bquote(R^2~"="~.(R2val))
  } else {
    annoCoef <- bquote(R~"="~.(Rval))
  }
  annoN <- bquote(N~"="~.(Nval))

  corPlot <- ggplot(plotTab, aes(x = x, y = y))

  if (!is.null(col) & is.null(shape)) {
    corPlot <- corPlot + geom_point(aes(fill = status), shape =21, size =3) +
      scale_fill_manual(values = dotCol)
  } else if (is.null(col) & !is.null(shape)) {
    corPlot <- corPlot + geom_point(aes(shape = statusShape), color = dotCol[1], size=3) +
      scale_shape_manual(values = dotShape)
  } else if (!is.null(col) & !is.null(shape)) {
    corPlot <- corPlot + geom_point(aes(shape = statusShape, color = status), size=3) +
      scale_shape_manual(values = dotShape) +
      scale_color_manual(values = dotCol)
  }
  else {
    corPlot <- corPlot + geom_point(fill = dotCol[1], shape =21, size=3)
  }

  corPlot <- corPlot +   geom_smooth(formula = y~x,method = "lm", se=FALSE, color = "grey50", linetype ="dashed" )

  if (annoPos == "right") {

    corPlot <- corPlot + annotate("text", x = max(plotTab$x), y = Inf, label = annoN,
                                  hjust=1, vjust =2, size = 5, parse = FALSE, col= textCol) +
      annotate("text", x = max(plotTab$x), y = Inf, label = annoP,
               hjust=1, vjust =4, size = 5, parse = FALSE, col= textCol) +
      annotate("text", x = max(plotTab$x), y = Inf, label = annoCoef,
               hjust=1, vjust =6, size = 5, parse = FALSE, col= textCol)

  } else if (annoPos== "left") {
    corPlot <- corPlot + annotate("text", x = min(plotTab$x), y = Inf, label = annoN,
                                  hjust=0, vjust =2, size = 5, parse = FALSE, col= textCol) +
      annotate("text", x = min(plotTab$x), y = Inf, label = annoP,
               hjust=0, vjust =4, size = 5, parse = FALSE, col= textCol) +
      annotate("text", x = min(plotTab$x), y = Inf, label = annoCoef,
               hjust=0, vjust =6, size = 5, parse = FALSE, col= textCol)
  }
  corPlot <- corPlot + ylab(y_lab) + xlab(x_lab) + ggtitle(title) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
    theme_full + theme(legend.position = legendPos,
                       plot.margin = margin(13,13,13,13))
  corPlot
}


#function for cox regression
com <- function(response, time, endpoint, scale =FALSE) {

  if (scale) {
    #calculate z-score
    response <- (response - mean(response, na.rm = TRUE))/sd(response, na.rm=TRUE)
  }
  surv <- coxph(Surv(time, endpoint) ~ response)


  tibble(p = summary(surv)[[7]][,5],
         HR = summary(surv)[[7]][,2],
         lower = summary(surv)[[8]][,3],
         higher = summary(surv)[[8]][,4])
}

# Function for Kaplan-Meier plot
km <- function(response, time, endpoint, titlePlot = "KM plot", pval = NULL,
               stat = "median", maxTime =NULL, showP = TRUE, showTable = FALSE,
               ylab = "Fraction", xlab = "Time (years)",
               table_ratio = c(0.7,0.3), yLabelAdjust = 0) {
  #function for km plot
  survS <- tibble(time = time,
                  endpoint = endpoint)

  if (!is.null(maxTime))
    survS <- mutate(survS, endpoint = ifelse(time > maxTime, FALSE, endpoint),
                    time = ifelse(time > maxTime, maxTime, time))

  if (stat == "maxstat") {
    ms <- maxstat.test(Surv(time, endpoint)  ~ response,
                       data = survS,
                       smethod = "LogRank",
                       minprop = 0.2,
                       maxprop = 0.8,
                       alpha = NULL)

    survS$group <- factor(ifelse(response >= ms$estimate, "high", "low"))
    p <- com(survS$group, survS$time, survS$endpoint)$p

  } else if (stat == "median") {
    med <- median(response, na.rm = TRUE)
    survS$group <- factor(ifelse(response >= med, "high", "low"))
    p <- com(survS$group, survS$time, survS$endpoint)$p

  } else if (stat == "binary") {
    survS$group <- factor(response)
    if (nlevels(survS$group) > 2) {
      sdf <- survdiff(Surv(survS$time,survS$endpoint) ~ survS$group)
      p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
    } else {
      p <- com(survS$group, survS$time, survS$endpoint)$p
    }
  }

  if (is.null(pval)) {
    if(p< 1e-16) {
      pAnno <- bquote(italic("P")~"< 1e-16")
    } else {
      pval <- formatNum(p, digits = 1)
      pAnno <- bquote(italic("P")~"="~.(pval))
    }

  } else {
     pval <- formatNum(pval, digits = 1)
     pAnno <- bquote(italic("P")~"="~.(pval))
  }

  if (!showP) pAnno <- ""

  colorPal <- colList[1:length(unique(survS$group))]
  p <- ggsurvplot(survfit(Surv(time, endpoint) ~ group, data = survS),
                  data = survS, pval = FALSE,  conf.int = FALSE, palette = colorPal,
                  legend = ifelse(showTable, "none","top"),
                  ylab = "Fraction", xlab = "Time (years)", title = titlePlot,
                  pval.coord = c(0,0.1), risk.table = showTable, legend.labs = sort(unique(survS$group)),
                  ggtheme = theme_half + theme(plot.title = element_text(hjust =0.5),
                                               panel.border = element_blank(),
                                               axis.title.y = element_text(vjust =yLabelAdjust)))
  if (!showTable) {
    p <- p$plot + annotate("text",label=pAnno, x = 0.1, y=0.1, hjust =0, size =5)
    return(p)
  } else {
    #construct a gtable
    pp <- p$plot + annotate("text",label=pAnno, x = 0.1, y=0.1, hjust =0, size=5)
    pt <- p$table + ylab("") + xlab("") + theme(plot.title = element_text(hjust=0, size =10))
    p <- plot_grid(pp,pt, rel_heights = table_ratio, nrow =2, align = "v")
    return(p)
  }
}



#Function to run multivariate Cox model on test table
runCox <- function(survTab, riskTab, time, endpoint) {
  survTab <- select(survTab, patientID, !!time, !!endpoint) %>%
    dplyr::rename(time = !!time, endpoint = !!endpoint) %>%
    filter(!is.na(time), !is.na(endpoint))
  testTab <- right_join(survTab, riskTab, by = "patientID") %>%
    select(-patientID)
  surv1 <- coxph(
    Surv(time, endpoint) ~
      .,
    data = testTab)
  return(surv1)
}

#Function to calculate C-index
calcConcord <- function(survTab, riskTab, time, endpoint) {
  survTab <- select(survTab, patientID, !!time, !!endpoint) %>%
    dplyr::rename(time = !!time, endpoint = !!endpoint) %>%
    filter(!is.na(time), !is.na(endpoint))
  testTab <- right_join(survTab, riskTab, by = "patientID") %>%
    select(-patientID)
  surv1 <- concordance(
    Surv(time, endpoint) ~
      .,
    data = testTab)
  return(surv1)
}

#Function to generate pretty scientific notation format for plot label
sciPretty <- function(n, digits = 2, bold = FALSE) {
  nForm <- strsplit(format(n, digits = digits, scientific = TRUE),split = "e")
  b <- nForm[[1]][1]
  i <- as.integer(nForm[[1]][2])
  #bquote(.(b)%*%10^.(i))
  if(bold) {
    sprintf("bold(%s%%*%%10^%s)",b,i)
  } else sprintf("%s%%*%%10^%s",b,i)
}

#function to remove highly correlated values and keep track of the removing
removeCorrelated <- function(x, cutoff = 0.6, method = "pearson", keep = NULL, record = TRUE) {

  if (!is.null(keep)) {
    #if specified some feature to keep, then reorder the matrix, to make sure they are not collapsed
    posNow <- grepl(paste(keep, collapse = "|"), colnames(x))
    posNew <- rev(c(colnames(x)[posNow],colnames(x)[!posNow]))
    x <- x[,posNew]
  }
  #input is a feature matrix, features are in columns
  if (method == "binary") {
    #use binary similarity if input is a binary matrix,
    #maybe also usefull is the input is a sparse matrix
    simiMat <- 1 - as.matrix(dist(t(x), method = "binary"))
  } else if (method == "pearson") {
    #otherwise, using pearson correlation
    simiMat <- cor(x)
  } else if (method == "euclidean") {
    simiMat <- 1 - as.matrix(dist(t(x), method = "euclidean"))
  } else if (method == "cosine") {
    # cosine similarity maybe prefered for sparse matrix
    cosineSimi <- function(x){
      x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2))))
    }
    simiMat <- cosineSimi(t(x))
  } else if (method == "canberra") {
    simiMat <- 1 - as.matrix(dist(t(x), method = "canberra"))/nrow(x)
  }

  #generate reduced matrix
  simiMat.ori <- simiMat
  simiMat[upper.tri(simiMat)] <- 0
  diag(simiMat) <- 0
  x.re <- x[,!apply(simiMat, 2, function(n) any(abs(n) >= cutoff))]

  if (record) {
    #a matrix keeping track of the removed features
    mapReduce <- simiMat.ori
    diag(mapReduce) <- 0
    mapList <- lapply(colnames(x.re), function(i) colnames(mapReduce)[mapReduce[i,]>=cutoff])
    names(mapList) <- colnames(x.re)
  } else mapList = NULL

  return(list(reduced = x.re,
              mapReduce = mapList))
}

# Run glmnet with gaussian familiy outcome
runGlm <- function(X, y, method = "lasso", repeats=20, folds = 3, testRatio = NULL, lambda = "lambda.1se") {
  modelList <- list()
  lambdaList <- c()
  r2Train <- c()
  r2Test <- c()
  coefMat  <- matrix(NA, ncol(X), repeats)
  rownames(coefMat) <- colnames(X)

  if (method == "lasso"){
    alpha = 1
  } else if (method == "ridge") {
    alpha = 0
  }

  for (i in seq(repeats)) {
    if (!is.null(testRatio)) {
      testIdx <- sample(seq_along(y), length(y)*testRatio)
      X.test <- X[testIdx,]
      X.train <- X[-testIdx,]
      y.test <- y[testIdx]
      y.train <- y[-testIdx]
    } else {
      X.train <- X
      y.train <- y
    }

    vecFold <- mltools::folds(y.train,nfolds = folds, stratified = TRUE)

    #train model
    res <- cv.glmnet(X.train,y.train, type.measure = "mse",
                     foldid = vecFold, alpha = alpha, standardize = FALSE,
                     intercept = TRUE, family = "gaussian")
    lambdaList <- c(lambdaList, res[[lambda]])

    #calculate variance explained for training model
    y.pred <- predict(res, s = lambda, newx = X.train)[,1]
    varExp <- cor(y.train,y.pred)^2
    r2Train <- c(r2Train, varExp)

    modelList[[i]] <- res
    coefMat[,i] <- coef(res, s = lambda)[-1]

    #test model if testRatio is speficied
    if(!is.null(testRatio)) {
      y.pred <- predict(res, s = lambda, newx = X.test)
      varExp <- cor(as.vector(y.test),as.vector(y.pred))^2
      r2Test <- c(r2Test, varExp)
    }
  }
  list(modelList = modelList, lambdaList = lambdaList, r2Train = r2Train, coefMat = coefMat,
       r2Test = r2Test)
}


#get number of selected features from a LASSO model
nFeature <- function(lassoRes) {
  coefMat <- lassoRes$coefMat
  return(colSums(coefMat != 0))
}



runCamera <- function(exprMat, design, gmtFile, id = NULL,
                      contrast = ncol(design),  method = "camera", pCut = 0.05,
                      ifFDR = FALSE, removePrefix = NULL, plotTitle = "", insideLegend = FALSE, setMap = NULL) {

  #prepare indices
  if (is.null(id)) id <- rownames(exprMat)

  if (is.character(gmtFile)) {
    idx <- limma::ids2indices(piano::loadGSC(gmtFile)$gsc, id)
  } else {
    idx <- limma::ids2indices(gmtFile,id)
  }

  #run camera for fry
  if (method == "camera") {
    res <- limma::camera(exprMat, idx, design, contrast)
  } else if (method == "fry") {
    res <- limma::fry(exprMat, idx, design, contrast)
  }

  #plot enrichment results as bar plot

  plotTab <- res %>% rownames_to_column("Name")
  if (!is.null(removePrefix)) plotTab <- mutate(plotTab, Name = str_remove(Name, removePrefix))

  if(!is.null(setMap)) {
    plotTab <- mutate(plotTab, newName = setMap[match(Name, setMap$setName),]$pathwayName) %>%
      mutate(Name = ifelse(is.na(newName), Name, newName))
  }

  plotTab <- plotTab %>%
    mutate(Direction= factor(Direction, levels =c("Down","Up"))) %>%
    arrange(desc(Direction),desc(PValue)) %>%
    mutate(Name = factor(Name, levels = Name))

  if (ifFDR) {
    plotTab <- dplyr::filter(plotTab, FDR <= pCut)
  } else {
    plotTab <- dplyr::filter(plotTab, PValue <= pCut)
  }

  if (nrow(plotTab) == 0) {
    print("No sets passed the criteria")
    return(list(enrichTab = res, enrichPlot = NULL))
  } else {
    p <- ggplot(data = plotTab, aes(x = Name, y = -log10(PValue), fill = Direction)) +
      geom_bar(position = "dodge", stat = "identity", width = 0.5) +
      scale_fill_manual(values = c(Up = colList[1], Down = colList[2])) +
      coord_flip() + xlab("") +
      ylab(expression(-log[10]*'('*italic(P)~value*')')) + ggtitle(plotTitle) +
      theme_full +
      theme(axis.text = element_text(size= 12))
    if (insideLegend) {
      p <- p + theme(legend.position = c(0.8,0.1))
    } else {
      p <- p + theme(legend.position = "right")
    }

    return(list(enrichTab = res, enrichPlot = p))
  }
}

#function to run fisher test for enrichment analysis
runFisher <- function (genes, reference, gmtFile, adj = "BH", verbose = FALSE, barCol = NULL, setName = "",
                       pCut = 0.05,ifFDR = FALSE, removePrefix = NULL, plotTitle = "", insideLegend = FALSE, setMap = NULL, topN = NULL) {

  genesets <- piano::loadGSC(gmtFile)$gsc
  tab = lapply(1:length(genesets), function(i) {
    if (verbose == TRUE) {
      cat("processing term", i, names(genesets)[i], "\n")
    }
    reference = reference[!reference %in% genes]
    RinSet = sum(reference %in% genesets[[i]])
    RninSet = length(reference) - RinSet
    GinSet = sum(genes %in% genesets[[i]])
    GninSet = length(genes) - GinSet
    fmat = matrix(c(GinSet, RinSet, GninSet, RninSet), nrow = 2,
                  ncol = 2, byrow = F)
    colnames(fmat) = c("inSet", "ninSet")
    rownames(fmat) = c("genes", "reference")
    fish = fisher.test(fmat, alternative = "greater")
    pval = fish$p.value
    inSet = RinSet + GinSet
    res = c(GinSet, inSet, pval)
    res
  })
  rtab = do.call("rbind", tab)
  rtab = data.frame(as.vector(names(genesets)), rtab)
  rtab = rtab[order(rtab[, 4]), ]
  colnames(rtab) = c("TermID", "genes", "all", "pval")
  padj = p.adjust(rtab[, 4], method = "BH")
  tab.out = data.frame(rtab, padj)

  plotTab <- tab.out %>% dplyr::rename(Name = TermID, PValue = pval, FDR = padj) %>% arrange(PValue)

  if (ifFDR) {
    plotTab <- dplyr::filter(plotTab, FDR <= pCut)
  } else {
    plotTab <- dplyr::filter(plotTab, PValue <= pCut)
  }

  if (nrow(plotTab) == 0) {
    print("No sets passed the criteria")
    return(list(enrichTab = tab.out, enrichPlot = NULL))
  }

  if (!is.null(topN)) {
    plotTab <- plotTab[1:min(nrow(plotTab),topN),]
  }

  if (!is.null(removePrefix)) plotTab <- mutate(plotTab, Name = str_remove(Name, removePrefix))
  if(!is.null(setMap)) {
    plotTab <- mutate(plotTab, newName = setMap[match(Name, setMap$setName),]$pathwayName) %>%
      mutate(Name = ifelse(is.na(newName), Name, newName))
  }

  plotTab <- plotTab %>%
    arrange(desc(PValue)) %>%
    mutate(Name = sprintf("%s(%s)",Name,genes)) %>%
    mutate(Name = factor(Name, levels = Name))



  if (is.null(barCol)) barCol <- colList[1]


  p <- ggplot(data = plotTab, aes(x = Name, y = -log10(PValue))) +
    geom_bar(position = "dodge", stat = "identity", width = 0.5, fill = barCol) +
    coord_flip() + xlab(setName) +
    ylab(expression(-log[10]*'('*italic(P)~value*')')) + ggtitle(plotTitle) +
    theme_full +
    theme(axis.text = element_text(size= 12))
  if (insideLegend) {
    p <- p + theme(legend.position = c(0.8,0.1))
  } else {
    p <- p + theme(legend.position = "right")
  }

  return(list(enrichTab = tab.out, enrichPlot = p))
}

#function to plot heatmap for each gene set (using complex heatmap)
plotSetHeatmap <- function(geneSigTab, setDir, setName, exprMat, colAnno, scale = TRUE, rowAnno = NULL, annoCol = NULL, highLight = NULL, plotName =NULL, italicize=TRUE) {
  geneList <- loadGSC(setDir)[["gsc"]][[setName]]
  sigGene <- dplyr::filter(geneSigTab, name %in% geneList) %>%
    arrange(desc(logFC))
  colAnno <- colAnno[order(colAnno[,1]),,drop = FALSE]
  colAnno <- colAnno[,rev(colnames(colAnno)),drop=FALSE]
  plotMat <- exprMat[sigGene[["name"]],rownames(colAnno)]

  if (scale) {
    #calculate z-score and sensor
    plotMat <- t(scale(t(plotMat)))
    plotMat[plotMat >= 4] <- 4
    plotMat[plotMat <= -4] <- -4
  }

  haCol <- ComplexHeatmap::HeatmapAnnotation(df = colAnno, col=annoCol, which = "column",
                                             annotation_name_gp = gpar(fontface = "bold"))
  if (!is.null(rowAnno)) {
    haRow <- ComplexHeatmap::HeatmapAnnotation(df = rowAnno[rownames(plotMat),,drop=FALSE], col=annoCol, which = "row",
                                               annotation_name_gp = gpar(fontface = "bold"))
  } else haRow <- NULL

  labelCol <- rep("black",nrow(plotMat))
  if (!is.null(highLight)) {
    labelCol[rownames(plotMat) %in% highLight] <-"red"
  }
  if (italicize) {
    labFont <- 3
  } else labFont <- 1

  plotName <- ifelse(is.null(plotName),setName,plotName)

  ComplexHeatmap::Heatmap(plotMat, col = colorRampPalette(c(colList[2],"white",colList[1]))(100),name = "z-score",
                          top_annotation = haCol, left_annotation = haRow, show_column_names = FALSE,
                          cluster_columns  = FALSE, cluster_rows = FALSE,
                          row_names_gp = gpar(col = labelCol, fontface = labFont),
                          column_title = plotName, column_title_gp = gpar(cex= 1.5, fontface = "bold")
  )
}


runFGSEA <- function(limmaRes, signatures, stat = "t", name = "symbol",minSize=15,
                     maxSize =1000,nperm = 10000, pcut = 0.1, ifFDR = TRUE, showFDR = FALSE) {
  geneRanks <- limmaRes[[stat]]
  names(geneRanks) <- limmaRes[[name]]
  fgseaRes <- fgsea(pathways = sigList,
                    stats = geneRanks,
                    minSize=minSize,
                    maxSize=maxSize,
                    nperm=nperm)
  if(ifFDR) {
    plotSets <- filter(fgseaRes, padj <= pcut)$pathway
  } else {
    plotSets <- filter(fgseaRes, pval <= pcut)$pathway
  }

  #plot gsea plots
  plots <- lapply(plotSets, function(setName) {
    pval <- formatNum(filter(fgseaRes, pathway == setName)$pval, digits = 1)
    FDR <- formatNum(filter(fgseaRes, pathway == setName)$padj, digits = 1)
    nes <- filter(fgseaRes, pathway == setName)$NES
    NES <- format(nes, digits = 1,nsmall = 1)
    #showValue <- ifelse(showFDR, sprintf("italic(P)~'value = %s\nFDR = %s\nNormalized ES=%s'", pval, FDR, NES),
    #                    sprintf("italic(P)~'value = %s\nNormalized ES=%s'", pval, NES))
    pp <- fgsea::plotEnrichment(signatures[[setName]],
                         geneRanks)
    pp <- pp + ggtitle(setName) +
      ylab("Enrichment score (ES)") + xlab("Rank") +
      theme(plot.title = element_text(face = "bold", size = 15),
            plot.margin = margin(2,1,2,1,"cm"),
            axis.text = element_text(size=15),
            axis.title = element_text(size=16)
            )
    if (nes >0) {
      showValue <- sprintf("atop('          '~italic(P)~'=%s','Normalized ES=%s')",pval,NES)
      pp <- pp + annotate("text", x = Inf, y = Inf,
                          label = showValue, parse= TRUE,
                          hjust=1, vjust =1.3, size = 5)
    } else {
      showValue <- sprintf("atop(italic(P)~'=%s            ','Normalized ES=%s')",pval,NES)
      pp <- pp + annotate("text", x = 0, y = -Inf,
                          label = showValue, parse=TRUE,
                          hjust=0, vjust =-0.5, size = 5)
    }
    pp
  })
  names(plots) <- plotSets
  return(list(table = fgseaRes, plots = plots))
}

#Function to format floats
formatNum <- function(i, limit = 0.01, digits =1, format="e") {
  r <- sapply(i, function(n) {
    if (n < limit) {
      formatC(n, digits = digits, format = format)
    } else {
      format(n, digits = digits)
    }
  })
  return(r)
}

sumToTidy <- function(seObject, rowID = "rowID", colID = "colID") {

  tidyTable <- lapply(assayNames(seObject),function(n) {
    valTab <- assays(seObject)[[n]] %>% data.frame() %>%
      rownames_to_column(rowID) %>%
      gather(key = !!colID, value = "val", -!!rowID) %>%
      mutate(assay = n)
  }) %>% bind_rows() %>%
    spread(key = assay, value = val)

  #append row annotations
  rowAnno <- rowData(seObject) %>% data.frame() %>%
    rownames_to_column(rowID)
  tidyTable <- left_join(tidyTable, rowAnno, by = rowID)

  #append column annotations
  colAnno <- colData(seObject) %>% data.frame() %>%
    rownames_to_column(colID)
  tidyTable <- left_join(tidyTable, colAnno, by = colID)


  return(as_tibble(tidyTable))
}

runGSEA <- function(inputTab,gmtFile,GSAmethod="gsea",nPerm=1000){
  require(piano)
  inGMT <- loadGSC(gmtFile,type="gmt")
  rankTab <- inputTab[order(inputTab[,1],decreasing = TRUE),,drop=FALSE] #re-rank by score

  if (GSAmethod == "gsea"){
    #readin geneset database
    #GSEA analysis
    res <- runGSA(geneLevelStats = rankTab,geneSetStat = GSAmethod,adjMethod = "fdr",
                  gsc=inGMT, signifMethod = 'geneSampling', nPerm = nPerm, verbose = FALSE)
    GSAsummaryTable(res)
  } else if (GSAmethod == "page"){
    res <- runGSA(geneLevelStats = rankTab,geneSetStat = GSAmethod,adjMethod = "fdr",
                  gsc=inGMT, signifMethod = 'nullDist', verbose = FALSE)
    GSAsummaryTable(res)
  }
}

plotEnrichmentBar <- function(resTab, pCut = 0.05, ifFDR = FALSE, setName = "", title="",
                              removePrefix = NULL, insideLegend = FALSE, setMap = NULL) {

    plotTab <- resTab

    if (ifFDR) {
      plotTab <- dplyr::filter(plotTab, `p adj (dist.dir.up)` <= pCut | `p adj (dist.dir.dn)` <= pCut)
    } else {
      plotTab <- dplyr::filter(plotTab, `p (dist.dir.up)` <= pCut | `p (dist.dir.dn)` <= pCut)
    }

    if (nrow(plotTab) == 0) {
      print("No sets passed the criteria")
      return(NULL)

    } else {
      #firstly, process the result table
      plotTab <- lapply(seq(nrow(plotTab)), function(i) {
        x <- plotTab[i,]
        statSign <- as.numeric(x[3])
        data.frame(Name = x[[1]], p = as.numeric(ifelse(statSign >= 0, x[[4]], x[[6]])),
                   geneNum = ifelse(statSign >= 0, x[[8]], x[[9]]),
                   Direction = ifelse(statSign > 0, "Up", "Down"), stringsAsFactors = FALSE)
      }) %>% bind_rows()

      if (!is.null(removePrefix)) plotTab <- mutate(plotTab, Name = str_remove(Name, removePrefix))
      if(!is.null(setMap)) {
        plotTab <- mutate(plotTab, newName = setMap[match(Name, setMap$setName),]$pathwayName) %>%
          mutate(Name = ifelse(is.na(newName), Name, newName))
      }

      plotTab$Name <- sprintf("%s (%s)",plotTab$Name,plotTab$geneNum)
      plotTab <- plotTab[with(plotTab,order(Direction, p, decreasing=TRUE)),]
      plotTab$Direction <- factor(plotTab$Direction, levels = c("Up","Down"))
      plotTab$Name <- factor(plotTab$Name, levels = plotTab$Name)
      #plot the barplot
      p <- ggplot(data=plotTab, aes(x=Name, y= -log10(p), fill=Direction)) +
        geom_bar(position="dodge",stat="identity", width = 0.5) +
        scale_fill_manual(values=c(Up = colList[1], Down = colList[2])) +
        coord_flip() + xlab(setName) +
        ylab(expression(-log[10]*'('*italic(P)*' value)')) +
        ggtitle(title) + theme_full + theme(plot.title = element_text(face = "bold", hjust =0.5),
                                        axis.title = element_text(size=15), legend.background = element_rect(fill = NA))

      if (insideLegend) {
        p <- p + theme(legend.position = c(0.8,0.1))
      } else {
        p <- p + theme(legend.position = "right")
      }
    }


  return(p)
}



plotVolcano <- function(pTab, fdrCut = 0.05, posCol = "red", negCol = "blue",
                        x_lab = "dm", plotTitle = "",ifLabel = FALSE, labelList=NULL,
                        colLabel = NULL, yLim = NULL) {
  plotTab <- pTab %>% mutate(ifSig = ifelse(adj.P.Val > fdrCut, "n.s.",
                                            ifelse(logFC > 0, "up","down"))) %>%
    mutate(ifSig = factor(ifSig, levels = c("up","down","n.s."))) %>%
    mutate(P.Value = ifelse(P.Value < 1e-16, 1e-16, P.Value))
  pCut <- -log10((filter(plotTab, ifSig != "n.s.") %>% arrange(desc(P.Value)))$P.Value[1])

  if (is.null(yLim)) yLim <- max(-log10(plotTab$P.Value))

  g <- ggplot(plotTab, aes(x=logFC, y=-log10(P.Value))) +
    geom_point(shape = 21, aes(fill = ifSig),size=3) +
    geom_hline(yintercept = pCut, linetype = "dashed") +
    annotate("text", x = -Inf, y = pCut, label = paste0(fdrCut*100,"% FDR"),
             size = 5, vjust = -1.2, hjust=-0.1) +
    scale_fill_manual(values = c(n.s. = "grey70",
                                 up = posCol, down = negCol),name = "") +
    ylim(0,yLim) +
    theme_full +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 15),
          axis.text = element_text(size=20),
          axis.title = element_text(size=20),
          plot.title = element_text(size=20)) +
    ylab(expression(-log[10]*'('*italic(P)~value*')')) +
    xlab(x_lab) + ggtitle(plotTitle)

  if (ifLabel) {
    if (is.null(labelList)) {
      labelTab <- plotTab
    } else {
      labelTab <- filter(plotTab, name %in% labelList, ifSig != "n.s.")
    }

    if (is.null(colLabel)) {
      g <- g + ggrepel::geom_label_repel(data = labelTab, aes(label = name),
                                         size=7, force = 5, col = "red")
    } else {
      g <- g+ggrepel::geom_label_repel(data = labelTab,
                                       aes_string(label = "name", col = colLabel),
                                       size=7, force = 5) +
        scale_color_manual(values = c(yes = "red",no = "black"))
    }
  }

  return(g)
}

plotBox <- function(plotTab, pValTabel = NULL, y_lab = "Protein expression") {
  plotList <- lapply(unique(plotTab$name), function(n) {
    eachTab <- filter(plotTab, !is.na(status), !is.na(count), name == n) %>%
      group_by(status) %>% mutate(n=n()) %>% ungroup() %>%
      mutate(group = sprintf("%s\n(N=%s)",status,n)) %>%
      arrange(status) %>% mutate(group = factor(group, levels = unique(group)))

    if (!is.null(pValTabel)) {
      pval <- filter(pValTabel, name == n)$P.Value
      if (pval < 10^-12) {
        pval <- 10^-12
        pval <- formatNum(pval, digits = 1, format="e")
        titleText <- bquote(italic(.(n))~" ("~italic("P")~"<"~.(pval)~")")
      } else {
        pval <- formatNum(pval, digits = 1, format="e")
        titleText <- bquote(italic(.(n))~" ("~italic("P")~"="~.(pval)~")")
      }
    } else {
      titleText <- bquote(italic(.(n)))
    }

    ggplot(eachTab, aes(x=group, y = count)) +
      geom_boxplot(width=0.3, aes(fill = group), outlier.shape = NA) +
      geom_beeswarm(col = "black", size =2.5,cex = 2, alpha=0.5) +
      ggtitle(titleText)+
      #ggtitle(sprintf("%s (p = %s)",geneName, formatNum(pval, digits = 1, format = "e"))) +
      ylab(y_lab) + xlab("") +
      scale_fill_manual(values = colList[3:5]) +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.1))+
      theme_full +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5),
            plot.margin = margin(10,10,10,10))

  })
  return(plotList)
}

plotPep <- function(pepCLL, protName, type = "count", stratifier = NULL, mutStatus = NULL, protExpr = NULL) {

  #get peptide expression
  pepTab <- assays(pepCLL[rowData(pepCLL)$symbol %in% protName,])[[type]] %>%
    data.frame() %>% rownames_to_column("pepID") %>%
    gather(key = "patID", value = "expr", -pepID) %>%
    mutate(sequence = rowData(pepCLL[pepID,])$pepSeq) %>%
    dplyr::select(patID, expr, sequence) %>%
    mutate(measure = "Peptide")

  #prepare protein expression annotation table
  if (!is.null(protExpr)) {
    protTab <- tibble(patID = names(protExpr), expr = protExpr, sequence = "Protein", measure = "Protein")
    pepTab <- bind_rows(pepTab, protTab)
  }

  # prepare plotting table
  plotTab <- pepTab %>%
    group_by(patID) %>% mutate(medVal = median(expr, na.rm=TRUE)) %>%
    ungroup() %>% arrange(desc(medVal)) %>%
    mutate(patID = factor(patID, levels = unique(patID))) %>%
    dplyr::select(-medVal)

  #add genomic
  if (!is.null(mutStatus)) {
    plotTab <- mutate(plotTab, mut = mutStatus[match(patID, names(mutStatus))]) %>%
      filter(!is.na(mut)) %>% mutate(mut = paste0(stratifier, " = ", mut))
  }

  p <- ggplot(plotTab, aes(x = patID, y = expr, col = sequence, group = sequence)) +
    geom_point() +  theme_full +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
    xlab("") + ylab("Expression") + ggtitle(protName)

  #if plot protein expression level
  if (!is.null(protExpr)) {
    p <- p + geom_line(aes(linetype = measure)) + scale_shape_manual(values = c(Peptide = "solid",Protein = "dotted"))
  } else {
    p <- p + geom_line()
  }

  if (!is.null(mutStatus)) {
    p <- p + facet_wrap(~mut, scales = "free_x",ncol =length(unique(mutStatus)), drop = TRUE)
    gp <- ggplotGrob(p)
    facet.columns <- gp$layout$l[grepl("panel", gp$layout$name)]
    x.var <- sapply(ggplot_build(p)$layout$panel_scales_x,
                    function(l) length(l$range$range))
    gp$widths[facet.columns] <- gp$widths[facet.columns] * x.var
    p <- gp
  }

  p
}

plotColList <- function() {
  plotTab <- tibble(ID = as.factor(seq(length(colList))), y =1)
  ggplot(plotTab, aes(x=ID, y=y, col = ID)) +
    geom_point(size=5) +
    scale_color_manual(values = colList)
}


shareLegend <- function(plotList, position = "left", ratio = c(1,0.1), ncol=2,
                        rel_widths = NULL) {
  #get legend
  legend <- get_legend(plotList[[1]])

  if (is.null(rel_widths)) rel_widths <- rep(1, length(plotList))

  #main plot
  plotList <- lapply(plotList, function(p) p+theme(legend.position = "none"))
  mainPlot <- plot_grid(plotlist = plotList, ncol=ncol,
                        rel_widths = rel_widths)

  #plot the whole plot
  if (position == "bottom") {
    plot_grid(mainPlot, legend, ncol=1, rel_heights  = ratio)
  } else if (position == "left") {
    plot_grid(mainPlot, legend, ncol=2, rel_widths  = ratio)

  }
}

