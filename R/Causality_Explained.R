#' Main function of Causality Explained
#'
#' @param CutOff which is the cutoff set to select IVs
#' @param DataPath folder containing the data (predictors, outcomes, information tables)
#' @param Outcomes outcomes to process
#' @param NBCores number of cores to run the pipeline in parallel
#' @param BonferroniCorrection if TRUE, apply a Bonferroni correction on predictor selection
#' @param MinNbrIVs The minimum number of significant IVs to use the trait as a predictor
#'
#' @return
#' @export
#'
#' @examples
Causality_Explained <- function(CutOff = 5e-08, DataPath, Outcomes, NBCores = 1, BonferroniCorrection = F, MinNbrIVs = 1){


  "%^%" <- function(x, n) with(eigen(x), vectors %*% (values^n * t(vectors)))

  # Get predictor data preprocessed
  load(paste0(DataPath,"AllData_UKBB_", CutOff))


  # Functions
  # modSummary(): For a given MR model, returns various indicators including adjusted predicted R2
  modSummary <- function(mod){
    getPredictedR2 <- function(mod){
      PRSS <- sum((weighted.residuals(mod) / (1 - lm.influence(mod)$hat))^2)
      predR2 <- 1 - PRSS / sum(anova(mod)$"Sum Sq")
      AdjPredR2 <- 1 - (1 - predR2) * ((length(mod$residuals) - 1) / (mod$df.residual -1))
      return(c(predR2 = predR2, AdjPredR2 = AdjPredR2))
    }
    cbind.data.frame(
      NBindepSNPs = length(mod$residuals),
      NBPhecodes = length(mod$coefficients),
      `P<0.05` = sum(summary(mod)$coefficients[, 4] < 0.05),
      `P<BF` = sum(summary(mod)$coefficients[, 4] < (0.05 / nrow(summary(mod)$coefficients))),
      R2 = summary(mod)$r.squared,
      AdjR2 = summary(mod)$adj.r.squared,
      t(getPredictedR2(mod))
    )
  }

  # PCmodel(): For a set of predictors beta coefficients and and outcome summary stats (beta and se) with the same IVs in rows, it performs the polytrait model with 3 PC selection procedures.
  PCmodel <- function(dataPCA, OutcomePCA, outcomePCA = outcome){
    dataPCA0 <- scale(sweep(dataPCA, 1, OutcomePCA$Sd ,"/"), scale = TRUE, center = TRUE)
    dataPCA <- scale(sweep(dataPCA, 1, OutcomePCA$Sd ,"/"), scale = FALSE, center = FALSE)
    resMetPCA0 <- svd(dataPCA0/sqrt(nrow(dataPCA0)))
    PCselect <- max(1, round(length(resMetPCA0$d) - sum(as.numeric(resMetPCA0$d^2 > 1) * (resMetPCA0$d^2 - 1)), digit = 0))
    resMetPCA <- svd(dataPCA/sqrt(nrow(dataPCA)))
    PCselect1 <- max(1, sum(resMetPCA$d^2 > 1))

    resModPCA <- lapply(c(PCselect, PCselect1, ncol(dataPCA)), function(NbPCs){ # 3 PC selection procedures
      if(NbPCs == 1)
        MetPCA <- data.frame(resMetPCA$u[, NbPCs] * resMetPCA$d[NbPCs])
      else
        MetPCA <- sweep(resMetPCA$u[, 1:NbPCs], 2, resMetPCA$d[1:NbPCs], "*")
      row.names(MetPCA) <- row.names(dataPCA)
      colnames(MetPCA) <- paste0("PC", 1:ncol(MetPCA))
      modPCA <- lm(as.formula(paste0("Beta ~ -1 + ", paste0("PC", 1:NbPCs, collapse = " + "))), data = cbind.data.frame(sweep(OutcomePCA, 1, OutcomePCA$Sd,"/"), MetPCA))
      res <- modSummary(modPCA)
      return(res)
    })
    resModPCA <- cbind.data.frame(Outcome = gsub("\\.txt", "", outcomePCA), PC = c("PCselect", "PCselect>1", "All"), do.call("rbind", resModPCA))
    return(resModPCA)
  }


  # Run MV MR
  resAllOutcomes <- mclapply(Outcomes, function(outcome){ # run the pipeline on all outcomes
    cat(paste(which(Outcomes == outcome), "\n"))

    Outcome <- fread(paste0(DataPath, outcome), data.table = FALSE)
    if(any(duplicated(Outcome$SNPid))){
      Duplicates <- unlist(lapply(Outcome$SNPid[duplicated(Outcome$SNPid)], function(snp) {
        dup <- Outcome[Outcome$SNPid == snp, ]
        dup <- setdiff(which(Outcome$SNPid == snp), which(Outcome$SNPid == snp & Outcome$N == max(dup$N))[1])
      }))
      Outcome <- Outcome[-Duplicates, ]
    }
    row.names(Outcome) <- Outcome$SNPid

    # Phecode category restriction
    phecodes <- fread(paste0(DataPath,"PheCodeExclusion_updated.txt"), data.table = FALSE) # File provided by Iain and Aine to filter out phecodes similar to the outcome
    phecodes <- phecodes$Phecode[phecodes$Trait == gsub("\\.txt", "", outcome)]
    phecodes <- unique(gsub("([0-9]+)(\\.[0-9]*)", "\\1", phecodes))
    if(length(phecodes) > 0)
      phecodes <- grep(paste0("^", phecodes, "\\.*.*", collapse = "|"), gsub("(PheCode_0*)([0-9\\.]+)(_Beta)", "\\2", grep("_Beta", colnames(UKBBData_AllSNPs), value = TRUE)), value = TRUE)
    phecodesDef <- fread(paste0(DataPath,"phenotype-information.txt"), data.table = FALSE)

    ########## Restriction to phecodes outside of cat
    row.names(UKBBData_AllSNPs) <- SNPs$SNPid
    UKBBData_AllSNPs <- UKBBData_AllSNPs[apply(abs(UKBBData_AllSNPs[, grep("Beta", colnames(UKBBData_AllSNPs))] / UKBBData_AllSNPs[, grep("Sd", colnames(UKBBData_AllSNPs))])> qnorm(CutOff/2, lower.tail = FALSE), 1, sum, na.rm = TRUE) > 0, ]
    Outcome <- Outcome[match(row.names(UKBBData_AllSNPs), Outcome$SNPid), c(6, 7)]

    ########## Restriction to phecodes outside of cat
    UKBBData_AllSNPsAll <- UKBBData_AllSNPs[, paste0(rep(gsub("_Beta", "", grep("_Beta", colnames(UKBBData_AllSNPs), value = TRUE)[ ! gsub("(PheCode_0*)([0-9\\.]+)(_Beta)", "\\2", grep("_Beta", colnames(UKBBData_AllSNPs), value = TRUE)) %in% phecodes]), each = 2), c("_Beta", "_Sd"))]

    ########## Restriction to phecodes with significant IVs and outside of cat
    UKBBData_AllSNPsO <- UKBBData_AllSNPs
    UKBBData_AllSNPsO <- UKBBData_AllSNPsO[, paste0(rep(gsub("_Beta", "", names(which(apply(abs(UKBBData_AllSNPsO[, grep("Beta", colnames(UKBBData_AllSNPsO))] /
          UKBBData_AllSNPsO[, grep("Sd", colnames(UKBBData_AllSNPsO))])> qnorm(CutOff/2, lower.tail = FALSE), 2, sum, na.rm = TRUE) > 0))), each = 2), c("_Beta", "_Sd"))]

    UKBBData_AllSNPsO <- UKBBData_AllSNPsO[, paste0(rep(gsub("_Beta", "", grep("_Beta", colnames(UKBBData_AllSNPsO), value = TRUE)[ ! gsub("(PheCode_0*)([0-9\\.]+)(_Beta)", "\\2", grep("_Beta", colnames(UKBBData_AllSNPsO), value = TRUE)) %in% phecodes]), each = 2), c("_Beta", "_Sd"))]
    UKBBData_AllSNPsO <- UKBBData_AllSNPsO[apply(abs(UKBBData_AllSNPsO[, grep("Beta", colnames(UKBBData_AllSNPsO))] / UKBBData_AllSNPsO[, grep("Sd", colnames(UKBBData_AllSNPsO))])> qnorm(CutOff/2, lower.tail = FALSE), 1, sum, na.rm = TRUE) > 0, ]

    # Remove missing values
    OutcomeO <- Outcome[match(row.names(UKBBData_AllSNPsO), row.names(Outcome)), ]
    UKBBData_AllSNPsO <- UKBBData_AllSNPsO[ ! is.na(OutcomeO$Beta), ]
    OutcomeO <- OutcomeO[ ! is.na(OutcomeO$Beta), ]
    UKBBData_AllSNPsAll <- UKBBData_AllSNPsAll[ ! is.na(Outcome$Beta), ]
    Outcome <- Outcome[ ! is.na(Outcome$Beta), ]

    GWlist <- abs(UKBBData_AllSNPsO[, grep("Beta", colnames(UKBBData_AllSNPsO))] / UKBBData_AllSNPsO[, grep("Sd", colnames(UKBBData_AllSNPsO))])> qnorm(CutOff/2, lower.tail = FALSE)
    GWlist <- melt(GWlist)
    GWlist <- GWlist[GWlist$value, 1:2]
    colnames(GWlist) <- c("SNP", "PheCode")

    # Remove phecodes with less than MinNbrIVs instruments
    UKBBData_AllSNPsO <- UKBBData_AllSNPsO[, colnames(UKBBData_AllSNPsO) %in% paste0(rep(gsub("Beta$", "", names(which(table(GWlist$PheCode) >= MinNbrIVs))), each = 2), c("Beta", "Sd"))]
    GWlist <- GWlist[GWlist$PheCode %in% names(which(table(GWlist$PheCode) >= MinNbrIVs)), ]


    ### Update data with final set of IVs
    UKBBData_AllSNPsAll<-UKBBData_AllSNPsAll[rownames(UKBBData_AllSNPsAll)%in%rownames(UKBBData_AllSNPsO),]
    Outcome<-Outcome[rownames(Outcome)%in%rownames(UKBBData_AllSNPsO),]

    #####################################################################
    # Effect of individual Phecode
    AllPhecodes <- na.omit(UKBBData_AllSNPsAll[, grep("Beta", colnames(UKBBData_AllSNPsAll), value = TRUE)])
    AllPhecodes <- scale(sweep(AllPhecodes, 1, Outcome$Sd ,"/"), scale = FALSE, center = FALSE)

    resModPhecodes <- cbind.data.frame(PheCode = gsub("^(PheCode_)(.*)(_Beta)$", "\\2", colnames(AllPhecodes)),
                                       do.call("rbind", lapply(colnames(AllPhecodes), function(phecode){
                                         mod <- lm(as.formula(paste0("Beta ~ -1 + ", phecode)), data = cbind.data.frame(sweep(Outcome, 1, Outcome$Sd,"/"), AllPhecodes))
                                         res <- cbind(Estimate = summary(mod)$coefficients[1, 1], Pvalue = summary(mod)$coefficients[1, 4], `Individual R2` = summary(mod)$r.squared, NbIVsPhecode = as.numeric(table(GWlist$PheCode)[phecode]))
                                         return(res)
                                       })))


    if(BonferroniCorrection == T){ #We filter Phecodes that are not sig in MR
      resModPhecodes <- resModPhecodes[resModPhecodes$Pvalue < (0.05/sum(!is.na(resModPhecodes$NbIVsPhecode))),]
    }

    # Get phecodes with GW IVs
    resModPhecodes$NbIVsPhecode[is.na(resModPhecodes$NbIVsPhecode)] <- 0
    PhecodesGW <- gsub("^(PheCode_)(.*)(_Beta)$", "\\2", unique(GWlist$PheCode))

    #####################################################################
    # Include phecodes 1 by 1

    cat(paste(which(Outcomes == outcome), " ------------- starting PC model----------\n"))

    FirstPheCode <- resModPhecodes$PheCode[resModPhecodes$PheCode %in% PhecodesGW][which.max(resModPhecodes$`Individual R2`[resModPhecodes$PheCode %in% PhecodesGW])]

    StepPhecodesPC <- resModPhecodes
    StepPhecodesPC$GW <- StepPhecodesPC$PheCode %in% PhecodesGW
    StepPhecodesPC <- StepPhecodesPC[order(StepPhecodesPC$GW, StepPhecodesPC$`Individual R2`, decreasing = TRUE), ]


    Ind <- c(2:sum(StepPhecodesPC$GW))
    Ind <- Ind[Ind <= sum(StepPhecodesPC$GW)]

    StepPhecodesPC <- cbind.data.frame(StepPhecodesPC, NBPhecodesTot = 1:nrow(StepPhecodesPC))
    StepPhecodesPC <- lapply(Ind, function(i){
      mod <- PCmodel(dataPCA = UKBBData_AllSNPsAll[, paste("PheCode", StepPhecodesPC$PheCode[1:i], "Beta", sep = "_")], OutcomePCA = Outcome, outcomePCA = outcome)
      return(cbind.data.frame(StepPhecodesPC[i, ], mod))
    })
    StepPhecodesPC <- do.call("rbind", StepPhecodesPC)

    mod1 <- lm(as.formula(paste0("Beta ~ -1 + ", paste0("PheCode_", FirstPheCode, "_Beta"))), data = cbind.data.frame(Outcome, UKBBData_AllSNPsAll), weights = 1/Sd^2)

    StepPhecodesPC <- rbind(
      cbind.data.frame(resModPhecodes[resModPhecodes$PheCode == FirstPheCode, ], GW = TRUE, NBPhecodesTot = 1, Outcome = StepPhecodesPC$Outcome[1], PC = c("All", "PCselect", "PCselect>1"), modSummary(mod1)),
      StepPhecodesPC
    )
    ####################################################################################

    res <- StepPhecodesPC
    return(res)
  }, mc.cores = NBCores
  )
}






