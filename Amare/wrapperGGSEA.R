#' @title a wrapper function for GGSEA.
#' @description \code{wrapperGGSEA} take expression data and 
#' @param DistM a distance matrix object containing the pairwise distances between genes.
#' @param PathwayL a list object containing the biological pathways to test for enrichment. 
#' @param TestType a type of enrichment test to use.
#' The enrichment test should be one of the function such as CalcNodeWeight, TreeDistBetweenness, SubTreeModule, MultidimDist, MFPT_distr, DistributionMaximum, EnrichFisherPathway, or SimpleThreshold.
##-------------working---------------------------------------------------------- 
wrapperGGSEA <-  function(DistM, PathwayL, TestType,desc = TRUE,
                          thres_quant = 0.05,...) {
  
  # Get the expression data
  # do the distance , correlation matrix, diff.distance 
  # Calculate the enrichment score using the selected test
  enrichmentScore <- switch(TestType,
                            "CalcNodeWeight" = CalcNodeWeight(DistM, PathwayL),
                            "TreeDistBetweenness" = TreeDistBetweenness(DistM, PathwayL),
                            "SubTreeModule" = SubTreeModule(DistM, PathwayL),
                            "MultidimDist" = MultidimDist(DistM,PathwayL),
                            "MFPT_distr" = MFPT_distr(DistM, PathwayL),
                            "DistributionMaximum" = DistributionMaximum(DistM, PathwayL),
                            "EnrichFisherPathway" = EnrichFisherPathway(DistM, PathwayL),
                            "SimpleThreshold" = SimpleThreshold(DistM, PathwayL, 
                                                                desc = desc, thres_quant = thres_quant)
  )
  # Combine the scores into a matrix with one row per pathway and one column per gene
  enrichmentScore
  
  # visualizations 
}

pathwayE <- wrapperGGSEA(DistM = distanceMatrix,
                         TestType = "TreeDistBetweenness", 
                         PathwayL = PathwaysTest)
head(pathwayE)


##------------------------------------------------------------------------------

wrapperGGSEA <-  function(DistM, PathwayL, TestType,desc = TRUE,
                          thres_quant = 0.05,...) {
  
  results <- sapply(names(PathwayL), function(pathwayName) {
    
    # extract genes in the current pathway and pass it to the DistM object,
    # by doing so you would make the test fast and save time
    pathwayGenes <- unique(unlist(PathwayL[["pathwayName"]]))
    DistM <- DistM[rownames(DistM) %in% pathwayGenes,  colnames(DistM) %in% pathwayGenes]
    
    # Calculate the enrichment score using the selected test
    enrichmentScore <- switch(TestType,
                              "CalcNodeWeight" = CalcNodeWeight(DistM, PathwayL[["pathwayName"]]),
                              "TreeDistBetweenness" = TreeDistBetweenness(DistM, PathwayL[["pathwayName"]]),
                              "SubTreeModule" = SubTreeModule(DistM, PathwayL[["pathwayName"]]),
                              "MultidimDist" = MultidimDist(DistM,PathwayL[["pathwayName"]]),
                              "MFPT_distr" = MFPT_distr(DistM, PathwayL[["pathwayName"]]),
                              "DistributionMaximum" = DistributionMaximum(DistM, PathwayL[["pathwayName"]]),
                              "EnrichFisherPathway" = EnrichFisherPathway(DistM, PathwayL[["pathwayName"]]),
                              "SimpleThreshold" = SimpleThreshold(DistM, PathwayL[["pathwayName"]], 
                                                                  desc = desc, thres_quant = thres_quant)
    )
    enrichmentScore
  })
  return(enrichmentScore)
}

ERes <- wrapperGGSEA(DistM = distanceMatrix,
                     TestType = "TreeDistBetweenness", 
                     PathwayL = PathwaysTest)


# alternative function however not finished 
wrapperGGSEA <- function(DistM, PathwayL, TestType, 
                         desc = TRUE, thres_quant = 0.05, ...) {
  
  # extract genes in the current pathway name and pass it to the DistM object
  pathwayGenes <- sapply(PathwayL, function(pathwayName) {
    unique(unlist(pathwayName))
  })
  
  # the use for loop to iterate in each pathway  and apply the different test types 
  for (i in seq_along(PathwayL)) {
    indices <-   pathwayGenes <- unique(unlist(PathwayL[[pathwayName]]))
    DistM <- DistM[rownames(DistM) %in% pathwayGenes,  colnames(DistM) %in% pathwayGenes]
    enrichmentScore <-  switch(TestType,
           "CalcNodeWeight" = {
             score <- CalcNodeWeight(DistM, PathwayL[[i]])
             #pval <- NA
           },
           "TreeDistBetweenness" = {
             score <- TreeDistBetweenness(DistM, PathwayL[[i]])
           },
           "SubTreeModule" = {
             score <- SubTreeModule(DistM, PathwayL[[i]])
           },
           "MultidimDist" = {
             score <- MultidimDist(DistM, PathwayL[[i]])
           },
           "MFPT_distr" = {
             score <- MFPT_distr(DistM, PathwayL[[i]])
           },
           "DistributionMaximum" = {
             score <- DistributionMaximum(DistM, PathwayL[[i]])
           },
           "EnrichFisherPathway" = {
             res <- EnrichFisherPathway(DistM, PathwayL[[i]])
           },
           "SimpleThreshold" = {
             score <- SimpleThreshold(DistM, PathwayL[[i]], desc = desc, thres_quant = thres_quant)
           }
    )
    
    # generate the data frame with four columns for pathway name, enrichment score, p_value and adj_pvalue 
    results <- data.frame(Pathway = character(length(PathwayL)),
                          EnrichmentScore = numeric(length(PathwayL)),
                          # PValue = numeric(length(PathwayL)),
                          # AdjPValue = numeric(length(PathwayL))
                          )
    results[i, "Pathway"] <- names(PathwayL)[i]
    results[i, "EnrichmentScore"] <- score
    # results[i, "PValue"] <- pval
    # results[i, "AdjPValue"] <- p.adjust(pval, method = "BH")
  }
  
  return(results)
}


ERes <- wrapperGGSEA(DistM = distanceMatrix,
                     TestType = "TreeDistBetweenness", 
                     PathwayL = PathwaysTest)





