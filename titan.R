#!/srv/gsfs0/projects/curtis/ruping/tools/R/bin/Rscript

## this is for WGS or WXS titan run

inputpar <- commandArgs(TRUE)
if (length(inputpar) < 7) stop("Wrong number of input parameters: 'path sampleName alleleCount tumorWig normalWig gcWig mapWig'")
#for (i in 1:7) {
#    if (!file.exists(inputpar[i]))
#          stop("File ",inputpar[i]," doesn't exist")
#  }


path <- inputpar[1]
sampleName <- inputpar[2]
alleleCount <- inputpar[3]
tumorWig <- inputpar[4]
normalWig <- inputpar[5]
gcWig <- inputpar[6]
mapWig <- inputpar[7]
exons <- inputpar[8]

library(TitanCNA)
library(HMMcopy)
library(doMC)

setwd(path)

#run titan
runTitan <- function(sampleName, snpFile, tumWig, normWig, gc, map, plp, plpe, normalc, exons) {

    #prepare data
    snpData <- loadAlleleCounts(snpFile)
    cnData <- correctReadDepth(tumWig, normWig, gc, map)
    if (length(inputpar) == 8){
      cnData <- correctReadDepth(tumWig, normWig, gc, map, targetedSequence = exons)
    }
    logR <- getPositionOverlap(snpData$chr, snpData$posn, cnData)
    snpData$logR <- log(2^logR) #transform the log ratio to natural logs
    snpData <- filterData(snpData, 1:22, minDepth = 10, maxDepth = 400, positionList = NULL)
    #prepare data

    registerDoMC(cores = 2)
    titancnaresults <- vector("list",2)
    
    for (j in 1:2) {
        numClusters <- j
        params <- loadDefaultParameters(copyNumber = 5, numberClonalClusters = numClusters)
        K <- length(params$genotypeParams$alphaKHyper)
        params$genotypeParams$alphaKHyper <- rep(1000, K)
        params$normalParams$n_0 <- normalc
        params$ploidyParams$phi_0 <- plp
        
        convergeParams <- runEMclonalCN(snpData, gParams = params$genotypeParams,
                                        nParams = params$normalParams,
                                        pParams = params$ploidyParams,
                                        sParams = params$cellPrevParams,
                                        maxiter = 10, maxiterUpdate = 500,
                                        useOutlierState = FALSE, txnExpLen = 1e9,
                                        txnZstrength = 1e9,
                                        normalEstimateMethod = "map",
                                        estimateS = TRUE, estimatePloidy = plpe)
        
        optimalPath <- viterbiClonalCN(snpData, convergeParams)
        if (length(unique(optimalPath)) == 1) next
        results <- outputTitanResults(snpData, convergeParams, optimalPath,
                                      filename = NULL, posteriorProbs = FALSE,
                                      subcloneProfiles = TRUE)
        ploidy <- tail(convergeParams$phi, 1)
        norm <- tail(convergeParams$n, 1)
        cellularity = 1 - convergeParams$s[, ncol(convergeParams$s)] # estimated cellular prevalence
        
        titancnaresults[[j]] <- list(S_DbwIndex=computeSDbwIndex(results)$S_DbwIndex,results=results,
                                     convergeParams=convergeParams)

        #generate segmentation files
        segmenttmp = titancna2seg(results, convergeParams)
        write.table(segmenttmp, file=paste(sampleName,"_nclones",numClusters,".TitanCNA.segments.txt",sep=""),
                    quote = F, row.names = F, sep = "\t")
        
        #make plots
        if (length(inputpar) == 7){
          for (chro in 1:22) {
            pdf(paste(sampleName,"_nclones",numClusters,"_chr", chro, ".TitanCNA.pdf",sep=""),width=11.5, height=8)
            if (is.null(titancnaresults[[j]])) next
            SD <- round(titancnaresults[[j]]$S_DbwIndex,3)
            nclones <- nrow(convergeParams$s)
            ploidy <- round(tail(convergeParams$phi, 1),2)
            meandepth <- round(mean(as.numeric(results$Depth)),2)
            npoints <- nrow(results)
            s <- round(convergeParams$s[1,ncol(convergeParams$s)],2)

            
            par(mfrow=c(2,1))
            par(mar=c(4,4,2,1))
            plotCNlogRByChr(results, chr = chro, ploidy = ploidy, ylim = c(-2, 2), cex=0.25,
                            main=paste(sampleName, " nc=", numClusters, sep=""),
                            xlab=paste("normC=", round(norm,3), " pl=", ploidy, " cellularity=", round(cellularity,3),
                              " SD=",SD," s=",s," nc=",nclones," np=",npoints," md=",meandepth,sep=""), cex.lab=0.8)
            par(mar=c(4,4,2,1))
            plotAllelicRatio(results, chr = chro, ylim = c(0, 1), cex = 0.25, xlab = paste("Chromosomes", chro, sep=" "), main = "", cex.lab=0.8)

            #par(mar=c(4,4,2,1))
            #plotClonalFrequency(results, chr = NULL, normal = norm, ylim = c(0, 1), cex = 0.25, xlab = "", main = "", cex.lab=.8)

            dev.off()
          }
        } else if (length(inputpar) == 8){
          
          pdf(paste(sampleName,"_nclones",numClusters,".TitanCNA.pdf",sep=""),width=11.5, height=8)
          if (is.null(titancnaresults[[j]])) next
          SD <- round(titancnaresults[[j]]$S_DbwIndex,3)
          nclones <- nrow(convergeParams$s)
          ploidy <- round(tail(convergeParams$phi, 1),2)
          meandepth <- round(mean(as.numeric(results$Depth)),2)
          npoints <- nrow(results)
          s <- round(convergeParams$s[1,ncol(convergeParams$s)],2)
          
          par(mfrow=c(2,1))
          par(mar=c(4,4,2,1))
          plotCNlogRByChr(results, chr = NULL, ploidy = ploidy, ylim = c(-2, 2), cex=0.25,
                          main=paste(sampleName, " nc=", numClusters, sep=""),
                          xlab=paste("normC=", round(norm,3), " pl=", ploidy, " cellularity=", round(cellularity,3),
                            " SD=",SD," s=",s," nc=",nclones," np=",npoints," md=",meandepth,sep=""), cex.lab=0.8)
          par(mar=c(4,4,2,1))
          plotAllelicRatio(results, chr = NULL, ylim = c(0, 1), cex = 0.25,xlab = "Chromosomes", main = "", cex.lab=0.8)

          #par(mar=c(4,4,2,1))
          #plotClonalFrequency(results, chr = NULL, normal = norm, ylim = c(0, 1), cex = 0.25, xlab = "", main = "", cex.lab=.8)
          dev.off()

        }
        
        #generate segmentation table (for adjusting somatic mutation frequencies)
        allstate <- paste(results$Chr,results$TITANstate,results$ClonalCluster)
        changepoints <- c(1,which(allstate[-1] != allstate[-length(allstate)])+1)      #get change points
        segments <- results[changepoints,c(1,2,2,6,7,8,9,10,11,12)]
        names(segments)[2:3] <- c("Position1","Position2")
        segments$Position2 <- results$Position[c(changepoints[-1]-1,length(allstate))]
        segments$NumMarker <- diff(c(changepoints,length(allstate)+1))
        segments[[2]] <- as.numeric(segments[[2]])
        segments[[3]] <- as.numeric(segments[[3]])
        segments[[4]] <- as.numeric(segments[[4]])
        segments[[5]] <- as.numeric(segments[[5]])
        segments[[6]] <- as.numeric(segments[[6]])
        segments[[7]] <- as.numeric(segments[[7]])
        segments[[10]] <- as.numeric(segments[[10]])
        segments$NumMarker <- diff(c(changepoints,length(allstate)+1))
        segments$normC <- norm
        segments$ploidy <- ploidy
        for (k in 1:nrow(segments)) {
            af <- as.numeric(results$AllelicRatio[changepoints[k]:(changepoints[k]+segments$NumMarker[k]-1)])
            segments$AllelicRatio[k] <- mean(pmax(af,1-af))    #smooth mirrored allelic ratio
            lr <- as.numeric(results$LogRatio[changepoints[k]:(changepoints[k]+segments$NumMarker[k]-1)])
            segments$LogRatio[k] <- mean(lr)                   #smooth depth log ratio
        }
        #write.table(segments, file=paste(sampleName,"_nclones",numClusters,".TitanCNA.result.txt",sep=""),
        #            quote = F, row.names = F, sep = "\t")
    }
    #save(titancnaresults,file=paste("./results/",sampleName,".TitanCNA.RData",sep=""))
}


titancna2seg <- function(titanresult,titanparams) {
  ploidy <- round(tail(titanparams$phi,1),3)
  n <- round(tail(titanparams$n,1),3)

  titanresult$Position <- as.integer(titanresult$Position)
  titanresult$LogRatio <- as.numeric(titanresult$LogRatio)
  titanresult$AllelicRatio <- as.numeric(titanresult$AllelicRatio)
  titanresult$CopyNumber <- as.numeric(titanresult$CopyNumber)
  titanresult$CellularPrevalence <- as.numeric(titanresult$CellularPrevalence)
  titanresult$ClonalCluster[is.na(titanresult$ClonalCluster)] <- 0
  
  cp2 <- c(which(titanresult$TITANstate[-1] != titanresult$TITANstate[-nrow(titanresult)] | 
                   titanresult$Chr[-1] != titanresult$Chr[-nrow(titanresult)] | 
                     titanresult$ClonalCluster[-1] !=  titanresult$ClonalCluster[-nrow(titanresult)]),
           nrow(titanresult))
  cp1 <- c(1,cp2[-length(cp2)]+1)

  cnv <- data.frame(chrom=titanresult$Chr[cp1],
                    loc.start=titanresult$Position[cp1],
                    loc.end=titanresult$Position[cp2],
                    num.mark=cp2-cp1+1,
                    seg.mean=titanresult$LogRatio[cp1],
                    copynumber=titanresult$CopyNumber[cp1],
                    allelicratio=titanresult$AllelicRatio[cp1],
                    LOHcall=titanresult$TITANcall[cp1],
                    cellularprevalence=titanresult$CellularPrevalence[cp1],
                    ploidy=ploidy,
                    normalproportion=n)
  for (j in 1:length(cp1)) {
    cnv$seg.mean[j] <- mean(titanresult$LogRatio[cp1[j]:cp2[j]])
    cnv$allelicratio[j] <- mean(0.5+abs(0.5-titanresult$AllelicRatio[cp1[j]:cp2[j]]))
  }
  cnv$logcopynumberratio <- log2(((cnv$copynumber - 2) * cnv$cellularprevalence + 2)/2)
  cnv$logcopynumberratio[is.na(cnv$logcopynumberratio)] <- 0
  ##cnv$seg.mean.2 <- cnv$seg.mean
  for (i in c(5,7,12)) {
    cnv[[i]] <- round(cnv[[i]],3)
  }
  return(cnv)
}



targetRegion = read.delim(exons, header=F)
targetRegion = data.frame(targetRegion[,1:3])
#targetRegion[,1] = gsub("chr","",targetRegion[,1])

runTitan(sampleName,alleleCount,tumorWig,normalWig,gcWig,mapWig,2,TRUE,0.5,targetRegion)
