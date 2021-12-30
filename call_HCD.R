################################################################################
################################    functions   ################################   
################################################################################

library(GENOVA)
makeBED <- function(idxes){
  cb <- cbind(idxes, idxes)
  prevLength = nrow(cb)
  newLength = nrow(cb) - 1
  while(prevLength != newLength){
    prevLength = nrow(cb)
    for(i in 1:newLength){
      if(i+1 > nrow(cb)){next()}
      if(cb[i,2] == (cb[i+1,1]-1)){
        cb[i,2] <- cb[i+1,2]
      }
      if(cb[i+1,1] %in% cb[i,1]:cb[i,2]){
        cb <- cb[-(i+1),]
      }
    }
    newLength = nrow(cb)
  }
  return(cb)
}

compartmentHammer <- function(CS_1, CS_2, chrom = 'chr1', arm_start = 0, arm_end = 1e6,
                              fragmentationSize = NULL, outlierQuantile = 0.995){
  
  # check the inputs
  CS_1 <- as.data.frame(CS_1)
  CS_2 <- as.data.frame(CS_2)
  
  if(is.null(fragmentationSize)){
    fragmentationSize <- median(diff(CS_1[,2])) * 5
  }
  # subset the CS's
  c1 <- CS_1[CS_1[,1] == chrom,]
  c2 <- CS_2[CS_2[,1] == chrom,]

  
  c1 <- c1[c1[,2] >= arm_start & c1[,3] <= arm_end,]
  c2 <- c2[c2[,2] >= arm_start & c2[,3] <= arm_end,]
  
  # filter outliers
  th <- quantile(c1[,4], c(1- outlierQuantile, outlierQuantile))
  c1[c1[,4] > th[2], 4] <- th[2]
  c1[c1[,4] < th[1], 4] <- th[1]
  
  th <- quantile(c2[,4], c(1- outlierQuantile, outlierQuantile))
  c2[c2[,4] > th[2], 4] <- th[2]
  c2[c2[,4] < th[1], 4] <- th[1]
  
  ## merge cs's
  MERGE <- merge(c1[,1:4], c2[,1:4], by = 1:3)
  MERGE <- MERGE[order(MERGE$start), ]
  INDF <- data.frame(pos1 = MERGE[,2],
                     pos2 = MERGE[,3],
                     c1  = MERGE[,4],
                     c2  = MERGE[,5],
                     absd   = abs(abs(MERGE[,4]) - abs(MERGE[,5])),
                     S   = as.integer(abs(apply(
                       cbind(sign(MERGE[,4]),sign(MERGE[,5])), 1, diff)/2)))

    if(cor(INDF$c1, INDF$c2) < 0){
    message("Compartment-scores seem flipped.
            Will now re-flip CS_2.")
    INDF$c2 <- -INDF$c2
    
  }
  
  # init hmm
  hmm <- depmix( list(c2~1, absd ~ 1),
                 family = list(gaussian(),gaussian()),
                 nstates = 2, data = INDF )
  # fit hmm
  hmmfit <- fit(hmm, verbose = F)
  post_probs <- posterior(hmmfit)
  
  # make beds
  LP = length(unique(post_probs$state))
  H = 0.5/LP
  bedlist = list()
  dlist = list()
  for(i in 1:LP){
    idxes <- which(post_probs$state == i)
    dlist[[i]] <- INDF[idxes, 'absd']
    
    bed <- makeBED(idxes)
    bed <- data.frame('seqnames' = chrom,
                      'start' = INDF[bed[,1],1],
                      'end' = INDF[bed[,2],2])
    bed <- bed[apply(bed[, 2:3], 1, diff) > fragmentationSize,]
    bedlist[[i]] <- bed
  }
  
  interestingState <- which.max(unlist(lapply(dlist, median)))
  if(interestingState != 1){
    tmp <- list()
    tmp[[1]] <- bedlist[[2]]
    tmp[[2]] <- bedlist[[1]]
    bedlist <- tmp
  }
  
  # annotate bedlist[[1]] with mean C1 and C2 values (i.e. hyper or hypo)
  bedlist[[1]]$d <- 0
  bedlist[[1]]$direction <- 'hyper'
  bedlist[[1]]$comp1 <- 'C'
  bedlist[[1]]$comp2 <- 'C'
  bedlist[[1]]$class <- 'meh'
  for(i in 1:nrow(bedlist[[1]])){
    region <- bedlist[[1]][i, ]
    data <- INDF[INDF$pos1 >= region$start & INDF$pos2 <= region$end, ]
    mc1 <- median(data$c1)
    mc2 <- median(data$c2)
    bedlist[[1]]$d[i] <- mc2 - mc1
    bedlist[[1]]$direction[i] <- ifelse( (abs(mc2) - abs(mc1)) > 0, 'hyper', 'hypo' )
    bedlist[[1]]$comp1[i] <- ifelse(mc1 < 0, 'B', 'A')
    bedlist[[1]]$comp2[i] <- ifelse(mc2 < 0, 'B', 'A')
    bedlist[[1]]$class[i] <- paste0(ifelse(mc1 < 0, 'B', 'A'),
                                    "_to_",
                                    ifelse(mc2 < 0, 'B', 'A'))
  }
  
  
  
  # output a list of bed-files per state and the model for other times.
  DIFF <- ifelse(abs(diff(unlist(lapply(dlist, mean)) )) > 0.25, T, F)
  if(!DIFF){
    message('The difference in |d| between the two states is low (< |.25|).
            Maybe consider this as a failed segmentation or no differences.')
  }
  return(list(BEDS = bedlist, DIFF = DIFF, ARM = ifelse(arm_start==0,'p','q')))
}

library(depmixS4)
################################################################################
################################    load data   ################################   
################################################################################
CS <- read.delim(h = F, '/DATA/references/human/hg19/ucsc.hg19.chrom.sizes')
cents <- read.delim('/DATA/references/human/hg19/cytobandAcen.bed', h = F)

ICE_40 = list.files('/DATA/projects/Hap1/Med12/hicpro/X201SC19114557-Z01-F001/', '40000_iced.matrix$', recursive = T, full.names = T)
ABS_40 = 'WT_40000_abs.bed'


contacts_40 = lapply(ICE_40[!grepl(ICE_40, pattern = 'days')], function(x){
  NAAM = reshape2::colsplit(dirname(x), pattern = "/", names = LETTERS)[,12]
  load_contacts(x, ABS_40, centromeres = cents,
                colour = MED12_PAL[reshape2::colsplit(NAAM,pattern = "_", names = LETTERS)[,1]],
                sample_name = NAAM)
})
names(contacts_40) = sapply(contacts_40, function(x) attr(x, 'samplename'))




K4mono_WT <- read.delim('4557_1_Hap1_K4mono_CCGTCC_S1_peaks.narrowPeak', h = F)
K4mono_MD <- read.delim('4590_6_3_K4mo_Med12_chipseq_CAGATC_S79_L007_peaks.narrowPeak', h = F)

################################################################################
###########################    compartment-score   #############################   
################################################################################

CSlist <- compartment_score(contacts_40, bed = K4mono_WT)

# calc CS
WT_40k_CS <- CSlist$compart_scores[,c(1,2,3,7)]
MD_40k_CS <- CSlist$compart_scores[,c(1,2,3,5)]


WT_40k_CS[,4][is.na(WT_40k_CS[,4])] <- 0
MD_40k_CS[,4][is.na(MD_40k_CS[,4])] <- 0

WT_40k_CS.bk <- WT_40k_CS
MD_40k_CS.bk <- MD_40k_CS

################################################################################
#################################    run HMM   #################################   
################################################################################

CHlist <- list()
for(C in levels(cents$V1)[-24]){
  for(A in c('p', 'q')){
    if(C %in% paste0('chr', c(9, 19, 15, 22))){ next()}
    message(paste(c(C, A), collapse = "_") )
    
    
    A_S <- 0
    A_E <- Inf
    if(A == 'p'){
      tmp = cents[cents$V1 == C,]
      A_S <- 0
      A_E <- tmp$V2
    } else {
      tmp = cents[cents$V1 == C,]
      A_S <- tmp$V3
      A_E <- CS[CS$V1 == C,2]
    }
    
    if (A_E/1e6 < 20) {
      next()
    }
    
    
    
    OUT <- compartmentHammer(CS_1 = WT_40k_CS, CS_2 = MD_40k_CS,
                             arm_start = A_S, arm_end = A_E,
                             chrom = C)
    CHlist[[ paste(c(C, A), collapse = "_")  ]] <- OUT$BEDS[[1]]
  }
}
dCDs <- unique(rbind_list(CHlist))

# output: [bed] class d Q
options(scipen = 9999999)
write.table('compartmentHMMout.bed',
            sep = "\t", quote = F, row.names = F, col.names = F,
            x = dCDs[, -c(6,7)])


HCD <- dCDs[dCDs$direction == 'hyper', ]

library(dplyr)
HCD$score = 0
for(C in unique(HCD$class)){
  tmp <- HCD[HCD$class == C , ]
  tmp$score <- as.numeric(cut(abs(tmp$d), quantile(abs(tmp$d), c(0, .25, .5, .75, 1)),
                              include.lowest = T))
  HCD[HCD$class == C , ] <- tmp
}

HCD$name <- HCD$class
HCD$strand <- "+"

################################################################################
##############################    filter classes   #############################   
################################################################################

twoClass_HCD  <- HCD[grepl(x = HCD$class, pattern = '_to_B') & HCD$direction == 'hyper',]

B <- RHWlib::df2gr(WT_40k_CS[WT_40k_CS$WT < 0, ])
B <- as.data.frame(reduce(B))
#library(RHWlib)
Bhcd <- as.data.frame(reduce(df2gr(HCD[grepl(x = HCD$class , pattern = 'B'),])))

options(scipen = 9999999)
write.table(twoClass_HCD[, c("seqnames", "start", "end", "class","score", "strand")],
            quote = F, sep = '\t', row.names = F, col.names = F,
            file = 'HCDs.bed')
write.table(B[, c("seqnames", "start", "end")],
            quote = F, sep = '\t', row.names = F, col.names = F,
            file = 'allB.bed')
write.table(Bhcd[, c("seqnames", "start", "end")],
            quote = F, sep = '\t', row.names = F, col.names = F,
            file = 'Bhcd.bed')

stop('Now intersect B with twoClass_HCD using bedtools intersect -v')
stop('bedtools intersect -a allB.bed -b Bhcd.bed -v > nHCD.bed; rm allB.bed; rm Bhcd.bed')


################################################################################
##############################    clean ranges   ###############################   
################################################################################
trimmedHCD <- NULL
treshD <- 1.5
fragmentationSize = (40e3)*5
for(i in 1:nrow(HCD)){
  hcd = HCD[i, ]
  
  score_WT <- mcols(WT_40k_CS[queryHits(findOverlaps(WT_40k_CS, df2gr(hcd)))])$score
  score_MD <- mcols(MD_40k_CS[queryHits(findOverlaps(MD_40k_CS, df2gr(hcd)))])$score
  
  D = abs(score_MD-score_WT)
  Treshhol = median(D) - (treshD*mad(D))

  dumpIDX <- which(D < Treshhol)
  
  keepIDX = which(!seq(1:length(D)) %in% dumpIDX)
  
  bedEntries <- seq(hcd$V2, hcd$V3, by = 40e3)
  
  if(length(unique(diff(keepIDX))) == 1){
    df = data.frame(hcd$V1, bedEntries[min(keepIDX)], bedEntries[max(keepIDX)])
    colnames(df) <- colnames(hcd[1:3])
    df <- df[(df[,3] - df[,2]) >= fragmentationSize, ]
    trimmedHCD = rbind(trimmedHCD, df)
  } else{
    df = data.frame(hcd$V1, bedEntries[ makeBED(keepIDX)[,1]], bedEntries[ makeBED(keepIDX)[,2]])
    colnames(df) <- colnames(hcd[1:3])
    df <- df[! (df[,3] - df[,2]) < fragmentationSize, ]
    trimmedHCD = rbind(trimmedHCD, df)
  }
  
}

write.table(trimmedHCD[complete.cases(trimmedHCD),],
            quote = F, sep = '\t', row.names = F, col.names = F,
            file = 'trimmedHCD.bed')