
getSequenceContext <- function(genome, position, chr, offsetL= 10, offsetR=50){
    sequence <- Biostrings::getSeq(genome, GRanges(seqnames=chr,ranges=IRanges(start=position-offsetL, end=position+offsetR)))
    context_string <- as.character(sequence)
    sequenceContext <- list(sequence=sequence,context_string=context_string)
}




attribute_sequence_contex_indel <- function(genome, in_CHROM, in_POS, 
						in_REF,in_ALT,
                                                 in_verbose = FALSE, 
                                                 in_offsetL= 10,
                                                 in_offsetR=50) {
  diff <- abs(nchar(in_REF)-nchar(in_ALT))
  out_dat <- list()
  out_dat$CHROM <- in_CHROM
  out_dat$POS <- in_POS
  out_dat$REF <- in_REF
  out_dat$ALT <- in_ALT

#  if (!is.character(in_REF)) {
#    print(out_dat)
#  }

  ref_length <- nchar(in_REF)
  alt_length <- nchar(in_ALT)
  if(ref_length > alt_length){
    out_dat$Type <- "Del"
    context_list <- getSequenceContext(genome = genome,
                                       position = in_POS,
                                       chr = in_CHROM,
                                       offsetL= in_offsetL,
                                       offsetR= in_offsetR)
    out_dat$SequenceContext <- context_list$context_string
   }else{
    out_dat$Type <- "Ins"
    context_list <- getSequenceContext(genome = genome,
                                       position = in_POS,
                                       chr = in_CHROM,
                                       offsetL= in_offsetL,
                                       offsetR= in_offsetR)
    out_dat$SequenceContext <- context_list$context_string
  }

  out_dat$Differnce <- diff
  out_dat$Change <- paste0(in_REF,">",in_ALT)
  return(out_dat)
}

attribution_of_indels <-function(genome, in_CHROM, in_POS,in_REF,in_ALT,
                                  in_verbose = FALSE,in_offsetL= 10,in_offsetR=50) { 

   in_dat_return <- attribute_sequence_contex_indel(genome=genome, 
						in_CHROM=in_CHROM, 
						in_POS=in_POS,
                                                in_REF=in_REF,
						in_ALT=in_ALT,
                                                in_offsetL=in_offsetL,
                                                in_offsetR=in_offsetR)
 
    sequence_string <- DNAString(in_dat_return$SequenceContext)
    if(in_dat_return$Differnce+11 > 61){
      motive_end_index <- 60
    }else{
    motive_end_index <- in_dat_return$Differnce+11
    }
    rest_start_index <- motive_end_index+1
    motive<-sequence_string[12:motive_end_index]
    motive_length <- length(motive)
    rest_string_R <-subseq(sequence_string, start= rest_start_index)
    rest_string_L <-subseq(sequence_string, start=1, end=11)
    match <- matchPattern(motive, rest_string_R)
    
    if(length(width(match))<5){
      rest_string_R_new <-xscat(rest_string_R, rep(motive, 5))
      match <- matchPattern(motive, rest_string_R_new)
    }
    if(in_dat_return$Type == "Del"){
      if(in_dat_return$Differnce == 1){
        if(identical(as.character(match),character(0))| start(match)[1]>1){
          if(as.character(motive)== "T" | as.character(motive) == "A"){
            in_dat_return$Motiv <- "T:A"
            in_dat_return$RepeatSize <- 0
            in_dat_return$IndelNumber <- "DEL_T_1_0"
          }else{
            in_dat_return$Motiv <- "C:G"
            in_dat_return$RepeatSize <- 0
            in_dat_return$IndelNumber <- "DEL_C_1_0"  
          }
        }else if(start(match)[1]==1 && start(match)[2] == 1+motive_length && 
                 start(match)[3] == 1+(2*motive_length) && 
                 start(match)[4] == 1+(3*motive_length) && 
                 start(match)[5] == 1+(4*motive_length)){
          if(as.character(motive)== "T" | as.character(motive) == "A"){
            in_dat_return$Motiv <- "T:A"
            in_dat_return$RepeatSize <- "5+"
            in_dat_return$IndelNumber <- "DEL_T_1_5+"
          }else{
            in_dat_return$Motiv <- "C:G"
            in_dat_return$RepeatSize <- "5+"
            in_dat_return$IndelNumber <- "DEL_C_1_5+"
          }
        }else if(start(match)[1]==1 && start(match)[2] == 1+motive_length && 
                 start(match)[3] == 1+(2*motive_length) && 
                 start(match)[4] == 1+(3*motive_length)){
          if(as.character(motive)== "T" | as.character(motive) == "A"){
            in_dat_return$Motiv <- "T:A"
            in_dat_return$RepeatSize <- 4
            in_dat_return$IndelNumber <- "DEL_T_1_4"
          }else{
            in_dat_return$Motiv <- "C:G"
            in_dat_return$RepeatSize <- 4
            in_dat_return$IndelNumber <- "DEL_C_1_4"
          }
        }else if(start(match)[1]==1 && start(match)[2] == 1+motive_length && 
                 start(match)[3] == 1+(2*motive_length)){
          if(as.character(motive)== "T" | as.character(motive) == "A"){
            in_dat_return$Motiv <- "T:A"
            in_dat_return$RepeatSize <- 3
            in_dat_return$IndelNumber <- "DEL_T_1_3"
          }else{
            in_dat_return$Motiv <- "C:G"
            in_dat_return$RepeatSize <- 3
            in_dat_return$IndelNumber <- "DEL_C_1_3"
          }
        }else if(start(match)[1]==1 && start(match)[2] == 1+motive_length){
          if(as.character(motive)== "T" | as.character(motive) == "A"){
            in_dat_return$Motiv <- "T:A"
            in_dat_return$RepeatSize <- 2
            in_dat_return$IndelNumber <- "DEL_T_1_2"
          }else{
            in_dat_return$Motiv <- "C:G"
            in_dat_return$RepeatSize <- 2
            in_dat_return$IndelNumber <- "DEL_C_1_2"
          }
        }else if(start(match)[1]==1){
          if(as.character(motive)== "T" | as.character(motive) == "A"){
            in_dat_return$Motiv <- "T:A"
            in_dat_return$RepeatSize <- 1
            in_dat_return$IndelNumber <- "DEL_T_1_1"
          }else{
            in_dat_return$Motiv <- "C:G"
            in_dat_return$RepeatSize <- 1
            in_dat_return$IndelNumber <- "DEL_C_1_1"
            }
        }else{
          in_dat_return$Motiv <- "NA"
          in_dat_return$RepeatSize <- "NA"
          in_dat_return$IndelNumber <- "NOT DEFINED"
          }
      }else if(in_dat_return$Differnce == 2){
        homologR1 <- motive[1:motive_length-1]
        homologL1 <- motive[2:motive_length]
        
        match_homologR1 <- matchPattern(homologR1, rest_string_R)
        match_homologL1 <- matchPattern(homologL1, rest_string_L)
        
        if(((length(as.character(match_homologL1)) != 0 |
             length(as.character(match_homologR1)) != 0) && 
                 (start(match)[1]!=1 | identical(as.character(match),
                                             character(0)))) &&
                 ((start(match_homologR1)[1]==1 && 
                   !is.na(start(match_homologR1)[1]==1))|
                 (end(match_homologL1)[1]==11 && 
                  !is.na(end(match_homologL1)[1]==11)))){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- "1bp"
          in_dat_return$IndelNumber <- "DEL_MH_2_1"
        }else if((identical(as.character(match),character(0)) && 
                  is.na(start(match)[1]!=1)) |  start(match)[1]!=1 |
                 identical(as.character(match_homologL1),character(0)) && 
             identical(as.character(match_homologR1),character(0))){
           in_dat_return$Motiv <- as.character(motive)
           in_dat_return$RepeatSize <- 0
           in_dat_return$IndelNumber <-"DEL_repeats_2_0"
        }else if(length(match) >= 5 && 
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length && 
                  start(match)[3] == 1+(2*motive_length) && 
                  start(match)[4] == 1+(3*motive_length) && 
                  start(match)[5] == 1+(4*motive_length))){
           in_dat_return$Motiv <- as.character(motive)
           in_dat_return$RepeatSize <- "5+"
           in_dat_return$IndelNumber <- "DEL_repeats_2_5+"
        }else if(length(match) >= 4 && 
                 (start(match)[1]==1 &&
                  start(match)[2] == 1+motive_length && 
                  start(match)[3] == 1+(2*motive_length) && 
             start(match)[4] == 1+(3*motive_length))){
           in_dat_return$Motiv <- as.character(motive)
           in_dat_return$RepeatSize <- 4
           in_dat_return$IndelNumber <- "DEL_repeats_2_4"
        }else if(length(match) >= 3 && 
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length &&
                  start(match)[3] == 1+(2*motive_length))){
           in_dat_return$Motiv <- as.character(motive)
           in_dat_return$RepeatSize <- 3
           in_dat_return$IndelNumber <- "DEL_repeats_2_3"
        }else if(length(match) >= 2 && 
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length)){
           in_dat_return$Motiv <- as.character(motive)
           in_dat_return$RepeatSize <- 2
           in_dat_return$IndelNumber <- "DEL_repeats_2_2"
        }else if(length(match) >= 1 && start(match)[1]==1){
           in_dat_return$Motiv <- as.character(motive)
           in_dat_return$RepeatSize <- 1
           in_dat_return$IndelNumber <- "DEL_repeats_2_1"
        }else{
           in_dat_return$Motiv <- "NA"
           in_dat_return$RepeatSize <- "NA"
           in_dat_return$IndelNumber <- "NOT DEFINED"
        }
      }else if(in_dat_return$Differnce == 3){
        homologR1 <- motive[1:motive_length-1]
        homologL1 <- motive[2:motive_length]
        homologR2 <- homologR1[1:length(homologR1)-1]
        homologL2 <- homologL1[2:length(homologL1)]
       
        match_homologR1 <- matchPattern(homologR1, rest_string_R)
        match_homologL1 <- matchPattern(homologL1, rest_string_L)
        match_homologR2 <- matchPattern(homologR2, rest_string_R)
        match_homologL2 <- matchPattern(homologL2, rest_string_L)
        
  
        if(((length(as.character(match_homologL1)) != 0 | 
             length(as.character(match_homologR1)) != 0) && 
             (start(match)[1]!=1 |
              identical(as.character(match),character(0)))) && 
                 ((start(match_homologR1)[1]==1 && 
                   !is.na(start(match_homologR1)[1]==1))|
                 (end(match_homologL1)[1]==11 && 
                  !is.na(end(match_homologL1)[1]==11)))){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- "2bp"
          in_dat_return$IndelNumber <- "DEL_MH_3_2" 
        }else if(((length(as.character(match_homologL2)) != 0 | 
                   length(as.character(match_homologR2)) != 0) && 
                 (start(match)[1]!=1 | 
                  identical(as.character(match),character(0)))) &&
                 (start(match_homologR2)[1]==1 && 
                  !is.na(start(match_homologR2)[1]==1))|
                 (end(match_homologL2)[1]==11 && 
                  !is.na(end(match_homologL2)[1]==11))){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- "1bp"
          in_dat_return$IndelNumber <- "DEL_MH_3_1" 
        }else if((identical(as.character(match),character(0)) && 
                  is.na(start(match)[1]!=1))| start(match)[1]!=1 | 
                 (identical(as.character(match_homologL1),character(0)) && 
                  identical(as.character(match_homologR1),character(0)) && 
                  identical(as.character(match_homologL2),character(0)) && 
                  identical(as.character(match_homologR2),character(0)))){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- 0
          in_dat_return$IndelNumber <- "DEL_repeats_3_0"
        }else if(length(match) >= 5 && 
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length && 
                  start(match)[3] == 1+(2*motive_length) &&
                  start(match)[4] == 1+(3*motive_length) && 
                  start(match)[5] == 1+(4*motive_length))){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- "5+"
          in_dat_return$IndelNumber <- "DEL_repeats_3_5+"
        }else if(length(match) >= 4 && 
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length && 
                  start(match)[3] == 1+(2*motive_length) &&
                  start(match)[4] == 1+(3*motive_length))){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- 4
          in_dat_return$IndelNumber <- "DEL_repeats_3_4"
        }else if(length(match) >= 3 && 
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length &&
                  start(match)[3] == 1+(2*motive_length))){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- 3
          in_dat_return$IndelNumber <- "DEL_repeats_3_3"
        }else if(length(match) >= 2 && 
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length)){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- 2
          in_dat_return$IndelNumber <- "DEL_repeats_3_2"
        }else if(length(match) >= 1 && start(match)[1]==1){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- 1
          in_dat_return$IndelNumber <- "DEL_repeats_3_1"
        }else{
          in_dat_return$Motiv <- "NA"
          in_dat_return$RepeatSize <- "NA"
          in_dat_return$IndelNumber <- "NOT DEFINED"
        }
      }else if(in_dat_return$Differnce == 4){
        homologR1 <- motive[1:motive_length-1]
        homologL1 <- motive[2:motive_length]
        homologR2 <- homologR1[1:length(homologR1)-1]
        homologL2 <- homologL1[2:length(homologL1)]
        homologR3 <- homologR2[1:length(homologR2)-1]
        homologL3 <- homologL2[2:length(homologL2)]
      
        
        match_homologR1 <- matchPattern(homologR1, rest_string_R)
        match_homologL1 <- matchPattern(homologL1, rest_string_L)
        match_homologR2 <- matchPattern(homologR2, rest_string_R)
        match_homologL2 <- matchPattern(homologL2, rest_string_L)
        match_homologR3 <- matchPattern(homologR3, rest_string_R)
        match_homologL3 <- matchPattern(homologL3, rest_string_L)
        
        
        if(((length(as.character(match_homologL1)) != 0 |
             length(as.character(match_homologR1)) != 0) &&
                 (start(match)[1]!=1 | 
                  identical(as.character(match),character(0)))) && 
                 ((start(match_homologR1)[1]==1 && 
                   !is.na(start(match_homologR1)[1]==1))|
                 (end(match_homologL1)[1]==11 && 
                  !is.na(end(match_homologL1)[1]==11)))){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- "3bp"
          in_dat_return$IndelNumber <- "DEL_MH_4_3"
        }else if(((length(as.character(match_homologL2)) != 0 | 
                   length(as.character(match_homologR2)) != 0) &&
                 (start(match)[1]!=1 | 
                  identical(as.character(match),character(0)))) && 
                 ((start(match_homologR2)[1]==1 &&
                   !is.na(start(match_homologR2)[1]==1))|
                 (end(match_homologL2)[1]==11 && 
                  !is.na(end(match_homologL2)[1]==11)))){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- "2bp"
          in_dat_return$IndelNumber <- "DEL_MH_4_2"
        }else if(((length(as.character(match_homologL3)) != 0 | 
                   length(as.character(match_homologR3)) != 0) &&
                 (start(match)[1]!=1 | 
                  identical(as.character(match),character(0)))) && 
                 ((start(match_homologR3)[1]==1 && 
                   !is.na(start(match_homologR3)[1]==1))|
                 (end(match_homologL3)[1]==11 && 
                  !is.na(end(match_homologL3)[1]==11)))){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- "1bp"
          in_dat_return$IndelNumber <- "DEL_MH_4_1"
        }else if((identical(as.character(match),character(0)) && 
                  is.na(start(match)[1]!=1)) |  start(match)[1]!=1 |
                 (identical(as.character(match_homologL1),character(0)) && 
                  identical(as.character(match_homologR1),character(0)) && 
                  identical(as.character(match_homologL2),character(0)) &&
                  identical(as.character(match_homologR2),character(0)) &&
                  identical(as.character(match_homologL3),character(0)) && 
                  identical(as.character(match_homologR3),character(0)))){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- 0
          in_dat_return$IndelNumber <- "DEL_repeats_4_0" 
        }else if(length(match) >= 5 && 
                 (start(match)[1] == 1 && 
                  start(match)[2] == 1+motive_length && 
                 start(match)[3] == 1+(2*motive_length) && 
                 start(match)[4] == 1+(3*motive_length) && 
                 start(match)[5] == 1+(4*motive_length))){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- "5+"
          in_dat_return$IndelNumber <- "DEL_repeats_4_5+" 
        }else if(length(match) >= 4 &&
                 (start(match)[1] == 1 &&
                  start(match)[2] == 1+motive_length && 
                  start(match)[3] == 1+(2*motive_length) &&
                  start(match)[4] == 1+(3*motive_length))){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- 4
          in_dat_return$IndelNumber <- "DEL_repeats_4_4" 
        }else if(length(match) >= 3 && 
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length && 
                  start(match)[3] == 1+(2*motive_length))){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- 3
          in_dat_return$IndelNumber <- "DEL_repeats_4_3" 
        }else if(length(match) >= 2 &&
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length)){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- 2
          in_dat_return$IndelNumber <- "DEL_repeats_4_2" 
        }else if(length(match) >= 1 && start(match)[1]==1){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- 1
          in_dat_return$IndelNumber <- "DEL_repeats_4_1" 
        }else{
          in_dat_return$Motiv <- "NA"
          in_dat_return$RepeatSize <- "NA"
          in_dat_return$IndelNumber <- "NOT DEFINED"
        }
      }else if(in_dat_return$Differnce >= 5){
        
           motive_first_5bp <- motive[1:5]
           start_last_5bp <- (motive_length-5)+1
           motive_last_5bp <- motive[start_last_5bp:motive_length]
           
           homologR1 <- motive_first_5bp
           homologL1 <- motive_last_5bp
           homologR2 <- homologR1[1:length(homologR1)-1]
           homologL2 <- homologL1[2:length(homologL1)]
           homologR3 <- homologR2[1:length(homologR2)-1]
           homologL3 <- homologL2[2:length(homologL2)]
           homologR4 <- homologR3[1:length(homologR3)-1]
           homologL4 <- homologL3[2:length(homologL3)]
           homologR5 <- homologR4[1:length(homologR4)-1]
           homologL5 <- homologL4[2:length(homologL4)]
           
           
           match_homologR1 <- matchPattern(homologR1, rest_string_R)
           match_homologL1 <- matchPattern(homologL1, rest_string_L)
           match_homologR2 <- matchPattern(homologR2, rest_string_R)
           match_homologL2 <- matchPattern(homologL2, rest_string_L)
           match_homologR3 <- matchPattern(homologR3, rest_string_R)
           match_homologL3 <- matchPattern(homologL3, rest_string_L)
           match_homologR4 <- matchPattern(homologR4, rest_string_R)
           match_homologL4 <- matchPattern(homologL4, rest_string_L)
           match_homologR5 <- matchPattern(homologR5, rest_string_R)
           match_homologL5 <- matchPattern(homologL5, rest_string_L)
           
         if(length(match) >= 5 && 
            (start(match)[1]==1 &&
             start(match)[2] == 1+motive_length &&
             start(match)[3] == 1+(2*motive_length) &&
             start(match)[4] == 1+(3*motive_length) &&
             start(match)[5] == 1+(4*motive_length))){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- "5+"
          in_dat_return$IndelNumber <- "DEL_repeats_5+_5+"
        }else if(length(match) >= 4 && 
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length && 
                  start(match)[3] == 1+(2*motive_length) && 
                  start(match)[4] == 1+(3*motive_length))){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- 4
          in_dat_return$IndelNumber <- "DEL_repeats_5+_4"
        }else if(length(match) >= 3 &&
                 (start(match)[1]==1 &&
                 start(match)[2] == 1+motive_length &&
                 start(match)[3] == 1+(2*motive_length))){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- 3
          in_dat_return$IndelNumber <- "DEL_repeats_5+_3"
        }else if(length(match) >= 2 && 
                 (start(match)[1]==1 && 
                 start(match)[2] == 1+motive_length)){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- 3
          in_dat_return$IndelNumber <- "DEL_repeats_5+_2"
        }else if(length(match) >= 1 && start(match)[1]==1){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- 1
          in_dat_return$IndelNumber <- "DEL_repeats_5+_1"
        }else if(((length(as.character(match_homologL1)) != 0 | 
                   length(as.character(match_homologR1)) != 0) && 
                 (start(match)[1]!=1 | 
                  identical(as.character(match),character(0)))) &&
                ((start(match_homologR1)[1]==1 && 
                  !is.na(start(match_homologR1)[1]==1))|
                 (end(match_homologL1)[1]==11 && 
                  !is.na(end(match_homologL1)[1]==11)))){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- "5bp"
          in_dat_return$IndelNumber <- "DEL_MH_5+_5+"
        }else if(((length(as.character(match_homologL2)) != 0 | 
                   length(as.character(match_homologR2)) != 0) && 
                  (start(match)[1]!=1 | 
                   identical(as.character(match),character(0)))) &&
                 ((start(match_homologR2)[1]==1 && 
                   !is.na(start(match_homologR2)[1]==1))|
                  (end(match_homologL2)[1]==11 && 
                   !is.na(end(match_homologL2)[1]==11)))){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- "4bp"
          in_dat_return$IndelNumber <- "DEL_MH_5+_4"
        }else if(((length(as.character(match_homologL3)) != 0 |
                   length(as.character(match_homologR3)) != 0) &&
                  (start(match)[1]!=1 | 
                   identical(as.character(match),character(0))))&&
                 ((start(match_homologR3)[1]==1 &&
                   !is.na(start(match_homologR3)[1]==1))|
                  (end(match_homologL3)[1]==11 &&
                   !is.na(end(match_homologL3)[1]==11)))){
            in_dat_return$Motiv <- as.character(motive)
            in_dat_return$RepeatSize <- "3bp"
            in_dat_return$IndelNumber <- "DEL_MH_5+_3"
        }else if(((length(as.character(match_homologL4)) != 0 |
                   length(as.character(match_homologR4)) != 0) && 
                  (start(match)[1]!=1 | 
                   identical(as.character(match),character(0))))&&
                 ((start(match_homologR4)[1]==1 && 
                   !is.na(start(match_homologR4)[1]==1))|
                  (end(match_homologL4)[1]==11 &&
                   !is.na(end(match_homologL4)[1]==11)))){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- "2bp"
          in_dat_return$IndelNumber <- "DEL_MH_5+_2"
        }else if(((length(as.character(match_homologL5)) != 0 |
                   length(as.character(match_homologR5)) != 0) &&
                  (start(match)[1]!=1 | 
                   identical(as.character(match),character(0)))) &&
                 ((start(match_homologR5)[1]==1 && 
                   !is.na(start(match_homologR5)[1]==1))|
                  (end(match_homologL5)[1]==11 &&
                   !is.na(end(match_homologL5)[1]==11)))){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- "1bp"
          in_dat_return$IndelNumber <- "DEL_MH_5+_1"
        }else if(((identical(as.character(match),character(0)) && 
                   is.na(start(match)[1]!=1)) |  start(match)[1]!=1 | 
                  identical(as.character(match_homologL1),character(0)) &&
                  identical(as.character(match_homologR1),character(0)) &&
                  identical(as.character(match_homologL2),character(0)) && 
                  identical(as.character(match_homologR2),character(0)) && 
                  identical(as.character(match_homologL3),character(0)) && 
                  identical(as.character(match_homologR3),character(0)) &&
                  identical(as.character(match_homologL4),character(0)) && 
                  identical(as.character(match_homologR4),character(0))) &&
                  identical(as.character(match_homologL5),character(0)) && 
                 identical(as.character(match_homologR5),character(0))
                 ){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- 0
          in_dat_return$IndelNumber <- "DEL_repeats_5+_0"
        }else{
          in_dat_return$Motiv <- "NA"
          in_dat_return$RepeatSize <- "NA"
          in_dat_return$IndelNumber <- "NOT DEFINED"
        }
      }
    }else{
      if(in_dat_return$Differnce == 1){
        if(start(match)[1]==1 && 
           start(match)[2] == 1+motive_length && 
           start(match)[3] == 1+(2*motive_length) && 
           start(match)[4] == 1+(3*motive_length) && 
           start(match)[5] == 1+(4*motive_length)){
          if(as.character(motive)== "T" | as.character(motive) == "A"){
            in_dat_return$Motiv <- "T:A"
            in_dat_return$RepeatSize <- "5+"
            in_dat_return$IndelNumber <- "INS_T_1_5+"
          }else{
            in_dat_return$Motiv <- "C:G"
            in_dat_return$RepeatSize <- "5+"
            in_dat_return$IndelNumber <- "INS_C_1_5+"
          }
        }else if(start(match)[1]==1 && 
                 start(match)[2] == 1+motive_length &&
                 start(match)[3] == 1+(2*motive_length) && 
                 start(match)[4] == 1+(3*motive_length)){
          if(as.character(motive)== "T" | as.character(motive) == "A"){
            in_dat_return$Motiv <- "T:A"
            in_dat_return$RepeatSize <- 4
            in_dat_return$IndelNumber <- "INS_T_1_4"
          }else{
            in_dat_return$Motiv <- "C:G"
            in_dat_return$RepeatSize <- 4
            in_dat_return$IndelNumber <- "INS_C_1_4"
          }
        }else if(start(match)[1]==1 && 
                 start(match)[2] == 1+motive_length && 
                 start(match)[3] == 1+(2*motive_length)){
          if(as.character(motive)== "T" | as.character(motive) == "A"){
            in_dat_return$Motiv <- "T:A"
            in_dat_return$RepeatSize <- 3
            in_dat_return$IndelNumber <- "INS_T_1_3"
          }else{
            in_dat_return$Motiv <- "C:G"
            in_dat_return$RepeatSize <- 3
            in_dat_return$IndelNumber <- "INS_C_1_3"
          }
        }else if(start(match)[1]==1 && start(match)[2] == 1+motive_length){
          if(as.character(motive)== "T" | as.character(motive) == "A"){
            in_dat_return$Motiv <- "T:A"
            in_dat_return$RepeatSize <- 2
            in_dat_return$IndelNumber <- "INS_T_1_2"
          }else{
            in_dat_return$Motiv <- "C:G"
            in_dat_return$RepeatSize <- 2
            in_dat_return$IndelNumber <- "INS_C_1_2"
          }
        }else if(start(match)[1]==1){
          if(as.character(motive)== "T" | as.character(motive) == "A"){
            in_dat_return$Motiv <- "T:A"
            in_dat_return$RepeatSize <- 1
            in_dat_return$IndelNumber <- "INS_T_1_1"
          }else{
            in_dat_return$Motiv <- "C:G"
            in_dat_return$RepeatSize <- 1
            in_dat_return$IndelNumber <- "INS_C_1_1"
          }
        }else{
          if(as.character(motive)== "T" | as.character(motive) == "A"){
            in_dat_return$Motiv <- "T:A"
            in_dat_return$RepeatSize <- 0
            in_dat_return$IndelNumber <- "INS_T_1_0"
          }else{
            in_dat_return$Motiv <- "C:G"
            in_dat_return$RepeatSize <- 0
            in_dat_return$IndelNumber <- "INS_C_1_0"  
          }
        }
      }else if(in_dat_return$Differnce == 2){
        if (identical(as.character(match),character(0)) | start(match)[1]!=1){
            in_dat_return$Motiv <- as.character(motive)
            in_dat_return$RepeatSize <- 0
            in_dat_return$IndelNumber <- "INS_repeats_2_0"
        }else if(length(match) >= 5 &&
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length && 
                  start(match)[3] == 1+(2*motive_length) &&
               start(match)[4] == 1+(3*motive_length) &&
               start(match)[5] == 1+(4*motive_length))){
            in_dat_return$Motiv <- as.character(motive)
            in_dat_return$RepeatSize <- "5+"
            in_dat_return$IndelNumber <- "INS_repeats_2_5+"
        }else if(length(match) >= 4 && 
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length &&
                  start(match)[3] == 1+(2*motive_length) &&
                start(match)[4] == 1+(3*motive_length))){
            in_dat_return$Motiv <- as.character(motive)
            in_dat_return$RepeatSize <- 4
            in_dat_return$IndelNumber <- "INS_repeats_2_4"
        }else if(length(match) >= 3 &&
                 (start(match)[1]==1 &&
                  start(match)[2] == 1+motive_length &&
                  start(match)[3] == 1+(2*motive_length))){
            in_dat_return$Motiv <- as.character(motive)
            in_dat_return$RepeatSize <- 3
            in_dat_return$IndelNumber <- "INS_repeats_2_3"
        }else if(length(match) >= 2 &&
                 (start(match)[1]==1 &&
                  start(match)[2] == 1+motive_length)){
            in_dat_return$Motiv <- as.character(motive)
            in_dat_return$RepeatSize <- 2
            in_dat_return$IndelNumber <- "INS_repeats_2_2"
        }else if(length(match) >= 1 && start(match)[1]==1){
            in_dat_return$Motiv <- as.character(motive)
            in_dat_return$RepeatSize <- 1
            in_dat_return$IndelNumber <- "INS_repeats_2_1"
        }else{
            in_dat_return$Motiv <- "NA"
            in_dat_return$RepeatSize <- "NA"
            in_dat_return$IndelNumber <- "NOT DEFINED"
        }
      }else if(in_dat_return$Differnce == 3){
        if (identical(as.character(match),character(0))| start(match)[1]!=1){
            in_dat_return$Motiv <- as.character(motive)
            in_dat_return$RepeatSize <- 0
            in_dat_return$IndelNumber <- "INS_repeats_3_0"
        }else if(length(match) >= 5 && 
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length &&
                  start(match)[3] == 1+(2*motive_length) &&
                  start(match)[4] == 1+(3*motive_length) && 
                  start(match)[5] == 1+(4*motive_length))){
            in_dat_return$Motiv <- as.character(motive)
            in_dat_return$RepeatSize <- "5+"
            in_dat_return$IndelNumber <- "INS_repeats_3_5+"
        }else if (length(match) >= 4 && 
                  (start(match)[1]==1 && 
                   start(match)[2] == 1+motive_length && 
                   start(match)[3] == 1+(2*motive_length) &&
                start(match)[4] == 1+(3*motive_length))){
            in_dat_return$Motiv <- as.character(motive)
            in_dat_return$RepeatSize <- 4
            in_dat_return$IndelNumber <- "INS_repeats_3_4"
        }else if(length(match) >= 3 && 
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length && 
                  start(match)[3] == 1+(2*motive_length))){
            in_dat_return$Motiv <- as.character(motive)
            in_dat_return$RepeatSize <- 3
            in_dat_return$IndelNumber <- "INS_repeats_3_3"
        }else if(length(match) >= 2 &&
                 (start(match)[1]==1 &&
                  start(match)[2] == 1+motive_length)){
            in_dat_return$Motiv <- as.character(motive)
            in_dat_return$RepeatSize <- 2
            in_dat_return$IndelNumber <- "INS_repeats_3_2"
        }else if(length(match) >= 1 && start(match)[1]==1){
            in_dat_return$Motiv <- as.character(motive)
            in_dat_return$RepeatSize <- 1
            in_dat_return$IndelNumber <- "INS_repeats_3_1"
        }else{
            in_dat_return$Motiv <- "NA"
            in_dat_return$RepeatSize <- "NA"
            in_dat_return$IndelNumber <- "NOT DEFINED"
        }
        }else if(in_dat_return$Differnce == 4){
        if (identical(as.character(match),character(0)) | start(match)[1]!=1){
            in_dat_return$Motiv <- as.character(motive)
            in_dat_return$RepeatSize <- 0
            in_dat_return$IndelNumber <- "INS_repeats_4_0"
        }else if(length(match) >= 5 && 
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length &&
                  start(match)[3] == 1+(2*motive_length) &&
               start(match)[4] == 1+(3*motive_length) &&
               start(match)[5] == 1+(4*motive_length))){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- "5+"
          in_dat_return$IndelNumber <- "INS_repeats_4_5+"
        }else if(length(match) >= 4 &&
                 (start(match)[1]==1 &&
                  start(match)[2] == 1+motive_length && 
                  start(match)[3] == 1+(2*motive_length) &&
                start(match)[4] == 1+(3*motive_length))){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- 4
          in_dat_return$IndelNumber <- "INS_repeats_4_4"
        }else if(length(match) >= 3 &&
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length &&
                  start(match)[3] == 1+(2*motive_length))){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- 3
          in_dat_return$IndelNumber <- "INS_repeats_4_3"
        }else if(length(match) >= 2 &&
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length)){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- 2
          in_dat_return$IndelNumber <- "INS_repeats_4_2"
        }else if(length(match) >= 1 && start(match)[1]==1){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- 1
          in_dat_return$IndelNumber <- "INS_repeats_4_1"
        }else{
          in_dat_return$Motiv <- "NA"
          in_dat_return$RepeatSize <- "NA"
          in_dat_return$IndelNumber <- "NOT DEFINED"
          }
      }else if(in_dat_return$Differnce >= 5){
        if (identical(as.character(match),character(0))| start(match)[1]!=1){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- 0
          in_dat_return$IndelNumber <- "INS_repeats_5+_0"
        }else if(length(match) >= 5 && 
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length && 
                  start(match)[3] == 1+(2*motive_length) &&
               start(match)[4] == 1+(3*motive_length) && 
               start(match)[5] == 1+(4*motive_length))){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- "5+"
          in_dat_return$IndelNumber <- "INS_repeats_5+_5+"
        }else if(length(match) >= 4 &&
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length && 
                  start(match)[3] == 1+(2*motive_length) &&
                start(match)[4] == 1+(3*motive_length))){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- 4
          in_dat_return$IndelNumber <- "INS_repeats_5+_4"
        }else if(length(match) >= 3 &&
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length && 
                  start(match)[3] == 1+(2*motive_length))){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- 3
          in_dat_return$IndelNumber <- "INS_repeats_5+_3"
        }else if(length(match) >= 2 && 
                 (start(match)[1]==1 && 
                  start(match)[2] == 1+motive_length)){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- 2
          in_dat_return$IndelNumber <- "INS_repeats_5+_2"
        }else if(length(match) >= 1 && start(match)[1]==1){
          in_dat_return$Motiv <- as.character(motive)
          in_dat_return$RepeatSize <- 1
          in_dat_return$IndelNumber <- "INS_repeats_5+_1"
        }else{
          in_dat_return$Motiv <- "NA"
          in_dat_return$RepeatSize <- "NA"
          in_dat_return$IndelNumber <- "NOT DEFINED"
          }
        }
      }   

      in_dat_return[['IndelNumber']]

}
