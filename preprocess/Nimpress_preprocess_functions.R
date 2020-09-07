######################################
## Nimpress preprocessing functions ##
######################################

## exit without error 
stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

## Not in function
'%!in%' <- function(x,y)!('%in%'(x,y));

## for output without colnames 
rename_cols <- function(df){
  for (i in seq_along(df)) {
    colnames(df)[i] = paste0('V', i)
  }
  return(df)
}


##################
## load bedfile ##
##################


bedfile_to_Granges <- function(ovlp){
  ovlp.df <- as.data.frame(ovlp)
  colnames(ovlp.df) <- c("chr", "start", "end")
  gr <- makeGRangesFromDataFrame(ovlp.df, keep.extra.columns = TRUE,starts.in.df.are.0based=TRUE)
}


###########################
## processing input file ##
###########################

check_gwas_file <- function(input){
  
  gwas_file <- read.table(input, sep=",", header=T, colClasses=c("character","character","numeric", "numeric"))
  
  ## remove blank rows
  blank <- which(gwas_file[,1] == "")
  na <- which(is.na(gwas_file[,1]))
  rm_blank <- c(blank, na)
  if(length(rm_blank) > 0){
    gwas_file <- gwas_file[-rm_blank,]
  }
  
  ## remove blank columns
  if(ncol(gwas_file) > 5){
    gwas_file <- gwas_file[,-(6:ncol(gwas_file))]
  }
  
  ## remove all leading and trailing blank spaces
  gwas_file <- as.data.frame(apply(gwas_file,2,function(x)gsub('\\s+', '',x)))
  
  ## check if all rsIDs are in propper format
  rsid_pattern <- "^rs\\d{1,}"
  sw <- grepl(rsid_pattern,gwas_file$rsID, ignore.case = T)
  sww <- which(sw !=TRUE)
  if(length(sww) > 0){
    stop(paste("Invalid format for", paste(gwas_file$rsID[sww],collapse=", ") , " on line ", sww, sep=""))
  }
  
  ## check ik all rsIDs are unique
  dups <- gwas_file$rsID[duplicated(gwas_file$rsID)]
  if(length(dups) > 0){
    stop(paste("The following rsIDs are duplicated ", paste(unique(dups),collapse=", "), sep=""))
  }
  
  ## check if valid nucleotides
  nva <- which(gwas_file$Risk_allele %!in% c("A", "T", "G", "C"))
  if(length(nva) > 0){
    stop(paste("Line ", paste(unique(nva), collapse=","), " does not contain a valid risk allele", sep=""))
  }
  
  ## Format Effect size
  or.ind <- grep("OR", colnames(gwas_file))
  beta.ind <- grep("Beta", colnames(gwas_file))
  
  if(length(or.ind) == 1){
    OR <- as.vector(gwas_file[,or.ind])
    gwas_file[,or.ind] <- log(as.numeric(sub('\\(.*', '', OR)))
    colnames(gwas_file)[or.ind] <- "Effect.size"
    message("Coverting OR to BETA")
  }else if(length(beta.ind) == 1){
    colnames(gwas_file)[beta.ind] <- "Effect.size"
  }
  
  message("All rsID ok")
  message("All Risk_allele ok")
  return(gwas_file)
}  



##########################################
## get rdID genomic location from dbSNP ##
##########################################

## extracts information from rsID dbSNP lookup: CHR, start, REF.Allele, ALT.Allele
format_dbSNP <- function(x){
  ## extract assembly, genome position and variant details 
  g1 <- gsub("^[A-Za-z]*=", "", x)
  ## remove everything after "|"
  g2 <- strsplit(gsub("\\|.*$", "", g1), ",")
  ## split chr from location
  g3 <- unlist(strsplit(g2[[1]], ":"))
  ## format into df
  odd <- seq_along(g3) %% 2 == 1
  NC_CHR <- g3[odd]
  nt <- gsub("^g\\.", "", g3[!odd])
  START <- gsub("[A-Z]>[A-Z]", "", nt)
  nt2 <- unlist(lapply(strsplit(nt,"[0-9]"),function(x) x[length(x)]))
  nt3 <- do.call(rbind, strsplit(nt2, ">"))
  colnames(nt3) <- c("REF", "ALT")
  df <- data.frame(NC_CHR,START, nt3)
  ## get the right genome build info based on assembly defined in main script 
  df2 <- unique(merge(assembly,df,by="NC_CHR"))
  df3 <- df2[,-1]
  
  ## if multiple alt alleles
  if(nrow(df3) > 1){
    df4 <- data.frame(CHR=unique(df3$CHR), START=unique(df3$START), REF=unique(df3$REF),ALT=paste(df3$ALT, collapse="|"), stringsAsFactors = F)
  }else{
    df4 <- df3
  }
  return(df4)
}

## this function lookes up information in the dbSNP database for each rsID
## Multi Alt Alleles are collapsed so it works better with downstream functions 

getrsID_info <- function(rsid_input){
  snp_term <- paste(rsid_input, "[RS]", sep="")
  r_search <- entrez_search(db="snp", term=snp_term)
  ## if no result returned
  if(length(r_search$id) == 0){
    final_snp <- cbind(rsid_input, NA,NA,NA,NA)
    colnames(final_snp) <- c("rsID", "CHR", "START", "REF.ALLELE", "ALT.ALLELE")
  }else{
    ## info associated with snp_id
    multi_summs <- entrez_summary(db="snp", id=r_search$id, version="2.0", retmode="xml")
    uid <- unique(extract_from_esummary(multi_summs, c("SNP_ID")))
    all_recs <- entrez_fetch(db="snp", id=uid, rettype="xml")
    tax_list <- XML::xmlToList(all_recs)
    
    if(tax_list$DocumentSummary$SNP_CLASS == "snv"){
      dbSNP_res <- format_dbSNP(tax_list$DocumentSummary$DOCSUM)
      ## multiple Alt Alleles are on different lines 
      final_snp <- cbind(rsid_input, dbSNP_res)
      colnames(final_snp) <- c("rsID", "CHR", "START", "REF.ALLELE", "ALT.ALLELE")
    ## if rsID doesn't represent SNV
    }else{
      final_snp <- cbind(rsid_input, NA,NA,NA,NA)
      colnames(final_snp) <- c("rsID", "CHR", "START", "REF.ALLELE", "ALT.ALLELE")
    }
  }
  return(final_snp)
}

#######################################################
## Function related to SNP strand and strandflipping ##
#######################################################

complement <- function(x) {
  switch (
    x,
    "A" = "T",
    "C" = "G",
    "T" = "A",
    "G" = "C",
    return(NA)
  )
}


###################################################
## Get coverage with blacklisted bad file region ##
###################################################

get_cov <- function(snp_info){
  snp <- GRanges(seqnames=as.numeric(snp_info$CHR), ranges=IRanges(start=as.numeric(snp_info$START), end=as.numeric(snp_info$START)+1),starts.in.df.are.0based=TRUE)
  hits <- findOverlaps(gr,snp)
  if(length(hits@from) > 0){
    bedcov = TRUE
  }else{
    bedcov = FALSE
  }
  return(bedcov)
}

###################
## Dlinkpipeline ##
###################

## function to check ldproxy res for coverage
ldproxy_check_cov <- function(ldpoxy_inter_df2){
  ## format chr to work with get_cov function
  ldpoxy_inter_df2$CHR <- gsub("chr", "", ldpoxy_inter_df2$CHR)
  ## loop to get first resul with coverage 
  res_ind <- vector()
  for(ldp in 1:nrow(ldpoxy_inter_df2)){
    gc <- get_cov(ldpoxy_inter_df2[ldp,])
    if(gc == FALSE){
      res_ind <- c(res_ind,ldp)
    }
  }
  ##if none of the returned rsIDs have coverage
  if(length(res_ind) == 0){
    ldpoxy_output <- LDproxy_NA_res(snp)
  }else{
    ldpoxy_output <- ldpoxy_inter_df2[res_ind,]
  }
  return(ldpoxy_output)
}

## return empty df (Quotes an issue?)
LDproxy_NA_res <- function(snp){
  ldpoxy_output <- cbind(snp, NA, NA, NA,NA,NA)
  colnames(ldpoxy_output) <- c("RSID_input", "RSID_Proxy", "CHR", "START", "REF.ALLELE", "ALT.ALLELE")
  return(ldpoxy_output)
}

## main function to get LDproxy results 
getLDproxy <- function(snp, pop, token, SNP_kept){
  ## run Query
  my_proxies <- LDproxy(snp, pop = pop, r2d = "r2", token = token, file = FALSE)
  ## error catching 
  if(grepl("error", my_proxies[1,1]) == TRUE){
    ldpoxy_output2 <- LDproxy_NA_res(snp)

  }else{
    ## extract only those with R2 >= 0.9 
    my_proxies_keep <- my_proxies[which(my_proxies$R2 >= 0.9),c(1,2,3)]
    ## remove those without rs number
    my_proxies_keep2 <- my_proxies_keep[grep("rs", my_proxies_keep$RS_Number),]
    ## remove those that are already in the dataset
    my_proxies_keep3 <- my_proxies_keep2[which(my_proxies_keep2$RS_Number %!in% SNP_kept),]
    if(nrow(my_proxies_keep3) == 0){
      ldpoxy_output2 <- LDproxy_NA_res(snp)
    }else{
      ## format coordinates
      crdf <- do.call(rbind,strsplit(my_proxies_keep3$Coord,":"))
      ## format alleles
      alsdf <- do.call(rbind,strsplit(gsub("\\(|\\)", "", my_proxies_keep3$Alleles), "/"))
      
      ## format results df
      ldpoxy_inter <- cbind(snp, my_proxies_keep3$RS_Number, crdf,alsdf)
      colnames(ldpoxy_inter) <- c("RSID_input", "RSID_Proxy", "CHR", "START", "REF", "ALT")
      ldpoxy_inter_df <- as.data.frame(ldpoxy_inter, stringAsFactors=F)
      
      ### remove those that aren't snps and aren't the 4 bases
      rm1 <- which(nchar(ldpoxy_inter_df$REF) != 1)
      rm2 <- which(nchar(ldpoxy_inter_df$ALT) != 1)
      rm3 <- which(ldpoxy_inter_df$REF %!in% c("A","T","G","C"))
      rm4 <- which(ldpoxy_inter_df$ALT %!in% c("A","T","G","C"))
      rm_all <- unique(c(rm1, rm2,rm3, rm4))
      
      ## if all removed
      if(length(rm_all) == nrow(ldpoxy_inter_df)){
        ldpoxy_output2 <- LDproxy_NA_res(snp)
        
      ## if some returned rsIDs are dodgy - remove
      }else if(length(rm_all) != 0){
        ldpoxy_inter_df2 <- ldpoxy_inter_df[-rm_all,]
        ## get ind of ldproxy without bedcov
        ldpoxy_output <- ldproxy_check_cov(ldpoxy_inter_df2)
        ## need to get info on ALT alleles via different function 
        lr <- lapply(ldpoxy_output[,"RSID_Proxy"], getrsID_info)
        ldpoxy_output2 <- cbind(snp, do.call(rbind, lr))
        colnames(ldpoxy_output2)[1:2] <- c("RSID_input", "RSID_Proxy")                        
        
      ## if all returned rsIDs are good    
      }else{
        ## check for bed coverage before output 
        ldpoxy_output <- ldproxy_check_cov(ldpoxy_inter_df)
        ## need to get info on ALT alleles via different function
        lr <- lapply(ldpoxy_output[,"RSID_Proxy"], getrsID_info)
        ldpoxy_output2 <- cbind(snp, do.call(rbind, lr))
        colnames(ldpoxy_output2)[1:2] <- c("RSID_input", "RSID_Proxy")       
      }
    }
  }
  ldpoxy_df <- as.data.frame(ldpoxy_output2, stringsAsFactors = FALSE)
  return(ldpoxy_df)
}  




