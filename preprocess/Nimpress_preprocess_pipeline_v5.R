####################################
##Automate data curation pipeline ##
####################################

## turn warnings off
oldw <- getOption("warn")
options(warn = -1)

## backend data files needed
assembly_file <- "../Suppl/GCF_000001405.13_GRCh37_assembly_report.txt"

pacman::p_load(docopt)

#############
## docopts ##
#############

'NIMPRESS preprocess
Usage:
  Nimpress_preprocess_pipeline_v3.R <input> [--bed=<bed> --outpath=<outpath> --out_prefix=<out_prefix> --title=<title> --description=<description> --citation=<citation> --genome_version=<genome version> --offset=<offset> ]  
  Nimpress_preprocess_pipeline_v3.R (-h | --help)
  Nimpress_preprocess_pipeline_v3.R --version
Options:
  -h --help                         Show this screen.
  --version                         Show version.
  --bed=<bed>                       Bed coverage file [default: NULL]
  --outpath=<outpath>               Path to output location [default: ./]
  --out_prefix=<out_prefix>         Prefix to outfile [default: Preprocessed]
  --title=<title>                   NIMPRESS file title [default: NIMPRESS preprocessed input]             
  --description=<description>       Description in NIMPRESS file [default: No descripton]
  --citation=<citation>             Data source [default: no citation]
  --genome_version=<genome version> Version of genome [default: GRCh37]
  --offset=<offset>                 Offset for NIMPRESS [default: 0.0]
Arguments:
    input     text file containing risk loci file in template format
    
' -> doc

arguments <- docopt(doc, version = 'NIMPRESS Preprocess v1')
#print(arguments)

dir.create(file.path(arguments$outpat, "Output"), showWarnings = FALSE)
dir.create(file.path(paste(arguments$outpat, "Output", sep="/"), "PLINK"), showWarnings = FALSE)
dir.create(file.path(paste(arguments$outpat, "Output", sep="/"), "NIMPRESS"), showWarnings = FALSE)


###############
## libraries ##
###############

pacman::p_load(rentrez,bedr,LDlinkR,stringr,GenomicRanges)

ela <- Sys.time() 

#arguments <- list()
#arguments$bed <- "/Users/ewilkie/Documents/Polygenic/DataCurationPipeline/mgrb_tier12.bed"
#arguments$input <- "./Data/Evans_Scott_2014_Table1.csv"

###################
## Initial setup ##
###################

## Get assembly info
## Instead downloaded file from: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13

ass <- read.table(assembly_file, header=F, sep="\t")
ass_sub <- ass[,c(3,7)]
assembly <- ass_sub[grep("NC_", ass_sub[,2]),]

## optional bed file
if(arguments$bed != "NULL"){
  ## use granges to check for bed in range
  bedfile <- arguments$bed
  ovlp <- read.table(bedfile, header = FALSE, stringsAsFactors = FALSE);
  colnames(ovlp) <- c("chr", "start", "end")
  gr <- makeGRangesFromDataFrame(ovlp, keep.extra.columns = TRUE)
}

###########################
## Format Loci risk file ##
###########################

test_file <- read.table(arguments$input, sep=",", header=T)

## remove blank rows
blank <- which(test_file[,1] == "")
na <- which(is.na(test_file[,1]))
rm_blank <- c(blank, na)
if(length(rm_blank) > 0){
 test_file <- test_file[-rm_blank,]
}

## rmeove blank columns
if(ncol(test_file) > 6){
  test_file <- test_file[,-(7:ncol(test_file))]
}

## remove all leading and trailing blank spaces
test_file <- as.data.frame(apply(test_file,2,function(x)gsub('\\s+', '',x)))

## Format Effect size

## if OR
or.ind <- grep("OR", colnames(test_file))
beta.ind <- grep("Beta", colnames(test_file))

if(length(or.ind) == 1){
  OR <- as.vector(test_file[,or.ind])
  test_file[,or.ind] <- log(as.numeric(sub('\\(.*', '', OR)))
  colnames(test_file)[or.ind] <- "Effect.size"
}else if(length(beta.ind) == 1){
  colnames(test_file)[beta.ind] <- "Effect.size"
}


## extract unique rsID
rsID_ind <- grep("rsID", colnames(test_file))
rsID <- test_file[,rsID_ind]
rsIDu <- as.vector(unique(rsID))

##########################################
## get rdID genomic location from dbSNP ##
##########################################

## this function lookes up information in the dbSNP database for each rsID
## extracts information for GRCh37: CHR, start, REF.Allele, ALT.Allele
## same results as a list it handle if multiple ALT alleles exist

getrsID_info <- function(rsid_input){
  
  ## return two lists
  final_snp <- list()
  
  snp_term <- paste(rsid_input, "[RS]", sep="")
  r_search <- entrez_search(db="snp", term=snp_term)
  multi_summs <- entrez_summary(db="snp", id=r_search$id)
  
  ## get the unique SNP_ID
  uid <- unique(extract_from_esummary(multi_summs, c("snp_id")))
  all_recs <- entrez_fetch(db="snp", id=uid, rettype="xml")
  tax_list <- XML::xmlToList(all_recs)
  
  ## extract assembly, genome position and variant details 
  g1 <- gsub("^[A-Za-z]*=", "", tax_list$DocumentSummary$DOCSUM)
  g2 <- strsplit(g1,"\\|")
  g3 <- as.data.frame(strsplit(g2[[1]][1],",")[[1]])
  colnames(g3) <- "ID"
  g4 <- apply(g3, 2, function(x) strsplit(x,":g\\."))
  g5 <- as.data.frame(do.call(rbind, g4$ID),stringsAsFactors=FALSE)
  colnames(g5) <- c("assembly", "START")
  g6 <- strsplit(g5$START,"[0-9]")
  g7 <- unlist(lapply(g6,function(x) x[length(x)]))
  g8 <- do.call(rbind, strsplit(g7, ">"))
  g9 <- cbind(g5,g8)
  rmg <- gsub("[^0-9]", "", g9$START)
  g9$START <- rmg
  if(!is.null(g9$g8)){
    g9 <- g9[-which(g9$g8 == "del"),]
  }
  if(nrow(g9) == 0 ){
    return(NA)   
  }else{
    g10 <- g9
    colnames(g10)[3:4] <- c("REF_Allele","ALT_Allele")
    
    ## subset to only NC
    g11 <- g10[grep("^NC",  g10[,1]),]
    
    ## extract CHR
    g12 <- gsub("NC_0+","",g11[,1])
    CHR <- gsub("\\.[0-9]*", "", g12)
    
    g13 <- data.frame(g11$assembly, CHR, g11$START, g11$REF_Allele, g11$ALT_Allele)
    
    ## Get right Assembly from SNP
    inter <- intersect(g13[,1], assembly[,2])
    g14 <- unique(g13[grep(inter,g13[,1]),])
    g15 <- g14[,-1]
    ## multiple Alt Alleles are on different lines 
    g16<- cbind(rsid_input, g15)
    colnames(g16) <- c("rsID", "CHR", "START", "REF.ALLELE", "ALT.ALLELE")
    
    final_snp$all.alleles <- g16
    
    g17 <-c(as.vector(unique(g16$CHR)), as.vector(unique(g16$START)), rsid_input)
    final_snp$rsID.loc <- g17
  
    return(final_snp)
  }
}


## list containing "rsID", "CHR", "START", "REF.ALLELE", "ALT.ALLELE" - for Risk allele check
rsID_genome <- list()
## list containing chr, start, rsid for coverage check
rsID_loc <- list()
for (rsid in 1: length(rsIDu)){
  res <- getrsID_info(rsIDu[rsid])
  if(!is.na(res)){
    rsID_loc[[rsid]] <- res$rsID.loc
    rsID_genome[[rsid]] <- res$all.allele
  }
}

rsID_loc_df <- do.call(rbind, rsID_loc)
colnames(rsID_loc_df) <- c("CHR","START","rsID")
urer <- rsID_loc_df

rsID_genome_df <- do.call(rbind, rsID_genome)


##################
## Get coverage ##
##################

get_cov <- function(urer,s){
  snp <- GRanges(seqnames=as.numeric(urer[s,1]), ranges=IRanges(start=as.numeric(urer[s,2]), end=as.numeric(urer[s,2])))
  hits <- findOverlaps(gr,snp)
  if(length(hits@from) > 0){
    cov = TRUE
  }else{
    cov = FALSE
  }
  return(cov)
}



if(arguments$bed == "NULL"){
  cov <- TRUE
  urercov <- cbind(urer, cov)
}else{
  urer <- rsID_loc_df
  cov <- vector()
  for (s in 1:nrow(urer)){
    nc <- get_cov(urer,s)
    cov <- c(cov, nc)
  
  }
  urercov <- cbind(urer, cov)
}

## this is for getLDproxy function
original_SNP <- urercov[,"rsID"]


###################
## Dlinkpipeline ##
###################

## coords are for GRCh37

## certain rsIDs can be perfectly linked. 
## in this case the LDproxy rsID turned out to be the same
##https://ldlink.nci.nih.gov/?var1=rs4948492&var2=rs4245597&pop=GBR&tab=ldpair
## could potentially also arise due to different issues, but the solution is to keep track of the new rsIDs and if that is already obtained, find another one, if no exist, drop the original rsID
## how to do this?

'%!in%' <- function(x,y)!('%in%'(x,y))

getLDproxy <- function(snp){
  ## for breaking
  new_rd <- TRUE
  ## run Query
  my_proxies <- LDproxy(snp, pop = "GBR", r2d = "r2", token = "cbe1b45bc8be", file = FALSE)
  ## extract only those with R2 >= 0.9
  my_proxies_keep <- my_proxies[which(my_proxies$R2 >= 0.9),c(1,2,3)]
  ## remove those without rs number
  my_proxies_keep2 <- my_proxies_keep[grep("rs", my_proxies_keep$RS_Number),]
  ## remove those that are already in the dataset
  my_proxies_keep3 <- my_proxies_keep2[-which(my_proxies_keep2$RS_Number %in% original_SNP),]
  ## if no proxies
  if(nrow(my_proxies_keep3) == 0){
    new_rd <- NA
  }else{
    ## split coords
    coord <- sub(".*:(\\d+).*", "\\1", my_proxies_keep3$Coord)
    coord_check <- paste(my_proxies_keep3$Coord , "-" , coord, sep="")
    for(i in 1:length(coord_check)){
      if(arguments$bed == "NULL"){
        getALT <- getrsID_info(my_proxies_keep3[i,1])
        if(!is.na(getALT)){
          ALTcol <- paste(sort(as.vector(getALT$all.alleles$ALT.ALLELE)), collapse=",")
          Alleles <- paste("(", unique(as.vector(getALT$all.alleles$REF.ALLELE)) , "/", ALTcol, ")", sep="")
          new_rd <- c(as.matrix(my_proxies_keep3[i,1:2]), Alleles)
          rsid_keep <- c(rsid_keep,as.vector(my_proxies_keep3[i,"RS_Number"]))
          break
        }else{
          next
        }
      }else{
        coord <- coord_check[i]
        snp_chr <- sub("chr","" , sub('\\:.*', '', coord))
        start <- as.numeric(sub("\\-.*", "", sub('.*\\:', '', coord)))
        snp <- GRanges(seqnames=snp_chr, ranges=IRanges(start=start, end=start))
        hits <- findOverlaps(gr,snp)
        if(length(hits@from) > 0){
          if(as.vector(my_proxies_keep3[i,"RS_Number"]) %!in% rsid_keep){
            ## need to get all alternative allelse fot this new rsID 
            getALT <- getrsID_info(as.vector(my_proxies_keep3[i,1]))
            if(!is.na(getALT)){
              ALTcol <- paste(sort(as.vector(getALT$all.alleles$ALT.ALLELE)), collapse=",")
              Alleles <- paste("(", unique(as.vector(getALT$all.alleles$REF.ALLELE)) , "/", ALTcol, ")", sep="")
              new_rd <- c(as.matrix(my_proxies_keep3[i,1:2]), Alleles)
              rsid_keep <- c(rsid_keep,as.vector(my_proxies_keep3[i,"RS_Number"]))
              break
            }else{
              next
            }
          }
        }else{
          new_rd <- NA
        } 
      }
    }
  }
  return(new_rd)
}  


rsid_keep <- vector()
LDres <- list()
for(n in 1:nrow(urercov)){
  if(urercov[n,4] == "FALSE"){
    LD <- getLDproxy(as.vector(urercov[n,3]))
    if (length(LD) == 3){
      coord <- strsplit(as.vector(LD[2]), ":")[[1]]
      j <- as.vector(LD[3])
      all <- regmatches(j, gregexpr("(?<=\\().*?(?=\\))", j, perl=T))[[1]]
      ref <- strsplit(all, "/")[[1]][1]
      alt <- strsplit(all, "/")[[1]][2]
      res <- c(as.vector(LD[1]), gsub("chr", "", coord[1]), coord[2], ref, alt)
      LDres[[n]] <- res
      rsid_keep <- c(rsid_keep, as.vector(LD[1]))
    }else if(is.na(LD)) {
      LDres[[n]] <- NA
    }
  }else{
    LDres[[n]] <- NA
  }
}


LDres_mat <- do.call(rbind, LDres)
if(ncol(LDres_mat) == 5){
  colnames(LDres_mat) <- c("ALT.rsID", "ALT.chr", "ALT.start","REF.new", "ALT.new")
}else{
  ## when there are no LDproxy results, either all cov == TRUE or when cov == FALSE, but no proxy exists 
  ## this is so code below doesn't break
  LDres_mat <- cbind(LDres_mat, LDres_mat,LDres_mat,LDres_mat,LDres_mat)
  colnames(LDres_mat) <- c("ALT.rsID", "ALT.chr", "ALT.start","REF.new", "ALT.new")
}

## combine original rsID and alternative
comb <- as.data.frame(cbind(urercov,LDres_mat))

#################################
#################################
## Combine results into output ##
#################################
#################################

##########################################################
## Check for strand FLIPPING AND DEFINE CORRECT ALLELES ##
##########################################################

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

check_multi_sub <- split(test_file, f = test_file$Subtype)

subs <- list()
for(type in 1:length(check_multi_sub)){
  out <- list()
  oc <- 0
  ## get risk allele
  test_f <- check_multi_sub[[type]]
  for(var in 1:nrow(comb)){
    oc <- oc + 1
    ## get rsID db alelles
    ref_dat <- merge(comb[var,], rsID_genome_df, by="rsID")
    ## match between original and final 
    relv <- test_f[which(as.vector(test_f$rsID) == as.vector(comb[var,"rsID"])),]
    ## if this rsID is not in this subtyope - skip
    if(nrow(relv) == 0){
      out[[oc]] <- NA
    }else{
      # if the SNP has coverage, LDproxy sub not needed
      if(unique(ref_dat$cov == TRUE)){
        ## determine which allele the risk alelle correpsonds to: ALT, REF or complement 
        ## need to incorporate multiple ALT.Alleles
        taa <- which(relv[,"Risk_allele"] %in% as.vector(ref_dat$ALT.ALLELE))
        tar <- which(relv[,"Risk_allele"] %in% as.vector(ref_dat$REF.ALLELE))
        
        ## alt_com
        taac <- which(relv[,"Risk_allele"] %in% sapply(as.vector(ref_dat$ALT.ALLELE), complement))
        tarc <- which(relv[,"Risk_allele"] %in% sapply(as.vector(ref_dat$REF.ALLELE), complement))
        
        ## if risk allele is ALT
        if(length(taa) == 1){
          o <- c(as.matrix(ref_dat[taa,c(1:3,12:13,13)]), as.matrix(relv[,3:ncol(relv)]))
          out[[oc]] <- o
          ## if risk allele is REF  
        }else if(length(tar) == 1){
          o <- c(as.matrix(ref_dat[tar,c(1:3,12:13,12)]), as.matrix(relv[,3:ncol(relv)]))
          out[[oc]] <- o
          ## if risk allele is complement of ALT  - strand flipping
        }else if(length(taac) == 1){
          ## change direction of odds ratio
          new.Effect.size <- as.numeric(as.vector(relv[,"Effect.size"])) * -1
          o <- c(as.matrix(ref_dat[taac,c(1:3,12:13,13)]), as.matrix(relv[,3]), new.Effect.size, as.matrix(relv[,c(5:ncol(relv))]))
          out[[oc]] <- o
          ## if not AT.ALLELE or REF.ALLELE -> strand filling has occured  
          ## check for complement of REF.ALLELE    
        }else if (length(tarc) == 1){
          ## change direction of odds ratio
          new.Effect.size <- as.numeric(as.vector(relv[,"Effect.size"])) * -1
          o <- c(as.matrix(ref_dat[tarc,c(1:3,12:13,12)]), as.matrix(relv[,3]), new.Effect.size, as.matrix(relv[,c(5:ncol(relv))]))
          out[[oc]] <- o
        }
        ## if multiple alt alleles change to include
        if(nrow(ref_dat) >1){
          alt <- paste(sort(as.vector(ref_dat$ALT.ALLELE)), collapse=",")
          out[[oc]][5] <- alt 
        }
        ## if the SNP doesn't have coverage, Need to do LDpoxy SUB
        }else if(unique(ref_dat$cov == FALSE)){
        ## ignore results which don't have an alternative
          if(!is.na(unique(as.vector(ref_dat$ALT.rsID)))){
            
            waa <- which(as.vector(relv[,"Risk_allele"]) %in% as.vector(ref_dat$ALT.ALLELE))
            war <- which(as.vector(relv[,"Risk_allele"]) %in% as.vector(ref_dat$REF.ALLELE))
            
            ## alt_com
            waac <- which(as.vector(relv[,"Risk_allele"]) %in% sapply(as.vector(ref_dat$ALT.ALLELE), complement))
            warc <- which(as.vector(relv[,"Risk_allele"]) %in% sapply(as.vector(ref_dat$REF.ALLELE), complement))
            
            ## if ALT.ALLELE
            if(length(waa) == 1){
              ## change to ALT.new 
              o <- c(as.matrix(ref_dat[waa,c(5:9,9)]), as.matrix(relv[,3:ncol(relv)]))
              out[[oc]] <- o
              ## if not ALT. ALLELE, then check for REF.Allele 
            }else if(length(war) == 1){
              ## change to REF.new 
              o <- c(as.matrix(ref_dat[war,c(5:9,8)]), as.matrix(relv[,3:ncol(relv)]))
              out[[oc]] <- o
              ## if not AT. ALLELE or REF.ALLELe -> strand filling has occured  
              ## check for complement of ALT.ALLELE  
            }else if (length(waac) == 1){
              ## change direction of odds ratio
              new.Effect.size <- as.numeric(as.vector(relv[,"Effect.size"])) * -1
              o <- c(as.matrix(ref_dat[waac,c(5:9,9)]), as.matrix(relv[,3]), new.Effect.size, as.matrix(relv[,c(5:ncol(relv))]))
              out[[oc]] <- o
              ## if not AT.ALLELE or REF.ALLELE -> strand flipping has occured  
              ## check for complement of REF.ALLELE    
            }else if (length(warc) == 1){
              ## change direction of odds ratio
              new.Effect.size <- as.numeric(as.vector(relv[,"Effect.size"])) * -1
              o <- c(as.matrix(ref_dat[warc,c(5:9,8)]), as.matrix(relv[,3]), new.Effect.size, as.matrix(relv[,c(5:ncol(relv))]))
              out[[oc]] <- o
            }
          ## ignore results which don't have an alternative
          }else{
          out[[oc]] <- NA
        }
      }## if multiple alt alleles change to include both as comma seperated
      if(nrow(ref_dat) > 1){
        alt <- paste(sort(as.vector(ref_dat$ALT.ALLELE)), collapse=",")
        out[[oc]][5] <- alt 
      }
    }
  }
  subs[[type]] <- do.call(rbind,out)
}

otdf <- do.call(rbind,subs)

na.ind <- which(is.na(otdf[,1]))
if (length(na.ind) == 0){
  colnames(otdf) <- c("rsID","CHR","START","REF.ALLELE","ALT.ALLELE","Risk_allele","Freq","Effect.size","P","Subtype")
  final <- as.data.frame(otdf,stringsAsFactors=FALSE)
}else{
  otdf <- otdf[-na.ind,]
  colnames(otdf) <- c("rsID","CHR","START","REF.ALLELE","ALT.ALLELE","Risk_allele","Freq","Effect.size","P","Subtype")
  final <- as.data.frame(otdf,stringsAsFactors=FALSE)
}

#######################################################################################
## If rsIDs after linkage turn out to be the same                                    ##
## flag this dataset and don't return results since currently don't know what to do  ##
#######################################################################################

sspl <- split(final, f = final$Subtype)

for(sub in 1:length(sspl)){

  allrsID <- as.vector(sspl[[sub]]$rsID)
  ursID <- unique(allrsID)
  
  if(length(allrsID) != length(ursID)){
    stop("Duplicated rsIDs in output")
  }
}


##################
## WRITE OUTPUT ##
##################

out <- split(final, f = final$Subtype)

for(type in 1:length(out)){

  #####################
  ## Nimpress output ##
  #####################
  ## file
  u <- gsub("\\s+", "_", unique(out[[type]]$Subtype))
  filen <- paste(arguments$outpath,"/Output/NIMPRESS/", arguments$out_prefix, "_", u, "_NIMPRESS_input.txt", sep="")
  ## if file exisits, delete content
  ## title
  write(arguments$title,file=filen, append=FALSE)
  ## description
  write(arguments$description,file=filen, append=TRUE)
  ## citation
  write(arguments$citation,file=filen, append=TRUE)
  ## genome version
  write(arguments$genome_version,file=filen, append=TRUE)
  ## ofset
  write(arguments$offset,file=filen, append=TRUE)
  
  ##data
  x <- out[[type]][,c("CHR","START", "REF.ALLELE", "Risk_allele", "Effect.size", "Freq")]
  
  ## recode NA to NaN 
  x$Freq <- "NaN"
  
  write.table(x, file=filen, sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)

  ###############################
  ## Scores for use with PLINK ##
  ###############################

  p <- gsub("\\s+", "_", unique(out[[type]]$Subtype))

  ## for locationID
  filep <- paste(arguments$outpath,"/Output/PLINK/", arguments$out_prefix, "_", p, "_PLINK_input.txt", sep="")
  
  out[[type]][which(out[[type]][,"REF.ALLELE"] == "-"),"REF.ALLELE"] <- "*"
  out[[type]][which(out[[type]][,"Risk_allele"] == "-"),"Risk_allele"] <- "*"
  
  z.rsID <- out[[type]][,c("rsID", "Risk_allele", "Effect.size", "P" )]
  write.table(z.rsID, file=filep, sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
}

ela <- Sys.time() - ela
print(ela)
