##########################################
## GWAS summary data curations pipeline ##
##########################################

## turn warnings off
# oldw <- getOption("warn")
# options(warn = -1)

suppressPackageStartupMessages(library(docopt))

## load processing functions files
source("functions.R")

#############
## docopts ##
#############

'makescore.R: Generate an input score file for nimpress from GWAS summary statistics.

Usage:
  makescore.R [options] <infile> <outdir>
  makescore.R (-h | --help)
  makescore.R --version

Arguments:
  <infile>                    Path to the file containing risk loci file in template format (see below).
  <outdir>                    Path to output folder. This will be created if it doesn\'t exist.

General options:
  -h --help                   Show this screen.
  --version                   Show version.

Output options:
  --description="<string>"    Description of the polygenic score, will be copied into the output file. [default: ""]
  --citation="<string>"       Citation of the polygenic score, will be copied into the output file. [default: ""]
  --offset=<offset>           Offset of the polygenic score. This value will be added to calculated scores (score = sum(dosage*weight) + offset). [default: 0.0]

Blacklist options:
  --blacklist=<file>          BED file of blacklist genome regions. Loci falling in these regions will either be imputed with proxy SNPs that are not in blacklist regions (if LDproxy_token and LDproxy_pop are supplied), or dropped (otherwise).

Proxy SNP imputation options:
  --LDproxy_token=<string>    LDproxy API token. Generate at https://ldlink.nci.nih.gov/?tab=apiaccess
  --LDproxy_pop=<string>      Background populations for LDproxy (eg. GBR)


Input file format
=================
  Tab-delimited, with a header and columns:
    rsID                      dbSNP SNP cluster ID for this score locus.
    effectAllele              The allele in this SNP associated with a differential trait phenotype or risk.
    effectAlleleFrequency     The population frequency (ideally, in the population of interest) of the effect allele.
    weight                    The contribution of a single copy of EffectAllele to the trait phenotype or risk.

  The final score is calculated as score = sum(dosage*weight) + offset, where dosage is the number of copies
  of the EffectAllele in each individual (dosage \\in {0, 1, 2}). Typically, weight is the reported Beta values
  (for continuous traits), or log odds (for binary traits / risk), of a GWAS.

  An example score file based on Wood et al (2014) can be found in example/wood_height_input.tsv.


Examples
========
1. Basic use, converts GWAS summary statistics to a polygenic score without any imputation:
   
   Rscript makescore.R example/wood_height_input_small.tsv example/wood_height_small_example1

2. Convert GWAS summary statistics, removing loci falling within blacklist regions:

   Rscript makescore.R --blacklist=example/wood_height_blacklist.bed example/wood_height_input_small.tsv example/wood_height_small_example2

3. Convert GWAS summary statistics, imputing loci falling within blacklist regions to linked loci outside of blacklist regions:

   Rscript makescore.R --blacklist=example/wood_height_blacklist.bed --LDproxy_token=<token> --LDproxy_pop=GBR example/wood_height_input_small.tsv example/wood_height_small_example3

' -> doc

arguments <- docopt(doc, version='nimpress makescore v1.0\n', strict=FALSE)

###############################################################################################
###############################################################################################
#####################                     Initial setup              ##########################
###############################################################################################
###############################################################################################

message("[1/12] Loading libraries")

suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(LDlinkR))
suppressPackageStartupMessages(library(rentrez))

######################
## Check input file ##
######################

message(paste("[2/12] Checking file: ", arguments$infile, sep="" ))

## check and format input file
gf_ok <- check_gwas_file(arguments$infile)

##################
## setup outdir ##
##################
message("[3/12] Setting up output dir...")

## get file name from input file
output_name1 <- gsub(".*/","", arguments$infile)
output_name2 <- gsub("\\.tsv","", output_name1)
if(is.null(arguments$outdir)){
  outdir <- "Nimpress_preprocess_Output"
  dir.create(file.path(getwd(),outdir), showWarnings = FALSE)
}else{
  outdir <- arguments$outdir
  dir.create(file.path(outdir), showWarnings = FALSE)
}

#######################################
## setup blacklisted genome bed file ##
#######################################

message("[4/12] Setting up blacklisted genome file...")
## this is for inbuilt blacklist file
if(!is.null(arguments$blacklist)){
  bedfile <- arguments$blacklist
  cust.bed <- fread(bedfile,sep="\t", header=TRUE,stringsAsFactors=FALSE)
  gr <- bedfile_to_Granges(cust.bed)
}else{
  bedfile = NULL
}

#############
## LDproxy ##
#############

message("[5/12] Testing LDproxy parameters... ")
## check if LDproxy is on
LDproxy_flag = "OFF"
if(!is.null(arguments$LDproxy_pop) & !is.null(arguments$LDproxy_token)){
  ## to check LDproxy pop input - get all available pops from package
  pop <- list_pop()
  bgpc <- pop$pop_code
  ## check valid background population
  if(arguments$LDproxy_pop %!in% bgpc){
    stop(paste(arguments$LDproxy_pop, " is not a valid background population. Select one from: ", paste(bgpc, collapse=",") , sep="" ))
  }
  ## test query
  testproxy <- LDproxy("rs456", arguments$LDproxy_pop, "r2", token = arguments$LDproxy_token)
  if(testproxy[1,1] != "  error: Invalid or expired API token. Please register using the API Access tab: https://ldlink.nci.nih.gov/?tab=apiaccess,"){
    message("LDproxy parameters ok")
    LDproxy_flag = "ON"
  }else{
    stop_quietly()
  }
}

if(LDproxy_flag == "ON" & is.null(bedfile)){
  stop("Parameters indicate LDproxy is desired, but no bedfile provided - please check you input parameters")
}

###################
## Genomic files ##
###################
 
message("[6/12] Setting up genomic environment...")

gv <- "GRCh37"
assembly_file <- "Suppl/GCF_000001405.13_GRCh37_assembly_report.txt"

ass <- read.table(assembly_file, header=F, sep="\t",stringsAsFactors = F)
ass_sub <- ass[,c(3,7)]
assembly <- ass_sub[grep("NC_", ass_sub[,2]),]
colnames(assembly) <- c("CHR", "NC_CHR")


#####################
#####################
## File processing ##
#####################
#####################

message(paste("[7/12] Processing file: ", arguments$infile," ...this may take a while...", sep="" ))

##########################################
## get rdID genomic location from dbSNP ##
##########################################

rsIDu <- as.vector(gf_ok$rsID)

rsID_loc <- list()
for (rsid in 1: length(rsIDu)){
  message(paste("Getting info for: ", rsIDu[rsid],sep=""))
  res <- getrsID_info(rsIDu[rsid])
  rsID_loc[[rsid]] <- res
}
rsID_loc_df <- do.call(rbind, rsID_loc)

## remove NA results returned
rs.info.na <- which(is.na(rsID_loc_df$CHR))
if(length(rs.info.na) != 0 ){
  rs.rm <- rsID_loc_df[rs.info.na,] 
  rs.rm.df <- cbind(rs.rm,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
  colnames(rs.rm.df) <- c("rsID", "CHR.x","START.x","REF.ALLELE.x","ALT.ALLELE.x","SEQ","strand","ambiguous","effectAllele","risk_type","flipped","bedcov","RSID_Proxy","CHR.y","START.y","REF.ALLELE.y","ALT.ALLELE.y")
  rsID.df <- rsID_loc_df[-rs.info.na,]
}else{
  rs.rm.df <- NULL
  rsID.df <- rsID_loc_df
}

##################################
## Check REF allele with genome ##
##################################

## create range of snps
snp_gr <- GRanges(seqnames=paste("chr", as.numeric(rsID.df$CHR), sep=""), ranges=IRanges(start=as.numeric(rsID.df$START)), starts.in.df.are.0based=F)

## results returned are in forward strand
seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, snp_gr)
seqs_df <- as.data.frame(seqs)
colnames(seqs_df) <- "SEQ"

SNP_info <- cbind(rsID.df, seqs_df)

## check for strand and ambigiouty between REF and ALT
strand <- list()
ambig <- list()

stp <- strsplit(SNP_info$ALT.ALLELE, "|")
for(rs in 1:nrow(SNP_info)){
  ## check whether REF is forwards or reverse strand (same as SEQ)
  if(SNP_info[rs,"REF.ALLELE"] == SNP_info[rs,"SEQ"]){
    strand[[rs]] <- "+"
  } else if(SNP_info[rs,"REF.ALLELE"] == complement(SNP_info[rs,"SEQ"])){
    strand[[rs]] <- "-"
  } else{
    strand[[rs]] <- NA
  }
  
  if(length(stp[[rs]]) == 1 && complement(SNP_info[rs,"REF.ALLELE"]) == stp[[rs]]){
    ambig[[rs]] <- "Y"
  }else if (length(stp[[rs]]) > 1 && complement(SNP_info[rs,"REF.ALLELE"]) %in% stp[[rs]]){
    ambig[[rs]] <- "D"
  }else{
    ambig[[rs]] <- "N"
  }
}

SNP_is <- cbind(SNP_info, strand=unlist(strand), ambiguous= unlist(ambig))

##############################
## Check for Ambiguous SNPS ##
##############################

message("[8/12] Defining ambigious SNPs...")

## Check whether Risk ALlele is REF or ALT or neither (wrong)/ Flipped 
rsm <- merge(SNP_is, gf_ok[,c("rsID","effectAllele")], by="rsID")
stp2 <- strsplit(rsm$ALT.ALLELE, "|")

risk_type <- list()
flipped <- list()

for(rll in 1:nrow(rsm)){
  ## if REF and ALT are ambiguous, then the risk_type and flipped are NA 
  if(rsm[rll, "ambiguous"] == "Y"){
    risk_type[[rll]] = NA
    flipped[[rll]] = NA
  }
  ## if ambigious flag is D or N - first check for ambigious then for check type and flipped
  else{
    if(rsm[rll,"effectAllele"] == complement(rsm[rll,"REF.ALLELE"]) & rsm[rll,"effectAllele"] %in% stp2[[rll]]){
      risk_type[[rll]] = NA
      flipped[[rll]] = NA
      rsm[rll,"ambiguous"] = "Y"
    }
    # Risk == REF and  == complemnt(ALT)
    # same as complment(risk) == complment(REF) & ALT
    else if(complement(rsm[rll,"effectAllele"]) == complement(rsm[rll,"REF.ALLELE"]) & complement(rsm[rll,"effectAllele"]) %in% stp2[[rll]]){
      risk_type[[rll]] = NA
      flipped[[rll]] = NA
      rsm[rll,"ambiguous"] = "Y"
    }
    else if(rsm[rll,"effectAllele"] == rsm[rll,"REF.ALLELE"]){
      risk_type[[rll]] = "REF"
      flipped[[rll]] = "N"
      rsm[rll,"ambiguous"] = "N"
    }else if (rsm[rll,"effectAllele"] %in% stp2[[rll]]){
      risk_type[[rll]] = "ALT"
      flipped[[rll]] = "N"
      rsm[rll,"ambiguous"] = "N"
      ## complement means flipped
    }else if(complement(rsm[rll,"effectAllele"]) == rsm[rll,"REF.ALLELE"]){
      risk_type[[rll]] = "REF"
      flipped[[rll]] = "Y"
      rsm[rll,"ambiguous"] = "N"
    }else if (complement(rsm[rll,"effectAllele"]) %in% stp2[[rll]]){
      risk_type[[rll]] = "ALT"
      flipped[[rll]] = "Y"
      rsm[rll,"ambiguous"] = "N"
    }
  }
}

SNP_fin <- cbind(rsm, risk_type=unlist(risk_type), flipped=unlist(flipped))


###################################################
## Get coverage of rsID with blacklisted regions ##
###################################################

message("[9/12] Checking input SNPS for coverage...")

if(is.null(bedfile)){
  bedcov <- FALSE
  urercov <- cbind(SNP_fin, bedcov)
}else{
  bedcov <- vector()
  for (s in 1:nrow(SNP_fin)){
    nc <- get_cov(SNP_fin[s,])
    bedcov <- c(bedcov, nc)
  }
  urercov <- cbind(SNP_fin, bedcov)
}

#############
## LDproxy ##
#############

message("[10/12] Coverage and LDproxy... this will take a while... ")

## if no bedfile all bedcov will equal FALSE, if no fall in region no sub needed
if(length(unique(urercov$bedcov)) == 1 && unique(urercov$bedcov) == FALSE){
  ## append empty columns for downstream flag adding 
  LDproxy_out <- cbind(urercov,NA,NA,NA,NA,NA)
  colnames(LDproxy_out) <- c("rsID", "CHR.x", "START.x", "REF.ALLELE.x", "ALT.ALLELE.x", "SEQ", "strand", "ambiguous", "effectAllele", "risk_type", "flipped", "bedcov", "RSID_Proxy", "CHR.y", "START.y", "REF.ALLELE.y", "ALT.ALLELE.y")
  
## if LDproxy not used and bedfile provided, remove those without coverage. 
}else if (LDproxy_flag == "OFF" && !is.null(bedfile)){
  message("Remove blacklisted regions enabled, but LDproxy disabled...")
  message("rsID without coverage will be removed")

  ## append empty columns for downstream flag adding 
  LDproxy_out <- cbind(urercov,NA,NA,NA,NA,NA)
  colnames(LDproxy_out) <- c("rsID", "CHR.x", "START.x", "REF.ALLELE.x", "ALT.ALLELE.x", "SEQ", "strand", "ambiguous", "effectAllele", "risk_type", "flipped", "bedcov", "RSID_Proxy", "CHR.y", "START.y", "REF.ALLELE.y", "ALT.ALLELE.y")
  
## do ldproxy and exlude those that fall in the bed regions  - do sub and resub
}else if (LDproxy_flag == "ON" && !is.null(bedfile)){
  
  ## only run LDproxy on unique results for those that don't have coverage or are ambigious
  run_ind <- which(urercov$bedcov == TRUE & urercov$ambiguous != "Y")
  ldproxy_input <- unique(urercov[run_ind,"rsID"])
  ## to keep track of which RSids are already in DF, without expanding loop 
  SNP_kept <- unique(urercov[-run_ind,"rsID"])
  
  if(length(ldproxy_input) > 0){
    ldproxy_ls <- list()
    for(s in 1:length(ldproxy_input)){
      ## need this or it will error
      snp <- ldproxy_input[s]
      # main function to get ldproxy results
      ldproxy_res <- getLDproxy(snp, arguments$LDproxy_pop, arguments$LDproxy_token, SNP_kept)
      ## if results is NA add it anyways for later removal 
      if(is.na(ldproxy_res$RSID_Proxy)[1]){
        ldproxy_res_keep <- ldproxy_res
      }else{
        ## check which not in dataset 
        wni <- which(ldproxy_res$RSID_Proxy %!in% SNP_kept)
        if(length(wni) == 0){
          ldproxy_res_keep <- LDproxy_NA_res(ldproxy_input[s])
        }else{
          ## keep one result only that is not already in dataset
          ldproxy_res_keep <- ldproxy_res[wni[1],]
          ## this is so that no duplicates appear in the data
          SNP_kept <- c(SNP_kept,ldproxy_res_keep$RSID_Proxy)
        }
      }
      
      ## combine ldproxy_res with orignal input - use all to get also those that return NA for ldproxy to later filter
      mres <- merge(urercov, ldproxy_res_keep, by.x="rsID", by.y="RSID_input")
      ldproxy_ls[[s]] <- mres
    }
    
    LDproxy_in <- do.call(rbind, ldproxy_ls)
    
    ## if not all results don't have coverage
    if(nrow(urercov) !=  length(run_ind)){
      nocov_padd <- cbind(urercov[-run_ind,],NA,NA,NA,NA,NA)
      colnames(nocov_padd) <- c("rsID", "CHR.x", "START.x", "REF.ALLELE.x", "ALT.ALLELE.x", "SEQ", "strand", "ambiguous", "effectAllele", "risk_type", "flipped", "bedcov", "RSID_Proxy", "CHR.y", "START.y", "REF.ALLELE.y", "ALT.ALLELE.y")
    
      # final output
      LDproxy_out <- rbind(nocov_padd, LDproxy_in)
    }else{
      LDproxy_out <- LDproxy_in
    }
  }else{
    nocov_padd <- cbind(urercov,NA,NA,NA,NA,NA)
    colnames(nocov_padd) <- c("rsID", "CHR.x", "START.x", "REF.ALLELE.x", "ALT.ALLELE.x", "SEQ", "strand", "ambiguous", "effectAllele", "risk_type", "flipped", "bedcov", "RSID_Proxy", "CHR.y", "START.y", "REF.ALLELE.y", "ALT.ALLELE.y")
    
    ## final output
    LDproxy_out <- nocov_padd
  }
}

###########################
## Add initially removed ##
###########################

all.res <- rbind(LDproxy_out, rs.rm.df)

#####################
## Set filter flag ##
#####################

message("[11/12] Generate intermediate file...")

## ldproxy flag. NA if bedcov = FALSE, Y if becov = TRUE & !is.na(RSID_Proxy), N if becov = TRUE & is.na(RSID_Proxy)
ldproxy_flag <- rep(NA, nrow(all.res))
ldproxy_flag[(all.res$bedcov == TRUE & !is.na(all.res$RSID_Proxy))] <- "Y"
ldproxy_flag[is.na(all.res$bedcov)] <- "N"
ldproxy_flag[(all.res$bedcov == TRUE & is.na(all.res$RSID_Proxy))] <- "N"
ldproxy_flag[(all.res$bedcov == FALSE & is.na(all.res$RSID_Proxy))] <- "N"
LDproxy_outl <- cbind(all.res, ldproxy_flag)

rm_flag <- rep(NA, nrow(LDproxy_outl))
rm_flag[LDproxy_outl$ambiguous == "Y"] <- "Y"
rm_flag[(LDproxy_outl$ambiguous == "N" & LDproxy_outl$bedcov == TRUE & LDproxy_outl$ldproxy_flag == "N")] <- "Y"
rm_flag[(LDproxy_outl$ambiguous == "N" & LDproxy_outl$bedcov == TRUE & LDproxy_outl$ldproxy_flag == "Y")] <- "N"
rm_flag[(LDproxy_outl$ambiguous == "N" & LDproxy_outl$bedcov == FALSE & LDproxy_outl$ldproxy_flag == "N")] <- "N"
rm_flag[is.na(LDproxy_outl$bedcov)] <- "Y"
LDproxy_outlf <- cbind(LDproxy_outl, rm_flag)

LDproxy_outlf$rm_flag[is.na(LDproxy_outlf$rm_flag)] <- "Y"

## merge with previous results
LDproxyf <- merge(LDproxy_outlf, gf_ok[,c("rsID","effectAlleleFrequency","weight")], by="rsID")

## change column names
## write this to outfile - as intermediate file
colnames(LDproxyf) <- c("INPUT.rsID", "dbSNP.CHR", "dbSNP.START", "dbSNP.REF.ALLELE", "dbSNP.ALT.ALLELE", paste(gv,".ALLELE", sep=""), "STRAND", "FLAG.AMBIGUOUS", "INPUT.RISK.ALLELE", "INPUT.RISK.TYPE", "FLAG.RISK.FLIPPED", "BED.COVERAGE", "LDPROXY.rsID", "LDPROXY.CHR", "LDPROXY.START", "LDPROXY.REF.ALLELE","LDPROXY.ALT.ALLELE", "FLAG.LDPROXY", "FLAG.RM", "INPUT.ALLELE.FREQ", "INPUT.BETA")

## order columns
interm <- LDproxyf[,c("INPUT.rsID","dbSNP.CHR","dbSNP.START","dbSNP.REF.ALLELE","dbSNP.ALT.ALLELE","GRCh37.ALLELE","STRAND","INPUT.RISK.ALLELE","FLAG.RM","FLAG.AMBIGUOUS", "BED.COVERAGE", "INPUT.RISK.TYPE","FLAG.RISK.FLIPPED","FLAG.LDPROXY","LDPROXY.rsID","LDPROXY.CHR", "LDPROXY.START","LDPROXY.REF.ALLELE","LDPROXY.ALT.ALLELE","INPUT.ALLELE.FREQ","INPUT.BETA")]

## order rows
intermo <- interm[order(interm$FLAG.RM, interm$FLAG.AMBIGUOUS, interm$BED.COVERAGE, interm$FLAG.LDPROXY, interm$INPUT.RISK.TYPE, interm$FLAG.RISK.FLIPPED),]

## write intermediate results to outfile
rm_output_file <- paste(outdir,"/", output_name2, "_Intermediate_results.csv", sep="")
write.table(intermo,rm_output_file, sep=",", row.names=F, quote=F)

#############################################
#############################################
## Kept risk loci - DEFINE CORRECT ALLELES ##
#############################################
#############################################

message("[12/12] Generating NIMPRESS input file...")

## seperte kept from removed
interm_keep <-  interm[interm$FLAG.RM == "N",]

if(nrow(interm_keep) == 0){
  stop("No rsIDs remaining after filtering. Please check intermediate file to determine cause")
}

##########################################################
## FOR BEDCOV == FALSE and not removed due to ambiguity ##
##########################################################

## dbSNP change
dbSNP_ale <- interm_keep[is.na(interm_keep$FLAG.LDPROXY) | interm_keep$FLAG.LDPROXY == "N",]
dnSNP_final.ls <- list()

if(nrow(dbSNP_ale) != 0){
  for(dbal in 1:nrow(dbSNP_ale)){
    if(dbSNP_ale[dbal,"FLAG.RISK.FLIPPED"] == "N"){
      INPUT.RISK.ALLELE <- dbSNP_ale[dbal,"INPUT.RISK.ALLELE"]
      dnSNP_final.ls[[dbal]] <- data.frame(dbSNP_ale[dbal,c("dbSNP.CHR", "dbSNP.START", "dbSNP.REF.ALLELE")], new.INPUT.RISK.ALLELE = INPUT.RISK.ALLELE, dbSNP_ale[dbal,c("INPUT.RISK.ALLELE", "FLAG.RISK.FLIPPED","INPUT.BETA", "INPUT.ALLELE.FREQ","FLAG.LDPROXY")])
    }else if(dbSNP_ale[dbal,"FLAG.RISK.FLIPPED"] == "Y"){
      ## flip risk allele
      INPUT.RISK.ALLELE <- complement(dbSNP_ale[dbal,"INPUT.RISK.ALLELE"])
      dnSNP_final.ls[[dbal]] <- data.frame(dbSNP_ale[dbal,c("dbSNP.CHR", "dbSNP.START", "dbSNP.REF.ALLELE")], new.INPUT.RISK.ALLELE = INPUT.RISK.ALLELE, dbSNP_ale[dbal,c("INPUT.RISK.ALLELE","FLAG.RISK.FLIPPED", "INPUT.BETA", "INPUT.ALLELE.FREQ","FLAG.LDPROXY")])
    }
  }
}

dnSNP_final.df <- do.call(rbind,dnSNP_final.ls)

#####################################
## FOR BEDCOV == TRUE with LDPROXY ##
#####################################

ldproxy_ale <- interm_keep[(interm_keep$FLAG.LDPROXY == "Y" & !is.na(interm_keep$FLAG.LDPROXY)),]
LDPROXY_final.ls <- list()

if( nrow(ldproxy_ale) != 0 ){
  for(alle in 1:nrow(ldproxy_ale)){
    ## if not flipped REF or ALT doesn't matter
   if(ldproxy_ale[alle,"INPUT.RISK.TYPE"] == "REF"){
     LDPROXY_final.ls[[alle]] <- ldproxy_ale[alle,c("LDPROXY.CHR", "LDPROXY.START", "LDPROXY.REF.ALLELE", "LDPROXY.REF.ALLELE","INPUT.RISK.TYPE","INPUT.BETA", "INPUT.ALLELE.FREQ", "FLAG.LDPROXY")]
   } else if(ldproxy_ale[alle,"INPUT.RISK.TYPE"] == "ALT"){
     LDPROXY_final.ls[[alle]] <- ldproxy_ale[alle,c("LDPROXY.CHR", "LDPROXY.START", "LDPROXY.REF.ALLELE", "LDPROXY.ALT.ALLELE","INPUT.RISK.TYPE","INPUT.BETA", "INPUT.ALLELE.FREQ", "FLAG.LDPROXY")]
   }
    colnames(LDPROXY_final.ls[[alle]]) <- c("LDPROXY.CHR","LDPROXY.START","LDPROXY.REF.ALLELE","new.RISK.ALLELE","INPUT.RISK.TYPE","INPUT.BETA","INPUT.ALLELE.FREQ","FLAG.LDPROXY")
  }
}

LDPROXY_final.df <- do.call(rbind,LDPROXY_final.ls)


#####################
## Combine results ##
#####################


if(is.null(dnSNP_final.df)){
  dnSNPfdf <- NULL
}else if(nrow(dnSNP_final.df) > 0){
  dnSNPfdf <- dnSNP_final.df[,c(1:4,7:8)]
  dnSNPfdf <- rename_cols(dnSNPfdf)
}

if(is.null(LDPROXY_final.df)){
  LDPROXYfdf <- NULL
}else if(nrow(LDPROXY_final.df) > 0){
  LDPROXYfdf <- LDPROXY_final.df[,c(1:4,6:7)]
  LDPROXYfdf <- rename_cols(LDPROXYfdf)
}

final <- rbind(LDPROXYfdf,dnSNPfdf)

## for nimpress compatibility
final[is.na(final[,6]),6] <- "NaN"

##################
## WRITE OUTPUT ##
##################

## setup output file
filen <- paste(outdir, "/", output_name2, ".score", sep="")

## title
write(output_name2,file=filen, append=FALSE)
## description
write(gsub("^\"", "", gsub("\"$", "", arguments$description)),file=filen, append=TRUE)
## citation
write(gsub("^\"", "", gsub("\"$", "", arguments$citation)),file=filen, append=TRUE)
## genome version
write(gv,file=filen, append=TRUE)
## ofset
write(arguments$offset,file=filen, append=TRUE)

final = final[order(final[,1], final[,2]),]
write.table(final, file=filen, sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)
 
message("Done")