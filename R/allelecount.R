################################################
# Allele counting functions - SNV, indel and SV
################################################
#' Run alleleCount
### modified by Shaghayegh Soudi
#' Dump allele counts from maf file
#'
#' Dump allele counts stored in the sample columns of the maf file. Output will go into a file
#' supplied as tumour_outfile and optionally normal_outfile. It will be a fully formatted
#' allele counts file as returned by alleleCounter.
#' @param maf_file The maf file to read in
#' @param tumour_outfile File to save the tumour counts to
#' @param normal_outfile Optional parameter specifying where the normal output should go
#' @param refence_genome Optional parameter specifying the reference genome build used
#' @param samplename Optional parameter specifying the samplename to be used for matching the right column in the maf
#' @author sd11
#' @export




dumpCounts.maf = function(maf_file, tumour_outfile,  normal_outfile=NA, samplename=NA, dummy_alt_allele=NA, dummy_ref_allele=NA) {
   
   for (maf_file in maf_file) { 
   maf = read.delim(maf_file, skip = 1, header = TRUE, sep = "\t")
   maf_SNV<-maf[maf$Variant_Type == "SNP",]
   counts_maf = maf_SNV[,c("t_ref_count","t_alt_count")]
   rownames(counts_maf)<-paste(maf_SNV$Chromosome,maf_SNV$Start_Position, sep = ":")
   allele.ref = as.character(maf_SNV[,"Reference_Allele"])
   allele.alt = as.character(maf_SNV[,"Tumor_Seq_Allele2"])


   count.ref = counts_maf[,1]
   count.alt = counts_maf[,2]
   #allele.ref = as.character(VariantAnnotation::ref(v))
   #allele.alt = unlist(lapply(VariantAnnotation::alt(v), function(x) { as.character(x[[1]]) }))
   
   output = array(0, c(length(allele.ref), 4))
   nucleotides = c("A", "C", "G", "T")
   # Propagate the alt allele counts
   nucleo.index = match(allele.alt, nucleotides)
  
   for (i in 1:nrow(output)) {
    output[i,nucleo.index[i]] = count.alt[i]
   } ## for i loop
  
  # Propagate the reference allele counts
  nucleo.index = match(allele.ref, nucleotides)
  for (i in 1:nrow(output)) {
    output[i,nucleo.index[i]] = count.ref[i]
  }  ## second for i loop

  if (nrow(maf_SNV)==0) {
    output = data.frame(matrix(ncol = 7, nrow = 0))
   } else {

    output<-cbind(as.character(maf_SNV$Chromosome),maf_SNV$Start_Position,output, rowSums(output))
    colnames(output) = c("#CHR","POS","Count_A","Count_C","Count_G","Count_T","Good_depth")
    output[,1]<-gsub("chr","",output[,1])  
    output<-data.frame(output)
   write.table(output, col.names=T, quote=F, row.names=F, file=tumour_outfile, sep="\t")

   }                      ### adjusted by Shaghayegh
   
   }

}

###############################################
###############################################



