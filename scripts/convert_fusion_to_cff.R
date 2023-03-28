library(dplyr)
library(data.table)
library(stringr)
cff_format <- c("chr1","pos1","strand1","chr2","pos2","strand2","library","sample_name",
                "sample_type","disease","tool",'split_cnt',"span_cnt","t_gene1","t_area1",
                "t_gene2","t_area2")

library(argparse)

opt = commandArgs(TRUE)

parser=ArgumentParser()
parser$add_argument("-f",'--forte_out_dir',type="character",default = NULL)
parser$add_argument("-m","--mapping_sample_types",type = "character",default = NULL)
parser$add_argument("-g","--filter_to_gene_list",type = "character",default = NULL)
parser$add_argument("-p","--prefix",type = "character", default = NULL)
parser$add_argument("-o",'--output_dir',type="character",default = getwd())
opt=parser$parse_args()



if(!is.null(opt$mapping_sample_types )){
  ### TO DO: Work out mapping sample types.... disease etc,,,
  stop("This has not been built yet")
}

my_samples_dirs <- list.dirs(paste0(opt$forte_out_dir,"/analysis"),recursive = F)

my_tools <- c("arriba","fusioncatcher","starfusion")

if(!all(my_tools %in% basename(list.dirs(my_samples_dirs[1],recursive = F)))){
  stop("A sample is missing an expected completed tool. Check your run")
}

#for testing convience
cff_format_df <- setNames(data.frame(matrix(ncol = length(cff_format), nrow = 0)), cff_format)

make_arriba <- function(sample_file){
  df <- as.data.frame(matrix(ncol = 0,nrow=nrow(sample_file)))

  df$t_gene1 <- sample_file[,1]
  df$t_gene2 <- sample_file$gene2
  df$chr1 <- str_split_fixed(sample_file$breakpoint1,":",2)[,1]
  df$pos1 <- str_split_fixed(sample_file$breakpoint1,":",2)[,2]
  df$chr2 <-  str_split_fixed(sample_file$breakpoint2,":",2)[,1]
  df$pos2 <- str_split_fixed(sample_file$breakpoint2,":",2)[,2]
  df$strand1 <- str_split_fixed(sample_file$`strand1(gene/fusion)`,"/",2)[,1]
  df$strand2 <- str_split_fixed(sample_file$`strand2(gene/fusion)`,"/",2)[,1]
  df$tool <- "arriba"
  df$split_cnt <- ifelse(!is.na(sample_file$split_reads1), sample_file$split_reads1, -1)
  df$span_cnt <- ifelse(!is.na(sample_file$discordant_mates), sample_file$discordant_mates, -1)
  df$t_area1 <- sample_file$site1
  df$t_area2 <- sample_file$site2

  return(df)
}

make_fusioncatcher <- function(sample_file){
  df <- as.data.frame(matrix(ncol = 0,nrow=nrow(sample_file)))

  df$t_gene1 <- sample_file[,1]
  df$t_gene2 <- sample_file[,2]
  df$chr1 <- str_split_fixed(sample_file[,"Fusion_point_for_gene_1(5end_fusion_partner)"],":",3)[,1]
  df$pos1 <- str_split_fixed(sample_file[,"Fusion_point_for_gene_1(5end_fusion_partner)"],":",3)[,2]
  df$chr2 <-  str_split_fixed(sample_file[,"Fusion_point_for_gene_2(3end_fusion_partner)"],":",3)[,1]
  df$pos2 <- str_split_fixed(sample_file[,"Fusion_point_for_gene_2(3end_fusion_partner)"],":",3)[,2]
  df$strand1 <- str_split_fixed(sample_file[,"Fusion_point_for_gene_1(5end_fusion_partner)"],":",3)[,3]
  df$strand2 <- str_split_fixed(sample_file[,"Fusion_point_for_gene_2(3end_fusion_partner)"],":",3)[,3]
  df$tool <- "fusioncatcher"
  df$split_cnt <- ifelse(!is.na(sample_file$Spanning_unique_reads), sample_file$Spanning_unique_reads, -1)
  df$span_cnt <- ifelse(!is.na(sample_file$Spanning_pairs), sample_file$Spanning_pairs, -1)
  df$t_area1 <- sample_file$Predicted_effect
  df$t_area2 <- sample_file$Predicted_effect


  return(df)


}

make_starfusion <- function(sample_file){

  df <- as.data.frame(matrix(ncol = 0,nrow=nrow(sample_file)))

  df$t_gene1 <- str_split_fixed(sample_file$LeftGene,"\\^",2)[,1]
  df$t_gene2 <- str_split_fixed(sample_file$RightGene,"\\^",2)[,1]
  df$chr1 <- str_replace(str_split_fixed(sample_file$LeftBreakpoint,":",3)[,1],"chr","")
  df$pos1 <- str_split_fixed(sample_file$LeftBreakpoint,":",3)[,2]
  df$chr2 <-    str_replace(str_split_fixed(sample_file$RightBreakpoint,":",3)[,1],"chr","")
  df$pos2 <-  str_split_fixed(sample_file$RightBreakpoint,":",3)[,2]
  df$strand1 <- str_split_fixed(sample_file$LeftBreakpoint,":",3)[,3]
  df$strand2 <- str_split_fixed(sample_file$RightBreakpoint,":",3)[,3]
  df$tool <- "starfusion"
  df$split_cnt <- ifelse(!is.na(sample_file$JunctionReadCount), sample_file$JunctionReadCount, -1)
  df$span_cnt <- ifelse(!is.na(sample_file$SpanningFragCount), sample_file$SpanningFragCount, -1)
  df$t_area1 <- sample_file$PROT_FUSION_TYPE
  df$t_area2 <- sample_file$PROT_FUSION_TYPE

  return(df)


}



make_sample_cff <- function(sample_dir){
  sample <- basename(sample_dir)
  print(paste0("Collapsing sample:", sample))
  sample_cff <- do.call(rbind,lapply(my_tools, function(tool){
    tool_cff <-  setNames(data.frame(matrix(ncol = length(cff_format), nrow = 0)), cff_format)
    if(tool == "arriba"){
      sample_file <- fread(paste0(sample_dir,"/arriba/",list.files(paste0(sample_dir,"/arriba/"), pattern = ".fusions.tsv")),data.table = F)

      if(nrow(sample_file) > 0){
        tool_cff <- rbind(tool_cff,make_arriba(sample_file))
      }
    } else if( tool == "fusioncatcher") {
      sample_file <- fread(paste0(sample_dir,"/fusioncatcher/",list.files(paste0(sample_dir,"/fusioncatcher/"), pattern = ".fusion-genes.hg19.txt")),data.table = F)

      if(nrow(sample_file) > 0){
        tool_cff <-rbind(tool_cff,make_fusioncatcher(sample_file))
      }
    } else if (tool == "starfusion"){
      sample_file <- fread(paste0(sample_dir,"/starfusion/",list.files(paste0(sample_dir,"/starfusion/"), pattern = ".abridged.coding_effect.tsv")),data.table = F)

      if(nrow(sample_file) > 0){
        tool_cff <-rbind(tool_cff,make_starfusion(sample_file))
      }
    }
    if(nrow(tool_cff) > 0 ){
      tool_cff$sample_name <- sample
      tool_cff$library <- "RNA"
      tool_cff$sample_type <- "Tumor"
      tool_cff$disease <- NA
    }
    return(tool_cff)

  }))
  return(sample_cff[,cff_format])
}



my_cff <- do.call(rbind,lapply(my_samples_dirs, make_sample_cff))
my_cff$strand1[my_cff$strand1 == "."] <- NA
my_cff$strand2[my_cff$strand2 == "."] <- NA

write.table(my_cff,paste0(opt$output_dir,"/",opt$prefix,".merged.cff"),sep = "\t",row.names = F,quote =  F,col.names =F)
