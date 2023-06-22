library(dplyr)
library(data.table)
library(stringr)
gtf <- rtracklayer::import('/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf')
gtf_df <- as.data.frame(gtf)


# Utilized gtf from igenomes for FORTE This corresponds to GRCh37 ensembl 75
# Add introns to gtf, convert to gff3
# bsub -R "rusage[mem=64]" -o add_introns_agat_%J.out singularity exec -B /juno/ \\
# -B /tmp -B /scratch/ docker://quay.io/biocontainers/agat:0.8.0--pl5262hdfd78af_0  \\
# /bin/bash -c "agat_sp_add_introns.pl -g /juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf\\
# -o genes.INTRONS.gff3"
# gff2bed < genes.INTRONS.gff3 > genes.INTRONS.agat.bed


total.introns.bed <- fread(file="/work/ccs/pintoa1/references/meta_fusion_bed_generation/genes.INTRONS.agat.bed", header = FALSE, stringsAsFactors = F, sep="\t", na.strings = "",data.table = F)
colnames(total.introns.bed) <- c("chr","start","end","gene_id","tmp","strand","gene_biotype","type","V9","description")
total.introns.bed$transcript_id <- gsub("\\;.*","",str_split_fixed(total.introns.bed$description,"transcript_id=",n=2)[,2])
total.introns.bed$gene_name <-gsub("\\;.*","",str_split_fixed(total.introns.bed$description,"gene_name=",n=2)[,2])


summarse_t_type <- total.introns.bed %>% group_by(transcript_id) %>%   summarise(total_cds = length(which(type == "UTR")), total_UTRS =length(which(type == "UTR")) , total_stops = length(which(type == "stop_codon")),total_starts = length(which(type == "start_codon")), both = all(c("stop_codon","start_codon") %in% type))
table(summarse_t_type$total_stops, summarse_t_type$total_starts,summarse_t_type$both)
# 
# , ,  = FALSE
# 
# 
#       0      1      2
# 0 141739    292      0
# 1    212      0      0
# 2      0      0      0
# 
# , ,  = TRUE
# 
# 
#        0      1      2
# 0      0      0      0
# 1      0  72385    181
# 2      0    314      0

table(summarse_t_type[summarse_t_type$total_UTRS > 0,c("total_stops" ,"total_starts")])
#               total_starts
# total_stops     0       1       2
# 0               28225   207     0
# 1               105     70813   181
# 2               0       314     0

### cannnot use start or stop codons for 28225 transcripts which have UTRS
table(summarse_t_type$total_cds > 0, summarse_t_type$total_UTRS)
# 
#        0      1      2      3      4      5      6      7      8      9     10     11     12     13     14     15     16     17     18     19     20     21     22     23     24     25     26     27     28
# FALSE 115278      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0
# TRUE       0  12104  38786  25936  10858   4750   2344   1364    834    617    397    328    276    199    195    144    102    102     83     82     81     46     28     32     15     16     18     12     13
# 
#         29     30     31     32     33     34     35     36     37     38     39     40     41     42     44     45     47     52     54     55     57     88
# FALSE      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0      0
# TRUE      14     13      6      7     10      6      5      5      1      2      1      2      2      1      1      1      1      1      1      1      1      1
### ALL TRANSCRIPTS WITH UTRS HAVE AT LEAST ONE CDS REGION, USING THAT FOR DETERMING UTR 3' and 5'


transcript_ids <- unique(total.introns.bed$transcript_id)
file.to_write <- "/work/ccs/pintoa1/references/meta_fusion_bed_generation/cleaned_metafusion_v75_gene.bed"

if(file.exists(file.to_write) ) {file.remove(file.to_write)}

#START CLOCK: THE INDEXING TAKES A LONG TIME, LIKE 5 HOURS
ptm <- proc.time()

# Index each transcript feature, incrementing when an intron is passed
## metafusion expects exon count 0 to (N(exons)-1)
## Forward strand: Exon 0 == Exon 1
### Reverse strand: Exon 0 == LAST EXON IN TRANSCRIPT
for (id in transcript_ids){
  transcript <- total.introns.bed[total.introns.bed$transcript_id == id,]
  # Remove exons if coding gene, since "exon" and "CDS" are duplicates of one another
  if ("CDS" %in% transcript$type){
    ## The following gene biotypes will have their "exon" values removed (assumption that CDS and exon are exactly the same for these biotypes)
    ### I found this is not always true... EX: "ENST00000390427" is a TR_V_gene which has differing exon and CDS start and end points.... (still following their logic tho)
    ####unique(gtf_df$gene_biotype[gtf_df$type == "CDS"])
    #[1] "protein_coding"         "polymorphic_pseudogene" "TR_V_gene"              "TR_D_gene"              "TR_J_gene"              "TR_C_gene"              "IG_C_gene"              "IG_J_gene"             
    # [9] "IG_D_gene"              "IG_V_gene"      
    
    transcript <- transcript[!transcript$type == "exon",]
  }
  # Order features by increasing bp 
  transcript <- transcript[order(transcript$start, decreasing = FALSE),]
  # Index features
  idx <- 0
  for (i in 1:nrow(transcript)){
    transcript$idx [i]<- idx
    if (transcript$type[i] == "intron"){
      idx <- idx + 1
    }
  }
  # REFORMAT TRANSCRIPT
  #Change strand info (+ --> f, - --> r)
  if (unique(transcript$strand) == "+"){
    transcript$strand <- 'f'
  } else if  (unique(transcript$strand) == "-"){
    transcript$strand <- 'r'
  } else {
    errorCondition("Strand info for this transcript is inconsistent")
  }
  #Add "chr" prefix to chromosomes
  transcript$chr <- sapply("chr", paste0,  transcript$chr)
  #Change CDS --> cds ### IF A TRANSCRIPT LACKS "CDS" THIS LINE WILL DO NOTHING, Changing exon values to UTRs later 
  if ("CDS" %in% unique(transcript$type)){transcript[transcript$type == "CDS",]$type <- "cds"}
  
  ## DETERMING UTR3 and UTR5
  
  
  ### INSTEAD OF START AND STOP, USE CDS LOCATIONS AND STRAND INFORMATION.....
  if ("UTR" %in% unique(transcript$type)){
    if( unique(transcript$strand) == "f"){
      #Forward strand 
      start_coding <- min(transcript[transcript$type == "cds","start"])
      stop_coding <-  max(transcript[transcript$type == "cds","end"])
      transcript$type[transcript$end <= start_coding &  transcript$type == "UTR"] <- "utr5"
      transcript$type[transcript$start >= stop_coding & transcript$type == "UTR"] <- "utr3"
    }else {
      start_coding <- max(transcript[transcript$type == "cds","end"])
      stop_coding <- min(transcript[transcript$type == "cds","start"])
      transcript$type[transcript$end <= start_coding &  transcript$type == "UTR"] <- "utr3"
      transcript$type[transcript$start >= stop_coding & transcript$type == "UTR"] <- "utr5"
    }
    
  }
  
  transcript <- transcript[,c("chr", "start", "end", "transcript_id", "type", "idx", "strand", "gene_name", "gene_id" )]
  write.table(transcript, file.to_write, append=TRUE, sep="\t", quote=F,  row.names=F, col.names=F)
}

time <- proc.time() - ptm
time
# 
# user    system   elapsed 
# 16657.116    32.227 16741.382 


new.bed <- fread(file.to_write,data.table = F)
colnames(new.bed) <- c("chr","start","end","transcript_id","type","idx","strand","gene_name","gene_id")

# 
# table(new.bed$type)
# 
# cds           exon           gene         intron Selenocysteine    start_codon     stop_codon     transcript           utr3           utr5 
# 791871         376580          63677        1091211            114          73353          73406         215117         142873         161162 
# 




#### Any exon that remains after teh cds change, is likely and untranslated region. change below

# Basically, subfeatures which are "exon" need to be changed (i.e. exon --> utr3/utr5)
#Forward strand
new.bed$type[new.bed$strand == "f" &  new.bed$type == "exon" ] <- "utr5"
#Reverse strand
new.bed$type[new.bed$strand == "r" &  new.bed$type == "exon"]<- "utr3"

table(new.bed$type)
# 
# cds           gene         intron Selenocysteine    start_codon     stop_codon     transcript           utr3           utr5 
# 791871          63677        1091211            114          73353          73406         215117         327592         353023 

expected_types <- c("cds","intron","utr3","utr5")
new.bed.ready <- new.bed[new.bed$type %in% c(expected_types),]

write.table(new.bed.ready, "/work/ccs/pintoa1/references/meta_fusion_bed_generation/v75_gene.bed",  sep="\t", quote=F,  row.names=F, col.names=F)
