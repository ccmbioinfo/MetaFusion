library(dplyr)
library(data.table)
library(stringr)


gene_info <- fread("/work/ccs/pintoa1/references/meta_fusion_bed_generation/Homo_sapiens_gene_info_may29th2023_ncbi.txt")
### 191478 genes
our_gtf <- rtracklayer::import("/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf")

our_gtf <- as.data.frame(our_gtf)



ensemblToGeneName <- unique(our_gtf[,c("gene_id","gene_name")])
length(unique(ensemblToGeneName$gene_id)) # 63677 gene ids
length(unique(ensemblToGeneName$gene_name))# 56638 gene names

gene_info$NCBI_accepted_Symbol <- gene_info$Symbol
gene_info$gene_id_ensg <- NA
### first identify rows with ensembl gene ids

gene_info$gene_id_ensg[grepl("Ensembl:",gene_info$dbXrefs)] <- regmatches(gene_info$dbXrefs[grepl("Ensembl:",gene_info$dbXrefs)],regexpr("ENSG[0-9]+",gene_info$dbXrefs[grepl("Ensembl:",gene_info$dbXrefs)] ))
table(is.na(gene_info$gene_id_ensg))
# 
# FALSE   TRUE 
# 36402 155076 
#36402 ensg ids exist in NCBI

## Are all of those ensgids in our gtf
clean_mapping <- merge(gene_info,ensemblToGeneName,by.x = "gene_id_ensg",by.y = "gene_id",all = F)
#33726 map cleanly
## however multple rows can exist for a single gene_id (some gene ids belong to two rows. I.E:
# clean_mapping[clean_mapping$gene_id_ensg == "ENSG00000004866",]
# gene_id_ensg #tax_id GeneID  Symbol LocusTag                               Synonyms                                                                      dbXrefs chromosome map_location
# 1: ENSG00000004866    9606   7982     ST7        - ETS7q|FAM4A|FAM4A1|HELG|RAY1|SEN4|TSG7 MIM:600833|HGNC:HGNC:11351|Ensembl:ENSG00000004866|AllianceGenome:HGNC:11351          7       7q31.2
# 2: ENSG00000004866    9606  93655 ST7-OT3        -                  NCRNA00026|ST7|ST7OT3            HGNC:HGNC:16045|Ensembl:ENSG00000004866|AllianceGenome:HGNC:16045          7       7q31.2
# description   type_of_gene Symbol_from_nomenclature_authority Full_name_from_nomenclature_authority Nomenclature_status
# 1: suppression of tumorigenicity 7 protein-coding                                ST7       suppression of tumorigenicity 7                   O
# 2:    ST7 overlapping transcript 3          ncRNA                            ST7-OT3          ST7 overlapping transcript 3                   O
# Other_designations Modification_date Feature_type Symbol_v75 NCBI_accepted_Symbol gene_name
# 1: suppressor of tumorigenicity 7 protein|family with sequence similarity 4, subfamily A, member 1|suppression of tumorigenicity 7 (breast)          20230329            -         NA                  ST7       ST7
# 2:          ST7 overlapping transcript 3 (non-coding RNA)|ST7 overlapping transcript 3 (non-protein coding)|suppression of tumorigenicity 7          20230329            -         NA              ST7-OT3       ST7
# 

missing_mapping <- gene_info[is.na(gene_info$gene_id_ensg) | !gene_info$gene_id_ensg %in% ensemblToGeneName$gene_id, ]
### 157752 genes missing mapping to v75

need_to_find <- unique(ensemblToGeneName[!ensemblToGeneName$gene_id %in% clean_mapping$gene_id_ensg,"gene_name"])

### match via gene name
### want to reacquire 25514 unqiue gene names that do not map from v75 to gene info via ENSG id
## problem genes are genes which do not exist in our gene info file but exist in v75 gtf.
### we will create an entry for problem genes
problem_genes <- c()
try_to_match_symbols_from_gene_info_to_v75 <- lapply(need_to_find, function(gene){

  syn_str <- paste("^",gene,"[\\|]|^",gene,"$","|[\\|]",gene,"[\\|]|[\\|]",gene,"$", sep="")
  
  if(any(missing_mapping$Symbol == gene)){
    tmp <- missing_mapping[missing_mapping$Symbol == gene,]
    tmp$gene_name <- gene
    tmp$gene_id_ensg <-paste0(ensemblToGeneName[ensemblToGeneName$gene_name == gene,"gene_id"],collapse = "|")
    return(tmp)
  } else if (any(grepl(syn_str,missing_mapping$Synonyms))){
    tmp <- missing_mapping[grepl(syn_str,missing_mapping$Synonyms),]
    tmp$Synonyms <- paste0(tmp$Synonyms,"|",tmp$Symbol)
    tmp$gene_name <- gene
    tmp$gene_id_ensg <-paste0(ensemblToGeneName[ensemblToGeneName$gene_name == gene,"gene_id"],collapse = "|")
    return(tmp)
  }else{
    problem_genes <<- c(problem_genes,gene)
    df <- data.frame(matrix(ncol = ncol(missing_mapping), nrow = 1))
    colnames(df) <- colnames(missing_mapping)
    df$gene_name <- gene
    df$gene_id_ensg <- paste0(ensemblToGeneName[ensemblToGeneName$gene_name == gene,"gene_id"],collapse = "|")
    return(df)
  }
  
  
  
})

try_to_match_symbols_from_gene_info_to_v75 <- do.call(rbind,try_to_match_symbols_from_gene_info_to_v75)
problem_genes <- unique(problem_genes)
#20307 genes do not have a match in the NCBI info file. These genes are added to the info file, but only have the v75 information in them


clean_mapping <- rbind(clean_mapping,try_to_match_symbols_from_gene_info_to_v75)
length(unique(ensemblToGeneName$gene_name)) ## 56638 unique gene NAMES trying to reaquire
length(unique(clean_mapping$gene_name))### all reaquired 56638 .... 59270 total rows

### ADD NCBISYMBOL TO SYNONMYS
### if no NCBI symbol/synonym populate
clean_mapping$Synonyms[is.na(clean_mapping$Synonyms)] <- "-"
mismatch_ncbi_v75 <- clean_mapping[!is.na(clean_mapping$NCBI_accepted_Symbol) & clean_mapping$NCBI_accepted_Symbol != clean_mapping$gene_name,]


### add NCBI symbol to synonyms
for( i in 1:nrow(mismatch_ncbi_v75)){
  gene <- mismatch_ncbi_v75$NCBI_accepted_Symbol[i]
  syn_str <- paste("^",gene,"[\\|]|^",gene,"$","|[\\|]",gene,"[\\|]|[\\|]",gene,"$", sep="")
  ### if 
  if(!grepl(syn_str,mismatch_ncbi_v75$Synonyms[i])){
    clean_mapping$Synonyms[clean_mapping$gene_id_ensg == mismatch_ncbi_v75$gene_id_ensg[i]] <- paste0(gene,"|",clean_mapping$Synonyms[clean_mapping$gene_id_ensg == mismatch_ncbi_v75$gene_id_ensg[i]])
  }
}

clean_mapping$Symbol <- clean_mapping$gene_name
clean_mapping$gene_name <- NULL
#### FROM HERE ON, GENE INFO IS IN METAFUSION ACCEPTED FORMAT TO CONVERT SYNONYMS TO V75 SYMBOL
write.table(clean_mapping,"/work/ccs/pintoa1/references/meta_fusion_bed_generation/gene_info.v75.june222023",row.names = F,sep = "\t")



star_fusion_ref <- rtracklayer::import("/work/taylorlab/cmopipeline/rnaseq_reference/GRCh37/starfusion/ctat_genome_lib_build_dir/ref_annot.gtf")
fusioncatcher_ref <- rtracklayer::import("/work/taylorlab/cmopipeline/rnaseq_reference/GRCh37/fusioncatcher/organism.gtf")

star_fusion_ref_df <- as.data.frame(star_fusion_ref)
fusioncatcher_ref_df <- as.data.frame(fusioncatcher_ref)
star_fusion_ref_name_to_id <- unique(star_fusion_ref_df[,c("gene_id","gene_name")])
## remove version 
colnames(star_fusion_ref_name_to_id) <- c("starfusion_gene_id","starfusion_gene_name")

star_fusion_ref_name_to_id$starfusion_gene_id_no_version <- str_split_fixed(star_fusion_ref_name_to_id$starfusion_gene_id,"\\.",n=2)[,1]


### match starfusion and v75 via ENSGid
v75_to_starfusion <- merge(ensemblToGeneName,star_fusion_ref_name_to_id,by.x = "gene_id",by.y ="starfusion_gene_id_no_version",all= F)
# 63232     
### if a gene id did not successfully map via geneid... 
missing_ensg_ids_starFusion <- star_fusion_ref_name_to_id[!star_fusion_ref_name_to_id$starfusion_gene_id %in% v75_to_starfusion$starfusion_gene_id,] ## 51
# starfusion_gene_id starfusion_gene_name starfusion_gene_id_no_version
# 2588077  ENSGR0000228572.2      LL0YNC03-29C1.1               ENSGR0000228572
# 2588081  ENSGR0000182378.8               PLCXD1               ENSGR0000182378
# 2588228  ENSGR0000178605.8               GTPBP6               ENSGR0000178605
###merge by gene name
missing_ensg_id_match_via_name <-  merge(ensemblToGeneName,missing_ensg_ids_starFusion,by.x = "gene_name",by.y ="starfusion_gene_name",all= F)## 49. all but 2 have a gene name match. 

still_missing <- missing_ensg_ids_starFusion[!missing_ensg_ids_starFusion$starfusion_gene_id %in% missing_ensg_id_match_via_name$starfusion_gene_id,]
# starfusion_gene_id starfusion_gene_name starfusion_gene_id_no_version
# 2815615         IGH.g@-ext             IGH@-ext                           IGH
# 2815616        IGH-.g@-ext            IGH-@-ext                          IGH-

### ignored IGH....
#check if all gene names match between gene ids
indirect_gene_name_match_starfusion_v75 <- v75_to_starfusion[v75_to_starfusion$gene_name != v75_to_starfusion$starfusion_gene_name,] 
### 498 gene names not matching 
### add these 498 to clean_mapping synonyms

for(i in 1:nrow(indirect_gene_name_match_starfusion_v75)){
  gene <- indirect_gene_name_match_starfusion_v75$starfusion_gene_name[i]
  syn_str <- paste("^",gene,"[\\|]|^",gene,"$","|[\\|]",gene,"[\\|]|[\\|]",gene,"$", sep="")
  
  if(any(!grepl(syn_str,clean_mapping$Synonyms[clean_mapping$Symbol == indirect_gene_name_match_starfusion_v75$gene_name[i]]))){
    clean_mapping$Synonyms[clean_mapping$Symbol == indirect_gene_name_match_starfusion_v75$gene_name[i]] <- ifelse( clean_mapping$Synonyms[clean_mapping$Symbol == indirect_gene_name_match_starfusion_v75$gene_name[i]] == '-',
                                                                                                                    gene ,
                                                                                                                    paste0(gene,"|",clean_mapping$Synonyms[clean_mapping$Symbol == indirect_gene_name_match_starfusion_v75$gene_name[i]]))
  }
}


#### any star fusion and v75 mismatches that do not already exist in synoyms have been added at this point. 
write.table(clean_mapping,"/work/ccs/pintoa1/references/meta_fusion_bed_generation/gene_info.v75.sf.june222023",row.names = F,sep = "\t")

clean_mapping <- fread("/work/ccs/pintoa1/references/meta_fusion_bed_generation/gene_info.v75.sf.june222023")
# 59270 mappings
#> length(unique(clean_mapping$Symbol))
# 56638 unique gene names in v75

### NOW DO IT FOR FUSIONCATCHER

### REDO FUSIONCATCHER MERGE WE ARE MISSING A CANCER GENE FROM ONCOKB
## we can reaquire itwith ENST/transcript_id
# cancer_genes[cancer_genes$`GRCh37 Isoform` %in% fusioncatcher_ref_df$transcript_id[fusioncatcher_ref_df$gene_name %in% problem_genes_fc],]
# Hugo Symbol Entrez Gene ID  GRCh37 Isoform GRCh37 RefSeq  GRCh38 Isoform GRCh38 RefSeq # of occurrence within resources (Column D-J)
# 1:      PTP4A1           7803 ENST00000370651   NM_003463.4 ENST00000626021   NM_003463.4                                             3
# OncoKB Annotated Is Oncogene Is Tumor Suppressor Gene MSK-IMPACT MSK-HEME FOUNDATION ONE FOUNDATION ONE HEME Vogelstein SANGER CGC(05/30/2017)
# 1:              Yes          No                       No        Yes      Yes             No                  No         No                     No
# Gene Aliases
# 1: PRL-1, PTPCAAX1
fusioncatcher_ref_name_to_id <- unique(fusioncatcher_ref_df[,c("gene_id","gene_name")]) #60620 unique mappings,59409 gene names
colnames(fusioncatcher_ref_name_to_id) <- c("FC_gene_id","FC_gene_name")

#### keep only genes which have differing gene names. attempt to map to v75
fusioncatcher_ref_name_to_id <- fusioncatcher_ref_name_to_id[!fusioncatcher_ref_name_to_id$FC_gene_name %in% clean_mapping$Symbol,]
###25027 gene names do not exist as Symbol

### keep problem genes which still have no mapping
problem_genes_fc <- c()
for(i in 1:nrow(fusioncatcher_ref_name_to_id)){
  gene <- fusioncatcher_ref_name_to_id$FC_gene_name[i]
  syn_str <- paste("^",gene,"[\\|]|^",gene,"$","|[\\|]",gene,"[\\|]|[\\|]",gene,"$", sep="")
  
  ## if gene doesnt exist in synonmys retain it
  if(!any(grepl(syn_str,clean_mapping$Synonyms))){
    problem_genes_fc <<- c(gene,problem_genes_fc)
  }
}

### 21139 gene names which do not exist within synonmys 
problem_genes_fc_df <-fusioncatcher_ref_name_to_id[fusioncatcher_ref_name_to_id$FC_gene_name %in% problem_genes_fc,]
### match starfusion and v75 via ENSGid
v75_to_fusioncatcher <- merge(ensemblToGeneName,problem_genes_fc_df,by.x = "gene_id",by.y ="FC_gene_id",all= F)
# 14443 matches between v75 and fusion catcher on GENE ID
### add fusioncatchers expectd gene name to synonmys 
for(i in 1:nrow(v75_to_fusioncatcher)){
  gene <- v75_to_fusioncatcher$FC_gene_name[i]
  syn_str <- paste("^",gene,"[\\|]|^",gene,"$","|[\\|]",gene,"[\\|]|[\\|]",gene,"$", sep="")
  
  if(any(!grepl(syn_str,clean_mapping$Synonyms[clean_mapping$Symbol == v75_to_fusioncatcher$gene_name[i]]))){
    clean_mapping$Synonyms[clean_mapping$Symbol == v75_to_fusioncatcher$gene_name[i]] <- ifelse( clean_mapping$Synonyms[clean_mapping$Symbol == v75_to_fusioncatcher$gene_name[i]] == '-',
                                                                                                             gene ,
                                                                                                             paste0(gene,"|",clean_mapping$Synonyms[clean_mapping$Symbol == v75_to_fusioncatcher$gene_name[i]]))
  }
}


missing_ensg_ids_fusioncatcher <- problem_genes_fc_df[!problem_genes_fc_df$FC_gene_id %in% v75_to_fusioncatcher$gene_id,] ## 6707 
#### there are some gene names which have more than 1 ENSG id. see if any of the other ENGs pairs matched for a gene name
## removed 10 because they already mapped EX:
# 863279 ENSG00000275631             U1
# 2894643 ENSG00000276932          Y_RNA
# 2894660 ENSG00000277428          Y_RNA
# 2894896 ENSG00000277374             U1
# 2895086 ENSG00000275405             U1
# 2895089 ENSG00000275987             U1
###  these gene names already exist within the v75_to_fusioncatcher dataframe as a different gene_id merge. 
missing_ensg_ids_fusioncatcher <- missing_ensg_ids_fusioncatcher[!missing_ensg_ids_fusioncatcher$FC_gene_name %in% v75_to_fusioncatcher$FC_gene_name,] #6696 unqiue genes
### match these via transcript id
ensemblTranscriptToGeneName <- unique(our_gtf[,c("transcript_id","gene_name")])
missing_ensg_ids_fusioncatcher_transcript <- unique(fusioncatcher_ref_df[fusioncatcher_ref_df$gene_name %in% missing_ensg_ids_fusioncatcher$FC_gene_name,c("transcript_id","gene_name")])
colnames(missing_ensg_ids_fusioncatcher_transcript) <- c("FC_transcript_id","FC_gene_name")
missing_ensg_ids_fusioncatcher_transcript <- missing_ensg_ids_fusioncatcher_transcript[!is.na(missing_ensg_ids_fusioncatcher_transcript$FC_transcript_id),]
v75_t_to_fc_t <-  merge(ensemblTranscriptToGeneName,missing_ensg_ids_fusioncatcher_transcript,by.x = "transcript_id",by.y ="FC_transcript_id",all = F) 
# length(unique(v75_t_to_fc_t$gene_name))
# [1] 309
# > length(unique(v75_t_to_fc_t$FC_gene_name))
# 354 gene names rescued by mapping via transcript id

### add these gene names to synonyms
for(i in 1:nrow(v75_t_to_fc_t)){
  gene <- v75_t_to_fc_t$FC_gene_name[i]

  clean_mapping$Synonyms[clean_mapping$Symbol == v75_t_to_fc_t$gene_name[i]] <- ifelse( clean_mapping$Synonyms[clean_mapping$Symbol == v75_t_to_fc_t$gene_name[i]] == '-',
                                                                                                 gene ,
                                                                                                 paste0(gene,"|",clean_mapping$Synonyms[clean_mapping$Symbol == v75_t_to_fc_t$gene_name[i]]))

}



still_missing_ensg_id_match_even_after_transcript <- missing_ensg_ids_fusioncatcher[!missing_ensg_ids_fusioncatcher$FC_gene_name %in% c(v75_t_to_fc_t$FC_gene_name),]
###6351 rows, 6342  gne namesSTILL MISSING EVEN AFTER MATCH VIA TRANSCRIPT 
# ## t FC_gene_id FC_gene_name
# 26  ENSG00000278267    MIR6859-1
# 152 ENSG00000279928   FO538757.1
# 159 ENSG00000279457       WASH9P
# 171 ENSG00000273874    MIR6859-2
# 437 ENSG00000278791   AC114498.2
# 642 ENSG00000285268   AL669831.6
table(fusioncatcher_ref_df$gene_biotype[fusioncatcher_ref_df$gene_name %in% still_missing_ensg_id_match_even_after_transcript$FC_gene_name & !duplicated(fusioncatcher_ref_df$gene_id)])

                     IG_pseudogene                          IG_V_gene                    IG_V_pseudogene                             lncRNA 
                                 1                                  7                                  5                               3445 
                             miRNA                           misc_RNA             polymorphic_pseudogene               processed_pseudogene 
                               431                                223                                  4                                364 
                    protein_coding                         pseudogene                           ribozyme                               rRNA 
                               260                                 16                                  1                                 24 
                   rRNA_pseudogene                             scaRNA                             snoRNA                              snRNA 
                                 7                                  2                                 14                                  1 
                              sRNA                                TEC                          TR_D_gene                          TR_J_gene 
                                 4                                957                                  1                                  6 
  transcribed_processed_pseudogene     transcribed_unitary_pseudogene transcribed_unprocessed_pseudogene    translated_processed_pseudogene 
                                21                                  9                                 89                                  1 
                unitary_pseudogene             unprocessed_pseudogene 
                                21                                437 
write.table(still_missing_ensg_id_match_even_after_transcript, "/work/ccs/pintoa1/references/meta_fusion_bed_generation/fusion_catcher_genes_missing_from_gene_info_metafusion_june262023",sep ="\t",quote = F)
write.table(clean_mapping,"/work/ccs/pintoa1/references/meta_fusion_bed_generation/gene_info.v75.sf.fc.june262023",sep = "\t",quote =  F)
