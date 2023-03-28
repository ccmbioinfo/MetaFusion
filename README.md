To run METAFUSION_forte:
Convert your forte results to a CFF file using the script `convert_fusion_to_cff.R` 

Command to run conversion: 
 ` Rscript convert_fusion_to_cff.R -f /path/to/my/forte/outdir/toplevel -p my_prefix -o /path/to/desired/output/directory`
  

This will output a file in your output directory: /path/to/desired/output/directory/my_prefix.merged.cff

Please ensure that you have a copy of the singularity image when running on HPC, as well as have downloaded the reference DB, described by:  https://github.com/ccmbioinfo/MetaFusion/wiki
You can then run METAFUSION, starting in the "scripts" directory of this repo. Using the RUN_MetaFusion_mskcc.sh 
script, you can simply provide arguments to the command line:
```` 
cd  /path/to/this/metafusion/repo/MetaFusion/scripts/
bsub -R "rusage[mem=16]" -o %J.out singularity exec -B /juno/ -B /path/to/this/metafusion/repo/MetaFusion/ -B /tmp -B /scratch/ metafusion.img bash RUN_MetaFusion_mskcc.sh my_prefix /path/to/this/metafusion/repo/ /path/to/desired/output/directory my_prefix.merged.cff truth_fusions.txt 
````

If you do not have truth fusions, you can run like this:
````
cd  /path/to/this/metafusion/repo/MetaFusion/scripts/
bsub -R "rusage[mem=16]" -o %J.out singularity exec -B /juno/ -B /path/to/this/metafusion/repo/MetaFusion/ -B /tmp -B /scratch/ metafusion.img bash  MetaFusion.sh --outdir /path/to/desired/output/directory/ \
                 --cff my_prefix.merged.cff \
                 --gene_bed /path/to/this/metafusion/repo/MetaFusion/reference_files/ens_known_genes.renamed.ENSG.bed  \
                 --genome_fasta /path/to/this/metafusion/repo/MetaFusion/reference_files/human_g1k_v37_decoy.fasta \
                 --gene_info /path/to/this/metafusion/repo/MetaFusion/reference_files/Homo_sapiens.gene_info  \
                 --num_tools=2  \
                 --recurrent_bedpe /path/to/this/metafusion/repo/MetaFusion/reference_files/blocklist_breakpoints.bedpe \
                 --scripts /path/to/this/metafusion/repo/MetaFusion/scripts
````

See the Wiki for documentation at https://github.com/ccmbioinfo/MetaFusion/wiki

Full text can be found here: https://www.biorxiv.org/content/10.1101/2020.09.17.302307v2

Manuscript published in Bioinformatics: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btab249/6263829

