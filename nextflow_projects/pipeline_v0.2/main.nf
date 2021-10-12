#!/usr/bin/env nextflow



/*
run:
sudo /home/zawada/tools/nextflow/nextflow run /home/zawada/projects/nextflow_projects/pipeline_v0.1/main.nf
in:
***/results_v01

run in ExomeServer:
sudo /home/zawada/tools/nextflow/nextflow run /home/zawada/projects/nextflow_projects/pipeline_v0.2/main.nf
or:
sudo /home/zawada/tools/nextflow/nextflow run /home/zawada/projects/nextflow_projects/pipeline_v0.1_cg_files/main_real_files_from_cg.nf -with-trace

Check after moving to server:
-if publishDir and storeDir works
*/


params.trimmomatic = "***/trimmomatic/Trimmomatic-0.38/trimmomatic-0.38.jar"
params.picard = "***/picardtools/picard-2.18.14/build/libs/picard.jar"
params.BWA = "***/bwa/bwa-0.7.17/bwa"
params.fastqc = "***/fastqc/FastQC/fastqc"
params.gatk4 = "***/gatk4/gatk-4.1.4.0/gatk"
params.mosdepth = "***/miniconda3/bin/mosdepth"
params.bedtools = "***/bedtools2/bin/./bedtools"
params.vcfTodb = "***/vcf_upload/vcf_to_db.py"

//params.reads = "***/fastq_test_data/fastq_NF_tests/*_R{1,2}.fastq" //3 pairs of normal size fastq
params.reads = "***/fastq/*_R{1,2}.fastq" //3 pairs of small size test fastq files
//params.reads = "***/fastq_test_data/fastq_CG/*_R{1,2}.fastq"
params.samplesConfig="***/current_run/samplesToRun_*.csv"

params.hg38 = "***/hg38/GRCh38.primary_assembly.genome.fa"
params.hg38VEP = "***/genomes/hg38/GRCh38.primary_assembly.genome.fa"
params.dbsnp38 = "***/dbsnp/hg38/GCF_000001405.38.chr_renamed.vcf.gz"
params.bedSureSelect = "***/hg38/sureselect_all_exon_v7_hs_hg38/S31285117_Regions_genes.bed"
params.bedSureSelect25bpExtend = "***/genomes/bed_hg38/agilent_bed/S31285117_agi_25bp_extended.bed"

params.hapmap = "***/hapmap_3.3.hg38_noChr.vcf.gz"
params.omni = "***/1000G_omni2.5.hg38_noChr.vcf.gz"
params.g1000G = "***/1000G_phase1.snps.high_confidence.hg38_noChr.vcf.gz"
params.mills = "***/Mills_and_1000G_gold_standard.indels.hg38_noChr.vcf.gz"

params.snpeff = "***/snpeff/snpEff/"
params.dbnsfp = "***/dbnsfp/dbNSFP4.0a.gz"
params.vep = "***/ensembl-vep/vep"

params.lowCpus = 10
params.highCpus = 30

///////////////////old version
/*
Channel
    .fromFilePairs(params.reads) //flat:true - consider this parameter to get output as flat list: 
    //[index, ***/fastq/index_R1.fastq, ***/fastq/index_R2.fastq]
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}"}
    //.transpose()  
    //.map{it -> [it[1].getSimpleName(),it[1]]}
    .set { read_pairs_fastq }
    //.println()


read_pairs_fastq.println()

Channel
    .fromPath(params.reads)
    .map {[file(it).getSimpleName(), file(it)]}
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set {read_single_fastq}
    

read_single_fastq.println()
*/
//////////////////////////version with sampleConfig file

//!!!!!!!!!!!!!!! check collectFile in nextflow. try to group them according to the family ID

Channel
    .fromPath(params.samplesConfig)
    .splitCsv(header:true)
    .map{ row-> tuple(row.order_id, tuple(file(row.read1), file(row.read2))) }
    .into { read_pairs_fastq; samples_for_single }


Channel
    .fromPath(params.samplesConfig)
    .splitCsv(header:true)
    .map{ row-> tuple(row.order_id, row.pipeline_run_id) }
    .set { ids_for_vcf_upload }


samples_for_single.transpose().map{it -> [it[1].getSimpleName(),it[1]]}.set{ read_single_fastq }

//read_single_fastq.println()
//ids_for_vcf_upload.println()





/*in the last step we are creating vcf file with all variants. 
To upload the variants we will use python script vcf_to_csv... . 
There we will put 2 parameters: file path from channel in the last rule (vep_annotated_vcf channel)
and pipeline_run_id from ids_for_vcf_upload channel created above. 
We need to join them (https://www.nextflow.io/docs/latest/operator.html#combining-operators --> use join)
according to sample_id from vep_annotated_vcf and order_id from ids_for_vcf_upload,
so that the correct pipeline_run_id will be lineked with correct file. 
Otherwise wrong variants will be uploaded to the patients (randomly mixed).
There's also orderID parameter necessary for vcf_upload, especially if we would have trio vcf file.
The bast way wolud be to pass it as parameter to python script as well.
Do we need gender info in sample config file? maybe to check if it's correct...
*/



process rawFastQC {
  tag "$raw_fastQC"
  

  //publishDir "***/results_v01/fastQC", mode: 'move'

  input: 
  tuple val(sample_id), file(read) from read_single_fastq

  output:
  file "${sample_id}_fastqc*" into raw_fastQC_report
  
  """
  sudo $params.fastqc $read -t $params.lowCpus
  """
}



process trimmomatic {
  tag "$trimmomatic"
  //storeDir "***/results_v01/trimmomatic"

  input:
  tuple val(sample_id), file(reads) from read_pairs_fastq

  output:
  tuple val(sample_id), file('*_trimmomatic_R?.fastq') into trimmed_paired_fastq, trimmed_fastq_for_QC_merged
  //tuple val(sample_id), path('*_trimmomatic_?U') into trimmed_unpaired_fastq
  tuple val(sample_id), file('*_stats') into stats_trimmomatic

  """
  java -jar $params.trimmomatic PE \
  -phred33 -threads 1 $reads \
  -baseout ${sample_id}_trimmomatic \
  -summary ${sample_id}_stats \
  ILLUMINACLIP:***/trimmomatic/Trimmomatic-0.38/adapters/TruSeq2-PE.fa:2:30:10\
   LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:32
  rm ${sample_id}_trimmomatic_?U 
  mv ${sample_id}_trimmomatic_1P ${sample_id}_trimmomatic_R1.fastq
  mv ${sample_id}_trimmomatic_2P ${sample_id}_trimmomatic_R2.fastq
  """
  //actually it's not necessary to remove trimmomatic_?U fiels, because work folder 
  //should be remove with all intermediate files at the end of the process
  //it's also not necessary to rename trimmomatic_?P files, because 
  //they are defined as channel, which is an input for thee next process. Name does't matter
}


trimmed_fastq_for_QC_merged
.transpose()
.map{it -> [it[1].getSimpleName(),it[1]]}
.set{trimmed_fastq_for_QC_final}
//Above transpose is necessary for fastqc, because fastqc is seperate for R1 and R2 as 2 processes
//Aligment is together for R1 and R2 as one process
//For Aligment: [father, [***/results_v01/work/5d/125de08931ff094528a23813056e86/father_trimmomatic_R1.fastq, ***/results_v01/work/5d/125de08931ff094528a23813056e86/father_trimmomatic_R2.fastq]]
//For fastqc: [father_trimmomatic_R1, ***/results_v01/work/5d/125de08931ff094528a23813056e86/father_trimmomatic_R1.fastq],[father_trimmomatic_R2, ***/results_v01/work/5d/125de08931ff094528a23813056e86/father_trimmomatic_R2.fastq]


process trimmedFastQC {
  tag "$trimmed_fastQC"

  //storeDir "***/results_v01/trimmomatic_fastQC"

  input:
  tuple val(sample_id), file(read) from trimmed_fastq_for_QC_final

  output:
  file "${sample_id}_fastqc*" into trimmed_fastQC_report

  """
  $params.fastqc $read -t $params.lowCpus
  """
}



process BWAMEM {
  tag "$BWAMEM"
  //publishDir "***/results_v01/aligner"

  input:
  tuple val(sample_id), file(reads) from trimmed_paired_fastq

  output:
  tuple val(sample_id), file('*.bam') into bam_file 
 


  """
  $params.BWA mem \
  -R "@RG\\tID:${sample_id}\\tPL:ILLUMINA\\tLB:${sample_id}\\tSM:${sample_id}" \
  -t$params.highCpus -M $params.hg38 \
  $reads | samtools sort - -m 100M\
  -o ${sample_id}.bam
  """
}

process MarkDuplicates {
  tag "$mark_duplicates"

  input:
  tuple val(sample_id), file(bam) from bam_file

  output:
  tuple val(sample_id), file('*_dedup.bam') into dedup_bam_for_recal, dedup_bam_for_applyBQSR

  """
  java -Djava.io.tmpdir=temp -jar $params.picard MarkDuplicates \
  I=$bam O=${sample_id}_dedup.bam REMOVE_DUPLICATES=true \
  METRICS_FILE=${sample_id}_dedupmetrics.txt 
  """
  //work folder should be removed when all files are generated, so no need of rm not dedup bamfile 
  //in original pipeline old bam file is removed 
}

process base_recalibrator {
  tag "$base_recalibrator"

  input:
  tuple val(sample_id), file(dedup_bam) from dedup_bam_for_recal

  output:
  tuple val(sample_id), file('*.recal_data.table') into base_recal_table

  """
  $params.gatk4 BaseRecalibrator \
  -R $params.hg38 \
  -I $dedup_bam \
  -O ${sample_id}.recal_data.table \
  --known-sites $params.dbsnp38
  """
}

//for applyBQSR script we need to use recal_data.table file and dedup_bam file.
// They come from diferent channels. In order to be sure that these both files are from
// the same patient, we need to join them to one tuple channel:
dedup_bam_for_applyBQSR.join(base_recal_table).set{dedupBam_recalTable_merged}


process apply_BQSR {
  tag "$apply_BQSR"
  //storeDir "***/results_v01/aligner"

  input:
  tuple val(sample_id), file(dedup_bam), file(recal_table) from dedupBam_recalTable_merged
  
  output:
  tuple val(sample_id), file('*.bqsr.bam') into bqsr_recal_bam_wgs_metr, bqsr_recal_bam_flagstat
  //tuple val(sample_id), file('*.bqsr.bai') into bqsr_recal_bai // -> not used, can be removed??
  tuple val(sample_id), file('*.bqsr.bam'), file('*.bqsr.bai') into mosdepth_bam_bai, haplotypeCall_bam_bai
  
  """
  $params.gatk4 ApplyBQSR \
  -R $params.hg38 \
  -I $dedup_bam \
  --bqsr-recal-file $recal_table \
  -O ${sample_id}.bqsr.bam
  """
  //there are 2 mv commands for this process in original pipeline.
  //It's removing bqsr from file names and replacing original bam files with recalibarted ones.
  //This process can be skipped, because file name doesn't have to be specified in next process.
  //Files are as output here and will be as input in next processes. 
}





process wgs_metrics {
  tag "$wgs_metrics"

  input:
  tuple val(sample_id), file(recal_bam) from bqsr_recal_bam_wgs_metr

  output:
  tuple val(sample_id), file('*_wgs_metrics.txt') into wgs_metrics_txt

  """
  java -jar $params.picard CollectWgsMetrics \
  I=$recal_bam \
  O=${sample_id}_wgs_metrics.txt \
  R=$params.hg38 \
  INTERVALS=${params.bedSureSelect}.interval
  """
}

process metrics_to_tsv {
  tag "$metrics_to_tsv"

  input:
  tuple val(sample_id), file(metrics_txt) from wgs_metrics_txt

  output:
  tuple val(sample_id), file('*_wgs_metrics.tsv') into wgs_metrics_tsv

  """
  grep -A 2 "## METR" $metrics_txt | grep -v "##" > ${sample_id}_wgs_metrics.tsv
  """
}

process exon_coverage {
  tag "$exon_coverage"

  input:
  tuple val(sample_id), file(recal_bam), file(recal_bai) from mosdepth_bam_bai
  //although bai file is not as input to command, mosdepth is looking for it in a folder
  //where bam file is located. Becuase bai is provided as recal_bai in input, symlink is created for 
  //bai as well as for recal_bam in exon_coverage execution folder. 

  output:
  tuple val(sample_id), file('*.mosdepth.global.dist.txt'), file('*.mosdepth.region.dist.txt') into mosdepth_txt
  tuple val(sample_id), file('*.regions.bed.gz'), file('*.thresholds.bed.gz') into mosdepth_bed

  """
  $params.mosdepth \
  --by $params.bedSureSelect \
  --thresholds 1,5,10,20,30,70,100 \
  $sample_id \
  $recal_bam -n
  """
  //In the original pipeline are *csi files removed and .gz files are unzipped.
  //work folder should be removed at the end of the pipeline so no need to remove now *csi
  //also we should check the reason why gz files are unzipped.
  //If they are uploadeed to database we should check why we need this data there
}

process samtools_flagstat {
  tag "$samtools_flagstat"
  
  input:
  tuple val(sample_id), file(recal_bam) from bqsr_recal_bam_flagstat

  output:
  tuple val(sample_id), file('*.flagstat') into stats_flagstat



  """
  samtools flagstat $recal_bam \
  > ${sample_id}.flagstat
  """
}

process haplotype_calling {
  tag "$haplotype_calling"
  input:
  tuple val(sample_id), file(recal_bam), file(recal_bai) from haplotypeCall_bam_bai
  //check comment in input, in exon_coverage process

  output:
  tuple val(sample_id), file('*.vcf') into vcf_raw_for_snp, vcf_raw_for_indel


  """
  $params.gatk4 HaplotypeCaller \
  --native-pair-hmm-threads $params.highCpus \
  -R $params.hg38 \
  -I $recal_bam \
  -O ${sample_id}.vcf
  """
}

process select_variants_snp {
  tag "$select_variants"

  input:
  tuple val(sample_id), file(vcf) from vcf_raw_for_snp

  output:
  tuple val(sample_id), file("*_snp_raw.vcf") into vcf_selected

  """
  $params.gatk4 SelectVariants \
  -R $params.hg38 \
  -V $vcf \
  -select-type SNP \
  -O ${sample_id}_snp_raw.vcf
  """
}

process filter_variants_snp {
  tag "$filter_variants_snp"

  input:
  tuple val(sample_id), file(vcf) from vcf_selected

  output:
  tuple val(sample_id), file("*_snp_filtered.vcf") into vcf_filtered_for_recal, vcf_filtered_for_VQSR

  """
  $params.gatk4 VariantFiltration \
  -R $params.hg38 \
  -V $vcf \
  --filter-expression "QD < 2.0" --filter-name "filter_QD" \
  --filter-expression "FS > 60.0" --filter-name "filter_FS" \
  --filter-expression "MQ < 40.0" --filter-name "filter_MQ" \
  --filter-expression "MQRankSum < -12.5" --filter-name "filter_MQRank" \
  --filter-expression "ReadPosRankSum < -8.0" --filter-name "filter_ReadPosRankSum" \
  -O ${sample_id}_snp_filtered.vcf
  """
}
//here we should add also --max-gaussians 4. In original pipeline it was also changed
process variant_recalibrator_snp {
  tag "$variant_recalibrator_snp"

  input:
  tuple val(sample_id), file(vcf) from vcf_filtered_for_recal

  output:
  tuple val(sample_id), file("*_snp_filtered.recal"), file("*_snp_filtered.tranches"), file("*_snp_filtered.recal.idx") into recal_tranches

  """
  $params.gatk4 VariantRecalibrator \
  -R $params.hg38 \
  -V $vcf \
  -O ${sample_id}_snp_filtered.recal \
  -tranches-file ${sample_id}_snp_filtered.tranches \
  --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $params.hapmap \
  --resource:omni,known=false,training=true,truth=true,prior=12.0 $params.omni \
  --resource:1000G,known=false,training=true,truth=false,prior=10.0 $params.g1000G \
  --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $params.dbsnp38 \
  -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP
  """
}

vcf_filtered_for_VQSR.join(recal_tranches).set{vcf_recal_tranches_for_VQSR}


process apply_vqsr_snp {
  tag "$apply_vqsr_snp"

  input:
  tuple val(sample_id), file(vcf), file(recal), file(tranches), file(recal_idx) from vcf_recal_tranches_for_VQSR

  output:
  tuple val(sample_id), file("*_snp_filtered_VQSR.vcf") into vqsr_vcf_snp

  """
  $params.gatk4 ApplyVQSR \
  -R $params.hg38 \
  -V $vcf \
  -tranches-file $tranches \
  -recal-file $recal \
  -O ${sample_id}_snp_filtered_VQSR.vcf \
  -ts-filter-level 99.5 -mode SNP
  """
}

/*
process bgzip_vcf_snp {
  tag "$bgzip_vcf_snp"

  input:
  tuple val(sample_id), file(vcf) from vqsr_vcf_snp

  output:
  tuple val(sample_id), file("*_snp_filtered_VQSR.vcf.gz") into bgzipped_snp_vcf

  """
  bgzip -c -@ 1 $vcf > ${sample_id}_snp_filtered_VQSR.vcf.gz
  """
}
*/
process tabix_vcf_snp {
  tag "$tabix_vcf_snp"

  input:
  tuple val(sample_id), file(vcf) from vqsr_vcf_snp

  output:
  tuple val(sample_id), file("*_snp_filtered_VQSR.vcf.gz"), file("*_snp_filtered_VQSR.vcf.gz.tbi") into tabix_snp_vcf

  """
  bgzip -c -@ 1 $vcf > ${sample_id}_snp_filtered_VQSR.vcf.gz
  tabix -p vcf ${sample_id}_snp_filtered_VQSR.vcf.gz

  """
}

process select_variants_indel {
  tag "$select_variants_snp"

  input:
  tuple val(sample_id), file(vcf) from vcf_raw_for_indel

  output:
  tuple val(sample_id), file("*_indel_raw.vcf") into vcf_selected_indel
  """
  $params.gatk4 SelectVariants \
  -R $params.hg38 \
  -V $vcf \
  -select-type INDEL \
  -O ${sample_id}_indel_raw.vcf
  """
}

process filter_variants_indel {
  tag "$filter_variants_indel"

  input:
  tuple val(sample_id), file(vcf) from vcf_selected_indel

  output:
  tuple val(sample_id), file("*_indel_filtered.vcf") into vcf_filtered_for_recal_indel, vcf_filtered_for_VQSR_indel

  """
  $params.gatk4 VariantFiltration \
  -R $params.hg38 \
  -V $vcf \
  --filter-expression "QD < 2.0" --filter-name "filter" \
  -O ${sample_id}_indel_filtered.vcf
  """
}

process variant_recalibrator_indel {
  tag "$variant_recalibrator_indel"

  input:
  tuple val(sample_id), file(vcf) from vcf_filtered_for_recal_indel

  output:
  tuple val(sample_id), file("*_indel_filtered.recal"), file("*_indel_filtered.tranches"), file("*_indel_filtered.recal.idx") into recal_tranches_indel

  """
  $params.gatk4 VariantRecalibrator \
  -R $params.hg38 \
  -V $vcf \
  -O ${sample_id}_indel_filtered.recal \
  -tranches-file ${sample_id}_indel_filtered.tranches \
  --max-gaussians 4 \
  --resource:mills,known=false,training=true,truth=true,prior=12.0 $params.mills \
  --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $params.dbsnp38 \
  -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -mode INDEL
  """
}

vcf_filtered_for_VQSR_indel.join(recal_tranches_indel).set{vcf_recal_tranches_for_VQSR_indel}


process apply_vqsr_indel {
  tag "$apply_vqsr_indel"

  input:
  tuple val(sample_id), file(vcf), file(recal), file(tranches), file(recal_idx) from vcf_recal_tranches_for_VQSR_indel

  output:
  tuple val(sample_id), file("*_indel_filtered_VQSR.vcf") into vqsr_vcf_indel

  """
  $params.gatk4 ApplyVQSR \
  -R $params.hg38 \
  -V $vcf \
  -tranches-file $tranches \
  -recal-file $recal \
  -O ${sample_id}_indel_filtered_VQSR.vcf  \
  -ts-filter-level 99.0 -mode INDEL  
  """
}

process tabix_vcf_indel {
  tag "$tabix_vcf_indel"

  input:
  tuple val(sample_id), file(vcf) from vqsr_vcf_indel

  output:
  tuple val(sample_id), file("*_indel_filtered_VQSR.vcf.gz"), file("*_indel_filtered_VQSR.vcf.gz.tbi") into tabix_indel_vcf

  """
  bgzip -c -@ 1 $vcf > ${sample_id}_indel_filtered_VQSR.vcf.gz
  tabix -p vcf ${sample_id}_indel_filtered_VQSR.vcf.gz

  """
}

tabix_snp_vcf.join(tabix_indel_vcf).set{vcfs_to_concat}


process vcf_concat {
  tag "$vcf_concat"
  echo true

  input:
  tuple val(sample_id), file(vcf_snp), file(snp_tbi), file(vcf_indel), file(indel_tbi) from vcfs_to_concat


  output:
  tuple val(sample_id), file("*_filtered.vcf.gz") into vcf_to_annotate

  """
  bcftools concat $vcf_snp $vcf_indel -a -o ${sample_id}_filtered.vcf.gz -Oz
  tabix -p vcf ${sample_id}_filtered.vcf.gz
  """
}

process snpeff_annotate {
  tag "$snpeff_annotate"

  input:
  tuple val(sample_id), file(vcf) from vcf_to_annotate

  output:
  tuple val(sample_id), file("*_SnpEff_HGVS.vcf") into snpeff_annotated_vcf

  """
  java -jar ${params.snpeff}snpEff.jar \
  -noLog -c ${params.snpeff}snpEff.config \
  -v -noStats -no-downstream -no-intergenic -no-upstream \
  -spliceRegionIntronMax 20 hg38 \
  $vcf > ${sample_id}_SnpEff_HGVS.vcf
  """
}

//Channel
//    .fromPath("***/dbnsfp/dbNSFP4.0a.gz.tbi")    
//    .set {dbNSFP4_index}




process dbnsfp_annotate {
  tag "$dbnsfp_annotate"

  input:
  tuple val(sample_id), file(vcf) from snpeff_annotated_vcf
  //file dbnsfp_idx from dbNSFP4_index

  output:
  tuple val(sample_id), file("*_Final.vcf") into dbnsfp_annotated_vcf

  """
  java -jar ${params.snpeff}SnpSift.jar dbnsfp -v -db $params.dbnsfp -f hg19_pos\\(1-based\\),MetaSVM_pred,MetaLR_pred,M-CAP_pred,FATHMM_pred,VEST4_score,Aloft_pred,rs_dbSNP151,SIFT_pred,MutationTaster_pred,Polyphen2_HVAR_pred,Polyphen2_HDIV_pred,CADD_phred,REVEL_score,REVEL_rankscore,1000Gp3_AC,1000Gp3_AF,1000Gp3_AFR_AC,1000Gp3_AFR_AF,1000Gp3_EUR_AC,1000Gp3_EUR_AF,1000Gp3_AMR_AC,1000Gp3_AMR_AF,1000Gp3_EAS_AC,1000Gp3_EAS_AF,1000Gp3_SAS_AC,1000Gp3_SAS_AF,TWINSUK_AC,TWINSUK_AF,ALSPAC_AC,ALSPAC_AF,ESP6500_AA_AC,ESP6500_AA_AF,ESP6500_EA_AC,ESP6500_EA_AF,ExAC_Adj_AC,ExAC_Adj_AF,ExAC_AFR_AC,ExAC_AFR_AF,ExAC_AMR_AC,ExAC_AMR_AF,ExAC_EAS_AC,ExAC_EAS_AF,ExAC_FIN_AC,ExAC_FIN_AF,ExAC_NFE_AC,ExAC_NFE_AF,ExAC_SAS_AC,ExAC_SAS_AF,gnomAD_exomes_AC,gnomAD_exomes_AF,gnomAD_exomes_AFR_AC,gnomAD_exomes_AFR_AF,gnomAD_exomes_AMR_AC,gnomAD_exomes_AMR_AF,gnomAD_exomes_ASJ_AC,gnomAD_exomes_ASJ_AF,gnomAD_exomes_EAS_AC,gnomAD_exomes_EAS_AF,gnomAD_exomes_FIN_AC,gnomAD_exomes_FIN_AF,gnomAD_exomes_NFE_AC,gnomAD_exomes_NFE_AF,gnomAD_exomes_SAS_AC,gnomAD_exomes_SAS_AF,clinvar_clnsig,clinvar_var_source,clinvar_review,Ensembl_geneid,Ensembl_transcriptid,Ensembl_proteinid -g hg38 $vcf > ${sample_id}_Final.vcf
  """
}

process vep_annotate {
  tag "$vep_annotate"

  input:
  tuple val(sample_id), file(vcf) from dbnsfp_annotated_vcf

  output:
  tuple val(sample_id), file("*_Final_vep.vcf") into vep_annotated_vcf

  """
  $params.vep \
  -i $vcf \
  -o ${sample_id}_Final_vep.vcf \
  --cache --dir_cache ***.vep \
  --dir ***.vep \
  --dir_plugins ***.vep \
  --force_overwrite --vcf \
  --fields "Allele,Consequence,Feature_type,Feature,IMPACT,SYMBOL" \
  --fasta $params.hg38VEP \
  --offline --fork 30
  """
  //in vep?annotate process I used I used reference genome from different location, 
  //which is a copy of the one used in other processes. It has something to do 
  //with indexing genome by vep. It should be investigated and the same genome should
  //be used for all processes. Maybe deleting indexed genome form nico 
  //and generating it new would work.  
}

process quality_filter_vcf {
  tag "$quality_filter_vcf"

  input:
  tuple val(sample_id), file(vcf) from vep_annotated_vcf

  output:
  tuple val(sample_id), file("*.recode.vcf") into qual_filtered_vcf, qual_filtered_vcf_for_clinsig

  """
  vcftools \
  --vcf $vcf \
  --out ${sample_id} \
  --minGQ 15 --minDP 10 --max-missing 1 \
  --remove-filtered-all \
  --recode --recode-INFO-all 
  """
}

process bed_filter_vcf {
  tag "$bed_filter_vcf"

  input:
  tuple val(sample_id), file(vcf) from qual_filtered_vcf

  output:
  tuple val(sample_id), file("*_bed_filtered.vcf") into bed_filtered_vcf

  """
  ${params.bedtools} intersect \
  -a $vcf \
  -b ${params.bedSureSelect25bpExtend} \
  -wa -header > ${sample_id}_bed_filtered.vcf 
  """
}

qual_filtered_vcf_for_clinsig.join(bed_filtered_vcf).into{vcf_for_append_clinsig;test}

test.println()


process append_clin_sig {
  tag "$append_clin_sig"

  input:
  tuple val(sample_id), file(vcf_qual_fitered), file(vcf_bed_fitered) from vcf_for_append_clinsig

  output:
  tuple val(sample_id), file("*_filtered_with_ClinSig.vcf") into filtered_with_ClinSig_vcf

  """
  grep -v '^##' $vcf_qual_fitered | grep -E 'dbNSFP_clinvar_clnsig|dbNSFP_rs_dbSNP151' >> $vcf_bed_fitered
  bcftools sort $vcf_bed_fitered -o ${sample_id}_sorted_with_dup.vcf
  awk '! a[\$0]++' ${sample_id}_sorted_with_dup.vcf > ${sample_id}_filtered_with_ClinSig.vcf
  """
  //In above commands we need single quots, because of dolar sign in awk
  //awk in this process removes dulplicated variants after appending variants wih ClinVar and rs_dbSNP entry
}

filtered_with_ClinSig_vcf.join(ids_for_vcf_upload).set{vcf_with_pipeline_id}

process vcf_to_db_upload {
  tag "$vcf_to_db_upload"

  maxForks 1

  input:
  tuple val(sample_id), file(vcf), val(pipeline_run_id) from vcf_with_pipeline_id

  """
  python3 ${params.vcfTodb} -o '$sample_id' -v '$vcf' -i '$pipeline_run_id'
  """
}



 

////////////////   !!!!!!!!!!!UPLOAD PROCESS ATTENTION!!!!!!!!!!!!!!!!    //////////////
/* For upload process run upload for one sample after the other, 
not all at once. Due to present "IF EXISTS" in sql upload script, running upload for all
samples in parallel might cause some errors. You can add "maxForks" in the process which shouldn't be parallelized







  /*
  """
  java -jar /${params.snpeff}/SnpSift.jar \
  dbnsfp -v -db ***/dbNSFP4.0a.gz \
  -f hg19_pos\\(1-based\\),MetaSVM_pred,MetaLR_pred,M-CAP_pred,FATHMM_pred,\
  VEST4_score,Aloft_pred,rs_dbSNP151,SIFT_pred,MutationTaster_pred,\
  Polyphen2_HVAR_pred,Polyphen2_HDIV_pred,CADD_phred,REVEL_score,REVEL_rankscore,\
  1000Gp3_AC,1000Gp3_AF,1000Gp3_AFR_AC,1000Gp3_AFR_AF,1000Gp3_EUR_AC,\
  1000Gp3_EUR_AF,1000Gp3_AMR_AC,1000Gp3_AMR_AF,1000Gp3_EAS_AC,1000Gp3_EAS_AF,\
  1000Gp3_SAS_AC,1000Gp3_SAS_AF,TWINSUK_AC,TWINSUK_AF,ALSPAC_AC,ALSPAC_AF,\
  ESP6500_AA_AC,ESP6500_AA_AF,ESP6500_EA_AC,ESP6500_EA_AF,ExAC_Adj_AC,ExAC_Adj_AF,\
  ExAC_AFR_AC,ExAC_AFR_AF,ExAC_AMR_AC,ExAC_AMR_AF,ExAC_EAS_AC,ExAC_EAS_AF,\
  ExAC_FIN_AC,ExAC_FIN_AF,ExAC_NFE_AC,ExAC_NFE_AF,ExAC_SAS_AC,ExAC_SAS_AF,\
  gnomAD_exomes_AC,gnomAD_exomes_AF,gnomAD_exomes_AFR_AC,gnomAD_exomes_AFR_AF,\
  gnomAD_exomes_AMR_AC,gnomAD_exomes_AMR_AF,gnomAD_exomes_ASJ_AC,gnomAD_exomes_ASJ_AF,\
  gnomAD_exomes_EAS_AC,gnomAD_exomes_EAS_AF,gnomAD_exomes_FIN_AC,gnomAD_exomes_FIN_AF,\
  gnomAD_exomes_NFE_AC,gnomAD_exomes_NFE_AF,gnomAD_exomes_SAS_AC,gnomAD_exomes_SAS_AF,\
  clinvar_clnsig,clinvar_var_source,clinvar_review,Ensembl_geneid,Ensembl_transcriptid,\
  Ensembl_proteinid -g hg38 \
  $vcf > ${sample_id}_Final.vcf
  """
  */











//trimmed_unpaired_fastq.println()
// ***/fastq/index_R1.fastq ***/fastq/index_R2.fastq


/* 

*************test process with channel form fake txt files**************

Channel
    .fromPath("***/fake_files/*.txt")
    .map {[file(it).getSimpleName(), file(it)]}
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set {test_files}

dedup_bam_file_BQSR.join(test_files).set{test_joined}

//test_files.println()
//dedup_bam_file_BQSR.println()

process test{
  echo true
  input:
  tuple val(sample_id), file(dedup_bam), file(text) from test_joined
  //tuple val(sample_id), file(bam) from dedup_bam_file_BQSR
  //tuple val(sample_id_2), file(test) from test_files

  
  //echo true is necessary for printf working
  """
  printf $sample_id
  printf "\n"
  printf $dedup_bam
  printf "\n"
  printf $text
  printf "\n"
  printf "________"
  """
}

*/