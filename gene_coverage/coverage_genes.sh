#!/bin/bash


can_exons_ref="***gencode.v32.annotation_only_canonical_exons_subtracted_UTR_gen_name.bed"

echo "Enter sample ID:"
read sample_id
#S30633A
bam_path="***/${sample_id}/results/bwa/index_dedup.bam"

#checking if bam file for selected sample exists, if not - aborting.
if [[ ! -f ${bam_path} ]] ; then
    echo "File ${bam_path} does not exist, aborting. Please check sample ID"
    exit
fi

echo -e "Input file: ${bam_path}\n"

#creating output directory. It consist sample ID and time when the script was executed
time_now=$(date '+%Y-%m-%d_%H-%M-%S')
o_path="***coverage_output/${sample_id}/${time_now}"

#creating directories with permissions for everyone
mkdir -p -m777 "***coverage_output/${sample_id}"
mkdir -p -m777 ${o_path}

echo -e "Do you want to generate coverage for all available genes? (y/n)"
read all_genes_q

if [[ $all_genes_q == "y" ]]
then
  echo "Including all genes...\n"
  gene_input=$(awk 'BEGIN { ORS = "," } { print }' ***gencode.v32.only_genes_names.txt)
  echo ${gene_input}
  
else
  echo "Enter coma seperated list of genes:"
  read gene_input
fi

#saving input genes to file
echo ${gene_input} > ${o_path}/genes_row.txt
#changing genes in a row to gens in column
awk -F, -v OFS="\n" '{$1=$1; print}' ${o_path}/genes_row.txt \
| awk '{gsub(/ /, "",$1); print$1}' > ${o_path}/genes_col.txt 
#selecting regions (exons) for entered genes and creating bed file  
awk 'NR==FNR{a[$1];next} ($5) in a' ${o_path}/genes_col.txt ${can_exons_ref} \
| awk '{gsub(/^chr/,""); print}' \
| awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$5":"$6}' > ${o_path}/regions.bed

echo -e "Regions of interest has been created.\n"

echo -e "Generating coverage...\n"

#runing coverage script
mosdepth --by ${o_path}/regions.bed --thresholds 1,5,10,20,30,70,100 ${o_path}/${sample_id} ${bam_path} -n
#unzip output
gunzip ${o_path}/*.gz
#in sample_id.thresholds.bed is number of bases covered 1x, 5x 10x ... per exon
#below command calculates percentage of covered bases per gene
awk -F '[\t:]' 'FNR==1 {next} BEGIN{OFS="\t"} BEGIN{print "#gene" "\t" "total" "\t" "1X" "\t" "5X" "\t" "10X" "\t" "20X" "\t" "30X" "\t" "70X" "\t" "100X"}{a[$4] += ($3-$2)}{b[$4] += $8}{c[$4] += $9}{d[$4] += $10}{e[$4] += $11}{f[$4] += $12}{g[$4] += $13}{h[$4] += $14} END{for (i in a) print i, a[i], b[i]/a[i], c[i]/a[i], d[i]/a[i], e[i]/a[i], f[i]/a[i], g[i]/a[i], h[i]/a[i]}' ${o_path}/${sample_id}.thresholds.bed > ${o_path}/${sample_id}.coverage_per_gene.tsv

echo -e "Genes with not comletely covered coding regions: \n"
#extracting genes where are regions without even one read
awk 'FNR==1 {next} BEGIN{OFS="\t"} BEGIN{print "#gene" "\t" "total" "\t" "1X" "\t" "5X" "\t" "10X" "\t" "20X" "\t" "30X" "\t" "70X" "\t" "100X"} $3<1 {print}' ${o_path}/${sample_id}.coverage_per_gene.tsv | tee ${o_path}/${sample_id}.not_covered_genes.tsv | cat
#below command calculates percentage of at least 1x covered bases per exon
awk 'FNR==1 {next} BEGIN{OFS="\t"} BEGIN{print "#gene" "\t" "total" "\t" "1X" "\t" "5X" "\t" "1X_%"} {$12=$3-$2} {if ($5!=0) $13=$5/$12; else $13=0} {print $4,$12,$5,$6,$13}' /${o_path}/${sample_id}.thresholds.bed > ${o_path}/${sample_id}.coverage_per_exon.tsv

echo -e "\nExons with not comletely covered coding regions: \n"
#extracting genes where are regions without even one read
awk 'FNR==1 {next} BEGIN{OFS="\t"} BEGIN{print "#gene" "\t" "total" "\t" "1X" "\t" "5X" "\t" "1X_%"} $5<1 {print}' ${o_path}/${sample_id}.coverage_per_exon.tsv | tee ${o_path}/${sample_id}.not_covered_exons.tsv | cat

echo -e "\nGenes not inclueded in the analysis (not found in ref file):"

awk 'NR==FNR {exclude[$1];next} !($1 in exclude)' ${o_path}/${sample_id}.coverage_per_gene.tsv ${o_path}/genes_col.txt | tee ${o_path}/not_included_genes.txt | cat

echo -e "\nDone. Data has been saved in ${o_path}.\n"

#ACAN, ACP5, AGPS, ALPL, ANKH, ANO5, ARHGAP31, ARSE, ATP6V0A2, B3GALT6, B3GAT3, B4GALT7, BMP1, BMP2, BMPR1B, CA2, CANT1, CASR, CC2D2A, CDH3, CDKN1C, CEP290, CHST14, CHST3, CHSY1, CLCN5, CLCN7, COL10A1, COL11A1, COL11A2, COL1A1, COL1A2, COL2A1, COL9A1, COL9A2, COL9A3, COMP, CRTAP, CTSK, CUL7, DDR2, DHCR24, DLL3, DLX3, DMP1, DYM, DYNC2H1, EBP, EIF2AK3, ENPP1, ESCO2, EVC, EVC2, EXT1, EXT2, FAM20C, FBLN1, FBN1, FBXW4, FERMT3, FGF10, FGF23, FGFR1, FGFR2, FGFR3, FKBP10, FLNA, FLNB, FMN1, GALNT3, GDF5, GLI3, GNAS, GNPAT, GORAB, GPC6, GREM1, HDAC4, HES7, HOXD13, HPGD, HSPG2, ICK, IDH1, IDH2, IFITM5, IFT122, IFT140, IFT80, IHH, KIF22, KIF7, LEMD3, LFNG, LIFR, LMBR1, LMNA, LRP4, LRP5, MAFB, MATN3, MESP2, MGP, MKS1, MMP13, MMP2, MMP9, MYCN, NEK1, NIPBL, NKX3-2, NOG, NOTCH2, NPR2, OBSL1, OSTM1, P3H1, PAPSS2, PCNT, PDE4D, PEX7, PHEX, PIGV, PITX1, PLOD2, PPIB, PPP3CA, PRKAR1A, PTH1R, PTHLH, PTPN11, PYCR1, RECQL4, ROR2, RPGRIP1L, RUNX2, SALL1, SALL4, SERPINF1, SERPINH1, SH3PXD2B, SHH, SHOX, SLC26A2, SLC34A3, SLC35D1, SLC39A13, SMAD4, SMARCAL1, SOST, SOX9, SP7, SULF1, TBCE, TBX15, TBX3, TBX5, TBX6, TBXAS1, TCIRG1, TGFB1, THPO, TMEM216, TMEM67, TNFRSF11A, TNFRSF11B, TNFSF11, TP63, TREM2, TRIP11, TRPS1, TRPV4, TYROBP, WDR35, WISP3, WNT3, WNT5A, WNT7A, ZMPSTE24 
