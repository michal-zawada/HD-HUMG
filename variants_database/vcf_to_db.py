import vcf
import os
import pandas as pd
import psycopg2
import pandas.io.sql as sqlio
import datetime, time
import subprocess
import csv
import sys
import argparse



#for nextflow remove orderID='index' in line ~331


ipath=os.path.join('***','WES')
opath=os.path.join('***','vcfs_with_variant_id')

parser=argparse.ArgumentParser()
parser.add_argument('-o', '--orders', action='store', dest='orderID', required=False, 
                    help="use --orders parameter with comma-seperated order_ids in quots, for example: --orders 'order_id1,order_id2'")
parser.add_argument('-v', '--vcf', action='store', dest='vcfFile', required=False)
parser.add_argument('-i', '--pipelineRunId', action='store', dest='pipelineRunID', required=False)


args=parser.parse_args()
sys.stdout.write('\nInput order(s): {}\n'.format(args.orderID))

#vcf_file="10k_S30633A_Final_vep.vcf"
#vcf_file="Final.vcf"
#vcf_file="sorted_with_ClinVar_25.vcf"
#orderID='S30633A'
""" < ---------------------------------------UNCOMMENT THIS BLOCK
orderID=args.orderID
i_vcf_path=args.vcfFile
#runID='190919_ST-K00352_0311_AHCKCCBBXY' #actually this we don't need with new db set up where we need pipelineRunID  
pipelineRunID=7 #S30633A has pipeline_run_id =7 in db #we could extract pipelineRunID also from db base on runID and orderID, because it's already there, but 
                #if we run pipeline twice for thje same sample from the same run, there might be for some time 2 pipelineRunIDs
pipelineRunID=int(args.pipelineRunID)
"""
updateDate='2020-11-20'

#i_vcf_path=(os.path.join(ipath,vcf_file))


orderID='index'
vcfVarID_path=os.path.join(opath,orderID)


pipelineRunID=14
i_vcf_path='***/test_2variants.vcf'
#populations table has to be one to many with variants if we have now update_date



#in vcf file there's empty column called "ID". 
#We are using it and add variant ID from our database if this variant exists 
#or add new variant_id if the variant doesn't exist.
#basically the way is first uploading variants table 
#and then getting varint_id from sql db

def dfFromSql(sqlString,db):
    """as arguments to this function put: 
    1) connection sql comand as string
    2)Dictionary with database info, for example:
        DATABASES = {'production':{'NAME': 'db_name','USER': 'user_name','PASSWORD': '**password**',
        'HOST': 'localhost','PORT': '',},}
        and then as db use: db = DATABASES['production'] 
    The function returns pandas df containing sql request"""
    try:
        connection = psycopg2.connect(user = db['USER'],
                                  password = db['PASSWORD'],
                                  host = db['HOST'],
                                  database = db['NAME'])


        cursor = connection.cursor()
        # Print PostgreSQL Connection properties
        print ( connection.get_dsn_parameters(),"\n")

        # Print PostgreSQL version
        cursor.execute("SELECT version();")
        
        record = cursor.fetchone()

        print("You are connected to - ", record,"\n")

        df=sqlio.read_sql_query(sqlString, connection) 
        #df=pd.read_sql_query(sqlString, connection)    
        return df

    except (Exception, psycopg2.Error) as error :
        print ("Error while connecting to PostgreSQL", error)
    finally:
        #closing database connection.
            if(connection):
                cursor.close()
                connection.close()
                print("PostgreSQL connection is closed")


def vcfTodb(i_vcf_path,db,order): ############ATTENTION! in order to work, it has to be at least one variant in variants table in db. It's due to INSERT ... WHERE NOT EXISTS ... #################


    print("Processing {}".format(i_vcf_path))
    vcf_reader = vcf.Reader(open(i_vcf_path))
    print(vcf_reader.samples)


    """
    f = open("snp_diff_param_than_transcript.txt", "a")
    sample=";".join(x for x in vcf_reader.samples)
    f.write("{}{}\n".format("Variants with dbsnp entry from 10k variants from sample ",sample))
    f.close()
    """

    n=0
    n_missense=0
    #list of populations frequencies from dbnsfp annotator
    dbnsfp="dbNSFP_ALSPAC_AC,dbNSFP_clinvar_var_source,dbNSFP_hg19_pos_1_based_,dbNSFP_ExAC_NFE_AF,dbNSFP_ExAC_SAS_AF,"\
            "dbNSFP_Ensembl_transcriptid,dbNSFP_ExAC_SAS_AC,dbNSFP_1000Gp3_AMR_AF,dbNSFP_ALSPAC_AF,dbNSFP_1000Gp3_AMR_AC,"\
            "dbNSFP_gnomAD_exomes_AMR_AF,dbNSFP_gnomAD_exomes_AMR_AC,dbNSFP_MetaSVM_pred,dbNSFP_1000Gp3_EAS_AC,dbNSFP_gnomAD_exomes_NFE_AF,"\
            "dbNSFP_FATHMM_pred,dbNSFP_ExAC_AFR_AF,dbNSFP_ExAC_AFR_AC,dbNSFP_gnomAD_exomes_AF,dbNSFP_gnomAD_exomes_SAS_AC,dbNSFP_gnomAD_exomes_AC,"\
            "dbNSFP_1000Gp3_EAS_AF,dbNSFP_gnomAD_exomes_SAS_AF,dbNSFP_Aloft_pred,dbNSFP_ExAC_FIN_AC,dbNSFP_ExAC_FIN_AF,dbNSFP_gnomAD_exomes_AFR_AF,"\
            "dbNSFP_rs_dbSNP151,dbNSFP_ESP6500_EA_AC,dbNSFP_clinvar_clnsig,dbNSFP_1000Gp3_AFR_AC,dbNSFP_ExAC_AMR_AF,dbNSFP_1000Gp3_AFR_AF,dbNSFP_MutationTaster_pred,"\
            "dbNSFP_MetaLR_pred,dbNSFP_ExAC_AMR_AC,dbNSFP_ExAC_NFE_AC,dbNSFP_REVEL_rankscore,dbNSFP_1000Gp3_SAS_AC,dbNSFP_ExAC_EAS_AC,dbNSFP_1000Gp3_SAS_AF,"\
            "dbNSFP_ExAC_EAS_AF,dbNSFP_ESP6500_EA_AF,dbNSFP_gnomAD_exomes_FIN_AF,dbNSFP_clinvar_review,dbNSFP_REVEL_score,dbNSFP_TWINSUK_AF,dbNSFP_ExAC_Adj_AC,"\
            "dbNSFP_TWINSUK_AC,dbNSFP_ExAC_Adj_AF,dbNSFP_gnomAD_exomes_AFR_AC,dbNSFP_Ensembl_proteinid,dbNSFP_Ensembl_geneid,dbNSFP_VEST4_score,"\
            "dbNSFP_gnomAD_exomes_NFE_AC,dbNSFP_1000Gp3_AC,dbNSFP_1000Gp3_AF,dbNSFP_gnomAD_exomes_EAS_AF,dbNSFP_gnomAD_exomes_EAS_AC,dbNSFP_CADD_phred,"\
            "dbNSFP_1000Gp3_EUR_AC,dbNSFP_Polyphen2_HDIV_pred,dbNSFP_gnomAD_exomes_ASJ_AC,dbNSFP_1000Gp3_EUR_AF,dbNSFP_gnomAD_exomes_ASJ_AF,dbNSFP_ESP6500_AA_AF,"\
            "dbNSFP_M_CAP_pred,dbNSFP_Polyphen2_HVAR_pred,dbNSFP_SIFT_pred,dbNSFP_ESP6500_AA_AC,dbNSFP_gnomAD_exomes_FIN_AC".split(",")

    qual_params_single=["AN","BaseQRankSum","DP","DS","END","ExcessHet","FS","InbreedingCoeff","MQ","MQRankSum","NEGATIVE_TRAIN_SITE","POSITIVE_TRAIN_SITE","QD","ReadPosRankSum","SOR","VQSLOD","culprit"]
    qual_params_list=["AC","AF","MLEAC","MLEAF"]


    snp_multi_param_list=[]

    se_records=[]
    vep_records =[]
    dbnsnp_freq=[]
    raw_variants=[]
    for record in vcf_reader:
        #record_list=[]
        if n < 50:             
            #print(record.INFO.keys())
            #try:
            #we can use eighter above try block or below if statemant to check if ANN is present for each variant
            #print(record.INFO)

    #creating list of variants for variants sql table
            var_gen={}
            var_gen["chrom"]=record.CHROM
            var_gen["pos"]=record.POS
            var_gen["ref"]=record.REF
            var_gen["alt"]=",".join([str(x) for x in record.ALT])

            raw_variants.append(var_gen)

            n+=1


###################### raw variants upload ######################

############ATTENTION! in order to work, it has to be at least one variant in variants table in db. 
#It's due to INSERT ... WHERE NOT EXISTS ... #################

    try:
        connection = psycopg2.connect(user = db['USER'],
                                    password = db['PASSWORD'],
                                    host = db['HOST'],
                                    database = db['NAME'])
        cursor = connection.cursor()
        # Print PostgreSQL Connection properties
        print ( connection.get_dsn_parameters(),"\n")

        # Print PostgreSQL version
        cursor.execute("SELECT version();")
        
        record = cursor.fetchone()

        print("You are connected to - ", record,"\n")


        
        begin_time_variants = datetime.datetime.now()
        #iterating through list of dictionaries containing raw variants. 
        #Each variant is checked in db if it exists
        for var_gen in raw_variants: 
            cursor.execute("INSERT INTO variants (chrom,pos,ref,alt) SELECT '{0}' AS chrom, {1} AS POS, '{2}' AS ref, '{3}' AS alt FROM variants WHERE NOT EXISTS(SELECT (chrom,pos,ref,alt) FROM variants WHERE chrom='{0}' AND pos={1} AND ref='{2}' AND alt='{3}') LIMIT 1;".format(var_gen["chrom"],var_gen["pos"],var_gen["ref"],var_gen["alt"]))
            #Warning: this is not safe if executed from multiple sessions at the same time (see caveats below). --> from stack overflow - discuss
            #cursor.execute("SELECT variant_id FROM variants WHERE chrom='{}' AND pos={} AND ref='{}' AND alt='{}';".format(var_gen["chrom"],var_gen["pos"],var_gen["ref"],var_gen["alt"]))            

            #record=cursor.fetchone()
            cursor.execute("SELECT variant_id FROM variants WHERE chrom='{}' AND pos={} AND ref='{}' AND alt='{}';".format(var_gen["chrom"],var_gen["pos"],var_gen["ref"],var_gen["alt"])) 


            record=cursor.fetchone()

            #here the variant_id is added to a variant dictionary in a list. Looks like below:
            #[{'chrom': '1', 'pos': 14464, 'ref': 'A', 'alt': 'T', 'variant_id': 44}, {'chrom': '1', 'pos': 14653, 'ref': 'C', 'alt': 'T', 'variant_id': 11},...]   
            record=record[0]            
            var_gen["variant_id"]=record

        #getting sample_id from db for our orderID and runID
        #cursor.execute("""SELECT sample_id FROM samples WHERE order_id='{0}' AND sequencing_run='{1}';""".format(orderID,runID))  
        #record=cursor.fetchone()
        #sampleID= record[0]
        #print(sampleID)



        connection.commit()
        print("Transaction completed successfully ")

    except (Exception, psycopg2.Error) as error :        
        """
        Shall we add below 2 lines of code? Does it work?
        if connnection:
            connection.rollback()
        """
        print ("Error while connecting to PostgreSQL", error)
    finally:
        #closing database connection.
            if(connection):
                cursor.close()
                connection.close()
                print("PostgreSQL connection is closed")

############### Generating vcf with varint_id form db #######################

        #below few lines are saving dicitionary with variants and variant_ids to vcf file and annotating original vcf file with variant_ids. 
    fieldnames=("chrom","pos","variant_id","ref","alt")
    with open('{}.vcf'.format(vcfVarID_path), 'w') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(raw_variants)
    os.system("""sed -i "1s/.*/"#CHROM'\t'POS'\t'ID'\t'REF'\t'ALT"/" {0}.vcf""".format(vcfVarID_path))
    os.system("""(head -n 1 {0}.vcf && tail -n +2 {0}.vcf | sort -k1,1V -k2,2n) > {0}_sorted.vcf""".format(vcfVarID_path))
    os.system("""bgzip -c {0}_sorted.vcf > {0}_sorted.vcf.gz""".format(vcfVarID_path))
    os.system("""tabix -p vcf {0}_sorted.vcf.gz""".format(vcfVarID_path))
    os.system("""bcftools annotate -a {0}_sorted.vcf.gz -c CHROM,POS,ID,REF,ALT {1} > {0}_withID.vcf""".format(vcfVarID_path,i_vcf_path))



    print("raw variants upload: ", datetime.datetime.now() - begin_time_variants)


    #Below 16 rows are to get variants_id from variants table, annotate vcf file with them and then fetch vep annotation info with corresponding variant_id. 
    #It's way faster selecting variant_id from variant table in order to upload vep annotations 
    n=0
    vep_records_with_ID=[]
    se_records_with_ID=[]
    dbnsfp_freq=[]
    detect=[]
    vcf_reader = vcf.Reader(open("{0}_withID.vcf".format(vcfVarID_path)))  
    for record in vcf_reader:
        if n < 50:
            if "CSQ" in record.INFO.keys():                
                for transcript in record.INFO["CSQ"]:
                    transcript_l = transcript.split("|")
                    transcript_l.insert(0,record.ID)
                    transcript_l.insert(7,updateDate)                        
                    vep_records_with_ID.append(tuple(transcript_l))

            if "ANN" in record.INFO.keys():                
                for transcript in record.INFO["ANN"]:                                        
                    transcript_l = transcript.split("|") 
                    if "/" in transcript_l[8]:
                        exons=transcript_l[8].split("/") #number of exon with the variant and total number of exons in this gene are as one value, so we split it
                        transcript_l[8]=exons[0]
                        transcript_l.insert(9,exons[1])
                    else:
                        transcript_l.insert(9,"")                    
                    if "/" in transcript_l[12]:
                        cDNA=transcript_l[12].split("/") #cDNA position of variant and number of the last cDNA in this gene
                        transcript_l[12]=cDNA[0]
                        transcript_l.insert(13,cDNA[1])
                    else:
                        transcript_l.insert(13,"")
                    
                    if "/" in transcript_l[14]:
                        CDS=transcript_l[14].split("/") #CDS position of variant and number of the last CDS in this gene
                        transcript_l[14]=CDS[0]
                        transcript_l.insert(15,CDS[1])
                    else:
                        transcript_l.insert(15,"")

                    if "/" in transcript_l[16]:
                        aminoAcids=transcript_l[16].split("/") #amino acid position of variant and number of the last amino acid in this gene
                        transcript_l[16]=aminoAcids[0]
                        transcript_l.insert(17,aminoAcids[1])
                    else:
                        transcript_l.insert(17,"")
                    transcript_l.insert(0,record.ID) # insert variant_id from our variants db
                    transcript_l.insert(21,updateDate)
                    se_records_with_ID.append(tuple(transcript_l))

            """        
            Creating list for dataframe for populations sql table   

            there are some problems with dbsnp annotation. If for a variant a parameter is missing (no value to annotate),
            then the parameter is not there. Not "prameter=''", it's just not there. This is filled here (with var_pop[i]=''), because the data 
            in db has to be consisent. Another problem is, that for a parameter there can be different number of values, 
            for example in silico tools give value per transcript. If there are different numbers of transcripts,
            then the values in in silico prameter are seperated witk comma, and for each variant the number of transcripts is different. 
            Even if we wold try to mach transcripts with in silico score, we are not sure which in silico score corresponds to which transcript,
            becuase even the number of transcripts sometimes is different than number of scores in in silico prarameter.
            That's why the values from parameter are joined into string and saved as one value in db. 
            """
            var_pop={}
            if (any(key in dbnsfp for key in record.INFO.keys())): 

                var_pop["variant_id"]=record.ID
                var_pop["pop_update"]=updateDate
                
                for i in dbnsfp:
                    for key in record.INFO.keys():                        
                        if i == key:                            
                            #joining values for a parameter into one string
                            # max 250 characters of a string
                            var_pop[i]=";".join([str(x) for x in record.INFO[key]])[:250] #!!! if not longer than 255 char
                    if i not in record.INFO.keys():
                        var_pop[i]='' #consider adding np.nan
                dbnsfp_freq.append(var_pop)


            var_detect={}
            var_detect["variant_id"]=record.ID
            var_detect["pipeline_run_id"]=pipelineRunID
            var_detect["qual"]=record.QUAL
            var_detect["filter"]=",".join([str(x) for x in record.FILTER])
            for i in qual_params_list:
                if i in record.INFO.keys(): 
                    var_detect[i]=";".join([str(x) for x in record.INFO[i]])
                if i not in record.INFO.keys():
                    var_detect[i]=""
            for i in qual_params_single:
                if i in record.INFO.keys():
                    var_detect[i]=record.INFO[i]
                if i not in record.INFO.keys():
                    var_detect[i]=""

            var_detect["gt_g"]=record.genotype(orderID).data[0] ############## check if it's also orderID written in vcf file after nextflow. It might change
            #here the ad_g0, ad_g1... , pl_g00 ... has to be replaced with the lines of code with join. For this the detections table has to be changed (ad_g col added and ad_g? cols removed),
            #db sql string needs to be corected 
            var_detect["ad_g0"]=record.genotype(orderID).data[1][0]                       
            print("AD: ", ','.join(str(x) for x in record.genotype(orderID).data[1]))  
            var_detect["ad_g1"]=record.genotype(orderID).data[1][1]
            var_detect["dp_g"]=record.genotype(orderID).data[2]
            var_detect["gq_g"]=record.genotype(orderID).data[3]
            var_detect["pl_g_00"]=record.genotype(orderID).data[4][0]
            print("PL: ", ','.join(str(x) for x in record.genotype(orderID).data[4]))         
            var_detect["pl_g_01"]=record.genotype(orderID).data[4][1]
            var_detect["pl_g_11"]=record.genotype(orderID).data[4][2]

            detect.append(var_detect)           




        n+=1
    
    sys.exit()        
##################### VEP upload ###############################    
    try:
        connection = psycopg2.connect(user = db['USER'],
                                    password = db['PASSWORD'],
                                    host = db['HOST'],
                                    database = db['NAME'])
        cursor = connection.cursor()
        # Print PostgreSQL Connection properties
        print ( connection.get_dsn_parameters(),"\n")
        # Print PostgreSQL version
        cursor.execute("SELECT version();")        
        record = cursor.fetchone()
        print("You are connected to - ", record,"\n")

        begin_time_vep = datetime.datetime.now()

        for transcript in vep_records_with_ID:        
            cursor.execute("""UPDATE vep_annotations SET allele_vep='{t[1]}', consequence_vep='{t[2]}', feature_type_vep='{t[3]}', impact_vep='{t[5]}', symbol_vep='{t[6]}', vep_update='{t[7]}' WHERE variant_id={t[0]} AND feature_vep='{t[4]}';INSERT INTO vep_annotations (variant_id, allele_vep, consequence_vep, feature_type_vep, feature_vep, impact_vep, symbol_vep, vep_update) SELECT {t[0]},'{t[1]}','{t[2]}', '{t[3]}', '{t[4]}', '{t[5]}', '{t[6]}', '{t[7]}' WHERE NOT EXISTS (SELECT (variant_id, feature_vep) FROM vep_annotations WHERE variant_id={t[0]} AND feature_vep='{t[4]}' LIMIT 1);""".format(t=transcript))   #'1', 14677, 'G', 'A', 'A', 'downstream_gene_variant', 'Transcript', 'ENST00000456328', 'MODIFIER', 'DDX11L1'
            
        ########## performance results 1k variants ############
        #WITH UPDATE in vep sql comand
            #when variants are there:
            #raw variants upload:  0:00:00.670117
            #vep upload:  0:00:08.593351

            #when variants are not there:
            #raw variants upload:  0:00:01.486114
            #vep upload:  0:00:08.416705

        #Without UPDATE in vep sql comand
            #When variants are there:
            #raw variants upload:  0:00:00.680382
            #vep upload:  0:00:02.558825
            #When variants are not there:
            #raw variants upload:  0:00:01.381642
            #vep upload:  0:00:04.480208

        

        #for transcript in vep_records:
        #with below INSERT ON CONFLICT ON CONSTRAINT (below line) the vep_annotation_id for the next transcript is increased even if transcript is not inserted (so probably inserted and deleted if exists). 
        #cursor.execute("""INSERT INTO vep_annotations (variant_id, allele_vep, consequence_vep, feature_type_vep, feature_vep, impact_vep, symbol_vep) VALUES ((SELECT variant_id FROM variants WHERE chrom='1' AND pos=14464 AND ref='A' AND alt='T'), 'T', 'downstream_gene_variant', 'Transcript', 'ENST00000450306', 'MODIFIER', 'DDX11L1') ON CONFLICT ON CONSTRAINT "vep_annotations_variant_id_feature_vep_key" DO UPDATE SET allele_vep='T', consequence_vep='downstream_gene_variant', feature_type_vep='Transcript', impact_vep='MODIFIER', symbol_vep='DDX11L3';""")
        #    cursor.execute("""UPDATE vep_annotations SET allele_vep='{4}', consequence_vep='{5}', feature_type_vep='{6}', impact_vep='{8}', symbol_vep='{9}' WHERE variant_id=(SELECT (variant_id) FROM variants WHERE chrom='{0}' AND pos={1} AND ref='{2}' AND alt='{3}' LIMIT 1) AND feature_vep='{7}';INSERT INTO vep_annotations (variant_id, allele_vep, consequence_vep, feature_type_vep, feature_vep, impact_vep, symbol_vep) SELECT (SELECT variant_id FROM variants WHERE chrom='{0}' AND pos={1} AND ref='{2}' AND alt='{3}'),'{4}','{5}', '{6}', '{7}', '{8}', '{9}' WHERE NOT EXISTS (SELECT (variant_id, feature_vep) FROM vep_annotations WHERE variant_id=(SELECT (variant_id) FROM variants WHERE chrom='{0}' AND pos={1} AND ref='{2}' AND alt='{3}' LIMIT 1) AND feature_vep='{7}' LIMIT 1);""".format(transcript[0],transcript[1],transcript[2],transcript[3],transcript[4],transcript[5],transcript[6],transcript[7],transcript[8],transcript[9]))   #'1', 14677, 'G', 'A', 'A', 'downstream_gene_variant', 'Transcript', 'ENST00000456328', 'MODIFIER', 'DDX11L1'
            #I changed gene name in vcf to check if UPDATE works properly and the gene was updatetd properly
            #We need to find a way to speed up uploading. We can do something like uploading only variants which are not there yet. To do it we need to check first which variants are there, upload remaining variants. 
            #We could also skip update for annotation, but then we have to write a SQL to update the values for annotation tables of all variants and run update once a while (once a month or so)
                      

        print("vep upload: ", datetime.datetime.now() - begin_time_vep)

        ########## performance results 1k variants ############
        #Without UPDATE in vep sql comand
            #when variants are there:
            #raw variants upload:  0:00:00.299635
            #vep upload:  0:00:03.833186

            #when variants are not there:
            #raw variants upload:  0:00:01.458215
            #vep upload:  0:00:13.931331


        connection.commit()
        print("Transaction completed successfully ")

    except (Exception, psycopg2.Error) as error :        
        """
        Shall we add below 2 lines of code? Does it work?
        if connnection:
            connection.rollback()
        """
        print ("Error while connecting to PostgreSQL", error)
    finally:
        #closing database connection.
            if(connection):
                cursor.close()
                connection.close()
                print("PostgreSQL connection is closed")

####################### snpEff upload ########################

    try:
        connection = psycopg2.connect(user = db['USER'],
                                    password = db['PASSWORD'],
                                    host = db['HOST'],
                                    database = db['NAME'])
        cursor = connection.cursor()
        # Print PostgreSQL Connection properties
        print ( connection.get_dsn_parameters(),"\n")
        # Print PostgreSQL version
        cursor.execute("SELECT version();")        
        record = cursor.fetchone()
        print("You are connected to - ", record,"\n")

        begin_time_se = datetime.datetime.now()
        

        for transcript in se_records_with_ID:        
            cursor.execute("""UPDATE se_annotations SET allele_se='{t[1]}', annotation_se='{t[2]}', annotation_imp_se='{t[3]}', gene_name='{t[4]}', gene_id='{t[5]}', feature_type_se=$${t[6]}$$, transcript_biotype='{t[8]}', rank_pos_se='{t[9]}', rank_length_se='{t[10]}', hgvs_c='{t[11]}', hgvs_p='{t[12]}', cdna_pos='{t[13]}', cdna_length='{t[14]}', cds_pos='{t[15]}', cds_length='{t[16]}', aa_pos='{t[17]}', aa_length='{t[18]}', distance='{t[19]}', errors='{t[20]}', se_update='{t[21]}' WHERE variant_id={t[0]} AND feature_id_se='{t[7]}';INSERT INTO se_annotations (variant_id, allele_se, annotation_se, annotation_imp_se, gene_name, gene_id, feature_type_se, feature_id_se, transcript_biotype, rank_pos_se, rank_length_se, hgvs_c, hgvs_p, cdna_pos, cdna_length, cds_pos, cds_length, aa_pos, aa_length, distance, errors, se_update) SELECT '{t[0]}', '{t[1]}', '{t[2]}', '{t[3]}', '{t[4]}', '{t[5]}', $${t[6]}$$, '{t[7]}', '{t[8]}', '{t[9]}', '{t[10]}', '{t[11]}', '{t[12]}', '{t[13]}', '{t[14]}', '{t[15]}', '{t[16]}', '{t[17]}', '{t[18]}', '{t[19]}', '{t[20]}', '{t[21]}'  WHERE NOT EXISTS (SELECT (variant_id, feature_id_se) FROM se_annotations WHERE variant_id={t[0]} AND feature_id_se='{t[7]}' LIMIT 1);""".format(t=transcript))   #'1', 14677, 'G', 'A', 'A', 'downstream_gene_variant', 'Transcript', 'ENST00000456328', 'MODIFIER', 'DDX11L1'


        print("se upload: ", datetime.datetime.now() - begin_time_se)
        #se upload when the data is there:  0:00:02.256056
        #se upload when the data is not there:  0:00:03.169911



        connection.commit()
        print("Transaction completed successfully ")

    except (Exception, psycopg2.Error) as error :        
        """
        Shall we add below 2 lines of code? Does it work?
        if connnection:
            connection.rollback()
        """
        print ("Error while connecting to PostgreSQL", error)
    finally:
        #closing database connection.
            if(connection):
                cursor.close()
                connection.close()
                print("PostgreSQL connection is closed")


####################### populations upload ########################

    try:
        connection = psycopg2.connect(user = db['USER'],
                                    password = db['PASSWORD'],
                                    host = db['HOST'],
                                    database = db['NAME'])
        cursor = connection.cursor()
        # Print PostgreSQL Connection properties
        print ( connection.get_dsn_parameters(),"\n")
        # Print PostgreSQL version
        cursor.execute("SELECT version();")        
        record = cursor.fetchone()
        print("You are connected to - ", record,"\n")

        begin_time_pop = datetime.datetime.now()
        

        for pop in dbnsfp_freq:
            #cursor.execute("""UPDATE populations SET metasvm_pred='{1}', metalr_pred='{2}', m_cap_pred='{3}', fathmm_pred='{4}', vest4_score='{5}', aloft_pred='{6}', rs_dbsnp151='{7}', sift_pred='{8}', mutationtaster_pred='{9}', polyphen2_hvar_pred='{10}', Polyphen2_hdiv_pred='{11}', cadd_phred='{12}', revel_score='{13}', revel_rankscore='{14}', g1000gp3_ac='{15}', g1000gp3_af='{16}', g1000gp3_afr_ac='{17}', g1000gp3_afr_af='{18}', g1000gp3_eur_ac='{19}', g1000gp3_eur_af='{20}', g1000gp3_amr_ac='{21}', g1000gp3_amr_af='{22}', g1000gp3_eas_ac='{23}', g1000gp3_eas_af='{24}', g1000gp3_sas_ac='{25}', g1000gp3_sas_af='{26}', twinsuk_ac='{27}', twinsuk_af='{28}', alspac_ac='{29}', alspac_af='{30}', esp6500_aa_ac='{31}', esp6500_aa_af='{32}', esp6500_ea_ac='{33}', esp6500_ea_af='{34}', exac_adj_ac='{35}', exac_adj_af='{36}', exac_afr_ac='{37}', exac_afr_af='{38}', exac_amr_ac='{39}', exac_amr_af='{40}', exac_eas_ac='{41}', exac_eas_af='{42}', exac_fin_ac='{43}', exac_fin_af='{44}', exac_nfe_ac='{45}', exac_nfe_af='{46}', exac_sas_ac='{47}', exac_sas_af='{48}', gnomad_exomes_ac='{49}', gnomad_exomes_af='{50}', gnomad_exomes_afr_ac='{51}', gnomad_exomes_afr_af='{52}', gnomad_exomes_amr_ac='{53}', gnomad_exomes_amr_af='{54}', gnomad_exomes_asj_ac='{55}', gnomad_exomes_asj_af='{56}', gnomad_exomes_eas_ac='{57}', gnomad_exomes_eas_af='{58}', gnomad_exomes_fin_ac='{59}', gnomad_exomes_fin_af='{60}', gnomad_exomes_nfe_ac='{61}', gnomad_exomes_nfe_af='{62}', gnomad_exomes_sas_ac='{63}', gnomad_exomes_sas_af='{64}', clinvar_clnsig='{65}', clinvar_var_source='{66}', clinvar_review='{67}', ensembl_geneid='{68}', ensembl_transcriptid='{69}', ensembl_proteinid='{70}',hg19_pos_1_based='{71}' WHERE variant_id={0}; INSERT INTO populations (variant_id, metasvm_pred, metalr_pred, m_cap_pred, fathmm_pred, vest4_score, aloft_pred, rs_dbsnp151, sift_pred, mutationtaster_pred, polyphen2_hvar_pred, Polyphen2_hdiv_pred, cadd_phred, revel_score, revel_rankscore, g1000gp3_ac, g1000gp3_af, g1000gp3_afr_ac, g1000gp3_afr_af, g1000gp3_eur_ac, g1000gp3_eur_af, g1000gp3_amr_ac, g1000gp3_amr_af, g1000gp3_eas_ac, g1000gp3_eas_af, g1000gp3_sas_ac, g1000gp3_sas_af, twinsuk_ac, twinsuk_af, alspac_ac, alspac_af, esp6500_aa_ac, esp6500_aa_af, esp6500_ea_ac, esp6500_ea_af, exac_adj_ac, exac_adj_af, exac_afr_ac, exac_afr_af, exac_amr_ac, exac_amr_af, exac_eas_ac, exac_eas_af, exac_fin_ac, exac_fin_af, exac_nfe_ac, exac_nfe_af, exac_sas_ac, exac_sas_af, gnomad_exomes_ac, gnomad_exomes_af, gnomad_exomes_afr_ac, gnomad_exomes_afr_af, gnomad_exomes_amr_ac, gnomad_exomes_amr_af, gnomad_exomes_asj_ac, gnomad_exomes_asj_af, gnomad_exomes_eas_ac, gnomad_exomes_eas_af, gnomad_exomes_fin_ac, gnomad_exomes_fin_af, gnomad_exomes_nfe_ac, gnomad_exomes_nfe_af, gnomad_exomes_sas_ac, gnomad_exomes_sas_af, clinvar_clnsig, clinvar_var_source, clinvar_review, ensembl_geneid, ensembl_transcriptid, ensembl_proteinid, hg19_pos_1_based) SELECT '{0}', '{1}', '{2}', '{3}', '{4}', '{5}', '{6}', '{7}', '{8}', '{9}', '{10}', '{11}', '{12}', '{13}', '{14}', '{15}', '{16}', '{17}', '{18}', '{19}', '{20}', '{21}', '{22}', '{23}', '{24}', '{25}', '{26}', '{27}', '{28}', '{29}', '{30}', '{31}', '{32}', '{33}', '{34}', '{35}', '{36}', '{37}', '{38}', '{39}', '{40}', '{41}', '{42}', '{43}', '{44}', '{45}', '{46}', '{47}', '{48}', '{49}', '{50}', '{51}', '{52}', '{53}', '{54}', '{55}', '{56}', '{57}', '{58}', '{59}', '{60}', '{61}', '{62}', '{63}', '{64}', '{65}', '{66}', '{67}', '{68}', '{69}', '{70}', '{71}' WHERE NOT EXISTS (SELECT (variant_id) FROM populations WHERE variant_id={0} LIMIT 1);""".format(pop["variant_id"],pop["dbNSFP_MetaSVM_pred"],pop["dbNSFP_MetaLR_pred"],pop["dbNSFP_M_CAP_pred"],pop["dbNSFP_FATHMM_pred"],pop["dbNSFP_VEST4_score"],pop["dbNSFP_Aloft_pred"],pop["dbNSFP_rs_dbSNP151"],pop["dbNSFP_SIFT_pred"],pop["dbNSFP_MutationTaster_pred"],pop["dbNSFP_Polyphen2_HVAR_pred"],pop["dbNSFP_Polyphen2_HDIV_pred"],pop["dbNSFP_CADD_phred"],pop["dbNSFP_REVEL_score"],pop["dbNSFP_REVEL_rankscore"],pop["dbNSFP_1000Gp3_AC"],pop["dbNSFP_1000Gp3_AF"],pop["dbNSFP_1000Gp3_AFR_AC"],pop["dbNSFP_1000Gp3_AFR_AF"],pop["dbNSFP_1000Gp3_EUR_AC"],pop["dbNSFP_1000Gp3_EUR_AF"],pop["dbNSFP_1000Gp3_AMR_AC"],pop["dbNSFP_1000Gp3_AMR_AF"],pop["dbNSFP_1000Gp3_EAS_AC"],pop["dbNSFP_1000Gp3_EAS_AF"],pop["dbNSFP_1000Gp3_SAS_AC"],pop["dbNSFP_1000Gp3_SAS_AF"],pop["dbNSFP_TWINSUK_AC"],pop["dbNSFP_TWINSUK_AF"],pop["dbNSFP_ALSPAC_AC"],pop["dbNSFP_ALSPAC_AF"],pop["dbNSFP_ESP6500_AA_AC"],pop["dbNSFP_ESP6500_AA_AF"],pop["dbNSFP_ESP6500_EA_AC"],pop["dbNSFP_ESP6500_EA_AF"],pop["dbNSFP_ExAC_Adj_AC"],pop["dbNSFP_ExAC_Adj_AF"],pop["dbNSFP_ExAC_AFR_AC"],pop["dbNSFP_ExAC_AFR_AF"],pop["dbNSFP_ExAC_AMR_AC"],pop["dbNSFP_ExAC_AMR_AF"],pop["dbNSFP_ExAC_EAS_AC"],pop["dbNSFP_ExAC_EAS_AF"],pop["dbNSFP_ExAC_FIN_AC"],pop["dbNSFP_ExAC_FIN_AF"],pop["dbNSFP_ExAC_NFE_AC"],pop["dbNSFP_ExAC_NFE_AF"],pop["dbNSFP_ExAC_SAS_AC"],pop["dbNSFP_ExAC_SAS_AF"],pop["dbNSFP_gnomAD_exomes_AC"],pop["dbNSFP_gnomAD_exomes_AF"],pop["dbNSFP_gnomAD_exomes_AFR_AC"],pop["dbNSFP_gnomAD_exomes_AFR_AF"],pop["dbNSFP_gnomAD_exomes_AMR_AC"],pop["dbNSFP_gnomAD_exomes_AMR_AF"],pop["dbNSFP_gnomAD_exomes_ASJ_AC"],pop["dbNSFP_gnomAD_exomes_ASJ_AF"],pop["dbNSFP_gnomAD_exomes_EAS_AC"],pop["dbNSFP_gnomAD_exomes_EAS_AF"],pop["dbNSFP_gnomAD_exomes_FIN_AC"],pop["dbNSFP_gnomAD_exomes_FIN_AF"],pop["dbNSFP_gnomAD_exomes_NFE_AC"],pop["dbNSFP_gnomAD_exomes_NFE_AF"],pop["dbNSFP_gnomAD_exomes_SAS_AC"],pop["dbNSFP_gnomAD_exomes_SAS_AF"],pop["dbNSFP_clinvar_clnsig"],pop["dbNSFP_clinvar_var_source"],pop["dbNSFP_clinvar_review"],pop["dbNSFP_Ensembl_geneid"],pop["dbNSFP_Ensembl_transcriptid"],pop["dbNSFP_Ensembl_proteinid"],pop["dbNSFP_hg19_pos_1_based_"]))
            cursor.execute("""UPDATE populations SET metasvm_pred='{dbNSFP_MetaSVM_pred}', metalr_pred='{dbNSFP_MetaLR_pred}', m_cap_pred='{dbNSFP_M_CAP_pred}', fathmm_pred='{dbNSFP_FATHMM_pred}', vest4_score='{dbNSFP_VEST4_score}', aloft_pred='{dbNSFP_Aloft_pred}', rs_dbsnp151='{dbNSFP_rs_dbSNP151}', sift_pred='{dbNSFP_SIFT_pred}', mutationtaster_pred='{dbNSFP_MutationTaster_pred}', polyphen2_hvar_pred='{dbNSFP_Polyphen2_HVAR_pred}', Polyphen2_hdiv_pred='{dbNSFP_Polyphen2_HDIV_pred}', cadd_phred='{dbNSFP_CADD_phred}', revel_score='{dbNSFP_REVEL_score}', revel_rankscore='{dbNSFP_REVEL_rankscore}', g1000gp3_ac='{dbNSFP_1000Gp3_AC}', g1000gp3_af='{dbNSFP_1000Gp3_AF}', g1000gp3_afr_ac='{dbNSFP_1000Gp3_AFR_AC}', g1000gp3_afr_af='{dbNSFP_1000Gp3_AFR_AF}', g1000gp3_eur_ac='{dbNSFP_1000Gp3_EUR_AC}', g1000gp3_eur_af='{dbNSFP_1000Gp3_EUR_AF}', g1000gp3_amr_ac='{dbNSFP_1000Gp3_AMR_AC}', g1000gp3_amr_af='{dbNSFP_1000Gp3_AMR_AF}', g1000gp3_eas_ac='{dbNSFP_1000Gp3_EAS_AC}', g1000gp3_eas_af='{dbNSFP_1000Gp3_EAS_AF}', g1000gp3_sas_ac='{dbNSFP_1000Gp3_SAS_AC}', g1000gp3_sas_af='{dbNSFP_1000Gp3_SAS_AF}', twinsuk_ac='{dbNSFP_TWINSUK_AC}', twinsuk_af='{dbNSFP_TWINSUK_AF}', alspac_ac='{dbNSFP_ALSPAC_AC}', alspac_af='{dbNSFP_ALSPAC_AF}', esp6500_aa_ac='{dbNSFP_ESP6500_AA_AC}', esp6500_aa_af='{dbNSFP_ESP6500_AA_AF}', esp6500_ea_ac='{dbNSFP_ESP6500_EA_AC}', esp6500_ea_af='{dbNSFP_ESP6500_EA_AF}', exac_adj_ac='{dbNSFP_ExAC_Adj_AC}', exac_adj_af='{dbNSFP_ExAC_Adj_AF}', exac_afr_ac='{dbNSFP_ExAC_AFR_AC}', exac_afr_af='{dbNSFP_ExAC_AFR_AF}', exac_amr_ac='{dbNSFP_ExAC_AMR_AC}', exac_amr_af='{dbNSFP_ExAC_AMR_AF}', exac_eas_ac='{dbNSFP_ExAC_EAS_AC}', exac_eas_af='{dbNSFP_ExAC_EAS_AF}', exac_fin_ac='{dbNSFP_ExAC_FIN_AC}', exac_fin_af='{dbNSFP_ExAC_FIN_AF}', exac_nfe_ac='{dbNSFP_ExAC_NFE_AC}', exac_nfe_af='{dbNSFP_ExAC_NFE_AF}', exac_sas_ac='{dbNSFP_ExAC_SAS_AC}', exac_sas_af='{dbNSFP_ExAC_SAS_AF}', gnomad_exomes_ac='{dbNSFP_gnomAD_exomes_AC}', gnomad_exomes_af='{dbNSFP_gnomAD_exomes_AF}', gnomad_exomes_afr_ac='{dbNSFP_gnomAD_exomes_AFR_AC}', gnomad_exomes_afr_af='{dbNSFP_gnomAD_exomes_AFR_AF}', gnomad_exomes_amr_ac='{dbNSFP_gnomAD_exomes_AMR_AC}', gnomad_exomes_amr_af='{dbNSFP_gnomAD_exomes_AMR_AF}', gnomad_exomes_asj_ac='{dbNSFP_gnomAD_exomes_ASJ_AC}', gnomad_exomes_asj_af='{dbNSFP_gnomAD_exomes_ASJ_AF}', gnomad_exomes_eas_ac='{dbNSFP_gnomAD_exomes_EAS_AC}', gnomad_exomes_eas_af='{dbNSFP_gnomAD_exomes_EAS_AF}', gnomad_exomes_fin_ac='{dbNSFP_gnomAD_exomes_FIN_AC}', gnomad_exomes_fin_af='{dbNSFP_gnomAD_exomes_FIN_AF}', gnomad_exomes_nfe_ac='{dbNSFP_gnomAD_exomes_NFE_AC}', gnomad_exomes_nfe_af='{dbNSFP_gnomAD_exomes_NFE_AF}', gnomad_exomes_sas_ac='{dbNSFP_gnomAD_exomes_SAS_AC}', gnomad_exomes_sas_af='{dbNSFP_gnomAD_exomes_SAS_AF}', clinvar_clnsig='{dbNSFP_clinvar_clnsig}', clinvar_var_source=$${dbNSFP_clinvar_var_source}$$, clinvar_review='{dbNSFP_clinvar_review}', ensembl_geneid='{dbNSFP_Ensembl_geneid}', ensembl_transcriptid='{dbNSFP_Ensembl_transcriptid}', ensembl_proteinid='{dbNSFP_Ensembl_proteinid}',hg19_pos_1_based='{dbNSFP_hg19_pos_1_based_}', pop_update='{pop_update}' WHERE variant_id={variant_id}; INSERT INTO populations (variant_id, metasvm_pred, metalr_pred, m_cap_pred, fathmm_pred, vest4_score, aloft_pred, rs_dbsnp151, sift_pred, mutationtaster_pred, polyphen2_hvar_pred, Polyphen2_hdiv_pred, cadd_phred, revel_score, revel_rankscore, g1000gp3_ac, g1000gp3_af, g1000gp3_afr_ac, g1000gp3_afr_af, g1000gp3_eur_ac, g1000gp3_eur_af, g1000gp3_amr_ac, g1000gp3_amr_af, g1000gp3_eas_ac, g1000gp3_eas_af, g1000gp3_sas_ac, g1000gp3_sas_af, twinsuk_ac, twinsuk_af, alspac_ac, alspac_af, esp6500_aa_ac, esp6500_aa_af, esp6500_ea_ac, esp6500_ea_af, exac_adj_ac, exac_adj_af, exac_afr_ac, exac_afr_af, exac_amr_ac, exac_amr_af, exac_eas_ac, exac_eas_af, exac_fin_ac, exac_fin_af, exac_nfe_ac, exac_nfe_af, exac_sas_ac, exac_sas_af, gnomad_exomes_ac, gnomad_exomes_af, gnomad_exomes_afr_ac, gnomad_exomes_afr_af, gnomad_exomes_amr_ac, gnomad_exomes_amr_af, gnomad_exomes_asj_ac, gnomad_exomes_asj_af, gnomad_exomes_eas_ac, gnomad_exomes_eas_af, gnomad_exomes_fin_ac, gnomad_exomes_fin_af, gnomad_exomes_nfe_ac, gnomad_exomes_nfe_af, gnomad_exomes_sas_ac, gnomad_exomes_sas_af, clinvar_clnsig, clinvar_var_source, clinvar_review, ensembl_geneid, ensembl_transcriptid, ensembl_proteinid, hg19_pos_1_based, pop_update) SELECT '{variant_id}', '{dbNSFP_MetaSVM_pred}', '{dbNSFP_MetaLR_pred}', '{dbNSFP_M_CAP_pred}', '{dbNSFP_FATHMM_pred}', '{dbNSFP_VEST4_score}', '{dbNSFP_Aloft_pred}', '{dbNSFP_rs_dbSNP151}', '{dbNSFP_SIFT_pred}', '{dbNSFP_MutationTaster_pred}', '{dbNSFP_Polyphen2_HVAR_pred}', '{dbNSFP_Polyphen2_HDIV_pred}', '{dbNSFP_CADD_phred}', '{dbNSFP_REVEL_score}', '{dbNSFP_REVEL_rankscore}', '{dbNSFP_1000Gp3_AC}', '{dbNSFP_1000Gp3_AF}', '{dbNSFP_1000Gp3_AFR_AC}', '{dbNSFP_1000Gp3_AFR_AF}', '{dbNSFP_1000Gp3_EUR_AC}', '{dbNSFP_1000Gp3_EUR_AF}', '{dbNSFP_1000Gp3_AMR_AC}', '{dbNSFP_1000Gp3_AMR_AF}', '{dbNSFP_1000Gp3_EAS_AC}', '{dbNSFP_1000Gp3_EAS_AF}', '{dbNSFP_1000Gp3_SAS_AC}', '{dbNSFP_1000Gp3_SAS_AF}', '{dbNSFP_TWINSUK_AC}', '{dbNSFP_TWINSUK_AF}', '{dbNSFP_ALSPAC_AC}', '{dbNSFP_ALSPAC_AF}', '{dbNSFP_ESP6500_AA_AC}', '{dbNSFP_ESP6500_AA_AF}', '{dbNSFP_ESP6500_EA_AC}', '{dbNSFP_ESP6500_EA_AF}', '{dbNSFP_ExAC_Adj_AC}', '{dbNSFP_ExAC_Adj_AF}', '{dbNSFP_ExAC_AFR_AC}', '{dbNSFP_ExAC_AFR_AF}', '{dbNSFP_ExAC_AMR_AC}', '{dbNSFP_ExAC_AMR_AF}', '{dbNSFP_ExAC_EAS_AC}', '{dbNSFP_ExAC_EAS_AF}', '{dbNSFP_ExAC_FIN_AC}', '{dbNSFP_ExAC_FIN_AF}', '{dbNSFP_ExAC_NFE_AC}', '{dbNSFP_ExAC_NFE_AF}', '{dbNSFP_ExAC_SAS_AC}', '{dbNSFP_ExAC_SAS_AF}', '{dbNSFP_gnomAD_exomes_AC}', '{dbNSFP_gnomAD_exomes_AF}', '{dbNSFP_gnomAD_exomes_AFR_AC}', '{dbNSFP_gnomAD_exomes_AFR_AF}', '{dbNSFP_gnomAD_exomes_AMR_AC}', '{dbNSFP_gnomAD_exomes_AMR_AF}', '{dbNSFP_gnomAD_exomes_ASJ_AC}', '{dbNSFP_gnomAD_exomes_ASJ_AF}', '{dbNSFP_gnomAD_exomes_EAS_AC}', '{dbNSFP_gnomAD_exomes_EAS_AF}', '{dbNSFP_gnomAD_exomes_FIN_AC}', '{dbNSFP_gnomAD_exomes_FIN_AF}', '{dbNSFP_gnomAD_exomes_NFE_AC}', '{dbNSFP_gnomAD_exomes_NFE_AF}', '{dbNSFP_gnomAD_exomes_SAS_AC}', '{dbNSFP_gnomAD_exomes_SAS_AF}', '{dbNSFP_clinvar_clnsig}', $${dbNSFP_clinvar_var_source}$$, '{dbNSFP_clinvar_review}', '{dbNSFP_Ensembl_geneid}', '{dbNSFP_Ensembl_transcriptid}', '{dbNSFP_Ensembl_proteinid}', '{dbNSFP_hg19_pos_1_based_}', '{pop_update}' WHERE NOT EXISTS (SELECT (variant_id) FROM populations WHERE variant_id={variant_id} LIMIT 1);""".format(**pop))
           
        
        print("pop upload: ", datetime.datetime.now() - begin_time_pop)
        #pop upload:  0:00:00.035055


        connection.commit()
        print("Transaction completed successfully ")

    except (Exception, psycopg2.Error) as error :        
        """
        Shall we add below 2 lines of code? Does it work?
        if connnection:
            connection.rollback()
        """
        print ("Error while connecting to PostgreSQL", error)
    finally:
        #closing database connection.
            if(connection):
                cursor.close()
                connection.close()
                print("PostgreSQL connection is closed")

####################### detections upload ########################

    try:
        connection = psycopg2.connect(user = db['USER'],
                                    password = db['PASSWORD'],
                                    host = db['HOST'],
                                    database = db['NAME'])
        cursor = connection.cursor()
        # Print PostgreSQL Connection properties
        print ( connection.get_dsn_parameters(),"\n")
        # Print PostgreSQL version
        cursor.execute("SELECT version();")        
        record = cursor.fetchone()
        print("You are connected to - ", record,"\n")

        begin_time_detect = datetime.datetime.now()
        

        for var_detect in detect:              
            #cursor.execute("""UPDATE populations SET metasvm_pred='{1}', metalr_pred='{2}', m_cap_pred='{3}', fathmm_pred='{4}', vest4_score='{5}', aloft_pred='{6}', rs_dbsnp151='{7}', sift_pred='{8}', mutationtaster_pred='{9}', polyphen2_hvar_pred='{10}', Polyphen2_hdiv_pred='{11}', cadd_phred='{12}', revel_score='{13}', revel_rankscore='{14}', g1000gp3_ac='{15}', g1000gp3_af='{16}', g1000gp3_afr_ac='{17}', g1000gp3_afr_af='{18}', g1000gp3_eur_ac='{19}', g1000gp3_eur_af='{20}', g1000gp3_amr_ac='{21}', g1000gp3_amr_af='{22}', g1000gp3_eas_ac='{23}', g1000gp3_eas_af='{24}', g1000gp3_sas_ac='{25}', g1000gp3_sas_af='{26}', twinsuk_ac='{27}', twinsuk_af='{28}', alspac_ac='{29}', alspac_af='{30}', esp6500_aa_ac='{31}', esp6500_aa_af='{32}', esp6500_ea_ac='{33}', esp6500_ea_af='{34}', exac_adj_ac='{35}', exac_adj_af='{36}', exac_afr_ac='{37}', exac_afr_af='{38}', exac_amr_ac='{39}', exac_amr_af='{40}', exac_eas_ac='{41}', exac_eas_af='{42}', exac_fin_ac='{43}', exac_fin_af='{44}', exac_nfe_ac='{45}', exac_nfe_af='{46}', exac_sas_ac='{47}', exac_sas_af='{48}', gnomad_exomes_ac='{49}', gnomad_exomes_af='{50}', gnomad_exomes_afr_ac='{51}', gnomad_exomes_afr_af='{52}', gnomad_exomes_amr_ac='{53}', gnomad_exomes_amr_af='{54}', gnomad_exomes_asj_ac='{55}', gnomad_exomes_asj_af='{56}', gnomad_exomes_eas_ac='{57}', gnomad_exomes_eas_af='{58}', gnomad_exomes_fin_ac='{59}', gnomad_exomes_fin_af='{60}', gnomad_exomes_nfe_ac='{61}', gnomad_exomes_nfe_af='{62}', gnomad_exomes_sas_ac='{63}', gnomad_exomes_sas_af='{64}', clinvar_clnsig='{65}', clinvar_var_source='{66}', clinvar_review='{67}', ensembl_geneid='{68}', ensembl_transcriptid='{69}', ensembl_proteinid='{70}',hg19_pos_1_based='{71}' WHERE variant_id={0}; INSERT INTO populations (variant_id, metasvm_pred, metalr_pred, m_cap_pred, fathmm_pred, vest4_score, aloft_pred, rs_dbsnp151, sift_pred, mutationtaster_pred, polyphen2_hvar_pred, Polyphen2_hdiv_pred, cadd_phred, revel_score, revel_rankscore, g1000gp3_ac, g1000gp3_af, g1000gp3_afr_ac, g1000gp3_afr_af, g1000gp3_eur_ac, g1000gp3_eur_af, g1000gp3_amr_ac, g1000gp3_amr_af, g1000gp3_eas_ac, g1000gp3_eas_af, g1000gp3_sas_ac, g1000gp3_sas_af, twinsuk_ac, twinsuk_af, alspac_ac, alspac_af, esp6500_aa_ac, esp6500_aa_af, esp6500_ea_ac, esp6500_ea_af, exac_adj_ac, exac_adj_af, exac_afr_ac, exac_afr_af, exac_amr_ac, exac_amr_af, exac_eas_ac, exac_eas_af, exac_fin_ac, exac_fin_af, exac_nfe_ac, exac_nfe_af, exac_sas_ac, exac_sas_af, gnomad_exomes_ac, gnomad_exomes_af, gnomad_exomes_afr_ac, gnomad_exomes_afr_af, gnomad_exomes_amr_ac, gnomad_exomes_amr_af, gnomad_exomes_asj_ac, gnomad_exomes_asj_af, gnomad_exomes_eas_ac, gnomad_exomes_eas_af, gnomad_exomes_fin_ac, gnomad_exomes_fin_af, gnomad_exomes_nfe_ac, gnomad_exomes_nfe_af, gnomad_exomes_sas_ac, gnomad_exomes_sas_af, clinvar_clnsig, clinvar_var_source, clinvar_review, ensembl_geneid, ensembl_transcriptid, ensembl_proteinid, hg19_pos_1_based) SELECT '{0}', '{1}', '{2}', '{3}', '{4}', '{5}', '{6}', '{7}', '{8}', '{9}', '{10}', '{11}', '{12}', '{13}', '{14}', '{15}', '{16}', '{17}', '{18}', '{19}', '{20}', '{21}', '{22}', '{23}', '{24}', '{25}', '{26}', '{27}', '{28}', '{29}', '{30}', '{31}', '{32}', '{33}', '{34}', '{35}', '{36}', '{37}', '{38}', '{39}', '{40}', '{41}', '{42}', '{43}', '{44}', '{45}', '{46}', '{47}', '{48}', '{49}', '{50}', '{51}', '{52}', '{53}', '{54}', '{55}', '{56}', '{57}', '{58}', '{59}', '{60}', '{61}', '{62}', '{63}', '{64}', '{65}', '{66}', '{67}', '{68}', '{69}', '{70}', '{71}' WHERE NOT EXISTS (SELECT (variant_id) FROM populations WHERE variant_id={0} LIMIT 1);""".format(pop["variant_id"],pop["dbNSFP_MetaSVM_pred"],pop["dbNSFP_MetaLR_pred"],pop["dbNSFP_M_CAP_pred"],pop["dbNSFP_FATHMM_pred"],pop["dbNSFP_VEST4_score"],pop["dbNSFP_Aloft_pred"],pop["dbNSFP_rs_dbSNP151"],pop["dbNSFP_SIFT_pred"],pop["dbNSFP_MutationTaster_pred"],pop["dbNSFP_Polyphen2_HVAR_pred"],pop["dbNSFP_Polyphen2_HDIV_pred"],pop["dbNSFP_CADD_phred"],pop["dbNSFP_REVEL_score"],pop["dbNSFP_REVEL_rankscore"],pop["dbNSFP_1000Gp3_AC"],pop["dbNSFP_1000Gp3_AF"],pop["dbNSFP_1000Gp3_AFR_AC"],pop["dbNSFP_1000Gp3_AFR_AF"],pop["dbNSFP_1000Gp3_EUR_AC"],pop["dbNSFP_1000Gp3_EUR_AF"],pop["dbNSFP_1000Gp3_AMR_AC"],pop["dbNSFP_1000Gp3_AMR_AF"],pop["dbNSFP_1000Gp3_EAS_AC"],pop["dbNSFP_1000Gp3_EAS_AF"],pop["dbNSFP_1000Gp3_SAS_AC"],pop["dbNSFP_1000Gp3_SAS_AF"],pop["dbNSFP_TWINSUK_AC"],pop["dbNSFP_TWINSUK_AF"],pop["dbNSFP_ALSPAC_AC"],pop["dbNSFP_ALSPAC_AF"],pop["dbNSFP_ESP6500_AA_AC"],pop["dbNSFP_ESP6500_AA_AF"],pop["dbNSFP_ESP6500_EA_AC"],pop["dbNSFP_ESP6500_EA_AF"],pop["dbNSFP_ExAC_Adj_AC"],pop["dbNSFP_ExAC_Adj_AF"],pop["dbNSFP_ExAC_AFR_AC"],pop["dbNSFP_ExAC_AFR_AF"],pop["dbNSFP_ExAC_AMR_AC"],pop["dbNSFP_ExAC_AMR_AF"],pop["dbNSFP_ExAC_EAS_AC"],pop["dbNSFP_ExAC_EAS_AF"],pop["dbNSFP_ExAC_FIN_AC"],pop["dbNSFP_ExAC_FIN_AF"],pop["dbNSFP_ExAC_NFE_AC"],pop["dbNSFP_ExAC_NFE_AF"],pop["dbNSFP_ExAC_SAS_AC"],pop["dbNSFP_ExAC_SAS_AF"],pop["dbNSFP_gnomAD_exomes_AC"],pop["dbNSFP_gnomAD_exomes_AF"],pop["dbNSFP_gnomAD_exomes_AFR_AC"],pop["dbNSFP_gnomAD_exomes_AFR_AF"],pop["dbNSFP_gnomAD_exomes_AMR_AC"],pop["dbNSFP_gnomAD_exomes_AMR_AF"],pop["dbNSFP_gnomAD_exomes_ASJ_AC"],pop["dbNSFP_gnomAD_exomes_ASJ_AF"],pop["dbNSFP_gnomAD_exomes_EAS_AC"],pop["dbNSFP_gnomAD_exomes_EAS_AF"],pop["dbNSFP_gnomAD_exomes_FIN_AC"],pop["dbNSFP_gnomAD_exomes_FIN_AF"],pop["dbNSFP_gnomAD_exomes_NFE_AC"],pop["dbNSFP_gnomAD_exomes_NFE_AF"],pop["dbNSFP_gnomAD_exomes_SAS_AC"],pop["dbNSFP_gnomAD_exomes_SAS_AF"],pop["dbNSFP_clinvar_clnsig"],pop["dbNSFP_clinvar_var_source"],pop["dbNSFP_clinvar_review"],pop["dbNSFP_Ensembl_geneid"],pop["dbNSFP_Ensembl_transcriptid"],pop["dbNSFP_Ensembl_proteinid"],pop["dbNSFP_hg19_pos_1_based_"]))
            cursor.execute("""UPDATE detections SET qual='{qual}', filter='{filter}', ac='{AC}', af='{AF}', an='{AN}', baseqranksum='{BaseQRankSum}', dp='{DP}', ds='{DS}', pos_end='{END}', excessHet='{ExcessHet}', fs='{FS}', inbreedingcoeff='{InbreedingCoeff}', mleac='{MLEAC}', mleaf='{MLEAF}', mq='{MQ}', mqranksum='{MQRankSum}', n_t_s='{NEGATIVE_TRAIN_SITE}', p_t_s='{POSITIVE_TRAIN_SITE}', qd='{QD}', readposranksum='{ReadPosRankSum}', vqslod='{VQSLOD}', culprint='{culprit}', sor='{SOR}', gt_g='{gt_g}', ad_g0='{ad_g0}', ad_g1='{ad_g1}', dp_g='{dp_g}', gq_g='{gq_g}', pl_g_00='{pl_g_00}', pl_g_01='{pl_g_01}', pl_g_11='{pl_g_11}' WHERE variant_id={variant_id} AND pipeline_run_id={pipeline_run_id}; INSERT INTO detections (variant_id, pipeline_run_id, qual, filter, ac, af, an, baseqranksum, dp, ds, pos_end, excessHet, fs, inbreedingcoeff, mleac, mleaf, mq, mqranksum, n_t_s, p_t_s, qd, readposranksum, vqslod, culprint, sor, gt_g, ad_g0, ad_g1, dp_g, gq_g, pl_g_00, pl_g_01, pl_g_11) SELECT '{variant_id}', '{pipeline_run_id}', '{qual}', '{filter}', '{AC}', '{AF}', '{AN}', '{BaseQRankSum}', '{DP}', '{DS}', '{END}', '{ExcessHet}', '{FS}', '{InbreedingCoeff}', '{MLEAC}', '{MLEAF}', '{MQ}', '{MQRankSum}', '{NEGATIVE_TRAIN_SITE}', '{POSITIVE_TRAIN_SITE}', '{QD}', '{ReadPosRankSum}', '{VQSLOD}', '{culprit}', '{SOR}', '{gt_g}', '{ad_g0}', '{ad_g1}', '{dp_g}', '{gq_g}', '{pl_g_00}', '{pl_g_01}', '{pl_g_11}' WHERE NOT EXISTS (SELECT (variant_id, pipeline_run_id) FROM detections WHERE variant_id={variant_id} AND pipeline_run_id={pipeline_run_id} LIMIT 1);""".format(**var_detect))
           
        
        print("detect upload: ", datetime.datetime.now() - begin_time_detect)
        #detect upload:  0:00:02.337389 for 1k variants


        connection.commit()
        print("Transaction completed successfully ")

    except (Exception, psycopg2.Error) as error :        
        """
        Shall we add below 2 lines of code? Does it work?
        if connnection:
            connection.rollback()
        """
        print ("Error while connecting to PostgreSQL", error)
    finally:
        #closing database connection.
            if(connection):
                cursor.close()
                connection.close()
                print("PostgreSQL connection is closed")


##############################################



DATABASES = {
    'test':{
        'NAME': 'test_variants_db',
        'USER': 'geneticsuser',
        'PASSWORD': '***',
        'HOST': 'localhost',
        'PORT': '',
    },
    'production':{
        'NAME': 'variants_db',
        'USER': 'geneticsuser',
        'PASSWORD': '***',
        'HOST': 'localhost',
        'PORT': '',
    },
}

# choose the database to use
db = DATABASES['test']


#i_vcfs=["Final.vcf"]
#i_vcfs=["missense.vcf"]
print("Following samples will be uploaded: \n {}".format(i_vcf_path))

vcfTodb(i_vcf_path,db,orderID)

###################################################################################################################################################################

"""
Testet on 1000 variants and empty db. 

Transaction completed successfully 
PostgreSQL connection is closed
raw variants upload:  0:00:01.432317
{'user': 'geneticsuser', 'dbname': 'variants_db', 'host': 'localhost', 'port': '5432', 'tty': '', 'options': '', 'sslmode': 'prefer', 'sslcompression': '1', 'krbsrvname': 'postgres', 'target_session_attrs': 'any'} 

You are connected to -  ('PostgreSQL 10.15 (Ubuntu 10.15-0ubuntu0.18.04.1) on x86_64-pc-linux-gnu, compiled by gcc (Ubuntu 7.5.0-3ubuntu1~18.04) 7.5.0, 64-bit',) 

vep upload:  0:00:08.446516
Transaction completed successfully 
PostgreSQL connection is closed
{'user': 'geneticsuser', 'dbname': 'variants_db', 'host': 'localhost', 'port': '5432', 'tty': '', 'options': '', 'sslmode': 'prefer', 'sslcompression': '1', 'krbsrvname': 'postgres', 'target_session_attrs': 'any'} 

You are connected to -  ('PostgreSQL 10.15 (Ubuntu 10.15-0ubuntu0.18.04.1) on x86_64-pc-linux-gnu, compiled by gcc (Ubuntu 7.5.0-3ubuntu1~18.04) 7.5.0, 64-bit',) 

se upload:  0:00:02.268176
Transaction completed successfully 
PostgreSQL connection is closed
{'user': 'geneticsuser', 'dbname': 'variants_db', 'host': 'localhost', 'port': '5432', 'tty': '', 'options': '', 'sslmode': 'prefer', 'sslcompression': '1', 'krbsrvname': 'postgres', 'target_session_attrs': 'any'} 

You are connected to -  ('PostgreSQL 10.15 (Ubuntu 10.15-0ubuntu0.18.04.1) on x86_64-pc-linux-gnu, compiled by gcc (Ubuntu 7.5.0-3ubuntu1~18.04) 7.5.0, 64-bit',) 

pop upload:  0:00:00.148213
Transaction completed successfully 
PostgreSQL connection is closed
{'user': 'geneticsuser', 'dbname': 'variants_db', 'host': 'localhost', 'port': '5432', 'tty': '', 'options': '', 'sslmode': 'prefer', 'sslcompression': '1', 'krbsrvname': 'postgres', 'target_session_attrs': 'any'} 

You are connected to -  ('PostgreSQL 10.15 (Ubuntu 10.15-0ubuntu0.18.04.1) on x86_64-pc-linux-gnu, compiled by gcc (Ubuntu 7.5.0-3ubuntu1~18.04) 7.5.0, 64-bit',) 

detect upload:  0:00:01.746906
Transaction completed successfully 
PostgreSQL connection is closed
[Finished in 17.4s]



2nd time the same:

Transaction completed successfully 
PostgreSQL connection is closed
raw variants upload:  0:00:02.882016
{'user': 'geneticsuser', 'dbname': 'variants_db', 'host': 'localhost', 'port': '5432', 'tty': '', 'options': '', 'sslmode': 'prefer', 'sslcompression': '1', 'krbsrvname': 'postgres', 'target_session_attrs': 'any'} 

You are connected to -  ('PostgreSQL 10.15 (Ubuntu 10.15-0ubuntu0.18.04.1) on x86_64-pc-linux-gnu, compiled by gcc (Ubuntu 7.5.0-3ubuntu1~18.04) 7.5.0, 64-bit',) 

vep upload:  0:00:12.387126
Transaction completed successfully 
PostgreSQL connection is closed
{'user': 'geneticsuser', 'dbname': 'variants_db', 'host': 'localhost', 'port': '5432', 'tty': '', 'options': '', 'sslmode': 'prefer', 'sslcompression': '1', 'krbsrvname': 'postgres', 'target_session_attrs': 'any'} 

You are connected to -  ('PostgreSQL 10.15 (Ubuntu 10.15-0ubuntu0.18.04.1) on x86_64-pc-linux-gnu, compiled by gcc (Ubuntu 7.5.0-3ubuntu1~18.04) 7.5.0, 64-bit',) 

se upload:  0:00:03.501990
Transaction completed successfully 
PostgreSQL connection is closed
{'user': 'geneticsuser', 'dbname': 'variants_db', 'host': 'localhost', 'port': '5432', 'tty': '', 'options': '', 'sslmode': 'prefer', 'sslcompression': '1', 'krbsrvname': 'postgres', 'target_session_attrs': 'any'} 

You are connected to -  ('PostgreSQL 10.15 (Ubuntu 10.15-0ubuntu0.18.04.1) on x86_64-pc-linux-gnu, compiled by gcc (Ubuntu 7.5.0-3ubuntu1~18.04) 7.5.0, 64-bit',) 

pop upload:  0:00:00.091386
Transaction completed successfully 
PostgreSQL connection is closed
{'user': 'geneticsuser', 'dbname': 'variants_db', 'host': 'localhost', 'port': '5432', 'tty': '', 'options': '', 'sslmode': 'prefer', 'sslcompression': '1', 'krbsrvname': 'postgres', 'target_session_attrs': 'any'} 

You are connected to -  ('PostgreSQL 10.15 (Ubuntu 10.15-0ubuntu0.18.04.1) on x86_64-pc-linux-gnu, compiled by gcc (Ubuntu 7.5.0-3ubuntu1~18.04) 7.5.0, 64-bit',) 

detect upload:  0:00:01.807185
Transaction completed successfully 
PostgreSQL connection is closed
[Finished in 24.3s]
"""