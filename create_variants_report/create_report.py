import psycopg2
import numpy as np
from openpyxl import load_workbook
import argparse
import pandas as pd
import pandas.io.sql as sqlio
import time




parser = argparse.ArgumentParser()
parser.add_argument('--input', action="store", dest='fileInput', required=True)
#parser.add_argument('--output', action="store", dest='fileOutput', default='output')
args = parser.parse_args()

print('Input sample: {}'.format(args.fileInput))
#print('Output file: {}'.format(args.fileOutput))

'''
This code can be used To parse more patients at once and generate reports for all of them at once:
parser = argparse.ArgumentParser()
parser.add_argument('-s', '--samples', action="store", dest='fileInput', required=True)
args = parser.parse_args()
print('Input sample(s): {}'.format(args.fileInput))
print([x.strip(' ') for x in args.fileInput.split(',')])
'''




def checkGenotype(x):
    if x=='0/1':
        return 'het'
    elif x=='1/1':
        return 'hom'

def getT1firstValue(x):
    if type(x)==str:
        return x.split(',')[0].split(' || ')[0]
    else:
        return x

def getVariantID(row):
    variant=[row['chromosome'],str(row['position']),row['ref'],row['alt']]
    variant_id = '-'.join(variant)
    return variant_id

def gnomADlink(row):
    if row['gnomAD_AC'] != 0 and len(row['variant_id'])<190:
        #variant=[row['chromosome'],str(row['position']),row['ref'],row['alt']]
        #variant_id = '-'.join(variant)       
        url='https://gnomad.broadinstitute.org/variant/{}?dataset=gnomad_r3'
        return '=HYPERLINK("%s", "%s")' % (url.format(row['variant_id']), row['variant_id'])       
    else:
        return row['variant_id']   

def addOMIMlink(x):
    url = "https://omim.org/entry/{}"
    if type(x)==str:
        return '=HYPERLINK("%s", "%s")' % (url.format(x), x)
    else:
        return x

def getPathogenic(x):   
    return (("Pathogenic" in x) or ("Likely_pathogenic" in x))

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


DATABASES = {
    'production':{
        'NAME': 'genetics',
        'USER': 'geneticsuser',
        'PASSWORD': '!s3q*n3Xt!',
        'HOST': 'localhost',
        'PORT': '',
    },
}

# choose the database to use
db = DATABASES['production']


#sample='A0210901'
#sample='S29752A'
sample=args.fileInput
#generating reports for many samples at once block code comes here. 
#for sample in args.samples:
    #sqlString.fromat() and run dfFromSql function

print('Generating report for {} sample...'.format(sample))

#get all variants from one patient
sqlString_patientVariants = ("SELECT * FROM varview_variant WHERE project_id='{}';").format(sample)
df = dfFromSql(sqlString_patientVariants,db)

#get all variants in house in order to calculate inhouse_MAF
sqlString_getAllVar = ("SELECT chromosome,position,ref,alt FROM varview_variant;")
df_inHouseVariants = dfFromSql(sqlString_getAllVar,db)

#get number of patients in house in order to calculate inhouse_MAF
sqlString_getNofPatientsIH = ("SELECT COUNT(*) FROM varview_project;")
nOfPatientsIH = dfFromSql(sqlString_getNofPatientsIH,db)['count'].iloc[0]

#create column variant_id. variant_id is a combiantion of chr, pos, ref and alt; for example:20_6543324_A_T
df_inHouseVariants['variant_id']=df_inHouseVariants.apply(lambda row: getVariantID(row), axis=1)

df_IHvariantsCounted=pd.DataFrame(df_inHouseVariants['variant_id'].value_counts()).reset_index()
df_IHvariantsCounted.columns=['variant_id','inhouse_counts']

df_IHvariantsCounted['inhouse_vaf']=df_IHvariantsCounted['inhouse_counts'].apply(lambda x: x/nOfPatientsIH)

#df_sortedCols = df[cols] 

df['gene']=df['snpeff_gene_symbol'].apply(getT1firstValue)

df['gnomAD_AF']=df['gnomAD_AF'].apply(getT1firstValue)

df['gnomAD_AC']=df['gnomAD_AC'].apply(getT1firstValue)

df['impact_T1']=df['snpeff_putative_impact']

#in column impact_T1 we use snpeff_putative_impact as first. If there's NA value it's filled with vep_impact value
df['impact_T1'].fillna(df['vep_impact'], inplace=True)

#in impact_T1 column we will have only one, the first one value 
df['impact_T1']=df['impact_T1'].apply(getT1firstValue)

df['index_genotype_T1']=df['index_genotype'].apply(checkGenotype)

df['vep_consequence_T1']=df['vep_consequence'].apply(getT1firstValue)

df['snpeff_hgvs_c_T1']=df['snpeff_hgvs_c'].apply(getT1firstValue)

df['snpeff_hgvs_p_T1']=df['snpeff_hgvs_p'].apply(getT1firstValue)

df.replace(['NA'], np.nan, inplace=True)

df=df.astype({'gnomAD_AC':'float64','gnomAD_AF':'float64'})

df['gnomAD_AF'].fillna(0.0, inplace=True)

df['gnomAD_AC'].fillna(0.0, inplace=True)
df=df.astype({'gnomAD_AC':'int64'})

#df['variant_id']=df[['chromosome','position','ref','alt']].astype(str).agg('-'.join,axis=1)

df['variant_id']=df.apply(lambda row: getVariantID(row), axis=1)

df=pd.merge(df,df_IHvariantsCounted,how='left', on=['variant_id'])

df['variant_id']=df.apply(lambda row: gnomADlink(row), axis=1)





cols=['chromosome',
 'position',
 'ref',
 'alt',
 'index_VAF',
 'index_DP',
 'gene',
 'index_genotype_T1',
 'impact_T1',
 'vep_consequence_T1',
 'snpeff_hgvs_c_T1',
 'snpeff_hgvs_p_T1',
 'omim_number_gene',
 'omim_phenotypes',
 'gnomAD_AC',
 'gnomAD_AF',
 'inhouse_counts',
 'inhouse_vaf',
 'variant_id',
 'clinvar_clnsig',
 'clinvar_var_source',
 'CADD_phred',
 'SIFT_pred',
 'MutationTaster_pred',
 'snpeff_refseq',
 'index_genotype',
 'snpeff_gene_symbol',
 'hg19_pos_1_based',
 'snpeff_sequence_ontology', 
 'snpeff_putative_impact', 
 'vep_feature',
 'snpeff_transcript_biotype',
 'vep_feature_type',
 'vep_consequence',
 'vep_impact',
 'snpeff_hgvs_c',
 'snpeff_hgvs_p', 
 'snpeff_hgvs_n', 
 'MetaSVM_pred',
 'Polyphen2_HVAR_pred',
 'Polyphen2_HDIV_pred',
 'MetaLR_pred',
 'M_CAP_pred',
 'FATHMM_pred',
 'VEST4_score',
 'Aloft_pred',
 'rs_dbSNP',
 'REVEL_score',
 'REVEL_rankscore',
 'clinpred_score',
 'k1000Gp3_AC',
 'k1000Gp3_AF',
 'k1000Gp3_AFR_AC',
 'k1000Gp3_AFR_AF',
 'k1000Gp3_EUR_AC',
 'k1000Gp3_EUR_AF',
 'k1000Gp3_AMR_AC',
 'k1000Gp3_AMR_AF',
 'k1000Gp3_EAS_AC',
 'k1000Gp3_EAS_AF',
 'k1000Gp3_SAS_AC',
 'k1000Gp3_SAS_AF',
 'TWINSUK_AC',
 'TWINSUK_AF',
 'ALSPAC_AC',
 'ALSPAC_AF',
 'ESP6500_AA_AC',
 'ESP6500_AA_AF',
 'ESP6500_EA_AC',
 'ESP6500_EA_AF',
 'ExAC_Adj_AC',
 'ExAC_Adj_AF',
 'ExAC_AFR_AC',
 'ExAC_AFR_AF',
 'ExAC_AMR_AC',
 'ExAC_AMR_AF',
 'ExAC_EAS_AC',
 'ExAC_EAS_AF',
 'ExAC_FIN_AC',
 'ExAC_FIN_AF',
 'ExAC_NFE_AC',
 'ExAC_NFE_AF',
 'ExAC_SAS_AC',
 'ExAC_SAS_AF',
 'gnomAD_exomes_AFR_AC',
 'gnomAD_exomes_AFR_AF',
 'gnomAD_exomes_AMR_AC',
 'gnomAD_exomes_AMR_AF',
 'gnomAD_exomes_ASJ_AC',
 'gnomAD_exomes_ASJ_AF',
 'gnomAD_exomes_EAS_AC',
 'gnomAD_exomes_EAS_AF',
 'gnomAD_exomes_FIN_AC',
 'gnomAD_exomes_FIN_AF',
 'gnomAD_exomes_NFE_AC',
 'gnomAD_exomes_NFE_AF',
 'gnomAD_exomes_SAS_AC',
 'gnomAD_exomes_SAS_AF',
 'clinvar_review',
 'Ensembl_geneid',
 'Ensembl_transcriptid',
 'Ensembl_proteinid', 
 'mother_VAF',
 'mother_DP',
 'mother_genotype',
 'father_VAF',
 'father_DP',
 'father_genotype',
 'inheritance',
 'project_id',
 'omim_comments',
 'omim_gene_name',
 'omim_gene_symbols', 
 'id']

df_sortedCols = df[cols]

df_sortedCols['omim_number_gene']=df_sortedCols['omim_number_gene'].apply(addOMIMlink)

df_hasClinVar=df_sortedCols[df_sortedCols['clinvar_clnsig'].notna()]

df_patClinVar = df_hasClinVar[df_hasClinVar['clinvar_clnsig'].apply(lambda x:(("Pathogenic" in x) or ("Likely_pathogenic" in x)))]

df_onlyOMIM=df_sortedCols[df_sortedCols['omim_number_gene'].notna()]

df_OnlyOMIM_MAF005 = df_onlyOMIM[df_onlyOMIM['gnomAD_AF']<0.05]

df_hom=df_OnlyOMIM_MAF005[df_OnlyOMIM_MAF005['index_VAF']>0.7]

df_OnlyOMIM_MAF001 = df_onlyOMIM[df_onlyOMIM['gnomAD_AF']<0.01]

df_het=df_OnlyOMIM_MAF001[(df_OnlyOMIM_MAF001['index_VAF']<=0.7)&(df_OnlyOMIM_MAF001['index_VAF']>=0.3)]

(df_OnlyOMIM_MAF005['gene'].value_counts()>=27)

df_CH=df_OnlyOMIM_MAF001[df_OnlyOMIM_MAF001['index_VAF']>0.2]

df_CH=df_CH.groupby('gene').filter(lambda x: len(x) > 1)

df_possMosaic=df_OnlyOMIM_MAF001[(df_OnlyOMIM_MAF001['index_VAF']<=0.3)&(df_OnlyOMIM_MAF001['index_VAF']>=0.1)]

len(df_possMosaic)

IF_genes='APC,MYH11,ACTA2,TMEM43,DSP,PKP2,DSG2,DSC2,BRCA1,BRCA2,SCN5A,RYR2,LMNA,MYBPC3,COL3A1,GLA,APOB,LDLR,MYH7,TPM1,MYBPC3,PRKAG2,TNNI3,MYL3,MYL2,ACTC1,RET,PCSK9,BMPR1A,SMAD4,TNNT2,TP53,TGFBR1,TGFBR2,TGFBR1,TGFBR2,SMAD3,KCNQ1,KCNH2,SCN5A,MLH1,MSH2,MSH6,PMS2,RYR1,CACNA1S,FBN1,TGFBR1,MEN1,RET,MUTYH,NF2,OTC,SDHD,SDHAF2,SDHC,SDHB,STK11,MUTYH,PTEN,RB1,TSC1,TSC2,VHL,WT1,ATP7B'.split(',')

df_IF_MAF005=df_OnlyOMIM_MAF005[df_OnlyOMIM_MAF005['gene'].apply(lambda x: x in IF_genes)]

timestr = time.strftime("%d%m%Y_%H%M%S")
path=r'***/{}_report_{}.xlsx'.format(sample,timestr)
writer = pd.ExcelWriter(path, engine='openpyxl')

df_sortedCols.to_excel(writer, sheet_name='all_variants', index=False)
df_patClinVar.to_excel(writer, sheet_name='ClinVar_pat', index=False)
df_onlyOMIM.to_excel(writer, sheet_name='Only_with_OMIM_nr', index=False)
df_hom.to_excel(writer, sheet_name='hom_MAF_0.05', index=False)
df_het.to_excel(writer, sheet_name='het_MAF_0.01', index=False)
df_CH.to_excel(writer, sheet_name='CH_MAF_0.05', index=False)
df_possMosaic.to_excel(writer, sheet_name='poss_mosaic', index=False)
df_IF_MAF005.to_excel(writer, sheet_name='df_IF_MAF005', index=False)

writer.close()
print('Report for {} patient with {} variants has been created in {}'.format(sample,len(df_sortedCols),path))
