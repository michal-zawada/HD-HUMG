import pandas as pd
import psycopg2
import pandas.io.sql as sqlio
import time, datetime
import sys
import os
import argparse
import subprocess

#To do:
#in csvForPipeline if te satus of selected sample for pipeline is status=0, then it shouldn't be accepted - user should get second warning and correct the chosen samples

#home/zawada/testdata/dkfz_files/201119_ST-K00265_0352_BHK7VNBBXY/Undetermined_1/fastq



timestr=time.strftime("%d%m%Y_%H%M%S")
fname="samplesToRun_{}.csv".format(timestr)

config_path=os.path.join('***','samples_to_nf')
archive_path=os.path.join(config_path,'archive')
currentRun_path=os.path.join(config_path,'current_run')
csv_path=os.path.join(currentRun_path,fname)
all_runs_path=os.path.join('***','dkfz_files')



parser=argparse.ArgumentParser()
parser.add_argument('-r', '--run', action='store', dest='inputRun', required=True, 
                    help="use --run parameter with run id in quots, for example: --run '201112_ST-K00246_0395_BHK7V2BBXY'")
args=parser.parse_args()
sys.stdout.write('\nInput run: {}\n'.format(args.inputRun))


run_id=args.inputRun
#run_id="201119_ST-K00265_0352_BHK7VNBBXY" #change parser.add_argument required=True 
dir_path=os.path.join(all_runs_path, run_id)
pipeline_version='v0.1'

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
        print ( connection.get_dsn_parameters())

        # Print PostgreSQL version
        cursor.execute("SELECT version();")
        
        record = cursor.fetchone()

        print("You are connected to - ", record)

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
                print("PostgreSQL connection is closed\n")




def cycleCountChanger(x):
    if x==101:
        x=100
    return x 

def seqReadChanger(x):
    if x=="PAIRED":
        x="pe"
    return x 


def getRunInfoFromFile(run_id):
    
    #Checking if given run exists aready in the db.
    seq_runs_str="""SELECT run_id FROM seq_runs;"""
    df_seq_runs=dfFromSql(seq_runs_str,db)
    if run_id in df_seq_runs['run_id'].values.tolist():
        print("{} run_id exists already in the db. Uploading terminated, all changes discarded.\n".format(run_id))
        sys.exit()

    selected_run_path=os.path.join(all_runs_path, run_id, "{}_meta.tsv".format(run_id))

    df=None
    if os.path.exists(selected_run_path):
        df=pd.read_table(selected_run_path, header=0)
    else:
        print("{} doesn't exist! Uploading terminated.".format(selected_run_path))
        sys.exit()


    cols=['FASTQ_FILE','READ', 'CENTER_NAME', 'RUN_ID', 'READ_COUNT', 'CYCLE_COUNT', 'SAMPLE_NAME','INSTRUMENT_MODEL',
       'SEQUENCING_READ_TYPE', 'ILSE_NO','LIB_PREP_KIT', 'AVERAGE_LANE_Q30']
    df=df[cols]
    df = df[df['SAMPLE_NAME'].notna()]
    df['CENTER_NAME']=df['CENTER_NAME'].apply(lambda x: x.lower())
    df['INSTRUMENT_MODEL']=df['INSTRUMENT_MODEL'].apply(lambda x: x.lower())
    df['CYCLE_COUNT']=df['CYCLE_COUNT'].apply(cycleCountChanger)
    df['SEQUENCING_READ_TYPE']=df['SEQUENCING_READ_TYPE'].apply(seqReadChanger)

    #print(df.to_string())

    #Selecting data for seq_run table upload
    cols_seqRun=["CENTER_NAME","RUN_ID","CYCLE_COUNT","INSTRUMENT_MODEL",
                    "SEQUENCING_READ_TYPE","ILSE_NO","LIB_PREP_KIT","AVERAGE_LANE_Q30"]
    df_seqRunToUpload=df[cols_seqRun]
    df_seqRunToUpload.columns=["seq_facility","run_id", "cycles_n", "seq_machine", "pe_sr", "submission", "lib_prep", "seq_qual_q30"]
    df_seqRunToUpload=df_seqRunToUpload.drop_duplicates()
    run_info=df_seqRunToUpload.to_dict('index') #index has nothing to do with index patient :p

    #Selecting data for content_seq_run table upload
    cols_contentSeqRun=["RUN_ID", "READ_COUNT", "SAMPLE_NAME"]
    df_contentSeqRunToUpload=df[cols_contentSeqRun]
    df_contentSeqRunToUpload.columns=["run_id","read_count","order_id"]
    df_contentSeqRunToUpload=df_contentSeqRunToUpload.drop_duplicates()
    run_content=df_contentSeqRunToUpload.to_dict('index') #index has nothing to do with index patient :p
    #Lane ratio is missing in file from dkfz, so it's set here to 0,25 as default. 
    for run in run_content.keys():
        run_content[run]["lane_ratio"]="0,25"
    sys.stdout.write("\nLane ratio was set to 0,25 as default. If you want to change it use content_seq_run_upload.py script.\n")

    #Selecting data for renaming the files
    cols_filesRename=['FASTQ_FILE','SAMPLE_NAME']
    df_filesRename=df[cols_filesRename]
    df_filesRename.columns=["fastq_names","order_id"]

    print("\n",df_filesRename.to_string(),"\n")
    old_newFileNames=df_filesRename.to_dict('index')

    cols_csvForPipeline=['SAMPLE_NAME','READ','FASTQ_FILE']
    df_csvForPipeline=df[cols_csvForPipeline]
    df_csvForPipeline.columns=["order_id","read","fastq_names"]    

    #order_ids=df_filesRename["order_id"].values.tolist()



    return (run_info,run_content,old_newFileNames,df_csvForPipeline)


def uploadRunInfoFromFile(runInfoDict,sql_str):
    #Uploading the data to db
    try:
        connection = psycopg2.connect(user = db['USER'],
                                    password = db['PASSWORD'],
                                    host = db['HOST'],
                                    database = db['NAME'])
        cursor = connection.cursor()
        # Print PostgreSQL Connection properties
        print ( connection.get_dsn_parameters())
        # Print PostgreSQL version
        cursor.execute("SELECT version();")        
        record = cursor.fetchone()
        print("You are connected to - ", record)        

        for run in runInfoDict.keys():                     
            cursor.execute(sql_str.format(**runInfoDict[run]))


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
                print("PostgreSQL connection is closed\n")


def renameFastq(dir_path, fileNamesDir):
    os.system("echo 'Compare sizes of the files with file names. Old file names:'")
    #print("ls -lh {}/*/fastq/*fastq.gz | cut -d' ' -f5,9".format(dir_path))
    os.system("ls -lh {}/*/fastq/*fastq.gz | cut -d' ' -f5,9".format(dir_path))
    renamed_dirs=[]        
    for f_name in fileNamesDir.keys():        
        fastqName = fileNamesDir[f_name]['fastq_names']
        orderID = fileNamesDir[f_name]['order_id']
        prefix=fastqName.split("_")[0]
        sample_old_DirPath=os.path.join(dir_path,prefix)
        sample_new_DirPath=os.path.join(dir_path,orderID)
        old_fastq_path = os.path.join(dir_path,orderID,"fastq/*")
        new_fastq_path = os.path.join(dir_path,orderID,"fastq",fastqName.replace(prefix,orderID))

        if sample_old_DirPath not in renamed_dirs: 
            os.system("mv {0} {1}".format(sample_old_DirPath,sample_new_DirPath))
            renamed_dirs.append(sample_old_DirPath)
        os.system("rename 's/{}/{}/' {}".format(prefix,orderID,old_fastq_path))
        #2 below commands didn't work with os.system, they work in commandline
        #os.system("for file in {0} ; do mv $file ${{file//{1}/{2}}} ; done".format(old_fastq_path,prefix,orderID))
        #os.system("for file in %s ; do mv $file ${file//%s/%s} ; done" % (old_fastq_path,prefix,orderID))
        #sys.stdout.write("Renamed {0} to {1}.\n".format(os.path.join(sample_old_DirPath,"fastq",fastqName),new_fastq_path))

    os.system("echo 'New file names:'")
    os.system("ls -lh {}/*/fastq/*fastq.gz | cut -d' ' -f5,9".format(dir_path))
    

##########creating csv file for nextflow. Pipeline_run_id should be downloaded from db and included in csv
#csvForPipeline function check first for which sample_run_id (from content_seq_runs table, depending on seq_run_id and order_id) pipeline was already ran. If it was ran for current one, gives Warning.
#If the user decides to run the pipeline anyway again for the sample_run_id which was already ran, then the sample_run_id recieves new pipeline_run_id and current pipeline_exec_date
#After uploading sample_ids to pipeline_runs table and giving by db new pipeline_run_ids, the pipeline_run_ids are downloaded form the db and saved in csv file with order and file names for nextflow.
#If the sample_run_id has been run again for a sample_run_id which is already there, the user should consider deleting one of the pipeline_run_id with the sample to avoid counting variants from a patient twice.

######WARNING! if we run pipeline again for the same sample, the status needs to be different than "0" for the previous ran sample, so we won't select the old pipeline_run_id to run again.   

def csvForPipeline(run_id,df_fileNames):
    
    df_fileNames["new_fastq_names"]=df_fileNames["order_id"]+"_"+df_fileNames["fastq_names"].apply(lambda x:x.split("_")[-1])
    
    df_fileNames["read"]="read"+df_fileNames["read"].astype(str)
    print(df_fileNames.to_string())   
    df_fileNames=df_fileNames.pivot(index="order_id", columns="read",values="new_fastq_names").reset_index()

    orders_from_tsv=df_fileNames["order_id"].values.tolist()


    orders_in_pipeline_runs="""SELECT pipeline_runs.pipeline_run_id, pipeline_runs.status, pipeline_runs.sample_run_id,pipeline_runs.pipeline_exec_date,content_seq_runs.order_id FROM pipeline_runs LEFT JOIN content_seq_runs ON pipeline_runs.sample_run_id=content_seq_runs.sample_run_id;"""
    df_ranOrders=dfFromSql(orders_in_pipeline_runs,db)
    ranOrders=df_ranOrders.to_dict('index')

    for sample in ranOrders.keys():        
        if ranOrders[sample]["order_id"] in orders_from_tsv:
            sys.stdout.write("\nWARNING!!! {} already analyzed by pipeline on {}, status {}!".format(ranOrders[sample]["order_id"],ranOrders[sample]["pipeline_exec_date"],ranOrders[sample]["status"]))    

    selectedOrders=[]
    sys.stdout.write("\nDo you want to select all below orders for next pipeline run? (y/n) ")
    for order in orders_from_tsv:
        sys.stdout.write("\n{}".format(order))

    user_input_all_samples=input("\n(y/n): ")
    #user_input_all_samples="y"
    if user_input_all_samples != "y":
        sys.stdout.write("\nSelect the samples for next pipeline run: (y/n) \n")        

        for order in orders_from_tsv:
            sys.stdout.write("{}".format(order))
            user_input_single_sample=input(" (y/n): ")            
            if user_input_single_sample == "y": 
                selectedOrders.append(order)
    else:
        selectedOrders=orders_from_tsv
        #selectedOrders=['P31671A', 'R32369A', 'S31468A']


    
    sys.stdout.write("\nYou have selected following samples: \n") 
    for selected in selectedOrders:
        sys.stdout.write("{}\n".format(selected))

    ordersForPipeline={}
    for i in range(0,len(selectedOrders)):
        ordersForPipeline[i]={}
        ordersForPipeline[i]["order_id"]=selectedOrders[i]
        ordersForPipeline[i]["run_id"]=run_id
        ordersForPipeline[i]["pipeline_version"]=pipeline_version
        ordersForPipeline[i]["status"]="0"

    print(ordersForPipeline)


    sqlStr_pipeline_runs="""INSERT INTO pipeline_runs (sample_run_id, pipeline_exec_date, pipeline_version, status) VALUES ((SELECT sample_run_id FROM content_seq_runs WHERE order_id='{order_id}' AND seq_run_id=(SELECT seq_run_id FROM seq_runs WHERE run_id='{run_id}')),DEFAULT, '{pipeline_version}', '{status}');"""
    uploadRunInfoFromFile(ordersForPipeline,sqlStr_pipeline_runs)

    ######WARNING! if we run pipeline again for the same sample, the status needs to be different than "0" for the previous ran sample, so we won't select the old pipeline_run_id to run again - in sqlStr_get_pipeline_runIds we are selecting pipeline_run_id where pipeline_runs.status='0'

    sqlStr_get_pipeline_runIds="""SELECT pipeline_runs.pipeline_run_id, content_seq_runs.order_id FROM pipeline_runs LEFT JOIN content_seq_runs ON pipeline_runs.sample_run_id=content_seq_runs.sample_run_id LEFT JOIN seq_runs ON content_seq_runs.seq_run_id=seq_runs.seq_run_id WHERE seq_runs.run_id='{}' AND pipeline_runs.status='0';""".format(run_id)
    df_pipelineIds=dfFromSql(sqlStr_get_pipeline_runIds,db)
    print(df_pipelineIds.to_string())



    df_fileNames=df_fileNames[df_fileNames['order_id'].isin(selectedOrders)]
    print(df_fileNames)
    df_forNF=pd.merge(df_fileNames, df_pipelineIds, on='order_id')
    print(df_forNF.to_string())
    df_forNF['read1']=df_forNF['read1'].apply(lambda x: os.path.join(dir_path,x.split("_")[0],'fastq',x))
    df_forNF['read2']=df_forNF['read2'].apply(lambda x: os.path.join(dir_path,x.split("_")[0],'fastq',x))

    # save above df to csv which will be used by nextflow, but before move all other files in folder for nextflow to archive

    os.system("mv {0}/* {1}".format(currentRun_path, archive_path))
    df_forNF.to_csv(csv_path,index=False, sep=",")

    #########################################################################################################################################
    #File with samples is saved in csv_path. Nextflow should be able to read it with the channel created there. 
    #Figure out how to pass pipeline_run_id from the last column in that csv file to parameter in nextflow, so it can be used in vcf_to_csv_with_sql5.py script as parameter to upload variants



sqlStr_seq_runs="""INSERT INTO seq_runs (run_id, seq_machine, pe_sr, cycles_n, seq_qual_q30, lib_prep, seq_facility, submission) SELECT '{run_id}', '{seq_machine}', '{pe_sr}', '{cycles_n}', '{seq_qual_q30}', '{lib_prep}', '{seq_facility}', '{submission}' WHERE NOT EXISTS (SELECT (run_id) FROM seq_runs WHERE run_id='{run_id}' LIMIT 1);"""
#Update has been removed from above sql command: UPDATE seq_runs SET seq_machine='{seq_machine}', pe_sr='{pe_sr}', cycles_n='{cycles_n}', seq_qual_q30='{seq_qual_q30}', lib_prep='{lib_prep}', seq_facility='{seq_facility}', submission='{submission}' WHERE run_id='{run_id}';
sqlStr_content_seq_runs="""INSERT INTO content_seq_runs (seq_run_id, order_id, read_count, lane_ratio) SELECT (SELECT seq_run_id FROM seq_runs WHERE run_id='{run_id}'), '{order_id}', '{read_count}', '{lane_ratio}' WHERE NOT EXISTS (SELECT (seq_run_id, order_id) FROM content_seq_runs WHERE seq_run_id=(SELECT seq_run_id FROM seq_runs WHERE run_id='{run_id}') AND order_id='{order_id}' LIMIT 1);"""
#Update has been removed from above sql command: UPDATE content_seq_runs SET read_count='{read_count}', lane_ratio='{lane_ratio}' WHERE seq_run_id=(SELECT seq_run_id FROM seq_runs WHERE run_id='{run_id}') AND order_id='{order_id}';

data_from_tsv=getRunInfoFromFile(run_id)


uploadRunInfoFromFile(data_from_tsv[0],sqlStr_seq_runs)
uploadRunInfoFromFile(data_from_tsv[1],sqlStr_content_seq_runs)
renameFastq(dir_path,data_from_tsv[2])
csvForPipeline(run_id,data_from_tsv[3])

##############write a script to rename the files







