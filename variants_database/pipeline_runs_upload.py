import pandas as pd
import psycopg2
import pandas.io.sql as sqlio
import time, datetime
import sys
import os
import argparse



#to-do:
#generating fine_names.csv in a similar way like seq_run_upload, or wven automaticaly base on given run - old file names from the file location and new file names form order_id in db.
#checking if given names are present in db for given run BEFORE renaming the files 
#in pipRunsUpload if te satus of selected sample for pipeline is status=0, then it shouldn't be accepted - user should get second warning and correct the chosen samples

parser=argparse.ArgumentParser()
parser.add_argument('-r', '--run', action='store', dest='inputRuns', required=False, 
                    help="use --run parameter with run id in quots, for example: --run '201112_ST-K00246_0395_BHK7V2BBXY'")
args=parser.parse_args()
sys.stdout.write('\nInput run(s): {}\n'.format(args.inputRuns))


pipeline_version='v0.1'
#run_id="190919_ST-K00352_0311_AHCKCCBBXY"
run_id="201119_ST-K00265_0352_BHK7VNBBXY"
dir_path=os.path.join('/media','data','zawada','fastq_test_data','dkfz_files', run_id)

timestr=time.strftime("%d%m%Y_%H%M%S")
fname="samplesToRun_{}.csv".format(timestr)

config_path=os.path.join('***','samples_to_nf')
archive_path=os.path.join(config_path,'archive')
currentRun_path=os.path.join(config_path,'current_run')
csv_path=os.path.join(currentRun_path,fname)

fnames_csv_path=os.path.join(config_path,'manual','file_names_current_run','file_names.csv')

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
        print ("Error while connecting tdir_patho PostgreSQL", error)
    finally:
        #closing database connection.
            if(connection):
                cursor.close()
                connection.close()
                print("PostgreSQL connection is closed")


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

def getFileNamesFromFile(f_path):
    
    df=None
    if os.path.exists(f_path):
        df_fileNames=pd.read_csv(f_path, header=0)
    else:
        print("{} doesn't exist! Uploading terminated.".format(f_path))
        sys.exit()
 
    return df_fileNames


def renameFastq(dir_path, fileNames_df):

    fileNamesDir=fileNames_df.to_dict('index') #index has nothing to do with index patient

    os.system("echo 'Compare sizes of the files with file names. Old file names:'")
    #print("ls -lh {}/*/fastq/*fastq.gz | cut -d' ' -f5,9".format(dir_path))
    os.system("ls -lh {}/*/fastq/*fastq.gz | cut -d' ' -f5,10".format(dir_path))
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
    os.system("ls -lh {}/*/fastq/*fastq.gz | cut -d' ' -f5,10".format(dir_path))

#pipRunsUpload function is similar to csvForPipeline but it has different arguments - only run_id, no df with orders. order_ids are taken from db base on run_id

######WARNING! if we run pipeline again for the same sample, the status needs to be different than "0" for the previous ran sample, so we won't select the old pipeline_run_id to run again. 

def pipRunsUpload(run_id,df_fileNames):

    df_fileNames["new_fastq_names"]=df_fileNames["order_id"]+"_"+df_fileNames["fastq_names"].apply(lambda x:x.split("_")[-1])    
    df_fileNames["read"]="read"+df_fileNames["read"].astype(str)  
    df_fileNames=df_fileNames.pivot(index="order_id", columns="read",values="new_fastq_names").reset_index() #df_fileNames come from csv file

    orders_in_content_seq_runs="""SELECT  content_seq_runs.order_id, seq_runs.run_id FROM content_seq_runs LEFT JOIN seq_runs ON content_seq_runs.seq_run_id=seq_runs.seq_run_id WHERE seq_runs.run_id='{}';""".format(run_id)
    df_seqOrders=dfFromSql(orders_in_content_seq_runs,db) #df_seqOrders is a df with order_ids only from selected run_id

    #checking if order_ids in csv file with new names (fnames_csv_path) are in given run_id in db.
    for order_id in df_fileNames['order_id'].values.tolist():
        if order_id not in df_seqOrders['order_id'].values.tolist():
            sys.stdout.write("\norder_id {} in csv file with order_ids to rename is not in db for given run {}".format(order_id,run_id))
            sys.stdout.write("\nPlease correct {} order_id in the file or in db!!".format(order_id))
            sys.stdout.write("\nfastq file names for given run {} has been already changed according to given order_ids. Please check them!!\nProcess terminated\n".format(run_id))
            sys.exit()


    orders_in_pipeline_runs="""SELECT pipeline_runs.pipeline_run_id, pipeline_runs.status, pipeline_runs.sample_run_id,pipeline_runs.pipeline_exec_date,content_seq_runs.order_id, seq_runs.run_id FROM pipeline_runs LEFT JOIN content_seq_runs ON pipeline_runs.sample_run_id=content_seq_runs.sample_run_id LEFT JOIN seq_runs ON content_seq_runs.seq_run_id=seq_runs.seq_run_id;"""
    df_ranOrders=dfFromSql(orders_in_pipeline_runs,db)
    ranOrders=df_ranOrders.to_dict('index') #ranOrders comes from pipeline_runs table in db. They are all order_ids, for which pipeline was already ran

    orders_in_selected_run=df_seqOrders['order_id'].values.tolist() #orders_in_selected_run come from db, base on run_id

    #checking if orders from selected run were already ran in pipeline
    for sample in ranOrders.keys():        
        if ranOrders[sample]["order_id"] in orders_in_selected_run:
            sys.stdout.write("\nWARNING!!! {} already analyzed by pipeline on {}, status {}!".format(ranOrders[sample]["order_id"],ranOrders[sample]["pipeline_exec_date"],ranOrders[sample]["status"]))    

    selectedOrders=[]
    sys.stdout.write("\nDo you want to select all below orders for next pipeline run? (y/n) ")
    for order in orders_in_selected_run:
        sys.stdout.write("\n{}".format(order))

    user_input_all_samples=input("\n(y/n): ")
    #user_input_all_samples="y"
    if user_input_all_samples != "y":
        sys.stdout.write("\nSelect the samples for next pipeline run: (y/n) \n")        

        for order in orders_in_selected_run:
            sys.stdout.write("{}".format(order))
            user_input_single_sample=input(" (y/n): ")            
            if user_input_single_sample == "y": 
                selectedOrders.append(order)
    else:
        selectedOrders=orders_in_selected_run
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

    sqlStr_pipeline_runs="""INSERT INTO pipeline_runs (sample_run_id, pipeline_exec_date, pipeline_version, status) VALUES ((SELECT sample_run_id FROM content_seq_runs WHERE order_id='{order_id}' AND seq_run_id=(SELECT seq_run_id FROM seq_runs WHERE run_id='{run_id}')),DEFAULT, '{pipeline_version}', '{status}');"""
    uploadRunInfoFromFile(ordersForPipeline,sqlStr_pipeline_runs)

    ######WARNING! if we run pipeline again for the same sample, the status needs to be different than "0" for the previous ran sample, so we won't select the old pipeline_run_id to run again. 

    sqlStr_get_pipeline_runIds="""SELECT pipeline_runs.pipeline_run_id, content_seq_runs.order_id FROM pipeline_runs LEFT JOIN content_seq_runs ON pipeline_runs.sample_run_id=content_seq_runs.sample_run_id LEFT JOIN seq_runs ON content_seq_runs.seq_run_id=seq_runs.seq_run_id WHERE seq_runs.run_id='{}' AND pipeline_runs.status='0';""".format(run_id)
    df_pipelineIds=dfFromSql(sqlStr_get_pipeline_runIds,db)

    df_fileNames=df_fileNames[df_fileNames['order_id'].isin(selectedOrders)]
    print(df_fileNames)
    df_forNF=pd.merge(df_fileNames, df_pipelineIds, on='order_id')
    print(df_forNF.to_string())
    df_forNF['read1']=df_forNF['read1'].apply(lambda x: os.path.join(dir_path,x.split("_")[0],'fastq',x))
    df_forNF['read2']=df_forNF['read2'].apply(lambda x: os.path.join(dir_path,x.split("_")[0],'fastq',x))
    # save above df to csv which will be used by nextflow, but before move all other files in folder for nextflow to archive

    os.system("mv {0}/* {1}".format(currentRun_path, archive_path))
    df_forNF.to_csv(csv_path,index=False, sep=",")

df_file_names = getFileNamesFromFile(fnames_csv_path)
renameFastq(dir_path,df_file_names)
pipRunsUpload(run_id,df_file_names)