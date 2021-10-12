import pandas as pd
import psycopg2
import pandas.io.sql as sqlio
import time, datetime
import sys
import os
import argparse

path=os.path.join('***','seq_runs_uploads')

timestr=time.strftime("%d%m%Y_%H%M%S")
fname="upload_seq_run_{}.xlsx".format(timestr)

fpath=os.path.join(path,fname)


parser=argparse.ArgumentParser()
parser.add_argument('-r', '--runs', action='store', dest='inputRuns', required=False, 
                    help="use --runs parameter with with comma-seperated sequencing run ids in quots, for example: --runs '201112_ST-K00246_0395_BHK7V2BBXY,190919_ST-K00352_0311_AHCKCCBBXY'")
args=parser.parse_args()
sys.stdout.write('\nInput run(s): {}\n'.format(args.inputRuns))

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


db_str_runs="""SELECT run_id,seq_machine,pe_sr,cycles_n,seq_qual_q30,lib_prep,seq_facility,submission FROM seq_runs;"""

df=dfFromSql(db_str_runs, db)

#if runs are given by user, then the table is created only with selected runs
if args.inputRuns:
    inputRuns=args.inputRuns.replace(" ","").split(",")       
    df_inputRuns=df[df['run_id'].isin(inputRuns)]
    for input_run in inputRuns:
        if input_run not in df['run_id'].values.tolist():
            sys.stdout.write("\n{0} run not found in DB. You can add it in the last row in the created table.\n".format(input_run))

    df=pd.concat([df.tail(3),df_inputRuns],sort=False) #I add 3 rows in the created excel file, so the user sees how should fill the table
else:
    df=df.tail(3)
df.drop_duplicates(inplace=True)


#df.fillna('', inplace=True)

df.to_excel(fpath, index=False)



run_charact = {'seq_machine':['miniseq','miseq','nextseq550','hiseq2500','hiseq4000','nextseq1000','nextseq2000','novaseq','unknown','other'],
                    'pe_sr':['pe','sr','unknown','other'],
                    'lib_prep':['agilent_ss','Agilent Low Input Exom-Seq Human v7','twist','illumina','unknown','other'],
                    'seq_facility':['dkfz','embl','uk_cardio','uk_patho','uk_hg','unknown','other']                    
                    }


sys.stdout.write("\nTable with run's data to edit has been generated in {}.\n".format(fpath))
sys.stdout.write("\nDon't change run_id!!!\n")
sys.stdout.write("\nPress 'y' after updating the table or 'n' if you want to discard the changes.\n")

user_input=input("\n(y/n): ")
if user_input != "y":
    print("Upload terminated, changes discarded.")
    sys.exit()


#checking if entered info are consistent
info_correct = False
while (info_correct == False):
    df=pd.read_excel(fpath)
    df.fillna('', inplace=True)
    seq_data=df.to_dict('index') #index has nothing to do with index patient :p        
    info_correct = True    
    for run in seq_data.keys():
        for feature in run_charact.keys():
            if seq_data[run][feature] not in run_charact[feature]:
                info_correct=False
                sys.stdout.write("\nChoosen feature {0} in row {1} is not available for {2}. Please choose one from available features: \n{3}.\n".format(seq_data[run][feature],(run+2),feature,', '.join(run_charact[feature])))

    if (info_correct == False):
        sys.stdout.write("\nPlease correct the file and press 'y' in order to upload the data or 'n' in order to terminate the process.\n")
        user_input=input("\n(y/n): ")
        if user_input != "y":
            print("Upload terminated, changes discarded.")
            sys.exit()
    else:
        sys.stdout.write("\nChanges accepted. Starting uploading...\n")



#Uploading the data to db
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


    

    for run in seq_data.keys():                     
        cursor.execute("""UPDATE seq_runs SET seq_machine='{seq_machine}', pe_sr='{pe_sr}', cycles_n='{cycles_n}', seq_qual_q30='{seq_qual_q30}', lib_prep='{lib_prep}', seq_facility='{seq_facility}', submission='{submission}' WHERE run_id='{run_id}';INSERT INTO seq_runs (run_id, seq_machine, pe_sr, cycles_n, seq_qual_q30, lib_prep, seq_facility, submission) SELECT '{run_id}', '{seq_machine}', '{pe_sr}', '{cycles_n}', '{seq_qual_q30}', '{lib_prep}', '{seq_facility}', '{submission}' WHERE NOT EXISTS (SELECT (run_id) FROM seq_runs WHERE run_id='{run_id}' LIMIT 1);""".format(**seq_data[run]))


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

