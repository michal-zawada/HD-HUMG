import pandas as pd
import psycopg2
import pandas.io.sql as sqlio
import time, datetime
import sys
import os
import argparse



#to-do:
#check if added order_ids to the run exist already in database, if not, then warning and corection by user before upload 




path=os.path.join('***','content_seq_runs_uploads')

timestr=time.strftime("%d%m%Y_%H%M%S")
fname="upload_content_seq_runs_{}.xlsx".format(timestr)

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
        'PASSWORD': '!s3q*n3Xt!',#'!s3q*n3Xt!',
        'HOST': 'localhost',
        'PORT': '',
    },
    'production':{
        'NAME': 'variants_db',
        'USER': 'geneticsuser',
        'PASSWORD': '!s3q*n3Xt!',#'!s3q*n3Xt!',
        'HOST': 'localhost',
        'PORT': '',
    },
}

# choose the database to use
db = DATABASES['test']

string_content="""SELECT seq_runs.run_id, content_seq_runs.order_id, content_seq_runs.read_count, content_seq_runs.lane_ratio FROM seq_runs LEFT JOIN content_seq_runs ON seq_runs.seq_run_id=content_seq_runs.seq_run_id"""

df_content=dfFromSql(string_content, db)

print(df_content.to_string())

if args.inputRuns:
    inputRuns=args.inputRuns.replace(" ","").split(",")       
    df_inputRuns=df_content[df_content['run_id'].isin(inputRuns)]
    for input_run in inputRuns:
        if input_run not in df_content['run_id'].values.tolist():
            sys.stdout.write("\n{0} run not found in DB. You can add it in the last row in the created table.\n".format(input_run))
    df_content=pd.concat([df_content.tail(3),df_inputRuns],sort=False) #I add 3 rows in the created excel file, so the user sees how should fill the table
else:
    df_content=df_content.tail(3)

df_content.drop_duplicates(inplace=True)

df_content.to_excel(fpath, index=False)

sys.stdout.write("\nTable with patient's data to edit has been generated in {}.\n".format(fpath))
sys.stdout.write("\nrun_id and order_id will not be changed!!!\n")
sys.stdout.write("\nPress 'y' after updating the table or 'n' if you want to discard the changes.\n")


user_input=input("\n(y/n): ")
if user_input != "y":
    print("Upload terminated, changes discarded.")
    sys.exit()



df=pd.read_excel(fpath)
df.fillna('', inplace=True)
seq_runs=df.to_dict('index') #index has nothing to do with index patient :p


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


    

    for run in seq_runs.keys():                     
        cursor.execute("""UPDATE content_seq_runs SET read_count='{read_count}', lane_ratio='{lane_ratio}' WHERE seq_run_id=(SELECT seq_run_id FROM seq_runs WHERE run_id='{run_id}') AND order_id='{order_id}';INSERT INTO content_seq_runs (seq_run_id, order_id, read_count, lane_ratio) SELECT (SELECT seq_run_id FROM seq_runs WHERE run_id='{run_id}'), '{order_id}', '{read_count}', '{lane_ratio}' WHERE NOT EXISTS (SELECT (seq_run_id, order_id) FROM content_seq_runs WHERE seq_run_id=(SELECT seq_run_id FROM seq_runs WHERE run_id='{run_id}') AND order_id='{order_id}' LIMIT 1);""".format(**seq_runs[run]))


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

