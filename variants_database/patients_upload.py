import pandas as pd
import psycopg2
import pandas.io.sql as sqlio
import time, datetime
import sys
import os
import argparse

"""
to do:
1. possibility to delete a patient using this script but only if the patient was not sequenced 
2. check if dob has correct date format: yyyy-mm-dd

"""


path=os.path.join('/media','data','WES','db_uploads','patients_uploads')

timestr=time.strftime("%d%m%Y_%H%M%S")
fname="update_patient_{}.xlsx".format(timestr)

fpath=os.path.join(path,fname)



parser=argparse.ArgumentParser()
parser.add_argument('-p', '--patients', action='store', dest='inputPatients', required=False, 
                    help="use --patients parameter with with comma-seperated patient_ids in quots, for example: --patients 'patient_id1,patient_id2'")
args=parser.parse_args()
sys.stdout.write('\nInput patient(s): {}\n'.format(args.inputPatients))


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


db_str_patients="""SELECT patients.patient_id, orders.order_id, families.family_id, orders.analysis_type, patients.aff_status, patients.dob, patients.gender, patients.idx_status,families.consanguineous_stat, families.region, hpos.hpo_number FROM families LEFT JOIN patients ON families.family_id=patients.family_id LEFT JOIN orders ON orders.patient_id=patients.patient_id LEFT JOIN phenotypes ON patients.patient_id=phenotypes.patient_id LEFT JOIN hpos ON phenotypes.hpo_id=hpos.hpo_id;"""
string_hpos="""SELECT hpo_number FROM hpos;"""

#getting dataframe with patients, hpos, analyisis etc.
df=dfFromSql(db_str_patients, db)

#getting dataframe with all available hpos from variants_db 
df_hpos=dfFromSql(string_hpos, db)
hpos_list=df_hpos['hpo_number'].values.tolist()
hpos_list.append('')




if len(df)>0:
    df.fillna('', inplace=True)
    df['hpo_number']=df[['patient_id', 'family_id', 'order_id', 'analysis_type', 'aff_status',
           'dob', 'gender', 'idx_status', 'consanguineous_stat', 'region',
           'hpo_number']].groupby(['order_id'])['hpo_number'].transform(lambda x: ','.join(x)).apply(lambda x: x.strip(","))
    df.drop_duplicates(inplace=True)


#if patients are given by user, then the table is created only with selected patients
if args.inputPatients:
    inputPatients=args.inputPatients.replace(" ","").split(",") 
    try:
        inputPatients=list(map(int, inputPatients))
    except:
        print("patient_id must be a number! Uploading terminated.")
        sys.exit()    
    df_inputPat=df[df['patient_id'].isin(inputPatients)]
    for input_pat in inputPatients:
        if input_pat not in df['patient_id'].values.tolist():
            sys.stdout.write("\n{0} patient not found in DB. You can add it in the last row in the created table.\n".format(input_pat))  
    df=pd.concat([df.tail(3),df_inputPat],sort=False)
else:
    df=df.tail(3) #I keep just 3 entries from db, so user can see how the table is filled
df.drop_duplicates(inplace=True)


#saving all ids in a dictionary to be sure that they are not changed after user update
ids_before_update={'patient_id':df['patient_id'].values.tolist(),'family_id':df['family_id'].values.tolist(),'order_id':df['order_id'].values.tolist()}

print("Ids before: \n",ids_before_update)

df.to_excel(fpath, index=False)

person_charact = {'consanguineous_stat':['yes','no','unknown'],
                    'aff_status':['affected','unaffected','unknown'],
                    'gender':['male','female','unknown'],
                    'idx_status':['index','mother','father','sister','brother','son','daughter','husband','wife','partner','grandmother','grandfather','cousin','uncle','aunt','unknown','other'],
                    'analysis_type':['wes_agi','wgs','other','unknown']
                    }
consistent_features=['hpo_number','dob','consanguineous_stat','aff_status','gender','idx_status']


sys.stdout.write("\nTable with patient's data to edit has been generated in {}.\n".format(fpath))
sys.stdout.write("\nDon't change patient_id, order_id and famliy_id!!!\n")
sys.stdout.write("\nPress 'y' after updating the table or 'n' if you want to discard the changes.\n")


user_input=input("\n(y/n): ")
if user_input != "y":
    print("Upload terminated, changes discarded.")
    sys.exit()



#checking if entered info are consistent
info_correct = False
info_consistent = False
while (info_correct == False) or (info_consistent == False):
    df=pd.read_excel(fpath)
    df.fillna('', inplace=True)
    patients=df.to_dict('index') #index has nothing to do with index patient :p    
    patient_2x=[]
    info_correct = True
    info_consistent = True
    for person in patients.keys():    
        patients[person]['hpo_number']=patients[person]['hpo_number'].replace(" ","").strip(",")
        #checking if info for patient which has 2 or more order_id (2 analysis) is the same in all rows except of order_id and analysis type
        if patients[person]['patient_id'] in patient_2x:
            for person2 in patients.keys():
                if patients[person]['patient_id']==patients[person2]['patient_id']:                    
                    for person_feature in consistent_features:
                        if patients[person][person_feature]!=patients[person2][person_feature]: 
                            info_consistent = False                         
                            sys.stdout.write("\nInconsistent {0} in the same patient_id, row {1} and {2}.\n".format(person_feature,(person+2),(person2+2)))                           
                            
        else:
            patient_2x.append(patients[person]['patient_id'])
        #checking if given info is valid in db, because all info need to be consistent in db. Acceptable features are given in person_charact dictionary
        for feature in person_charact.keys():
            if patients[person][feature] not in person_charact[feature]:
                info_correct=False
                sys.stdout.write("\nChoosen feature {0} in row {1} is not available for {2}. Please choose one from available features: \n{3}.\n".format(patients[person][feature],(person+2),feature,', '.join(person_charact[feature])))
        #checking if given hpos are present in db
        for hpo in patients[person]['hpo_number'].split(","):            
            if hpo not in hpos_list:
                info_correct=False
                sys.stdout.write("\n {0} hpo given in row {1} is not present in set of hpos in db. Please correct it or replace with other.\n".format(hpo,(person+2)))
        #order_id and family_id are necessary to proceed. Checking if they are there
        for feature_id in ['family_id','order_id']:
            if patients[person][feature_id]=='':
                info_correct=False
                sys.stdout.write("\n {0} in row {1} cannot be empty. Please correct it.\n".format(feature_id,(person+2)))
    #check if ids in rows before update are the same like after - order_id, patient_id, family_id
    for id_b_update in ids_before_update.keys():
        for i in range(len(ids_before_update[id_b_update])):
            if ids_before_update[id_b_update][i]!=df[id_b_update].values.tolist()[i]:
                info_consistent = False
                sys.stdout.write("\nForbidden action detected: {0} {1} has been changed. Terminating update. All changes will be discarded.\n".format(id_b_update, ids_before_update[id_b_update][i]))
                sys.exit()

    if (info_correct == False) or (info_consistent==False):
        sys.stdout.write("\nPlease correct the file and press 'y' in order to upload the data or 'n' in order to terminate the process.\n")
        user_input=input("\n(y/n): ")
        if user_input != "y":
            print("Upload terminated, changes discarded.")
            sys.exit()
    else:
        sys.stdout.write("\nChanges accepted. Starting uploading...\n")


print(df.to_string())

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


    
    #here we are going to update patient's info and and if patient doesn't exists we will create it. 
    #patient_id, family_id and order_id cannot be updated (changed) this way, because we need to konw which data to update, but new patient with new order_id can be created. 
    for person in patients.keys():             
        cursor.execute("""UPDATE families SET consanguineous_stat='{consanguineous_stat}', region='{region}' WHERE family_id={family_id};INSERT INTO families (family_id, consanguineous_stat, region) SELECT {family_id}, '{consanguineous_stat}', '{region}' WHERE NOT EXISTS (SELECT (family_id) FROM families WHERE family_id={family_id} LIMIT 1);""".format(**patients[person]))
        cursor.execute("""UPDATE patients SET aff_status='{aff_status}', dob='{dob}', gender='{gender}', idx_status='{idx_status}' WHERE patient_id={patient_id}; INSERT INTO patients (patient_id, family_id, aff_status, dob, gender, idx_status) SELECT {patient_id}, {family_id}, '{aff_status}', '{dob}', '{gender}', '{idx_status}'  WHERE NOT EXISTS (SELECT (patient_id) FROM patients WHERE patient_id={patient_id} LIMIT 1);""".format(**patients[person]))
        cursor.execute("""UPDATE orders SET analysis_type='{analysis_type}' WHERE order_id='{order_id}'; INSERT INTO orders (order_id, patient_id, analysis_type) SELECT '{order_id}', {patient_id}, '{analysis_type}' WHERE NOT EXISTS (SELECT (order_id) FROM orders WHERE order_id='{order_id}' LIMIT 1);""".format(**patients[person]))
    
    for person in patients.keys():         
        #we need to delete first all hpos from patients, because maybe the list of hpos was changed 
        cursor.execute("""DELETE FROM phenotypes WHERE patient_id = {patient_id}""".format(**patients[person]))
        #The Phe patient's list of HPOs could have been changed. It#s n:m connsction so we need to delete all hpo entries for the patients and enter new one.   
        for hpo in patients[person]['hpo_number'].split(","):
            #we need this if statement here to avoid entering empty string to phenotype table
            if hpo!='':
                print('{patient_id},{hpo_number}'.format(patient_id=patients[person]['patient_id'],hpo_number=hpo))                
                cursor.execute("""INSERT INTO phenotypes (patient_id, hpo_id) VALUES ({patient_id},(SELECT hpo_id FROM hpos WHERE hpo_number='{hpo_number}'));""".format(patient_id=patients[person]['patient_id'],hpo_number=hpo))


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



