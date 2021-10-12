import os
import pandas as pd
import psycopg2
import pandas.io.sql as sqlio
from sqlalchemy import create_engine


#before runing this script, all values from 3 tables: pheno_genes, gene_hpo and hpos has to be deleted! 
#to-do: change the script in the way that will be updating the hpos, not tht we have to delete the old and upload new

def uploadInfoFromDict(infoDict,sql_str):
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

        for info in infoDict.keys():                     
            cursor.execute(sql_str.format(**infoDict[info]))


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

engine=create_engine('postgresql://geneticsuser:***@localhost/test_variants_db') ################## change test to production db

col_names=["hpo_number","hpo_label","gene_id","gene_symbol","Additional_Info","G-D_source","disease_ID"]

df=pd.read_table("phenotype_to_genes.txt", names=col_names, header=0)

df_uniqueHPOs=df[['hpo_number','hpo_label']].drop_duplicates().reset_index(drop=True)

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
    
    
    
    df_uniqueHPOs.to_sql('hpos',engine,if_exists='append',index=False)

    #df=sqlio.read_sql_query(sqlString, connection) 
    #df=pd.read_sql_query(sqlString, connection)
    


except (Exception, psycopg2.Error) as error :
    print ("Error while connecting to PostgreSQL", error)
finally:
    #closing database connection.
        if(connection):
            cursor.close()
            connection.close()
            print("PostgreSQL connection is closed")

df[['hpo_number','hpo_label']].drop_duplicates().reset_index(drop=True)

df_uniqueGenes=df[['gene_id','gene_symbol']].drop_duplicates().reset_index(drop=True)

df_uniqueGenes.columns=['gene_entrez_id','gene_symbol']

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
    
    
    
    df_uniqueGenes.to_sql('pheno_genes',engine,if_exists='append',index=False)

    #df=sqlio.read_sql_query(sqlString, connection) 
    #df=pd.read_sql_query(sqlString, connection)
    


except (Exception, psycopg2.Error) as error :
    print ("Error while connecting to PostgreSQL", error)
finally:
    #closing database connection.
        if(connection):
            cursor.close()
            connection.close()
            print("PostgreSQL connection is closed")

df_gene_hpo=df[['gene_id','hpo_number']].drop_duplicates().reset_index(drop=True)

gene_hpo_dict=df_gene_hpo.to_dict('index')

sqlStr_gene_hpo="""INSERT INTO gene_hpo (pheno_gene_id, hpo_id) SELECT (SELECT pheno_gene_id FROM pheno_genes WHERE gene_entrez_id = '{gene_id}'), (SELECT hpo_id FROM hpos WHERE hpo_number = '{hpo_number}');"""

#SELECT * FROM pheno_genes LEFT JOIN gene_hpo ON pheno_genes.pheno_gene_id=gene_hpo.pheno_gene_id LEFT JOIN hpos ON gene_hpo.hpo_id=hpos.hpo_id WHERE gene_symbol = 'SKI';

#HP:0004322,HP:0040064,HP:0002857,HP:0002650,HP:0001638,HP:0001388

uploadInfoFromDict(gene_hpo_dict,sqlStr_gene_hpo)