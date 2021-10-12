/*
su -l postgres
psql



*/
CREATE DATABASE test_variants_db;

\list 
\c variants_db
\dt

/*
Giving previliges to a user. We should consider which should be given, so far with 
bellow commands all previliges are given
 */
GRANT CONNECT ON DATABASE test_variants_db TO geneticsuser;
GRANT USAGE ON SCHEMA public TO geneticsuser;
GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA public TO geneticsuser;
GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA public TO geneticsuser;


GRANT CONNECT ON DATABASE variants_db TO geneticsuser;
GRANT USAGE ON SCHEMA public TO geneticsuser;
GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA public TO geneticsuser;
GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA public TO geneticsuser;

________________________

DROP TABLE IF EXISTS variants;
CREATE TABLE variants (
	variant_id bigserial PRIMARY KEY,
	chrom VARCHAR (255) NOT NULL,
	pos int NOT NULL,
	ref VARCHAR (255) NOT NULL,
	alt VARCHAR (255) NOT NULL,
	UNIQUE (chrom,pos,ref,alt)
);

ALTER TABLE variants DROP CONSTRAINT variants_variant_id_chrom_pos_ref_alt_key;
ALTER TABLE variants ADD UNIQUE (chrom,pos,ref,alt);

INSERT INTO variants (chrom,pos,ref,alt) VALUES ('1',5432165,'A','T') ON CONFLICT ON CONSTRAINT "variants_chrom_pos_ref_alt_key" DO NOTHING; 
/*Below command gives error if variant exists*/
INSERT INTO variants (chrom,pos,ref,alt) VALUES ('1',5432165,'A','T');

INSERT INTO variants (chrom,pos,ref,alt)
SELECT '1' AS chrom, 5432164 AS pos, 'A' AS ref, 'T' AS alt FROM variants
WHERE NOT EXISTS(
            SELECT (chrom,pos,ref,alt) FROM variants WHERE chrom='1' AND pos=5432164 AND ref='A' AND alt='T'
    )
LIMIT 1;

/*
With INSERT ON CONFLICT ON CONSTRAINT the variant_id is increased even if variant is not inserted.
With INSERT with SELECT the variant_id is not increased if the variant is not inserted, but not sure 
if the table is big selecting will take a lot of time for every variant. So far INSERT with SELECT 
is implemented in vcf_to_csv_sql.py script.
*/



DROP TABLE IF EXISTS se_annotations;
CREATE TABLE se_annotations (
	se_annotation_id bigserial PRIMARY KEY,
	variant_id bigint NOT NULL,
	allele_se VARCHAR (255),
	annotation_se VARCHAR (255),
	annotation_imp_se VARCHAR (255),
	gene_name VARCHAR (255),
	gene_id VARCHAR (255),
	feature_type_se VARCHAR (255),
	feature_id_se VARCHAR (255),
	transcript_biotype VARCHAR (255),
	rank_pos_se VARCHAR (255),
	rank_length_se VARCHAR (255),
	hgvs_c VARCHAR (255),
	hgvs_p VARCHAR (255),
	cdna_pos VARCHAR (255),
	cdna_length VARCHAR (255),
	cds_pos VARCHAR (255),
	cds_length VARCHAR (255),
	aa_pos VARCHAR (255),
	aa_length VARCHAR (255),
	distance VARCHAR (255),
	errors VARCHAR (255),
	se_update date,
	UNIQUE (variant_id, feature_id_se, se_update),
	FOREIGN KEY (variant_id)
		REFERENCES variants (variant_id) ON UPDATE CASCADE ON DELETE CASCADE
);

CREATE TABLE vep_annotations (
	vep_annotation_id bigserial PRIMARY KEY,
	variant_id bigint NOT NULL,
	allele_vep VARCHAR (255),
	consequence_vep VARCHAR (255),
	feature_type_vep VARCHAR (255),
	feature_vep VARCHAR (255),
	impact_vep VARCHAR (255),
	symbol_vep VARCHAR (255),
	vep_update date,
	UNIQUE (variant_id, feature_vep,vep_update),
	FOREIGN KEY (variant_id)
		REFERENCES variants (variant_id) ON UPDATE CASCADE ON DELETE CASCADE
);

ALTER TABLE vep_annotations ADD UNIQUE (variant_id,feature_vep);
ALTER TABLE vep_annotations ALTER COLUMN feature_vep SET NOT NULL;



INSERT INTO vep_annotations (variant_id, allele_vep, consequence_vep, feature_type_vep, feature_vep, impact_vep, symbol_vep) VALUES ((SELECT variant_id FROM variants WHERE chrom='1' AND pos=14464 AND ref='A' AND alt='T'), 'T', 'downstream_gene_variant', 'Transcript', 'ENST00000450305', 'MODIFIER', 'DDX11L1');
['1', 14464, 'A', 'T', 'T', 'downstream_gene_variant', 'Transcript', 'ENST00000450305', 'MODIFIER', 'DDX11L1']

/*Similar situation like for variants table.
With INSERT ON CONFLICT ON CONSTRAINT (below line) the vep_annotation_id for the next transcript is increased even if transcript is not inserted.*/
INSERT INTO vep_annotations (variant_id, allele_vep, consequence_vep, feature_type_vep, feature_vep, impact_vep, symbol_vep) VALUES ((SELECT variant_id FROM variants WHERE chrom='1' AND pos=14464 AND ref='A' AND alt='T'), 'T', 'downstream_gene_variant', 'Transcript', 'ENST00000450306', 'MODIFIER', 'DDX11L1') ON CONFLICT ON CONSTRAINT "vep_annotations_variant_id_feature_vep_key" DO UPDATE SET allele_vep='T', consequence_vep='downstream_gene_variant', feature_type_vep='Transcript', impact_vep='MODIFIER', symbol_vep='DDX11L3';


UPDATE vep_annotations SET allele_vep='T', consequence_vep='downstream_gene_variant', feature_type_vep='Transcript', impact_vep='MODIFIER', symbol_vep='DDX11L3' WHERE vep_annotation_id = 15;
#Below command works, but it's super long. We have to always select variants from variants table to get variant_id. Maybe since we know variant_id (because raw variants are uploaded as first) we can write it in vcf file in ID column?
INSERT INTO vep_annotations (variant_id, allele_vep, consequence_vep, feature_type_vep, feature_vep, impact_vep, symbol_vep) SELECT (SELECT variant_id FROM variants WHERE chrom='1' AND pos=14464 AND ref='A' AND alt='T'),'T','downstream_gene_variant', 'Transcript', 'ENST00000450308', 'MODIFIER', 'DDX11L1' WHERE NOT EXISTS (SELECT (variant_id, feature_vep) FROM vep_annotations WHERE variant_id=(SELECT (variant_id) FROM variants WHERE chrom='1' AND pos=14464 AND ref='A' AND alt='T' LIMIT 1) AND feature_vep='ENST00000450308' LIMIT 1);

#Implemented solution is uploading first raw varints to variants table and then downloading variant_ids. Saving variant_ids and raw variants in another vcf and annotating the orgiginal vcf with variant_ids. 
#For uploading annotations (here vep_annotation) is used annotated vcf with variant_ids from our database'- below comand 
#"""UPDATE vep_annotations SET allele_vep='{1}', consequence_vep='{2}', feature_type_vep='{3}', impact_vep='{5}', symbol_vep='{6}' WHERE variant_id={0} AND feature_vep='{4}';INSERT INTO vep_annotations (variant_id, allele_vep, consequence_vep, feature_type_vep, feature_vep, impact_vep, symbol_vep) SELECT {0},'{1}','{2}', '{3}', '{4}', '{5}', '{6}' WHERE NOT EXISTS (SELECT (variant_id, feature_vep) FROM vep_annotations WHERE variant_id={0} AND feature_vep='{4}' LIMIT 1);""".format(transcript[0],transcript[1],transcript[2],transcript[3],transcript[4],transcript[5],transcript[6])  
#In above comand is also UPDATE part, which takes some time. We could skip it if from time to time we would update the info in the whole db, for example when the gene name changes. 




DROP TABLE IF EXISTS classifications;
CREATE TABLE classifications (
	classification_id bigserial PRIMARY KEY,
	variant_id bigint NOT NULL,
	variant_class VARCHAR (255),
	class_datetime TIMESTAMP DEFAULT NOW(),
	comment VARCHAR (255),
	person VARCHAR (5),
	FOREIGN KEY (variant_id)
		REFERENCES variants (variant_id) ON UPDATE CASCADE ON DELETE CASCADE

);


CREATE TABLE populations (
	variant_id bigint PRIMARY KEY,
	metasvm_pred VARCHAR (255),
	metalr_pred VARCHAR (255),
	m_cap_pred VARCHAR (255),
	fathmm_pred VARCHAR (255),
	vest4_score VARCHAR (255),
	aloft_pred VARCHAR (255),
	rs_dbsnp151 VARCHAR (255),
	sift_pred VARCHAR (255),
	mutationtaster_pred VARCHAR (255),
	polyphen2_hvar_pred VARCHAR (255),
	Polyphen2_hdiv_pred VARCHAR (255),
	cadd_phred VARCHAR (255),
	revel_score VARCHAR (255),
	revel_rankscore VARCHAR (255),
	g1000gp3_ac VARCHAR (255),
	g1000gp3_af VARCHAR (255),
	g1000gp3_afr_ac VARCHAR (255),
	g1000gp3_afr_af VARCHAR (255),
	g1000gp3_eur_ac VARCHAR (255),
	g1000gp3_eur_af VARCHAR (255),
	g1000gp3_amr_ac VARCHAR (255),
	g1000gp3_amr_af VARCHAR (255),
	g1000gp3_eas_ac VARCHAR (255),
	g1000gp3_eas_af VARCHAR (255),
	g1000gp3_sas_ac VARCHAR (255),
	g1000gp3_sas_af VARCHAR (255),
	twinsuk_ac VARCHAR (255),
	twinsuk_af VARCHAR (255),
	alspac_ac VARCHAR (255),
	alspac_af VARCHAR (255),
	esp6500_aa_ac VARCHAR (255),
	esp6500_aa_af VARCHAR (255),
	esp6500_ea_ac VARCHAR (255),
	esp6500_ea_af VARCHAR (255),
	exac_adj_ac VARCHAR (255),
	exac_adj_af VARCHAR (255),
	exac_afr_ac VARCHAR (255),
	exac_afr_af VARCHAR (255),
	exac_amr_ac VARCHAR (255),
	exac_amr_af VARCHAR (255),
	exac_eas_ac VARCHAR (255),
	exac_eas_af VARCHAR (255),
	exac_fin_ac VARCHAR (255),
	exac_fin_af VARCHAR (255),
	exac_nfe_ac VARCHAR (255),
	exac_nfe_af VARCHAR (255),
	exac_sas_ac VARCHAR (255),
	exac_sas_af VARCHAR (255),
	gnomad_exomes_ac VARCHAR (255),
	gnomad_exomes_af VARCHAR (255),
	gnomad_exomes_afr_ac VARCHAR (255),
	gnomad_exomes_afr_af VARCHAR (255),
	gnomad_exomes_amr_ac VARCHAR (255),
	gnomad_exomes_amr_af VARCHAR (255),
	gnomad_exomes_asj_ac VARCHAR (255),
	gnomad_exomes_asj_af VARCHAR (255),
	gnomad_exomes_eas_ac VARCHAR (255),
	gnomad_exomes_eas_af VARCHAR (255),
	gnomad_exomes_fin_ac VARCHAR (255),
	gnomad_exomes_fin_af VARCHAR (255),
	gnomad_exomes_nfe_ac VARCHAR (255),
	gnomad_exomes_nfe_af VARCHAR (255),
	gnomad_exomes_sas_ac VARCHAR (255),
	gnomad_exomes_sas_af VARCHAR (255),
	clinvar_clnsig VARCHAR (255),
	clinvar_var_source VARCHAR (255),
	clinvar_review VARCHAR (255),
	ensembl_geneid VARCHAR (255),
	ensembl_transcriptid VARCHAR (255),
	ensembl_proteinid VARCHAR (255),
	hg19_pos_1_based VARCHAR (255),
	pop_update date,
	UNIQUE (variant_id, pop_update),
	FOREIGN KEY (variant_id)
		REFERENCES variants (variant_id) ON UPDATE CASCADE ON DELETE CASCADE

);

#CONSTRAINT fk_variant_id FOREIGN KEY (variant_id)   #### decided with Aaron and Martin to add pop_update_date column, so now it has to be one to many connection with variants table
#	REFERENCES variants (variant_id) ON UPDATE CASCADE ON DELETE CASCADE

CREATE TABLE variants_counts (
	variant_id bigint PRIMARY KEY,
	hom_inhouse int,
	het_inhouse int,
	CONSTRAINT fk_variant_id FOREIGN KEY (variant_id)
		REFERENCES variants (variant_id)
);

CREATE TABLE families (
	family_id int PRIMARY KEY,
	consanguineous_stat VARCHAR (20),
	region VARCHAR (255)
);

CREATE TABLE hpos (
	hpo_id serial PRIMARY KEY,
	hpo_number VARCHAR (20),
	hpo_label VARCHAR (255)
);

CREATE TABLE patients (
	patient_id int PRIMARY KEY,
	family_id int NOT NULL,
	aff_status VARCHAR (15) NOT NULL, 
	dob date NOT NULL,
	gender VARCHAR (25) NOT NULL,
	idx_status VARCHAR (25) NOT NULL,
	FOREIGN KEY (family_id)
		REFERENCES families (family_id) ON UPDATE CASCADE ON DELETE CASCADE
);

CREATE TABLE phenotypes (
	patient_id int,
	hpo_id int,
	PRIMARY KEY (patient_id, hpo_id),
	FOREIGN KEY (patient_id)
		REFERENCES patients (patient_id) ON UPDATE CASCADE ON DELETE CASCADE,
	FOREIGN KEY (hpo_id)
		REFERENCES hpos (hpo_id) ON UPDATE CASCADE ON DELETE CASCADE
);

ALTER TABLE patients ADD CONSTRAINT patients_family_id_fkey FOREIGN KEY (family_id) REFERENCES families (family_id) ON UPDATE CASCADE ON DELETE CASCADE;

#INSERT INTO patients (patient_id,family_id,aff_status,dob,gender,idx_status,region) VALUES (12345,23456,'affected','2011-03-23','male','index','Germany');

#INSERT INTO patients (patient_id,family_id,aff_status,dob,gender,idx_status) VALUES (13332,21334,'unaffected','1995-03-23','female','mother');

CREATE TABLE orders (
	order_id VARCHAR (20) PRIMARY KEY, 
	patient_id int NOT NULL,
	analysis_type VARCHAR (25) NOT NULL,	
	FOREIGN KEY (patient_id)
		REFERENCES patients (patient_id) ON UPDATE CASCADE ON DELETE CASCADE
);
###############now is int as order_id here, but we should check with projodis if it can stay like this


#INSERT INTO orders (order_id,patient_id,analysis_type) VALUES ('S30633A',12345,'WES_agi');
#orders cannot be integer because in our order numbers are letters



INSERT INTO pipeline_runs (sample_run_id,pipeline_exec_date,pipeline_version,status) VALUES (5,DEFAULT,'v0,1',2);
INSERT INTO pipeline_runs (seq_run_id, pipeline_exec_date) VALUES ((SELECT seq_run_ID FROM seq_runs WHERE order_id = 'S30633A' AND run_id = '190919_ST-K00352_0311_AHCKCCBBXY'),DEFAULT);
UPDATE pipeline_runs SET pipeline_exec_date = (SELECT pipeline_exec_date FROM pipeline_runs WHERE pipeline_run_id = 1) WHERE sample_run_id = 4;



CREATE TABLE seq_runs (
	seq_run_id serial PRIMARY KEY,
	run_id VARCHAR (255) NOT NULL,
	seq_machine VARCHAR (255),
	pe_sr VARCHAR (10),
	cycles_n VARCHAR (10),
	seq_qual_q30 VARCHAR (20),
	lib_prep VARCHAR (255),
	seq_facility VARCHAR (255),	
	submission VARCHAR (10),
	UNIQUE(run_id)

);

CREATE TABLE content_seq_runs (
	sample_run_id serial PRIMARY KEY,
	seq_run_id int,
	order_id VARCHAR (20),	
	read_count VARCHAR (20),
	lane_ratio VARCHAR (10),	
	FOREIGN KEY (order_id)
		REFERENCES orders (order_id) ON UPDATE CASCADE ON DELETE CASCADE,
	FOREIGN KEY (seq_run_id)
		REFERENCES seq_runs (seq_run_id) ON UPDATE CASCADE ON DELETE CASCADE
);

CREATE TABLE pipeline_runs (
	pipeline_run_id serial PRIMARY KEY,
	sample_run_id int NOT NULL,
	pipeline_exec_date TIMESTAMP DEFAULT NOW(),
	pipeline_version VARCHAR (25),
	status VARCHAR (10),
	coverage_mean VARCHAR (10),
	coverage_1x VARCHAR (10),
	coverage_5x VARCHAR (10),
	coverage_10x VARCHAR (10),
	coverage_20x VARCHAR (10),
	coverage_50x VARCHAR (10),
	coverage_70x VARCHAR (10),
	coverage_100x VARCHAR (10),
	FOREIGN KEY (sample_run_id)
		REFERENCES content_seq_runs (sample_run_id) ON UPDATE CASCADE ON DELETE CASCADE
);


SELECT patients.patient_id, families.family_id, orders.order_id, orders.analysis_type, patients.aff_status, patients.dob, patients.gender, patients.idx_status,families.consanguineous_stat, families.region,patients.dob FROM families INNER JOIN patients ON families.family_id=patients.family_id INNER JOIN orders ON orders.patient_id=patients.patient_id;
SELECT patients.patient_id, families.family_id, orders.order_id, orders.analysis_type, patients.aff_status, patients.dob, patients.gender, patients.idx_status,families.consanguineous_stat, families.region,patients.dob, hpos.hpo_number FROM families INNER JOIN patients ON families.family_id=patients.family_id INNER JOIN orders ON orders.patient_id=patients.patient_id LEFT JOIN phenotypes ON patients.patient_id=phenotypes.patient_id LEFT JOIN hpos ON phenotypes.hpo_id=hpos.hpo_id;

CREATE TABLE detections (
	variant_id bigint,
	pipeline_run_id int,
	qual VARCHAR (25),
	filter VARCHAR (255),
	ac VARCHAR (25),  
	af VARCHAR (25), 
	an VARCHAR (25), 
	baseqranksum VARCHAR (25), 
	dp VARCHAR (25), 
	ds VARCHAR (25),
	pos_end VARCHAR (255),
	excessHet VARCHAR (25), 
	fs VARCHAR(25), 
	inbreedingcoeff VARCHAR (255),
	mleac VARCHAR (25),
	mleaf VARCHAR (25), 
	mq VARCHAR (25), 
	mqranksum VARCHAR (25),
	n_t_s VARCHAR (25),
	p_t_s VARCHAR (25),
	qd VARCHAR (25), 
	readposranksum VARCHAR (25),
	vqslod VARCHAR (25),
	culprint VARCHAR (25), 
	sor VARCHAR (25), 
	gt_g VARCHAR (10),
	ad_g0 VARCHAR (10),
	ad_g1 VARCHAR (10),
	dp_g VARCHAR (10), 
	gq_g VARCHAR (10), 
	pl_g_00 VARCHAR (10), 
	pl_g_01 VARCHAR (10),
	pl_g_11 VARCHAR (10),
	PRIMARY KEY (variant_id, pipeline_run_id),
	FOREIGN KEY (variant_id)
		REFERENCES variants (variant_id) ON UPDATE CASCADE ON DELETE CASCADE,
	FOREIGN KEY (pipeline_run_id)
		REFERENCES pipeline_runs (pipeline_run_id) ON UPDATE CASCADE ON DELETE CASCADE
);

CREATE TABLE pheno_genes ( 
pheno_gene_id serial PRIMARY KEY, 
gene_entrez_id VARCHAR (20), 
gene_symbol VARCHAR (20) 
); 

  

CREATE TABLE gene_hpo ( 
pheno_gene_id int,  
hpo_id int, 
PRIMARY KEY (pheno_gene_id, hpo_id),
FOREIGN KEY (pheno_gene_id) 
REFERENCES pheno_genes (pheno_gene_id) ON UPDATE CASCADE ON DELETE CASCADE, 
FOREIGN KEY (hpo_id) 
REFERENCES hpos (hpo_id) ON UPDATE CASCADE ON DELETE CASCADE 
); 


GRANT CONNECT ON DATABASE test_variants_db TO geneticsuser; ###################change name of db here from test to real
GRANT USAGE ON SCHEMA public TO geneticsuser;
GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA public TO geneticsuser;
GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA public TO geneticsuser;
#shall we include original samples names in sample table?



Table variants_counts {
  variant_id bigint
  hom_inhouse int
  het_inhouse int
}
//we should discuss if we need this table
Ref: variants_counts.variant_id - variants.variant_id

 


CREATE TABLE [IF NOT EXISTS] table_name (
   column1 datatype(length) column_contraint,
   column2 datatype(length) column_contraint,
   column3 datatype(length) column_contraint,
   table_constraints
);

CREATE TABLE accounts (
	user_id serial PRIMARY KEY,
	username VARCHAR ( 50 ) UNIQUE NOT NULL,
	password VARCHAR ( 50 ) NOT NULL,
	email VARCHAR ( 255 ) UNIQUE NOT NULL,
	created_on TIMESTAMP NOT NULL,
        last_login TIMESTAMP 
);