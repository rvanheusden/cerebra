"""
This program takes in two sets of vcfs, single-cell and bulk (peripheral
blood, ie. germline) and filters out the common variants found in both sc 
and bulk. Creates a new directory, filteredOut, that contans the filtered vcfs. 
"""

import os
import sys
import shutil
import warnings
import pandas as pd
import click
import VCF
warnings.simplefilter(action='ignore', category=FutureWarning)


def create_final_outdir():
	""" creates a finalized out dir that has all of the filtered vcfs as
		well as the ones we didnt have germline for """

	filterDir = cwd + 'filteredOut/'
	filterDir_list = os.listdir(filterDir)

	filteredCells = []
	for f in filterDir_list:
		cell = f.strip('_unique.vcf')
		filteredCells.append(cell)

	epiDir = cwd + 'scVCF/'
	epiDir_list = os.listdir(epiDir)

	epiCells = []
	for f in epiDir_list:
		cell = f.strip('.vcf')
		epiCells.append(cell)

	# get cells in epiCells but NOT filteredCells
	nonFilteredCells = set(epiCells) - set(filteredCells)

	nonFilteredCells_list = []
	for cell in nonFilteredCells:
		f = cell + '.vcf' 
		nonFilteredCells_list.append(f)

	cmd1 = 'sudo mkdir -p ' + cwd + 'scVCF_filtered_all/'
	cmd2 = 'sudo chmod -R 777 ' + cwd + 'scVCF_filtered_all/'
	os.system(cmd1)
	os.system(cmd2)

	# copy over the non-filtered cells
	outPATH = cwd + 'scVCF_filtered_all/'
	for file in nonFilteredCells_list:
		src = epiDir + file
		dst = outPATH + file
		shutil.copyfile(src, dst)

	# copy over all the filtered cells
	for file in filterDir_list:
		f = file.strip('_unique.vcf')
		f = f + '.vcf'
		src = filterDir + file
		dst = outPATH + f
		shutil.copyfile(src, dst)



def write_vcf(df, outStr_):
	""" routine for writing VCF files, from an existing dataframe. 
		essentially just adding in this horrible vcf header. """
	with open(cwd + 'vcfheader.txt', 'r') as f:
		header = f.read()
	
		df['QUAL'] = 60.28
		df['FILTER'] = '.'
		df['INFO'] = 'AC=2;AF=1.00;AN=2;DP=7;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=3.00;QD=30.14;SOR=2.303'

		output_VCF = outStr_
		with open(output_VCF, 'w') as vcf:
			vcf.write(header)

		df.to_csv(output_VCF, sep="\t", mode='a', index=False)

	# replace this damn CHROM field
	with open(output_VCF) as f:
		newText = f.read().replace('CHROM', '#CHROM')

	with open(output_VCF, "w") as f:
		f.write(newText)



def get_unique_vcf_entries(patient):
	""" do the germline filter, and return a dataframe with only the 
		UNIQUE entries for a given patient """
	germlinePATH = cwd + 'bulkVCF/' + patient + 
	tumorExomePATH = cwd + 'tumorExome/' + patient + '.vcf'
	
	try:
		germline_df = VCF.dataframe(germlinePATH)
		tumorExome_df = VCF.dataframe(tumorExomePATH)
	except FileNotFoundError:
		print('FILE NOT FOUND: %s' % germlinePATH)
		return
    
	germline_df_trimmed = germline_df[['CHROM', 'POS', 'ID', 'REF', 'ALT']]
	tumorExome_df_trimmed = tumorExome_df[['CHROM', 'POS', 'ID', 'REF', 'ALT']]
    
	# get whats SHARED between patient and cell 
	germline_tumorExome_concat = pd.concat([germline_df_trimmed, tumorExome_df_trimmed])
	rowsToKeep = germline_tumorExome_concat.duplicated()
	germline_tumorExome_shared = germline_tumorExome_concat[rowsToKeep]
	germline_tumorExome_shared = germline_tumorExome_shared.reset_index(drop=True)

	# now go back to the original cell df, pull out whats UNIQUE 
	tumorExome_tumorExome_concat = pd.concat([tumorExome_df_trimmed, germline_tumorExome_shared])
	tumorExome_tumorExome_concat_noDups = tumorExome_tumorExome_concat.drop_duplicates(keep=False)
	tumorExome_tumorExome_concat_noDups = tumorExome_tumorExome_concat_noDups.reset_index(drop=True)
    
	return(tumorExome_tumorExome_concat_noDups)



""" launch """
@click.command()
@click.option('--test', default = False)
@click.option('--wrkdir', default = '/Users/lincoln.harris/code/cerebra/cerebra/wrkdir/', prompt='s3 import directory', required=True)



def germline_filter_tumor_exome(test, wrkdir):
	""" given a set of tumor exome (bulk) vcfs and bulk peripheral blood vcfs, this
		program subtracts out the mutations common to tumor exome and bulkVCF. """
	global patientMetadata
	global cwd

	cwd = wrkdir

	# read in patient metadata
	patientMetadata = pd.read_csv(cwd + 'metadata_all_cells_4.10.19.csv')

	# get a list of all the single-cell VCF files
	vcfDir = cwd + 'tumorExome/'
	tumorExome_list = os.listdir(vcfDir)

	# get list of bulk VCF files
	bulkVCF_dir = cwd + 'bulkVCF/'
	bulkVCF_list = os.listdir(bulkVCF_dir)

	patientsRun = []

	cmd1 = 'sudo mkdir -p ' + cwd + 'filteredOut'
	cmd2 = 'sudo chmod -R 777 ' + cwd + 'filteredOut/' 
	os.system(cmd1)
	os.system(cmd2)

	# outer loop -- by PATIENT
	for item in bulkVCF_list:
		currSample = item.strip('.vcf')
		currPatient = currSample.split('_')[0]
		suffix1 = currSample.split('_')[1]
		try:
			suffix2 = currSample.split('_')[2]
		except IndexError:
			suffix2 = ''
	
		if suffix2 != '' and currPatient not in patientsRun:
			print('TUMOR EXOME FOUND, for %s' % currPatient)

			tumorExome_unique = get_unique_vcf_entries(currPatient)
			outStr = cwd + 'filteredOut/' + currPatient + '_unique.vcf'
			write_vcf(tumorExome_unique, outStr)
			
			patientsRun.append(currPatient)

	if test:
		cmd = 'cp -r ' + cwd + 'filteredOut/ ' + cwd + 'test/germline_filter/' 
		os.system(cmd)
	else:
		create_final_outdir()
