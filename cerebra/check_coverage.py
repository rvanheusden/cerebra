""" this tool compares VCF records to gVCF for a given cell, to determine coverage
	to a given loci. keep in mind the loci should be SMALL, ie. individual SNPs / 
	small indels. it is not intended for whole exon or whole transcript queries """

import pandas as pd
import numpy as np
from . import VCF
import sys
import multiprocessing as mp
import os
import click
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning) # fuck this message


def get_filenames():
	""" get file names given path """
	files = []
	for file in os.listdir(cwd + "scVCF_filtered_all/"):
		if file.endswith(".vcf"):
			fullPath = cwd + 'scVCF_filtered_all/' + file 
			files.append(fullPath)
    
	return files



def build_outfile_line(outCode, depth, cell_name):
	""" generic func for creating a single line pd df with cellName, 
		coverage (bool) and depth """
	colNames = ['cellName', 'coverage_bool', 'depth']

	if outCode == 1: # no records found
		toAddRow = pd.DataFrame([[cell_name, 0, 0]], columns=colNames)
	elif outCode == 2: # single record found
		toAddRow = pd.DataFrame([[cell_name, 1, depth]], columns=colNames)
	else: # multiple records found
		toAddRow = pd.DataFrame([[cell_name, 1, depth]], columns=colNames)

	return(toAddRow)



def get_depth(df, cellName_):
	""" given a dataframe containing multiple records, reports 
		depth for every record within that df """
	outCode_ = 0 # were gonna send this to buildOutputFile()

	if len(df.index) == 0: 		# no records found
		outCode_ = 1
		toAddRow_ = build_outfile_line(outCode_, 0, cellName_)
	elif len(df.index) == 1: 	# single record found
		outCode_ = 2
		infoStr = df['INFO']
		infoStr = str(infoStr)
		DP = infoStr.split('DP')[1].split(';')[0].strip('=')
		toAddRow_ = build_outfile_line(outCode_, DP, cellName_)
	else:						 # multiple records found
		outCode_ = 3						
		infoDF = df['INFO']
		DP_vec = []

		for i in range(0, len(infoDF.index)-1):
			line = infoDF.iloc[i]
			line = str(line)
			try:
				DP = line.split('DP')[1].split(';')[0].strip('=')
				DP_vec.append(DP)
			except IndexError:
				continue

		toAddRow_ = build_outfile_line(outCode_, DP_vec, cellName_)

	return(toAddRow_)



def get_GOI_record(record, *args):
	""" defines a list of records corresponding to the GOI """
	chrom = 'chr' + str(args[0])
	start = int(args[1])
	end = int(args[2])

	if record['CHROM'] == chrom:
		if end >= record['POS'] >= start:
			return 1
		else:
			return 0
	else:
		return 0



def run_batch(file):
	""" implements BATCH MODE. for every cell, call subroutines to search
		for ROI, get depth, and output """
	try:
		cellName = file.strip(cwd + 'scVCF_filtered_all/')
		cellName = cellName.strip('.vcf')
		
		vcf_path = file
		gvcf_path = cwd + 'gVCF/' + cellName
		gvcf_path_strip = gvcf_path.rstrip() + '.g.vcf'

		vcf = VCF.dataframe(vcf_path)
		gvcf = VCF.dataframe(gvcf_path_strip)

		# get a list of the records we actually care about
		toKeepList_v = vcf.apply(get_GOI_record, axis=1, args=(chrom_, start_ ,end_))
		toKeepList_g = gvcf.apply(get_GOI_record, axis=1, args=(chrom_, start_, end_))

		# subset by relevant records
		vcf_GOI = vcf[np.array(toKeepList_v, dtype=bool)]
		gvcf_GOI = gvcf[np.array(toKeepList_g, dtype=bool)]

		# get depth of coverage, for relevant records
		outputRow_v = get_depth(vcf_GOI, cellName)
		outputRow_g = get_depth(gvcf_GOI, cellName)

		# make the combined row, with both vcf and gvcf fields filled in
		outputRow_comb = pd.DataFrame(columns=colNames) # colNames is a global
		outputRow_comb['cellName'] = outputRow_v['cellName']
		outputRow_comb['coverage_bool_vcf'] = outputRow_v['coverage_bool']
		outputRow_comb['depth_vcf'] = outputRow_v['depth']
		outputRow_comb['coverage_bool_gvcf'] = outputRow_g['coverage_bool']
		outputRow_comb['depth_gvcf'] = outputRow_g['depth']
	
	except:
		outputRow_comb = pd.DataFrame(columns=colNames) # just an empty row
		# fill in this row with something 
	return(outputRow_comb)



def init_pool(filenames, outfile_name):
	""" function to set up the thread pool """

	# init outFile
	outputDF_init = pd.DataFrame(columns=colNames) 	
			
	print('creating pool')

	p = mp.Pool(processes=nThreads)

	print('running...')
	outputRows = p.map(run_batch, filenames)
	p.close()
	p.join()
	print('done!')
	print(' ')

	# join all of the rows into single df
	cmd = 'sudo mkdir -p ' + cwd + 'coverage/'
	cmd1 = 'sudo chmod -R 777 ' + cwd + 'coverage/'
	os.system(cmd)
	os.system(cmd1)

	outputDF = outputDF_init.append(outputRows)
	outputDF.to_csv(cwd + 'coverage/' + outfile_name, index=False)



""" get cmdline input """
@click.command()
@click.option('--chrom', default = 7, prompt='chromosome', required=True, type=str)
@click.option('--start_pos', default = 55191820, prompt='start position', required=True, type=str)
@click.option('--end_pos', default = 55191822, prompt='end position', required=True, type=str)
@click.option('--nthreads', default = 4, prompt='number of threads', required=True, type=int)
@click.option('--wrkdir', default = '/Users/lincoln.harris/code/cerebra/cerebra/wrkdir/', 
	prompt='s3 import directory', required=True)
@click.option('--batch_mode', default = 0, prompt='batch mode', required=True, type=int)
@click.option('--outfile', default = 'egfr_L858R_coverageByCell.csv', 
	prompt='name of output file (leave blank for batch mode)', required=True)
@click.option('--test', default = False)



def check_coverage(chrom, start_pos, end_pos, nthreads, wrkdir, batch_mode, outfile, test):
	""" check coverage to a given ROI """
	global cellName
	global vcf_s3_path
	global gvcf_s3_path
	global chrom_
	global start_
	global end_
	global colNames
	global cwd
	global nThreads

	cwd = wrkdir
	nThreads = nthreads
	
	print(' ')
	print('this tool should be used for loci specific coverage queries.')
	print('it is NOT intended for calculating coverage at the exon/transcript level.')

	fNames = get_filenames()
	colNames = ['cellName', 'coverage_bool_vcf', 'depth_vcf', 'coverage_bool_gvcf', 'depth_gvcf']

	if batch_mode:
		cov_df = pd.read_csv(cwd + '../coverageBatch.csv')

		for i in range(0, len(cov_df.index)):
			currRow = cov_df.iloc[i]
			chrom_ = currRow['chrom']
			start_ = currRow['start_pos']
			end_ = currRow['end_pos']
			outfile_ = currRow['outfile']

			init_pool(fNames, outfile_)
	else:
		chrom_ = chrom
		start_ = start_pos
		end_ = end_pos
		outfile_ = outfile

		init_pool(fNames, outfile_)
