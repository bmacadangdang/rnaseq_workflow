import os, glob, csv
import pandas as pd
import numpy as np
from pathlib import Path

raw_data_folder = '00_Raw'
qc_folder = '01_QC'
aligned_folder = '02_Aligned'
count_folder = '03_Counts'
counts_filtered_folder = '04_Counts_Filtered'
normalized_folder = '05_Results_normalized'

genome_file = 'bacteroides_thetaiotamicron_VPI5482_ref_genome.fa'
gtf_file = 'bacteroides_thetaiotamicron_VPI5482.gtf'
gff_file = 'bacteroides_thetaiotamicron_VPI5482_gene_annotations.gff'

def quality_filtering():
	Path(qc_folder).mkdir(parents=True, exist_ok=True)
	files = glob.glob('%s/*.fastq.gz' % (raw_data_folder))
	for f in files:
		os.system('trim_galore --gzip -o %s -j 4 %s' % (qc_folder, f))

def align_reads_bwa():
	Path(aligned_folder).mkdir(parents=True, exist_ok=True)
	files = glob.glob('%s/*.fq.gz' % (qc_folder))
	os.system('bwa index %s' % (genome_file))

	already_aligned = glob.glob('%s/*.bam' % (aligned_folder))
	already_aligned = [f.split('/')[-1].split('-')[0] for f in already_aligned]

	for f in files:
		filename = f.split('/')[-1].split('.')[0]
		if filename not in already_aligned:
			sam_file = '%s/%s.sam' % (aligned_folder, filename)
			print('Aligning %s using bwa' % (filename))
			os.system('bwa mem -v 1 -t 4 %s %s > %s' % (genome_file, f, sam_file))
			print('Converting sam to bam')
			bam_file = '%s/%s.bam' % (aligned_folder, filename)
			os.system('samtools view --threads 4 -b -o %s %s' % (bam_file, sam_file))
			os.system('rm %s' % (sam_file))
			print('Sorting bam file')
			sorted_bam = '%s/%s-sorted.bam' % (aligned_folder, filename) 
			os.system('samtools sort -m 4G --threads 2 -o %s %s' % (sorted_bam, bam_file))
			os.system('rm %s' % (bam_file))
			print('')

def count_genes():
	Path(count_folder).mkdir(parents=True, exist_ok=True)
	files = glob.glob('%s/*-sorted.bam' % (aligned_folder))
	for f in files:
		filename = '_'.join(f.split('/')[-1].split('.')[0].split('_')[:-1])
		count_file = '%s/%s_counts.txt' % (count_folder, filename)
		os.system('featureCounts -T 4 -t gene -s 2 -a %s -o %s %s' % (gtf_file, count_file, f))

def remove_rrna():
	Path(counts_filtered_folder).mkdir(parents=True, exist_ok=True)
	rrna_genes = []
	with open(gff_file, 'r') as f:
		reader = csv.reader(f, delimiter='\t')
		for line in reader:
			if line[0][0] != '#':
				if line[2] == 'rRNA':
					parts = line[8].split(';')
					for p in parts:
						parts2 = p.split('=')
						if parts2[0] == 'ID':
							rrna_genes.append(parts2[1].split('-')[1])
	
	files = glob.glob('%s/*_counts.txt' % (count_folder))
	for f in files:
		filename = f.split('/')[-1].split('.')[0]
		counts_filtered_file = '%s/%s_filtered.txt' % (counts_filtered_folder, filename)
		with open(f, 'r') as counts_file, open(counts_filtered_file, 'w') as filtered_file:
			reader = csv.reader(counts_file, delimiter='\t')
			writer = csv.writer(filtered_file, delimiter='\t')
			next(reader)
			header = next(reader)
			writer.writerow(header)
			for line in reader:
				gene = line[0]
				if gene not in rrna_genes:
					writer.writerow(line)

def normalize_counts():
	Path(normalized_folder).mkdir(parents=True, exist_ok=True)
	first = True
	genes_lengths = {}
	with open(gff_file, 'r') as f:
		reader = csv.reader(f, delimiter='\t')
		for line in reader:
			if line[0][0] != '#':
				if line[1] == 'RefSeq':
					gene = line[8].split(';')[0][8:]
					gene_length = int(line[4]) - int(line[3]) + 1
					genes_lengths[gene] = gene_length

	files = glob.glob('%s/*_counts_filtered.txt' % (counts_filtered_folder))
	for f in files:
		filename = '_'.join(f.split('/')[-1].split('_')[:-2])
		data = pd.read_csv(f, sep='\t')
		gene_ids = data['Geneid']
		counts = np.array(data[data.columns[-1]])
		if first:
			all_counts = pd.DataFrame({filename: counts}, index=gene_ids)
			first = False
		else:
			gene_counts = pd.DataFrame({filename : counts}, index=gene_ids)
			all_counts = all_counts.merge(gene_counts, 'outer', left_index=True, right_index=True)

	all_counts.to_csv('%s/all_counts.txt' % (normalized_folder), sep='\t')

	rpkm = pd.read_csv('%s/all_counts.txt' % (normalized_folder), sep='\t', index_col='Geneid')
	for col in rpkm.columns:
		total_count = np.sum(rpkm[col])
		factor = total_count / 1000000.
		rpkm[col] = rpkm[col] / factor

	for gene_id in rpkm.index:
		rpkm.loc[gene_id] = rpkm.loc[gene_id] / genes_lengths[gene_id]

	rpkm.to_csv('%s/rpkm.txt' % (normalized_folder), sep='\t')

	tpm = pd.read_csv('%s/all_counts.txt' % (normalized_folder), sep='\t', index_col='Geneid')
	for gene_id in tpm.index:
		tpm.loc[gene_id] = tpm.loc[gene_id] / genes_lengths[gene_id]

	for col in tpm.columns:
		total_counts = np.sum(tpm[col])
		factor = total_counts / 1000000.
		tpm[col] = tpm[col] / factor

	tpm.to_csv('%s/tpm.txt' % (normalized_folder), sep='\t')




if __name__ == '__main__':
	#quality_filtering()
	#align_reads_bwa()
	#count_genes()
	#remove_rrna()
	normalize_counts()