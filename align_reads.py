import os, glob
from pathlib import Path

raw_data_folder = '00_Raw'
qc_folder = '01_QC'
aligned_folder = '02_Aligned'
genome_file = 'bacteroides_thetaiotamicron_VPI5482_ref_genome.fa'
gtf_file = 'bacteroides_thetaiotamicron_VPI5482.gtf'

def quality_filtering():
	Path(qc_folder).mkdir(parents=True, exist_ok=True)
	files = glob.glob('%s/*.fastq.gz' % (raw_data_folder))
	for f in files:
		os.system('trim_galore --gzip -o %s -j 4 %s' % (qc_folder, f))

def align_reads_bwa():
	Path(aligned_folder).mkdir(parents=True, exist_ok=True)
	files = glob.glob('%s/*.fq.gz' % (qc_folder))
	os.system('bwa index %s' % (genome_file))
	for f in files:
		filename = f.split('/')[-1].split('.')[0]
		sam_file = '%s/%s.sam' % (aligned_folder, filename)
		print('Aligning %s using bwa' % (filename))
		os.system('bwa mem -v 1 %s %s > %s' % (genome_file, f, sam_file))
		print('Converting sam to bam')
		bam_file = '%s/%s.bam' % (aligned_folder, filename)
		os.system('samtools view --threads 4 -b -o %s %s' % (bam_file, sam_file))
		os.system('rm %s' % (sam_file))
		print('Sorting bam file')
		sorted_bam = '%s/%s-sorted.bam' % (aligned_folder, filename) 
		os.system('samtools sort -m 7G --threads 4 -o %s %s' % (bam_file, sorted_bam))
		os.system('rm %s' % (bam_file))
		print('')

if __name__ == '__main__':
	#quality_filtering()
	align_reads_bwa()