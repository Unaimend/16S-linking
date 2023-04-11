#All reads
read_folder="/work_beegfs/sukem127/clean/REAL_16s/reads"
sample_list="/work_beegfs/sukem127/clean/REAL_16s/reads/samples"
forward_reads_name="2M_finished_forward_reads.fastq.gz"
reverse_reads_name="2M_finished_reverse_reads.fastq.gz"
pooled_forward_reads=read_folder+"/"+forward_reads_name
pooled_reverse_reads=read_folder+"/"+reverse_reads_name
bin_folder="/work_beegfs/sukem127/clean/REAL_16s/fastaBins"
out_folder="/work_beegfs/sukem127/clean/16ss"
assembly_threads=16

rule	assemble16s:
	input:
		pooled_forward_reads,
		pooled_reverse_reads
	output:
		"phyloFlash/LIB.all.final.phyloFlash.dbhits.fa"
			
	conda:
		"16s_recovery"
	shell:
		"""
			phyloFlash.pl -lib LIB -almosteverything -read1 {pooled_forward_reads} -read2 {pooled_reverse_reads} -CPUs={assembly_threads}
			# this mkdir should only be executed if phyFlash returned succesfully. How do we do this?
			# mkdir phyloFlash
			mv LIB.* phyloFlash/
		"""

rule	create_16s_index:
	input:
		"phyloFlash/LIB.all.final.phyloFlash.dbhits.fa",
	output:
		"db_16s/16s_db.ndb"
	conda:
		"16s_recovery"

	shell:
		"""
			makeblastdb -dbtype nucl -in phyloFlash/LIB.all.final.phyloFlash.dbhits.fa -out 16s_db 
			mv 16s_db* db_16s
		"""

	
rule	create_bin_yaml:
	input:
		bin_folder
	output:
		"bin_list.yaml"
	log:
		out = "logs/bin_list_creation.log"
	script:
		"scripts/python/convert_folder_to_bin_yaml.py"


rule	create_bin_txt:
	input:
		bin_folder
	output:
		"BinNames.txt"
	shell:	
		"""
		find {bin_folder}/*.fa -printf "%f\n"  > BinNames.txt
		"""	

		

def read_bins():
	import yaml

	with open("bin_list.yaml", "r") as f:
		dict = yaml.safe_load(f)
		samples = dict["samples"]
		return list(samples.values())


#use checkpoints for paralel.. https://stackoverflow.com/questions/57925104/snakemake-running-single-jobs-in-parallel-from-all-files-in-folder
#TALK WITH LENA AGAIN: PROPER RENAMIN_G
#TODO LATER IN THE DOC. WE DO CREATE CONDUCTED BINS I>E> RENAME THEM IN THE BELOW FASHION AND MERGE THEM. SO MAYBE THIS IS UNNCECESSARY
rule rename_read_ids:
	input:
		#Maybe use glob wildcards https://carpentries-incubator.github.io/snakemake-novice-bioinformatics/05-expansion/index.html
		expand("{folder}/{t}", t = read_bins(), folder=bin_folder)
	output:
		directory("renamed_bins"),
		expand("renamed_bins/{t}", t = read_bins(), folder=bin_folder)
	conda:
		"test.yaml"
	script:
		"scripts/python/rename_fasta_header.py"


rule pool_bins:
	input:
		expand("{out_folder}/renamed_bins/{t}", t = read_bins(), out_folder=out_folder)
	output:
		"{out_folder}/renamed_bins/pooled.fa"
	shell:
		"""	
		cat {input} > "{out_folder}/renamed_bins/pooled.fa"
		"""


rule all_16s_to_all_bins: 
	input:
		out_folder + "/renamed_bins/pooled.fa"
	output:
		"blast_16s_all_bins.outfmt6"
	conda:
		"test.yaml"
	shell:
		"""
		blastn \
		  -query {input} \
		  -db db_16s/16s_db \
		  -max_target_seqs 50 \
		  -outfmt 6 \
		  -evalue 0.000001 \
		  -num_threads 20 \
		  -out  blast_16s_all_bins.outfmt6
		"""			
#Filter Blast results (min ident 97% and min coverage 100bp)
rule filter_blast_and_format:
	input: 
		"blast_16s_all_bins.outfmt6"
	output:
		"16snames.txt"
	shell:
		"""	
		cat {input} | \
		awk '$3 > 97' | awk '$4 > 100' | cut -f1,2,12 > filtered-BlastOut.tab && \
		grep ">"  phyloFlash/LIB.all.final.phyloFlash.dbhits.fa | cut -d ' ' -f 1 | tr -d '>' > 16snames.txt

		"""

rule create_blast_matrix:
	input:
		out_folder + "/BinNames.txt",
		out_folder + "/filtered-BlastOut.tab",
		out_folder + "/16snames.txt"
	output:
		"matrix-Blast.tab"
	shell: 
		"""
		# has to be chmodded chmod +x
		scripts/bash/reformat-blast-results.sh {input[1]} {input[2]} {input[0]} && rm temp-pasting*
		"""

# Build index over the bins
rule bbsplit_bin_index:
	input:
		bin_folder
	output:
		directory(out_folder + "/bin_index/ref/index/")
	conda:
		"bbmap.yaml"
	shell:
		"""
		bbsplit.sh ref={input} path={out_folder}/bin_index
		"""


def read_bins2():
	with open(str(sample_list), "r") as f:
		return f.read().splitlines()


# Maps reads from all the samples to all bins
rule map_reads_to_bins:
	input:
		samples={sample_list},
		read_path={read_folder}
	output: 
		#EXPAND OUTPUT AND SEE 
	conda:
		"bbmap.yaml"
	resources:
		slurm_extra="--mail-type=ALL --mail-user=stu203329@mail.uni-kiel.de --array 0-51"
	shell:
		"""
		#TODO MANUALLY EDIT THE ARRAY SIZE PLEASE
		#TODO IS THIS A BAD IDEA?  DELETING LOGS WITHOUT USER PERM
		#ES GIBT --array wir koennte da die lines count passen
		#TODO ADD pATH FOR SLURM ERRORS AS pARAMETER
		sbatch "scripts/bash/slurm/run_bbsplit.sh" {sample_list} {read_folder} {out_folder}
		#TODO WHY DOES THIS NOT WORK
		mv {out_folder}/splited-* {out_folder}/slurm_out/reads_to_bins/
		"""

rule map_mapped_bin_reads_to_16s:
	input:
		expand("{folder}/{t}", t = read_bins(), folder=bin_folder)
	output: 
		#EXPAND OUTPUT AND SEE 
	resources:
		slurm_extra="--mail-type=ALL --mail-user=stu203329@mail.uni-kiel.de --array 0-180"
	shell:
		"""
		#TODO MANUALLY EDIT THE ARRAY SIZE PLEASE
		#TODO IS THIS A BAD IDEA?  DELETING LOGS WITHOUT USER PERM
		#ES GIBT --array wir koennte da die lines count passen
		#TODO ADD pATH FOR SLURM ERRORS AS pARAMETER
		sbatch "scripts/bash/slurm/map_mapped_bin_reads_to_16s.sh" \
			{input}
		#TODO WHY DOES THIS NOT WORK
		#mv {out_folder}/splited-* {out_folder}/slurm_out/reads_to_bins/
		"""
#TODO WILL NET 
rule summarize_bin_16s_align_statistics:
	#Not way to specifiy a dir as input?
	input:
	
	params:
		path =  "scafstats-statsfiles-unpaired-all/"
	conda:
		"test.yaml"
	shell:
		"python3 scripts/python/make-splitmapping-stats.py"



rule create_16s_abund:
	output: 
		directory("16SrRNA-cov-statsfiles"),
		directory("16SrRNA-rpkm-statsfiles"),
		directory("16SrRNA-scafstats-statsfiles"),
		directory("16SrRNA-statsfiles")
	params:
		samples=read_bins2(),
		read_path={read_folder}
	conda:
		"bbmap.yaml"
	resources:
		slurm_extra="--mail-type=ALL --mail-user=stu203329@mail.uni-kiel.de --array 1-52" #52
	shell:
		""" #TODO THIS IS WRONG DO MNOT USE sbatch THIS CIRCUMVENTS SNAKEMAKES MECHANISM
		sbatch "scripts/bash/slurm/create_16s_abund_over_sampls.sh" {params.read_path} {params.samples} 
		"""

rule summarize_16S_abundance_statistics:
	input:
		expand("{folder}/{t}.scafstats", t = read_bins2(), folder="16SrRNA-scafstats-statsfiles")
	output:
		"16S-abundances.tab"
	shell:
		"""
			"scripts/bash/make-16S-mapping-stats.sh"
		"""		
	


rule create_ref_over_pooled_bins:
	input:
		#### POOLED & RENAMED BINS BINS
		#TODO BBMAP DIES WHEN USING LEADING /
		"renamed_bins/pooled.fa"
	conda:
		"bbmap.yaml"
	shell: "bbmap.sh ref={input}"

	


rule create_bin_abund:
	#output: 
	#	directory("bin-statsfiles"),
	#	directory("bin-cov-statsfiles"),
	#	directory("bin-scafstats-statsfiles"),
	#	directory("bin-rpkm-statsfiles")
	params:
		samples=read_bins2(),
		read_path={read_folder}
	conda:
		"bbmap.yaml"
	resources:
		slurm_extra="--mail-type=ALL --mail-user=stu203329@mail.uni-kiel.de --array 1-52 -e  /zfshome/sukem127/e_%A_%a.log  -o  /zfshome/sukem127/o_%A_%a.log" #52
	shell:
		""" #TODO THIS IS WRONG DO MNOT USE sbatch THIS CIRCUMVENTS SNAKEMAKES MECHANISM
		sbatch "scripts/bash/slurm/create_bin_abund_over_sample.sh" {params.read_path} {params.samples} 
		"""



rule summarize_bin_abundance_statistics:
	input:
		expand("{folder}/{t}.scafstats", t = read_bins2(), folder="bin-scafstats-statsfiles")
	output:
		"bin-abundances.tab"
	shell:
		"""
			"scripts/bash/make-bin-mapping-stats.sh"
		"""		

rule create_bin_metadata:
	input: 
		expand("{folder}/{t}", t = read_bins(), folder=bin_folder)
	output: "BinMetaData.csv"
	shell:
		"""
			for item in {input}; do
				echo $item
				BinFileName=${{item##*/}}
				BinName=${{BinFileName%.*}}
				echo $BinName
				fgrep -v ">" ${{item}} | wc -m | sed -e "s/^/${{BinName}}\t/" >> BinMetaData.csv
			done
		"""		
