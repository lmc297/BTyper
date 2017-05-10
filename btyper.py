#!/usr/bin/env python

# BTyper v1.0:  Bacillus spp. Typer
# October 10, 2015
# Created by Laura Carroll
# lmc297@cornell.edu
# Usage:  python btyper.py -i <path/to/input/bacillus.fasta> -o <path/to/desired/output_directory/> -t <type_of_input_data (sequence, PE reads, SE reads)>
# Optional arguments:  [-draft_genome]

# import required packages
import argparse, sys, os, re, glob, collections
from Bio.Blast.Applications import NcbiblastnCommandline, NcbitblastnCommandline
from Bio import SeqIO
from Bio.Blast import NCBIXML


# parse arguments
parser = argparse.ArgumentParser(usage='python btyper.py -i </path/to/input/file.extension> -o </path/to/desired/output_directory/> -t [input data format (seq, pe, se, sra, or sra-get)] [-other options]')
parser.add_argument('-i','--input', nargs='+', help='Enter the path to the Bacillus cereus group sequence data you would like to input, or enter an SRA accession number',required=True)
parser.add_argument('-o', '--output', nargs='+', help='Specify a path to your desired output directory',required=True)
parser.add_argument('-t', '--type', nargs='+', help='Specify your type of data: seq for genomes or contigs in fasta format, pe for paired-end Illumina reads, se for single-end Illumina reads, sra for a sra file, or sra-get with an SRA accession number',required=True)
parser.add_argument('--draft_genome', action='store_true', default=False, help='Optional argument for use with contigs in fasta format; concatenates draft genome contigs into pseudochromosome')
parser.add_argument('-v', '--virulence', nargs='?', default=True, help='Optional argument, True or False; perform virulence gene typing; default=True')
parser.add_argument('-m', '--mlst', nargs='?', default=True, help='Optional argument, True or False; perform MLST using Bacillus cereus MLST database; default=True')
parser.add_argument('-r', '--rpoB', nargs='?', default=True, help='Optional argument, True or False; perform rpoB typing; default=True')
parser.add_argument('-p', '--panC', nargs='?', default=True, help='Optional argument, True or False; perform panC typing; default=True')
parser.add_argument('-s', nargs='?', default=False, help='Optional argument, True or False; BLAST 16s DNA sequence; not recommended for inferring species or pathogenicity; default=False')
parser.add_argument('--spades_m', nargs='?', default=250, help='Optional argument for use with Illumina reads; integer; set SPAdes memory limit -m/--memory option in Gb; default is 250 Gb, the default for SPAdes')
parser.add_argument('--spades_t', nargs='?', default=16, help='Optional argument for use with Illumina reads; integer; set number of threads to use for SPAdes -t/--threads option; default is 16, the default for SPAdes')
parser.add_argument('--spades_k', nargs='?', default=77, help='Optional argument for use with Illumina reads; comma-separated list of integers; set k-mer sizes to use with SPAdes -k option; default is 77')
parser.add_argument('-v_db', '--virulence_database', nargs='?', default="aa", help='Optional argument for use with -v/--virulence option; specify virulence database to be used: nuc for nucleotide database or aa for amino acid database; default=aa')
parser.add_argument('-nuc_p', '--nucleotide_p', nargs='?',default=75, help='Optional argument for use with -v/--virulence option and nucleotide database -v_db nuc option; integer between 0 and 100; minimum percent nucleotide identity for virulence gene detection; default=75')
parser.add_argument('-nuc_q','--nucleotide_q', nargs='?',default=90, help='Optional argument for use with -v/--virulence option and nucleotide database -v_db nuc option; integer between 0 and 100; minimum percent coverage for virulence gene detection; default=90')
parser.add_argument('-aa_p','--amino_acid_p',nargs='?',default=50, help='Optional argument for use with -v/--virulence option and amino acid database -v_db aa option; integer between 0 and 100; minimum percent amino acid identity for virulence gene detection; default=50')
parser.add_argument('-aa_q','--amino_acid_q', nargs='?',default=70, help='Optional argument for use with -v/--virulence option and amino acid database -v_db aa option; integer between 0 and 100; minimum percent coverage for virulence gene detection; default=70')
parser.add_argument('-e', '--evalue', nargs='?', default=1e-5, help='Optional argument; float >= to 0; maximum blast e-value for a hit to be saved; default=1e-5')
parser.add_argument('--version', action="version", version='%(prog)s 1.0.0', help="Print version")

args=parser.parse_args()


# path to program btyper.py, as well as databases
btyper_path = os.path.realpath(__file__)


# split arguments provided by user
def arg_splitter(argin):
	return str(argin).strip("[']")

iarg=arg_splitter(args.input)
oarg=arg_splitter(args.output)
targ=arg_splitter(args.type)
varg=arg_splitter(args.virulence)
marg=arg_splitter(args.mlst)
rarg=arg_splitter(args.rpoB)
parg=arg_splitter(args.panC)
sarg=arg_splitter(args.s)
spadesm=arg_splitter(args.spades_m)
spadest=arg_splitter(args.spades_t)
spadesk=arg_splitter(args.spades_k)
dgarg=arg_splitter(args.draft_genome)
vdbarg=arg_splitter(args.virulence_database)
nucparg=arg_splitter(args.nucleotide_p)
nucqarg=arg_splitter(args.nucleotide_q)
aaparg=arg_splitter(args.amino_acid_p)
aaqarg=arg_splitter(args.amino_acid_q)
earg=arg_splitter(args.evalue)

# append to output directory argument, if necessary
if oarg[-1]!="/":
	oarg=oarg+"/"

# between_sections: prints a blank line between sections in final output file, if desried
def between_sections(finalfile_string):
	finalfile=open(finalfile_string,"a")
	print >> finalfile, ""
	finalfile.close()

# check if file exists and isn't empty
def file_check(input_file):
	try:
		if os.stat(input_file).st_size > 0:
			print input_file+" exits. Continuing..."
			return iarg
		else:
			print "Your input file {0} is empty! Please use a different input file.".format(input_file)
		sys.exit()
	except OSError:
		print "Your input file {0} does not exist. Please make sure your path and file name are correct, or specify a different input file.".format(input_file)
		sys.exit()

# dbparse: creates dictionary mapping virulence gene/mlst/rpoB/panC/16s id names from a given database to their respective sequences
def dbparse(db_path):
	argseq={}
	infile=open(db_path,"r")
	for record in SeqIO.parse(infile, "fasta"):
		seqid=str(record.description).strip()
		newid=re.sub('[^0-9a-zA-Z]+', '_', seqid)
		seqseq=str(record.seq).strip()
		seqlen=len(seqseq)
		if newid not in argseq.keys():
			argseq[newid]=seqseq
	infile.close()
	return argseq

# seqparse: creates dictionary that maps query sequence IDs to sequences
def seqparse(seq_path):
	newseq={}
	newinfile=open(seq_path,"r")
	for record in SeqIO.parse(newinfile, "fasta"):
		seqid=str(record.id).strip()
		newid=re.sub('[^0-9a-zA-Z]+', '_', seqid)
		seqseq=str(record.seq).strip()
		reduced_seq=seqseq.replace("N","")
		if newid not in newseq.keys():
			newseq[newid]=seqseq
	newinfile.close()
	return newseq

# make_blast_xml: takes sequences, makes blast xml file for each
def make_blast_xml(newseq,argdict,query_path,task,shorttask,evalue_thresh,pident_thresh,qcov_thresh):
	# for each sequence in the newseq dictionary
	for key in newseq.keys():
		# count the sequences
		counter=1
		# call query sequence queryseq
		queryseq=newseq[key]
		# get the length of the queryseq (genelen)
		genelen=len(queryseq)
		# create a temporary query file with the gene sequence and id
		queryfile=open(oarg+key.strip()+"_querytemp.fasta","a")
		print >> queryfile, ">"+key.strip()
		print >> queryfile, queryseq
		queryfile.close()
		# check if final results directory exists
		if os.path.isdir(oarg+"btyper_final_results"):
			print "Final results directory exists."
		# if not, make one
		else:
			os.system("mkdir "+oarg+"btyper_final_results")
		# check if isolatefiles directory exists within final results directory
		if not os.path.isdir(oarg+"btyper_final_results/isolatefiles"):
			os.system("mkdir "+oarg+"btyper_final_results/isolatefiles")
		rdir_root=oarg+"btyper_final_results/isolatefiles/"+key#.split(".fasta")[0]
		#print rdir_root
		# make a results directory for each sequence
		if not os.path.isdir(rdir_root+"_results"):
			os.system("mkdir "+rdir_root+"_results")
		rdir=rdir_root+"_results/"
		# create a database out of the temporary query file with the individual sequence
		if os.path.isfile(oarg+key.strip()+"_querytemp.fasta.nsq"):
			print "Database already exists."
		else:
			os.system("makeblastdb -in "+oarg+key.strip()+"_querytemp.fasta -dbtype nucl")
		out=rdir_root+"_blast_results.xml"
		# create final file named for each isolate if necessary
		finalfile_string=key.strip()+"_final_results.txt"
		if not os.path.isfile(oarg+"btyper_final_results/"+finalfile_string):
			finalfile=open(oarg+"btyper_final_results/"+finalfile_string,"a")
			print >> finalfile, "BTyper Results for "+key
			print >> finalfile, ""
			finalfile.close()
		if not os.path.isdir(oarg+"btyper_final_results/genefiles"):
			os.system("mkdir "+oarg+"btyper_final_results/genefiles")
		# open the final file for the isolate
		finalfile=open(oarg+"btyper_final_results/"+finalfile_string,"a")
		# print task-specific header to final output file
		if shorttask=="virulence":
			# print the task
			print >> finalfile, task.strip()
			print >> finalfile, ""
			table_header=["Hit #", "Virulence Gene Name","E-Value","Percent (%) Identity","Percent (%) Coverage"]
			print >>finalfile, "\t".join([t for t in table_header])
		elif shorttask=="rpoB":
			# print the task
			print >> finalfile, task.strip()
			print >> finalfile, ""
			print >> finalfile, "Predicted rpoB Allelic Type"+"\t"+"Percent (%) Identity"+"\t"+"Percent (%) Coverage"
		elif shorttask=="panC":
			print >> finalfile, task.strip()
			print >> finalfile, ""
			table_header=["panC Group Name","Closest Typed Strain","Percent (%) Identity","Percent (%) Coverage"]
                        print >>finalfile, "\t".join([t for t in table_header])
		elif shorttask=="16s":
			print >> finalfile, task.strip()
			print >> finalfile, ""
			print >> finalfile, "Predicted 16s Type"+"\t"+"Percent (%) Identity"+"\t"+"Percent (%) Coverage"
		# close the final file for the isolate
		finalfile.close()
		# run blast with query=database of sequences, db = querytemp.fasta with single isolate sequence, out = outdir, outfmt = xml
		if query_path==btyper_path+"seq_virulence_db/virul_aa_db.fasta":
			cline=NcbitblastnCommandline(query=str(query_path),db=oarg+key.strip()+"_querytemp.fasta",out=str(out),outfmt='"5"')
		else:
			cline=NcbiblastnCommandline(query=str(query_path),db=oarg+key.strip()+"_querytemp.fasta",out=str(out),outfmt='"5"')
		stdout, stderr = cline()
		#remove temporary temporary fasta files
		os.system("rm "+oarg+key.strip()+"_querytemp.fasta*")
		# parse blast xml files according to task
		if shorttask=="virulence":
			# parse_blast_xml_seq: parses blast xml file for virulence task (getting a list of all genes present in a sequence)
			parse_blast_xml_seq(xml=out,rdir=rdir,rdir_root=rdir_root,key=key,counter=counter,shorttask=shorttask,finalfile_string=finalfile_string,evalue_thresh=evalue_thresh,pident_thresh=pident_thresh,qcov_thresh=qcov_thresh)
		else:
			# parse_panC: parses blast xml for non-virulence tasks (getting the "top hit", rather than a list of all genes)
			parse_panC(xml=out,rdir=rdir,rdir_root=rdir_root,key=key,counter=counter,shorttask=shorttask,finalfile_string=finalfile_string,evalue_thresh=evalue_thresh,pident_thresh=pident_thresh,qcov_thresh=qcov_thresh)

# dictionary_add: create a new entry for a dictionary, or append if it already exists as a key
def dictionary_add(dict_name,key,blast):
	if key not in dict_name.keys():
		dict_name[key]=[]
		dict_name[key].append(blast)
	else:
		dict_name[key].append(blast)

# print_virulence_gene_report: prints isolatefiles and genefiles
def print_virulence_gene_report(filtres,maxgene, genequery,key,rdir,shorttask):
	print >> filtres, "\t".join([str(m).strip() for m in maxgene])
	genequery=str(genequery).strip()
	if shorttask=="virulence":
	###genequery=re.sub('[^0-9a-zA-Z]+', '_', genequery)
		genequery=genequery.split("|")[0]
	else:
		genequery=re.sub('[^0-9a-zA-Z]+', '_', genequery)
	genequery=genequery.strip()
	if shorttask=="virulence":
		genequery_short=genequery
	elif shorttask=="panC":
		genequery_short="panC"
	elif shorttask=="mlst":
		genequery_short=genequery.split("_")[0]
	elif shorttask=="rpoB":
		genequery_short="rpoB"
	elif shorttask=="16s":
		genequery_short="16s"
	# print to genefile (genefile=1 per gene, with all isolates' sequences)
	if shorttask!="virulence":
		genefile=open(oarg+"btyper_final_results/genefiles/"+genequery_short+"_genefile.fasta","a")
		print >> genefile, ">"+key+"___"+genequery
		print >> genefile, str(maxgene[14]).strip()
		genefile.close()
	# print to isolatefile (isolatefile=1 per isolate, with all genes detected per isolate)
	isolatefile=open(rdir+key+"_"+shorttask+"_sequences.fasta","a")
	print >> isolatefile, ">"+genequery+"___"+key
	print >> isolatefile, str(maxgene[14]).strip()
	isolatefile.close()	

# print_final_report: print to the final report file i.e. BTyper's final output file
def print_final_report(finalfile_string, shorttask, maxgene):
	finalout=open(oarg+"btyper_final_results/"+finalfile_string,"a")
	if shorttask=="virulence":
		gene_row=[maxgene[0],maxgene[3],maxgene[5],maxgene[8],maxgene[7],maxgene[6],maxgene[14]]
		print >> finalout, "\t".join([str(gr).strip() for gr in gene_row])
		#finalout.close()
		#prune_alleles(finalpath=finalfile_string)
		#finalfile=open(oarg+"btyper_final_results/"+finalfile_string,"a")	
	elif shorttask=="rpoB":
		rpotype=maxgene[3]
		rpotype=str(rpotype).replace("_","|")
		if float(maxgene[8])<100:
			pid=str(maxgene[8]).strip()+"*"
		else:
			pid=str(maxgene[8]).strip()
		if float(maxgene[7])<100:
			qid=str(maxgene[7]).strip()+"*"
		else:
			qid=str(maxgene[7]).strip()
		if "*" in pid or "*" in qid:
			rpotype=rpotype.strip()+"*"
		gene_row=[rpotype,pid,qid]
		print >> finalout, "\t".join([str(gr).strip() for gr in gene_row])
		if "*" in rpotype:
			print >> finalout, "*No rpoB allele in the current database matches with 100% identity and coverage."
	elif shorttask=="panC":
		pan=str(maxgene[3]).strip()
		panpid=maxgene[8]
		panqid=maxgene[7]
		if float(panpid.strip())<75:
			pan="None*___None"
		if float(panpid.strip())<90 and float(panpid.strip())>75:
			pan="?*___?" 
		pan1=pan.split("___")[0]
		pan2=pan.split("___")[1]
                print >> finalout, pan1.strip()+"\t"+pan2.strip()+"\t"+panpid.strip()+"\t"+panqid.strip()
		if pan.strip()=="None*___None":
			print >> finalout, "*No panC gene could be detected in your sequence at > 75% identity. Your panC allele may not be significantly associated with the Bacillus cereus group."
		if pan.strip()=="?*___?":
			print >> finalout, "*A panC clade could not be determined for your isolate."
	elif shorttask=="16s":
		stype=maxgene[3]
		if float(maxgene[8])<97:
			pid=str(maxgene[8]).strip()+"*"
		else:
			pid=str(maxgene[8]).strip()
		if float(maxgene[7])<100:
			qid=str(maxgene[7]).strip()+"*"
		else:
			qid=str(maxgene[7]).strip()
		if "*" in pid or "*" in qid:
			stype=stype.strip()+"*"
		gene_row=[stype,pid,qid]
		print >> finalout, "\t".join([str(gr).strip() for gr in gene_row])
		if "*" in stype:
			print >> finalout, "*No 16s gene in the current database matches with >97% identity and/or 100% coverage. Your isolate may not be a member of the Bacillus cereus group."
	finalout.close()

# get_st: parse mlst output for final report file; gives ST from database of STs and ATs
def get_st(mlst_infile, st_file, finalfile_string,mlst_genes):
	mlstdict={}	
	lowconf=[]
	for line in mlst_infile:
		splits=line.split("\t")
		if targ:
			if "alignment_title" not in line:
				alleleid=splits[3]
				geneid=str(alleleid).split("_")[0]
				at=str(alleleid).split("_")[1]
				mlstdict[str(geneid).strip()]=str(at).strip()
				pid=splits[8]
				qid=splits[7]	
				if float(pid.strip())<100 or float(qid.strip())<100:
					lowconf.append(geneid.strip())
	
	mlst_infile.close()
	stdict={}
	stfile=open(st_file,"r")
	for line in stfile:
		splits=line.split("\t")
		st=splits[0]
		all_at=splits[1:8]
		dictionary_add(stdict,str(st).strip(),[a.strip() for a in all_at])
	stfile.close()
	st_order=[s for i in stdict["ST"] for s in i]
	myst=[]
	for s in st_order:
		if s in mlstdict.keys():
			myst.append(mlstdict[s])
		else:
			myst.append("?")

	
	finalfile=open(finalfile_string,"a")
	print >> finalfile, "Predicted MLST Profile:"
	print >> finalfile, ""
	print >> finalfile, "ST"+"\t"+"\t".join([str(mg).split(".fas")[0] for mg in mlst_genes])
	if "?" in myst:
		finalmm=[]
		for mm,tt in zip(myst,mlst_genes):
			testmlst=tt.split(".fas")[0]
			if testmlst.strip() in lowconf:
				finalmm.append(mm+"*")
			else:
				finalmm.append(mm)
		finalfile=open(finalfile_string,"a")
		if "*" in str(finalmm):	
			print >> finalfile, "?*"+"\t"+"\t".join([str(m).strip() for m in finalmm])
			print >> finalfile, "*No allele in the current database matches with 100% identity and coverage."
		else:
			print >> finalfile, "?"+"\t"+"\t".join([str(m).strip() for m in finalmm])
		print >> finalfile, ""
		finalfile.close()
	else:
		confirmst=0
		for stkey in stdict.keys():
			st_match2=[s for i in stdict[stkey] for s in i]
			if myst==st_match2:
				confirmst=1
				finalmm=[]
				print "Found matching ST."
				for mm,tt in zip(myst,mlst_genes):
					testmlst=tt.split(".fas")[0]
					if testmlst.strip() in lowconf:
						finalmm.append(mm+"*")
					else:
						finalmm.append(mm)
						
				if "*" in str(finalmm):
					stkey=stkey.strip()+"*"
				finalfile=open(finalfile_string,"a")
				print >> finalfile, stkey.strip()+"\t"+"\t".join([str(m).strip() for m in finalmm])
				if "*" in stkey:
					print >> finalfile, "*No allele in the current database matches with 100% identity and coverage."
				print >> finalfile, ""
				finalfile.close()
		if confirmst==0:
			finalmm=[]
			for mm,tt in zip(myst,mlst_genes):
				testmlst=tt.split(".fas")[0]
				if testmlst.strip() in lowconf:
					finalmm.append(mm+"*")
				else:
					finalmm.append(mm)
			if "*" in str(finalmm):
				stkey="?*"
			else:
				stkey="?"
			finalfile=open(finalfile_string,"a")
			print >> finalfile, stkey.strip()+"\t"+"\t".join([str(m).strip() for m in finalmm])
			if "*" in stkey:
				print >> finalfile, "*No allele in the current database matches with 100% identity and coverage."
			print >> finalfile, ""
			finalfile.close()
			
			
		

# parse_blast_xml_seq: parse the blast xml output into results file
def parse_blast_xml_seq(xml,rdir,rdir_root,key,counter,shorttask,finalfile_string,evalue_thresh,qcov_thresh,pident_thresh):
	result_handle=open(xml)
	blast_records=NCBIXML.parse(result_handle)
	# open file in individual isolate results directory
	filtres=open(rdir+key+"_"+shorttask+"_results.txt","a")
	# print the header to the results file
	headerlist=["hit_number","alignment_title","query_id",shorttask+"_gene_query","alignment_length","evalue","blast_bitscore","query_coverage","percent_idenity",shorttask+"_gene_start",shorttask+"_gene_end","genome_sequence_start","genome_sequence_end",shorttask+"_gene_sequence","genome_sequence","match_sequence"]
	print >> filtres, "\t".join([h for h in headerlist])
	# loop through each record in the blast output
	for record in blast_records:
		for alignment in record.alignments:
			# create a dictionary for HSPs
			hspdict={}
			# for each hsp in the alignment
			for hsp in alignment.hsps:
				#query coverage (qcov) can be calculated as one of the following:
				#qcov=float(float(len(hsp.sbjct))/float(len(hsp.query)))
				qcov=float(hsp.align_length)/float(record.query_length)*100		
				pident=float(hsp.identities)/float(hsp.align_length)*100
				# add information for each HSP if the e-value is lower than the threshold
				if hsp.expect<=evalue_thresh and qcov>=qcov_thresh and pident>=pident_thresh:
					genefacts=[counter,alignment.title,record.query_id,record.query,hsp.align_length,hsp.expect,hsp.bits,qcov,pident,hsp.query_start,hsp.query_end,hsp.sbjct_start,hsp.sbjct_end,hsp.query,hsp.sbjct,hsp.match]
					dictionary_add(hspdict,str(record.query).strip(),"\t".join([str(g).strip() for g in genefacts]))
					# add the start and end positions of the hsp to hspdict
			# loop through each hsp in the alignment
			for hkey in hspdict.keys():
				# get values (start and end locations of each hsp, bit score)
				hsps=hspdict[hkey]
				# if more than 1 hsp exists (i.e. there is a chance that overlapping hsps might exist)
				if len(hsps)>1:
					# create list to check for overlapping HSPs
					allnums=[]
					# loop through each HSP
					for h in hsps:
						hsplits=h.split("\t")
						# get genome start position
						gstart=hsplits[11]
						gstart=int(gstart.strip())
						# get genome end position
						gend=hsplits[12]
						gend=int(gend.strip())
						if gstart<gend:
							# append sequence of base numbers covered by HSP to list
							allnums.append(range(gstart,gend+1,1))
						else:
							allnums.append(range(gend, gstart+1,1))
					# get the intersection of genome positions when more than one HSP exists in the genome
					inset=set.intersection(*map(set,allnums))
					# if an intersection exists (i.e. there are overlapping HSPs)
					if inset:
						# loop through HSPs, keeping track of HSP with highest bit score
						maxbits=0
						for h in hsps:
							hsplits=h.split("\t")
							# get bit score of current HSP
							testbits=hsplits[6]
							# if bit score of current HSP is larger than current maximum bit score:
							if float(testbits.strip())>maxbits:
								# the current bit score becomes the max
								maxbits=float(testbits.strip())
								# the current dictionary information becomes the best gene information
								maxgene=hsplits
								# the current gene query becomes the best gene query
								genequery=maxgene[3]
						# print the best-scoring HSP to the isolatefiles and genefiles
						counter=counter+1
						print_virulence_gene_report(filtres=filtres, maxgene=maxgene, genequery=genequery, key=key, rdir=rdir, shorttask=shorttask)
						print_final_report(finalfile_string=finalfile_string, shorttask=shorttask, maxgene=maxgene)
					# if no intersection exists, i.e. there are multiple HSPs, but they are found in different regions of the genome:
					else:
						# loop through each nonoverlapping HSP
						for h in hsps:
							maxgene=h.split("\t")
							genequery=maxgene[3]
							# print each HSP to the appropriate isolatefile/genefile
							counter=counter+1
							print_virulence_gene_report(filtres=filtres, maxgene=maxgene, genequery=genequery, key=key, rdir=rdir, shorttask=shorttask)
							print_final_report(finalfile_string=finalfile_string, shorttask=shorttask, maxgene=maxgene)

				# if there is only 1 HSP
				else:
					for h in hsps:
						maxgene=h.split("\t")
						genequery=maxgene[3]
						#print it to the appropriate isolatefile/genefile
						counter=counter+1
						print_virulence_gene_report(filtres=filtres, maxgene=maxgene, genequery=genequery, key=key, rdir=rdir, shorttask=shorttask)
						print_final_report(finalfile_string=finalfile_string, shorttask=shorttask, maxgene=maxgene)
	filtres.close()
	prune_alleles(finalpath=finalfile_string)

# parse_panC: parse blast output for panC task 
def parse_panC(xml,rdir,rdir_root,key,counter,shorttask,finalfile_string,evalue_thresh,qcov_thresh,pident_thresh):
	counter=0
        result_handle=open(xml)
        blast_records=NCBIXML.parse(result_handle)
        # open file in individual isolate results directory
        filtres=open(rdir+key+"_"+shorttask+"_results.txt","a")
        # print the header to the results file
        headerlist=["hit_number","alignment_title","query_id",shorttask+"_gene_query","alignment_length","evalue","blast_bitscore","query_coverage","percent_idenity",shorttask+"_gene_start",shorttask+"_gene_end","genome_sequence_start","genome_sequence_end",shorttask+"_gene_sequence","genome_sequence","match_sequence"]
        print >> filtres, "\t".join([h for h in headerlist])
        # loop through each record in the blast output
	maxp=0 
        for record in blast_records:
                for alignment in record.alignments:
                        # for each hsp in the alignment
                        for hsp in alignment.hsps:
                                #qcov=float(float(len(hsp.sbjct))/float(len(hsp.query)))
                                qcov=float(hsp.align_length)/float(record.query_length)*100 
                                pident=float(hsp.identities)/float(hsp.align_length)*100
				testbits=float(hsp.bits)
                                # add information for each HSP if the e-value is lower than the threshold
                                if hsp.expect<=evalue_thresh and qcov>=qcov_thresh and testbits>=pident_thresh and testbits>maxp:
					maxp=testbits
                                        genefacts=[counter,alignment.title,record.query_id,record.query,hsp.align_length,hsp.expect,hsp.bits,qcov,pident,hsp.query_start,hsp.query_end,hsp.sbjct_start,hsp.sbjct_end,hsp.query,hsp.sbjct,hsp.match]
        genefacts="\t".join([str(g).strip() for g in genefacts])
	maxgene=genefacts.split("\t")
	genequery=maxgene[3]
        #print it to the appropriate isolatefile/genefile
        counter=counter+1
        print_virulence_gene_report(filtres=filtres, maxgene=maxgene, genequery=genequery, key=key, rdir=rdir, shorttask=shorttask)
	print_final_report(finalfile_string=finalfile_string, shorttask=shorttask, maxgene=maxgene)
        filtres.close()

# assemble_reads: assemble genome from reads
def assemble_reads_pe(forward,reverse):
	os.system("spades.py -k "+spadesk+" --careful -1 "+forward+" -2 "+reverse+" -o "+oarg+"spades_assembly -t "+spadest+" -m "+spadesm)
def assemble_reads_se(reads):
	os.system("spades.py -k "+spadesk+" --careful -s "+reads+" -o "+oarg+"spades_assembly -t "+spadest+" -m "+spadesm)

def fastq_check(reads):
	if reads.endswith(".fastq"):
		os.system("gzip reads")
		reads=reads.strip()+".gz"
	if reads.endswith(".fastq.gz"):
		return()
	else:
		print "It looks like your file "+reads.strip()+" does not end with either '.fastq.gz' or '.fastq'. Please make sure that you are supplying reads in fastq or fastq.gz format, and, if necessary, please rename your file with the appropriate extension."
		sys.exit()

def sra_check(sra):
	if sra.endswith(".sra"):
		return()
	else:
		print "It looks like your file "+reads.strip()+" does not end with '.sra'. Please make sure that you are supplying sequencing data in sra format, and, if necessary, please rename your file with the extension '.sra'."

def prune_alleles(finalpath):
	f=oarg+"btyper_final_results/"+finalpath	
	print "Pruning alleles for "+f+"..."
	os.system("mv "+f+" "+f+"_temporary.txt")
	old=open(f+"_temporary.txt","r")
	new=open(f,"a")
	bitdict={}
	for o in old:
		if "BTyper Results for " in o:
			print >> new, o.strip()
			print >> new, ""
		if "Predicted Virulence " in o and ":" in o:
			print >> new, o.strip()
			print >> new, ""
		if "Hit #" in o and "Virulence Gene Name" in o:
			print >> new, o.strip()
		if "|" in o and "___" not in o and "rpoB|AT" not in o:
			print >> new, o.strip()
		elif "|" in o and "___" in o and "rpoB|AT" not in o:	
			splits=o.split("\t")
			gname=splits[1]	
			gpref=gname.split("___")[0]
			gsuff=gname.split("___")[1]
			testbits=splits[5]
			if gpref.strip() not in bitdict:
				bitdict[gpref.strip()]=[]
				bitdict[gpref.strip()].append(gsuff.strip())
				for s in splits[2:]:
					bitdict[gpref.strip()].append(str(s).strip())	
			else:
				maxlist=bitdict[gpref.strip()]
				maxbits=maxlist[4]
				if float(testbits.strip())>float(maxbits.strip()):
					bitdict[gpref.strip()]=[]
					bitdict[gpref.strip()].append(gsuff.strip())
					for s in splits[2:]:
						bitdict[gpref.strip()].append(str(s).strip())	
	old.close()
	counter=0
	for key in bitdict.keys():
		counter=counter+1
		finalg=bitdict[key]
		print >> new, str(counter)+"\t"+"\t".join([str(fg).strip() for fg in finalg[:-2]])
		shortgene=finalg[0].split("|")[0]
		genefile=open(oarg+"btyper_final_results/genefiles/"+shortgene.strip()+"_genefile.fasta","a")
		isoname=f.split("/")[-1]
		print >> genefile, ">"+isoname.strip()+"___"+shortgene
		print >> genefile, finalg[-1].strip()
		genefile.close()
	new.close()
	os.system("rm "+f+"_temporary.txt")
				

# run program based on user input
if targ!="seq":
	os.system("mkdir "+oarg+"spades_assembly")
	if targ=="pe":
		spacecount=iarg.count(" ")
		if spacecount>1:
			print "It looks like there may be whitespace in one or more of your file names. Please rename your file(s), and try again."
			sys.exit()
		elif spacecount==0:
			print "It looks like you've chosen the paired-end reads (-t pe) option but have only supplied one file of reads. Please supply one forward reads file and one reverse reads file separated by a space, or select a different option for -t (se for single-end reads, sra for files in sra format, or sra-get to download sequencing data in SRA format)."
			sys.exit()
		else:
			forward_reads=iarg.split(" ")[0]
			reverse_reads=iarg.split(" ")[1]
			forward_reads=re.sub("[,']","",forward_reads)
			reverse_reads=re.sub("[,']","",reverse_reads)
			fastq_check(forward_reads)
			fastq_check(reverse_reads)
			prefix=forward_reads.split(".fastq.gz")[0]
			if "/" in prefix:
				prefix=prefix.split("/")[-1]
			assemble_reads_pe(forward=forward_reads,reverse=reverse_reads)
	elif targ=="se":
		sreads=iarg.strip()
		fastq_check(sreads)
		prefix=sreads.split(".fastq.gz")[0]
		if "/" in prefix:
			prefix=prefix.split("/")[-1]
		assemble_reads_se(reads=sreads)
	elif targ=="sra-get" or targ=="sra":
		sra=iarg.strip()
		if targ=="sra-get":
			os.system("prefetch -v "+sra)
			os.system("fastq-dump --outdir "+oarg+" --split-files --gzip ~/ncbi/public/sra/"+sra+".sra")
		elif targ=="sra":
			sra=iarg.split(".sra")[0]
			os.system("fastq-dump --outdir "+oarg+" --split-files --gzip "+iarg.strip())
		forward_reads=oarg+sra.strip()+"_1.fastq.gz"
		fastq_check(forward_reads)
		prefix=forward_reads.split(".fastq.gz")[0]
		if "/" in prefix:
			prefix=prefix.split("/")[-1]
		if os.path.exists(oarg+sra+"_2.fastq.gz"):
			reverse_reads=oarg+sra+"_2.fastq.gz"
			fastq_check(reverse_reads)
			assemble_reads_pe(forward=forward_reads,reverse=reverse_reads)
		else:
			assemble_reads_se(reads=forward_reads)
	os.system("cp "+oarg+"spades_assembly/contigs.fasta "+oarg+prefix.strip()+"_spades_assembly.fasta")
	contigs=open(oarg+prefix.strip()+"_spades_assembly.fasta","r")
	pseudo=open(oarg+prefix.strip()+"_pseudochrom.fasta","a")
	print >> pseudo, ">"+prefix.strip()
	print >> pseudo, "".join(line.strip() for line in contigs if ">" not in line)
	contigs.close()
	pseudo.close()
	targ="seq"
	iarg=oarg+prefix.strip()+"_pseudochrom.fasta"




# if the user inputs a fasta file
if targ=="seq":
	print "-t seq has been selected. Expecting a fasta file."
	# check if the file is a fasta
	file_check(iarg)
	# if genome is a draft genome
	if dgarg=="True":
		prefix=iarg.split(".")[0]
		if "/" in prefix:
			prefix=prefix.split("/")[-1]
		contigs=open(iarg,"r")
		pseudo=open(oarg+prefix.strip()+"_pseudochrom.fasta","a")
		print >> pseudo, ">"+prefix.strip()
		print >> pseudo, "".join(line.strip() for line in contigs if ">" not in line)
		contigs.close()
		pseudo.close()
		iarg=oarg+prefix.strip()+"_pseudochrom.fasta"
	# create a dictionary from the user-supplied fasta, mapping ids to sequences 
	dictionaries=seqparse(iarg)

	# if performing virulence typing
	if varg=="True":
		# define database
		query_paths=[]
		if vdbarg=="nuc":
			query_paths.append(btyper_path+"seq_virulence_db/bc_virul_nt.fasta")
		if vdbarg=="aa":
			query_paths.append(btyper_path+"seq_virulence_db/virul_aa_db.fasta")
		for query_path in query_paths:
			# map db ids to sequences
			mydb=dbparse(query_path)
			#task="Predicted Virulence Genes:"
			shorttask="virulence"
			# define minimum e-value for virulence genes
			evalue_thresh=float(earg)
			if query_path==btyper_path+"seq_virulence_db/bc_virul_nt.fasta":
				task="Predicted Virulence Genes:"
				# define minimum percent identity for virulence genes
				pident_thresh=int(nucparg)
				# define minimum query coverage for virulence genes
				qcov_thresh=int(nucqarg)
			elif query_path==btyper_path+"seq_virulence_db/virul_aa_db.fasta":
				task="Predicted Virulence Proteins:"
				# define minimum percent identity for virulence genes
				pident_thresh=int(aaparg)
				# define minimum query coverage for virulence genes
				qcov_thresh=int(aaqarg)
			# run virulence typing
			make_blast_xml(newseq=dictionaries,argdict=mydb,query_path=query_path,task=task,shorttask=shorttask,evalue_thresh=evalue_thresh,pident_thresh=pident_thresh,qcov_thresh=qcov_thresh)
			# print spacing between under final results
			#prune_alleles(finalpath=oarg+"btyper_final_results/")
			for f in glob.glob(oarg+"btyper_final_results/*_final_results.txt"):
				between_sections(finalfile_string=f)
			#prune_alleles(oarg+"btyper_final_results")
	
	# if performing panC typing:
	if parg=="True":
		# define database
		query_path=btyper_path+"seq_panC_db/panC.fasta"
		mydb=dbparse(query_path)
		task="Predicted panC Clade Designation:"
		shorttask="panC"
		evalue_thresh=float(earg)
		pident_thresh=0
		qcov_thresh=0
		# run panC typing, if panC gene is detected
		try:
			make_blast_xml(newseq=dictionaries,argdict=mydb,query_path=query_path,task=task,shorttask=shorttask,evalue_thresh=evalue_thresh,pident_thresh=pident_thresh,qcov_thresh=qcov_thresh)
			for f in glob.glob(oarg+"btyper_final_results/*_final_results.txt"):
				between_sections(finalfile_string=f)
		except UnboundLocalError:
			print "No sequences found for "+shorttask
	
	# if performing mlst:	
	if marg=="True":
		mlst_genes=["glp.fas","gmk.fas","ilv.fas","pta.fas","pur.fas","pyc.fas","tpi.fas"]
		# get best-matching AT for each MLST gene
		for mlst in mlst_genes:
			query_path=btyper_path+"seq_mlst_db/"+mlst
			mydb=dbparse(query_path)
			task="Predicted MLST Profile:"
			shorttask="mlst"
			evalue_thresh=float(earg)
			pident_thresh=0
			qcov_thresh=0
			# run AT for each gene, if sequence deteced
			try:
				make_blast_xml(newseq=dictionaries,argdict=mydb,query_path=query_path,task=task,shorttask=shorttask,evalue_thresh=evalue_thresh,pident_thresh=pident_thresh,qcov_thresh=qcov_thresh)
			except UnboundLocalError:
				print "No sequences found for "+mlst
		# loop through isolatefiles
		for root, dirs, files in os.walk(oarg+"btyper_final_results/isolatefiles/"):
			for d in dirs:
				dirroot=d.split("_results")[0]
				# open mlst results file, and get ST from ATs
				try:	
					newf=open(oarg+"btyper_final_results/isolatefiles/"+d+"/"+dirroot.strip()+"_mlst_results.txt","r")
					finalfile_string=oarg+"btyper_final_results/"+dirroot.strip()+"_final_results.txt"
					ff=open(finalfile_string,"r")
					flines=ff.readlines()	
					if not any("Predicted MLST Profile" in fl.strip() for fl in flines):
						get_st(mlst_infile=newf,st_file=btyper_path+"seq_mlst_db/b_cereus_mlst_db.txt",finalfile_string=finalfile_string,mlst_genes=mlst_genes)
				except IOError:
					print "No sequences found for "+shorttask

	# if performing rpoB typing:
	if rarg=="True":
		# define database
		query_path=btyper_path+"seq_rpoB_db/rpobdatabase08122015.fa"
		mydb=dbparse(query_path)
		task="Predicted rpoB Allelic Type:"
		shorttask="rpoB"
		evalue_thresh=float(earg)
		pident_thresh=0
		qcov_thresh=0
		# perform rpoB typing, if gene is present
		try:
			make_blast_xml(newseq=dictionaries,argdict=mydb,query_path=query_path,task=task,shorttask=shorttask,evalue_thresh=evalue_thresh,pident_thresh=pident_thresh,qcov_thresh=qcov_thresh)
			for f in glob.glob(oarg+"btyper_final_results/*_final_results.txt"):
				between_sections(finalfile_string=f)
		except UnboundLocalError:
			print "No sequences found for "+shorttask
	
	# if performing 16s typing:
	if sarg=="True":
		# define database
		query_path=btyper_path+"seq_16s_db/b_cereus_group_16s_db.fasta"
		mydb=dbparse(query_path)
		task="Predicted 16s Type"
		shorttask="16s"
		evalue_thresh=float(earg)
		pident_thresh=0
		qcov_thresh=0
		# perform 16s typing, if gene is present
		try:
			make_blast_xml(newseq=dictionaries,argdict=mydb,query_path=query_path,task=task,shorttask=shorttask,evalue_thresh=evalue_thresh,pident_thresh=pident_thresh,qcov_thresh=qcov_thresh)
			for f in glob.glob(oarg+"btyper_final_results/*_final_results.txt"):
				between_sections(finalfile_string=f)
		except UnboundLocalError:
			print "No sequences found for "+shorttask
print "Typing complete...how neat is that?"
print ""
print "Thank you for using BTyper! For more fun, take your output files to BMiner, BTyper's companion application for data aggregation and visualization."
print ""
print "To cite BTyper and/or BMiner, please use the following:"
print "Carroll, Laura M., Jasna Kovac, Rachel A. Miller, Martin Wiedmann. 2017. Rapid, high-throughput identification of anthrax-causing and emetic Bacillus cereus group genome assemblies using BTyper, a computational tool for virulence-based classification of Bacillus cereus group isolates using nucleotide sequencing data. Submitted to Applied and Environmental Microbiology."


