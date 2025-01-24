#!/usr/bin/env python

import re 
import argparse
import os
from operator import itemgetter
from Bio import SeqIO
import gzip
import subprocess
import time

pssm_dict = {380458: 'C', 380456: 'C', 380461: 'C', 380454: 'C', 380470: 'C', 380465: 'C', 380466: 'C', 380455: 'C', 334202: 'C', 380457: 'E', 275325: 'PT', 376382: 'cAT', 379688: 'DH', 214837: 'DH', 239531: 'DH', 224941: 'DH', 214838: 'AT', 223408: 'AT', 214836: 'KS', 238429: 'KS', 238425: 'KS', 273826: 'KS', 238383: 'KS', 215161: 'KS', 274452: 'KS', 223381: 'KS', 238428: 'KS', 238430: 'KS', 238427: 'KS', 238426: 'KS', 374347: 'SAT', 273136: 'PPT', 223807: 'PPT', 369777: 'MT', 369778: 'MT', 372616: 'MT', 379312: 'MT', 316372: 'MT', 100107: 'MT', 214839: 'MT', 273787: 'TR', 187546: 'TR', 214835: 'TE', 366397: 'TE', 223730: 'TE', 366166: 'TE', 341253: 'A', 366135: 'A', 214834: 'ACP', 177047: 'ACP', 373139: 'ACP', 376348: 'ACP', 375752: 'ACP', 214833: 'KR', 187653: 'KR', 227314: 'ER', 214840: 'ER', 176179: 'ER', 176231: 'ER', 176645: 'ER'}
'''
pssm_dict taken rules.json from Synthaser github.
Gilchrist, C. L., & Chooi, Y. H. (2021).
Synthaser: a CD-Search enabled Python toolkit for analysing domain architecture of fungal secondary metabolite megasynth (et) ases.
Fungal Biology and Biotechnology, 8(1), 1-19.

Each PSSM corresponds to a PKS/NRPS domain.'''

def andor(a, b):
    '''Fxn to be able to return if at least one is True'''
    if a and b:
        return True
    elif a or b:
        return True
    else:
        return False

def set_variables(args):
    '''
    Setting the variables required.
        Fasta file, output directory always necessary.
        Input file (of CDD blast results) //OR// path to blast database required to get blast results
            * if blast database is given, BLAST will be run on that database. This will increase run time.
        antiSMASH results or 'True' to run antiSMASH; optional
            * if True is given, antiSMASH will be run. This will increase runtime. 
        collapsed default is True to collapse results that are named the same and overlap. This helps with viewing in final output. 
        num_matches and distance_required are default set to 3 and 10,000bp. This can be changed if you want more or less stringent results. 
            * num_matches requires there to be at least a group of 3 BGC domain hits that are named differently at least once (no KS,KS,KS. but KS,AT,KS would pass.)
            * distance_required requires the group of 3 to be at least within 10,000 bp of each other. 
                ** This orders them by start position, then compares the end of 1 to the start of 2 and the start of 1 to the start of 2 and so on. Each hit only needs to be within 10kbp of the one before it, so the BGCs can easily be larger than 10kbp
    '''

    if not args.fasta:
        parser.error(f'ERROR: Please include a fasta file for your genome of choice with -f or --fasta.')
    if not args.output_directory:
        parser.error(f'ERROR: Please indicate an output directory for BiGEST results with -o or --output_directory')
    fasta_file = args.fasta
    output_directory = args.output_directory
    num_matches_required = args.num_matches
    distance_required = args.distance
    os.makedirs(output_directory,exist_ok=True)
    fasta_name, file_extension = os.path.splitext(os.path.basename(fasta_file))
    if fasta_file.endswith(".gz"):
        fasta_name, file_extension = os.path.splitext(os.path.basename(fasta_name))
    if str(args.antismash_genbank) == 'True':
        if args.gff3 != None:
            if os.path.exists(args.gff3):
                gff3 = args.gff3
            else:
                raise ValueError("ERROR: Path for gff3 does not exist. Please check the path. If you would not like to include a gff3 file for the antiSMASH run, do not use the -g or --gff3 flag.")
        else:
            gff3 = None
        antismash_gbk = run_antismash(fasta_file,output_directory,fasta_name,gff3)
    else:
        if args.gff3 != None:
            print('WARNING: gff3 provided will be ignored: antiSMASH has already been run.')
        antismash_gbk = args.antismash_genbank

    output_combined_gbk = os.path.join(output_directory, f"combined_{fasta_name}.gbk")
    if str(args.collapsed) == 'True':
        collapsed = True
    else:
        collapsed = False
    
    if not args.input:
        if not args.db:
            raise ValueError('ERROR: Input blast text file (--input) or blast DB (--db) must be given.')
        else:
            blast_data = run_blast(args.db,fasta_file,fasta_name,output_directory)
    else:
        blast_data = args.input

    return fasta_file,output_directory,fasta_name,num_matches_required,distance_required,antismash_gbk,output_combined_gbk,collapsed,blast_data

def run_blast(db,fasta_file,fasta_name,output_directory):
    '''
    Runs blast for the program. This will unzip and zip files as required and print the standard out and standard error to the terminal to allow for troubleshooting. 
    '''
    output_file = os.path.join(output_directory,f'{fasta_name}_rpstblastn.txt')
    start_time = time.time()
    if fasta_file.endswith(".gz"):
        print("Unzipping file. File will be zipped again after blast run")
        result = os.system(f'gunzip {fasta_file}')
        if result == 0:
            filename = fasta_file.split('.gz')[0]
            cmd = f'rpstblastn -num_threads 6 -evalue 1e-5 -query {filename} -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore stitle gaps qseq sseq sacc slen" >> {output_file}'
            result = subprocess.run(cmd,stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell=True, universal_newlines=True)
            print(result.stdout)
            print(result.stderr)
            print('Blast done. Zipping fasta file again.')
            result = os.system(f'gzip {filename}')
            if result!=0:
                print("WARNING: Error in zipping fasta file. Script will continue")
        else:
            raise ValueError("ERROR: Error in unzipping file. Please unzip file before running.")
    else:
        cmd = f'rpstblastn -num_threads 6 -evalue 1e-5 -query {fasta_file} -db {db} -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore stitle gaps qseq sseq sacc slen" >> {output_file}'
        result = subprocess.run(cmd,stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell=True, universal_newlines=True)
        print(result.stdout)
        print(result.stderr)
        print('Blast done.')
    if os.stat(output_file).st_size == 0:
        raise ValueError("ERROR: No CDD blast results. Please check the above blast stdout/stderr for any potential blast errors.")
    elapsed_time = round((time.time()-start_time),2)
    if elapsed_time >= 3600:
        elapsed_time = round(elapsed_time / 3600,2)
        text = f'Blast search took {elapsed_time} hours'
    elif elapsed_time >= 60: 
        elapsed_time = round(elapsed_time / 60,2)
        text = f'Blast search took {elapsed_time} minutes'
    else:
        text = f'Blast search took {elapsed_time} seconds'
    print(text)
    return output_file

def run_antismash(fasta_file,output_directory,fasta_name,gff3):
    '''This will run antiSMASH. It will print the standard out and standard error to the terminal for troubleshooting.'''
    start_time = time.time()
    output = os.path.join(output_directory,f'antismash_output',f'{fasta_name}')
    os.makedirs(output,exist_ok=True)
    if gff3 == None:
        cmd = f'antismash {fasta_file} --genefinding-tool glimmerhmm --fullhmmer --taxon fungi --output-dir {output} --output-basename {fasta_name}'
    else:
        cmd = f'antismash {fasta_file} --gene-finding-tool {gff3} --taxon fungi --output-dir {output} --output-basename {fasta_name}'
    result = subprocess.run(cmd,stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell=True, universal_newlines=True)
    print(result.stderr)
    print(result.stdout)
    antismash_gbk = os.path.join(output,f'{fasta_name}.gbk')
    try: 
        os.stat(antismash_gbk).st_size
        antismash_gbk=antismash_gbk
    except FileNotFoundError:
        print(f'WARNING: No antiSMASH results found at {antismash_gbk}.\n Please troubleshoot with the above antiSMASH output if you would still like to include it. \nThe BiGEST search will continue without antiSMASH results.')
        antismash_gbk = None
    elapsed_time = round((time.time()-start_time),2)
    if elapsed_time >= 3600:
        elapsed_time = round(elapsed_time / 3600,2)
        text = f'antiSMASH search took {elapsed_time} hours'
    elif elapsed_time >= 60: 
        elapsed_time = round(elapsed_time / 60,2)
        text = f'antiSMASH search took {elapsed_time} minutes'
    else:
        text = f'antiSMASH search took {elapsed_time} seconds'
    print(text)
    return antismash_gbk

########## processing blast output ####################

def access_data(blast_file_name,coverage_min,antismash_gbk,fasta_name):
    '''Get info from the blast.txt output and save the information in the class BiGEST and/or antiSMASH
    This will filter the BiGEST results based on coverage requirements. 
    This will filter antiSMASH results to have only NRPS and PKS results in the aSModule, aSDomain, and protoclusters.
    '''
    contigs_dict = dict()
    antismash_dict = dict()
    with open(blast_file_name, "r") as f: 
        for line in f:
            splits = line.split("\t")
            contig = splits[0]
            try:
                pssm_value = int(splits[14].split(":")[1].strip("\n"))
            except IndexError:
                pass
            if pssm_value in pssm_dict:
                if contig not in contigs_dict:
                    if check_coverage(line,coverage_min) == True:
                        contigs_dict[contig] = BiGEST(contig) ## create object, add to list
                        contigs_dict[contig].get_info(line,fasta_name) ## get important info from that object. 
                else:
                    if check_coverage(line,coverage_min)== True:
                        contigs_dict[contig].get_info(line,fasta_name) ## contig already seen, get import info from that object
    if antismash_gbk!=None:
        with open(antismash_gbk,"r") as my_file:
            for record in SeqIO.parse(my_file,"genbank"):
                name = record.description.split(' ')[0]
                for feature in record.features:
                    if 'aSModule'in feature.type:
                        if 'nrps' in ','.join(feature.qualifiers['domains']):
                            if name not in antismash_dict:
                                antismash_dict[name] = antiSMASH(name)
                                antismash_dict[name].get_anti_mod_info(feature,record)
                            else:
                                antismash_dict[name].get_anti_mod_info(feature,record)
                    if 'aSDomain' in feature.type:
                        if 'nrps' in feature.qualifiers['aSTool'][0]:
                            if name not in antismash_dict:
                                antismash_dict[name] = antiSMASH(name)
                                antismash_dict[name].get_anti_domain_info(feature,record)
                            else:
                                antismash_dict[name].get_anti_domain_info(feature,record)
                    if 'protocluster' in feature.type:
                        if 'NRPS' in feature.qualifiers['category'] or 'PKS' in feature.qualifiers['category']:
                            #print(name)
                            if name not in antismash_dict:
                                antismash_dict[name] = antiSMASH(name)
                                antismash_dict[name].get_anti_protocluster_info(feature,record)
                            else:
                                antismash_dict[name].get_anti_protocluster_info(feature,record)
                    
    ### keep only groups of {num_matches_required} within {distance} and remove the rest of the options
        keys_for_delete = set()
        for key,value in antismash_dict.items():
            if len(value.domain_info) < num_matches_required:
                keys_for_delete.add(key)
            else:
                value = value.are_they_close_enough()
                if value == None:
                    keys_for_delete.add(key)
        for k in keys_for_delete:
            del antismash_dict[k]
    return contigs_dict,antismash_dict

def check_coverage(line,coverage_min):  
    '''
    Used in access_data()
    This checks the BLAST hit coverage of the subject to ensure it covers a reasonable amount of the reference domain.
    '''
    line = line.split('\t')
    sseqid = line[1]
    sseqlength = int(line[15])
    if sseqlength == 0:
        return True
    length = int(line[3]) - int(line[11]) #length - gaps
    coverage = round(length/sseqlength,2)
    if coverage >= float(coverage_min):
        return True
    else:
        return False

class BiGEST:

    def __init__(self,name):
        '''this holds the information of a contig (record) in a fasta file. Information is only saved if the blast output is part of a BGC'''

        ##the data##
        self.file_name = fasta_name # the name of the file
        self.contig = name ## record.description
        self.info = [] ## for BiGEST
        self.updated_info = [] ## [start,end,bsr,line]
        self.full_start = '' ## 10K before where the BGC starts
        self.full_end = '' ## 10K after the BGC startsx

    ## the methods ##

    def get_info(self,line,fasta_name):
        '''Save info as needed'''
        split_lines = line.split('\t')
        start,end = int(split_lines[4]),int(split_lines[5])
        bsg = self.find_BSG(split_lines)
        ## 1 based start and end
        self.info.append([start,end,bsg,line])
        self.file_name = fasta_name
        return 

    def find_BSG(self,splits):
        ''' this will take the line and compare to the synthaser database to give the overall bisynthetic gene.'''
        pssm_value = int(splits[14].split(":")[1].strip("\n"))
        if pssm_value in pssm_dict:
            return pssm_dict[pssm_value]

    def are_they_close_enough(self,info):
        '''This will sort the information by start position, 
        then loop through and find the BGC groups that are close enough, with enough matches, 
        that are not all the same (ie KS,KS,KS)'''
        sorts = sorted(info,key = itemgetter(0)) ## info is [(start,end,BSG,line),(start,end,BSG,line)] for each contig
        current_group = [sorts[0]]
        current_end = max(current_group[0][0],current_group[0][1])
        groups_of_BGCs = [sorts[0][2]]
        if len(sorts)< num_matches_required:
            return None
        i = 0
        for obj in sorts[1:]:
            i += 1
            start = min(obj[0],obj[1])
            if abs(current_end - start) <= distance_required: ## they are within the boundary
                current_group.append(obj)
                if groups_of_BGCs[-1] != obj[2]:
                    groups_of_BGCs.append(obj[2]) ## pks
                current_end = max(obj[0],obj[1])
            else: ## they are too far apart
                if len(groups_of_BGCs)>= num_matches_required:
                    return current_group,i
                else:
                    groups_of_BGCs = [obj[2]]
                    current_group=[obj]
                    current_end = max(obj[0],obj[1])
        if len(groups_of_BGCs)>= num_matches_required:
            return current_group,i
        else:
            return None

    def joining(self):
        '''This joins together the domain hits that are exactly the same with exactly the same positions'''
        # updated_info is [(start,end,BSG,line),(start,end,BSG,line)] for each contig
        for item in self.updated_info:
            if item[0] > item[1]:
                item[0],item[1] = item[1],item[0]
                if len(item)==4:
                    item.append('-') 
            else:
                if len(item)==4:
                    item.append('+')
            ### now updated info has  [(start,end,BSG,line,strand),(start,end,BSG,line,strand)] for each contig
            ### from here, the start is always a lower number than the end. 
        sorted_objects = sorted(self.updated_info,key=itemgetter(2,0)) ## sort by BSG and then start location
        split_groups = []
        current_group = [sorted_objects[0]]
        current_min_start,current_max_end = current_group[0][0],current_group[0][1]
        more_descriptive = current_group[0][3].split('\t')[10].split(',')[1].strip() ## this is the stitle information from the orignial BLAST hit
        current_BSG_to_compare = f"{current_group[0][2]}_{more_descriptive}" ## I want to compare if the names are exactly the same for the joining
        current_strand = current_group[0][4]
        for obj in sorted_objects[1:]:
            more_descriptive = obj[3].split('\t')[10].split(',')[1].strip()
            next_BSG_to_compare = f"{obj[2]}_{more_descriptive}"
            next_min_start, next_max_end =  obj[0],obj[1]
            if current_BSG_to_compare == next_BSG_to_compare and andor(current_min_start==next_min_start,current_max_end==next_max_end) and obj[4] == current_strand:
                current_group.append(obj)
            else:
                split_groups.append(current_group)
                current_group = [obj]
                current_min_start,current_max_end = current_group[0][0],current_group[0][1]
                this_bgc = current_group[0][3].split('\t')[10].split(',')[1].strip()
                current_BSG_to_compare = f"{current_group[0][2]}_{this_bgc}"
        split_groups.append(current_group)
        return split_groups

    def joining_collapsed(self):
        '''This only happens if collpased == True at the start of the script.
        This will join together domain hits that overlap and are named the same by CDD BLAST'''
        for item in self.updated_info:
            if item[0] > item[1]:
                item[0],item[1] = item[1],item[0]
                if len(item)==4:
                    item.append('-')
            else:
                if len(item)==4:
                    item.append('+')
        sorted_objects = sorted(self.updated_info,key=itemgetter(2,0,3)) ## sort by BSG and then location, then by direction
        split_groups = []
        current_group = [sorted_objects[0]]
        current_min_start,current_max_end = current_group[0][0],current_group[0][1]
        current_strand = current_group[0][4]
        for obj in sorted_objects[1:]:
            if obj[2] == current_group[-1][2] and max(current_min_start, obj[0]) <= min(current_max_end,obj[1]) and obj[4] == current_strand:
                current_min_start=min(current_min_start,obj[0])
                current_max_end= max(current_max_end,obj[1])
                current_strand = current_group[0][4]
                current_group.append(obj)
            else:
                split_groups.append(current_group)
                current_group = [obj]
                current_min_start,current_max_end = current_group[0][0],current_group[0][1]
                current_strand = current_group[0][4]
        split_groups.append(current_group)
        return split_groups

    def get_full_ends(self,fasta_file):
        '''To get the start and end of what section of the contig the BGC is in. This helps with viewing to not have the whole contig, 
        especailly if that contig is incredibly long. This should then only display a few kbp instead of all million bases. 
        The start and end, in reference to the contig, will be noted in the genbank notes at the beginning of the record.'''
        distance_list = self.updated_info
        starts,ends = [],[]
        for item in distance_list:
            if item[0]>item[1]:
                start,end = item[1],item[0]
            else:
                start,end = item[0],item[1]
            starts.append(start)
            ends.append(end)
        seq_len = len(get_sequence(fasta_file,self.contig))
        object.full_start = min(starts) - 10000 if min(starts) >= 10000 else 0 ## add 10k on either side so that there's some buffer in the GBK
        object.full_end = max(ends)+10000 if seq_len >= (max(ends)+10000) else seq_len
        return

class antiSMASH:

    def __init__(self,name):
        '''this holds the information of a contig (record) in a fasta file, the antiSMASH results.
        Information is only saved if the blast output is part of a BGC'''

        ##the data##
        self.file_name = fasta_name # the name of the file
        self.contig = name ## record.description
        self.mod_info = []
        self.domain_info = []
        self.proto_info = []
        self.full_start = '' 
        self.full_end = ''

    def get_anti_mod_info(self,feature,record):
        '''This will get the module information. This adds in only the NRPS/PKS module hits to the contig's mod_info'''
        type = feature.qualifiers['type'][0] ##'pks or nrps'
        starts,ends,locations = list(),list(),list()
        typer = feature.type ##'asmodule'
        qualifiers = feature.qualifiers
        matches = re.findall(r'\d+', str(feature.location))
        for iter,item in enumerate(matches):
            if iter%2 == 0: ## the starts, at the even values
                starts.append(item)
            else: ## the ends, at the odd values. 
                ends.append(item)
        for iter,item in enumerate(starts):
            locations.append((item,ends[iter]))
        direction = str(feature.location).split('(')[1].split(')')[0]  
    
        ## antismash_start is needed because antismash will separate contigs if they are far enough apart and have multiple genbank results files per contig.
        ## By keeping track of an antismash start and end per contig, then we can make sure they're all adjusted correctly.
        antismash_start,antismash_end,is_cut = self.get_antismash_starts(record,len(record.seq))

        if self.full_start == '' or self.full_start > antismash_start:
            self.full_start = antismash_start
        else:
            if self.full_start > antismash_start:
                self.full_start = antismash_start
        if self.full_end == '':
            self.full_end = antismash_end
        else:
            if self.full_end < antismash_end:
                self.full_end == antismash_end
        if is_cut == False:
            antismash_start = 0

        self.mod_info.append([locations,direction,typer,qualifiers,type,antismash_start])     
        return 

    def get_anti_domain_info(self,feature,record):
        '''Pulling out domain hits from the antiSMASH records. This will only pull the NRPS/PKS domain hits to the contigs domain_info'''
        antismash_start,antismash_end,is_cut = self.get_antismash_starts(record,len(record.seq))
        if self.full_start == '' or self.full_start > antismash_start:
            self.full_start = antismash_start
        else:
            if self.full_start > antismash_start:
                self.full_start = antismash_start
        if self.full_end == '':
            self.full_end = antismash_end
        else:
            if self.full_end < antismash_end:
                self.full_end = antismash_end
        if is_cut == False:
            antismash_start = 0
        direction = str(feature.location).split('(')[1].split(')')[0]
        domain_name = feature.qualifiers['domain_id'][0].split('_')
        new_domain_name = []
        for splits in domain_name:
            splits = ''.join(splits.split(' '))
            if 'nrpspks' in splits or 'ctg' in splits or splits.isdigit() == True or 'input.path' in splits:
                pass
            else:
                if splits == 'Condensation':
                    splits = 'C'
                new_domain_name.append(splits.split('.')[0])
        new_domain_name = '_'.join(new_domain_name)
        new_domain_name = new_domain_name.upper()
        starts,ends,locations = list(),list(),list()
        typer = feature.type
        qualifiers = feature.qualifiers
        matches = re.findall(r'\d+', str(feature.location))
        for iter,item in enumerate(matches):
            if iter%2 == 0: ## the starts, at the even values
                starts.append(item)
            else: ## the ends, at the odd values. 
                ends.append(item)
        for iter,item in enumerate(starts):
            locations.append((item,ends[iter]))
        self.domain_info.append([locations,direction,typer,qualifiers,new_domain_name,antismash_start])
        return

    def get_anti_protocluster_info(self,feature,record):
        antismash_start,antismash_end,is_cut = self.get_antismash_starts(record,len(record.seq))

        if self.full_start == '' or self.full_start > antismash_start:
            self.full_start = antismash_start
        else:
            if self.full_start > antismash_start:
                self.full_start = antismash_start
        if self.full_end == '':
            self.full_end = antismash_end
        else:
            if self.full_end < antismash_end:
                self.full_end == antismash_end
        if is_cut == False:
            antismash_start = 0
        
        type = feature.qualifiers['category'][0] ##'pks or nrps'
        starts,ends,locations = list(),list(),list()
        typer = feature.type ##'protocluster'
        qualifiers = feature.qualifiers
        matches = re.findall(r'\d+', str(feature.location))
        for iter,item in enumerate(matches):
            if iter%2 == 0: ## the starts, at the even values
                starts.append(item)
            else: ## the ends, at the odd values. 
                ends.append(item)
        for iter,item in enumerate(starts):
            locations.append((item,ends[iter]))
        direction = str(feature.location).split('(')[1].split(')')[0]
        self.proto_info.append([locations,direction,typer,qualifiers,type,antismash_start])  

        return

    def get_antismash_starts(self,record,lenseq):
        '''antiSMASH already has cut a lot of the sections of contigs out for easier viewing. This is to pull that start,end information 
        so that we can adjust the start,end information that is pulled out when reading the genbank file.'''
        is_cut = True
        try:
            antismash_start = int(record.annotations['structured_comment']['antiSMASH-Data']['Orig. start']) # like three dictionaries within each other
        except KeyError:
            is_cut = False
        try:
            antismash_end = int(record.annotations['structured_comment']['antiSMASH-Data']['Orig. end'])
        except KeyError:
            is_cut = False

        if is_cut == False:
            locations = []
            for feature in record.features:
                if 'aSModule'in feature.type:
                    if 'nrps' in ','.join(feature.qualifiers['domains']):
                        matches = re.findall(r'\d+', str(feature.location))
                        for x in matches:
                            locations.append(int(x))
                if 'aSDomain' in feature.type:
                    if 'nrps' in feature.qualifiers['aSTool'][0]:
                        matches = re.findall(r'\d+', str(feature.location))
                        for x in matches:
                            locations.append(int(x))
                if 'protocluster' in feature.type:
                    if 'NRPS' in feature.qualifiers['category'] or 'PKS' in feature.qualifiers['category']:
                        matches = re.findall(r'\d+', str(feature.location))
                        for x in matches:
                            locations.append(int(x))
            mini = min(locations)
            maxi = max(locations)
            return mini,maxi,is_cut
        else:
            return antismash_start,antismash_end,is_cut

    def are_they_close_enough(self):
        '''I am applying the same standards to antiSMASH as I am to BiGEST. They must be within 10kbbp and in groups of 3 or more. (adjustable in parser)'''
        all_groupings = []
        to_sort = []
        for y in self.domain_info:
            s,e = self.get_start_end(y[0])
            to_sort.append((s,e,y))
        sorts = sorted(to_sort, key= itemgetter(0))## info is [(start,end,PKS,line),(start,end,PKS,line)] for each contig
        domain_info_to_use = []
        for s,e,y in sorts:
            domain_info_to_use.append(y)
        sorts = domain_info_to_use
        current_group = [sorts[0]]
        current_end = self.get_start_end(current_group[0][0])[1]
        if len(sorts) < num_matches_required:
            return None
        i = 0
        for obj in sorts[1:]:
            i += 1
            start = self.get_start_end(obj[0])[0]
            if abs(current_end - start) <= distance_required: ## they are within the boundary
                current_group.append(obj)
                current_end = self.get_start_end(obj[0])[1]
            else: ## they are too far apart
                if len(current_group)>= num_matches_required:
                    all_groupings.append(current_group)
                    current_group = [obj]
                    current_end = self.get_start_end(obj[0])[1]
                else:
                    current_group=[obj]
                    current_end = self.get_start_end(obj[0])[1]
        if len(current_group)>= num_matches_required:
            all_groupings.append(current_group)
        new = []
        for i in all_groupings:
            new.extend(i)
        unsorted_original_info = self.domain_info
        if all_groupings !=[]:
            to_remove = []
            for i in unsorted_original_info:
                if i not in new:
                    to_remove.append(i)
            for i in to_remove:
                unsorted_original_info.remove(i)
            self.domain_info = unsorted_original_info
            return new
        else:
            return None

    def get_start_end(self,locations):
        '''and another start and end grabber'''
        #locations,direction,typer,new_domain_name,qualifiers,antismash_start = self.domain_info
        nums = []
        for s,e in locations:
            nums.append(int(s))
            nums.append(int(e))
        start = min(nums)
        end = max(nums)
        return start,end 

class tbl:
    def __init__(self,finder_type,contig_name,antismash_start,this_antismash_start,locations,strand,result,featuretype,min_if_combined):

        ##the data##
        self.finder_type = finder_type
        self.contig = contig_name
        self.antismash_start = antismash_start ## this will be None for BiGEST
        self.strand = strand
        self.this_antismash_start = this_antismash_start ## none for BiGEST. Some contigs can have multiple entries in AS. 
        self.result = result ## this will be the type, domain name (for antismash) or the (bsg,result_name) --> like (KS,LCL_NRPS)
        self.featuretype = featuretype ## useful for antismash, module or domain, and 'still BiGEST' lol
        self.is_combined = min_if_combined
        self.locations = self.locations_to_int(locations)

    def adjust_antismash_starts(self, locations):
        if self.is_combined != None:## both combined and antiSMASH match
            new_locations = []
            for start,end in locations:
                start = int(start) + int(self.this_antismash_start) -int(self.is_combined)
                end = int(end) + int(self.this_antismash_start) -int(self.is_combined)
                new_locations.append((start,end))
            return new_locations
        else:
            return locations
    
    def locations_to_int(self,locations):
        new_locations = list()
        for location in locations:
            new_locations.append(tuple(map(int, location)))
        return new_locations

def filling_tbl_info_BiGEST(object,tbl_info,total_min):
    for start,end,bsg,line,strandz in object.updated_info:
        start = start-total_min-1
        end = end - total_min-1
        result = line.split()[11].split(",")[0]
        tbl_info.append(tbl('BiGEST',object.contig,None,None,[(start,end)],strandz,(bsg,result),'still BiGEST',None))
    return tbl_info

def get_sequence(fasta_file, contig_name,start=None,end=None):
    '''To pull sequence from fasta file required to write BiGEST output.'''
    if fasta_file.endswith(".gz"):
        with gzip.open(fasta_file,'rt') as fasta:
            for record in SeqIO.parse(fasta, "fasta"):
                if record.id == contig_name:
                    if start==None:
                        sequence = str(record.seq)
                    else:
                        sequence = str(record.seq)
                        sequence = sequence[start:end]
                    return sequence
    else:
        with open(fasta_file,'r') as fasta:
            for record in SeqIO.parse(fasta, "fasta"):
                if record.id == contig_name:
                    if start==None:
                        sequence = str(record.seq)
                    else:
                        sequence = str(record.seq)
                        sequence = sequence[start:end]
                    return sequence

def determine_number_of_BGC_groups(object,index_at_cut,num_matches_required):
    '''There could be multiple BGC groups on one contig. If there is a break in the distance (over 10kbp), that doesn't mean
    that is the only BGC groups found on that contig. This allows you to search across the entire contig to find them all.'''
    object.info = sorted(object.info,key = itemgetter(0))
    lists = list()
    distance_list = object.are_they_close_enough(object.info) ## will return None if no matches greater than 8 within 10000 bp
    while True:
        if len(object.info[index_at_cut:])< num_matches_required: ## there aren't enough BGC genes found within one section.
            break
        if distance_list == None: ### there aren't any aren't and BGCs found (either at all or within 10,000bp)
            break
        distance_list = object.are_they_close_enough(object.info[index_at_cut:]) ## do it again, starting from 0. 
        if distance_list:
            distance_list, indexed = distance_list ## split the output
            index_at_cut = index_at_cut + indexed ## increase the index at cut
            lists.append(distance_list) ### save the 3+ match 
    return lists

def write_BiGEST_only_gbk(object,fasta_file,output_directory,collapsed,distance_required):
    '''Writing out results for BiGEST only. This is the only file created if antiSMASH is not run/given. 
    
    Additionally, this fxn will create 
        a text file of just the BGC matches from the blast results
        a BED file for BiGEST results
        a GFF3 for BiGEST results 
    '''
    object.get_full_ends(fasta_file)
    output_gbk = os.path.join(output_directory, f"{object.file_name}_BiGEST_bgc_hits.gbk")
    lines = ''
    with open(output_gbk,"a") as output:
        lines = gbk_header(object.contig,object.file_name,object.full_start,object.full_end,'BiGEST')
        output.write(lines)
        if collapsed == True:
            joined_data = object.joining_collapsed()
        else:
            joined_data = object.joining() ## returns a 3 nested list: Outside list is the entire dataset. Second list is the groupings of "join", last is the "line" info originally 
        object.new_full_start = object.full_start
        fillings_data = write_fillings(object,joined_data)
        output.write(fillings_data)
        footer_data = gbk_footer(get_sequence(fasta_file,object.contig,object.full_start,object.full_end))
        output.write(footer_data)
    write_out_new_rpstblastn_txt_file_just_bgc_matches(object.updated_info,object.file_name)
    write_BED_BIG(object.contig,output_directory,object,fasta_name)
    write_gff3(object.contig,output_directory,object,fasta_name,joined_data,distance_required)
    return output_gbk

def get_protocluster(object):
    protocluster_filling = ''
    PKS = ['AT'',PT','KS','DH','ACP','KR','ER','MT','TE']
    NRPS = ['C','A','E','TE']
    for g in object.bgc_groups:
        mini,maxi = 999999999999999999,0 
        this_clust = set()
        for start,end,bsr,line,strand in g:
            mini = min((int(start)-object.new_full_start),int(end)-object.new_full_start,mini)
            maxi = max((int(start)-object.new_full_start),int(end)-object.new_full_start,maxi)
            this_clust.add(bsr)

        compare = {'nrps':len(list(set(NRPS).intersection(this_clust))),'pks':len(list(set(PKS).intersection(this_clust)))}
        max_value = max(compare.values())
        res = [key for key in compare if compare[key] == max_value]
        if len(res) == 2:
            string = 'PKS or NRPS'
        else:
            string = res[0]

        protocluster_filling += f"  BiGEST_cluster        {mini}..{maxi}\n"
        protocluster_filling += f"                   /label={string} hypothetical cluster\n"
        protocluster_filling += f"                   /note=\"This is the hypothetical cluster type concluded from domains included. All BiGEST domain types are {this_clust}\"\n"
    return protocluster_filling

def write_fillings(object,joined_data,combined=False):
    '''The information from the blast results need to be formatted for the genbank. This is the bulk of the information that goes in the genbank.'''
    lines = ''
    for grouping in joined_data:
        if len(grouping)==1: ## they don't need to be joined
            for start,end,bsr,line,strand in grouping:
                h = line.split('\t')
                q_start, q_end, sseqid, stitle, pident, evalue, gaps, qseq, sseq = int(h[4]), int(h[5]), h[1], h[10], h[2], h[8], h[11], h[12], h[13]
                if strand == '-':
                    new_start,new_end = q_end-object.new_full_start,q_start-object.new_full_start
                    location = f"complement({new_start}..{new_end})"
                else:
                    comp = False
                    new_start,new_end = q_start-object.new_full_start,q_end-object.new_full_start
                    location = f"{new_start}..{new_end}"
                label = stitle.split(',')[1].strip()
                lines += f"  BiGEST_domain         {location}\n"
                lines += f"                   /label={bsr}\t{label}\n"
                lines += f"                   /note=\"Evalue:{evalue}\tIdentity:{pident}%\tGaps:{gaps}\tQuery_seq:{qseq}\tSubject_Seq:{sseq}\t{stitle}\"\n"
        else: ## they do need to be joined. 

            sorted_joining_list = sorted(grouping,key=itemgetter(0))
            new_string = ''
            other_lines = ''
            for start,end,bsr,line,strand in sorted_joining_list:
                h = line.split('\t')
                q_start, q_end, sseqid, stitle, pident, evalue, gaps, qseq, sseq = int(h[4]), int(h[5]), h[1], h[10], h[2], h[8], h[11], h[12], h[13]
                label = stitle.split(',')[1].strip()
                other_lines+= f"                   /label={bsr}\t{label}\n"
                other_lines += f"                   /note=\"Evalue:{evalue}\tIdentity:{pident}%\tGaps:{gaps}\tQuery_seq:{qseq}\tSubject_Seq:{sseq}\t{stitle}\"\n"
                if strand == '+':
                    new_start, new_end = q_start-object.new_full_start,q_end-object.new_full_start
                    new_string+=f'{new_start}..{new_end},'
                    comp = False
                else:
                    new_start,new_end = q_end-object.new_full_start,q_start-object.new_full_start
                    new_string+=f'{new_start}..{new_end},'
                    comp= True
            new_string = new_string[0:-1] ## take off the last comma
            if comp == True:
                loc = f'complement(join({new_string}))'
            else:
                loc = f'join({new_string})'

            lines += f"  BiGEST_domain         {loc}\n"
            lines += other_lines
    lines += get_protocluster(object)
    return lines

def gbk_header(fasta_name,file_name,start,end,type):
    '''The GBK header will have information on the locus, where the original start and end are for the contig (calculated in get_start_end for both classes), 
    and start the formatting for the features in the genbank. '''
    lines = ''
    lines += f"LOCUS  {fasta_name}: {file_name} \t {end-start} bp\n"
    lines += f"AUTHORS     Lisa Adriani MS., Taehyung Kwon, PhD., Blake Hovde, PhD.\n"
    lines += f'COMMENT     ##COMMENT START##\n'
    lines += f'            BiGEST (Biosynthetic Gene Eukaryotic Search Tool) Results\n'
    if type == 'combined':
        lines += f'            This combines the results from BiGEST and antiSMASH into one genbank file to be viewed together\n'
    elif type == 'antismash':
        lines += f'            This is an antiSMASH record that did not match any BiGEST results.\n'
    elif type == 'BiGEST':
        pass
    elif type == 'BiGEST+antismash':
        lines += f'            This combines the results from BiGEST and antiSMASH into one genbank file to be viewed together\n'
    lines += f'            NOTE: This is a single region extracted from a larger record!\n'
    lines += f'            Orig. start  :: {start+1}\n'
    lines += f'            Orig. end :: {end+1}\n'
    lines += f'            ##COMMENT END##\n'
    lines += "FEATURES             Location/Qualifiers\n"
    return lines 

def gbk_footer(seq):
    '''The footer will format the sequence correctly for a genbank file.'''
    lines = ''
    lines += "ORIGIN\n"
    for i in range(0,len(seq), 60):
        lines += f"{i+1:9} {seq[i:i+60]}\n"
    lines += "//\n"
    return lines

def write_out_new_rpstblastn_txt_file_just_bgc_matches(info,fasta_name):
    with open(os.path.join(output_directory,f'{fasta_name}.BiGEST_bgc_hits.rpstblastn.txt'),'a+') as output:
        for obj in info:
            output.write(obj[3])
    return

def write_combined_output(successful_contigs,contigs_dict,antismash_dict, output_directory,fasta_file,collapsed,fasta_name):
    '''This will loop through all of the BiGEST hits on each of the contigs of the fasta file. If antiSMASH is given, it will include these results within each contig, with different feature names.
    If there are no BiGEST hits on the contig and there are antiSMASH results, those will be included in the "combined" genbank. '''
    ### write the gbk ## 
    output_gbk_name = f"combined_{fasta_name}.gbk" 
    output_gbk = os.path.join(output_directory,output_gbk_name)
    tbl_info = []
    written_non_BiGEST = False
    with open(output_gbk,'w') as output:
        for contig in successful_contigs:
            try:
                total_min = min(contigs_dict[contig].full_start,antismash_dict[contig].full_start)
                total_max = max(contigs_dict[contig].full_end,antismash_dict[contig].full_end)
                type_of_gbk_entry = 'combined'
            except KeyError:
                total_min = contigs_dict[contig].full_start
                total_max = contigs_dict[contig].full_end
                type_of_gbk_entry = 'BiGEST'
            if type_of_gbk_entry !='BiGEST':
                written_non_BiGEST = True
            contigs_dict[contig].new_full_start = total_min
            output.write(gbk_header(contig,fasta_name,total_min,total_max,type_of_gbk_entry))

            if collapsed == True:
                filling_data = write_fillings(contigs_dict[contig],contigs_dict[contig].joining_collapsed(),combined=True)
            else:
                filling_data = write_fillings(contigs_dict[contig],contigs_dict[contig].joining(),combined=True)# returns a 3 nested list: Outside list is the entire dataset. Second list is the groupings of "join", last is the "line" info originally 
            tbl_info = filling_tbl_info_BiGEST(contigs_dict[contig],tbl_info,total_min)
            output.write(filling_data)
            ## if there is antismash, get those written in as well
            if type_of_gbk_entry == 'combined' or type_of_gbk_entry == 'BiGEST+antismash':
                if antismash_dict[contig].mod_info != []:
                    for locations,direction,typer,qualifiers,type,this_antismash_start in antismash_dict[contig].mod_info:
                        output.write(write_anti_filling(typer,qualifiers,locations,direction,antismash_dict[contig].full_start,this_antismash_start,None,total_min,antismash_dict[contig]))
                        write_BED_anti(contig,output_directory,locations,direction,qualifiers,fasta_name,this_antismash_start,None)
                        tbl_info.append(tbl('antiSMASH',contig,antismash_dict[contig].full_start,this_antismash_start,locations,direction,type.upper(),'module',total_min))
                    #write_gff3_anti(contig,output_directory,antismash_dict[contig],fasta_name,'mod')
                if antismash_dict[contig].domain_info != []:
                    for locations,direction,typer,qualifiers,new_domain_name,this_antismash_start in antismash_dict[contig].domain_info:
                        output.write(write_anti_filling(typer,qualifiers,locations,direction,antismash_dict[contig].full_start,this_antismash_start,new_domain_name,total_min,antismash_dict[contig]))
                        write_BED_anti(contig,output_directory,locations,direction,qualifiers,fasta_name,this_antismash_start,new_domain_name)
                        tbl_info.append(tbl('antiSMASH',contig,antismash_dict[contig].full_start,this_antismash_start,locations,direction,new_domain_name,'domain', total_min))
                    write_gff3_anti(contig,output_directory,antismash_dict[contig],fasta_name,'domain')
                if antismash_dict[contig].proto_info != []:
                    for locations,direction,typer,qualifiers,new_domain_name,this_antismash_start in antismash_dict[contig].proto_info:
                        output.write(write_anti_filling(typer,qualifiers,locations,direction,antismash_dict[contig].full_start,this_antismash_start,None,total_min,antismash_dict[contig]))
            output.write(gbk_footer(get_sequence(fasta_file,contig,total_min,total_max)))

        ### antismash only ## 
        for contig_name,information in antismash_dict.items():
            if contig_name not in successful_contigs:
                written_non_BiGEST = True
                typ = 'antismash'
                total_min,total_max = information.full_start,information.full_end
                fillings = ''
                if antismash_dict[contig_name].mod_info != []:
                    for locations,direction,typer,qualifiers,type,this_antismash_start in antismash_dict[contig_name].mod_info:
                        fillings += write_anti_filling(typer,qualifiers,locations,direction,antismash_dict[contig_name].full_start,this_antismash_start,None,total_min = total_min)
                        write_BED_anti(contig_name,output_directory,locations,direction,qualifiers,fasta_name,this_antismash_start,None)
                        tbl_info.append(tbl('antiSMASH',contig_name,antismash_dict[contig_name].full_start,this_antismash_start,locations,direction,type.upper(),'module',None))
                    #write_gff3_anti(contig,output_directory,antismash_dict[contig_name],fasta_name,'mod')
                if antismash_dict[contig_name].domain_info != [] and len(antismash_dict[contig_name].domain_info)>1:
                    for locations,direction,typer,qualifiers,new_domain_name,this_antismash_start in antismash_dict[contig_name].domain_info:
                        fillings += write_anti_filling(typer,qualifiers,locations,direction,antismash_dict[contig_name].full_start,this_antismash_start,new_domain_name,total_min = total_min)
                        write_BED_anti(contig_name,output_directory,locations,direction,qualifiers,fasta_name,this_antismash_start,new_domain_name)
                        tbl_info.append(tbl('antiSMASH',contig_name,antismash_dict[contig_name].full_start,this_antismash_start,locations,direction,new_domain_name,'domain',None))
                    write_gff3_anti(contig_name,output_directory,antismash_dict[contig_name],fasta_name,'domain')
                if antismash_dict[contig_name].proto_info != []:
                    for locations,direction,typer,qualifiers,new_domain_name,this_antismash_start in antismash_dict[contig_name].proto_info:
                        fillings += (write_anti_filling(typer,qualifiers,locations,direction,antismash_dict[contig_name].full_start,this_antismash_start,None,total_min,antismash_dict[contig_name]))
                if fillings != '':
                    output.write(gbk_header(contig_name,fasta_name,total_min,total_max,typ))
                    output.write(fillings)
                    output.write(gbk_footer(get_sequence(fasta_file,contig_name,total_min,total_max)))
    output_table_name = f"combined_table_fmt_{fasta_name}.txt"
    output_table = os.path.join(output_directory,output_table_name)
    tbl_info = check_for_double(tbl_info)
    with open(output_table,'w') as output:
        for one in tbl_info:
            if one.finder_type == 'antiSMASH':
                for start,end in one.locations:
                    if one.featuretype == 'module':
                        output.write(f'{one.contig}\tantiSMASH\t\t\t{one.result} Module\t{start}\t{end}\t{one.strand}\n')
                                    ## contig featuretype domain_name.upper() start , end strand
                    else: ## featureType is asDomain
                        if one.result.startswith('PKS'):
                            bsg = one.result.split('PKS_')[1]
                        else:
                            bsg = ''
                        output.write(f'{one.contig}\tantiSMASH\t{bsg}\t{one.result}\t\t{start}\t{end}\t{one.strand}\n')

            if one.finder_type == 'BiGEST':
                output.write(f'{one.contig}\tBiGEST\t{one.result[0]}\t{one.result[1]}\t\t{one.locations[0][0]}\t{one.locations[0][1]}\t{one.strand}\n')

    if os.path.getsize(output_table) == 0:
        os.remove(output_table)
    if written_non_BiGEST == False:
        os.remove(output_gbk)
        print(f'No non-BiGEST hits were found. combined file is deleted.')

    return output_gbk

def check_for_double(all_tbls):
    for one in all_tbls:
        if one.finder_type == 'antiSMASH':
            one.locations = one.adjust_antismash_starts(one.locations)
    sorted_list = sorted(all_tbls, key=lambda x: (x.contig, int(x.locations[0][0])))
    previous_lowest,previous_highest = None,None
    counter = 1
    current_contig_name = None
    current_contig_grouping = None
    for one in sorted_list:
        current_highest = max(max(one.locations, key=lambda x: max(x))) # get the typle with the max value, get the max value in the tuple. 
        current_lowest = min(min(one.locations, key=lambda x: min(x)))
        one.og_contig = one.contig
        if previous_lowest != None:
            if abs(previous_lowest-current_lowest) > 10000 and abs(current_highest - previous_highest)>10000:
                if current_contig_grouping == one.og_contig:
                    counter+=1
                    current_contig_name = f'{one.contig}_{counter}'
                    one.contig = current_contig_name
                else:
                    counter = 1
                    current_contig_name = f'{one.contig}_{counter}'
                    one.contig = current_contig_name
            else:
                if current_contig_name!=None:
                    if current_contig_grouping == one.contig:
                        one.contig = current_contig_name
                    else:
                        counter = 1
                        current_contig_name = f'{one.contig}_{counter}'
                        one.contig = current_contig_name
        else:
            current_contig_name = f'{one.contig}_{counter}'
            one.contig = current_contig_name
            
        current_contig_grouping = one.og_contig
        previous_lowest = current_lowest
        previous_highest= current_highest
    return all_tbls

def write_BED_anti(contig_name,output_directory,locations,direction,qualifiers,fasta_name,this_antismash_start,fixed_domain_name=None):
    '''Writing a BED file for the antiSMASH results for easy additional analysis'''
    bed_file = os.path.join(output_directory,f'{fasta_name}.BED')
    if fixed_domain_name == None:
        for key,value in qualifiers.items():
            if key == 'label':
                fixed_domain_name = value
    if fixed_domain_name == None:
        fixed_domain_name = "aSModule"
    with open(bed_file,'a+') as b:
        for s,e in locations:
            ## bed is 0 based for start, 1 based for end
            s = int(s) + int(this_antismash_start)
            e = int(e) + 1 + int(this_antismash_start)
            b.write(f'{contig_name}\t{s}\t{e}\tantiSMASH:{fixed_domain_name}\t\t{direction}\n')
    return

def write_BED_BIG(contig_name,output_directory,contigs_obj,fasta_name):
    '''Writing a BED file for BiGEST results for easy additional analysis.'''
    bed_file = os.path.join(output_directory,f'{fasta_name}.BED')
    with open(bed_file,'a+') as b:
        for start,end,bsr,line,strand in contigs_obj.updated_info:
            ## bed is 0 based for start, 1 based for end
            start = int(start) -1
            end = int(end)
            stitle = line.split('\t')[10].split(',')[1].strip()
            label = f'{bsr}_{stitle}'
            b.write(f'{contig_name}\t{start}\t{end}\tBiGEST:{label}\t\t{strand}\n')
    return 

def write_anti_filling(typer,qualifiers,locations,direction,antismash_start,this_antismash_start,fixed_domain_name = None,total_min=None,object=None):
    '''Formatting the antiSMASH information for the genbank file. This information is formatted differently from BiGEST and needs to be edited.'''
    lines = ''
    new_locations = list()
    for start,end in locations:
        if total_min!=None:
            new_locations.append((int(start)+int(this_antismash_start)-int(total_min),int(end)+int(this_antismash_start) - int(total_min)))
        else:
            new_locations.append((int(start),int(end)))

    sorted_locations = sorted(new_locations,key=itemgetter(0))
    new_string = ''
    if len(sorted_locations)>1:
        for start,end in sorted_locations:
            new_string+=f'{start+1}..{end+1},'
        new_string = new_string[0:-1] ## take off the last comma
        if direction == '-':
            loc = f'complement(join({new_string}))'
        else:
            loc = f'join({new_string})'
    else:
        for start,end in sorted_locations:
            if direction == '-':
                loc = f'complement({int(start) + 1}..{int(end) +1})'
            else:
                loc = f'{int(start)+1}..{int(end)+1}'
    if typer == 'protocluster':
        typer = 'aSprotocluster'
    if len(f"  {typer}      ") <= 20:
        diff = 20 - len(f"  {typer}      ")
        extra_spaces = ' '*diff
        lines += f'  {typer}      {extra_spaces}{loc}\n'
    else:
        lines += f'  {typer}      {loc}\n'
    
    for key,value in qualifiers.items():
        value = '\t'.join(value)
        if key == 'label' and fixed_domain_name != None:
            value = fixed_domain_name
        lines += f'                   /{key}="{value}"\n'
    return lines

def write_gff3(contig_name,output_dir,contigs_obj,fasta_name,joined_data,distance_required):
    '''Writing a GFF3 for BiGEST results. This is for additional analysis.'''
    gff3_file = os.path.join(output_dir,f'{fasta_name}.gff3')
    seqid = contig_name
    source = 'BiGEST'
    with open(gff3_file,'a+') as g:
        if os.stat(gff3_file).st_size == 0:
            g.write('##gff-version 3\n')
        CDSs = ''
        label = ''
        totalmini, totalmaxi = None,None
        data = sorted(contigs_obj.updated_info,key=itemgetter(0)) ## sort by BSG and then location, then by direction
        for start,end,bsr,line,strand in data:
            #print(bsr)
            start = int(start) +1
            end = int(end)
            ## first time through
            if totalmini == None:
                totalmini,totalmaxi = start,end
                mini,maxi = start,end
                lasts,laste,lastbsr  = start,end,bsr
            else:
                ## make sure that they're the same. If not, write in the CDS's and restart 
                if lastbsr != bsr:
                    type = 'CDS'
                    score = '.'
                    strand = strand
                    phase = '.'
                    attributes = label[:-1]
                    CDSs += f'{seqid}\t{source}\t{type}\t{mini}\t{maxi}\t{score}\t{strand}\t{phase}\t{attributes}\n'
                    label = ''
                    mini,maxi = None,None
                else:
                    if max(start,lasts) <= min(end,laste):
                        pass
                    else:
                        type = 'CDS'
                        score = '.'
                        strand = strand
                        phase = '.'
                        attributes = label[:-1]
                        CDSs += f'{seqid}\t{source}\t{type}\t{mini}\t{maxi}\t{score}\t{strand}\t{phase}\t{attributes}\n'
                        label = ''
                        mini,maxi = None,None
                if abs(laste-start)>= distance_required:
                    g.write(f'{seqid}\t{source}\tgene\t{totalmini}\t{totalmaxi}\t.\t.\t.\tID:{contig_name}\n')
                    g.write(CDSs)
                    CDSs = ''
                    totalmini,totalmaxi = None,None
                else:
                    if mini == None:
                        mini,maxi = start,end
                    else:
                        if start < mini:
                            mini = start
                        if end > maxi:
                            maxi = end
                        if start < totalmini:
                            totalmini = start
                        if end > totalmaxi:
                            totalmaxi = end
                lasts,laste = start,end
                lastbsr = bsr
                stitle = line.split("\t")[10].split(',')[1].strip()
                lastlabel = f'{bsr}_{stitle}_'
                label += lastlabel

        type = 'gene'
        score = '.'
        phase = '.'
        strand = '.'
        attributes = f'ID:{contig_name}'
        g.write(f'{seqid}\t{source}\t{type}\t{totalmini}\t{totalmaxi}\t{score}\t{strand}\t{phase}\t{attributes}\n')
        g.write(CDSs)
    return

def write_gff3_anti(contig,output_directory,obj,fasta_name,information_type):
    '''Writing a GFF3 file for antiSMASH results for additional downstream analysis.'''
    gff3_file = os.path.join(output_directory,f'{fasta_name}.gff3')
    seqid = contig
    source = 'antiSMASH'
    with open(gff3_file,'a+') as g:
        if os.stat(gff3_file).st_size == 0:
            g.write('##gff-version 3\n')
        totalmini,totalmaxi = None,None
        CDSs = ''
        if information_type == 'mod':
            info = []
            fixed_domain_name = None
        else:
            fixed_domain_name = True
            info = obj.domain_info
        laste,lasts = None,None
        for locations,direction,typer,qualifiers,type,this_antismash_start in info:
            label = ''
            mini = None
            maxi = None
            for s,e in locations:
                s = int(s) +1 + int(this_antismash_start) 
                e = int(e) + int(this_antismash_start) 
                if totalmini == None:
                    totalmini,totalmaxi = s,e
                    mini,maxi = s,e 
                    laste,lasts = e,s
                elif mini == None:
                    mini,maxi = s,e
                if abs(laste-s)>= distance_required:
                    g.write(f'{seqid}\t{source}\tgene\t{totalmini}\t{totalmaxi}\t.\t.\t.\tID:{contig}\n')
                    g.write(CDSs)
                    CDSs = ''
                    totalmini,totalmaxi = s,e
                    mini,maxi = s,e
                else:
                    if s<mini:
                        mini = s 
                    if e > maxi:
                        maxi = e 
                    if s<totalmini:
                        totalmini = s 
                    if e > totalmaxi:
                        totalmaxi = e 
                laste, lasts = e,s
            for key,value in qualifiers.items():
                value = '\t'.join(value)
                if key == 'label' and fixed_domain_name != None:
                    value = type
                label += f'{key}={value}_'
            attributes = label[:-1]
            type = 'CDS'
            score = '.'
            strand = direction
            phase = '.'
            attributes = label[:-1]
            CDSs += f'{seqid}\t{source}\tCDS:{information_type}\t{mini}\t{maxi}\t{score}\t{strand}\t{phase}\t{attributes}\n'
            g.write(CDSs)
            CDSs = ''
        type = 'gene'
        score = '.'
        phase = '.'
        strand = '.'
        attributes = f'ID:{seqid}'

        g.write(f'{seqid}\t{source}\t{type}\t{totalmini}\t{totalmaxi}\t{score}\t{strand}\t{phase}\t{attributes}\n')
        g.write(CDSs)

    return

def visualize(output_directory,fasta_name):
    output_table_name = f"combined_table_fmt_{fasta_name}.txt"
    r_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "visualize.R")
    cmd = f'Rscript {r_script}  {output_directory}/{output_table_name} {output_directory}'
    result = subprocess.run(cmd,stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell=True, universal_newlines=True)
    print(result.stdout)
    print(result.stderr)

    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='BLAST-bassed sequence processing with Biosynthetic Gene Cluster Identification.')
    parser.add_argument('-i', '--input', help='Input txt file with blast output. Remember to follow the directions. If input blast results is not provided, rpsblastn will be run when blast_db provided.')
    parser.add_argument('-d', '--db', help='Path to CDD or subsetted-CDD blast database if blast needs to be run')
    parser.add_argument('-o', '--output_directory',required = True, help = "Output directory for PDF/GBK files")
    parser.add_argument('-f,','--fasta', required = True, help = "Path to original fasta file used for blast")
    parser.add_argument('--distance', default = 10000, type = int, help = 'The distance apart you will allow the end of one match to the start of the next match. 10,000bp default')
    parser.add_argument('--num_matches',default = 3, type = int, help = 'The number of matches required in order for a BGC to be "found". Default 3.')
    parser.add_argument('-a','--antismash_genbank', help = "The GenBank output from antiSMASH if you'd like integrated results. If you still need to run antiSMASH, set this equal to 'True'. ",default = None)
    parser.add_argument('-g','--gff3',help='If you have a gff3 to be run with antiSMASH, it can be included here. default = None',default=None)
    parser.add_argument('-c','--collapsed',help="If there are matches that overlap with the same classification, they will be collapsed in the final output. Default True. ",default = True)
    parser.add_argument('--coverage_min',help = 'Minimum coverage of the CDD subject match to be counted as a "match". Default is 0.0 (0%).', default = 0.0,type = float)
    args = parser.parse_args()

    start_time = time.time()
    fasta_file,output_directory,fasta_name,num_matches_required,distance_required,antismash_gbk,output_combined_gbk,collapsed,blast_data = set_variables(args)

    print(f'Starting on {fasta_name}')

    ## pull data into dictionaries
    contigs_dict,antismash_dict = access_data(blast_data,args.coverage_min,antismash_gbk,fasta_name) ## find the data that has BGC clusters
    successful_contigs = []

    ### Here's where the magic happens ######
    if len(contigs_dict) == 0:
        print("No matches found with BiGEST for this genome. Exiting program.")
    else:
        for object in contigs_dict.values():
            ### if there are multiple BGC groups on one contig, the information needs to be handled separately instead of combining them.
            BGC_GROUPS = determine_number_of_BGC_groups(object,0,num_matches_required) ## this will return a list containing a list of each BGC group found (>3,within 10kbp)
            object.bgc_groups = BGC_GROUPS
            if len(BGC_GROUPS)==1: ## only one group of BGC was found. 
                object.updated_info, indexed =  object.are_they_close_enough(object.info)
                file_name = write_BiGEST_only_gbk(object,fasta_file,output_directory,collapsed,distance_required)
                successful_contigs.append(object.contig)
                print(f'Output file created for {object.contig}')
            elif len(BGC_GROUPS)>1:
                object.updated_info = [item for cluster in BGC_GROUPS for item in cluster] ## loop through all bgc groups and append to updated info.
                file_name = write_BiGEST_only_gbk(object,fasta_file,output_directory,collapsed,distance_required)
                successful_contigs.append(object.contig)
                print(f'More than one group of {num_matches_required} found in {object.contig}. \n {str([len(item) for item in BGC_GROUPS])}')
                print(f'Output file created for {object.contig}')
            else:
                print(f'{object.contig} has no matches fulfilling requirements')
        if successful_contigs != []: ## there were BiGEST matches 
            combined_gbk = write_combined_output(successful_contigs,contigs_dict,antismash_dict,output_directory,fasta_file,collapsed,fasta_name)
            visualize(output_directory,fasta_name)
    print(f'Done with {fasta_name}')
    end_time = time.time()
    elapsed_time = round((time.time()-start_time),2)
    print(f'ELAPSED: {elapsed_time}')
