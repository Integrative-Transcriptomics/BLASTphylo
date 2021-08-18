# /***
# * BLASTphylo
# * processing of the data
# *
# * @author jennifer mueller
# *
# ***/


# packages
import sys
import os
import pandas as pd
import numpy as np
import tempfile
import subprocess
import json
from io import StringIO


# packages for Blast + Taxonomy mapping
from Bio.Blast.Applications import NcbiblastpCommandline as Blastp
from Bio.Blast.Applications import NcbiblastxCommandline as Blastx
from Bio.Blast.Applications import NcbiblastnCommandline as Blastn

from ete3 import NCBITaxa  # NCBI taxonomy (localy stored, require ~300MB)
#ncbi.update_taxonomy_database() # update actual NCBI taxonomy version

from ete3 import Tree  # handle trees


# packages for Phylogeny calculation
from server.external_tools import NCBIblastdbcmdCommandline as Blastdbcmd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MafftCommandline as Mafft
from server.external_tools import FastTreeCommandline as FastTree


# global variables 
count_clades = 0
normalization_dic = {}



########################################################################################################################
#                                                                                                                  BLAST
########################################################################################################################
''' run and parse the blast result dependent on the server input 
    prot                    amino acid sequence or direct Blast result as csv file
    prot_file_type:         0= aa sequence, 1= csv file
    blast_type:             blastp or blastx dependent on the input (by now only blastp)   --> TODO: change database if blastn should be implemented
    eValue:                 significant threshold for the blast search
    min_align_ident         minimal identity between query and subject sequence in the alignment
    min_query_cover:        minimal coverage with the query
    min_subject_cover:      minimal coverage with the subject (hit)
    entrez_query:           restrict the blast search to 'roots' of the given taxonomy
    out_dir:                path to output directory
'''
def run_blast(prot, prot_file_type, blast_type, eValue, min_align_ident, min_query_cover, min_subject_cover, entrez_query, out_dir):
    header = ['qacc', 'sacc', 'qstart', 'qend', 'sstart', 'send', 'slen', 'nident', 'evalue', 'pident', 'staxids', 'qcovhsp', 'sseq']

    header_basic = ['qacc', 'sacc', 'pident', 'alen', 'mm', 'g', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bit', 'm']
    if prot_file_type == "1": # input was csv file
        try:
            preFilter = pd.read_csv(StringIO(prot), sep='\t')
        except:
            preFilter = pd.read_csv(StringIO(prot), sep=',')

        if len(preFilter.columns) == len(header):
            preFilter.columns = header
        else:
            preFilter.columns = header_basic

        ''''# filter for alignment identity, query coverage, subject coverage, evalue
        subjectAlignedLength = abs(preFilter['send']-preFilter['sstart'])
        subjectCoverage = (subjectAlignedLength/preFilter['slen']) > (int(min_subject_cover)/100)

        alignIdent = preFilter['pident'] > (float(min_align_ident))

        result = preFilter[(preFilter['evalue'] < float(eValue)) & (preFilter['qcovhsp'] > (float(min_query_cover))) & subjectCoverage & alignIdent]

        print(result.shape)'''

    elif (prot_file_type == '2') or (prot_file_type == '3'):
        preFilter = pd.read_csv(prot, sep='\t')
        preFilter.columns = header

    else: # input was a protein
        blast_output_columns = '6 qacc sacc qstart qend sstart send slen nident evalue pident staxids qcovhsp sseq'

        # switch beteen the blasttypes
        if blast_type == 'blastp':
            blastp_cline = Blastp(cmd=blast_type, remote=True, query=prot, db='nr', evalue=eValue, max_hsps=1, num_alignments=1000,
                              qcov_hsp_perc=min_query_cover, entrez_query='\'' + entrez_query + '\'',
                              outfmt=blast_output_columns, out='blast_result.csv')
        elif blast_type == 'blastx':
            blastp_cline = Blastx(cmd=blast_type, remote=True, query=prot, db='nr', evalue=eValue, max_hsps=1, num_alignments=1000,
                              qcov_hsp_perc=min_query_cover, entrez_query='\'' + entrez_query + '\'',
                              outfmt=blast_output_columns)
        elif blast_type == 'blastn':
            blastp_cline = Blastn(cmd=blast_type, remote=True, query=prot, db='nt', evalue=eValue, max_hsps=1, num_alignments=1000,
                              qcov_hsp_perc=min_query_cover, entrez_query='\'' + entrez_query + '\'',
                              outfmt=blast_output_columns)
        elif blast_type == 'diamond': # implemented diamond but not tested --> miss local stored database
            blastp_cline = 'used diamond with default setting for blastp'
            # set up temporary output file
            temp_diamond_file = tempfile.NamedTemporaryFile(prefix=output_dir)
            output_diamond = temp_diamond_file.name

            # set up command for diamond
            '''
                need to set up a database --> pathToDatabase
                protein needs to be a file
                -p = threads 
                can automatically add the full sequences for query and subject (full_qseq, full_sseq)
            '''
            command_arguments_diamond = ['./diamond', 'blastp', '-q', prot, '-d', 'pathToDatabase', '-o', output_diamond,
                                         '-f', blast_output_columns, '-e', str(eValue), '--id', str(min_align_ident),
                                         '--query-cover', str(min_query_cover), '--subject-cover', str(min_subject_cover),
                                         '--max-target-seqs', '1000', '-p', '3']
            diamond = subprocess.Popen(command_arguments_diamond, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            diamond.wait()

        print(blastp_cline)
        stdout, stderr = blastp_cline()
        #print(stdout)
        blast_data = StringIO(stdout)
        print(stderr)
        preFilter = pd.read_csv('blast_result.csv', sep='\t', names=header)

    # filter for alignment identity, query coverage, subject coverage, evalue
    # evalue and query coverage are necessary given that BLAST stop when the first subject sequences exceed the threshold
    subjectAlignedLength = abs(preFilter['send']-preFilter['sstart'])
    subjectCoverage = (subjectAlignedLength/preFilter['slen']) > (int(min_subject_cover)/100)

    alignIdent = preFilter['pident'] > (float(min_align_ident))

    result = preFilter[(preFilter['evalue'] < float(eValue)) & (preFilter['qcovhsp'] > (float(min_query_cover))) & subjectCoverage & alignIdent]
    print(result.shape)

    return result



########################################################################################################################
#                                                                                                Filtering and Hit count
########################################################################################################################
'''
    Scoring function for the best hit of each taxa (taxa-based phylogeny)
    row:        actual hit 
    score:      previous best score

'''
def best_hit_scoring_function(row, score):

    subjectCoverage = abs(row['send']-row['sstart'])/row['slen']
    row_score = 0.25*(float(row['pident'])/100) + 0.25*float(row['evalue']) + 0.2*subjectCoverage + 0.2*(float(row['qcovhsp'])/100)

    if score > row_score:
        return score
    else:
        return row_score


'''Filter blast result for hits similar to the species in the tree
    blast_record:      result of the blast search
    tree_taxIDs:       set of all taxonomic IDs in the given taxonomy
'''
def filter_blast_result(blast_record, tree_tax_IDs):
    filtered_result = {}
    sequence_dic = {}  # taxid encoded seq dic --> for taxa-based phylogeny
    uniqueAccs = {}  # accession ID encoded seq dic  --> for unique sequence-based phylogeny
    filtered_blast = blast_record[blast_record['staxids'].notna()]  # remove all rows with no taxid information
    ncbi = NCBITaxa()

    number_of_queries_np = filtered_blast['qacc'].unique()
    number_of_queries = number_of_queries_np.tolist()

    for index, row in filtered_blast.iterrows():  # taxids from the hit table with number of occurrence

        if isinstance(row['staxids'], float):  # column only contain one element
            taxIDlist = [int(row['staxids'])]
        elif isinstance(row['staxids'], int):  # column only contain one element
            taxIDlist = [row['staxids']]
        else:
            taxIDlist = row['staxids'].split(';')

        for start_taxa in taxIDlist:
            if int(start_taxa) in tree_tax_IDs:
                uniqueAccs[row['sacc']] = [row['sseq'], len(taxIDlist)]

                try: # taxon was already in filtered_result
                    filtered_result[int(start_taxa)][number_of_queries.index(row['qacc'])] += 1

                    new_score = best_hit_scoring_function(row, sequence_dic[int(start_taxa)][2])
                    if new_score != sequence_dic[int(start_taxa)][2]:
                        sequence_dic[int(start_taxa)] = [row['sacc'], row['sseq'], new_score]
                except KeyError:
                    filtered_result[int(start_taxa)] = [0]*len(number_of_queries)
                    filtered_result[int(start_taxa)][number_of_queries.index(row['qacc'])] = 1
                    sequence_dic[int(start_taxa)] = [row['sacc'], row['sseq'], best_hit_scoring_function(row, 0)]

            else:
                try: # dependent on the ncbi taxonomy database version the taxon can not be present in the database
                    taxID_lineage = ncbi.get_lineage(start_taxa)
                except:
                    continue

                uniqueAccs[row['sacc']] = [row['sseq'], 0]
                for line_taxa in taxID_lineage:  # all higher level taxa of these taxids
                    if line_taxa in tree_tax_IDs:
                        uniqueAccs[row['sacc']][1] += 1

                        try: # taxon was already in filtered_result
                            filtered_result[line_taxa][number_of_queries.index(row['qacc'])] += 1

                            new_score = best_hit_scoring_function(row, sequence_dic[int(line_taxa)][2])
                            if new_score != sequence_dic[int(start_taxa)][2]:
                                sequence_dic[int(line_taxa)] = [row['sacc'], row['sseq'], new_score]
                        except KeyError:
                            filtered_result[line_taxa] = [0] * len(number_of_queries)
                            filtered_result[line_taxa][number_of_queries.index(row['qacc'])] = 1
                            sequence_dic[int(line_taxa)] = [row['sacc'], row['sseq'], best_hit_scoring_function(row, 0)]

                if uniqueAccs[row['sacc']][1] == 0: # remove the sequence when it was not present
                    uniqueAccs.pop(row['sacc'])

    return filtered_result, sequence_dic, uniqueAccs, number_of_queries



########################################################################################################################
#                                                                                                              Phylogeny
########################################################################################################################
'''
    extract best hit for each taxID
    seqs:           dictionary with sequence array per taxID
    output_file:    path to fasta for aligned sequences
'''
def extract_best_hit(seqs, output_file):
    best_seqs = []
    [best_seqs.append(SeqRecord(Seq(seqs[taxID][1]), id=str(taxID), description='')) for taxID in seqs.keys()]
    SeqIO.write(best_seqs, output_file, "fasta")
    print('Taxa-based phylogeny number of seqs. ' + str(len(best_seqs)))

'''
    generate Fasta from dict
    seqs:           dictionary with sequence array per taxID
    output_file:    path to fasta for aligned sequences
'''
def generate_fasta(seqs, output_file):
    fasta = []
    [fasta.append(SeqRecord(Seq(seqs[key][0]), id=key, description='')) for key in seqs.keys()]
    SeqIO.write(fasta, output_file, "fasta")
    print('unique sequence-based phylogeny number of seqs. ' + str(len(fasta)))


'''
   transfer Mafft alignment in taxid version 
   mafft_out:           path to mafft output
   taxid_accession:     dictionary with subject aligned seqs and their accessions
'''
def generate_fasttree_input(mafft_out, taxid_accession):
    new_mafft = []
    records = SeqIO.parse(mafft_out, 'fasta')
    for record in records:
        try:
            [new_mafft.append(SeqRecord(record.seq, id=str(taxId), description='')) for taxId in taxid_accession[record.id[:-2]]]
        except KeyError:
            continue
    SeqIO.write(new_mafft, mafft_out.split('.')[0] + '_NEW.fasta', 'fasta')

    return None


''' calculate phylogeny for given hit sequences
    subject_seqs:         dictionary with subject aligned seqs and their accessions
    fasta:                placeholder for phylogeny calculation with full subject sequences 
    output_dir:           path to output directory
    from_aligned_seq:     boolean, true: aligned subject sequences, false: full subject sequences
    uniqueSeqs:           boolean, true: accession encoded seqs, false: taxid encoded seqs
    based on: #https://biopython.org/wiki/Phylo
'''
def calculate_phylogeny(subject_seqs, fasta, output_dir, from_aligned_seq, uniqueSeqs):

    # generation of tmp file
    # temp fasta for MAFFT
    temp_fasta_file = tempfile.NamedTemporaryFile(prefix=output_dir)
    output_seqs = temp_fasta_file.name

    # temp tree for FastTree
    temp_tree_file = tempfile.NamedTemporaryFile(prefix=output_dir)
    output_tree = temp_tree_file.name


    if from_aligned_seq == 'True':
        if uniqueSeqs:
            print('Generate FastA file')
            generate_fasta(subject_seqs, output_seqs)
        else:
            print('Extract best Hit')
            extract_best_hit(subject_seqs, output_seqs)
    else:
        output_seqs = fasta
        print('Input was fasta file')

        ### Code for full sequence alignments if a local protein database is stored, TODO: need to be changed (path to database) and tested
        # list_of_accessions = [subject_seqs[key][0][0] for key in subject_seqs.keys()]
        # blastdbcmd_cline = Blastdbcmd(db="/share/references/library/refseq_protein/refseq_protein", entry=','.join(list_of_accessions), out=output_seqs)  #db path = default path to refseq database on Romanov
        # blastdbcmd_stdout, blastdbcmd_stderr = blastdbcmd_cline()

    # calculate MSA
    '''mafft_file = "mafft_aligned_unique.fasta" if uniqueSeqs else "mafft_aligned.fasta"
    output_fasta = output_dir + mafft_file

    mafft_cline = Mafft(input=output_seqs, amino=True, clustalout=False,                                             # old version: MSA file generation
                        parttree=True, thread=4)  # treeout --> guidetree in output                                     # is needed if full sequence alignments should be calculated
    print(mafft_cline)
    stdout_Mafft, stderr_Mafft = mafft_cline()  # stdout_Mafft = MSA
    with open(output_fasta, "w") as handle:
        handle.write(stdout_Mafft)
    handle.close()

    if from_aligned_seq == 'False':  # need to convert the Mafft full seq output to a similar shape as for aligned seqs
        acc_taxids = {}
        for key in subject_seqs:
            try:
                acc_taxids[subject_seqs[key][0]].append(key)
            except KeyError:
                acc_taxids[subject_seqs[key][0]] = [key]
        generate_fasttree_input(output_fasta, acc_taxids)
        output_fasta = output_dir + mafft_file.split('.')[0] + '_NEW.fasta'

    # generate phylogenetic tree
    fasttree_file = "fasttree_unique.tree" if uniqueSeqs else "fasttree.tree"
    output_tree = output_dir + fasttree_file

    fasttree_cline = FastTree(input=output_fasta, quiet=True,                                                        # old version: input is MSA file path
                              out=output_tree, fastest=True)  # fastest --> faster calculation by focus on best hit



    print(fasttree_cline)
    stdout_FastTree, stderr_FastTree = fasttree_cline()'''

    # alternative calculation to directly pipe the MAFFT output in FastTree
    command_arguments_mafft = ['mafft', '--thread', '-1', output_seqs]
    command_arguments_fasttree = [ 'fasttree', '-fastest', '-quiet', '-out', output_tree]                               ########################################### ATTENTION: check your FastTree command
    mafft = subprocess.Popen(command_arguments_mafft, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #print(mafft.stdout)
    fasttree = subprocess.Popen(command_arguments_fasttree, stdin=mafft.stdout, stdout=subprocess.PIPE)
    #print(fasttree.stdout)
    mafft.wait()
    fasttree.wait()

    print('Start Tree calculation and processing data')
    tree, tree_tax_ids, _ = read_tree_input(output_tree, '2', True)
    global count_clades

    # close temp files
    temp_tree_file.close()
    temp_fasta_file.close()

    if uniqueSeqs:  # unique sequence-based phylogeny
        count_clades = 0
        return unique_phylogeny_data(tree, subject_seqs)
    else: # taxa-based phylogeny
        count_clades = 0
        d3Tree, parentnode_info, mosthitAcc = phylogeny_data(tree, subject_seqs, tree_tax_ids)
        newickTree = tree.write(format=3)
        temp_fasta_file.close()
        return d3Tree, newickTree, parentnode_info, mosthitAcc

'''
    Taxon name has to be in a newick string valid form
    taxa:   string with scientific name
'''
def newick_valid_name(taxa):
    nonValidChar = [' ', '[', ']', '(', ')', ':', '/', '|', '\'']
    newName = taxa
    for nonValid in nonValidChar:
        newName = newName.replace(nonValid, '_')
    return newName

''' replace newick string with taxa IDs to newick string with scientific names
    newick:     taxa ID newick string (output of FastTree)
    treeIDs:    list of all tree taxa IDs
'''
def replace_ID_Name(newick, treeIDs):
    ncbi = NCBITaxa()
    taxaNames = ncbi.get_taxid_translator(list(treeIDs))
    newTree = newick
    for taxa in list(treeIDs):
        taxaName = newick_valid_name(taxaNames[taxa])
        newTree = newTree.replace(str(taxa)+':', taxaName+':')
    return newTree


''' 1. Phylogeny: taxa based 
    dependent of the visualisation different additional information are needed for the further calculation
    d3 Visualisation: hold all information
    PhyloTree Visualisation: only contain the tree information

    tree:           FastTree newickString as tree like datastructure
    subject_seqs:   taxa encoded dic with hit Accession and sequences
    treeIDs:        list of all taxa IDs in the tree
'''
def phylogeny_data(tree, subject_seqs, treeIDs):
    ncbi = NCBITaxa()

    pseudo_filtered_blast = {}  # instead of the hit counts the parent taxid is saved
    parentnode_info = {}


    for key in list(treeIDs):
        defaultvalue = 1
        try: # extract phylum information
            key_lineage = ncbi.get_rank(ncbi.get_lineage(key))
            parent_of_key = [key for key in key_lineage if
                             key_lineage[key] == 'phylum']
            pseudo_filtered_blast[key] = [parent_of_key[0]]
            parentID, parentName = translate_node(str(parent_of_key[0]), 0)
            parentnode_info[parentID] = parentName  # parentnode_info contain parent sci name + taxid for visualisation
        except:
            pseudo_filtered_blast[key] = [defaultvalue]

    # extract all taxonomic IDs with unique accessions for the best hit
    accvalues = [subject_seqs[key][0] for key in subject_seqs.keys()]
    uniqueSeqs = [key for key in subject_seqs.keys() if accvalues.count(subject_seqs[key][0]) == 1]

    # generate cladogram dictionary
    newTree = wrapper_transfer_own_tree(tree, pseudo_filtered_blast, '2', 10, 1)
    return newTree, parentnode_info, uniqueSeqs



''' 2. Phylogeny: unique hit based
    tree:           see 1. Phylogeny
    uniqueAccs:     hit accession encoded dic with seqs + #taxa connected to hit
'''
def unique_phylogeny_data(tree, uniqueAccs):
    unique_count = {key: uniqueAccs[key][1] for key in uniqueAccs.keys()}
    d3_tree = unique_transfer_tree(tree, unique_count, 20)
    newickTree = tree.write(format=3)

    return d3_tree, newickTree



########################################################################################################################
#                                                                                                 Taxonomy/Tree handling
########################################################################################################################
'''read tree input 
    tree_input:     tree as newick string, file or string of sientific names/IDs to generate a tree from the NCBI taxonomy
    ncbi_boolean:   boolean to seperate the two input types (0 = NCBI, 1 = own taxonomy) 
                    NCBI example: Staphylococcus,Staphylococcus aureus|subtree
    needTaxIDs:     boolean
'''
def read_tree_input(tree_input, ncbi_boolean, needTaxIDs):
    tree = None
    roots = ''  # ncbi taxonomy can begin with more than one 'root'
    tree_taxIDs = set()

    if ncbi_boolean == '0': # NCBI taxonomy
        roots, tree_taxIDs = translate_string_in_taxa(tree_input)

    elif ncbi_boolean == '1': # own taxonomic phylogeny
        tree = Tree(tree_input, format=8)
        tree_taxIDs = translate_nodes([node.name for node in tree.traverse('preorder')])
        roots = 'txid' + translate_node(tree.name.replace('_', ' '), 0)[1] + '[ORGN]'

    elif ncbi_boolean == '2' and needTaxIDs: # phylogeny calculation
        tree = Tree(tree_input)
        tree_taxIDs = translate_nodes([node.name for node in tree.traverse('preorder')])
        roots = 'txid' + translate_node(tree.name.replace('_', ' '), 0)[1] + '[ORGN]'

    elif ncbi_boolean == '2' and not needTaxIDs: #phylogeny data generation
        tree = Tree(tree_input)
        tree_taxIDs = set()

    else:
        sys.stderr.write('no tree given')

    return tree, tree_taxIDs, roots

def translate_string_in_taxa(taxonomy):
    tree_taxa = set()
    regular_exp = ''
    root = ''
    not_option = ''
    previous_symbol = ''
    ncbi = NCBITaxa()

    if taxonomy[len(taxonomy)-1] == ',':  # taxonomy was generated with Search option --> last character is ','
        taxonomy = taxonomy[:-1]
    print(taxonomy)

    i = 0
    while i < len(taxonomy):
        #print(root)
        if taxonomy[i] == '|':
            if taxonomy[i+1] == '!':   # input was taxa|!(....)
                previous_symbol = '!'
                #print(root)
                j = i+3
                while taxonomy[j] != ')':     # iterate over inner nodes
                    not_option += taxonomy[j]
                    j = j + 1
                #print(not_option)
                not_regular, not_taxa = translate_string_in_taxa(not_option)  # evaluate subexpression

                rootID, rootName = translate_node(root, 0)
                root_descendents =  ncbi.get_descendant_taxa(rootID, intermediate_nodes=True)
                root_descendents.append(rootID)

                if 'OR' in not_regular:  # for complex subexpressions brackets around subexpression | root without subexpression
                    regular_exp = regular_exp + '(txid' + str(rootID) + '[ORGN] NOT (' + not_regular + '))'
                else:
                    regular_exp = regular_exp + '(txid' + str(rootID) + '[ORGN] NOT ' + not_regular + ')'
                tree_taxa.update(set(root_descendents).difference(not_taxa))                             # root taxa without subexpression taxa

                root = ''
                not_option = ''
                previous_symbol = ''
                i = j+1

            elif taxonomy[i+1:i+8] == 'subtree':   # input was taxa|subtree
                #print('extract subtree')
                rootID, rootName = translate_node(root, 0)
                tree_taxa.update(set(ncbi.get_descendant_taxa(rootID, intermediate_nodes=True)))
                tree_taxa.update(set([rootID]))
                regular_exp = regular_exp + 'txid' + str(rootID) + '[ORGN]'              # important have to be ID
                root = ''
                not_option = ''
                previous_symbol = ''
                i = i+8

        elif taxonomy[i] == ',':  # check next list element
            if len(root) > 0:  # input was  taxa,
                rootID, rootName = translate_node(root, 0)
                regular_exp = regular_exp + 'txid' + str(rootID) + '[ORGN] OR '
                tree_taxa.update(set([rootID]))
                root = ''
            else:
                regular_exp = regular_exp + ' OR '
            i += 1
        else:
            if previous_symbol == '!':
                not_option += taxonomy[i]
            else:
                root += taxonomy[i]
            i += 1

    if len(root) > 0:   # add last element of the list
        rootID, rootName = translate_node(root, 0)
        regular_exp = regular_exp + 'txid' + str(rootID) + '[ORGN]'
        tree_taxa.update(set([rootID]))

    return regular_exp, tree_taxa

'''Translate all nodes to NCBI taxonomy ID
    tree_nodes:     list of all tree nodes  
'''
def translate_nodes(tree_nodes):
    ncbi = NCBITaxa()  # setup NCBI database
    translated_set = set()
    count_not_valid_ids = 0
    for node in tree_nodes:
        if node != '':
            if node.isnumeric():
                translated_set.add(int(node))
            else: # need to translate the scientific name
                try:
                    translated_node = ncbi.get_name_translator([node.replace('_', ' ')])
                    translated_set.add(translated_node[node.replace('_', ' ')][0])
                except KeyError:
                    count_not_valid_ids += 1
    print('During translation of the tree nodes ' + str(count_not_valid_ids) + ' could not be translated')
    return translated_set


########################################################################################################################
#                                                                                                   d3-compatible format
########################################################################################################################
'''Transfer the given tree in a d3 compatible hierarchical format
    tree:                   input tree 
    filtered_blast_result:  filtered blast results with found hits 
    ncbi_boolean:           dependent on the input tree the nodes could have a mixed labeling 
    number_of_queries:      number of input queries 
'''
def transfer_tree_in_d3(tree, filtered_blast_result, ncbi_boolean, number_of_queries):
    d3_tree = {}
    if ncbi_boolean == '0':
        d3_tree = wrapper_transfer_ncbi_tree(filtered_blast_result, ncbi_boolean, 30, number_of_queries)
    elif ncbi_boolean == '1':
        d3_tree = wrapper_transfer_own_tree(tree, filtered_blast_result, ncbi_boolean, 30, number_of_queries)
    else:
        sys.stderr.write('check the input again')
    return d3_tree


''' generate tree and transfer it into d3 format
    filtered_blast_result:  result of the blast filtering is used to generate the minimal tree to connact the nodes
    ncbi_boolean:           see transfer_tree_in_d3     
    branch_length:          default branch length (unique sequence-based phylogeny need more space)
    number_of_queries:      number of input queries  
'''
def wrapper_transfer_ncbi_tree(filtered_blast_result, ncbi_boolean, branch_length, number_of_queries):

    ncbi = NCBITaxa()  # setup NCBI database

    # get rank information for all taxa in filtered blast result
    taxa_keys = filtered_blast_result.keys()
    tree = ncbi.get_topology(taxa_keys)
    treeTaxIDs = translate_nodes([node.name for node in tree.traverse('preorder')])
    treeRanks = ncbi.get_rank(list(treeTaxIDs))

    if len(taxa_keys) == 1: # tree only contains root
        root_value = [[filtered_blast_result[int(tree.name)][i], 0] for i in range(number_of_queries)]
        root_value.append(treeRanks[int(tree.name)])
        return {'size': [15, branch_length], 'name': tree.sci_name, 'value': root_value, 'leaf_counter': 1}

    else: # iterate recursive over children
        root_children = list(
            filter(None, [transfer_tree(child, filtered_blast_result, treeRanks, ncbi_boolean, branch_length, number_of_queries) for child in tree.get_children()]))

        root_value = [[] for i in range(number_of_queries)]  # contain the hits for the root and the sum of all children hits
        try:
            [root_value[i].append(filtered_blast_result[int(tree.name)][i]) for i in range(number_of_queries)]
        except KeyError:
            [root_value[i].append(0) for i in range(number_of_queries)]

        subtree_hits = [0]*number_of_queries
        root_size = 0
        leaf_counter = 0
        for child in root_children:
            for i in range(number_of_queries):
                subtree_hits[i] += sum(child['value'][i][:2])
            root_size += child['size'][0]
            leaf_counter += child['leaf_counter']


        if ncbi_boolean == '2': # taxa-based phylogeny
            root_value.append(int(tree.name))
        else:
            [root_value[i].append(subtree_hits[i]) for i in range(number_of_queries)]
            root_value.append(treeRanks[int(tree.name)])

        return {'size': [root_size, branch_length], 'name': tree.sci_name, 'children': root_children,
                'value': root_value, 'leaf_counter': leaf_counter}


''' transfer own tree with possible mixed labeling in d3 format
    tree:                   given tree structure
    filtered_blast_result:  dictionary of hits  
    ncbi_boolean:           see transfer_tree_in_d3
    branch_length:          default branch length
    number_of_queries:      number of input queries
'''
def wrapper_transfer_own_tree(tree, filtered_blast_result, ncbi_boolean, branch_length, number_of_queries):
    global count_clades
    ncbi = NCBITaxa()

    # get rank for all taxa in tree
    taxa_keys = filtered_blast_result.keys()
    treeTaxIDs = translate_nodes([node.name for node in tree.traverse('preorder')])
    treeRanks = ncbi.get_rank(list(treeTaxIDs))
    rootID, rootName = translate_node(tree.name.replace('_', ' '), count_clades)

    if len(taxa_keys) == 1: # filtered blast only contains root
        root_value = [[filtered_blast_result[rootID][i], 0] for i in range(number_of_queries)]
        root_value.append(treeRanks[rootID])
        return {'size': [15, branch_length], 'name': tree.sci_name, 'value': root_value, 'leaf_counter':1}

    else: # iterate recursive over children
        root_children = list(
            filter(None, [transfer_tree(child, filtered_blast_result, treeRanks, ncbi_boolean, branch_length, number_of_queries) for child in tree.get_children()]))

        root_value = [[] for i in range(number_of_queries)]  # contain the hits for the root and the sum of all children hits
        try:
            [root_value[i].append(filtered_blast_result[rootID][i]) for i in range(number_of_queries)]
        except KeyError:
            [root_value[i].append(0) for i in range(number_of_queries)]


        subtree_hits = [0]*number_of_queries
        root_size = 0
        leaf_counter = 0
        list_of_child_ids = []  # phylogeny: all child ids to calculate the last common ancestor
        count_clade_members = 0  # phylogeny: addition to cladeName to distinguish between clades
        for child in root_children:
            if ncbi_boolean != '2':
                for i in range(number_of_queries):
                    subtree_hits[i] += sum(child['value'][i][:2])
            root_size += child['size'][0]
            leaf_counter += child['leaf_counter']
            if ncbi_boolean == '2':
                count_clade_members += child['value'][2]
                if 'clade' not in child['name']:
                    list_of_child_ids.append(child['value'][1])

        # get the last common ancestor (LCA) for the inner nodes in the taxa-based phylogeny
        # same LCA for different branches are distinguished by numbers
        if ncbi_boolean == '2':
            root_value = root_value[0]
            if len(list_of_child_ids) == 1:
                try:
                    rootName = ncbi.get_taxid_translator(list_of_child_ids)[list_of_child_ids[0]] + '_' + str(count_clade_members) + '.' + str(count_clades)
                    count_clades += 1
                    rootID = list_of_child_ids[0]
                except:
                    print('No name for:' + str(list_of_child_ids[0]))

            elif len(list_of_child_ids) > 0:
                try:
                    children_tree = ncbi.get_topology(list_of_child_ids)
                    rootName = children_tree.sci_name + '_' + str(count_clade_members) + '.' + str(count_clades)
                    count_clades += 1
                    rootID = int(children_tree.name)
                except:
                    print('No topology for :' + ','.join([str(item) for item in list_of_child_ids]))
            tree.name = newick_valid_name(rootName)
            root_value.append(rootID)
            root_value.append(count_clade_members)
        else:
            [root_value[i].append(subtree_hits[i]) for i in range(number_of_queries)]
            root_value.append(treeRanks[rootID])

        return {'size': [root_size, branch_length], 'name': rootName, 'children': root_children,
                'value': root_value, 'leaf_counter': leaf_counter}


''' recursively transfer the subtrees in d3 compatible format
    tree:                   subtree
    filtered_blast_result:  dictionary of hits
    treeRanks:              dictionary of all ranks for the taxa in the tree
    ncbi_boolean:           dependent on the tree format the node information is different
    branch_length:          default branch length
    number_of_queries:      number of input queries
'''
def transfer_tree(tree, filtered_blast_result, treeRanks, ncbi_boolean, branch_length, number_of_queries):
    nodeID = None
    nodeName = None
    global count_clades
    ncbi = NCBITaxa()

    if ncbi_boolean == '0':  # input was a ncbi taxonomy
        nodeID = int(tree.name)
        nodeName = tree.sci_name

    else:  # input was own tree
        nodeID, nodeName = translate_node(tree.name.replace('_', ' '), count_clades)

    if len(tree.get_children()) == 0:  # hit a leaf
        if nodeID in filtered_blast_result.keys():
            if ncbi_boolean == '2': # taxa-based phylogeny
                tree.name = newick_valid_name(nodeName)
                return {'size': [15, branch_length], 'name': nodeName, 'value': [filtered_blast_result[nodeID], nodeID, 1], 'leaf_counter': 1}
            else:
                value = [[filtered_blast_result[nodeID][i], 0] for i in range(number_of_queries)]
                value.append(treeRanks[nodeID])

                return {'size': [15, branch_length], 'name': nodeName, 'value': value, 'leaf_counter': 1}

    else:  # inner node
        tree_children = list(
            filter(None, [transfer_tree(child, filtered_blast_result, treeRanks, ncbi_boolean, branch_length, number_of_queries) for child in tree.get_children()]))

        child_value = [[] for i in range(number_of_queries)]

        try:
            [child_value[i].append(filtered_blast_result[nodeID][i]) for i in range(number_of_queries)]
        except KeyError:
            [child_value[i].append(0) for i in range(number_of_queries)]

        subtree_hits = [0]*number_of_queries
        tree_size = 0
        leaf_counter = 0
        list_of_child_ids = []                       # phylogeny: all child ids to calculate the last common ancestor
        count_clade_members = 0                      # phylogeny: addition to cladeName to distinguish between clades
        for child in tree_children:
            if ncbi_boolean != '2':
                for i in range(number_of_queries):
                    subtree_hits[i] += sum(child['value'][i][:2])
            tree_size += child['size'][0]
            leaf_counter += child['leaf_counter']

            if ncbi_boolean == '2':
                count_clade_members += child['value'][2]
                if 'clade' not in child['name']:
                    list_of_child_ids.append(child['value'][1])

        # get the last common ancestor (LCA) for the inner nodes in the taxa-based phylogeny
        # same LCA for different branches are distinguished by numbers
        if ncbi_boolean == '2':
            child_value = child_value[0]
            if len(list_of_child_ids) == 1:
                try:
                    nodeName = ncbi.get_taxid_translator(list_of_child_ids)[list_of_child_ids[0]] + '_' + str(count_clade_members) + '.' + count_clades
                    count_clades += 1
                    nodeID = list_of_child_ids[0]
                except:
                    print('No name for:' + str(list_of_child_ids[0]))

            elif len(list_of_child_ids) > 0:
                try:
                    children_tree = ncbi.get_topology(list_of_child_ids)
                    nodeName = children_tree.sci_name + '_' + str(count_clade_members) + '.' + str(count_clades)
                    count_clades += 1
                    nodeID = int(children_tree.name)
                except:
                    print('No topology for :' + ','.join([str(item)  for item in list_of_child_ids]))
            tree.name = newick_valid_name(nodeName)
            child_value.append(nodeID)
            child_value.append(count_clade_members)

            if sum(child_value[:2]) != 0:
                return {'size': [tree_size + 10, branch_length], 'name': nodeName, 'children': tree_children,
                        'value': child_value, 'leaf_counter': leaf_counter}

        else:
            [child_value[i].append(subtree_hits[i]) for i in range(number_of_queries)]
            child_value.append(treeRanks[nodeID])

            if sum([sum(child_value[i][:2]) for i in range(number_of_queries)]) != 0:
                return {'size': [tree_size+10, branch_length], 'name': nodeName, 'children': tree_children,
                    'value': child_value, 'leaf_counter': leaf_counter}


'''recursively transfer the tree in d3 compatible format + other data for hit tree
    unique_count:   dic with hit accesion + counted taxa
'''
def unique_transfer_tree(tree, unique_count, branch_length):
    global count_clades
    if len(tree.name) == 0:
        tree.name = 'clade' + str(count_clades)
        count_clades += 1

    if len(tree.get_children()) == 0: # hit a leaf
        return {'size': [branch_length, branch_length], 'name': tree.name, 'value': [unique_count[tree.name], 0]}
    else: # inner node
        tree_children = list(filter(None, [unique_transfer_tree(child, unique_count, branch_length) for child in tree.get_children()]))
        tree_size = sum([child['size'][0] for child in tree_children])
        return {'size': [tree_size+branch_length, branch_length], 'name': tree.name,
                'children': tree_children, 'value': [sum([child['value'][0] for child in tree_children]), 0]}


''' Translate taxid and return sci_name and taxid
    node_label: given label of the node (taxid or sci_name)
    count:      for the phylogeny calculation count the number of unlabeled clades
'''
def translate_node(node_label, count):
    ncbi = NCBITaxa()  # setup NCBI database
    global count_clades

    if node_label == '':
        count_clades += 1
        return count, 'clade' + str(count)
    else:
        try:
            if node_label.isnumeric(): # scientific name missing
                nodeName = ncbi.get_taxid_translator([node_label])
                return int(node_label), nodeName[int(node_label)]
            else: # taxonomic id missing
                nodeId = ncbi.get_name_translator([node_label])
                return nodeId[node_label][0], node_label

        except KeyError: # taxon not present in stored NCBI database
            print('no valid key')
            count_clades += 1
            return count, 'clade' + str(count)


########################################################################################################################
#                                                                                         Pipeline for taxonomic Mapping
########################################################################################################################
'''
prot_data:                  aa sequence or blast result as xml file     
prot_file_type:             0= aa sequence, 1= xml file     
tree_data:                  newick string/file (user defined) or list of taxa
tree_menu:                  0= ncbi taxonomy, 1= user defined taxonomy
blast_type:                 default = blastp  (possibility to expend to blastx)
eValue:                     E-value for blast search threshold
min_align_ident:            minimal match between query and hit (percentage)
min_query_cover:            minimal query coverage (percentage)
min_hit_cover:              minimal hit coverage (percentage)
out_dir:	   
'''
def run_blastphylo(prot_data, prot_file_type, tree_data, tree_menu, blast_type, eValue, min_align_ident, min_query_cover,  min_hit_cover, out_dir):
    d3_tree = {}
    sequence_dic = {}
    uniqueAccs = {}
    number_of_queries = []

    # read tree
    #try:
    print('Start Tree calculation')
    tree, taxid_tree_set, entrez_query = read_tree_input(tree_data, tree_menu, True)
    print('complete')
    '''except:
        sys.stderr.write('Tree calculation failed')
        return d3_tree, sequence_dic, uniqueAccs, len(number_of_queries)
    
    # run BLAST
    try:'''
    print('Start Blast run')
    blast_result = run_blast(prot_data, prot_file_type, blast_type, eValue, min_align_ident, min_query_cover, min_hit_cover, entrez_query, out_dir)
    print('complete')
    '''except:
        sys.stderr.write('Blast run failed')
        return d3_tree, sequence_dic, uniqueAccs, len(number_of_queries)'''

    # filtering
    if blast_result.size > 0:
        #try:
        print('Start filtering')
        filtered_blast, sequence_dic, uniqueAccs, number_of_queries = filter_blast_result(blast_result, taxid_tree_set)
        print(number_of_queries)
        print('complete')
        '''except:
            sys.stderr.write('Filtering failed')
            return d3_tree, sequence_dic, uniqueAccs, len(number_of_queries)'''
    else:
        sys.stderr.write('Blast prefiltering lead to no hits')
        return d3_tree, sequence_dic, uniqueAccs, number_of_queries

    if filtered_blast:
        #try:
        print('Start conversion')
        d3_tree = transfer_tree_in_d3(tree, filtered_blast, tree_menu, len(number_of_queries))
        #generate_blastphylo_output(d3_tree, out_dir, len(number_of_queries))
        print('complete')
        return d3_tree, sequence_dic, uniqueAccs, number_of_queries
        ''''except:
            sys.stderr.write('Transfer in d3 compatible tree/newick failed')
            return d3_tree, sequence_dic, uniqueAccs, len(number_of_queries)'''
    else:
        sys.stderr.write('None of the BLAST hits were present in the tree')
        return d3_tree, sequence_dic, uniqueAccs, number_of_queries


########################################################################################################################
#                                                                                                    Output/Export files
########################################################################################################################
'''Generate Table like output of the d3 tree 
   d3_tree:             tree in d3 compatible format
   number_of_queries:   number of input queries
'''
def generate_tree_output(d3_tree, number_of_queries):
    global  normalization_dic

    if number_of_queries == 1:
        ncbi_nodes = normalization_dic[d3_tree['name']][1]
        percent_nodes = round((d3_tree['leaf_counter']/ncbi_nodes)*100, 2)
        try:
            return d3_tree['name'].replace(' ', '_') + ',' + d3_tree['value'][1].replace(' ', '_') + ',' + ','.join([str(val) for val in d3_tree['value'][0][:2]]) + ',' + \
                   str(d3_tree['leaf_counter']) + ',' + str(ncbi_nodes) + ',' + str(percent_nodes) +'\n' + '\n'.join([generate_tree_output(child, number_of_queries) for child in d3_tree['children']])
        except KeyError:
            return d3_tree['name'].replace(' ', '_') + ',' + d3_tree['value'][1].replace(' ', '_') + ',' + ','.join([str(val) for val in d3_tree['value'][0][:2]]) + ',' + \
                   str(d3_tree['leaf_counter']) + ',' + str(ncbi_nodes) + ',' + str(percent_nodes)
    else:
        ncbi_nodes = normalization_dic[d3_tree['name']][1]
        percent_nodes = round((d3_tree['leaf_counter']/ncbi_nodes)*100, 2)
        try:
            child_string = ','.join([str(val) for row in d3_tree['value'][:number_of_queries] for val in row])
            return d3_tree['name'].replace(' ', '_') + ',' +  d3_tree['value'][number_of_queries].replace(' ', '_')  + ',' + child_string + ',' +\
                   str(d3_tree['leaf_counter']) + ',' + str(ncbi_nodes) + ',' + str(percent_nodes) + '\n' + '\n'.join([generate_tree_output(child, number_of_queries) for child in d3_tree['children']])
        except KeyError:
            child_string = ','.join([str(val) for row in d3_tree['value'][:number_of_queries] for val in row])
            return d3_tree['name'].replace(' ', '_') + ',' + d3_tree['value'][number_of_queries].replace(' ', '_') + ',' + child_string + ',' + \
                   str(d3_tree['leaf_counter']) + ',' + str(ncbi_nodes) + ',' + str(percent_nodes)


''' Generate csv file with header for the d3-compatible tree
    tree:               d3-compatible tree structure (dictionary)
    outdir:             path to folder 
    number_of_queries:  number of input queries
'''
def generate_blastphylo_output(tree, number_of_queries):
    global normalization_dic

    if len(number_of_queries) == 1:
        query_string = '#hits,#subtree_hits,#nodesTree,#nodesNCBI,%nodes'
    else:
        query_header = [query+'#hits,'+ query+'#subtree' for query in number_of_queries]
        query_string = ','.join(query_header) + ',#nodesTree,#nodesNCBI,%nodes'

    actual_dir = os.getcwd()
    json_file = actual_dir[:-6] + '/src/data/ncbi_normalisation.json'
    with open(json_file, 'r') as f:
        normalization_dic = json.load(f)

    table_tree = 'Sci_Name,rank,' + query_string + '\n' + generate_tree_output(tree, len(number_of_queries))
    return table_tree


''' Generate Newick string from dictionary like tree
    tree:       d3-compatible tree structure (dictionary)
'''
def generate_Newick_from_dict(tree):
    if 'children' in tree.keys():
        return '(' + ''.join([generate_Newick_from_dict(child) for child in tree['children']])[:-1] + ')' + newick_valid_name(tree['name']) + ','
    else:
        return newick_valid_name(tree['name']) + ','
