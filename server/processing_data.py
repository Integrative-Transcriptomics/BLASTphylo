# /***
# * PhyloBlast
# * processing of the data
# *
# * @author jennifer mueller
# *
# ***/


# packages
import sys
import pandas as pd
import numpy as np


# packages for Blast + Taxonomy mapping
from Bio.Blast.Applications import NcbiblastpCommandline as Blastp
from ete3 import NCBITaxa  # NCBI taxonomy (localy stored, require ~300MB)
#ncbi.update_taxonomy_database() # update actual NCBI taxonomy version

from ete3 import Tree  # handle trees


# packages for Phylogeny calculation
from external_tools import NCBIblastdbcmdCommandline as Blastdbcmd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MafftCommandline as Mafft
from external_tools import FastTreeCommandline as FastTree


# global variables 
count_clades = 0

# functions prot_data, prot_file_type, blast_type, eValue, min_align_ident, min_query_cover, min_hit_cover, entrez_query, out_dir
######################################################################################################################## BLAST
''' run and parse the blast result dependent on the server input 
    prot                    amino acid sequence or direct Blast result as xml file
    prot_file_type:         0= aa sequence, 1= xml file
    blast_type:             blastp or blastx dependent on the input (by now only blastp)
    eValue:                 significant threshold for the blast search
    min_align_ident         minimal identity between query and subject sequence in the alignment
    min_query_cover:        minimal coverage with the query
    min_subject_cover:      minimal coverage with the subject (hit)
    entrez_query:           restrict the blast search to 'roots' of the given taxonomy
    out_dir:                path to output directory

    output: Tabular only WP identifier lesser hits as with XML file ?
    tab sep: # Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score, % positives
'''
def run_blast(prot, prot_file_type, blast_type, eValue, min_align_ident, min_query_cover, min_subject_cover, entrez_query, out_dir):
    header = ['qacc', 'sacc', 'qstart', 'qend', 'sstart', 'send', 'slen', 'nident', 'evalue', 'pident', 'staxids', 'qcovhsp', 'sseq']

    header_basic = ['qacc', 'sacc', 'pident', 'alen', 'mm', 'g', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bit', 'm']

    if prot_file_type == "1":
        try:
            preFilter = pd.read_csv(prot, sep='\t')
        except:
            preFilter = pd.read_csv(prot, sep=',')

        if len(preFilter.columns) == len(header):
            preFilter.columns = header
        else:
            preFilter.columns = header_basic

        subjectAlignedLength = abs(preFilter['send']-preFilter['sstart'])
        #print(subjectAlignedLength)
        subjectCoverage = (subjectAlignedLength/preFilter['slen']) > (int(min_subject_cover)/100)
        #print(subjectCoverage)
        alignIdent = preFilter['pident'] > (float(min_align_ident))
        #print(alignIdent)
        #print(preFilter['pident'])
        result = preFilter[(preFilter['evalue'] < float(eValue)) & (preFilter['qcovhsp'] > (float(min_query_cover))) & subjectCoverage & alignIdent]

        print(result.shape)

    else:
        blast_out_path = out_dir + 'blast_result.csv'
        blastp_cline = Blastp(cmd=blast_type, remote=True, query=prot, db='nr', evalue=eValue, max_hsps=1, num_alignments=1000,
                              qcov_hsp_perc=min_query_cover, entrez_query='\'' + entrez_query + '\'',
                              outfmt='6 qacc sacc qstart qend sstart send slen nident evalue pident staxids qcovhsp sseq', out=blast_out_path)
        print(blastp_cline)
        stdout, stderr = blastp_cline()
        print(stderr)
        preFilter = pd.read_csv(blast_out_path, sep='\t', names=header)

        subjectAlignedLength = abs(preFilter['send']-preFilter['sstart'])
        #print(subjectAlignedLength)
        subjectCoverage = (subjectAlignedLength/preFilter['slen']) > (float(min_subject_cover))
        #print(subjectCoverage)
        alignIdent = preFilter['pident'] > (int(min_align_ident)/100)

        result = preFilter[(preFilter['evalue'] < float(eValue)) & (preFilter['qcovhsp'] > (float(min_query_cover))) & subjectCoverage & alignIdent]
        print(result.shape)

    return result


######################################################################################################################## Filtering and Hit count
'''Filter blast result for hits similar to the species in the tree
    blast_record:      result of the blast search
    tree_taxIDs:       set of all taxonomic IDs in the given taxonomy
'''
def filter_blast_result(blast_record, tree_tax_IDs):
    filtered_result = {}
    sequence_dic = {}  # taxid encoded seq dic
    uniqueAccs = {}  # accession ID encoded seq dic
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

                try:
                    filtered_result[int(start_taxa)][number_of_queries.index(row['qacc'])] += 1
                    sequence_dic[int(start_taxa)][1].append(row['sseq'])
                    sequence_dic[int(start_taxa)][0].append(row['sacc'])
                except KeyError:
                    filtered_result[int(start_taxa)] = [0]*len(number_of_queries)
                    filtered_result[int(start_taxa)][number_of_queries.index(row['qacc'])] = 1
                    sequence_dic[int(start_taxa)] = [[row['sacc']], [row['sseq']]]

            else:
                try:
                    taxID_lineage = ncbi.get_lineage(start_taxa)
                except:
                    continue
                uniqueAccs[row['sacc']] = [row['sseq'], 0]
                for line_taxa in taxID_lineage:  # all higher level taxa of these taxids
                    if line_taxa in tree_tax_IDs:
                        uniqueAccs[row['sacc']][1] += 1
                        try:
                            filtered_result[line_taxa][number_of_queries.index(row['qacc'])] += 1
                            sequence_dic[line_taxa][1].append(row['sseq'])
                            sequence_dic[line_taxa][0].append(row['sacc'])
                        except KeyError:
                            filtered_result[line_taxa] = [0] * len(number_of_queries)
                            filtered_result[line_taxa][number_of_queries.index(row['qacc'])] = 1
                            sequence_dic[line_taxa] = [[row['sacc']], [row['sseq']]]
                if uniqueAccs[row['sacc']][1] == 0:
                    uniqueAccs.pop(row['sacc'])

    return filtered_result, sequence_dic, uniqueAccs, number_of_queries



######################################################################################################################## sequence-based Phylogeny
'''
    extract best hit for each taxID
    seqs:           dictionary with sequence array per taxID
    output_file:    path to fasta for aligned sequences
'''
def extract_best_hit(seqs, output_file):
    best_seqs = []
    [best_seqs.append(SeqRecord(Seq(seqs[taxID][1][0]), id=str(taxID), description='')) for taxID in seqs.keys()]
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
    output_dir:           path to output directory
    from_aligned_seq:     boolean, true: aligned subject sequences, false: full subject sequences
    uniqueSeqs:           boolean, true: accession encoded seqs, false: taxid encoded seqs
    based on: #https://biopython.org/wiki/Phylo
'''
def calculate_phylogeny(subject_seqs, fasta, output_dir, from_aligned_seq, uniqueSeqs):
    output_seqs = output_dir + "blast_subject_seqs.fasta"

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
        # list_of_accessions = [subject_seqs[key][0][0] for key in subject_seqs.keys()]
        # blastdbcmd_cline = Blastdbcmd(db="/share/references/library/refseq_protein/refseq_protein", entry=','.join(list_of_accessions), out=output_seqs)  #db path = default path to refseq database on Romanov
        # blastdbcmd_stdout, blastdbcmd_stderr = blastdbcmd_cline()

    # calculate MSA
    mafft_file = "mafft_aligned_unique.fasta" if uniqueSeqs else "mafft_aligned.fasta"
    output_fasta = output_dir + mafft_file

    mafft_cline = Mafft(input=output_seqs, amino=True, clustalout=False,
                        parttree=True, thread=4)  # treeout --> guidetree in output
    print(mafft_cline)
    stdout_Mafft, stderr_Mafft = mafft_cline()  # stdout_Mafft = MSA
    with open(output_fasta, "w") as handle:
        handle.write(stdout_Mafft)
    handle.close()

    if from_aligned_seq == 'False':  # need to convert the Mafft full seq output to a similar shape as for aligned seqs
        acc_taxids = {}
        for key in subject_seqs:
            try:
                acc_taxids[subject_seqs[key][0][0]].append(key)
            except KeyError:
                acc_taxids[subject_seqs[key][0][0]] = [key]
        generate_fasttree_input(output_fasta, acc_taxids)
        output_fasta = output_dir + mafft_file.split('.')[0] + '_NEW.fasta'

    # generate phylogenetic tree
    fasttree_file = "fasttree_unique.tree" if uniqueSeqs else "fasttree.tree"
    output_tree = output_dir + fasttree_file

    fasttree_cline = FastTree(input=output_fasta, quiet=True,
                              out=output_tree, fastest=True)  # fastest --> faster calculation by focus on best hit
    print(fasttree_cline)
    stdout_FastTree, stderr_FastTree = fasttree_cline()

    print('Start Tree calculation and processing data')
    tree, tree_tax_ids, _ = read_tree_input(output_tree, '2', True)
    global count_clades
    if uniqueSeqs:

        count_clades = 0

        return unique_phylogeny_data(tree, subject_seqs, output_tree)
    else:
        count_clades = 0
        d3Tree, parentnode_info, mosthitAcc = phylogeny_data(tree, subject_seqs, tree_tax_ids, output_tree, True)
        #newickTree = phylogeny_data(tree, subject_seqs, tree_tax_ids, output_tree, False)
        newickTree = tree.write(format=3)
        #print(newickTree)
        return d3Tree, newickTree, parentnode_info, mosthitAcc

'''
    Taxon name have to be in a newick string valid form
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
    treefile:       path to FastTree output
    d3_version:     boolean, True: d3 ouput needed, False: Phylotree output
'''
def phylogeny_data(tree, subject_seqs, treeIDs, treefile, d3_version):
    ncbi = NCBITaxa()

    pseudo_filtered_blast = {}  # instead of the hit counts the parent taxid is saved
    parentnode_info = {}

    if d3_version:
        for key in list(treeIDs):
            defaultvalue = 1
            try:
                key_lineage = ncbi.get_rank(ncbi.get_lineage(key))
                parent_of_key = [key for key in key_lineage if
                                 key_lineage[key] == 'phylum']
                pseudo_filtered_blast[key] = [parent_of_key[0]]
                parentID, parentName = translate_node(str(parent_of_key[0]), 0)
                parentnode_info[parentID] = parentName  # parentnode_info contain parent sci name + taxid for visualisation
            except:
                pseudo_filtered_blast[key] = [defaultvalue]

        accvalues = [subject_seqs[key][0][0] for key in subject_seqs.keys()]
        #uniqueSeqs = [key for key in subject_seqs.keys() if len(subject_seqs[key][0]) == 1]

        uniqueSeqs = [key for key in subject_seqs.keys() if accvalues.count(subject_seqs[key][0][0]) == 1]
        newTree = wrapper_transfer_own_tree(tree, pseudo_filtered_blast, '2', 10, 1)
        return newTree, parentnode_info, uniqueSeqs

    '''else:
        with open(treefile, 'r') as treehandle:
            newickTree = treehandle.read()
        treehandle.close()
        newTree = replace_ID_Name(newickTree, treeIDs)
        return newTree'''


''' 2. Phylogeny: unique hit based
    tree, tree file:    see 1. Phylogeny
    uniqueAccs:         hit accession encoded dic with seqs + #taxa connected to hit
'''
def unique_phylogeny_data(tree, uniqueAccs, treefile):
    unique_count = {key: uniqueAccs[key][1] for key in uniqueAccs.keys()}
    d3_tree = unique_transfer_tree(tree, unique_count, 20)
    newickTree = tree.write(format=3)
    '''with open(treefile, 'r') as treehandle:
        newickTree = treehandle.read()
    treehandle.close() '''

    return d3_tree, newickTree

######################################################################################################################## Tree handling
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
        ncbi_taxa = []
        subtrees = tree_input.split(',')
        for node in subtrees:  # user could select the complete subtree
            ncbi = NCBITaxa()
            descendants = node.split('|')
            if len(descendants) > 1:
                ncbi_taxa.append(descendants[0])
                tree_taxIDs.update(ncbi.get_descendant_taxa(descendants[0],
                                                            intermediate_nodes=True))  # output of descendant is a list of ids
            else:
                ncbi_taxa.append(descendants[0])
        rootIDs = translate_nodes(ncbi_taxa)
        root = ['txid' + str(taxon) + '[ORGN]' for taxon in list(rootIDs)]
        roots = ' OR '.join(root)
        tree_taxIDs = tree_taxIDs.union(translate_nodes(ncbi_taxa))

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
            else:
                try:
                    translated_node = ncbi.get_name_translator([node.replace('_', ' ')])
                    translated_set.add(translated_node[node.replace('_', ' ')][0])
                except KeyError:
                    count_not_valid_ids += 1
    print('During translation of the tree nodes ' + str(count_not_valid_ids) + ' could not be translated')
    return translated_set


################################ d3 compatible format #############################################
'''Transfer the given tree in a d3 compatible hierarchical format
    tree:                   input tree 
    filtered_blast_result:  filtered blast results with found hits 
    ncbi_boolean:           dependent on the input tree the nodes could have a mixed labeling 
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
'''
def wrapper_transfer_ncbi_tree(filtered_blast_result, ncbi_boolean, branch_length, number_of_queries):

    ncbi = NCBITaxa()  # setup NCBI database
    taxa_keys = filtered_blast_result.keys()
    tree = ncbi.get_topology(taxa_keys)
    treeTaxIDs = translate_nodes([node.name for node in tree.traverse('preorder')])
    treeRanks = ncbi.get_rank(list(treeTaxIDs))

    if len(taxa_keys) == 1:
        root_value = [[filtered_blast_result[int(tree.name)][i], 0] for i in range(number_of_queries)]
        root_value.append(treeRanks[int(tree.name)])
        return {'size': [15, branch_length], 'name': tree.sci_name, 'value': root_value}

    else:
        root_children = list(
            filter(None, [transfer_tree(child, filtered_blast_result, treeRanks, ncbi_boolean, branch_length, number_of_queries) for child in tree.get_children()]))

        root_value = [[] for i in range(number_of_queries)]  # contain the hits for the root and the sum of all children hits
        try:
            [root_value[i].append(filtered_blast_result[int(tree.name)][i]) for i in range(number_of_queries)]
        except KeyError:
            [root_value[i].append(0) for i in range(number_of_queries)]

        subtree_hits = [0]*number_of_queries
        root_size = 0
        for child in root_children:
            for i in range(number_of_queries):
                subtree_hits[i] += sum(child['value'][i][:2])
            root_size += child['size'][0]


        if ncbi_boolean == '2':
            root_value.append(int(tree.name))
        else:
            [root_value[i].append(subtree_hits[i]) for i in range(number_of_queries)]
            root_value.append(treeRanks[int(tree.name)])

        return {'size': [root_size, branch_length], 'name': tree.sci_name, 'children': root_children,
                'value': root_value}


''' transfer own tree with possible mixed labeling in d3 format
    tree:                  given tree structure
    filtered_blast_result: dictionary of hits  
    ncbi_boolean:           see transfer_tree_in_d3
'''
def wrapper_transfer_own_tree(tree, filtered_blast_result, ncbi_boolean, branch_length, number_of_queries):
    global count_clades
    ncbi = NCBITaxa()
    taxa_keys = filtered_blast_result.keys()
    treeRanks = ncbi.get_rank(taxa_keys)
    rootID, rootName = translate_node(tree.name.replace('_', ' '), count_clades)

    if len(taxa_keys) == 1:
        root_value = [[filtered_blast_result[rootID][i], 0] for i in range(number_of_queries)]
        root_value.append(treeRanks[rootID])
        return {'size': [15, branch_length], 'name': tree.sci_name, 'value': root_value}

    else:
        root_children = list(
            filter(None, [transfer_tree(child, filtered_blast_result, treeRanks, ncbi_boolean, branch_length, number_of_queries) for child in tree.get_children()]))

        root_value = [[] for i in range(number_of_queries)]  # contain the hits for the root and the sum of all children hits
        try:
            [root_value[i].append(filtered_blast_result[rootID][i]) for i in range(number_of_queries)]
        except KeyError:
            [root_value[i].append(0) for i in range(number_of_queries)]


        subtree_hits = [0]*number_of_queries
        root_size = 0
        list_of_child_ids = []  # phylogeny: all child ids to calculate the last common ancestor
        count_clade_members = 0  # phylogeny: addition to cladeName to distinguish between clades
        for child in root_children:
            if ncbi_boolean != '2':
                for i in range(number_of_queries):
                    subtree_hits[i] += sum(child['value'][i][:2])
            root_size += child['size'][0]
            if ncbi_boolean == '2':
                count_clade_members += child['value'][2]
                if 'clade' not in child['name']:
                    list_of_child_ids.append(child['value'][1])

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
                'value': root_value}


''' recursively transfer the subtrees
    tree:                   subtree
    filtered_blast_result:  dictionary of hits
    ncbi_boolean:           dependent on the tree format the node information is different
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
                return {'size': [15, branch_length], 'name': nodeName, 'value': [filtered_blast_result[nodeID], nodeID, 1]}
            else:
                value = [[filtered_blast_result[nodeID][i], 0] for i in range(number_of_queries)]
                value.append(treeRanks[nodeID])

                return {'size': [15, branch_length], 'name': nodeName, 'value': value}

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
        list_of_child_ids = []                       # phylogeny: all child ids to calculate the last common ancestor
        count_clade_members = 0                      # phylogeny: addition to cladeName to distinguish between clades
        for child in tree_children:
            if ncbi_boolean != '2':
                for i in range(number_of_queries):
                    subtree_hits[i] += sum(child['value'][i][:2])
            tree_size += child['size'][0]

            if ncbi_boolean == '2':
                count_clade_members += child['value'][2]
                if 'clade' not in child['name']:
                    list_of_child_ids.append(child['value'][1])

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
                        'value': child_value}

        else:
            [child_value[i].append(subtree_hits[i]) for i in range(number_of_queries)]
            child_value.append(treeRanks[nodeID])

            if sum([sum(child_value[i][:2]) for i in range(number_of_queries)]) != 0:
                return {'size': [tree_size+10, branch_length], 'name': nodeName, 'children': tree_children,
                    'value': child_value}


'''recursively transfer the tree in d3 compatible format + other data for hit tree
    unique_count:   dic with hit accesion + counted taxa
'''
def unique_transfer_tree(tree, unique_count, branch_length):
    global count_clades
    if len(tree.name) == 0:
        tree.name = 'clade' + str(count_clades)
        count_clades += 1

    if len(tree.get_children()) == 0:
        return {'size': [branch_length, branch_length], 'name': tree.name, 'value': [unique_count[tree.name], 0]}
    else:
        tree_children = list(filter(None, [unique_transfer_tree(child, unique_count, branch_length) for child in tree.get_children()]))
        tree_size = sum([child['size'][0] for child in tree_children])
        return {'size': [tree_size+branch_length, branch_length], 'name': tree.name,
                'children': tree_children, 'value': [sum([child['value'][0] for child in tree_children]), 0]}


''' Translate taxid and return sci_name and taxid
    node_label: given label of sthe node (taxid or sci_name)
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
            if node_label.isnumeric():
                nodeName = ncbi.get_taxid_translator([node_label])
                return int(node_label), nodeName[int(node_label)]
            else:
                nodeId = ncbi.get_name_translator([node_label])
                return nodeId[node_label][0], node_label
        except KeyError:
            print('no valid key')
            count_clades += 1
            return count, 'clade' + str(count)


######################################################################################################################## Pipeline
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
def run_phyloblast(prot_data, prot_file_type, tree_data, tree_menu, blast_type, eValue, min_align_ident, min_query_cover,  min_hit_cover, out_dir):
    d3_tree = {}
    sequence_dic = {}
    uniqueAccs = {}

    # read tree
    #try:
    print('Start Tree calculation')
    tree, taxid_tree_set, entrez_query = read_tree_input(tree_data, tree_menu, True)
    print('complete')
    '''except:
        sys.stderr.write('Tree calculation failed')
        return d3_tree, sequence_dic, uniqueAccs'''
    
    # run BLAST
    #try:
    print('Start Blast run')
    blast_result = run_blast(prot_data, prot_file_type, blast_type, eValue, min_align_ident, min_query_cover, min_hit_cover, entrez_query, out_dir)
    print('complete')
    '''except:
        sys.stderr.write('Blast run failed')
        return d3_tree, sequence_dic, uniqueAccs '''

    # filtering
    if blast_result.size > 0:
        try:
            print('Start filtering')
            filtered_blast, sequence_dic, uniqueAccs, number_of_queries = filter_blast_result(blast_result, taxid_tree_set)
            print('complete')
        except:
            sys.stderr.write('Filtering failed')
            return d3_tree, sequence_dic, uniqueAccs
    else:
        sys.stderr.write('Blast prefiltering lead to no hits')
        return d3_tree, sequence_dic, uniqueAccs

    if filtered_blast:
        #try:
        print('Start conversion')
        d3_tree = transfer_tree_in_d3(tree, filtered_blast, tree_menu, len(number_of_queries))
        generate_phyloblast_output(d3_tree, out_dir, len(number_of_queries))
        return d3_tree, sequence_dic, uniqueAccs
        ##except:
         #   sys.stderr.write('Transfer in d3 compatible tree/newick failed')
         #   return d3_tree, sequence_dic, uniqueAccs
    else:
        sys.stderr.write('Non of the BLAST hits were present in the tree')
        return d3_tree, sequence_dic, uniqueAccs


######################################################################################################################## Output/Export files
'''Generate Table like output of the d3 tree 
   d3_tree: tree in d3 compatible format
'''
def generate_tree_output(d3_tree, number_of_queries):

    if number_of_queries == 1:
        try:
            return d3_tree['name'].replace(' ', '_') + ',' + d3_tree['value'][1].replace(' ', '_') + ',' + ','.join([str(val) for val in d3_tree['value'][0][:2]]) + '\n' + '\n'.join([generate_tree_output(child, number_of_queries) for child in d3_tree['children']])
        except KeyError:
            return d3_tree['name'].replace(' ', '_') + ',' + d3_tree['value'][1].replace(' ', '_') + ',' + ','.join([str(val) for val in d3_tree['value'][0][:2]])
    else:
        try:
            child_string = ','.join([str(val) for row in d3_tree['value'][:number_of_queries] for val in row])
            return d3_tree['name'].replace(' ', '_') + ',' + d3_tree['value'][number_of_queries].replace(' ', '_') + ',' + child_string + '\n' + '\n'.join([generate_tree_output(child, number_of_queries) for child in d3_tree['children']])
        except KeyError:
            child_string = ','.join([str(val) for row in d3_tree['value'][:number_of_queries] for val in row])
            return d3_tree['name'].replace(' ', '_') + ',' + d3_tree['value'][number_of_queries].replace(' ', '_') + ',' + child_string


def generate_phyloblast_output(tree, outdir, number_of_queries):
    if number_of_queries == 1:
        query_string = '#hits,#subtree_hits'
    else:
        query_header = ['q'+str(i)+'#hits,q'+str(i)+'#subtree' for i in range(number_of_queries)]
        query_string = ','.join(query_header)

    table_tree = 'Sci_Name,rank,' + query_string + '\n' + generate_tree_output(tree, number_of_queries)
    tree_out = outdir + 'taxonomicMapping.csv'
    with open(tree_out, 'w') as w:
        w.write(table_tree)

''' Generate newick string from dictionary
    d3_tree:    tree in d3 compatible format
'''
def generate_newick_from_dic(d3_tree):
    if 'children' in d3_tree.keys():
        return (','.join([generate_newick_from_dic(child) for child in de_tree['children']])) + d3_tree['name']
    else:
        return d3_tree['name']

