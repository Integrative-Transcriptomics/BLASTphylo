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


# packages for Blast + Taxonomy mapping
from Bio.Blast.Applications import NcbiblastpCommandline as Blastp
from ete3 import NCBITaxa  # NCBI taxonomy (localy stored, require ~300MB)
# ncbi.update_taxonomy_database() # update actual NCBI taxonomy version

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
branch_length = 30 

# functions
######################################################################################################################## BLAST
''' run and parse the blast result dependent on the server input 
    prot         amino acid sequence or direct Blast result as xml file
    prot_file_type:     0= aa sequence, 1= xml file
    blast_type:         blastp or blastx dependent on the input (by now only blastp)
    eValue:             significant threshold for the blast search
    min_identity:       minimal identity with the query
    taxa:               root of the given taxonomy
    out_dir:            path to output directory

    output: Tabular only WP identifier lesser hits as with XML file ?
    tab sep: # Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score, % positives
'''
def run_blast(prot, prot_file_type, blast_type, eValue, min_query_cover, out_dir):
    header = ['qacc', 'sacc', 'qstart', 'qend', 'evalue', 'pident', 'staxids', 'qcovs', 'sseq']

    if prot_file_type == "1":
        result = pd.read_csv(prot, sep='\t', names=header)
        # print(result.shape)
    else:
        blast_out_path = out_dir + 'blast_result.csv'
        blastp_cline = Blastp(cmd=blast_type, remote=True, query=prot, db='refseq_protein', evalue=eValue,
                              qcov_hsp_perc=min_query_cover,
                              outfmt='6 qacc sacc qstart qend evalue pident staxids qcovs sseq', out=blast_out_path)
        print(blastp_cline)
        stdout, stderr = blastp_cline()
        print(stderr)
        result = pd.read_csv(blast_out_path, sep='\t', names=header)

    return result


######################################################################################################################## Filtering and Hit count
'''Filter blast result for hits similar to the species in the tree
    blast_record:      result of the blast search
    tree_taxIDs:       set of all taxonomic IDs in the given taxonomy
    min_query_ident:   minimal query ident of the blast result to be counted as hit 
'''
def filter_blast_result(blast_record, tree_tax_IDs, min_query_ident):
    filtered_result = {}
    sequence_dic = {}
    filtered_blast = blast_record[blast_record['pident'] > (int(min_query_ident) / 100)]

    ncbi = NCBITaxa()
    for index, row in filtered_blast.iterrows():  # taxids from the hit table with number of occurence
        for start_taxa in row['staxids'].split(';'):
            if int(start_taxa) in tree_tax_IDs:
                try:
                    filtered_result[int(start_taxa)] += 1
                    sequence_dic[int(start_taxa)][1].append(row['sseq'])
                    sequence_dic[int(start_taxa)][0].append(row['sacc'])
                except KeyError:
                    filtered_result[int(start_taxa)] = 1
                    sequence_dic[int(start_taxa)] = [[row['sacc']], [row['sseq']]]

            else:
                taxID_lineage = ncbi.get_lineage(start_taxa)
                for line_taxa in taxID_lineage:  # all higher level taxa of these taxids
                    if line_taxa in tree_tax_IDs:

                        try:
                            filtered_result[line_taxa] += 1
                            sequence_dic[line_taxa][1].append(row['sseq'])
                            sequence_dic[line_taxa][0].append(row['sacc'])
                        except KeyError:
                            filtered_result[line_taxa] = 1
                            sequence_dic[line_taxa] = [[row['sacc']], [row['sseq']]]
    return filtered_result, sequence_dic



######################################################################################################################## sequence-based Phylogeny
'''
    calculate consensus sequence for each taxID (majority vote)
    seqs:           dictionary with sequence array per taxID
    output_file:    path to fasta for aligned sequences
'''
def generate_consensus(seqs, output_file):
    consensus_seqs = []
    for taxID in seqs.keys():
        consensus_seqs.append(SeqRecord(Seq(seqs[taxID][1][0]), id=str(taxID), description=''))
    SeqIO.write(consensus_seqs, output_file, "fasta")

'''
   transfer Mafft alignment in taxid version 
   mafft_out:           path to mafft output
   taxid_accession:     dictionary with subject aligned seqs and their accessions
'''
def generate_fasttree_input(mafft_out, taxid_accession):
    new_mafft = []
    records = SeqIO.parse(mafft_out, 'fasta')
    count = 0
    for record in records:
        try:
            [new_mafft.append(SeqRecord(record.seq, id=str(taxId), description='')) for taxId in taxid_accession[record.id[:-2]]]
        except KeyError:
            count += 1
    SeqIO.write(new_mafft, mafft_out.split('.')[0] + '_NEW.fasta', 'fasta')

    return None

''' calculate phylogeny for given hit sequences
    subject_seqs:         dictionary with subject aligned seqs and their accessions
    output_dir:           path to output directory
    from_aligned_seq:     boolean, true: aligned subject sequences, false: full subject sequences

    based on: #https://biopython.org/wiki/Phylo
'''


def calculate_phylogeny(subject_seqs, fasta, output_dir, from_aligned_seq):
    output_seqs = output_dir + "blast_subject_seqs.fasta"

    if from_aligned_seq == 'True':
        print('Calculate consensus')
        generate_consensus(subject_seqs, output_seqs)
    else:
        output_seqs = fasta
        print('Input was fasta file')
        # list_of_accessions = [subject_seqs[key][0][0] for key in subject_seqs.keys()]
        # blastdbcmd_cline = Blastdbcmd(db="/share/references/library/refseq_protein/refseq_protein", entry=','.join(list_of_accessions), out=output_seqs)  #db path = default path to refseq database on Romanov
        # blastdbcmd_stdout, blastdbcmd_stderr = blastdbcmd_cline()

    # calculate MSA
    output_fasta = output_dir + "mafft_aligned.fasta"

    mafft_cline = Mafft(input=output_seqs, amino=True, clustalout=False,
                        parttree=True)  # treeout --> guidetree in output
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
        output_fasta = output_dir + "mafft_aligned_NEW.fasta"

    # generate phylogenetic tree
    output_tree = output_dir + "fasttree.tree"

    fasttree_cline = FastTree(input=output_fasta,
                              out=output_tree, fastest=True)  # fastest --> faster calculation by focus on best hit
    print(fasttree_cline)
    stdout_FastTree, stderr_FastTree = fasttree_cline()

    tree, tree_tax_ids, _ = read_tree_input(output_tree, '2')

    pseudo_filtered_blast = {}  # instead of the hit counts the parent taxid is saved
    parentnode_info = {}
    hitAcc_info = {}

    ncbi = NCBITaxa()
    for key in list(tree_tax_ids):
        try:
            key_lineage = ncbi.get_rank(ncbi.get_lineage(key))
            parent_of_key = [key for key in key_lineage if
                             key_lineage[key] == 'phylum']  # last entry in lineage is the taxid itself
            pseudo_filtered_blast[key] = parent_of_key[0]
            parentID, parentName = translate_node(str(parent_of_key[0]), 0)
            parentnode_info[parentID] = parentName  # parentnode_info contain parent sci name + taxid for visualisation
        except:
            pseudo_filtered_blast[key] = 1
        hitAcc_info[key] = subject_seqs[key][0][0]  # subject_seqs[key][0][0] == accession number of the used hit

    d3_phylogeny = wrapper_transfer_own_tree(tree, pseudo_filtered_blast, '2')

    return d3_phylogeny, parentnode_info, hitAcc_info


######################################################################################################################## Tree handling
'''read tree input 
    tree_input:     tree as newick string, file or string of sientific names/IDs to generate a tree from the NCBI taxonomy
    ncbi_boolean:   boolean to seperate the two input types (0 = NCBI, 1 = own taxonomy) 
                    NCBI example: Staphylococcus,Staphylococcus aureus|subtree
'''
def read_tree_input(tree_input, ncbi_boolean):
    tree = None
    root = None
    tree_taxIDs = set()
    if ncbi_boolean == '0':
        ncbi_taxa = []
        subtrees = tree_input.split(',')
        root = subtrees[0].split('|')[0]
        for node in subtrees:  # user could select the complete subtree
            ncbi = NCBITaxa()
            descendants = node.split('|')
            if len(descendants) > 1:
                ncbi_taxa.append(descendants[0])
                tree_taxIDs.update(ncbi.get_descendant_taxa(descendants[0],
                                                            intermediate_nodes=True))  # output of descendant is a list of ids
            else:
                tree_taxIDs.append(descendants[0])

        tree_taxIDs = tree_taxIDs.union(translate_nodes(ncbi_taxa))

    elif ncbi_boolean == '1':
        tree = Tree(tree_input, format=1)
        root = tree.name.replace('_', ' ')
        tree_taxIDs = translate_nodes([node.name for node in tree.traverse('preorder')])

    elif ncbi_boolean == '2':
        tree = Tree(tree_input)
        root = tree.name.replace('_', ' ')
        tree_taxIDs = translate_nodes([node.name for node in tree.traverse('preorder')])

    else:
        sys.stderr.write('no tree given')

    return tree, tree_taxIDs, translate_node(root, 0)[1]


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
def transfer_tree_in_d3(tree, filtered_blast_result, ncbi_boolean):
    d3_tree = None
    if ncbi_boolean == '0':
        d3_tree = wrapper_transfer_ncbi_tree(filtered_blast_result, ncbi_boolean)
    elif ncbi_boolean == '1':
        d3_tree = wrapper_transfer_own_tree(tree, filtered_blast_result, ncbi_boolean)
    else:
        sys.stderr.write('check the input again')
    return d3_tree


''' generate tree and transfer it into d3 format
    filtered_blast_result:  result of the blast filtering is used to generate the minimal tree to connact the nodes
    ncbi_boolean:           see transfer_tree_in_d3       
'''
def wrapper_transfer_ncbi_tree(filtered_blast_result, ncbi_boolean):
    global branch_length
    ncbi = NCBITaxa()  # setup NCBI database
    tree = ncbi.get_topology(filtered_blast_result.keys())

    root_children = list(
        filter(None, [transfer_tree(child, filtered_blast_result, ncbi_boolean) for child in tree.get_children()]))

    root_value = []  # contain the hits for the root and the sum of all children hits
    try:
        root_value.append(filtered_blast_result[int(tree.name)])
    except KeyError:
        root_value.append(0)

    subtree_hits = 0
    root_size = 0
    for child in root_children:
        subtree_hits += sum(child['value'])
        root_size += child['size'][1]
    root_value.append(subtree_hits)
    return {'size': [branch_length, root_size], 'name': tree.sci_name, 'children': root_children,
            'value': root_value}


''' transfer own tree with possible mixed labeling in d3 format
    tree:                  given tree structure
    filtered_blast_result: dictionary of hits  
    ncbi_boolean:           see transfer_tree_in_d3
'''
def wrapper_transfer_own_tree(tree, filtered_blast_result, ncbi_boolean):
    global count_clades
    global branch_length

    rootID, rootName = translate_node(tree.name.replace('_', ' '), count_clades)
    root_children = list(
        filter(None, [transfer_tree(child, filtered_blast_result, ncbi_boolean) for child in tree.get_children()]))

    root_value = []  # contain the hits for the root and the sum of all children hits

    try:
        root_value.append(filtered_blast_result[rootID])
    except KeyError:
        root_value.append(0)


    subtree_hits = 0
    root_size = 0
    for child in root_children:
        subtree_hits += sum(child['value'])
        root_size += child['size'][1]
    if ncbi_boolean == '2':
        root_value.append(rootID)
    else:
        root_value.append(subtree_hits)

    return {'size': [branch_length, root_size], 'name': rootName, 'children': root_children,
            'value': root_value}


''' recusively transfer the subtrees
    tree:                   subtree
    filtered_blast_result:  dictionary of hits
    ncbi_boolean:           dependent on the tree format the node information is different
'''
def transfer_tree(tree, filtered_blast_result, ncbi_boolean):
    nodeID = None
    nodeName = None
    global count_clades
    global branch_length

    if ncbi_boolean == '0':  # input was a ncbi taxonomy
        nodeID = int(tree.name)
        nodeName = tree.sci_name

    else:  # input was own taxonomic phylogeny / seq-based phylogeny
        nodeID, nodeName = translate_node(tree.name.replace('_', ' '), count_clades)


    if len(tree.get_children()) == 0:  # hit a leaf
        if nodeID in filtered_blast_result.keys():
            if ncbi_boolean == '2':
                return {'size': [branch_length, 10], 'name': nodeName, 'value': [filtered_blast_result[nodeID], nodeID]}
            else:
                return {'size': [branch_length, 10], 'name': nodeName, 'value': [filtered_blast_result[nodeID], 0]}

    else:  # inner node
        tree_children = list(
            filter(None, [transfer_tree(child, filtered_blast_result, ncbi_boolean) for child in tree.get_children()]))
        child_value = []

        try:
            child_value.append(filtered_blast_result[nodeID])
        except KeyError:
            child_value.append(0)

        subtree_hits = 0
        tree_size = 0
        for child in tree_children:
            subtree_hits += sum(child['value'])
            tree_size += child['size'][1]
        if ncbi_boolean == '2':
            child_value.append(nodeID)
        else:
            child_value.append(subtree_hits)

        if sum(child_value) != 0:
            return {'size': [branch_length, tree_size+10], 'name': nodeName, 'children': tree_children,
                    'value': child_value}


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
            if node_label.isnumeric():
                return int(node_label), ncbi.get_taxid_translator([node_label])[int(node_label)]
            else:
                return ncbi.get_name_translator([node_label])[node_label][0], node_label
        except KeyError:
            # print('no valid key')
            count_clades += 1
            return count, 'clade' + str(count)


######################################################################################################################## Pipeline
'''
prot_data:          aa sequence or blast result as xml file     
prot_file_type:     0= aa sequence, 1= xml file     
tree_data:          newick string/file (user defined) or list of taxa
tree_menu:          0= ncbi taxonomy, 1= user defined taxonomy
blast_type:         default = blastp  (possibility to expend to blastx)
eValue:             E-value for blast search threshold
min_query_cover:    minimal query coverage (percentage)
min_identity:       minimal match between query and target (percentage)
out_dir:	   
'''
def run_phyloblast(prot_data, prot_file_type, tree_data, tree_menu, blast_type, eValue, min_query_cover, min_identity, out_dir):
    d3_tree = None
    sequence_dic = None

    # read tree
    try:
        print('Start Tree calculation')
        tree, taxid_tree_set, root_taxa = read_tree_input(tree_data, tree_menu)
        print('complete')
    except:
        sys.stderr.write('Tree calculation failed')
        return None, None

    # run BLAST
    try:
        print('Start Blast run')
        blast_result = run_blast(prot_data, prot_file_type, blast_type, eValue, min_query_cover, out_dir)
        print('complete')
    except:
        sys.stderr.write('Blast run failed')
        return None, None

    # filtering
    try:
        print('Start filtering')
        filtered_blast, sequence_dic = filter_blast_result(blast_result, taxid_tree_set, min_identity)
        print('complete')
    except:
        sys.stderr.write('Filtering failed')
        return None, None

    if filtered_blast:
        try:
            d3_tree = transfer_tree_in_d3(tree, filtered_blast, tree_menu)
            return d3_tree, sequence_dic
        except:
            sys.stderr.write('Transfer in d3 compatible tree/newick failed')
            return None, None
