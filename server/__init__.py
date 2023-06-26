##/***
# * BLASTphylo
# * @authors Jennifer Müller, Susanne Zabel
# *
# ***/

# packages
import sys
import json
import os
import shutil, re
import copy

# check menu inputs and preprocess them
import pandas as pd
import tempfile
from io import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ete3 import Tree
from server.processing_data import translate_nodes, run_blastphylo, calculate_phylogeny, generate_blastphylo_output, generate_Newick_from_dict

# flask server packages
from flask import Flask, render_template, redirect, url_for, request, jsonify, make_response



app = Flask(__name__, static_folder='../build', static_url_path='/')
here = os.path.dirname(__file__)
print('here', here)
# general functions and settings
# create tmp output folder for output file
flask_tmp_dir = here + '/flask_tmp/'
#flask_tmp_dir = tempfile.TemporaryDirectory()
# shutil.rmtree(flask_tmp_dir, ignore_errors=True)
app.config['NCBI_TAXONOMY_SEARCHBAR_ENTRIES'] = os.path.join(here, 'data', 'searchbar_entries.json')
app.config['TEST_EXAMPLE_MAPPING'] = os.path.join(here, 'test_example/mapping_example/blast_result_sada.csv')
app.config['TEST_EXAMPLE_COMPARISON'] = os.path.join(here, 'test_example/comparison_example/blast_result_mpsA_B.csv')
try:                          # during first installation need to generate flask_tmp folder and download NCBI database
    os.mkdir(flask_tmp_dir)
except:
    print('flask_tmp folder is already generated')


# generate fasta file from textfield input and check for amino acid sequence (can handle any number of queries)
def generate_fasta_from_input(textfieldinput, outdir, blast_type):
    fastas = []
    valid_pro_seq = True
    query_header = []
    regAA = "^[ACDEFGHIKLMNPQRSTVWY]*$"
    regNuc = "^[ACTG]*$"

    split_input = textfieldinput.split('>')
    print('splitinput', split_input)

    for fasta in split_input:
        print('fasta', fasta)
        if len(fasta) > 0:
            if '\n' in fasta:
                split_fasta = fasta.split('\n')
                seq = ''.join(split_fasta[1:]).replace('\r', '')
                if seq =='':
                    return False, query_header
                if blast_type == 'blastp':  # check for valid protein sequence
                    valid_pro_seq = re.search(regAA, seq)
                elif blast_type == 'blastx' or blast_type == 'blastn': # check for valid nucleotid sequence
                    valid_pro_seq = re.search(regNuc, seq)

                if valid_pro_seq:
                    query_header.append(split_fasta[0])
                    fastas.append(SeqRecord(Seq(seq), id=split_fasta[0], description=''))
                else:
                    print(valid_pro_seq)
                    return valid_pro_seq, query_header
            else:
                return False, query_header
    SeqIO.write(fastas, outdir, "fasta")
    return valid_pro_seq, query_header


# global parameter for the run and for the export of the data
hit_seqs = {}
accs_seqs = {}
d3_tree = {}
queries = []
taxa_newick = ''
unique_newick = ''


######################################################################################################################## connection to the front end
# start the flask server
@app.route('/')
def index():
    return app.send_static_file('index.html')


# accession for data export
@app.route('/server/exportData', methods=['POST', 'GET'])
def exportData():
    global d3_tree
    global queries
    global taxa_newick
    global unique_newick

    if request.method == 'POST':
        print(request.form['datatype'])
        datatype = request.form['datatype']

        if len(d3_tree) > 0:
            if 'mapping.csv' in datatype:
                table_tree = generate_blastphylo_output(d3_tree, queries)
                return {'data': table_tree, 'data_type': 'table'}
            elif 'newick_taxonomic' in datatype:
                d3_newick = generate_Newick_from_dict(d3_tree)
                d3_newick = d3_newick[:-1] + ';'
                return {'data': d3_newick, 'data_type': 'newick'}
            elif 'taxa_based' in datatype:
                return {'data': taxa_newick, 'data_type': 'newick'}
            elif 'unique' in datatype:
                return {'data': unique_newick, 'data_type': 'newick'}
            else:
                print('All data are ready for the download')
                return {'data': '', 'data_type': 'error'}
        else:
            print('No hits found')
            return {'data': '', 'data_type': 'error'}
    else:
        return None


# search for taxon
@app.route('/server/searchNcbiTaxa', methods=['POST'])
def searchNcbiTaxa():
    if request.method == 'POST':
        searchquery = request.form['searchquery']
        ncbi_taxa_file = app.config['NCBI_TAXONOMY_SEARCHBAR_ENTRIES']

        with open(ncbi_taxa_file, 'r') as f:
            ncbi_taxa = json.load(f)

        matching_taxa = [taxon for taxon in ncbi_taxa if searchquery.lower() in taxon['name'].lower()]

        return {'result': matching_taxa}
    else:
        return {'result': None}


# taxa-based phylogeny
@app.route('/server/phylogeny', methods=['POST', 'GET'])
def phylogeny():

    global hit_seqs
    global taxa_newick
    errors = []
    d3Tree = None
    newick_phylogeny=None,
    phylum_info = None,

    print('Start taxa-based phylogeny calculation')
    data, errors = calculate_phylogeny(hit_seqs, None, flask_tmp_dir, 'True',  False)
    if len(errors)>0:
        return {'tree': None, 'newick': None, 'extraInfo': None,
           'actualTree': None, 'errors': errors}
    else:
        d3Tree, newick_phylogeny, phylum_info, acc_info = data
        taxa_newick = newick_phylogeny

        return {'tree': d3Tree, 'newick': newick_phylogeny[:-1], 'extraInfo': [phylum_info, acc_info],
           'actualTree': copy.deepcopy(d3Tree), 'errors': errors}
    #return {'tree': d3Tree, 'newick': newick_phylogeny[:-1], 'extraInfo': [phylum_info, acc_info],
    #        'actualTree': copy.deepcopy(d3Tree)}

# unique sequence-based phylogeny
@app.route('/server/phylogenyUnique', methods=['POST', 'GET'])
def phylogenyUnique():
    global accs_seqs
    global unique_newick
    print('Start unique sequence-based phylogeny calculation')
    data, errors = calculate_phylogeny(accs_seqs, None, flask_tmp_dir, 'True',  True)
    if len(errors)>0:
        return {'tree': None, 'newick': None, 'extraInfo': [0],
                'actualTree': None, 'errors': errors}
    else:
        d3_phylogeny, newick_phylogeny = data
        unique_newick = newick_phylogeny

        return {'tree': d3_phylogeny, 'newick': newick_phylogeny[:-1], 'extraInfo': [0],
            'actualTree': copy.deepcopy(d3_phylogeny), 'errors': errors}


# data generation for taxonomic mapping
@app.route('/server/menu', methods=['POST'])
def menu():
    # remove all old files + set global variable

    if request.method == 'POST':
        global hit_seqs
        global accs_seqs
        global d3_tree
        global queries
        hit_seqs = {}
        accs_seqs = {}
        d3_tree = {}
        queries = []
        taxa_newick = ''
        unique_newick = ''
        tree_data = None
        tree_menu_selection = None
        error = []
        protein = None
        print('Parameter of this run:')

        # Blast Search
        blasttype = request.form['blasttype']
        print('BLAST search:   ' + blasttype)
        protein_seq = request.form['protein']
        if len(protein_seq) > 0:
            if '>' not in protein_seq:
                error.append({'message': 'Protein sequence(s) does not contain > identifier.'})
                protein = None
            # generation of tmp file
            temp_prot_file = tempfile.NamedTemporaryFile(prefix=flask_tmp_dir)
            valid_pro_seq, query_header = generate_fasta_from_input(protein_seq, temp_prot_file.name, blasttype)
            queries = query_header
            if valid_pro_seq:
                protein = temp_prot_file.name
            else:
                error.append({'message': 'One of the queries contains no or irregular letters.'})
                protein = None


        protein_file_type = ''
        fasta_filename = ''
        if len(protein_seq) == 0:
            protein_file_type = request.form['protein_file_type']
            fasta_filename = request.form['fasta_filename']
            if protein_file_type == '0': # fasta like file
                protein_data = request.files['fasta_file'].read()

                # generation of tmp file
                temp_prot_file = tempfile.NamedTemporaryFile(prefix=flask_tmp_dir)
                temp_prot_file.write(protein_data)
                temp_prot_file.read()
                protein = temp_prot_file.name
                print('Sequence from file')
            elif protein_file_type == '1': # blast result csv file
                protein_data = request.files['fasta_file'].read().decode('utf-8')
                protein = protein_data
                blastResult = pd.read_csv(StringIO(protein), sep='\t')
                if len(blastResult.columns) != 13:
                    blastResult = pd.read_csv(StringIO(protein), sep=',')
                queries_np = blastResult.iloc[:,0].unique()
                queries = queries_np.tolist()
                if len(blastResult.columns) != 13:
                    error.append({'message': 'Uploaded BLAST result is incorrectly formatted. Check the help page for more details.'})
                    protein = None
                print('Data input: CSV table')
            elif protein_file_type == '2':
                if fasta_filename == '2':
                    protein = app.config['TEST_EXAMPLE_MAPPING']
                else :
                    protein = app.config['TEST_EXAMPLE_COMPARISON']



        # tree
        tree_menu_selection = request.form['tree_menu_selection']
        print(tree_menu_selection)
        # ncbi taxonomy
        taxaIDs = None
        if tree_menu_selection == '0':
            print('tree menu selection')
            tree_data = request.form['tree_data']

            # check if given scientific names and taxon IDs are valid
            tree_taxa = tree_data.split(',')
            taxa = [taxon.split('|')[0] for taxon in tree_taxa]
            try:
                print('try')
                subtree_strings = [taxon.split('|')[1] for taxon in tree_taxa]
                # check if string behind | is subtree
                for i in subtree_strings:
                    if i != 'subtree':
                        error.append({'message': 'String behind | is not subtree.'})
            except:
                print('except')
            taxaIDs = translate_nodes(taxa)

        # user defined taxonomy
        elif tree_menu_selection == '1':
            tree_data = request.form['tree_data']
            if len(tree_data) == 0:
                tree_file = request.files['newick_file'].read().decode('utf-8')
                tree_data = tree_file.replace("\n", "")
            try:
                tree = Tree(tree_data, format=8)
            except:
                error.append({'message': 'Uploaded tree file contains a BLASTphylo incompatible format. Newick string needs labels for all nodes.'})
                tree_data = None


        # filter parameter
        eValue = request.form['eValue']
        print('E-value: ' + eValue)
        min_align_identity = request.form['align_ident']
        print('min. alignment identity: ' +  min_align_identity)
        min_query_cover = request.form['query_cover']
        print('min. query coverage: ' + min_query_cover)
        min_hit_cover = request.form['hit_cover']
        print('min. hit coverage: ' + min_hit_cover)

        if tree_menu_selection == '0' and len(taxaIDs) == 0:
            error.append({'message': 'Given taxa are not present in NCBI taxonomy.'})

        if len(error) > 0:
            return {'error': error}
        else:
            # start processing of the data
            print('\nStart BLASTphylo')
            try:
                d3_tree, d3_no_hit_tree, hit_seqs, accs_seqs, queries, blasterrors = run_blastphylo(protein, protein_file_type, tree_data, tree_menu_selection, blasttype, eValue, min_align_identity, min_query_cover, min_hit_cover, flask_tmp_dir)
                if fasta_filename == '2':
                    queries = ['sada']
                elif fasta_filename == '3':
                    queries = ['mpasA', 'mpsB']

                # remove temporary filesvalue
                if protein_file_type == '0':
                    temp_prot_file.close()
                if len(d3_tree) > 0:
                    print('Root of the tree: ' + d3_tree['name'])
                    return {'tree': d3_tree, 'noHitTree': d3_no_hit_tree, 'error': None, 'queries': queries, 'actualTree': d3_tree}
                else:
                    return {'tree': None, 'noHitTree': None, 'error': blasterrors[0], 'queries': queries, 'actualTree': None}
            except:
                if protein_file_type == '0':
                    temp_prot_file.close()
                return {'tree': None, 'noHitTree': None, 'error': 'undefined', 'queries': queries, 'actualTree': None}
    else:
        return None



# start flask server
if __name__ == '__main__':
    app.config['TEMPLATES_AUTO_RELOAD'] = True
    app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0
    app.run()

