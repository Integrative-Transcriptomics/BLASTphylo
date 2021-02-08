##/***
# * PhyloBlast
# * @author jennifer mueller
# *
# ***/

# packages
import sys
import json
import os
import shutil, re

# check menu inputs and preprocess them
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ete3 import Tree
import processing_data

# flask server packages
from flask import Flask, render_template, redirect, url_for, request, jsonify, make_response
from werkzeug.utils import secure_filename
from flask_cors import CORS, cross_origin
import logging

# create tmp output folder for output file

actual_dir = os.getcwd()
flask_tmp_dir = actual_dir + '/flask_tmp'
shutil.rmtree(flask_tmp_dir, ignore_errors=True)
os.mkdir(flask_tmp_dir)

UPLOAD_FOLDER = '/server/flask_tmp'

# start server
app = Flask(__name__)
app.config['CORS_HEADERS'] = 'Content-Type'
CORS(app, resources={r"*": {"origins": "*"}})

hit_seqs = {}
accs_seqs = {}

@app.route('/server/phylogeny', methods=['POST', 'GET'])
def phylogeny():
    global hit_seqs
    output_tree = "flask_tmp/fasttree.tree"
    if os.path.isfile(output_tree):
        tree, tree_tax_ids,_ = processing_data.read_tree_input(output_tree, '2', True)
        d3Tree, phylum_info, acc_info = processing_data.phylogeny_data(tree, hit_seqs,tree_tax_ids, output_tree, True)
        newick_phylogeny = processing_data.phylogeny_data(tree, hit_seqs, tree_tax_ids, output_tree, False)
    else:
        print('Start Phylogeny calculation')
        d3Tree, newick_phylogeny, phylum_info, acc_info = processing_data.calculate_phylogeny(hit_seqs, None, 'flask_tmp/', 'True',  False)

    return {'tree': d3Tree, 'newick': newick_phylogeny[:-1], 'extraInfo': [phylum_info, acc_info]}

@app.route('/server/phylogenyUnique', methods=['POST', 'GET'])
def phylogenyUnique():
    global accs_seqs
    output_tree = "flask_tmp/fasttree_unique.tree"
    if os.path.isfile(output_tree):
        tree, tree_tax_ids, _ = processing_data.read_tree_input(output_tree, '2', True)
        d3_phylogeny, newick_phylogeny = processing_data.unique_phylogeny_data(tree, accs_seqs, output_tree)
    else:
        print('Start Phylogeny calculation')
        d3_phylogeny, newick_phylogeny = processing_data.calculate_phylogeny(accs_seqs, None, 'flask_tmp/', 'True',  True)

    return {'tree': d3_phylogeny, 'newick': newick_phylogeny[:-1]}

@app.route('/server/phyloblast', methods=['POST', 'GET'])
def phyloblast():
    # remove all old files
    [os.remove(os.path.join('flask_tmp/', f)) for f in os.listdir('flask_tmp')]

    if request.method == 'POST':
        global hit_seqs
        global accs_seqs
        d3_tree = None
        tree_data = None
        tree_menu_selection = None
        error = []
        print('Parameter of this run:')
        regAA = "^[ACDEFGHIKLMNPQRSTVWY]*$"

        # Blast Search
        protein_seq = request.form['protein']
        valid_pro_seq = re.search(regAA, protein_seq)
        if valid_pro_seq:
            protein_data = SeqRecord(Seq(protein_seq), id='query 1', description='protein sequence')
            SeqIO.write(protein_data, "flask_tmp/protein.fasta", "fasta")  # generate fasta file
            protein = 'flask_tmp/protein.fasta'
        else:
            error.append({'message': 'Protein sequence contain irregular amino acids'})
            protein = None

        protein_file_type = ''
        if len(protein_seq) == 0:
            protein_file_type = request.form['protein_file_type']
            if protein_file_type == '0': # fasta like file
                protein_data = request.files['fasta_file'].save('flask_tmp/protein.fasta')
                protein = 'flask_tmp/protein.fasta'
                print('Sequence from file')
            elif protein_file_type == '1': # blast result xml file
                protein_data = request.files['fasta_file'].save('flask_tmp/blast_result.csv')
                protein = 'flask_tmp/blast_result.csv'
                try:
                    blastResult = pd.read_csv(protein, sep='\t')
                except:
                    blastResult = pd.read_csv(protein, sep=',')
                if len(blastResult.columns) != 13:
                    error.append({'message': 'Uploaded BLAST result has to much/less columns. Check the help page for '
                                             'more details.'})
                    protein = None

        print(protein)

        # tree
        tree_menu_selection = request.form['tree_menu_selection']

        # ncbi taxonomy
        taxaIDs = None
        if tree_menu_selection == '0':
            tree_data = request.form['tree_data']

            # check if given scientific names and taxon IDs are valid
            tree_taxa = tree_data.split(',')
            taxa = [taxon.split('|')[0] for taxon in tree_taxa]
            taxaIDs = processing_data.translate_nodes(taxa)

        # user defined taxonomy
        elif tree_menu_selection == '1':
            tree_data = request.form['tree_data']
            if len(tree_data) == 0:
                tree_file = request.files['newick_file'].read().decode('utf-8')
                tree_data = tree_file.replace("\n", "")
            try:
                print(tree_data)
                tree = Tree(tree_data)
                print(tree)
            except:
                error.append({'message': 'Uploaded tree file contain a no valid newick string'})
                tree_data = None

        print(tree_data)

        # filter parameter
        eValue = request.form['eValue']
        print('E-value: ' + eValue)
        min_query_identity = request.form['query_ident']
        print('min. query identity: ' +  min_query_identity)
        min_query_cover = request.form['query_cover']
        print('min. query coverage: ' + min_query_cover)
        min_hit_identity = request.form['hit_ident']
        print('min. hit identity: ' +  min_hit_identity)
        min_hit_cover = request.form['hit_cover']
        print('min. hit coverage: ' + min_hit_cover)

        if tree_menu_selection == '0' and len(taxaIDs) == 0:
            error.append({'message': 'Given taxa are not present in NCBI taxonomy'})
        print(error)

        if len(error) > 0:
            return {'error': error}
        else:
            # start processing of the data
            # run_phyloblast(prot_data, prot_file_type, tree_data, tree_menu, blast_type, eValue, min_query_cover, min_identity, protein_length):
            print('\nStart PhyloBlast')
            try:
                d3_tree, hit_seqs, accs_seqs = processing_data.run_phyloblast(protein, protein_file_type, tree_data, tree_menu_selection, 'blastp', eValue, min_query_cover, min_query_identity, min_hit_cover, min_hit_identity, 'flask_tmp/')
                #print(d3_tree)
                if len(d3_tree) > 0:
                    processing_data.generate_phyloblast_output(d3_tree, 'flask_tmp/')  # generate mapping output as csv

                    print('Root of the tree: ' + d3_tree['name'])
                    return {'tree': d3_tree, 'error': None}
                else:
                    return {'tree': None, 'error': None}
            except:
                return {'tree': None, 'error': None}
    else:
        return None
      





