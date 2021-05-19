##/***
# * BLASTphylo
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
import logging, click
import time


######################################################################################################################## general functions and settings
# create tmp output folder for output file
actual_dir = os.getcwd()
flask_tmp_dir = actual_dir + '/flask_tmp'
shutil.rmtree(flask_tmp_dir, ignore_errors=True)

try:                          # during first installation need to generate flask_tmp folder and download NCBI database
    os.mkdir(flask_tmp_dir)
    print('Download NCBI database')
    from ete2 import NCBITaxa
    ncbi = NCBITaxa()
except:
    print('flask_tmp folder is already generated')
UPLOAD_FOLDER = '/server/flask_tmp'


# generate fasta file from textfield input and check for amino acid sequence (can handle any number of queries)
def generate_fasta_from_input(textfieldinput, outdir, blast_type):
    fastas = []
    valid_pro_seq = True
    regAA = "^[ACDEFGHIKLMNPQRSTVWY]*$"
    regNuc = "^[ACTG]*$"

    split_input = textfieldinput.split('>')

    for fasta in split_input:
        if len(fasta) > 0:
            split_fasta = fasta.split('\n')
            seq = ''.join(split_fasta[1:]).replace('\r', '')

            if blast_type == 'blastp':  # check for valid protein sequence
                valid_pro_seq = re.search(regAA, seq)
            elif blast_typ == 'blastx' or blast_type == 'blastn': # check for valid nucleotid sequence
                valid_pro_seq = re.search(regNuc, seq)

            if valid_pro_seq:
                fastas.append(SeqRecord(Seq(seq), id=split_fasta[0], description=''))
            else:
                return valid_pro_seq
    SeqIO.write(fastas, outdir, "fasta")
    return valid_pro_seq

# start server
app = Flask(__name__)

# global parameter for the run
hit_seqs = {}
accs_seqs = {}
d3_tree = {}


######################################################################################################################## connection to the front end
# taxa-based phylogeny
@app.route('/server/phylogeny', methods=['POST', 'GET'])
def phylogeny():
    global hit_seqs
    output_tree = "flask_tmp/fasttree.tree"
    if os.path.isfile(output_tree): # generate only data for the visualization
        tree, tree_tax_ids,_ = processing_data.read_tree_input(output_tree, '2', True)
        d3Tree, phylum_info, acc_info = processing_data.phylogeny_data(tree, hit_seqs,tree_tax_ids)
        newick_phylogeny = tree.write(format=3)

    else: # perform complete calculation
        print('Start Phylogeny calculation')
        d3Tree, newick_phylogeny, phylum_info, acc_info = processing_data.calculate_phylogeny(hit_seqs, None, 'flask_tmp/', 'True',  False)
    return {'tree': d3Tree, 'newick': newick_phylogeny[:-1], 'extraInfo': [phylum_info, acc_info]}

# unique sequence-based phylogeny
@app.route('/server/phylogenyUnique', methods=['POST', 'GET'])
def phylogenyUnique():
    global accs_seqs
    output_tree = "flask_tmp/fasttree_unique.tree"
    if os.path.isfile(output_tree): # generate only data for the visualization
        tree, tree_tax_ids, _ = processing_data.read_tree_input(output_tree, '2', True)
        d3_phylogeny, newick_phylogeny = processing_data.unique_phylogeny_data(tree, accs_seqs)

    else: # perform complete calculation
        print('Start Phylogeny calculation')
        d3_phylogeny, newick_phylogeny = processing_data.calculate_phylogeny(accs_seqs, None, 'flask_tmp/', 'True',  True)

    return {'tree': d3_phylogeny, 'newick': newick_phylogeny[:-1]}


# data generation for taxonomic mapping
@app.route('/server/menu', methods=['POST', 'GET'])
def menu():
    # remove all old files
    [os.remove(os.path.join('flask_tmp/', f)) for f in os.listdir('flask_tmp')]

    if request.method == 'POST':
        global hit_seqs
        global accs_seqs
        global d3_tree
        tree_data = None
        tree_menu_selection = None
        error = []
        print('Parameter of this run:')

        # Blast Search
        blasttype = request.form['blasttype']
        print('BLAST search:   ' + blasttype)
        protein_seq = request.form['protein']
        if len(protein_seq) > 0:
            if '>' not in protein_seq:
                error.append({'message': 'Protein sequence(s) does not contain > identifier'})
                protein = None
            valid_pro_seq = generate_fasta_from_input(protein_seq, "flask_tmp/protein.fasta", blasttype)
            if valid_pro_seq:
                protein = 'flask_tmp/protein.fasta'
            else:
                error.append({'message': 'Protein sequence(s) contain irregular amino acids'})
                protein = None

        protein_file_type = ''
        if len(protein_seq) == 0:
            protein_file_type = request.form['protein_file_type']
            if protein_file_type == '0': # fasta like file
                protein_data = request.files['fasta_file'].save('flask_tmp/protein.fasta')
                protein = 'flask_tmp/protein.fasta'
                print('Sequence from file')
            elif protein_file_type == '1': # blast result csv file
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
                tree = Tree(tree_data, format=8)
                #print(tree)
            except:
                error.append({'message': 'Uploaded tree file contain a BLASTphylo incompatible format. Newick string need labels for all nodes.'})
                tree_data = None

        print(tree_data)

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
            error.append({'message': 'Given taxa are not present in NCBI taxonomy'})
        print(error)

        if len(error) > 0:
            return {'error': error}
        else:
            # start processing of the data
            print('\nStart PhyloBlast')
            try:
                d3_tree, hit_seqs, accs_seqs = processing_data.run_blastphylo(protein, protein_file_type, tree_data, tree_menu_selection, blasttype, eValue, min_align_identity, min_query_cover, min_hit_cover, 'flask_tmp/')
                #print(d3_tree)
                if len(d3_tree) > 0:
                    print('Root of the tree: ' + d3_tree['name'])
                    return {'tree': d3_tree, 'error': None}
                else:
                    return {'tree': None, 'error': None}
            except:
                return {'tree': None, 'error': None}
    else:
        return None
      





