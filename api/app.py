##/***
# * PhyloBlast
# * @author jennifer mueller
# *
# ***/

# packages
from __future__ import print_function
import sys
import json
import os
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import processing_data
from flask import Flask, render_template, redirect, url_for, request, jsonify, make_response
from werkzeug.utils import secure_filename
from flask_cors import CORS, cross_origin
import logging

# create tmp output folder for output file

actual_dir = os.getcwd()
flask_tmp_dir = actual_dir + '/flask_tmp'
shutil.rmtree(flask_tmp_dir, ignore_errors=True)
os.mkdir(flask_tmp_dir)

UPLOAD_FOLDER = '/api/flask_tmp'
ALLOWED_EXTENSIONS = set(['fasta', 'fastq', 'csv'])

# start server
app = Flask(__name__)
app.config['CORS_HEADERS'] = 'Content-Type'
CORS(app, resources={r"*": {"origins": "*"}})




@app.route('/api/phyloblast', methods=['POST', 'GET'])
def phyloblast():

    if request.method == 'POST':
        d3_tree = None
        tree_data = None
        tree_menu_selection = None
        print('Parameter of this run:')

        # Blast Search
        protein_seq = request.form['protein']
        protein_data = SeqRecord(Seq(protein_seq), id='query 1', description='protein sequence')
        SeqIO.write(protein_data, "flask_tmp/protein.fasta", "fasta")  # generate fasta file
        protein = 'flask_tmp/protein.fasta'

        protein_file_type = ''
        if len(protein_seq) == 0:
            protein_file_type = request.form['protein_file_type']
            if protein_file_type == '0': # fasta like file
                protein_data = request.files['fasta_file'].save('flask_tmp/protein.fasta')
                protein = 'flask_tmp/protein.fasta'
                print('Seq from file')
            elif protein_file_type == '1': # blast result xml file
                protein_data = request.files['fasta_file'].save('flask_tmp/blast_result.csv')
                protein = 'flask_tmp/blast_result.csv'
        print(protein)

        # tree
        tree_menu_selection = request.form['tree_menu_selection']

        # ncbi taxonomy
        if tree_menu_selection == '0':
            tree_data = request.form['tree_data']
        # user defined taxonomy
        elif tree_menu_selection == '1':
            tree_data = request.form['tree_data']
            if len(tree_data) == 0:
                tree_file = request.files['newick_file'].read().decode('utf-8')
                tree_data = tree_file.replace("\n", "")
        print(tree_data)

        # filter parameter
        eValue = request.form['eValue']
        print('E-value: ' + eValue)
        min_query_cover = request.form['query_cover']
        print('min. query coverage: ' + min_query_cover)
        min_query_identity = request.form['query_ident']
        print('min. query identity: ' +  min_query_identity)


        # start processing of the data
        # run_phyloblast(prot_data, prot_file_type, tree_data, tree_menu, blast_type, eValue, min_query_cover, min_identity, protein_length):

        d3_tree, seq_dic =  processing_data.run_phyloblast(protein, protein_file_type, tree_data, tree_menu_selection, 'blastp', eValue, min_query_cover, min_query_identity, 'flask_tmp/')

        print('Root of the tree: ' + d3_tree['name'])
        return d3_tree
    else:
        return None
      





