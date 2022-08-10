# BLASTphylo

BLASTphylo is an interactive web-tool, which applies BLAST to a given protein/gene sequence, maps the resulting hits onto a given taxonomic tree, and calculates a phylogeny, all in an automated pipeline. Taxonomic and phylogenetic trees are visualized interactively and allow the user to manipulate the tree on demand, e.g. collapsing the taxonomy to a certain taxonomic rank. If used locally please follow the installation instructions below. <br>
<br>
BLASTphylo is available as a web-tool an can be accessed here: [Visit BLASTphylo](https://tuevis.cs.uni-tuebingen.de/blastphylo/)

## Required Tools

BLASTphylo calculations require local installations of [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279690/) (suggest: v2.60+), [MAFFT](https://mafft.cbrc.jp/alignment/software/index.html) (suggest: v7.310) and [FastTree](http://www.microbesonline.org/fasttree/) (suggest: v2.1.10, linux recommand installation with `sudo apt-get install fasttree`) . <br>
Check your FastTree installation with `fasttree -help` if the command is not found check the installation again <br>
with `FastTree -help`. Depending on your FastTree installation open the **/server/external_tools.py** file and <br>
change the **cmd="fasttree"** parameter of `def __init__(self, cmd="fasttree", **kwargs):` into your FastTree <br>
command. 


## Installation Manual

This project was bootstrapped with [Create React App](https://github.com/facebook/create-react-app) and requires [npm](https://www.npmjs.com/get-npm) (suggest: node v14.15.0, npm v6.14.8) and [python](https://www.python.org/) (suggest: v3.6.9) for the installation

### Step 1: Create a venv in the parent directory of BLASTphylo

For example use:
`cd <path/to/BLASTphylo>` <br>
`virtualenv --python=/usr/bin/python3.6 <your-venv>`

### Step2: Installation of the required back-end python packages

 `source <your-venv>/bin/activate` <br>
 `pip3 install -r server/requirements.txt` <br>
 `pip install .`


### Step 3: Installations for the front-end

 `npm install`  <br>
command will end with `found 3 vulnerabilities (2 moderate, 1 high)`

### Step 4: Start BLASTphylo

Change into the project directory and run
 `npm run blastphylo`. Start the Flask server. 

Make sure that all required modules are installed before running the app.
Open [http://localhost:3000](http://localhost:3000) to open BLASTphylo in the browser. Enable the pop-up-window function for this app. 

## Available Examples

The project directory folder **server/test_examples** contains several example data which can be used as an input for BLASTphylo. <br>


## Help
For more details about the different features of BLASTphylo switch
to the **help page** in the top right corner of the application.
