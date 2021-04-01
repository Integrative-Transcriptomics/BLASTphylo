# BLASTphylo

BLASTphylo is an interactive web-tool, which applies a BLASTp search for a given protein sequence and maps the result on a given taxonomic tree. In addition, a phylogeny calculation of the remaining BLAST hits can be performed and visualised.


## Required Tools

BLASTphylo calculations require local installations of [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279690/) (suggest: v2.60+), [MAFFT](https://mafft.cbrc.jp/alignment/software/index.html) PartTree (suggest: v7.310) and [FastTree](http://www.microbesonline.org/fasttree/) Double-precision (suggest: v2.1.10) . <br>
Check your FastTree installation with `fasttree -help` if the command is not found check the installation again <br>
with `FastTree -help`. Depending on your FastTree installation open the **/server/external_tools.py** file and <br>
change the **cmd="fasttree"** parameter of `def __init__(self, cmd="fasttree", **kwargs):` into your FastTree <br>
command. 


## Installation Manual

This project was bootstrapped with [Create React App](https://github.com/facebook/create-react-app) and require [npm](https://www.npmjs.com/get-npm) (suggest: node v14.15.0, npm v6.14.8) and [python](https://www.python.org/) (suggest: v3.6.9) for the installation

### Step 1: Installation of the required back-end python packages

Change in the `server` directory and run 
#### `pip install requirements.txt`  


### Step 2: Installation of the front-end

Change in the project directory and run 
#### `npm install` 

### Step 3: Start BLASTphylo

Change in the project directory and run
#### `npm run dev`
Make sure that all needed modules are installed before running the app.
Runs the app and the server in the background.
Open [http://localhost:3000](http://localhost:3000) to view BLASTphylos 
*menu* page in the browser. <br>
Enable the pop-up-window function for this app. 

## Available Examples

In the project directory folder **test examples** contain a example BLAST
result, which was used to implement BLASTphylo.



## Help
For more details about the different features of BLASTphylo switch
to the **help page** in the top right corner of the *menu* page.
