# BLASTphylo

This project was bootstrapped with [Create React App](https://github.com/facebook/create-react-app).

## Required Tools

BLASTphylo calculations require local installations of [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279690/), [MAFFT]() PartTree 
and [FastTree](http://www.microbesonline.org/fasttree/) Double-precision. <br>
Check your FastTree installation with `fasttree -help` if the command is not found check the installation again <br>
with `FastTree -help`. Depending on your FastTree installation open the **/server/external_tools.py** file and <br>
change the **cmd="fasttree"** parameter of `def __init__(self, cmd="fasttree", **kwargs):` into your FastTree <br>
command. 


## Available Scripts

In the project directory, you can run: 

### `npm ci`

Runs a clean install of all needed modules

### `npm run dev`
Make sure that all needed modules are installed before running the app.
Runs the app and the server in the background.
Open [http://localhost:3000](http://localhost:3000) to view BLASTphylos 
*menu* page in the browser. <br>
Enable the pop-up-window function for this app. 

## Available Examples

In the project directory folder **test examples** contain two example BLAST
results, which are used to implement BLASTphylo.

Example 1: mapping based <br>
Example 2: phylogeny example 

## Help
For more details about the different features of BLASTphylo switch
to the **help page** in the top right corner of the *menu* page.