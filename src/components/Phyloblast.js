// used libraries
import React, {Component} from 'react';
import * as d3v6 from 'd3v6';
import {Link, BrowserRouter, Route, Switch } from "react-router-dom";


// own components and style sheets
import './phyloblastStyle.css';
import {ExportTrees} from '../App.js'
import Phylogeny from './Phylogeny.js';


// own visualisations
import {chart, hitBars, startTreevis, collapseTree, publicationReady} from '../visualisations/phyloblast.js';

class Phyloblast extends Component {
    constructor(props) {
        super(props);

        this.state = {hitSelect: "2",
                      rankSelect: 'class'};

        this.handleHits = this.handleHits.bind(this);
        this.handleRanks = this.handleRanks.bind(this);

        d3v6.select('#treeVis').style('border', '2px solid #69a2c9')
                                .style('border-radius', '5px');

        // generation of the tree visualisation
        var tree = startTreevis(this.props.phyloData);
        if (tree !== 0){
            chart(tree, null, 6, 6, true);
        }else{
          d3v6.select('#phyloblastAlert').remove();
          d3v6.select('#tree').append('div').attr('id', 'phyloblastAlert')
                    .text('Found 0 hits. Return to the main page and try another phylogentic tree');
        }


   }

   handleHits(event){
     this.setState({hitSelect: event.target.value});

     if(event.target.value !== '2'){
        hitBars();
     }else{
        d3v6.select('#hitbars')
        .remove();
     }
   }

   handleRanks(event){
      if(document.getElementById('public_ready')){
              this.setState({rankSelect: event.target.value});
              collapseTree(this.props.phyloData);
      }
   }

   render(){
    var taxonomyLevel = ['life', 'domain', 'superkingdom', 'kingdom', 'clade', 'phylum', 'class', 'order', 'family', 'genus', 'species group','species', 'strain'];
    var MakeItem = function(X){
        return <option value={X}>{X}</option>;
    };


    return(
       <div id="phyloblast">

          <div id="phyloblastMenu">
              <h2>taxonomic Analysis</h2>
              <BrowserRouter>
              <ul>
              <Link to='/phylogeny' target='_blank'>taxa seqs. based phylogeny </Link>
              <Link to='/phylogenyUnique' target='_blank'>unique seqs. based phylogeny</Link>
              </ul>
              <Switch>
                    <Route exact path="/phylogny" render={(props) => <Phylogeny {...props} changeComp={this.changeComp} />}/>
                    <Route exact path="/phylogenyUnique" component={Phylogeny}/>
              </Switch>
              </BrowserRouter>
              <br></br>
              <label>Collapse to:</label>
              <select id="collapse_menu" name="collapse_menu" value={this.state.rankSelect} onChange={this.handleRanks}>
              {taxonomyLevel.map(MakeItem)}
              </select>

              <label>Show barchart for:</label>
               <select id="tree_menu" name="tree_menu" value={this.state.hitSelect} onChange={this.handleHits} >
                <option value="2">none</option>
                <option value="0">node hits</option>
                <option value="1">subtree hits</option>
               </select>
               <button id="public_ready" title="Result will be a static tree. Reload the page to start a new interpolation of the basic tree" onClick={publicationReady}>publication ready</button>
               <ExportTrees />

          </div>
     </div>
    );
    }
}

export default Phyloblast;
