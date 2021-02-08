// used libraries
import React, {Component} from 'react';
import * as d3v6 from 'd3v6';
import * as d3 from 'd3';
import * as $ from "jquery";
import * as _ from "lodash";
import "../../node_modules/phylotree/src/main";

// own components and style sheets
import './phylogenyStyle.css';
import {ExportTrees} from '../App.js'

// other style sheets
import '../../node_modules/phylotree/phylotree.css'

// own visualisations
import {chart, showClades, startTreevis, publicationReady} from '../visualisations/phyloblast.js';


class Phylogeny extends Component {
    constructor(props) {
        // remove old vis and tooltips
        d3v6.select('#tree_vis').remove();
        d3v6.select('#hitbars').remove();
        d3v6.select('#tooltip').remove();
        super(props);

        var treeData = this.props.data.tree;
        if( Object.keys(this.props.data).length === 3){
            this.state = {tree: treeData, extra: this.props.data.extraInfo, actualTree: null, counter: 0};
        }else{
            this.state = {tree: treeData, extra: [1], actualTree: null, counter: 0};
        }
        console.log(this.state);

        // own functions
        this.handleHomeClick = this.handleHomeClick.bind(this);
        this.handleHelpClick = this.handleHelpClick.bind(this);
        this.distTree = this.distTree.bind(this);
        this.d3Tree = this.d3Tree.bind(this);
        this.showAdditional = this.showAdditional.bind(this);


    }

    handleHomeClick(){
        this.props.changeComp('actual', 'menu');
    }

    handleHelpClick(){
        this.props.changeComp('actual', 'help');
    }


    distTree(){
        d3v6.selectAll('#tree_vis').remove(); // remove old vis & new empty svg
        d3v6.select('#tree').append('svg').attr('id', 'tree_vis');
        var libTree = d3.layout.phylotree()
          // create a tree layout object
          .svg(d3.select("#tree_vis"))
          .options({'align-tips': true,
                     'left-offset': 20,
                     'selectable': true});
        // render to this SVG element
            libTree(this.props.data.newick)
              // parse the Newick into a d3 hierarchy object with additional fields
              .layout();
            // layout and render the tree
            // for syntax highlighting
           // hljs.initHighlightingOnLoad();

        this.setState({actualTree: libTree});
        if(document.getElementById('clade_info')){
            showClades(this.state.extra[0], this.state.extra[1], libTree);
        }


    }

    d3Tree(){
        var tree = startTreevis(this.state.tree);
        //console.log(tree)
        if (tree !== 0){

            chart(tree, this.state.extra, null, null, true);
            if (this.state.counter !== 0 && document.getElementById('clade_info')){
                showClades(this.state.extra[0], this.state.extra[1], this.state.actualTree);
            }

        }else{
          d3v6.select('#tree').append('div').text('Found 0 hits. Return to the main page and try another phylogentic tree');
        }
    }

    showAdditional(event){
        //console.log(this.state.extra);
        showClades(this.state.extra[0], this.state.extra[1], this.state.actualTree);
    }

    render(){
        d3v6.select('#treeVis').style('border', '2px solid #5e66b4')
                                .style('border-radius', '5px');

        if (this.state.counter === 0){
            this.d3Tree();
            this.setState({counter: 1});
        }

        let additionalInformation;
        let labeling;
        if (this.state.extra.length === 2){
            additionalInformation = <input id='clade_info' type='checkbox' onChange={this.showAdditional}/>;
            labeling = <label>show additional clade information</label>;
        }



        return(
                    <div id="phylogeny">
                        <div id='phylogenyMenu'>
                            <h2>phylogenetic analysis</h2>
                            <li ><button id="distTree" onClick={this.distTree} >distance focused tree </button></li>
                            <li ><button id="d3Tree" onClick={this.d3Tree}>clade focused tree </button></li>
                            <br></br><br></br>
                           <button id="public_ready" title="Result will be a static tree. Reload the page to start a new interpolation of the basic tree" onClick={publicationReady}>publication ready</button>
                           <ExportTrees />
                           <br></br>
                           <br></br>
                           {additionalInformation}
                           {labeling}
                        </div>
                    </div>

        );
    }
}

export default Phylogeny;

