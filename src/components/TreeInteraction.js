// used libraries
import React, {Component} from 'react';
import * as d3v6 from 'd3v6';
import * as d3 from 'd3';
import * as $ from "jquery";
import * as _ from "lodash";
import "../../node_modules/phylotree/src/main";
import {ButtonToolbar, ButtonGroup, Button, DropdownButton, Dropdown, Tooltip, OverlayTrigger} from 'react-bootstrap';
import {AiOutlineNodeCollapse, AiOutlineBarChart} from 'react-icons/ai';
//import {FaBezierCurve} from 'react-icons/fa';
//import {VscListTree} from 'react-icons/vsc';

// own visualisations
import { hitBars, collapseTree, chart, startTreevis, showClades} from '../visualisations/phyloblast2.js';


// component for all direct tree interactions and frequently used functions
class TreeInteraction extends Component{
     constructor(props){
        super(props);

        if(this.props.calculationMethod === 'taxa'){
            this.state = {hitSelect: "2",
                      rankSelect: 'class'};
        }else{
            // remove old vis and tooltips
            d3v6.select('#tree_vis').remove();
            d3v6.select('#hitbars').remove();
            d3v6.select('#tooltip').remove();

            var treeData = this.props.data.tree;
            if( Object.keys(this.props.data).length === 3){ // taxa-based phylogeny
                this.state = {tree: treeData, extra: this.props.data.extraInfo, actualTree: treeData, counter: 0};
            }else{ // sequence based phylogeny
                this.state = {tree: treeData, extra: [1], actualTree: treeData, counter: 0};
            }
        }

        // functions to handle the interactions and tooltips
        this.handleHits = this.handleHits.bind(this);
        this.handleRanks = this.handleRanks.bind(this);
        this.distTree = this.distTree.bind(this);
        this.d3Tree = this.d3Tree.bind(this);
        this.showTooltip = this.showTooltip.bind(this);
    }

    // taxonomic analysis: tree interactions
    // handle hit barchart
    handleHits(event){
        this.setState({hitSelect: event});
        if(event !== '2'){
            hitBars(event);
        }else{
            d3v6.select('#hitbars')
                .remove();
        }
    }

    // handle 'collapse to' interaction
    handleRanks(event){
        this.setState({rankSelect: event});
        var copy = {...this.props.phyloData};
        collapseTree(copy, event);
    }

    // phylogenetic analysis: types of visualisations
    // visualise the distance-focused visualisation
    distTree(){

        // set scroll bar to the top of the visualisation
        document.getElementById("treeVis").scrollTop = 0;
        d3v6.selectAll('#tree_vis').remove(); // remove old vis & new empty svg
        d3v6.select('#tree').append('svg').attr('id', 'tree_vis');
        var libTree = d3.layout.phylotree()
            .svg(d3.select("#tree_vis"))
            .options({
              brush: false,
              'align-tips': true,
              zoom: false,
              selectable: false,
              collapsible: true,
              hide: true,
              reroot: false,
              transitions: true,
              "internal-names": true,
              "draw-size-bubbles": true,
        });
        // render to this SVG element
            libTree(this.props.data.newick)
              .layout();

        if(document.getElementById('infoSelection')){
            showClades(this.state.extra[0], this.state.extra[1], libTree);
        }
    }

    // visualise the clade-focused visualisation
    d3Tree(){
        const treeVis = {...this.state.actualTree};
        chart(treeVis, this.state.extra, null, null, true);
        if(document.getElementById('infoSelection')){
            showClades(this.state.extra[0], this.state.extra[1], null, null);
        }
    }

    // hide the tooltip in 5 seconds
    showTooltip(event){
        setTimeout(function(){d3v6.select('#' + event.id).style("visibility","hidden");}, 5000);
    }


    render(){
        if(this.props.calculationMethod === 'taxa'){ // tree interaction for taxonomic mapping
            var taxonomyLevel = ['life', 'domain', 'superkingdom', 'kingdom', 'clade', 'phylum', 'class', 'order', 'family', 'genus', 'species group','species', 'strain'];
            var MakeItem = function(X){
                    return <Dropdown.Item eventKey={X}>{X}</Dropdown.Item>;
            };

            const renderCollapseTooltip = (props) => (
                    <Tooltip id="collapse_tooltip" {... props}>
                         collapse all nodes below selected taxonomic rank
                    </Tooltip>
            );
            const renderBarchartTooltip = (props) => (
                    <Tooltip id="barchart_tooltip" {... props}>
                         visualise node/subtree hits for actual state of the tree
                    </Tooltip>
            );


            return(
                <div id='treeInteraction'>
                    <ButtonToolbar aria-label='Toolbar with button groups'>
                        <ButtonGroup className='mr-2' aria-label='First group'>
                            <OverlayTrigger placement='top'  overlay={renderCollapseTooltip} onEnter={this.showTooltip}>
                            <DropdownButton onSelect={this.handleRanks} as={ButtonGroup} title={<AiOutlineNodeCollapse size={25}/>} id='collapse_menu'>
                                {taxonomyLevel.map(MakeItem)}
                            </DropdownButton>
                            </OverlayTrigger>
                            <OverlayTrigger placement='top' overlay={renderBarchartTooltip} onEnter={this.showTooltip}>
                            <DropdownButton onSelect={this.handleHits} eventKey={this.state.hitSelect}
                            as={ButtonGroup} title={<AiOutlineBarChart size={25}/>} id='tree_menu'>
                                <Dropdown.Item eventKey='2'>none</Dropdown.Item>
                                <Dropdown.Item eventKey='0'>node hits</Dropdown.Item>
                                <Dropdown.Item eventKey='1'>subtree hits</Dropdown.Item>
                            </DropdownButton>
                            </OverlayTrigger>
                        </ButtonGroup>
                    </ButtonToolbar>
                </div>
            );
        } else{ // tree interaction for phylogenetic analysis
            if (this.state.counter === 0){
                const tree = startTreevis(this.state.tree);
                this.setState({actualTree: tree});
                if (tree !== 0){

                    this.d3Tree();
                }else{
                  d3v6.select('#tree').append('div').text('Found 0 hits. Return to the main page and try another phylogentic tree');
                }
                this.setState({counter: 1});
            }

            const renderDistanceTooltip = (props) => (
                    <Tooltip id="distance_tooltip" {... props}>
                         visualise sequence distances
                    </Tooltip>
            );
            const renderCladeTooltip = (props) => (
                    <Tooltip id="clade_tooltip" {... props}>
                         visualise groups in a straight tree
                    </Tooltip>
            );
            return(
                <div id='treeInteraction'>
                <ButtonToolbar aria-label='Toolbar with button groups'>
                    <OverlayTrigger placement='top' overlay={renderCladeTooltip} onEnter={this.showTooltip}>
                    <Button onClick={this.d3Tree}>Cladogram</Button>
                    </OverlayTrigger>
                    <OverlayTrigger placement='top' overlay={renderDistanceTooltip} onEnter={this.showTooltip}>
                    <Button onClick={this.distTree}>Phylogram</Button>
                    </OverlayTrigger>
                </ButtonToolbar>
                </div>
            );
        }
    }
}


export default TreeInteraction;

// Icon distTree: {<VscListTree size={25}/>} Icon d3Tree: {<FaBezierCurve size={25}/>}