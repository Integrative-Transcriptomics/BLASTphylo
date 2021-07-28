// used libraries
import React, {Component} from 'react';
import * as d3v6 from 'd3v6';
import * as d3 from 'd3';
import * as $ from "jquery";
import * as _ from "lodash";
import "../../node_modules/phylotree/src/main";
import {ButtonToolbar, ButtonGroup, Button, DropdownButton, Dropdown, Tooltip, OverlayTrigger, Form} from 'react-bootstrap';
import {AiOutlineNodeCollapse, AiOutlineBarChart, AiOutlineArrowLeft} from 'react-icons/ai';

// own visualisations
import { hitBars, collapseTree, chart, startTreevis, showClades, stackBars} from '../visualisations/phyloblast2.js';


// constant for taxonomic ranks
var taxonomyLevel = ['life', 'domain', 'superkingdom', 'kingdom', 'clade', 'phylum', 'class', 'order', 'family', 'genus', 'species group','species', 'strain'];


// component for all direct tree interactions and frequently used functions
class TreeInteraction extends Component{
     constructor(props){
        super(props);
        console.log(props)
        this.state = {hitSelect: '-', rankSelect: 'class', counter: 0};

        // functions to handle the interactions and tooltips
        this.handleHits = this.handleHits.bind(this);
        this.handleRanks = this.handleRanks.bind(this);
        this.distTree = this.distTree.bind(this);
        this.d3Tree = this.d3Tree.bind(this);
    }


    // taxonomic analysis: tree interactions
    // handle hit barchart
    handleHits(event){
        this.setState({hitSelect: event});

        if(!(taxonomyLevel.includes(this.props.data.tree['value'][2]))){
            hitBars(event);
        }else{
            stackBars(event);
        }
    }

    // handle 'collapse to' interaction
    handleRanks(event){
        if(!(document.getElementById('public_ready_taxa').disabled)){
            document.getElementById('collapse_menu').disabled = false;
            this.setState({rankSelect: event});
            collapseTree(event);
            console.log(event);

            // no white space on the strain level --> disable publication ready
            if(event === 'strain'){
                document.getElementById('public_ready_taxa').disabled = true;
            }else{
                document.getElementById('public_ready_taxa').disabled = false;
            }
        }
    }

    // phylogenetic analysis: types of visualisations
    // visualise the distance-focused visualisation (phylogram)
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
              hide: false,
              reroot: false,
              transitions: true,
              "internal-names": true,
              "draw-size-bubbles": true
        });

        // render to this SVG element
            libTree(this.props.data.newick)
              .layout();

        // update additional information if nodes will be collapsed
        if(document.getElementById('infoSelection')){
            //console.log('run showClades')
            showClades(this.props.data.extraInfo[0], this.props.data.extraInfo[1], null, libTree);
            // Select the node that will be observed for mutations
            const targetNode = document.getElementById('tree_vis');

            // Options for the observer (which mutations to observe)
            const config = { attributes: true };
            const self = this;
            // Callback function to execute when mutations are observed
            const callback = function(mutationsList, observer) {
                // Use traditional 'for loops' for IE 11
                for(const mutation of mutationsList) {
                    if (mutation.type === 'attributes' && mutation.attributeName === 'height') {
                        showClades(self.props.data.extraInfo[0], self.props.data.extraInfo[1], null, null);
                    }
                }
            };

            // Create an observer instance linked to the callback function
            const observer = new MutationObserver(callback);

            // Start observing the target node for configured mutations
            observer.observe(targetNode, config);
        }

        // change settings for phylogram
        document.getElementById('treeVis').style.height = '80vh';
        document.getElementById('public_ready_phylo').disabled = true;

    }

    // visualise the clade-focused visualisation (cladogram)
    d3Tree(){
        d3v6.selectAll('#tree_vis').remove();
        const rank = this.state.rankSelect;
        const treeCopy = {...this.props.data.actualTree};
        //console.log(treeCopy)
        if(this.state.counter === 0){
            chart(treeCopy, this.props.data.extraInfo, null, false, true, false);
        }else{
            chart({'size': treeCopy['size']}, this.props.data.extraInfo, null, false, true, false);
        }
        if(document.getElementById('infoSelection')){
            showClades(this.props.data.extraInfo[0], this.props.data.extraInfo[1], null, null);
        }

        // change settings for cladogram
        document.getElementById('treeVis').style.height = '80vh'
        if(this.state.counter > 0){
            document.getElementById('public_ready_phylo').disabled = false;
        }
    }

    // hide the tooltip in 5 seconds
    showTooltip(event){
        setTimeout(function(){d3v6.select('#' + event.id).style("visibility","hidden");}, 5000);
    }

    componentDidMount(){
        if (this.props.calculationMethod.includes('phylo')){
            d3v6.select('#clade_vis').remove();
            const tree = startTreevis(this.props.data.actualTree);
            //console.log(tree);
            if (tree !== 0){
                this.d3Tree();
                this.setState({counter: 1})
            }else{
                d3v6.select('#tree').append('div').text('Found 0 hits. Return to the main page and try another phylogentic tree');
            }
        }

    }
    componentDidUpdate(prevProps){
        //console.log(this.props)
        //console.log(prevProps)
        if (this.props.calculationMethod.includes('phylo')){
            d3v6.select('#clade_vis').remove();
            d3v6.select('#hitbars').remove();
            d3v6.select('#tooltip').remove();

            const tree = startTreevis(this.props.data.actualTree);
            //console.log(tree);
            this.d3Tree();
        }
    }

    render(){

        if(this.props.calculationMethod === 'taxa'){ // tree interaction for taxonomic mapping
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
                    visualize node/subtree hits for actual state of the tree
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
                                <Dropdown.Item eventKey='-'>none</Dropdown.Item>
                                <Dropdown.Item eventKey='node hits'>node hits</Dropdown.Item>
                                <Dropdown.Item eventKey='subtree hits'>subtree hits</Dropdown.Item>
                            </DropdownButton>
                            </OverlayTrigger>
                            <Form id='selectedValues' inline>
                                <Form className='selectedValues' id='taxoRank' eventKey={this.state.rankSelect}><i>taxonomic Rank: {this.state.rankSelect}</i></Form>
                                <Form className='selectedValues' id='barChart' eventKey={this.state.hitSelect}><i>bar chart: {this.state.hitSelect}</i></Form>
                            </Form>
                        </ButtonGroup>
                    </ButtonToolbar>

                </div>
            );
        } else{ // tree interaction for phylogenetic analysis
            const renderDistanceTooltip = (props) => (
                    <Tooltip id="distance_tooltip" {... props}>
                         visualize sequence distances
                    </Tooltip>
            );
            const renderCladeTooltip = (props) => (
                    <Tooltip id="clade_tooltip" {... props}>
                         visualize groups in a straight tree
                    </Tooltip>
            );
            return(
                <div id='treeInteraction'>
                <ButtonToolbar aria-label='Toolbar with button groups'>
                    <OverlayTrigger placement='top' overlay={renderCladeTooltip} onEnter={this.showTooltip}>
                    <Button id='cladogram' eventKey='cladogramButton' onClick={this.d3Tree}>Cladogram</Button>
                    </OverlayTrigger>
                    <OverlayTrigger placement='top' overlay={renderDistanceTooltip} onEnter={this.showTooltip}>
                    <Button id='phylogram' eventKey='phylogramButton' onClick={this.distTree}>Phylogram</Button>
                    </OverlayTrigger>
                </ButtonToolbar>
                </div>
            );
        }
    }
}



export default TreeInteraction;

