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
import { hitBars, collapseTree, chart, startTreevis, showClades, stackBars, publicationReady} from '../visualisations/phyloblast2.js';


// constant for taxonomic ranks
var taxonomyLevel = ['life', 'domain', 'superkingdom', 'kingdom', 'clade', 'phylum', 'class', 'order', 'family', 'genus', 'species group','species', 'strain'];


// component for all direct tree interactions and frequently used functions
class TreeInteraction extends Component{
     constructor(props){
        super(props);
        console.log('Tree interaction props', props)
        this.state = {hitSelect: '-',
                      rankSelect: 'class',
                      counter: 0,
                      publicationReady: false,
                      noHitTree: false};

        // functions to handle the interactions and tooltips
        this.handleHits = this.handleHits.bind(this);
        this.handleRanks = this.handleRanks.bind(this);
        this.distTree = this.distTree.bind(this);
        this.d3Tree = this.d3Tree.bind(this);
        this.handlePublicationReady = this.handlePublicationReady.bind(this);
        this.handleNoHit = this.handleNoHit.bind(this);
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
    handleReturn(event){


        var taxonomyLevel = ['life', 'domain', 'superkingdom', 'kingdom', 'clade', 'phylum', 'class', 'order', 'family', 'genus', 'species group','species', 'strain'];
        //document.getElementById('returnButton').style.display = 'none';
        const rank = taxonomyLevel.indexOf(this.state.rankSelect);
        console.log('handleReturn', rank);

        d3v6.select('#tree_vis').remove();
        const treeCopy = {...this.props.data.actualTree};

        if (this.state.noHitTree==true){
            chart(this.props.data.actualTree, null, rank, true, false, false);
        }
        // update visualizations
        // add if else depending on taxonomy or phylogeny
        if(this.props.calculationMethod=='taxa'){
            chart({'size': treeCopy['size']}, null, rank, false, true, true);
        }else{
            chart({'size': treeCopy['size']}, this.props.data.extraInfo, rank, false, true, true);}

        if(!(taxonomyLevel.includes(treeCopy['value'][2]))){
                hitBars(this.state.hitSelect);
        }else{
                stackBars(this.state.hitSelect);
        }

        // change visualisation settings
        document.getElementById('treeVis').style.height = '80vh';
        //document.getElementById('collapse_menu').disabled = false;
    }
    // handle 'collapse to' interaction
    handleRanks(event){
        if(!this.state.publicationReady){
            document.getElementById('collapse_menu').disabled = false;
            this.setState({rankSelect: event});
            console.log(this.state.rankSelect);
            collapseTree(event);
            console.log(event);

            // no white space on the strain level --> disable publication ready
//            if(event === 'strain'){
//                document.getElementById('public_ready_taxa').disabled = true;
//            }else{
//                document.getElementById('public_ready_taxa').disabled = false;
//            }
        }
    }

    handleNoHit(){
        if (!this.state.noHitTree){
            d3v6.select('#tree_vis').remove();
            console.log(this.props.data);
            chart(this.props.data.noHitTree, false, this.state.rankSelect, true, false, false);
            this.setState({noHitTree: !this.state.noHitTree});
        }
        else {
            this.handleReturn();
            this.setState({noHitTree: !this.state.noHitTree});
        }
    }

    handlePublicationReady(){
        console.log('print publicationReady state', this.state.publicationReady);
        console.log(document.getElementById("publicationReady"));
        if (!this.state.publicationReady){
            publicationReady();
            this.setState({publicationReady: !this.state.publicationReady});
        }
        else {
            console.log(this.state.publicationReady);
            console.log(document.getElementById("publicationReady"));
            console.log("return from publication ready")
            this.handleReturn();
            this.setState({publicationReady: !this.state.publicationReady});
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
        document.getElementById('publicationReadyPhylo').disabled = true;

    }

    // visualise the clade-focused visualisation (cladogram)
    d3Tree(){
        d3v6.selectAll('#tree_vis').remove();
        const rank = this.state.rankSelect;
        document.getElementById('publicationReadyPhylo').disabled = false;

        //chart(treeCopy, this.props.data.extraInfo, null, true, true, true);
        if(this.state.counter === 0){
            const treeCopy = {...this.props.data.actualTree};
            console.log(treeCopy['size'])
            chart(treeCopy, this.props.data.extraInfo, null, true, true, false);
        }else{
            const treeCopy = {...this.props.data.actualTree};
            var sizeOriginal = this.props.data.tree['size'][0];
            treeCopy['size'][0] = sizeOriginal;
            console.log(treeCopy['size'])
            chart(treeCopy, this.props.data.extraInfo, null, false, true, true);
        }
        if(document.getElementById('infoSelection')){
            showClades(this.props.data.extraInfo[0], this.props.data.extraInfo[1], null, null);
        }

        // change settings for cladogram
        document.getElementById('treeVis').style.height = '80vh'
        if(document.getElementById('public_ready_phylo')){
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
        if ((this.props.calculationMethod !== prevProps.calculationMethod) &&
        (this.props.calculationMethod.includes('phylo'))){
            d3v6.select('#clade_vis').remove();
            d3v6.select('#hitbars').remove();
            d3v6.select('#tooltip').remove();

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
                    collapse tree to selected taxonomic rank
                </Tooltip>
            );
            const renderBarchartTooltip = (props) => (
                <Tooltip id="barchart_tooltip" {... props}>
                    select counts displayed in barchart next to the tree
                </Tooltip>
            );

            return(
                <div id='treeInteraction'>
                    <ButtonToolbar aria-label='Toolbar with button groups'>
                        <ButtonGroup className='mr-2' aria-label='First group'>
                            <OverlayTrigger placement='top'  overlay={renderCollapseTooltip} onEnter={this.showTooltip}>
                            <DropdownButton onSelect={this.handleRanks} as={ButtonGroup} eventKey={this.state.rankSelect} title={<span> <AiOutlineNodeCollapse size={30}/> {this.state.rankSelect} </span>} id='collapse_menu'>
                                {taxonomyLevel.map(MakeItem)}
                            </DropdownButton>
                            </OverlayTrigger>
                            <OverlayTrigger placement='top' overlay={renderBarchartTooltip} onEnter={this.showTooltip}>
                            <DropdownButton onSelect={this.handleHits} eventKey={this.state.hitSelect}
                            as={ButtonGroup} title={<span> <AiOutlineBarChart size={30}/> {this.state.hitSelect} </span>} id='tree_menu'>
                                <Dropdown.Item eventKey='-'>none</Dropdown.Item>
                                <Dropdown.Item eventKey='node hits'>node hits</Dropdown.Item>
                                <Dropdown.Item eventKey='subtree hits'>subtree hits</Dropdown.Item>
                                <Dropdown.Item eventKey='covered taxa'>covered taxa</Dropdown.Item>
                            </DropdownButton>
                            </OverlayTrigger>
                        </ButtonGroup>
                        <ButtonGroup className='mr-2' aria-label='Second group'>
                            <div class="d-flex align-items-center flex-row" as={ButtonGroup}>
                                        <Form>
                                            <Form.Check
                                            type="switch"
                                            id="publicationReady"
                                            label="PUBLICATION READY"
                                            onChange={this.handlePublicationReady}
                                            />
                                        </Form>
                            </div>
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
                <ButtonToolbar id="phylogenybar" aria-label='Toolbar with button groups'>
                        <Form inline>
                            <Form.Check inline type={'radio'} name='treeVisualizationOptions' id='cladogram' label={'CLADOGRAM'} onChange={this.d3Tree} defaultChecked/>
                            <Form.Check inline type={'radio'} name='treeVisualizationOptions' id='phylogram' label={'PHYLOGRAM'} onChange={this.distTree} />
                        </Form>
                        <ButtonGroup className='mr-2' aria-label='Second group'>
                            <div class="d-flex align-items-center flex-row" as={ButtonGroup}>
                                        <Form>
                                            <Form.Check
                                            type="switch"
                                            id="publicationReadyPhylo"
                                            label="PUBLICATION READY"
                                            onChange={this.handlePublicationReady}
                                            />
                                        </Form>
                            </div>
                        </ButtonGroup>
                </ButtonToolbar>
                </div>
            );
        }
    }
}



export default TreeInteraction;

/*<ButtonGroup className='mr-2' aria-label='Third group'>
                            <div class="d-flex align-items-center flex-row" as={ButtonGroup}>
                                        <Form>
                                            <Form.Check
                                            type="switch"
                                            id="noHitToggle"
                                            label="SHOW TAXA WITHOUT HITS"
                                            onChange={this.handleNoHit}
                                            />
                                        </Form>
                            </div>
                        </ButtonGroup>*/