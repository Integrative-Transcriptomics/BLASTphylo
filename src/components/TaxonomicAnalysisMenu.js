// used libraries
import React, {Component} from 'react';
import * as d3v6 from 'd3v6';
import {Nav, Accordion, Card, Tooltip, OverlayTrigger, Button, Form} from 'react-bootstrap';
import {BsBoxArrowUpRight} from "react-icons/bs";
import {AiOutlineArrowLeft} from 'react-icons/ai';

// own components and style sheets
import './phyloblastStyle.css';
import {ExportTrees} from '../App.js';


// own visualisations
import {chart, startTreevis, publicationReady, hitBars, stackBars} from '../visualisations/phyloblast2.js';


class TaxonomicAnalysisMenu extends Component{
     constructor(props){
        super(props);

        // generation of the tree visualisation
        var treeData = startTreevis(this.props.phyloData);
        this.state = {hitSelect: "2",
                      rankSelect: 'class',
                      tree: treeData};

        if (treeData !== 0){
            d3v6.select('#visualisation').style('border', '2px solid #69a2c9')
                                .style('border-radius', '5px');
            document.getElementById('treeVis').style.height = '80vh';
            chart(treeData, null, 6, true, true);
        }else{ // if tree was empty show a message to make clear that the calculation worked but no hits were found
            d3v6.select('#phyloblastAlert').remove();
            d3v6.select('#tree').append('div').attr('id', 'phyloblastAlert')
                    .text('Found 0 hits. Return to the main page and try another taxonomy');
        }

        // functions
        this.showTooltip = this.showTooltip.bind(this);
        this.handleReturn = this.handleReturn.bind(this);

    }

    // hide any tooltip in 5 seconds
    showTooltip(event){
        setTimeout(function(){d3v6.select('#' + event.id).style("visibility","hidden");}, 5000);
    }

    // restore old tree state if publication ready was clicked
    handleReturn(event){
        var taxonomyLevel = ['life', 'domain', 'superkingdom', 'kingdom', 'clade', 'phylum', 'class', 'order', 'family', 'genus', 'species group','species', 'strain'];
        document.getElementById('returnButton').style.display = 'none';
        const rank = taxonomyLevel.indexOf(document.getElementById('taxoRank').attributes.eventkey.value);
        d3v6.select('#tree_vis').remove();
        const treeCopy = {...this.state.tree};

        // update visualizations
        chart({'size': treeCopy['size']}, null, rank, false, true, true);
        if(!(taxonomyLevel.includes(treeCopy['value'][2]))){
                hitBars(document.getElementById('barChart').attributes.eventkey.value);
        }else{
                stackBars(document.getElementById('barChart').attributes.eventkey.value);
        }

        // change visualisation settings
        document.getElementById('treeVis').style.height = '80vh';
        document.getElementById('public_ready_taxa').disabled = false;
        document.getElementById('collapse_menu').disabled = false;
    }


    render(){


        return(
            <div id='phyloblast'>
                <div id='phyloblastMenu'>
                    <h2>taxonomic analysis</h2>
                    <Accordion defaultActiveKey='2'>
                        <Card>
                            <Accordion.Toggle as={Card.Header} eventKey='0'>
                                phylogeny calculation
                            </Accordion.Toggle>
                            <Accordion.Collapse eventKey='0'>
                                <Card.Body>
                                    <Form inline>
                                        Links will start the phylogeny calculation in a new tab. A detailed description of the difference between the phylogenies
                                        can be found in the <Nav.Link style={{margin: '0px'}} href='/help' target='_blank'>help page</Nav.Link>
                                        <br/><br/><br/>
                                     </Form>
                                      <Nav className="flex-column">
                                        <Nav.Link href="/phylogeny" target='_blank'>taxa-based phylogeny <BsBoxArrowUpRight size={20}/></Nav.Link>
                                        <Nav.Link href="/phylogenyUnique" target='_blank'>unique sequence-based phylogeny <BsBoxArrowUpRight size={20} /></Nav.Link>
                                      </Nav>
                                </Card.Body>
                            </Accordion.Collapse>
                        </Card>
                        <Card>
                            <Accordion.Toggle as={Card.Header} eventKey='1'>
                                export tree visualization
                            </Accordion.Toggle>
                            <Accordion.Collapse eventKey='1'>
                                <Card.Body>
                                    Different export options for the taxonomic tree including the bar chart. The <b>publication ready</b> option will remove the
                                    white spaces between the nodes and visualizes the full tree as static visualization. <br/>
                                    The <b>unlock publication ready</b> button unlocks the static visualization. <br/><br/>
                                    <Nav className='mr-auto'>
                                        <button id="public_ready_taxa"  onClick={publicationReady}>publication ready</button>
                                        <ExportTrees />
                                        <button eventKey='returnButton' id='returnButton' onClick={this.handleReturn} style={{display: 'none'}}>unlock publication ready</button>
                                    </Nav>
                                </Card.Body>
                            </Accordion.Collapse>
                        </Card>
                    </Accordion>
                    <br></br>
                </div>
            </div>

        );
    }
}

export default TaxonomicAnalysisMenu;