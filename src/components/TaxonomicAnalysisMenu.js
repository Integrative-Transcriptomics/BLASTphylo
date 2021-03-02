// used libraries
import React, {Component} from 'react';
import * as d3v6 from 'd3v6';
import {Nav, Accordion, Card, Tooltip, OverlayTrigger} from 'react-bootstrap';
import {BsBoxArrowUpRight} from "react-icons/bs";

// own components and style sheets
import './phyloblastStyle.css';
import {ExportTrees} from '../App.js';


// own visualisations
import {chart, startTreevis, publicationReady} from '../visualisations/phyloblast2.js';


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
            chart(treeData, null, 6, 6, true);
        }else{
            d3v6.select('#phyloblastAlert').remove();
            d3v6.select('#tree').append('div').attr('id', 'phyloblastAlert')
                    .text('Found 0 hits. Return to the main page and try another taxonomy');
        }

        this.handlePublicationReady = this.handlePublicationReady.bind(this);

    }

    handlePublicationReady(){
        var treeCopy = {...this.state.tree};
        publicationReady(treeCopy);
    }

    render(){

            const renderPublicReadyTooltip = (props) => (
                <Tooltip id='public_ready_tooltip' {... props}>
                   Result will be a static tree.
                </Tooltip>
            )

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
                                      <Nav className="flex-column">
                                        <Nav.Link href="/phylogeny" target='_blank'>taxa sequence based phylogeny <BsBoxArrowUpRight size={20}/></Nav.Link>
                                        <Nav.Link href="/phylogenyUnique" target='_blank'>unique sequence based phylogeny <BsBoxArrowUpRight size={20} /></Nav.Link>
                                      </Nav>
                                </Card.Body>
                            </Accordion.Collapse>
                        </Card>
                        <Card>
                            <Accordion.Toggle as={Card.Header} eventKey='1'>
                                export tree visualisation
                            </Accordion.Toggle>
                            <Accordion.Collapse eventKey='1'>
                                <Card.Body>
                                    <Nav className='mr-auto'>
                                        <OverlayTrigger placement='bottom' delay={{show:150, hide:50}} overlay={renderPublicReadyTooltip}>
                                        <button id="public_ready"  onClick={this.handlePublicationReady}>publication ready</button>
                                        </OverlayTrigger>
                                        <ExportTrees />
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