// used libraries
import React, {Component} from 'react';
import * as d3v6 from 'd3v6';
import {Nav, Accordion, Card, Tooltip, OverlayTrigger, Button} from 'react-bootstrap';
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
            chart(treeData, null, 6, 6, true);
        }else{
            d3v6.select('#phyloblastAlert').remove();
            d3v6.select('#tree').append('div').attr('id', 'phyloblastAlert')
                    .text('Found 0 hits. Return to the main page and try another taxonomy');
        }

        this.showTooltip = this.showTooltip.bind(this);
        this.handleReturn = this.handleReturn.bind(this);

    }

    // hide the tooltip in 5 seconds
    showTooltip(event){
        setTimeout(function(){d3v6.select('#' + event.id).style("visibility","hidden");}, 5000);
    }

    // restore old tree state if publication ready was clicked
    handleReturn(event){
        var taxonomyLevel = ['life', 'domain', 'superkingdom', 'kingdom', 'clade', 'phylum', 'class', 'order', 'family', 'genus', 'species group','species', 'strain'];
        document.getElementById('returnButton').style.display = 'none';
        // old version const rank = taxonomyLevel.indexOf(this.state.rankSelect);
        const rank = taxonomyLevel.indexOf(document.getElementById('taxoRank').attributes.eventkey.value);
        console.log(['Return button', this.state.rankSelect])
        d3v6.select('#tree_vis').remove();
        const treeCopy = {...this.state.tree};
        console.log(treeCopy)
        chart({'size': treeCopy['size']}, null, rank, rank, true, true);
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
        const renderReturnTooltip = (props) => (
            <Tooltip id="return_tooltip" {... props}>
                return to dynamic visualisation
            </Tooltip>
        );

        const returnButton = () => (
                                <OverlayTrigger placement='top' overlay={renderReturnTooltip} onEnter={this.showTooltip}>
                                <button eventKey='returnButton' id='returnButton' onClick={this.handleReturn} style={{display: 'none'}}>unlock publication ready</button>
                    </OverlayTrigger>
        );

        const renderPublicReadyTooltip = (props) => (
                <Tooltip id='public_ready_tooltip' {... props}>
                   Result will be a static tree.
                </Tooltip>
        );

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
                                        <button id="public_ready_taxa"  onClick={publicationReady}>publication ready</button>
                                        </OverlayTrigger>
                                        <ExportTrees />
                                        <OverlayTrigger placement='top' overlay={renderReturnTooltip} onEnter={this.showTooltip}>
                                        <button eventKey='returnButton' id='returnButton' onClick={this.handleReturn} style={{display: 'none'}}>unlock publication ready</button>
                                        </OverlayTrigger>
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