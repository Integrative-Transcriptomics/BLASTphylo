// used libraries
import React, {Component} from 'react';
import * as d3v6 from 'd3v6';
import {CSVLink} from 'react-csv';
import axios from 'axios';
import {Nav, Accordion, Card, Tooltip, OverlayTrigger, Button, Form} from 'react-bootstrap';
import {BsBoxArrowUpRight} from "react-icons/bs";
import {AiOutlineArrowLeft} from 'react-icons/ai';

// own components and style sheets
import './phyloblastStyle.css';
import {ExportTreeImage, ExportCsvData, exportSVG, exportJPEG} from './HandleDataExports.js';


// own visualisations
import {chart, startTreevis, publicationReady, hitBars, stackBars} from '../visualisations/phyloblast2.js';


class TaxonomicAnalysisMenu extends Component{
     constructor(props){
        super(props);
        console.log(props)
        // generation of the tree visualisation
        //var treeData = startTreevis(this.props.data.tree, this.props.queries);
        //console.log(treeData['size'])
        this.state = {hitSelect: "-",
                      rankSelect: 'class',
                      treeExportFormat: 'svgSelect',
                      dataToDownload: ''};

        // functions
        this.showTooltip = this.showTooltip.bind(this);
        this.handleReturn = this.handleReturn.bind(this);
        this.exportTree = this.exportTree.bind(this);
        this.handleTreeExportFormat = this.handleTreeExportFormat.bind(this);
        this.exportTaxonomicMapping = this.exportTaxonomicMapping.bind(this);
    }


    componentDidMount(){
        d3v6.select('#tree_vis').remove();
        d3v6.select('#clade_vis').remove();
        if (this.props.data.tree !== 0){
            d3v6.select('#visualisation').style('border', '2px solid #fcb42d')
                                .style('border-radius', '5px');
            document.getElementById('treeVis').style.height = '80vh';

            chart(startTreevis(this.props.data.tree, this.props.queries), null, 6, true, true, true);

        }else{ // if tree was empty show a message to make clear that the calculation worked but no hits were found
            d3v6.select('#phyloblastAlert').remove();
            d3v6.select('#tree').append('div').attr('id', 'phyloblastAlert')
                    .text('Found 0 hits. Return to the main page and try another taxonomy');
        }
    }


    // hide any tooltip in 5 seconds
    showTooltip(event){
        setTimeout(function(){d3v6.select('#' + event.id).style("visibility","hidden");}, 5000);
    }

    // restore old tree state if publication ready was clicked
    handleReturn(event){
        var taxonomyLevel = ['life', 'domain', 'superkingdom', 'kingdom', 'clade', 'phylum', 'class', 'order', 'family', 'genus', 'species group','species', 'strain'];
        //document.getElementById('returnButton').style.display = 'none';
        const rank = this.state.rankSelect;
        d3v6.select('#tree_vis').remove();
        const treeCopy = {...this.props.data.actualTree};
        console.log(treeCopy['size'])

        // update visualizations
        chart({'size': treeCopy['size']}, null, rank, false, true, true);


        if(!(taxonomyLevel.includes(treeCopy['value'][2]))){
                hitBars(this.state.hitSelect);
        }else{
                stackBars(this.state.hitSelect);
        }

        // change visualisation settings
        document.getElementById('treeVis').style.height = '80vh';
        document.getElementById('collapse_menu').disabled = false;
    }

    exportTree(){
        if (this.state.treeExportFormat === "svgSelect") {
            exportSVG();
        }
        if (this.state.treeExportFormat === "jpegSelect") {
            exportJPEG();
        }
    }

    handleTreeExportFormat(event){
        console.log(event.target.id);
        this.setState({treeExportFormat: event.target.id});
    }

    exportTaxonomicMapping(){
        var data_to_download = 'Hello';
        this.setState({ dataToDownload: data_to_download }, () => {
         // click the CSVLink component to trigger the CSV download
        this.csvLink.link.click()
      })
    }

    exportNewickString(){
    }


    render(){

        return(
            <div id='phyloblast'>
                <div id='phyloblastMenu'>
                    <Accordion>
                        <Card>
                            <Accordion.Toggle as={Card.Header} eventKey='0'>
                                DOWNLOAD
                            </Accordion.Toggle>
                            <Accordion.Collapse eventKey='0'>
                                <Card.Body>
                                    <div class="d-flex align-items-center flex-row">
                                        <button class="btn btn-primary" onClick={this.exportTree}>EXPORT VIEW</button>

                                        <Form inline>
                                            <Form.Check inline type={'radio'} name='treeExportOptions' id='svgSelect' label={'SVG'} onChange={this.handleTreeExportFormat} defaultChecked/>
                                            <Form.Check inline type={'radio'} name='treeExportOptions' id='jpegSelect' label={'JPEG'} onChange={this.handleTreeExportFormat} />
                                        </Form>

                                    </div>

                                    <ExportCsvData dataName='taxonomic Mapping' filename='taxonomic_mapping.csv' />
                                    <ExportCsvData dataName='Newick string' filename='newick_taxonomic_mapping.txt' />

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
                                    /*<Nav className='mr-auto'>
                                        <button id="public_ready_taxa"  onClick={publicationReady}>publication ready</button>
                                        <ExportTreeImage />
                                        <ExportCsvData dataName='taxonomic Mapping' filename='taxonomic_mapping.csv' />
                                        <ExportCsvData dataName='Newick string' filename='newick_taxonomic_mapping.txt' />
                                        <button eventKey='returnButton' id='returnButton' onClick={this.handleReturn} style={{display: 'none'}}>unlock publication ready</button>
                                    </Nav>*/
/*<Accordion>
                <Accordion.Item eventKey='0'>
                    <Accordion.Header>
                        DOWNLOAD
                    </Accordion.Header>
                    <Accordion.Body>
                        <Nav className='mr-auto'>
                            <button id="public_ready_taxa"  onClick={publicationReady}>publication ready</button>
                            <ExportTreeImage />
                            <ExportCsvData dataName='taxonomic Mapping' filename='taxonomic_mapping.csv' />
                            <ExportCsvData dataName='Newick string' filename='newick_taxonomic_mapping.txt' />
                            <button eventKey='returnButton' id='returnButton' onClick={this.handleReturn} style={{display: 'none'}}>unlock publication ready</button>
                        </Nav>
                    </Accordion.Body>
                </Accordion.Item>
            </Accordion>*/