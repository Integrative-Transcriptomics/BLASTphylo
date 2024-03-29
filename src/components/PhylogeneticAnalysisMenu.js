// used libraries
import React, {Component, useState} from 'react';
import * as d3v6 from 'd3v6';
import * as d3 from 'd3';
import * as $ from "jquery";
import * as _ from "lodash";
import "../../node_modules/phylotree/src/main";
import {Nav, Accordion, Card, Tooltip, OverlayTrigger, Form, Popover} from 'react-bootstrap';
import {BiHelpCircle} from 'react-icons/bi';

// other style sheets
import '../../node_modules/phylotree/phylotree.css'

// own components and style sheets
import './phylogenyStyle.css';
import {ExportTreeImage, ExportCsvData} from './HandleDataExports.js'


// own visualisations
import {showClades, publicationReady, chart} from '../visualisations/phyloblast2.js';


class PhylogeneticAnalysisMenu extends Component{
    constructor(props){
        super(props);
        //console.log(props)
        // set state for visualisations
        var treeData = this.props.data.tree;
        this.state = { ownInfo: null, counter: 0};



        // own functions
        this.showAdditional = this.showAdditional.bind(this);
        this.handleUploadAdditional = this.handleUploadAdditional.bind(this);
        this.convertCSVtoJson = this.convertCSVtoJson.bind(this);
        this.showTooltip = this.showTooltip.bind(this);
        this.handleReturn = this.handleReturn.bind(this);
        this.handleTabChange = this.handleTabChange.bind(this);

    }

    /****componentDidMount(){
        if (this.props.phylogenyState === 'taxa'){
            console.log(this.props.phylogenyState)
            document.getElementById('phylogenyTab').addEventListener('click', this.handleTabChange);
        }else if (this.props.phylogenyState === 'seq'){
            document.getElementById('phylogenyUniqueTab').addEventListener('click', this.handleTabChange);
        }
    }***/

    handleTabChange(){
        console.log('handle Tab change')
        this.handleReturn('Tabchange');
    }

    componentDidUpdate(prevProps){
        if( (this.props.calculationMethod !== prevProps.calculationMethod) &&
        (document.getElementById('infoSelection'))){
            showClades(this.props.data.extraInfo[0], this.props.data.extraInfo[1], null, null);
        }

    }


    // show additional information for the taxa-based phylogeny
    showAdditional(event){
        showClades(this.props.data.extraInfo[0], this.props.data.extraInfo[1], this.state.ownInfo, null);
    }

    // convert the given additional information in a json object for better handling in d3
    // based on: https://stackoverflow.com/questions/27979002/convert-csv-data-into-json-format-using-javascript
    convertCSVtoJson(csv){
        var result = {};
        const self = this;
        csv.onload = function(){

            var lines = this.result.split('\n');
            const headers = lines[0].split(',');

            for (let i = 1; i < lines.length; i++) {
                if (!lines[i])
                    continue
                const obj = {};
                const currentline = lines[i].split(',');

                for (let j = 1; j < headers.length; j++) {
                    obj[headers[j]] = currentline[j];
                }
                result[currentline[0]] = obj;
            }
            self.setState({ownInfo: result});
        }
    }

    // handle the upload of the additional information
    handleUploadAdditional(event){
        const reader = new FileReader();
        reader.readAsText(event.target.files[0])
        this.convertCSVtoJson(reader);
        document.getElementById('uploadOwnAdditional').style.display = "none";
        document.getElementById('ownCheckbox').style.display = "block";
    }

    // hide a tooltip in 5 seconds
    showTooltip(event){
        setTimeout(function(){d3v6.select('#' + event.id).style("visibility","hidden");}, 5000);
    }

    // restore old tree state if publication ready was clicked
    handleReturn(event){
        var taxonomyLevel = ['life', 'domain', 'superkingdom', 'kingdom', 'clade', 'phylum', 'class', 'order', 'family', 'genus', 'species group','species', 'strain'];
        document.getElementById('returnButton').style.display = 'none';
        const rank = taxonomyLevel.indexOf(this.state.rankSelect);
        d3v6.select('#tree_vis').remove();
        const treeCopy = {...this.props.data.tree};
        console.log(treeCopy['size'])
        // update the visualizations
        chart({'size': treeCopy['size']}, this.props.data.extraInfo, rank, false, true, true);
        /***if(event === 'Tabchange'){
            chart(treeCopy, this.props.data.extraInfo, rank, false, true, true);
        }else{
            chart({'size': treeCopy['size']}, this.props.data.extraInfo, rank, false, true, true);
        }***/
        if(document.getElementById('infoSelection')){
            showClades(this.props.data.extraInfo[0], this.props.data.extraInfo[1], null, null);
        }

        // change visualisation settings
        document.getElementById('treeVis').style.height = '80vh';
        document.getElementById('public_ready_phylo').disabled = false;
    }

    render(){
        d3v6.select('#visualisation').style('border', '2px solid #5e66b4')
                                .style('border-radius', '5px');

        let additionalInformation;
        if (this.props.phylogenyState === 'taxa'){
            additionalInformation = additionalCladeInformation(this.showAdditional, this.handleUploadAdditional);
        }

        let newickFilename;

        if(this.props.phylogenyState === 'seq'){
            newickFilename = 'unique_sequence_based_phylogeny.txt';
        }else{
            newickFilename = 'taxa_based_phylogeny.txt';
        }

        return(
            <div id='phylogeny'>
                <div id='phylogenyMenu'>
                    <h2>                        </h2>
                    <Accordion defaultActiveKey='1'>
                        {additionalInformation}
                        <Card>
                            <Accordion.Toggle as={Card.Header} eventKey='1'>
                                download tree information
                            </Accordion.Toggle>
                            <Accordion.Collapse eventKey='1'>
                                <Card.Body>
                                    Different export options for the tree visualization. The <b>publication ready</b> option will remove the
                                    white spaces between the nodes of the cladogram and visualizes the full tree as static visualization.<br/>
                                    The <b>unlock publication ready</b> button unlocks the static visualization. <br/><br/>
                                    <Nav className='mr-auto'>
                                        <button id="public_ready_phylo" onClick={publicationReady}>publication ready</button>
                                        <ExportTreeImage />
                                        <ExportCsvData dataName='Newick string' filename={newickFilename} />
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


// Handle the display of  the additional clade information menu if the taxa-based phylogeny was calculated
function additionalCladeInformation(showAdditional, handleUploadAdditional){
    return(
        <Card>
            <Accordion.Toggle as={Card.Header} eventKey='0'>
                visualize additional clade information
            </Accordion.Toggle>
            <Accordion.Collapse eventKey='0'>
                <Card.Body>
                    Visualize additional information for each node as heat map on the right of the tree. The heat map will be included in
                    the exported visualization. <br/><br/>
                    <Form id='infoSelection'>
                        <div key={'inline-checkbox'} className='mb-3'>
                            <Form.Group inline id='ownCheckbox' style={{display: 'none'}}>
                                    <Form.Check inline type={'checkbox'} id='ownAdditional' label={'own information'} onChange={showAdditional} />
                            </Form.Group>
                            <Form.Check inline type={'checkbox'}  id='uniqueCheck' label={'unique seqs.'}  onChange={showAdditional}/>
                            <Form.Check inline type={'checkbox'}  id='phylumCheck' label={'phylum rank'} onChange={showAdditional}/>
                            <br /><br />
                            <Form inline id='uploadOwnAdditional'>
                                <Form.Label>upload additional information as csv</Form.Label>
                                <OverlayTrigger trigger='click' placement='right' overlay={
                                    <Popover id='popover-basic' >
                                        <Popover.Title as='h3'><strong>Additional information</strong></Popover.Title>
                                        <Popover.Content>Upload csv file with additional information about the taxa.
                                            Each column will treated as <b>boolean attribute</b>, whereby 0 indicates false and any value greater than 0 as true.
                                            The first line will be treated as header and the naming convention has to be consistent with the tree labels.</Popover.Content>
                                    </Popover>}>
                                    <BiHelpCircle style={{color: 'blue', 'margin': '0px 10px 0px 5px'}}/>
                                </OverlayTrigger>
                                <Form.File id='additional_file' name='additional_file'
                                onChange={handleUploadAdditional} accept='.csv' />
                            </Form>
                        </div>
                    </Form>
                </Card.Body>
            </Accordion.Collapse>
        </Card>
    );

}



export default PhylogeneticAnalysisMenu;