// used libraries
import React, {Component} from 'react';
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
import {ExportTrees} from '../App.js'


// own visualisations
import {showClades, publicationReady, chart} from '../visualisations/phyloblast2.js';


class PhylogeneticAnalysisMenu extends Component{
    constructor(props){
        super(props);


        // set state for visualisations
        var treeData = this.props.data.tree;
        if( Object.keys(this.props.data).length === 3){
            this.state = {tree: treeData, extra: this.props.data.extraInfo, ownInfo: null, actualTree: treeData, counter: 0};
        }else{
            this.state = {tree: treeData, extra: [1], actualTree: treeData, counter: 0};
        }

        // own functions
        this.showAdditional = this.showAdditional.bind(this);
        this.handleUploadAdditional = this.handleUploadAdditional.bind(this);
        this.convertCSVtoJson = this.convertCSVtoJson.bind(this);
        this.showTooltip = this.showTooltip.bind(this);
        this.handleReturn = this.handleReturn.bind(this);

    }

    // show additional information for the taxa-based phylogeny
    showAdditional(event){
        showClades(this.state.extra[0], this.state.extra[1], this.state.ownInfo, null);
    }

    // based on: https://stackoverflow.com/questions/27979002/convert-csv-data-into-json-format-using-javascript
    convertCSVtoJson(csv){
        var result = {};
        const self = this;
        csv.onload = function(){
            //console.log(this.result);

            var lines = this.result.split('\n');
            //console.log(lines)
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
            //console.log(result)
            self.setState({ownInfo: result});
        }
    }

    handleUploadAdditional(event){
        const reader = new FileReader();
        reader.readAsText(event.target.files[0])
        this.convertCSVtoJson(reader);
        document.getElementById('uploadOwnAdditional').style.display = "none";
        document.getElementById('ownCheckbox').style.display = "block";
    }

        // hide the tooltip in 5 seconds
    showTooltip(event){
        setTimeout(function(){d3v6.select('#' + event.id).style("visibility","hidden");}, 5000);
    }

    // restore old tree state if publication ready was clicked
    handleReturn(event){
        var taxonomyLevel = ['life', 'domain', 'superkingdom', 'kingdom', 'clade', 'phylum', 'class', 'order', 'family', 'genus', 'species group','species', 'strain'];
        document.getElementById('returnButton').style.display = 'none';
        const rank = taxonomyLevel.indexOf(this.state.rankSelect);
        console.log(['Return button', this.state.rankSelect])
        d3v6.select('#tree_vis').remove();
        const treeCopy = {...this.state.tree};
        console.log(treeCopy)
        chart({'size': treeCopy['size']}, this.state.extra, rank, false, true, true);
        if(document.getElementById('infoSelection')){
            showClades(this.state.extra[0], this.state.extra[1], null, null);
        }

        // change visualisation settings
        document.getElementById('treeVis').style.height = '80vh';
        document.getElementById('public_ready_phylo').disabled = false;
    }

    render(){
        d3v6.select('#visualisation').style('border', '2px solid #5e66b4')
                                .style('border-radius', '5px');

        let additionalInformation;
        if (this.state.extra.length === 2){
            additionalInformation = additionalCladeInformation(this.showAdditional, this.handleUploadAdditional);
        }

        let phylogenyType;
        if(window.location.href.includes('Unique')){
            phylogenyType = 'unique sequence based';
        }else{
            phylogenyType = 'taxa based';
        }

        return(
            <div id='phylogeny'>
                <div id='phylogenyMenu'>
                    <h2>phylogenetic analysis - <i>{phylogenyType}</i></h2>
                    <Accordion >
                        {additionalInformation}
                        <Card>
                            <Accordion.Toggle as={Card.Header} eventKey='1'>
                                export tree visualization
                            </Accordion.Toggle>
                            <Accordion.Collapse eventKey='1'>
                                <Card.Body>
                                    Different export options for the tree visualization. The <b>publication ready</b> option will remove the
                                    white spaces between the nodes of the cladogram and visualizes the full tree as static visualization.<br/>
                                    The <b>unlock publication ready</b> button unlocks the static visualization. <br/><br/>
                                    <Nav className='mr-auto'>
                                        <button id="public_ready_phylo" onClick={publicationReady}>publication ready</button>
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