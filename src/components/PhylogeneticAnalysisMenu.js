// used libraries
import React, {Component} from 'react';
import * as d3v6 from 'd3v6';
import * as d3 from 'd3';
import * as $ from "jquery";
import * as _ from "lodash";
import "../../node_modules/phylotree/src/main";
import {Nav, Accordion, Card, Tooltip, OverlayTrigger, Form} from 'react-bootstrap';

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

        var treeData = this.props.data.tree;
        if( Object.keys(this.props.data).length === 3){
            this.state = {tree: treeData, extra: this.props.data.extraInfo, actualTree: treeData, counter: 0};
        }else{
            this.state = {tree: treeData, extra: [1], actualTree: treeData, counter: 0};
        }

        // own functions
        this.showAdditional = this.showAdditional.bind(this);
    }

    // show additional information for the taxa-based phylogeny
    showAdditional(event){
        showClades(this.state.extra[0], this.state.extra[1], null);
    }




    render(){
        d3v6.select('#visualisation').style('border', '2px solid #5e66b4')
                                .style('border-radius', '5px');

        let additionalInformation;
        if (this.state.extra.length === 2){
            additionalInformation = additionalCladeInformation(this.showAdditional);
        }

        const renderPublicReadyTooltip = (props) => (
            <Tooltip id='public_ready_tooltip' {... props}>
                Result will be a static tree.
                Reload the page to start a new interpolation of the basic tree
            </Tooltip>
        )

        return(
            <div id='phylogeny'>
                <div id='phylogenyMenu'>
                    <h2>phylogenetic analysis</h2>
                    <Accordion >
                        {additionalInformation}
                        <Card>
                            <Accordion.Toggle as={Card.Header} eventKey='1'>
                                export tree visualisation
                            </Accordion.Toggle>
                            <Accordion.Collapse eventKey='1'>
                                <Card.Body>
                                    <Nav className='mr-auto'>
                                        <OverlayTrigger placement='bottom' delay={{show:150, hide:50}} overlay={renderPublicReadyTooltip}>
                                        <button id="public_ready" onClick={publicationReady}>publication ready</button>
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



function additionalCladeInformation(showAdditional){
    return(
        <Card>
            <Accordion.Toggle as={Card.Header} eventKey='0'>
                visualise additional clade information
            </Accordion.Toggle>
            <Accordion.Collapse eventKey='0'>
                <Card.Body>
                    <Form id='infoSelection'>
                        <div key={'inline-checkbox'} className='mb-3'>
                            <Form.Check inline type={'checkbox'}  id='uniqueCheck' label={'unique seqs.'}  onChange={showAdditional}/>
                            <Form.Check inline type={'checkbox'}  id='phylumCheck' label={'phylum rank'} onChange={showAdditional}/>
                        </div>
                    </Form>
                </Card.Body>
            </Accordion.Collapse>
        </Card>
    );

}

export default PhylogeneticAnalysisMenu;