// used libraries
import React, {Component} from 'react';
import * as d3v6 from 'd3v6';
import * as d3 from 'd3';
import * as $ from "jquery";
import * as _ from "lodash";
import "../../node_modules/phylotree/src/main";
import {Nav, Accordion, Card, Tooltip, OverlayTrigger, Form} from 'react-bootstrap';
import {BsBoxArrowUpRight} from "react-icons/bs";

// own components and style sheets
import './phylogenyStyle.css';
import {ExportTrees} from '../App.js'


// own visualisations
import {chart, showClades, startTreevis, publicationReady} from '../visualisations/phyloblast2.js';


class PhylogeneticAnalysisMenu extends Component{
     constructor(props){
        super(props);
        console.log(this.props)

        var treeData = this.props.data.tree;
        if( Object.keys(this.props.data).length === 3){
            this.state = {tree: treeData, extra: this.props.data.extraInfo, actualTree: null, counter: 0};
        }else{
            this.state = {tree: treeData, extra: [1], actualTree: null, counter: 0};
        }
        console.log(this.state);

        // own functions
        this.handleMenuClick = this.handleMenuClick.bind(this);
        this.handleHelpClick = this.handleHelpClick.bind(this);
        this.showAdditional = this.showAdditional.bind(this);
    }

    handleMenuClick(){
        this.props.changeComp('actual', 'menu');
    }

    handleHelpClick(){
        this.props.changeComp('actual', 'help');
    }

    showAdditional(event){
        //console.log(event.target.value)
        //console.log(this.state.extra);
        if(event.target.value === '3'){
            d3v6.select('#clade_vis').remove();
        }else{
            showClades(this.state.extra[0], this.state.extra[1], this.state.actualTree, event.target.value);
        }
    }

    render(){
        d3v6.select('#visualisation').style('border', '2px solid #5e66b4')
                                .style('border-radius', '5px');


        let additionalInformation;
        let labeling;
        if (this.state.extra.length === 2){
            additionalInformation = additionalCladeInformation();
        }
        d3v6.selectAll("input").on("change", this.showAdditional);

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
                                        <button id="public_ready"  onClick={publicationReady}>publication ready</button>
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

function additionalCladeInformation(){
    return(
        <Card>
            <Accordion.Toggle as={Card.Header} eventKey='0'>
                visualise additional clade information
            </Accordion.Toggle>
            <Accordion.Collapse eventKey='0'>
                <Card.Body>
                    <Form id='infoSelection'>
                        <label>
                            unique seqs.
                            <input type='radio' name='valueSel' value='0' />
                        </label>
                        <label>
                            phylum rank
                            <input type='radio' name='valueSel' value='1' />
                        </label>
                        <label>
                            both
                            <input type='radio' name='valueSel' value='2' />
                        </label>
                        <label>
                            none
                            <input type='radio' name='valueSel' value='3' />
                        </label>
                    </Form>
                </Card.Body>
            </Accordion.Collapse>
        </Card>
    );

}

export default PhylogeneticAnalysisMenu;