// used libraries
import React, {Component} from 'react';
import * as d3v6 from 'd3v6';
import "bootstrap/dist/css/bootstrap.min.css";
import {Table, Form, OverlayTrigger, Tooltip, Nav} from 'react-bootstrap';
import {IoMdArrowDropleft, IoMdArrowDropright} from 'react-icons/io';
import {BiHelpCircle} from 'react-icons/bi';
import {AiFillHome} from 'react-icons/ai';


// own components and style sheets
import HandlePhylogenyServer from './HandlePhylogenyServer.js';
import TaxonomicAnalysisMenu from './TaxonomicAnalysisMenu.js';
import PhylogeneticAnalysisMenu from './PhylogeneticAnalysisMenu.js';
import TreeInteraction from './TreeInteraction.js';




class Tabhandling extends Component {

    constructor(props){
        super(props);
        console.log(props);
        this.state = {taxonomicMapData: this.props.data,
                       taxaPhyloData: null,
                       seqPhyloData:null,
                       isActualComponent: 'taxonomy',
                       phylogenyState: null}

        // functions to handle view update, reset, help page and alert messages
        this.handleTabSwitch = this.handleTabSwitch.bind(this);
        this.saveActualData = this.saveActualData.bind(this);
        this.getActualData = this.getActualData.bind(this);
    }

    handleTabSwitch(event){
        console.log(event);
        if(event === 'phylogeny'){
            this.setState({isActualComponent: 'handlePhylogeny'});
            this.setState({phylogenyState: event});
            /****if(this.state.taxaPhyloData){
                //console.log(this.state.taxaPhyloData)
                this.setState({isActualComponent: event});
                this.setState({phylogenyState: null});
            }else{
                this.setState({isActualComponent: 'handlePhylogeny'});
                this.setState({phylogenyState: event});
            }****/
        }else if(event === 'phylogenyUnique'){
            this.setState({isActualComponent: 'handlePhylogeny'});
            this.setState({phylogenyState: event});
            /***if(this.state.seqPhyloData){
                //console.log(this.state.seqPhyloData)
                this.setState({isActualComponent: event});
                this.setState({phylogenyState: null});
            }else{
                this.setState({isActualComponent: 'handlePhylogeny'});
                this.setState({phylogenyState: event});
            }***/
        }else if(event === 'taxonomy'){
            this.setState({isActualComponent: event});
        }else{
            console.log('Tabhandling - switch Tab: Origin does not match Tabs');
        }
    }

    saveActualData(originData, data){
        //console.log(data)
        //console.log(originData)
        if (originData === 'server/phylogeny'){
            this.setState({taxaPhyloData: data});
            this.setState({isActualComponent: this.state.phylogenyState});
        }else if (originData === 'server/phylogenyUnique'){
            this.setState({seqPhyloData: data});
            this.setState({isActualComponent: this.state.phylogenyState});
        }else{
            console.log('Tabhandling - save data: Origin does not match Tabs');
        }
    }

    getActualData(origin){
        if(origin === 'phylogenyTab'){
            return this.state.taxaPhylodata;
        }else if (origin === 'phylogenyUniqueTab'){
            return this.state.seqPhylodata;
        }else{
            console.log('Tabhandling: Origin does not match Tabs');
        }
    }



    render() {


        var isActualComponent = this.state.isActualComponent;
        const taxonomicMapData = this.state.taxonomicMapData;
        var userMenu = <div />;  // sub div: menus for the visualizations (export visualization, etc)
        console.log(this.state)
        let actualComponent;
        d3v6.select('#phyloblastAlert').remove();
        d3v6.select('#tree_vis').remove();


        // switch view based on the isActualComponent state
        if (isActualComponent === 'taxonomy'){ // taxonomic Mapping
            const copy = {...taxonomicMapData[0]};
            actualComponent = <TaxonomicAnalysisMenu phyloData={copy} queries={taxonomicMapData[1]}  />;
            userMenu = <TreeInteraction data={{tree: copy}} queries={taxonomicMapData[1]} calculationMethod={'taxa'}/>;

        }else if (isActualComponent === 'handlePhylogeny'){ // handle data processing of the phylogenetic analysis
            actualComponent = <HandlePhylogenyServer data={this.state.phylogenyState} saveActualData={this.saveActualData} />;
        }else if (isActualComponent === 'phylogeny'){ // visualize the phylogeny
            //console.log(this.state.taxaPhyloData)
            actualComponent = <PhylogeneticAnalysisMenu data={this.state.taxaPhyloData} phylogenyState={'taxa'}/>;
            userMenu = <TreeInteraction data={this.state.taxaPhyloData} calculationMethod={'phylo'}  />;
        }else if (isActualComponent === 'phylogenyUnique'){ // visualize the phylogeny
            //console.log(this.state.seqPhyloData)
            actualComponent = <PhylogeneticAnalysisMenu data={this.state.seqPhyloData} phylogenyState={'seq'} />;
            userMenu = <TreeInteraction data={this.state.seqPhyloData} calculationMethod={'phyloUnique'} />;
        }else {
            d3v6.select('#visualisation').style('border', 'hidden');
            console.log('Tabhandling: No tab from selection')
        }
        console.log(this.state.isActualComponent)


            return (
                    <div className='Tabhandling'>
                        <Nav variant='tabs' defaultActiveKey='taxonomy' onSelect={this.handleTabSwitch}>
                            <Nav.Item>
                                <Nav.Link id='taxonomyTab' eventKey='taxonomy'>taxonomic Mapping</Nav.Link>
                            </Nav.Item>
                            <Nav.Item>
                                <Nav.Link id='phylogenyTab' eventKey='phylogeny'>taxa-centered Phylogeny</Nav.Link>
                            </Nav.Item>
                            <Nav.Item>
                                <Nav.Link id='phylogenyUniqueTab' eventKey='phylogenyUnique'>sequence-centered Phylogeny</Nav.Link>
                            </Nav.Item>
                        </Nav>
                        <div className='App-body'>
                            {actualComponent}
                            <div id="visualisation">
                                <div id='treeInteraction'>
                                {userMenu}
                                </div>
                                <div id='treeVis'>
                                    <div id="tree"></div>
                                    <div id="additionalInfo"></div>
                                </div>
                            </div>
                        </div>
                    </div>

            );
    }
}

export default Tabhandling;
