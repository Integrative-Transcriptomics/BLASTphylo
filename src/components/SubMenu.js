// used libraries
import React, {Component} from 'react';
import * as d3v6 from 'd3v6';
import {Nav, Navbar, NavItem, NavDropdown, Dropdown, DropdownButton, Form, FormControl, Button} from 'react-bootstrap';
import {Link, BrowserRouter, Route, Switch } from "react-router-dom";


// own components and style sheets
import './submenuStyle.css';
import {ExportTrees} from '../App.js'
import Phylogeny from './Phylogeny.js';


// own visualisations
import {chart, hitBars, startTreevis, collapseTree, publicationReady} from '../visualisations/phyloblast.js';


class SubMenu extends Component{
     constructor(props){
        super(props);
        console.log(this.props)

        this.state = {hitSelect: "2",
                      rankSelect: 'class'};

        this.handleHits = this.handleHits.bind(this);
        this.handleRanks = this.handleRanks.bind(this);

        d3v6.select('#treeVis').style('border', '2px solid #69a2c9')
                                .style('border-radius', '5px');

        // generation of the tree visualisation
        var tree = startTreevis(this.props.phyloData);
        if (tree !== 0){
            chart(tree, null, 6, 6, true);
        }else{
            d3v6.select('#phyloblastAlert').remove();
            d3v6.select('#tree').append('div').attr('id', 'phyloblastAlert')
                    .text('Found 0 hits. Return to the main page and try another phylogentic tree');
        }
    }

    handleHits(event){
        this.setState({hitSelect: event});
        console.log(event)
        if(event !== '2'){
            hitBars();
        }else{
            d3v6.select('#hitbars')
                .remove();
        }
    }

    handleRanks(event){
        if(document.getElementById('public_ready')){
              this.setState({rankSelect: event.target.value});
              collapseTree(this.props.phyloData);
        }
    }

    render(){
            var taxonomyLevel = ['life', 'domain', 'superkingdom', 'kingdom', 'clade', 'phylum', 'class', 'order', 'family', 'genus', 'species group','species', 'strain'];
            var MakeItem = function(X){
                    return <NavDropdown.Item eventKey={X}>{X}</NavDropdown.Item>;
            };
        return(
            <div id='phyloblast'>
                <div id='SubMenu'>
                    <h2>taxonomic Analysis</h2>
                    <Navbar id='menubar' collapseOnSelect expand="lg" bg="lg" >
                      <Navbar.Brand  id='menutitle' >menu</Navbar.Brand>
                      <Navbar.Toggle id='menubutton' aria-controls="basic-navbar-nav" />
                      <Navbar.Collapse id="basic-navbar-nav">
                        <Nav className="mr-auto">phylogeny calculation
                          <Nav.Link href="/phylogeny" target='_blank'>taxa sequence based phylogeny</Nav.Link>
                          <Nav.Link href="/phylogenyUnique" target='_blank'>unique sequence based phylogeny</Nav.Link>
                        </Nav>
                        <Nav  className='mr-auto'>tree interactions
                          <NavDropdown onSelect={this.handleRanks} title="collapse to" id="collasible-nav-dropdown">
                               {taxonomyLevel.map(MakeItem)}
                          </NavDropdown>
                          <NavDropdown onSelect={this.handleHits} title="show barchart of" id="collasible-nav-dropdown">
                            <NavDropdown.Item eventKey='2'>none</NavDropdown.Item>
                            <NavDropdown.Item eventKey='0'>node hits</NavDropdown.Item>
                            <NavDropdown.Item eventKey='1'>subtree hits</NavDropdown.Item>
                          </NavDropdown>
                        </Nav>
                        <Nav id='lastChapter' className='mr-auto'>export tree
                            <Form inline>
                                <button id="public_ready" title="Result will be a static tree. Reload the page to start a new interpolation of the basic tree" onClick={publicationReady}>publication ready</button>
                                <ExportTrees />
                            </Form>
                        </Nav>
                      </Navbar.Collapse>
                    </Navbar>
                    <Navbar id='treeInteractions' collapseOnSelect expand='lg' bg='lg'>
                     <Navbar.Brand id='treeMenu'>tree interactions</Navbar.Brand>
                     <Navbar.Toggle id='treeMenuButton' aria-controls='basic-navbar-nav' />
                     <Navbar.Collapse  id='basic-nabar-nav'>
                        <Nav className='mr-auto'>
                          <NavDropdown onSelect={this.handleRanks} title="collapse to" id="collasible-nav-dropdown">
                               {taxonomyLevel.map(MakeItem)}
                          </NavDropdown>
                          <NavDropdown onSelect={this.handleHits} title="show barchart of" id="collasible-nav-dropdown">
                            <NavDropdown.Item eventKey='2'>none</NavDropdown.Item>
                            <NavDropdown.Item eventKey='0'>node hits</NavDropdown.Item>
                            <NavDropdown.Item eventKey='1'>subtree hits</NavDropdown.Item>
                          </NavDropdown>
                        </Nav>
                     </Navbar.Collapse>
                    </Navbar>
                </div>
            </div>
        );
    }
}


export default SubMenu;