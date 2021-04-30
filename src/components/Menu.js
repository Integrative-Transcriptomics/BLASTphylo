// used libraries
import React, {Component} from 'react';
import axios from 'axios';
import * as d3v6 from 'd3v6';
import {Spinner, Button, Form, Popover, OverlayTrigger} from 'react-bootstrap';
import {BiHelpCircle} from 'react-icons/bi';


// own style sheets
import './menuStyle.css';

// constant dictionary for help messages
const helpMessages = {
    'prot': ['Protein sequence', 'Enter query sequence without white spaces in the text area. Only PROTEIN sequences'],
    'protFile': ['Protein files', 'Either a protein sequence or already calculated BLAST result with the columns: qacc sacc qstart qend sstart send slen nident evalue pident staxids qcovs sseq'],
    'NCBI': ['NCBI taxonomy', 'Enter comma-separated list of scientific names or taxonomic IDs. Addition of \'|subtree\' will select complete subtree. Check https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi for more information' ],
    'ownTaxa': ['Newick string', 'Enter a newick string in the text area or upload a text file which contain a newick string.'],
    'eValue': ['E-value', 'Common filter parameter for BLAST searches'],
    'query': ['Query coverage', 'Result of the BLAST search will be filtered for the minimal query coverage (= alignment length/query length)'],
    'hit': ['Hit or Subject coverage', 'Result of the BLAST search will be filtered for the minimal hit coverage (= alignment length/subject length)'],
    'alignment': ['Alignment identity', 'Result of the BLAST search will be filtered for the minimal alignment identity (= matches/ alignment length)']
};

class Menu extends Component {

    constructor(props) {
        super(props);
        this.protFileInput = React.createRef();
        this.treeFileInput = React.createRef();
        this.state = { protein: '',
                       protein_file_type: '0',
                       tree_data: '',
                       tree_menu_selection: '0',
                       eValue: '0.05',
                       query_cover: '80',
                       align_ident: '80',
                       hit_cover: '50'}

        // functions
        this.handleSubmit = this.handleSubmit.bind(this);
        this.handleChange = this.handleChange.bind(this);

   }


    // handle click on Submit button: extract actual state of the parameters and send them to the back end
    //dependent on the return value the taxonomic Mapping will be visualized
    handleSubmit(){
        const formData = new FormData();
        var error = [];
         if (this.state.protein === ''){
            var protFileContent = this.protFileInput.current.files[0];
            if(typeof protFileContent === 'undefined'){
                error.push({'message': 'Protein is undefined. Please, define a protein'});
            }else{
                formData.append("fasta_file", protFileContent);
                formData.append("fasta_filename", protFileContent.name);
            }
         }
        if(this.state.tree_menu_selection === '1'){
            var treeFileContent = this.treeFileInput.current.files[0];
            console.log(treeFileContent)
            if((typeof treeFileContent === 'undefined') && (this.state.tree_data.length === 0)){
                error.push({'message': 'Taxonomy is undefined. Please, define a taxonomy'});
            }else if(this.state.tree_data.length === 0){
                formData.append('newick_file', treeFileContent);
                formData.append('newick_filename', treeFileContent.name);
            }
        }

        if(error.length === 0){
            for (var key in this.state){
                console.log([key, this.state[key]]);
                formData.append(String(key), this.state[key]);
            }
           /***for(var pair of formData.entries()) {
           console.log(pair[0] + ':  ' + pair[1]);
           }***/

           // switch to loading button
           if ((this.state.protein !== '') || (this.state.protein_file_type === '0') ){
               document.getElementById('submit').style.display = "none";
               document.getElementById('loadingButton').style.display = "block";
           }

           // send data to back end
           var self = this;
            axios.post("server/menu", formData)
             .then(function (response) {
                 console.log(response.data);
                 if(response.data.error === null){
                     if(document.getElementById('alert')){
                        document.getElementById('alert').remove();
                     }
                     self.props.changeComp('data', response.data.tree);
                     self.props.changeComp('actual', 'phyloblast');
                 }else{
                     self.props.changeComp('error', response.data.error);
                     self.props.changeComp('actual', 'menu');
                 }
            })
            .catch(error => {
                console.log(error);
            })
        }
    }

    // Handle changes in the parameter settings
    handleChange(event) {

        if (event.target.name === "fasta_seq"){
          this.setState({protein: event.target.value});
        }
        else if (event.target.name === "file_type"){
          this.setState({protein_file_type: event.target.value});
        }
        else if (event.target.name === "tree_menu"){
          this.setState({tree_menu_selection: event.target.value});
          this.taxaChoice(event);
        }
        else if (event.target.name === "taxa"){
          this.setState({tree_data: event.target.value});
        }
        else if (event.target.name === "newick_string"){
          this.setState({tree_data: event.target.value});
          this.setState({tree_menu_selection: '1'});
        }
        else if (event.target.name === "eValue"){
          this.setState({eValue: event.target.value});
        }
        else if (event.target.name === "align_ident"){
          this.setState({align_ident: event.target.value});
        }
        else if (event.target.name === "query_cover"){
          this.setState({query_cover: event.target.value});
        }
        else if (event.target.name === "hit_ident"){
          this.setState({hit_ident: event.target.value});
        }
        else if (event.target.name === "hit_cover"){
          this.setState({hit_cover: event.target.value});
        }
        else if (event.target.name === "newick_file"){
          this.setState({tree_menu_selection: '1'});
        }
        else {
          console.log('no matching state present')
        }
        //console.log([event.target.name, event.target.value]);
   }

    // handle the switch between the two taxonomy input options
    taxaChoice(event){
        if (event.target.value === '0'){
            document.getElementById('own taxonomy').style.display = "none";
            document.getElementById('NCBI taxonomy').style.display = "block";

        } else {
            document.getElementById('NCBI taxonomy').style.display = "none";
            document.getElementById('own taxonomy').style.display = "block";
        }
   }

    render() {
        // remove old visualisations and reduce size of default divs
        d3v6.select('#tree_vis').remove();
        d3v6.select('#hitbars').remove();
        d3v6.select('#clade_vis').remove();


        // function to generate the help message pop up
        var MakeItem = function(X){
                        return(
                <Popover id='popover-basic' >
                    <Popover.Title as='h3'><strong>{X[0]}</strong></Popover.Title>
                    <Popover.Content>{X[1]}</Popover.Content>
                </Popover>);
        };

        return(

            <div id="menu" >
            <Form as='fieldset' id='BlastpSearch'>
            <legend>Blastp Search</legend>
            <Form.Group >
                <Form inline>
                    <Form.Label>Enter protein sequence:</Form.Label>
                    <OverlayTrigger trigger='click' placement='right' overlay={MakeItem(helpMessages['prot'])}>
                        <BiHelpCircle style={{color: 'blue'}}/>
                    </OverlayTrigger>
                </Form>
                <Form.Control as='textarea' rows={5} cols={50} id='fasta_seq' name='fasta_seq' onChange={this.handleChange} />
            </Form.Group>
            <Form.Group>
                <Form inline>
                    <Form.Label>Or, upload</Form.Label>
                    <Form.Control as='select' id='file_type' name='file_type' onChange={this.handleChange}>
                        <option value="0">Fasta file</option>
                        <option value="1">BLAST result</option>
                    </Form.Control>
                    <OverlayTrigger trigger='click' placement='right' overlay={MakeItem(helpMessages['protFile'])}>
                        <BiHelpCircle style={{color: 'blue', 'margin': '0px 10px 0px 5px'}}/>
                    </OverlayTrigger>
                </Form>
                <Form.File ref={this.protFileInput} id='fasta_file' name='fasta_file'
                onChange={this.handleChange} accept='.fasta,.fastq,.csv' />
            </Form.Group>
            </Form>
            <br />
            <Form as='fieldset' id='Taxonomy'>
            <legend>Taxonomy</legend>
            <Form.Group>
                <Form inline>
                    <Form.Label>Select input type:</Form.Label>
                    <Form.Control as='select' id='tree_menu' name='tree_menu' onChange={this.handleChange}>
                        <option value='0'>NCBI taxonomy</option>
                        <option value='1'>own taxonomic phylogeny</option>
                    </Form.Control>
                </Form>
            </Form.Group>
            <Form.Group id='NCBI taxonomy'>
                <Form inline>
                    <Form.Label>Enter list of taxa as scientific names of taxonomic IDs:</Form.Label>
                    <OverlayTrigger trigger='click' placement='right' overlay={MakeItem(helpMessages['NCBI'])}>
                        <BiHelpCircle style={{color: 'blue'}}/>
                    </OverlayTrigger>
                </Form>
                <Form.Control as='textarea' rows={5} cols={50} id='taxa' name='taxa'
                placeholder='Staphylococcus,Staphylococcus aureus|subtree' onChange={this.handleChange} />
            </Form.Group>
            <Form.Group id='own taxonomy' style={{display: 'none'}}>
                <Form inline>
                    <Form.Label>Enter taxonomy as newick string:</Form.Label>
                    <OverlayTrigger trigger='click' placement='right' overlay={MakeItem(helpMessages['ownTaxa'])}>
                        <BiHelpCircle style={{color: 'blue'}}/>
                    </OverlayTrigger>
                </Form>
                <Form.Control as='textarea' rows={5} cols={50} id='newick_string' name='newick_string'
                placeholder='(A,(C,D)E)F;' onChange={this.handleChange} />
                <Form inline>
                        <Form.Label>Or, upload file:</Form.Label>
                        <Form.File ref={this.treeFileInput} id='newick_file' name='newick_file'
                        onChange={this.handleChange} accept='.txt,.dnd' />
                </Form>
            </Form.Group>
            </Form>
            <br />
            <Form as='fieldset' id='FilterConditions'>
            <legend>Filter conditions</legend>
            <Form.Group inline>
                <Form.Label>E-value:</Form.Label>
                <input type="number" id="eValue"  name='eValue' step='0.01' min='0' max='100'
                value={this.state.eValue} size='6' onChange={this.handleChange} />
            </Form.Group>
            <Form.Group inline>
                <Form.Label>min. alignment identity:</Form.Label>
                <OverlayTrigger trigger='click' placement='right' overlay={MakeItem(helpMessages['alignment'])}>
                    <BiHelpCircle style={{color: 'blue', 'margin': '0px 10px 0px 0px'}}/>
                </OverlayTrigger>
                <input type="number" id="align_ident" name="align_ident" step="1" min='5' max='100'
                    value={this.state.align_ident} size='6' onChange={this.handleChange} />
            </Form.Group>
            <Form.Group inline>
                <Form.Label>min. query coverage:</Form.Label>
                <OverlayTrigger trigger='click' placement='right' overlay={MakeItem(helpMessages['query'])}>
                    <BiHelpCircle style={{color: 'blue', 'margin': '0px 10px 0px 0px'}}/>
                </OverlayTrigger>
                <input type="number" id="query_cover" name="query_cover" step="1" min='5' max='100'
                     value={this.state.query_cover} size='6' onChange={this.handleChange} />
            </Form.Group>
            <Form.Group inline>
                <Form.Label>min. hit coverage:</Form.Label>
                <OverlayTrigger trigger='click' placement='right' overlay={MakeItem(helpMessages['hit'])}>
                    <BiHelpCircle style={{color: 'blue', 'margin': '0px 10px 0px 0px'}}/>
                </OverlayTrigger>
                <input type="number" id="hit_cover" name="hit_cover" step="1" min='5' max='100'
                    value={this.state.hit_cover} size='6' onChange={this.handleChange} />
            </Form.Group>
            </Form>
            <br />
            <Button id='submit' type="submit" value="Submit" onClick={this.handleSubmit}>Submit </Button>
            <div id='loadingButton' style={{display: 'none'}}>
                <Button variant="primary" disabled>
                    <Spinner
                      as="span"
                      animation="border"
                      size="sm"
                      role="status"
                      aria-hidden="true"
                    />
                    Loading...
                </Button>
            </div>
          </div>
        );
    }
}

export default Menu;