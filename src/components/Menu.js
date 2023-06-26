// used libraries
import React, {Component, useState} from 'react';
import axios from 'axios';
import * as d3v6 from 'd3v6';
import {Spinner, Button, Form, Popover, OverlayTrigger, Dropdown, FormControl} from 'react-bootstrap';
import {BiHelpCircle} from 'react-icons/bi';
import {ReactSearchAutocomplete } from 'react-search-autocomplete';

// own style sheets
import './menuStyle.css';
import {checkBrowser} from '../App.js';
import SearchBar from './SearchBar.js'

// constant dictionary for help messages
const helpMessages = {
    'blasttype': ['BLAST search', 'Select the BLAST search you want to use.'],
    'prot': ['Sequence', 'Enter query sequence without white spaces in the text area'],
    'protFile': ['FastA files', 'Upload either a Fasta sequence file or an already calculated BLAST result with the columns: qacc sacc qstart qend sstart send slen nident evalue pident staxids qcovs sseq'],
    'NCBI': ['NCBI taxonomy', 'Enter a comma-separated list of scientific names or taxonomic IDs without spaces. Adding of \'|subtree\' will select the complete subtree. Check https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi for more information' ],
    'ownTaxa': ['Newick string', 'Enter a newick string in the text area or upload a text file which contains a newick string.'],
    'eValue': ['E-value', 'Common filter parameter for BLAST searches'],
    'query': ['Query coverage', 'Results of the BLAST search will be filtered for the minimal query coverage (= alignment length/query length)'],
    'hit': ['Hit or Subject coverage', 'Results of the BLAST search will be filtered for the minimal hit coverage (= alignment length/subject length)'],
    'alignment': ['Alignment identity', 'Results of the BLAST search will be filtered for the minimal alignment identity (= matches/alignment length)']
};

class Menu extends Component {

    constructor(props) {
        super(props);
        this.protFileInput = React.createRef();
        this.treeFileInput = React.createRef();
        this.state = { blasttype: 'blastp',
                       protein: '',
                       protein_file_type: '0',
                       example: '',
                       tree_data: '',
                       tree_menu_selection: '0',
                       eValue: '0.05',
                       query_cover: '80',
                       align_ident: '80',
                       hit_cover: '50'}


        // functions
        this.handleSubmit = this.handleSubmit.bind(this);
        this.handleChange = this.handleChange.bind(this);
        this.inputField = this.inputField.bind(this);
        //this.handleChangeSearchbar = this.handleChangeSearchbar(this);

   }



    // handle click on Submit button: extract actual state of the parameters and send them to the back end
    //dependent on the return value the taxonomic Mapping will be visualized
    handleSubmit(){
        const formData = new FormData();
        var error = [];
         if (this.state.protein_file_type === '2') {
            formData.append('fasta_file', null);
            formData.append('fasta_filename', this.state.example);
         }
         else if (this.state.protein === ''){
            var protFileContent = this.protFileInput.current.files[0];
            if(typeof protFileContent === 'undefined'){
                error.push({'message': 'Sequence is undefined. Please, enter a nucleotide or protein sequence'});
            }else{
                formData.append("fasta_file", protFileContent);
                formData.append("fasta_filename", protFileContent.name);
            }
         }
        if(this.state.tree_menu_selection === '1'){
            var treeFileContent = this.treeFileInput.current.files[0];
            if((typeof treeFileContent === 'undefined') && (this.state.tree_data.length === 0)){
                error.push({'message': 'Taxonomy is undefined. Please, define a taxonomy'});
            }else if(this.state.tree_data.length === 0){
                formData.append('newick_file', treeFileContent);
                formData.append('newick_filename', treeFileContent.name);
            }
        }
        if((this.state.tree_menu_selection === '0') && (this.state.tree_data.length === 0)){
            error.push({'message': 'Taxonomy is undefined. Please, define a taxonomy'});
        }
        if(error.length === 0){
            for (var key in this.state){
                //console.log([key, this.state[key]]);
                formData.append(String(key), this.state[key]);
            }
           for(var pair of formData.entries()) {
                console.log(pair[0] + ':  ' + pair[1]);
           }

           // switch to loading button
           document.getElementById('submit').style.display = "none";
           document.getElementById('loadingButton').style.display = "block";

           // send data to back end
           var self = this;
           axios.post('server/menu', formData)
             .then(function (response) {
                 console.log('responsedata', response.data);
                 if(response.data.error === null){
                     if(document.getElementById('alert')){
                        document.getElementById('alert').remove();
                     }
                     self.props.changeComp('data', response.data);
                     self.props.changeComp('actual', 'tabhandling');
                 }else{
                     if(response.data.error === 'noHits'){
                        alert('No significant hits found within the taxonomy specified.')}
                     else if(response.data.error === 'treeD3failed'){
                        alert('Failed to generate D3 tree.')}
                     else if(response.data.error === 'undefined'){
                        alert('An undefined error occured.')}
                     else if(response.data.error === 'wrongColumns'){
                        alert('Your uploaded BLAST result contains invalid values.')}
                     else if(response.data.error === 'blastfailed'){
                        alert('BLAST run failed due to server timeout. We kindly ask the user to run BLAST oneself and upload te result to BLASTphylo.')}
                     else {
                        console.log('here');
                        let error_message = '';
                        for (let i = 0; i < response.data.error.length; i++) {
                            if (error_message === '') {
                                error_message = error_message.concat(response.data.error[i]['message']);
                            } else {
                                error_message = error_message.concat('\n').concat(response.data.error[i]['message']);
                            }
                        }

                        console.log(response.data.error[0]['message']);
                        alert(error_message);
                        }

                     //self.props.sendErrors(response.data.error);
                     self.props.changeComp('actual', 'menu');
                     //document.getElementById('submit').style.display = "none";
                     document.getElementById('loadingButton').style.display = "none";
                     document.getElementById('submit').style.display = 'initial';
                 }
            })
            .catch(error => {
                alert('BLAST run failed due to server timeout (max. 5 min). We kindly ask the user to run BLAST himself and upload the result to BLASTphylo.');
                self.props.changeComp('actual', 'menu');
                document.getElementById('loadingButton').style.display = "none";
                document.getElementById('submit').style.display = 'initial';
            }
           )
        }else{
            this.props.sendErrors(error);
        }
    }


    handleChangeSearchbar = (value) => {
        this.setState({tree_data: value});
        console.log('change');
    }

    // Handle changes in the parameter settings
    handleChange(event) {
        if (event.target.name === 'blast'){
            this.setState({blasttype: event.target.id})
        }
        else if (event.target.name === "fasta_seq"){
            this.setState({protein: event.target.value});
        }
        else if (event.target.name === "file_type"){
             this.setState({protein_file_type: event.target.value});
        }
        else if (event.target.name === "example"){
            this.setState({example: event.target.value})
            document.getElementById('taxa').value = 'Bacteria|subtree';
            this.setState({tree_data: '2|subtree'});
        }
        else if (event.target.name === "tree_menu"){
            this.setState({tree_menu_selection: event.target.value});
            this.taxaChoice(event);
        }
        else if (event.target.name === "taxa"){
            this.setState({tree_data: event.target.value});
            console.log('change');
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
   inputField = () => {
   // check browser and switch newline symbol
         const actualBrowser = checkBrowser();
         let placeholderFasta;
         if(actualBrowser === 'safari'){
            placeholderFasta = '>query1\nMEMEFNENNIDLETIIRDEVNKYLSRDI\nGDLPATQQAPLELREKYEKMEVPNKGRDIYEV';
         }else{
            placeholderFasta = '>query1\nMEMEFNENNIDLETIIRDEVNKYLSRDI\nGDLPATQQAPLELREKYEKMEVPNKGRDIYEV';
         }
         let input_file;
         input_file=<Form.File ref={this.protFileInput} id='fasta_file' name='fasta_file'
                onChange={this.handleChange} accept='.fasta,.fastq,.csv,.fna' />;
    switch (this.state.protein_file_type) {
      case "0":
        return (
          <Form.Group>
            <Form inline>
                    <Form.Label>Enter nucleotide/protein sequence:</Form.Label>
                </Form>
                <Form.Control as='textarea' placeholder={placeholderFasta} rows={5} cols={50} id='fasta_seq' name='fasta_seq' onChange={this.handleChange} />
                <Form inline>
                    <Form.Label>Or upload a Fasta file:</Form.Label>
                </Form>
                <Form.File ref={this.protFileInput} id='fasta_file' name='fasta_file'
                onChange={this.handleChange} accept='.fasta,.fastq,.fna' />
          </Form.Group>
        );
      case "1":
        return (
          <Form.Group>
             <Form inline>
                    <Form.Label>Upload BLAST result table:</Form.Label>
             </Form>
             <Form.File ref={this.protFileInput} id='blast_file' name='blast_file'
                onChange={this.handleChange} accept='.csv' />
          </Form.Group>
        );
      case "2":
        return (<Form.Group>
             <Form inline>
                    <Form.Label>Load example data:</Form.Label>
             </Form>
             <Form id='example_data' name='example_data'>
                        <Form.Check name='example' type='radio' value='2' id='example0' label={'One query sequence'} onChange={this.handleChange} />
                        <Form.Check name='example' type='radio' value='3' id='example1' label={'Comparing two query sequences'} onChange={this.handleChange} />
                </Form>
          </Form.Group>
        );
      default:
        return null;
    }
    };


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

         // check browser and switch newline symbol
         const actualBrowser = checkBrowser();
         let placeholderFasta;
         if(actualBrowser === 'safari'){
            placeholderFasta = '>query1&#x0a; MEMEFNENNIDLETIIRDEVNKYLSRDI&#x0a; GDLPATQQAPLELREKYEKMEVPNKGRDIYEV';
         }else{
            placeholderFasta = '>query1\nMEMEFNENNIDLETIIRDEVNKYLSRDI\nGDLPATQQAPLELREKYEKMEVPNKGRDIYEV';
         }

        return(

            <div id="menu" >
            <Form as='fieldset'>
            <legend>About</legend>
            <Form inline>
            <Form.Label> BLASTphylo is a web tool to visualize the occurrence of a protein in the taxonomy for different taxonomic ranks and to visualize the phylogenetic relationships of these proteins as an interactive tree.
            For more information on the input formats or potential user interactions with BLASTphylo's visualizations, please visit the help page found by clicking on the '?' (upper right corner).
            </Form.Label>
            <br/>
            <Form.Label class="text-warning"> Please note that the current version of BLASTphylo is limited to bacterial proteins and genes.</Form.Label>
            </Form>
            </Form>
            <br/>
            <Form as='fieldset' id='BlastpSearch'>
            <legend>Blast Search</legend>
            <Form.Group >
                <Form inline>
                    <Form.Label>Select your BLAST search:</Form.Label>
                    <OverlayTrigger trigger='click' placement='right' overlay={MakeItem(helpMessages['blasttype'])}>
                        <BiHelpCircle style={{color: 'blue'}}/>
                    </OverlayTrigger>
                </Form>
                <Form id='blasttypes' inline>
                        <Form.Check inline type={'radio'} name='blast' id='blastp' label={'blastp'} defaultChecked onChange={this.handleChange} />
                        <Form.Check inline type={'radio'} name='blast' id='blastn' label={'blastn'} onChange={this.handleChange} />
                        <Form.Check inline type={'radio'} name='blast' id='blastx' label={'blastx'} onChange={this.handleChange} />
                </Form>
                <br/>
                <Form inline>
                    <Form.Label>Select your input type or choose from several demo examples:</Form.Label>
                </Form>
                <Form id='input_type' name='input_type' inline>
                        <Form.Check inline type={'radio'} name='file_type' value='0' id='fasta' label={'Query sequence'} defaultChecked onChange={this.handleChange} />
                        <Form.Check inline type={'radio'} name='file_type' value='1' id='blastresult' label={'BLAST result'} onChange={this.handleChange} />
                        <Form.Check inline type={'radio'} name='file_type' value='2' id='demo' label={'Demo'} onChange={this.handleChange} />
                </Form>
                {this.inputField()}
                <br/>

            </Form.Group>
            </Form>
            <br />
            <Form as='fieldset' id='Taxonomy'>
            <legend>Taxonomy</legend>
            <Form.Group>
                <Form inline>
                    <Form.Label>Use:</Form.Label>
                    <Form.Control as='select' id='tree_menu' name='tree_menu' onChange={this.handleChange}>
                        <option value='0'>NCBI taxonomy</option>
                        <option value='1'>own taxonomic phylogeny</option>
                    </Form.Control>
                </Form>
            </Form.Group>
            <Form.Group id='NCBI taxonomy'>
                <Form inline>
                    <Form.Label>Specify taxon name(s) or taxID(s) to restrict the search to part of the taxonomy:</Form.Label>
                    <OverlayTrigger trigger='click' placement='right' overlay={MakeItem(helpMessages['NCBI'])}>
                        <BiHelpCircle style={{color: 'blue'}}/>
                    </OverlayTrigger>
                </Form>
                <Form.Control as='textarea' rows={5} cols={50} id='taxa' name='taxa'
                placeholder='Example formats: &#10;Bacterial subtree: Bacteria|subtree or 2|subtree &#10;Specific taxon: Staphylococcus or 1279 &#10;Multiple taxa: 1279|subtree,Proteobacteria|subtree &#10;Subtree of specific taxon: Staphylococcus|subtree or 1279|subtree &#10;Exclude part of a subtree: 1279|!(1280|subtree)' onChange={this.handleChange} />
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
                <Form.Group>
                <Form inline>
                    <Form.Label>Optional: Search for taxon name or taxId to be used in the field above:</Form.Label>
                </Form>
                </Form.Group>
                <SearchBar handleChange={this.handleChangeSearchbar}/>
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