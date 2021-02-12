// used libraries
import React, {Component} from 'react';
import axios from 'axios';
import * as d3v6 from 'd3v6';
import {Spinner, Button} from 'react-bootstrap';

// own style sheets
import './menuStyle.css';

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
                   query_ident: '80',
                   hit_cover: '50',
                   hit_ident: '50'}

    this.handleSubmit = this.handleSubmit.bind(this);
    this.handleChange = this.handleChange.bind(this);

   }

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
        if(typeof treeFileContent === 'undefined'){
            error.push({'message': 'Taxonomy is undefined. Please, define a taxonomy'});
        }else{
            formData.append('newick_file', treeFileContent);
            formData.append('newick_filename', treeFileContent.name);
        }
    }

    if(error.length === 0){
        for (var key in this.state){
            formData.append(String(key), this.state[key]);
        }
       // for(var pair of formData.entries()) {
       // console.log(pair[0] + ':  ' + pair[1]);
       // }

       // switch to loading button
       if ((this.state.protein !== '') || (this.state.protein_file_type === '0') ){
           document.getElementById('submit').style.display = "none";
           document.getElementById('loadingButton').style.display = "block";
       }


       var self = this;
        axios.post("server/phyloblast", formData)
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
    }else{
        this.props.changeComp('error', error);
        this.props.changeComp('actual', 'menu');
    }


   }

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
    else if (event.target.name === "query_ident"){
      this.setState({query_ident: event.target.value});
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
   }


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
    d3v6.select('#tree_vis').remove();

    return(

        <div id="menu" >
        <fieldset>
		<legend>BLASTp search</legend>
        <label>
        Enter protein sequence:</label>
        <br></br>
        <textarea id="fasta_seq" name="fasta_seq" cols="50" rows="5" onChange={this.handleChange}></textarea>
        <br></br>
        <label>Or, upload </label>
        <select id="file_type" name="file_type" onChange={this.handleChange}>
	         <option value="0">Fasta file</option>
             <option value="1">BLAST result</option>
        </select>
        <input  ref={this.protFileInput}  type="file" id="fasta_file" name="fasta_file" accept=".fasta,.fastq, .csv"/>
        </fieldset>
        <br></br>
        <br></br>
	    <fieldset>
		<legend>Taxonomy</legend>
        <label>Select input type</label>
		    <select  id="tree_menu" name="tree_menu" onChange={this.handleChange}>
			    <option value='0'>NCBI taxonomy</option>
			    <option value='1'>own taxonomic phylogeny</option>
		    </select>
        <div id="NCBI taxonomy">
		  <label>Enter list of taxa as scientific names or taxIDs</label><br></br>
		  <textarea id="taxa" name="taxa" placeholder="Staphylococcus,Staphylococcus aureus|subtree"
			cols="50" rows="5" onChange={this.handleChange}></textarea>
	     </div>
        <div id="own taxonomy" style={{display: 'none'}}>
		  <label>Enter newick string:</label><br></br>
  		  <textarea id="newick_string" name="newick_string" placeholder="(A:0.1,(C:0.3,D:0.4)E:0.5)F;"
			cols="50" rows="5" onChange={this.handleChange}></textarea><br></br>
        	   <label>Or, upload dnd file:</label>
		   <input ref={this.treeFileInput} type="file" id="newick_file" name="newick_file" onChange={this.handleChange} accept=".txt, .dnd" /><br></br>
	     </div>
	     </fieldset>
         <br></br>
        <fieldset>
		<legend>Filter conditions</legend>
         <label>E-value:</label>
	     <input type="number" id="eValue"  name='eValue' step='0.01' min='0' max='100' value={this.state.eValue} size='4' onChange={this.handleChange} /><br></br><br></br>
	     <strong>Query filter:</strong><br></br>
	     <label>min. query identity:</label>
	     <input type="number" id="query_ident" name="query_ident" step="1" min='5' max='100' value={this.state.query_ident} size='4' onChange={this.handleChange} /><br></br>
	     <label>min. query coverage:</label>
	     <input type="number" id="query_cover" name="query_cover" step="1" min='5' max='100' value={this.state.query_cover} size='4' onChange={this.handleChange} />
         <br></br><br></br>
         <strong>Hit filter:</strong><br></br>
	     <label>min. hit identity:</label>
	     <input type="number" id="hit_ident" name="hit_ident" step="1" min='5' max='100' value={this.state.hit_ident} size='4' onChange={this.handleChange} /><br></br>
	     <label>min. hit coverage:</label>
	     <input type="number" id="hit_cover" name="hit_cover" step="1" min='5' max='100' value={this.state.hit_cover} size='4' onChange={this.handleChange} />
           <br></br>
            </fieldset>
            <br></br>
           <button id='submit' type="submit" value="Submit" onClick={this.handleSubmit}>Submit </button>
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
           <br></br>
      </div>
    );
  }
}

export default Menu;