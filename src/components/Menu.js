// used libraries
import React, {Component} from 'react';
import axios from 'axios';
import * as d3v6 from 'd3v6';

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
     // send data to flask: https://stackoverflow.com/questions/59126998/how-to-pass-data-from-react-form-flask-backend-react-component-does-it-ha
     const formData = new FormData();
     if (this.state.protein === ''){
        formData.append("fasta_file", this.protFileInput.current.files[0]);
        formData.append("fasta_filename", this.protFileInput.current.files[0].name);
     }
    if(this.state.tree_menu_selection === '1'){

        formData.append('newick_file', this.protFileInput.current.files[0]);
        formData.append('newick_filename', this.protFileInput.current.files[0].name);
      alert(
      `Selected tree file - ${this.treeFileInput.current.files[0].name}`);
    }
    for (var key in this.state){
        formData.append(String(key), this.state[key]);
    }
   // for(var pair of formData.entries()) {
   // console.log(pair[0] + ':  ' + pair[1]);
   // }

   var self = this;
    axios.post("server/phyloblast", formData)
     .then(function (response) {
         console.log(response.data);
         if(response.data.error === null){
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

   handleChange(event) {
   // handle files: https://medium.com/excited-developers/file-upload-with-react-flask-e115e6f2bf99

    if (event.target.name === "fasta_seq"){
      this.setState({protein: event.target.value});
    }
    else if (event.target.name === "file_type"){
      this.setState({protein_file_type: event.target.value});
    }
    else if (event.target.name === "tree_menu"){
      this.setState({tree_menu_selection: event.target.value});
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
        <br></br>
        <label>Or, upload </label>
        <select id="file_type" name="file_type" onChange={this.handleChange}>
	     <option value="0">Fasta file</option>
             <option value="1">BLAST result</option>
        </select>
        <input  ref={this.protFileInput}  type="file" id="fasta_file" name="fasta_file"/>
        </fieldset>
        <br></br>
        <br></br>
	    <fieldset>
		<legend>Taxonomy</legend>
        <label>Select input type</label>
		    <select  id="tree_menu" name="tree_menu" onChange={this.taxaChoice}>
			    <option value='0'>NCBI taxonomy</option>
			    <option value='1'>own taxonomic phylogeny</option>
		    </select>
             <br></br>
             <br></br>
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
		   <input ref={this.treeFileInput} type="file" id="newick_file" name="newick_file" onChange={this.handleChange}/><br></br>
	     </div>
	     </fieldset>
         <br></br>
        <fieldset>
		<legend>Filter conditions</legend>
         <label>E-value:</label>
	     <input type="number" id="eValue"  name='eValue' step='0.01' min='0' max='100' value={this.state.eValue} size='3' onChange={this.handleChange} /><br></br><br></br>
	     <strong>Query filter:</strong><br></br>
	     <label>min. query identity:</label>
	     <input type="number" id="query_ident" name="query_ident" step="1" min='50' max='100' value={this.state.query_ident} size='2' onChange={this.handleChange} /><br></br><br></br>
	     <label>min. query coverage:</label>
	     <input type="number" id="query_cover" name="query_cover" step="1" min='50' max='100' value={this.state.query_cover} size='2' onChange={this.handleChange} />
         <br></br><br></br>
         <strong>Hit filter:</strong><br></br>
	     <label>min. hit identity:</label>
	     <input type="number" id="hit_ident" name="hit_ident" step="1" min='50' max='100' value={this.state.hit_ident} size='2' onChange={this.handleChange} /><br></br><br></br>
	     <label>min. hit coverage:</label>
	     <input type="number" id="hit_cover" name="hit_cover" step="1" min='50' max='100' value={this.state.hit_cover} size='2' onChange={this.handleChange} />
           <br></br>
            </fieldset>
            <br></br>
           <button id='submit' type="submit" value="Submit" onClick={this.handleSubmit}>Submit </button>
           <br></br>
      </div>
    );
  }
}

export default Menu;