import React from 'react';
//import NewWindow from 'react-new-window'
import axios from 'axios';
import * as d3 from 'd3';


// stylesheets 
import './baseStyle.css';
import './menuStyle.css';


// own javascripts 
import {chart, hitBars, showClades, startTreevis} from './phyloblast.js';


class App extends React.Component {
  
  constructor(props){
    super(props);
    //console.log(this.props)
    this.changeComp = this.changeComp.bind(this);
    this.handleHomeClick = this.handleHomeClick.bind(this);
    this.handleHelpClick = this.handleHelpClick.bind(this);
   

    this.state = {isActualComponent: 'menu',
                   data: null};
  }	

  changeComp(stateName, stateValue){
    if(stateName === 'actual'){
    	this.setState({isActualComponent: stateValue});
    } else if (stateName === 'data'){
        this.setState({data: stateValue});
    } else {
      console.log('Wrong stateName: ' + stateName + 'was given') 
    }
  }
      
  handleHomeClick(){
    this.setState({isActualComponent: 'menu'});
  } 

  
  handleHelpClick(){
    this.setState({isActualComponent: 'help'});	
  }

  render() {

   const isActualComponent = this.state.isActualComponent;
   const data = this.state.data;
   //console.log(data)
   let actualComponent;
   if (isActualComponent === 'help'){
     actualComponent = <Help />;
   }
   else if (isActualComponent === 'wait'){
     actualComponent = <Wait />;
   }
   else if (isActualComponent === 'phyloblast'){
     actualComponent = <Phyloblast phyloData={data} changeComp={this.changeComp}/>;
   }else {
     actualComponent = <Menu isActualComponent={isActualComponent} changeComp={this.changeComp}/>;
   }
   console.log(this.state)
   
  return (
    <div className="App">
      <header className="App-header">
      <nav>
	  <h1>PhyloBlast</h1>
 	  <ul>
        <li>
	  <button onClick={this.handleHomeClick} >Home </button>
        </li>
      	<li>
          <button onClick={this.handleHelpClick} >Help </button>
        </li>
        <li>
        </li>
	  </ul>
      </nav>
      </header>
     {actualComponent}
     <div id="treeVis">
	        <div id="tree">
	        </div>
            <div id="additionalInfo"></div>
     </div>
    </div>
  );
 }
}

class Menu extends React.Component {
  
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
                   query_ident: '80'}

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
    axios.post("api/phyloblast", formData)
     .then(function (response) {
         console.log(response.data);
         self.props.changeComp('data', response.data)
         self.props.changeComp('actual', 'phyloblast');
        
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


    return( 
        
        <div  >
        <p> Enter Blast search:</p>
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
        <br></br>
        <br></br>
	<p>  Taxonomy</p>
        <label>Select input type</label><br></br>
		    <select  id="tree_menu" name="tree_menu" onChange={this.taxaChoice}>	
			    <option value='0'>NCBI taxonomy</option>
			    <option value='1'>own taxonomic phylogeny</option>
		    </select>
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
         <br></br>
        <p>  Filter conditions</p>
         <label>E-value:</label>
	     <input type="number" id="eValue"  name='eValue' step='0.01' min='0' max='100' value={this.state.eValue} size='3' onChange={this.handleChange} /><br></br><br></br>
	     <label>min. query identity:</label>
	     <input type="number" id="query_ident" name="query_ident" step="1" min='50' max='100' value={this.state.query_ident} size='2' onChange={this.handleChange} /><br></br><br></br>
	     <label>min. query coverage:</label>
	     <input type="number" id="query_cover" name="query_cover" step="1" min='50' max='100' value={this.state.query_cover} size='2' onChange={this.handleChange} />
           <br></br>
            <br></br>
           <button type="submit" value="Submit" onClick={this.handleSubmit}>Submit </button>
           <br></br>
      </div>
    );
  }
}


class Phyloblast extends React.Component {
    constructor(props) {
    super(props);
    this.state = {actualData: null,
                  svgHeigth: 0}

    var tree = startTreevis(this.props.phyloData);
    console.log(tree)
    if (tree !== 0){
        var margin = ({top: 10 , right: 10, bottom: 40, left: 20});
        var barSVG = 300 ;
        const branchLength = 30;
        chart(tree, null, branchLength, barSVG, margin);

    }else{
      d3.select('#tree').append('div').text('Found 0 hits. Return to the main page and try another phylogentic tree');
    }

   }

   render(){
    /***
    var example_tree = "(((EELA:0.150276,CONGERA:0.213019):0.230956,(EELB:0.263487,CONGERB:0.202633):0.246917):0.094785,((CAVEFISH:0.451027,(GOLDFISH:0.340495,ZEBRAFISH:0.390163):0.220565):0.067778,((((((NSAM:0.008113,NARG:0.014065):0.052991,SPUN:0.061003,(SMIC:0.027806,SDIA:0.015298,SXAN:0.046873):0.046977):0.009822,(NAUR:0.081298,(SSPI:0.023876,STIE:0.013652):0.058179):0.091775):0.073346,(MVIO:0.012271,MBER:0.039798):0.178835):0.147992,((BFNKILLIFISH:0.317455,(ONIL:0.029217,XCAU:0.084388):0.201166):0.055908,THORNYHEAD:0.252481):0.061905):0.157214,LAMPFISH:0.717196,((SCABBARDA:0.189684,SCABBARDB:0.362015):0.282263,((VIPERFISH:0.318217,BLACKDRAGON:0.109912):0.123642,LOOSEJAW:0.397100):0.287152):0.140663):0.206729):0.222485,(COELACANTH:0.558103,((CLAWEDFROG:0.441842,SALAMANDER:0.299607):0.135307,((CHAMELEON:0.771665,((PIGEON:0.150909,CHICKEN:0.172733):0.082163,ZEBRAFINCH:0.099172):0.272338):0.014055,((BOVINE:0.167569,DOLPHIN:0.157450):0.104783,ELEPHANT:0.166557):0.367205):0.050892):0.114731):0.295021)";
    var tree = d3.layout.phylotree()
        .svg(d3.select('#tree'));

    tree(example_tree)
      .layout();
    console.log(tree)
    ***/
    //d3.select(svg).remove();
    //console.log(svg)
    return(
       <div id="vis">
          <div id="menu">
          <label>Show barchart for:</label>
           <select id="tree_menu" name="tree_menu" >
            <option value="2">none</option>
            <option value="0">node hits</option>
            <option value="1">subtree hits</option>
           </select>
           <br></br>
          <input id='phylogeny' type='checkbox' align='left' /><label>Calculate phylogeny</label>
        </div>
     </div>
    );

    }
}

function Wait(){
    return(
       <h1> Please wait. Data processing can take some minutes dependent on the workload of the Blast server. </h1>
    );

}

function Help() {
    return(
       <h1> Help </h1>
    );
}

export default App;
