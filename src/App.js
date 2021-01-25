// used libraries
import React, {Component} from 'react';
import * as d3v6 from 'd3v6';
import domtoimage from 'dom-to-image';
import "bootstrap/dist/css/bootstrap.min.css";
//import "./visualisations/defaultFirefox.css";

// own components and style sheets
import './baseStyle.css';
import Menu from './components/Menu.js';
import Help from './components/Help.js';
import Phyloblast from './components/Phyloblast.js';
import Phylogeny from './components/Phylogeny.js';
import HandlePhylogenyServer from './components/HandlePhylogenyServer.js';


/*** #################################### IMPORTANT
    1.
        update of window.location.href for home/help links has to be updated if tool become online version
    2.
        at the moment tool is adapted to linux like file systems most likely don't work with Windows
    3.
        export of figures is limited by node number especially unique sequence based phylogeny

***/

class App extends Component {

  constructor(props){
    super(props);
    //console.log(this.props)

    // check if new Tab for phylogeny was "opened"
    var actualWebsite = window.location.href.split('/');
    if(actualWebsite[actualWebsite.length-1].length === 0){
        this.state = {isActualComponent: 'menu',
                   data: null,
                   error: null};
    }else{
        this.state = {isActualComponent: 'handlePhylogeny',
                   data: actualWebsite[actualWebsite.length-1],
                   error: null};
    }

    this.changeComp = this.changeComp.bind(this);
    this.handleHomeClick = this.handleHomeClick.bind(this);
    this.handleHelpClick = this.handleHelpClick.bind(this);
    this.handleAlert = this.handleAlert.bind(this);
  }	

  changeComp(stateName, stateValue){
    if(stateName === 'actual'){
    	this.setState({isActualComponent: stateValue});
    } else if (stateName === 'data'){
        this.setState({data: stateValue});
    } else if (stateName === 'error'){
        this.setState({error: stateValue});
    } else {
      console.log('Wrong stateName: ' + stateName + 'was given') 
    }
  }
      
  handleHomeClick(){
    window.location.href = 'http://localhost:3000/';
    this.setState({isActualComponent: 'menu'});
  }
  
  handleHelpClick(){
    window.location.href = 'http://localhost:3000/';
    this.setState({isActualComponent: 'help'});	
  }

  handleAlert(){
    document.getElementById("alert").style.display='none';
  }

  render() {
   var alerts = null;
   const isActualComponent = this.state.isActualComponent;
   const data = this.state.data;

   //console.log(data)
   let actualComponent;
   d3v6.select('#treeVis').style('border', 'none');
   d3v6.select('#phyloblastAlert').remove();

   if (isActualComponent === 'help'){
     actualComponent = <Help />;
   }
   else if (isActualComponent === 'wait'){
     actualComponent = <Wait />;
   }
   else if (isActualComponent === 'phyloblast'){
     const copy = {...data};
     actualComponent = <Phyloblast phyloData={copy} changeComp={this.changeComp}/>;

   }else if (isActualComponent === 'handlePhylogeny'){
     actualComponent = <HandlePhylogenyServer data={data} changeComp={this.changeComp} />;
   }else if (isActualComponent === 'phylogeny'){
     actualComponent = <Phylogeny data={data} changeComp={this.changeComp} />;
   } else {
     actualComponent = <Menu isActualComponent={isActualComponent} changeComp={this.changeComp}/>;
   }
   console.log(this.state.isActualComponent)

   if(this.state.error !== null){
      alerts = this.state.error.map((alertmessage) =>
                <li><strong>Warning: </strong>{alertmessage.message}</li>);

      return (
        <div className="App">
          <header className="App-header">
          <nav>
          <h1 id='title'>BLASTPhylo</h1>
          <ul>
            <li id='link1'>
               <button id='menuLink' onClick={this.handleHomeClick} >home </button>
            </li>
            <li id='link2'>
              <button id='helpLink' onClick={this.handleHelpClick} >help </button>
            </li>
          </ul>
          </nav>
          </header>
          <div id="alert">
              <span class="closebtn" onClick={this.handleAlert}>&times;</span>
              <ul>
                {alerts}
              </ul>
          </div>
         {actualComponent}
         <br></br>
         <br></br>
         <div id="treeVis">
                <div id="tree">
                </div>
                <div id="additionalInfo"></div>
         </div>
        </div>
      );
  }else{
      return (
    <div className="App">
      <header className="App-header">
      <nav>
	  <h1 id='title'>BLASTPhylo</h1>
 	  <ul>
        <li id='link1'>
	       <button id='menuLink' onClick={this.handleHomeClick} >home </button>
        </li>
      	<li id='link2'>
          <button id='helpLink' onClick={this.handleHelpClick} >help </button>
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
}


function Wait(){
    return(
       <h1> Please wait. Data processing can take some minutes dependent on the workload of the Blast server. </h1>
    );

}



class ExportTrees extends Component{
    constructor(props) {
        super(props);
        this.exportJPEG = this.exportJPEG.bind(this);
    }



    exportJPEG(){
        var element = document.getElementById('treeVis');
        var borderStyle = element.style.border;
        element.style.border = 'none';

        // labeling of the tree
        var figureName = '';
        if(borderStyle.includes('#69a2c9')){
            figureName = 'taxamapping.jpeg';
        }else{
            figureName = 'phylogeny.jpeg';
        }

        domtoimage.toJpeg(element, { quality: 0.95, bgcolor: 'white',
                                            style:{overflow:'visible'} })
            .then(function (dataUrl) {
                var link = document.createElement('a');
                link.download = figureName;
                link.href = dataUrl;
                link.click();
                link.remove();
                element.style.border = borderStyle;
        });


    }

   render(){
        return(<button id="export_pdf" onClick={this.exportJPEG}>export tree as jpeg</button>);
   }
}


export {ExportTrees};
export default App;
