// used libraries
import React, {Component} from 'react';
import html2canvas from 'html2canvas';



// own components and style sheets
import './baseStyle.css';
import Menu from './components/Menu.js';
import Phyloblast from './components/Phyloblast.js';
import Phylogeny from './components/Phylogeny.js';
import Help from './components/Help.js';

class App extends Component {
  
  constructor(props){
    super(props);
    //console.log(this.props)
    this.changeComp = this.changeComp.bind(this);
    this.handleHomeClick = this.handleHomeClick.bind(this);
    this.handleHelpClick = this.handleHelpClick.bind(this);
    this.handleAlert = this.handleAlert.bind(this);

    this.state = {isActualComponent: 'menu',
                   data: null,
                   error: null};
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
    this.setState({isActualComponent: 'menu'});
  }
  
  handleHelpClick(){
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
   if (isActualComponent === 'help'){
     actualComponent = <Help />;
   }
   else if (isActualComponent === 'wait'){
     actualComponent = <Wait />;
   }
   else if (isActualComponent === 'phyloblast'){
     const copy = {...data};
     actualComponent = <Phyloblast phyloData={copy} changeComp={this.changeComp}/>;

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
        // change size to full treesize + bars/clade info
        var element = document.getElementById('treeVis');
        var styleOrig = element.getBoundingClientRect();
        var treefigure = document.getElementById('tree_vis').getBoundingClientRect();
        element.style.width = treefigure.width+350;
        element.style.height = treefigure.height;

        html2canvas(element).then(function (canvas) {
            javascript:void(window.open().location = canvas.toDataURL("image/jpeg", 5.0));
        });


        // resize window to original size
        element.style.width = styleOrig.width;
        element.style.height = styleOrig.height;
    }

   render(){
        return(<button id="export_pdf" onClick={this.exportJPEG}>export tree as jpeg</button>);
   }
}


export {ExportTrees};
export default App;
