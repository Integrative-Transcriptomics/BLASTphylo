// used libraries
import React, {Component} from 'react';
import * as d3v6 from 'd3v6';
import domtoimage from 'dom-to-image';
import "bootstrap/dist/css/bootstrap.min.css";
import {BsBoxArrowUpRight} from "react-icons/bs";
//import "./visualisations/defaultFirefox.css";

// own components and style sheets
import './baseStyle.css';
import Menu from './components/Menu.js';
import Help from './components/Help.js';
import HandlePhylogenyServer from './components/HandlePhylogenyServer.js';
import TaxonomicAnalysisMenu from './components/TaxonomicAnalysisMenu.js';
import PhylogeneticAnalysisMenu from './components/PhylogeneticAnalysisMenu.js';
import TreeInteraction from './components/TreeInteraction.js';

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

        // check if new Tab for phylogeny was "opened"
        var actualWebpage = window.location.href.split('/');
        if(actualWebpage[actualWebpage.length-1].includes('phylogeny')){
            this.state = {isActualComponent: 'handlePhylogeny',
                   data: actualWebpage[actualWebpage.length-1],
                   error: null};
        }else if(actualWebpage[actualWebpage.length-1] === 'help'){
            this.state = {isActualComponent: 'help',
                   data: null,
                   error: null};
        }else{
            this.state = {isActualComponent: 'menu',
                   data: null,
                   error: null};
        }

        // functions to handle view update, reset, help page and alert messages
        this.changeComp = this.changeComp.bind(this);
        this.handleMenuClick = this.handleMenuClick.bind(this);
        this.handleAlert = this.handleAlert.bind(this);
    }

    // change view independent of component
    changeComp(stateName, stateValue){
        if(stateName === 'actual'){
            this.setState({isActualComponent: stateValue});
        }else if (stateName === 'data'){
            this.setState({data: stateValue});
        } else if (stateName === 'error'){
            this.setState({error: stateValue});
        } else {
          console.log('Wrong stateName: ' + stateName + 'was given')
        }
    }

    // reset view to the menu page
    handleMenuClick(){
        if(window.location.href !== 'http://localhost:3000/'){
            window.location.href = 'http://localhost:3000/';
        }
        this.setState({isActualComponent: 'menu'});
    }



    // close alert messages
    handleAlert(){
        document.getElementById("alert").style.display='none';
    }



    render() {
        var alerts = null;

        const isActualComponent = this.state.isActualComponent;
        const data = this.state.data;
        var userMenu = <div />;

        //console.log(data)
        let actualComponent;
        d3v6.select('#phyloblastAlert').remove();

        // switch view dependent of the component state
        if (isActualComponent === 'help'){
            d3v6.select('#visualisation').style('border', 'hidden');
            actualComponent = <Help />;
        } else if (isActualComponent === 'phyloblast'){
            //this.getTaxonomicMapping();
            const copy = {...data};
            actualComponent = <TaxonomicAnalysisMenu phyloData={copy} changeComp={this.changeComp} />;
            userMenu = <TreeInteraction phyloData={copy} calculationMethod={'taxa'}/>;

        }else if (isActualComponent === 'handlePhylogeny'){
            actualComponent = <HandlePhylogenyServer data={data} changeComp={this.changeComp} />;
        }else if (isActualComponent === 'phylogeny'){
            actualComponent = <PhylogeneticAnalysisMenu data={data} changeComp={this.changeComp} />;
            userMenu = <TreeInteraction data={data} calculationMethod={'phylo'} />;
        } else {
            d3v6.select('#visualisation').style('border', 'hidden');
            actualComponent = <Menu isActualComponent={isActualComponent} changeComp={this.changeComp} />;
        }
        console.log(this.state.isActualComponent)

        if(this.state.error !== null){ // errors occurred during input
            alerts = this.state.error.map((alertmessage) =>
                <li><strong>Warning: </strong>{alertmessage.message}</li>);

            return (
                <div className="App">
                    <header className="App-header">
                    <nav>
                    <h1 id='title'>BLASTPhylo</h1>
                    <ul>
                        <li id='link1'>
                            <button id='menuLink' onClick={this.handleMenuClick} >menu </button>
                        </li>
                        <li id='link2'>
                            <a id='helpLink' href='/help' target='_blank' >help </a>
                        </li>
                    </ul>
                    </nav>
                    </header>
                    <div className='App-body'>
                        <div id="alert">
                            <span class="closebtn" onClick={this.handleAlert}>&times;</span>
                            <ul>
                                {alerts}
                            </ul>
                        </div>
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
        }else{ // calculation can run with success
            return (
                <div className="App">
                    <header className="App-header">
                    <nav>
                    <h1 id='title'>BLASTPhylo</h1>
                    <ul>
                        <li id='link1'>
                            <button id='menuLink' onClick={this.handleMenuClick} >menu </button>
                        </li>
                        <li id='link2'>
                            <a id='helpLink' href='/help' target='_blank' >help </a>
                        </li>
                    </ul>
                    </nav>
                    </header>
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
}

// export the view of the treeVis div (phlogeny or taxonomic mapping) as jpeg
class ExportTrees extends Component{
    constructor(props) {
        super(props);
        this.exportJPEG = this.exportJPEG.bind(this);
        this.exportSVG = this.exportSVG.bind(this);
    }



    exportJPEG(){
        var element = document.getElementById('treeVis');
        console.log(element)
        var borderStyle = element.style.border;
        element.style.border = 'none';
        // labeling of the tree
        var figureName = '';
        if(window.location.href.includes('phylogeny')){
            figureName = 'phylogeny';
        }else{
            figureName = 'taxonomicmapping';
        }

        domtoimage.toJpeg(element, { quality: 1, bgcolor: 'white',
                                            style:{overflow:'visible'} })
            .then(function (dataUrl) {
                var link = document.createElement('a');
                link.download = figureName;
                link.href = dataUrl;
                link.click();
                link.remove();
                element.style.border = borderStyle;
        });
        element.style.height = '80vh';
    }

        exportSVG(){
        var element = document.getElementById('treeVis');
        console.log(element)
        var borderStyle = element.style.border;
        element.style.border = 'none';
        // labeling of the tree
        var figureName = '';
        if(window.location.href.includes('phylogeny')){
            figureName = 'phylogeny';
        }else{
            figureName = 'taxonomicmapping';
        }

        domtoimage.toSvg(element, { quality: 1, bgcolor: 'white',
                                            style:{overflow:'visible'} })
            .then(function (dataUrl) {
                var link = document.createElement('a');
                link.download = figureName;
                link.href = dataUrl;
                link.click();
                link.remove();
                element.style.border = borderStyle;
        });
        element.style.height = '80vh';
    }

   render(){
        return(<div>
               <button id="export_svg" onClick={this.exportSVG}>export tree as svg <BsBoxArrowUpRight size={20}/> </button>
               <button id="export_pdf" onClick={this.exportJPEG}>export tree as jpeg <BsBoxArrowUpRight size={20}/> </button>
               </div>);
   }
}


export {ExportTrees};
export default App;
