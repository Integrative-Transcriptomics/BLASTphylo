// used libraries
import React, {Component} from 'react';
import * as d3v6 from 'd3v6';
import "bootstrap/dist/css/bootstrap.min.css";
import {Table, Form, OverlayTrigger, Tooltip} from 'react-bootstrap';
import {IoMdArrowDropleft, IoMdArrowDropright} from 'react-icons/io';
import {BiHelpCircle} from 'react-icons/bi';
import {AiFillHome} from 'react-icons/ai';

//import "./visualisations/defaultFirefox.css";

// own components and style sheets
import './baseStyle.css';
import Menu from './components/Menu.js';
import Help from './components/Help.js';
import HandlePhylogenyServer from './components/HandlePhylogenyServer.js';
import TaxonomicAnalysisMenu from './components/TaxonomicAnalysisMenu.js';
import PhylogeneticAnalysisMenu from './components/PhylogeneticAnalysisMenu.js';
import TreeInteraction from './components/TreeInteraction.js';
import headerTree from './visualisations/header_tree.png';


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

    // switch between the componets and allow data transfer
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

        const isActualComponent = this.state.isActualComponent;   // main div: menu or visualizations
        const data = this.state.data;
        var userMenu = <div />;  // sub div: menus for the visualizations (export visualization, etc)

        //console.log(data)
        let actualComponent;
        d3v6.select('#phyloblastAlert').remove();

        //Navbar tooltips
        const generalTooltip = function(name){ return(<Tooltip id='tooltip-bottom'>
                                {name}
                            </Tooltip>)};


        // switch view based on the isActualComponent state
        if (isActualComponent === 'help'){
            d3v6.select('#visualisation').style('border', 'hidden');
            actualComponent = <Help />;
        } else if (isActualComponent === 'phyloblast'){ // taxonomic Mapping
            const copy = {...data};
            actualComponent = <TaxonomicAnalysisMenu phyloData={copy} changeComp={this.changeComp} />;
            userMenu = <TreeInteraction phyloData={copy} calculationMethod={'taxa'}/>;

        }else if (isActualComponent === 'handlePhylogeny'){ // handle data processing of the phylogenetic analysis
            actualComponent = <HandlePhylogenyServer data={data} changeComp={this.changeComp} />;
        }else if (isActualComponent === 'phylogeny'){ // visualize the phylogeny
            actualComponent = <PhylogeneticAnalysisMenu data={data} changeComp={this.changeComp} />;
            userMenu = <TreeInteraction data={data} calculationMethod={'phylo'} />;
        } else {
            d3v6.select('#visualisation').style('border', 'hidden');
            actualComponent = <Menu isActualComponent={isActualComponent} changeComp={this.changeComp} />;
        }
        //console.log(this.state.isActualComponent)

        if(this.state.error !== null){ // errors occurred during input
            alerts = this.state.error.map((alertmessage) =>
                <li><strong>Warning: </strong>{alertmessage.message}</li>);

            return (
                <div className="App">
                    <header className="App-header">
                    <nav >
                        <li id='headerName' >
                            <h1 id='title'>BLASTphylo</h1>
                            <h1 id='subtitle'>taxonomy <IoMdArrowDropleft /> blast result <IoMdArrowDropright /> phylogeny</h1>
                        </li>
                        <img src={headerTree} alt='headerPicture' id='headerPicture'/>
                        <ul  style={{width:'80%'}} id='subMenu'>
                            <li id='link1' >
                                    <button id='menuLink' onClick={this.handleMenuClick} >
                                        <OverlayTrigger key='tooltip_home' placement='bottom' overlay={generalTooltip('home')}>
                                        <AiFillHome style={{color: '#002060'}} size={25}/>
                                        </OverlayTrigger>
                                    </button>
                            </li>
                            <li id='link2'>
                                    <a id='helpLink' href='/help' target='_blank' >
                                        <OverlayTrigger key='tooltip_help' placement='bottom'  overlay={generalTooltip('help')}>
                                        <BiHelpCircle style={{color: '#002060'}} size={25}/>
                                        </OverlayTrigger>
                                    </a>
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
                    <nav >
                        <li id='headerName' >
                            <h1 id='title'>BLASTphylo</h1>
                            <h1 id='subtitle'>taxonomy <IoMdArrowDropleft /> blast result <IoMdArrowDropright /> phylogeny</h1>
                        </li>
                        <img src={headerTree} alt='headerPicture' id='headerPicture'/>
                        <ul  style={{width:'80%'}} id='subMenu'>
                            <li id='link1' >
                                    <button id='menuLink' onClick={this.handleMenuClick} >
                                        <OverlayTrigger key='tooltip_home' placement='bottom' overlay={generalTooltip('home')}>
                                        <AiFillHome style={{color: '#002060'}} size={25}/>
                                        </OverlayTrigger>
                                    </button>
                            </li>
                            <li id='link2'>
                                    <a id='helpLink' href='/help' target='_blank' >
                                        <OverlayTrigger key='tooltip_help' placement='bottom'  overlay={generalTooltip('help')}>
                                        <BiHelpCircle style={{color: '#002060'}} size={25}/>
                                        </OverlayTrigger>
                                    </a>
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



// check the actual browser the application is running in
//  based on: https://stackoverflow.com/questions/9847580/how-to-detect-safari-chrome-ie-firefox-and-opera-browser (8.6.2021)
function checkBrowser(){
    // Opera 8.0+
    //var isOpera = !!window.opr || !!window.opera || navigator.userAgent.indexOf(' OPR/') >= 0;

    // Firefox 1.0+
    //var isFirefox = typeof InstallTrigger !== 'undefined';

    // Safari 3.0+ "[object HTMLElementConstructor]"
    var isSafari = /constructor/i.test(window.HTMLElement) || (function (p) { return p.toString() === "[object SafariRemoteNotification]"; })(!window['safari'] || (typeof safari !== 'undefined' && window['safari'].pushNotification));

    // Internet Explorer 6-11
    //var isIE = /*@cc_on!@*/false || !!document.documentMode;

    // Edge 20+
    //var isEdge = !isIE && !!window.StyleMedia;

    // Chrome 1 - 79
    //var isChrome = !!window.chrome && (!!window.chrome.webstore || !!window.chrome.runtime);

    // Edge (based on chromium) detection
    //var isEdgeChromium = isChrome && (navigator.userAgent.indexOf("Edg") != -1);

    // Blink engine detection
    //var isBlink = (isChrome || isOpera) && !!window.CSS;

    if(isSafari){
        return 'safari'
    }else{
        return 'other'
    }
}


export {checkBrowser};
export default App;
