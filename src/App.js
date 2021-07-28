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
import DefaultScreen from './components/DefaultScreen.js';
import Tabhandling from './components/Tabhandling.js';
import HandlePhylogenyServer from './components/HandlePhylogenyServer.js';
import TaxonomicAnalysisMenu from './components/TaxonomicAnalysisMenu.js';
import PhylogeneticAnalysisMenu from './components/PhylogeneticAnalysisMenu.js';
import TreeInteraction from './components/TreeInteraction.js';
import headerTree from './visualisations/header_tree.png';



class App extends Component {

    constructor(props){
        super(props);


        this.state = {isActualComponent: 'defaultScreen',
                   data: null,
                   error: null};


        // functions to handle view update, reset, help page and alert messages
        this.changeComp = this.changeComp.bind(this);
        this.handleMenuClick = this.handleMenuClick.bind(this);
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
        this.setState({isActualComponent: 'menu'});
    }



    render() {

        const isActualComponent = this.state.isActualComponent;   // main div: menu or visualizations
        const data = this.state.data;

        //console.log(data)
        let actualComponent;

        //Navbar tooltips
        const generalTooltip = function(name){ return(<Tooltip id='tooltip-bottom'>
                                {name}
                            </Tooltip>)};


        // switch view based on the isActualComponent state
        if (isActualComponent === 'help'){
            actualComponent = <Help />;
        } else if (isActualComponent === 'tabhandling'){ // taxonomic Mapping
            actualComponent = <Tabhandling data={data} changeComp={this.changeComp} />;
        } else {
            actualComponent = <DefaultScreen isActualComponent={isActualComponent} changeComp={this.changeComp} />;
        }
        //console.log(this.state.isActualComponent)


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
                                    <a id='helpLink' href='/help' >
                                        <OverlayTrigger key='tooltip_help' placement='bottom'  overlay={generalTooltip('help')}>
                                        <BiHelpCircle style={{color: '#002060'}} size={25}/>
                                        </OverlayTrigger>
                                    </a>
                            </li>
                        </ul>
                    </nav>
                    </header>
                    {actualComponent}
                </div>
        );
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
