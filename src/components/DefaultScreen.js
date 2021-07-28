// used libraries
import React, {Component} from 'react';
import * as d3v6 from 'd3v6';
import "bootstrap/dist/css/bootstrap.min.css";
import {Table, Form, OverlayTrigger, Tooltip} from 'react-bootstrap';
import {IoMdArrowDropleft, IoMdArrowDropright} from 'react-icons/io';
import {BiHelpCircle} from 'react-icons/bi';
import {AiFillHome} from 'react-icons/ai';


// own components and style sheets
import Menu from './Menu.js';
import Help from './Help.js';

class DefaultScreen extends Component {

    constructor(props){
        super(props);
        console.log(props);
        this.state = {error: null}


        // functions to handle view update, reset, help page and alert messages
        this.handleAlert = this.handleAlert.bind(this);
        this.handleErrorMessages = this.handleErrorMessages.bind(this);
    }

    // set given alert messages from Menuinput
    handleErrorMessages(messages){
        console.log(messages);
        this.setState({error: messages});
    }


    // close alert messages
    handleAlert(){
        document.getElementById("alert").style.display='none';
    }

    render() {
        var alerts = null;
        var errormessage = <div />;
        d3v6.select('#phyloblastAlert').remove();

        const actualComponent = <Menu changeComp={this.props.changeComp} sendErrors={this.handleErrorMessages}/>;


        if(this.state.error !== null){ // errors occurred during input
            alerts = this.state.error.map((alertmessage) =>
            <li><strong>Warning: </strong>{alertmessage.message}</li>);

            errormessage = <div style={{color: 'red'}}><b>Error occurred. Check error messages above.</b></div>;

            return (

                    <div className='App-body'>
                        <div id="alert">
                            <span class="closebtn" onClick={this.handleAlert}>&times;</span>
                            <ul>
                                {alerts}
                            </ul>
                        </div>
                        {actualComponent}
                        {errormessage}
                    </div>

            );
        }else{ // calculation can run with success
            return (
                    <div className='App-body'>
                        {actualComponent}
                    </div>

            );
        }
    }
}

export default DefaultScreen;
