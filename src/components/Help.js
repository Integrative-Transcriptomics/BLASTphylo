// used libraries
import React, {Component} from 'react';
import * as d3v6 from 'd3v6';


class Help extends Component {

    render() {
        d3v6.select('#tree_vis').remove();
        return(
           <h1> Help </h1>
        );
    }
}

export default Help;