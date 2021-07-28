// used libraries
import React, {Component} from 'react';
import axios from 'axios';
import * as d3v6 from 'd3v6';
import {Spinner, Table} from 'react-bootstrap';


class HandlePhylogenyServer extends Component{
    constructor(props){
        super(props);
        console.log(props)
        d3v6.select('#tree_vis').remove();
        d3v6.select('#clade_vis').remove();
        d3v6.select('#visualisation').style('border', '2px solid white')
                                .style('border-radius', '5px');

    }

    // dependent on the given link information taxa-based or unique sequence-based phylogeny will calculated in the back end
    componentDidMount(){
        var self = this;
        var path = 'server/' + this.props.data;
        axios.post(path, null)
         .then(function (response) {
             //console.log(response.data);
             self.props.saveActualData(path, response.data);
         })
         .catch(error => {
            console.log('Error occurred')
            //console.log(error);
         })
   }

   // visualize a rough run time expectation
   render(){
        return(
            <div style={{margin:"20px"}}>
                <p>Calculation of the phylogeny can take up to <b>10 min </b> dependent
            on the number of hits, size of the taxonomic tree and GC content</p><br/>
                <Table>
                 <thead>
                    <tr>
                        <th># sequences</th><th>run time</th>
                    </tr>
                 </thead>
                 <tbody>
                    <tr>
                        <td>&lt; 50</td><td>&lt; 1 sec</td>
                    </tr>
                    <tr>
                        <td>50-250</td><td>&lt; 1 min</td>
                    </tr>
                    <tr>
                        <td>250-500</td><td>&lt; 2 min</td>
                    </tr>
                    <tr>
                        <td>&gt; 500</td><td>&gt; 3 min</td>
                    </tr>
                 </tbody>
                </Table>
                <div className="d-flex justify-content-center">
                    <Spinner animation="border" role="status">
                        <span className="sr-only">Loading...</span>
                    </Spinner>
                </div>
            </div>
        );
   }
}

export default HandlePhylogenyServer;