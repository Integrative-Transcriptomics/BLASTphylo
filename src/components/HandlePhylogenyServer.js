// used libraries
import React, {Component} from 'react';
import axios from 'axios';
import {Spinner, Table} from 'react-bootstrap';


class HandlePhylogenyServer extends Component{
    constructor(props){
        super(props);
        console.log(this.props)
        this.calculatePhylogeny = this.calculatePhylogeny.bind(this);

    }

    calculatePhylogeny(){
        var self = this;
        var path = 'server/' + this.props.data;
        axios.post(path, null)
         .then(function (response) {
             //console.log(response.data);
             self.props.changeComp('data', response.data);
             self.props.changeComp('previous', 'handlePhylogeny'); // set state for help page content
             self.props.changeComp('actual', 'phylogeny');
         })
         .catch(error => {
            console.log('Error occurred')
            //console.log(error);
         })
   }
   render(){
        this.calculatePhylogeny();
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