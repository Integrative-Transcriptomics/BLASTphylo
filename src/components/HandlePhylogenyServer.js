// used libraries
import React, {Component} from 'react';
import axios from 'axios';
import Spinner from 'react-bootstrap/Spinner';


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
            on the number of hits and size of the taxonomic tree</p>
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