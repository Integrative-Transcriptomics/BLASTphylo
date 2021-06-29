// used libraries
import React, {Component} from 'react';
import axios from 'axios';
import SearchField from 'react-search-field';
import {Dropdown, Form, Spinner} from 'react-bootstrap';
import * as d3v6 from 'd3v6';



class SearchBar extends Component{
    constructor(props){
        super(props);
        this.state = {taxon: '',
                      matching_taxa: [],
                      htmlSize: 0}


        //own functions
        this.searchBacteria = this.searchBacteria.bind(this)
        this.handleInputChange = this.handleInputChange.bind(this)
        this.handleSelection = this.handleSelection.bind(this)
    }

    // access server and search for matching pattern
    searchBacteria(event){
        var self = this;
        var path = 'server/searchNcbiTaxa'
        document.getElementById('searchSpinner').style.opacity = 1;

        if(event !== ''){
            const formData = new FormData();
            formData.append('searchquery', event)

            axios.post(path, formData)
                .then(function (response){
                    console.log(response);
                    self.setState({matching_taxa:response.data.result});
                    self.setState({htmlSize: 10});
                })
                .catch(error => {
                    console.log(error);
                })
        }
        document.getElementById('searchSpinner').style.opacity = 0;


    }

    handleInputChange(value){
        //console.log(value)
        this.setState({taxon: value});
    }

    handleSelection(event){
        console.log(event);

        var oldValue = document.getElementById('taxa').value;
        document.getElementById('taxa').value = oldValue.concat(event.target.value).concat(',');

    }

    // visualize a rough run time expectation
    render(){
        var matching_taxa = this.state.matching_taxa;

        var MakeItem = function(X){
            return <option value={X['menuInput']}>{X['name']}</option>
        };

        if (matching_taxa.length !== 0){
            document.getElementById('SearchResult').style.opacity = 1;
        }

        return(
            <div id='Searchbar'>
                <div id='SearchbarContent'>
                <div id='SearchField'>
                        <SearchField
                            placeholder='Search for bacteria'
                            onEnter={this.searchBacteria}
                            onChange={this.handleInputChange}
                            searchText={this.state.taxon}
                       />
                        <Spinner
                            id='searchSpinner'
                            as="span"
                            animation="border"
                            size="sm"
                            role="status"
                            aria-hidden="true"
                        />
               </div>
               <Form.Control as='select' htmlSize={this.state.htmlSize} id='SearchResult' onClick={this.handleSelection} >
                    {matching_taxa.map(MakeItem)}
               </Form.Control>
               </div>
            </div>
        );
   }
}

export default SearchBar;