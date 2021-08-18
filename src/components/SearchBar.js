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
                      htmlSize: 0,
                      spinnerStyle: 0}


        //own functions
        this.searchBacteria = this.searchBacteria.bind(this)
        this.handleInputChange = this.handleInputChange.bind(this)
        this.handleSelection = this.handleSelection.bind(this)
    }

    // access server and search for matching pattern
    searchBacteria(event){
        this.setState({spinnerStyle: 1});
        console.log(event)
        var self = this;
        var path = 'server/searchNcbiTaxa'

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
        this.setState({spinnerStyle: 0});


    }

    handleInputChange(value){
        //console.log(value)
        this.setState({taxon: value});
    }

    handleSelection(event){
        console.log(event);


        if(this.state.taxon.length > 0){
            var oldValue = document.getElementById('taxa').value;
            document.getElementById('taxa').value = oldValue.concat(event.target.value).concat(',');
        }
    }

    // visualize a rough run time expectation
    render(){
        var spinnerStyle = this.state.spinnerStyle;
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
                            onSearchClick={this.searchBacteria}
                            onChange={this.handleInputChange}
                            searchText={this.state.taxon}
                       />
                       <Spinner style={{opacity:spinnerStyle}}
                            id='searchSpinner'
                            animation="border"
                            size="sm"
                            role="status"
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