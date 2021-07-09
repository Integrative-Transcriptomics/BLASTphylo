// used libraries
import React, {Component} from 'react';
import {CSVLink} from 'react-csv';
import axios from 'axios';
import * as d3v6 from 'd3v6';
import domtoimage from 'dom-to-image';

import {BsBoxArrowUpRight} from "react-icons/bs";
import {AiOutlineDownload} from 'react-icons/ai';



// export the view of the treeVis div (phlogeny or taxonomic mapping) as jpeg or SVG
class ExportTreeImage extends Component{
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
        return(<div id='treeImageExport'>
               <button id="export_svg" onClick={this.exportSVG}>export tree as svg <BsBoxArrowUpRight size={20}/> </button>
               <button id="export_pdf" onClick={this.exportJPEG}>export tree as jpeg <BsBoxArrowUpRight size={20}/> </button>
               </div>);
   }
}



// export data as CSV file
class ExportCsvData extends Component {
    constructor(props){
        super(props);
        this.state = {dataString: '',
                      downloadLink: '',
                      setDownloadLink: ''};

        // functions
        this.accessionOfData = this.accessionOfData.bind(this);
        this.makeTextFile = this.makeTextFile.bind(this);
    }

    // get data from server
    accessionOfData(){
        var path = 'server/exportData';
        const formData = new FormData();
        formData.append("datatype", this.props.filename);
        console.log(this.props)

        return axios.post(path, formData)
         .then(function (response) {
             console.log(response.data);
             return response.data

         })
         .catch(error => {
            console.log('Error occurred')
            //console.log(error);
         })

    }

    makeTextFile(dataString){
        //console.log(dataString)
        const data = new Blob([dataString], {type: 'text/plain'});

        if (this.state.downloadLink !== '') window.URL.revokeObjectURL(this.state.downloadLink);
        this.setState({downloadLink: window.URL.createObjectURL(data)});
    }

    componentDidMount(){
        this.accessionOfData().then(data => {
            //console.log(data);
            if(data.data_type === 'table'){
                this.setState({dataString: data['data']});
            }else if(data.data_type === 'newick'){
                this.makeTextFile(data['data']);
            }else{
                this.setState({dataString: 'No data present'})
            }
        });

    }


    render(){
        const exportedData = this.props.dataName;
        const filenameExportedData = this.props.filename;

        if(filenameExportedData.includes('.csv')){

            return(
                <button><CSVLink id='csvDataExport' data={this.state.dataString} filename={filenameExportedData} enclosingCharacter={'\''}>
                    download {exportedData} <AiOutlineDownload  size={20}/>
                </CSVLink></button>);
        }else{
            return(
                <button><a id='txtDataExport' download={filenameExportedData} href={this.state.downloadLink}>
                    download {exportedData} <AiOutlineDownload size={20}/>
                </a></button>);

        }
    }

}

export {ExportTreeImage, ExportCsvData};