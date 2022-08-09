// used libraries
import React, {Component} from 'react';
import * as d3 from 'd3';
import phylogeny from '../visualisations/phylogeny.png';
import phylogeny2 from '../visualisations/phylogeny2.png';

//import taxonomy from '../visualisations/taxonomic_map_label.png';
import taxonomy from '../visualisations/taxonomy.png';

class Help extends Component {
    render() {

        return(
            <div id='helpInfo'>
                <div id='generalInfo'>
                    <h2 class='helpHeader'>About</h2>
                    <p> BLASTphylo is an interactive web-tool, which applies BLAST to a given protein/gene sequence,
                        maps the resulting hits onto a given taxonomic tree, and calculates a phylogeny, all in an automated pipeline.
                        Taxonomic and phylogenetic trees are visualized interactively and allow the user to manipulate the tree on demand, e.g. collapsing the taxonomy to a certain taxonomic rank.
                        If used locally installations of
                        <a class='menulinks' href='https://www.ncbi.nlm.nih.gov/books/NBK279690/' rel='nofollow'> BLAST</a>,
                        <a class='menulinks' href='https://mafft.cbrc.jp/alignment/software/index.html' rel='nofollow'> MAFFT </a>
                        and
                        <a class='menulinks' href='http://www.microbesonline.org/fasttree/' rel='nofollow'> FastTree </a> are required.
                    </p>
                </div>
                <br />
                <HowToDoHelpInfo />
                <br />
                <TaxonomyHelpInfo />
                <br />
                <PhylogenyHelpInfo />
            </div>
        );
    }
}

function HowToDoHelpInfo(){
    return(
        <div id='howtoHelpInfo'>
            <h2 class='helpHeader'>How to get started:</h2>
            <ol>
                <li>Select your BLAST type: BLASTn, BLASTp, or BLASTx</li>
                <li>Enter one/two sequence(s) or upload a FASTA file in the following format<br/>
                    <p style={{'margin-left':'20px'}}>
                        >query1 <br/>
                        AAATTTGMMM<br />
                        >query2 <br />
                        TTTGCGPPP<br />
                    </p>
                Alternatively, upload a pre-computed BLAST result in .csv format using the following columns:
                <div class='subList'>qacc sacc qstart qend sstart send slen nident evalue pident staxids qcovs sseq<br />
                </div>
                <p style={{color:"#fcb42d"}}>We highly recommend to use a pre-computed BLAST result as an input when using the web-version of BLASTphylo for performance reasons!</p>
                </li>
                <li>Define your taxonomic tree:
                    <ul style={{margin:'5px'}} class='sublist'>
                        <li>Enter a comma-separated list of taxonomic IDs or scientific names (Note: do not use spaces between list items)</li>
                        <li>The keyword <b>|subtree</b> selects the complete subtree rooted at the specified taxon <b>taxon|subtree</b></li>
                        <li>The keyword <b>|!(....)</b> removes all taxa in brackets from the subtree rooted at the taxon <b>taxon|!(....)</b><br />
                        <p style={{margin:'8px'}} > <i>example:</i> Staphylococcus|!(1280|subtree) = Staphylococcus except for the Staphylococcus aureus (txid: 1280) subtree </p>
                        </li>
                        <li>search bar: Start typing to search for bacteria. Select a taxon from the list to add them to the taxonomy.</li>
                    </ul>
                </li>
                <li>Check parameter used to filter the BLAST result:
                    <ul class='subList'>
                        <li><b>E-value: </b> number of expected hits of similar quality (score) that could be found just by chance.</li>
                        <li><b>alignment identity:  </b>minimal identity between query and subject sequence within the alignment</li>
                        <li><b>query coverage:  </b>alignment has to cover at least <b>x</b>% of the query sequence</li>
                        <li><b>subject coverage:    </b>alignment has to cover at least <b>x</b>% of the HSP for the subject sequence</li>
                    </ul>
                </li>
                <li>Click on <b>Submit</b>. You will be redirected to the taxonomic mapping as soon as the calculation is completed.</li>
            </ol>
        </div>
    );
}



function MenuHelpInfo(){
    return(
        <div id='menuHelpInfo'>
            <h2 class='helpHeader'>Menu</h2>
            <h4 class='helpHeader'>BLASTp search</h4>
            <p> BLASTphylo allows two different input formats for the BLASTp search: A <b>protein sequence</b> or an already
                calculated BLAST result as <b>.csv file</b>. <br/>
                The protein sequence can be entered in the text area or uploaded as a FASTA file. If the input was a nucleotide
                sequence it will be interpreted as protein sequence. <br/>
                Direct uploads of BLAST results have to fulfil specific properties for further calculations. These are that
                the BLAST result contains the following ordered columns: <b>qacc, sacc, qstart, qend, sstart, send, slen,
                nident, evalue, pident, staxids, qcovhsp, sseq</b>.
            </p>
            <h4 class='helpHeader'>Taxonomy</h4>
            <p> The taxonomic tree can be defined as a Newick string or as a regular expression based on the NCBI taxonomy. <br/>
                To enter a Newick string or upload a text file switch the <b>input type</b> to <i>'own taxonomic phylogeny'</i>.
                With this option, the given tree structure will be conserved during the calculation. <br/>
                For the regular expression of the NCBI taxonomy, you can enter a <b>comma-separated list</b> of taxonomic IDs or
                scientific names. With the addition of <b>|subtree</b> the hole subtree of this taxon can be selected. <br/>
                For example would the input <i>staphylococcus|subtree,bacillus|subtree</i> extract all taxa below the roots staphylococcus and bacillus.
            </p>
            <h4 class='helpHeader'>Filter parameter</h4>
            <p>
            Similar to the NCBI BLASTp search BLASTphylo allows the filtering of the result by a cut of value for <b>the E-value,
            the query and subject (hit) coverage, and the alignment identity</b>. The default values for these parameters are:
            0.05 (E-value), 80% (query coverage and alignment identity) and 50% (hit coverage).
            </p>
        </div>
    );
}

function TaxonomyHelpInfo(){
    return(
        <div id='taxonomyHelpInfo'>
            <h2 class='helpHeader'>Taxonomic Mapping</h2>
            <p> This visualization shows the subset of taxa of the specified taxonomy for which BLAST found at least one homologous sequence to the query sequence. It offers the following features:
                <ul>
                <li><b>(1)</b> Auto-collapse the tree to the specified taxonomic rank.</li>
                <li><b>(2)</b> Choose one of the following count information which are displayed in a barplot next to the tree:</li>
                <ul class='subList'>
                        <li>node hits: number of hits found for the respective leaf taxon</li>
                        <li>subtree hits: number of hits found within the subtree rooted at the respective leaf taxon</li>
                        <li>covered taxa: percentage of strains within the subtree rooted at the respective leaf taxon for which a homologous sequence was found</li>
                    </ul>
                <li><b>(3)</b> Toggle for an optimized layout.</li>
                <li><b>(4)</b> Hovering over nodes provides additional information.</li>
                </ul>
            </p>
            <p style={{"text-align":"center"}}>
            <img src={taxonomy} alt='taxonomyPicture' width="800" height='auto'/>
            </p>
        </div>
    );
}

function PhylogenyHelpInfo(){
    return(
        <div id='phylogenyHelpInfo'>
            <h2 class='helpHeader'>Phylogenetic Analysis</h2>
            <p> MAFFT is used to compute the multiple sequence alignment. The phylogenetic tree is calculated using FastTree.
                BLASTphylo includes two different types of phylogeny calculations: taxa-centered and unique sequence-centered.
                The taxa-centered phylogeny uses the best scoring hit sequence for each taxon, whereas the unique sequence-centered
                phylogeny uses all BLAST hits.
                <br/>
                <br/>
                The following visualization shows a taxa-centered phylogeny. Leaf nodes are labeled by taxa.
            </p>
            <p style={{"text-align":"center"}}>
            <img src={phylogeny2} alt='phylogenyPicture2' width="500" height='auto'/>
            </p>
                <br/>
                <p>
                The following visualization shows a sequence-centered phylogeny. Leaf nodes are labeled by protein ID.
                </p>
            <p style={{"text-align":"center"}}>
            <img src={phylogeny} alt='phylogenyPicture' width="500" height='auto'/>
            </p>
        </div>
    );
}


export default Help;