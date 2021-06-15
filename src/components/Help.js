// used libraries
import React, {Component} from 'react';
import * as d3 from 'd3';
import phylogeny from '../visualisations/phylogeny_label.png';
import taxonomy from '../visualisations/taxonomic_map_label.png';

class Help extends Component {
    render() {

        return(
            <div id='helpInfo'>
                <div id='generalInfo'>
                    <h2 class='helpHeader'>About</h2>
                    <p> BLASTphylo is an interactive web-tool, which applies a BLASTp search for a given protein sequence and
                        maps the result on a given taxonomic tree. In addition, a phylogeny calculation of the remaining BLAST
                        hits can be performed and visualized.
                        For the calculations, local installations of
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
                <li>Selcet your BLAST type: BLASTn, BLASTp, or BLASTx</li>
                <li>Enter one/two sequence(s) or upload a FASTA file <br/>
                <div class='subList'><i>required format:</i><br />
                    <p style={{'margin-left':'20px'}}>
                        >query1 <br/>
                        AAATTTGMMM....<br />
                        ....CCCPLLLLLL<br />
                        >query2 <br />
                        TTTGCGPPP...<br />
                        ....AACCCMMLL<br />
                    </p>
                </div>
                </li>
                <li>Define your taxonomic tree:
                    <ul style={{margin:'5px'}} class='sublist'>
                        <li>Enter list of taxonomic IDs or scientific names</li>
                        <li><b>|subtree</b> select complete subtree </li>
                        <li><b>taxa|!(....)</b> remove all taxon in brackets from the subtree below the taxa <br />
                        <p style={{margin:'8px'}} > <i>example:</i> Staphylococcus|!(1280|subtree) = Staphylococcus except for the Staphylococcus aureus (txid: 1280) subtree </p>
                        </li>
                        <li><b>search bar:</b> Start typing to search for bacteria. Click on the listed bacteria to add them to the taxonomy</li>
                    </ul>
                </li>
                <li>Check setting of the filter parameter. They will influence the used sequences for the taxa-bassed phylogeny
                    <ul class='subList'>
                        <li><b>E-value: </b>expect value for random hit</li>
                        <li><b>alignment identity:  </b>minimal identity between query and subject sequencen within the alignment</li>
                        <li><b>query coverage:  </b>alignment has to cover at least <b>x</b>% of the query sequence</li>
                        <li><b>subject coverage:    </b>alignment has to cover at least <b>x</b>% of the HSP for the subject sequence</li>
                    </ul>
                </li>
                <li>Click on <b>Submit</b>. You will be redircted to the taxonomic mappind when the calculation is completed.</li>
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
            <p> Visualization of the taxonomic mapping for the BLAST result, which only includes taxa with minimal one hit
                or subtree hit. Hits are directly attached to BLAST hits, whereas subtree hits represent the sum of all hits
                below a taxon. Hover over the nodes to see the actual hit values and the taxonomic rank of the node <b>(5)</b>.
                For an easier comparison of the values, a <b>bar chart (3)</b> can be visualized.
                Furthermore, the visualized taxonomic level can be selected with the <b>collapse to (3)</b> option. The actual state
                of the <b> bar chart </b> and <b> collapse to </b> option is written next to the buttons <b>(4)</b>. <br/>
            </p>
            <img src={taxonomy} alt='taxonomyPicture' />
            <p> With the <b>phylogeny calculation (1)</b> you can start two different types of phylogeny calculations.
                One the one hand a taxa-based phylogeny and on the other hand unique sequence-based phylogeny. <br/>
                The taxa-based phylogeny extracts the best hit for each taxon, whereas the unique sequence-based phylogeny uses
                all hit (subject) sequences.
                Finally, you can export the visualization with the <b>export tree visualization (2)</b> option. The export
                includes a <b>publication-ready</b> option to remove large white spaces between nodes.
            </p>
        </div>
    );
}

function PhylogenyHelpInfo(){
    return(
        <div id='phylogenyHelpInfo'>
            <h2 class='helpHeader'>Phylogenetic Analysis</h2>
            <p> Visualizations of the phylogeny calculation based on the taxonomic mapping. For the MSA calculation
                MAFFT and for the tree construction, FastTree was used. <br/>
                BLASTphylo includes two different types of phylogeny calculations: taxa-based (shown below) and unique sequence-based.
                The taxa-based phylogeny extracted the best hit sequence for each taxon, whereas the unique sequence-based
                phylogeny, independent of the taxa, uses all BLAST hits.
            </p>
            <img src={phylogeny} alt='phylogenyPicture' />
            <p> Both phylogeny calculations can be visualized as <b>clade-based (3)</b> or <b>distance-based (4)</b> trees.
                The distance-based visualization allows the comparison of the phylogenetic distance between the sequences and the
                clade-based focus on the tree structure for an easier comparison of direct neighborhoods. <br/>
                In the taxa-based phylogeny, additional information <b>(1)</b> about the uniqueness and phylum level of
                the taxa can be visualized. For the unique sequence-based phylogeny, the node size represents the number of
                taxa which are attached to this hit. Similar to the taxonomic mapping, the tree can be exported as a jpeg or SVG file <b>(2)</b>.
            </p>
        </div>
    );
}


export default Help;