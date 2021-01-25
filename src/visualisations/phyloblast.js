/***
# * PhyloBlast
# * @author jennifer mueller
# * 
# * based on: https://observablehq.com/@d3v6/collapsible-tree
              https://observablehq.com/@bumbeishvili/horizontal-collapsible-d3v6-flextree
# ***/

import * as d3v6 from 'd3v6';

// global variables
var clicked_nodes = {}
var treeVis = null;
var extraInfo = null;


// ****************** General Functions ***************************
// Function to remove "null" roots
function removeNullroot(tree){
  while (tree.name === null && tree['children'].length === 1 ){
	tree = tree['children'][0];
  }
  return tree;
}

// check if the dataset have any hits
function startTreevis(tree){
   if(Object.keys(tree).length === 0){
	return 0
   }else{
        return removeNullroot(tree);
   }
}


// ****************** Tooltip: show Nodevalues ***************************
const tooltip_x_pos = 10;
const tooltip_y_pos = 35;

function showTooltip(node, barhover, extraDataLength){
  var tooltip = d3v6.select("#tree")
	.append("div")
	.attr('id', 'tooltip')
	.style("position","absolute")
	.style("visibility","hidden")
	.style("border-radius", "10px")
	.style("padding", "5px")
	.style("color", "white")
	.style("background-color", "rgba(35, 29, 51, 0.7)");

  var text = null;
  if(extraDataLength === 1){
    text = '<b>' + node.data.name + '</b> <br /> taxa: ' + String(node.data.value[0]) + '<br />';
  }else if (extraDataLength === 2){
    text = '<b>' + node.data.name + '</b>';
  }else{
    text = '<b>' + node.data.name + '</b> <br /> hits: ' + String(node.data.value[0]) + '<br />' +
	 'subtree hits: ' + String(node.data.value[1]);
  }

  tooltip.style("visibility", "visible")
	    .html(text)
	    .style("left", (window.event.pageX + tooltip_x_pos) + "px")
	    .style("top", (window.event.pageY + tooltip_y_pos) + "px");

  // change styling of the circle to indicate on click action
  if(barhover){
	  if(node.children || node._children){
	 	 d3v6.select('#nodecircle' + node.id).attr('fill', 'crimson')
					    .attr('r', 5)
 	  }
  }
}

function moveTooltip(node, barhover){
  d3v6.select('#tooltip').style("left", (window.event.pageX + tooltip_x_pos) + "px")
         .style("top", (window.event.pageY + tooltip_y_pos) + "px");

}

function hideTooltip(node, barhover, oldCircleSize){
  d3v6.selectAll('#tooltip').remove();

  // change styling of the circle back to normal
  if(barhover){
	  if(node.children || node._children){
	  	d3v6.select('#nodecircle' + node.id).attr('r', oldCircleSize)
	 	 .attr("fill", function(d){ if(d.parent === null){return '#377ba8';} // violet: #9038c7
                                        else if(d._children || d.children){return "#555";}
                                        else{return "#999";}})
 	 }
  }

}


// ****************** Treevisualisation: basic interactive tree ***************************

// general variables
const flextree = require('d3-flextree').flextree;
const taxonomyLevel = ['life', 'domain', 'superkingdom', 'kingdom', 'clade', 'phylum', 'class', 'order', 'family', 'genus', 'species group','species', 'strain'];
var treeData = null;
var treeHeight = 0;
var svgWidth = 300;
var margin = ({top: 20 , right: 30, bottom: 40, left: 60});
var leftmostnode = null;
var rightmostnode = null;


// Collapse the node and all it's children
function collapse(d, branchLength) {
    if (d.children) {
        d._children = d.children;
        const newData = {'name':d.data.name, 'children':d.data.children, 'value':d.data.value,
        '_size':d.data.size, 'size':[15, branchLength]};
        d._children.forEach(collapse);
        d.data = newData
        d.children = null;

    }
    return d
}

// check if children have lower taxonomic level than in select
function checkTaxonomicLevel(children, taxonLevel){
    for (var i = 0; i < children.length; i++){
        if(taxonomyLevel.indexOf(children[i].data.value[2]) > taxonLevel){
            return true;
        }
    }
    return false;

}

// main visualisation of the tree
function chart(data, extraData, taxonomicLevel, previousTaxonomicLevel, onclickInteraction) {
    var branchLength = data['size'][1]
    const duration = 750;
    let i = 0;
    treeVis = data;
    extraInfo = extraData;

    if(document.getElementById('tree').getElementsByTagName('svg').length >= 1){
        d3v6.select('#tree_vis').remove();
    }

    const flexLayout = flextree();
    let root = d3v6.hierarchy(data, function(d) {
            return d.children;
    });
    var rootHeight =  root.data.size[0];

    if(extraData != null){ // data = phylogeny
        var maxTaxa = d3v6.max(root.descendants(), function(d){return d.data.value[0];});
        var circleSize = d3v6.scaleSqrt().domain([0, Math.log2(maxTaxa)]).range([2, 10]);
    }


    // define the maximal space for the barchart
    var max_name_length = 0
    root.descendants().forEach(function(d){
       var data_name = d.data.name;
       if((data_name.length) > max_name_length){ //4.96 == 0.31 em in pixel, 10px max height
          max_name_length = data_name.length;
       }

        // collapse tree
        if((d.data.name in clicked_nodes) && (clicked_nodes[d.data.name] === 0)){  // node was collapsed
            d = collapse(d, branchLength);
        }

        if (d.data.value[2]  && taxonomicLevel !== (taxonomyLevel.length-1)){
            if (((taxonomyLevel.indexOf(d.data.value[2]) >= taxonomicLevel) ||    //  nodes has lower rank  OR no rank AND not clicked
             (d.data.value[2] === 'no rank')) && !(d.data.name in clicked_nodes)){
                d = collapse(d, branchLength);
             }
             else if (d.children && !(d.data.name === root.data.name)){  // inner node and one of the children has lower rank
                var childWithLowerTaxon = checkTaxonomicLevel(d.children, taxonomicLevel);
                if(childWithLowerTaxon && !(d.data.name in clicked_nodes)){
                    d = collapse(d, branchLength);
                }
             }
        }

    });

    // publication ready tree visualisation have no onclickInteractions
    // --> resize the node to their actual child size
    if(!(onclickInteraction)){
        root.descendants().reverse().forEach(function(d){
            if(d.children){
                d.data.size[0] = d3v6.sum(d.children, function(child){return child.data.size[0];});
                d.data.size[0] = d.data.size[0] + 10;
            }
        });
    }

    //              max tree depth                 max label             spaces
    var width = (root.height * branchLength) +  (max_name_length*6) + 6 + margin.left ;
    treeHeight = rootHeight + margin.top + margin.bottom ;


    // append the svg object to the body of the page
    var svg = d3v6.select('#tree')
            .append('svg')
            .attr('id', 'tree_vis')
      		.style("font", "10px sans-serif")
     		.style("user-select", "none")

    var g = svg
        .attr("width", width)
        .attr("height", treeHeight)
        .append("g")
        .attr('id', 'group')
        .attr('transform',`translate(` + margin.left + `,${treeHeight/2})`)

    // Collapse after second level
    root.x0 = 0;
    root.y0 = 0;
  
    // center the scrolling bar to the root
    document.getElementById("treeVis").scrollTop = treeHeight/2-window.innerHeight/2;

    update(root);


    function update(source) {
         
        // Assigns the x and y position to the nodes
        treeData = flexLayout(root);

        leftmostnode = treeData;
	    rightmostnode = treeData;
	    treeData.eachBefore(node => {
	        if(node.x < leftmostnode.x) leftmostnode = node;
	        if(node.x > rightmostnode.x) rightmostnode = node;
	    });

        // Switch x and y coordinates for horizontal layout
        treeData.each(d=>{
          const x = d.x;
          d.x = d.y;
          d.y = x;
        })

        // Compute the new tree layout.
        var nodes = treeData.descendants(),
            links = treeData.descendants().slice(1);

        // ****************** Nodes section ***************************

        // Update the nodes...
        var node = g.selectAll('g.treenode')
            .data(nodes, d=> d.id || (d.id = ++i));

        // Enter any new modes at the parent's previous position.
        var nodeEnter = node.enter().append('g')
            .attr('class', 'treenode')
            .attr("transform", function(d) {
                return "translate(" + source.x0 + "," + source.y0 + ")";
            });

        // Add Circle for the nodes
        nodeEnter.append('circle')
            .attr('class', 'node')
	        .attr('id', function(d){return 'nodecircle' + d.id;})
            .attr('r', d => extraData != null && Object.keys(extraData).length === 1 ? circleSize(Math.log2(d.data.value[0])) : 4)
            .attr("fill", function(d,i){ if( i === 0){return '#377ba8';} // violet: #9038c7
                                        else if(d._children || d.children){return "#555";}
                                        else{return "#999";}})
            .attr("stroke-width", 10)
            .attr("cursor", "pointer")
            .on('click', click)
            .on("mouseover", function(d, i){ if (extraData != null){showTooltip(i, true, Object.keys(extraData).length);}
                                           else{ showTooltip(i, true, 0);}})
            .on("mouseout", function(d,i){ if (extraData != null && Object.keys(extraData).length === 1){hideTooltip(i, true, circleSize(Math.log2(i.data.value[0])));}
                                           else{ hideTooltip(i, true, 4);}})
            .on("mousemove", function(d, i){moveTooltip(i, true);});

        // Add labels for the nodes
        nodeEnter.append('text')
            .attr('id', function(d){ return 'nodetext' + d.id;})
            .attr('dy', '0.31em')
            .attr("x", d => d.children ? -6 : 12)
            .attr("text-anchor", d => d.children  ? "end" : "start")
            .attr('fill', d =>  d.children ? 'transparent' : 'black')
            .text(d =>  d.children ? "p" : d.data.name)
            .attr('cursor', 'pointer')
            .on('click', click)
            .on("mouseover", function(d, i){if (extraData != null){showTooltip(i, true, Object.keys(extraData).length);}
                                           else{ showTooltip(i, true, 0);}})
            .on("mouseout", function(d,i){ if (extraData != null && Object.keys(extraData).length === 1){hideTooltip(i, true, circleSize(Math.log2(i.data.value[0])));}
                                           else{ hideTooltip(i, true, 4);}})
            .on("mousemove", function(d, i){moveTooltip(i, true);});

        // UPDATE
        nodeEnter.merge(node)
            .transition()
            .duration(duration)
            .attr("transform", function(event, i, arr) { 
                const d = d3v6.select(this).datum();
                return "translate(" + d.x + "," + d.y + ")";
            });


        // Remove any exiting nodes
        node.exit().transition()
            .duration(duration)
            .attr("transform", function(event, i, arr) {
                return "translate(" + source.x + "," + source.y + ")";
            })
            .remove();



        // ****************** links section ***************************

        // Update the links...
        var link = g.selectAll('path.link')
            .data(links, function(d) {
                return d.id;
            });

        // Enter any new links at the parent's previous position.
        var linkEnter = link.enter().insert('path', "g")
            .attr("class", "link")
            .attr('d', function(d) {
                var o = {
                    x: source.x0,
                    y: source.y0
                }
                return diagonal(o, o)
            });

        // UPDATE
        var linkUpdate = linkEnter.merge(link)
            .attr("fill", "none")
            .attr("stroke", "#555")
            .attr("stroke-opacity", 0.4)
            .attr("stroke-width", 1.5)

        // Transition back to the parent element position
        linkUpdate.transition()
            .duration(duration)
            .attr('d', function(d) {
                return diagonal(d, d.parent)
            });

        // Remove any exiting links
        link.exit().transition()
            .duration(duration)
            .attr('d', function(event, i, arr) {
                var o = {
                    x: source.x,
                    y: source.y
                }
                return diagonal(o, o)
            })
            .remove();

        // Store the old positions for transition.
        nodes.forEach(function(d) {
            d.x0 = d.x;
            d.y0 = d.y;
        });
	

        // Creates a curved (diagonal) path from parent to the child nodes
        function diagonal(s, d) {
            const path = `M ${s.x} ${s.y}
            C ${(s.x + d.x) / 2} ${s.y},
              ${(s.x + d.x) / 2} ${d.y},
              ${d.x} ${d.y}`

            return path
        }

	function click(event, d) {
            if(onclickInteraction){
                if (d.children) {
                    clicked_nodes[d.data.name] = 0;  // node is collapsed
                    d3v6.select('#nodetext' + d.id).attr('fill', 'black')
                                .attr("x",  6)
                                .attr("text-anchor", "start")
                                .text(d.data.name);
                    d._children = d.children;
                    d.children = null;
                    d.data._size = d.data.size;
                    d.data.size = [15, branchLength];

                } else {
                    clicked_nodes[d.data.name] = 1;  // node is not collapsed
                    d3v6.select('#nodetext' + d.id).attr('fill', d =>  d.children || d._children ? 'transparent' : 'black')
                                             .text(d =>  d.children || d._children ? 'p' : d.data.name);
                    d.children = d._children;
                    d._children = null;
                    if(d.data._size){
                        var oldSize = d.data.size;
                        d.data.size = d.data._size;
                        d.data._size = oldSize;
                    }
                }
                update(d);


                if (extraData === null){   // taxonomic mapping
                    hitBars();
                } else if (extraData.length === 2) {    // phylogeny calculation (taxa-based)
                    showClades(extraData[0], extraData[1], null);
                }

            }

    }

    }

}




// ******************************** Barcharts   ***************************

// show hit values based on the actual state of the tree (phyloblast.html)
function hitBars(){
  var value = document.getElementById('tree_menu').value;
  var barheight = 7;

  var hitValue = null;
  if (value === '0'){  // show node bars
	hitValue = 0;
  } else if (value === '1'){
    hitValue = 1;
  } else {
         d3v6.select('#hitbars')
        .remove();
        return
  }

  d3v6.select('#hitbars')
        .remove();

  var nodes = treeData.descendants();
  var max_hit = d3v6.max(treeData.leaves(), function(d){return d.data.value[hitValue];});

  var scaleX = d3v6.scaleLinear().domain([0, max_hit]).range([0, svgWidth]);

  var hitbars = d3v6.select('#additionalInfo')
           .append('svg')
           .attr('id', 'hitbars')
           .attr("width", svgWidth+margin.left+margin.right)
           .attr("height", treeHeight)
           .style("font", "10px sans-serif")
           .style("user-select", "none");
      
      hitbars.append('g')
             .attr('transform', 'translate(' + margin.left + ',' + (treeHeight/2-barheight+leftmostnode.y) + ')')
             .call(d3v6.axisTop(scaleX));

      hitbars.append('g')
           .attr('transform', `translate(` + margin.left + `,${treeHeight/2})`)
           .selectAll('.bars')
           .data(nodes)
           .enter()
           .append('rect')
           .attr('class', 'bars')
           .attr('transform', function(d) { return 'translate(0,' + (d.y-4) + ')';})
           //.attr("x", 0 )
	   //.attr("y", function(d){return d.y;})
   	   .attr('fill',  d => d.children  ? "transparent" : "#377ba8")	
	   .attr("width", function(d){ return scaleX(d.data.value[hitValue]);})
	   .attr("height", barheight)
        .on("mouseover", function(d, i){showTooltip(i, false, null);})
        .on("mouseout", function(d, i){hideTooltip(i, false, 4);})
        .on("mousemove", function(d, i){moveTooltip(i, false);});

     if((rightmostnode.y-leftmostnode.y) >= (window.innerHeight*0.8)){
         hitbars.append('g')
                 .attr('transform', 'translate(' + margin.left + ',' + (treeHeight/2+rightmostnode.y+barheight) + ')')
                 .call(d3v6.axisBottom(scaleX));
     }
}

// translate node names into valid newick string names
function validNodename(taxa){
    var nonValidChar = [' ', '[', ']', '(', ')', ':', '/', '|', '\\'];
    var res = taxa;
    for (var i = 0; i < nonValidChar.length; i++){
        res = res.replaceAll(nonValidChar[i], '_');
    }

    return res;
}

// show parent taxa as color-coded bars  (phylogeny.html)
function showClades(taxData, accData, libTree){
   var checked = document.getElementById('clade_info').checked;
   
   // color-encoding of the taxids   Problem to much ids/accs as colors
   var taxids = Object.keys(taxData);
   var colorsParent = d3v6.scaleOrdinal().domain(taxids).range(d3v6.schemeCategory10);
  //console.log(taxids)

//[...new Set(d3v6.schemeDark2.concat(d3v6.schemePaired).concat(d3v6.schemeCategory10).concat(d3v6.schemeTableau10))])

   if(checked){
      //remove old taxids 
      d3v6.select('#clade_vis')
        .remove();
      
      var rectSize = 8;

      var nodes = treeData.descendants();
      var dataInfo = null;
      var actualTreeHeight = treeHeight/2;

      var phylotreePresent = (document.getElementsByClassName('phylotree-container').length === 1);
      if(phylotreePresent){
            var libTreenodes = libTree.get_nodes();
            dataInfo = {};                       // position of the nodes in distance tree
            for (var i = 0; i < libTreenodes.length; i++){
                dataInfo[libTreenodes[i].name] = libTreenodes[i];
            }
            var treeContainer = document.getElementsByClassName('phylotree-container')[0];
            var offset = treeContainer.transform.baseVal[0].matrix;
            actualTreeHeight = offset.f;
            //console.log(dataInfo)
      }
  
      var clades = d3v6.select('#additionalInfo')
           .append('svg')
           .attr('id', 'clade_vis')
           .attr("width", svgWidth+margin.left+margin.right)
           .attr("height", treeHeight)
           .style("font", "10px sans-serif")
           .style("user-select", "none");

      clades.append('g')
           .attr('transform', 'translate(' + margin.left + ',' + actualTreeHeight + ')')
           .selectAll('.accession')
           .data(nodes)
           .enter()
	       .append('rect')
           .attr('class', 'accession')
           .attr('transform',function(d,i) { var taxaName = validNodename(d.data.name);
                                if(phylotreePresent && dataInfo[taxaName] !== undefined){
                                return 'translate(0,' + (dataInfo[taxaName].screen_y-(rectSize/2)) + ')';}
                                else{return 'translate(0,' + (d.y-(rectSize/2)) + ')'; }})
           .attr('stroke', d => accData.includes(d.data.value[1]) ? 'black' : 'none')
   	       .attr('fill', function(d){if(accData.includes(d.data.value[1])) {
				return 'blue';}
                                else{
                                return 'transparent';}})
	       .attr("width", rectSize*2)
	       .attr("height", rectSize);

	  clades.append('g')
           .attr('transform', 'translate(' + margin.left + ',' + actualTreeHeight + ')')
           .selectAll('.clades')
           .data(nodes)
           .enter()
           .append('rect')
           .attr('class', 'clades')
           .attr('transform', function(d,i) {var taxaName = validNodename(d.data.name);
                                if(phylotreePresent && dataInfo[taxaName] !== undefined){
                                return 'translate(' + (rectSize*2+10)+ ',' + (dataInfo[taxaName].screen_y-(rectSize/2)) + ')';}
                                else{return 'translate(' + (rectSize*2+10)+ ',' + (d.y-(rectSize/2)) + ')'; }})
           .attr('stroke', d => d.data.value[0] === 0 || d.children ? 'none' : 'black')
   	       .attr('fill', function(d){if(d.data.value[0] === 0 || d.children) {
				return 'transparent';}
                                else{
                                return colorsParent(d.data.value[0]);}})
	       .attr("width", rectSize*2)
	       .attr("height", rectSize);

      var legend =  clades.selectAll('.legend')
                            .data(taxids)
                            .enter()
                            .append('g')
                            .attr('class', 'legend')
                            .attr('transform', 'translate(' + (8*rectSize) + ','+(treeHeight/2+leftmostnode.y) + ')');

	  legend.append('rect')
	        .attr('transform', function(d,i){return 'translate('+ (8*rectSize) + ',' + i*(rectSize+15) + ')';})
	        .attr('fill', function(d){return colorsParent(d);})
	        .attr('stroke', 'black')
	        .attr('width', rectSize*1.5)
	        .attr('height', rectSize*1.5)

	  legend.append('text')
	        .attr('transform', function(d,i){return 'translate('+ (8*rectSize) + ',' + i*(rectSize+15) + ')';})
	        .attr('x',rectSize*1.5+3)
	        .attr('y', rectSize*1.5/2)
	        .attr('dy', '0.31em')
            .attr("text-anchor", 'start')
            .attr('fill', 'black')
            .text(d => taxData[d]);

      var labels = ['unique', 'phylum'];

      clades.append('g')
            .attr('transform', 'translate(' + margin.left + ',' + (treeHeight/2-rectSize*3+leftmostnode.y) + ')')
            .selectAll('.labels')
            .data(labels)
            .enter()
            .append('text')
            .attr('class', 'labels')
            .attr('dy', '0.31em')
            .attr("text-anchor", "middle")
            .attr("dominant-baseline", "central")
            .attr('transform', function(d,i){return 'rotate(-45,' + i*(rectSize+15) + ', 0 )';})
            .attr('x', function(d,i){ return i*(rectSize+20);})
            .attr('y', 0)
            .text(function(d){return d;});



  }     
  else{
      d3v6.select('#clade_vis')
        .remove();

  }    
   
}

// collapse full tree to given taxonomic level
var previousTaxonomicLevel = null;
function collapseTree(tree){

    var taxonLevel = taxonomyLevel.indexOf(document.getElementById('collapse_menu').value);
    console.log(document.getElementById('collapse_menu').value)
    d3v6.select('#tree_vis').remove();
    chart(tree, null, taxonLevel, previousTaxonomicLevel, true);
    hitBars();
    previousTaxonomicLevel = taxonLevel;
}


// remove all withe spaces between nodes of the actual tree
function publicationReady(){
    d3v6.select('#tree_vis').remove();

    chart(treeVis, extraInfo, previousTaxonomicLevel, previousTaxonomicLevel, false);

    if (extraInfo === null){
        hitBars();
    } else if (document.getElementById('clade_info')){
        showClades(extraInfo[0], extraInfo[1], null);
    }
    d3v6.select('#public_ready').remove();
}

// ****************** Eventlistener ***************************
/***if (document.getElementById('tree_menu') != null){
    document.getElementById('tree_menu').addEventListener('change', hitBars);
} else if (document.getElementById('clade_info') != null){
    document.getElementById('clade_info').addEventListener('change', showClades);
}
***/

export {chart, hitBars, showClades,  startTreevis, collapseTree, publicationReady};
