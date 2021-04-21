/***
# * PhyloBlast
# * @author jennifer mueller
# * 
# * based on: https://observablehq.com/@d3v6/collapsible-tree
              https://observablehq.com/@bumbeishvili/horizontal-collapsible-d3v6-flextree
# ***/

import * as d3v6 from 'd3v6';
import "../../node_modules/phylotree/src/main";

// global variables
var clicked_nodes = {};
var treeVis = null;
var extraInfo = null;
var ownAddInfo = null;
const flextree = require('d3-flextree').flextree;
const taxonomyLevel = ['life', 'domain', 'superkingdom', 'kingdom', 'clade', 'phylum', 'class', 'order', 'family', 'genus', 'species group','species', 'strain'];
var treeData = null;
var treeHeight = 0;
var svgWidth = 300;
var margin = ({top: 60 , right: 50, bottom: 40, left: 60});
var leftmostnode = null;
var rightmostnode = null;
var hitSelection = '-';
var phyloTreetree = null;


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
const tooltip_x_pos = 20;
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

    var text = null;   // define the content of the tooltip
    if(extraDataLength === 1){
        text = '<b>' + node.data.name + '</b> <br /> taxa: ' + String(node.data.value[0]) + '<br />';
    }else if (extraDataLength === 2){
        text = '<b>' + node.data.name + '</b> <br /> nodes: ' + String(node.data.value[2]) + '<br />';
    }else if (node.key){
        text = '<b>' + node.key + '</b> <br /> nodes: ' + String(node) + '<br />';
    }else{
        if(taxonomyLevel.includes(node.data.value[2]) || node.data.value[2] === 'no rank'){
            text = '<b>' + node.data.name + '</b> <br /> query1 hits: ' + String(node.data.value[0]) + '<br />' +
	                'query2 hits: ' + String(node.data.value[1]) + '<br /> rank: ' + node.data.value[2] + '<br />';
	    }else{
	        text = '<b>' + node.data.name + '</b> <br /> hits: ' + String(node.data.value[0][0]) + '<br />' +
	                'subtree hits: ' + String(node.data.value[0][1]) + '<br /> rank: ' + node.data.value[1] + '<br />';
	    }
    }

    tooltip.style("visibility", "visible")
            .html(text)
            .style("left", (window.event.pageX + tooltip_x_pos) + "px")
            .style("top", (window.event.pageY - tooltip_y_pos) + "px");

    // change styling of the circle to indicate on click action
    if(barhover){
	    if(node.children || node._children){
	 	    d3v6.select('#nodecircle' + node.id).attr('fill', 'crimson')
			                                    .attr('r', 5);
 	    }
    }
}

function moveTooltip(node, barhover){
    d3v6.select('#tooltip')
            .style("left", (window.event.pageX + tooltip_x_pos) + "px")
            .style("top", (window.event.pageY - tooltip_y_pos) + "px");

}

function hideTooltip(node, barhover, oldCircleSize){
    d3v6.selectAll('#tooltip').remove();

    // change styling of the circle back to normal
    if(barhover){
	    if(node.children || node._children){
	  	    d3v6.select('#nodecircle' + node.id).attr('r', oldCircleSize)
	 	        .attr("fill", function(d){ if(d.parent === null){return '#377ba8';}
                                           else if(d._children || d.children){return "#555";}
                                           else{return "#999";}});
 	    }
    }
}


// ****************** Treevisualization: basic interactive tree + necessary functions ***************************
// Collapse the node and all it's children
function collapse(d, b, markChilds) {
    if (d.children) {
        d._children = d.children;
        d._children.forEach(function(child){collapse(child, b, markChilds);});
        d.children = null;

        if(markChilds){
            clicked_nodes[d.data.name] = 2;
        }
    }
    return d
}

// collapse the nodes and all it's children after a 'collapse to' interaction
function collapseClickedNodesAfterSelection(d, b){
    if (d.children) {
        d._children = d.children;
        d._children.forEach(function(child){collapse(child, b, true);});
        d.children = null;
    }
    return d
}

// expand the nodes and all it's children after a 'collapse to' interaction if their taxonomic rank is higher
function expand(d, b){
    if (d._children) {
        d.children = d._children;
        d._children.forEach(function(child){expand(child, b);});
        d._children = null;
        delete clicked_nodes[d.data.name];
    }
    return d
}

// check if children have lower taxonomic level than in select
/***function checkTaxonomicLevel(children, taxonLevel){
    for (var i = 0; i < children.length; i++){
        if(taxonomyLevel.indexOf(children[i].data.value[taxaIndex]) > taxonLevel){
            return true;
        }
    }
    return false;
}***/




// main visualisation of the tree
function chart(data, extraData, taxonomicLevel, previousTaxonomicLevel, onclickInteraction, returnFromStatic) {
    //console.log(clicked_nodes)

    var branchLength = data['size'][1];
    const duration = 750;
    let i = 0;
    if(Object.keys(data).length > 1){
        treeVis = data;
    }
    extraInfo = extraData;
    var taxaIndex = 1;
    if(taxonomyLevel.includes(treeVis['value'][2])){
        taxaIndex = 2;
    }
    if(document.getElementById('tree').getElementsByTagName('svg').length >= 1){
        d3v6.select('#tree_vis').remove();
    }

    const flexLayout = flextree();
    let root = d3v6.hierarchy(treeVis, function(d) {
            return d.children;
    });
    //console.log(root)

    if((root.data.name in clicked_nodes) && (clicked_nodes[root.data.name] === 0)){  // root was collapsed before collapse to interaction
        delete clicked_nodes[root.data.name];
    }
    var rootHeight = root.data.size[0];

    if(returnFromStatic && root.data._size){
        console.log('resize root')
        root.data.size = root.data._size.slice();
        delete root.data._size;
    }

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

        // resize tree if back (<--) button was clicked
        if(returnFromStatic && d.data._size){
            //console.log(d.data)
            //var oldSize = d.data.size.slice();
            d.data.size = d.data._size.slice();
            //d.data._size = oldSize;
            //console.log(d.data)
            delete d.data._size;
        }


        // collapse tree
        var indexOfNode = taxonomyLevel.length-1;
        if(d.data.value[taxaIndex]){
            indexOfNode = taxonomyLevel.indexOf(d.data.value[taxaIndex]);
        }

        if((data_name in clicked_nodes)) { // node was collapsed
        //(!(isNaN(d.data.value[taxaIndex])) && (d.data.value[taxaIndex] <= 30) )){                 // node contain less than 30 children
            if ((clicked_nodes[data_name] === 2) && (indexOfNode < taxonomicLevel)
                && (d.data.value[taxaIndex] !== 'no rank')){
                //console.log(['expand', data_name])
                //console.log(d.data.value[taxaIndex])
                d = expand(d, branchLength);
            }else if((clicked_nodes[data_name] === 2) && taxonomicLevel !== (taxonomyLevel.length-1)){
                //console.log(['col1', data_name])
                //console.log(d.data.value[taxaIndex])
                d = collapse(d, branchLength, false);
            }else if (clicked_nodes[data_name] === 0){
                d = collapseClickedNodesAfterSelection(d, branchLength);
                //console.log(['col2', data_name])
               // console.log(d.data.value[taxaIndex])
            }
        }

        if(!(extraInfo)){ // only for taxonomic mapping
            if (d.data.value[taxaIndex]  && (taxonomicLevel !== (taxonomyLevel.length-1)) && // strain level --> all nodes are visible
                !(data_name in clicked_nodes)){

                if (((indexOfNode >= taxonomicLevel) ||    //  nodes has lower rank  OR no rank AND not clicked
                 (d.data.value[taxaIndex] === 'no rank')) && !(data_name in clicked_nodes)){
                    //console.log(['colTaxa', data_name])
                    //console.log(d.data.value[taxaIndex])
                    d = collapse(d, branchLength, true);
                 }
                 /***else if (d.children && !(data_name === root.data.name)){  // inner node and one of the children has lower rank
                    var childWithLowerTaxon = checkTaxonomicLevel(d.children, taxonomicLevel);
                    if(childWithLowerTaxon && !(data_name in clicked_nodes)){
                        d = collapse(d, branchLength, true);
                    }
                 }***/
            }
        }
    });

    // publication ready tree visualisation have no onclickInteractions
    // --> resize the node to their actual child size
    if(!(onclickInteraction)){
        console.log('run resize of tree')
        root.descendants().reverse().forEach(function(d){
            if(d.children){
                d.data._size = d.data.size.slice();
                d.data.size[0] = d3v6.sum(d.children, function(child){return child.data.size[0];});
                d.data.size[0] = d.data.size[0] + 10;
            }else if(d._children){
                d.data._size = d.data.size.slice();
                d.data.size[0] = 15; //d3v6.sum(d.children, function(child){return child.data.size[0];});
                //d.data.size[0] = d.data.size[0] + 10;
            }
        });
        rootHeight = root.data.size[0];
        treeHeight = rootHeight + margin.top + margin.bottom;
        document.getElementById('treeVis').style.height = String(treeHeight).concat('px');

    }else if(returnFromStatic && !(extraInfo)){
        root.data.size[0] = d3v6.sum(root.data.children, function(child){return child.size[0];});
        rootHeight = root.data.size[0]
        treeHeight = rootHeight + margin.top + margin.bottom ;
    }else{
        rootHeight = root.data.size[0]
        treeHeight = rootHeight + margin.top + margin.bottom ;
    }

    //              max tree depth                 max label             spaces
    var width = (root.height * branchLength) +  (max_name_length*6) + 6 + margin.left ;



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
    document.getElementById("treeVis").scrollTop = treeHeight/2-(window.innerHeight*0.8)/2;

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
            .on('click',click)
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
                                .attr("x",  12)
                                .attr("text-anchor", "start")
                                .text(d.data.name);
                    d._children = d.children;
                    d.children = null;
                    //d.data._size = d.data.size;
                    //d.data.size = [15, branchLength];

                    // if nodes is going to be collapsed --> enable publication ready
                    if(window.location.href.includes('phylo')){
                        document.getElementById('public_ready_phylo').disabled = false;
                    }else{
                        document.getElementById('public_ready_taxa').disabled = false;
                    }

                }else if(!(d.children) && !(d._children) && extraData != null && Object.keys(extraData).length === 1){
                    window.open('https://www.ncbi.nlm.nih.gov/protein/'.concat(d.data.name), '_blank');
                } else {
                    clicked_nodes[d.data.name] = 1;  // node is not collapsed
                    d3v6.select('#nodetext' + d.id).attr('fill', d =>  d.children || d._children ? 'transparent' : 'black')
                                             .text(d =>  d.children || d._children ? 'p' : d.data.name);
                    d.children = d._children;
                    d._children = null;
                    /***if(d.data._size){
                        var oldSize = d.data.size;
                        d.data.size = d.data._size;
                        d.data._size = oldSize;
                    }***/
                }
                update(d);

                console.log(hitSelection)
                if (extraData === null){   // taxonomic mapping
                    if(taxonomyLevel.includes(d.data.value[2]) || d.data.value[2] === 'no rank'){
                        stackBars(hitSelection);
                        console.log('run stack')
                    }else{
                        hitBars(hitSelection);
                    }
                } else if (extraData.length === 2) {    // phylogeny calculation (taxa-based)
                    showClades(extraData[0], extraData[1], ownAddInfo, null);
                }

            }

        }

    }

}




// ******************************** Barcharts   ***************************

// show hit values based on the actual state of the tree (phyloblast.html)
function hitBars(value){
    hitSelection = value;
    var barheight = 7;

    var hitValue = null;
    if (value === 'node hits'){  // show node bars
	    hitValue = 0;
    } else if (value === 'subtree hits'){
        hitValue = 1;
    } else {
        d3v6.select('#hitbars')
             .remove();
        return
    }

    d3v6.select('#hitbars')
        .remove();

    var nodes = treeData.descendants();
    var max_hit = d3v6.max(treeData.leaves(), function(d){return d.data.value[0][hitValue];});
    var ticksStep = 0;
    if((max_hit >=  1000) || (max_hit <= 10)){
        ticksStep = 5;
    }else {
        ticksStep = 10;
    }
    var scaleX = d3v6.scaleLinear().domain([0, max_hit]).range([0, svgWidth]);

    var hitbars = d3v6.select('#additionalInfo')
           .append('svg')
           .attr('id', 'hitbars')
           .attr("width", svgWidth+margin.left+margin.right)
           .attr("height", treeHeight)
           .style("font", "10px sans-serif")
           .style("overflow","visible")
           .style("user-select", "none");
      
    hitbars.append('g')
             .attr('transform', 'translate(' + margin.left + ',' + (treeHeight/2-barheight+leftmostnode.y) + ')')
             .call(d3v6.axisTop(scaleX)
                    .ticks(ticksStep, 'f'));

    hitbars.append('g')
            .attr('transform', `translate(` + margin.left + `,${treeHeight/2})`)
            .selectAll('.bars')
            .data(nodes)
            .enter()
            .append('rect')
            .attr('class', 'bars')
            .attr('transform', function(d) { return 'translate(0,' + (d.y-4) + ')';})
   	        .attr('fill',  d => d.children  ? "transparent" : "#377ba8")
	        .attr("width", function(d){ return scaleX(d.data.value[0][hitValue]);})
	        .attr("height", barheight)
            .on("mouseover", function(d, i){showTooltip(i, false, null);})
            .on("mouseout", function(d, i){hideTooltip(i, false, 4);})
            .on("mousemove", function(d, i){moveTooltip(i, false);});

    if((rightmostnode.y-leftmostnode.y) >= (window.innerHeight*0.8)){
        hitbars.append('g')
                 .attr('transform', 'translate(' + margin.left + ',' + (treeHeight/2+rightmostnode.y+barheight) + ')')
                 .call(d3v6.axisBottom(scaleX)
                        .ticks(ticksStep, 'f'));
    }
}

function stackBars(value){
    hitSelection = value;
    var barheight = 7;

    var hitValue = null;
    if (value === 'node hits'){  // show node bars
	    hitValue = 0;
    } else if (value === 'subtree hits'){
        hitValue = 1;
    } else {
        d3v6.select('#hitbars')
            .remove();
        return
    }
    d3v6.select('#hitbars')
        .remove();

    var nodes = treeData.descendants();
    var stacks = [];
    var max_hit = d3v6.max(treeData.leaves(), function(d){return d.data.value[0][hitValue] + d.data.value[1][hitValue];});
    nodes.forEach(function(d){
        stacks.push({'node': d, 'query1': d.data.value[0][hitValue], 'query2': d.data.value[1][hitValue]});
    });
    //console.log(stacks)

    // define ticks of the axis
    var ticksStep = 0;
    if((max_hit >=  1000) || (max_hit <= 10)){
        ticksStep = 5;
    }else {
        ticksStep = 10;
    }

    var scaleX = d3v6.scaleLinear().domain([0, max_hit]).range([0, svgWidth-50]);

    var hitbars = d3v6.select('#additionalInfo')
           .append('svg')
           .attr('id', 'hitbars')
           .attr("width", svgWidth+margin.left+margin.right)
           .attr("height", treeHeight)
           .style("font", "10px sans-serif")
           .style("overflow","visible")
           .style("user-select", "none");

    // x-Axis
    hitbars.append('g')
             .attr('transform', 'translate(' + margin.left + ',' + (treeHeight/2-barheight+leftmostnode.y) + ')')
             .call(d3v6.axisTop(scaleX)
                    .ticks(ticksStep, 'f'));

    // colors of the segment, domain parameter need to be adapted for more than 2 queries
    var colors = d3v6.scaleOrdinal().domain(['query1', 'query2']).range(d3v6.schemeTableau10);

    // stack each value of the subgroup above each other
    var stack = d3v6.stack().keys(['query1', 'query2']);
                  //  .order(d3v6.stackOrderNone)
                   // .offset(d3v6.stackOffsetExpand);
    var stackedData = stack(stacks);
    //console.log(stackedData)

    // generate the bars
    hitbars.append("g")
        .attr('transform', `translate(` + margin.left + `,${treeHeight/2})`)
        .selectAll(".stackData")
        .data(stackedData)
        .enter()
        .append("g")
        // loop over profiles
        .attr("class", 'stackData')
        .attr("fill", function(d,i){return colors(d.key);})
        .selectAll("rect")
        // loop over taxa
        .data(function(d){return d;})
        .enter()
        .append("rect")
        .attr('class', d => d.key)
        .attr('opacity', function(d,i){if(stacks[i].node.children){return 0;}
                                      else{return 1;}})
        .attr('x', function(d){return scaleX(d[0]);})
        .attr('y', function(d,i) {return stacks[i].node.y-4;})
        .attr("width", function(d,i){return scaleX(d[1]) - scaleX(d[0]);})
        .attr("height", barheight);

    var legend = hitbars.selectAll(".legend")
        .data(['query1', 'query2'])
        .enter()
        .append("g")
        .attr("class", "legend")
        .attr('transform', 'translate(' + margin.left + ',' + (treeHeight/2-barheight+leftmostnode.y) + ')');


   legend.append("rect")
        .attr("x", svgWidth-40)
        .attr("y", function(d,i){return i*19;})
        .attr("width", barheight*2)
        .attr("height", barheight*2)
        .style("fill", function(d){return colors(d); });

   legend.append("text")
        .attr("x", svgWidth-20)
        .attr("y", function(d,i){return i*19+barheight;})
        .attr("dy", ".35em")
        .style("text-anchor", "start")
        .text(function(d, i){ return d; });
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
var actualVisualisation = 'clade'
function showClades(taxData, accData, ownAddData, libTree){
    // trace the click state of the additional information
    var uniqueCheck = document.getElementById("uniqueCheck").checked;
    var phylumCheck = document.getElementById("phylumCheck").checked;
    var ownAdditional = false;
    var numberOfAdditionalFeatures = 0;

    if (document.getElementById('ownCheckbox').style.display === 'block'){
        ownAdditional = document.getElementById('ownAdditional').checked;
        if(ownAddData){
            ownAddInfo = ownAddData;
        }else{
            ownAddData = ownAddInfo;
        }
        if(ownAdditional){
            var firstKey = Object.keys(ownAddData)[0]
            numberOfAdditionalFeatures = Object.keys(ownAddData[firstKey]).length;
        }
    }

    // add event Listener to button to trace the actual visualisation
    document.getElementById('cladogram').addEventListener('click', function(){
    actualVisualisation = 'clade';});
    document.getElementById('phylogram').addEventListener('click', function(){
    actualVisualisation = 'phylo';});

    if(libTree){
        phyloTreetree = libTree;
        actualVisualisation = 'phylo';
    }

    if(uniqueCheck || phylumCheck || ownAdditional){
        // color-encoding of the taxids   Problem to much ids/accs as colors
        var taxids = Object.keys(taxData);
        var colorsParent = d3v6.scaleOrdinal().domain(taxids).range(d3v6.schemeCategory10);
        //console.log(taxids)

        //[...new Set(d3v6.schemeDark2.concat(d3v6.schemePaired).concat(d3v6.schemeCategory10).concat(d3v6.schemeTableau10))])


        //remove old taxids
        d3v6.select('#clade_vis')
            .remove();

        var rectSize = 8;

        var nodes = treeData.descendants();
        var dataInfo = null;
        var actualTreeHeight = treeHeight/2;
        var phylotreePresent = actualVisualisation === 'phylo' ? true : false;
        if(phylotreePresent){
            var libTreenodes = phyloTreetree.get_nodes();
            dataInfo = {};                       // position of the nodes in distance tree
            for (var i = 0; i < libTreenodes.length; i++){
                dataInfo[libTreenodes[i].name] = libTreenodes[i];
            }
            var treeContainer = document.getElementsByClassName('phylotree-container')[0];
            var offset = treeContainer.transform.baseVal[0].matrix;
            actualTreeHeight = offset.f;
            console.log(dataInfo)
        }
        var clades = d3v6.select('#additionalInfo')
               .append('svg')
               .attr('id', 'clade_vis')
               .attr("width", ((numberOfAdditionalFeatures+2)*(rectSize*2+10))+105+margin.right)
               .attr("height", treeHeight)
               .style("font", "10px sans-serif")
               .style("overflow","visible")
               .style('background', '#F0F0F0')
               .style("user-select", "none");

        var counter = 0;
        var firstElementSpace = 0;   // first element should be as close as possible to the tree

        // additional information was uploaded iterate over the columns of the input csv
        if(ownAdditional){

            var colorCounter = 0;
            for(var key in ownAddData[firstKey]){
                if(counter === 0){
                    firstElementSpace = 0;
                }else{
                    firstElementSpace = counter*(rectSize*2+10);
                }
                clades.append('g')
                        .attr('transform', 'translate(' + 5 + ',' + actualTreeHeight + ')')
                        .selectAll('.ownInfo')
                        .data(nodes)
                        .enter()
                        .append('rect')
                        .attr('class', 'ownInfo' + key)
                        .attr('transform',function(d,i) { var taxaName = validNodename(d.data.name);
                                        if(phylotreePresent && dataInfo[taxaName] !== undefined && !(dataInfo[taxaName].hidden)){
                                        return 'translate(' + (firstElementSpace)+ ','  + (dataInfo[taxaName].screen_y-(rectSize/2)) + ')';}
                                        else{return 'translate(' + (firstElementSpace)+ ',' + (d.y-(rectSize/2)) + ')'; }})
                        .attr('stroke', function(d){var taxaName = validNodename(d.data.name);
                                if((taxaName in ownAddData) && !(d.children) && (ownAddData[taxaName][key] !== '0')){
                                    if(phylotreePresent && (!(taxaName in dataInfo) || dataInfo[taxaName].hidden )){
                                         return 'none';
                                    }else{return 'black';}
                                }else{
                                    return 'none';}})
                        .attr('fill', function(d){var taxaName = validNodename(d.data.name);
                                if((taxaName in ownAddData) && !(d.children) && (ownAddData[taxaName][key] !== '0')) {

                                    if(phylotreePresent && (!(taxaName in dataInfo) || dataInfo[taxaName].hidden )){
                                        return 'transparent';
                                    }else{return d3v6.schemeDark2[colorCounter];}
                                }else{return 'transparent';}})
                        .attr("width", rectSize*2)
                        .attr("height", rectSize);
                counter = counter + 1;
                colorCounter = colorCounter + 1;
                if(colorCounter === 8){
                    colorCounter = 0;
                }
            }

        }

        if(uniqueCheck){
            if(ownAdditional){
                firstElementSpace = firstElementSpace + rectSize*2+10;
            }/***else{
                firstElementSpace = 0;
            }***/

            clades.append('g')
                   .attr('transform', 'translate(' + 5 + ',' + actualTreeHeight + ')')
                   .selectAll('.accession')
                   .data(nodes)
                   .enter()
                   .append('rect')
                   .attr('class', 'accession')
                   .attr('transform',function(d,i) { var taxaName = validNodename(d.data.name);
                                        if(phylotreePresent && dataInfo[taxaName] !== undefined && !(dataInfo[taxaName].hidden)){
                                        return 'translate(' + (firstElementSpace)+ ','  + (dataInfo[taxaName].screen_y-(rectSize/2)) + ')';}
                                        else{return 'translate(' + (firstElementSpace)+ ',' + (d.y-(rectSize/2)) + ')'; }})
                   .attr('stroke', function(d){var taxaName = validNodename(d.data.name);
                        if(accData.includes(d.data.value[1]) && !(d.children)){
                            if(phylotreePresent && (!(taxaName in dataInfo) || dataInfo[taxaName].hidden )){
                                return 'none';
                            }else{return 'black';}
                        }else{
                            return 'none';}})
                   .attr('fill', function(d){var taxaName = validNodename(d.data.name);
                        if(accData.includes(d.data.value[1]) && !(d.children)) {
                            if(phylotreePresent && (!(taxaName in dataInfo) || dataInfo[taxaName].hidden )){
                                return 'transparent';
                            }else{return d3v6.schemeSet1[3];}
                        }else{return 'transparent';}})
                   .attr("width", rectSize*2)
                   .attr("height", rectSize);
        }

        if(phylumCheck){
            if(ownAdditional && !(uniqueCheck)){
                firstElementSpace = firstElementSpace + rectSize*2+10;
            }

            clades.append('g')
                   .attr('transform', 'translate(' + 5 + ',' + actualTreeHeight + ')')
                   .selectAll('.clades')
                   .data(nodes)
                   .enter()
                   .append('rect')
                   .attr('class', 'clades')
                   .attr('transform', function(d,i) {var taxaName = validNodename(d.data.name);
                                        if(phylotreePresent && dataInfo[taxaName] !== undefined){
                                        return 'translate(' + (firstElementSpace+rectSize*2+10)+ ',' + (dataInfo[taxaName].screen_y-(rectSize/2)) + ')';}
                                        else{return 'translate(' + (firstElementSpace+rectSize*2+10)+ ',' + (d.y-(rectSize/2)) + ')'; }})
                   .attr('stroke', function(d){var taxaName = validNodename(d.data.name);
                        if(d.data.value[0] in taxData) {
                            if(phylotreePresent && (!(taxaName in dataInfo) || dataInfo[taxaName].hidden )){
                                         return 'none';
                            }else{return 'black';}}
                        else{
                            return 'none';}})
                   .attr('fill', function(d){var taxaName = validNodename(d.data.name);
                        if(d.data.value[0] in taxData) {
                            if(phylotreePresent && (!(taxaName in dataInfo) || dataInfo[taxaName].hidden )){
                                         return 'none';
                            }else{ return colorsParent(d.data.value[0]);}}
                        else{
                            return 'none';}})
                   .attr("width", rectSize*2)
                   .attr("height", rectSize);
        }

        var spaceTillLegend = 0;
        if(numberOfAdditionalFeatures === 0){
            spaceTillLegend = (rectSize*8)+5;
        }else{
            spaceTillLegend =(numberOfAdditionalFeatures+2)*(rectSize*2+10)+15;
        }

        if(phylumCheck){
            var legend =  clades.selectAll('.legend')
                                    .data(taxids)
                                    .enter()
                                    .append('g')
                                    .attr('class', 'legend')
                                    .attr('transform', 'translate(' + (spaceTillLegend) + ','+(treeHeight/2+leftmostnode.y) + ')');

            legend.append('rect')
                    .attr('transform', function(d,i){return 'translate(0,' + i*(rectSize+15) + ')';})
                    .attr('fill', function(d){return colorsParent(d);})
                    .attr('stroke', 'black')
                    .attr('width', rectSize*1.5)
                    .attr('height', rectSize*1.5)

            legend.append('text')
                    .attr('transform', function(d,i){return 'translate(0,' + i*(rectSize+15) + ')';})
                    .attr('x',rectSize*1.5+3)
                    .attr('y', rectSize*1.5/2)
                    .attr('dy', '0.31em')
                    .attr("text-anchor", 'start')
                    .attr('fill', 'black')
                    .text(d => taxData[d]);
        }

        var addLabels = [];
        if(ownAdditional){
            addLabels = Object.keys(ownAddData[firstKey]);
        }
        var labels = [];
        if(uniqueCheck && phylumCheck){
            labels = addLabels.concat(['unique', 'phylum']);
        }else if(uniqueCheck){
            labels = addLabels.concat(['unique']);
        }else if(phylumCheck){
            labels = addLabels.concat(['', 'phylum']);
        }else{
            labels = addLabels;
        }

        var labelSpace = rectSize*2;
        if (phylotreePresent){
            labelSpace = rectSize*3;
        }
        clades.append('g')
                .attr('transform', 'translate(' + 5 + ',' + (treeHeight/2-(labelSpace)+leftmostnode.y) + ')')
                .selectAll('.labels')
                .data(labels)
                .enter()
                .append('text')
                .attr('class', 'labels')
                .attr('dy', '0.2em')
                .attr("text-anchor", "start")
                .attr("dominant-baseline", "central")
                .attr('transform', function(d,i){return 'rotate(-90,' + (i*(rectSize*2+10)+5) + ', 0)';})
                .attr('x', function(d,i){ return (i*(rectSize*2+10)+5);})
                .attr('y', 0)
                .text(function(d){return d.substr(0,6);});

    }else{
        d3v6.select('#clade_vis')
            .remove();
    }


   
}

// collapse full tree to given taxonomic level
var previousTaxonomicLevel = 'class';
function collapseTree(rank){

    console.log('run collapse')
    var taxonLevel = taxonomyLevel.indexOf(rank);
    d3v6.select('#tree_vis').remove();
    chart(treeVis, null, taxonLevel, previousTaxonomicLevel, true, false);
    if(taxonomyLevel.includes(treeVis['value'][1]) || (treeVis['value'][1] === 'no rank')){
        hitBars(hitSelection);
    }else{
        stackBars(hitSelection);
    }

    previousTaxonomicLevel = taxonLevel;

}


// remove all white spaces between nodes of the actual tree
function publicationReady(){
    if((previousTaxonomicLevel !== taxonomyLevel.length -1) || (Object.keys(clicked_nodes).length > 0)){
        document.getElementById('returnButton').style.display = 'block'; // enable reset to dynamic visualisation
        d3v6.select('#tree_vis').remove();
        console.log(treeVis)
        chart(treeVis, extraInfo, previousTaxonomicLevel, previousTaxonomicLevel, false, false);
        if (extraInfo === null){
            if(taxonomyLevel.includes(treeVis.value[1])){
                hitBars(hitSelection);
            }else{
                stackBars(hitSelection);
            }
        } else if (document.getElementById('clade_vis')){
            showClades(extraInfo[0], extraInfo[1], ownAddInfo, null);
        }

        // prevent double-click on publication ready
        if(window.location.href.includes('phylo')){
            document.getElementById('public_ready_phylo').disabled = true;
        }else{
            document.getElementById('public_ready_taxa').disabled = true;
            document.getElementById('collapse_menu').disabled = true; // disable rank selection
        }
    }


}


export {chart, hitBars, showClades,  startTreevis, collapseTree, publicationReady, stackBars};
