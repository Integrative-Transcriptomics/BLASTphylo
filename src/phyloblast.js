/***
# * PhyloBlast
# * @author jennifer mueller
# * 
# * based on: https://observablehq.com/@d3/collapsible-tree
              https://observablehq.com/@bumbeishvili/horizontal-collapsible-d3-flextree
# ***/

import * as d3 from 'd3';



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
   if  (tree === null){
	return 0
   }else{
        return removeNullroot(tree);
   }
}







// ****************** Tooltip: show Nodevalues ***************************
function showTooltip(node, barhover){

  var tooltip = d3.select("#tree")
	.append("div")
	.attr('id', 'tooltip')
	.style("position","absolute")
	.style("visibility","hidden")
	.style("border-radius", "10px")
	.style("padding", "5px")
	.style("color", "white")
	.style("background-color", "rgba(35, 29, 51, 0.7)");

  var tooltip_x_pos = 10;
  var tooltip_y_pos = 10;
  tooltip.style("visibility", "visible")
	.html('<b>' + node.data.name + '</b> <br />' + 'hits: ' + String(node.data.value[0]) + '<br />' +
	 'subtree hits: ' + String(node.data.value[1]))
	.style("left", (window.event.pageX + tooltip_x_pos) + "px")
	.style("top", (window.event.pageY + tooltip_y_pos) + "px")

  // change styling of the circle to indicate on click action
  if(barhover){
	  if(node.children || node._children){
	 	 d3.select('#nodecircle' + node.id).attr('fill', 'crimson')
					    .attr('r', 5)
 	  }
  }
}



function hideTooltip(node, barhover){
  d3.select('#tooltip').remove();

  // change styling of the circle back to normal
  if(barhover){
	  if(node.children || node._children){
	  	d3.select('#nodecircle' + node.id).attr('r', 3.5)
	 	 .attr("fill", node._children || node.children ? "#555" : "#999")
 	 }
  }

}
const flextree = require('d3-flextree').flextree;

// ****************** Treevisualisation: basic interactive tree ***************************
// main visualisation of the tree
function chart(data, extraData, branchLength, svgWidth, margin) {
    var treeData = null;
    const duration = 750;
    let i = 0;

    // Reverse size parameters, in order to maintain order in horizontal layout
    loopOverHierarchy(data,d=>{
      if(Array.isArray(d.size)){
        if(!d._size) d._size = d.size.slice();
        d.size = d._size.slice().reverse();
        d._size = [10, branchLength]
      }
    })


    const flexLayout = flextree();
    let root = flexLayout.hierarchy(data, function(d) {
            return d.children;
    });

    var number_of_leaves = root.leaves().length;
    // define the maximal space for the barchart
    var max_name_length = 0
    root.descendants().forEach(function(d){
       var data_name = d.data.name;
       if((data_name.length*4.96) > max_name_length){ //4.96 == 0.31 em in pixel, 10px max height
          max_name_length = data_name.length;
       }

        // collapse tree if more than 50 leaves should be displayed
        if (number_of_leaves >= 50 && d.depth >= (root.height/2)){
            collapse(d);
        }
    });



    var width = (root.height * branchLength) + (max_name_length*4.96);
    var treeHeight = root.data.size[0] + margin.top + margin.bottom ;
    //this.setState({svgHeight: treeHeight});

    // append the svg object to the body of the page
    var svg = d3.select('#tree')
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
  



    update(root);



    // Collapse the node and all it's children
    function collapse(d) {
        if (d.children) {
            d._children = d.children
            var size = d.data._size;
            d.data._size = d.data.size;
            d.data.size = size;
            d._children.forEach(collapse)
            d.children = null
        }
    }

    function loopOverHierarchy(d,callback){
      callback(d);
      if(d.children) d.children.forEach(c=>loopOverHierarchy(c,callback))
      if(d._children) d._children.forEach(c=>loopOverHierarchy(c,callback))
    }

    function update(source) {
         
        // Assigns the x and y position to the nodes
        treeData = flexLayout(root);
	    //this.setState({actualData: treeData});


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
        var node = g.selectAll('g.node')
            .data(nodes, d=> d.id || (d.id = ++i));

        // Enter any new modes at the parent's previous position.
        var nodeEnter = node.enter().append('g')
            .attr('class', 'node')
            .attr("transform", function(d) {
                return "translate(" + source.x0 + "," + source.y0 + ")";
            });

        // Add Circle for the nodes
        nodeEnter.append('circle')
            .attr('class', 'node')
	        .attr('id', function(d){return 'nodecircle' + d.id;})
            .attr('r', 3.5)
            .attr("fill", d => d._children || d.children ? "#555" : "#999")
            .attr("stroke-width", 10)
            .attr("cursor", "pointer")
            .on('click', click)
            .on("mouseover", function(d, i){showTooltip(i, true);})
            .on("mouseout", function(d, i){hideTooltip(i, true);});

        // Add labels for the nodes
        nodeEnter.append('text')
            .attr('id', function(d){ return 'nodetext' + d.id;})
            .attr('dy', '0.31em')
            .attr("x", d => d.children ? -6 : 6)
            .attr("text-anchor", d => d.children  ? "end" : "start")
            .attr('fill', d =>  d.children ? 'transparent' : 'black')
            .text(d =>  d.children ? "p" : d.data.name)
            .attr('cursor', 'pointer')
            .on('click', click)
            .on("mouseover", function(d, i){showTooltip(i, true);})
            .on("mouseout", function(d, i){hideTooltip(i, true);});

        // UPDATE
        nodeEnter.merge(node)
            .transition()
            .duration(duration)
            .attr("transform", function(event, i, arr) { 
                const d = d3.select(this).datum();
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


    }
            // Toggle children on click.
    function click(event, d) {
            if (d.children) {
                d3.select('#nodetext' + d.id).attr('fill', 'black')
            				.attr("x",  6)
            				.attr("text-anchor", "start")
                            .text(d.data.name);
                d._children = d.children;
                d.children = null;
                var oldSize = d.data._size;
                d.data._size = d.data.size;
                d.data.size = oldSize;

            } else {
                d3.select('#nodetext' + d.id).attr('fill', d =>  d.children || d._children ? 'transparent' : 'black')
                                         .text(d =>  d.children || d._children ? 'p' : d.data.name);
                d.children = d._children;
                d._children = null;
                var oldSize = d.data.size;
                d.data.size = d.data._size;
                d.data._size = oldSize;
            }
            update(d);

            if (extraData == null){
              hitBars(treeData, treeHeight, svgWidth, margin);
            } else {
              showClades(extraData[0], extraData[1], treeData, treeHeight, svgWidth, margin);
            }
    }


  //return  svg.node();
}




// ******************************** Barcharts   ***************************


// show hit values based on the actual state of the tree (phyloblast.html)
function hitBars(treeData, treeHeight, svgWidth, margin){
  var value = document.getElementById('tree_menu').value;
  var hitValue = null;
  if (value === '0'){  // show node bars
	hitValue = 0;
  } else if (value === '1'){
    hitValue = 1;
  } else {
         d3.select('#hitbars')
        .remove();
        return
  }

  d3.select('#hitbars')
        .remove();

  var nodes = treeData.descendants();
  var max_hit = d3.max(treeData.leaves(), function(d){return d.data.value[hitValue];});

  var scaleX = d3.scaleLinear().domain([0, max_hit]).range([0, svgWidth]);

  var hitbars = d3.select('#additionalInfo')
           .append('svg')
           .attr('id', 'hitbars')
           .attr("width", svgWidth+margin.left+margin.right)
           .attr("height", treeHeight)
           .style("font", "10px sans-serif")
           .style("user-select", "none");
      
      hitbars.append('g')
             .attr('transform', 'translate(' + margin.left + ',' + (margin.top+10) + ')')
             .call(d3.axisTop(scaleX));

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
	   .attr("height", 7)
        .on("mouseover", function(d, i){showTooltip(i, false);})
        .on("mouseout", function(d, i){hideTooltip(i, false);});

     hitbars.append('g')
             .attr('transform', 'translate(' + margin.left + ',' + (treeHeight-margin.bottom+20) + ')')
             .call(d3.axisBottom(scaleX));

}

// show parent taxa as color-coded bars  (phylogeny.html)
function showClades(taxData, accData, treeData, treeHeight, svgWidth, margin){

   var checked = document.getElementById('clade_info').checked
   
   // color-encoding of the taxids   Problem to much ids/accs as colors
   var taxids = Object.keys(taxData);
   var colorsParent = d3.scaleOrdinal().domain(taxids).range(d3.schemeCategory10);
  //console.log(taxids)
  // color-encoding for the hit accessions
  /***var accs = new Set()
  for (var key in accData){
     accs.add(accData[key]);
  }
  //console.log(accs) 
  var colorsAccs = d3.scaleOrdinal().domain([...accs]).range(d3.schemeDark2);
    ***/
//[...new Set(d3.schemeDark2.concat(d3.schemePaired).concat(d3.schemeCategory10).concat(d3.schemeTableau10))])

   if(checked){
      //remove old taxids 
      d3.select('#clade_vis')
        .remove();
      
      var rectSize = 8;

      var nodes = treeData.descendants();
  
     var clades = d3.select('#additionalInfo')
           .append('svg')
           .attr('id', 'clade_vis')
           .attr("width", svgWidth+margin.left+margin.right)
           .attr("height", treeHeight)
           .style("font", "10px sans-serif")
           .style("user-select", "none");


      clades.append('g')
           .attr('transform', `translate(` + margin.left + `,${treeHeight/2})`)
           .selectAll('.clades')
           .data(nodes)
           .enter()
           .append('rect')
           .attr('class', 'clades')
           .attr('transform', function(d) { return 'translate(0,' + (d.y-4) + ')';})
           //.attr("x", 0 )
	   //.attr("y", function(d){return d.y;})
   	   .attr('fill', function(d){if(d.data.value[0] === 0 || d.children) {
				return 'transparent';}
                                else{ 
                                return colorsParent(d.data.value[0]);}})	
	   .attr("width", rectSize)
	   .attr("height", rectSize);
     /***
      clades.append('g')
           .attr('transform', `translate(` + margin.left + `,${treeHeight/2})`)
           .selectAll('.accession')
           .data(nodes)
           .enter()
	   .append('rect')
           .attr('class', 'accession')
           .attr('transform', function(d) { return 'translate(' + (rectSize+10) +',' + (d.y-4) + ')';})
           //.attr("x", 0 )
	   //.attr("y", function(d){return d.y;})
   	   .attr('fill', function(d){if(d.data.value[1] in Object.keys(accData) || d.children) {
				return 'transparent';}
                                else{ 
                                return colorsAccs(accData[d.data.value[1]]);}})
	   .attr("width", rectSize)
	   .attr("height", rectSize);
    

      
      var legend = d3.select('#clade_vis')
         .selectAll('.legend')
         .data(taxids)
         .enter()
         .append('g')
         .attr('class', 'legend')
         .attr("transform", function(d,i){return "translate(10," + i*15 + ")"; });

     legend.append("rect")
	.attr("x", margin.left + 20)
	.attr("y", margin.top - 10)	
	.attr("width", 10)
	.attr("height", 10)
	.style("fill", function(d){ return colors(d);} );

     legend.append("text")
	.attr("x", margin.left + 35 )
	.attr("y", margin.top - 5)
	.attr("dy", ".35em")
	.style("text-anchor", "start")
        .attr('fill', 'black')
	.text(function(d, i){ return dataset[1][d][0]; })
     ***/ 
  }     
  else{
      d3.select('#clade_vis')
        .remove();

  }    
   
}

// ****************** Eventlistener ***************************
/***if (document.getElementById('tree_menu') != null){
    document.getElementById('tree_menu').addEventListener('change', hitBars);
} else if (document.getElementById('clade_info') != null){
    document.getElementById('clade_info').addEventListener('change', showClades);
}
***/

export {chart, hitBars, showClades,  startTreevis};
