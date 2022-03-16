"use strict";

/**
 * List of groups to apply a styles to nodes
 */
var node_groups = {
  legend : 'node_legend',
  reactant : 'reactant', 
  product : 'product', 
  main : 'mainspecies'
};

var threshold_types = {
    simple : 'simple',
    cumulated : 'cumulated'
};

/**
 * Manages species network
 */
var mol_network = {
  network_events : false,
  container_area : "species_network",
  network : null,
  species : null,
  nodes : null,
  edges : null,
  height : 600,
  width : 700,
  current_species : null,
  species_history : [],   // navigation history in graph
  current_threshold : null,
  current_threshold_type : threshold_types.simple,
  dataset_id : 0,
  theme_color : "rgba(132, 197, 187, 0.61)",

    
  removeGraph : function(){
    if (this.network !== null)
      this.network.destroy();  
    this.network = null;
    this.network_events = false;
    document.getElementById(mol_network.container_area).innerHTML = null;
  },  
    
  updateGraph : function(species_name, rate_min){
    this.current_species = species_name;
    
    if ( this.species_history.length > 0 ){
      if( species_name !== this.species_history[this.species_history.length-1] )
         this.species_history.push(species_name);
    }else
      this.species_history.push(species_name);
    
 
    this.current_threshold = rate_min;
    if(this.species != null && this.species[species_name] !== undefined ){
      this.buildDatasets(species_name, rate_min, this.current_threshold_type);
      this.enableBorders("1");
      
      var do_plot = true;
      if( !(threshold_types.simple === this.current_threshold_type && rate_min == 0) &&
        !(threshold_types.cumulated === this.current_threshold_type && rate_min == 1) ){
        this.plotNetwork();
      }
      
    }
  },    
  
  /**
   * Manages navigation list
   * returns previous species, false if none
   */
  goBackSpecies : function(){
    var result;

    if ( this.species_history.length <= 1  )
      return false;
        
    this.species_history.pop();           //current species
    result = this.species_history.pop();  //previous species
    return result;
  },
  
  /**
   * take a png snapshot of current graph
   */
  takeSnapshot : function (){
    var img    = this.network.canvas.frame.canvas.toDataURL("image/png");
    window.open(img,'_blank');
  },
  
  enableBorders : function(width){
    document.getElementById(this.container_area).style.borderWidth = width+'px';
  },
  
  plotNetwork : function(){
    var container = document.getElementById(this.container_area);
   
    var data= {
      nodes: this.nodes,
      edges: this.edges
    };

    var theme_font_size = 15;
    
    var options = {
      locale: 'en',
      width: '600px',
      height: '500px',
      interaction : {
        zoomView : false,
        hover : true,
        selectConnectedEdges : false,
        hoverConnectedEdges : false
      },
      manipulation: {
        enabled: false,
        initiallyActive: false
      },
      edges : {
        width : 1.2,
        physics : true,
        shadow : {
          enabled : false,
          size : 10
        },
        arrows : {
          to : {
            scaleFactor : 0.5,
            enabled : true
          }
        },
        font : {
          size : theme_font_size,
          align : 'horizontal',
          strokeWidth: 10 // px
          
        },
        color : {
          hover : this.theme_color   
        }
          
      },
      nodes:{
        shape: 'dot', // has to be dot to have a scalable size        
        physics : false,
        size : 50,
        fixed : {
          x : false,
          y : false
        }
      },    
      groups: {
        reactant: {
          shape: 'box',
          font: {
            color : 'black',
            size: theme_font_size
          },          
          color: {
            border: 'white',
            background: 'white',
            highlight: {
              border: 'yellow',
              background: 'orange'
            },
            hover: {
              border : this.theme_color,
              background : this.theme_color
            }
          }
        },
        node_legend: {
          borderWidth : 1,
          shape: 'box',
          font : {
            color: 'black',
            size: theme_font_size
          },
          color: {
            border: 'white',
            background: 'white',
            highlight: {
              border: 'yellow',
              background: 'orange'
            },
            hover: {
              border : this.theme_color,
              background : this.theme_color
            }
          }
        },   
        product: {
          shape: 'box',
          font : {
            color: 'black',
            size: theme_font_size
          },
          color: {
            border: 'white',
            background: 'white',
            highlight: {
              border: 'yellow',
              background: 'orange'
            },
            hover: {
              border : this.theme_color,
              background : this.theme_color
            }
          }
        },     
        mainspecies: {
          shape: 'box',
          font : {
            color: 'black',
            size: theme_font_size
          },
          color: {
            border: 'black',
            background: 'white',
            highlight: {
              border: 'yellow',
              background: 'orange'
            },
            hover: {
              border : this.theme_color,
              background : this.theme_color
            }
          }
        }
      }
    };
    
    // initialize plot    
    if(this.network == null){
      this.network = new vis.Network(container, data, options);      
    }else{
      this.network.setData(data);
    }
    this.network.focus(1);  
  },
  

  
  buildNode : function(id, reaction, label, group, x, y){
    return [{ id : id, 
              label : label, 
              group : group, 
              x : x, 
              y : y,  
              allowedToMoveX : false,  
              allowedToMoveY : false,
              metadata : {
                reaction : reaction.reaction, 
                rate: reaction.rate, 
                reactants : reaction.reactants,                
                abundances : reaction.reactants_abundances,
                proportion : reaction.proportion                
                //abundanceUnits : reaction.reactants_abundance_units,
                //rateUnit : reaction.unit,                
              }      
            }];
  },
  
  buildEdge : function(id, from, to, label){
    return [{ id : id, 
              from : from, 
              to : to, 
              label : label}];
  },  
  
  removeCurrentSpecies : function(species){
    var result = [];
    for(var i = 0, len = species.length; i < len;i++){
      if ( species[i] !== this.current_species)
        result.push(species[i]);
    }
    return result;
  },
  
  getListByMinimumThreshold : function(reactions, rate_min){
    var result = [];

    for(var i = 0, len = reactions.length; i < len; i++){
      if(get_percentage(reactions[i].proportion) >= rate_min)
        result.push(reactions[i]);
    }
    
    return result;      
  },
  
  getListByCumulatedThreshold: function(reactions, rate_min){
    var result = [];  
    var total = 0;   
    var i = 0;      
    var end = false;
    var rate = 0;
    
    while( end === false && i < reactions.length){
      rate = get_percentage(reactions[i].proportion);
      if(total+rate > rate_min){
        result.push(reactions[i]); 
        end = true;
      }else{      
        total += rate;
        result.push(reactions[i]);      
        i++;
      }
    }
    
    return result;
  },
  
  getListByThreshold : function(reactions, rate_min){    
    if(this.current_threshold_type === threshold_types.simple){
      return this.getListByMinimumThreshold(reactions, rate_min);
    }
    
    if(this.current_threshold_type === threshold_types.cumulated){
      return this.getListByCumulatedThreshold(reactions, rate_min);
    }
    
    //error case, unknown threshold
    alert("invalid threshold type selected");
    return this.getListByMinimumThreshold(reactions, rate_min);   
    
  },
  
  /**
   * sorts species by decreasing mass value
   */
  sortSpeciesList : function(list){  
    var result = [list[0]];
    var is_added = false;  
          
    for(var i = 1; i < list.length; i++){   
      is_added = false;
      for(var j = 0; j < result.length; j++){
        if(list[i].mass > result[j].mass){
          result.splice(j, 0, list[i]);
          is_added = true;
          break;
        }
      }
      if(is_added == false)
        result.push(list[i]);
    }
    return result;
  },  
  
  buildDatasets : function(species_name, rate_min){
    var sp = this.species[species_name];
  
    this.nodes = new vis.DataSet();
    this.edges = new vis.DataSet();
    //central node
    this.nodes.add(this.buildNode(1, [], species_name, node_groups.main, 0, 0));

    var current_id= 2;   
    var reactions = this.getListByThreshold(sp.creation, rate_min); 
    var reactions_count = reactions.length;     
    var spacing = this.height / (reactions_count+1);
    
    if(reactions_count == 1)
      spacing = 0;

    var init_pos = -(reactions_count/2)*spacing + (spacing / 2);    
    
    for(var i = 0, len = reactions_count; i < len; i++){   
      //sorted reactants              
      var reactants = this.sortSpeciesList(reactions[i].reactants);       
      this.nodes.add(this.buildNode(current_id, reactions[i], reactants[0].name, node_groups.reactant, -200, init_pos+spacing*i));   //node
      this.edges.add(this.buildEdge(current_id, current_id, 1, reactants[1].name));  
      current_id++;
      this.nodes.add(this.buildNode(current_id, reactions[i], reactions[i].proportion, node_groups.legend, -280, init_pos+spacing*i));//legend
      current_id++;  
    }

    var reactions = this.getListByThreshold(sp.destruction, rate_min);
    var reactions_count = reactions.length;
    
    var spacing = this.height / (reactions_count+1);
    
    if(reactions_count == 1)
      spacing = 0;

    var init_pos = -(reactions_count/2)*spacing + (spacing / 2);   
    
    for(var i = 0, len = reactions_count ; i < len; i++){ 
      //products 
      var products = this.sortSpeciesList(this.removeCurrentSpecies(reactions[i].products));   
      var reactants = this.sortSpeciesList(this.removeCurrentSpecies(reactions[i].reactants));     
      
      //checks if edge label is different from central node label
      var edge_label = reactants[0].name;
      if(edge_label === species_name){
        edge_label = reactants[1].name;
      }
      
      this.nodes.add(this.buildNode(current_id, reactions[i], products[0].name, node_groups.product, 200, init_pos+spacing*i));    //node
      this.edges.add(this.buildEdge(current_id, 1, current_id, edge_label));    
      current_id++;
      this.nodes.add(this.buildNode(current_id, reactions[i], reactions[i].proportion, node_groups.legend, 280, init_pos+spacing*i)); //legend
      current_id++;
      
    }
  } 
}

//
// Plot
//
var plot2d = {
  plotter : {
    container_area : "plot_main_pos_var",  
    plot : null, 
    log_x_enabled : false,
    log_y_enabled : false,
    dataset : [],
    log_dataset : [],
    log_x_dataset : [],
    log_y_dataset : [],
    
    options : {
      legend : {            
          container : page.mainPosVarGraphLegend,
          noColumns : 2
      },        
      yaxes: [ { 
        show : false
      }, {
        // align if we are to the right
        alignTicksWithAxis: 1,
        position: "right",
        show : false
      }],
      xaxis : {
        show : false
      }
    },
    
    logFunction : function (v) { return Math.log(v); },
    
    toggleLog : function(x_enabled, y_enabled){
      this.log_x_enabled = x_enabled;
      this.log_y_enabled = y_enabled;
    },
       
    addToDataset : function(title, data){
      
      var log_data = [];
      var log_x_data = [];
      var log_y_data = [];
      var log_x, log_y = null;
      for(var i =0; i < data.length; i++ ){
        log_x = Math.log(data[i][0]);
        log_y = Math.log(data[i][1]);
        log_data.push([log_x, log_y]);
        log_x_data.push([log_x, data[i][1]]);
        log_y_data.push([data[i][0], log_y]);
      }
      
      this.pushInDataset(title, this.dataset, data);
      this.pushInDataset(title, this.log_dataset, log_data);
      this.pushInDataset(title, this.log_x_dataset, log_x_data);
      this.pushInDataset(title, this.log_y_dataset, log_y_data);
    },
    
    pushInDataset : function(title, dataset, data){
      var yaxis = '';
      if(dataset.length > 0){
        yaxis = 2;
        
        dataset.splice(1, 1);
        dataset.push({label : title, data : data, yaxis : yaxis});      

      }else{
        dataset.push({label : title, data : data, yaxis : 1});
      }          
    },
    
    clearDataset : function(){
      this.dataset = [];
      this.log_dataset = [];
      this.log_x_dataset = [];
      this.log_y_dataset = [];
    },
    
    doPlot : function (){
      var self = this;
      $(function() {
        if(self.log_x_enabled === true && self.log_y_enabled === true){
          self.plot = $.plot("#"+self.container_area, self.log_dataset, self.options);
        } else if(self.log_x_enabled === true && self.log_y_enabled === false){
          self.plot = $.plot("#"+self.container_area, self.log_x_dataset, self.options);
        } else if(self.log_x_enabled === false && self.log_y_enabled === true){
          self.plot = $.plot("#"+self.container_area, self.log_y_dataset, self.options);
        } else if(self.log_x_enabled === false && self.log_y_enabled === false){
          self.plot = $.plot("#"+self.container_area, self.dataset, self.options);
        }
      });
    },
    
    highlight : function(point){
      if(this.plot !== null){
        for(var i = 0;i < this.dataset.length; i++){
          this.plot.highlight(i, point);
        }
      }
    },
    
    unhighlight : function(){
      if(this.plot !== null)
        this.plot.unhighlight();
    }
  },
  
  plotTemperature : function(main_pos_vars, temperatures, step){
    var exponentialMainPosVars = []
    var plotted_main_pos_vars = [];
    
    // display only 1/main_pos_var_step of values
    exponentialMainPosVars[0] = main_pos_vars[0].toExponential(4);   
    for(var i = step, len = main_pos_vars.length; i < len; i = i+step){ 
      exponentialMainPosVars[i] = main_pos_vars[i].toExponential(4);   
      plotted_main_pos_vars.push([exponentialMainPosVars[i], temperatures[i]]);
    }

    this.plotter.addToDataset("Temperature", plotted_main_pos_vars);
    this.plotter.doPlot(); 
    //this.plotter.highlight(0);    
  },  
  
  plotDensity : function(main_pos_vars, densities, step, species){
    var plotted_n = [];
    // display only 1/main_pos_var_step of values
    for(var i = step, len = main_pos_vars.length; i < len; i = i + step){ 
      plotted_n.push(([main_pos_vars[i], densities[i].toExponential(2)]));
    }        
    
    this.plotter.addToDataset("n("+ species + ")", plotted_n);
    this.plotter.doPlot();      
  },
  
  toggleLog : function(x_enabled, y_enabled){
    this.plotter.toggleLog(page.isLogXPlot.checked, page.isLogYPlot.checked);
    this.plotter.doPlot();    
  },
  
  highlight : function(point){
    this.plotter.unhighlight();
    this.plotter.highlight(page.mainPosVarInput.value/config.main_pos_var_step);  
  }
  
};




//
// Spinner
//
var opts = {
  lines: 5, // The number of lines to draw
  length: 36, // The length of each line
  width: 7, // The line thickness
  radius: 4, // The radius of the inner circle
  corners: 1, // Corner roundness (0..1)
  rotate: 56, // The rotation offset
  direction: 1, // 1: clockwise, -1: counterclockwise
  color: '#000', // #rgb or #rrggbb or array of colors
  speed: 1.6, // Rounds per second
  trail: 93, // Afterglow percentage
  shadow: false, // Whether to render a shadow
  hwaccel: false, // Whether to use hardware acceleration
  className: 'spinner', // The CSS class to assign to the spinner
  zIndex: 2e9, // The z-index (defaults to 2000000000)
  top: '30%', // Top position relative to parent
  left: '40%' // Left position relative to parent
};

