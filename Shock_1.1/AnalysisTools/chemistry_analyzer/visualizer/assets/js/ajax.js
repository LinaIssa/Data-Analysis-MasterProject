
"use strict";

// 
// Configuration variables
//
var config = {
  service_url : "http://localhost:12345/chemistry-analyzer",
  main_pos_var_step : 1
}

//
// Ajax functions calling server
// 
var ajax = {
  //allSpecies : [],
  sessionId : null,   // session id, created by client and submitted to server 
                                      // ok while local application
                                      
   
  getData : function(species, main_pos_var){
    var spinner = null;
    var self = this;
    jQuery.ajax({
        type: "GET",
        url: config.service_url+"/session/"+this.sessionId+"/mainposvar/"+main_pos_var+"/species/"+species,
        dataType: 'json',
        beforeSend : function(){
          if (page.isNewFile){
            mol_network.removeGraph();
            page.showSpinner();
            page.isNewFile = false;
          }
        },
        success: function(data){        
          mol_network.species = data;
          mol_network.updateGraph(species, page.thresholdInput.value/100, mol_network.current_threshold_type);  
          
          page.enableButtons();
          page.toogleBackButton(function(){if(mol_network.species_history.length >= 2) return true; return false});        
          page.updateSpeciesDisplay(species, data[species]["abundance"], data[species]["abundance_unit"]);
          page.updateInfosDisplay(mol_network.nodes, data[species]["abundance"]);        
          
          if(mol_network.network_events === false){
            mol_network.network.on('hoverNode', function (hovered) {  
              if(mol_network.nodes.get(hovered.node).group === node_groups.product || mol_network.nodes.get(hovered.node).group === node_groups.reactant)              
                page.initReactionsTable();                
                if( !(threshold_types.simple === mol_network.current_threshold_type && mol_network.current_threshold == 0) 
                    && !(threshold_types.cumulated === mol_network.current_threshold_type && mol_network.current_threshold == 1)){
                  //all reactions displayed in reactions table but not on graph if threshold == 0 or 100
                  page.highlightReaction('reaction-'+hovered.node);
                }
            });
            mol_network.network.on('select', function (properties) {               
              var species = null;
              if(properties.nodes.length == 1 ){ // click on a node
                species = mol_network.nodes.get(properties.nodes[0]).label;   
              }else if(properties.edges.length == 1){ // click on an edge
                species = mol_network.edges.get(properties.edges[0]).label;                        
              }
              if(species != null && page.allSpecies[species] != undefined){ //no request for grain/photon ... 
                page.goBackSpecies(species);      
                self.getData(species, page.mainPosVarInput.value);                  
                self.getAbundances(species);       
              }
            });        

          }  
          mol_network.network_events = true;
          page.clearSavedFiles();
          
          page.hideSpinner();
        },
        error: function(jqXHR, textStatus, errorThrown ){
          alert("Error while querying server");    
        }
    });
  },

  getGasConditions : function(main_pos_var){
    jQuery.ajax({
        type: "GET",
        url: config.service_url+"/session/"+this.sessionId+"/gasconditions/mainposvar/"+main_pos_var,
        dataType: 'json',
        success: function(data){
          var result = "<ul>";
          var others = [];
          var temperatures = [];
          var densities = [];        
          
          for(var cond in data){ 
            if(cond.indexOf("nd") > -1 && cond !== "nd(e-)")  {
              densities.push(cond);
            }else if (cond.indexOf("Td") > -1){
              temperatures.push(cond);
            }else{          
              others.push(cond);
            }
          }
          temperatures.sort();
          densities.sort();

          page.updateGasConditionsDisplay(datalist_to_html(others, data)); 
          page.updateGrainConditionsDisplay(datalist_to_html(temperatures, data), datalist_to_html(densities, data));    
        },
        error: function(jqXHR, textStatus, errorThrown ){
          alert("Error while querying server");    
        }
    });
  },

  getMainPosVars : function(){
    var self = this;
    jQuery.ajax({
        type: "GET",
        url: config.service_url+"/session/"+this.sessionId+"/allmainposvars",
        dataType: 'json',
        beforeSend : function(){
          if (page.isNewFile){
            mol_network.removeGraph();
            page.showSpinner();
            page.isNewFile = false;
          }
        },
        success: function(data){
          page.updateMainPosVarInput(data['main_pos_vars']);
          plot2d.plotTemperature(data['main_pos_vars'], data['temperatures'], config.main_pos_var_step);
          page.hideSpinner();
        },
        error: function(jqXHR, textStatus, errorThrown ){
          alert("Error while querying server");    
        }
    });
  },
  
  saveNetwork : function(filename){
    var self = this;
    jQuery.ajax({
        type: "GET",
        url: config.service_url+"/session/"+this.sessionId+"/network",
        dataType: 'json',
        success: function(data){
          self.saveData(filename, data, page.displayNetworkFile); 
        },
        error: function(jqXHR, textStatus, errorThrown ){
          alert("Error while querying server");    
        }
    });
  },

  getSpecies : function(){
    var self = this;
    jQuery.ajax({
        type: "GET",
        url: config.service_url+"/session/"+this.sessionId+"/allspecies",
        dataType: 'json',
        success: function(data){
          page.updateSpeciesList(data.sort())
          page.triggerEvent(page.speciesList, new Event('change'));
        },
        error: function(jqXHR, textStatus, errorThrown ){
          alert("Error while querying server");    
        }
    });
  },

  saveData : function(filename, data, success_callback){
    jQuery.ajax({
        type: "POST",
        url: config.service_url+"/savefile",
        dataType: 'json',
        data: { 'filename' : filename, 'string_data' : data},
        success: function(data){                
          success_callback(config.service_url, filename);      
        },
        error: function(jqXHR, textStatus, errorThrown ){
          alert("Error while writing file");    
        }
    });
  },
  
  checkAndSaveData : function(filename, savefunction){
    var self = this;
    var svfunction = savefunction;
    jQuery.ajax({
        type: "GET",
        url: config.service_url+"/existfile/"+filename,
        dataType: 'json',
        success: function(response){
          // checks if file already exists
          var answer = false;
          if(response['status'] === true){
            answer = confirm ("File already exists. Overwrite ?");
          }
          //writes if not
          if(response['status'] === false || answer === true){
            svfunction();
          }
        },
        error: function(jqXHR, textStatus, errorThrown ){
          alert("Error while checking file status");    
        }
    });
  },  
  
  saveRatesMainPosVarData : function(species, threshold, filename){
    jQuery.ajax({
        type: "GET",
        url: config.service_url+"/session/"+this.sessionId+"/saverates/species/"+species+"/threshold/"+threshold+"/filename/"+filename,
        dataType: 'json',
        beforeSend : function(){
          page.showSpinner();
        },
        success: function(data){
          page.hideSpinner();
          page.displayFormationDestructionFiles(config.service_url, filename, data['formation'], data['destruction']);
        },
        error: function(jqXHR, textStatus, errorThrown ){
          alert("Error while writing file");    
        }
    });
  },

  getAbundances : function(species){
    var self = this;
    jQuery.ajax({
        type: "GET",
        url: config.service_url+"/session/"+this.sessionId+"/abundances/species/"+species,
        dataType: 'json',
        success: function(data){        
          plot2d.plotDensity(page.allMainPosVars, data, config.main_pos_var_step, species);
        },
        error: function(jqXHR, textStatus, errorThrown ){
          alert("Error while querying server");    
        }
    });
  },
  
  setMainPositionVar : function(){
    var self = this;
    jQuery.ajax({
        type: "GET",
        url: config.service_url+"/session/"+this.sessionId+"/positionvarname",
        dataType: 'json',
        success: function(data){        
          page.setMainPositionVarName(data['name'], data['unit']);
        },
        error: function(jqXHR, textStatus, errorThrown ){
          alert("Error while querying server");    
        }
    });
  },

  getFiles : function(){
    jQuery.ajax({
        type: "GET",
        url: config.service_url+"/allfiles",
        dataType: 'json',
        success: function(data){
          page.updateFilesList(data);
        },
        error: function(jqXHR, textStatus, errorThrown ){
          alert("Error while querying server");    
        }
    });
  },

  loadFile : function(filename){
    var self = this;
    if( this.sessionId != null ){
      jQuery.ajax({
          type: "GET",
          url: config.service_url+"/session/"+this.sessionId+"/changefile/file/"+filename,
          dataType: 'json',
          success: function(data, textStatus, xhr){
            self.getMainPosVars();
            self.getSpecies();
            page.setModelName(filename);
          },
          complete: function( xhr, settings){
            if(xhr.status == "204"){
              alert("File not found");
            }
          },
          error: function(xhr, textStatus, errorThrown ){
            alert("Error while querying server");    
          }
      });
    }else{
      alert("No session id");
    }
  },
  
  runSession : function(exec){
    var self = this;
    var toexec = exec;
    jQuery.ajax({
        type: "GET",
        url: config.service_url+"/opensession",
        dataType: 'json',
        success: function(data, textStatus, xhr){
          self.sessionId = data["sessionid"];          
          toexec();
        },
        error: function(xhr, textStatus, errorThrown ){
          alert("Unable to open a session. Is the server running ? ");    
        }
    });
  }
}
