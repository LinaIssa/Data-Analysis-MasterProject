"use strict";

/**
 * finds index of main_pos_vars value closest to user input among all main_pos_var value
 */
function get_closest_main_pos_var(user_value){
  var i = 0;
  var found = false;
  var result = 0;  
  
  //value >= main_pos_var_max
  if(user_value >= parseFloat(page.allMainPosVars[page.allMainPosVars.length-1])){
    found = true;
    result = page.allMainPosVars.length-1;
  } 
  
  //value <= main_pos_var_min
  else if(user_value <= parseFloat(page.allMainPosVars[0])){
    found = true;
    result = 0;
  }
  
  // search closest main_pos_var
  else {
    while(i < page.allMainPosVars.length && found === false){
      var compare_value = parseFloat(page.allMainPosVars[i]);      
      
      if( compare_value >= user_value){ 
        if( Math.abs(user_value-compare_value) >  Math.abs( user_value-parseFloat(page.allMainPosVars[i-config.main_pos_var_step]) ) ) {
          result = i-config.main_pos_var_step;
        } else {
          result = i;                
        }
        found = true;
      }else if ( compare_value === user_value){
        result = i;
        found = true;
      }
      i++;
    }
  }  
  return result;
}

/**
 * updates display of gas conditions
 */
function update_gas_conditions(){
  ajax.getData(page.speciesList.value, page.mainPosVarInput.value);
  ajax.getGasConditions(page.mainPosVarInput.value);         
}

/**
 * updates species graph display
 */
function update_graph(threshold){
  page.updateGraphVisibility(page.thresholdInput.value);
  mol_network.updateGraph(mol_network.current_species, threshold/100);
}


//
// Page initialization
//
var page_events = function() {
  
  ajax.getFiles();    

  // update graph from local data
  page.thresholdInput.addEventListener('input', function(){
    page.displayThresholdValue();
    
    update_graph(page.thresholdInput.value);
    page.updateInfosDisplay(mol_network.nodes, mol_network.species[page.speciesList.value]["abundance"]); 
  });

  page.speciesList.addEventListener('change', function(){    
    update_gas_conditions();
    ajax.getAbundances(page.speciesList.value);
  });
  
  page.simpleThreshold.addEventListener('click', function(){
    mol_network.current_threshold_type = threshold_types.simple;
    update_graph(page.thresholdInput.value);
  });
  
  page.cumulatedThreshold.addEventListener('click', function(){
    mol_network.current_threshold_type = threshold_types.cumulated;
    update_graph(page.thresholdInput.value);
  });

  // call server and update graph
  page.mainPosVarInput.addEventListener('change', function(){    
    page.mainPosVarManualInput.value = page.allMainPosVars[page.mainPosVarInput.value].toExponential(4);
    if(page.speciesList.value !== ""){    
      update_gas_conditions();
    }  
  });  
  
  // call server and update graph
  page.mainPosVarInput.addEventListener('input', function(){  
    plot2d.highlight(page.mainPosVarInput.value/config.main_pos_var_step);    
  });  
  
  //user enters main_pos_var value manually
  page.mainPosVarManualInput.onkeypress = function(e){
    
    if (!e) e = window.event;
    var keyCode = e.keyCode || e.which;
    var user_value = parseFloat(page.mainPosVarManualInput.value);
    var searched_main_pos_var = 0;
    if (keyCode === 13){// "return" key
      searched_main_pos_var = get_closest_main_pos_var(user_value);
      // select value in range, if not in range will return closest one
      page.mainPosVarInput.value = searched_main_pos_var; 
      page.mainPosVarManualInput.value = page.allMainPosVars[page.mainPosVarInput.value];
      //update value for range element            
      update_gas_conditions(); 
      return false;
    }
  }  
    page.mainPosVarInput.addEventListener('input', function(){
    page.mainPosVarManualInput.value = page.allMainPosVars[page.mainPosVarInput.value];
  });  
  
  page.filesList.addEventListener('change', function(){    
    ajax.loadFile(page.filesList.value); 
    ajax.setMainPositionVar();
    page.isNewFile = true;   
  });
  
  page.isLogXPlot.addEventListener('change', function(){
    plot2d.toggleLog(page.isLogXPlot.checked, page.isLogYPlot.checked );
  });   
  
  page.isLogYPlot.addEventListener('change', function(){
    plot2d.toggleLog(page.isLogXPlot.checked, page.isLogYPlot.checked);
  });   
   
  //checkboxes to display gas conditions  
  page.setInitialGasCondition();
  page.setInitialGrainConditions();  

  page.hasGasConditions.addEventListener('change', function(){
    page.displayGasConditions();
  });      

  page.hasGrainConditions.addEventListener('change', function(){
    page.displayGrainConditions();
  });
 
  page.displayGrainConditions();
  page.setInitialReactionInformations();  
  
  page.hasReactionInformations.addEventListener('change', function(){
    page.displayReactionInformations();
  });  
  
  //buttons
  page.exportRatesDisplayedButton.addEventListener('click', function(){
    if(mol_network !== undefined){
      var filename = page.getModelName().replace('.hdf5', '')+"_"+mol_network.current_species+"_rates.txt";
      var data = formatter.dataAsText(mol_network.nodes, mol_network.species[mol_network.current_species]["abundance"]);
      ajax.checkAndSaveData(filename, function(){ajax.saveData(filename, data, page.displayReactionRatesFile);});
    }
  });  
  
  page.exportPictureButton.addEventListener('click', function(){    
    if(mol_network !== undefined){
      mol_network.takeSnapshot();
    }
  });
  
  page.exportRatesMainPosVarTxtButton.addEventListener('click', function(){
    if(mol_network !== undefined){
      page.clearSavedFiles();
      var filename = mol_network.current_species+"_"+page.filesList.value+".txt";
      ajax.checkAndSaveData("formation_"+filename, function(){ajax.saveRatesMainPosVarData(mol_network.current_species, mol_network.current_threshold, filename);});      
    }
  });
  
  page.exportNetworkButton.addEventListener('click', function(){
    if(mol_network !== undefined){
      var filename = page.getModelName().replace('.hdf5', '')+"_network.txt";
      ajax.saveNetwork(filename);
    }
  });  

  page.backButton.addEventListener('click', function(){
    page.clearSavedFiles();
    var backSpecies = mol_network.goBackSpecies();
    if (backSpecies !== false){
      ajax.getData(backSpecies, page.mainPosVarInput.value);
      ajax.getAbundances(backSpecies);
      page.goBackSpecies(backSpecies);
    }    
  });
  
  page.hasControlPanel.addEventListener('click', function(){
    page.switchLeftPanel(page.hasControlPanel.id);
  });  

  page.hasPhysicalConditionPanel.addEventListener('click', function(){
    page.switchLeftPanel(page.hasPhysicalConditionPanel.id);
  });
  
  page.switchLeftPanel(page.hasControlPanel.id);
}

$( document ).ready(function() { 
  //get an id for the session and load the page
  ajax.runSession(page_events);
});
