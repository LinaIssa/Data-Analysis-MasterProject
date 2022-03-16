"use strict";

/**
 * pad a string to reach desired size
 */
function pad(pad, str, padLeft) {
  if (typeof str === 'undefined') 
    return pad;
  if (padLeft) {
    return (pad + str).slice(-pad.length);
  } else {
    return (str + pad).substring(0, pad.length);
  }
}

var formatter = {

  dataAsHtml : function(nodes, abundance){
    var reactions = {formation : '', destruction : ''};
    var node = null;
    var ids = nodes.getIds();
    var str = "";
    for (var i = 0, len = ids.length; i < len; i++ ){
      str = "";
      node = nodes.get(ids[i]);          

      if(node.metadata.reaction !== undefined && node.group !== node_groups.legend){  
        str += '<tr id="reaction-'+ids[i]+'">';
        str += '<td class="left_align">'+this.formatReactionHtml(node.metadata.reaction)+"</td>";
        str += "<td>"+str_to_exp_number(node.metadata.abundances[0])+"</td>";
        str += "<td>"+str_to_exp_number(node.metadata.abundances[1])+"</td>";
        str += "<td>"+str_to_exp_number(get_constant(node.metadata.rate, node.metadata.abundances))+"</td>";  
        str += "<td>"+str_to_exp_number(node.metadata.rate)+"</td>";      
        str += "<td>"+node.metadata.proportion+"</td>";         
        str += "</tr>";
      }
      
      if(node.group === node_groups.product)
        reactions.destruction = reactions.destruction + str;
        
      if(node.group === node_groups.reactant)
        reactions.formation = reactions.formation + str;
    }

    var colnbr = 7;
    var result = '<table id="all_reactions_table">';
    result = result + '<tr style="border-bottom : 1px solid #BDBDBD"><th colspan="'+colnbr+'"><div>formation reactions</div></th></tr>';   
    result = result + '<tr style="border-bottom : 1px solid #BDBDBD"><th class="left_align">A + B &#8594; C + D</th><th>n(A)<br/>(cm-3)</th><th>n(B)<br/>(cm-3)</th><th><span title="reaction constant">k</span></th><th>rate<br/>(cm-3 s-1)</th><th>proportion</th></tr>';
    result = result + reactions.formation;
    result = result + '<tr><th colspan="'+colnbr+'" >&nbsp;</th></tr>'; 
    result = result + '<tr><th colspan="'+colnbr+'" >&nbsp;</th></tr>'; 
    result = result + '<tr><th colspan="'+colnbr+'" >&nbsp;</th></tr>'; 
    result = result + '<tr style="border-bottom : 1px solid #BDBDBD"><th colspan="'+colnbr+'"><div>destruction reactions</div></th></tr>';    
    result = result + '<tr style="border-bottom : 1px solid #BDBDBD"><th class="left_align">A + B &#8594; C + D</th><th>n(A)<br/>(cm-3)</th><th>n(B)<br/>(cm-3)</th><th><span title="reaction constant">k</span></th><th>rate<br/>(cm-3 s-1)</th><th>proportion</th></tr>';   
    result = result + reactions.destruction;
    result = result + '</table>';  
    
    return result;
    
  },
  
  /**
   * Replaces all occurences of target by replace in string
   */ 
  replaceAll : function(string, target, replace){  
    return string.split(target).join(replace);    
  },
  
  /**
   * Add <sup> markup around + characters in species names
   */
  formatReactionHtml : function(reaction){        
    var result = " "+reaction;
  
    result = result.replace("e-", "e<sup>-</sup>");    
    
    var reg=new RegExp("_(..)", "g");
    result = result.replace(reg, "<sup>$1</sup>");    

    // numbers after element name
    var reg=new RegExp("([^0-9\s])([0-9]+)", "g");
    result = result.replace(reg, "$1<sub>$2</sub>");
    
    // numbers before element name
    var reg=new RegExp("(\s)([0-9]+)", "g");
    result = result.replace(reg, "<sup>$1</sup>$2");        
    
    result = this.replaceAll(result, " + ", "#");
    result = this.replaceAll(result, '+', "<sup>+</sup>");
    result = this.replaceAll(result, "#", " + ");
    result = this.replaceAll(result, " > ", " &#8594; ");
    


    
    return result;    
  },
  
  formatReaction : function(reaction){
    var result = "";
    var sep = " + ";

    var parts = reaction.split(' > ');
    
    var reactants = parts[0].split(sep);
    for( var i =0; i<reactants.length; i++){
      reactants[i] = pad(new Array(8).join(" "), reactants[i].trim() );
    }
    
    
    var products = parts[1].split(sep);
    for( var i =0; i<products.length; i++){
      products[i] = pad( new Array(8).join(" "), products[i].trim() );
    }   
    
    result = reactants.join(sep) + ' > ' + products.join(sep);  
    return pad(new Array(35).join(" "),result);
  },
    
  dataAsText : function(nodes, abundance){
    var reactions = {formation : '', destruction : ''};
    var node = null;
    var ids = nodes.getIds();
    var str = "";
    var pad_short = new Array(20).join(" ");
    var separator = " ";
    for (var i = 0, len = ids.length; i < len; i++ ){
      str = "";
      node = nodes.get(ids[i]);          
      if(node.metadata.reaction !== undefined && node.group !== node_groups.legend){      
        str += this.formatReaction(node.metadata.reaction)+separator;
        str += pad( pad_short, str_to_exp_number(node.metadata.abundances[0])+separator, true);
        str += pad( pad_short, str_to_exp_number(node.metadata.abundances[1])+separator, true);
        str += pad( pad_short, str_to_exp_number(get_constant(node.metadata.rate, node.metadata.abundances))+separator, true);
        str += pad( pad_short, str_to_exp_number(node.metadata.rate)+separator, true);
        str += pad( pad_short, str_to_exp_number(node.metadata.rate/abundance), true);
        str += pad( pad_short, node.metadata.proportion, true);
        str += "\n";
      }
      
      if(node.group === node_groups.product)
        reactions.destruction = reactions.destruction + str;
        
      if(node.group === node_groups.reactant)
        reactions.formation = reactions.formation + str;
    }
    
    
    var result = 'Formation reactions \n';
    
    result += pad( new Array(36).join(' '), 'A + B -> C + D'+separator )+ pad( pad_short, 'n(A) (cm-3)'+separator, true );
    result += pad( pad_short, 'n(B) (cm-3)'+separator, true ) + pad( pad_short, 'constant'+separator, true )+ pad( pad_short, 'rate (cm-3 s-1)'+separator, true ) + pad( pad_short, 'freq (s-1)', true ) + pad( pad_short, 'proportion', true );
    result += '\n'; 
    result += reactions.formation;
    result += '\n'; 
    result += 'Destruction reactions\n';    
    result += pad( new Array(36).join(' '), 'A + B -> C + D'+separator )+ pad( pad_short, 'n(A) (cm-3)'+separator, true );
    result += pad( pad_short, 'n(B) (cm-3)'+separator, true ) + pad( pad_short, 'constant'+separator, true )+ pad( pad_short, 'rate (cm-3 s-1)'+separator, true ) + pad( pad_short, 'freq (s-1)', true )+ pad( pad_short, 'proportion', true );
    result += '\n'; 
    result += reactions.destruction;
    result += '\n';  
    
    return result;    
    
  },
  data_as_json : function(){}
  
}

//
// Utility functions
//
function get_constant(rate, abundances){
  var result = parseFloat(rate);

  if (abundances[0] !== "N/A" && abundances[1] !== "N/A"){
    result =  (parseFloat(rate)/(parseFloat(abundances[0])*parseFloat(abundances[1]))).toExponential(2);  
  } else if (abundances[0] !== "N/A" && abundances[1] === "N/A"){
    result = (parseFloat(rate)/(parseFloat(abundances[0]))).toExponential(2);  
  } else if (abundances[0] === "N/A" && abundances[1] !== "N/A"){
    result = (parseFloat(rate)/(parseFloat(abundances[1]))).toExponential(2);  
  }
    
  return result;
}

function picture(canvas){
  var img    = canvas.toDataURL("image/png");
  document.write('<img src="'+img+'"/>');  
}

function map_to_string(map){
  var result = "";
  for (var key in map) {
    result = result + map[key];
  }  
  return result;
}

function get_percentage(rate){
  var value = parseFloat(rate.replace('%', ''));
  return value/100;
}

function str_to_exp_number(value){
    if ( value !== "N/A" )  
      return  parseFloat(value).toExponential(2);
      
    return "";
}

function datalist_to_html(sorted_data, data){
  var result = "";
  for(var i = 0, len = sorted_data.length; i< len; i++){
    result = result + "<p>";
    result = result +"<span class='title'>"+sorted_data[i]+"</span>";
    result = result +"<span class='value' title='"+data[sorted_data[i]].value.toExponential(4)+"'>"+data[sorted_data[i]].value.toExponential(1)+"</span>";
    result = result + "<span class='unit'>"+data[sorted_data[i]].unit+"</span>";
    result = result + "</p>";  
  }
  return  result + "";  
}

//
// Global variables
//
var page = {  
  thresholdInput : document.getElementById("threshold"),
  simpleThreshold : document.getElementById('minvalue_thresholdtype'),
  cumulatedThreshold : document.getElementById('cumulated_thresholdtype'),
  mainPosVarInput : document.getElementById("main_pos_var"),
  mainPosVarName : document.getElementById("main_pos_var_name"),
  mainPosVarManualInput : document.getElementById("main_pos_var_input_value"),
  controlPanel : document.getElementById('interface'),
  hasControlPanel : document.getElementById('with_control_panel'),
  controlPanelItem : document.getElementById('with_control_panel_item'),
  hasPhysicalConditionPanel : document.getElementById('with_physical_conditions'),
  physicalConditionPanelItem : document.getElementById('with_physical_conditions_item'),
  physicalConditionPanel : document.getElementById('gas_conditions'),
  mainPosVarGraphLegend : document.getElementById("plot_main_pos_var_legend"),
  speciesList : document.getElementById("species_list"),
  filesList : document.getElementById("files_list"),
  isNewFile : false,
  hasGasConditions : document.getElementById('with_gas_conditions'),
  gasConditionsDisplay : document.getElementById("gas_conditions_display"),  
  hasGrainConditions : document.getElementById('with_grain_conditions'),
  grainTemperaturesDisplay : document.getElementById("grain_temperatures_display"),
  grainDensitiesDisplay : document.getElementById("grain_densities_display"),
  hasReactionInformations : document.getElementById('with_reaction_informations'),
  infosDisplay : document.getElementById("infos"),
  exportPictureButton : document.getElementById("picture_export"),
  exportRatesDisplayedButton : document.getElementById("rates_displayed_export"),
  exportRatesMainPosVarTxtButton : document.getElementById("rates_main_pos_vars_txt_export"),
  exportNetworkButton : document.getElementById("network_export"),
  spinnerArea : document.getElementById('spinner_area'),
  graphDisabler : document.getElementById('graph_disabler'),
  backButton : document.getElementById("go_back"),
  allReactionsDisplay : document.getElementById("all_reactions_display"),
  savedFiles : document.getElementById("saved_files"),
  currentModel : document.getElementById("current_model_name"),
  isLogXPlot : document.getElementById('with_plot_x_log'),
  isLogYPlot : document.getElementById('with_plot_y_log'),
  spinner : null,
  allSpecies : [],
  allMainPosVars : [],
  
  updateMainPosVarInput : function(main_pos_vars){
    this.allMainPosVars = main_pos_vars;
    var str = '<option value="0">'+main_pos_vars[0]+'</option>';  
    this.mainPosVarInput.step = config.av_step;
    for(var i = config.main_pos_var_step, len = main_pos_vars.length; i < len; i = i+config.main_pos_var_step){ // display only 1/main_pos_var_step of values
      str = str + '<option value="'+i+'">'+main_pos_vars[i]+'</option>';
      page.mainPosVarInput.max = i; 
    }
    this.mainPosVarInput.innerHTML = str;       
  },
  
  updateGraphVisibility : function(threshold){
    if(threshold == 0 && this.simpleThreshold.checked === true){
      this.disableGraph();
    } else if(threshold == 100 && this.cumulatedThreshold.checked === true){
      this.disableGraph();
    } else{
      this.enableGraph();  
    }
  },
  
  disableGraph : function(){
    this.graphDisabler.style.visibility = "visible";
  },
  
  enableGraph : function(){
    this.graphDisabler.style.visibility = "hidden";
  },
  
  getThresholdType : function(){
    if(simpleThreshold.checked === true)
      return 
  },
  
  goBackSpecies : function(species){
    this.speciesList.selectedIndex = this.allSpecies[species];  
  },
  
  setModelName : function(name){
    this.currentModel.innerHTML = name;
  },
  
  setMainPositionVarName : function(name, unit){
    this.mainPosVarName.innerHTML = name + ' ( ' + unit + ' ) ';
  },
  
  getModelName : function(){
    return this.currentModel.innerHTML;
  },
  
  switchLeftPanel :function(source){
    if(source === this.hasControlPanel.id ){
      this.displayControlPanel(true);
      this.displayPhysicalConditionPanel(false);
    }else if(source === this.hasPhysicalConditionPanel.id){
      this.displayControlPanel(false);
      this.displayPhysicalConditionPanel(true);
    }

  },

  displayControlPanel : function(is_visible){
    if(is_visible === false){
      this.controlPanel.style.display = "none";
      this.controlPanelItem.className = "pure-menu-item";
    }else{
      this.controlPanel.style.display = "";
      this.controlPanelItem.className = "pure-menu-item pure-menu-selected";
    }
    
  },
  
  displayPhysicalConditionPanel : function(is_visible){
    if(this.physicalConditionPanel !== null){
      if(is_visible === false){
        this.physicalConditionPanel.style.display = "none";
        this.physicalConditionPanelItem.className = "pure-menu-item";
      }else{
        this.physicalConditionPanel.style.display = "";
        this.physicalConditionPanelItem.className = "pure-menu-item pure-menu-selected";
      }
    }
  },  
 
  displayGasConditions : function(){
    if(this.hasGasConditions.checked === false){
      this.gasConditionsDisplay.style.display  = "none";
    }else{
      this.gasConditionsDisplay.style.display = "";
    }   
  },
  
  displayGrainConditions : function(){
    if(this.hasGrainConditions.checked === false){
      this.grainTemperaturesDisplay.style.display  = "none";
      this.grainDensitiesDisplay.style.display  = "none";  
    }else{
      this.grainTemperaturesDisplay.style.display = "";
      this.grainDensitiesDisplay.style.display  = "";  
    }      
  },
  
  displayReactionInformations : function(){
    if(this.hasReactionInformations.checked == false){
      this.infosDisplay.style.display  = "none";
    }else{
      this.infosDisplay.style.display = "";
    }  
  },
  
  /**
   * updates list of available files
   * @param data array of file names
   */
  updateFilesList : function(data){
    var files_string = '<option value="" selected>-- None --</option>';
    for(var i = 0, len = data.length; i < len; i = i+1){ 
      files_string = files_string + '<option value="'+data[i]+'">'+data[i]+'</option>';
    }
    this.filesList.innerHTML = files_string;
  },
  
  /**
   * updates list of available species
   * @param data array of species names
   */
  updateSpeciesList : function(data){
    var species_string = "";//'<option selected>-- None --</option>';
    var selected = "";
    for(var i = 0, len = data.length; i < len; i++){
      if(data[i] === 'H2')
        selected += 'selected="selected"';
      species_string = species_string+'<option value="'+data[i]+'"'+ selected +'>'+data[i]+'</option>';
      page.allSpecies[data[i]] = i;
      if(selected != "")
        selected = "";      
    }
    this.speciesList.innerHTML = species_string;   
  },
  
  
  setInitialGasCondition : function(){
    if (this.hasGasConditions.checked == false)
      this.gasConditionsDisplay.style.display  = "none";
  },
  
  setInitialGrainConditions : function(){
    if (this.hasGrainConditions.checked == false){
      this.grainTemperaturesDisplay.style.display  = "none";  
      this.grainDensitiesDisplay.style.display  = "none";   
    }   
  },
  
  setInitialReactionInformations : function(){
    if (this.hasReactionInformations.checked == false)
      this.infosDisplay.style.display  = "none";  
  },
  
  updateGrainConditionsDisplay : function(temperatures, densities){  
    this.grainTemperaturesDisplay.innerHTML =  "<h2>Grain temperatures</h2>"+temperatures;
    this.grainDensitiesDisplay.innerHTML =  "<h2>Grain densities</h2>"+densities;
  },
  
  updateGasConditionsDisplay : function(content){
    this.gasConditionsDisplay.innerHTML = "<h2>Gas conditions</h2>"+content;
  },
  
  enableButtons : function(){
    this.exportPictureButton.disabled = "";
    this.exportRatesDisplayedButton.disabled = "";
    this.exportRatesMainPosVarTxtButton.disabled = "";
    this.exportNetworkButton.disabled = "";
  },
  
  toogleBackButton : function(status){
    if(status)
      this.backButton.disabled = "";
    else
      this.backButton.disabled = "disabled";
  },
  
  updateInfosDisplay : function(nodes, abundance){    
    this.infosDisplay.innerHTML = formatter.dataAsHtml(nodes, abundance);
  },
  
  updateSpeciesDisplay : function(species, abundance, abundance_unit){ 
    var html = "<p><strong>n("+species+") = </strong>";
    html = html + parseFloat(abundance).toExponential(2)+" "+abundance_unit+"</p>" ;
    document.getElementById("species_display").innerHTML = html;
  },
  
  highlightReaction : function(element_id){
    var current = document.getElementById(element_id);
    if(current !== null) 
      current.style.backgroundColor = mol_network.theme_color;  
  },
  
  initReactionsTable : function(){
    var table = document.getElementById("all_reactions_table");
    for( var i = 0, len = table.rows.length; i < len; i++)
      table.rows[i].style.backgroundColor = "white";
  },
  
  displayThresholdValue : function(){
    document.getElementById("threshold_value_display").innerHTML = this.thresholdInput.value+" %";
  },
  
  displayReactionRatesFile : function(url, filename){
    page.savedFiles.innerHTML = page.createLink(url +'/static/' , filename, 'Reaction rates');
  },
  
  displayFormationDestructionFiles : function(url, filename, formation_file, destruction_file){
    var content = '<p>' + page.createLink(url +'/static/' , formation_file, "Formation data") + '</p>';
    content += '<p>' + page.createLink(url +'/static/' , destruction_file, "Destruction data") + '</p>';
    page.savedFiles.innerHTML = content;
  },
  
  displayNetworkFile : function(url, filename){
    page.savedFiles.innerHTML = page.createLink(url +'/static/' , filename, "Network"); 
  },
  
  createLink : function(url, filename, linkText){
    return  '<a href="' + url + filename + '" target="_blank">' + linkText + '</a> ';
  },

  clearSavedFiles : function(){
    this.savedFiles.innerHTML = "";
  },
  
  showSpinner : function(){
    this.spinnerArea.style.visibility = "visible";
    this.spinner = new Spinner(opts).spin(this.spinnerArea);      
   
  },
  
  hideSpinner : function(){
    this.spinner.stop();
    this.spinnerArea.style.visibility = "hidden";
  },
  
  triggerEvent : function(element, event){
    if (document.createEvent) {
      element.dispatchEvent(event);
    } else {
      element.fireEvent("on" + event.eventType, event);
    }
  }  

}
