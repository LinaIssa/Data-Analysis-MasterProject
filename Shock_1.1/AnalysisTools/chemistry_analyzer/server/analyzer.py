#-*- coding: utf-8 -*-
import sys
import os.path
import json
import collections
import time
import hdf as hdf_lib


"""Object providing an interface to extract data from a hdf5 file
"""
class Analyzer(object):
  def __init__(self, data_directory):
    
    # check data directory exists
    if os.path.isdir(data_directory+"") is False : 
      raise Exception("Directory %s does not exist" % data_directory)
    else : 
      self.data_directory = data_directory
    
    self.__complete_directory_path()
    
    # hdf object
    self.hdf = None #hdf_lib.PdrHDF(hdf_path+"/"+hdf_filename)
    
    # list of all main_pos_vars contained in the file
    self.mainPosVars = None
  
  """
    append trailing / if necesary  to data directory
  """
  def __complete_directory_path(self):      
    if self.data_directory[len(self.data_directory)-1] != "/" :
      self.data_directory = self.data_directory + "/"
      
  """loads a hdf file
  :param hdf_path: path of file
  """
  def loadHdf(self, hdf_file):    
    self.hdf = hdf_lib.PdrHDF(self.data_directory  + hdf_file)
    self.__checkHdf()

  """returns a dictionary containing data to be displayed by analyzer web page
  
  :param reaction: string describing the reaction
  :param rate: rate of this reaction
  :param proportion: part of this reaction in the formation/destruction of the species at given main_pos_var
  :param species: considered species
  :param main_pos_var: considered main_pos_var value
  :param unit: rate unit
  
  :returns:  dictionary of data  
  """
  def __build_data_dictionary(self, reaction, rate, proportion, species, main_pos_var, unit):
    result = {}      
    reaction_string = reaction[4].split("|")[1]  

    parts = reaction_string.split(">")  
    
    products = parts[1].strip().split(" + ")  
    reactants = parts[0].strip().split(" + ")  
    
    result["products"] = []
    result["reactants"] = []
    
    for product in products :
      result["products"].append({"name" : product, "mass" : self.hdf.getMass(product)})
      
    for reactant in reactants :
      result["reactants"].append({"name" : reactant, "mass" :  self.hdf.getMass(reactant)})   
    
    
    result["reactants_abundances"] = []
    
    for sp in result["reactants"] :
      try : 
        result["reactants_abundances"].append(self.hdf.getAbundance(sp['name'], main_pos_var))
      except Exception  :
        result["reactants_abundances"].append("N/A")
        
    result["reactants_abundance_units"] = []
    
    for sp in result["reactants"] :
      try : 
        result["reactants_abundance_units"].append(self.hdf.getAbundanceUnit(sp['name']))
      except Exception  :
        result["reactants_abundance_units"].append("N/A")
      
    result["rate"] = rate
    result["proportion"] =  "{0:.4g} %".format(proportion) #  format(str(proportion))+"%"
    result["reaction"] = reaction_string
    
    result["unit"] = unit
      
    return result
    
  """check that hdf attribute is not empty
  
  """
  def __checkHdf(self):
    if self.hdf is None :
      raise Exception("No HDF5 data loaded")
      
  
  """ returns a dictionary of rate for all reactions at all main_pos_vars for a given species
  :param sp: considered species
  :param threshold: threshold above which we keep the reaction
  :param reactions: list of reactions
  :param main_pos_vars: list of main_pos_var values
  :param rates: list of rates
  :param rateGetter: function returning cumulated rate values at a a given main_pos_var, to calculate reaction proportion (lambda expr)
  
  :returns: dictionary
  """
  def __getAllReactions(self, sp, threshold, reactions, main_pos_vars, rates, rateGetter):
    result = {}
    total_rates = {}
    
    result["reactions"] = []
    result["units"] = []
    result["main_pos_vars"] = {}
    
    reaction_buffer = []
    
    for reaction in reactions :       
      current_reaction = reaction[4]
      values = []
      for i in range(0, len(main_pos_vars)) : 
        rate = rates[i][int(reaction[2])] 
        
        if i not in total_rates : 
          total_rates[i] = rateGetter(sp, i)

        proportion =  rate/total_rates[i]
        
        
        if proportion > threshold : 
          reaction_buffer.append(reaction)
          break
              
    for reaction in reaction_buffer :
      current_reaction = reaction[4]
      values = []
      for i in range(0, len(main_pos_vars)) : 
        rate = rates[i][int(reaction[2])] 
        
        if i not in total_rates : 
          total_rates[i] = rateGetter(sp, i)

        proportion =  rate/total_rates[i]      
        if i not in result["main_pos_vars"] :
          result["main_pos_vars"][i] = {}          
                    
        result["main_pos_vars"][i][current_reaction] = {"value" : rate, "proportion" : proportion}
        if current_reaction not in result['reactions']:
          result['reactions'].append(current_reaction)
          result['units'].append(reaction[6])             
    
    return result
  
  """returns a dictionary of all reactions (formation and destruction, all main_pos_vars) and rates for a given species at a given threshold
  :param sp: considered species
  :param threshold: considered threshold
  
  :returns: dictionary
  
  """
  def getAllReactions(self, sp, threshold):
    
    main_pos_vars = self.hdf.getTemperaturesMainPosVars()['main_pos_vars']
    result = {}   
    
    #pass function to obtain rates as a parameter
    getter = lambda sp, i : self.hdf.getTotalFormationRate(sp, i)    
    result["formation"] = self.__getAllReactions(sp, threshold, self.hdf.getFormationBySpecies(sp), main_pos_vars, self.hdf.getFormationRates(sp), getter)
    getter = lambda sp, i : self.hdf.getTotalDestructionRate(sp, i)    
    result["destruction"] = self.__getAllReactions(sp, threshold, self.hdf.getDestructionBySpecies(sp), main_pos_vars, self.hdf.getDestructionRates(sp), getter)
    
    return result

  """returns a dictionary of reactions data, indexed on the reaction proportion
  
  :param total_rate: some of all rates for the reactions at the given main_pos_vars and species
  :param sp: considered species
  :param main_pos_vars: considered main_pos_vars value
  :param reactions: list of reactions
  
  :returns: dictionary
  """
  def __getReactions(self, total_rate, sp, main_pos_var, reactions):    
    result = {}    
    i = 0  
    for reaction in reactions : # dict of metadata      
      rate = self.hdf.getData(reaction, main_pos_var)
      proportion = (rate/total_rate)*100      
      #~ if proportion > 0 : 
      entry = self.__build_data_dictionary(reaction, rate, proportion, sp, main_pos_var, reaction[6])
      key = "%s#%s"%(proportion,i) # add a incremented counter for equal proportions
      result[key] = entry
      i = i + 1
     
    # sort results after removing the counter
    result = collections.OrderedDict(sorted(result.items(), key=lambda x:float(x[0].split('#')[0]), reverse=True))
    return result.values()
    
  """returns a dictionary of reactions indexed by species for a given main_pos_var
  
  :param sp: considered species
  :param main_pos_var_index: index of the considered main_pos_var value in value list
  
  :returns: dictionary
  
  """
  def getReactions(self, sp, main_pos_var_index): 
    result = {}
    result[sp] = {}
    result[sp]['abundances'] = self.hdf.getAbundances(sp) # all densities
    result[sp]['abundance'] = self.hdf.getAbundance(sp, main_pos_var_index) # density for given main_pos_var value
    result[sp]['abundance_unit'] = self.hdf.getAbundanceUnit(sp)
    result[sp]['creation'] = self.__getReactions(self.hdf.getTotalFormationRate(sp, main_pos_var_index), sp, main_pos_var_index, self.hdf.getFormationBySpecies(sp))
    result[sp]['destruction'] = self.__getReactions(self.hdf.getTotalDestructionRate(sp, main_pos_var_index), sp, main_pos_var_index, self.hdf.getDestructionBySpecies(sp)) 
    return result


  """ returns gas conditions info at given main_pos_var index (index of main_pos_var value in main_pos_var values list)
  
  :param main_pos_var_index: main_pos_var value index
  
  :returns: list  
  """
  def getGasConditions(self, main_pos_var_index):
    return self.hdf.getGasConditions(main_pos_var_index)
        
  """ returns a list of all species in hdf file
  
  :returns: list
  
  """
  def getSpecies(self):  
    return self.hdf.getAllSpecies()
    
  """ returns a list of all main_pos_vars and temperatures values in hdf file
  
  :returns: dict[main_pos_var, temperatures]
  """
  def getTemperaturesMainPosVars(self):
    if self.mainPosVars is None :
      self.mainPosVars = self.hdf.getTemperaturesMainPosVars()
      
    return self.mainPosVars
    
  """ returns abundances of species in the hdf file
  
  :param sp: considered species
  
  :returns: list
  
  """
  def getAbundances(self, sp):
    return self.hdf.getAbundances(sp)
    
  """ return Network table
  
  """
  def getNetwork(self):
    return self.hdf.getNetwork()  
    
  """ return Network table
  
  """
  def getPositionVarName(self):
    return self.hdf.getPositionVarName() 


