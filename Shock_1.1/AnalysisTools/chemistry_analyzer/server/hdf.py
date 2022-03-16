'''
Notes :

PS: Table Metadata :
 Metadata contains :
PathHDF5 | DataSet | Col.Index | PubID | HName  | DataType | Unit | SKOS | UCD | utype | Description | ObjParent | Group |
   0     |    1    |     2     |   3   |   4    |    5     |  6   |  7   |  8  |  9    |   10        |    11     |  12   |
   
 Metadata_ObjType contains :
 PubID | Hname | SKOS | Description

'''
import sys
import string as string
import h5py as h5
import numpy as numpy
import time


__author__ = "Nicolas Moreau"
__copyright__ = "Copyright 2012"
__credits__ = ["Nicolas Moreau"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Nicolas Moreau"
__email__ = "nicolas.moreau@obspm.fr"
__status__ = "Development"  

 
class PdrHDF(object):
  """
      Access to a pdr HDF5 file
  """
  def __init__(self, filename):
    
    self.chemistry_path = '/Local quantities/Chemistry'
    self.conditions_path = '/Local quantities/Gas conditions'
    self.position_path = '/Local quantities/Positions'
    self.masses_path = '/Parameters/Species mass'
    self.network_path = '/Parameters/Network'
    
    #hdf file    
    self.f =  h5.File(filename, 'r')    

    self.destructions = None
    self.formations = None
    self.gasConditions = None
    self.densities = None
    self.masses = None
    self.network = None
    self.positionVar = None
      
  def close(self):
    """
        close hdf file
    """
    self.f.close()
    
  def __buildSet(self):
    
    self.formations = {}
    self.destructions = {}
    self.gasConditions = {}
    self.densities = {}
    self.masses = {}
    self.metadata = self.f['/Metadata/Metadata']
    
    
    i = 0
    # build list of reactions and quantities
    for value in self.metadata :
      if value[1].startswith('Destruction rate'):
        if value[1] not in self.destructions :
          self.destructions[value[1]] = []
        self.destructions[value[1]].append(value)

      elif value[1].startswith('Formation rate'):
        if value[1] not in self.formations :
          self.formations[value[1]] = []
        self.formations[value[1]].append(value)  
        
      elif value[1] == "Gas conditions" : 
        self.gasConditions[str(value[3])] = value.tolist()   
        
      elif value[1] == "Positions" : 
        if self.positionVar is None : 
          self.positionVar = {'name' : value[4], 'unit' : value[6], 'short_name' : value[3]}
        
      elif value[1].startswith("Densities"):         
        self.densities[str(value[4])] = value.tolist() 
        
      elif "mass_of" in value[3] :
        self.masses[value[1].replace("Mass ", "")] = float(self.f[value[0]+'/'+value[1]][0])     
        
    self.getNetwork()
    
  def getAllSpecies(self):
    return self.f[self.chemistry_path].keys()  
    
  def getGetConditionsValues(self, main_pos_vars, quantity_name):
    if self.gasConditions is None :
      self.__buildSet()
      
    result = []
    for key in sorted(self.gasConditions.keys()) :
      if key == quantity_name :        
        value = self.gasConditions[key]
        for i in range(0,len(main_pos_vars)) :
          result.append(self.f[str(value[0])+'/'+str(value[1])][i].tolist()[int(value[2])])
          
    return result

    
  def getTemperaturesMainPosVars(self):
    """
      returns values that stay constant for a given file :
        main_pos_var and temperatures
        
      returns dict[main_pos_vars, temperatures]
    
    """
    result = {}
    main_pos_vars = []
    
    for value in self.f[self.position_path]:
      main_pos_vars.append(value[0])    
      
    result['temperatures'] = self.getGetConditionsValues(main_pos_vars, 'tgas')
    result['main_pos_vars'] = main_pos_vars
    return result
  
  def getFormationBySpecies(self, species):
    if self.formations is None : 
      self.__buildSet()
      
    # field can be missing in shock hdf file
    try :
      return self.formations['Formation rates '+species]
    except KeyError :
      return []
    
  def getDestructionBySpecies(self, species):
    if self.destructions is None :
      self.__buildSet()
      
    # field can be missing in shock hdf file
    try :
      return self.destructions['Destruction rates '+species]  
    except KeyError:      
      return []
    
  def getData(self, meta, main_pos_var):
    data = self.f[meta[0]+'/'+meta[1]]
    return data[main_pos_var][int(meta[2])]
    
  def getGasConditions(self, main_pos_var):
    if self.gasConditions is None :
      self.__buildSet()
      
    result = {}    
    for key in sorted(self.gasConditions.keys()) :
      value = self.gasConditions[key]
      entry = {}
      entry['name'] = value[4]
      
      if value[6] != "no unit" :
        entry['unit'] = value[6]
      else : 
        entry['unit'] = ''
        
      # numpy has its own list type, integer stored as string
      entry['value'] = self.f[str(value[0])+'/'+str(value[1])][main_pos_var].tolist()[int(value[2])]
      result[value[4]] = entry
    return result
    
  """
    concatenate content of Network table into a string and return it
  """
  def getNetwork(self):
    if self.network is None : 
      result = []
      try :
        for line in self.f[self.network_path]['Network'][0] :
          result.append(line+"\n")
        self.network = ''.join(result)
      except KeyError as e:
        print "Network dataset is not defined in hdf5 file"
        print e
        self.network = 'Network is not available in the HDF5 file'
      
    return self.network
    
  def getAbundance(self, species, main_pos_var):
    return self.f[self.chemistry_path+'/'+species+'/Densities '+species][main_pos_var][0]
    
  def getAbundances(self, species):
    result = []
    for value in self.f[self.chemistry_path+'/'+species+'/Densities '+species][:] :
      result.append(value[0])
    return result
    
  # if abundance is always cm-3, this should be removed to save memory and iterations
  def getAbundanceUnit(self, species):
    #return "cm-3"
    if self.densities is None :
      self.__buildSet() 
      
    #return self.densities["dens_spec_"+species.lower().replace('+', 'p')][6]
    return self.densities["n("+species+")"][6]  
  
        
  def getTotalDestructionRate(self, species, main_pos_var):
    return self.f[self.chemistry_path+'/'+species+'/Total dest rates '+species][main_pos_var][0]
    
  def getTotalFormationRate(self, species, main_pos_var):
    return self.f[self.chemistry_path+'/'+species+'/Total form rates '+species][main_pos_var][0]
    
  def getTotalFormationRates(self, species, main_pos_var):
    return self.f[self.chemistry_path+'/'+species+'/Total form rates '+species][main_pos_var][0]    
    
  def getFormationRates(self, species):
    return self.f[self.chemistry_path+'/'+species+'/Formation rates '+species]
    
  def getDestructionRates(self, species):
    return self.f[self.chemistry_path+'/'+species+'/Destruction rates '+species]
    
  def getMass(self, species):
    try :
      return self.masses[species]
    except KeyError :   
      return 0

  def getMasses(self):
    return self.masses


  def getPositionVarName(self):
    if self.positionVar is None :
      self.__buildSet() 
    
    return self.positionVar
