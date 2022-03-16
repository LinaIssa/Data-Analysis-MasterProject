#-*- coding: utf-8 -*-
import os
import json
import collections
import uuid
from bottle import route, run, get, request, view, default_app, BaseTemplate, response, abort, static_file, BaseRequest
from analyzer import Analyzer



""" enable cross domain ajax requests

"""
def enable_cors(fn):
  def _enable_cors(*args, **kwargs):
    # set CORS headers
    response.headers['Access-Control-Allow-Origin'] = '*'
    response.headers['Access-Control-Allow-Methods'] = 'GET, POST, PUT, OPTIONS'
    response.headers['Access-Control-Allow-Headers'] = 'Origin, Accept, Content-Type, X-Requested-With, X-CSRF-Token'

    if request.method != 'OPTIONS':
        # actual request; reply with the actual response
        return fn(*args, **kwargs)

  return _enable_cors
  
"""  creates directory if it does not exist

:param directory: name of a directory to create if it does not exist
:returns: name of directory 
"""
def check_directory(directory):
  if os.path.exists(directory) is False : 
    try : 
      os.mkdir(directory)
    except Exception as e: 
      print("error while creating '%s' directory"%directory)
      print(e)
      
  return directory

"""create a file a write a string in it

:param filename: name of the file to create, if it exists it will be replaced
:param string_data: content to write in the file

:returns: True if the written file content is identical to the string parameter, or False if not

"""
def save_data(filename, string_data):
  f = open(save_directory+ "/" + filename, "w")
  f.write(string_data)
  f.close()
  
  f = open(save_directory+ "/" + filename, "r")
  content = f.read()
  f.close()
  
  #check content of written file
  if content == string_data :
    return True
    
  return False

"""writes a file containing rates for all main_pos_vars for a species

:param main_pos_var:  name of position variable to navigate data (
                      AV or position )
:param data:          a dictionnary contaning data
                      data['reactions'] list of reactions
                      data['main_pos_vars'][main_pos_var_value][reaction_name][rate] values
            
:returns: a string containing the formatted data

"""
def build_main_pos_var_export(main_pos_var, data):
  col_index = 1
  string = "#[%s] %s     "% (col_index, main_pos_var)  
  
  for reaction in data['reactions'] :
    col_index += 1
    string += "[%s] %s      "%(col_index, reaction)      
  string += "\n#[1]      "
  
  col_index = 1
  for unit in data['units'] :
    col_index += 1
    string += "[%s] %s      "%(col_index, unit)      
  string += "\n#"
  
  
  for i in range(col_index):
    string += ("[%s]"%str(i+1)).rjust(20)
  string += "\n"
  for key  in collections.OrderedDict(sorted(data["main_pos_vars"].items())) :
    string += ("%2.5e"%key).rjust(20)
    for reaction in data['reactions']:
      if reaction in data['main_pos_vars'][key] :
        string += ("%2.10e"%data['main_pos_vars'][key][reaction]['value']).rjust(20)
      else :
        string += "_".rjust(20)
        
    string += "\n"
    
  return string


"""
 Manages Analyzer object instances
 Each instance contains one hdf file
"""
class AnalyzerManager():
  
  def __init__(self):
    self.analyzers = {}
    self.currentAnalyzer = None
  
  """ returns Analyzer corresponding to the given client id
   :param client_id: id of the client process
   :returns: Azalyzer object
  """
  def getAnalyzer(self, client_id):
    return self.analyzers[client_id]
    
  """ creates a new Analyzer object for the given client id
  :param client_id: id of the client process
  :param filename: name of the file to open   
  """
  def addAnalyzer(self, client_id, filename):
    self.analyzers[client_id] = self.instantiateAnalyzer(filename)
    
  """returns a Analyzer object for the hdf5 file filename

  :param filename: name of hdf5 file
  :returns: Analyzer object
  """
  def instantiateAnalyzer(self, filename):
    analyzer = Analyzer(data_directory)
    analyzer.loadHdf(filename)
    return  analyzer


"""
    urls
"""
base_url = "/chemistry-analyzer"

""" writes locally rates for all reactions of given species at all main_pos_vars if at least one rate is above threshold and
returns a dictionary containing the name of the generated files
:param session: id of the client process
:param species: species name
:param threshold: threshold value 
:param filename: suffix for file name

:returns: a dictionary of file names

"""
@route( base_url+'/session/<session>/saverates/species/<species>/threshold/<threshold>/filename/<filename>', name='ratespecies')
@enable_cors
def rate_species(session, species, threshold, filename):
  global analyzerManager
  main_pos_var = analyzerManager.getAnalyzer(session).getPositionVarName()
  data = analyzerManager.getAnalyzer(session).getAllReactions(species, float(threshold))
  save_data("formation_"+filename, build_main_pos_var_export(main_pos_var['short_name'], data['formation']))
  save_data("destruction_"+filename, build_main_pos_var_export(main_pos_var['short_name'], data['destruction']))
  return json.dumps({"formation":"formation_"+filename, "destruction":"destruction_"+filename})
  
""" generates and returns an id for this session to identify connection with client

:returns: random uuid

"""
@route( base_url+'/opensession', name='opensession')
@enable_cors
def open_session():    
  return json.dumps({'sessionid' : uuid.uuid4().hex})


""" returns list of reactions for a species at a given main_pos_var
:param session: id of the client process
:param main_pos_var: main_pos_var value
:param species: species name

:returns: a dictionary of reactions

"""
@route( base_url+'/session/<session>/mainposvar/<main_pos_var:int>/species/<species>', name='mainposvarspecies')
@enable_cors
def main_pos_var_species(session, main_pos_var, species):  
  global analyzerManager
  return analyzerManager.getAnalyzer(session).getReactions(species, main_pos_var)
  
""" returns abundances of the given species at all main_pos_vars
:param session: id of the client process
:param species: species name

:returns: a list of abundance at each main_pos_var

"""
@route( base_url+'/session/<session>/abundances/species/<species>', name='abundancesspecies')
@enable_cors
def abundances(session, species):  
  global analyzerManager
  return json.dumps(analyzerManager.getAnalyzer(session).getAbundances(species))
  
""" returns network in given hdf file
:returns: string

"""
@route( base_url+'/session/<session>/network', name='network')
@enable_cors
def network(session):  
  global analyzerManager
  return json.dumps(analyzerManager.getAnalyzer(session).getNetwork())
  
""" returns network in given hdf file
:returns: string

"""
@route( base_url+'/session/<session>/positionvarname', name='positionvarname')
@enable_cors
def positionvarname(session):  
  global analyzerManager
  return json.dumps(analyzerManager.getAnalyzer(session).getPositionVarName())
  
""" returns gas conditions at given main_pos_var value
:param session: id of the client process
:param main_pos_var: main_pos_var value

:returns: a dictionary containing gas conditions
"""
@route( base_url+'/session/<session>/gasconditions/mainposvar/<main_pos_var:int>', name='gasconditions')
@enable_cors
def gas_conditions(session, main_pos_var):  
  global analyzerManager  
  return json.dumps(analyzerManager.getAnalyzer(session).getGasConditions(main_pos_var))

""" returns a list of all main_pos_var positions
:param session: id of the client process

:returns: a list of values
"""  
@route( base_url+'/session/<session>/allmainposvars', name='main_pos_vars')
@enable_cors
def mainposvars(session):  
  global analyzerManager
  return json.dumps(analyzerManager.getAnalyzer(session).getTemperaturesMainPosVars())
  

""" returns a list of all species
:param session: id of the client process

:returns: a list of species
"""  
@route( base_url+'/session/<session>/allspecies', name='species')
@enable_cors
def species(session):  
  global analyzerManager
  return json.dumps(analyzerManager.getAnalyzer(session).getSpecies())
  

""" returns a list of files contained in the data directory (file name must contain _c string)

:returns: a list of file names
"""  
@route( base_url+'/allfiles', name='files')
@enable_cors
def files(): 
  result = []
  for entry in os.listdir(data_directory):
    if "_c" in entry :
      result.append(entry)
  return json.dumps(result)
  

""" loads a new file
:param session: id of the client process
:param filename: name of the file to load

:returns: 200 status if file is loaded, 204 if not found

"""  
@route( base_url+'/session/<session>/changefile/file/<filename>', name='changefile')
@enable_cors
def change_file(session, filename): 
  global analyzerManager
  
  try : 
    analyzerManager.addAnalyzer(session, filename)
    
    if analyzerManager.getAnalyzer(session).hdf is not None : 
      response.status = 200
      return {}
    else : 
      response.status = 500  
  except IOError as e :
    response.status = 204
    
""" returns the content of a file into the save directory
:param filename: name of the file to read

:returns: file content, 404 error if not found
"""
@route(base_url+'/session/<session>/positionvar')
def positionvar():
  global analyzerManager
  return json.dumps(analyzerManager.getAnalyzer(session).getMainPositionVariable())


""" saves a string into a file in save directory
:param filename: name of the file to write
:param string_data: content to write

:returns: an empty list if the file has been written, HTTP 500 error if a problem occured

"""
@route( base_url+'/savefile', name='savefile', method="POST")
@enable_cors
def save_file(): 

  filename = request.forms.get('filename')
  string_data = request.forms.get('string_data')  
  write_status = save_data(filename, string_data)
  
  #check content of written file
  if write_status is True :
    response.status = 200
    return []
  else : 
    abort(500, "Server interal error : File was not correctly written.")
    
    
""" checks if a file already exists in the save directory
:param filename: name of the file to test

:returns: True if the file exists, False if not
"""
@route( base_url+'/existfile/<filename>', name='existfile', method="GET")
@enable_cors
def file_exists(filename):   
  result = {'status' : os.path.isfile(save_directory+'/'+filename) }
  return result
  
  
     
""" returns the content of a file into the save directory
:param filename: name of the file to read

:returns: file content, 404 error if not found
"""
@route(base_url+'/static/<filename:path>')
def static_content(filename):
    return static_file(filename, root='../save/')
    

""" home page

:returns: string message
"""
@route( base_url , name='home', method="GET")
@enable_cors
def home():   
  result = "chemistry analyzer server default page" 
  return result
  
def runserver(port, debug):
  app = default_app()
  BaseTemplate.defaults['get_url'] = app.get_url
  run(host='localhost', port=port, debug=debug)

analyzerManager = AnalyzerManager()

# increase max post size to 8 Mb to save big files
BaseRequest.MEMFILE_MAX = 8 * 1024 * 1024 
data_directory = check_directory("../data")
save_directory = check_directory("../save")
runserver(12345, True)
