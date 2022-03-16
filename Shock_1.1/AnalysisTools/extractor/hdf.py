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
import votable as vot
import numpy as numpy
import config as cfg
import convert as conv
import util


__author__ = "Nicolas Moreau"
__copyright__ = "Copyright 2012"
__credits__ = ["Nicolas Moreau"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Nicolas Moreau"
__email__ = "nicolas.moreau@obspm.fr"
__status__ = "Development"  


class PointerFactory(object):
  """
  creates a Pointer object
  """
  def __init__(self):
    self.pointer = Pointer()
      
  def withPath(self, path):
    self.pointer.path = path
    return self

  def withDataset(self, dataset):
    self.pointer.dataset = dataset
    return self        
      
  def withName(self, name):
    self.pointer.name = name
    return self
      
  def withColumn(self, column):
    self.pointer.column = column
    return self
      
  def withUnit(self, unit):
    self.pointer.unit = unit
    return self
      
  def withUcd(self, ucd):
    self.pointer.ucd = ucd
    return self
      
  def withSkos(self, skos):
    self.pointer.skos = skos
    return self
      
  def withUtype(self, utype):
    self.pointer.utype = utype
    return self
      
  def withDescription(self, description):
    description = description.strip()
    if description != '' and description != "none":
      self.pointer.description = description.strip()
    return self
      
  def withObjParent(self, parent):
    self.pointer.objParent = parent
    return self
      
  def withGroup(self, group):
    self.pointer.group = group
    return self
      
  def getPointer(self):
    return self.pointer
    
  def getErrorPointer(self, quantity):
    self.pointer = Pointer()
    self.pointer.name = quantity
    self.pointer.unit = "N/A"
    self.pointer.description = "error"
    return  self.pointer
    
    
class Pointer(object):
    """
        Points on data in a dataset
    """
    def __init__(self):
        #dataset path
        self.path = None
        #dataset name
        self.dataset = None
        #quantity human readable name
        self.name = None
        #column position in dataset
        self.column = None
        self.ucd = None
        self.skos = None
        self.utype = None
        self.unit = None
        self.description = None
        self.objParent = None
        self.group = None

  
class PdrHDF(object):
  """
      Access to a pdr HDF5 file
  """
  def __init__(self, filename):
    #hdf file
    self.f =  h5.File(filename, 'r')
    #create a pivoted array with the metadata
    # ------------------------------------------------------------------------------------
    #The old versions of numpy do not treat automatically the fact that the metadata
    #table is composed of string of variable lengths. Beside updating numpy, another
    #solution consists in setting the type to "object" when creating the numpy array
    #see http://stackoverflow.com/questions/16037824/how-to-convert-numpy-object-array-into-str-unicode-array
    #note: setting the type to unicode or post-converting object type via astype('U')
    #      do not work either with the old version of numpy
    # ------------------------------------------------------------------------------------
    self.metadata = numpy.array(self.f['/Metadata/Metadata'],dtype=object)
    #separator in path
    self.separator = "#"
    #indexes all the data
    self.index, self.quantities = self.__buildIndex()

      
  def close(self):
    """
        close hdf file
    """
    self.f.close()
      
  def getByGroup(self, group):
    """
        returns the object located at group position
    """
    try : 
      return self.f[group]
    except KeyError:
      print "group : " +group+ " not found"       
          
  def getDefaultColumns(self):
    result = []
    for key in  self.index['Positions']:
      result.append(key)
    return result

  def __buildIndex(self):
    """
        creates the index
    """        
    datasets = {} # map of Pointer objects
    quantities = [] # list of quantities with their complete path in tree
    
    for line in self.metadata:
      quantities.append(line[0].replace('/', self.separator)+self.separator+line[1]+self.separator+line[4])  
          
      if line[1] not in datasets : 
        datasets[line[1]] = {}   
      pf = PointerFactory()
      pf.withPath(line[0]).withDataset(line[1])
      pf.withName(line[4]).withColumn(line[2])
      pf.withUnit(line[6]).withUcd(line[8])
      pf.withSkos(line[7]).withUtype(line[9])
      pf.withDescription(line[10]).withObjParent(line[11])
      pf.withGroup(line[12]).getPointer()
      
      # pointer is identified by the human readable name of a quantity
      datasets[line[4]] = pf.getPointer()
    quantities.sort()
    return datasets, quantities
      
  def extractQuantity(self, p):  
    """
        extract a quantity according to the pointer p
    """
    ds = self.f[p.path+'/'+p.dataset]
    #1D array, returned directly
    if len(ds.shape) == 1:
      return ds
    else :        
      return ds[:,int(p.column)]
 
class Writer(object):
  """
      write extracted data in a file
  """
  def __init__(self, filename, header=None, separator=',', isdouble=False):
    #output file name
    self.filename = filename
    #a list of pointers
    self.pointers = []
    #a list of columns
    self.columns = []
    #separator between values
    self.separator = separator
    #
    self.header = header
    # precision is .6 if False, .15 if True
    self.isDouble = isdouble
      
      
      
  def addColumn(self, pointer, column):
    """
        add a column of data and its pointer
    """
    self.pointers.append(pointer)
    self.columns.append(column)
      
  def clearColumns(self):
    self.pointers = []
    self.columns = []
    
  def __getEmptyData(self, length):
    result = []
    for i in range(0, length):
      result.append('Not_In_HDF5')
      
    return result
    
  def __getArraySize(self):
    length = 0
    for i in range(0, len(self.columns) ):  
      if len(self.columns[i]) > 0:
        length = len(self.columns[i])
        break    
    return length
      
  def writePlotConfigFile(self):
    """
        write all data
    """        
    try:
      path = self.filename.split('.')
      
      f = open(path[0]+'_plot_xml.xml', 'w')
      result = '<?xml version="1.0" encoding="UTF-8"?>\n'
      result += '<descriptor readHeader="true">\n'
      result += '\t<data>\n'
      result += '\t\t<separator>'+self.separator+'</separator>\n'
      result += '\t\t<line label="data" ID="data">\n'
      for idx, pointer in enumerate(self.pointers) : 
        result += '\t\t\t<quantity>\n'
        result += '\t\t\t\t<name>'+pointer.name+'</name>\n'
        result += '\t\t\t\t<unit>'+pointer.unit+'</unit>\n'
        result += '\t\t\t</quantity>\n'
        
      result += '\t\t</line>\n'
      result += '\t</data>\n'
      result += '</descriptor>\n'

      f.write(result)
      f.close()
        
    except Exception as e:
      raise e

      
  def writeText(self):
    """
        write all data
    """        
    try:
      precision = "float"
      #~ if self.isDouble:
        #~ precision = "double"
      self.__columnsAreWritable()
      f = open(self.filename, 'w')
      header = ''
      data = []
            
      widths = []
      length = self.__getArraySize()
       
      for i in range(0, len(self.columns)) :
        if len(self.columns[i]) == 0:
          self.columns[i] = self.__getEmptyData(length)
          
      for i in range(0, length ): 
        values = [column[i] for column in self.columns]   
        towrite = []
        for value in values :
          if isinstance(value,numpy.string_) is False and isinstance(value, basestring) is False:
            value = conv.getValue(value, precision)          
          widths.append(len(value))
          towrite.append(value)          
        data.append(self.separator.join(towrite))
        
              
      if length > 1 :
        # common case, several lines, several columns
        if self.header is not None :
          header += "".join(self.header)
        line1 = []
        line2 = []
        line3 = []
        for idx, pointer in enumerate(self.pointers) : 
          line1.append('['+str(idx+1)+']'+pointer.name)
          line2.append('['+str(idx+1)+']'+pointer.unit)
          line3.append(string.rjust('['+str(idx+1)+']', widths[idx]))
        header = "#"+self.separator.join(line1) + "\n#" +self.separator.join(line2) + "\n#" + self.separator.join(line3)[1:] + "\n"      
        
        f.write(header)
        f.write("\n".join(data))
      else : 
        # file with only one line and several columns are plotted as
        # one column with several lines
        f.write("#%s%s       %s    %s\n"%("index".rjust(6), "value".rjust(35), 'unit'.rjust(12), "quantity" ))
        for idx, value in enumerate(data[0].split(self.separator)):
          f.write("%s%s   #   %s     %s\n"%(str(idx).rjust(7), str(value).rjust(35), self.pointers[idx].unit.rjust(11), self.pointers[idx].name ))

      f.close()
    except Exception as e:
      raise e
       
        
  def writeVoTable(self):
    try:
      self.__columnsAreWritable()
      table = vot.Votable()      
      datatype='float'
      if self.isDouble : 
        datatype = 'double'

      for i in range(0,len(self.columns)):
        builder = vot.FieldBuilder()
        field = self.pointers[i]
        link = vot.LinkBuilder().withContentRole('type').withHref(field.skos).getLink()
        table.addField(builder.withName(field.name).withUnit(field.unit).withUcd(field.ucd).withUtype(field.utype).withDatatype(datatype).withLink(link).getField())
        table.addColumn(self.columns[i])
      table.toFile(self.filename)           

    except Exception as e:
      raise e   
   
  def __columnsAreWritable(self):
    """
        check data validity
    """
    if len(self.columns) != len(self.pointers):
        raise Exception("columns and pointers do not match")
        
    if len(self.columns) > 0 :
      size = self.__getArraySize()
      for i in range(1, len(self.columns)):
        if len(self.columns[i]) > 0 and len(self.columns[i]) != size:
          raise Exception('all the columns do not have the same size')
      return True        
    raise Exception("no data")
        
class Exporter(object):
  @staticmethod
  def exportAsText(hdf5, header, columns, outputfile, separator=",", isdouble=False): 
      """
      export selected data into text file
      """      
      factory = PointerFactory()
      w = Writer(outputfile, header, separator, isdouble)

      for i in range(0,len(columns)):
        try : 
          column = columns[i]
          w.addColumn(hdf5.index[column], hdf5.extractQuantity(hdf5.index[column]))
        except Exception as e :
          pointer = factory.getErrorPointer(column)
          w.addColumn(pointer, [])
      
      
      try : 
        w.writeText()
      except Exception as e :
        raise e
      
  @staticmethod
  def exportAsVotable(hdf5, header, columns, outputfile, isdouble=False): 
      """
      export selected data into text file
      """        
      w = Writer(outputfile, header, isdouble=isdouble)
      rejected = []
      for i in range(0,len(columns)):
        try : 
          column = columns[i]
          w.addColumn(hdf5.index[column], hdf5.extractQuantity(hdf5.index[column]))          
        except Exception as e :
          rejected.append(column)
          
      if len(rejected) > 0 :
        Exporter.exportErrors(outputfile, rejected)
        
      try : 
        w.writeVoTable()
      except Exception as e :
        raise e
        
  @staticmethod
  def exportErrors(outputfile, errors):
    parts = outputfile.split('.')
    filepath = parts[0]+'_errors.'+parts[1]
    f = open(filepath, 'w')
    f.write("unknown columns : \n")
    for error in errors : 
      f.write(error+'\n')
    f.close()
    


