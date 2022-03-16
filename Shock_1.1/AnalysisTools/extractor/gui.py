#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = "Nicolas Moreau"
__copyright__ = "Copyright 2012"
__credits__ = ["Nicolas Moreau"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Nicolas Moreau"
__email__ = "nicolas.moreau@obspm.fr"
__status__ = "Development"  

import sys, os, time, os.path
from PyQt4 import QtGui
from PyQt4 import QtCore

import hdf
import config as cfg
import util as u
import sampy as s


class Samp(object):
  """
  class managing samp connection and sending messages
  """
  def __init__(self, applicationName, applicationDescription = ''):
  
    self.sentTableCount = 0
    """ nuimber of table sent """
    self.client = s.SAMPIntegratedClient()
    """ samp client """
    self.applicationName = applicationName
    self.applicationDescription = applicationDescription

  
  def sendTable(self, filepath) :
    """
    send a votable
    filepath : url of the table
    """
    if self.client.isConnected() is False :  
      try : 
        self.__connectHub()
      except s.SAMPHubError as e: 
        raise e        
        
    if self.client.isConnected():  
      self.sentTableCount += 1
      self.client.notifyAll({"samp.mtype": "table.load.votable",
                             "samp.params": {"url": "file:"+filepath, "name":self.applicationName+" table "+str(self.sentTableCount)}})    


  def isConnected(self) :
    """
      returns connection status
    """
    return self.client.isConnected()
    
  def closeConnection(self):
    """
      closes connection to hub
    """
    self.client.disconnect()
    self.sentTableCount = 0
      

  def __hubNotification(self, private_key, sender_id, mtype, params, extra):
    """ 
    dispatch messages received from samp hub
    """
    if mtype == 'samp.hub.event.shutdown':
      self.closeConnection()     

  def __connectHub(self):   
    """
    establish a connection with samp hub
    """
    try:            
      self.client.connect()
      self.client.declareMetadata(metadata = {"samp.name":self.applicationName,
                 "samp.description.text":self.applicationDescription,
                 "cli1.version":"1"})
      self.client.bindReceiveNotification("samp.hub.event.shutdown", self.__hubNotification)
      self.client.bindReceiveNotification("samp.hub.disconnect", self.__hubNotification)           
      self.client.declareSubscriptions()
    except s.SAMPHubError as e: 
      raise e


class ExtractorWindow(QtGui.QMainWindow):
      
  def __init__(self):
    super(ExtractorWindow, self).__init__()
    
    #self.defaultColumn = None   
    """ column automcatically exported with data for plotting """
    
    #~ self.defaultColumnAdded = False
    """ True when the default column has already been appended to the current selection"""
    
     
    self.hdfData = None
    """ PdrHDF object containing PDR HDF5 data """
    
    self.datasets = []
    """ list of dataset names (string) """
    
    self.treeWidget = None
    """ Widget displaying tree to navigate quantities """
    
    self.quantities = []
    """ List of all quantities  """
    
    self.rootItem = None
    """ Root element in the data tree  """
    
    self.separator = cfg.separator
    """ File menu in the menubar """
    
    self.exportedData = []
    """ Name of quantities to extract  """    
    
    self.samp = None
    """ Samp connector  """
       
    self.EXPORT_AS_TEXT = "txt"    
    self.EXPORT_AS_VOTABLE = "votable"
    
    
    """  buttons  """
    self.exportTextButton = None
    self.exportVOTableButton = None
    self.saveScriptButton = None
    self.sendTableButton = None
    self.searchBarButton = None
    self.saveScriptAction = None
    self.applyScriptAction = None
    
    self.quantityList = None
    """ QListWidget containing quantities to extract  """    
    
    self.searchBar = None  
    """ search bar widget  """
    
    self.openedFilePath = None
        
    self.initUI()
    
  def _getSamp(self):
    if self.samp is None : 
      self.samp =  Samp(cfg.applicationName, cfg.applicationDescription)
      
    return self.samp
    
  def closeEvent(self, event):
    """
      redefine closeEvent method, close samp connection if it exists
    """
    if self.samp is not None and self._getSamp().isConnected():
      self._getSamp().closeConnection()
      
    return super(ExtractorWindow, self).closeEvent(event)
    
      
  def initUI(self):           
    """
      init GUI
    """    

    self.openedFilePath = self.__getDefaultDataDirectory()
    self.setCentralWidget(QtGui.QWidget())
    mainLayout = QtGui.QGridLayout()       
    self.initMenuBar()
    
    searchGroup = QtGui.QGroupBox()
    searchGroup.setFlat(True)
    searchGroup.setStyleSheet('QGroupBox{border : none; padding-top : 10px;}')
    searchGroup.setFlat(True)
    searchLayout = QtGui.QVBoxLayout()
    
    searchBarLayout = QtGui.QHBoxLayout()
    self.searchBar = QtGui.QLineEdit()
    self.searchBar.setDisabled(True)
    self.searchBar.setAttribute(QtCore.Qt.WA_MacShowFocusRect, False)
    self.searchBar.returnPressed.connect(self.__addSearchBarQuantity)
    self.searchBarButton = QtGui.QPushButton('Confirm')
    self.searchBarButton.clicked.connect(self.__addSearchBarQuantity)
    self.searchBarButton.setEnabled(False)
    searchBarLayout.addWidget(self.searchBar)
    searchBarLayout.addWidget(self.searchBarButton)    
    searchLayout.addLayout(searchBarLayout) 
       
    self.treeWidget = self.getTree()
    self.treeWidget.setAttribute(QtCore.Qt.WA_MacShowFocusRect, False)
    searchLayout.addWidget(self.treeWidget)
    searchLayout.setStretch(1, 1)
    
    selectionGroup = QtGui.QGroupBox()
    selectionGroup.setFlat(True)
    selectionGroup.setStyleSheet('QGroupBox{border : none; padding-top : 10px;}')
    selectionLayout = QtGui.QVBoxLayout()
    
    self.quantityList = QtGui.QListWidget()
    self.quantityList.itemClicked.connect(self.__removeQuantity)
    self.quantityList.setAttribute(QtCore.Qt.WA_MacShowFocusRect, False)
    

    quantityButtonLayout = QtGui.QHBoxLayout()
    
    self.removeAllButton = QtGui.QPushButton("Remove All")
    self.removeAllButton.clicked.connect(self.__removeAllQuantities)
    self.removeAllButton.setEnabled(False)
    quantityButtonLayout.addWidget(self.removeAllButton)
    selectionLayout.addLayout(quantityButtonLayout)
    selectionLayout.addWidget(self.quantityList)
        
    
    searchGroup.setLayout(searchLayout)    
    selectionGroup.setLayout(selectionLayout)
    
    exportAreaLayout = QtGui.QHBoxLayout()
    exportAreaLayout.addStretch(0)

    self.exportTextButton = QtGui.QPushButton("Export as Text")
    self.exportTextButton.setEnabled(False)
    self.exportTextButton.clicked.connect(self.__exportAsText)
    exportAreaLayout.addWidget(self.exportTextButton, QtCore.Qt.AlignJustify)
    
    self.exportVOTableButton = QtGui.QPushButton("Export as VOTable")
    self.exportVOTableButton.setEnabled(False)
    self.exportVOTableButton.clicked.connect(self.__exportAsVoTable)
    exportAreaLayout.addWidget(self.exportVOTableButton, QtCore.Qt.AlignJustify)
   
    self.sendTableButton = QtGui.QPushButton("Send Table")
    self.sendTableButton.setEnabled(False)
    self.sendTableButton.clicked.connect(self.__sendSampTable)
    exportAreaLayout.addWidget(self.sendTableButton, QtCore.Qt.AlignJustify)
    exportAreaLayout.addStretch(9)
    
    mainLayout.addWidget(searchGroup, 0, 0, 1, 1)
    mainLayout.addWidget(selectionGroup, 0, 2, 1, 2)
    mainLayout.addLayout(exportAreaLayout, 1, 0, 1, 3 )    
  
    self.centralWidget().setLayout(mainLayout)
      
    self.setGeometry( 300, 300, 850, 600)
    self.setWindowTitle('PDR Extractor')    
    self.show()
    
    
  def initMenuBar(self):

    fileMenu = self.menuBar().addMenu("&File")  # Alt+F as shortcut
    scriptMenu = self.menuBar().addMenu("&Script")
    outputMenu = self.menuBar().addMenu("Output")

    openAction = QtGui.QAction("&Open HDF5", self)  
    openAction.triggered.connect(self.__parseHdf)
    fileMenu.addAction(openAction)
    
    self.saveScriptAction = QtGui.QAction("Save Script", self)  
    self.saveScriptAction.setEnabled(False)
    self.saveScriptAction.triggered.connect(self.__saveScript)
    scriptMenu.addAction(self.saveScriptAction)
        
    self.applyScriptAction = QtGui.QAction("Apply Script", self)  
    self.applyScriptAction.setEnabled(False)
    self.applyScriptAction.triggered.connect(self.__askOpenScriptFile)
    scriptMenu.addAction(self.applyScriptAction)
    
    separatorAction = QtGui.QAction('Separator', self)
    separatorAction.triggered.connect(self.__askChangeSeparator)
    outputMenu.addAction(separatorAction)  
       
    
  def getTree(self):
       
    treeWidget = QtGui.QTreeWidget()
    treeWidget.setSelectionMode(QtGui.QAbstractItemView.MultiSelection)
    treeWidget.itemSelectionChanged.connect(self.__addTreeQuantity)
    treeWidget.setColumnCount(1)
    items = []
    
    self.rootItem = treeWidget.invisibleRootItem()
    
    treeWidget.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding) 
    treeWidget.header().setStretchLastSection(False)
    treeWidget.header().setResizeMode(0, QtGui.QHeaderView.ResizeToContents)
    treeWidget.setHeaderLabels([''])

    return treeWidget  
    
    
  def __enableButtons(self):
    self.exportTextButton.setEnabled(True)
    self.exportVOTableButton.setEnabled(True)
    self.sendTableButton.setEnabled(True)
    self.searchBarButton.setEnabled(True)
    self.removeAllButton.setEnabled(True)
    self.saveScriptAction.setEnabled(True)
    self.applyScriptAction.setEnabled(True)
    self.searchBar.setEnabled(True)
    
  def __parseHdf(self):
    """
    read hdf file
    """

    filename = QtGui.QFileDialog.getOpenFileName(self, "Choose a PDR file", self.__getDataDirectory(), "HDF5 Files (*.hdf5)").toUtf8()
    filename = str(filename)
      
    if filename != "" :       
      self.openedFilePath = os.path.dirname(filename)
      self.hdfData = hdf.PdrHDF(filename)
      self.datasets = sorted(self.hdfData.index.keys())    
      self.loadNewFile(os.path.basename(filename))
      
  def __initArrays(self):
    self.datasets = []
    self.quantities = []   
    self.exportedData = []

    
  def loadNewFile(self, fileName):   
    """  
        update the GUI after a new file has been loaded
    """        

    nodes = {}    
    path = []
    if self.treeWidget is not None :
      self.treeWidget.clear()
      
    self.__initArrays()      

    for element in self.hdfData.quantities : 
      currentElement = self.rootItem
      path = element.split(self.hdfData.separator)
      path = filter(None, path) #removes empty elements after split
      self.quantities.append(path[len(path)-1])
      indexPath = ''
      for part in path :
        indexPath = indexPath + '/' + part
        
        if indexPath in nodes :
          child = nodes[indexPath]
        else :
          child = None
       
        if child is None :  
          item = QtGui.QTreeWidgetItem(currentElement)
          item.setText(0, part)
          item.setData(1, QtCore.Qt.UserRole, QtCore.QVariant(part))
          
          # add tooltip on tree elements
          if part in self.hdfData.index :
            pointer = self.hdfData.index[part]
            try :
              if pointer.description is not None : 
                # clean description content
                description = pointer.description.replace("<![CDATA[", "").replace("]]>", "")
                item.setToolTip(0, description)

            except AttributeError :
              pass
            
          currentElement = item
          nodes[indexPath] = item
        else : 
          currentElement = child   
          
    self.treeWidget.setHeaderLabels([fileName])
    self.__initCompleter()
    self.__enableButtons()
    
  def __initCompleter(self):
    """    
      Init search bar widget with quantity list
    """
    completer = QtGui.QCompleter(self.quantities, self)
    completer.setCaseSensitivity(QtCore.Qt.CaseSensitive)
    self.searchBar.setCompleter(completer)
   
  def __removeQuantity(self, item):
    """
      remove selected quantity from quantity list
    """
    self.quantityList.takeItem(self.quantityList.row(item))

  def __removeAllQuantities(self):
    """
      remove all selected quantities
    """
    self.quantityList.clear()
    self.treeWidget.clearSelection()
    
  def __addSearchBarQuantity(self):
    """
      add a quantity from the search bar
    """
    value = self.searchBar.text()
    if value in self.quantities : 
      if  len(self.quantityList.findItems(value, QtCore.Qt.MatchExactly)) == 0 :
        self.quantityList.addItem(value)
    else :
      QtGui.QMessageBox.warning(self, "Warning", "Invalid quantity")
  
  def __addTreeQuantity(self):
    """
      add the selected quantity in the tree into the selected quantity list
    """
    for item in self.treeWidget.selectedItems():
      value = item.data(1, QtCore.Qt.UserRole).toString()
      if len(self.quantityList.findItems(value, QtCore.Qt.MatchExactly)) == 0 and item.childCount() == 0 :        
        self.quantityList.addItem(value)
            
  def __setExportedData(self):
    """    
      build a list of quantities to extract
    """
    self.exportedData = []
    
    for i in range(0, self.quantityList.count()):
      if str(self.quantityList.item(i).text()) != '' :
        self.exportedData.append(str(self.quantityList.item(i).text()))    
      
  def __exportAsText(self): 
    """
    export selected data into text file
    """
    self.__setExportedData()
    self.__export(self.EXPORT_AS_TEXT)

  def __exportAsVoTable(self): 
    """
    export selected data into votable file
    """
    self.__setExportedData()
    self.__export(self.EXPORT_AS_VOTABLE)
    
    
  def __getDataDirectory(self):   
    """
    returns a directory where to look for and save data file
    """
    if self.openedFilePath is not None : 
      return  self.openedFilePath    
    return os.getcwd()   

  def __getDefaultDataDirectory(self):   
    """
    returns a directory where to look for and save data file
    """
    if os.path.isdir(os.path.join("..", cfg.pdrDataDirectory)):
      return os.path.join("..", cfg.pdrDataDirectory)
    return os.getcwd()   
  
      
  def __export(self, mode): 
    """
    export selected data into text file
    """        
    isDouble = False
    
    path = self.__getDataDirectory()
    
    try :             
      if mode == self.EXPORT_AS_TEXT :  
        filename = str(QtGui.QFileDialog.getSaveFileName(self, "Export data to", path, ""))
        hdf.Exporter.exportAsText(self.hdfData, None, self.exportedData, filename, self.separator, isDouble)
      elif mode == self.EXPORT_AS_VOTABLE :  
        filename = str(QtGui.QFileDialog.getSaveFileName(self, "Export data to", path, "XML Files (*.xml)"))
        hdf.Exporter.exportAsVotable(self.hdfData, None, self.exportedData, filename, isDouble)
      else : 
        QtGui.QMessageBox.warning(self, "Warning", "Unknown export mode")
        
      QtGui.QMessageBox.information(self, "Information", "Data have been exported")

    except Exception as e :
      QtGui.QMessageBox.warning(self, "Warning", str(e))
      
  def __saveScript(self):
    """
    export selected coumns as script file
    """
    userfilename = QtGui.QFileDialog.getSaveFileName(self, "Export script to", self.__getDataDirectory(), cfg.scriptFileDescription + " ( "+ "*."+cfg.scriptFileExtension +" )")
    filename = self.__setExtension( str(userfilename) , cfg.scriptFileExtension )
    f = open(filename,'w')
    
    for i in range(0, self.quantityList.count()):
      if str(self.quantityList.item(i).text()) != '' :
        f.write(str(self.quantityList.item(i).text())+"\n")       

    f.close()   
    QtGui.QMessageBox.information(self, "Information", "Script has been saved") 
      
  def __setExtension(self, filename, extension):
    if filename.endswith(extension) :
      return filename
    else :
      parts = filename.split(".")
      if len(parts) > 1 :
        return parts[0] + "." + extension
      else :  
        return filename + "." + extension
        
  def __askOpenScriptFile(self):
    """Returns an opened file in read mode."""
    filename = str(QtGui.QFileDialog.getOpenFileName(self, "Select a template", os.getcwd(), cfg.scriptFileDescription + " ( "+ "*."+cfg.scriptFileExtension +" )"))

    if os.path.isfile(filename):
        self.__execScript(filename)
    else:
      QtGui.QMessageBox.warning(self, "Warning", filename+" is not a file")
      
  def __askChangeSeparator(self):
    result = QtGui.QInputDialog.getText(self, "Set column separator", "New separator : ", text=self.separator)
    if len(result) > 0:
      self.separator = str(result[0])
    else :
      QtGui.QMessageBox.warning(self, "Warning", "Separator can not be empty")
      
  def __sendSampTable(self):  
    #~ self.__appendMeshColumns()
    #~ isDouble = True
    #~ if self.precision == cfg.floatPrecision:
    self.__setExportedData()
    isDouble = False
    filename = os.path.expanduser("~/.extractor.xml")

    try :
      hdf.Exporter.exportAsVotable(self.hdfData, None, self.exportedData, filename, isDouble)
      self._getSamp().sendTable(filename)
    except Exception as e:
      QtGui.QMessageBox.warning(self, "Warning", str(e))


  def __execScript(self, filename):      
    """
    run an export script
    """
    self.exportedData = u.ScriptReader.getColumns(filename)
    #self.precision, self.exportedDataHeader = u.ScriptReader.getHeader(filename)
    self.__export(self.EXPORT_AS_TEXT)
    del self.exportedData[:]



def launch():    
  app = QtGui.QApplication(sys.argv)
  app.setApplicationName("PDR Extractor")
  ex = ExtractorWindow()
  sys.exit(app.exec_())
