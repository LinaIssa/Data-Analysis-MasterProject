# -*- coding: utf-8 -*-
#!/usr/bin/python
import sys
import getopt
import hdf
import config as cfg
import util as u

__author__ = "Nicolas Moreau"
__copyright__ = "Copyright 2012"
__credits__ = ["Nicolas Moreau"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Nicolas Moreau"
__email__ = "nicolas.moreau@obspm.fr"
__status__ = "Development"   
             
        
def display_help():
    """
        display help of the program
    """
    print '-h, --help : display this help\n\
-f, --file : HDF5 file to read\n\
-t, --template : script to apply (.esf file)\n\
-o, --output : name of output file\n\
-s, --separator : separator between columns, default is ,'

def error(message):
    """
        print message and exit the program
    """
    print message
    sys.exit()


def main(argv):    
    try:                                
        opts, args = getopt.getopt(argv, "hf:t:o:s:", ["help", "file=","template=", "output=", 'separator='])
    except getopt.GetoptError:          
        display_help()                         
        sys.exit(2)

    if len(opts) == 0:
        import gui
        gui.launch()        
    else:
        hdffile = None
        scriptfile = None
        outputfile = None
        separator = cfg.separator
        
        for opt, arg in opts:       
            if opt in ("-h", "--help"):      
                display_help()                     
                sys.exit()     
                
            if opt in ("-f", "--file"):
                hdffile = arg
                
            if opt in ("-t","--template"):
                scriptfile = arg
                
            if opt in ("-o","--output"):
                outputfile = arg
                
            if opt in ("-s","--separator"):
                separator = arg
                
        if hdffile is None :         
            error("HDF5 file is missing")            
        
        elif scriptfile is None : 
            error("Template file is missing")
            
        elif outputfile is None : 
            error("Output file name is missing")
        
        else: 
            precision, header = u.ScriptReader.getHeader(scriptfile)
            if precision == cfg.doublePrecision:
                isDouble = True
            else: isDouble = False
            
            hdf.Exporter.exportAsText(hdf.PdrHDF(hdffile), header, u.ScriptReader.getColumns(scriptfile), outputfile, separator, isDouble)           


if __name__ == "__main__":
    main(sys.argv[1:])
