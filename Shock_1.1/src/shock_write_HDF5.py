# Copyright (C) 2019 - Observatoire de Paris
#                      CNRS (France)
#                      Contacts : antoine.gusdorf@lra.ens.fr
#                                 sylvie.cabrit@obspm.fr
#                                 benjamin.godard@obspm.fr
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#                                     PXDR_WRITE.py
#
#  Description : shock_write_hdf5.py is a script that exports some quantities from the ASCII
#                files produced by the hock code and store them in a HDF5 file. It adds 
#                metadata to quantities as units, ucd, utypes, ... that help users to 
#                understand the quantities computed by the shock code and facilitate the
#                integration of shock results in VO tools. The script is launched from 
#                the shock code itself ... see mhd_vode.f90

# ========================================================================================
# IMPORT ALL THE LIBRARIES NEEDED
# ========================================================================================
import numpy as np
import os
import shutil
import h5py
import sys
import datetime

# ========================================================================================
# DEFINITION OF A FUNCTION THAT COMPUTES MEAN VALUES
# ========================================================================================
def comp_mean(x,y):
   mean = 0.0
   npt  = len(x)
   for i in range(0,len(x)-1) :
       mean = mean + 0.5 * (y[i+1] + y[i]) * (x[i+1] - x[i])
   mean = mean / (x[npt-1]-x[0])
   return mean

# ========================================================================================
# START THE SCRIPT - MESSAGE ON SCREEN
# ========================================================================================
print " "
print "==================================================="
print "            START RUNNING PYTHON SCRIPT            "
print "==================================================="

# ========================================================================================
# STORE THE ARGUMENTS OF THE PYTHON SCRIPT
#    - 1st argument ...... path of the directory where ASCII files are stored
#    - 2rd argument ...... name of the output directory
#    - 3rd argument ...... name of the model
#    - 4th argument ...... ifaf iteration number (used in PDR model so far)
#    - 5th argument ...... code version (used in PDR model so far)
#    - 6th argument ...... flag to write the stat file
# Example to run the script manually:
# python2.7 shock_write_HDF5.py output/ASCII_files_dat/ output/ modelname
# ========================================================================================
ascii_dir    = sys.argv[1]     # Name of the directory where ASCII files are stored
out_dir      = sys.argv[2]     # Name of the output directory
out_file     = sys.argv[3]     # Name of the PDR model
#ifaf_nb      = sys.argv[4]     # ifaf iteration number
#code_version = sys.argv[5]     # Version of the code 
write_stat   = sys.argv[4]     # Flag: determine if .stat file must be written

pos_file     = ascii_dir+"Positions.dat"        # ASCII file containing positions
par_file     = ascii_dir+"InputParameters.dat"  # ASCII file containing input parameters
inf_file     = ascii_dir+"Informations.dat"     # ASCII file containing informations on models
net_file     = ascii_dir+"Network.dat"          # ASCII file containing chemical network

# ========================================================================================
# CREATE THE ID AND THE NAMES OF THE HDF5 AND STAT FILES TO PRODUCE
#    - the ID is based on the date to get unique file names
#    - the hdf5 name cannot contain any '.' -> replace by '_'
# ========================================================================================
modelID  = out_file # + "_" + ifaf_nb
modelID  = modelID.replace(".","_")
date     = datetime.datetime.now()

hdf5name = out_dir + modelID + ".hdf5"
hdf5file = h5py.File(hdf5name,'w')  # Open the HDF5 file in write mode

if (write_stat == "1"):
   statname = out_dir + modelID + ".stat"
   statfile = open(statname, 'w')      # Open the stat file in write mode
   statfile.write(str(date)    + "\n")
   #statfile.write(code_version + "\n")
   statfile.write(modelID      + "\n")
   statfile.write("\n")

# ========================================================================================
# GET THE POSITION TABLE FROM THE CORRESPONDING ASCII FILE
#    - since lines begin with "val" in the ascii file, use "usecols=(1,)"
#      in the function np.loadtxt to prevent reading the first empty column
# ========================================================================================
pos_data = np.loadtxt(pos_file,comments="@",delimiter="val=",usecols=(0,))
npt_Av   = len(pos_data)

# ========================================================================================
# SET THE LIST OF METADATA
#    - all columns in all datasets in the HDF5 file are described in the Metadata dataset
# ========================================================================================
#PathHDF5 | DataSet | Col.Index | PubID | HName  | DataType | Unit | SKOS | UCD | utype | Description | ObjParent | Group |
#   0     |    1    |     2     |   3   |   4    |    5     |  6   |  7   |  8  |  9    |     10      |    11     |  12   |
MetaList = [[],[],[],[],[],[],[],[],[],[],[],[],[]]

# ========================================================================================
# LIST THE NAMES OF ASCII FILES IN THE OUTPUT DIRECTORY
# ========================================================================================
import glob
ListAscii = glob.glob(ascii_dir +"*.dat")
ndataset    = len(ListAscii)

# ========================================================================================
# LOOP ON ASCII FILE: READ METADATA AND DATA AND WRITE EVERYTHING IN HDF5 AND STATS FILES
# ========================================================================================
for i in range(ndataset) :
   print "Transform ASCII file : "+ListAscii[i].strip()
   ASCIIFile = open(ListAscii[i].strip(), 'r')

   # -------------------------------------------------------------------------------------
   # Initialization of booleans to check if all metadata have been read
   # -------------------------------------------------------------------------------------
   path_found   = False
   dset_found   = False
   PubID_found  = False
   Hname_found  = False
   dtype_found  = False
   unit_found   = False
   SKOS_found   = False
   UCD_found    = False
   utype_found  = False
   Desc_found   = False
   Group_found  = False
   ObjPar_found = False

   for line in ASCIIFile :
      # ----------------------------------------------------------------------------------
      # Read the metadata associated to each column in each ASCII files
      #    - split the lines with the separator @
      #    - be careful, the lines in ASCIIfile starts with @, so there is one more
      #      element than it should which corresponds to the first "empty" column
      # ----------------------------------------------------------------------------------
      if line.startswith("@path=")  == True :
         path_found     = True
         line_path      = line
         line_path_el   = line_path.split("@")
      if line.startswith("@dset=")  == True :
         dset_found     = True
         line_dset      = line
         line_dset_el   = line_dset.split("@")
      if line.startswith("@PubID=") == True :
         PubID_found    = True
         line_PubID     = line
         line_PubID_el  = line_PubID.split("@")
      if line.startswith("@Name=") == True :
         Hname_found    = True
         line_name      = line
         line_name_el   = line_name.split("@")
      if line.startswith("@dtype") == True :
         dtype_found    = True
         line_dtype     = line
         line_dtype_el  = line_dtype.split("@")
      if line.startswith("@unit=") == True :
         unit_found     = True
         line_unit      = line
         line_unit_el   = line_unit.split("@")
      if line.startswith("@SKOS=") == True :
         SKOS_found     = True
         line_SKOS      = line
         line_SKOS_el   = line_SKOS.split("@")
      if line.startswith("@UCD=") == True :
         UCD_found      = True
         line_UCD       = line
         line_UCD_el    = line_UCD.split("@")
      if line.startswith("@utype=") == True :
         utype_found    = True
         line_utype     = line
         line_utype_el  = line_utype.split("@")
      if line.startswith("@Desc=") == True :
         Desc_found     = True
         line_Desc      = line
         line_Desc_el   = line_Desc.split("@")
      if line.startswith("@Group=") == True :
         Group_found    = True
         line_Group     = line
         line_Group_el  = line_Group.split("@")
      if line.startswith("@ObjPar=") == True :
         ObjPar_found   = True
         line_ObjPar    = line
         line_ObjPar_el = line_ObjPar.split("@")

   # -------------------------------------------------------------------------------------
   # Nb of columns in the ASCII File. Remove 1 because lines start with @ (see above)
   # -------------------------------------------------------------------------------------
   ncol = len(line_unit_el) - 1

   # -------------------------------------------------------------------------------------
   # Fill the list of metadata - check if they have been read or give default values
   # note: the PubID must only contain lower case and no capital letters
   # note: save the index at which we continue the metadata list
   # -------------------------------------------------------------------------------------
   index_list = len(MetaList[0])
   for j in range(1, ncol+1) :
      if (path_found == True) :
         MetaList[0].append(line_path_el[j].replace("path=","").replace("\n",""))
      else :
         MetaList[0].append("/Default")
      if (dset_found == True) :
         MetaList[1].append(line_dset_el[j].replace("dset=","").replace("\n",""))
      else :
         MetaList[1].append("Missing Dataset")
      MetaList[2].append(str(j-1))
      if (PubID_found == True) :
         MetaList[3].append(line_PubID_el[j].replace("PubID=","").replace("\n","").lower())
      else :
         MetaList[3].append("Missing PubID")
      if (Hname_found == True) :
         MetaList[4].append(line_name_el[j].replace("Name=","").replace("\n",""))
      else :
         MetaList[4].append("Missing Human name")
      if (dtype_found == True) :
         MetaList[5].append(line_dtype_el[j].replace("dtype=","").replace("\n",""))
      else :
         MetaList[5].append("To Fill with a datatype")
      if (unit_found == True) :
         MetaList[6].append(line_unit_el[j].replace("unit=","").replace("\n",""))
      else :
         MetaList[6].append("To Fill with a unit")
      if (SKOS_found == True) :
         MetaList[7].append(line_SKOS_el[j].replace("SKOS=","").replace("\n",""))
      else :
         MetaList[7].append("To Fill with SKOS")
      if (UCD_found == True) :
         MetaList[8].append(line_UCD_el[j].replace("UCD=","").replace("\n",""))
      else :
         MetaList[8].append("To Fill with UCD")
      if (utype_found == True) :
         MetaList[9].append(line_utype_el[j].replace("utype=","").replace("\n",""))
      else :
         MetaList[9].append("To Fill with utype")
      if (Desc_found == True) :
         MetaList[10].append(line_Desc_el[j].replace("Desc=","").replace("\n",""))
      else :
         MetaList[10].append("No description")
      if (Group_found == True) :
         MetaList[11].append(line_Group_el[j].replace("Group=","").replace("\n",""))
      else :
         MetaList[11].append("No Group")
      if (ObjPar_found == True) :
         MetaList[12].append(line_ObjPar_el[j].replace("ObjPar=","").replace("\n",""))
      else :
         MetaList[12].append("No Parent Object")

   # -------------------------------------------------------------------------------------
   # Check if all the columns in a given ascii file have the same path and dataset. This
   # is mandatory for the production of a coherent tree corresponding to the ascii files
   # -------------------------------------------------------------------------------------
   for j in range(1, ncol+1) :
      if (line_path_el[j].replace("path=","").strip() != line_path_el[1].replace("path=","").strip()) :
         print "*** Error: the columns in the output ascii files have different paths"
         print "***    - file   = ", ListAscii[i].strip()
         print "***    - path 1 = ", line_path_el[1]
         print "***    - path j = ", line_path_el[j]
         exit()
      if (line_dset_el[j].replace("dset=","").strip() != line_dset_el[1].replace("dset=","").strip()) :
         print "*** Error: the columns in the output ascii files have different datasets"
         print "***    - file   = ", ListAscii[i].strip()
         print "***    - dset 1 = ", line_dset_el[1]
         print "***    - dset j = ", line_dset_el[j]
         exit()

   # -------------------------------------------------------------------------------------
   # Create groups and tree in the HDF5 file
   # -------------------------------------------------------------------------------------
   path   = line_path_el[1].replace("path=","").strip()
   datset = line_dset_el[1].replace("dset=","").strip()
   group  = hdf5file.require_group(path)

   # -------------------------------------------------------------------------------------
   # Read the data in the ascii file - differentiate the parameter file, info file or
   # network file (which contains strings) from all the other files which contain reals.
   # -------------------------------------------------------------------------------------
   if (ListAscii[i].strip() == par_file or ListAscii[i].strip() == inf_file  or ListAscii[i].strip() == net_file) :
      Data = np.loadtxt(ListAscii[i].strip(), dtype="string", comments="@", delimiter="val=")
   else :
      Data = np.loadtxt(ListAscii[i].strip(), comments="@", delimiter="val=", ndmin=2)

   # -------------------------------------------------------------------------------------
   # Reshape the data if necessary - this step is necessary when only one line is present
   # in data. If not done, it will appear as one column in HDF5. Since all data have to
   # be stored in NUMPY arrays, scalar have to be stored in a (1,1) matrix
   # -------------------------------------------------------------------------------------
   if Data.ndim == 0 :
      Data = Data.reshape(1,1)
   if Data.ndim == 1 and Data.shape[0] != npt_Av :
      Data = Data.reshape(1,Data.size)

   # -------------------------------------------------------------------------------------
   # Store Data in the HDF5 dataset and add atributes to the dataset
   # -------------------------------------------------------------------------------------
   dset = group.create_dataset(datset,data=Data)
   dset.attrs["CLASS"] = "TABLE"

   # -------------------------------------------------------------------------------------
   # Compute the mean values of the data and write the results in the stat file
   # -------------------------------------------------------------------------------------
   if (write_stat == "1"):
      for j in range(1, ncol+1) :
         # Quantity function of AV (CAREFUL NOT ROBUST TEST)
         if Data.shape[0] == npt_Av :
            column = Data.take([j-1],axis=1)
            column_mean = comp_mean(pos_data,column[:,0])
            statfile.write("min  # %s # %E # %s # %s # %s\n"   % (MetaList[3][index_list+j-1], column.min(), MetaList[6][index_list+j-1], MetaList[12][index_list+j-1], modelID+"_"+MetaList[12][index_list+j-1]))
            statfile.write("mean # %s # %E # %s # %s # %s\n"   % (MetaList[3][index_list+j-1], column_mean,  MetaList[6][index_list+j-1], MetaList[12][index_list+j-1], modelID+"_"+MetaList[12][index_list+j-1]))
            statfile.write("max  # %s # %E # %s # %s # %s\n"   % (MetaList[3][index_list+j-1], column.max(), MetaList[6][index_list+j-1], MetaList[12][index_list+j-1], modelID+"_"+MetaList[12][index_list+j-1]))

         # Quantity on a line that is not Input Parameter or Information
         elif Data.shape[0] == 1 and ListAscii[i].strip() != par_file and ListAscii[i] != inf_file and ListAscii[i] != net_file :
            column_mean = Data.take([j-1],axis=1)[[0]]
            statfile.write("value # %s # %E # %s # %s # %s\n"  % (MetaList[3][index_list+j-1], column_mean,  MetaList[6][index_list+j-1], MetaList[12][index_list+j-1], modelID+"_"+MetaList[12][index_list+j-1]))

         elif Data.shape[0] != npt_Av and Data.shape[0] != 1 :
            print "Case not understood"
            print "Data.shape[0] = " + str(Data.shape[0]) + " " + MetaList[3][index_list+j-1]

      ASCIIFile.close()

# ========================================================================================
# DEFINE THE PATH TOWARDS THE METADATA DATASET
# ========================================================================================
MetaDataPath = hdf5file.create_group("/Metadata")

# ========================================================================================
# WRITE METADATA TABLES IN THE HDF5 FILE
#    - convert python array_like object to numpy array
#    - write the transposed table in the hdf5 file
#
# HDF5 supports types which have no direct NumPy equivalent. We use here a variable-length
# (VL) type to store arbitrary-length vectors of a given base type. This is essential to
# limit the size of the HDF5 file, which otherwise, stores a table of string whose size
# seem to be that of the largest string in Metalist.
#    - definition of dt type
#    - use dt to store the Metadata table.
#
# Use additional HDF5 compression options. Level 9 of gzip compression reduces the size of
# the output hdf5 files by about 20%
# ========================================================================================
dt = h5py.special_dtype(vlen=str)
MetaTable = np.array(MetaList)
MetaTable = MetaTable.T
dset = MetaDataPath.create_dataset("Metadata",data=MetaTable,dtype=dt,compression="gzip", compression_opts=9)
dset.attrs["CLASS"] = "TABLE"

# ========================================================================================
# SET THE INFORMATION STORED IN METATABLE_OBJ_TYPE
#    - There are three types of objects: grains, meshcell, and cloud
# ========================================================================================
# PubID | Hname | SKOS | Description | ObjParent | Number_of_Object |
#   0   |   1   |   2  |      3      |     4     |         5        |
MetaTable_ObjType = [[],[],[],[],[],[]]

MetaTable_ObjType[0].append("grain")
MetaTable_ObjType[1].append("Grain")
MetaTable_ObjType[2].append("NO SKOS")
MetaTable_ObjType[3].append("Interstellar grain")
MetaTable_ObjType[4].append("meshcell")
MetaTable_ObjType[5].append("17")

MetaTable_ObjType[0].append("meshcell")
MetaTable_ObjType[1].append("Mesh Cell")
MetaTable_ObjType[2].append("http://purl.org/astronomy/vocab/DataObjectTypes/AdaptiveMeshCell")
MetaTable_ObjType[3].append("Adaptive Mesh cell")
MetaTable_ObjType[4].append("cloud")
MetaTable_ObjType[5].append(str(npt_Av))

MetaTable_ObjType[0].append("cloud")
MetaTable_ObjType[1].append("Cloud")
MetaTable_ObjType[2].append("http://purl.org/astronomy/vocab/AstronomicalObjects/MolecularCloud")
MetaTable_ObjType[3].append("Interstellar cloud")
MetaTable_ObjType[4].append("")
MetaTable_ObjType[5].append("1")

# ========================================================================================
# WRITE METADATA OBJECT TYPES TABLES IN THE HDF5 FILE
#    - convert python array_like object to numpy array
#    - write the transposed table in the hdf5 file
# ========================================================================================
MetaTable_ObjType_arr = np.array(MetaTable_ObjType)
MetaTable_ObjType_arr = MetaTable_ObjType_arr.T
dset = MetaDataPath.create_dataset("MetaData_ObjectType",data=MetaTable_ObjType_arr)
dset.attrs["CLASS"] = "TABLE"

# ========================================================================================
# WRITE METADATA OBJECT TYPES TABLES IN THE STAT FILE
# ========================================================================================
if (write_stat == "1"):
   for j in range (0,3) :
      statfile.write( "value # number_of_objects # %s # dimensionless # %s # %s \n" % \
      (MetaTable_ObjType[5][j], MetaTable_ObjType[0][j], modelID+"_"+MetaTable_ObjType[0][j]) )

# ========================================================================================
# CLOSE THE HDF5 AND STAT FILES AND END THE PYTHON SCRIPT
# ========================================================================================
hdf5file.close() # close the HDF5 file
if (write_stat == "1"):
   statfile.close() # close the stat file

print('removing the {} Ascii files folder ...'.format(ascii_dir))
shutil.rmtree(ascii_dir)
print('done')

print " "
print "==================================================="
print "               END OF PYTHON SCRIPT                "
print "==================================================="
