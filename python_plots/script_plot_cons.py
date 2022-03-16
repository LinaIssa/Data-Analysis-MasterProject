# ========================================================================================
# 0 - IMPORT LIBRARIES / DEFINE CONSTANTS
# ========================================================================================
from scipy import array
import matplotlib.pyplot as plt
import numpy as np
import argparse


# ----------------------------------------------------------------------------------------
# Add arguments
# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('-type', type=str,   help="output format:"+'\n'+"file (pdf) or interacticve plot (itv)")
parser.add_argument('-axis', type=str,   help="x axis choice:"+'\n'+"time (time) or distance (dist)")
parser.add_argument('-xmin', type=float, help="min x value")
parser.add_argument('-xmax', type=float, help="max x value")
args = parser.parse_args()

if args.type:
   if args.type == 'pdf':
      scptype='pdf'
      tsize=5
      msize=3
   elif args.type == 'itv':
      scptype='itv'
      tsize=10
      msize=7
else:
   scptype='itv'
   tsize=10
   msize=7

if args.axis:
   if args.axis == 'time':
      xaxis='timeN'
      xlabel='time (yr)'
      xmin = 1.0E-7
      xmax = 1.0E+9
   elif args.axis == 'dist':
      xaxis='distance'
      xlabel='distance (cm)'
      xmin = 1.0E+6
      xmax = 1.0E+22
else:
   xaxis='timeN'
   xlabel='time (yr)'
   xmin = 1.0E-7
   xmax = 1.0E+9

if args.xmin:
   xmin = args.xmin

if args.xmax:
   xmax = args.xmax
      
nameout='plot_cons.pdf'



# ========================================================================================
# A - FUNCTIONS DEFINITIONS
# ========================================================================================
def test(v,E,e):
   if (v[0]=='-'):
      v=-999.
   if (E[0]=='-'):
      E=-999.
   if (e[0]=='-'):
      e=-999.
   return v, E, e

def test_1(x):
   if (x[0]=='-'):
      x=-999.
   return x

def sumtab(l_x,v_x,E_x,e_x,r_x,l_y,v_y,E_y,e_y,r_y):
   l_sum = []
   v_sum = []
   E_sum = []
   e_sum = []
   r_sum = []
   for i in range(0, len(v_x)):
      if   l_x[i] == '<' or  l_y[i] == '<':
         l = '<'
         l_sum.append(l)
      elif l_x[i] == '>' or  l_y[i] == '>':
         l = '<'
         l_sum.append(l)
      else:
         l = ' '
         l_sum.append(l)
      if v_x[i] < 0 or v_y[i] < 0:
         v = -999.
      else:
         v = v_x[i] + v_y[i]
      E = np.sqrt(E_x[i]*E_x[i] + E_y[i]*E_y[i])
      e = np.sqrt(e_x[i]*e_x[i] + e_y[i]*e_y[i])
      r = r_x[i] + r_y[i]
      v_sum.append(v)
      E_sum.append(E)
      e_sum.append(e)
      r_sum.append(r)
   v_sum = array(v_sum)
   E_sum = array(E_sum)
   e_sum = array(e_sum)
   return l_sum, v_sum, E_sum, e_sum, r_sum

def dvdtab(l_x,v_x,E_x,e_x,r_x,l_y,v_y,E_y,e_y,r_y):
   l_dvd = []
   v_dvd = []
   E_dvd = []
   e_dvd = []
   r_dvd = []
   for i in range(0, len(v_x)):
      if   l_x[i] == '<' or  l_y[i] == '>':
         l = '<'
         l_dvd.append(l)
      elif l_x[i] == '>' or  l_y[i] == '<':
         l = '<'
         l_dvd.append(l)
      else:
         l = ' '
         l_dvd.append(l)
      if v_x[i] <= 0 or v_y[i] <= 0:
         v = -999.
         E = -999.
         e = -999.
      else:
         v = v_x[i] / v_y[i]
         E = np.sqrt( v * (E_x[i]*E_x[i]/v_x[i]/v_x[i] + E_y[i]*E_y[i]/v_y[i]/v_y[i]) ) 
         e = np.sqrt( v * (e_x[i]*e_x[i]/v_x[i]/v_x[i] + e_y[i]*e_y[i]/v_y[i]/v_y[i]) ) 
      r = r_x[i] + r_y[i]
      v_dvd.append(v)
      E_dvd.append(E)
      e_dvd.append(e)
      r_dvd.append(r)
   v_dvd = array(v_dvd)
   E_dvd = array(E_dvd)
   e_dvd = array(e_dvd)
   return l_dvd, v_dvd, E_dvd, e_dvd, r_dvd

def dvdbg09(v_x,v_y):
   v_dvd = []
   for i in range(0, len(v_x)):
      if v_x[i] < 0 or v_y[i] < 0:
         v = -999.
      else:
         v = v_x[i] / v_y[i]
      v_dvd.append(v)
   v_dvd = array(v_dvd)
   return v_dvd

def filter(l_x,v_x,l_y,v_y,cx,cy,xmin,xmax,ymin,ymax):
   x_filt = []
   y_filt = []
   for i in range(0, len(v_x)):
      if l_x[i] == cx and l_y[i] == cy and v_x[i] > xmin and v_y[i] > ymin and v_x[i] < xmax and v_y[i] < ymax:
         x_filt.append(v_x[i])
         y_filt.append(v_y[i])
   x_filt = array(x_filt)
   y_filt = array(y_filt)
   return x_filt, y_filt

def checkValue(value):
    # Check if value should be a float or flagged as missing
    if value == "---":
       value = np.ma.masked
    else:
       value = float(value)
    return value


# ========================================================================================
# B1 - READ INPUT CHEMICAL FILE
# ========================================================================================

# -------------------------------------------
# 1 - open file                              
# 2 - skip first line                        
# 3 - read header                            
# 4 - read data                              
# -------------------------------------------
fchem = open('mhd_speci.out','r')
dum=fchem.readline()
chem_names = fchem.readline().split()
chem_lines = fchem.readlines()

# -------------------------------------------
# 1 - Create a data dictionary, containing a 
#     list of values for each variable       
# 2 - Add entry to the dictionary for each   
#     column                                 
# 3 - Loop through each value: append to each 
#     column
# -------------------------------------------
chem_data = {}

for chem_name in chem_names:
   chem_data[chem_name] = np.ma.zeros(len(chem_lines), 'f', fill_value = -999.999)

for (line_count, line) in enumerate(chem_lines):
    items = line.split()
    for (col_count, col_name) in enumerate(chem_names):
        value = items[col_count]
        chem_data[col_name][line_count] = checkValue(value)

#print chem_names

# ========================================================================================
# B2 - READ INPUT PHYSICAL FILE
# ========================================================================================

# -------------------------------------------
# 1 - open file                              
# 2 - read header                            
# 3 - read data                              
# -------------------------------------------
fphys = open('mhd_phys.out','r')
phys_names = fphys.readline().split()
phys_lines = fphys.readlines()

# -------------------------------------------
# 1 - Create a data dictionary, containing a 
#     list of values for each variable       
# 2 - Add entry to the dictionary for each   
#     column                                 
# 3 - Loop through each value: append to each 
#     column
# -------------------------------------------
phys_data = {}

for phys_name in phys_names:
   phys_data[phys_name] = np.ma.zeros(len(phys_lines), 'f', fill_value = -999.999)

for (line_count, line) in enumerate(phys_lines):
    items = line.split()
    for (col_count, col_name) in enumerate(phys_names):
        value = items[col_count]
        phys_data[col_name][line_count] = checkValue(value)

#print phys_names

# ========================================================================================
# B3 - READ INPUT H2 FILE
# ========================================================================================

# -------------------------------------------
# 1 - open file                              
# 2 - read header                            
# 3 - read data                              
# -------------------------------------------
fhmol = open('H2_lev.out','r')
hmol_names = fhmol.readline().split()
hmol_lines = fhmol.readlines()

# -------------------------------------------
# 1 - Create a data dictionary, containing a 
#     list of values for each variable       
# 2 - Add entry to the dictionary for each   
#     column                                 
# 3 - Loop through each value: append to each 
#     column
# -------------------------------------------
hmol_data = {}

for hmol_name in hmol_names:
   hmol_data[hmol_name] = np.ma.zeros(len(hmol_lines), 'f', fill_value = -999.999)

for (line_count, line) in enumerate(hmol_lines):
    items = line.split()
    for (col_count, col_name) in enumerate(hmol_names):
        value = items[col_count]
        hmol_data[col_name][line_count] = checkValue(value)

#print hmol_names




# ========================================================================================
# C - PLOT A VARIABLE VS ANOTHER
# ========================================================================================

# -------------------------------------------
# subplot 1
# -------------------------------------------
share_ax =plt.subplot(2,2,1)

plt.plot(chem_data[xaxis], np.fabs(chem_data["sum(i)-sum(neg)"]),'-', c='blue',  alpha=1.0, markeredgecolor='none',markersize=msize,label='|n$_{\mathrm{i}}$-n$_{\mathrm{e}}$|')
plt.legend(loc='upper right', borderaxespad=0., fontsize=tsize)

plt.xlabel(xlabel, fontsize=tsize)
plt.ylabel('density (cm$^3$)', fontsize=tsize)

ax = plt.gca()
ax.set_xlim([xmin,xmax])
ax.set_ylim([1e-20,1e-15])
ax.tick_params(axis='both', which='major', labelsize=tsize)
ax.tick_params(axis='both', which='minor', labelsize=tsize)
ax.set_xscale('log')
ax.set_yscale('log')

plt.grid(True)

# -------------------------------------------
# subplot 2
# -------------------------------------------
ax2=plt.subplot(2,2,2,sharex=share_ax)

plt.plot(chem_data[xaxis], chem_data["elab_H"],  '-', c='black',   alpha=1.0, markeredgecolor='none',markersize=msize,label='[H]' )
plt.plot(chem_data[xaxis], chem_data["elab_He"], '-', c='cyan',    alpha=1.0, markeredgecolor='none',markersize=msize,label='[He]')
plt.plot(chem_data[xaxis], chem_data["elab_C"],  '-', c='red',     alpha=1.0, markeredgecolor='none',markersize=msize,label='[C]' )
plt.plot(chem_data[xaxis], chem_data["elab_N"],  '-', c='green',   alpha=1.0, markeredgecolor='none',markersize=msize,label='[N]' )
plt.plot(chem_data[xaxis], chem_data["elab_O"],  '-', c='blue',    alpha=1.0, markeredgecolor='none',markersize=msize,label='[O]' )
plt.plot(chem_data[xaxis], chem_data["elab_Si"], '-', c='magenta', alpha=1.0, markeredgecolor='none',markersize=msize,label='[Si]')
plt.plot(chem_data[xaxis], chem_data["elab_S"],  '-', c='orange',  alpha=1.0, markeredgecolor='none',markersize=msize,label='[S]' )
plt.plot(chem_data[xaxis], chem_data["elab_G"],  '-', c='grey',    alpha=1.0, markeredgecolor='none',markersize=msize,label='[G]' )
plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left', borderaxespad=0., fontsize=tsize)

plt.xlabel(xlabel, fontsize=tsize)
plt.ylabel('abundance', fontsize=tsize)

ax = plt.gca()
ax.set_xlim([xmin,xmax])
#ax.set_ylim([1e2,1e8])
ax.tick_params(axis='both', which='major', labelsize=tsize)
ax.tick_params(axis='both', which='minor', labelsize=tsize)
ax.set_xscale('log')
ax.set_yscale('log')

plt.grid(True)

# -------------------------------------------
# subplot 3
# -------------------------------------------
ax3=plt.subplot(2,2,3,sharex=share_ax)

plt.plot(hmol_data[xaxis], (hmol_data["sum_lev"]/hmol_data["H2tot"]-1)*1.0e5, '-', c='blue',    alpha=1.0, markeredgecolor='none',markersize=msize,label='somme-1e5')
plt.legend(loc='lower right', borderaxespad=0., fontsize=tsize)

plt.xlabel(xlabel, fontsize=tsize)
plt.ylabel('H2* / H2 * 1e5', fontsize=tsize)

ax = plt.gca()
ax.set_xlim([xmin,xmax])
ax.tick_params(axis='both', which='major', labelsize=tsize)
ax.tick_params(axis='both', which='minor', labelsize=tsize)
ax.set_xscale('log')
#ax.set_yscale('log')

plt.grid(True)


# -------------------------------------------
# subplot 4
# -------------------------------------------
ax4=plt.subplot(2,2,4,sharex=share_ax)

plt.plot(chem_data[xaxis], chem_data["ab_PAH"], '-', c='blue',  alpha=1.0, markeredgecolor='none',markersize=msize,label='PAH')
plt.legend(loc='lower right', borderaxespad=0., fontsize=tsize)

plt.xlabel(xlabel, fontsize=tsize)
plt.ylabel('abundance', fontsize=tsize)

ax = plt.gca()
ax.set_xlim([xmin,xmax])
ax.set_ylim([3e-9,3e-6])
ax.tick_params(axis='both', which='major', labelsize=tsize)
ax.tick_params(axis='both', which='minor', labelsize=tsize)
ax.set_xscale('log')
ax.set_yscale('log')

plt.grid(True)



if   scptype == 'pdf':
   plt.savefig(nameout)
elif scptype == 'itv':
   plt.show()
