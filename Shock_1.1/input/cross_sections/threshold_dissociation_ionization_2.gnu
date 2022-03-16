set terminal postscript portrait enhanced color solid 12
set output "| ps2pdf - threshold_dissociation_ionization_2.pdf"


set for [i=1:5] linetype i dt i


set multiplot
set size 1.0,0.30

set noxlabel
set nologscale x
set format x " "
set xrange [1:23]
set noxtics

set ylabel "wavelength (A)"
set yrange [0:5000]
set ytics 1000

y1_H      =-1000 ; y2_H      =   91.2 ;
y1_C      =-1000 ; y2_C      =  110 ;
y1_N      =-1000 ; y2_N      =   85 ;
y1_O      =-1000 ; y2_O      =   91.04 ;
y1_Si     =-1000 ; y2_Si     =  152 ;
y1_S      =-1000 ; y2_S      =  120 ;
y1_Cl     =-1000 ; y2_Cl     =   96 ;
y1_Fe     =-1000 ; y2_Fe     =  158 ;
y1_H2     =  274 ; y2_H2     =   82 ;
y1_H2p    =  290 ; y2_H2p    =-1000 ;
y1_H3p    =  283 ; y2_H3p    =-1000 ;
y1_CH     =  358 ; y2_CH     =  117 ;
y1_CHp    =  303 ; y2_CHp    =-1000 ;
y1_CH2    =  417 ; y2_CH2    =  119 ;
y1_CH2p   =  204 ; y2_CH2p   =-1000 ;
y1_CH3    =  271 ; y2_CH3    =  126 ;
y1_CH4    =  277 ; y2_CH4    =   98 ;
y1_CH4p   = 1031 ; y2_CH4p   =-1000 ;
y1_C2     =  193 ; y2_C2     =  102 ;
y1_C2H    =  253 ; y2_C2H    =  109 ;
y1_C2H2   =  217 ; y2_C2H2   =  109 ;
y1_C2H4   =  258 ; y2_C2H4   =  118 ;
y1_C3     =  268 ; y2_C3     =  102 ;
y1_l_C3H  =  379 ; y2_l_C3H  =  144 ;
y1_c_C3H  =  289 ; y2_c_C3H  =  129 ;
y1_HC3H   =  400 ; y2_HC3H   =  138 ;
y1_l_C3H2 =  320 ; y2_l_C3H2 =  119 ;
y1_c_C3H2 =  284 ; y2_c_C3H2 =  136 ;
y1_l_C4   =  263 ; y2_l_C4   =  116 ;
y1_l_C4H  =  267 ; y2_l_C4H  =  129 ;
y1_OH     =  279 ; y2_OH     =   95 ;
y1_OHp    =  247 ; y2_OHp    =-1000 ;
y1_H2O    =  242 ; y2_H2O    =   98 ;
y1_O2     =  242 ; y2_O2     =  103 ;
y1_O2p    =  186 ; y2_O2p    =-1000 ;
y1_HO2    =  476 ; y2_HO2    =  109 ;
y1_H2O2   =  556 ; y2_H2O2   =  117 ;
y1_O3     = 1180 ; y2_O3     =   99 ;
y1_CO     =  110 ; y2_CO     =   88 ;
y1_COp    =  149 ; y2_COp    =-1000 ;
y1_CO2    =  227 ; y2_CO2    =   90 ;
y1_HCO    = 2037 ; y2_HCO    =  152 ;
y1_HCOp   =  145 ; y2_HCOp   =-1000 ;
y1_H2CO   =  361 ; y2_H2CO   =  114 ;
y1_NH     =  362 ; y2_NH     =   92 ;
y1_NHp    =  281 ; y2_NHp    =-1000 ;
y1_NH2    =  314 ; y2_NH2    =  111 ;
y1_NH3    =  301 ; y2_NH3    =  122 ;
y1_N2     =  127 ; y2_N2     =   79 ;
y1_NO     =  191 ; y2_NO     =  134 ;
y1_NO2    =  398 ; y2_NO2    =  127 ;
y1_N2O    =  735 ; y2_N2O    =   96 ;
y1_CN     =  160 ; y2_CN     =   87 ;
y1_HCN    =  243 ; y2_HCN    =   91 ;
y1_CH3OH  =  280 ; y2_CH3OH  =  113 ;
y1_SH     =  345 ; y2_SH     =  119 ;
y1_SHp    =  320 ; y2_SHp    =-1000 ;
y1_H2S    =  318 ; y2_H2S    =  119 ;
y1_CS     =  169 ; y2_CS     =  110 ;
y1_CS2    =  278 ; y2_CS2    =  123 ;
y1_OCS    =  280 ; y2_OCS    =  111 ;
y1_S2     =  283 ; y2_S2     =  133 ;
y1_SO     =  233 ; y2_SO     =  120 ;
y1_SO2    =  219 ; y2_SO2    =  100 ;
y1_SiH    =  417 ; y2_SiH    =  155 ;
y1_SiHp   =  391 ; y2_SiHp   =-1000 ;
y1_SiO    =  150 ; y2_SiO    =  108 ;
y1_HCl    =  279 ; y2_HCl    =   97 ;
y1_HClp   =  267 ; y2_HClp   =-1000 ;

h_line_01 = 6562.909442 ; 
h_line_02 = 6562.867336 ; 
h_line_03 = 6562.851769 ; 
h_line_04 = 6562.771533 ; 
h_line_05 = 6562.751807 ; 
h_line_06 = 6562.724827 ; 
h_line_07 = 6562.709702 ; 
h_line_08 = 1215.673645 ; 
h_line_09 = 1215.668237 ; 
h_line_10 = 1025.722966 ; 
h_line_11 = 1025.722855 ; 
h_line_12 = 1025.721827 ; 
h_line_13 = 1025.721825 ; 
h_line_14 = 1025.721447 ; 



set parametric ;
set trange [1:35]
set pointsize 0.7


set origin 0.0,0.66
set nolabel

set label "H     " font ",08" textcolor rgb "black" at first  1 , -800 rotate left
set label "C     " font ",08" textcolor rgb "black" at first  2 , -800 rotate left
set label "N     " font ",08" textcolor rgb "black" at first  3 , -800 rotate left
set label "O     " font ",08" textcolor rgb "black" at first  4 , -800 rotate left
set label "Si    " font ",08" textcolor rgb "black" at first  5 , -800 rotate left
set label "S     " font ",08" textcolor rgb "black" at first  6 , -800 rotate left
set label "Cl    " font ",08" textcolor rgb "black" at first  7 , -800 rotate left
set label "Fe    " font ",08" textcolor rgb "black" at first  8 , -800 rotate left
set label "H2    " font ",08" textcolor rgb "black" at first  9 , -800 rotate left
set label "H2+   " font ",08" textcolor rgb "black" at first 10 , -800 rotate left
set label "H3+   " font ",08" textcolor rgb "black" at first 11 , -800 rotate left
set label "CH    " font ",08" textcolor rgb "black" at first 12 , -800 rotate left
set label "CH+   " font ",08" textcolor rgb "black" at first 13 , -800 rotate left
set label "CH2   " font ",08" textcolor rgb "black" at first 14 , -800 rotate left
set label "CH2+  " font ",08" textcolor rgb "black" at first 15 , -800 rotate left
set label "CH3   " font ",08" textcolor rgb "black" at first 16 , -800 rotate left
set label "CH4   " font ",08" textcolor rgb "black" at first 17 , -800 rotate left
set label "CH4+  " font ",08" textcolor rgb "black" at first 18 , -800 rotate left
set label "C2    " font ",08" textcolor rgb "black" at first 19 , -800 rotate left
set label "C2H   " font ",08" textcolor rgb "black" at first 20 , -800 rotate left
set label "C2H2  " font ",08" textcolor rgb "black" at first 21 , -800 rotate left
set label "C2H4  " font ",08" textcolor rgb "black" at first 22 , -800 rotate left
set label "C3    " font ",08" textcolor rgb "black" at first 23 , -800 rotate left

set nokey
plot \
 1 , 10*y1_H      with points pt(5) lc rgb "#006400" lw(1),\
 2 , 10*y1_C      with points pt(5) lc rgb "#006400" lw(1),\
 3 , 10*y1_N      with points pt(5) lc rgb "#006400" lw(1),\
 4 , 10*y1_O      with points pt(5) lc rgb "#006400" lw(1),\
 5 , 10*y1_Si     with points pt(5) lc rgb "#006400" lw(1),\
 6 , 10*y1_S      with points pt(5) lc rgb "#006400" lw(1),\
 7 , 10*y1_Cl     with points pt(5) lc rgb "#006400" lw(1),\
 8 , 10*y1_Fe     with points pt(5) lc rgb "#006400" lw(1),\
 9 , 10*y1_H2     with points pt(5) lc rgb "#006400" lw(1),\
10 , 10*y1_H2p    with points pt(5) lc rgb "#006400" lw(1),\
11 , 10*y1_H3p    with points pt(5) lc rgb "#006400" lw(1),\
12 , 10*y1_CH     with points pt(5) lc rgb "#006400" lw(1),\
13 , 10*y1_CHp    with points pt(5) lc rgb "#006400" lw(1),\
14 , 10*y1_CH2    with points pt(5) lc rgb "#006400" lw(1),\
15 , 10*y1_CH2p   with points pt(5) lc rgb "#006400" lw(1),\
16 , 10*y1_CH3    with points pt(5) lc rgb "#006400" lw(1),\
17 , 10*y1_CH4    with points pt(5) lc rgb "#006400" lw(1),\
18 , 10*y1_CH4p   with points pt(5) lc rgb "#006400" lw(1),\
19 , 10*y1_C2     with points pt(5) lc rgb "#006400" lw(1),\
20 , 10*y1_C2H    with points pt(5) lc rgb "#006400" lw(1),\
21 , 10*y1_C2H2   with points pt(5) lc rgb "#006400" lw(1),\
22 , 10*y1_C2H4   with points pt(5) lc rgb "#006400" lw(1),\
23 , 10*y1_C3     with points pt(5) lc rgb "#006400" lw(1),\
 1 , 10*y2_H      with points pt(7) lc rgb "#8B0000" lw(1),\
 2 , 10*y2_C      with points pt(7) lc rgb "#8B0000" lw(1),\
 3 , 10*y2_N      with points pt(7) lc rgb "#8B0000" lw(1),\
 4 , 10*y2_O      with points pt(7) lc rgb "#8B0000" lw(1),\
 5 , 10*y2_Si     with points pt(7) lc rgb "#8B0000" lw(1),\
 6 , 10*y2_S      with points pt(7) lc rgb "#8B0000" lw(1),\
 7 , 10*y2_Cl     with points pt(7) lc rgb "#8B0000" lw(1),\
 8 , 10*y2_Fe     with points pt(7) lc rgb "#8B0000" lw(1),\
 9 , 10*y2_H2     with points pt(7) lc rgb "#8B0000" lw(1),\
10 , 10*y2_H2p    with points pt(7) lc rgb "#8B0000" lw(1),\
11 , 10*y2_H3p    with points pt(7) lc rgb "#8B0000" lw(1),\
12 , 10*y2_CH     with points pt(7) lc rgb "#8B0000" lw(1),\
13 , 10*y2_CHp    with points pt(7) lc rgb "#8B0000" lw(1),\
14 , 10*y2_CH2    with points pt(7) lc rgb "#8B0000" lw(1),\
15 , 10*y2_CH2p   with points pt(7) lc rgb "#8B0000" lw(1),\
16 , 10*y2_CH3    with points pt(7) lc rgb "#8B0000" lw(1),\
17 , 10*y2_CH4    with points pt(7) lc rgb "#8B0000" lw(1),\
18 , 10*y2_CH4p   with points pt(7) lc rgb "#8B0000" lw(1),\
19 , 10*y2_C2     with points pt(7) lc rgb "#8B0000" lw(1),\
20 , 10*y2_C2H    with points pt(7) lc rgb "#8B0000" lw(1),\
21 , 10*y2_C2H2   with points pt(7) lc rgb "#8B0000" lw(1),\
22 , 10*y2_C2H4   with points pt(7) lc rgb "#8B0000" lw(1),\
23 , 10*y2_C3     with points pt(7) lc rgb "#8B0000" lw(1),\
t, h_line_01 with lines notitle 'cou' lt(4) lc rgb "blue" lw(2),\
t, h_line_08 with lines notitle 'cou' lt(4) lc rgb "blue" lw(2),\
t, h_line_11 with lines notitle 'cou' lt(4) lc rgb "blue" lw(2)



set origin 0.0,0.33
set nolabel

set label "l-C3H " font ",08" textcolor rgb "black" at first  1 , -800 rotate left
set label "c-C3H " font ",08" textcolor rgb "black" at first  2 , -800 rotate left
set label "HC3H  " font ",08" textcolor rgb "black" at first  3 , -800 rotate left
set label "l-C3H2" font ",08" textcolor rgb "black" at first  4 , -800 rotate left
set label "c-C3H2" font ",08" textcolor rgb "black" at first  5 , -800 rotate left
set label "l-C4  " font ",08" textcolor rgb "black" at first  6 , -800 rotate left
set label "l-C4H " font ",08" textcolor rgb "black" at first  7 , -800 rotate left
set label "OH    " font ",08" textcolor rgb "black" at first  8 , -800 rotate left
set label "OH+   " font ",08" textcolor rgb "black" at first  9 , -800 rotate left
set label "H2O   " font ",08" textcolor rgb "black" at first 10 , -800 rotate left
set label "O2    " font ",08" textcolor rgb "black" at first 11 , -800 rotate left
set label "O2+   " font ",08" textcolor rgb "black" at first 12 , -800 rotate left
set label "HO2   " font ",08" textcolor rgb "black" at first 13 , -800 rotate left
set label "H2O2  " font ",08" textcolor rgb "black" at first 14 , -800 rotate left
set label "O3    " font ",08" textcolor rgb "black" at first 15 , -800 rotate left
set label "CO    " font ",08" textcolor rgb "black" at first 16 , -800 rotate left
set label "CO+   " font ",08" textcolor rgb "black" at first 17 , -800 rotate left
set label "CO2   " font ",08" textcolor rgb "black" at first 18 , -800 rotate left
set label "HCO   " font ",08" textcolor rgb "black" at first 19 , -800 rotate left
set label "HCO+  " font ",08" textcolor rgb "black" at first 20 , -800 rotate left
set label "H2CO  " font ",08" textcolor rgb "black" at first 21 , -800 rotate left
set label "NH    " font ",08" textcolor rgb "black" at first 22 , -800 rotate left
set label "NH+   " font ",08" textcolor rgb "black" at first 23 , -800 rotate left

plot \
 1 , 10*y1_l_C3H  with points pt(5) lc rgb "#006400" lw(1),\
 2 , 10*y1_c_C3H  with points pt(5) lc rgb "#006400" lw(1),\
 3 , 10*y1_HC3H   with points pt(5) lc rgb "#006400" lw(1),\
 4 , 10*y1_l_C3H2 with points pt(5) lc rgb "#006400" lw(1),\
 5 , 10*y1_c_C3H2 with points pt(5) lc rgb "#006400" lw(1),\
 6 , 10*y1_l_C4   with points pt(5) lc rgb "#006400" lw(1),\
 7 , 10*y1_l_C4H  with points pt(5) lc rgb "#006400" lw(1),\
 8 , 10*y1_OH     with points pt(5) lc rgb "#006400" lw(1),\
 9 , 10*y1_OHp    with points pt(5) lc rgb "#006400" lw(1),\
10 , 10*y1_H2O    with points pt(5) lc rgb "#006400" lw(1),\
11 , 10*y1_O2     with points pt(5) lc rgb "#006400" lw(1),\
12 , 10*y1_O2p    with points pt(5) lc rgb "#006400" lw(1),\
13 , 10*y1_HO2    with points pt(5) lc rgb "#006400" lw(2),\
14 , 10*y1_H2O2   with points pt(5) lc rgb "#006400" lw(2),\
15 , 10*y1_O3     with points pt(5) lc rgb "#006400" lw(2),\
16 , 10*y1_CO     with points pt(5) lc rgb "#006400" lw(2),\
17 , 10*y1_COp    with points pt(5) lc rgb "#006400" lw(2),\
18 , 10*y1_CO2    with points pt(5) lc rgb "#006400" lw(2),\
19 , 10*y1_HCO    with points pt(5) lc rgb "#006400" lw(2),\
20 , 10*y1_HCOp   with points pt(5) lc rgb "#006400" lw(2),\
21 , 10*y1_H2CO   with points pt(5) lc rgb "#006400" lw(2),\
22 , 10*y1_NH     with points pt(5) lc rgb "#006400" lw(2),\
23 , 10*y1_NHp    with points pt(5) lc rgb "#006400" lw(2),\
 1 , 10*y2_l_C3H  with points pt(7) lc rgb "#8B0000" lw(1),\
 2 , 10*y2_c_C3H  with points pt(7) lc rgb "#8B0000" lw(1),\
 3 , 10*y2_HC3H   with points pt(7) lc rgb "#8B0000" lw(1),\
 4 , 10*y2_l_C3H2 with points pt(7) lc rgb "#8B0000" lw(1),\
 5 , 10*y2_c_C3H2 with points pt(7) lc rgb "#8B0000" lw(1),\
 6 , 10*y2_l_C4   with points pt(7) lc rgb "#8B0000" lw(1),\
 7 , 10*y2_l_C4H  with points pt(7) lc rgb "#8B0000" lw(1),\
 8 , 10*y2_OH     with points pt(7) lc rgb "#8B0000" lw(1),\
 9 , 10*y2_OHp    with points pt(7) lc rgb "#8B0000" lw(1),\
10 , 10*y2_H2O    with points pt(7) lc rgb "#8B0000" lw(1),\
11 , 10*y2_O2     with points pt(7) lc rgb "#8B0000" lw(1),\
12 , 10*y2_O2p    with points pt(7) lc rgb "#8B0000" lw(1),\
13 , 10*y2_HO2    with points pt(7) lc rgb "#8B0000" lw(2),\
14 , 10*y2_H2O2   with points pt(7) lc rgb "#8B0000" lw(2),\
15 , 10*y2_O3     with points pt(7) lc rgb "#8B0000" lw(2),\
16 , 10*y2_CO     with points pt(7) lc rgb "#8B0000" lw(2),\
17 , 10*y2_COp    with points pt(7) lc rgb "#8B0000" lw(2),\
18 , 10*y2_CO2    with points pt(7) lc rgb "#8B0000" lw(2),\
19 , 10*y2_HCO    with points pt(7) lc rgb "#8B0000" lw(2),\
20 , 10*y2_HCOp   with points pt(7) lc rgb "#8B0000" lw(2),\
21 , 10*y2_H2CO   with points pt(7) lc rgb "#8B0000" lw(2),\
22 , 10*y2_NH     with points pt(7) lc rgb "#8B0000" lw(2),\
23 , 10*y2_NHp    with points pt(7) lc rgb "#8B0000" lw(2),\
t, h_line_01 with lines notitle 'cou' lt(4) lc rgb "blue" lw(2),\
t, h_line_08 with lines notitle 'cou' lt(4) lc rgb "blue" lw(2),\
t, h_line_11 with lines notitle 'cou' lt(4) lc rgb "blue" lw(2)



set origin 0.0,0.0
set nolabel

set label "NH2   " font ",08" textcolor rgb "black" at first  1 , -800 rotate left
set label "NH3   " font ",08" textcolor rgb "black" at first  2 , -800 rotate left
set label "N2    " font ",08" textcolor rgb "black" at first  3 , -800 rotate left
set label "NO    " font ",08" textcolor rgb "black" at first  4 , -800 rotate left
set label "NO2   " font ",08" textcolor rgb "black" at first  5 , -800 rotate left
set label "N2O   " font ",08" textcolor rgb "black" at first  6 , -800 rotate left
set label "CN    " font ",08" textcolor rgb "black" at first  7 , -800 rotate left
set label "HCN   " font ",08" textcolor rgb "black" at first  8 , -800 rotate left
set label "CH3OH " font ",08" textcolor rgb "black" at first  9 , -800 rotate left
set label "SH    " font ",08" textcolor rgb "black" at first 10 , -800 rotate left
set label "SH+   " font ",08" textcolor rgb "black" at first 11 , -800 rotate left
set label "H2S   " font ",08" textcolor rgb "black" at first 12 , -800 rotate left
set label "CS    " font ",08" textcolor rgb "black" at first 13 , -800 rotate left
set label "CS2   " font ",08" textcolor rgb "black" at first 14 , -800 rotate left
set label "OCS   " font ",08" textcolor rgb "black" at first 15 , -800 rotate left
set label "S2    " font ",08" textcolor rgb "black" at first 16 , -800 rotate left
set label "SO    " font ",08" textcolor rgb "black" at first 17 , -800 rotate left
set label "SO2   " font ",08" textcolor rgb "black" at first 18 , -800 rotate left
set label "SiH   " font ",08" textcolor rgb "black" at first 19 , -800 rotate left
set label "SiH+  " font ",08" textcolor rgb "black" at first 20 , -800 rotate left
set label "SiO   " font ",08" textcolor rgb "black" at first 21 , -800 rotate left
set label "HCl   " font ",08" textcolor rgb "black" at first 22 , -800 rotate left
set label "HCl+  " font ",08" textcolor rgb "black" at first 23 , -800 rotate left

plot \
 1 , 10*y1_NH2   with points pt(5) lc rgb "#006400" lw(2),\
 2 , 10*y1_NH3   with points pt(5) lc rgb "#006400" lw(2),\
 3 , 10*y1_N2    with points pt(5) lc rgb "#006400" lw(2),\
 4 , 10*y1_NO    with points pt(5) lc rgb "#006400" lw(2),\
 5 , 10*y1_NO2   with points pt(5) lc rgb "#006400" lw(2),\
 6 , 10*y1_N2O   with points pt(5) lc rgb "#006400" lw(2),\
 7 , 10*y1_CN    with points pt(5) lc rgb "#006400" lw(2),\
 8 , 10*y1_HCN   with points pt(5) lc rgb "#006400" lw(2),\
 9 , 10*y1_CH3OH with points pt(5) lc rgb "#006400" lw(2),\
10 , 10*y1_SH    with points pt(5) lc rgb "#006400" lw(2),\
11 , 10*y1_SHp   with points pt(5) lc rgb "#006400" lw(2),\
12 , 10*y1_H2S   with points pt(5) lc rgb "#006400" lw(2),\
13 , 10*y1_CS    with points pt(5) lc rgb "#006400" lw(2),\
14 , 10*y1_CS2   with points pt(5) lc rgb "#006400" lw(2),\
15 , 10*y1_OCS   with points pt(5) lc rgb "#006400" lw(2),\
16 , 10*y1_S2    with points pt(5) lc rgb "#006400" lw(2),\
17 , 10*y1_SO    with points pt(5) lc rgb "#006400" lw(2),\
18 , 10*y1_SO2   with points pt(5) lc rgb "#006400" lw(2),\
19 , 10*y1_SiH   with points pt(5) lc rgb "#006400" lw(2),\
20 , 10*y1_SiHp  with points pt(5) lc rgb "#006400" lw(2),\
21 , 10*y1_SiO   with points pt(5) lc rgb "#006400" lw(2),\
22 , 10*y1_HCl   with points pt(5) lc rgb "#006400" lw(2),\
23 , 10*y1_HClp  with points pt(5) lc rgb "#006400" lw(2),\
 1 , 10*y2_NH2   with points pt(7) lc rgb "#8B0000" lw(2),\
 2 , 10*y2_NH3   with points pt(7) lc rgb "#8B0000" lw(2),\
 3 , 10*y2_N2    with points pt(7) lc rgb "#8B0000" lw(2),\
 4 , 10*y2_NO    with points pt(7) lc rgb "#8B0000" lw(2),\
 5 , 10*y2_NO2   with points pt(7) lc rgb "#8B0000" lw(2),\
 6 , 10*y2_N2O   with points pt(7) lc rgb "#8B0000" lw(2),\
 7 , 10*y2_CN    with points pt(7) lc rgb "#8B0000" lw(2),\
 8 , 10*y2_HCN   with points pt(7) lc rgb "#8B0000" lw(2),\
 9 , 10*y2_CH3OH with points pt(7) lc rgb "#8B0000" lw(2),\
10 , 10*y2_SH    with points pt(7) lc rgb "#8B0000" lw(2),\
11 , 10*y2_SHp   with points pt(7) lc rgb "#8B0000" lw(2),\
12 , 10*y2_H2S   with points pt(7) lc rgb "#8B0000" lw(2),\
13 , 10*y2_CS    with points pt(7) lc rgb "#8B0000" lw(2),\
14 , 10*y2_CS2   with points pt(7) lc rgb "#8B0000" lw(2),\
15 , 10*y2_OCS   with points pt(7) lc rgb "#8B0000" lw(2),\
16 , 10*y2_S2    with points pt(7) lc rgb "#8B0000" lw(2),\
17 , 10*y2_SO    with points pt(7) lc rgb "#8B0000" lw(2),\
18 , 10*y2_SO2   with points pt(7) lc rgb "#8B0000" lw(2),\
19 , 10*y2_SiH   with points pt(7) lc rgb "#8B0000" lw(2),\
20 , 10*y2_SiHp  with points pt(7) lc rgb "#8B0000" lw(2),\
21 , 10*y2_SiO   with points pt(7) lc rgb "#8B0000" lw(2),\
22 , 10*y2_HCl   with points pt(7) lc rgb "#8B0000" lw(2),\
23 , 10*y2_HClp  with points pt(7) lc rgb "#8B0000" lw(2),\
t, h_line_01 with lines notitle 'cou' lt(4) lc rgb "blue" lw(2),\
t, h_line_08 with lines notitle 'cou' lt(4) lc rgb "blue" lw(2),\
t, h_line_11 with lines notitle 'cou' lt(4) lc rgb "blue" lw(2)
