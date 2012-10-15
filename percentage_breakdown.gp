# percentage_breakdown.gp
# Started by Matthew Demas
# on 27 September 2012
##############################################################################
# to run: 
#  1. make directory 'figures' below this script
#  2. type "gnuplot percentage_breakdown.gp"
#  3. then run (to generate publication quality images)
#   for i in `ls ./figures/*.eps | gawk -F'.' '{print $1}'`;
#   do convert -geometry 3600 -density 600 -units pixelsperinch \
#       ./figures/$i.eps ./figures $i.png;
#   done
##############################################################################
# - Three plots showing the histogram of the breakdown for the 1st
#   coordination sphere
# - No carbon
# - make publication quality
# 1. Na, Mg, K, Ca, Zn 
# 2. Fe2+, Fe3+
# 3. All elements for:
#    (have cutoff 0.3)
#    a. Zn
#    b. Fe2+
#    c. Fe3+
##############################################################################

# file/output options
set terminal postscript color enhanced "Helvetica" 40 eps solid \
    size 12in, 7.5in lw 5
#set terminal wxt
set datafile separator "|"

# general graph formatting options
set ytics nomirror 
set xtics nomirror
set tics out
set format x "%.1f"
set pointsize 2.5
set xlabel "Bond Valence Weight"
set ylabel "Number of Sites"
set key spacing 1.5

########################################
# 1st Plot
########################################
set output "./figures/NaMgKCaZn.eps"
# set thresholds for Zn(II), Fe(II/III)
set style arrow 5 nohead front lc rgb 'black'
set arrow from 0.3,graph 0 to 0.3,graph 0.90 as 5
set label "weight threshold\nfor Zn" at 0.3,graph 0.97 center
# set thresholds for Na,Mg,K,Ca
set style arrow 5 nohead front lw 2 lt 4 
set arrow from 0.5,graph 0 to 0.5,graph 0.3  as 5
set label "weight threshold\nfor Na, Mg, K, Ca" at 0.5,graph 0.37 center
plot \
'../allmetals_percentage_query_file.csv' u 1:8 w linespoints t 'Na',\
'../allmetals_percentage_query_file.csv' u 1:6 w linespoints t 'Mg',\
'../allmetals_percentage_query_file.csv' u 1:7 w linespoints t 'K',\
'../allmetals_percentage_query_file.csv' u 1:3 w linespoints t 'Ca',\
'../allmetals_percentage_query_file.csv' u 1:2 w linespoints t 'Zn'
unset arrow
unset label
#pause -1


########################################
# 2nd Plot
########################################
set output "./figures/Fe2Fe3.eps"
plot \
'../allmetals_percentage_query_file.csv' u 1:4 w linespoints pt 2 t 'Fe(II)',\
'../allmetals_percentage_query_file.csv' u 1:5 w linespoints pt 5 t 'Fe(III)'


########################################
# 3rd Plot
########################################
#set terminal postscript color enhanced "Helvetica" 40 eps solid \
#    size 12in, 7.5in lw 5
####################
# 3rd Plot part "A"
# Plot for all Zinc 2 interactions except for Carbon
####################
set output "./figures/Zn_ele.eps"
plot \
'../zinc2_percentage_query_file.csv' u 1:4 w linespoints  lc rgb "blue"    t 'N',\
'../zinc2_percentage_query_file.csv' u 1:5 w linespoints  lc rgb "red"     t 'O',\
'../zinc2_percentage_query_file.csv' u 1:6 w linespoints  lc rgb "#22CCCC" t 'F' ,\
'../zinc2_percentage_query_file.csv' u 1:7 w linespoints  lc rgb "#C73F17" t 'P',\
'../zinc2_percentage_query_file.csv' u 1:8 w linespoints  lc rgb "#DDBB40" t 'S',\
'../zinc2_percentage_query_file.csv' u 1:9 w linespoints  lc rgb "#20EE20" t 'Cl',\
'../zinc2_percentage_query_file.csv' u 1:10 w linespoints lc rgb "#6B8E23" t 'Se',\
'../zinc2_percentage_query_file.csv' u 1:11 w linespoints lc rgb "#B22222" t 'Br',\
'../zinc2_percentage_query_file.csv' u 1:12 w linespoints lc rgb "#4B0082" t 'I'

####################
# 3rd Plot part "B"
# Plot for all Iron 2 interactions except for Carbon
####################
set output "./figures/Fe2_ele.eps"
set ytics 0,50
set yrange[0:120]
plot \
'../iron2_percentage_query_file_no-carbon.csv' u 1:3 w linespoints  lc rgb "blue"     t 'N',\
'../iron2_percentage_query_file_no-carbon.csv' u 1:4 w linespoints  lc rgb "red"      t 'O',\
'../iron2_percentage_query_file_no-carbon.csv' u 1:5 w linespoints  lc rgb "#22CCCC"  t 'F' ,\
'../iron2_percentage_query_file_no-carbon.csv' u 1:6 w linespoints  lc rgb "#FF8000"  t 'P',\
'../iron2_percentage_query_file_no-carbon.csv' u 1:7 w linespoints  lc rgb "#DDBB40"  t 'S',\
'../iron2_percentage_query_file_no-carbon.csv' u 1:8 w linespoints  lc rgb "#20EE20"  t 'Cl',\
'../iron2_percentage_query_file_no-carbon.csv' u 1:9 w linespoints lc rgb "#6B8E23"  t 'Se',\
'../iron2_percentage_query_file_no-carbon.csv' u 1:10 w linespoints lc rgb "#B22222"  t 'Br',\
'../iron2_percentage_query_file_no-carbon.csv' u 1:11 w linespoints lc rgb "#4B0082"  t 'I'

####################
# 3rd Plot part "C"
# Plot for all Iron 3 interactions except for Carbon
####################
set output "./figures/Fe3_ele.eps"
set yrange [0:200]
plot \
'../iron3_percentage_query_file.csv' u 1:4 w linespoints  lc rgb "blue"     t 'N',\
'../iron3_percentage_query_file.csv' u 1:5 w linespoints  lc rgb "red"      t 'O',\
'../iron3_percentage_query_file.csv' u 1:6 w linespoints  lc rgb "#22CCCC"  t 'F' ,\
'../iron3_percentage_query_file.csv' u 1:7 w linespoints  lc rgb "#FF8000"  t 'P',\
'../iron3_percentage_query_file.csv' u 1:8 w linespoints  lc rgb "#DDBB40"  t 'S',\
'../iron3_percentage_query_file.csv' u 1:9 w linespoints  lc rgb "#20EE20"  t 'Cl',\
'../iron3_percentage_query_file.csv' u 1:10 w linespoints lc rgb "#6B8E23"  t 'Se',\
'../iron3_percentage_query_file.csv' u 1:11 w linespoints lc rgb "#B22222"  t 'Br',\
'../iron3_percentage_query_file.csv' u 1:12 w linespoints lc rgb "#4B0082"  t 'I'
