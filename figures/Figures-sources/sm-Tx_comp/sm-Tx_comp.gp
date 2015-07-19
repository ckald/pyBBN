reset

set xlabel "x"
set ylabel "T / T_{/Symbol n}"

set term post eps enh dl 2. lw 2.5 font "Helvetica, 25"
set ytics scale 2.
set xtics scale 1.5

#set key bottom left
set style data l

set log x

set out "Paper/figures/sm-Tx_comp.eps"

set yrange [1.:1.41]
set xrange [0.5:40.]

plot    'Data/sm/semi_sm_upd_x.dat' u 1: 3  t "",\
	'Digitized_data/sm_temp.csv' u 1: 2  t ""
	
     

set out 
!gv -scale=2 Paper/figures/sm-Tx_comp.eps