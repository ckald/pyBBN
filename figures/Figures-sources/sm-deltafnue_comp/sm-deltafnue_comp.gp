reset

set xlabel "x"
set ylabel "{/Symbol d} f_{{/Symbol n}_{e}} / f_{(eq)}"

set term post eps enh dl 2. lw 2.5 font "Helvetica, 25"
set ytics scale 2.
set xtics scale 1.5

set key left
set style data l

set label 1 at 3, 0.005 "y=3"
set label 2 at 2, 0.011 "y=5"
set label 3 at 1.3, 0.017 "y=7"

set log x

set out "Paper/figures/sm-deltafnue_comp.eps"

set yrange [0.:0.025]
set xrange [0.1:7.]


plot    'Data/sm/semi_sm_upd_x.dat' u 1:(($13)*(exp(3)+1)-1)  t "",\
	'Data/sm/semi_sm_upd_x.dat' u 1:(($14)*(exp(5)+1)-1)  t "",\
	'Data/sm/semi_sm_upd_x.dat' u 1:(($15)*(exp(7)+1)-1)  t "",\
	'Digitized_data/sm-dfnue3-semi.csv' u 1:($2) t "",\
	'Digitized_data/sm-dfnue5-semi.csv' u 1:($2) t "",\
	'Digitized_data/sm-dfnue7-semi.csv' u 1:($2) t ""#,\
	#'Data/sm/semi_sm_upd_20p_x.dat' u 1:(($13)*(exp(3)+1)-1)  t "",\
	#'Data/sm/semi_sm_upd_20p_x.dat' u 1:(($14)*(exp(5)+1)-1)  t "",\
	#'Data/sm/semi_sm_upd_20p_x.dat' u 1:(($15)*(exp(7)+1)-1)  t ""
     

set out 
!gv -scale=2 Paper/figures/sm-deltafnue_comp.eps