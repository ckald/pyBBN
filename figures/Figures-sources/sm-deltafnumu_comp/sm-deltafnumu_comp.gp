reset

set xlabel "x"
set ylabel "{/Symbol d} f_{{/Symbol n}_{/Symbol m}}/ f^{eq}_{/Symbol n}"

set term post eps enh dl 2. lw 2.5 font "Helvetica, 25"
set ytics scale 2.
set xtics scale 1.5

set key  left
set style data l

set log x

set label 1 at 2, 0.002 "y=3"
set label 2 at 1.5, 0.005 "y=5"
set label 3 at 1, 0.008 "y=7"

set out "Paper/figures/sm-deltafnumu_comp.eps"

set yrange [0.:0.01]
set xrange [0.1:7.]

plot    'Data/sm/semi_sm_upd_x.dat' u 1:(($27)*(exp(3)+1)-1)  t "",\
	'Data/sm/semi_sm_upd_x.dat' u 1:(($28)*(exp(5)+1)-1)  t "",\
	'Data/sm/semi_sm_upd_x.dat' u 1:(($29)*(exp(7)+1)-1)  t "",\
	'Digitized_data/sm-dfnumu3-semi.csv' u 1:($2) t "",\
	'Digitized_data/sm-dfnumu5-semi.csv' u 1:($2) t "",\
	'Digitized_data/sm-dfnumu7-semi.csv' u 1:($2) t ""#,\
     	#'Data/sm/semi_sm_upd_20p_x.dat' u 1:(($27)*(exp(3)+1)-1)  t "",\
	#'Data/sm/semi_sm_upd_20p_x.dat' u 1:(($28)*(exp(5)+1)-1)  t "",\
	#'Data/sm/semi_sm_upd_20p_x.dat' u 1:(($29)*(exp(7)+1)-1)  t ""

#plot  'Data/sm/semi_sm_upd_x.dat' u 1:(($15)*(exp(7)+1)-1)  t "y=7",\
 #     'Digitized_data/sm-dfnue7-semi.csv' u 1:($2) t "semi",\
  #    	'Data/sm/semi_sm_upd_20p_x.dat' u 1:(($15)*(exp(7)+1)-1)  t "double"
     

set out 
!gv -scale=2 Paper/figures/sm-deltafnumu_comp.eps
