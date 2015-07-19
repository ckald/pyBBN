reset
#set grid

set xlabel "y"
set ylabel "{/Symbol d}f_{{/Symbol n}_e} y^2"

set term post eps enh dl 2. lw 2.5 font "Helvetica, 25"
set ytics scale 2.
set xtics scale 1.5

set key  right
set style data l

set label 1 at 2, 0.012 "M_{/Symbol t}=1MeV"
set label 2 at 2.2, 0.033 "M_{/Symbol t}=3MeV"
set label 3 at 2.6, 0.052 "M_{/Symbol t}=7MeV"
set label 4 at 3.5, 0.072 "M_{/Symbol t}=20MeV"

set out "Paper/figures/masstau_fnue_comp.eps"

set yrange [0.:0.08]
set xrange [0.:10.]
#set xrange [2.5:3.5]

plot    'Data/Tests_etc/masstau_1mev_upd.dat' u 1: ($4-1./(exp($1)+1))*($1)*($1) t "",\
	'Data/Tests_etc/masstau_3mev_upd.dat' u 1: ($4-1./(exp($1)+1))*($1)*($1) t "",\
	'Data/Tests_etc/masstau_7mev_upd.dat' u 1: ($4-1./(exp($1)+1))*($1)*($1) t "",\
	'Data/Tests_etc/masstau_20mev_upd.dat' u 1: ($4-1./(exp($1)+1))*($1)*($1) t "",\
	'Digitized_data/df_1_imp.csv' u 1:2 t "",\
	'Digitized_data/df_3_imp.csv' u 1:2 t "",\
	'Digitized_data/df_7_imp.csv' u 1:2 t "",\
	'Digitized_data/df_20_imp.csv' u 1:2 t ""

#plot 	'Data/Tests_etc/masstau_7mev_upd.dat' u 1: ($4-1./(exp($1)+1))*($1)*($1) t "",\
	#'Digitized_data/df_7_imp.csv' u 1:2 t "M=7 MeV"
#'Data/corr_masstau_3mev_decay2test.dat' u 1:($4-1./(exp($1)+1))*($1)*($1)  t "nu_{e}, M=3 MeV,decay2,corr",\
 #'Data/semi_sm_masstau_3mev.dat' u 1:($4-1./(exp($1)+1))*($1)*($1)  t "nu_{e}, M=3 MeV",\     
 #'Data/semi_sm_masstau_7mev_test.dat' u 1:($4-1./(exp($1)+1))*($1)*($1)  t "nu_{e}, M=7 MeV"
      #  'Data/semi_sm_x.dat' u 1: ($4-1./(exp($1)+1))*($1)*($1) t "nu_{e}, SM",\
     # 'Data/semi_sm_masstau_7mev_test_x.dat' u 1:($9/0.575727-1.) t "nu_{e}, M=7 MeV",\
     # 'Data/semi_sm_masstau_20mev_test_x.dat' u 1:($9/0.575727-1.) t "nu_{e},M=20 MeV",\
      # 'Data/semi_sm_masstauTdec100_20mev_test_x.dat' u 1:($9/0.575727-1.) t "nu_{e},Tdec=100MeV M=20 MeV"
     

set out
!gv -scale=2 Paper/figures/masstau_fnue_comp.eps