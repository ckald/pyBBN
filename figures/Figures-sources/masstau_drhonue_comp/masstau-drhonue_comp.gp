reset
set log x

set xlabel "x"
set ylabel "{/Symbol d}{/Symbol r}_{{/Symbol n}_e}/{/Symbol r}_{(eq)}"

set term post eps enh dl 2. lw 2.5 font "Helvetica, 25"
set ytics scale 2.
set xtics scale 1.5

set key  left
set style data l

set label 1 at 2, 0.02 "M_{/Symbol t}=0MeV"
set label 2 at 2.4, 0.09 "M_{/Symbol t}=3MeV"
set label 3 at 2.8, 0.18 "M_{/Symbol t}=7MeV"
set label 4 at 3.2, 0.255 "M_{/Symbol t}=20MeV"


set out "Paper/figures/masstau_drhonue_comp.eps"

set yrange [0.:0.3]
set xrange [0.05:20.]

plot  'Data/sm/semi_sm_upd_x.dat' u 1:($9/0.575727-1.) t "",\
      'Digitized_data/drho0_1_imp.csv' u 1: 2 t "",\
      'Data/Tests_etc/masstau_3mev_upd_x.dat' u 1:($9/0.575727-1.) t "",\
      'Digitized_data/drho3_1_imp.csv' u 1:2 t "",\
      'Data/Tests_etc/masstau_7mev_upd_x.dat' u 1:($9/0.575727-1.) t "",\
      'Digitized_data/drho7_1_imp.csv' u 1:2 t "",\
      'Data/Tests_etc/masstau_20mev_upd2_x.dat' u 1:($9/0.575727-1.) t "",\
      'Digitized_data/drho20_1_imp.csv' u 1:2 t ""
     

set out
!gv -scale=2 Paper/figures/masstau_drhonue_comp.eps
 
