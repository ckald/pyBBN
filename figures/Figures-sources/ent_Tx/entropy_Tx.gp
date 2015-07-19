reset
set log x

set xlabel "x"
set ylabel "xT"

set term post eps enh dl 2. lw 2.5 font "Helvetica, 25"
set ytics scale 2.
set xtics scale 1.5

set key left
set style data l

set label 1 at 1.0, 1.05 "M_s=100MeV,{/Symbol t}=0.5s"
set label 2 at 0.1,1.6 "M_s=1030MeV,{/Symbol t}=0.1s"
set label 3 at 0.5, 1.9 "M_s=580MeV,{/Symbol t}=1s"


set out "Paper/figures/ent_Tx.ps"

set yrange [0.95:2]
set xrange [0.05:8.]

plot 'Data/numsm/Entropy/entropy_580MeV_1.01s_th2e_10p_upd2_x.dat' u 1: 3 t "",\
     'Data/numsm/Entropy/entropy_1030MeV_0.102s_th2tau_10p_upd2_x.dat' u 1: 3 t "",\
     'Data/numsm/Entropy/entropy_100MeV_0.5s_th2mu_10p_upd2_x.dat' u 1: 3 t "",\
     'Data/entropy_2_upd_01/Tx1030.dat' u 2: 3 t "",\
     'Data/entropy_2_upd_01/Tx100.dat' u 2: 3 t "",\
     'Data/entropy_2_upd_01/Tx580.dat' u 2: 3 t ""

#plot 'Data/numsm/Entropy/entropy_100MeV_0.5s_th2mu_10p_upd2_x.dat' u 1: 3 t "M_s=580 MeV, t=1s", 'Data/entropy_2_upd_01/Tx100.dat' u 2: 3 t ""

set out
!gv -scale=2 Paper/figures/ent_Tx.ps
 
