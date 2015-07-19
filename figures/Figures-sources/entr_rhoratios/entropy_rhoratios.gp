reset
set log x

set xlabel "x"
set ylabel "{/Symbol r}_s / {/Symbol r}_{SM}"

set term post eps enh dl 2. lw 2.5 font "Helvetica, 25"
set ytics scale 2.
set xtics scale 1.5

set key right
set style data l

set label 1 at 0.08, 0.5 "M_s=100MeV,{/Symbol t}=0.5s"
set label 2 at 0.031, 1.88 "M_s=1030MeV,{/Symbol t}=0.1s"
set label 3 at 0.07, 2.7 "M_s=580MeV,{/Symbol t}=1s"

set out "Paper/figures/entr_rhoratios.ps"

set xrange [0.03:5.]

plot 'Data/numsm/Entropy/entropy_580MeV_1.01s_th2e_10p_upd2_x.dat' u 1: ($30)/($5+$7+3*$9) t "",\
     'Data/numsm/Entropy/entropy_1030MeV_0.102s_th2tau_10p_upd2_x.dat' u 1: ($30)/($5+$7+3*$9)  t "",\
     'Data/numsm/Entropy/entropy_100MeV_0.5s_th2mu_10p_upd2_x.dat' u 1: ($30)/($5+$7+3*$9) t "",\
     'Data/entropy_2_upd_01/ratioRho1030.dat' u 2: 3 t "",\
     'Data/entropy_2_upd_01/ratioRho100.dat' u 2: 3 t "",\
     'Data/entropy_2_upd_01/ratioRho580.dat' u 2: 3 t ""

#plot  'Data/numsm/Entropy/entropy_1030MeV_0.102s_th2tau_10p_upd2_x.dat' u 1: ($30)/($5+$7+3*$9) t "M_s=580 MeV",  'Data/entropy_2_upd_01/ratioRho1030.dat' u 2: 3 t ""


set out
!gv -scale=2 Paper/figures/entr_rhoratios.ps
 
