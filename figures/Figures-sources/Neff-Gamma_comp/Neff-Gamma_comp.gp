reset
set log x

set xlabel "{/Symbol G}, s^{-1}"
set ylabel "N_{/Symbol n}"

set term post eps enh dl 2. lw 2.5 font "Helvetica, 25"
set ytics scale 2.
set xtics scale 1.5

set style data l

set label 1 at 5, 0.6 "Present paper result"
set label 2 at 0.5,0.5 "Hannestad"
set label 3 at 3., 2.5 "Kawasaki et al."

set out "Paper/figures/Neff_Gamma_comp.eps"

set yrange [0.:3.]
set xrange [0.2:100.]
#set xrange [7:10]

plot  'Data/Tests_etc/kks4_summary.dat' u 1: 2 t "",\
      'Digitized_data/hanne_g.csv' t "",\
      'Digitized_data/kks_tau.csv' t ""

set out
!gv -scale=2 Paper/figures/Neff_Gamma_comp.eps
  
