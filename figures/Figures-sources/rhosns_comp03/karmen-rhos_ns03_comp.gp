reset
#set grid
set log x

set xlabel "x"
set ylabel "{/Symbol r}_s /(M_s n_s)"

set term post eps enh dl 2. lw 2.5 font "Helvetica, 25"
set ytics scale 2.
set xtics scale 1.5

set key right
set style data l

#set term pop
set term post eps enh

set out "Paper/figures/rhosns_comp03.eps"

#set yrange [0.7:1.7]
set xrange [0.027:1.3]

plot  'Data/Tests_etc/03sec_nudec12_x.dat' u ($1): ($30/$31/$1/33.9) t "Present paper result",\
     'Digitized_data/ratio_rhos_ns03_semikoz.dat' u 1:2  t "Dolgov et al."

set out
!gv -scale=2 Paper/figures/rhosns_comp03.eps
 
 
