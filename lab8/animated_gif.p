set terminal gif animate delay 20
set output "multiplot_animated.gif"
n=199
do for[i=0:(n-1)]{
    set multiplot layout 1,1
    plot[][0:0.03] 'psi_tutta.dat' every :::i::i w l lw 5
    unset multiplot
}