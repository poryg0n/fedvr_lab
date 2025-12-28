  set term pngcairo              # with png --- expected result 
  set xlabel "x (a.u.)" enhanced

# set yrange [-5.5:5.5] 

  set grid
  set size ratio 0.3182 1,1


 set label 1 left at graph 0.3,0.7 "N_p = 5" font ",11"  rotate by 0
 set label 2 left at graph 0.3,0.6 "N_s = 4" font ",11"  rotate by 0

 set label 4 left at graph 0.4,0.1 "xrange=[-5:5]" font ",11"  rotate by 0

# -------- Modification zone 

#  set xrange [-5.00:5.00] 
#  set xrange [-5.00:-0.00] 


   rn = ARG1 
   s = ARG2 
   i = ARG3 
   n = ARG4 

   nnbr = 5

#  set yrange [-5.00:5.00] 
#  set xrange [-5.00:-4.00] 



 
#  set output "data/run_".rn."/fedvr_sect".s."_auto.png"
#  set ylabel "{/Symbol p}^".s."_i(r)" enhanced

   set label 3 left at graph 0.4,0.8 "Sector ".s font ",11"  rotate by 0
 
#   plot "data/run_".rn."/pppic_".s."_".i."_.dat" u 7:9      with lines lw 2 lc 1 title '{/Symbol p}^'.s.'_'.i.'(1)', \
#        "data/run_".rn."/ppic_".s."_".i.".dat"   u 7:9      with lines lw 2 lc 2 title '{/Symbol p}^'.s.'_'.i.'(2)', \
#        "data/run_".rn."/qquad_pts_wts_.dat"    u 6:($13)  with linesp dt (10,5) lt 6 lc 8 title 'r^s_j'
#
#
#   set ylabel "{/Symbol p}^".s."_i(r)" enhanced
#
#   set label 3 left at graph 0.4,0.8 "Sector ".s font ",11"  rotate by 0
# 
#
#
#   set output "data/run_".rn."/fedvr_sect".s."_comp_func_and_deriv.png"
#   plot "data/run_".rn."/ppic_".s."_".i.".dat"  u 7:($9)   with lines lw 2 lc 4 title '{/Symbol p}^'.s.'_'.i.'(1)', \
#        "data/run_".rn."/ppic_".s."_".i.".dat"  u 7:($11)  with lines lw 2 lc 5 title 'd{/Symbol p}^'.s.'_'.i.'(1)/dx', \
#        "data/run_".rn."/ppi_".s."_".i.".dat"   u 5:($11)  with linesp dt(10,5) lt 6 lc 3 title 'ref', \
#        "data/run_".rn."/quad_pts_wts.dat"      u 6:($13)  with linesp dt (10,5) lt 6 lc 8 title 'r^s_j'
#
#
##   plot "data/run_".rn."/ppi_".s."_1.dat" u 6:7 with lines lw 2 lc 1 title '{/Symbol p}^'.s.'_1', \
##        "data/run_".rn."/ppi_".s."_2.dat" u 6:7 with lines lw 2 lc 2 title '{/Symbol p}^'.s.'_2', \
##        "data/run_".rn."/ppi_".s."_3.dat" u 6:7 with lines lw 2 lc 3 title '{/Symbol p}^'.s.'_3', \
##        "data/run_".rn."/ppi_".s."_4.dat" u 6:7 with lines lw 2 lc 4 title '{/Symbol p}^'.s.'_4', \
##        "data/run_".rn."/ppi_".s."_5.dat" u 6:7 with lines lw 2 lc 5 title '{/Symbol p}^'.s.'_5', \
##        "data/run_".rn."/quad_pts_wts.dat" u 4:8 with linesp  dt (10,5) lt 6 lc 8 title 'r^s_j'



    set output "data/run_".rn."/fedvr_sect".s."_g_pis_all.png"
    set ylabel "{/Symbol P}_n(r)" enhanced
 
   set label 3 left at graph 0.4,0.8 "Sector ".s font ",11"  rotate by 0
   set key top left
  
    plot "data/run_".rn."/ppic_".s."_2.dat" u 7:($9) with lines lw 2 lc 1 title sprintf("{/Symbol P}_{%g}", 1+(nnbr-1)*(s-1)), \
         "data/run_".rn."/ppic_".s."_3.dat" u 7:($9) with lines lw 2 lc 2 title sprintf("{/Symbol P}_{%g}", 2+(nnbr-1)*(s-1)), \
         "data/run_".rn."/ppic_".s."_4.dat" u 7:($9) with lines lw 2 lc 3 title sprintf("{/Symbol P}_{%g}", 3+(nnbr-1)*(s-1)), \
         "data/run_".rn."/ppic_".s."_5.dat" u 7:($9) with lines lw 2 lc 4 title sprintf("{/Symbol P}_{%g}", 4+(nnbr-1)*(s-1)), \
         "data/run_".rn."/ppi_".s."_".i.".dat" u 5:($8) lc 3 lt 7 lw 3 title sprintf("ref_{%g}", i+(nnbr-1)*(s-1)-1), \
         "data/run_".rn."/quad_pts_wts.dat" u 6:($13) with linesp  dt (10,5) lt 6 lc 8 title 'r^s_j'


    set output "data/run_".rn."/fedvr_sect".s."_dg_pis_all.png"
    set ylabel "d{/Symbol P}_n(x)/dx" enhanced
 
   set label 3 left at graph 0.4,0.8 "Sector ".s font ",11"  rotate by 0
  
    plot "data/run_".rn."/ppic_".s."_2.dat" u 7:($11) with lines lw 2 lc 1 title sprintf("d{/Symbol P}_{%g}/dx", 1+(nnbr-1)*(s-1)), \
         "data/run_".rn."/ppic_".s."_3.dat" u 7:($11) with lines lw 2 lc 2 title sprintf("d{/Symbol P}_{%g}/dx", 2+(nnbr-1)*(s-1)), \
         "data/run_".rn."/ppic_".s."_4.dat" u 7:($11) with lines lw 2 lc 3 title sprintf("d{/Symbol P}_{%g}/dx", 3+(nnbr-1)*(s-1)), \
         "data/run_".rn."/ppic_".s."_5.dat" u 7:($11) with lines lw 2 lc 4 title sprintf("d{/Symbol P}_{%g}/dx", 4+(nnbr-1)*(s-1)), \
         "data/run_".rn."/ppi_".s."_".i.".dat" u 5:($12) lc 3 lt 7 title sprintf("ref_{%g}", i+(nnbr-1)*(s-1)-1), \
         "data/run_".rn."/quad_pts_wts.dat" u 6:($13) with linesp  dt (10,5) lt 6 lc 8 title 'r^s_j'



#   set output "data/run_".rn."/fedvr_sect".s."_dg_pis_all.png"
#   set ylabel "{/Symbol P}_n(r)" enhanced
#
#  set label 3 left at graph 0.4,0.8 "Sector ".s font ",11"  rotate by 0
# 
#   plot "data/run_".rn."/ppi_".s."_2.dat" u 5:($11) with lines lw 2 lc 1 title sprintf("{/Symbol P}_{%g}", 2+(nnbr-1)*(s-1)), \
#        "data/run_".rn."/ppi_".s."_3.dat" u 5:($11) with lines lw 2 lc 2 title sprintf("{/Symbol P}_{%g}", 3+(nnbr-1)*(s-1)), \
#        "data/run_".rn."/ppi_".s."_4.dat" u 5:($11) with lines lw 2 lc 3 title sprintf("{/Symbol P}_{%g}", 4+(nnbr-1)*(s-1)), \
#        "data/run_".rn."/ppi_".s."_5.dat" u 5:($11) with lines lw 2 lc 4 title sprintf("{/Symbol P}_{%g}", 5+(nnbr-1)*(s-1)), \
#        "data/run_".rn."/quad_pts_wts.dat" u 6:($13) with linesp  dt (10,5) lt 6 lc 8 title 'r^s_j'



    unset label 3
    set ylabel "{/Symbol Y}(x,t)" enhanced
    set output "data/run_".rn."/fedvr_wavefunc_.png"
    plot "analytic_wfc.dat"  u 1:($2)   with lines lw 2 lc 6 title '{/Symbol Y}_0(x)', \
         "analytic_wfc.dat"  u 1:($2)   pt 7 lc 8 title 'numerical', \
         "data/run_".rn."/quad_pts_wts.dat"      u 6:($12)  with linesp dt (10,5) lt 6 lc 8 title 'r^s_j'


     set label 1 left at graph 0.6,0.7 "N_p = 5" font ",11"  rotate by 0
     set label 2 left at graph 0.6,0.6 "N_s = 4" font ",11"  rotate by 0

     set label 4 left at graph 0.6,0.1 "xrange=[-5:5]" font ",11"  rotate by 0
 
    set ylabel "d{/Symbol Y}(x)/dx" enhanced
    set output "data/run_".rn."/fedvr_deriv_wavefunc_.png"
    plot "analytic_dwfc.dat"  u 1:($2)   with lines lw 2 lc 6 title 'd{/Symbol Y}_0(x)/dx', \
         "analytic_dwfc.dat"  u 1:($2)   pt 7 lc 8 title 'ref', \
         "analytic_dwfc.dat"  u 1:($2)   pt 7 lc 7 title '', \
         "data/run_".rn."/quad_pts_wts.dat"      u 6:($12)  with linesp dt (10,5) lt 6 lc 8 title 'r^s_j'


    set ylabel "d{/Symbol Y}(x)/dx" enhanced
    set output "data/run_".rn."/fedvr_deriv_wavefunc_fix.png"
    plot "fort.81"  u 1:($2)   with lines lw 2 lc 6 title 'd{/Symbol Y}_0(x)/dx', \
         "fort.74"  u 1:($2)   pt 7 lc 8 title 'ref', \
         "data/run_".rn."/quad_pts_wts.dat"      u 6:($12)  with linesp dt (10,5) lt 6 lc 8 title 'r^s_j'
 
 
 
#   set output "data/run_".rn."/fedvr_sect".s."_dpi_fin_diff.png"
#   set ylabel "d{/Symbol P}_n(r)/dr (fin. diff.)" enhanced
#
#   set label 3 left at graph 0.4,0.8 "Sector ".s font ",11"  rotate by 0
# 
#   plot "data/run_".rn."/dpi_".s."_1.dat" u 4:5 with lines lw 2 lc 1 title sprintf("d{/Symbol P}_{%g}/dr", 1+(nnbr-1)*(s-1)), \
#        "data/run_".rn."/dpi_".s."_2.dat" u 4:5 with lines lw 2 lc 2 title sprintf("d{/Symbol P}_{%g}/dr", 2+(nnbr-1)*(s-1)), \
#        "data/run_".rn."/dpi_".s."_3.dat" u 4:5 with lines lw 2 lc 3 title sprintf("d{/Symbol P}_{%g}/dr", 3+(nnbr-1)*(s-1)), \
#        "data/run_".rn."/dpi_".s."_4.dat" u 4:5 with lines lw 2 lc 4 title sprintf("d{/Symbol P}_{%g}/dr", 4+(nnbr-1)*(s-1)), \
#        "data/run_".rn."/dpi_".s."_5.dat" u 4:5 with lines lw 2 lc 5 title sprintf("d{/Symbol P}_{%g}/dr", 5+(nnbr-1)*(s-1)), \
#        "data/run_".rn."/quad_pts_wts.dat" u 4:8 with linesp  dt (10,5) lt 6 lc 8 title 'r^s_j'
#
#
##   plot "data/run_".rn."/dpi_".s."_1.dat" u 4:7  lc 1 title sprintf('d{/Symbol P}_{%g}/dr', 1+(nnbr-1)*(s-1)), \
##        "data/run_".rn."/dpi_".s."_2.dat" u 4:7  lc 2 title sprintf('d{/Symbol P}_{%g}/dr', 2+(nnbr-1)*(s-1)), \
##        "data/run_".rn."/dpi_".s."_3.dat" u 4:7  lc 3 title sprintf('d{/Symbol P}_{%g}/dr', 3+(nnbr-1)*(s-1)), \
##        "data/run_".rn."/dpi_".s."_4.dat" u 4:7  lc 4 title sprintf('d{/Symbol P}_{%g}/dr', 4+(nnbr-1)*(s-1)), \
##        "data/run_".rn."/dpi_".s."_5.dat" u 4:7  lc 5 title sprintf('d{/Symbol P}_{%g}/dr', 5+(nnbr-1)*(s-1)), \
##        "data/run_".rn."/dp2_".s."_".i.".dat" u 5:10 title sprintf('ref')
##
###       "data/run_".rn."/dp2_".s.".dat" u 5:9 with linesp dt (1,1) lt 7 lc 8 title sprintf('ref'), \
#
#
#
#
#
#
#
#
#
##  set yrange [-.50:0.50] 
#   set output "data/run_".rn."/wavfun_analyt_real.png"
#   set ylabel "Re(d{/Symbol Y}(x,t)/dx)" enhanced
#
#   set label 3 left at graph 0.4,0.8 "{/Symbol S}_{s=1}^".s." {/Symbol S}_{n=1}^".n font ",11"  rotate by 0
#
#    
#   plot "data/run_".rn."/wavfun_and_deriv_cont.dat" u 2:3 with lines lw 1 lc 6 title '{/Symbol Y}(x)', \
#        "data/run_".rn."/wavfun_and_deriv_cont.dat" u 2:5 with lines lw 1 lc 8 title 'd{/Symbol Y}(x)/dx', \
#        "data/run_".rn."/wavfun_and_deriv_disc.dat" u 2:5 lt 1 lc 7 title 'ref (sum)',            \
#        "data/run_".rn."/quad_pts_wts.dat" u 4:8 with linesp  dt (10,5) lt 6 lc 8 title 'r^s_j'
#
##  set yrange [-.50:0.50] 
#   set output "data/run_".rn."/wavfun_analyt_imag.png"
#   set ylabel "Im(d{/Symbol Y}(x,t)/dx)" enhanced
#
#   set label 3 left at graph 0.4,0.8 "{/Symbol S}_{s=1}^".s." {/Symbol S}_{n=1}^".n font ",11"  rotate by 0
#
#    
#   plot "data/run_".rn."/wavfun_and_deriv_cont.dat" u 2:4 with lines lw 2 lc 6 title '{/Symbol Y}(x)', \
#        "data/run_".rn."/wavfun_and_deriv_cont.dat" u 2:6 with lines lw 2 lc 8 title 'd{/Symbol Y}(x)/dx', \
#        "data/run_".rn."/wavfun_and_deriv_disc.dat" u 2:6 lt 5 lc 7 title 'ref (sum)',            \
#        "data/run_".rn."/quad_pts_wts.dat" u 4:8 with linesp  dt (10,5) lt 6 lc 8 title 'r^s_j', \
#
#
#
#
##  set yrange [-.50:0.50] 
#   set output "data/run_".rn."/wavfun_analyt0_real.png"
#   set ylabel "Re(d{/Symbol Y}(x,t)/dx)" enhanced
#
#   set label 3 left at graph 0.4,0.8 "{/Symbol S}_{s=1}^".s." {/Symbol S}_{n=1}^".n font ",11"  rotate by 0
#
#    
#   plot "data/run_".rn."/wavfun_and_deriv_cont.dat" u 2:3 with lines lw 1 lc 6 title '{/Symbol Y}_0(x)', \
#        "data/run_".rn."/wavfun_and_deriv_cont.dat" u 2:8 with lines lw 1 lc 7 title 'num',
#        "data/run_".rn."/quad_pts_wts.dat" u 4:8 with linesp  dt (10,5) lt 6 lc 8 title 'r^s_j'
#
##  plot "fort.889" u 4:5 with lines lw 2 lc 1, \
##       "fort.889" u 4:6 with lines lw 2 lc 2
#
#
##   set output "fedvr_dpi_".s."_".i.".png"
##   set ylabel "d{/Symbol P}_{".n."}(r)/dr" enhanced
##   set yrange [-2.00:2.00] 
#    set xrange [-5.00:-4.00] 
## 
## # set arrow nohead
## 
## set label 3 left at graph 0.4,0.9 "s=".s font ",11"  rotate by 0
## set label 4 left at graph 0.4,0.8 "i=".i font ",11"  rotate by 0
## set label 5 left at graph 0.4,0.7 "n=".n font ",11"  rotate by 0
## 
##   plot "pi_".s."_".i.".dat" u 6:8 with lines  lw 2 lc 2 title '{/Symbol P}_{'.n.'}', \
##        "dpi_".s."_".i.".dat" u 4:7 with lines  lw 3 lc 6 title 'd{/Symbol P}_{'.n.'}/dr', \
##        "quad_pts_wts.dat" u 4:8 with linesp  dt (10,5) lt 6 lc 8 title 'r^s_j'
###       "dpi_".s."_".i.".dat" u 4:8 with lines  lw 1 lc 4 title 'fourier', \
###       "dpi_".s."_".i.".dat" u 4:5 with lines  dt (8,10) lw 3 lc 1 title '(fin. diff.)', \
###       "dpi_".s."_".i.".dat" u 4:6 with lines  lw 2 lc 5 title '', \
###       "dpi2_".s."_".i.".dat" u 4:7 with linesp  dt (10,2) lt 7 lw 1 lc 8 title 'ref', \
###       "dpi2_2_5.dat" u 4:7 with linesp  dt (10,2) lt 7 lw 1 lc 8 title '', \
###       "dpi_3_1.dat" u 4:6 with lines  lw 2 lc 5 title '', \
####      "pi_1_5.dat" u 6:8 with lines dt (10,3) lt 7 lw 2 lc 8 title '{/Symbol P}_5', \
##
####
####   set output "fedvr_dg_pi_n_5_.png"
####   set ylabel "d{/Symbol p}^s_i(r)/dr" enhanced
####   set yrange [-2.00:2.00] 
####   set xrange [-5.00:1.66] 
#####  set log y
#### 
#### # set arrow nohead
#### 
#### set label 3 left at graph 0.4,0.8 "s=1" font ",11"  rotate by 0
#### set label 4 left at graph 0.4,0.7 "i=1" font ",11"  rotate by 0
#### 
####   plot "dpi_1_1.dat" u 4:5 with lines  lw 1 lc 2 title '(fin.diff.)', \
####        "dpi_1_1.dat" u 4:6 with lines  lw 1 lc 7 title '(analyt.)', \
####        "dpi_1_1.dat" u 4:7 with lines  lw 1 lc 6 title 'd{/Symbol P}_5/dr', \
####        "quad_pts_wts.dat" u 4:8 with linesp  dt (10,5) lt 6 lc 8 title 'r^s_j'
#####       "dpi_2_1.dat" u 4:6 with lines  lw 1 lc 6 title '{/Symbol p}^2_1', \
#####       "dpi_2_1.dat" u 4:7 with lines  lw 2 lc 8 title 'd{/Symbol P}_5/dr', \
#####       "dpi_1_1.dat" u 6:8 with lines  dt (8,10) lw 1 lc 7 title '(fin. diff.)', \
#####       "dpi2_1_1.dat" u 4:7 with linesp  dt (10,2) lt 6 lw 2 lc 8 title 'ref', \
#### #      "pi_1_5.dat" u 6:8 with lines dt (10,3) lt 7 lw 2 lc 8 title '{/Symbol P}_5', \
