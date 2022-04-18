# Usage: awk -v spread=1 -f consolidate.awk times.txt
# or:    awk -f consolidate.awk times.txt  (for not spread results)
function line() {
   for (i=1; i<14*8; i++) printf "-"; print ""
}
function header() {
   line()
   h11="Time of BLAS operation"; h12="Time of RMD-operation"
   h21="(BLAS-time; ps/flop)";   h22="(multiple of BLAS-time)"
   split("                   Blas    N.of    Matrix", h20)
   split("Computer-CPU-Speed library threads size  dgemv dgemm dsyrk dtrsm dgemv dgemm dsyrk", h3);
   h31 = "dtrsm"
   printf(".\t\t\t.\t.\t.\t%s\t\t%s\n.\t\t\t", h11, h12)
   for (i in h20) printf "%s\t", h20[i]
   printf "%s\t\t%s\n", h21, h22
   for (i in h3) printf "%s\t", h3[i];
   printf "%s\n", h31
}   

BEGIN {
   OFS="\t"
   if (spread) {
      k=7
      print "SPREAD"
   }
   else {
      k=6
      print "REPEATED"
   }
   computer["mimir" ] = "Xeon E5-2650 2.0 GHz";
   computer["jotunn"] = "Xeon E5-2640 2.6 GHz";
   computer["ubuntu"] = "Ubuntu Core-i5 3.4 GHz";
   computer["tg"    ] = "Mac Core-i5 2.9 GHz"
   header()   
   while (getline < "times.txt" == 1) {
      dgemm = $k; dgemm_rmd = $(k+2); getline < "times.txt"
      dgemv = $k; dgemv_rmd = $(k+2); getline < "times.txt"
      dsyrk = $k; dsyrk_rmd = $(k+2); getline < "times.txt"
      dtrsm = $k; dtrsm_rmd = $(k+2)
      if ($1 != last1) line();
      print computer[$1],$2,$3,$4,dgemv,dgemm,dsyrk,dtrsm,dgemv_rmd,dgemm_rmd,dsyrk_rmd,dtrsm_rmd
      last1 = $1
   }
   line()
}
