# BUILD_MORE_TESTS.AWK
#
# USAGE
#   awk -v ofile=<file> [-v idir=<idir>] -f build_test_more.awk
#
# DESCRIPTION
#   Generates a test program, <file>, which carries out more tests of a routine
#   <oper>_rmd than those done with test_<oper>_rmd.f90. The file part of <file>
#   should be test_more_<oper>.f90.
#
#   A template test program is read from <idir>/test-more-template.txt and
#   control parameters are read from <idir>/test-more-param.txt. Default for
#   <idir> is the current directory.
#
#   The additional tests are that:
#      (a) all arrays retained from blas call remain unchanged
#      (b) not selected adjoints always remain unchanged
#      (c) selected adjoints always have the same values
#      (d) added to adjoints are indeed added to
#      (e) places inbetween incx positions remain uchanged
#
# EXAMPLE
#   (Generate "test_more_sdot.f90" in subdirectory more)
#   awk -v ofile=more/test_more_sdot.f90 -v idir=../tools -f build_test_more.awk

function error(s) {print "Error: " s; erroroccurred = 1; exit}

function append(x, y   ,nx, i) { # append y to x
  for (i in x) ++nx
  for (i in y) x[++nx] = y[i]
}

function printv(s, x  ,i) {print s,"="; for (i in x) print i ":" x[i]}

function concat2(c, x, y) { # c := [x y]
  append(c, x)
  append(c, y)
}

function concat4(c, x, y, z, t) { # c := [x y z]
  append(c, x)
  append(c, y)
  append(c, z)
  append(c, t)
}

function intersect(c, x, y  ,i,iny,j) { # c := x intersect y
  for (i in y) iny[y[i]] = 1
  for (i in x) if (x[i] in iny) c[++j] = x[i]
}

function setdiff(c, x, y  ,i,iny,j) {
  for (i in y) iny[y[i]] = 1
  for (i in x) if (!(x[i] in iny)) c[++j] = x[i]
}

BEGIN {
  if (!ofile) error("-v ofile=<file> not specified on command line")
  if (!idir) idir = "."
  OPER = ofile
  sub(".*/?test_more_", "", OPER);
  sub(".f90$", "", OPER);
  if (substr(OPER, 1, 1) != "s" || match(OPER, "[^a-z0-9]") != 0)  \
    error("Wrong name for ofile")
  parfile = idir "/test-more-param.txt"
  templatefile =  idir "/test-more-template.txt"
  while (getline < parfile == 1) {
    if ($0 == OPER) break
  }
  if ($0 != OPER) error("Operation " OPER " not found in param file")
  while (getline < parfile == 1) {
    sub(": *", ":")
    ksp = split($0, line, ":"); key = line[1]; val = line[2]
    if (ksp != 2) break
    if (key=="workspace")   split(val, wrk)
    else if (key=="sel")    split(val, sel)
    else if (key=="sels")   split(val, sels)
    else if (key=="always") split(val, alw)
    else if (key=="scalar") split(val, sca)
    else if (key=="vector") split(val, vec)
    else if (key=="matrix") split(val, mat)
    else if (key=="packed") split(val, pck)
    else if (key=="opt")    split(val, opt)
    else if (key=="inc")    split(val, inc)
    else if (key=="update") split(val, upd)
    else if (key=="destroyed") split(val, des)
    else if (key=="wrksel")    wrksel = val
    else if (key=="callrmd")   callrmd = val
    else if (key=="callrmds")  callrmds = val
    else error("Unknown key (" key ") for " OPER)
  }
  if (!callrmd) error("No callrmd-line for " OPER)
  nsela = 0
  for (i in sel) nsela++
  nsels = 0
  for (i in sels) nsels++
  concat2(para, sel, alw)
  append(para, sels)
  intersect(scaa, sca, para)
  intersect(veca, vec, para)
  intersect(mata, mat, para)
  intersect(pcka, pck, para)
  concat4(allpar, sca, vec, mat, pck)
  setdiff(pard, allpar, para)
  setdiff(par, pard, des)
  for (i in sca) isscalar[sca[i]] = 1
  for (i in vec) isvector[vec[i]] = 1
  exit
}

function addincrement(par,line  ,e,parinc,increm) {
  if (isvector[par]) {
    par1 = substr(par,1,1)
    e = substr(par,length(par),1)
    if (e == "i") par = substr(par, 1, length(par)-1)
    parinc = "inc" par1
    increm1 = "(1:n*" parinc ":" parinc ")" 
    increm2 = "(2:n*" parinc ":" parinc ")"
    gsub("<vec1>", par1, line)
    gsub("<inc1>", increm1, line)
    gsub("<inc2>", increm2, line)
  }
  else {
    gsub("<inc1>", "", line)
    gsub("<inc2>", "", line)
  }
  return line
}

function scan(template, par  ,i,j,parj,line,lines,blk,n) {
  if ($0 ~ template && $0 !~ /^ *!!/) {
    for (j in par) {
      parj = par[j]
      line = $0;
      gsub(template, parj, line)
      gsub("<j>", j, line)
      gsub("<wrksel>", wrksel, line)
      if (isscalar[parj]) gsub("<all>", "", line)
      else                gsub("<all>", "all", line)
      if (parj ~ /^trans/) gsub("<optlist>", "nt", line)
      if (parj == "side")  gsub("<optlist>", "lr", line)
      if (parj == "diag")  gsub("<optlist>", "nu", line)
      if (parj == "uplo")  gsub("<optlist>", "ul", line)
      line = addincrement(parj,line)
      blk = substr(line, 1, match(line, "[^ ]")-1)
      n = split(line, lines, ";")
      print lines[1] > ofile
      for (i=2; i<=n; i++) print (blk lines[i]) > ofile
    }
  }
}

END {
  if (erroroccurred) exit 1
  skip = 0
  while(getline < templatefile == 1) {
    if ($0 ~ /<ifnselnonzero>/ && nsela==0) skip = 1
    if ($0 ~ /<ifnselsnonzero>/ && nsels==0) skip = 1
    if ($0 ~ /<endif>/) {skip = 0; continue}
    if (skip) continue
    gsub("<callrmd>", callrmd)
    gsub("<callrmds>", callrmds)
    gsub("<nsel_array>", nsela)
    gsub("<nsel_scalar>", nsels)
    gsub("<OPER>", toupper(OPER))
    gsub("<oper>", OPER)
    if ($0 !~ /<[a-z]*>/ && $0 !~ /^ *!!/) print > ofile
    scan("<sca>", sca)
    scan("<vec>", vec)
    scan("<wrk>", wrk)
    scan("<mat>", mat)
    scan("<pck>", pck)
    scan("<scaa>", scaa)
    scan("<veca>", veca)
    scan("<mata>", mata)
    scan("<pcka>", pcka)
    scan("<par>", par)
    scan("<pard>",pard)
    scan("<para>", para)
    scan("<upd>", upd)
    scan("<opt>", opt)
    scan("<inc>", inc)
    scan("<sel>", sel)
    scan("<sels>", sels)
  }
}
