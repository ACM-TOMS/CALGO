function main() {
   print
   print "   Example 1: Simple Substitution"
   print "   ------------------------------"
   br = "brown"
   print "   The quick " br " fox."
   print
   print "   Example 2: Substitution inside a String"
   print "   ---------------------------------------"
   r = "row"
   print "   The quick b" r "n fox."
   print
   print "   Example 3: Expression Substitution"
   print "   ----------------------------------"
   a = 4
   b = 3
   print "   The quick " 2*a + b " foxes."
   print
   print "   Example 4: Macros References inside a Macro"
   print "   -------------------------------------------"
   M["fox"] = "$[q] $[b] $[f]"
   M["q"] = "quick"
   M["b"] = "brown"
   M["f"] = "fox"
   print "   The " eval(M["fox"]) "."
   print
   print "   Example 5: Array Reference Substitution"
   print "   ---------------------------------------"
   x[7] = "brown"
   b = 3
   print "   The quick " x[2*b+1] " fox."
   print
   print "   Example 6: Function Reference Substitution"
   print "   ------------------------------------------"
   print "   The quick " color(1,2) " fox."
   print
   print "   Example 7: Substitution of Special Characters"
   print "   ---------------------------------------------"
   print "\#  The \$ quick \\ brown $# fox. $$"
}
function color(i,j) {
   print "   The lazy dog."
   if (i == j)
      return "blue"
   else
      return "brown"
}

function eval(inp	,isplb,irb,out,name) {

   splb = SP "["
   out = ""

   while( isplb = index(inp, splb) ) {
      irb = index(inp, "]")
      if ( irb == 0 ) {
         out = out substr(inp,1,isplb+1)
         inp = substr( inp, isplb+2 )
      } else {
         name = substr( inp, isplb+2, irb-isplb-2 )
         sub( /^ +/, "", name )
         sub( / +$/, "", name )
         out = out substr(inp,1,isplb-1) eval(M[name])
         inp = substr( inp, irb+1 )
      }
   }

   out = out inp

   return out
}
BEGIN {
   SP = "$"
   main()
   exit
}
