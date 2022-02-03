#/bin/bash

#converts old html stuff like <br> into correct <br/>
mv $1 tmp
cat tmp | sed -e 's/<\([bB][rR][^/]*\)>/<\1\/>/g' > $1
rm tmp


#replace <area ...> by <area .../>
mv $1 tmp
cat tmp | sed -e 's/<\([aA][rR][eE][aA][^/]*\)>/<\1\/>/g' | sed -e 's/\(<[aA][rR][eE][aA].*NOHREF[^/]*\)\/>/\1=""\/>/g'> $1
rm tmp


#replace <img ...> by <img .../>
mv $1 tmp
cat tmp | sed -e 's/<\([iI][mM][gG].*\)>/<\1\/>/g' > $1
rm tmp


#replace <something //> by <something /> -> for when it was run twice...
mv $1 tmp
cat tmp | sed -e 's/\(<.*\/\)\/>/\1>/g' > $1
rm tmp


#put <SCRIPT> tags to lowercase...
mv $1 tmp
cat tmp | sed -e 's/<[sS][cC][rR][iI][pP][tT]\(.*>\)/<script \1/g' | sed -e 's/<\/[sS][cC][rR][iI][pP][tT]\(.*>\)/<\/script \1/g' > $1
rm tmp

