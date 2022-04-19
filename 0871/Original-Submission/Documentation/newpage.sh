
#if( $$ != 3 

html=`echo "$1.html" | tr A-Z a-z`
bml=`echo "$1.bml" | tr A-Z a-z` 

echo -n "" > $bml

echo -e  "<?xml version=\"1.0\"?>\n"  >> $bml
echo -e  "<document extend=\"style.bml\">\n"  >> $bml

echo -e  "  <pageTitle>$1</pageTitle>"  >> $bml
echo -e  "  <pageLink>"  >> $bml
echo -e  "    <define name=\"$html\" scope=\"global\"/>"  >> $bml
echo -e  "  </pageLink>\n"  >> $bml

echo -e  "  <import document=\"prec_sitemap.bml\"/>"  >> $bml

echo -e  "  <pageNavigation>"  >> $bml
echo -e  "    <maplink link=\"#seclink\" title=\"sectitle\"/>"  >> $bml
echo -e  "  </pageNavigation>\n\n"  >> $bml

echo -e  "  <content>\n"  >> $bml
echo -e  "    <section id=\"seclink\" title=\"Section sectitle\">\n"  >> $bml
echo -e  "    </section>\n"  >> $bml
echo -e  "  </content>\n"  >> $bml

echo -e  "</document>"  >> $bml
