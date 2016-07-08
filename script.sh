mv _images/ images
mv _static/ static
mv _sources/ sources
sed -i "s/_static/static/" *.html
sed -i "s/_static/static/" static/websupport.js
sed -i "s/_sources/sources/" static/searchtools.js
sed -i "s/_sources/sources/" *.html
sed -i "s/_images/images/" inoutput.html
sed -i "s/_images/images/" inoutput.html
sed -i "s/_images/images/" handson.html
sed -i "s/_images/images/" handson.html

