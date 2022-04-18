*****************************
 Welcome to the JSR Toolbox! 
*****************************

The Joint Spectral Radius of a set of matrices characterizes the maximal asymptotic rate 
of growth of a product of matrices taken in this set, when the length of the product increases.
It is known to be very hard to compute. 

In recent years, many different methods have been proposed to approximate it.
These methods have different advantages, depending on the application considered, the type of matrices 
considered, the desired accuracy or running time, etc.  

The goal of this toolbox is to provide the practitioner with the best available methods and 
an easy-to-use blackbox method: jsr.m, and propose an easy tool for the researcher to compare
different methods.  

Everyone is invited to participate to the project by including its own method, which fits the toolbox standards.
If you want your method incorporated in this toolbox (with due acknowledgement), 
send us an email at jsr.louvain@gmail.com.

Please note that some methods require SeDuMi (by default) *version 1.3*, 
it can be found at :
http://perso.uclouvain.be/raphael.jungers/sites/default/files/sedumi.zip


************************
Get started
************************

Installation: see file Install.txt 

Tutorial: Launch demo1_JSR, demo2_JSR and demo3_JSR in Matlab

Please note that some methods require SeDuMi (by default) *version 1.3*, 
it can be found at :
http://perso.uclouvain.be/raphael.jungers/sites/default/files/sedumi.zip

************************
Version
************************

The current version is v1.2 (beta).

See ChangeLog.txt for changes from past release.

Please report any bug, question or suggestion to jsr.louvain@gmail.com.

*************************
License & Credentials
*************************

See License.txt for the terms and conditions on the use and spread of this toolbox.

*************************
AUTHORS
*************************

V1 implemented by Guilllaume Vankeerberghen, Julien Hendrickx and Raphaël Jungers, 
and also makes use of several algorithms implemented by Chia-Tche Chang and Vincent Blondel.

V1.1 implemented by Damien Scieur, Julien Hendrickx and Raphaël Jungers.

V1.2 implemented by Damien Scieur, Léopold Cambier, Julien Hendrickx and Raphaël Jungers.


*************************
CITATION
*************************
If this software was useful to you, please cite the following article in your publication

G. Vankeerberghen, J. M. Hendrickx, R. M. Jungers, 
JSR: a Toolbox to Compute the Joint Spectral Radius
HSCC'14, April 15-17, 2014, Berlin, Germany
http://dx.doi.org/10.1145/2562059.2562124
