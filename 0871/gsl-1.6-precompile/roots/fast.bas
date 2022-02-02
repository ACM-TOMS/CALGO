#!/home/wschrep/bin/wsbasic

#define LIBNAME	       "libgslroots.so"
#define LIBLOADER_PATH "/home/wschrep/various/libloader-release"

function libworking( files )
begin
	LIBRARIES = "-L/home/wschrep/ArithmosRelease/Libraries/Arithmos/lib/ -lArithmos -lMpIeee -lgmp -lContinuedFraction -lSpecHardware -lSpecialValue -lRational -lEasyval"
	print "linking files = ",files," ..."
	d=run("g++ -shared -fPIC -O2 -Wall -o "+LIBNAME+" "+files+ " " + LIBRARIES + " 2>&1")
	d=run("cp "+LIBNAME+" "+LIBLOADER_PATH)
	ok = run(LIBLOADER_PATH+"/test "+LIBNAME +" 2>/dev/null")

	#print "test run='",ok,"'"
	if( ok == "" ) return 0
	else return 1
end

working = run( "cat working.log" )

if( libworking( working ) )
	print "ok new lib generated and copied"
else
	print "something went wrong with a thing called justice."

