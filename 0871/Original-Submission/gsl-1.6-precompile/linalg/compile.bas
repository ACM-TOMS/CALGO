#!/home/wschrep/bin/wsbasic

#define LIBNAME	       "libgsllinalg.so"
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


working    = ""
testfiles  = ""

workingold = ""
first = 1

#localpath=run("pwd")
#dummy=run("export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"+pwd)

print("--< COMPILING shared library "+LIBNAME+" >--")
print("LD_LIBRARY_PATH="+run("echo $LD_LIBRARY_PATH"))
print("Make sure this is included in the path : export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"+LIBLOADER_PATH)


while( (workingold != working) or first ) 
begin
	first=0
	workingold=working
	
	foreach F in run("ls -1")
	begin
		if( ( run( "echo "+ F +" | grep [.]o " ) != "" ) and 
		    ( run( "echo "+ working + " | grep "+F ) == "" )
 		  )
			
		begin
			print ("\n--------------------------\n\nlinking "+F+"...")
			testfiles = working+F+" "
			if( libworking(testfiles) )
				begin
					print F+" ADDED to "+LIBNAME+"."
					working=working+F+" "
				end
			else
				begin
					print F+" FAILED!"
				end
		
		end
	end
end

write( "working.log", working )


print "making the largest possible .so which works..."
if( libworking( working ) )
	print "ALL DONE! Created working shared library '"+LIBNAME+"' ."
else
	print "DOH!, FAILED loading shared library '"+LIBNAME+"' !"


