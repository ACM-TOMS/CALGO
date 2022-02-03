#this shell script exercises bertini_real, by decomposing a curve and several surfaces, with various user-options.  this script requires no user input.


#first, ensure that matlab is available as a system command.
if hash matlab 2>/dev/null; then
		rm -rf examples_run
    cp -r examples examples_run
#if it is, then run the examples.

#first, we'll run a curve decomposition, on a curve requiring randomization as written.
	cd examples_run/twisted_cubic
	echo "decomposing twisted cubic"
	../../sources/bertini_real-1.4.0/bertini_real

#decompose a torus using projection pi = <1,0,0>,<0,1,0>
	cd ../torus
	echo "decomposing torus"
	../../sources/bertini_real-1.4.0/bertini_real -pi user_provided_projection
	../../sources/bertini_real-1.4.0/sampler

#decompose the whitney umbrella using a random projection, but the unit sphere
	cd ../whitney_umbrella
	echo "decomposing whitney umbrella"
	../../sources/bertini_real-1.4.0/bertini_real -sphere user_provided_sphere

#decompose a five-variable surface, with multiple components
	cd ../nested_spheres
	echo "decomposing nested spheres"
	../../sources/bertini_real-1.4.0/bertini_real -d 2 -c -1

#return to the top directory
	cd ../../

#success
	echo "examples complete"

else
	echo "Unable to run examples, as Matlab is not installed.  Please also ensure you have the symbolic toolbox installed."
fi
