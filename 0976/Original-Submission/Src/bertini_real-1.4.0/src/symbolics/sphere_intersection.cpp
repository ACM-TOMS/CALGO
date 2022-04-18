#include "symbolics/sphere_intersection.hpp"



void create_sphere_system(boost::filesystem::path input_file, boost::filesystem::path output_file,
                          comp_mp sphere_radius,
                          vec_mp sphere_center,
                          const WitnessSet & W)
{
	
#ifdef functionentry_output
	std::cout << "create_sphere_system" << std::endl;
#endif
	
	
	
	
	
	
	// a bit of error checking
    
	
	if (W.num_var_names()==0) {
		throw std::logic_error("trying to create a sphere system, but witness set does not have the variable names.");
	}
	
	
	// got here, so ok to continue.
	
	
	int *declarations = NULL;
	
	partition_parse(&declarations, input_file, "func_input", "config", 0); // the 0 means not self conjugate.
	free(declarations);
	
	FILE *OUT = safe_fopen_write(output_file);
	
	
	
	// copy the original input file as parsed above.
	FILE *IN = safe_fopen_read("func_input");
	fprintf(OUT,"INPUT\n\n");
	copyfile(IN,OUT);
	fclose(IN);
	
	
	
	// now put in the sphere equations
	
	int rand_index = rand(); // what if there are multiple sphere equations???  put in this random number to be safe.
	
	fprintf(OUT,"function sphere_%d;\n",rand_index);
    
	fprintf(OUT,"sphere_%d = ",rand_index);
	for (int jj=1; jj<W.num_variables(); jj++) { // start at 1 to omit the homogenizing variable
		fprintf(OUT," (%s-(", W.name(jj).c_str());
		mpf_out_str(OUT,10,0,sphere_center->coord[jj-1].r);
		fprintf(OUT,"+I*");
		mpf_out_str(OUT,10,0,sphere_center->coord[jj-1].i);
		fprintf(OUT,"))^2");
		
		if (jj!=W.num_variables()-1) {
			fprintf(OUT," + ");
		}
		
	}
    
	fprintf(OUT," - (");
	mpf_out_str(OUT,10,0,sphere_radius->r);
	fprintf(OUT,")^2 ");
    
	fprintf(OUT, ";\n\nEND;\n\n\n\n\n");
	
    for (unsigned int ii=0; ii<W.num_patches(); ii++) {
		std::stringstream linname;
		linname << "patch_" << ii;
		write_vector_as_constants(W.patch(ii), linname.str(), OUT);
		linname.clear();  linname.str("");
	}
    
    
	fclose(OUT);
	
	
	
}





