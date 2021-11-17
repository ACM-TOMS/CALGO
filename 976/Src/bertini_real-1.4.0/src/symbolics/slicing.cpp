#include "symbolics/slicing.hpp"




void create_sliced_system(boost::filesystem::path input_file, boost::filesystem::path output_file,
                          vec_mp * linears, int num_to_add,
                          const WitnessSet & W)
{
#ifdef functionentry_output
	std::cout << "create_sliced_system" << std::endl;
#endif
	
	
	
	
	if (W.num_var_names()==0) {
		std::cout << "trying to create a sliced system, but witness set does not have the variable names." << std::endl;
		deliberate_segfault();
	}
	int *declarations = NULL;
	
	partition_parse(&declarations, input_file, "func_input", "config", 0); // the 0 means not self conjugate.
	free(declarations);
	
	FILE *OUT = safe_fopen_write(output_file);
	FILE *IN = safe_fopen_read("func_input");
	
	
	
	
	fprintf(OUT,"INPUT\n\n");
	copyfile(IN,OUT);
	fclose(IN);
	
	std::vector< int > indicators;
	for (int ii=0; ii<num_to_add; ii++) {
		indicators.push_back(rand());
	}
	
	for (int ii=0; ii<num_to_add; ii++) {
		std::stringstream linname;
		linname << "supp_lin_" << indicators[ii];
		write_vector_as_constants(linears[ii], linname.str(), OUT);
		
		linname.clear();  linname.str("");
	}
	
	for (int ii=0; ii<num_to_add; ii++) {
		fprintf(OUT,"function supp_lin_%d;\n",indicators[ii]);
	}
	for (int ii=0; ii<num_to_add; ii++) {
		fprintf(OUT,"supp_lin_%d = supp_lin_%d_1",indicators[ii],indicators[ii]);
		for (int jj=1; jj<W.num_variables(); jj++) {
			fprintf(OUT," + %s*supp_lin_%d_%d",W.name(jj).c_str(), indicators[ii],jj+1);
		}
		fprintf(OUT, ";\n\n");
	}
	fprintf(OUT,"END;\n\n\n\n");
	
    
    for (unsigned int ii=0; ii<W.num_patches(); ii++) {
		std::stringstream linname;
		linname << "patch_" << ii;
		write_vector_as_constants(W.patch(ii), linname.str(), OUT);
		linname.clear();  linname.str("");
	}
    
	fclose(OUT);
	
	
	
	
}










