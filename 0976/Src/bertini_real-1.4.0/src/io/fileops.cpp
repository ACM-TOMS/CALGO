#include "io/fileops.hpp"


void WaitOnGeneratedFile(const std::string & name)
{
	
	FILE *IN = fopen(name.c_str(),"r");
	int num_waits = 0;
	while (IN==NULL){
		sleep(1); 
		if (num_waits==0)
			printf("waiting for %s\n",name.c_str());
		else if (num_waits%15==0)
			printf("waiting for file '%s' (it's been %d seconds)\n",name.c_str(),num_waits);
		IN = fopen(name.c_str(),"r");
		num_waits++;
	}
	fclose(IN);
}



int partition_parse(int **declarations,
										boost::filesystem::path input_filename,
										boost::filesystem::path functions_filename,
										boost::filesystem::path config_filename,
										int not_sc_flag)
{
	
	FILE *IN = safe_fopen_read(input_filename);
	
	int retval = partitionParse(declarations, IN,
															const_cast<char *> (functions_filename.c_str()),
															const_cast<char *> (config_filename.c_str()),
															not_sc_flag); // the 0 means not self conjugate.
	fclose(IN);
	
	return retval;
}



void rename_bertini_files_dotbak(){
	
	rename("arr.out","arr.out.bak");
	rename("num.out","num.out.bak");
	rename("deg.out","deg.out.bak");
	rename("config","config.bak");
	rename("preproc_data","preproc_data.bak");
}

void restore_bertini_files_dotbak(){
	rename("config.bak","config");
	rename("arr.out.bak","arr.out");
	rename("num.out.bak","num.out");
	rename("deg.out.bak","deg.out");
	rename("preproc_data.bak","preproc_data");
	
}

//will remove all files which do not start with a . from directory.
void purge_previous_directory(char *directoryName)
{
	DIR *dp;
	struct dirent *ep;
	char tempname[1000];
	
	
	dp = opendir (directoryName);
	if (dp != NULL)
	{
		while ( (ep = readdir (dp)) )
			if (ep->d_name[0] != '.')
	    {
				sprintf(tempname,"%s/%s",directoryName, ep->d_name);
				remove(tempname);
			}
		
		
		(void) closedir (dp);
	}
	return;
}



FILE *safe_fopen_read(boost::filesystem::path filename)
{
	
	FILE* IN;
	
	if (boost::filesystem::is_directory(filename)){
		throw std::runtime_error("trying to open directory " + filename.string() + " as a file!");
	}
	
	if (!boost::filesystem::exists(filename)){
		throw std::runtime_error("unable to find specified file to read: " + filename.string());
	}
	
	
	IN = fopen(filename.c_str(),"r");

	if (IN == NULL) {
		throw std::runtime_error("failed to open file, pointer is NULL: " + filename.string());
	}
	
	return IN;
	
}


FILE *safe_fopen_write(boost::filesystem::path filename)
{
	
	FILE* OUT;
	
	
	OUT = fopen(filename.c_str(),"w");
	
	
	if (OUT == NULL) {
		throw std::runtime_error("unable to open specified file to write: " + filename.string());
	}
	
	return OUT;
	
}


FILE *safe_fopen_append(boost::filesystem::path filename)
{
	FILE* OUT;

	struct stat stat_p;
	stat (filename.c_str(), &stat_p);
	
	if (S_ISDIR(stat_p.st_mode)){
		throw std::runtime_error("trying to open directory as a file: " + filename.string());
	}
	
	
	OUT = fopen(filename.c_str(),"a");
	
	
	if (OUT == NULL) {
		throw std::runtime_error("unable to open specified file to write: " + filename.string());
	}
	
	return OUT;
}


void copyfile(boost::filesystem::path input_file, boost::filesystem::path OUTfile)
{
	
	FILE *IN  = safe_fopen_read(input_file);
	FILE *OUT = safe_fopen_write(OUTfile);
	
	copyfile(IN,OUT);
	
	fclose(IN);
	fclose(OUT);
}



void copyfile(FILE *IN,FILE *OUT)
{
	char ch;
	while ((ch = fgetc(IN)) != EOF)
		fprintf(OUT, "%c", ch);

}






void deliberate_segfault()
{
	printf("the following segfault is deliberate\n");
	int *faulty = NULL;
	faulty[-10] = faulty[10]+faulty[0];
	
}


void br_exit(int errorCode)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: exits Bertini_real - either standard or using MPI           *
 \***************************************************************/
{	
	std::cout << "bertini_real quitting\a" << std::endl;
	
#ifdef debug_compile
	deliberate_segfault();
#endif
	
	MPI_Abort(MPI_COMM_WORLD, errorCode);
}


















bool parseInteger( std::string const& text, int& results )
{
    std::istringstream parser( text );
	parser >> results;
	return !(parser.fail());// >> std::ws && parser.peek() == EOF;
}





int getInteger()
{
	
	std::cin.ignore( std::numeric_limits<std::streampos>::max(), '\n' ) ;
	std::cin.clear() ; // to be safe, clear error flags which may be set
	
	int results;
	std::string line;
	while (1){
		while ( ! std::getline( std::cin, line )  )
		{
			
			std::cin.ignore( std::numeric_limits<std::streampos>::max(), '\n' ) ;
			std::cin.clear() ; // to be safe, clear error flags which may be set
			
		}
		
		if (! parseInteger( line, results )) {
			std::cout << "Only 'numeric' value(s) allowed:" << std::endl;
		}
		else{
			break;
		}
		
	}
	return results;
}


//  use this function to get input from the user and ensure it is an integer (actually fails to detect non-integer numeric inputs
int get_int_choice(std::string display_string,int min_value,int max_value){
	
	std::cout << display_string;
	
	int userinput = min_value - 1;
	
	std::string tmpstr;
	
	while (1)
	{
		userinput = getInteger();
		if (userinput <= max_value && userinput >= min_value){
			break;
		}
		else{
			std::cout << "value out of bounds." << std::endl;
		}
    }
	return userinput;
}



int get_int_choice(std::string display_string,const std::set<int> & valid_values)
{
	if (valid_values.size()==0) {
		std::cout << "trying to get int choice from empty set of valid values..." << std::endl;
		return 0;
	}
	
	
	std::cout << display_string;
	
	int userinput = *(valid_values.begin()) - 1;
	

	while (1)
	{
		userinput = getInteger();
		if (valid_values.find(userinput)==valid_values.end()){
			std::cout << "value out of bounds." << std::endl;
		}
		else{
			break;
		}
    }
	return userinput;
}




















