Do not delete either the folder •00_Dependencies• or its contents. 
It contains additional files (like vertex and fragment shaders, 
model files, and -- after successful installation -- binary 
libraries and include files) which are used by:

- the files 
	../01_Library/cagd/cagd.pro, 
	../01_Library/Makefile.Windows.Debug.{x86|x64}, and
	../01_Library/Makefile.Windows.Release.{x86|x64}
  under Windows platforms;

- all test applications provided inside the folders 
  •../03_Examples/Example_0x•, where x = 1,2,3,4,5 (on all platforms).
  

The subfolders •Include• and •Lib• are important only on Windows 
platforms (for more details see subsections 1.c and 2.1.c of the 
readme-file •../01_Library/readme.txt•).
