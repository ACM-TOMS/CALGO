Example1
--------
This folder contains a makefile for the two runs associated with
"Example1" as described in the user's manual in pages 4-7. 

Type 
   make example1
in the Examples directory to run these tests.

Example2
--------
This folder contains a makefile for the run associated with
"Example2" as described in the user's manual in pages 7-8. 

Type 
   make example2
in the Examples directory to run this test.

Fun
---
This folder contains a makefile for the 13 runs associated with
"Fun" as described in the user's manual in pages 12-13. 

Type 
   make fun
in the Examples directory to run these tests.

Fun100
------
This folder contains a makefile for the five runs associated with
"Example1" as described in the user's manual in pages 10-11. 

Type 
   make fun100
in the Examples directory to run these tests.

Fun1000
-------
This folder contains a makefile for the five runs associated with
"Example1" as described in the user's manual in pages 4-7. 

Type 
   make fun1000
in the Examples directory to run these tests.

To run all the tests
--------------------
Type 
  make
in the Examples directory to run all the above tests

Expected results
----------------
The expected results files are all stored in the Results
directory in each of the subdirectories.

Cleaning the directories
------------------------
To remove all files associated with each subdirectory, use
  cd <subdir> && make clean
from the Examples directory where <subdir>is one of Example1,
Example2, Fun, Fun100, Fun1000

To clean all the directories use
  make clean
from the Examples directory.

