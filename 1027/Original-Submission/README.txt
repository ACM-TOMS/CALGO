1. The user manual may be made by issuing the command

   make

   This will generate a file called nomad4.pdf. 

2. In order to install the nomad library, etc under Linux Ubuntu 22.04, I had
   to install petsc-dev in order to satisfy cmake. It required the directory

   /usr/lib/x86_64-linux-gnu/openmpi/include/

   which was missing (/usr/lib/x86_64-linux-gnu/openmpi/lib did exist)

Tim Hopkins
June 30th 2022
