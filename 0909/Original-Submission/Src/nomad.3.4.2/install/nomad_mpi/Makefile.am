CXX = mpic++
CXXFLAGS = -O3
INCLUDES = -I$(top_srcdir)/nomad_mpi
METASOURCES = AUTO
AM_CXXFLAGS = -DGCC_X -DUSE_MPI
noinst_HEADERS = Barrier.hpp Cache.hpp Cache_File_Point.hpp Cache_Point.hpp \
	Cache_Search.hpp Clock.hpp Direction.hpp Directions.hpp Display.hpp Double.hpp \
	Eval_Point.hpp Evaluator.hpp Evaluator_Control.hpp Exception.hpp Extended_Poll.hpp \
	Filter_Point.hpp LH_Search.hpp L_Curve.hpp Mads.hpp Mesh.hpp Multi_Obj_Evaluator.hpp \
	Parameter_Entries.hpp Parameter_Entry.hpp Parameters.hpp Pareto_Front.hpp Pareto_Point.hpp \
	Phase_One_Evaluator.hpp Phase_One_Search.hpp Point.hpp Priority_Eval_Point.hpp \
	Random_Pickup.hpp Search.hpp Set_Element.hpp Signature.hpp Signature_Element.hpp \
	Slave.hpp Speculative_Search.hpp Stats.hpp Uncopyable.hpp VNS_Search.hpp \
	Variable_Group.hpp defines.hpp nomad.hpp utils.hpp
lib_LIBRARIES = libnomad.MPI.a
libnomad_MPI_a_SOURCES = Barrier.cpp Cache.cpp Cache_File_Point.cpp \
	Cache_Point.cpp Cache_Search.cpp Clock.cpp Direction.cpp Directions.cpp Display.cpp \
	Double.cpp Eval_Point.cpp Evaluator.cpp Evaluator_Control.cpp Exception.cpp \
	Extended_Poll.cpp LH_Search.cpp L_Curve.cpp Mads.cpp Mesh.cpp Multi_Obj_Evaluator.cpp \
	Parameter_Entries.cpp Parameter_Entry.cpp Parameters.cpp Pareto_Front.cpp Pareto_Point.cpp \
	Phase_One_Evaluator.cpp Phase_One_Search.cpp Point.cpp Priority_Eval_Point.cpp \
	Random_Pickup.cpp Signature.cpp Slave.cpp Speculative_Search.cpp Stats.cpp VNS_Search.cpp \
	Variable_Group.cpp utils.cpp
bin_PROGRAMS = nomad.MPI
nomad_MPI_SOURCES = nomad.cpp
nomad_MPI_LDADD = $(top_builddir)/nomad_mpi/libnomad.MPI.a
