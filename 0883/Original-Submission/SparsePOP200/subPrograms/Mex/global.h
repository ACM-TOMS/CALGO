/* -------------------------------------------------------------

This file is a component of SparsePOP
Copyright (C) 2007 SparsePOP Project

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA

------------------------------------------------------------- */
#ifndef _global_
#define _global_
// Header files
#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <list>
#include <ctime>
#include <algorithm>
#include <numeric>
#include <set>
#include <sstream>
#include <functional>
#include <cassert>

// macros
#define EPS (1.0e-15)
#define MAX (1.0e8)
#define MIN (-1.0e8)
#define YES  1
#define NO  -1

//def. of typeCones
enum Cones{
	EQU=-1,		// equality constraint
	INE=1,		// inequality constriant or objective function	
	SOC=2,		// second-order cone constraint	
	SDP=3,		// positive semidefinite cone constraint	
};

#endif
