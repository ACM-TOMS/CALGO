/*
 This file is part of EASAL. 

 EASAL is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 EASAL is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
/*
 * ThreadShare.h
 *
 *  Created on: May 20, 2016
 *      Author: rprabhu 
 */
#ifndef THREADSHARE_H_
#define THREADSHARE_H_

#include "AtlasBuilder.h"
#include "Atlas.h"
#include "SaveLoader.h"
namespace ThreadShare {
extern SaveLoader *save_loader;
extern AtlasBuilder *atlas_builder;
extern Atlas *atlas_view;
}

#endif
