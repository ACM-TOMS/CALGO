/*
 This file is part of CayMos. 

 CayMos is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 CayMos is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package ccs.graph;

import java.io.Serializable;

public class Vertex implements Serializable{

	// now vertex does not have position
	// index should only be used as name?
	// position retrieved by graph, using vertex or v.index as key?

	public int index;

	public Vertex(int index) {
		this.index = index;
	}

	public String toString() {
		return "v" + (index);
	}
}
