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

package ccs;

import java.util.ArrayList;

import ccs.graph.Vertex;

public class TraceModel {
	private TraceModel() {
		tracingVertices = new ArrayList<Vertex>();
		traces = new ArrayList<ArrayList<Trace>>();
	}

	private static TraceModel me = new TraceModel();

	public static TraceModel getInstance() {
		return me;
	}

	private ArrayList<Vertex> tracingVertices;

	private ArrayList<ArrayList<Trace>> traces;

	public ArrayList<Vertex> getTracingVertices() {
		return tracingVertices;
	}

	public ArrayList<ArrayList<Trace>> getTraces() {
		return traces;
	}

	public void clear() {
		tracingVertices.clear();
		traces.clear();
	}
}
