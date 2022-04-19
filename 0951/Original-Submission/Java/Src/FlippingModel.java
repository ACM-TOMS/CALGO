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
import java.util.Collection;

import ccs.graph.Realization;
import ccs.graph.RealizationType;
import ccs.graph.TDLinkage;
import ccs.graph.Vertex;

public class FlippingModel {
	private FlippingModel() {
		toFlip = new ArrayList<Vertex>();
	}

	private static FlippingModel me = new FlippingModel();

	public static FlippingModel getInstance() {
		return me;
	}

	private TDLinkageModel tdModel = TDLinkageModel.getInstance();

	/**
	 * Vertices to flip orientations
	 */
	private ArrayList<Vertex> toFlip;

	public void clear() {
		toFlip.clear();
	}

	public Collection<Vertex> getVerteicsToFlip() {
		return toFlip;
	}

	public void doFlip() {
		/*
		 * assert (getStatus() == Status.generated || getStatus() ==
		 * Status.picking);
		 */
		RealizationType t = tdModel.getForwardSolutionType();
		for (Vertex v : toFlip)
			t.flipOrientation(tdModel.isStepVertex(v));

		Realization g = tdModel.tryRealize(t);
		if (g != null) {
			tdModel.setPoints(g);
		} else {
			// TODO: say not possible
		}
		clear();
	}

	/**
	 * Add v to flip list if it is not already there. If v is in list, remove v
	 */
	public void flipVertex(Vertex v) {
		if (toFlip.contains(v)) {
			toFlip.remove(v);
		} else {
			int index = tdModel.isStepVertex(v);
			if (index >= 0) {
				toFlip.add(v);
			}
		}
	}

}
