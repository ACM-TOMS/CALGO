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

import java.util.AbstractList;
import java.util.ArrayList;

public class ComponentCompleteCayleySampler implements NodeSampler<ArrayList<Double>> {
	int dimension;
	int samplesPerNode;
	TDLinkage t;

	public ComponentCompleteCayleySampler(TDLinkage t, int spn) {
		assert (spn > 0);
		this.t = t;
		dimension = t.getNumOfConstructStep();
		samplesPerNode = spn;
	}

	// ???
	public int getSamplePointDimension() {
		return dimension;
	}

	public AbstractList<SamplePoint<ArrayList<Double>>> sample(OrientedInterval n) {
		ArrayList<SamplePoint<ArrayList<Double>>> list = new ArrayList<SamplePoint<ArrayList<Double>>>();

		for (int i = 0; i < samplesPerNode; ++i) {
			ArrayList<Double> l = new ArrayList<Double>();

			double percentage = (double) i / samplesPerNode;
			double cayley = n.sampleCayleyAt(percentage);
			Realization g = t.tryRealize(cayley, n.getSolutionType());
			assert (g != null);
			for (Edge e : t.completeCayleyVector()) {
				double le = g.length(e);
				l.add(le);
			}
			list.add(new SamplePoint<ArrayList<Double>>(l));
		}

		return list;
	}

}