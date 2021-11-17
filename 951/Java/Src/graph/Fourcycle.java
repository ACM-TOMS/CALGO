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
import java.util.HashSet;

import ccs.Debug;

/** computes and stores the canonical edges of a TreeDecomp t
 */
public class Fourcycle {

	/** 
	 * @param t	the associated TreeDecomp graph 
	 */
	public Fourcycle(TDLinkage t) {
		this.t = t;
	}

	TDLinkage t;
	ArrayList<Edge> completeCayleyVector;
	// ArrayList<Pair<Cluster>> l; // adjacent cluster pairs in four-cycles;
	boolean isLow = false;

	/** computes the canonical base non-edges
	 * @return whether the associated TreeDecomp has low Cayley complexity
	 */
	public boolean init() {
		ArrayList<Pair<Cluster>> validBasePairs = new ArrayList<Pair<Cluster>>();
		completeCayleyVector = new ArrayList<Edge>();
		completeCayleyVector.add(t.getBaseNonedge());

		HashSet<Cluster> constructedClusters = new HashSet<Cluster>();

		for (int i = 0; i < t.getNumOfConstructStep(); ++i) {
			ConstructionStep s = t.getConstructionStep(i);

			if (i > 0) {
				// (1) U and W
				Vertex v1 = s.v1(), v2 = s.v2();
				HashSet<Cluster> cv1 = new HashSet<Cluster>();
				cv1.addAll(t.getClusters(v1));
				cv1.retainAll(constructedClusters);
				HashSet<Cluster> cv2 = new HashSet<Cluster>();
				cv2.addAll(t.getClusters(v2));
				cv2.retainAll(constructedClusters);

				// (2) find (c1,c2): adjacent base cluster pair
				boolean isLowStep = false;
				Pair<Cluster> pair = null;
				outer: for (Cluster c1 : cv1) {
					for (Cluster c2 : cv2) {
						assert (!c1.equals(c2));
						pair = new Pair<Cluster>(c1, c2);
						if (validBasePairs.contains(pair)) {
							isLowStep = true;
							break outer; // TODO: pick the pair constructed
											// earliest
						}
					}
				}
				if (!isLowStep) {
					Debug.warnMsg("Step " + s.v() + " makes t not low");
					return false;
				}
				assert (pair != null);
				Vertex v = Cluster.sharedVertex(pair.getFirst(),
						pair.getSecond());
				assert (v != null);
				completeCayleyVector.add(new Edge(v, s.v()));

				// (3)
				for (Cluster c1 : cv1) {
					for (Cluster c2 : cv2) {
						assert (!c1.equals(c2));
						if (Cluster.sharedVertex(c1, c2) != null) {
							validBasePairs.add(new Pair<Cluster>(s.c1(), c1));
							validBasePairs.add(new Pair<Cluster>(s.c2(), c2));
						}
					}
				}
			}
			validBasePairs.add(s.clusterPair);
			constructedClusters.add(s.c1());
			constructedClusters.add(s.c2());
		}
		isLow = true;
		return true;

	}

	public AbstractList<Edge> getCompleteCayleyVector() {
		assert (isLow);
		return completeCayleyVector;
	}

}
