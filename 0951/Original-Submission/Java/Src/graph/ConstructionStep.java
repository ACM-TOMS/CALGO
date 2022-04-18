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

import java.util.ArrayList;
import java.util.HashSet;

import ccs.Debug;

public class ConstructionStep {
	Vertex stepVertex;
	ClusterPair clusterPair;
	Edge baseVertices;

	public Vertex v1() {
		return baseVertices.v1();
	}

	public Vertex v2() {
		return baseVertices.v2();
	}

	public Vertex v() {
		return stepVertex;
	}

	public Cluster c1() {
		return clusterPair.c1();
	}

	public Cluster c2() {
		return clusterPair.c2();
	}

	public ConstructionStep(Vertex stepVertex, ClusterPair clusterPair,
			Edge baseVertices) {
		super();
		this.stepVertex = stepVertex;
		this.clusterPair = clusterPair;
		this.baseVertices = baseVertices;
	}

	public ConstructionStep(Vertex stepVertex, Cluster c1, Cluster c2,
			Vertex v1, Vertex v2) {
		this(stepVertex, new ClusterPair(c1, c2), new Edge(v1, v2));
	}

	@Override
	public String toString() {
		return v() + " < " + " (" + v1() + "@" + c1() + ", " + v2() + "@"
				+ c2() + ")";
	}

	public static ArrayList<ConstructionStep> generateConstructionSequence(
			TDLinkage t) {
		ArrayList<ConstructionStep> constructionSequence = new ArrayList<ConstructionStep>();

		ArrayList<Vertex> constructedV = new ArrayList<Vertex>(), remainingV = new ArrayList<Vertex>();
		remainingV.addAll(t.getSharedVertices());
		HashSet<Cluster> constructedC = new HashSet<Cluster>(), remainingC = new HashSet<Cluster>();
		remainingC.addAll(t.getClusters());

		remainingV.remove(t.getBaseNonedge().v1());
		remainingV.remove(t.getBaseNonedge().v2());
		constructedV.add(t.getBaseNonedge().v1());
		constructedV.add(t.getBaseNonedge().v2());

		// TODO: First attempt: Brute Force
		while (remainingV.size() > 0) {
			for (int i = 0; i < constructedV.size(); ++i)
				for (int j = 0; j < i; ++j) {
					Vertex v1 = constructedV.get(i);
					Vertex v2 = constructedV.get(j);
					for (Cluster c1 : t.getClusters(v1))
						for (Cluster c2 : t.getClusters(v2)) {
							if (c1 == c2 || constructedC.contains(c1)
									|| constructedC.contains(c2))
								continue;

							Vertex v = Cluster.sharedVertex(c1, c2);
							if (v == null)
								continue;
							assert (remainingV.contains(v));

							ConstructionStep step;
							if (v1.index < v2.index)
								step = new ConstructionStep(v, c1, c2, v1, v2);
							else
								step = new ConstructionStep(v, c2, c1, v2, v1);
							constructionSequence.add(step);
							Edge e1 = new Edge(v1, v), e2 = new Edge(v2, v);

							constructedV.add(v);
							remainingV.remove(v);
							HashSet<Vertex> clusterVSet = c1
									.unionSharedVertices(c2);
							for (Vertex sv : clusterVSet) {
								if (!constructedV.contains(sv)) {
									constructedV.add(sv);
									remainingV.remove(sv);
								}
							}
							constructedC.add(c1);
							constructedC.add(c2);
							// if (DEBUG)
							// System.out.print(v + " constructed \t");
						}
				}
		}
		return constructionSequence;
	}

}