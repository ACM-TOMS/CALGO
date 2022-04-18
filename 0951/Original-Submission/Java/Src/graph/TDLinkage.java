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
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;

import ccs.Debug;

public class TDLinkage {
	final static double ACCURACY = 0.0001;

	// adjacent-list, indexed by v
	private LinkedList<Cluster> clusters;

	public Collection<Cluster> getClusters() {
		return clusters;
	}

	private ArrayList<Vertex> sharedVertices;

	public Collection<Vertex> getSharedVertices() {
		return sharedVertices;
	}

	private HashMap<Vertex, LinkedList<Cluster>> vClusters;

	// Usually Vertex v0, v1;
	private Edge baseNonEdge;

	private ArrayList<ConstructionStep> constructionSequence;

	// 0:not initialized, 1:true, -1:false
	private int isTriangleFree = 0;
	private int is1Path = 0;
	private Fourcycle fourcycleInfo;

	public Realization graph;

	public CayleyConfigSpace cayleyConfigSpace;

	// Initialize: clusters, vClusters, sharedVertices <- decomposition of g
	// does not do base non-edge.
	// may not be 1-dof t-d
	public TDLinkage(Realization g) {
		this.graph = g;
		vClusters = new HashMap<Vertex, LinkedList<Cluster>>();
		for (int i = 0; i < g.getSize(); ++i) {
			Vertex v = g.getVertex(i);
			vClusters.put(v, new LinkedList<Cluster>());
		}
		sharedVertices = new ArrayList<Vertex>();
		this.decompose();
	}

	private void decompose() {
		this.clusters = Util.getDecompose(graph);

		for (Cluster c : clusters) {
			for (Vertex v : c.getVertices()) {
				LinkedList<Cluster> vclusters = getClusters(v);
				assert (!vclusters.contains(c));
				vclusters.add(c);
				if (vclusters.size() > 1 && !sharedVertices.contains(v))
					sharedVertices.add(v);
			}
		}
		for (Vertex v : sharedVertices) {
			for (Cluster c : getClusters(v))
				c.addSharedVertex(v);
		}
	}

	public boolean isTreeDecomposable() {
		return clusters.size() == 1;
	}

	// set the base non-edge
	public boolean is1DofTreeDecomposable() {
		if (baseNonEdge != null)
			return true;
		if (isTreeDecomposable()) {
			Debug.msg("Not 1-dof: is tree-decomposable");
			return false;
		}

		// TODO: remove last level v?
		// brute force: for all possible base non-edge
		for (int i = 0; i < graph.getSize(); ++i) {
			for (int j = 0; j < i; ++j) {
				Vertex v1 = graph.getVertex(i), v2 = graph.getVertex(j);
				if (!graph.isAdjacent(v1, v2)) {
					Debug.msg("Trying base non-edge: (" + v1 + ", " + v2 + ")");
					graph.addEdge(v1, v2);
					if (Util.getDecompose(graph).size() == 1) {
						graph.removeEdge(v1, v2);
						this.baseNonEdge = new Edge(v2, v1); // index ordered
						assert (baseNonEdge.v1().index < baseNonEdge.v2().index);
						Debug.msg("new base nonedge:" + baseNonEdge);
						return true;
					}
					graph.removeEdge(v1, v2);
				}
			}
		}
		return false;
	}

	// ======================= Base Non Edge =====================

	public double getBaseNoneedgeLength() {
		return graph.length(baseNonEdge);
	}

	public Edge getBaseNonedge() {
		return baseNonEdge;
	}

	public boolean setBaseNonedge(Edge base) {
		if (graph.isAdjacent(base))
			return false;
		graph.addEdge(base);
		if (Util.getDecompose(graph).size() != 1) {
			graph.removeEdge(base);
			return false;
		} else {
			graph.removeEdge(base);
			if (base.v1().index > base.v2().index)
				base.swap();
			this.baseNonEdge = base;
			this.generateConstructionSequence();
			this.is1Path = 0;
			Debug.msg("base non-edge reset to:" + baseNonEdge);
			return true;
		}
	}

	// ======================= CLUSTER =========================

	public LinkedList<Cluster> getClusters(Vertex v) {
		return vClusters.get(v);
	}

	public Cluster getClusterBetween(Vertex v1, Vertex v2) {
		for (Cluster c1 : getClusters(v1))
			for (Cluster c2 : getClusters(v2))
				if (c1 == c2)
					return c1;
		return null;
	}

	public int getCDeg(Vertex v) {
		int d = getClusters(v).size();
		return d;
	}

	public int indexOf(Cluster c) {
		return clusters.indexOf(c);
	}

	// ========================= LEVEL, PATH ================

	public boolean isLastLevel(Vertex v) {
		int d = getCDeg(v);
		if (d != 2)
			return false;
		Cluster c1 = getClusters(v).getFirst(), c2 = getClusters(v).getLast();
		if (c1.getNumOfSharedVertices() > 2 || c2.getNumOfSharedVertices() > 2)
			return false;
		return true;
	}

	public boolean is1Path() {
		if (is1Path != 0)
			return is1Path == 1;

		int count = 0;
		is1Path = -1;
		for (Vertex v : sharedVertices) {
			if (v == baseNonEdge.v1() || v == baseNonEdge.v2())
				continue;
			if (isLastLevel(v)) {
				if (count == 0)
					++count;
				else
					return false;
			}
		}
		is1Path = 1;
		return true;
	}

	public boolean isTriangleFree() {
		if (isTriangleFree != 0)
			return isTriangleFree == 1;

		isTriangleFree = -1;
		for (Cluster c : clusters)
			if (c.getSize() > 2) {
				assert (c.getSize() >= 3);
				if (Realization.DEBUG) {
					System.out.println("triangle in cluster: " + c);
				}
				return false;
			}
		isTriangleFree = 1;
		return true;
	}

	// ====================== LOW =============================
	public boolean isLow() {
		boolean b = is1DofTreeDecomposable();
		assert (b);
		if (!this.constructionSequenceGenerated()) {
			generateConstructionSequence();
			this.normalizeBaseNonEdge();
		}

		if (fourcycleInfo == null) {
			fourcycleInfo = new Fourcycle(this);
			return fourcycleInfo.init();
		} else {
			return fourcycleInfo.isLow;
		}

	}

	public AbstractList<Edge> completeCayleyVector() {
		boolean b = isLow();
		assert (b);
		return fourcycleInfo.getCompleteCayleyVector();
	}

	// =================== CONSTRUCTION SEQUENCE ==================
	public int getNumOfConstructStep() {
		if (!constructionSequenceGenerated()) {
			this.generateConstructionSequence();
			this.normalizeBaseNonEdge();
		}
		return constructionSequence.size();
	}

	public ConstructionStep getConstructionStep(int i) {
		assert (i < getNumOfConstructStep());
		return constructionSequence.get(i);
	}

	public boolean constructionSequenceGenerated() {
		return this.constructionSequence != null;
	}

	// sort shared vertices & generate cluster pairs
	// TODO: !!! sometimes I don't wanna normalize the base non edge !!!!
	public void generateConstructionSequence() {
		// Debug.msg("\n ----------\n graph:" + this);
		boolean b = is1DofTreeDecomposable();
		assert (b);
		// Debug.warnMsg(this.baseNonEdge);

		constructionSequence = ConstructionStep
				.generateConstructionSequence(this);

		Debug.msg(constructionSequence.toString());
		/*
		 * this.normalizeBaseNonEdge(); Debug.msg("after normalization:" +
		 * this.constructionSequence);
		 */
	}

	public void normalizeBaseNonEdge() {
		if (this.getNumOfConstructStep() < 2)
			return;
		ConstructionStep s1 = constructionSequence.get(0);
		ConstructionStep s2 = constructionSequence.get(1);
		assert (s1.baseVertices.equals(baseNonEdge));
		if (!baseNonEdge.equals(s2.baseVertices)) {
			if (Realization.DEBUG)
				System.out.println(baseNonEdge + "->" + s2.baseVertices);

			baseNonEdge = s2.baseVertices;

			assert (baseNonEdge.v1().index < baseNonEdge.v2().index);
			if (!s1.c1().containsSharedVertex(baseNonEdge.v1())) {
				assert (s1.c1().containsSharedVertex(baseNonEdge.v2()));
				assert (s1.c2().containsSharedVertex(baseNonEdge.v1()));
				s1.clusterPair.swap();
			} else {
				assert (s1.c2().containsSharedVertex(baseNonEdge.v2()));
			}
			s1.baseVertices = baseNonEdge;
		}
		Debug.msg("after normalization:" + this.constructionSequence);
	}

	// ===================== orientations ======================
	private int orientationOf(ConstructionStep s, Realization g) {
		// If this cross product is positive, then v1->v is clockwise from
		// v1->v2; if negative, it is counterclockwise.
		Vertex v1 = s.v1(), v2 = s.v2();
		assert (v1.index < v2.index); // constrain the order of base points
		// Debug.warnMsg("g:"+g+",v1:"+v1+",v2:"+v2+",s:"+s);
		return g.orientationOf(v1, v2, s.v());
	}

	// for a subgraph?
	public RealizationType getForwardSolutionType(Realization g, int numOfSteps) {
		if (!this.constructionSequenceGenerated()) {
			this.generateConstructionSequence();
			this.normalizeBaseNonEdge();
		}
		RealizationType forward = new RealizationType(getNumOfConstructStep());
		for (int i = 0; i < numOfSteps; ++i) {
			ConstructionStep s = constructionSequence.get(i);
			forward.setOrientation(i, orientationOf(s, g));
		}
		return forward;
	}

	public RealizationType getForwardSolutionType(Realization g) {
		return getForwardSolutionType(g, getNumOfConstructStep());
	}

	public RealizationType getForwardSolutionType() {
		return getForwardSolutionType(graph);
	}

	// ======================= extreme ==============================

	// Does not actually connect the extreme edge
	// Returns a t-d w/ extreme non-edge as base
	public TDLinkage getExtremeGraph(int index) {
		if (!constructionSequenceGenerated())
			generateConstructionSequence();
		assert (index < getNumOfConstructStep());
		ArrayList<Vertex> vset = new ArrayList<Vertex>();
		vset.add(baseNonEdge.v1());
		vset.add(baseNonEdge.v2());
		for (int i = 0; i < index; ++i) {
			ConstructionStep s = constructionSequence.get(i);
			// System.out.println("Step "+s);
			HashSet<Vertex> clusterV = s.c1().unionVertices(s.c2());
			clusterV.remove(s.v1());
			clusterV.remove(s.v2());
			// clusterV.remove(s.v());
			for (Vertex v : clusterV)
				vset.add(v);
		}
		Debug.msg("extreme's v: " + vset);
		Realization eg = graph.inducedSubgaph(vset);
		TDLinkage t = new TDLinkage(eg);
		ConstructionStep s = constructionSequence.get(index);
		// System.out.println(s.getV1()+","+s.getV2());
		Edge base = new Edge(s.v1(), s.v2());

		// TODO:
		t.generateConstructionSequence();
		// System.out.println("set base to: "+base);
		// assert (t.setBaseNonedge(base)); //!!!!
		// Debug.warnMsg(index + "th extreme's base is " + t.getBaseNonedge()
		// + ", should be " + base);
		if (!t.setBaseNonedge(base))
			t.baseNonEdge = null;
		// Debug.warnMsg(index + "th extreme's base is now " +
		// t.getBaseNonedge());

		// System.out.println(t.baseNonEdge);
		return t;
	}

	/**
	 * @return the achievable interval for the base pair of vertices of the
	 *         given construction step
	 */
	Interval extremeInterval(ConstructionStep step) {
		assert (step != null);
		double l1, l2;
		l1 = graph.distance(step.v1(), step.v());
		l2 = graph.distance(step.v2(), step.v());
		double lower = Math.abs(l1 - l2);
		double upper = l1 + l2;
		return new Interval(lower, upper);
	}

	// =========================== CCS ===============================
	public Double barLength(Edge e) {
		// System.out.println("get edge length:"+e);
		return graph.initLength(e);
	}

	public Double barLength(Vertex v1, Vertex v2) {
		return graph.initDistance(v1, v2);
	}

	public ArrayList<Realization> solve(double length) {
		Point2D p1 = graph.getPoint(baseNonEdge.v1());
		Point2D p2 = new Point2D(p1.x() + length, p1.y());
		return solve(new EdgePos(p1, p2));
	}

	/**
	 * Given lf, solve the possible realizations w/ any solution type. Does not
	 * modify "this".
	 * 
	 * @param newBasePos
	 *            new position of base non-edge
	 * @return list of all possible realizations
	 */
	public ArrayList<Realization> solve(EdgePos newBasePos) {
		ArrayList<Realization> solutions = new ArrayList<Realization>();
		if (!constructionSequenceGenerated()) {
			generateConstructionSequence();
			this.normalizeBaseNonEdge();
		}

		Realization g = graph.clone();
		g.setPoint(baseNonEdge.v1(), newBasePos.p1());
		g.setPoint(baseNonEdge.v2(), newBasePos.p2());

		this.solve(solutions, 0, g);

		return solutions;
	}

	// TODO: CLEAR UP ======================================

	/**
	 * Called by solve(EdgePos). Recursively solve curGraph, does not restrict
	 * solution type
	 * 
	 * @param solutions
	 *            storing the resulted realizations
	 * @param stepIndex
	 *            the step to start solving
	 * @param curGraph
	 * @return whether results exist
	 */
	private boolean solve(ArrayList<Realization> solutions, int stepIndex,
			Realization curGraph) {

		if (stepIndex >= this.getNumOfConstructStep()) {
			// successfully obtained one solution!
			solutions.add(curGraph.clone());
			return true;
		}

		ConstructionStep s = constructionSequence.get(stepIndex);
		Vertex v1 = s.v1(), v2 = s.v2(), v = s.v();
		assert (v1.index < v2.index);

		// compute (from old position) e1, e2
		double r1 = barLength(v1, v);// graph.distance(v1, v);
		double r2 = barLength(v2, v);// graph.distance(v2, v);

		Point2D p1 = curGraph.getPoint(v1), p2 = curGraph.getPoint(v2);

		Point2D pv = Util.solve(p1, p2, r1, r2, 1);
		if (pv == null)
			return false;
		curGraph.setPoint(v, pv);

		// transform cluster
		HashSet<Vertex> vC1 = s.c1().getVertices();
		vC1.remove(v);
		vC1.remove(v1);
		if (!vC1.isEmpty())
			curGraph.transformVertices(vC1, s.c1(), new Edge(v1, v),// graph.getPoints(v1,
																	// v),
					curGraph.getPoints(v1, v));

		HashSet<Vertex> vC2 = s.c2().getVertices();
		vC2.remove(v);
		vC2.remove(v2);
		if (!vC2.isEmpty())
			curGraph.transformVertices(vC2, s.c2(), new Edge(v2, v), // graph.getPoints(v2,
																		// v),
					curGraph.getPoints(v2, v));

		boolean b = this.solve(solutions, stepIndex + 1, curGraph);

		Point2D pv2 = Util.solve(p1, p2, r1, r2, -1);

		if (Math.abs(pv2.x() - pv.x()) > ACCURACY
				|| Math.abs(pv2.y() - pv.y()) > ACCURACY) {
			// if (Math.abs(Xd - Xc) > ACCURACY || Math.abs(Yd - Yc) > ACCURACY)
			// {

			curGraph.setPoint(v, pv2);
			// curGraph.setPoint(v, Xd, Yd);

			// if (!vC1.isEmpty())
			// curGraph.transformVertices(vC1, graph.getPoints(v1, v),
			// curGraph.getPoints(v1, v));
			if (!vC1.isEmpty())
				curGraph.transformVertices(vC1, s.c1(), new Edge(v1, v),
						curGraph.getPoints(v1, v));
			// if (!vC2.isEmpty())
			// curGraph.transformVertices(vC2, graph.getPoints(v2, v),
			// curGraph.getPoints(v2, v));
			if (!vC2.isEmpty())
				curGraph.transformVertices(vC2, s.c2(), new Edge(v2, v),
						curGraph.getPoints(v2, v));

			boolean b2 = this.solve(solutions, stepIndex + 1, curGraph);
			b = b || b2;
		}

		// if (Graph.DEBUG)
		// System.out.println("r1,r2,r3: " + Vertex.distance(u1, u) + ","
		// + Vertex.distance(u2, u) + "," + Vertex.distance(u1, u2));

		return b;
	}

	/**
	 * Given an lf, solve for realization with a specific solution type
	 * 
	 * @param length
	 *            length of new base non-edge
	 * @param forward
	 *            specified solution type
	 * @return the realization, or null if none exists
	 */
	public Realization solve(double length, RealizationType forward) {
		Point2D p1 = graph.getPoint(baseNonEdge.v1());
		Point2D p2 = new Point2D(p1.x() + length, p1.y());
		return solve(new EdgePos(p1, p2), forward);
	}

	// solve w/ sol type
	private Realization solve(EdgePos newBasePos, RealizationType forward) {
		if (!constructionSequenceGenerated()) {
			generateConstructionSequence();
			this.normalizeBaseNonEdge();
		}
		Realization g = graph.clone();

		g.setPoint(baseNonEdge.v1(), newBasePos.p1());
		g.setPoint(baseNonEdge.v2(), newBasePos.p2());

		for (int i = 0; i < this.getNumOfConstructStep(); ++i) {
			ConstructionStep s = this.getConstructionStep(i);
			Vertex v1 = s.v1(), v2 = s.v2(), v = s.v();
			// if (Graph.DEBUG)
			// Debug.msg("solving " + v + " from construction step " + s);
			assert (v1.index < v2.index);

			// compute (from old position) e1, e2
			double r1 = barLength(v1, v);// graph.distance(v1, v);
			double r2 = barLength(v2, v);// graph.distance(v2, v);

			Point2D p1 = g.getPoint(v1), p2 = g.getPoint(v2);

			int orient = forward.getOrientation(i);

			Point2D pv = Util.solve(p1, p2, r1, r2, orient);
			if (pv == null)
				return null;
			g.setPoint(v, pv);

			// transform cluster
			HashSet<Vertex> vC1 = s.c1().getVertices();
			vC1.remove(v);
			vC1.remove(v1);
			if (!vC1.isEmpty())
				g.transformVertices(vC1, s.c1(), new Edge(v1, v),
						g.getPoints(v1, v));
			// if (!vC1.isEmpty())
			// g.transformVertices(vC1,
			// // new EdgePos(g.getPoint(v1), graph.getPoint(v)),
			// graph.getPoints(v1, v), // ??? right?
			// g.getPoints(v1, v));

			HashSet<Vertex> vC2 = s.c2().getVertices();
			vC2.remove(v);
			vC2.remove(v2);
			if (!vC2.isEmpty())
				g.transformVertices(vC2, s.c2(), new Edge(v2, v),
						g.getPoints(v2, v));
			// if (!vC2.isEmpty())
			// g.transformVertices(vC2, graph.getPoints(v2, v),
			// g.getPoints(v2, v));
		}
		// System.out.println("old graph:" + graph.printPoints());
		// System.out.println("new graph:" + g.printPoints());
		return g;
	}

	public Realization tryRealize(RealizationType t) {
		Realization g = tryRealize(getBaseNoneedgeLength(), t);
		if (g != null)
			return g;

		boolean b = ccsGenerated();
		assert (b);
		OrientedCayleyConfigSpace occs = cayleyConfigSpace.getOrientedCCS(t);
		if (occs != null) {
			Interval it = occs.intervals.get(0);
			g = tryRealize((it.lower + it.upper) / 2, t);
		}

		return g;
	}

	// add accuracy
	public ArrayList<Realization> tryRealize(double length) {
		ArrayList<Realization> solutions = solve(length);
		if (solutions.isEmpty()) {
			solutions = solve(length + ACCURACY);
			if (solutions.isEmpty())
				solutions = solve(length - ACCURACY);
		}
		return solutions;

	}

	// TODO: efficient implementation
	public Realization tryRealize(double length, RealizationType forward) {
		Realization g = solve(length, forward);
		if (g == null) {
			g = solve(length + ACCURACY, forward);
			if (g == null)
				g = solve(length - ACCURACY, forward);
		}
		return g;

		// ArrayList<ListGraph> l = tryRealize(length);
		// ArrayList<ListGraph> solutions = new ArrayList<ListGraph>();
		// for (ListGraph g : l) {
		// if (getForwardSolutionType(g).compatible(forward)) {
		// solutions.add(g);
		// }
		// }
		// // if (solutions.size() > 1) {
		// // System.out.println("solutions for " + length);
		// // for (ListGraph g : solutions) {
		// // System.out.println(g.printPoints());
		// // }
		// // }
		// if (solutions.isEmpty()) {
		// l = tryRealize(length + accuracy);
		// for (ListGraph g : l) {
		// if (getForwardSolutionType(g).compatible(forward)) {
		// solutions.add(g);
		// }
		// }
		// if (solutions.isEmpty()) {
		// l = tryRealize(length - accuracy);
		// for (ListGraph g : l) {
		// if (getForwardSolutionType(g).compatible(forward)) {
		// solutions.add(g);
		// }
		// }
		// }
		// }
	}

	// add accuracy
	public boolean realizable(double length) {
		ArrayList<Realization> solutions = tryRealize(length);
		return !solutions.isEmpty();
	}

	public boolean realizable(double length, RealizationType forward) {
		Realization g = tryRealize(length, forward);
		return g != null;
	}

	/**
	 * ??????????????? Pretty sure that length should have solution. Iteratively
	 * try values near length
	 * 
	 * @param direction
	 *            -1: try values before length; +1: try values after length
	 * @return the nearest value to length having solution
	 */
	private double iterateTry(RealizationType t, double length, int direction) {
		final double STEP = (direction > 0 ? ACCURACY * 20 : -ACCURACY * 20);

		double tryLength = length + STEP;
		while (!realizable(tryLength, t))
			tryLength += STEP;
		return tryLength;
	}

	public void genCayleyConfigSpace() {
		cayleyConfigSpace = CayleyConfigSpace.generateCCS(this);
	}

	// assume candidates is sorted
	OrientedCayleyConfigSpace genOrientedCCSFromCandidates(ArrayList<Double> candidates,
			RealizationType type) {
		OrientedCayleyConfigSpace ccs = new OrientedCayleyConfigSpace(type);
		Double lastEndpoint = null;

		// TODO: numerical error exists -> unrealizable !
		double first = candidates.get(0);
		double second = candidates.get(1);

		// if (!realizable(first))
		// first += errorAdjust;
		if (realizable(first, type)) {
			if (realizable((first + second) / 2, type)) { // first is endpoint
				lastEndpoint = first;
				Debug.msg(first + " is interval start");
			} else { // first is isolated point
				ccs.appendInterval(first, first);
				Debug.msg(first + "is isolated point");
			}
		} else
			Debug.msg(first + " not realizable");

		for (int i = 1; i < candidates.size() - 1; ++i) {
			double cur = candidates.get(i);
			Debug.msg("Processing Candidate: " + cur);
			// if (!realizable(cur)) {
			// cur -= errorAdjust;
			// if (!realizable(cur))
			// cur += 2 * errorAdjust;
			// }
			if (!realizable(cur, type)) {
				Debug.msg(cur + " not realizable");

				// ad-hoc modification for assertion error???
				/*
				 * if (lastEndpoint != null) { // current candidate should be an
				 * end point... double tryCur = iterateTry(type, cur, -1);
				 * assert (tryCur > lastEndpoint); assert (cur - tryCur < 1);
				 * Debug.warnMsg(tryCur + " is realizable !"); } else {
				 */
				assert (lastEndpoint == null);
				continue;
				// }
			}
			double prev = candidates.get(i - 1);
			double next = candidates.get(i + 1);
			double p = (prev + cur) / 2;
			double n = (cur + next) / 2;
			boolean P = realizable(p, type);
			boolean N = realizable(n, type);
			if (!P && !N) {
				ccs.appendInterval(cur, cur);
				Debug.msg(cur + " is isolated point");
			} else if (P && N) {
				assert (lastEndpoint != null);
				Debug.msg(cur + " in middle of interval");
			} else if (P && !N) {
				// ad-hoc modification for assertion error????
				/*
				 * if (lastEndpoint == null) { // prev candidate should be a
				 * start point... double tryPrev = iterateTry(type, prev, +1);
				 * Debug.warnMsg(tryPrev + " is realizable as interval start!");
				 * assert (tryPrev < p); assert (tryPrev - prev < 1);
				 * lastEndpoint = prev; }
				 */
				assert (lastEndpoint != null);
				ccs.appendInterval(lastEndpoint, cur);
				lastEndpoint = null;
				Debug.msg(cur + " is interval end");
			} else {
				assert (lastEndpoint == null);
				assert (!P && N);
				lastEndpoint = cur;
				Debug.msg(cur + " is interval start");
			}
		}

		double last = candidates.get(candidates.size() - 1);
		double lastbutone = candidates.get(candidates.size() - 2);
		// if (!realizable(last))
		// last -= errorAdjust;
		if (realizable(last, type)) {
			if (realizable((lastbutone + last) / 2, type)) { // last is endpoint
				assert (lastEndpoint != null);
				ccs.appendInterval(lastEndpoint, last);
				Debug.msg(last + " is interval end");
			} else {
				ccs.appendInterval(last, last);
				Debug.msg(last + " is isolated point");
			}
		} else
			Debug.msg(last + " not realizable");

		return ccs;
	}

	public boolean ccsGenerated() {
		return cayleyConfigSpace != null;
	}

	// ========================= continuous path ==============
	public ContinuousMotionPath findPath(double l1, double l2) {
		ContinuousMotionPath result = null;
		ArrayList<Realization> solutions1 = tryRealize(l1);
		Realization g = tryRealize(l1, getForwardSolutionType());
		if (g != null)
			solutions1.add(0, g);
		ArrayList<Realization> solutions2 = tryRealize(l2);
		g = tryRealize(l2, getForwardSolutionType());
		if (g != null)
			solutions2.add(0, g);
		for (Realization g1 : solutions1)
			for (Realization g2 : solutions2) {
				ContinuousMotionPath path = findPath(g1, g2);
				if (path != null) {
					// System.out.println("path length:" + path.size() + "; " +
					// result);
					if (result == null || result.size() > path.size())
						result = path;
				}
				// return path;
			}
		return result;
	}

	public ContinuousMotionPath findPath(Realization g1, Realization g2) {
		return ContinuousMotionPath.findPath(this, g1, g2);
	}

	public ConnectedComponent findPathFrom(Realization g1) {
		return ConnectedComponent.findComponent(this, g1);
	}

	public int isStepVertex(Vertex v) {
		for (int i = 0; i < this.getNumOfConstructStep(); ++i) {
			ConstructionStep s = constructionSequence.get(i);
			if (s.stepVertex == v)
				return i;
		}
		return -1;
	}

	/**
	 * Treat the length of canonical base non-edges as a vector, and compute the
	 * vector distance between g1 and g2
	 */
	public double cayleyDistance(Realization g1, Realization g2) {
		if (g1.equals(g2))
			return 0;

		// sqrt( sum[(ai - bi)^2] )
		double sum = 0;
		for (Edge e : completeCayleyVector()) {
			double l1 = g1.length(e);
			double l2 = g2.length(e);
			sum += (l1 - l2) * (l1 - l2);
		}
		return Math.sqrt(sum);
	}

	public String toString() {
		return "f:" + this.getBaseNonedge() + "\n" + graph.toString();
	}
}