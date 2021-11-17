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

import java.awt.Color;
import java.util.LinkedList;
import java.util.Random;

import ccs.Debug;

public class Util {

	static int border = 20;

	private static Point2D randomPoint() {
		int x = (int) (Math.random() * (800 - 2 * border)) + border;
		int y = (int) (Math.random() * (600 - 2 * border)) + border;
		return new Point2D(x, y);
	}

	public static Realization generateListGraph(int n) {
		Realization g = new Realization();
		for (int i = 0; i < n; ++i) {
			g.addVertex(randomPoint());
		}

		for (int i = 0; i < g.getSize(); ++i) {
			Vertex v = g.getVertex(i);
			int m = (int) (Math.random() * n);
			for (int j = 0; j < m; ++j) {
				Vertex u;
				do {
					u = g.getRandomVertex();
				} while (v == u);
				if (!g.isAdjacent(v, u)) {
					if (Realization.DEBUG)
						System.out.println("add edge: (" + v + ", " + u + ")");
					g.addEdge(v, u);
				}
			}
		}
		return g;
	}

	// Generate an arbitrary 1-dof Henneberg-I graph
	public static Realization generateHennebergGraph(int n, boolean isTriangleFree) {
		// base non-edge:
		Realization g = new Realization();
		for (int i = 0; i < Math.min(2, n); ++i)
			g.addVertex(randomPoint());
		if (n <= 2)
			return g;

		for (int i = 0; i < n - 2; ++i) {
			Vertex v1, v2;
			boolean looping = true;
			do {
				v1 = g.getRandomVertex();
				v2 = g.getRandomVertex();

				looping = (v1 == v2);
				if (isTriangleFree)
					looping = (looping || g.isAdjacent(v1, v2));
			} while (looping);

			Vertex v = g.addVertex(randomPoint());

			g.addEdge(v1, v);
			g.addEdge(v, v2);

			if (Realization.DEBUG)
				System.out.println(v + " on (" + v1 + "," + v2 + ")");
		}
		return g;
	}

	// does not modify g
	public static LinkedList<Cluster> getDecompose(Realization g) {
		// generate a list of clusters

		// TODO: use Fudo's method
		// first attempt: brute-force: merge to get all clusters
		LinkedList<Cluster> clusters = new LinkedList<Cluster>();

		// Initialize: for each edge: add a cluster
		// !!! for each isolated vertex: add a cluster
		int n = g.getSize();
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < i; ++j) {
				Vertex v1 = g.getVertex(i), v2 = g.getVertex(j);
				if (g.isAdjacent(v1, v2)) {
					Cluster c = new Cluster(g);
					c.addVertex(v1, g.getPoint(v1));
					c.addVertex(v2, g.getPoint(v2));
					clusters.add(c);
				}
			}
		for (int i = 0; i < n; ++i) {
			Vertex v = g.getVertex(i);
			if (g.getNeighbors(v).isEmpty()) {
				Cluster c = new Cluster(g);
				c.addVertex(v, g.getPoint(v));
				clusters.add(c);
			}
		}
		// Merge step
		Cluster.merge(clusters);
		return clusters;
	}

	// public static double twoDigits(double val) {
	// val *= 100;
	// int n = (int) val;
	// return n / 100.0;
	// }

	public static Color randomColor(int seed) {
		Random rand = new Random();
		rand.setSeed(seed);
		float r = rand.nextFloat();
		float g = rand.nextFloat();
		float b = rand.nextFloat();
		return new Color(r, g, b);
	}

	// distance from point v to line (v1,v2)
	public static double distance(Point2D v, Point2D v1, Point2D v2) {
		double a = v2.y() - v1.y();
		double b = v1.x() - v2.x();
		double c = (v2.x() - v1.x()) * v1.y() - (v2.y() - v1.y()) * v1.x();
		double d = Math.abs(a * v.x() + b * v.y() + c)
				/ Math.sqrt(a * a + b * b);
		return d;
	}

	// public static int orientationOf(Vertex v1, Vertex v2, Vertex v) {
	// return orientationOf(v1.pos, v2.pos, v.pos);
	// }

	// ???? necessary?
	// public static boolean forwardSolutionTypeCompatible(TreeDecomp superT,
	// TreeDecomp subT) {
	// assert (superT.constructionSequenceGenerated());
	// assert (subT.graph.superGraph == superT.graph);
	// int n = subT.getNumOfConstructStep();
	// for (int i = 0; i < n; ++i) {
	// ConstructionStep step = superT.constructionSequence.get(i);
	// Vertex v1 = step.getV1(), v2 = step.getV2(), v = step.getV();
	// Vertex sv1 = subT.graph.getVertexFromSuperGraphVertex(v1), sv2 =
	// subT.graph
	// .getVertexFromSuperGraphVertex(v2), sv = subT.graph
	// .getVertexFromSuperGraphVertex(v);
	// if (Util.orientationOf(sv1, sv2, sv) != Util.orientationOf(v1, v2,
	// v))
	// return false;
	// }
	// return true;
	// }

	// ???? put it here or not
	// public static boolean forwardSolutionTypeSame(TreeDecomp t, ListGraph g1,
	// ListGraph g2) {
	// assert (t.constructionSequenceGenerated());
	// assert (g1.getSize() == g2.getSize());
	// for (ConstructionStep step : t.constructionSequence) {
	// int i1 = step.getV1().index, i2 = step.getV2().index, i3 = step
	// .getV().index;
	// Vertex v1 = g1.getVertex(i1), v2 = g1.getVertex(i2), v = g1
	// .getVertex(i3);
	// Vertex u1 = g2.getVertex(i1), u2 = g2.getVertex(i2), u = g2
	// .getVertex(i3);
	// if (Util.orientationOf(v1, v2, v) != Util.orientationOf(u1, u2, u))
	// return false;
	// }
	// return true;
	// }

	/**
	 * @return whether value is between end1 and end2
	 */
	public static boolean between(double value, double end1, double end2) {
		if (end1 <= end2)
			return end1 <= value && value <= end2;
		else
			return end1 >= value && value >= end2;
	}

	/**
	 * Given coordinates of two base vertices, two bar lengths and local
	 * orientation, solve for a third point constructed using ruler-and-compass
	 * 
	 * @param orient
	 *            +1: v1->v clockwise from v1->v2; -1: counterclockwise; 0:
	 *            collinear | don't care
	 * @return null if no solution exists
	 */
	public static Point2D solve(Point2D p1, Point2D p2, double r1, double r2,
			int orient) {

		// TODO:
		double r3 = p1.distance(p2);
		assert (r3 != 0);

		// a = (r0^2 - r1^2 + d^2 ) / (2 d)
		double a = (r1 * r1 - r2 * r2 + r3 * r3) / (2 * r3);
		// double a = ((2 * pp1.x() * pp1.x() - 2 * pp1.x() * pp.x() + 2 *
		// pp1.y()
		// * pp1.y() - 2 * pp1.y() * pp.y() + 2 * pp2.x() * pp.x() + 2
		// * pp2.y() * pp.y() - 2 * pp1.x() * pp2.x() - 2 * pp1.y()
		// * pp2.y()) / (2 * r3));

		if (r1 * r1 < a * a) {
			// if (Graph.DEBUG) System.out.println("quit for r1 < a:" + a);
			return null;
		}
		// h^2 = r0^2 - a^2
		double h = Math.sqrt(r1 * r1 - a * a);
		// if (Graph.DEBUG) System.out.println("h,a: " + h + "," + a);

		// P2 = P0 + a ( P1 - P0 ) / d
		Point2D p3 = p1.add(p2.minus(p1).multiply(a).divide(r3));

		// if (Graph.DEBUG) System.out.println("p2x,p2y: " + p2x + "," +
		// p2y);

		double Xc = p3.x() + h * (p2.y() - p1.y()) / r3;
		double Yc = p3.y() - h * (p2.x() - p1.x()) / r3;
		// System.out.println("Xc,Yc: " + Xc + "," + Yc); // orientation +1
		double Xd = p3.x() - h * (p2.y() - p1.y()) / r3;
		double Yd = p3.y() + h * (p2.x() - p1.x()) / r3;
		// System.out.println("Xd,Yd: " + Xd + "," + Yd); // orientation -1

		Point2D pv = new Point2D(Xc, Yc);

		int o = Point2D.orientationOf(p1, p2, pv);
		if (o == 0 || orient == 0 || o == orient) {
			if (orient == 0)
				Debug.warnMsg("???? orient == 0 ????");
		} else {
			pv = new Point2D(Xd, Yd);
			o = Point2D.orientationOf(p1, p2, pv);
		}

		return pv;
	}
}
