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

import java.awt.Graphics;
import java.awt.geom.AffineTransform;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;

public class Cluster {
	private Realization graph;

	private HashSet<Vertex> vertices;
	private HashSet<Vertex> sharedVertices;

	// TODO: always store the inital points
	// private HashMap<Vertex, Point2D> initPoints;

	public Cluster(Realization g) {
		graph = g;

		vertices = new HashSet<Vertex>();
		sharedVertices = new HashSet<Vertex>();
		// initPoints = new HashMap<Vertex, Point2D>();
	}

	public void addVertex(Vertex v, Point2D p) {
		vertices.add(v);
		// initPoints.put(v, p);
	}

	public String toString() {
		return vertices.toString();
	}

	// for 1-dof t-d only
	public static Vertex sharedVertex(Cluster c1, Cluster c2) {
		HashSet<Vertex> set = new HashSet<Vertex>();
		set.addAll(c1.sharedVertices);
		set.retainAll(c2.sharedVertices);
		assert (set.size() <= 1);
		if (set.size() > 0)
			return set.iterator().next();
		return null;
	}

	// for 1-dof t-d only
	public static ArrayList<Vertex> sharedVertices(ArrayList<Cluster> clusters) {
		HashSet<Vertex> vset = new HashSet<Vertex>();
		for (int i = 0; i < clusters.size(); ++i) {
			for (int j = 0; j < i; ++j) {
				Cluster c1 = clusters.get(i), c2 = clusters.get(j);
				HashSet<Vertex> set = new HashSet<Vertex>();
				set.addAll(c1.sharedVertices);
				set.retainAll(c2.sharedVertices);
				assert (set.size() <= 1);
				if (set.size() > 0) {
					Vertex v = set.iterator().next();
					vset.add(v);
				}
			}
		}
		ArrayList<Vertex> vertices = new ArrayList<Vertex>();
		vertices.addAll(vset);
		return vertices;
	}

	public static boolean mergable(Cluster c1, Cluster c2, Cluster c3) {
		assert (!c1.equals(c2) && !c1.equals(c3));
		HashSet<Vertex> set = new HashSet<Vertex>();
		set.addAll(c1.vertices);
		set.retainAll(c2.vertices);
		if (set.size() != 1)
			return false;
		Vertex v1 = set.iterator().next();
		set.clear();
		set.addAll(c1.vertices);
		set.retainAll(c3.vertices);
		if (set.size() != 1)
			return false;
		Vertex v2 = set.iterator().next();
		if (v1.equals(v2))
			return false;
		set.clear();
		set.addAll(c2.vertices);
		set.retainAll(c3.vertices);
		if (set.size() != 1)
			return false;
		Vertex v3 = set.iterator().next();
		if (v1.equals(v3) || v2.equals(v3))
			return false;
		return true;
	}

	public static Cluster merge(Cluster c1, Cluster c2, Cluster c3) {
		assert (c1.graph == c2.graph && c1.graph == c3.graph);
		Cluster c = new Cluster(c1.graph);
		c.vertices.addAll(c1.vertices);
		c.vertices.addAll(c2.vertices);
		c.vertices.addAll(c3.vertices);
		// c.initPoints.putAll(c1.initPoints);
		// c.initPoints.putAll(c2.initPoints);
		// c.initPoints.putAll(c3.initPoints);
		return c;
	}

	public static void merge(LinkedList<Cluster> clusters) {
		boolean modified = true;
		while (modified) {
			modified = false;
			forloop: for (int i = 0; i < clusters.size(); ++i) {
				for (int j = 0; j < i; ++j) {
					for (int k = 0; k < j; ++k) {
						Cluster c1 = clusters.get(i), c2 = clusters.get(j), c3 = clusters
								.get(k);
						if (Cluster.mergable(c1, c2, c3)) {
							Cluster c = Cluster.merge(c1, c2, c3);
							clusters.remove(c1);
							clusters.remove(c2);
							clusters.remove(c3);
							clusters.add(c);
							// if (DEBUG)
							// System.out.println(c1 + " " + c2 + " " + c3
							// + " -> " + c);
							modified = true;
							break forloop;
						}
					}
				}
			}
		}
	}

	public void addSharedVertex(Vertex v) {
		assert (vertices.contains(v));
		this.sharedVertices.add(v);
	}

	public boolean containsSharedVertex(Vertex v) {
		return sharedVertices.contains(v);
	}

	public int getNumOfSharedVertices() {
		return sharedVertices.size();
	}

	public int getSize() {
		return vertices.size();
	}

	public HashSet<Vertex> getVertices() {
		HashSet<Vertex> set = new HashSet<Vertex>();
		set.addAll(vertices);
		return set;
	}

	public HashSet<Vertex> unionSharedVertices(Cluster that) {
		HashSet<Vertex> set = new HashSet<Vertex>();
		set.addAll(this.sharedVertices);
		set.addAll(that.sharedVertices);
		return set;
	}

	public HashSet<Vertex> unionVertices(Cluster that) {
		HashSet<Vertex> set = new HashSet<Vertex>();
		set.addAll(this.vertices);
		set.addAll(that.vertices);
		return set;
	}

	// public double getDistance(Edge e) {
	// return initPoints.get(e.v1()).distance(initPoints.get(e.v2()));
	// }

	public HashMap<Vertex, Point2D> getTransform(Edge e, EdgePos newPos) {
		// generate a transformation from oldpos -> newpos
		// TODO: incorrect?

		// TODO: first try to move one point together?
		// EdgePos oldPos = new EdgePos(initPoints.get(e.v1()), initPoints.get(e
		// .v2()));
		EdgePos oldPos = new EdgePos(graph.getInitPoint(e.v1()),
				graph.getInitPoint(e.v2()));

		// System.out.println(initPoints);
		// System.out.println(oldPos);
		assert (Math.abs(oldPos.length() - newPos.length()) < TDLinkage.ACCURACY * 10);
		//System.out.println("transform cluster:" + this);
		//System.out.println("old pos/new pos: " + oldPos + "/" + newPos);

		AffineTransform r = new AffineTransform();
		double theta = newPos.angle(oldPos);
		r.rotate(-theta, oldPos.x1(), oldPos.y1());
		//System.out.println("ccw angle: " + theta / Math.PI + "PI");

		AffineTransform t = new AffineTransform();
		t.translate(newPos.x1() - oldPos.x1(), newPos.y1() - oldPos.y1());
		//System.out.println("translate (" + (newPos.x1() - oldPos.x1()) + ","
		//		+ (newPos.y1() - oldPos.y1()) + ")");

		t.concatenate(r);

		HashMap<Vertex, Point2D> points = new HashMap<Vertex, Point2D>();
		for (Vertex v : getVertices()) {
			// System.out.print(v + "'s old pos:" + getPoint(v) + "; ");
			//Point2D p = initPoints.get(v).transform(t);// transformVertex
			Point2D p = graph.getInitPoint(v).transform(t);
			// System.out.println(" new pos:" + getPoint(v));
			points.put(v, p);
		}

		return points;
	}

}