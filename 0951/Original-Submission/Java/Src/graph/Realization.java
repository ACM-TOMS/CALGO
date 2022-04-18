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

import java.awt.geom.AffineTransform;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;

public class Realization implements Serializable {
	public static final boolean DEBUG = true;

	// TODO: resize? screen to mouse?
	// static int width = 800, height = 600;

	// vertices
	private ArrayList<Vertex> vertices;

	// point coordinates for each vertex
	private HashMap<Vertex, Point2D> points;
	private HashMap<Vertex, Point2D> initPoints;

	// adjacent-list, indexed by v
	private HashMap<Vertex, LinkedList<Vertex>> vNeighbors;

	public Realization() {
		this.vertices = new ArrayList<Vertex>();
		this.points = new HashMap<Vertex, Point2D>();
		this.initPoints = new HashMap<Vertex, Point2D>();
		this.vNeighbors = new HashMap<Vertex, LinkedList<Vertex>>();
	}

	public int getSize() {
		return vertices.size();
	}

	public LinkedList<Vertex> getNeighbors(Vertex v) {
		return vNeighbors.get(v);
	}

	public Point2D getInitPoint(Vertex v) {
		return initPoints.get(v);
	}

	public Point2D getPoint(Vertex v) {
		assert (points.containsKey(v));
		return points.get(v);
	}

	public EdgePos getPoints(Edge e) {
		return new EdgePos(getPoint(e.v1()), getPoint(e.v2()));
	}

	public EdgePos getPoints(Vertex v1, Vertex v2) {
		return getPoints(new Edge(v1, v2));
	}

	public void setPoint(Vertex v, Point2D newloc) {
		assert (points.containsKey(v));
		points.put(v, newloc);
	}

	public void setPoint(Vertex v, double x, double y) {
		setPoint(v, new Point2D(x, y));
	}

	public void resetInitPoints() {
		for (Vertex v : vertices) {
			setInitPoint(v, getPoint(v));
		}
	}

	public void setInitPoint(Vertex v, Point2D newloc) {
		assert (initPoints.containsKey(v));
		initPoints.put(v, newloc);
	}

	public void setInitPoint(Vertex v, double x, double y) {
		setInitPoint(v, new Point2D(x, y));
	}

	public Vertex getRandomVertex() {
		return vertices.get((int) (Math.random() * vertices.size()));
	}

	// insert existing vertex
	public void insertVertex(Vertex v, Point2D p) {
		if (vertices.contains(v))
			return;
		vertices.add(v);
		points.put(v, p);
		initPoints.put(v, p);
		vNeighbors.put(v, new LinkedList<Vertex>());
	}

	public void insertVertex(Vertex v, double x, double y) {
		Point2D p = new Point2D(x, y);
		insertVertex(v, p);
	}

	// create a new vertex
	public Vertex addVertex(Point2D p) {
		int index = vertices.size();
		Vertex v = new Vertex(index);
		insertVertex(v, p);
		return v;
	}

	public Vertex addVertex(double x, double y) {
		Point2D p = new Point2D(x, y);
		return this.addVertex(p);
	}

	public Vertex removeVertex(Vertex v) {
		LinkedList<Vertex> neighbors = this.getNeighbors(v);
		for (Vertex u : neighbors) {
			getNeighbors(u).remove(v);
		}
		neighbors.clear();
		vertices.remove(v);
		vNeighbors.remove(v);

		points.remove(v);
		initPoints.remove(v);

		for (int i = 0; i < vertices.size(); ++i)
			vertices.get(i).index = i;

		return v;
	}

	public void addEdge(Vertex v1, Vertex v2) {
		if (v1 == v2 || isAdjacent(v1, v2))
			return;
		assert (vertices.contains(v1));
		assert (vertices.contains(v2));
		getNeighbors(v1).add(v2);
		getNeighbors(v2).add(v1);
	}

	public void addEdge(Edge e) {
		addEdge(e.v1(), e.v2());
	}

	public void removeEdge(Vertex v1, Vertex v2) {
		assert (this.isAdjacent(v1, v2));
		getNeighbors(v1).remove(v2);
		getNeighbors(v2).remove(v1);
	}

	public void removeEdge(Edge e) {
		removeEdge(e.v1(), e.v2());
	}

	public boolean isAdjacent(Vertex v1, Vertex v2) {
		if (v1 == v2)
			return false;
		assert (vertices.contains(v1));
		assert (vertices.contains(v2));
		boolean r = (getNeighbors(v1).contains(v2));
		assert (r == getNeighbors(v2).contains(v1));
		return r;
	}

	public boolean isAdjacent(Edge edge) {
		return isAdjacent(edge.v1(), edge.v2());
	}

	@Override
	public String toString() {
		String s = new String();
		// int n = vertices.size();
		for (Vertex v : vertices) {
			s += v + ": ";
			for (Vertex u : vNeighbors.get(v))
				s += u + " ";
			s += "\n";
		}
		return s;
	}

	public String pointsToString() {
		String s = "";
		for (Vertex v : vertices) {
			s += v.toString() + getPoint(v) + "\t";
		}
		return s;
	}

	// returns a subgraph induced by vertices in vset
	// share vertices, copied points
	public Realization inducedSubgaph(ArrayList<Vertex> vset) {
		int n = vset.size();
		Realization g = new Realization();
		for (Vertex v : vset)
			g.insertVertex(v, getPoint(v));
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < i; ++j) {
				Vertex v1 = vset.get(i), v2 = vset.get(j);
				if (this.isAdjacent(v1, v2))
					g.addEdge(v1, v2);
			}
		}
		return g;
	}

	public Realization clone() {
		return this.inducedSubgaph(vertices);
		// Graph g = new Graph();
		// for (Vertex v : vertices)
		// g.insertVertex(v, getPoint(v));
		// for (int i = 0; i < getSize(); ++i)
		// for (int j = 0; j < i; ++j) {
		// Vertex v1 = getVertex(i), v2 = getVertex(j);
		// if (isAdjacent(v1, v2))
		// g.addEdge(v1, v2);
		// }
		// return g;
	}

	public HashMap<Vertex, Point2D> getPoints() {
		HashMap<Vertex, Point2D> map = new HashMap<Vertex, Point2D>();
		map.putAll(points);
		return map;
	}

	public void setPoints(Realization g) {
		assert (g.getSize() == this.getSize());
		for (int i = 0; i < getSize(); ++i) {
			Point2D newp = g.getPoint(g.vertices.get(i));
			this.setPoint(vertices.get(i), newp);
		}
	}

	public void setPoints(HashMap<Vertex, Point2D> map) {
		assert (map.size() == points.size());
		for (Vertex v : map.keySet()) {
			assert (points.containsKey(v));
			points.put(v, map.get(v));
		}
	}

	public Vertex getVertex(int i) {
		assert (i < vertices.size());
		return vertices.get(i);
	}

	public Vertex getVertex(Point2D p, double precision) {
		for (Vertex v : vertices) {
			if (getPoint(v).distance(p) < precision) {
				return v;
			}
		}
		return null;
	}

	public double initDistance(Vertex v1, Vertex v2) {
		Point2D p1 = getInitPoint(v1), p2 = getInitPoint(v2);
		return p1.distance(p2);
	}

	public double initLength(Edge e) {
		return initDistance(e.v1(), e.v2());
	}

	public double distance(Vertex v1, Vertex v2) {
		Point2D p1 = getPoint(v1), p2 = getPoint(v2);
		return p1.distance(p2);
	}

	public double length(Edge e) {
		return distance(e.v1(), e.v2());
	}

	public int orientationOf(Vertex v1, Vertex v2, Vertex v) {
		return Point2D.orientationOf(getPoint(v1), getPoint(v2), getPoint(v));
	}

	public void transformVertex(Vertex v, AffineTransform transform) {
		Point2D newp = getPoint(v).transform(transform);
		setPoint(v, newp);
	}

	public void transformVertices(AffineTransform transform) {
		for (Vertex v : vertices) {
			Point2D newp = getPoint(v).transform(transform);
			setPoint(v, newp);
		}
	}

	public void transformVertices(Collection<Vertex> list, Cluster c, Edge e,
			EdgePos newPos) {
		/*HashMap<Vertex, Point2D> newpoints = c.getTransform(e, newPos);*/

		for (Vertex v : list) {
			
			Vertex v1 = e.v1(), v2 = e.v2();
			Point2D p1 = newPos.p1(), p2 = newPos.p2();
			double l1 = initDistance(v1, v), l2 = initDistance(v2, v);
			int orient = Point2D.orientationOf(getInitPoint(v1),
					getInitPoint(v2), getInitPoint(v));
			Point2D newp = Util.solve(p1, p2, l1, l2, orient);
			this.setPoint(v, newp);
			
			// System.out.print(v + "'s old pos:" + getPoint(v) + "; ");
			// transformVertex(v, t);
			/*this.setPoint(v, newpoints.get(v));*/
			// System.out.println(" new pos:" + getPoint(v));
		}
	}

	public void writeToStream(ObjectOutputStream o) throws IOException {
		int n = getSize();
		o.writeInt(n); // write size

		// write positions
		for (int i = 0; i < n; ++i) {
			o.writeObject(getPoint(getVertex(i)));
		}

		// write adjacency matrix
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < i; ++j) {
				if (isAdjacent(getVertex(i), getVertex(j)))
					o.writeInt(1);
				else
					o.writeInt(0);
			}
		}
	}

	public static Realization readFromStream(ObjectInputStream in)
			throws IOException, ClassNotFoundException {
		Realization g = new Realization();
		int n = in.readInt(); // read size

		// read vertices
		for (int i = 0; i < n; ++i) {
			Point2D p = (Point2D) in.readObject();
			g.addVertex(p);
		}

		// read adjacency matrix
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < i; ++j) {
				int adjacent = in.readInt();
				if (adjacent == 1)
					g.addEdge(g.getVertex(i), g.getVertex(j));
			}
		}

		return g;
	}
}
