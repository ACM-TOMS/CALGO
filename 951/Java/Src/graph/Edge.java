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

public class Edge extends Pair<Vertex> {

	Edge(Vertex o1, Vertex o2) {
		super(o1, o2);
	}

	public Vertex v1() {
		return getFirst();
	}

	public Vertex v2() {
		return getSecond();
	}

	@Override
	public String toString() {
		return "(" + v1() + ", " + v2() + ")";
	}

	@Override
	public boolean equals(Object that) {
		if (!(that instanceof Edge))
			return false;
		Edge e = (Edge) that;
		return (v1() == e.v1() && v2() == e.v2())
				|| (v1() == e.v2() && v2() == e.v1());
	}

	@Override
	public int hashCode() {
		if (v1().index < v2().index)
			return (v1().hashCode() + ":" + v2().hashCode()).hashCode();
		else
			return (v2().hashCode() + ":" + v1().hashCode()).hashCode();
	}
}

class EdgePos extends Pair<Point2D> {

	public EdgePos(Point2D p1, Point2D p2) {
		super(p1, p2);
	}

	public Point2D p1() {
		return getFirst();
	}

	public Point2D p2() {
		return getSecond();
	}

	public double x1() {
		return getFirst().x();
	}

	public double x2() {
		return getSecond().x();
	}

	public double y1() {
		return getFirst().y();
	}

	public double y2() {
		return getSecond().y();
	}

	@Override
	public String toString() {
		return "(" + p1() + ", " + p2() + ")";
	}

	public double length() {
		return getFirst().distance(getSecond());
	}

	// return ccw angle from e1 to e2
	public double angle(EdgePos that) {
		return new Vector2D(p1(), p2())
				.angle(new Vector2D(that.p1(), that.p2()));
	}
}

class ClusterPair extends Pair<Cluster> {

	public ClusterPair(Cluster o1, Cluster o2) {
		super(o1, o2);
	}

	public Cluster c1() {
		return getFirst();
	}

	public Cluster c2() {
		return getSecond();
	}
}