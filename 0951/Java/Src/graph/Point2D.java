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
import java.awt.geom.Point2D.Double;
import java.io.Serializable;

//import org.jscience.mathematics.number.Real;

/* Immutable Object */

public class Point2D implements Serializable {

	public static Point2D orig = new Point2D(0, 0);

	private Double p;

	public Point2D(double x, double y) {
		p = new Double(x, y);
	}

	public Point2D(Point2D point) {
		this(point.x(), point.y());
	}

	public double x() {
		return p.getX();
	}

	public double y() {
		return p.getY();
	}

	public Vector2D minus(Point2D that) {
		return new Vector2D(that, this);
	}

	public Point2D add(Vector2D v) {
		return new Point2D(x() + v.x(), y() + v.y());
	}

	public double distance(Point2D that) {
		return p.distance(that.p);
	}

	public double distance(double x, double y) {
		return p.distance(x, y);
	}

	/**
	 * @return +1: v1->v is clockwise from v1->v2; -1: counterclockwise; 0:
	 *         collinear
	 */
	public static int orientationOf(Point2D v1, Point2D v2, Point2D v) {
		// (v - v1) X (v2 - v1) = (x - x1)(y2 - y1) - (x2 - x1)(y - y1)
		// If this cross product is positive, then v1->v is clockwise from
		// v1->v2; if negative, it is counterclockwise.
		double cross = (v.x() - v1.x()) * (v2.y() - v1.y()) - (v2.x() - v1.x())
				* (v.y() - v1.y());
		// ??? add accuracy ???
		final double distanceError = 2;
		if (Math.abs(cross) < TDLinkage.ACCURACY
				|| Util.distance(v, v1, v2) < distanceError) {
			cross = 0;
		}
		int r = (int) Math.signum(cross);
		return r;
	}

	public Point2D transform(AffineTransform t) {
		Double dst = (Double) t.transform(p, null);
		return new Point2D(dst.getX(), dst.getY());
	}

	public String toString() {
		return "(" + x() + "," + y() + ")";
	}
}
