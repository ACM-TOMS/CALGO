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

import ccs.graph.Realization;
import ccs.graph.Point2D;

public class BoundingBox {

	double x_min, x_max, y_min, y_max;

	public double getWidth() {
		return x_max - x_min;
	}

	public double getHeight() {
		return y_max - y_min;
	}

	public BoundingBox(double x1, double x2, double y1, double y2) {
		x_min = x1;
		x_max = x2;
		y_min = y1;
		y_max = y2;
	}

	public BoundingBox(Realization g) {
		x_min = Double.POSITIVE_INFINITY;
		x_max = Double.NEGATIVE_INFINITY;
		y_min = Double.POSITIVE_INFINITY;
		y_max = Double.NEGATIVE_INFINITY;

		for (int i = 0; i < g.getSize(); ++i) {
			Point2D p = g.getPoint(g.getVertex(i));
			if (p.x() < x_min)
				x_min = p.x();
			if (p.x() > x_max)
				x_max = p.x();
			if (p.y() < y_min)
				y_min = p.y();
			if (p.y() > y_max)
				y_max = p.y();
		}
	}

	public String toString() {
		return "[" + x_min + "," + x_max + ";" + y_min + "," + y_max + "]";
	}
}
