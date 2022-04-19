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

public class OrientedCayleyConfigSpace implements CCSInterface{
	RealizationType forwardSolutionType;
	ArrayList<Interval> intervals;

	// private void sort() {
	// Collections.sort(intervals);
	// }
	
	public RealizationType getSolutionType(){
		return forwardSolutionType;
	}
	
	public Interval contains(double val) {
		for (Interval i : intervals) {
			if (i.contains(val))
				return i;
		}
		return null;
	}

	public ArrayList<Interval> getIntervals() {
		return intervals;
	}

	public double getMin() {
		if (intervals.isEmpty())
			return java.lang.Double.NEGATIVE_INFINITY;
		else
			return intervals.get(0).lower;
	}

	public boolean isEmpty() {
		return intervals.isEmpty();
	}

	public double getMax() {
		if (intervals.isEmpty())
			return java.lang.Double.POSITIVE_INFINITY;
		else
			return intervals.get(intervals.size() - 1).upper;
	}

	public OrientedCayleyConfigSpace(RealizationType type) {
		this.forwardSolutionType = type;
		intervals = new ArrayList<Interval>();
	}

	// append interval at the end
	public void appendInterval(double lower, double upper) {
		assert (lower > getMin());
		intervals.add(new Interval(lower, upper));
	}

	// TODO: arraylist + binary search?
	public Interval getContainingInterval(double value) {
		for (Interval interval : intervals) {
			if (interval.contains(value))
				return interval;
		}
		return null;
	}

	@Override
	public String toString() {
		String s = "";
		for (Interval interval : intervals) {
			s += interval.toString();
		}
		return s;
	}
}