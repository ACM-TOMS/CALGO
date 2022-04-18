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

import java.util.List;

public class ContinuousMotionPath extends ContinuousMotion {

	public ContinuousMotionPath(TDLinkage t) {
		super(t);

		// tempNodeSearchMap = new ArrayList<Double>();
	}

	// accumulated length -> node
	/*
	 * private ArrayList<Double> tempNodeSearchMap; // private double[]
	 * nodeSearchMap; // private double totalLength = 0;
	 */

	// TODO: WHAT IS NEAREST??
	/*
	 * public Pair<MotionPath> getNearestPoint(MotionPath that) { double minDiff
	 * = Math.abs(that.endNode().lf - this.endNode().lf); Node endthis =
	 * this.endNode(), endthat = that.endNode(); double thislf =
	 * this.endNode().lf, thatlf = that.endNode().lf; for (Node n1 : this.nodes)
	 * { for (Node n2 : that.nodes) { Interval i1 = (n1.direction == 1 ? new
	 * Interval(n1.getStart(), n1.getEnd()) : new Interval(n1.getEnd(),
	 * n1.getStart())); Interval i2 = (n2.direction == 1 ? new
	 * Interval(n2.getStart(), n2.getEnd()) : new Interval(n2.getEnd(),
	 * n2.getStart())); if (Interval.hasIntersection(i1, i2)) { // TODO: what to
	 * do? return null; } double diff = (i1.upper < i2.lower ? i2.lower -
	 * i1.upper : i1.lower - i2.upper); if (diff < minDiff) { minDiff = diff;
	 * endthis = n1; endthat = n2; if (i1.upper < i2.lower) { thislf = i1.upper;
	 * thatlf = i2.lower; } else { thislf = i1.lower; thatlf = i2.upper; } } } }
	 * MotionPath subthis = new MotionPath(t); for (Node n : this.nodes) { if (n
	 * == endthis) { subthis.addEndNode(thislf, n.interval, n.o, n.direction);
	 * break; } subthis.nodes.add(n); } MotionPath subthat = new MotionPath(t);
	 * for (Node n : that.nodes) { if (n == endthat) {
	 * subthat.addEndNode(thatlf, n.interval, n.o, n.direction); break; }
	 * subthat.nodes.add(n); } return new Pair<MotionPath>(subthis, subthat); }
	 */

	// assume the entire path is length 1
	// NOT Pretty Sure
	/*
	 * public Graph getRealizationAt(double position) { assert (position <= 1 &&
	 * position >= 0); double posLength = position * totalLength;
	 * 
	 * int index = Arrays.binarySearch(nodeSearchMap, posLength); if (index < 0)
	 * { index = -(index + 1); } double remainder; if (index > 0) remainder =
	 * posLength - nodeSearchMap[index - 1]; else remainder = posLength; //
	 * System.out.println("path: " + this); //
	 * System.out.print("nodeSearchMap: "); // for (double d : nodeSearchMap) //
	 * System.out.print(d + ","); // System.out.println(); //
	 * System.out.println("posLength: " + posLength + " index: " + index // +
	 * " reaminder: " + remainder);
	 * 
	 * Node node = nodes.get(index); double cayley =
	 * node.getCayleyAt(remainder); return t.tryRealize(cayley,
	 * node.getSolutionType()); }
	 */

	/*
	 * public void addSingletonNode(double l1, double l2, Interval interval,
	 * OrientedCCS o) { assert (!ended); assert (nodes.isEmpty()); ended = true;
	 * SingletonNode node = new SingletonNode(l1, l2, interval, o);
	 * nodes.add(node);
	 * 
	 * double l = node.getLength(); totalLength = l; tempNodeSearchMap.add(l);
	 * nodeSearchMap = new double[tempNodeSearchMap.size()]; for (int i = 0; i <
	 * nodeSearchMap.length; ++i) nodeSearchMap[i] = tempNodeSearchMap.get(i);
	 * tempNodeSearchMap = null; }
	 * 
	 * public void addEndNode(double lf, Interval interval, OrientedCCS o, int
	 * direction) { assert (!ended); boolean isStart = (nodes.isEmpty() ? true :
	 * false); EndNode node = new EndNode(lf, interval, o, direction, isStart);
	 * nodes.add(node);
	 * 
	 * double l = node.getLength(); totalLength += l;
	 * tempNodeSearchMap.add(totalLength);
	 * 
	 * if (nodes.size() > 1) { ended = true;
	 * 
	 * nodeSearchMap = new double[tempNodeSearchMap.size()]; for (int i = 0; i <
	 * nodeSearchMap.length; ++i) nodeSearchMap[i] = tempNodeSearchMap.get(i);
	 * tempNodeSearchMap = null; } }
	 * 
	 * public void addNode(Interval interval, OrientedCCS o, int direction) {
	 * assert (!ended); Node node = new Node(interval, o, direction);
	 * nodes.add(node);
	 * 
	 * double l = node.getLength(); totalLength += l;
	 * tempNodeSearchMap.add(totalLength); }
	 */

	/*
	 * 
	 * public Graph getEndRealization(int which) { final int start = 0, end = 1;
	 * 
	 * Node n; if (nodes.size() == 1) { n = (SingletonNode) nodes.get(0); //
	 * double l; // if (which == start) { // l = n.getStart();// (n.direction ==
	 * 1 ? n.l1 : n.l2); // } else // l = n.getEnd();// (n.direction == 1 ? n.l2
	 * : n.l2); // ListGraph g = t.tryRealize(l, n.getSolutionType()); // return
	 * g; } else n = (which == start ? startNode() : endNode()); double l =
	 * (which == start ? n.getStart() : n.getEnd()); Graph g = t.tryRealize(l,
	 * n.getSolutionType()); return g; }
	 */

	public static ContinuousMotionPath findPath(TDLinkage t, Realization g1, Realization g2) {
		Edge baseNonEdge = t.getBaseNonedge();
		double l1 = g1.length(baseNonEdge), l2 = g2.length(baseNonEdge);
		RealizationType startType = t.getForwardSolutionType(g1);
		RealizationType targetType = t.getForwardSolutionType(g2);
		OrientedCayleyConfigSpace startO = t.cayleyConfigSpace.getOrientedCCS(startType);
		OrientedCayleyConfigSpace targetO = t.cayleyConfigSpace.getOrientedCCS(targetType);
		Interval startI = startO.getContainingInterval(l1);
		Interval targetI = targetO.getContainingInterval(l2);

		ConnectedComponent component = ConnectedComponent
				.findComponent(t, g1);
		int endIndex = -1;
		for (int i = 0; i < component.size(); ++i) {
			OrientedInterval n = component.get(i);
			if (n.contains(l2, targetType)) {
				endIndex = i;
				break;
			}
		}
		if (endIndex >= 0) {
			List<OrientedInterval> nodesOnPath = component.subList(0, endIndex + 1);
			ContinuousMotionPath p = new ContinuousMotionPath(t);
			p.addAll(nodesOnPath);

			OrientedInterval start = p.startNode();
			OrientedInterval newStart = new OrientedInterval(l1, start.getEnd(), start.getInterval(),
					start.getSolutionType());
			p.set(0, newStart);
			OrientedInterval end = p.endNode();
			OrientedInterval newEnd = new OrientedInterval(end.getStart(), l2, end.getInterval(),
					end.getSolutionType());
			p.set(p.size() - 1, newEnd);

			return p;
		} else {
			return null; // ????
		}

	}

	@Override
	public OrientedInterval startNode() {
		assert (!isEmpty());
		return this.get(0);
	}

	@Override
	public OrientedInterval endNode() {
		assert (!isEmpty());
		return this.get(size() - 1);
	}


	/*
	 * // only start & end need to have lf // inner node: only need (type?), o,
	 * interval, direction, i think
	 * 
	 * class SingletonNode extends EndNode { double l1; double l2;
	 * 
	 * // l1 is the start & l2 is the end. Direction comes from that. // if that
	 * is the case, why do we need direction? public SingletonNode(double l1,
	 * double l2, Interval interval, OrientedCCS o) { super(0, interval, o, (l2
	 * > l1 ? 1 : 0), true); // both start & end this.l1 = l1; this.l2 = l2; }
	 * 
	 * @Override public double getStart() { return l1; }
	 * 
	 * @Override public double getEnd() { return l2; } }
	 * 
	 * class EndNode extends Node { double lf; boolean isStart;
	 * 
	 * public EndNode(double lf, Interval interval, OrientedCCS o, int
	 * direction, boolean isStart) { super(interval, o, direction); this.lf =
	 * lf; this.isStart = isStart; }
	 * 
	 * @Override public double getStart() { if (isStart) return lf; else return
	 * (direction == 1 ? interval.lower : interval.upper); }
	 * 
	 * @Override public double getEnd() { if (!isStart) return lf; else return
	 * (direction == 0 ? interval.lower : interval.upper); }
	 * 
	 * }
	 * 
	 * class Node { // containing space & interval Interval interval;
	 * OrientedCCS o;
	 * 
	 * // 0: next is lower; +1: next is upper; int direction;
	 * 
	 * public double getStart() { return (direction == 1 ? interval.lower :
	 * interval.upper); }
	 * 
	 * public double getEnd() { return (direction == 1 ? interval.upper :
	 * interval.lower); }
	 * 
	 * public SolutionType getSolutionType() { return o.forwardSolutionType; }
	 * 
	 * public Node(Interval interval, OrientedCCS o, int direction) {
	 * this.interval = interval; this.o = o; this.direction = direction; }
	 * 
	 * public String toString() { return getStart() + "~" + getEnd() + "\t" +
	 * getSolutionType().toString() + "@" + interval + "->" + direction; }
	 * 
	 * public double getLength() { return Math.abs(getEnd() - getStart()); }
	 * 
	 * public double getCayleyAt(double position) { assert (position >= 0 &&
	 * position <= getLength()); double lf = (direction == 0 ? getStart() -
	 * position : getStart() + position); // double lf = (direction == 0 ?
	 * interval.lower + position // : interval.upper - position); return lf; } }
	 */

}
