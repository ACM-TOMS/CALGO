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
import java.util.LinkedList;
import java.util.ListIterator;

import ccs.Debug;
import ccs.TDLinkageModel;

public class ConnectedComponent extends ContinuousMotion {

	public ConnectedComponent(TDLinkage t) {
		super(t);
	}

	@Override
	public OrientedInterval startNode() {
		assert (!isEmpty());
		return get(0);
	}

	@Override
	public OrientedInterval endNode() {
		return startNode();
	}

	public static ConnectedComponent findComponent(TDLinkage t, Realization g1) {
		double startl = g1.length(t.getBaseNonedge());
		RealizationType startType = t.getForwardSolutionType(g1);
		OrientedCayleyConfigSpace startO = t.cayleyConfigSpace.getOrientedCCS(startType);
		Interval startI = startO.getContainingInterval(startl);

		assert (startI != null);

		ConnectedComponent component = new ConnectedComponent(t);

		// Using the entire interval as start node
		double endl = startI.toArray()[0];
		startl = startI.toArray()[1];
		OrientedInterval startNode = new OrientedInterval(startl, endl, startI, startType);

		component.add(startNode);

		RealizationType type = startType;

		while (true) {
			Debug.msg("from end point: " + endl);

			// Flip the orientation of extreme construction step
			Realization g = t.tryRealize(endl, type);
			RealizationType extremeType = t.getForwardSolutionType(g);
			if (!extremeType.isExtreme()) {
				Debug.warnMsg(extremeType + "");
				Debug.warnMsg(endl+"");
				for (int i = 0; i < TDLinkageModel.getInstance()
						.getNumOfConstructStep(); ++i)
					Debug.msg(i
							+ ":"
							+ TDLinkageModel.getInstance()
									.getConstructionStep(i));

			}
			assert (extremeType.isExtreme());
			int zeroIndex = extremeType.indexOfZero();
			type.flipOrientation(zeroIndex);

			OrientedCayleyConfigSpace o = t.cayleyConfigSpace.getOrientedCCS(type);
			assert (o != null);
			Interval curI = o.getContainingInterval(endl);
			assert (curI != null);

			if (o == startO && curI == startI) {
				Debug.msg("BACK TO START. ");
				// Only contains the start node once.
				/* addNode(startl, startI, startO, 0); */
				return component;
			}

			// continue finding
			startl = endl;
			if (Math.abs(endl - curI.lower) < TDLinkage.ACCURACY) {
				endl = curI.upper;
			} else {
				assert (Math.abs(endl - curI.upper) < TDLinkage.ACCURACY);
				endl = curI.lower;
			}
			OrientedInterval n = new OrientedInterval(startl, endl, curI, type);
			component.add(n);
		}
	}

	/**
	 * @return a list of all connected components of t
	 */
	public static ArrayList<ConnectedComponent> findAllComponents(TDLinkage t) {
		ArrayList<ConnectedComponent> list = new ArrayList<ConnectedComponent>();

		// 1. list every intervals from every Occs.
		LinkedList<TwoTuple<RealizationType, Interval>> intervalList = new LinkedList<TwoTuple<RealizationType, Interval>>();
		for (OrientedCayleyConfigSpace oc : t.cayleyConfigSpace.getOrientedCCSs()) {
			RealizationType type = oc.getSolutionType().clone();
			for (Interval i : oc.getIntervals()) {
				intervalList.add(new TwoTuple<RealizationType, Interval>(type, i));
			}
		}

		// 2. for each interval:
		while (!intervalList.isEmpty()) {
			// (1) gen component.
			TwoTuple<RealizationType, Interval> orientedInt = intervalList
					.getFirst();
			RealizationType type = orientedInt.getFirst();
			Interval inv = orientedInt.getSecond();
			double lf = (inv.lower + inv.upper) / 2;
			Realization realization = t.solve(lf, type);
			assert (realization != null);

			Debug.msg("remove first oriented interval:" + orientedInt);

			ConnectedComponent component = findComponent(t, realization);
			list.add(component);

			// (2) remove the intervals passed by the component from the list of
			// intervals
			for (OrientedInterval n : component) {
				TwoTuple<RealizationType, Interval> nTuple = new TwoTuple<RealizationType, Interval>(
						n.getSolutionType().clone(), n.getInterval());
				Debug.msg("interval list:" + intervalList);
				Debug.msg("interval to remove:" + nTuple);
				boolean b = intervalList.remove(nTuple);
				assert (b);
			}
		}
		Debug.msg("total # of components:"+list.size(), 0);
		return list;
	}

	/**
	 * TODO: Equal components can have different order in nodes.
	 */
	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof ContinuousMotion))
			return false;
		ContinuousMotion that = (ContinuousMotion) obj;
		if (size() != that.size())
			return false;

		OrientedInterval thisFirst = this.get(0);
		if (!that.contains(thisFirst))
			return false;
		ListIterator<OrientedInterval> iterThis = this.listIterator(), iterThat = that
				.listIterator(that.indexOf(thisFirst));
		while (iterThis.hasNext()) {
			if (!iterThat.hasNext())
				iterThat = that.listIterator();
			if (!iterThis.next().equals(iterThat.next()))
				return false;
		}
		return true;
	}
}
