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

import java.util.ArrayList;
import java.util.Collection;

import ccs.graph.CayleyConfigSpace;
import ccs.graph.CCSInterface;
import ccs.graph.Interval;
import ccs.graph.OrientedCayleyConfigSpace;
import ccs.graph.Pair;
import ccs.graph.RealizationType;

public class CCSModel {
	private static CCSModel me;

	// private TreeDecompModel tdModel;
	// private ControlPanel control;

	private CCSModel() {
		typesToDraw = new ArrayList<RealizationType>();
	}

	public static CCSModel getInstance() {
		if (me == null)
			me = new CCSModel();
		return me;
	}

	private CayleyConfigSpace ccs;
	private OrientedCayleyConfigSpace occs;
	private boolean isOriented;
	private double curVal = -1;

	private ArrayList<RealizationType> typesToDraw; // o-ccs along a path

	ArrayList<ConfigPanel> markedPoints = new ArrayList<ConfigPanel>();

	public Collection<ConfigPanel> getMarkedPoints() {
		return markedPoints;
	}

	/**
	 * Set current ccs, oriented/nonoriented according to control specification.
	 */
	public void refresh() {
		TDLinkageModel tdModel = TDLinkageModel.getInstance();
		ControlPanel control = ControlPanel.getInstance();

		CayleyConfigSpace nccs = tdModel.getCayleyConfigurationSpace();
		if (control.isOrientedCCS()) {
			isOriented = true;
			occs = nccs.getOrientedCCS(tdModel.getForwardSolutionType());
		} else {
			isOriented = false;
			occs = null;
		}
		ccs = nccs;
		typesToDraw.clear();
	}

	public Pair<CCSInterface> getCCS() {
		return new Pair<CCSInterface>(ccs, occs);
	}

	public RealizationType getSolutionType() {
		if (isOriented)
			return occs.getSolutionType();
		else
			return null;
	}

	public Collection<RealizationType> getTypesToDraw() {
		return typesToDraw;
	}

	/**
	 * 
	 * @return the interval in ccs containing val. If val is not in current ccs,
	 *         return null
	 */
	public Interval contains(double val) {
		ArrayList<Interval> union;
		if (isOriented)
			union = occs.getIntervals();
		else
			union = ccs.getIntervals();
		for (Interval i : union)
			if (i.contains(val))
				return i;
		return null;
	}

	public boolean isGenerated() {
		return ccs != null;
	}

	public boolean isOriented() {
		return isOriented;
	}

	public double getCurrent() {
		return curVal;
	}

	public void setCurrent(double val) {
		assert (contains(val) != null);
		curVal = val;
	}

	/**
	 * Reset current Cayley config from current realization
	 */
	public void refreshCurrent() {
		double length = TDLinkageModel.getInstance().getBaseNoneedgeLength();
		setCurrent(length);
	}

	public double getMin() {
		if (isOriented)
			return occs.getMin();
		else
			return ccs.getMin();
	}

	public double getMax() {
		if (isOriented)
			return occs.getMax();
		else
			return ccs.getMax();
	}

	public void addMarkedPoint(ConfigPanel c) {
		if (!markedPoints.contains(c)) {
			markedPoints.add(c);
		}
	}

	public void removeMarkedPoint(ConfigPanel c) {
		if (markedPoints.contains(c)) {
			markedPoints.remove(c);
		}
	}

	public void addType(RealizationType t) {
		//assert (!isOriented);
		assert (((CayleyConfigSpace) ccs).getOrientedCCS(t) != null);
		if (!typesToDraw.contains(t))
			typesToDraw.add(t);
	}

	public void removeType(RealizationType t) {
		typesToDraw.remove(t);
	}

	public void clearTypes() {
		typesToDraw.clear();
	}

	/*
	 * Refresh CCS according to current orientation. Only necessary if
	 * isOriented ??
	 * 
	 * // public void refreshCCS() { // init(main.td.ccs); // // if (isOriented)
	 * { // // SolutionType type = tdModel.getForwardSolutionType(); // //
	 * addType(type); // initCCS(); // } // repaint(); // }
	 */

}
