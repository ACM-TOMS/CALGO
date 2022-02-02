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
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import ccs.Debug;

public class CayleyConfigSpace implements CCSInterface {
	// final static int width = 1000, height = 180;

	private ArrayList<RealizationType> forwardSolutionTypes;

	// decomposed, indexed by forward solution type
	private ArrayList<OrientedCayleyConfigSpace> orientedCCS;

	private ArrayList<Interval> unionIntervals;
	
	private int numOrientedIntervals = 0;

	private double min = Double.POSITIVE_INFINITY,
			max = Double.NEGATIVE_INFINITY;

	public CayleyConfigSpace() {
		forwardSolutionTypes = new ArrayList<RealizationType>();
		orientedCCS = new ArrayList<OrientedCayleyConfigSpace>();
	}

	public Collection<OrientedCayleyConfigSpace> getOrientedCCSs() {
		return orientedCCS;
	}

	public OrientedCayleyConfigSpace getOrientedCCS(RealizationType type) {
		int index;
		for (index = 0; index < forwardSolutionTypes.size(); ++index) {
			RealizationType t = forwardSolutionTypes.get(index);
			if (t.equals(type))
				break;
		}
		if (index >= forwardSolutionTypes.size())
			return null;
		// int index = forwardSolutionTypes.indexOf(type);
		return orientedCCS.get(index);
	}

	public void addOrientedCCS(OrientedCayleyConfigSpace o) {
		if (Realization.DEBUG)
			for (RealizationType t : forwardSolutionTypes)
				assert (!t.equals(o.forwardSolutionType));
		forwardSolutionTypes.add(o.forwardSolutionType);
		orientedCCS.add(o);
		double min = o.getMin(), max = o.getMax();
		this.min = (min < this.min ? min : this.min);
		this.max = (max > this.max ? max : this.max);
		numOrientedIntervals += o.intervals.size();
	}

	public int getNumOfSolutionTypes() {
		return forwardSolutionTypes.size();
	}

	public double getMin() {
		return min;
	}

	public double getMax() {
		return max;
	}

	@Override
	public String toString() {
		String s = "";
		for (int i = 0; i < forwardSolutionTypes.size(); ++i) {
			s += forwardSolutionTypes.get(i).toString() + ":\n";
			s += orientedCCS.get(i).toString() + "\n";
		}
		return s;
	}

	public void generateUnionIntervals() {
		unionIntervals = new ArrayList<Interval>();
		ArrayList<Interval> l = new ArrayList<Interval>();
		for (int i = 0; i < forwardSolutionTypes.size(); ++i) {
			l.addAll(orientedCCS.get(i).intervals);
		}
		if (l.isEmpty())
			return;
		// System.out.println(l);
		Collections.sort(l);
		
		Debug.msg("intervals to union: " + l);
		unionIntervals.add(l.get(0));
		for (int i = 0; i < l.size(); ++i) {
			int lastIndex = unionIntervals.size() - 1;
			Interval cur = l.get(i);
			Interval last = unionIntervals.get(lastIndex);
			if (cur.hasIntersection(last)) {
				unionIntervals.set(lastIndex, cur.union(last));
			} else {
				unionIntervals.add(cur);
			}
		}
	}

	public ArrayList<Interval> getIntervals() {
		return unionIntervals;
	}

	public String printUnionIntervals() {
		return unionIntervals.toString();
	}

	public static CayleyConfigSpace generateCCS(TDLinkage td) {
		boolean b = td.isLow();
		assert (b);
		b = td.constructionSequenceGenerated();
		assert (b);
		// for low: reverse solve for each extreme graph
		// TODO: generate oriented for each fwd sol type
		// ArrayList<java.lang.Double> candidate = new
		// ArrayList<java.lang.Double>();
		// SolutionType forward = getForwardSolutionType();
		HashMap<RealizationType, ArrayList<Double>> candidates = new HashMap<RealizationType, ArrayList<Double>>();
		ArrayList<RealizationType> types = RealizationType.generateSolutionTypes(td
				.getNumOfConstructStep());
		for (RealizationType type : types) {
			candidates.put(type, new ArrayList<Double>());
		}

		for (int i = 0; i < td.getNumOfConstructStep(); ++i) {
			TDLinkage t = td.getExtremeGraph(i);
			// Debug.warnMsg(i+"th extreme graph:"+t);
			Interval interval = td.extremeInterval(td.getConstructionStep(i));
			Debug.msg(i + "th extreme Interval:" + interval);

			for (double length : interval.toArray()) {
				// added Accuracy ...
				// ??? May return duplicate result, Why? -- come from
				// reflection?
				ArrayList<Realization> solutions = t.tryRealize(length);
				// Debug.warnMsg("extreme l:"+length+"; extreme solution size:"+solutions.size());

				for (Realization g : solutions) {
					// System.out.println(g.printPoints());
					double distance = g.length(td.getBaseNonedge());
					RealizationType extremeSolType = td.getForwardSolutionType(g,
							i);
					for (RealizationType type : types) {
						if (type.compatible(extremeSolType)) {
							ArrayList<Double> typeCandidates = candidates
									.get(type);
							boolean addit = true;
							for (double d : typeCandidates) {
								if (Math.abs(d - distance) < TDLinkage.ACCURACY) {
									addit = false;
									break;
								}
							}
							if (addit) {
								Debug.msg("for type " + type
										+ ", add solution of " + length
										+ " to candidates: " + distance);
								typeCandidates.add(distance);
							}
						}
					}
				}
			}

		}
		CayleyConfigSpace ccs = new CayleyConfigSpace();
		for (RealizationType type : types) {
			ArrayList<Double> typeCandidates = candidates.get(type);
			if (typeCandidates.isEmpty())
				continue;
			Collections.sort(typeCandidates);
			Debug.msg("for type " + type + ", sorted candidates:"
					+ typeCandidates);

			// Check candidate points & set up intervals
			OrientedCayleyConfigSpace orientedCCS = td.genOrientedCCSFromCandidates(
					typeCandidates, type);
			if (orientedCCS.isEmpty())
				continue;
			Debug.msg("oriented ccs:" + orientedCCS + "\n");
			ccs.addOrientedCCS(orientedCCS);
		}
		Debug.msg(ccs);
		ccs.generateUnionIntervals();
		Debug.msg("overall: " + ccs.printUnionIntervals(), 1);
		Debug.msg("# of oriented intervals: " + ccs.numOrientedIntervals, 1);
		return ccs;
	}
}
