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

import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;

import ccs.Debug;

/**
 * Interface for MotionPath and ConnectedComponent
 */
public abstract class ContinuousMotion implements List<OrientedInterval> {
	protected TDLinkage tdLinkage;

	protected ArrayList<OrientedInterval> orientedIntervals;

	public ContinuousMotion(TDLinkage t) {
		this.tdLinkage = t;
		orientedIntervals = new ArrayList<OrientedInterval>();
	}

	public abstract OrientedInterval startNode();

	public abstract OrientedInterval endNode();

	public AbstractList<RealizationType> getSolutionTypes() {
		ArrayList<RealizationType> l = new ArrayList<RealizationType>();
		for (OrientedInterval n : this) {
			l.add(n.getSolutionType());
		}
		return l;
	}

	public String toString() {
		return orientedIntervals.toString();
	}

	/*
	 * public AbstractList<Node> getNodes() { // return nodes;
	 */

	public int size() {
		return orientedIntervals.size();
	}

	public OrientedInterval set(int index, OrientedInterval n) {
		return orientedIntervals.set(index, n);
	}

	public boolean add(OrientedInterval n) {
		// Debug.msg("adding node " + n);
		return orientedIntervals.add(n);
	}

	public OrientedInterval get(int index) {
		return orientedIntervals.get(index);
	}

	public void clear() {
		orientedIntervals.clear();
	}

	public boolean contains(Object arg0) {
		return orientedIntervals.contains(arg0);
	}

	public boolean containsAll(Collection<?> arg0) {
		return orientedIntervals.containsAll(arg0);
	}

	public int indexOf(Object arg0) {
		return orientedIntervals.indexOf(arg0);
	}

	public boolean isEmpty() {
		return orientedIntervals.isEmpty();
	}

	public Iterator<OrientedInterval> iterator() {
		return orientedIntervals.iterator();
	}

	/**
	 * Not supported.
	 * 
	 * @throws UnsupportedOperationException
	 */
	public int lastIndexOf(Object arg0) {
		throw new UnsupportedOperationException();
	}

	public ListIterator<OrientedInterval> listIterator() {
		return orientedIntervals.listIterator();
	}

	public ListIterator<OrientedInterval> listIterator(int arg0) {
		return orientedIntervals.listIterator(arg0);
	}

	/**
	 * @param fromIndex
	 *            included
	 * @param toIndex
	 *            not included
	 * @see java.util.List#subList(int, int)
	 */
	public List<OrientedInterval> subList(int fromIndex, int toIndex) {
		return orientedIntervals.subList(fromIndex, toIndex);
	}

	public Object[] toArray() {
		return orientedIntervals.toArray();
	}

	public <T> T[] toArray(T[] arg0) {
		return orientedIntervals.toArray(arg0);
	}

	public void add(int arg0, OrientedInterval arg1) {
		throw new UnsupportedOperationException();
	}

	public boolean addAll(Collection<? extends OrientedInterval> arg0) {
		for (OrientedInterval n : arg0)
			if (!add(n))
				return false;
		return true;
	}

	public boolean addAll(int arg0, Collection<? extends OrientedInterval> arg1) {
		throw new UnsupportedOperationException();
	}

	public boolean remove(Object arg0) {
		throw new UnsupportedOperationException();
	}

	public OrientedInterval remove(int arg0) {
		throw new UnsupportedOperationException();
	}

	public boolean removeAll(Collection<?> arg0) {
		throw new UnsupportedOperationException();
	}

	public boolean retainAll(Collection<?> arg0) {
		throw new UnsupportedOperationException();
	}

	/**
	 * Using a NodeSampler to sample every node.
	 * 
	 * @param <E>
	 * 
	 * @param s
	 *            a sampler defining sampling method for each node
	 * @return a sequence of sample points obtained from each node
	 */
	public <E> ContinuousMotionSamples<E> getRealizations(NodeSampler<E> s) {
		ContinuousMotionSamples<E> samples = new ContinuousMotionSamples<E>(
				this);
		for (OrientedInterval n : this) {
			samples.addAll(s.sample(n));
		}
		return samples;
	}

	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof ContinuousMotion))
			return false;
		ContinuousMotion that = (ContinuousMotion) obj;
		if (size() != that.size())
			return false;
		ListIterator<OrientedInterval> iterThis = this.listIterator(), iterThat = that
				.listIterator();
		while (iterThis.hasNext()) {
			if (!iterThis.next().equals(iterThat.next()))
				return false;
		}
		return true;
	}

}

/**
 * Contains one consecutive interval on a continuous motion.
 * 
 */
class OrientedInterval {
	/**
	 * Start Cayley configuration
	 */
	private double start;

	/**
	 * End Cayley configuration
	 */
	private double end;

	private Interval interval;

	private RealizationType realizationType;

	public OrientedInterval(double start, double end, Interval interval, RealizationType type) {
		this.start = start;
		this.end = end;
		this.interval = interval;
		this.realizationType = type.clone();
	}

	public double getStart() {
		return start;
	}

	public double getEnd() {
		return end;
	}

	public Interval getInterval() {
		return interval;
	}

	public RealizationType getSolutionType() {
		return realizationType;
	}

	public double getLength() {
		return Math.abs(getStart() - getEnd());
	}

	public boolean contains(double lf, RealizationType type) {
		if (!type.equals(getSolutionType()))
			return false;
		return Util.between(lf, getStart(), getEnd());
	}

	protected double sampleCayleyAt(double percentage) {
		assert (percentage <= 1 && percentage >= 0);
		double l = getLength() * percentage;
		return (getStart() < getEnd() ? getStart() + l : getStart() - l);
	}

	@Override
	public String toString() {
		return getStart() + "~" + getEnd() + "\t"
				+ getSolutionType().toString();
	}

	@Override
	public boolean equals(Object obj) {
		// TODO Auto-generated method stub
		if (!(obj instanceof OrientedInterval))
			return false;
		OrientedInterval that = (OrientedInterval) obj;
		return start == that.start && end == that.end
				&& interval.equals(that.interval) && realizationType.equals(that.realizationType);
	}

}
