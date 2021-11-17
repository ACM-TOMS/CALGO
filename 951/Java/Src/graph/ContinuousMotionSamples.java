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
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;

public class ContinuousMotionSamples<E> implements List<SamplePoint<E>> {
	// private int samplePointDimension;

	private ContinuousMotion motion;

	private ArrayList<SamplePoint<E>> samples;

	public ContinuousMotionSamples( ContinuousMotion motion) {
		// samplePointDimension = n;
		this.motion = motion;
		samples = new ArrayList<SamplePoint<E>>();
	}

	public ContinuousMotion getMotion() {
		return motion;
	}
	
	public List<SamplePoint<E>> getList(){
		return this;
	}

	/*
	 * public int getDimension() { return samplePointDimension; }
	 */

	public boolean add(SamplePoint<E> arg0) {
		return samples.add(arg0);
	}

	public void add(int arg0, SamplePoint<E> arg1) {
		samples.add(arg0, arg1);
	}

	public boolean addAll(Collection<? extends SamplePoint<E>> arg0) {
		return samples.addAll(arg0);
	}

	public boolean addAll(int arg0, Collection<? extends SamplePoint<E>> arg1) {
		return samples.addAll(arg0, arg1);
	}

	public void clear() {
		samples.clear();
	}

	public boolean contains(Object arg0) {
		return samples.contains(arg0);
	}

	public boolean containsAll(Collection<?> arg0) {
		return samples.containsAll(arg0);
	}

	public SamplePoint<E> get(int arg0) {
		return samples.get(arg0);
	}

	public int indexOf(Object arg0) {
		return samples.indexOf(arg0);
	}

	public boolean isEmpty() {
		return samples.isEmpty();
	}

	public Iterator<SamplePoint<E>> iterator() {
		return samples.iterator();
	}

	public int lastIndexOf(Object arg0) {
		return samples.lastIndexOf(arg0);
	}

	public ListIterator<SamplePoint<E>> listIterator() {
		return samples.listIterator();
	}

	public ListIterator<SamplePoint<E>> listIterator(int arg0) {
		return samples.listIterator(arg0);
	}

	public boolean remove(Object arg0) {
		return samples.remove(arg0);
	}

	public SamplePoint<E> remove(int arg0) {
		return samples.remove(arg0);
	}

	public boolean removeAll(Collection<?> arg0) {
		return samples.removeAll(arg0);
	}

	public boolean retainAll(Collection<?> arg0) {
		return samples.retainAll(arg0);
	}

	public SamplePoint<E> set(int arg0, SamplePoint<E> arg1) {
		return samples.set(arg0, arg1);
	}

	public int size() {
		return samples.size();
	}

	public List<SamplePoint<E>> subList(int arg0, int arg1) {
		return samples.subList(arg0, arg1);
	}

	public Object[] toArray() {
		return samples.toArray();
	}

	public <T> T[] toArray(T[] arg0) {
		return samples.toArray(arg0);
	}
}
