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

public class Pair<T> {
	// not ordered
	private T o1;
	private T o2;

	public Pair(T o1, T o2) {
		this.o1 = o1;
		this.o2 = o2;
	}

	public T getFirst() {
		return o1;
	}

	public T getSecond() {
		return o2;
	}

	public void setFirst(T o) {
		o1 = o;
	}

	public void setSecond(T o) {
		o2 = o;
	}

	static boolean same(Object o1, Object o2) {
		return o1 == null ? o2 == null : o1.equals(o2);
	}

	@SuppressWarnings("unchecked")
	public boolean equals(Object obj) {
		if (!(obj instanceof Pair<?>))
			return false;
		Pair<T> p = (Pair<T>) obj;
		// order does not matter
		return same(p.o1, this.o1) && same(p.o2, this.o2)
				|| same(p.o1, this.o2) && same(p.o2, this.o1);
	}

	public String toString() {
		return "Pair{" + o1 + ", " + o2 + "}";
	}

	public void swap() {
		T temp = o1;
		o1 = o2;
		o2 = temp;
	}
}

