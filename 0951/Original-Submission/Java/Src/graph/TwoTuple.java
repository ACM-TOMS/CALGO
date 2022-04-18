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

public class TwoTuple<T, V> {
	public T o1;
	public V o2;

	public TwoTuple(T t, V v) {
		o1 = t;
		o2 = v;
	}

	public T getFirst() {
		return o1;
	}

	public V getSecond() {
		return o2;
	}

	public void setFirst(T o) {
		o1 = o;
	}

	public void setSecond(V o) {
		o2 = o;
	}

	static boolean same(Object o1, Object o2) {
		return o1 == null ? o2 == null : o1.equals(o2);
	}

	@SuppressWarnings("unchecked")
	public boolean equals(Object obj) {
		if (!(obj instanceof TwoTuple<?, ?>))
			return false;
		TwoTuple<T, V> p = (TwoTuple<T, V>) obj;
		// order matters
		return same(p.o1, this.o1) && same(p.o2, this.o2);
	}

	public String toString() {
		return "{" + o1 + ", " + o2 + "}";
	}

}
