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

import java.awt.Color;
import java.util.ArrayList;

public class RealizationType {
	int orientations[];

	public RealizationType(int size) {
		orientations = new int[size];
	}

	public int getSize() {
		return orientations.length;
	}

	public int getOrientation(int index) {
		return orientations[index];
	}

	public void setOrientation(int index, int value) {
		assert (value == 1 || value == -1 || value == 0);
		orientations[index] = value;
	}

	@Override
	public RealizationType clone() {
		RealizationType t = new RealizationType(this.getSize());
		t.orientations = orientations.clone();
		return t;
	}

	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof RealizationType))
			return false;

		RealizationType that = (RealizationType) obj;
		if (that.getSize() != this.getSize())
			return false;
		for (int i = 0; i < getSize(); ++i) {
			int o1 = that.getOrientation(i), o2 = this.getOrientation(i);
			if (o1 != o2)
				return false;
		}
		return true;
	}

	@Override
	public String toString() {
		String s = "";
		for (int i : orientations) {
			s += i + " ";
		}
		return s;
	}

	public boolean compatible(RealizationType that) {
		if (this.equals(that))
			return true;
		for (int i = 0; i < Math.min(this.getSize(), that.getSize()); ++i) {
			int o1 = that.getOrientation(i), o2 = this.getOrientation(i);
			if (o1 != 0 && o2 != 0 && o1 != o2)
				return false;
		}
		return true;
	}

	// ??? is it the case that every solution type realizable for some |f|?
	public static ArrayList<RealizationType> generateSolutionTypes(int length) {
		ArrayList<RealizationType> solutionTypes = new ArrayList<RealizationType>();
		double upper = Math.pow(2, length);
		for (int i = 0; i < upper; ++i) {
			RealizationType type = new RealizationType(length);
			for (int j = 0; j < length; ++j)
				type.setOrientation(j, -1);
			String bin = Integer.toBinaryString(i);
			// System.out.println(bin);
			for (int j = 0; j < bin.length(); ++j) {
				if (bin.charAt(bin.length() - j - 1) == '1')
					type.setOrientation(j, 1);
			}
			solutionTypes.add(type);
		}
		return solutionTypes;
	}

	protected boolean isExtreme() {
		int count = 0;
		for (int o : orientations) {
			if (o == 0)
				count++;
		}
		return count >= 1;
	}

	protected int indexOfZero() {
		for (int i = 0; i < getSize(); ++i)
			if (getOrientation(i) == 0)
				return i;
		return -1;
	}

	public void flipOrientation(int index) {
		int o = orientations[index];
		orientations[index] = -o;
	}

	public int getEncoding() {
		int code = 0;
		for (int o : orientations) {
			switch (o) {
			case 1:
				code = code * 3 + 1;
				break;
			case -1:
				code = code * 3 + 2;
				break;
			case 0:
				code = code * 3 + 0;
				break;
			}
		}
		return code;
	}

	public Color getColor() {
		Color c = Util.randomColor(getEncoding()).darker();
		return c;
	}
}
