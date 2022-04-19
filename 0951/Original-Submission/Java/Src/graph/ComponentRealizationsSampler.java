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

public class ComponentRealizationsSampler implements NodeSampler<Realization> {
	private int samplesPerNode;
	private TDLinkage t;

	public ComponentRealizationsSampler(TDLinkage t, int spn) {
		assert (spn > 0);
		this.t = t;
		samplesPerNode = spn;
	}
	
	final double ACCURACY = 2;

	public AbstractList<SamplePoint<Realization>> sample(OrientedInterval n) {
		ArrayList<SamplePoint<Realization>> list = new ArrayList<SamplePoint<Realization>>();
		
		int spn = Math.max(samplesPerNode, (int)(n.getLength()/ACCURACY));

		for (int i = 0; i < spn; ++i) {
			double percentage = (double) i / spn;
			double cayley = n.sampleCayleyAt(percentage);

			Realization g = t.tryRealize(cayley, n.getSolutionType());
			assert (g != null);
			/*// ad-hoc modification
			double pp = percentage;
			while(g==null){
				pp -= 0.0001;
				cayley = n.sampleCayleyAt(pp);
				g = t.tryRealize(cayley,n.getSolutionType());
			}*/
			list.add(new SamplePoint<Realization>(g));
			
		}

		return list;
	}

}
