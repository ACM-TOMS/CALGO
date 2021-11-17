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

import javax.swing.JSpinner;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import ccs.graph.Realization;
import ccs.graph.RealizationType;

public class CCSSpinnerListener implements ChangeListener {

	TDLinkageModel tdModel = TDLinkageModel.getInstance();
	CCSModel ccsModel = CCSModel.getInstance();

	public void stateChanged(ChangeEvent e) {
		JSpinner s = ((JSpinner) (e.getSource()));
		if (!s.isVisible())
			return;
		Double val = (Double) s.getValue();

		// only realize those displayed on ccsPanel
		if (ccsModel.contains(val) == null)
			return;

		ccsModel.setCurrent(val);

		RealizationType forward = tdModel.getForwardSolutionType();
		Realization g = tdModel.tryRealize(val, forward);
		if (g != null)
			tdModel.setPoints(g);
		else if (!ccsModel.isOriented()) {
			ArrayList<Realization> gs = tdModel.tryRealize(val);
			if (!gs.isEmpty())
				tdModel.setPoints(gs.get(0));
		}
		// main.graphPanel.repaint();
		// TODO: ccs Panel refresh
		
		GPanel.getInstance().repaint();
		FloatingPanels.getInstance().repaintAll();
		CCSPanel.getInstance().repaint();
	}

}
