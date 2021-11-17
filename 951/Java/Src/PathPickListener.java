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

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.AbstractButton;

import ccs.graph.Realization;

public class PathPickListener implements ActionListener {

	private ControlPanel control = ControlPanel.getInstance();

	private TDLinkageModel tdModel = TDLinkageModel.getInstance();
	private CCSModel ccsModel = CCSModel.getInstance();
	private MotionModel motionModel = MotionModel.getInstance();

	public void actionPerformed(ActionEvent e) {
		AbstractButton b = (AbstractButton) e.getSource();
		int picking = control.getPickingEnd(b);
		assert (picking == 1 || picking == 2);

		// assert (status.isPicking());
		if (control.getPathType() == 2) {
			Realization real = tdModel.getGraphClone();
			if (picking == 1)
				motionModel.setStartRealization(real);
			else
				motionModel.setEndRealization(real);
			// main.graphPanel.repaint();
		}
		// } else if (status == PathStatus.cayleyPicking) {
		double val = ccsModel.getCurrent();
		if (picking == 1)
			motionModel.setStartCayleyConfig(val);
		else
			motionModel.setEndCayleyConfig(val);

		control.checkGenPathBtn();
		
		GPanel.getInstance().repaint();
		CCSPanel.getInstance().repaint();
		
		// main.ccsInfo.ccsPanel.repaint();
	}

}
