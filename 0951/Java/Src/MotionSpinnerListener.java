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

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import ccs.graph.ConnectedComponent;
import ccs.graph.ContinuousMotion;
import ccs.graph.Realization;
import ccs.graph.SamplePoint;

public class MotionSpinnerListener implements ChangeListener {

	/**
	 * 0: pathSpinner; 1:p1; 2:p2
	 */
	private int which = 0;

	/**
	 * @param n
	 *            Indicates which spinner. 0: pathSpinner; 1:p1; 2:p2
	 */
	public MotionSpinnerListener(int n) {
		which = n;
	}

	private TDLinkageModel tdModel = TDLinkageModel.getInstance();
	private MotionModel motionModel = MotionModel.getInstance();

	public void stateChanged(ChangeEvent e) {
		DialSpinner<SamplePoint<Realization>> spinner = ((DialSpinner<SamplePoint<Realization>>) (e
				.getSource()));

		if (!spinner.isVisible())
			return;

		ContinuousMotion p;
		if (which == 0)
			p = motionModel.getMotion();
		else if (which == 1)
			p = motionModel.getM1();
		else
			p = motionModel.getM2();
		// double percentage = spinner.getPercentage();

		// ContinuousMotionSamples<Graph> pSamples = motionModel
		// .getMotionSamples();

		// int index = (int) percentage * pSamples.size();
		// SamplePoint<Graph> point = pSamples.get(index);
		SamplePoint<Realization> point = spinner.getValue();
		Realization g = point.getValue();

		if (g != null) {
			switch (which) {
			case 0:
				tdModel.setPoints(g);
				break;
			/*
			 * case 1: motionModel.getStartRealization().setPoints(g); break;
			 * case 2: motionModel.getEndRealization().setPoints(g); break;
			 */
			}

		}

		GPanel.getInstance().repaint();
		CCSModel.getInstance().refreshCurrent();
		CCSPanel.getInstance().repaint();
		FloatingPanels.getInstance().repaintAll();
	}

}
