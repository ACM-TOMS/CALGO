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

import ccs.graph.Edge;

/**
 * Organizing floating panels
 */
public class FloatingPanels {
	private FloatingPanels() {
		curve3Dpanel = new Curve3DContainer();
		curve3Dpanel.setVisible(false);

	}

	private static FloatingPanels me = new FloatingPanels();

	public static FloatingPanels getInstance() {
		return me;
	}

	private Curve3DContainer curve3Dpanel;
	private EdgeLengthContainer edgeLengthPanel;

	public boolean addCurve3D(Curve3DMotion c) {
		return curve3Dpanel.addCurve(c);
	}

	public boolean removeCurve3D(Curve3D c) {
		return curve3Dpanel.removeCurve(c);
	}

	public void clearCurves() {
		curve3Dpanel.clear();
	}

	public void showCurve3DPanel() {
		curve3Dpanel.setVisible(true);
		curve3Dpanel.setAlwaysOnTop(true);
		curve3Dpanel.createDialog();
	}

	public void closeCurve3DPanel() {
		curve3Dpanel.closeDialog();
		curve3Dpanel.setVisible(false);
	}

	public void repaintAll() {
		if (curve3Dpanel.isVisible()) {
			// Debug.warnMsg("painting 3d panel");
			curve3Dpanel.repaint();
		}

		if (edgeLengthPanel.isVisible())
			edgeLengthPanel.repaint();
	}

	public Edge[] getProjectionEdges() {
		return edgeLengthPanel.getProjectionEdges();
	}

	public void showEdgeLengthPanel() {
		edgeLengthPanel.createDialog();
		edgeLengthPanel.setAlwaysOnTop(true);
	}

	public void closeEdgeLengthPanel() {
		edgeLengthPanel.closeDialog();
	}

	public int indexOf(Edge e) {
		return edgeLengthPanel.indexOf(e);
	}

	public void refresh() {
		if (TDLinkageModel.getInstance().isLow()) {
			edgeLengthPanel = new EdgeLengthContainer(TDLinkageModel
					.getInstance().getTd());
			edgeLengthPanel.setVisible(false);
		}
		// edgeLengthPanel.refresh(TreeDecompModel.getInstance().getTd());
	}
}
