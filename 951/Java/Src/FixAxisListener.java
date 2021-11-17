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
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

import ccs.graph.Point2D;
import ccs.graph.Vertex;

public class FixAxisListener implements ActionListener, MouseListener {
	private FixAxisListener() {
	}

	private static FixAxisListener me = new FixAxisListener();

	public static FixAxisListener getInstance() {
		return me;
	}

	private boolean isPickingAxis = false;

	private Vertex o, x;

	public Vertex getO() {
		return o;
	}

	private static final int clickPrecision = 15;

	public void mouseClicked(MouseEvent e) {
		if (!isPickingAxis)
			return;

		Point2D clickP = new Point2D(e.getX(), e.getY());
		TDLinkageModel tdModel = TDLinkageModel.getInstance();

		if (o == null) {
			o = tdModel.getVertex(clickP, clickPrecision);
		} else {
			x = tdModel.getVertex(clickP, clickPrecision);
			tdModel.setAxis(o, x);
			o = null;
			x = null;
			isPickingAxis = false;
		}
		GPanel.getInstance().repaint();
	}

	public void mouseEntered(MouseEvent arg0) {
		// TODO Auto-generated method stub

	}

	public void mouseExited(MouseEvent arg0) {
		// TODO Auto-generated method stub

	}

	public void mousePressed(MouseEvent arg0) {
		// TODO Auto-generated method stub

	}

	public void mouseReleased(MouseEvent arg0) {
		// TODO Auto-generated method stub

	}

	public void actionPerformed(ActionEvent e) {
		// control.fixAxisButton is clicked
		isPickingAxis = true;
	}

}
