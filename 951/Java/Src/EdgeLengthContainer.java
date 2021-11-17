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

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.ArrayList;
import java.util.Arrays;

import javax.swing.JButton;
import javax.swing.JDialog;

import ccs.graph.Edge;
import ccs.graph.TDLinkage;

public class EdgeLengthContainer extends JDialog implements MouseListener {
	// Main main;
	// Graph graph;
	BoundingBox box;

	final int WIDTH = 160;
	final int HEIGHT_INT = 15;
	int HEIGHT;

	// private TreeDecompModel tdModel = TreeDecompModel.getInstance();
	private TDLinkage td;
	private int n;

	private ArrayList<Edge> projectionEdges;
	private Edge confirmedProjectionEdges[];

	public Edge[] getProjectionEdges() {
		return confirmedProjectionEdges;
	}

	private JButton confirmBtn;

	public void refresh(TDLinkage t) {
		td = t;
		n = td.getNumOfConstructStep();
		HEIGHT = HEIGHT_INT * (n + 1) + 2 * starty;

		setVisible(true);
		setPreferredSize(new Dimension(WIDTH, HEIGHT));
		this.setVisible(false);
	}

	public EdgeLengthContainer(TDLinkage t) {
		// main = m;
		// graph = g;
		if(t.getNumOfConstructStep()<3)
			return;
		projectionEdges = new ArrayList<Edge>();
		confirmedProjectionEdges = new Edge[3];
		for (int i = 0; i < 3; ++i) {
			confirmedProjectionEdges[i] = t.completeCayleyVector().get(i);
			projectionEdges.add(confirmedProjectionEdges[i]);
		}

		addMouseListener(this);

		refresh(t);
	}

	/**
	 * Confirm the set of projection edges.
	 */
	private void confirm() {
		assert (projectionEdges.size() == 3);
		for (int i = 0; i < 3; ++i) {
			confirmedProjectionEdges[i] = projectionEdges.get(i);
		}
		Debug.msg("confirmed");
	}

	final private int starty = 25;
	final private int startx = 10;

	/**
	 * @return index of edge e in the list of projection edges
	 */
	public int indexOf(Edge e) {
		for (int i = 0; i < 3; ++i) {
			if (e.equals(confirmedProjectionEdges[i]))
				return i;
		}
		return -1;
	}

	public void paint(Graphics g) {
		super.paint(g);
		Graphics2D g2d = (Graphics2D) g;

		int i = 0;
		for (Edge e : td.completeCayleyVector()) {
			double distance = td.graph.length(e);
			int endx = (int) (distance / 3) + startx;

			int axisy = HEIGHT_INT * (i + 1) + starty;
			g2d.setStroke(MyStrokes.solid);

			switch (projectionEdges.indexOf(e)) {
			case 0:
				g2d.setColor(Color.red);
				break;
			case 1:
				g2d.setColor(Color.blue);
				break;
			case 2:
				g2d.setColor(Color.green);
				break;
			default:
				g2d.setColor(Color.DARK_GRAY);
				break;
			}

			g2d.drawLine(startx, axisy, endx, axisy);
			g2d.setFont(new Font("SansSerif", Font.PLAIN, 10));
			g2d.drawString(
					e.toString() + ":  " + String.format("%.2f", distance),
					startx, axisy - 3);
			++i;
		}
		g2d.setColor(Color.black);
		g2d.setFont(new Font("SansSerif", Font.BOLD, 15));
		int axisy = HEIGHT_INT * (i + 1) + starty / 2 * 3;
		g2d.drawString("OK", WIDTH / 2 - 10, axisy);

	}

	public void createDialog() {
		setSize(new Dimension(WIDTH, HEIGHT));
		setLocationRelativeTo(GPanel.getInstance());
		setVisible(true);
		repaint();
	}

	public void closeDialog() {
		setVisible(false);
		dispose();
	}

	public void mouseClicked(MouseEvent e) {
		int my = e.getY();
		int edgeIndex = (my - starty) / HEIGHT_INT; // TODO: ???
		Debug.msg("mouseY:" + my + "; edge:" + edgeIndex);
		if (edgeIndex >= n) {
			if (projectionEdges.size() == 3)
				confirm();
			return;
		}
		if (edgeIndex < 0)
			return;
		pick(edgeIndex);
	}

	/**
	 * Add / remove the given edge to projection edge set
	 * 
	 * @param edgeIndex
	 *            index of the edge
	 * @return 0: added; 1: removed; -1: can't add because projection edge set
	 *         is full
	 */
	private int pick(int edgeIndex) {

		Edge e = td.completeCayleyVector().get(edgeIndex);
		if (projectionEdges.contains(e)) {
			projectionEdges.remove(e);
			// confirmBtn.setEnabled(false);
			repaint();
			return 1;
		} else if (projectionEdges.size() < 3) {
			projectionEdges.add(e);
			// if (projectionEdges.size() == 3)
			// confirmBtn.setEnabled(true);
			repaint();
			return 0;
		} else
			return -1;
	}

	public void mouseEntered(MouseEvent arg0) {
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

}
