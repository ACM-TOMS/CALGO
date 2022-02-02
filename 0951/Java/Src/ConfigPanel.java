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
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.geom.AffineTransform;

import javax.swing.JLabel;
import javax.swing.JPanel;

import ccs.graph.Edge;
import ccs.graph.Realization;
import ccs.graph.Point2D;
import ccs.graph.RealizationType;
import ccs.graph.Util;
import ccs.graph.Vertex;

public class ConfigPanel extends JPanel implements MouseMotionListener,
		MouseListener {
	Main main;
	Realization graph;

	BoundingBox box;

	double lf;

	int myid = -1;
	JLabel idLabel;

	// TODO:
	static int id = 0;
	
	public Realization getGraph(){
		return graph;
	}

	public ConfigPanel(Main m, Realization g) {
		this(m, g, 0.4, true);
	}

	private TDLinkageModel tdModel = TDLinkageModel.getInstance();

	public ConfigPanel(Main m, Realization g, double ratio, boolean isThumb) {
		this.main = m;
		lf = g.length(tdModel.getBaseNonedge());

		if (isThumb) {
			id++;
			myid = id;
			idLabel = new JLabel("(" + myid + ")");
			this.add(idLabel);
		}

		box = new BoundingBox(g);
		// System.out.println(box);

		AffineTransform t = new AffineTransform();
		t.scale(ratio, ratio);
		t.translate(-box.x_min + 25, -box.y_min + 25);
		graph = g.clone(); // ???
		graph.transformVertices(t);

		this.setBackground(Color.gray);
		this.setPreferredSize(new Dimension(
				(int) (box.getWidth() * ratio) + 40,
				(int) (box.getHeight() * ratio) + 40)); // 0.4 * GraphPanel
		// this.setPreferredSize(new Dimension(300, 220)); //0.4 * GraphPanel
		this.addMouseMotionListener(this);
		this.addMouseListener(this);
	}

	final int pointsize = 6;

	@Override
	public void paint(Graphics g) {
		// super.paint(g);
		if (myid > 0)
			idLabel.paint(g);

		Graphics2D g2d = (Graphics2D) g;

		for (int i = 0; i < graph.getSize(); ++i) {
			Vertex v = graph.getVertex(i);
			for (Vertex u : graph.getNeighbors(v)) {
				Color c = Color.red;
				if (tdModel.isTd()) {
					RealizationType t = tdModel.getForwardSolutionType(graph);
					c = t.getColor();
				}
				g2d.setColor(c);
				Point2D pv = graph.getPoint(v), pu = graph.getPoint(u);
				g2d.setStroke(MyStrokes.solid_thin);
				g2d.drawLine((int) pv.x(), (int) pv.y(), (int) pu.x(),
						(int) pu.y());
			}
		}

		if (tdModel.isTd()) {
			g2d.setStroke(MyStrokes.dashed);
			g2d.setColor(Color.blue);
			Edge base = tdModel.getBaseNonedge();
			Point2D pv = graph.getPoint(base.v1()), pu = graph.getPoint(base
					.v2());
			g2d.drawLine((int) pv.x(), (int) pv.y(), (int) pu.x(), (int) pu.y());
			if (myid > 0)
				g2d.drawString(((double) ((int) (lf * 100))) / 100 + "",
						(float) ((pv.x() + pu.x()) / 2),
						(float) ((pv.y() + pu.y()) / 2));
		}

		for (int i = 0; i < graph.getSize(); ++i) {

			Vertex v = graph.getVertex(i);
			Color c = Util.randomColor(i);
			g2d.setColor(c);

			Point2D p = graph.getPoint(v);
			g2d.fillOval((int) p.x() - pointsize / 2, (int) p.y() - pointsize
					/ 2, pointsize, pointsize);

			g2d.setColor(Color.black);
			g2d.setFont(new Font("SansSerif", Font.BOLD, 18));
			g2d.drawString(v.toString(), (float) p.x() + pointsize,
					(float) p.y() + pointsize);
		}
	}

	int pressX, pressY;

	public void mouseDragged(MouseEvent e) {
		this.setLocation(e.getXOnScreen() - pressX, e.getYOnScreen() - pressY);
		// this.getParent().validate();
		this.getParent().repaint();
	}

	public void mouseMoved(MouseEvent arg0) {
	}

	public void mouseClicked(MouseEvent e) {
		if (e.getClickCount() == 2) {
			GPanel.getInstance().removeThumbNail(this);
		}

	}

	public void mouseEntered(MouseEvent e) {
	}

	public void mouseExited(MouseEvent e) {
	}

	public void mousePressed(MouseEvent e) {
		pressX = e.getX();
		pressY = e.getY();
		// System.out.println("press @ ("+ pressX +","+pressY+")");
	}

	public void mouseReleased(MouseEvent e) {
	}
}
