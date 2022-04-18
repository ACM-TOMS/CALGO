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
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.Point;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.util.ArrayList;

import javax.swing.JDialog;
import javax.swing.WindowConstants;

import ccs.graph.TwoTuple;

/**
 * Connected component curve(s) projected on 3D.
 * 
 */
public class Curve3DContainer extends JDialog implements MouseListener,
		MouseMotionListener {

	private int width, height;
	private int mx, my; // the most recently recorded mouse coordinates

	/**
	 * List of curves to be shown on this panel
	 */
	private ArrayList<Curve3DMotion> curves;

	private Image backbuffer;
	// private Image curBuffer;

	Graphics backg;

	int azimuth = 35, elevation = 30;

	// private Main main;
	// BoundingBox box;

	private final int WIDTH = 400;
	private final int HEIGHT = 400;

	// Point3D[] vertices;
	// Edge[] edges;

	public Curve3DContainer() {
		curves = new ArrayList<Curve3DMotion>();
		init();
	}

	/**
	 * Redraw the curves
	 */
	public void refresh() {
		backg = backbuffer.getGraphics();
		drawCurves(backg);
	}

	public void init() {
		// width = getSize().width;
		// height = getSize().height;
		width = WIDTH;
		height = HEIGHT;
		this.setVisible(true);
		this.setPreferredSize(new Dimension(WIDTH, HEIGHT));
		backbuffer = createImage(WIDTH, HEIGHT);
		setDefaultCloseOperation(WindowConstants.DO_NOTHING_ON_CLOSE);

		// Debug.warnMsg("backbuffer:" + backbuffer);

		addMouseListener(this);
		addMouseMotionListener(this);

		this.setVisible(false);
	}

	double theta, phi;
	float cosT, sinT, cosP, sinP, cosTcosP, cosTsinP, sinTcosP, sinTsinP;

	int scaleFactor;
	float near = 3; // distance from eye to near plane
	float nearToObj = 1.5f; // distance from near plane to center of object

	private Point project(Point3D p) {
		int x0 = p.x;
		int y0 = p.y;
		int z0 = p.z;

		// compute an orthographic projection
		float x1 = cosT * x0 + sinT * z0;
		float y1 = -sinTsinP * x0 + cosP * y0 + cosTsinP * z0;

		// now adjust things to get a perspective projection
		float z1 = cosTcosP * z0 - sinTcosP * x0 - sinP * y0;
		x1 = x1 * near / (z1 + near + nearToObj);
		y1 = y1 * near / (z1 + near + nearToObj);

		// the 0.5 is to round off when converting to int
		return new Point((int) (width / 2 + scaleFactor * x1 + 0.5),
				(int) (height / 2 - scaleFactor * y1 + 0.5));
	}

	/**
	 * Draw the curves on g
	 */
	private void drawCurves(Graphics g) {
		// compute coefficients for the projection
		theta = Math.PI * azimuth / 180.0;
		phi = Math.PI * elevation / 180.0;
		cosT = (float) Math.cos(theta);
		sinT = (float) Math.sin(theta);
		cosP = (float) Math.cos(phi);
		sinP = (float) Math.sin(phi);
		cosTcosP = cosT * cosP;
		cosTsinP = cosT * sinP;
		sinTcosP = sinT * cosP;
		sinTsinP = sinT * sinP;

		scaleFactor = width / 4;

		// Background: white
		g.setColor(Color.white);
		g.fillRect(0, 0, width, height);

		Graphics2D g2d = (Graphics2D) g;
		g2d.setStroke(MyStrokes.solid_thin);

		// Draw 3 axis
		Point o = project(new Point3D(0, 0, 0));
		Point x = project(new Point3D(1000, 0, 0));
		Point y = project(new Point3D(0, 1000, 0));
		Point z = project(new Point3D(0, 0, 1000));
		g.setColor(Color.red);
		g.drawLine(o.x, o.y, x.x, x.y);
		g.setColor(Color.blue);
		g.drawLine(o.x, o.y, y.x, y.y);
		g.setColor(Color.green);
		g.drawLine(o.x, o.y, z.x, z.y);

		for (Curve3DMotion curve : curves) {
			// project vertices onto the 2D viewport
			Point[] points;
			points = new Point[curve.size()];
			int j;
			for (j = 0; j < curve.size(); ++j) {
				points[j] = project(curve.getVertex(j));
			}

			// draw the curve
			// g.setColor(Color.white);
			for (j = 0; j < curve.edgeSize(); ++j) {
				g.setColor(curve.getColor(j));
				g.drawLine(points[curve.getEdge(j).a].x,
						points[curve.getEdge(j).a].y,
						points[curve.getEdge(j).b].x,
						points[curve.getEdge(j).b].y);
			}

			// draw the focal point, if there is one
			for (TwoTuple<Point3D, Color> tuple : curve.getFocalPoints()) {
				Point3D p = tuple.getFirst();
				if (p != null) {
					g.setColor(tuple.getSecond());
					draw3DPoint(p, g);
				}
			}
		}
	}

	public void mouseEntered(MouseEvent e) {
	}

	public void mouseExited(MouseEvent e) {
	}

	public void mouseClicked(MouseEvent e) {
	}

	public void mousePressed(MouseEvent e) {
		mx = e.getX();
		my = e.getY();
		e.consume();
	}

	public void mouseReleased(MouseEvent e) {
	}

	public void mouseMoved(MouseEvent e) {
	}

	public void mouseDragged(MouseEvent e) {
		// get the latest mouse position
		int new_mx = e.getX();
		int new_my = e.getY();

		// adjust angles according to the distance travelled by the mouse
		// since the last event
		azimuth -= new_mx - mx;
		elevation += new_my - my;

		// update the backbuffer
		drawCurves(backg);

		// update our data
		mx = new_mx;
		my = new_my;

		repaint();
		e.consume();
	}

	private void draw3DPoint(Point3D p3d, Graphics g) {
		Point p = project(p3d);
		g.fillOval(p.x - 5, p.y - 5, 10, 10);
	}

	private void drawCurrent(Graphics g) {
		g.setColor(Color.red);

		for (Curve3D curve : curves) {
			if (curve instanceof Curve3DMotion) {
				Curve3DMotion c = (Curve3DMotion) curve;
				draw3DPoint(c.getCurrentPoint(), g);
				// g.fillOval(p.x, p.y, 10, 10);
			}
		}
	}

	public void paint(Graphics g) {
		super.paint(g);

		g.drawImage(backbuffer, 0, 0, this);

		// Draw the current point
		drawCurrent(g);

	}

	/**
	 * Show a dialog
	 */
	public void createDialog() {
		Debug.msg(this.curves.get(0).toString());
		this.refresh();

		/*
		 * GPanel.getInstance().add(this); this.setVisible(true);
		 * GPanel.getInstance().validate(); GPanel.getInstance().repaint();
		 * 
		 * dialog = new JDialog();
		 */

		/*
		 * JButton closeBtn = new JButton("Close");
		 * closeBtn.addActionListener(new ActionListener() { public void
		 * actionPerformed(ActionEvent arg0) { dialog.setVisible(false);
		 * dialog.dispose(); } });
		 */

		// JPanel contentPane = new JPanel();
		// contentPane.add(this);
		// contentPane.add(closeBtn);
		// contentPane.setOpaque(true);
		// dialog.setContentPane(this);

		// dialog.getContentPane().add(this);
		JDialog dialog = this;
		dialog.setSize(new Dimension(WIDTH + 10, HEIGHT + 10));
		dialog.setLocationRelativeTo(GPanel.getInstance());
		dialog.setVisible(true);
		repaint();
	}

	public void closeDialog() {
		JDialog dialog = this;
		dialog.setVisible(false);
		dialog.dispose();
	}

	public boolean addCurve(Curve3DMotion c) {
		return curves.add(c);
	}

	public boolean removeCurve(Curve3D c) {
		return curves.remove(c);
	}

	public void clear() {
		curves.clear();
	}
}
