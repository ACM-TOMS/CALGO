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
import java.awt.Font;
import java.awt.GradientPaint;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.geom.AffineTransform;
import java.util.ArrayList;
import java.util.HashSet;

import javax.swing.JPanel;

import ccs.graph.Cluster;
import ccs.graph.ConstructionStep;
import ccs.graph.Edge;
import ccs.graph.Realization;
import ccs.graph.Point2D;
import ccs.graph.RealizationType;
import ccs.graph.Util;
import ccs.graph.Vertex;

public class GPanel extends JPanel implements MouseListener,
		MouseMotionListener {

	private static GPanel me = new GPanel();

	public static GPanel getInstance() {
		return me;
	}

	private GPanel() {
		addMouseListener(this);
		addMouseMotionListener(this);

		addMouseListener(FixAxisListener.getInstance());

		thumbnails = new ArrayList<ConfigPanel>();
	}

	private Main main;

	public void setMain(Main main) {
		this.main = main;
	}

	private ControlPanel control = ControlPanel.getInstance();

	private TDLinkageModel tdModel = TDLinkageModel.getInstance();
	private CCSModel ccsModel = CCSModel.getInstance();
	private MotionModel motionModel = MotionModel.getInstance();
	private TraceModel traceModel = TraceModel.getInstance();

	/**
	 * The vertex being dragged
	 */
	private Vertex selectedV = null;

	/**
	 * The start vertex for add edge
	 */
	private Vertex edgeStartV = null;

	/**
	 * Currently shown construction step
	 */
	public int showStep = -1;

	/**
	 * List of thumbnail panels
	 */
	private ArrayList<ConfigPanel> thumbnails;

	public void clearThumbnails() {
		for (ConfigPanel c : thumbnails) {
			this.remove(c);
		}
		// ConfigPanel.id = 0;
	}

	private boolean isHidingMainRealization = false;

	/**
	 * Panels for showing path start / end
	 */
	private ConfigPanel pathStart, pathEnd;

	public void clearPathEnds() {
		pathStart = pathEnd = null;
	}

	public void clearDrawingBuffer() {
		edgeStartV = null;
		repaint();
	}

	public void clear() {
		clearThumbnails();
		clearPathEnds();
	}

	// for canonical base non-edges
	// private boolean showCanonical = false;
	private EdgeLengthContainer canonicalEdgePanel;

	/**
	 * Drawing size of vertices
	 */
	private final int pointsize = 16;

	/**
	 * Draw the base non-edge using dash line
	 * 
	 * @param graph
	 *            the realization to draw
	 */
	private void drawBase(Realization graph, Graphics2D g2d) {
		g2d.setStroke(MyStrokes.dashed);
		g2d.setColor(Color.blue);
		assert (tdModel.isTd());
		Edge base = tdModel.getBaseNonedge();
		assert (base != null);
		drawEdge(graph, base, g2d);
		if (!control.isDrawing()) {
			drawEdgeLength(graph, base, g2d);
		}
	}

	private void drawEdgeLength(Realization graph, Edge e, Graphics2D g2d) {
		drawEdgeLength(graph, e.v1(), e.v2(), g2d);
	}

	private void drawEdgeLength(Realization graph, Vertex v, Vertex u, Graphics2D g2d) {
		Point2D pv = graph.getPoint(v), pu = graph.getPoint(u);
		g2d.setFont(new Font("SansSerif", Font.PLAIN, 10));
		g2d.drawString(String.format("%.3f", graph.distance(v, u)),
				(float) ((pv.x() + pu.x()) / 2),
				(float) ((pv.y() + pu.y()) / 2));
	}

	private final boolean gradientStyle = false;

	/**
	 * Draw a dot. Called by drawVertex.
	 * 
	 * @param p
	 *            the point to draw
	 */
	private void drawPoint(Point2D p, Graphics2D g2d) {

		// if (Graph.DEBUG) {
		// g2d.setFont(new Font("SansSerif", Font.PLAIN, 12));
		// g2d.drawString(String.format("(%.2f,%.2f)", p.x(), p.y()),
		// (float) p.x(), (float) p.y());
		// }

		// For gradient effect
		if (gradientStyle) {
			GradientPaint gp = new GradientPaint((int) p.x() - pointsize / 2,
					(int) p.y() - pointsize / 2, Color.white, (int) p.x()
							+ pointsize / 2, (int) p.y() + pointsize / 2,
					g2d.getColor());
			g2d.setPaint(gp);
		}

		g2d.fillOval((int) p.x() - pointsize / 2, (int) p.y() - pointsize / 2,
				pointsize, pointsize);

	}

	private void drawVertex(Realization graph, Vertex v, Graphics2D g2d) {
		Point2D p = graph.getPoint(v);
		drawPoint(p, g2d);
		g2d.setColor(Color.black);
		g2d.setFont(new Font("SansSerif", Font.BOLD, 24));
		g2d.drawString(v.toString(), (float) p.x() + pointsize, (float) p.y()
				+ pointsize);
	}

	private void drawEdge(Realization graph, Vertex v, Vertex u, Graphics2D g2d) {
		Point2D pv = graph.getPoint(v), pu = graph.getPoint(u);
		if (control.isDrawing()) {
			drawEdgeLength(graph, v, u, g2d);
		}

		Color c = g2d.getColor();
		// Point2D mid = pv.add(pu.minus(pv).divide(2));
		if (gradientStyle) {
			GradientPaint gp = new GradientPaint((int) pv.x(), (int) pv.y(),
					Color.white, (int) pu.x(), (int) pu.y(), c, true);
			g2d.setPaint(gp);
		}
		g2d.drawLine((int) pv.x(), (int) pv.y(), (int) pu.x(), (int) pu.y());
	}

	private void drawEdge(Realization graph, Edge e, Graphics2D g2d) {
		drawEdge(graph, e.v1(), e.v2(), g2d);
	}

	private void drawGraph(Realization graph, Graphics2D g2d) {
		g2d.setStroke(MyStrokes.solid);
		for (int i = 0; i < graph.getSize(); ++i) {
			Vertex v = graph.getVertex(i);
			for (Vertex u : graph.getNeighbors(v)) {
				if (control.isShowingCluster()) {
					// set color according to clusters
					Cluster c = tdModel.getClusterBetween(v, u);
					if (c != null) {
						g2d.setColor(Util.randomColor(tdModel.indexOf(c)));
					}
				} else {
					// set color according to solution type
					Color c = Color.red;
					if (tdModel.isTd()) {
						RealizationType t = tdModel.getForwardSolutionType(graph);
						c = Util.randomColor(t.getEncoding());
					}
					g2d.setColor(c);
				}

				drawEdge(graph, v, u, g2d);
			}
		}

		if (tdModel.isTd())
			drawBase(graph, g2d);

		// For drawing: pink the chosen vertex
		for (int i = 0; i < graph.getSize(); ++i) {
			Vertex v = graph.getVertex(i);
			if (v == edgeStartV)
				g2d.setColor(Color.PINK);
			else {
				Color c = Util.randomColor(v.index);
				g2d.setColor(c);
			}
			drawVertex(graph, v, g2d);
		}
	}

	private void drawStep(Graphics2D g2d) {
		Debug.msg("drawing steps to " + showStep);
		g2d.setColor(Color.GREEN);
		for (int stepNum = 0; stepNum <= showStep; ++stepNum) {
			ConstructionStep s = tdModel.getConstructionStep(stepNum);
			HashSet<Vertex> vset = s.c1().unionVertices(s.c2());
			ArrayList<Vertex> vs = new ArrayList<Vertex>();
			vs.addAll(vset);
			for (int i = 0; i < vs.size(); ++i) {
				Vertex v = vs.get(i);
				for (int j = 0; j < i; ++j) {
					Vertex u = vs.get(j);
					if (tdModel.isAdjacent(v, u)) {
						drawEdge(tdModel.getGraph(), v, u, g2d);
					}
				}
				drawVertex(tdModel.getGraph(), v, g2d);
			}
		}
		drawBase(tdModel.getGraph(), g2d);
	}

	// TODO: change this.
	private void drawTrace(Graphics2D g2d) {
		// /g2d.setColor(Color.lightGray);
		for (int index = 0; index < traceModel.getTracingVertices().size(); ++index) {
			Vertex tracingVertex = traceModel.getTracingVertices().get(index);

			g2d.setColor(Util.randomColor(tracingVertex.index));
			Point2D pv = tdModel.getPoint(tracingVertex);

			boolean traceCreated = false;
			double lf = tdModel.getBaseNoneedgeLength();
			RealizationType t = tdModel.getForwardSolutionType();
			for (Trace trace : traceModel.getTraces().get(index)) {
				for (int i = 1; i < trace.tracePoints.size(); ++i) {
					Point2D p1 = trace.tracePoints.get(i - 1).point, p2 = trace.tracePoints
							.get(i).point;
					g2d.setStroke(MyStrokes.solid_thin);
					g2d.drawLine((int) p1.x(), (int) p1.y(), (int) p2.x(),
							(int) p2.y());
				}
				if (trace.canAddPoint(lf, t)) {
					trace.addPoint(pv, lf);
					traceCreated = true;
				}
			}
			if (!traceCreated) {
				Trace trace = new Trace(tracingVertex, tdModel.getTd().cayleyConfigSpace
						.getOrientedCCS(t).getContainingInterval(lf), t,
						tdModel.getTd());
				trace.addPoint(pv, lf);
				traceModel.getTraces().get(index).add(trace);
			}
		}

		// for (Point2D p : trace) {
		//
		// g2d.fillOval((int) p.x() - min_distance, (int) p.y()
		// - min_distance, min_distance * 2, min_distance * 2);
		//
		// if (p.distance(pv) < min_distance)
		// addit = false;
		// }
		// if (addit)
		// trace.add(pv);
	}

	/**
	 * End of tracing paths: clear the temporary motion panels
	 */
	public void clearMotionEndPanels() {
		pathStart = null;
		pathEnd = null;
	}

	@Override
	public void paint(Graphics g) {
		super.paint(g);
		Graphics2D g2d = (Graphics2D) g;

		if (showStep >= 0) {
			drawStep(g2d);
			return;
		}

		// draw trace
		if (control.isTracingVertices()) {// && tracingVertex != null) {
			drawTrace(g2d);
		}

		// Draw start realization
		Realization start = motionModel.getStartRealization();
		if (start != null) {
			if (pathStart == null) {
				pathStart = new ConfigPanel(main, start, 0.5, false);
				pathStart.paint(g.create(20, 20, getWidth(), getHeight()));
			} else {
				pathStart.paint(g.create(20, 20, getWidth(), getHeight()));
			}
		}

		// Draw end realization
		Realization end = motionModel.getEndRealization();
		if (end != null) {
			if (pathEnd == null) {
				pathEnd = new ConfigPanel(main, end, 0.5, false);
				pathEnd.paint(g.create(
						(int) (getWidth() - pathEnd.box.getWidth()), 20,
						(int) pathEnd.box.getWidth(),
						(int) pathEnd.box.getHeight()));

			} else {
				pathEnd.paint(g.create(
						(int) (getWidth() - pathEnd.box.getWidth()), 20,
						(int) pathEnd.box.getWidth(),
						(int) pathEnd.box.getHeight()));
			}
		}

		// Draw main realization
		g2d.setColor(Color.black);

		// No path & 2 components: draw realization for each component
		if (control.isTracingMotion() == 3) {
			// TODO: mechanism for scale / fix graphs
			Realization g1 = motionModel.getSpinnerModel1().getValue().getValue();
			final double ratio = 0.6;
			BoundingBox box = new BoundingBox(g1);
			AffineTransform t = new AffineTransform();
			t.translate(-box.getWidth() / 6, box.getWidth() / 2);
			t.scale(ratio, ratio);
			Realization g11 = g1.clone(); // ???
			g11.transformVertices(t);
			drawGraph(g11, g2d);

			Realization g2 = motionModel.getSpinnerModel2().getValue().getValue();
			box = new BoundingBox(g2);
			t = new AffineTransform();
			t.translate(box.getWidth() * 2 / 3, box.getWidth() / 2);
			t.scale(ratio, ratio);
			Realization g22 = g2.clone(); // ???
			g22.transformVertices(t);
			drawGraph(g22, g2d);
		} else if (!isHidingMainRealization)
			// Draw the main realization
			drawGraph(tdModel.getGraph(), g2d);

		// Emphasize vertices to be flipped
		if (control.isFlipping()) {
			for (Vertex v : FlippingModel.getInstance().getVerteicsToFlip()) {
				g2d.setColor(Color.RED);
				drawVertex(tdModel.getGraph(), v, g2d);
			}
		}

		// // emphasize the Vertex being traced
		// if (tracing && tracingVertex != null) {
		// g2d.setColor(Color.red);
		// drawVertex(graph, tracingVertex, g2d);
		// }

		// draw length of canonimcal edges
		if (control.isShowingCompleteCayley() && control.isTracingMotion() != 3) {
			drawCompleteCayley(g2d);
		}

		// when fixing axis, emphasize the origin picked
		if (FixAxisListener.getInstance().getO() != null) {
			g2d.setColor(Color.YELLOW);
			this.drawVertex(tdModel.getGraph(), FixAxisListener.getInstance()
					.getO(), g2d);
		}
	}

	private void drawCompleteCayley(Graphics2D g2d) {
		for (Edge e : tdModel.getTd().completeCayleyVector()) {
			g2d.setStroke(MyStrokes.dashedMedium);
			switch (FloatingPanels.getInstance().indexOf(e)) {
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
				g2d.setStroke(MyStrokes.dashed);
				g2d.setColor(Color.GRAY);
				break;
			}
			drawEdge(tdModel.getGraph(), e, g2d);
		}
	}

	// =================================

	public void removeThumbNail(ConfigPanel c) {
		remove(c);
		ccsModel.removeMarkedPoint(c);
		validate();
		repaint();
	}

	// ========================= LISTENER ================================

	private static final int clickPrecision = 15;

	public void mouseClicked(MouseEvent e) {
		GPanel graphPanel = (GPanel) e.getSource();
		Point2D clickP = new Point2D(e.getX(), e.getY());

		// Hide the main graph if Alt+Click
		if (e.isAltDown()) {
			isHidingMainRealization = !isHidingMainRealization;
			repaint();
			// System.out.println("hide graph");
			return;
		}

		// Create a thumbnail when clicked three times
		if (e.getClickCount() == 3 && ccsModel.isGenerated()) {
			ConfigPanel c = new ConfigPanel(main, tdModel.getGraph());
			graphPanel.add(c);
			ccsModel.addMarkedPoint(c);
			c.setVisible(true);
			graphPanel.validate();
			graphPanel.repaint();
			return;
		}

		// flip orientations
		if (control.isFlipping()) {
			Vertex v = tdModel.getVertex(clickP, clickPrecision);
			if (v != null)
				FlippingModel.getInstance().flipVertex(v);
			repaint();
			return;
		}

		// pick tracing vertex
		if (control.isTracingVertices()) { // && graphPanel.tracingVertex ==
											// null) {
			// assert (graphPanel.getStatus() == Status.generated);
			Vertex v = tdModel.getVertex(clickP, clickPrecision);
			if (v != null && !traceModel.getTracingVertices().contains(v)) {
				traceModel.getTracingVertices().add(v);
				traceModel.getTraces().add(new ArrayList<Trace>());
				control.setTracingVertex(v);
			}
			return;
		}

		// ================= draw =================

		Vertex v;
		switch (control.getDrawOption()) {
		case 0:
			break;
		case 1:
			tdModel.addVertex(clickP);
			break;
		case 3:
			v = tdModel.getVertex(clickP, clickPrecision);
			if (v != null)
				tdModel.removeVertex(v);
			break;
		case 2:
			v = tdModel.getVertex(clickP, clickPrecision);
			if (v == null)
				return;
			if (graphPanel.edgeStartV == null) {
				// pick the start vertex
				graphPanel.edgeStartV = v;
			} else {
				// determine the end vertex
				tdModel.addEdge(graphPanel.edgeStartV, v);
				graphPanel.edgeStartV = null;
			}
			break;
		case 4:
			v = tdModel.getVertex(clickP, clickPrecision);
			if (v == null)
				return;
			if (graphPanel.edgeStartV == null) {
				// pick the start vertex
				graphPanel.edgeStartV = v;
			} else {
				// determine the end vertex
				if (tdModel.isAdjacent(graphPanel.edgeStartV, v)) {
					tdModel.removeEdge(graphPanel.edgeStartV, v);
					graphPanel.edgeStartV = null;
				}
			}
			break;
		}
		graphPanel.repaint();
	}

	public void mouseEntered(MouseEvent e) {
	}

	public void mouseExited(MouseEvent e) {
	}

	public void mousePressed(MouseEvent e) {
		if (!control.isDrawing())
			// TODO: move the entire graph
			return;

		if (control.getDrawOption() != 0)
			return;

		// move a vertex
		Point2D pt = new Point2D(e.getX(), e.getY());
		assert (selectedV == null);
		selectedV = tdModel.getVertex(pt, clickPrecision);
	}

	public void mouseReleased(MouseEvent e) {
		// Drop a vertex
		if (selectedV != null) {
			selectedV = null;
		}
	}

	public void mouseDragged(MouseEvent e) {
		// move a vertex
		if (selectedV != null) {
			int x = e.getX(), y = e.getY();
			tdModel.setPoint(selectedV, x, y);
			// p.main.graph.setInitPoint(p.selectedV, x, y);
			repaint();
		}
	}

	public void mouseMoved(MouseEvent arg0) {
	}

}
