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

import java.awt.geom.AffineTransform;
import java.awt.geom.NoninvertibleTransformException;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;

import ccs.graph.CayleyConfigSpace;
import ccs.graph.Cluster;
import ccs.graph.ConnectedComponent;
import ccs.graph.ConstructionStep;
import ccs.graph.ContinuousMotion;
import ccs.graph.Edge;
import ccs.graph.Realization;
import ccs.graph.Point2D;
import ccs.graph.RealizationType;
import ccs.graph.TDLinkage;
import ccs.graph.Vector2D;
import ccs.graph.Vertex;

public class TDLinkageModel {
	private Realization graph = new Realization();
	private TDLinkage td; // when upadte the graph, set td to null;

	private static TDLinkageModel me;

	// private ControlPanel control = ControlPanel.getInstance();

	private TDLinkageModel() {
	}

	public static TDLinkageModel getInstance() {
		if (me == null)
			me = new TDLinkageModel();
		return me;
	}

	public Realization getGraph() {
		// return graph rotated to specified axis.
		if (origin != null && xAxis != null) {
			Realization gg = graph.clone();

			Point2D curXP = gg.getPoint(xAxis);
			Point2D curOP = gg.getPoint(origin);
			Vector2D rotateV = curXP.minus(curOP);

			AffineTransform r = new AffineTransform();
			r.rotate(rotateV.x(), rotateV.y());
			try {
				r.invert();
			} catch (NoninvertibleTransformException e) {
				e.printStackTrace();
			}
			gg.transformVertices(r);

			curOP = gg.getPoint(origin);
			// Debug.warnMsg("original op:" + originP + ", current op:" +
			// curOP);
			AffineTransform t = new AffineTransform();
			Vector2D translateV = originP.minus(curOP);
			t.translate(translateV.x(), translateV.y());

			gg.transformVertices(t);
			return gg;
		}

		return graph;
	}

	public TDLinkage getTd() {
		return me.td;
	}

	public void clearTd() {
		me.td = null;
		me.origin = null;
		me.originP = null;
		me.xAxis = null;
	}

	/**
	 * @return whether current graph is 1-dof td
	 */
	public boolean isTd() {
		return me.td != null;
	}

	public boolean isLow() {
		assert (isTd());
		return me.td.isLow();
	}

	public boolean is1Path() {
		assert (isTd());
		return me.td.is1Path();
	}

	public boolean isTriangleFree() {
		assert (isTd());
		return me.td.isTriangleFree();
	}

	public void refresh() {
		me.clearTd();
		// Step 1. test if it is td
		TDLinkage t = new TDLinkage(me.graph);
		if (t.is1DofTreeDecomposable()) {
			me.td = t;
		} else {
			me.td = null;
			return;
		}

		// Step 2. test 1-path / t-free / low
		// if (!me.td.constructionSequenceGenerated())
		me.td.generateConstructionSequence();
		me.td.normalizeBaseNonEdge();

		//Debug.msg(me.td+"");
		FloatingPanels.getInstance().refresh();
	}

	public void writeToStream(ObjectOutputStream o) throws IOException {
		graph.writeToStream(o);
	}

	public void readFromStream(ObjectInputStream in) throws IOException {
		try {
			graph = Realization.readFromStream(in);
		} catch (ClassNotFoundException c) {
			Debug.warnMsg("ListGraph class not found.");
			c.printStackTrace();
			return;
		}
		refresh();
	}

	public CayleyConfigSpace getCayleyConfigurationSpace() {
		if (!td.ccsGenerated())
			td.genCayleyConfigSpace();
		return td.cayleyConfigSpace;
	}

	public RealizationType getForwardSolutionType() {
		return td.getForwardSolutionType();
	}

	public RealizationType getForwardSolutionType(Realization g) {
		return td.getForwardSolutionType(g);
	}

	public Realization tryRealize(Double val, RealizationType forward) {
		return td.tryRealize(val, forward);
	}

	public ArrayList<Realization> tryRealize(Double val) {
		return td.tryRealize(val);
	}

	public Realization tryRealize(RealizationType t) {
		return td.tryRealize(t);
	}

	public void setPoints(Realization g) {
		graph.setPoints(g);
	}

	public double getBaseNoneedgeLength() {
		return td.getBaseNoneedgeLength();
	}

	public Realization getGraphClone() {
		return graph.clone();
	}

	public Edge getBaseNonedge() {
		return td.getBaseNonedge();
	}

	public Cluster getClusterBetween(Vertex v, Vertex u) {
		return td.getClusterBetween(v, u);
	}

	public int indexOf(Cluster c) {
		return td.indexOf(c);
	}

	public ConstructionStep getConstructionStep(int stepNum) {
		return td.getConstructionStep(stepNum);
	}

	public boolean isAdjacent(Vertex v, Vertex u) {
		return graph.isAdjacent(v, u);
	}

	public Vertex getVertex(Point2D clickP, int clickprecision) {
		return getGraph().getVertex(clickP, clickprecision);
	}

	public int isStepVertex(Vertex v) {
		return td.isStepVertex(v);
	}

	public void addVertex(Point2D clickP) {
		graph.addVertex(clickP);
	}

	public void removeVertex(Vertex v) {
		graph.removeVertex(v);
	}

	public void addEdge(Vertex v, Vertex u) {
		graph.addEdge(v, u);
	}

	public void removeEdge(Vertex v, Vertex u) {
		graph.removeEdge(v, u);
	}

	public void setPoint(Vertex selectedV, int x, int y) {
		graph.setPoint(selectedV, x, y);
	}

	public Point2D getPoint(Vertex v) {
		return getGraph().getPoint(v);
	}

	public ContinuousMotion findPath(Realization startRealization,
			Realization endRealization) {
		return td.findPath(startRealization, endRealization);
	}

	public ContinuousMotion findPath(Double startCayleyConfig,
			Double endCayleyConfig) {
		return td.findPath(startCayleyConfig, endCayleyConfig);
	}

	public ConnectedComponent findComponent() {
		return td.findPathFrom(graph);
	}

	public ConnectedComponent findComponent(Realization g) {
		return td.findPathFrom(g);
	}

	public int getNumOfConstructStep() {
		return td.getNumOfConstructStep();
	}

	public Iterable<Edge> getCompleteCayleyVector() {
		return td.completeCayleyVector();
	}

	private Vertex origin;
	private Point2D originP;
	private Vertex xAxis;

	public void setAxis(Vertex o, Vertex x) {
		origin = o;
		originP = graph.getPoint(o);
		xAxis = x;
	}

	// TODO: remove this.
	Realization transformG(Realization g) {
		if (origin != null && xAxis != null) {
			Realization gg = g.clone();

			Point2D curXP = gg.getPoint(xAxis);
			Point2D curOP = gg.getPoint(origin);
			Vector2D rotateV = curXP.minus(curOP);

			AffineTransform r = new AffineTransform();
			r.rotate(rotateV.x(), rotateV.y());
			try {
				r.invert();
			} catch (NoninvertibleTransformException e) {
				e.printStackTrace();
			}
			gg.transformVertices(r);

			curOP = gg.getPoint(origin);
			// Debug.warnMsg("original op:" + originP + ", current op:" +
			// curOP);
			AffineTransform t = new AffineTransform();
			Vector2D translateV = originP.minus(curOP);
			t.translate(translateV.x(), translateV.y());

			gg.transformVertices(t);
			return gg;
		}
		return g;
	}
}
