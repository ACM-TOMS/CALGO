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
import java.awt.Component;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.Iterator;

import javax.swing.AbstractButton;
import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JToggleButton;
import javax.swing.filechooser.FileNameExtensionFilter;

import ccs.graph.ConnectedComponent;
import ccs.graph.ContinuousMotion;
import ccs.graph.Edge;
import ccs.graph.Realization;
import ccs.graph.SamplePoint;
import ccs.graph.Vertex;

public class ControlPanel extends JPanel {
	private ControlPanel() {
	}

	private static ControlPanel me;;

	public static ControlPanel getInstance() {
		if (me == null)
			me = new ControlPanel();
		return me;
	}

	private Main main;

	public void setMain(Main m) {
		main = m;
	}

	private TDLinkageModel tdModel = TDLinkageModel.getInstance();
	private CCSModel ccsModel = CCSModel.getInstance();
	private MotionModel motionModel = MotionModel.getInstance();
	private CCSPanel ccsPanel = CCSPanel.getInstance();

	/**
	 * Set all components to be invisible
	 */
	private void hideAll() {
		for (Component c : getComponents()) {
			c.setVisible(false);
		}
	}

	/**
	 * Status: drawing graph
	 */
	private void setDrawing() {
		hideAll();
		tdModel.clearTd();
		showDrawControl();
		// controlPanel.repaint();
	}

	/**
	 * Status: have drawn or loaded a graph
	 */
	private void setDrawn() {
		hideAll();
		tdModel.refresh();
		if (!tdModel.isTd()) {
			tdLabel.setText("1-dof tree-decomposable: NO");

			drawBtn.setVisible(true);
			showFileControl();
			tdLabel.setVisible(true);
			return;
		}
		tdLabel.setText("1-dof tree-decomposable: YES");
		onePathLabel.setText(tdModel.is1Path() ? "1-path: YES" : "1-path: NO");
		tfreeLabel.setText(tdModel.isTriangleFree() ? "triangle-free: YES"
				: "triangle-free: NO");
		lowLabel.setText(tdModel.isLow() ? "low Cayley complexity: YES"
				: "low Cayley complexity: NO");

		drawBtn.setVisible(true);
		showFileControl();
		showTdControl();
		if (tdModel.isLow()) {
			showCCSControl();
			ccsSpinner.setVisible(false);
		}

		GPanel.getInstance().clear();
		GPanel.getInstance().repaint();
	}

	public void init() {
		drawBtn = new JToggleButton("Draw");
		drawBtn.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (drawBtn.isSelected())
					setDrawing();
				else
					setDrawn();

			}
		});
		add(drawBtn);

		initFileControl();

		initDrawControl();

		initTdControl();

		initCCSControl();

		initPathControl();

		hideAll();
		showFileControl();
		showDrawControl();
	}

	private JToggleButton drawBtn;

	// JButton generateBtn;

	public boolean isDrawing() {
		return drawBtn.isSelected();
	}

	private JButton saveBtn, loadBtn;
	private JFileChooser fc;

	private void initFileControl() {
		fc = new JFileChooser();
		FileNameExtensionFilter filter = new FileNameExtensionFilter(
				"CCS Graph files", "ccs");
		fc.setFileFilter(filter);

		saveBtn = new JButton("Save");
		saveBtn.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent event) {
				// does not add .ccs extension automatically
				int returnVal = fc.showSaveDialog(main);
				if (returnVal == JFileChooser.APPROVE_OPTION) {
					File file = fc.getSelectedFile();
					try {
						FileOutputStream fileOut = new FileOutputStream(file);
						ObjectOutputStream out = new ObjectOutputStream(fileOut);
						tdModel.writeToStream(out);
						out.close();
						fileOut.close();
					} catch (IOException i) {
						i.printStackTrace();
					}
				}
			}
		});
		add(saveBtn);

		loadBtn = new JButton("Load");
		loadBtn.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent event) {
				int returnVal = fc.showDialog(main, "Load");
				if (returnVal == JFileChooser.APPROVE_OPTION) {
					File file = fc.getSelectedFile();
					try {
						FileInputStream fileIn = new FileInputStream(file);
						ObjectInputStream in = new ObjectInputStream(fileIn);
						tdModel.readFromStream(in);
						in.close();
						fileIn.close();
					} catch (IOException i) {
						i.printStackTrace();
						return;
					}
					setDrawn();
				}
			}
		});
		add(loadBtn);
	}

	private void showFileControl() {
		saveBtn.setVisible(true);
		loadBtn.setVisible(true);
	}

	private JPanel drawOptionPanel;
	private JToggleButton addVertexBtn, addEdgeBtn, removeVertexBtn,
			removeEdgeBtn;

	private void initDrawControl() {
		drawOptionPanel = new JPanel();
		drawOptionPanel.setPreferredSize(new Dimension(220, 80));
		drawOptionPanel
				.setBorder(BorderFactory.createLineBorder(Color.gray, 2));

		addVertexBtn = new JToggleButton("add v");
		addEdgeBtn = new JToggleButton("add e");
		removeVertexBtn = new JToggleButton("remove v");
		removeEdgeBtn = new JToggleButton("remove e");

		final JToggleButton drawOptionBtns[] = { addVertexBtn, addEdgeBtn,
				removeVertexBtn, removeEdgeBtn };

		DrawingListener listener = new DrawingListener(drawOptionBtns);
		for (JToggleButton b : drawOptionBtns) {
			b.addActionListener(listener);
			drawOptionPanel.add(b);
		}

		add(drawOptionPanel);
	}

	private void showDrawControl() {
		drawBtn.setVisible(true);
		drawOptionPanel.setVisible(true);
	}

	/**
	 * @return which draw option button is seletec. 1:add v; 2:add e; 3:remove
	 *         v; 4:remove e; 0: nothing
	 */
	public int getDrawOption() {
		if (addVertexBtn.isSelected())
			return 1;
		else if (addEdgeBtn.isSelected())
			return 2;
		else if (removeVertexBtn.isSelected())
			return 3;
		else if (removeEdgeBtn.isSelected())
			return 4;
		return 0;
	}

	private JLabel tdLabel;
	private JLabel lowLabel, onePathLabel, tfreeLabel;

	/**
	 * Toggle showing clusters
	 */
	private JToggleButton clusterBtn;

	public boolean isShowingCluster() {
		return clusterBtn.isSelected();
	}

	/**
	 * Demonstrate construction steps.
	 */
	private JToggleButton constructBtn;

	/**
	 * Choose a specific bar / non-edge to fix on axis
	 */
	private JButton fixAxisBtn;

	private void initTdControl() {
		tdLabel = new JLabel("1-dof tree-decomposable");
		add(tdLabel);

		onePathLabel = new JLabel("1-Path");
		add(onePathLabel);

		tfreeLabel = new JLabel("Triangle-free");
		add(tfreeLabel);

		lowLabel = new JLabel("low Cayley Complexity");
		add(lowLabel);

		clusterBtn = new JToggleButton("view clusters");
		clusterBtn.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				GPanel.getInstance().repaint();
			}
		});
		add(clusterBtn);

		constructBtn = new JToggleButton("Demonstrate construction");
		constructBtn.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (constructBtn.isSelected()) {
					constructBtn.setEnabled(false);
					Thread t = new Thread() {
						public void run() {
							for (int i = 0; i < tdModel.getTd()
									.getNumOfConstructStep(); ++i) {
								GPanel.getInstance().showStep = i;
								Debug.msg("step " + i);
								GPanel.getInstance().paintImmediately(
										GPanel.getInstance().getBounds());
								try {
									synchronized (this) {
										this.wait(1000);
									}
								} catch (InterruptedException e) {
									e.printStackTrace();
								}
							}
							constructBtn.setSelected(false);
							constructBtn.setEnabled(true);
							GPanel.getInstance().repaint();
							GPanel.getInstance().showStep = -1;
						}
					};
					t.run();
				}
			}
		});
		add(constructBtn);

		fixAxisBtn = new JButton("fix x-axis");
		fixAxisBtn.addActionListener(FixAxisListener.getInstance());
		add(fixAxisBtn);
	}

	private void showTdControl() {
		tdLabel.setVisible(true);
		lowLabel.setVisible(true);
		onePathLabel.setVisible(true);
		tfreeLabel.setVisible(true);
		clusterBtn.setVisible(true);
		constructBtn.setVisible(true);
		fixAxisBtn.setVisible(true);
	}

	private JPanel ccsControlPanel;

	/**
	 * Select from oriented/nonoriented CCS
	 */
	private JRadioButton orientedRadioBtn, nonorientedRadioBtn;

	/**
	 * Generate ccs.
	 */
	private JButton ccsBtn;

	/**
	 * Flip orientations if pressed.
	 */
	private JToggleButton orientationBtn;

	public boolean isFlipping() {
		return orientationBtn.isSelected();
	}

	/**
	 * Navigates ccs using Cayley configurations.
	 */
	private DoubleSpinner ccsSpinner;

	/**
	 * Show the complete cayley configuration space, represented by lengths of
	 * canonical base non-edges
	 */
	private JToggleButton showCompleteCayleyBtn;

	/**
	 * Status: generated the CCS / come back from path tracing
	 */
	private void setGenerated() {
		hideAll();

		drawBtn.setVisible(true);
		showFileControl();
		showTdControl();

		showCCSControl();
		refreshCCSSpinner();

		showPathControl();

		motionModel.clear();

		ccsPanel.repaint();

		GPanel.getInstance().clearPathEnds();
		GPanel.getInstance().repaint();

		tracingMotion = 0;

	}

	private void initCCSControl() {
		ccsBtn = new JButton("Cayley configuration space");
		ccsBtn.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				ccsModel.refresh();
				setGenerated();
			}
		});

		orientedRadioBtn = new JRadioButton("oriented");
		nonorientedRadioBtn = new JRadioButton("non-oriented", true);
		ButtonGroup ccsRadioGroup = new ButtonGroup();
		ccsRadioGroup.add(nonorientedRadioBtn);
		ccsRadioGroup.add(orientedRadioBtn);

		ccsSpinner = new DoubleSpinner();
		ccsSpinner.setStepSize(1);
		ccsSpinner.addChangeListener(new CCSSpinnerListener());

		ccsControlPanel = new JPanel();
		ccsControlPanel.setPreferredSize(new Dimension(220, 160));
		ccsControlPanel
				.setBorder(BorderFactory.createLineBorder(Color.gray, 2));
		ccsControlPanel.add(ccsBtn);
		ccsControlPanel.add(nonorientedRadioBtn);
		ccsControlPanel.add(orientedRadioBtn);
		ccsControlPanel.add(ccsSpinner);

		add(ccsControlPanel);

		orientationBtn = new JToggleButton("flip orientations");
		orientationBtn.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				boolean selected = orientationBtn.isSelected();
				if (!selected) {
					FlippingModel.getInstance().doFlip();
					if (ccsModel.isOriented()) {
						ccsModel.refresh();
					}
					refreshCCSSpinnerCur();
					GPanel.getInstance().repaint();
					CCSPanel.getInstance().repaint();
				}
			}
		});
		ccsControlPanel.add(orientationBtn);

		showCompleteCayleyBtn = new JToggleButton("show complete ccs");
		showCompleteCayleyBtn.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				boolean selected = showCompleteCayleyBtn.isSelected();
				FloatingPanels fp = FloatingPanels.getInstance();
				if (selected)
					fp.showEdgeLengthPanel();
				else
					fp.closeEdgeLengthPanel();

				// GPanel.getInstance().setShowCompleteCayley(selected);
			}
		});
		ccsControlPanel.add(showCompleteCayleyBtn);
	}

	private void refreshCCSSpinner() {
		double min = ccsModel.getMin(), max = ccsModel.getMax();
		ccsSpinner.setMinimum(min);
		ccsSpinner.setMaximum(max);
		refreshCCSSpinnerCur();
	}

	private void refreshCCSSpinnerCur() {
		ccsSpinner.setValue(tdModel.getBaseNoneedgeLength());
	}

	private void showCCSControl() {
		ccsControlPanel.setVisible(true);
		for (Component c : ccsControlPanel.getComponents()) {
			c.setVisible(true);
			c.setEnabled(true);
		}
	}

	public boolean isOrientedCCS() {
		return (orientedRadioBtn.isSelected());
	}

	public boolean isShowingCompleteCayley() {
		return showCompleteCayleyBtn.isSelected();
	}

	private JPanel pathPanel;

	/**
	 * Generate connected component from current realization
	 */
	private JToggleButton genComponentBtn;

	/**
	 * Show the projection of current connected components in 3D
	 */
	private JToggleButton showComponentCurveBtn;

	/**
	 * Show the projection of all connected components in 3D
	 */
	private JToggleButton showComponentsCurveBtn;

	public boolean isShowingComponentCurve() {
		return showComponentCurveBtn.isSelected()
				|| showComponentsCurveBtn.isSelected();
	}

	/**
	 * Start picking path from two cayleys / two realizations
	 */
	private JToggleButton cayleyPathBtn, realizationPathBtn;
	/**
	 * Pick start / end
	 */
	private JButton startBtn, endBtn;

	private JButton genPathBtn;

	private DialSpinner<SamplePoint<Realization>> pathSpinner;
	private DialSpinner<SamplePoint<Realization>> pathSpinner1, pathSpinner2;

	/**
	 * TODO: navigate all components. shown when generated component
	 */
	private DialSpinner<ConnectedComponent> componentSpinner;

	/**
	 * Choose vertices to be traced when pressed.
	 */
	private JToggleButton traceVertexBtn;

	public boolean isTracingVertices() {
		return traceVertexBtn.isSelected();
	}

	/**
	 * Add v to caption of traceVertexBtn
	 */
	public void setTracingVertex(Vertex v) {
		if (traceVertexBtn.getText() == "trace a vertex")
			traceVertexBtn.setText("tracing");
		traceVertexBtn.setText(traceVertexBtn.getText() + " " + v);
	}

	private boolean isTracingComponent() {
		return genComponentBtn.isSelected();
	}

	/**
	 * Called when genComponentBtn/genPathBtn pressed: set the status to
	 * tracing.
	 */
	private void setTracing() {

		if (isTracingComponent()) {
			// Case 1: tracing a component
			pathSpinner.setVisible(true);
			traceVertexBtn.setVisible(true);
			ccsSpinner.setEnabled(false);
			showComponentCurveBtn.setVisible(true);
			pathSpinner.setCyclic(true);
			cayleyPathBtn.setVisible(false);
			realizationPathBtn.setVisible(false);
			componentSpinner.setVisible(true);
			tracingMotion = 1;
		} else {
			// Tried to find a path
			startBtn.setVisible(false);
			endBtn.setVisible(false);
			genPathBtn.setVisible(false);

			if (motionModel.getMotion() != null) {
				// Case 2: found a path
				pathSpinner.setVisible(true);
				traceVertexBtn.setVisible(true);
				ccsSpinner.setEnabled(false);
				showComponentCurveBtn.setVisible(true);
				tracingMotion = 2;
			} else {
				// Case 3: no path
				// TODO: no path: display the 2 nearest points on the ccsPanel
				motionModel.gen2Components();
				pathSpinner1.setVisible(true);
				pathSpinner2.setVisible(true);
				showComponentCurveBtn.setVisible(true);
				tracingMotion = 3;
			}
		}
		repaint();
		ccsPanel.repaint();
	}

	/**
	 * 1: tracing component; 2: tracing path; 3: no path, tracing two
	 * components; 0: not tracing
	 */
	private int tracingMotion = 0;

	/**
	 * @return 1: tracing component; 2: tracing path; 3: no path, tracing two
	 *         components; 0: not tracing
	 */
	public int isTracingMotion() {
		return tracingMotion;
	}

	private void initPathControl() {
		pathPanel = new JPanel();
		pathPanel.setPreferredSize(new Dimension(220, 180));
		pathPanel.setBorder(BorderFactory.createLineBorder(Color.gray, 2));

		genComponentBtn = new JToggleButton("generate connected component");
		genComponentBtn.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (!genComponentBtn.isSelected()) {
					// toggled off
					pathSpinner.setCyclic(false);
					ccsModel.clearTypes();
					setGenerated();
				} else {
					motionModel.genComponent();
					setTracing();
				}

			}
		});

		cayleyPathBtn = new JToggleButton("find path between two cayleys");
		cayleyPathBtn.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (!cayleyPathBtn.isSelected()) {
					ccsModel.clearTypes();
					setGenerated();
				} else {
					startBtn.setVisible(true);
					endBtn.setVisible(true);
					realizationPathBtn.setVisible(false);
					genComponentBtn.setVisible(false);
				}
			}
		});

		realizationPathBtn = new JToggleButton(
				"find path between two realizations");
		realizationPathBtn.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (!realizationPathBtn.isSelected()) {
					ccsModel.clearTypes();
					setGenerated();
				} else {
					startBtn.setVisible(true);
					endBtn.setVisible(true);
					cayleyPathBtn.setVisible(false);
					genComponentBtn.setVisible(false);
				}
			}
		});

		startBtn = new JButton("pick start");
		startBtn.addActionListener(new PathPickListener());

		endBtn = new JButton("pick end");
		endBtn.addActionListener(new PathPickListener());

		genPathBtn = new JButton("generate path");
		genPathBtn.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				MotionModel.getInstance().genPath();
				setTracing();
			}
		});

		showComponentCurveBtn = new JToggleButton("show curve of component");
		showComponentCurveBtn.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				FloatingPanels fp = FloatingPanels.getInstance();

				if (showComponentCurveBtn.isSelected()) {
					if (isTracingMotion() == 1) {
						Curve3DMotion c = new Curve3DMotion(0);
						c.refresh(fp.getProjectionEdges());

						fp.addCurve3D(c);
					} else if (isTracingMotion() == 2) {
						Curve3DMotion c = new Curve3DMotion(0);
						c.refresh(fp.getProjectionEdges());
						// TODO: need multiple focal points
						c.setFocalPoint(motionModel.getStartRealization(),
								Color.red);
						c.setFocalPoint(motionModel.getEndRealization(),
								Color.blue);
						fp.addCurve3D(c);
					} else {
						assert (isTracingMotion() == 3);
						Curve3DMotion c1 = new Curve3DMotion(1);
						Curve3DMotion c2 = new Curve3DMotion(2);
						Edge es[] = fp.getProjectionEdges();
						c1.refresh(es);
						c1.setFocalPoint(motionModel.getNearestG1(), Color.pink);
						c2.refresh(es);
						c2.setFocalPoint(motionModel.getNearestG2(), Color.cyan);
						fp.addCurve3D(c1);
						fp.addCurve3D(c2);
					}
					fp.showCurve3DPanel();
				} else {
					fp.closeCurve3DPanel();
					fp.clearCurves();
				}
			}
		});

		// TODO:

		// showComponentsCurveBtn = new
		// JToggleButton("show curves of all components");
		// showComponentsCurveBtn.addActionListener(new ActionListener() {
		// public void actionPerformed(ActionEvent e) {
		// FloatingPanels fp = FloatingPanels.getInstance();
		//
		// if (showComponentsCurveBtn.isSelected()) {
		// //if (isTracingMotion() == 1) {
		// Curve3DMotion c = new Curve3DMotion(0);
		// c.refresh(fp.getProjectionEdges());
		//
		// fp.addCurve3D(c);
		// //} else {
		// assert (isTracingMotion() == 3);
		// Curve3DMotion c1 = new Curve3DMotion(1);
		// Curve3DMotion c2 = new Curve3DMotion(2);
		// Edge es[] = fp.getProjectionEdges();
		// c1.refresh(es);
		// c1.setFocalPoint(motionModel.getNearestG1());
		// c2.refresh(es);
		// c2.setFocalPoint(motionModel.getNearestG2());
		// fp.addCurve3D(c1);
		// fp.addCurve3D(c2);
		// }
		// fp.showCurve3DPanel();
		// } else {
		// fp.closeCurve3DPanel();
		// fp.clearCurves();
		// }
		// }
		// });

		pathSpinner = new DialSpinner<SamplePoint<Realization>>(
				motionModel.getSpinnerModel());
		pathSpinner.addChangeListener(new MotionSpinnerListener(0));

		pathSpinner1 = new DialSpinner<SamplePoint<Realization>>(
				motionModel.getSpinnerModel1());
		pathSpinner1.setCyclic(true);
		pathSpinner1.addChangeListener(new MotionSpinnerListener(1));

		pathSpinner2 = new DialSpinner<SamplePoint<Realization>>(
				motionModel.getSpinnerModel2());
		pathSpinner2.setCyclic(true);
		pathSpinner2.addChangeListener(new MotionSpinnerListener(2));

		componentSpinner = new DialSpinner<ConnectedComponent>(
				motionModel.getComponentSpinnerModel());
		componentSpinner.setCyclic(true);
		componentSpinner.addChangeListener(motionModel);

		traceVertexBtn = new JToggleButton("trace a vertex");
		traceVertexBtn.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				JToggleButton b = (JToggleButton) e.getSource();
				if (b.isSelected())
					;// main.graphPanel.setTracing(true);
				else {
					// main.graphPanel.setTracing(false);
					b.setText("trace a vertex");
				}
			}
		});

		pathPanel.add(genComponentBtn);
		pathPanel.add(realizationPathBtn);
		pathPanel.add(cayleyPathBtn);
		pathPanel.add(startBtn);
		pathPanel.add(endBtn);
		pathPanel.add(genPathBtn);
		pathPanel.add(pathSpinner);
		pathPanel.add(pathSpinner1);
		pathPanel.add(pathSpinner2);
		pathPanel.add(showComponentCurveBtn);
		pathPanel.add(traceVertexBtn);
		pathPanel.add(componentSpinner);

		add(pathPanel);
	}

	/*
	 * Let the user pick 3 base non-edges for component curve projection
	 * 
	 * private Edge[] pickEdges() { // TODO: let the user pick the edges.
	 * Iterator<Edge> ei = tdModel.getCanonicalBaseNonedges().iterator(); Edge
	 * e1 = ei.next(); Edge e2 = ei.next(); Edge e3 = ei.next(); Edge[] es = {
	 * e1, e2, e3 }; return es; }
	 */

	/**
	 * Show path control after gen CCS.
	 */
	private void showPathControl() {
		pathPanel.setVisible(true);

		for (Component c : pathPanel.getComponents()) {
			c.setVisible(false);
			if (c instanceof AbstractButton)
				((AbstractButton) c).setSelected(false);
			c.setEnabled(true);
		}

		realizationPathBtn.setVisible(true);
		if (!isOrientedCCS()) {
			genComponentBtn.setVisible(true);
			cayleyPathBtn.setVisible(true);
		}

	}

	/**
	 * @param b
	 *            the button clicked
	 * @return 1: start button; 2: end button; 0: dunno
	 */
	public int getPickingEnd(AbstractButton b) {
		if (b == startBtn)
			return 1;
		else if (b == endBtn)
			return 2;
		else
			return 0;
	}

	/**
	 * @return 1: path between cayleys; 2: path between realizations; 0: not
	 *         picking
	 */
	public int getPathType() {
		if (cayleyPathBtn.isSelected())
			return 1;
		else if (realizationPathBtn.isSelected())
			return 2;
		else
			return 0;
	}

	/**
	 * Enable genPathBtn if both start & end are picked
	 */
	public void checkGenPathBtn() {
		if (MotionModel.getInstance().getStartCayleyConfig() != null
				&& MotionModel.getInstance().getEndCayleyConfig() != null)
			genPathBtn.setVisible(true);
	}

}
