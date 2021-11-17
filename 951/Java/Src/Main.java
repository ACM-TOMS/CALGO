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

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JOptionPane;

import ccs.graph.Realization;
import ccs.graph.RealizationType;

public class Main extends JFrame {

	TDLinkageModel tdModel = TDLinkageModel.getInstance();

	// panels
	private GPanel gpanel;
	private ControlPanel controls;
	private CCSPanel ccsPanel;

	public static void main(String args[]) {
		Main main = new Main();
		main.init();
	}

	// GPanel graphPanel;

	/*
	 * 
	 * case generated: graph.resetInitPoints();
	 * 
	 * drawBtn.setVisible(true); generateBtn.setVisible(true);
	 * 
	 * td.generateCCS(); propagateStatus(s);
	 * 
	 * spacePanel.repaint(); break;
	 * 
	 * case picking: drawBtn.setVisible(false); generateBtn.setVisible(false);
	 * propagateStatus(s); break;
	 * 
	 * case tracing: //assert (oldStatus == Status.picking); propagateStatus(s);
	 * break;
	 * 
	 * case tracingNearest: propagateStatus(s); break; } }
	 */

	public void init() {

		// ================= Panels ====================
		gpanel = GPanel.getInstance();
		gpanel.setPreferredSize(new Dimension(750, 550));
		gpanel.setBackground(Color.white);
		gpanel.setMain(this);
		this.add(gpanel, BorderLayout.CENTER);

		controls = ControlPanel.getInstance();
		controls.init();
		controls.setMain(this);
		controls.setPreferredSize(new Dimension(260, 550));
		this.add(controls, BorderLayout.LINE_END);

		ccsPanel = CCSPanel.getInstance();
		ccsPanel.setPreferredSize(new Dimension(1005, 768 - 550 - 80));
		this.add(ccsPanel, BorderLayout.PAGE_END);

		this.setSize(1024, 768);
		this.setBackground(Color.lightGray);
		this.setVisible(true);
		// this.setStatus(Status.drawn);

		// test();
	}

	private void test() {
		final JButton test = new JButton("test");
		test.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				// TODO Auto-generated method stub
				String s = JOptionPane.showInputDialog(null, "soltype", "", 1);
				double length = Double.parseDouble(JOptionPane.showInputDialog(
						null, "Enter base.length : ", "", 1));
				String a[] = s.split(" ");
				RealizationType t = new RealizationType(a.length);
				for (int i = 0; i < a.length; ++i) {
					t.setOrientation(i, Integer.parseInt(a[i]));
				}
				Realization g = TDLinkageModel.getInstance().getTd()
						.solve(length, t);
				if (g != null) {
					TDLinkageModel.getInstance().setPoints(g);
					GPanel.getInstance().repaint();
				} else {
					test.setText("null");
				}
			}

		});
		gpanel.add(test);
	}
	/*
	 * private void test() { JPanel testPanel = new JPanel();
	 * 
	 * JButton test = new JButton("test solve"); test.addActionListener(new
	 * ActionListener() { public void actionPerformed(ActionEvent e) { double
	 * length = Double.parseDouble(JOptionPane.showInputDialog( null,
	 * "Enter base.length : ", "", 1)); ArrayList<Graph> solutions =
	 * td.solve(length); ((JButton) (e.getSource())).setText(solutions.size() +
	 * "");
	 * 
	 * // TreeDecomp temp = td; Graph temp = td.graph; assert (graph ==
	 * td.graph); // td = null; for (int i = 0; i < solutions.size(); ++i) {
	 * graph = solutions.get(i); // td.graph = graph;
	 * graphPanel.paintImmediately(graphPanel.getBounds()); try { synchronized
	 * (this) { this.wait(1000); } } catch (InterruptedException ee) {
	 * ee.printStackTrace(); } } // td.graph = temp; graph = temp;
	 * graphPanel.repaint(); } }); testPanel.add(test);
	 * 
	 * JButton testf = new JButton("test solve in fwd");
	 * testf.addActionListener(new ActionListener() { public void
	 * actionPerformed(ActionEvent e) { double length =
	 * Double.parseDouble(JOptionPane.showInputDialog( null,
	 * "Enter base.length : ", "", 1)); Graph g = td.tryRealize(length,
	 * td.getForwardSolutionType()); ((JButton) (e.getSource())).setText((g !=
	 * null) + "");
	 * 
	 * if (g != null) { // TreeDecomp temp = td; Graph temp = td.graph; assert
	 * (graph == td.graph); // td = null; // for (int i = 0; i <
	 * solutions.size(); ++i) { graph = g; // td.graph = graph;
	 * graphPanel.paintImmediately(graphPanel.getBounds()); try { synchronized
	 * (this) { this.wait(1000); } } catch (InterruptedException ee) {
	 * ee.printStackTrace(); } // } // td.graph = temp; graph = temp;
	 * graphPanel.repaint(); } } }); testPanel.add(testf);
	 * 
	 * // JButton test2 = new JButton("test orient"); //
	 * test2.addActionListener(new ActionListener() { // public void
	 * actionPerformed(ActionEvent e) { // td.generateConstructionSequence(); //
	 * for (ConstructionStep s : td.constructionSequence) //
	 * System.out.println(td.orientationOf(s) + ":" + s); // } // }); //
	 * controlPanel.add(test2);
	 * 
	 * // JButton testSetBase = new JButton("test set base"); //
	 * testSetBase.addActionListener(new ActionListener() { // public void
	 * actionPerformed(ActionEvent e) { // int i1 =
	 * Integer.parseInt(JOptionPane.showInputDialog(null, // "Enter edge.v1 : ",
	 * "", 1)); // int i2 = Integer.parseInt(JOptionPane.showInputDialog(null,
	 * // "Enter edge.v2 : ", "", 1)); // ((JButton) (e.getSource())).setText(""
	 * // + td.setBaseNonedge(new Edge(graph.getVertex(i1), graph //
	 * .getVertex(i2)))); // graphPanel.repaint(); // } // });
	 * 
	 * // JButton testEG = new JButton("Text EG"); //
	 * testEG.addActionListener(new ActionListener() { // public void
	 * actionPerformed(ActionEvent e) { // int i =
	 * Integer.parseInt(JOptionPane.showInputDialog(null, //
	 * "?th extrem graph: ", "", 1)); // TreeDecomp t = td.getExtremeGraph(i);
	 * // // ListGraph temp = graph; // graph = t.graph; // td = t; //
	 * graphPanel.repaint(); // // graph = temp; // } // }); //
	 * controlPanel.add(testEG); controlPanel.add(testPanel); }
	 */

}
