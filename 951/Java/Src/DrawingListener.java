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

import javax.swing.BorderFactory;
import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.JToggleButton;

/**
 * When click a draw option button, automatically unselect other draw option
 * buttons, and clear drawing buffer of GPanel
 * 
 */
public class DrawingListener implements ActionListener {

	private JToggleButton drawOptionBtns[];

	public DrawingListener(JToggleButton drawOptionBtns[]) {
		this.drawOptionBtns = drawOptionBtns;
	}

	public void actionPerformed(ActionEvent e) {
		for (JToggleButton b : drawOptionBtns) {
			if (b != e.getSource())
				b.setSelected(false);
		}

		GPanel.getInstance().clearDrawingBuffer();
	}
}