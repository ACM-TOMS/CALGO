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
import java.util.List;

import javax.swing.JComponent;
import javax.swing.JSpinner;
import javax.swing.SpinnerListModel;
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

/**
 * A JSpinner using DialEditor and DialListModel.
 * 
 * @param <E>
 *            type of DialListModel's list element
 */
public class DialSpinner<E> extends JSpinner {
	/*
	 * public DialSpinner(DialListModel<E> model) { super(model);
	 * this.setEditor(new DialEditor(this)); // this.setBorder(null); }
	 */

	public DialSpinner() {
		super(new PercentageListModel<E>());
		this.setEditor(new DialEditor(this));
	}

	public DialSpinner(PercentageListModel<E> model) {
		super(model);
		this.setEditor(new DialEditor(this));
	}

	public void setList(List<E> l) {
		getModel().setList(l);
	}

	/*
	 * public Double getDoubleValue() { return (Double) getValue(); }
	 */

	public double getPercentage() {
		return ((DialListModel<E>) getModel()).getPercentage();
	}

	@Override
	public DialListModel<E> getModel() {
		return (DialListModel<E>) super.getModel();
	}

	@Override
	public void setModel(SpinnerModel model) {
		if (!(model instanceof DialListModel<?>))
			return;
		super.setModel(model);
	}

	/*
	 * public void setMinimum(double minimum) { ((SpinnerNumberModel)
	 * this.getModel()).setMinimum(minimum); }
	 * 
	 * public void setMaximum(double maximum) { ((SpinnerNumberModel)
	 * this.getModel()).setMaximum(maximum); }
	 * 
	 * public void setStepSize(double stepSize) { ((SpinnerNumberModel)
	 * this.getModel()).setStepSize(stepSize); }
	 */

	@Override
	public E getValue() {
		// TODO Auto-generated method stub
		return (E) super.getValue();
	}

	public void fire() {
		// TODO Auto-generated method stub
		this.fireStateChanged();
	}

	public boolean isCyclic() {
		return getModel().isCyclic();
	}

	public void setCyclic(boolean cyclic) {
		getModel().setCyclic(cyclic);
	}

	public void setInitValue() {
		getModel().setInitValue();
	}
}

/**
 * Acting similarly as SpinnerListModel. Can get the current percentage in the
 * list, as well as change the list.
 * 
 * @param <E>
 *            the type of list
 */
interface DialListModel<E> extends SpinnerModel {
	/**
	 * @return current percentage in the list
	 */
	double getPercentage();

	/**
	 * Change the list used in model.
	 * 
	 * @param l
	 *            the new list to use.
	 */
	void setList(List<E> l);

	boolean isCyclic();

	void setCyclic(boolean cyclic);

	public void setInitValue();
}

class PercentageListModel<E> implements DialListModel<E> {
	// private List<E> list;
	private int curIndex = 0;
	private SpinnerListModel model = new SpinnerListModel();
	private boolean isCyclic = false;
	private E firstValue, lastValue;

	public boolean isCyclic() {
		return isCyclic;
	}

	public void setCyclic(boolean cyclic) {
		isCyclic = cyclic;
		/*
		 * if (!isCyclic() && cyclic) { isCyclic = cyclic; } else if (isCyclic()
		 * && !cyclic) { }
		 */
	}

	public PercentageListModel() {
		model = new SpinnerListModel();
	}

	public PercentageListModel(List<E> l) {
		model = new SpinnerListModel(l);
		firstValue = l.get(0);
		lastValue = l.get(l.size() - 1);
	}

	public void setList(List<E> l) {
		model.setList(l);
		firstValue = l.get(0);
		lastValue = l.get(l.size() - 1);
	}

	public double getPercentage() {
		return (double) curIndex / size();
	}

	public void addChangeListener(ChangeListener l) {
		model.addChangeListener(l);
	}

	public void removeChangeListener(ChangeListener l) {
		model.removeChangeListener(l);
	}

	public int size() {
		return model.getList().size();
	}

	public Object getNextValue() {
		Object value = model.getNextValue();
		if (value == null) {
			if (isCyclic()) {
				this.setValue(firstValue);
				//Debug.warnMsg("return to first value");
				return firstValue;
			} else {
				// does not change curIndex
				return null;
			}
		} else {
			curIndex++;
			return value;
		}
	}

	public Object getPreviousValue() {
		Object value = model.getPreviousValue();
		if (value == null) {
			if (isCyclic()) {
				this.setValue(lastValue);
				//Debug.warnMsg("return to last value");
				return lastValue;
			} else {
				// does not change curIndex
				return null;
			}
		} else {
			curIndex--;
			return value;
		}
	}

	public E getValue() {
		return (E) model.getValue();
	}

	/**
	 * Changes the curIndex together with value.
	 * 
	 * @see javax.swing.SpinnerModel#setValue(java.lang.Object)
	 */
	public void setValue(Object value) {
		model.setValue(value);
		curIndex = model.getList().indexOf(value);
	}

	public void setInitValue() {
		model.setValue(firstValue);
		curIndex = 0;
	}
}

class DialEditor extends JComponent implements ChangeListener {
	double percentDone = 0;

	public DialEditor(JSpinner spinner) {
		DialListModel<?> model = (DialListModel<?>) (spinner.getModel());
		percentDone = model.getPercentage();
		spinner.addChangeListener(this);

		Dimension size = new Dimension(35, 35);
		setMinimumSize(size);
		setPreferredSize(size);
	}

	public void paint(Graphics g) {
		int startAngle = 90;
		int doneAngle = (int) (percentDone * 360);

		g.setColor(getBackground());
		g.fillArc(3, 3, getSize().width - 8, getSize().height - 8, 0, 360);

		g.setColor(Color.MAGENTA);
		g.fillArc(3, 3, getSize().width - 8, getSize().height - 8, startAngle,
				doneAngle);
		// System.out.println(startAngle+" "+doneAngle);

		g.setColor(Color.black);
		g.drawArc(3, 3, getSize().width - 8, getSize().height - 8, 0, 360);
	}

	public void stateChanged(ChangeEvent e) {
		JSpinner spinner = (JSpinner) (e.getSource());
		DialListModel<?> model = (DialListModel<?>) (spinner.getModel());
		percentDone = model.getPercentage();
		// System.out.println(percentDone);
		repaint();
	}
}

// ============================

class DoubleSpinner extends JSpinner {
	final Dimension editorDimension = new Dimension(80, 15);

	public DoubleSpinner(double value, double minimum, double maximum,
			double stepSize) {
		super(new SpinnerNumberModel(value, minimum, maximum, stepSize));
		this.getEditor().setMinimumSize(editorDimension);
		this.getEditor().setPreferredSize(editorDimension);
	}

	public DoubleSpinner() {
		super(new SpinnerNumberModel(0, 0, 100, 0.1));
		this.getEditor().setMinimumSize(editorDimension);
		this.getEditor().setPreferredSize(editorDimension);
	}

	public Double getDoubleValue() {
		return (Double) getValue();
	}

	public void setMinimum(double minimum) {
		((SpinnerNumberModel) this.getModel()).setMinimum(minimum);
	}

	public void setMaximum(double maximum) {
		((SpinnerNumberModel) this.getModel()).setMaximum(maximum);
	}

	public void setStepSize(double stepSize) {
		((SpinnerNumberModel) this.getModel()).setStepSize(stepSize);
	}
}
