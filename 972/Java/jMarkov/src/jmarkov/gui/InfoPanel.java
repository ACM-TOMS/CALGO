/*
 * Created on Jul 27, 2003
 *
 * 
 */
package jmarkov.gui;

import java.awt.Component;

import javax.swing.JPanel;

import jmarkov.MarkovProcess;

/**
 * This class is the superclass of all panels that give info about the Markov
 * Chain.
 * 
 * @author German
 */
public abstract class InfoPanel extends JPanel {

	
	/**
	 * The Markov Process associated with this GUI.
	 */
	protected MarkovProcess mp = null;

	/**
	 *  
	 */
	public InfoPanel() {
		super();
	}

	/**
	 * Sets the MarkovProcess and shows the corresponding information. This method
	 * calls the method updateMP, which must be implemented by the subclasses.
	 * 
	 * @param mp
	 *          The MarkovProcess to show.
	 */
	public final void setMP(MarkovProcess mp) {
		// if (mp==this.mp)return;
		this.mp = mp;
		updateMP();
		Component[] comps = this.getComponents();
		for (int i = 0, n = comps.length; i < n; i++) {
			try {
				((InfoPanel) comps[i]).setMP(mp);
			} catch (Exception e) {
			}
		}
	}

	/**
	 * Unloads the Markov Process.
	 */
	public final void unloadMP() {
		mp = null;
		updateMP();
		Component[] comps = this.getComponents();
		for (int i = 0, n = comps.length; i < n; i++) {
			try {
				((InfoPanel) comps[i]).unloadMP();
			} catch (Exception e) {
			}
		}
	}

	/**
	 * This method is called whenever the SimpleMarkovProcess associated with this Panel
	 * has changed. The classes must implement this to show the relevant
	 * information. Note that <code>mp</code> might be null.
	 */
	protected abstract void updateMP();

	/**
	 * Return the SimpleMarkovProcess associted with this GUI.
	 * 
	 * @return Th MP.
	 */
	public MarkovProcess getMP() {
		return mp;
	}

}