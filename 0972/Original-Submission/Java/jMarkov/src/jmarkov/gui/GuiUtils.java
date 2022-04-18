/*
 * Created on 16/07/2003
 */

package jmarkov.gui;

import java.awt.Container;
import java.awt.Window;
import java.net.URL;
import java.util.Hashtable;

import javax.swing.Action;
import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JTabbedPane;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.text.JTextComponent;

/**
 * GUI utilities for actions
 * 
 * @author Germán Riaño. Universidad de los Andes.
 * 
 */
public class GuiUtils {

    static Hashtable<String, Action> actions;

    /**
     * 
     */
    public GuiUtils() {
        super();
    }

    /**
     * Changes the current look and feel for the given Window.
     * 
     * @param lnfName
     *            Look/And/Feel name
     * @param win
     *            Window
     */
    public static void changeLook(String lnfName, Window win) {
        try {
            UIManager.setLookAndFeel(lnfName);
        } catch (Exception e) {
        }
        ;
        SwingUtilities.updateComponentTreeUI(win);
        win.pack();
    }

    /**
     * Utility Function to enable/disable the tab associated with this component
     * 
     * @param comp
     *            The component
     * @param enabl
     *            True or false.
     */
    public static void setTabEnabled(JComponent comp, boolean enabl) {
        Container papi = comp.getParent();
        Container nene = comp;
        int idx = -1;
        while (papi != null) {
            if (papi instanceof JTabbedPane) {
                idx = ((JTabbedPane) papi).indexOfComponent(nene);
                ((JTabbedPane) papi).setEnabledAt(idx, enabl);
                try {
                    if (!enabl)
                        ((TextPanel) nene).clear();
                } catch (Exception e) {
                }
            }
            nene = papi;
            papi = papi.getParent();
        }
    }

    /**
     * 
     * The following two methods allow us to find an action provided by the
     * editor kit by its name.
     * @param textComponent 
     */
    public static void createActionTable(JTextComponent textComponent) {
        actions = new Hashtable<String, Action>();
        Action[] actionsArray = textComponent.getActions();
        for (int i = 0; i < actionsArray.length; i++) {
            Action a = actionsArray[i];
            actions.put((String) a.getValue(Action.NAME), a);
        }
    }

    /**
     * Creates a fancy Label with border ans title around.
     * @param title
     * @param text
     * @return The Label
     */
    public static JLabel fancyLabel(String title, String text) {
        JLabel lbl = new JLabel(text);
        lbl.setBorder(BorderFactory.createCompoundBorder(BorderFactory
                .createTitledBorder(title), BorderFactory.createEmptyBorder(5,
                5, 5, 5)));
        return lbl;
    }

    /**
     * Returns the Action object with this name
     * @param name
     * @return Action object.
     */
    public static Action getActionByName(String name) {
        return (actions.get(name));
    }

    /*
     * Utility function to get icons
     */
    static ImageIcon createIcon(String location) {
        URL url = MarkovGUI.class.getResource(location);
        return new ImageIcon(url);
    }

    /**
     * Aux method: deletes text from a button if icon is available
     * @param a Action
     * @return Button object
     *
     */
    public static JButton getButton(Action a) {
        JButton but = new JButton();
        but.setAction(a);
        if (but.getIcon() != null)
            but.setText("");
        return but;
    }

}
