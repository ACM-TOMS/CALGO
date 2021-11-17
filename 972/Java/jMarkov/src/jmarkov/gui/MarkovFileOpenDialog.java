/**
 * 
 */
package jmarkov.gui;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.BoxLayout;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.WindowConstants;
import javax.swing.filechooser.FileFilter;

import jmarkov.DebugReporter;

/**
 * @author Germán Riaño. Universidad de los Andes.
 * 
 */
public class MarkovFileOpenDialog extends JDialog implements ActionListener {

    private static final long serialVersionUID = 19692005;

    private JPanel pnlMain = null;
    private JPanel pnlMP = null;
    private JPanel pnlEventFile = null;
    private JPanel pnlStateFile = null;

    private JTextField txtMarkovProcessFile = null;
    private JTextField txtEventFile = null;
    private JTextField txtStateFile = null;

    private JLabel lblMarkovProcessFile = null;
    private JLabel lblStateFile = null;
    private JLabel lblEventFile = null;

    private JButton btnMarkovProcess = null;
    private JButton btnStateFile = null;
    private JButton btnEventFile = null;
    private JLabel lblIntro = null; // @jve:decl-index=0:visual-constraint="567,12"
    private JPanel pnlControls = null;
    private JButton btnOK = null;
    private JButton btnCancel = null;

    private JButton[] btnsArray = null;
    private JTextField[] txtfldArray = null;
    private static String fileNames[] = null;

    /**
     * The result types returned by the File open dialog.
     * 
     * @author Germán Riaño. Universidad de los Andes.
     */
    public enum ResultTypes {
        /** The user hit OK */
        OK,
        /** The user hit cancel */
        CANCEL,
        /** The user closed the dialog hitting the (X) */
        CLOSE
    };

    private ResultTypes dialogResult = ResultTypes.CANCEL;

    private DebugReporter reporter = null;

    /**
     * This is the default constructor
     */
    public MarkovFileOpenDialog() {
        super();
        setModal(true);
        this.setDefaultCloseOperation(WindowConstants.HIDE_ON_CLOSE);
        initialize();
    }

    /**
     * Constructor pointing to a given Reporter and GUI.
     * @param reporter
     * @param owner
     */
    public MarkovFileOpenDialog(DebugReporter reporter, MarkovGUI owner) {
        super(owner, "JMarkov - Choose File", true);
        this.reporter = reporter;
        this.setDefaultCloseOperation(WindowConstants.HIDE_ON_CLOSE);
        if (fileNames == null) {
            fileNames = new String[] { "MyMarkovProcessEvent.class",
                    "MyMarkovProcessState.class", "MyMarkovProcess.class", };
        }
        initialize();
    }

    /**
     * Constructor with given file candidates (usually the last thre files used).
     * @param files
     * @param reporter
     * @param owner
     */
    public MarkovFileOpenDialog(String[] files, DebugReporter reporter,
            MarkovGUI owner) {
        this(reporter, owner);
        fileNames = files;
        for (int i = 0; i < 3; i++) {
            txtfldArray[i].setText(fileNames[i]);
        }
    }

    /**
     * This method initializes this dialog
     * 
     */
    private void initialize() {
        this.setSize(558, 298);
        this.setContentPane(getPnlMain());
        this.setLocationByPlatform(true);
        btnsArray = new JButton[] { btnEventFile, btnStateFile,
                btnMarkovProcess };
        txtfldArray = new JTextField[] { txtEventFile, txtStateFile,
                txtMarkovProcessFile };
        if (fileNames != null) {
            for (int i = 0; i < 3; i++) {
                txtfldArray[i].setText(fileNames[i]);
            }
        }
        SwingUtilities.updateComponentTreeUI(this);
        this.pack();
    }

    /*
     * Interface procedures
     * 
     */

    /**
     * Gets the result (OK or Cancel)
     * @return The result according to what the user hit.
     */
    public ResultTypes getDialogResult() {
        return dialogResult;
    }

    /**
     * Returns a String array with the names of the three files.
     * 
     * @return a String array with the names of the three files
     */
    public String[] getFileNames() {
        fileNames = new String[3];
        for (int i = 0; i < 3; i++) {
            fileNames[i] = txtfldArray[i].getText();
        }
        return fileNames;
    }

    /*
     * Button Events.
     */

    JFileChooser fc = null;

    /**
     * Open the File selection dialog
     * 
     */
	public void actionPerformed(ActionEvent e) {
        JButton btn = ((JButton) e.getSource());
        int idx = 0;
        for (idx = 0; btn != btnsArray[idx] && idx < 3; idx++)
            ;
        if (idx == 3)
            return;// button not found.
        String fullName = "";
        File file = null;
        if (fc == null) {
            fc = new JFileChooser();
            fc.setFileFilter(new FileFilter() {
                @Override
                public boolean accept(File f) {
                    String path = f.getAbsolutePath();
                    return (path.endsWith("class") || f.isDirectory());
                }

                @Override
                public String getDescription() {
                    return "Class Files (*.class) ";
                }
            });
        }
        File currFile = new File(txtfldArray[idx].getText());
        String path = currFile.getParent();
        String fileName = currFile.getName();
        String userDir = System.getProperty("user.dir");
        if (path == "" || path == null) {
            path = userDir;
        }
        fc.setSelectedFile(new File(path, fileName));
        int returnVal = fc.showOpenDialog(this);
        if (returnVal == JFileChooser.APPROVE_OPTION) {
            try {
                file = fc.getSelectedFile();
                if (file.getParent() == userDir) {
                    fullName = file.getName();
                } else {
                    fullName = file.getCanonicalPath();
                }
                txtfldArray[idx].setText(fullName);
            } catch (Exception ex) {
                reporter.debug(0, "Error openning file " + fullName
                        + "\nException:" + ex);
            }
        }
    }

    /**
     * This method initializes jContentPane
     * 
     * @return javax.swing.JPanel
     */
    private JPanel getPnlMain() {
        if (pnlMain == null) {
            lblIntro = new JLabel();
            lblIntro
                    .setText("Please select the file where the new classes reside:");
            lblIntro
                    .setFont(new java.awt.Font("Dialog", java.awt.Font.BOLD, 14));
            lblIntro.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
            lblIntro
                    .setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
            pnlMain = new JPanel();
            pnlMain.setLayout(new BoxLayout(pnlMain, BoxLayout.Y_AXIS));
            pnlMain.add(lblIntro, null);
            pnlMain.add(getPnlMP(), null);
            pnlMain.add(getPnlStateFile(), null);
            pnlMain.add(getPnlEventFile(), null);
            pnlMain.add(getPnlControls(), null);
        }
        return pnlMain;
    }

    /**
     * This method initializes jPanel
     * 
     * @return javax.swing.JPanel
     */
    private JPanel getPnlMP() {
        if (pnlMP == null) {
            BorderLayout borderLayout5 = new BorderLayout();
            borderLayout5.setHgap(9);
            borderLayout5.setVgap(9);
            lblMarkovProcessFile = new JLabel();
            lblMarkovProcessFile.setText("SimpleMarkovProcess File:");
            lblMarkovProcessFile
                    .setHorizontalAlignment(javax.swing.SwingConstants.LEADING);
            pnlMP = new JPanel();
            pnlMP.setLayout(borderLayout5);
            pnlMP.setBorder(javax.swing.BorderFactory.createEmptyBorder(5, 5,
                    5, 5));
            pnlMP.add(getTxtMarkovProcessFile(), java.awt.BorderLayout.CENTER);
            pnlMP.add(getBtnMarkovProcess(), java.awt.BorderLayout.EAST);
            pnlMP.add(lblMarkovProcessFile, java.awt.BorderLayout.NORTH);
        }
        return pnlMP;
    }

    /**
     * This method initializes jTextField
     * 
     * @return javax.swing.JTextField
     */
    private JTextField getTxtMarkovProcessFile() {
        if (txtMarkovProcessFile == null) {
            txtMarkovProcessFile = new JTextField();
            txtMarkovProcessFile.setText("MyMarkovProcess.class");
            txtMarkovProcessFile
                    .setHorizontalAlignment(SwingConstants.LEFT);
        }
        return txtMarkovProcessFile;
    }

    /**
     * This method initializes jTextField1
     * 
     * @return javax.swing.JTextField
     */
    private JTextField getTxtEventFile() {
        if (txtEventFile == null) {
            txtEventFile = new JTextField();
            txtEventFile.setText("MyEventProcessEvent.class");
        }
        return txtEventFile;
    }

    /**
     * This method initializes jButton
     * 
     * @return javax.swing.JButton
     */
    private JButton getBtnMarkovProcess() {
        if (btnMarkovProcess == null) {
            btnMarkovProcess = new JButton();
            btnMarkovProcess.setText("Browse ...");
            btnMarkovProcess.setIcon(new ImageIcon(getClass().getResource(
                    "/jmarkov/gui/images/open.gif")));
            btnMarkovProcess.addActionListener(this);
        }
        return btnMarkovProcess;
    }

    /**
     * This method initializes jTextField2
     * 
     * @return javax.swing.JTextField
     */
    private JTextField getTxtStateFile() {
        if (txtStateFile == null) {
            txtStateFile = new JTextField();
            txtStateFile.setText("MyMarkovProcessState.class");
        }
        return txtStateFile;
    }

    /**
     * This method initializes jPanel1
     * 
     * @return javax.swing.JPanel
     */
    private JPanel getPnlStateFile() {
        if (pnlStateFile == null) {
            BorderLayout borderLayout6 = new BorderLayout();
            borderLayout6.setHgap(9);
            borderLayout6.setVgap(9);
            lblStateFile = new JLabel();
            lblStateFile.setText("State File:");
            lblStateFile
                    .setHorizontalAlignment(javax.swing.SwingConstants.LEADING);
            lblStateFile
                    .setHorizontalTextPosition(javax.swing.SwingConstants.RIGHT);

            pnlStateFile = new JPanel();
            pnlStateFile.add(getTxtStateFile(), java.awt.BorderLayout.CENTER);
            pnlStateFile.setBorder(javax.swing.BorderFactory.createEmptyBorder(
                    5, 5, 5, 5));
            pnlStateFile.setLayout(borderLayout6);
            pnlStateFile.add(getTxtStateFile(), java.awt.BorderLayout.CENTER);
            pnlStateFile.add(getBtnStateFile(), java.awt.BorderLayout.EAST);
            pnlStateFile.add(lblStateFile, java.awt.BorderLayout.NORTH);
        }
        return pnlStateFile;
    }

    /**
     * This method initializes jPanel2
     * 
     * @return javax.swing.JPanel
     */
    private JPanel getPnlEventFile() {
        if (pnlEventFile == null) {
            BorderLayout borderLayout7 = new BorderLayout();
            borderLayout7.setHgap(9);
            borderLayout7.setVgap(9);
            lblEventFile = new JLabel();
            lblEventFile.setText("Event File:");
            lblEventFile
                    .setHorizontalAlignment(javax.swing.SwingConstants.LEADING);
            pnlEventFile = new JPanel();
            pnlEventFile.add(getTxtEventFile(), java.awt.BorderLayout.CENTER);
            pnlEventFile.setBorder(javax.swing.BorderFactory.createEmptyBorder(
                    5, 5, 5, 5));
            pnlEventFile.setLayout(borderLayout7);
            pnlEventFile.add(getTxtEventFile(), java.awt.BorderLayout.CENTER);
            pnlEventFile.add(getBtnEventFile(), java.awt.BorderLayout.EAST);
            pnlEventFile.add(lblEventFile, java.awt.BorderLayout.NORTH);
        }
        return pnlEventFile;
    }

    /**
     * This method initializes jButton
     * 
     * @return javax.swing.JButton
     */
    private JButton getBtnStateFile() {
        if (btnStateFile == null) {
            btnStateFile = new JButton();
            btnStateFile.setText("Browse ...");
            btnStateFile.setIcon(new ImageIcon(getClass().getResource(
                    "/jmarkov/gui/images/open.gif")));
            btnStateFile.addActionListener(this);
        }
        return btnStateFile;
    }

    /**
     * This method initializes jButton1
     * 
     * @return javax.swing.JButton
     */
    private JButton getBtnEventFile() {
        if (btnEventFile == null) {
            btnEventFile = new JButton();
            btnEventFile.setText("Browse ...");
            btnEventFile.setIcon(new ImageIcon(getClass().getResource(
                    "/jmarkov/gui/images/open.gif")));
            btnEventFile.addActionListener(this);
        }
        return btnEventFile;
    }

    /**
     * This method initializes jPanel
     * 
     * @return javax.swing.JPanel
     */
    private JPanel getPnlControls() {
        if (pnlControls == null) {
            pnlControls = new JPanel();
            pnlControls.add(getBtnOK(), null);
            pnlControls.add(getBtnCancel(), null);
        }
        return pnlControls;
    }

    /**
     * This method initializes jButton
     * 
     * @return javax.swing.JButton
     */
    private JButton getBtnOK() {
        if (btnOK == null) {
            btnOK = new JButton();
            btnOK.setText("OK");
            btnOK.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
                    dialogResult = ResultTypes.OK;
                    setVisible(false);
                }
            });
        }
        return btnOK;
    }

    /**
     * This method initializes jButton1
     * 
     * @return javax.swing.JButton
     */
    private JButton getBtnCancel() {
        if (btnCancel == null) {
            btnCancel = new JButton();
            btnCancel.setText("Cancel");
            btnCancel.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
                    dialogResult = ResultTypes.CANCEL;
                    setVisible(false);
                }
            });
        }
        return btnCancel;
    }

} // @jve:decl-index=0:visual-constraint="35,19"
