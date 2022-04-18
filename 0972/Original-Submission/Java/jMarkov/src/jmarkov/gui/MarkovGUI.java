/*
 * Created on 8/07/2003
 */
package jmarkov.gui;

import static jmarkov.MarkovProcess.Status.GENERATED;
import static jmarkov.MarkovProcess.Status.IDLE;
import static jmarkov.MarkovProcess.Status.NoModel;
import static jmarkov.MarkovProcess.Status.RUNNING;
import static jmarkov.MarkovProcess.Status.SUSPENDED;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.HeadlessException;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.io.File;
import java.io.PrintWriter;
import java.net.URL;
import java.text.DateFormat;
import java.util.Date;
import java.util.StringTokenizer;

import javax.swing.AbstractAction;
import javax.swing.AbstractButton;
import javax.swing.Action;
import javax.swing.ButtonGroup;
import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.JSpinner;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.JToolBar;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;
import javax.swing.Timer;
import javax.swing.UIManager;
import javax.swing.WindowConstants;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.text.DefaultEditorKit;

import jmarkov.DebugReporter;
import jmarkov.MarkovProcess;
import jmarkov.MarkovProcess.Status;

/**
 * @author Germán Riaño. Universidad de los Andes.
 */
public class MarkovGUI extends JFrame implements ChangeListener,
        ComponentListener, ActionListener {

    private static final long serialVersionUID = 1969;

    MarkovProcess<?, ?> mp = null;

    MarkovClassLoader cld = null; // used to load classes
                                    // dynamically

    // Main GUIs:
    private JPanel mainContentPane = null; // @jve:decl-index=0:visual-constraint="30,186"
    private JSplitPane splitPane = null;;
    private JPanel statusBar = null;
    private JToolBar toolbar = null;
    private TextPanel outTextPanel = null;
    private TextPanel dbgTextPanel = null;
    private JTabbedPane tbpTop = null;
    // On the status bar:
    private JLabel lblStatus = null;
    private JProgressBar progBar = null;
    // menu stuff
    private JMenuBar menuBar = null;
    private JMenu controlMenu = null;
    private JMenu editMenu = null;
    private JMenu optionsMenu = null;
    // mainTabPane tabbed panelcomponents
    private JPanel mainTabPane = null;
    private InfoPanel pnlBrowse;
    private InfoPanel pnlStates;
    private InfoPanel pnlRates;
    private InfoPanel pnlMOPs;
    private InfoPanel pnlEvents;
    private InfoPanel pnlTransient;
    private TextPanel pnlDesc;
    private JLabel lblClass;
    // actions
    private Action accGo;
    private Action accStep;
    private Action accPause;
    private Action accStop;
    private Action accReStart;
    private Action accLoad;
    private Action accReload;
    private Action accUnload;
    private Action accExit;

    private JSpinner spnDebug;

    // Other elements
    private Timer timer = new Timer(500, this); // used for progress
                                                // bar
    private JMenu fileMenu = null;
    private DebugReporter reporter = null;
    private MarkovFileOpenDialog mfod = null;
    private String lastClassesLoaded[] = null;
    private boolean allowedToRun = true; // wheather the model can
                                            // continue

    // Constructors

    /**
     * Creates GUI with thr given SimpleMarkovProcess.
     * @param mp
     */
    public MarkovGUI(MarkovProcess mp) {
        this();
        loadMP(mp); // loads the MP instance
    }

    /**
     * Creates a default GUI.
     * @throws java.awt.HeadlessException
     */
    public MarkovGUI() throws HeadlessException {
        super("JMarkov -- No Model Loaded");
        initialize();
    }

    /**
     * This method initiallizes alll GUI Components
     */
    private void initialize() {
        initActions();
        mainContentPane = getMainContentPane();
        this.setContentPane(getMainContentPane());
        this.setLocationByPlatform(true);
        System.setProperty("java.awt.Window.locationByPlatform", "true");
        String os = System.getProperty("os.name");
        if (os.startsWith("Windows"))
            GuiUtils.changeLook(
                    "com.sun.java.swing.plaf.windows.WindowsLookAndFeel", this);
        this.setSize(800, 500);
        this.setPreferredSize(new java.awt.Dimension(800, 500));
        this.setJMenuBar(getTheMenuBar());
        this.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        outTextPanel.captureStandardOut();
        // dbgTextPanel.captureStandardErr();
        updateStatus();
        pack();
        setVisible(true);
    }

    // interface

    /**
     * @return Returns true if allowedToRun.
     */
    public boolean isAllowedToRun() {
        return allowedToRun;
    }

    // GUI Elements builders

    private JPanel getMainContentPane() {
        if (mainContentPane == null) {
            mainContentPane = new JPanel();
            mainContentPane.setLayout(new BorderLayout());
            outTextPanel = new TextPanel("Output");
            mainContentPane.add(getToolBar(), BorderLayout.NORTH);
            // Create a split pane for the change log and the output
            // area.
            mainContentPane.add(getSplitPane(), BorderLayout.CENTER);

            // Create the status area.
            mainContentPane.add(getStatusBar(), BorderLayout.SOUTH);
            // set up menus
            menuBar = getTheMenuBar();
            setJMenuBar(menuBar);
        }
        return mainContentPane;
    }

    /**
     * Sets up the status Bar
     */
    private JPanel getStatusBar() {
        if (statusBar == null) {
            statusBar = new JPanel();
            statusBar.setLayout(new GridLayout(1, 2, 10, 10));
            lblStatus = new JLabel();
            progBar = getProgressBar();

            statusBar.add(lblStatus);
            statusBar.add(progBar);
        }
        return statusBar;
    }

    /**
     * Sets up and retuns the Progress bar
     * @return th Progress bar
     */
    private JProgressBar getProgressBar() {
        if (progBar == null) {
            progBar = new JProgressBar();
            progBar.setStringPainted(true);
            progBar.setString("");
        }
        return progBar;
    }

    /* sets upp all the menus */
    private JMenuBar getTheMenuBar() {
        // Set up the menu bar.
        if (menuBar == null) {
            GuiUtils.createActionTable(outTextPanel.txtArea);
            editMenu = getEditMenu();
            controlMenu = getControlMenu();
            optionsMenu = getOptionsMenu();
            menuBar = new JMenuBar();
            menuBar.add(getFileMenu());
            menuBar.add(editMenu);
            menuBar.add(controlMenu);
            menuBar.add(optionsMenu);
        }
        return menuBar;
    }

    private JTabbedPane getUpperTabPane() {
        // Tabbed pane
        if (tbpTop == null) {
            tbpTop = new JTabbedPane();
            lblClass = GuiUtils.fancyLabel("Class:", "");
            mainTabPane = getMainTabPane();

            tbpTop.addTab("Main", null, mainTabPane, "Main Information");

            pnlBrowse = new SparseMatrixPanel(true);
            pnlBrowse.addComponentListener(this);
            tbpTop.addTab("Browse", null, pnlBrowse, "Browse the model");

            pnlStates = new TextPanel("States");
            pnlStates.addComponentListener(this);
            tbpTop.addTab("States", null, pnlStates, "States information");

            pnlRates = new MatrixPanel();
            pnlRates.addComponentListener(this);
            tbpTop.addTab("Rates", null, pnlRates,
                    "Transition Rates information");

            pnlMOPs = new TextPanel("Measures of Performance");
            pnlMOPs.addComponentListener(this);

            pnlTransient = new TextPanel("Transient Bahavior. "
                    + "Not Implemented yet.");

            tbpTop.addTab("MOPs", null, pnlMOPs, "Measures of performance");
            pnlEvents = new TextPanel("Events Rates");
            pnlEvents.addComponentListener(this);
            tbpTop.addTab("Events", null, pnlEvents, "Events Rates");
            tbpTop
                    .addTab("Transient", null, pnlTransient,
                            "Transient Analisys");
            tbpTop.addTab("Output", null, outTextPanel, "Output");
            tbpTop.setSelectedIndex(0);
            tbpTop.setPreferredSize(new java.awt.Dimension(500, 200));
        }
        return tbpTop;
    }

    private JPanel getMainTabPane() {
        if (mainTabPane == null) {
            mainTabPane = new JPanel();
            mainTabPane.setLayout(new BorderLayout());
            mainTabPane.add(lblClass, BorderLayout.NORTH);
            pnlDesc = new TextPanel("Description:");
            // pnlDesc = new TextPanel("Description:", "", 6, 110);
            mainTabPane.add(pnlDesc, BorderLayout.CENTER);
        }
        return mainTabPane;
    }

    private JSplitPane getSplitPane() {
        if (splitPane == null) {
            dbgTextPanel = new TextPanel("Debug Info");
            tbpTop = getUpperTabPane();
            splitPane = new JSplitPane();
            splitPane.setOrientation(JSplitPane.VERTICAL_SPLIT);
            splitPane.setTopComponent(tbpTop);
            splitPane.setBottomComponent(dbgTextPanel);
            splitPane.setOneTouchExpandable(true);
            splitPane.setResizeWeight(1);
        }
        return splitPane;
    }

    private JMenu getControlMenu() {
        controlMenu = new JMenu("Control");
        JMenuItem menuItem = null;

        Action[] actions = { accGo, accPause, accReStart, accStop };
        for (int i = 0; i < actions.length; i++) {
            menuItem = new JMenuItem(actions[i]);
            controlMenu.add(menuItem);
        }
        return controlMenu;
    }

    private JMenu getEditMenu() {
        editMenu = new JMenu("Edit");
        // These actions come from the default editor kit.
        // Get the ones we want and stick them in the menu.
        // menu.add(getActionByName(DefaultEditorKit.cutAction));
        editMenu.add(GuiUtils.getActionByName(DefaultEditorKit.copyAction));
        editMenu.addSeparator();
        editMenu
                .add(GuiUtils.getActionByName(DefaultEditorKit.selectAllAction));
        editMenu.add(GuiUtils
                .getActionByName(DefaultEditorKit.selectParagraphAction));
        return editMenu;
    }

    /**
     * Sets up and returns the tool bar.
     * @return the Tool bar
     */
    private JToolBar getToolBar() {
        if (toolbar == null) {
            toolbar = new JToolBar();

            /*
             * toolbar.add(getBtnGo()); btnPause = new
             * JButton("PAUSE"); btnPause.setAction(accPause);
             * toolbar.add(btnPause);
             */
            for (Commands cmd : Commands.values()) {
                toolbar.add(cmd.getButton());
            }
            Commands.PAUSE.getButton().setVisible(false);

            toolbar.addSeparator();
            toolbar.add(new JLabel("Debug Level:"));
            spnDebug = getSpinner();
            toolbar.add(spnDebug);

            toolbar.setFloatable(false);
        }
        return toolbar;
    }

    private JSpinner getSpinner() {
        if (spnDebug == null) {
            // spnDebug = new JSpinner(new SpinnerListModel(new
            // String[]
            // {"0","1","2","3","4"}));
            spnDebug = new JSpinner(new SpinnerNumberModel(0, 0, 4, 1));
            spnDebug.addChangeListener(this);
            spnDebug.setMaximumSize(new Dimension(50, 20));

        }
        return spnDebug;
    }

    private JMenu getOptionsMenu() {
        if (optionsMenu == null) {
            optionsMenu = new JMenu("Options");
            JMenu lookItem = new JMenu("Interface look");
            final JFrame frame = this;
            optionsMenu.add(lookItem);
            JRadioButtonMenuItem menuItem = null;
            UIManager.LookAndFeelInfo looks[] = UIManager
                    .getInstalledLookAndFeels();
            ButtonGroup group = new ButtonGroup();
            for (int i = 0; i < looks.length; i++) {
                menuItem = new JRadioButtonMenuItem(looks[i].getName());
                group.add(menuItem);
                lookItem.add(menuItem);
                menuItem.setActionCommand(looks[i].getClassName());
                menuItem.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        GuiUtils.changeLook(e.getActionCommand(), frame);
                        if (mfod == null)
                            GuiUtils.changeLook(e.getActionCommand(), mfod);
                    }
                });
                if (UIManager.getLookAndFeel().getName() == looks[i].getName())
                    menuItem.setSelected(true);
            }
        }
        return optionsMenu;
    }

    /**
     * Called to ensure that all the GUI shows correctly according to
     * Model status
     */
    public void updateStatus() {
        Status status = NoModel;
        if (mp != null)
            status = mp.getStatus();
        boolean goVisible = (status != RUNNING);
        Commands.GO.getButton().setVisible(goVisible);
        Commands.PAUSE.getButton().setVisible(!goVisible);

        accGo.setEnabled(goVisible && (status == IDLE || status == SUSPENDED));
        accPause.setEnabled((!goVisible) && (status == RUNNING));

        accStop.setEnabled(status == RUNNING);
        accStep.setEnabled(status == SUSPENDED || status == IDLE);
        accReStart.setEnabled((status != IDLE) && (status != RUNNING)
                && (status != NoModel));
        accLoad.setEnabled(mp == null);
        accReload.setEnabled(lastClassesLoaded != null);
        accUnload.setEnabled(mp != null);
        String msg = "No model loaded";
        if (mp != null)
            msg = "   " + mp.getStatusMsg();
        lblStatus.setText(msg);
        // setTabEnabled(pnlDesc, status != -1);
        setTabEnabled(pnlDesc, true);
        setTabEnabled(pnlStates, status == GENERATED);
        setTabEnabled(pnlMOPs, status == GENERATED);
        setTabEnabled(pnlRates, status == GENERATED);
        // setTabEnabled(pnlRates, true);
        setTabEnabled(pnlEvents, status == GENERATED);
        // TODO: change when transient is implemented.
        // setTabEnabled(pnlTransient, status ==
        // SimpleMarkovProcess.GENERATED);
        setTabEnabled(pnlTransient, false);
        if (status == IDLE) {
            // progBar.setIndeterminate(false);
        }
        if (status == RUNNING) {
            if (!timer.isRunning())
                timer.start();
        } else {
            progBar.setIndeterminate(false);
            progBar.setValue(progBar.getMinimum());
        }

        if (status == GENERATED) {
            progBar.setIndeterminate(false);
            progBar.setMaximum(mp.getNumStates());
        }
        if (status != NoModel)
            spnDebug.setValue(new Integer(mp.getDebugLevel()));
    }

    /*
     * COMMANDS
     */

    /**
     * Closes the program.
     */
    public void cmdExit() {
        timer.stop();
        System.exit(0);
    }

    /**
     * Generates the currently loaded model.
     */
    public void cmdGo() {
        Commands.GO.getButton().setVisible(false);
        Commands.PAUSE.getButton().setVisible(true);
        accGo.setEnabled(false);
        accPause.setEnabled(true);
        allowedToRun = true;
        timer.start();
        progBar.setIndeterminate(true);
        oneStep = false;
        runModel();
    }

    /**
     * Pauses generation of the current model.
     */
    public void cmdPause() {
        allowedToRun = false;
        mp.pause();
        progBar.setIndeterminate(false);
        timer.stop();
    }

    /**
     * Reloads the model. If the model has been recompiled it loads
     * the most recent version.
     */
    public void cmdReLoad() {
        timer.stop();
        progBar.setIndeterminate(false);
        progBar.setValue(progBar.getMinimum());
        dbgTextPanel.clear();
        pnlDesc.clear();
        lblClass.setText("");
        setTabEnabled(lblClass, true);
        tbpTop.setSelectedIndex(0);
        unLoadMP();
        loadMP(lastClassesLoaded);
    }

    /**
     * Re-starts the generation of the model.
     */
    public void cmdReStart() {
        timer.stop();
        progBar.setIndeterminate(false);
        progBar.setValue(progBar.getMinimum());
        mp.reset();
        dbgTextPanel.clear();
        tbpTop.setSelectedComponent(mainTabPane);
    }

    /**
     * Runs one step of model generation.
     */
    public void cmdStep() {
        timer.start();
        progBar.setIndeterminate(true);
        allowedToRun = true;
        oneStep = true;
        runModel();
    }

    /**
     * Activates the file load dialog.
     */
    private void loadMPFromFileDialog() {
        // Class[] classesArray = new Class[] { Event.class,
        // State.class,
        // SimpleMarkovProcess.class };
        if (mfod == null)
            mfod = new MarkovFileOpenDialog(reporter, this);
        mfod.setVisible(true);
        if (mfod.getDialogResult() == MarkovFileOpenDialog.ResultTypes.CANCEL) {
            reporter.debug(0, "Action Canceled");
            return;
        }
        String[] files = mfod.getFileNames();
        loadMP(files);
    }

    /**
     * Loads a model from the given classes
     * @param files The name of Event, State and Model classes, in
     *        that order.
     */
    public void loadMP(String files[]) {
        Object objects[] = new Object[3];
        File file = null;
        System.gc();
        if (cld == null)
            cld = new MarkovClassLoader();
        for (int i = 0; i < 3; i++) {
            try {
                if (files[i] != "") {
                    file = new File(files[i]);
                    // instantiates only MP class, which is i==2.
                    objects[i] = cld.loadFromFile(file, i == 2);
                }
            } catch (Exception e) {
                reporter.debug(0, e.getMessage());
                JOptionPane.showMessageDialog(this, e.getMessage(),
                        "JMarkov -- Error Opening model from files",
                        JOptionPane.ERROR_MESSAGE);
            }
        }
        if (objects[2] instanceof MarkovProcess) {
            mp = (MarkovProcess) objects[2];
        }
        if (mp != null) {
            loadMP(mp);
            lastClassesLoaded = files;
        }
    }

    /**
     * Loads the given SimpleMarkovProcess in the GUI
     * @param mp the markov Process to load.
     */
    public void loadMP(MarkovProcess<?, ?> mp) {
        this.mp = mp;
        int dbgLevel = mp.getDebugLevel();
        mp.setDebugReporter(this.getDebugReporter());
        reporter.setDebugLevel(dbgLevel); // recover previous Debug
                                            // Level

        String files[] = new String[3];
        files[0] = filePathForClass(mp.getEventClass());
        reporter.debug(1, "Event Class file loaded:" + files[0]);
        files[1] = filePathForClass(mp.getStateClass());
        reporter.debug(1, "State Class file loaded:" + files[1]);
        files[2] = filePathForClass(mp.getClass());
        reporter.debug(1, "Model Class file loaded:" + files[2]);
        if (mfod == null)
            mfod = new MarkovFileOpenDialog(files, reporter, this);
        if (files[2] != "")
            lastClassesLoaded = files;
        String className = mp.getClass().getName();
        this.setTitle("JMarkov -- " + className);
        try {
            File file = new File(files[2]);
            Date date = new Date(file.lastModified());
            lblClass.setText(mp.getClass().getName() + " ("
                    + file.getCanonicalPath() + " -- "
                    + DateFormat.getInstance().format(date) + ")");
        } catch (Exception e) {
            lblClass.setText(className);
        }
        pnlDesc.setText(mp.description());
        spnDebug.setValue(new Integer(mp.getDebugLevel()));
        tbpTop.setSelectedIndex(0);
        updateStatus();
    }

    /**
     * @param cls
     * @return filename where this class is located. returns "" if not
     *         found.
     */
    private String filePathForClass(Class cls) {
        String fileName = "";
        String canonicalName = cls.getCanonicalName();
        String simpleName = cls.getSimpleName();
        String sep = System.getProperty("file.separator");
        StringTokenizer st = new StringTokenizer(canonicalName, ".");
        String className = (st.hasMoreTokens()) ? st.nextToken() : "";
        while (st.hasMoreTokens()) {
            className += sep + st.nextToken();
        }
        try {
            URL url = cls.getClassLoader().getResource(className + ".class");
            File file = new File(url.toURI());
            String path = file.getParent() + sep;
            fileName = path + simpleName + ".class";
        } catch (Exception e) {
        }
        return fileName;
    }

    /**
     * Unload the SimpleMarkovProcess
     */
    public void unLoadMP() {
        if (mp != null)
            mp.clearMOPs();
        mp = null;
        pnlStates.unloadMP();
        pnlRates.unloadMP();
        pnlMOPs.unloadMP();
        pnlEvents.unloadMP();
        pnlTransient.unloadMP();
        pnlDesc.unloadMP();
        setTitle("JMarkov -- No model loaded");
        cld = null;
        System.gc();
    }

    /**
     * Calls the appropriate function so the Panel gets updated
     * @param aPanel The States, MOPs, etc panel.
     */
    public void updateTxtPanel(TextPanel aPanel) {
        if (mp != null) {
            final TextPanel panel = aPanel; // needs to be final to be
                                            // used
            final SwingWorker worker = new SwingWorker() {
                @Override
                public Object construct() {
                    TextPanel txtPanel = panel;
                    PrintWriter wtr = null;
                    wtr = new PrintWriter(txtPanel.getStream(), true);
                    txtPanel.setText("");
                    if (txtPanel == pnlStates) {
                        mp.printStates(wtr);
                    }
                    if (txtPanel == pnlMOPs) {
                        mp.printMOPs(wtr);
                    }
                    if (txtPanel == pnlEvents) {
                        mp.printEventsRates(wtr);
                    }
                    return null; // return value not used by this
                                    // program
                }

                // Runs on the event-dispatching thread.
                @Override
                public void finished() {
                    updateStatus();
                }
            };
            worker.start();
        }
    }

    /**
     * @return The GUI reporter for debug info.
     */
    public DebugReporter getDebugReporter() {
        if (reporter == null) {
            // load it with auto-flush = true
            reporter = new DebugReporter(new PrintWriter(dbgTextPanel
                    .getStream(), true));

            dbgTextPanel.setAutoShow(true);
        }
        return reporter;
    }

    /**
     * Runs the model. If oneStep is true, it runs just one step,
     * otherwise runs it all.
     */
    boolean oneStep = false;

    /**
     * Runs the model
     */
    public void runModel() {
        final SwingWorker runWorker = new SwingWorker() {
            @Override
            public Object construct() {
                if (oneStep) {
                    mp.goStep();
                } else
                    mp.go();
                return null;
            }

            @Override
            public void finished() {
                updateStatus();
            }
        };
        Runnable modelRunner = new Runnable() {
            public void run() {
                runWorker.start(); // required for SwingWorker 3
            }
        };
        SwingUtilities.invokeLater(modelRunner);
    }

    /**
     * Initialize all actions. Sets up icons.
     */
    private void initActions() {
        for (Commands cmd : Commands.values()) {
            cmd.setAction(new ButtonOrMenuUserAction(cmd.getName(), cmd
                    .getIcon(), cmd.getDescription()));
        }
        accGo = Commands.GO.getAction();
        accPause = Commands.PAUSE.getAction();
        accStep = Commands.STEP.getAction();
        accStop = Commands.STOP.getAction();
        accReStart = Commands.RESTART.getAction();
        accLoad = Commands.LOAD.getAction();
        accUnload = Commands.UNLOAD.getAction();
        accReload = Commands.RELOAD.getAction();
        accExit = Commands.EXIT.getAction();

    }

    /*
     * EVENTS HANDLERS
     */

    /*
     * (non-Javadoc)
     * @see javax.swing.event.ChangeListener#stateChanged(javax.swing.event.ChangeEvent)
     */
    public void stateChanged(ChangeEvent e) {
        int level = ((SpinnerNumberModel) ((JSpinner) e.getSource()).getModel())
                .getNumber().intValue();
        mp.setDebugLevel(level);
    }

    /*
     * @see java.awt.event.ComponentListener#componentHidden(java.awt.event.ComponentEvent)
     */
    public void componentHidden(ComponentEvent e) {

    }

    /*
     * @see java.awt.event.ComponentListener#componentMoved(java.awt.event.ComponentEvent)
     */
    public void componentMoved(ComponentEvent e) {
    }

    /*
     * @see java.awt.event.ComponentListener#componentResized(java.awt.event.ComponentEvent)
     */
    public void componentResized(ComponentEvent e) {
    }

    /*
     * @see java.awt.event.ComponentListener#componentShown(java.awt.event.ComponentEvent)
     */
    public void componentShown(ComponentEvent e) {
        if (e.getComponent() instanceof InfoPanel) {
            InfoPanel infoPanel = (InfoPanel) e.getComponent();
            infoPanel.setMP(mp);
        }
        if (e.getComponent() instanceof TextPanel) {
            TextPanel txtpanel = (TextPanel) e.getComponent();
            String currentText = null;
            currentText = txtpanel.getText();
            if ((currentText == null) || currentText.equals("")) {
                txtpanel.setText("Wait...");
                updateTxtPanel(txtpanel);
            }
        }
    }

    /**
     * Utility Function to enable/disable the tab associated with this
     * component
     * @param comp The component
     * @param enabl Troe or False.
     */
    public void setTabEnabled(JComponent comp, boolean enabl) {
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

    /*
     * (non-Javadoc) Called only by timer
     * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
     */
    public void actionPerformed(ActionEvent e) {
        if (e.getSource() instanceof AbstractButton) {
        } else {
            // Tmier action
            Status st = mp.getStatus();
            progBar.setValue((int) mp.getProgress());
            if (st != RUNNING) {
                timer.stop();
                if (st == GENERATED)
                    progBar.setValue(progBar.getMinimum());
            }
        }

    }

    /**
     * This method initializes jMenu
     * @return javax.swing.JMenu
     */
    private JMenu getFileMenu() {
        if (fileMenu == null) {
            fileMenu = new JMenu();
            fileMenu.setText("File");
            JMenuItem menuItem = null;
            Action[] actions = { accLoad, accUnload, accReload, accExit };
            for (Action ac : actions) {
                menuItem = new JMenuItem(ac);
                fileMenu.add(menuItem);
            }
        }
        return fileMenu;
    }

    // this class is internal to the Enumeration
    class ButtonOrMenuUserAction extends AbstractAction {

        // must serve some purpose?
        private static final long serialVersionUID = 1969;

        /**
         * Constructor.
         * @param name The name
         * @param icon
         * @param desc Description
         * @param mnemonic
         */
        public ButtonOrMenuUserAction(String name, Icon icon, String desc,
                Integer mnemonic) {
            super(name, icon);
            putValue(SHORT_DESCRIPTION, desc);
            putValue(MNEMONIC_KEY, mnemonic);
        }

        /**
         * @param name
         * @param icon
         * @param desc
         */
        public ButtonOrMenuUserAction(String name, Icon icon, String desc) {
            this(name, icon, desc, new Integer(0));
        }

        /*
         * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
         */
        public void actionPerformed(ActionEvent e) {

            // A button action
            Action acc = ((AbstractButton) e.getSource()).getAction();
            if (acc == accGo) {
                cmdGo();
            }
            if (acc == accPause) {
                cmdPause();
            }
            if (acc == accStep) {
                cmdStep();
            }
            if (acc == accStop) {
                cmdStop();
            }
            if (acc == accReStart) {
                cmdReStart();
            }
            if (acc == accLoad) {
                cmdLoad();
            }
            if (acc == accReload) {
                cmdReLoad();
            }
            if (acc == accUnload) {
                cmdUnload();
            }
            if (acc == accExit) {
                cmdExit();
            }
            updateStatus();
        }

        /**
         * 
         */
        public void cmdUnload() {
            timer.stop();
            progBar.setIndeterminate(false);
            progBar.setValue(progBar.getMinimum());
            dbgTextPanel.clear();
            pnlDesc.clear();
            lblClass.setText("No class loaded");
            setTabEnabled(lblClass, true);
            tbpTop.setSelectedIndex(0);
            unLoadMP();
        }

        /**
         * 
         */
        public void cmdLoad() {
            timer.stop();
            progBar.setIndeterminate(false);
            progBar.setValue(progBar.getMinimum());
            dbgTextPanel.clear();
            loadMPFromFileDialog();
        }

        /**
         * 
         */
        public void cmdStop() {
            timer.stop();
            progBar.setIndeterminate(false);
            progBar.setValue(progBar.getMinimum());
            allowedToRun = false;
            mp.pause();
            mp.reset();
            tbpTop.setSelectedComponent(mainTabPane);
        }
    };// end ButtonOrMenuUserAction class

    /**
     * Runs the user interface. If no argumentas are given it runs
     * with no associated model
     * @param a
     */
    public static void main(String[] a) {
        MarkovGUI theGUI = new MarkovGUI();
        theGUI.setVisible(true);
    }

    /**
     * Command enumeration encapsulates the functionallity of all
     * functions.
     * @author Germán Riaño. Universidad de los Andes. (C) 2005
     */
    enum Commands {
        /** Runs model. */
        GO("GO", "Play16.gif", "Generates the model"), //
        /** Puses Execution */
        PAUSE("PAUSE", "Pause16.gif", "Pauses Execution of the model"), //
        /** Runs one step */
        STEP("STEP", "StepForward16.gif", "Generates one step of the model"), //
        /** Re-starts execution */
        RESTART("Re-Start", "Rewind16.gif", "Starts over execution"), //
        /** Stops execution */
        STOP("Stop", "Stop16.gif", "Stops Execution"), //
        /** Load a new model */
        LOAD("Load", "Import16.gif", "Loads the problem"), //
        /** Unloads current model */
        UNLOAD("Unload", "Export16.gif", "Unloads the problem"), //
        /** Re-loads last model. */
        RELOAD("Reload", "Refresh16.gif", "Reloads the problem"), //
        /** Exits JMarkov. */
        EXIT("Exit", "Remove16.gif", "Exit JMarkov");

        private Action action;
        private JButton button;
        private ImageIcon icon = null;
        private String name;
        private String description;

        private Commands(String name, String iconFile, String description) {
            try {
                icon = GuiUtils.createIcon("/jmarkov/gui/images/" + iconFile);
            } catch (Exception e) {
            }
            this.name = name;
            this.description = description;
        }

        /**
         * @return The Button
         */
        public JButton getButton() {
            return button;
        }

        /**
         * @return The associated Action Object
         */
        public Action getAction() {
            return action;
        }

        /**
         * @return The Name
         */
        public String getName() {
            return name;
        }

        /**
         * @return the Description
         */
        public String getDescription() {
            return description;
        }

        /**
         * @return associated Icon
         */
        public ImageIcon getIcon() {
            return icon;
        }

        /**
         * Sets a new Action
         * @param ac action
         */
        public void setAction(Action ac) {
            action = ac;
            button = GuiUtils.getButton(action);
        }

    }// end enumeration
} // @jve:decl-index=0:visual-constraint="10,10"

