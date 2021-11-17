/*
 * Created on 9/07/2003
 */
package jmarkov.gui;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PipedInputStream;
import java.io.PipedOutputStream;
import java.io.PrintStream;
import java.io.Writer;

import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextArea;
import javax.swing.JToolBar;
import javax.swing.filechooser.FileFilter;

/**
 * This class is an extended TextArea, that includes a toollbar for commands
 * such as saving, selecting and cleaning the text. Its main purpose is to
 * report information to the user and thus it does not allow editing.
 * 
 * @author Germán Riaño. Universidad de los Andes.
 */
public class TextPanel extends InfoPanel implements ActionListener {
    private static final long serialVersionUID = 1969;

    JTextArea txtArea = null;

    JButton btnSave = null;

    JButton btnClear = null;

    JButton btnCopy = null;

    // TODO poner boton buscar...
    static JFileChooser fc = null;

    String title = "";

    boolean autoShow = false;

    /**
     * Creates a TextPanel with this title and text.
     * 
     * @param title
     *            The title of the panel
     */
    public TextPanel(String title) {
        this(title, "", 0, 0);
    }

    @SuppressWarnings("unused")
    private TextPanel(String title, int rows, int cols) {
        this(title, "", rows, cols);
    }

    /**
     * Creates a TextPanel with this title and text, and the given number of
     * initial cols and rows.
     * 
     * @param title
     * @param text
     * @param rows
     * @param cols
     */
    private TextPanel(String title, String text, int rows, int cols) {
        super();
        this.title = title;
        // String NavLocation = "toolbarButtonGraphics/media/";
        String ImgLoc = "/jmarkov/gui/images/";
        ImageIcon icnSave = null;
        ImageIcon icnClear = null;
        // ImageIcon icnSaveAs = null;
        ImageIcon icnCopy = null;
        // ImageIcon icnReStart = null;
        try {
            // icnSave = new ImageIcon(ImgLoc + "Save16.gif");
            icnSave = GuiUtils.createIcon(ImgLoc + "Save16.gif");
            icnClear = GuiUtils.createIcon(ImgLoc + "Delete16.gif");
            // icnSaveAs = GuiUtils.createIcon(this, ImgLoc + "SaveAs16.gif");
            icnCopy = GuiUtils.createIcon(ImgLoc + "Copy16.gif");
            // icnReStart = Utils.createIcon(ImgLoc + "Rewind16.gif");
        } catch (NullPointerException e) {
        }

        JToolBar toolbar = new JToolBar(JToolBar.VERTICAL);
        btnSave = new JButton(icnSave);
        btnSave.setToolTipText("Saves the File");
        toolbar.add(btnSave);
        btnSave.addActionListener(this);

        btnClear = new JButton(icnClear);
        btnClear.setToolTipText("Clears the log");
        btnClear.addActionListener(this);
        toolbar.add(btnClear);

        btnCopy = new JButton(icnCopy);
        btnCopy.setToolTipText("Copy selected text to the clipboard");
        btnCopy.addActionListener(this);
        toolbar.add(btnCopy);

        toolbar.setFloatable(false);
        // txtArea= new JTextPane();
        txtArea = new JTextArea(text, rows, cols);
        txtArea.setEditable(false);
        // txtArea.setFont(new Font("Monospaced",Font.PLAIN,12));
        txtArea.setFont(new Font("Courier New", Font.PLAIN, 12));
        JScrollPane areaScrollPane = new JScrollPane(txtArea);
        // areaScrollPane.setVerticalScrollBarPolicy(
        // JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
        // areaScrollPane.setPreferredSize(new Dimension(750, 250));
        setLayout(new BorderLayout(5, 5));
        // setLayout(new BorderLayout());
        add(toolbar, BorderLayout.EAST);
        add(areaScrollPane, BorderLayout.CENTER);
        setBorder(BorderFactory.createCompoundBorder(BorderFactory
                .createTitledBorder(title), BorderFactory.createEmptyBorder(5,
                5, 5, 5)));
        setText(text);
    }

    /**
     * Sets the text (deleting whatever was there)
     * 
     * @param s
     *            The text to write
     */
    public void setText(String s) {
        setText(s, true);
    }

    /**
     * @param s
     *            Text to set
     * @param showText
     *            whether it should scroll to meke text visible.
     */
    public void setText(String s, boolean showText) {
        txtArea.setText(s);
        if (showText)
            showTxt();
    }

    /**
     * Appends the string and scrolls to make it visible.
     * 
     * @param s
     */
    public void append(String s) {
        append(s, autoShow);
    }

    /**
     * Appends text.
     * 
     * @param s
     *            The string to add.
     * @param showText
     *            Whether it should scroll to show the text.
     */
    public void append(String s, boolean showText) {
        txtArea.append(s);
        if (showText)
            showTxt();
    }

    /**
     * @return The current Text on the control
     */
    public String getText() {
        return txtArea.getText();
    }

    /**
     * Ensures that the last line of text writen is visible.
     */
    public void showTxt() {
        Container papi = this.getParent();
        Container nene = this;
        while (papi != null) {
            if (papi instanceof JTabbedPane) {
                ((JTabbedPane) papi).setSelectedComponent(nene);
            }
            nene = papi;
            papi = papi.getParent();
        }

        txtArea.setCaretPosition(txtArea.getDocument().getLength());
        txtArea.getCaret().setSelectionVisible(true);
        txtArea.revalidate();
        txtArea.repaint();
        Thread.yield();

        // TODO make gain focus, visible
    }

    /**
     * @return The OutputStream assciated with the control
     */
    public OutputStream getStream() {
        return new StreamCapturer();
    }

    /**
     * Forces this control to capture std. out.
     */
    public void captureStandardOut() {
        StreamCapturer sc = new StreamCapturer(System.out);
        System.setOut(new PrintStream(sc));
    }

    /**
     * Forces this control to capture std. err
     */
    public void captureStandardErr() {
        StreamCapturer sc = new StreamCapturer(System.err);
        System.setErr(new PrintStream(sc));
    }

    /**
     * Captures the given PrintStream
     * 
     * @param out
     */
    public void captureStream(PrintStream out) {
        StreamCapturer sc = new StreamCapturer(out);
        // TODO: esto debe estar mal..
        System.setOut(new PrintStream(sc));
    }

    void oldCaptureStandarIO() {
        PrintStream originalOut = System.out;
        try {
            PipedOutputStream wtr = new PipedOutputStream();
            PipedInputStream rdr = new PipedInputStream(wtr);
            System.setOut(new PrintStream(wtr));
            // txtArea.read(new BufferedReader(new
            // InputStreamReader(rdr)),null);
            new PipeThread(new BufferedReader(new InputStreamReader(rdr)))
                    .start();
        } catch (IOException ioe) {
            originalOut.println(ioe);
        }
        ;
    }

    /*
     * (non-Javadoc)
     * 
     * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
     */
    public void actionPerformed(ActionEvent e) {
        if (e.getSource() == btnClear) {
            clear();
        }
        if (e.getSource() == btnSave) {
            saveFile();
        }
        if (e.getSource() == btnCopy) {
            if (txtArea.getSelectionEnd() == txtArea.getSelectionStart())
                txtArea.selectAll();
            txtArea.copy();
        }
    }

    /**
     * Clears the current text.
     */
    public void clear() {
        txtArea.setText("");
    }

    /**
     * Calls the save dialog.
     */
    public void saveFile() {
        if (fc == null) {
            fc = new JFileChooser();
        }
        fc.setFileFilter(new FileFilter() {
            @Override
            public boolean accept(File f) {
                String path = f.getAbsolutePath();
                return (path.endsWith("txt") || f.isDirectory());
            }

            @Override
            public String getDescription() {
                return "Text Files (*.txt) ";
            }
        });
        String path = System.getProperty("user.dir")
                + System.getProperty("file.separator") + title + ".txt";
        fc.setSelectedFile(new File(path));
        int returnVal = fc.showSaveDialog(this);
        if (returnVal == JFileChooser.APPROVE_OPTION) {
            File file = fc.getSelectedFile();
            // debug("Opening: " + file.getName() + "." + "\n");
            try {
                Writer out = new BufferedWriter(new FileWriter(file));
                txtArea.write(out);
                out.close();
            } catch (IOException e) {
                System.out.println("Error writing fle " + file);
            }
        }
    }

    class StreamCapturer extends OutputStream {
        boolean keepold = false;

        PrintStream oldStream = null;

        /**
         * Constructor
         */
        public StreamCapturer() {
            super();
        }

        /**
         * @param ps
         */
        public StreamCapturer(PrintStream ps) {
            this(ps, false);
        }

        /**
         * @param ps
         *            The PrintStream to capture
         * @param keepold
         *            true if old should be kept.
         */
        public StreamCapturer(PrintStream ps, boolean keepold) {
            super();
            this.keepold = keepold;
            oldStream = ps;
        }

        /**
         * 
         * @param s 
         * @see java.io.PrintStream#print(java.lang.String)
         */

        public void print(String s) {
            if (keepold)
                oldStream.print(s);
            append(s);
        }

        /**
         * (non-Javadoc)
         * @param s 
         * 
         * @see java.io.PrintStream#print(java.lang.String)
         */

        public void println(String s) {
            if (keepold)
                oldStream.println(s);
            append(s + "\n");
        }

        /**
         * (non-Javadoc)
         * 
         * @see java.io.PrintStream#println()
         */

        public void println() {
            if (keepold)
                oldStream.println();
            append("\n");
        }

        /**
         * (non-Javadoc)
         * @param b 
         * 
         * @see java.io.PrintStream#print(boolean)
         */

        public void print(boolean b) {
            if (keepold)
                oldStream.print(b);
            append(new Boolean(b).toString());
        }

        /**
         * (non-Javadoc)
         * @param c 
         * 
         * @see java.io.PrintStream#print(char)
         */

        public void print(char c) {
            if (keepold)
                oldStream.print(c);
            append(new String(new char[] { c }));
        }

        /**
         * (non-Javadoc)
         * @param s 
         * 
         * @see java.io.PrintStream#print(char[])
         */

        public void print(char[] s) {
            if (keepold)
                oldStream.print(s);
            append(new String(s).toString());
        }

        /**
         * (non-Javadoc)
         * @param d 
         * 
         * @see java.io.PrintStream#print(double)
         */

        public void print(double d) {
            if (keepold)
                oldStream.print(d);
            append(new Double(d).toString());
        }

        /**
         * (non-Javadoc)
         * @param f 
         * 
         * @see java.io.PrintStream#print(float)
         */

        public void print(float f) {
            if (keepold)
                oldStream.print(f);
            append(new Float(f).toString());
        }

        /**
         * (non-Javadoc)
         * @param i 
         * 
         * @see java.io.PrintStream#print(int)
         */

        public void print(int i) {
            if (keepold)
                oldStream.print(i);
            append(new Integer(i).toString());
        }

        /**
         * (non-Javadoc)
         * @param l 
         * 
         * @see java.io.PrintStream#print(long)
         */

        public void print(long l) {
            if (keepold)
                oldStream.print(l);
            append(new Long(l).toString());
        }

        /**
         * (non-Javadoc)
         * @param obj 
         * 
         * @see java.io.PrintStream#print(java.lang.Object)
         */

        public void print(Object obj) {
            if (keepold)
                oldStream.print(obj);
            append(obj.toString());
        }

        /**
         * (non-Javadoc)
         * @param b 
         * 
         * @see java.io.PrintStream#print(boolean)
         */

        public void println(boolean b) {
            if (keepold)
                oldStream.println(b);
            append(new Boolean(b).toString() + "\n");
        }

        /**
         * (non-Javadoc)
         * @param c 
         * 
         * @see java.io.PrintStream#print(char)
         */

        public void println(char c) {
            if (keepold)
                oldStream.println(c);
            append(new String(new char[] { c }) + "\n");
        }

        /**
         * (non-Javadoc)
         * @param s 
         * 
         * @see java.io.PrintStream#println(char[])
         */

        public void println(char[] s) {
            if (keepold)
                oldStream.println(s);
            append(new String(s).toString() + "\n");
        }

        /**
         * (non-Javadoc)
         * @param d 
         * 
         * @see java.io.PrintStream#println(double)
         */

        public void println(double d) {
            if (keepold)
                oldStream.println(d);
            append(new Double(d).toString() + "\n");
        }

        /**
         * (non-Javadoc)
         * @param f 
         * 
         * @see java.io.PrintStream#println(float)
         */

        public void println(float f) {
            if (keepold)
                oldStream.println(f);
            append(new Float(f).toString() + "\n");
        }

        /**
         * (non-Javadoc)
         * @param i 
         * 
         * @see java.io.PrintStream#println(int)
         */

        public void println(int i) {
            if (keepold)
                oldStream.println(i);
            append(new Integer(i).toString() + "\n");
        }

        /**
         * (non-Javadoc)
         * @param l 
         * 
         * @see java.io.PrintStream#println(long)
         */

        public void println(long l) {
            if (keepold)
                oldStream.println(l);
            append(new Long(l).toString() + "\n");
        }

        /**
         * (non-Javadoc)
         * @param obj 
         * 
         * @see java.io.PrintStream#println(java.lang.Object)
         */

        public void println(Object obj) {
            if (keepold)
                oldStream.println(obj);
            append(obj.toString() + "\n");
        }

        /**
         * 
         * @see java.io.OutputStream#write(byte[], int, int)
         */

        @Override
        public void write(byte[] buf, int off, int len) {
            if (keepold)
                oldStream.write(buf, off, len);
            append(new String(buf, off, len));
        }

        /**
         * (non-Javadoc)
         * 
         * @see java.io.OutputStream#write(int)
         */

        @Override
        public void write(int b) {
            if (keepold)
                oldStream.write(b);
            append(new String(new byte[] { (byte) b }));
        }

    }

    class PipeThread extends Thread {

        private BufferedReader in = null;

        /**
         * @param in
         *            The reader.
         */
        public PipeThread(BufferedReader in) {
            super("PipeThread");
            this.in = in;
        }

        @Override
        public void run() {
            while (in != null) {
                try {
                    String input;
                    while ((input = in.readLine()) != null) {
                        append(input);
                        // out.flush();
                    }
                    yield(); // let other things run
                } catch (IOException e) {
                    System.err.println("PipedThread run: " + e);
                }
            }
        }
    }

    /**
     * @return true if showText() is called whenever append is called.
     */
    public boolean isAutoShow() {
        return autoShow;
    }

    /**
     * Sets whether append should call showText.
     * 
     * @param show
     *            True if showText() is called automatically.
     */
    public void setAutoShow(boolean show) {
        autoShow = show;
    }

    /*
     * (non-Javadoc)
     * 
     * @see jmarkov.gui.InfoPanel#updateMP()
     */
    @Override
    public void updateMP() {
        if (mp == null)
            setText("");
    }
}
