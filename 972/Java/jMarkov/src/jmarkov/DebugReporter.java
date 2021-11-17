/*
 * Created on 9/07/2003
 */
package jmarkov;

import java.io.PrintWriter;

/**
 * A debug reporter is used to report debug Information from a program. It has
 * an internal integer number, called the Debuglevel, where 0 means no debug
 * information will be reporter and 5 a lot of verbose information will be
 * reported. The information is reported to standard output, to a PrintWriter
 * (which can be associated with a file), or to a TextPanel which can be
 * included in a graphic user interface.
 * 
 * @author German Riano. Universidad de los Andes.
 * 
 */
public class DebugReporter {
    // Current debug level
    private int curLevel = 1;
    private PrintWriter dbgWt = null;

    /**
     * Creates a debug reporter that will report to standard I/O.
     * 
     * @param initDebugLevel
     *            Initial debug level
     */
    public DebugReporter(int initDebugLevel) {
        super();
        curLevel = initDebugLevel;
    }

    /**
     * Creates a debug reporter that will send its output to the given
     * PrintWriter.
     * 
     * @param dbgWt
     *            the PrintWriter where the debug information will be sent.
     */
    public DebugReporter(PrintWriter dbgWt) {
        super();
        this.dbgWt = new PrintWriter(dbgWt, true);
    }

    /**
     * Sets the debug level, where level=0 means no debug info, level = 5
     * verbose info.
     * 
     * @param curLevel
     *            The curLevel to set to.
     */
    public final void setCurLevel(int curLevel) {
        this.curLevel = curLevel;
    }

    /**
     * @return Returns the curLevel.
     */
    public final int getCurLevel() {
        return curLevel;
    }

    /**
     * Reports this debug information. Newline and indent are true.
     * 
     * @see #debug(int, String, boolean, boolean)
     * @param level
     *            Use level=0 for very important things, level=5 less important.
     * @param s
     *            The String to report
     */
    public void debug(int level, String s) {
        debug(level, s, true, true);
    }

    /**
     * Reports this debug information. Info is indented if newline is selected.
     * 
     * @param level
     *            Use level=0 for very important things, level=5 less important.
     * @param s
     *            The String to report
     * @param newline
     *            whether newline is added at the end.
     * @see #debug(int, String, boolean, boolean)
     */
    public void debug(int level, String s, boolean newline) {
        debug(level, s, newline, newline);
    }

    /**
     * Reports this debug information.
     * 
     * @param level
     *            Use level=0 for very important things, level=5 less important.
     * @param s
     *            The String to report
     * @param newline
     *            whether newline is added at the end.
     * @param indent
     *            whwether information should go indented according to debug
     *            level.
     */
    public void debug(int level, String s, boolean newline, boolean indent) {
        if (level <= curLevel) {
            s = (indent) ? blank(2 * level) + s : s;
            if (dbgWt == null)
                dbgWt = new PrintWriter(System.out, true);
            if (newline)
                dbgWt.println(s);
            else
                dbgWt.print(s);
            Thread.yield();
        }
    }

    /**
     * @return current debug level, where level=0 means no debug info and level =
     *         5 verbose info.
     */
    public synchronized int getDebugLevel() {
        return curLevel;
    }

    /**
     * Sets the debug level, where level=0 means no debug info, level = 5
     * verbose info.
     * 
     * @param level
     *            debug level
     */
    public synchronized void setDebugLevel(int level) {
        curLevel = level;
    }

    /*
     * Returns a blank string with b characters.
     */
    private String blank(int b) {
        String stg = "";
        for (int i = 0; i < b; i++)
            stg += ' ';
        return stg;
    }

}
