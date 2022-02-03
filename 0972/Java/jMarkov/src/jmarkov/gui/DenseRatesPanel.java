/*
 * Created on Jul 27, 2003
 *
 */
package jmarkov.gui;

import java.awt.Component;
import java.awt.Dimension;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.ScrollPaneConstants;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;

import jmarkov.basic.State;

/**
 * This class presents the rate matrix through a table.
 * @author Germán Riaño. Universidad de los Andes
 */
public class DenseRatesPanel extends InfoPanel {

    private static final long serialVersionUID = 3278485673228567288L;

    JTable table = null;
    JTable rowHeaders = null; // @jve:decl-index=0:visual-constraint="37,71"
    String[] columnNames = new String[] { "No data yet ..." };
    JScrollPane scrollPane = null;

    /**
     * Default constructor.
     */
    public DenseRatesPanel() {
        super();
        MatrixTableModel tableModel = new MatrixTableModel();
        table = new JTable(tableModel);
        // table.setPreferredScrollableViewportSize(new Dimension(70,
        // 70));
        // Create the scroll pane and add the table to it.
        scrollPane = new JScrollPane(table,
                ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS,
                ScrollPaneConstants.HORIZONTAL_SCROLLBAR_ALWAYS);
        table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        // Add the scroll pane to this window.
        // setLayout(new BorderLayout(5,5));
        // add(scrollPane, BorderLayout.CENTER);
        setLayout(new BoxLayout(this, BoxLayout.LINE_AXIS));
        add(scrollPane);
    }

    @SuppressWarnings("unchecked")
    @Override
    public void updateMP() {
        if (mp == null)
            return;
        // Loads the data
        // this.matrix = mp.getRates();
        State[] states = mp.getStates().toStateArray();
        int N = states.length;
        // This is needed to create the first column.
        // By default JTable has column headers, but no row's.
        State[][] firstRowData = new State[N][1];
        for (int i = 0; i < N; i++) {
            firstRowData[i][0] = states[i];
        }

        StatesTableHeader stateRenderer = new StatesTableHeader(states);
        // Configure row headers
        rowHeaders = new JTable(firstRowData, new String[] { "" });
        TableColumn col0 = rowHeaders.getColumnModel().getColumn(0);
        col0.setHeaderRenderer(stateRenderer);
        col0.setCellRenderer(stateRenderer);
        int width = 100;
        int height = table.getPreferredSize().height;
        rowHeaders.setPreferredScrollableViewportSize(new Dimension(width,
                height));
        scrollPane.setRowHeaderView(rowHeaders);

        columnNames = new String[N];
        for (int i = 0; i < N; i++) {
            columnNames[i] = states[i].label();
        }
        MatrixTableModel tblModel = new MatrixTableModel();
        table.setModel(tblModel);
        for (int i = 0; i < N; i++) {
            TableColumn col = table.getColumnModel().getColumn(i);
            col.setMinWidth(50);
            col.setPreferredWidth(100);
            col.setHeaderRenderer(stateRenderer);
        }
    }

    /**
     * This model represents the data in a Matrix.
     * @author Germán Riaño. Universidad de los Andes.
     */
    class MatrixTableModel extends AbstractTableModel {

        /**        */
        private static final long serialVersionUID = -2060267097161157594L;

        /**
         * Constructor
         */
        public MatrixTableModel() {
        }

		public int getColumnCount() {
            return columnNames.length;
        }

		@SuppressWarnings("unchecked")
        public int getRowCount() {
            return (mp != null) ? (int) mp.getNumStates() : 1;
        }

        @Override
        public String getColumnName(int col) {
            return columnNames[col];
        }

		@SuppressWarnings("unchecked")
        public Object getValueAt(int row, int col) {
            /*
             * Object ob = new Object(); if (col==0) ob =
             * columnNames[row]; else ob = new
             * Double(matrix[row][col-1]); return ob;
             */
            double val = 0.0;
            if (mp != null) {
                State sts[] = mp.getStates().toStateArray();
                val = mp.getFinalRate(sts[row], sts[col]);
            }
            return (val > 0) ? new Double(val) : "";
        }

        private Class getColumnClassToDelete(int c) {
            return getValueAt(0, c).getClass();
        }

        public boolean isCellEditable(int row, int col) {
            return false;
        }
    }

    /**
     * Used in columns and row headers.
     * @author Germán Riaño. Universidad de los Andes
     */
    class StatesTableHeader extends Object implements TableCellRenderer {
        State[] states;

        /**
         * This class is used to render the states's GUI
         * @param states The model States
         */
        public StatesTableHeader(State[] states) {
            super();
            this.states = states;
        }

        /*
         * @see javax.swing.table.TableCellRenderer#getTableCellRendererComponent(javax.swing.JTable,
         *      java.lang.Object, boolean, boolean, int, int)
         */
		public Component getTableCellRendererComponent(JTable table,
                Object value, boolean isSelected, boolean hasFocus, int row,
                int column) {
            // JComponent comp = new JLabel(value.toString(),
            // JLabel.CENTER);
            JComponent comp = new JButton(value.toString());
            if (states != null) {
                if (value instanceof State) // row header
                    comp.setToolTipText(((State) value).description());
                if (value instanceof String)// col header
                    comp.setToolTipText(states[column].description());
            }
            return comp;
        }
    }

}