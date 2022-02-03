package jphase.GUI;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Rectangle;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.JDialog;

/**
 * Frame to specify the parameters of a new PH variable
 * @author Juan F. Pérez
 * @version 1.0
 */
public class ParamFrame extends JDialog{
    /**
     * Class Version
     */
    private static final long serialVersionUID = 1L;


    /**
     * Main Panel
     */
    private JPanel principalPanel;

    /**
     * Alert label 
     */
    JLabel alertLabel;

    /**
     * Parameter names
     */
    JLabel[] paramNames;

    /**
     * Parameter Values
     */
    JTextField[] paramValues;

    /**
     * Accept Button
     */
    JButton yesButton;

    /**
     * Cancel Button
     */
    JButton noButton;

    /**
     * True if all the parameters are completely specified
     */
    public boolean res = false;

    /**
     * Resulting Phase variable
     */
    //private PhaseVar var = null;

    /**
     * Frame width
     */
    private int width = 360;

    /**
     * Frame height
     */
    private int height = 220;

    /**
     * Minimum allowed frame width 
     */
    private int minWidth = 360;

    /**
     * Minimum allowed frame height
     */
    //private int minHeight = 220;



    /**
     * Main frame constructor 
     * @param varType Variable definition in String version
     */
    public ParamFrame(String varType){
        this.setTitle("JPhase - New Variable Parameters");
        this.setResizable(false);
        this.setModal(true);

        principalPanel = new JPanel();
        principalPanel.setLayout(null);
        principalPanel.setPreferredSize(new java.awt.Dimension(width, height));
        principalPanel.setOpaque(true);

        if(varType.equals("Expo")){
            alertLabel= new JLabel("Parameters of the "+varType+" Variable");
            alertLabel.setBounds(new Rectangle(width/2-110, 10, 260, 20));
            alertLabel.setFont(new java.awt.Font("@Arial Unicode MS", 1, 12));
            alertLabel.setOpaque(true);
            principalPanel.add(alertLabel, null);
            
            paramNames = new JLabel[1] ;
            paramNames[0] = new JLabel("lambda");

            paramValues = new JTextField[1] ;
            paramValues[0] = new JTextField();

            for(int i = 0; i < paramNames.length; i++){
                paramNames[i].setBounds(new Rectangle(width/2 - 100, 40*(i+1), 80, 20));
                principalPanel.add(paramNames[i],null );

                paramValues[i].setBackground(Color.WHITE);
                paramValues[i].setBounds(new Rectangle(width/2 + 20, 40*(i+1), 80, 20));
                principalPanel.add(paramValues[i],null );
            }


        }else if(varType.equals("Erlang")){


            alertLabel= new JLabel("Parameters of the "+varType+" Variable");
            alertLabel.setBounds(new Rectangle(width/2-110, 10, 260, 20));
            alertLabel.setFont(new java.awt.Font("@Arial Unicode MS", 1, 12));
            alertLabel.setOpaque(true);
            principalPanel.add(alertLabel, null);
            
            paramNames = new JLabel[2];
            paramNames[0] = new JLabel("lambda");
            paramNames[1] = new JLabel("n");

            paramValues = new JTextField[2];
            paramValues[0] = new JTextField();
            paramValues[1] = new JTextField();

            for(int i = 0; i < paramNames.length; i++){
                paramNames[i].setBounds(new Rectangle(width/2 - 100, 40*(i+1), 80, 20));
                principalPanel.add(paramNames[i],null );

                paramValues[i].setBackground(Color.WHITE);
                paramValues[i].setBounds(new Rectangle(width/2 + 20, 40*(i+1), 80, 20));
                principalPanel.add(paramValues[i],null );
            }


        }else if(varType.equals("HyperExponential")){
            InputFrame newHypExFrame = new InputFrame(
                    "New HyperExponential Variable",
            "Enter the number of phases");
            newHypExFrame.setVisible(true);
            newHypExFrame.setFocusable(true);

            if(newHypExFrame.getRes()){
                int n = Integer.parseInt(newHypExFrame.getValue());
                width = minWidth;
                height = 120 + 25*(n+1);
                principalPanel.setPreferredSize(new java.awt.Dimension(width, height));

                paramNames = new JLabel[4+n];
                paramNames[0] = new JLabel("Number Of Phases");
                paramNames[0].setBounds(new Rectangle(width/2 - 70, 40, 60, 20));
                paramNames[1] = new JLabel("Phase");
                paramNames[1].setBounds(new Rectangle(width/2 - 140, 70, 40, 20));
                paramNames[2] = new JLabel("Rate");
                paramNames[2].setBounds(new Rectangle(width/2 - 80, 70, 80, 20));
                paramNames[3] = new JLabel("Initial probability");
                paramNames[3].setBounds(new Rectangle(width/2 + 20, 70, 120, 20));

                principalPanel.add(paramNames[0],null );
                principalPanel.add(paramNames[1],null );
                principalPanel.add(paramNames[2],null );
                principalPanel.add(paramNames[3],null );

                for(int i = 0; i < n; i++){
                    paramNames[4+i] = new JLabel(""+i);
                    paramNames[4+i].setBounds(new Rectangle(
                            width/2 - 140, 70 + 25*(i+1), 40, 20));
                    principalPanel.add(paramNames[4+i],null );

                }

                paramValues = new JTextField[2*n + 1];
                paramValues[0] = new JTextField(""+n);
                paramValues[0].setEditable(false);
                paramValues[0].setBounds(new Rectangle(width/2 + 10, 40, 60, 20));
                principalPanel.add(paramValues[0],null );

                for(int i = 0; i < n; i++){
                    paramValues[i+1] = new JTextField("");
                    paramValues[i+1].setBounds(new Rectangle(
                            width/2 - 80, 70 + 25*(i+1), 80, 20));
                    principalPanel.add(paramValues[i+1],null );
                    
                    paramValues[n+i+1] = new JTextField("");
                    paramValues[n+i+1].setBounds(new Rectangle(
                            width/2 + 20, 70 + 25*(i+1), 80, 20));
                    principalPanel.add(paramValues[n+i+1],null );
                }
            }
        }else if(varType.equals("Coxian")){

            InputFrame newCoxFrame = new InputFrame(
                    "JPhase: New Coxian Variable",
            "Enter the number of phases in the Coxian");
            newCoxFrame.setVisible(true);
            newCoxFrame.setFocusable(true);

            if(newCoxFrame.getRes()){
                int n = Integer.parseInt(newCoxFrame.getValue());
                width = minWidth;
                height = 120 + 25*(n+1);
                principalPanel.setPreferredSize(new java.awt.Dimension(width, height));

                paramNames = new JLabel[4+n];
                paramNames[0] = new JLabel("n");
                paramNames[0].setBounds(new Rectangle(width/2 - 70, 40, 60, 20));
                paramNames[1] = new JLabel("Phase");
                paramNames[1].setBounds(new Rectangle(width/2 - 140, 70, 40, 20));
                paramNames[2] = new JLabel("Rate");
                paramNames[2].setBounds(new Rectangle(width/2 - 80, 70, 80, 20));
                paramNames[3] = new JLabel("Non-Absorption Prob");
                paramNames[3].setBounds(new Rectangle(width/2 + 20, 70, 120, 20));

                principalPanel.add(paramNames[0],null );
                principalPanel.add(paramNames[1],null );
                principalPanel.add(paramNames[2],null );
                principalPanel.add(paramNames[3],null );

                for(int i = 0; i < n; i++){
                    paramNames[4+i] = new JLabel(""+i);
                    paramNames[4+i].setBounds(new Rectangle(
                            width/2 - 140, 70 + 25*(i+1), 40, 20));
                    principalPanel.add(paramNames[4+i],null );

                }

                paramValues = new JTextField[2*n];
                paramValues[0] = new JTextField(""+n);
                paramValues[0].setEditable(false);
                paramValues[0].setBounds(new Rectangle(width/2 + 10, 40, 60, 20));
                principalPanel.add(paramValues[0],null );

                for(int i = 0; i < n; i++){
                    paramValues[i+1] = new JTextField("");
                    paramValues[i+1].setBounds(new Rectangle(
                            width/2 - 80, 70 + 25*(i+1), 80, 20));
                    principalPanel.add(paramValues[i+1],null );
                    if(i<n-1){
                        paramValues[n+i+1] = new JTextField("");
                        paramValues[n+i+1].setBounds(new Rectangle(
                                width/2 + 20, 70 + 25*(i+1), 80, 20));
                        principalPanel.add(paramValues[n+i+1],null );
                    }
                }


            }

        }else if(varType.equals("General Phase")){
        
        }else if(varType.equals("HyperErlang")){
            InputFrame newHypExFrame = new InputFrame(
                    "New HyperErlang Variable",
            "Enter the number of branches");
            newHypExFrame.setVisible(true);
            newHypExFrame.setFocusable(true);

            if(newHypExFrame.getRes()){
                int n = Integer.parseInt(newHypExFrame.getValue());
                width = minWidth+  80;
                height = 120 + 25*(n+1);
                principalPanel.setPreferredSize(new java.awt.Dimension(width, height));

                paramNames = new JLabel[5+n];
                paramNames[0] = new JLabel("Number Of Branches");
                paramNames[0].setBounds(new Rectangle(width/2 - 70, 40, 60, 20));
                paramNames[1] = new JLabel("Phase");
                paramNames[1].setBounds(new Rectangle(width/2 - 140, 70, 40, 20));
                paramNames[2] = new JLabel("Rate");
                paramNames[2].setBounds(new Rectangle(width/2 - 80, 70, 80, 20));
                paramNames[3] = new JLabel("Initial probability");
                paramNames[3].setBounds(new Rectangle(width/2 + 20, 70, 120, 20));
                paramNames[4] = new JLabel("Phases");
                paramNames[4].setBounds(new Rectangle(width/2 + 120, 70, 40, 20));

                principalPanel.add(paramNames[0],null );
                principalPanel.add(paramNames[1],null );
                principalPanel.add(paramNames[2],null );
                principalPanel.add(paramNames[3],null );
                principalPanel.add(paramNames[4],null );

                for(int i = 0; i < n; i++){
                    paramNames[5+i] = new JLabel(""+i);
                    paramNames[5+i].setBounds(new Rectangle(
                            width/2 - 140, 70 + 25*(i+1), 40, 20));
                    principalPanel.add(paramNames[5+i],null );

                }

                paramValues = new JTextField[3*n + 1];
                paramValues[0] = new JTextField(""+n);
                paramValues[0].setEditable(false);
                paramValues[0].setBounds(new Rectangle(width/2 + 10, 40, 60, 20));
                principalPanel.add(paramValues[0],null );

                for(int i = 0; i < n; i++){
                    paramValues[i+1] = new JTextField("");
                    paramValues[i+1].setBounds(new Rectangle(
                            width/2 - 80, 70 + 25*(i+1), 80, 20));
                    principalPanel.add(paramValues[i+1],null );
                    
                    paramValues[n+i+1] = new JTextField("");
                    paramValues[n+i+1].setBounds(new Rectangle(
                            width/2 + 20, 70 + 25*(i+1), 80, 20));
                    principalPanel.add(paramValues[n+i+1],null );
                    
                    paramValues[2*n+i+1] = new JTextField("");
                    paramValues[2*n+i+1].setBounds(new Rectangle(
                            width/2 + 120, 70 + 25*(i+1), 80, 20));
                    principalPanel.add(paramValues[2*n+i+1],null );
                }
            }
        }else{
            System.out.println("Non-known distribution");       
        }

        yesButton = new JButton("Enter");
        yesButton.setBounds(new Rectangle(width/2 - 100, height - 40, 80, 25));
        yesButton.addActionListener(
                new ParamFrame_yesButton_actionAdapter(this));
        yesButton.setOpaque(true);
        principalPanel.add(yesButton, null);

        noButton = new JButton("Cancel");
        noButton.setBounds(new Rectangle(width/2 + 20, height - 40, 80, 25));
        noButton.addActionListener(
                new ParamFrame_noButton_actionAdapter(this));
        principalPanel.add(noButton, null);

        this.getContentPane().add(principalPanel, BorderLayout.CENTER);
        this.centrarFrame();
        pack();

        // Screen size 
        Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
        Dimension frameSize = this.getSize();
        if (frameSize.height > screenSize.height) frameSize.height = screenSize.height;
        if (frameSize.width > screenSize.width) frameSize.width = screenSize.width;
        this.setLocation((screenSize.width - frameSize.width) / 2, (screenSize.height - frameSize.height) / 2);
    }

    /**
     * 
     * @return True if all parameter are fully specified
     */
    public boolean getRes(){
        return this.res;
    }

    /**
     * @param e Event
     */
    void yesButton_actionPerformed(ActionEvent e) {
        if(this.paramValues!=null){
            for(int i = 0; i < this.paramValues.length; i++)
                if(this.paramValues[i].getText().equals("")){
                    JOptionPane.showMessageDialog(null,  
                            "You must enter ALL the Parameters", 
                            "JPhase Alert", 
                            JOptionPane.INFORMATION_MESSAGE);
                }else{
                    this.res = true;
                    this.setVisible(false);
                }
        }
    }


    /**
     * @param e Event
     */
    void noButton_actionPerformed(ActionEvent e) {
        this.setVisible(false);
    }




    /**
     * Centers the main frame 
     */ 
    private void centrarFrame() {
        Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
        Dimension frameSize = this.getSize();
        if (frameSize.height > screenSize.height) 
            frameSize.height = screenSize.height;
        if (frameSize.width > screenSize.width) 
            frameSize.width = screenSize.width;
        this.setLocation((screenSize.width - width) / 2 , 
                (screenSize.height - height)/2);
    }

    /**
     * Yes button action listener 
     * @author Juan F. Perez
     *
     */
    class ParamFrame_yesButton_actionAdapter
    implements java.awt.event.ActionListener {

        /**
         * 
         */
        ParamFrame adaptee;

        /**
         * 
         * @param adaptee
         */
        ParamFrame_yesButton_actionAdapter(ParamFrame adaptee) {
            this.adaptee = adaptee;
        }

        public void actionPerformed(ActionEvent e) {
            adaptee.yesButton_actionPerformed(e);
        }
    }

    /**
     * No button action listener 
     * @author Juan F. Perez
     *
     */
    class ParamFrame_noButton_actionAdapter
    implements java.awt.event.ActionListener {

        /**
         * 
         */
        ParamFrame adaptee;

        /**
         * 
         * @param adaptee
         */
        ParamFrame_noButton_actionAdapter(ParamFrame adaptee) {
            this.adaptee = adaptee;
        }

        public void actionPerformed(ActionEvent e) {
            adaptee.noButton_actionPerformed(e);
        }
    }
}

