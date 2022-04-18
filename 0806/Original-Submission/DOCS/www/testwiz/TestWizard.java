/*******************************************************************************
   TestWizard class
*******************************************************************************/
import java.awt.*;

//******************************************************************************
// TestWizard class provides the main interface where user specifies the test
// parameters and obtains the command line translations.
//******************************************************************************

public class TestWizard
       extends Frame
{
   //---------------------------------------------------------------------------
   // Menu items under "File" menu.
   protected MenuItem  exitMI;
   // Menu items under "Test" menu.
   protected MenuItem  eqdistrMI;
   protected MenuItem  serialMI;
   protected MenuItem  gapMI;
   protected MenuItem  permMI;
   protected MenuItem  runsMI;
   protected MenuItem  collisionsMI;
   protected MenuItem  couponMI;
   protected MenuItem  maxtMI;
   protected MenuItem  pokerMI;
   protected MenuItem  sumMI;
   // Menu items under "Generator" menu.
   protected MenuItem  lcgMI;
   protected MenuItem  lfgMI;
   protected MenuItem  lcg64MI;
   protected MenuItem  cmrgMI;
   protected MenuItem  pmlcgMI;
   // Menu items under "Options" menu.
   protected MenuItem  prefMI;

   // Parameter panels for various tests and generators.
   protected ParmGrpPanel  testParmPanel;
   protected ParmGrpPanel  genParmPanel;
   // Button for invoking command line translation.
   protected Button  translButton;

   // Output screen.
   protected OutputFrame  outputScrn;
   // Preference screen.
   protected PrefDialog   prefScrn;

   // Memory limits on tests and generators.
   protected long  maxTestMemory;
   protected long  maxGenMemory;

   // Official name for this application.
   public final static String  APP_NAME = "Test Wizard";
   // Conversion factor for kilobytes.
   public final static int  KILO_BYTE = 1024;

   //***************************************************************************
   // The constructor.
   //---------------------------------------------------------------------------
   public TestWizard() {
      super(APP_NAME);

      // Change the layout manager to GridBagLayout.
      setLayout(new GridBagLayout());
      // Add menu bar.
      setMenuBar(createMenuBar());

      // Initialize components.
      testParmPanel = null;
      genParmPanel  = null;

      translButton  = new Button ("Translate");

      outputScrn = null;
      prefScrn = null;

      // Default memory limits to 50Mbytes.
      maxTestMemory = 50000;
      maxGenMemory = 5000;
   }

   //---------------------------------------------------------------------------
   // finalize() prints out a message. Used for debugging purpose only.
   //---------------------------------------------------------------------------
   protected void finalize() throws Throwable {
      System.out.println("TestWizard::finalize()");
      super.finalize();
   }

   //---------------------------------------------------------------------------
   // handleEvent() catches "close window" event, then invokes collapse().
   //---------------------------------------------------------------------------
   public boolean handleEvent(Event eve) {
      if (eve.id == Event.WINDOW_DESTROY) {
         collapse();
         return true;
      }

      return super.handleEvent(eve);
   }

   //---------------------------------------------------------------------------
   // action() catches events on the button and the menu bar.
   //---------------------------------------------------------------------------
   public boolean action(Event eve, Object arg) {
      // Translate the parameters if the button is clicked.
      if (eve.target == translButton) {
         translate();
         return true;
      }

      // Display the "Equidistribution Test" panel for that menu choice.
      if (eve.target == eqdistrMI) {
         // Create and add a new "Equidistribution Test" panel only if none is
         // displayed.
         if (!(testParmPanel instanceof EquidistPanel)) {
            testParmPanel = new EquidistPanel(this);
            adjustPanel();
         }
         return true;
      }

      // Display the "Serial Test" panel for that menu choice.
      if (eve.target == serialMI) {
         // Create and add a new "Serial Test" panel only if none is displayed.
         if (!(testParmPanel instanceof SerialPanel)) {
            testParmPanel = new SerialPanel(this);
            adjustPanel();
         }
         return true;
      }

      // Display the "Gap Test" panel for that menu choice.
      if (eve.target == gapMI) {
         // Create and add a new "Gap Test" panel only if none is displayed.
         if (!(testParmPanel instanceof GapPanel)) {
            testParmPanel = new GapPanel(this);
            adjustPanel();
         }
         return true;
      }

      // Display the "Permutation Test" panel for that menu choice.
      if (eve.target == permMI) {
         // Create and add a new "Permutation Test" panel only if none is
         // displayed.
         if (!(testParmPanel instanceof PermPanel)) {
            testParmPanel = new PermPanel(this);
            adjustPanel();
         }
         return true;
      }

      // Display the "Runs Test" panel for that menu choice.
      if (eve.target == runsMI) {
         // Create and add a new "Runs Test" panel only if none is displayed.
         if (!(testParmPanel instanceof RunsPanel)) {
            testParmPanel = new RunsPanel(this);
            adjustPanel();
         }
         return true;
      }


      // Display the "Collisions Test" panel for that menu choice.
      if (eve.target == collisionsMI) {
         // Create and add a new "Collisions Test" panel only if none is displayed.
         if (!(testParmPanel instanceof CollisionsPanel)) {
            testParmPanel = new CollisionsPanel(this);
            adjustPanel();
         }
         return true;
      }

      // Display the "Coupon Test" panel for that menu choice.
      if (eve.target == couponMI) {
         // Create and add a new "Coupon Test" panel only if none is displayed.
         if (!(testParmPanel instanceof CouponPanel)) {
            testParmPanel = new CouponPanel(this);
            adjustPanel();
         }
         return true;
      }

      // Display the "Maxt Test" panel for that menu choice.
      if (eve.target == maxtMI) {
         // Create and add a new "Maxt Test" panel only if none is displayed.
         if (!(testParmPanel instanceof MaxtPanel)) {
            testParmPanel = new MaxtPanel(this);
            adjustPanel();
         }
         return true;
      }

      // Display the "Poker Test" panel for that menu choice.
      if (eve.target == pokerMI) {
         // Create and add a new "Poker Test" panel only if none is displayed.
         if (!(testParmPanel instanceof PokerPanel)) {
            testParmPanel = new PokerPanel(this);
            adjustPanel();
         }
         return true;
      }

      // Display the "Sum Test" panel for that menu choice.
      if (eve.target == sumMI) {
         // Create and add a new "Sum Test" panel only if none is displayed.
         if (!(testParmPanel instanceof SumPanel)) {
            testParmPanel = new SumPanel(this);
            adjustPanel();
         }
         return true;
      }

      // Display the "LCG" panel for that menu choice.
      if (eve.target == lcgMI) {
         // Create and add a new "LCG" panel only if none is displayed.
         if (!(genParmPanel instanceof LcgPanel)) {
            genParmPanel = new LcgPanel(this);
            adjustPanel();
         }
         return true;
      }

      // Display the "LFG" panel for that menu choice.
      if (eve.target == lfgMI) {
         // Create and add a new "LFG" panel only if none is displayed.
         if (!(genParmPanel instanceof LfgPanel)) {
            genParmPanel = new LfgPanel(this);
            adjustPanel();
         }
         return true;
      }


      // Display the "LCG64" panel for that menu choice.
      if (eve.target == lcg64MI) {
         // Create and add a new "LCG64" panel only if none is displayed.
         if (!(genParmPanel instanceof lcg64Panel)) {
            genParmPanel = new lcg64Panel(this);
            adjustPanel();
         }
         return true;
      }


      // Display the "CMRG" panel for that menu choice.
      if (eve.target == cmrgMI) {
         // Create and add a new "CMRG" panel only if none is displayed.
         if (!(genParmPanel instanceof cmrgPanel)) {
            genParmPanel = new cmrgPanel(this);
            adjustPanel();
         }
         return true;
      }


      // Display the "PMLCG" panel for that menu choice.
      if (eve.target == pmlcgMI) {
         // Create and add a new "PMLCG" panel only if none is displayed.
         if (!(genParmPanel instanceof pmlcgPanel)) {
            genParmPanel = new pmlcgPanel(this);
            adjustPanel();
         }
         return true;
      }

      // Pop up the preference screen for that menu choice.
      if (eve.target == prefMI) {
         // Create the preference screen.
         prefScrn = new PrefDialog(this);

         // Dummy statement to reduce flickering when the dialog pops
         // up in an applet.
         prefScrn.resize(300,70);
         // Pack and display the preference screen.
         prefScrn.pack();
         prefScrn.show();

         return true;
      }

      // Destroy TestWizard for that menu choice.
      if (eve.target == exitMI) {
         collapse();
         return true;
      }

      return super.action(eve, arg);
   }

   //***************************************************************************
   // createMenuBar() creates the menu bar.
   //---------------------------------------------------------------------------
   protected MenuBar createMenuBar() {
      MenuBar  mBar = new MenuBar();
      Menu     fileMenu = new Menu("File");
      Menu     testMenu = new Menu("Test");
      Menu     genMenu  = new Menu("Generator");
      Menu     optionMenu = new Menu("Options");

      // Create and add items under menu "Exit".
      exitMI = new MenuItem("Exit");
      fileMenu.add(exitMI);

      // Create and add items under menu "Test".
      eqdistrMI = new MenuItem("Equidistribution");
      serialMI = new MenuItem("Serial");
      gapMI = new MenuItem("Gap");
      permMI = new MenuItem("Permutation");
      runsMI = new MenuItem("Runs");
      collisionsMI = new MenuItem("Collisions");
      couponMI = new MenuItem("Coupon");
      maxtMI = new MenuItem("Maxt");
      pokerMI = new MenuItem("Poker");
      sumMI = new MenuItem("Sum");
      testMenu.add(eqdistrMI);
      testMenu.add(serialMI);
      testMenu.add(gapMI);
      testMenu.add(permMI);
      testMenu.add(runsMI);
      testMenu.add(collisionsMI);
      testMenu.add(couponMI);
      testMenu.add(maxtMI);
      testMenu.add(pokerMI);
      testMenu.add(sumMI);

      // Create and add items under menu "Generator".
      lcgMI = new MenuItem("48 bit Linear Congruential");
      lfgMI = new MenuItem("Lagged Fibonacci");
      lcg64MI = new MenuItem("64 Linear Congruential");
      cmrgMI = new MenuItem("Combined Mulitple Recursive");
      pmlcgMI = new MenuItem("Prime Modulus Linear Congruential");
      genMenu.add(lcgMI);
      genMenu.add(lfgMI);
      genMenu.add(lcg64MI);
      genMenu.add(cmrgMI);
      genMenu.add(pmlcgMI);

      // Create and add items under menu "Options".
      prefMI = new MenuItem("Preferences...");
      optionMenu.add(prefMI);

      // Add menus onto menu bar.
      mBar.add(fileMenu);
      mBar.add(testMenu);
      mBar.add(genMenu);
      mBar.add(optionMenu);

      return mBar;
   }

   //---------------------------------------------------------------------------
   // adjustPanel() adjusts the display to accomodate the new set of parameter
   // panels.
   //---------------------------------------------------------------------------
   protected void adjustPanel() {
      GridBagLayout  layoutMgr = (GridBagLayout) getLayout();
      GridBagConstraints  constraint = new GridBagConstraints();

      // All other constraints use their default values.
      constraint.gridx = 0;
      constraint.fill  = GridBagConstraints.HORIZONTAL;
      constraint.insets  = new Insets(20,10,10,10);
      constraint.weightx = 1;

      // Remove all components, so we have a clean slate.
      removeAll();

      // Add in the test parameter panel if there is one.
      if (testParmPanel != null) {
         layoutMgr.setConstraints(testParmPanel, constraint);
         add(testParmPanel);
      }

      // Add in the generator parameter panel if there is one.
      if (genParmPanel != null) {
         layoutMgr.setConstraints(genParmPanel, constraint);
         add(genParmPanel);
      }

      // Modify constraints for the button.
      constraint.anchor = GridBagConstraints.EAST;
      constraint.fill = GridBagConstraints.NONE;
      constraint.weighty = 1;
      // Add in the button.
      layoutMgr.setConstraints(translButton, constraint);
      add(translButton);

      // Enable the button if both test and generator panel are shown.
      if ((testParmPanel!=null) && (genParmPanel!=null)) translButton.enable();
      else translButton.disable();

      // Adjust the frame size.
      pack();
   }

   //---------------------------------------------------------------------------
   // translate() translates the user-specified parameters into the appropriate
   // command line.
   //---------------------------------------------------------------------------
   protected void translate() {
      String  message;
      String  command;
      String  testParms;
      String  genParms;

      // Check parameters for errors and concatenate error messages.
      message = testParmPanel.checkParms()+genParmPanel.checkParms();
      // Obtain the translated test parameters.
      testParms = testParmPanel.getParms();
      // Obtain the translated generator parameters.
      genParms = genParmPanel.getParms();

      // Formulate the command line.
      if ((testParms.length()==0) || (genParms.length()==0)) command = "";
      else command = testParmPanel.getCmdId() + "." +
                     genParmPanel.getCmdId() + " " +
                     genParms + " " + testParms;

      // Make sure the output screen exists.
      if (outputScrn == null) outputScrn = new OutputFrame(this);
      // Update and display the output screen.
      outputScrn.setText(message, command);
      outputScrn.show();
   }

   //---------------------------------------------------------------------------
   // collapse() destroys TestWizard and its children, and releases the GUI
   // resources. 
   //---------------------------------------------------------------------------
   protected void collapse() {
      // Destroys output screen and releases its GUI resources.
      if (outputScrn != null) {
         outputScrn.hide();
         outputScrn.dispose();
         outputScrn = null;
      }

      // Destroys preference screen and releases its GUI resources.
      if (prefScrn != null) {
         prefScrn.hide();
         prefScrn.dispose();
         prefScrn = null;
      }

      // Self-destruct.
      hide();
      dispose();
   }

   //---------------------------------------------------------------------------
   // Entry point for running TestWizard standalone.
   //---------------------------------------------------------------------------
   public static void main(String arg[]) {
      TestWizard  twiz = new TestWizard();

      twiz.setFont(new Font("Courier", Font.PLAIN, 16));
      twiz.resize(500,200);
      twiz.show();
   }
}
