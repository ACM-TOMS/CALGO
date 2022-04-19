/*******************************************************************************
   OutputFrame class
*******************************************************************************/
import java.awt.*;

//******************************************************************************
// OutputFrame class defines a frame that includes a message text area and a 
// command line text field for output.
//******************************************************************************

public class OutputFrame
       extends Frame
{
   //---------------------------------------------------------------------------
   // Message area.
   protected TextArea  msgArea;
   // Command line text field.
   protected TextField  cmdLnField;
   // Button for closing the window.
   protected Button  closeButton;

   // Owner of this frame.
   protected TestWizard  wiz;

   //***************************************************************************
   // The constructor constructs OutputFrame with owner <wizard>.
   //---------------------------------------------------------------------------
   public OutputFrame(TestWizard wizard) {
      super(TestWizard.APP_NAME+": Output");
      // Synchronize with the owner's font.
      setFont(wizard.getFont());

      // Add the components.
      add("Center", createMsgPanel());
      add("South", createCmdPanel());

      // The size of <msgArea> limits the minimal pack size.
      pack();

      // Remember the owner.
      wiz = wizard;
   }

   //---------------------------------------------------------------------------
   // insets() customizes this frame's insets.
   //---------------------------------------------------------------------------
   public Insets insets() {
      Insets  oldVal = super.insets();

      // Increase from default value by 5.
      return new Insets(oldVal.top + 5, oldVal.left + 5,
                        oldVal.bottom + 5, oldVal.right + 5);
   }

   //---------------------------------------------------------------------------
   // handleEvent() catches "close window" event, then hides this frame.
   //---------------------------------------------------------------------------
   public boolean handleEvent(Event eve) {
      if (eve.id == Event.WINDOW_DESTROY) {
         hide();
         return true;
      }

      return super.handleEvent(eve);
   }

   //---------------------------------------------------------------------------
   // action() catches events on the button, then hides this frame.
   //---------------------------------------------------------------------------
   public boolean action(Event event, Object arg) {
      // Check if there is action on the button.
      if (event.target == closeButton) {
         hide();
         return true;
      }

      return super.action(event, arg);
   }

   //---------------------------------------------------------------------------
   // setText() copies <message> into the message text area and <command> into
   // the command line text field.
   //---------------------------------------------------------------------------
   public void setText(String message, String command) {
      msgArea.setText(message);
      cmdLnField.setText(command);
   }

   //***************************************************************************
   // createMsgPanel() creates a panel that contains the message area.
   //---------------------------------------------------------------------------
   protected Panel createMsgPanel() {
      Panel  msgPanel;
      Panel  namePanel;

      // Create a panel that contains the title of the message area.
      namePanel = new Panel();
      namePanel.setLayout(new BorderLayout());
      namePanel.add("West", new Label("Messages"));

      // Create the message text area.
      msgArea = new TextArea(15,50);

      // Create a panel that will contain the message area and the title.
      msgPanel = new Panel();
      msgPanel.setLayout(new BorderLayout());
      // Add in the components.
      msgPanel.add("North", namePanel);
      msgPanel.add("Center", msgArea);

      return msgPanel;
   }

   //---------------------------------------------------------------------------
   // createCmdPanel() creates a panel that contains the command line text field
   // and the "Close" button.
   //---------------------------------------------------------------------------
   protected Panel createCmdPanel() {
      Panel  cmdPanel;
      Panel  namePanel;
      Panel  buttonPanel;

      // Create a panel that contains the title of the command line.
      namePanel = new Panel();
      namePanel.setLayout(new BorderLayout());
      namePanel.add("West", new Label("Command Line"));

      // Create the command line text field.
      cmdLnField = new TextField();

      // Create a panel that contains the "Close" button.
      buttonPanel = new Panel();
      buttonPanel.setLayout(new GridLayout(1,3));
      // Add dummy labels to space the button right.
      buttonPanel.add(new Label());
      buttonPanel.add(new Label());
      // Create and add in the button.
      closeButton = new Button("Close");
      buttonPanel.add(closeButton);

      // Create a panel that will contain the command line text field and the
      // button.
      cmdPanel = new Panel();
      cmdPanel.setLayout(new BorderLayout());
      // Add in the components.
      cmdPanel.add("North", namePanel);
      cmdPanel.add("Center", cmdLnField);
      cmdPanel.add("South", buttonPanel);

      return cmdPanel;
   }
}
