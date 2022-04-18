/*******************************************************************************
   PrefDialog class
*******************************************************************************/
import java.awt.*;

//******************************************************************************
// PrefDialog class defines a dialog window where user can adjust preferences.
//******************************************************************************

public class PrefDialog
       extends Dialog
{
   //---------------------------------------------------------------------------
   // Text fields for adjusting maximum test memory and maximum generator
   // memory.
   protected IntParmField  maxTestMemField;
   protected IntParmField  maxGenMemField;

   // Usual buttons.
   protected Button  okButton;
   protected Button  cancelButton;

   // Owner of this dialog.
   protected TestWizard  wiz;

   //***************************************************************************
   // The constructor constructs PrefDialog with owner <wizard>.
   //---------------------------------------------------------------------------
   public PrefDialog(TestWizard wizard) {
      // Create a modal dialog.
      super(wizard, TestWizard.APP_NAME+": Preferences", true);
      // Change layout manager to BoarderLayout.
      setLayout(new BorderLayout(0,20));

      // Remember the owner.
      wiz = wizard;

      // Add the components.
      add("North", createInputPanel());
      add("South", createButtonPanel());
   }

   //---------------------------------------------------------------------------
   // insets() customizes this dialog's insets.
   //---------------------------------------------------------------------------
   public Insets insets() {
      Insets  oldVal = super.insets();

      // Increase from default value by 5.
      return new Insets(oldVal.top + 10, oldVal.left + 10,
                        oldVal.bottom + 10, oldVal.right + 10);
   }

   //---------------------------------------------------------------------------
   // action() catches events on the buttons, then makes sure any adjustment on
   // the preferences are valid and close the dialog window.
   //---------------------------------------------------------------------------
   public boolean action(Event eve, Object arg) {
      if (eve.target == okButton) {
         // Check for errors before accepting the new preference settings.
         if (maxTestMemField.checkParm().length() == 0) {
            wiz.maxTestMemory = maxTestMemField.getParmValue();
         }
         if (maxGenMemField.checkParm().length() == 0) {
            wiz.maxGenMemory = maxGenMemField.getParmValue();
         }
         // Close and destroy the dialog window.
         hide();
         dispose();
         return true;
      }

      if (eve.target == cancelButton) {
         // Close and destroy the dialog window.
         hide();
         dispose();
         return true;
      }

      return super.action(eve, arg);
   }

   //***************************************************************************
   // createInputPanel() creates a panel that contains the text fields for
   // adjusting preferences
   //---------------------------------------------------------------------------
   protected Panel createInputPanel() {
      Panel  inputPanel;

      // Create preference text fields.
      maxTestMemField = new IntParmField("     Maximum Test Memory (KBytes)",
                                         1, Integer.MAX_VALUE);
      maxGenMemField  = new IntParmField("Maximum Generator Memory (KBytes)",
                                         1, Integer.MAX_VALUE);
      // Remove the "Random" button and display current preference values.
      maxTestMemField.no_default();
      maxTestMemField.setParm(wiz.maxTestMemory);
      maxGenMemField.no_default();
      maxGenMemField.setParm(wiz.maxGenMemory);

      // Create the panel that will contain the preference fields.
      inputPanel = new Panel();
      inputPanel.setLayout(new GridLayout(2,1));
      // Add in the preference fields.
      inputPanel.add(maxTestMemField);
      inputPanel.add(maxGenMemField);

      return inputPanel;
   }

   //---------------------------------------------------------------------------
   // createButtonPanel() creates a panel that contains the buttons.
   //---------------------------------------------------------------------------
   protected Panel createButtonPanel() {
      Panel  buttonPanel;

      // Create the buttons.
      okButton = new Button("OK");
      cancelButton = new Button("Cancel");

      // Create the panel that will contain the buttons.
      buttonPanel = new Panel();
      buttonPanel.setLayout(new GridLayout(1,5,5,0));
      // Add in dummy labels for spacing.
      buttonPanel.add(new Label());
      buttonPanel.add(new Label());
      // Add in the buttons.
      buttonPanel.add(okButton);
      buttonPanel.add(cancelButton);

      return buttonPanel;
   }
}
