/*******************************************************************************
   ParmField class
*******************************************************************************/
import java.awt.*;

//******************************************************************************
// ParmField class defines the basic behavior and GUI structure of a parameter
// field, which includes a label, a text field, and a button called "Random".
//******************************************************************************

public abstract class ParmField
                extends Panel
{
   //---------------------------------------------------------------------------
   // Label of the parameter.
   protected Label  parmLabel;
   // Text field where user can enter value for the parameter.
   protected TextField  entryField;
   // "Default" button.
   protected Button  defaultButton;

   // Ratio of (default text field length) / (parameter label length).
   protected final static double  RATIO = 1.0;

   //***************************************************************************
   // The constructor constructs ParmField with title <parmName>.
   //---------------------------------------------------------------------------
   public ParmField(String parmName) {
      // Change layout manager to BorderLayout.
      setLayout(new BorderLayout());

      // Create the label, the text field, and the button.
      parmLabel = new Label(parmName);
      // Default text field length is set to (parameter label length) * <RATIO>.
      entryField = new TextField(15);
      defaultButton = new Button("Default");

      // Add the components.
      add("West", parmLabel);
      add("Center", entryField);
      add("East", defaultButton);
   }

   //---------------------------------------------------------------------------
   // action() catches events on the button, then invokes default().
   //---------------------------------------------------------------------------
   public boolean action(Event event, Object arg) {
      // Check for action on the button.
      if (event.target == defaultButton) {
         setDefault();
         return true;
      }

      return super.action(event, arg);
   }

   //---------------------------------------------------------------------------
   // no_default() removes the button.
   //---------------------------------------------------------------------------
   public void no_default() {
      remove(defaultButton);
      defaultButton = null;
   }

   //---------------------------------------------------------------------------
   // getParmTitle() returns the name of the parameter.
   //---------------------------------------------------------------------------
   public String getParmTitle() {
      return parmLabel.getText().trim();
   }

   //---------------------------------------------------------------------------
   // setDefault() assigns a default value to the parameter.
   //---------------------------------------------------------------------------
   public abstract void setDefault();

   //---------------------------------------------------------------------------
   // checkParm() checks the validity of the parameter value, and returns a
   // description of the error.  Empty string indicates no error.
   //---------------------------------------------------------------------------
   public abstract String checkParm();

   //---------------------------------------------------------------------------
   // getParmString() returns the parameter as a string.
   // NOTE: Before calling getParmString(), call checkParm() first.
   //---------------------------------------------------------------------------
   public String getParmString() {
      return entryField.getText();
   }
}
