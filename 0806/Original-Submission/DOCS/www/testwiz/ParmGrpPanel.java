/*******************************************************************************
   ParmGrpPanel class
*******************************************************************************/
import java.awt.*;

//******************************************************************************
// ParmGrpPanel class defines the basic behavior and GUI structure of a group
// of parmeter fields embedded in a panel.
//******************************************************************************

public abstract class ParmGrpPanel
                extends Panel
{
   //---------------------------------------------------------------------------
   // Panel that contains the group of parameter fields.
   protected Panel  parmPanel;
   // Command line id that is associated with this group of parameters.
   protected String  cmdId;

   //***************************************************************************
   // The constructor constructs ParmGrpPanel with title <title>, that will
   // hold <numField> parameter fields.
   //---------------------------------------------------------------------------
   public ParmGrpPanel(String title, int numField) {
      // Create a panel that contains the title of the parameter group.
      Panel  namePanel = new Panel();
      namePanel.add(new Label(title));

      // Create the panel that will contain <numField> parameter fields.
      parmPanel = new Panel();
      parmPanel.setLayout(new GridLayout(numField, 1));

      // Change layout manager to BorderLayout.
      setLayout(new BorderLayout(0,10));
      
      // Add the components.
      add("North", namePanel);
      add("South", parmPanel);
   }

   //---------------------------------------------------------------------------
   // addField() adds the component <field>.
   //---------------------------------------------------------------------------
   public void addField(Component field) {
      parmPanel.add(field);
   }

   //---------------------------------------------------------------------------
   // getCmdId() returns the command line id associated with this group of
   // parameters.
   //---------------------------------------------------------------------------
   public String getCmdId() {
      return cmdId;
   }

   //---------------------------------------------------------------------------
   // checkParms() checks the validity of the parameter values, and returns a
   // description of the error.
   //
   // No error: checkParms() returns empty string
   // Warning : checkParms() returns warning message (non-empty string)
   //           getParms() returns command line segment (non-empty string)
   // Error   : checkParms() returns error message (non-empty string)
   //           getParms() returns empty string
   //---------------------------------------------------------------------------
   public abstract String checkParms();

   //---------------------------------------------------------------------------
   // getParms() translates the user-specified parameters into a command line
   // segment, and returns it in a string.
   // NOTE: Before calling getParms(), call checkParms() first.
   //---------------------------------------------------------------------------
   public abstract String getParms();
}
