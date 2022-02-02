/*******************************************************************************
   SerialPanel class
*******************************************************************************/


//******************************************************************************
// SerialPanel class defines a panel for entering and translating parameters of
// "Serial Test".
//******************************************************************************

public class SerialPanel
       extends ParmGrpPanel
{
   //---------------------------------------------------------------------------
   // Parameter field for number of divisions.
   protected IntParmField  numDivnField;
   // Parameter field for number of samples.
   protected IntParmField  numSampField;

   // The command line translation of the parameters.
   protected String  parmStr;
   // Owner of this panel.
   protected TestWizard  wiz;

   //***************************************************************************
   // The constructor constructs SerialPanel with owner <wizard>.
   //---------------------------------------------------------------------------
   public SerialPanel(TestWizard wizard) {
      // Set up a group panel that will contain two parameter fields.
      super("Serial Test", 2);
      // Set the command line id.
      super.cmdId = "serial";

      // Create the parameter fields.
      numDivnField = new IntParmField("Generate integers in [0,?]",
                                      2, Integer.MAX_VALUE,64);
      numSampField = new IntParmField("    Number of Pairs",
                                      1, Integer.MAX_VALUE,500000);

      // Initialize internal variables.
      parmStr = "";
      wiz = wizard;

      // Add the components.
      addField(numDivnField);
      addField(numSampField);
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
   public String checkParms() {
      String  errorMsg;
      long  numDivn;
      long  numSamp;
      long  limit;

      // Check parameters for errors and concatenate error messages.
      errorMsg = numDivnField.checkParm() + numSampField.checkParm();
      // If error exists, set command line translation to empty string and
      // return the error message.
      if (errorMsg.length() > 0) {
         parmStr = "";
         return errorMsg;
      }

      // Translate the parameters into command line segment and cache it.
      parmStr = numDivnField.getParmString() + " " +
                numSampField.getParmString();

      // Get the numerical values of the parameters.
      numDivn = numDivnField.getParmValue();
      numSamp = numSampField.getParmValue();

      // Create a warning message if the sampling requirement "minimum of 5
      // in the least likely bin" is not met.
      if (numSamp/(numDivn*numDivn) < 5) {
         errorMsg += "WARNING: sampling requirement not met\n" +
                     "         make sure " +
                     numSampField.getParmTitle() + " >= " +
                     numDivnField.getParmTitle() + "^2 * 5\n";
      }

      limit = wiz.maxTestMemory * wiz.KILO_BYTE / 4;
      // Create a warning message if memory used by the test will exceed
      // the memory constraint.
      if (numDivn*numDivn > limit) {
         errorMsg += "WARNING: memory constraint exceeded\n" +
                     "         make sure " +
                     numDivnField.getParmTitle() + "^2 <= " +
                     limit + "\n";
      }

      return errorMsg;
   }

   //---------------------------------------------------------------------------
   // getParms() translates the user-specified parameters into a command line
   // segment, and returns it in a string.
   // NOTE: Before calling getParms(), call checkParms() first.
   //---------------------------------------------------------------------------
   public String getParms() {
      // Return the cached translation.
      return parmStr;
   }
}
