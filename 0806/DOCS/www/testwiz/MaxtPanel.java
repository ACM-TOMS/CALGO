/*******************************************************************************
   MaxtPanel class
*******************************************************************************/


//******************************************************************************
// MaxtPanel class defines a panel for entering and translating parameters of
// "Maxt Test".
//******************************************************************************

public class MaxtPanel
       extends ParmGrpPanel
{
   //---------------------------------------------------------------------------
   // Parameter field for block length.
   protected IntParmField  blockLengthField;
   // Parameter field for number of blocks.
   protected IntParmField  numSampField;

   // The command line translation of the parameters.
   protected String  parmStr;
   // Owner of this panel.
   protected TestWizard  wiz;

   //***************************************************************************
   // The constructor constructs PermPanel with owner <wizard>.
   //---------------------------------------------------------------------------
   public MaxtPanel(TestWizard wizard) {
      // Set up a group panel that will contain two parameter fields.
      super("Maxt Test", 2);
      // Set the command line id.
      super.cmdId = "maxt";

      // Create the parameter fields.
      numSampField   = new IntParmField("Numbers of blocks",
                                        2, Integer.MAX_VALUE,100000);
      blockLengthField = new IntParmField("Length of each block",
                                        2, Integer.MAX_VALUE,10);
 
      // Initialize internal variables.
      parmStr = "";
      wiz = wizard;

      // Add the components.
      addField(numSampField);
      addField(blockLengthField);
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
      long  numSamp;
      long  blockLength;

      String  errorMsg;
      long  limit;

      // Check parameters for errors and concatenate error messages.
      errorMsg = numSampField.checkParm() + blockLengthField.checkParm();
      // If error exists, set command line translation to empty string and
      // return the error message.
      if (errorMsg.length() > 0) {
         parmStr = "";
         return errorMsg;
      }

      // Translate the parameters into command line segment and cache it.
      parmStr = numSampField.getParmString() + " " +
                blockLengthField.getParmString();

      // Get the numerical values of the parameters.
      numSamp = numSampField.getParmValue();
      blockLength = blockLengthField.getParmValue();

      limit = wiz.maxTestMemory * wiz.KILO_BYTE / 4;
      // Create a warning message if memory used by the test will exceed
      // the memory constraint.
      if (numSamp > limit) {
         errorMsg += "WARNING: memory constraint exceeded\n" +
                     "         make sure " + "\n";
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
