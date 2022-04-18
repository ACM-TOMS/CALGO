/*******************************************************************************
   RunsPanel class
*******************************************************************************/


//******************************************************************************
// RunsPanel class defines a panel for entering and translating parameters of
// "Runs Test".
//******************************************************************************

public class RunsPanel
       extends ParmGrpPanel
{
   //---------------------------------------------------------------------------
   // Parameter field for maximum run length.
   protected IntParmField  maxRunLenField;
   // Parameter field for number of samples.
   protected IntParmField  numSampField;

   // The command line translation of the parameters.
   protected String  parmStr;
   // Owner of this panel.
   protected TestWizard  wiz;

   //***************************************************************************
   // The constructor constructs PermPanel with owner <wizard>.
   //---------------------------------------------------------------------------
   public RunsPanel(TestWizard wizard) {
      // Set up a group panel that will contain two parameter fields.
      super("Runs Test", 2);
      // Set the command line id.
      super.cmdId = "runs";

      // Create the parameter fields.
      maxRunLenField = new IntParmField("Maximum Run Length",
                                        2, Integer.MAX_VALUE,7);
      numSampField   = new IntParmField("    Number of Runs",
                                        1, Integer.MAX_VALUE,600000);

      // Initialize internal variables.
      parmStr = "";
      wiz = wizard;

      // Add the components.
      addField(maxRunLenField);
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
      long  maxRunLen;
      long  numSamp;

      String  errorMsg;
      double  binProb;
      long  limit;

      // Check parameters for errors and concatenate error messages.
      errorMsg = maxRunLenField.checkParm() + numSampField.checkParm();
      // If error exists, set command line translation to empty string and
      // return the error message.
      if (errorMsg.length() > 0) {
         parmStr = "";
         return errorMsg;
      }

      // Translate the parameters into command line segment and cache it.
      parmStr = maxRunLenField.getParmString() + " " +
                numSampField.getParmString();

      // Get the numerical values of the parameters.
      maxRunLen = maxRunLenField.getParmValue();
      numSamp = numSampField.getParmValue();

      // <maxRunLen>! is calculated by using Sterling approximation.
      binProb = 1 / (Math.sqrt(2*Math.PI*maxRunLen) *
                     Math.pow(maxRunLen/Math.E, maxRunLen)) / (maxRunLen+1);
      // Create a warning message if the sampling requirement "minimum of 5
      // in the least likely bin" is not met.
      if (numSamp*binProb < 5) {
         errorMsg += "WARNING: sampling requirement not met\n" +
                     "         make sure the least likely bin will have " +
                     "at least 5 samples on the average\n";
      }

      limit = wiz.maxTestMemory * wiz.KILO_BYTE / 4;
      // Create a warning message if memory used by the test will exceed
      // the memory constraint.
      if (maxRunLen+1 > limit) {
         errorMsg += "WARNING: memory constraint exceeded\n" +
                     "         make sure " +
                     maxRunLenField.getParmTitle() + " <= " +
                     (limit-1) + "\n";
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
