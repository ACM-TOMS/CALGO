/*******************************************************************************
   PermPanel class
*******************************************************************************/


//******************************************************************************
// PermPanel class defines a panel for entering and translating parameters of
// "Permutation Test".
//******************************************************************************

public class PermPanel
       extends ParmGrpPanel
{
   //---------------------------------------------------------------------------
   // Parameter field for group length.
   protected IntParmField  grpLenField;
   // Parameter field for number of samples.
   protected IntParmField  numSampField;

   // The command line translation of the parameters.
   protected String  parmStr;
   // Owner of this panel.
   protected TestWizard  wiz;

   //***************************************************************************
   // The constructor constructs PermPanel with owner <wizard>.
   //---------------------------------------------------------------------------
   public PermPanel(TestWizard wizard) {
      // Set up a group panel that will contain two parameter fields.
      super("Permutation Test", 2);
      // Set the command line id.
      super.cmdId = "perm";

      // Create the parameter fields.
      grpLenField  = new IntParmField("Length of Each Subsequence",
                                      2, Integer.MAX_VALUE,5);
      numSampField = new IntParmField("    Number of Subsequences",
                                      1, Integer.MAX_VALUE,200000);

      // Initialize internal variables.
      parmStr = "";
      wiz = wizard;

      // Add the components.
      addField(grpLenField);
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
      long  grpLen;
      long  numSamp;

      String  errorMsg;
      double  numBins;
      long  limit;

      // Check parameters for errors and concatenate error messages.
      errorMsg = grpLenField.checkParm() + numSampField.checkParm();
      // If error exists, set command line translation to empty string and
      // return the error message.
      if (errorMsg.length() > 0) {
         parmStr = "";
         return errorMsg;
      }

      // Translate the parameters into command line segment and cache it.
      parmStr = grpLenField.getParmString() + " " +
                numSampField.getParmString();

      // Get the numerical values of the parameters.
      grpLen = grpLenField.getParmValue();
      numSamp = numSampField.getParmValue();

      // Use Sterling approximation to calculate the number of bins.
      numBins = Math.sqrt(2*Math.PI*grpLen) * Math.pow(grpLen/Math.E, grpLen);

      // Create a warning message if the sampling requirement "minimum of 5
      // in the least likely bin" is not met.
      if (numSamp/numBins < 5) {
         errorMsg += "WARNING: sampling requirement not met\n" +
                     "         make sure " +
                     numSampField.getParmTitle() + " >= (" +
                     grpLenField.getParmTitle() + ")! * 5\n";
      }

      limit = wiz.maxTestMemory * wiz.KILO_BYTE / 4;
      // Create a warning message if memory used by the test will exceed
      // the memory constraint.
      if (numBins > limit) {
         errorMsg += "WARNING: memory constraint exceeded\n" +
                     "         make sure (" +
                     grpLenField.getParmTitle() + ")! <= " +
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
