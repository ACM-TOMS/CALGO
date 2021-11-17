/*******************************************************************************
   GapPanel class
*******************************************************************************/


//******************************************************************************
// GapPanel class defines a panel for entering and translating parameters of
// "Gap Test".
//******************************************************************************

public class GapPanel
       extends ParmGrpPanel
{
   //---------------------------------------------------------------------------
   // Parameter field for maximum gap length.
   protected IntParmField  maxGapLenField;
   // Parameter field for number of samples.
   protected IntParmField  numSampField;
   // Parameter fields for the lower and upper bound.
   protected DecParmField  lowerBoundField;
   protected DecParmField  upperBoundField;

   // The command line translation of the parameters.
   protected String  parmStr;
   // Owner of this panel.
   protected TestWizard  wiz;

   //***************************************************************************
   // The constructor constructs PermPanel with owner <wizard>.
   //---------------------------------------------------------------------------
   public GapPanel(TestWizard wizard) {
      // Set up a group panel that will contain four parameter fields.
      super("Gap Test", 4);
      // Set the command line id.
      super.cmdId = "gap";

      // Create the parameter fields.
      maxGapLenField  = new IntParmField("Maximum Gap Length",
                                        2, Integer.MAX_VALUE,20);
      numSampField    = new IntParmField("    Number of Gaps",
                                        1, Integer.MAX_VALUE,100000);
      lowerBoundField = new DecParmField("       Lower Bound", 0.0, 1.0, 0.5);
      upperBoundField = new DecParmField("       Upper Bound", 0.0, 1.0, 0.6);

      // Initialize internal variables.
      parmStr = "";
      wiz = wizard;

      // Add the components.
      addField(maxGapLenField);
      addField(lowerBoundField);
      addField(upperBoundField);
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
      long    maxGapLen;
      double  lowerBound;
      double  upperBound;
      long    numSamp;

      String  errorMsg;
      double  prob;
      double  binProb;
      long    limit;

      // Check parameters for errors and concatenate error messages.
      errorMsg = maxGapLenField.checkParm() + lowerBoundField.checkParm() +
                 upperBoundField.checkParm() + numSampField.checkParm();
      // If error exists, set command line translation to empty string and
      // return the error message.
      if (errorMsg.length() > 0) {
         parmStr = "";
         return errorMsg;
      }

      // Get the numerical values of the parameters.
      maxGapLen = maxGapLenField.getParmValue();
      lowerBound = lowerBoundField.getParmValue();
      upperBound = upperBoundField.getParmValue();
      numSamp = numSampField.getParmValue();

      // Set command line translation to empty string and return an error
      // message if the lower bound is not smaller than the upper bound.
      if (lowerBound >= upperBound) {
         parmStr = "";
         return "ERROR: make sure " + lowerBoundField.getParmTitle() +
                " < " + upperBoundField.getParmTitle();
      }

      // Translate the parameters into command line segment and cache it.
      parmStr = maxGapLenField.getParmString() + " " +
                lowerBoundField.getParmString() + " " +
                upperBoundField.getParmString() + " " +
                numSampField.getParmString();

      // Find the probability of the least likely bin.
      prob = upperBound - lowerBound;
      if (prob < 0.5) {
         binProb = Math.pow((1-prob), maxGapLen) * prob;
      } else {
         binProb = Math.pow((1-prob), maxGapLen+1);
      }
      // Create a warning message if the sampling requirement "minimum of 5
      // in the least likely bin" is not met.
      if (binProb * numSamp < 5) {
         errorMsg += "WARNING: sampling requirement not met\n" +
                     "         make sure the least likely bin will have " +
                     "at least 5 samples on the average\n";
      }

      // Create a warning message if memory used by the test will exceed
      // the memory constraint.
      limit = wiz.maxTestMemory * wiz.KILO_BYTE / 4;
      if (maxGapLen+2 > limit) {
         errorMsg += "WARNING: memory constraint exceeded\n" +
                     "         make sure " +
                     maxGapLenField.getParmTitle() + " <= " +
                     (limit-2) + "\n";
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
