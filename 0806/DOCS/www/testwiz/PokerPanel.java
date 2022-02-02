/*******************************************************************************
   PokerPanel class
*******************************************************************************/


//******************************************************************************
// PokerPanel class defines a panel for entering and translating parameters of
// "Poker Test".
//******************************************************************************

public class PokerPanel
       extends ParmGrpPanel
{
   //---------------------------------------------------------------------------
   // Parameter field for block length.
   protected IntParmField  blockLengthField;
   // Parameter field for number of samples.
   protected IntParmField  numSampField;
   // Parameter field for range of integers to generate.
   protected IntParmField  maxNumField;

   // The command line translation of the parameters.
   protected String  parmStr;
   // Owner of this panel.
   protected TestWizard  wiz;

   //***************************************************************************
   // The constructor constructs PermPanel with owner <wizard>.
   //---------------------------------------------------------------------------
   public PokerPanel(TestWizard wizard) {
      // Set up a group panel that will contain three parameter fields.
      super("Poker Test", 3);
      // Set the command line id.
      super.cmdId = "poker";

      // Create the parameter fields.
      numSampField   = new IntParmField("Numbers of blocks",
                                        2, Integer.MAX_VALUE,100000);
      blockLengthField = new IntParmField("Length of each block",
                                        2, Integer.MAX_VALUE,10);
      maxNumField = new IntParmField("Generate integers in [0,?).",
                                        2, Integer.MAX_VALUE,10);

      // Initialize internal variables.
      parmStr = "";
      wiz = wizard;

      // Add the components.
      addField(numSampField);
      addField(blockLengthField);
      addField(maxNumField);
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
      long  maxNum;
      long  blockLength;

      String  errorMsg;
      long  limit;

      // Check parameters for errors and concatenate error messages.
      errorMsg = numSampField.checkParm() + blockLengthField.checkParm() + 
	maxNumField.checkParm();
      // If error exists, set command line translation to empty string and
      // return the error message.
      if (errorMsg.length() > 0) {
         parmStr = "";
         return errorMsg;
      }

      // Translate the parameters into command line segment and cache it.
      parmStr = numSampField.getParmString() + " " +
                blockLengthField.getParmString() + " " +
                maxNumField.getParmString();

      // Get the numerical values of the parameters.
      numSamp = numSampField.getParmValue();
      maxNum = maxNumField.getParmValue();
      blockLength = blockLengthField.getParmValue();

      limit = wiz.maxTestMemory * wiz.KILO_BYTE / 4;
      // Create a warning message if memory used by the test will exceed
      // the memory constraint.
      if (blockLength > limit) {
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
