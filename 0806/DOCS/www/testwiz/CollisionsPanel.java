/*******************************************************************************
   CollisionsPanel class
*******************************************************************************/


//******************************************************************************
// CollisionsPanel class defines a panel for entering and translating parameters of
// "Collisions Test".
//******************************************************************************

public class CollisionsPanel
       extends ParmGrpPanel
{
   //---------------------------------------------------------------------------
   // Parameter field for number of digits.
   protected IntParmField  logmdField;
   // Parameter field for log of base of each digit.
   protected IntParmField  logdField;
   // Parameter field for number of samples.
   protected IntParmField  numSampField;

   // The command line translation of the parameters.
   protected String  parmStr;
   // Owner of this panel.
   protected TestWizard  wiz;

   //***************************************************************************
   // The constructor constructs PermPanel with owner <wizard>.
   //---------------------------------------------------------------------------
   public CollisionsPanel(TestWizard wizard) {
      // Set up a group panel that will contain three parameter fields.
      super("Collisions Test", 3);
      // Set the command line id.
      super.cmdId = "collisions";

      // Create the parameter fields.
      numSampField   = new IntParmField("Form ? numbers per test by ...",
                                        2, Integer.MAX_VALUE,25000);
      logdField = new IntParmField("... concatenating ? bits from ...",
                                        1, Integer.MAX_VALUE,5);
      logmdField = new IntParmField("... ? random numbers each.",
                                        1, Integer.MAX_VALUE,4);

      // Initialize internal variables.
      parmStr = "";
      wiz = wizard;

      // Add the components.
      addField(numSampField);
      addField(logdField);
      addField(logmdField);
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
      long  logd;
      long  logmd;

      String  errorMsg;
      long  limit;

      // Check parameters for errors and concatenate error messages.
      errorMsg = numSampField.checkParm() + logdField.checkParm() + 
	logmdField.checkParm();
      // If error exists, set command line translation to empty string and
      // return the error message.
      if (errorMsg.length() > 0) {
         parmStr = "";
         return errorMsg;
      }

      // Translate the parameters into command line segment and cache it.
      parmStr = numSampField.getParmString() + " " +
                logmdField.getParmString() + " " +
                logdField.getParmString();

      // Get the numerical values of the parameters.
      numSamp = numSampField.getParmValue();
      logmd = logmdField.getParmValue();
      logd = logdField.getParmValue();

      limit = wiz.maxTestMemory * wiz.KILO_BYTE / 4;
      // Create a warning message if memory used by the test will exceed
      // the memory constraint.
      if ((numSamp+(1<<(logmd*logd))) > limit) {
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
