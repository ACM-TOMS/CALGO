/*******************************************************************************
   CouponPanel class
*******************************************************************************/


//******************************************************************************
// CouponPanel class defines a panel for entering and translating parameters of
// "Coupon Test".
//******************************************************************************

public class CouponPanel
       extends ParmGrpPanel
{
   //---------------------------------------------------------------------------
   // Parameter field for maximum sequence length.
   protected IntParmField  maxSeqField;
   // Parameter field for largest integer+1.
   protected IntParmField  maxNumField;
   // Parameter field for number of samples.
   protected IntParmField  numSampField;

   // The command line translation of the parameters.
   protected String  parmStr;
   // Owner of this panel.
   protected TestWizard  wiz;

   //***************************************************************************
   // The constructor constructs PermPanel with owner <wizard>.
   //---------------------------------------------------------------------------
   public CouponPanel(TestWizard wizard) {
      // Set up a group panel that will contain three parameter fields.
      super("Coupon Test", 3);
      // Set the command line id.
      super.cmdId = "coupon";

      // Create the parameter fields.
      numSampField   = new IntParmField("Numbers of complete sets",
                                        2, Integer.MAX_VALUE,20000);
      maxSeqField = new IntParmField("Lump together sets larger than:",
                                        2, Integer.MAX_VALUE,30);
      maxNumField = new IntParmField("Generate integers in [0,?).",
                                        2, Integer.MAX_VALUE,10);

      // Initialize internal variables.
      parmStr = "";
      wiz = wizard;

      // Add the components.
      addField(numSampField);
      addField(maxSeqField);
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
      long  maxSeq;

      String  errorMsg;
      long  limit;

      // Check parameters for errors and concatenate error messages.
      errorMsg = numSampField.checkParm() + maxSeqField.checkParm() + 
	maxNumField.checkParm();
      // If error exists, set command line translation to empty string and
      // return the error message.
      if (errorMsg.length() > 0) {
         parmStr = "";
         return errorMsg;
      }

      // Translate the parameters into command line segment and cache it.
      parmStr = numSampField.getParmString() + " " +
                maxSeqField.getParmString() + " " +
                maxNumField.getParmString();

      // Get the numerical values of the parameters.
      numSamp = numSampField.getParmValue();
      maxNum = maxNumField.getParmValue();
      maxSeq = maxSeqField.getParmValue();

      limit = wiz.maxTestMemory * wiz.KILO_BYTE / 4;
      // Create a warning message if memory used by the test will exceed
      // the memory constraint.
      if (maxSeq > limit) {
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
