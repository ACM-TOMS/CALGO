/*******************************************************************************
   IntParmField class
*******************************************************************************/


//******************************************************************************
// IntParmField class defines a parameter field for entering integers.
//******************************************************************************

public class IntParmField
       extends ParmField
{
   //---------------------------------------------------------------------------
   // Minimum and maximum on allowed values.
   protected long  minValue;
   protected long  maxValue;

   // Default value.
   protected long  defaultValue;

   // Value of the parameter.
   protected long  parmVal;

   //***************************************************************************
   // The constructor constructs IntParmField with title <parmName>, minimum
   // allowed value <min>, maximum allowed value <max>, and default value.
   //---------------------------------------------------------------------------
   public IntParmField(String parmName, long min, long max, long def) {
      super(parmName);

      minValue = min;
      maxValue = max;
      defaultValue = def;
   }

   //---------------------------------------------------------------------------
   // The constructor constructs IntParmField with title <parmName>, minimum
   // allowed value <min>, maximum allowed value <max>, and default value.
   //---------------------------------------------------------------------------
   public IntParmField(String parmName, long min, long max) {
      this(parmName, min, max, 0);
   }

   //---------------------------------------------------------------------------
   // Assigns the default value to the parameter.
   //---------------------------------------------------------------------------
   public void setDefault() {
      super.entryField.setText(Long.toString(defaultValue));
      // Update the internal variable.
      parmVal = defaultValue;
   }

   //---------------------------------------------------------------------------
   // checkParm() checks the validity of the parameter value, and returns a
   // description of the error.  Empty string indicates no error.
   //---------------------------------------------------------------------------
   public String checkParm() {
      try {
         // Convert the parameter string into a long integer.
         parmVal = Long.parseLong(super.entryField.getText());
      } catch (NumberFormatException e) {
         // The parameter string cannot be interpreted as a long integer.
         return "ERROR: " + getParmTitle() + " must be an integer\n";
      }

      // Check if the parameter value is outside the max/min limits.
      if ((parmVal < minValue) || (parmVal > maxValue)) {
         return "ERROR: make sure " + minValue + " <= " +
                getParmTitle() + " <= " + maxValue + "\n";
      }

      return "";
   }

   //---------------------------------------------------------------------------
   // setParm() sets the parameter to <theValue>.
   //---------------------------------------------------------------------------
   public void setParm(long theValue) {
      super.entryField.setText(Long.toString(theValue));
      parmVal = theValue;
   }

   //---------------------------------------------------------------------------
   // getParmValue() returns the parameter as a long integer.
   // NOTE: Before calling getParmValue(), call checkParm() first.
   //---------------------------------------------------------------------------
   public long getParmValue() {
      return parmVal;
   }
}
