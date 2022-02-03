/*******************************************************************************
   DecParmField class
*******************************************************************************/


//******************************************************************************
// DecParmField class defines a parameter field for entering decimal numbers.
//******************************************************************************

public class DecParmField
       extends ParmField
{
   //---------------------------------------------------------------------------
   // Minimum and maximum on allowed values.
   protected double  minValue;
   protected double  maxValue;

   // Default value.
   protected double  defaultValue;

   // Value of the parameter.
   protected double  parmVal;

   //***************************************************************************
   // The constructor constructs DecParmField with title <parmName>, minimum
   // allowed value <min>, maximum allowed value <max> and default value.
   //---------------------------------------------------------------------------
  public DecParmField(String parmName, double min, double max,double def) {
      super(parmName);

      minValue = min;
      maxValue = max;

      defaultValue = def;
   }

   //---------------------------------------------------------------------------
   // The constructor constructs DecParmField with title <parmName>, minimum
   // allowed value <min>, and maximum allowed value <max>.
   //---------------------------------------------------------------------------
   public DecParmField(String parmName, double min, double max) {
      this(parmName, min, max, 0.0);
   }

   //---------------------------------------------------------------------------
   // Assigns the default value to the parameter.
   //---------------------------------------------------------------------------
   public void setDefault() {
      super.entryField.setText(Double.toString(defaultValue));
      // Update the internal variable.
      parmVal = defaultValue;
   }

   //---------------------------------------------------------------------------
   // checkParm() checks the validity of the parameter value, and returns a
   // description of the error.  Empty string indicates no error.
   //---------------------------------------------------------------------------
   public String checkParm() {
      try {
         // Convert the parameter string into a double.
         parmVal = Double.valueOf(super.entryField.getText())
                                 .doubleValue();
      } catch (NumberFormatException e) {
         // The parameter string cannot be interpreted as a double.
         return "ERROR: " + getParmTitle() + " must be a real\n";
      }

      // Check if the parameter value is outside the max/min limits.
      if ((parmVal < minValue) || (parmVal >= maxValue)) {
         return "ERROR: make sure " + minValue + " <= " +
                getParmTitle() + " < " + maxValue + "\n";
      }

      return "";
   }

   //---------------------------------------------------------------------------
   // setParm() sets the parameter to <theValue>.
   //---------------------------------------------------------------------------
   public void setParm(double theValue) {
      // Update the screen.
      super.entryField.setText(Double.toString(theValue));
      // Update the internal variable.
      parmVal = theValue;
   }

   //---------------------------------------------------------------------------
   // getParmValue() returns the parameter as a double.
   // NOTE: Before calling getParmValue(), call checkParm() first.
   //---------------------------------------------------------------------------
   public double getParmValue() {
      return parmVal;
   }
}
