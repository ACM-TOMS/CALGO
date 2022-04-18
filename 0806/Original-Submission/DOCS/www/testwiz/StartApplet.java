/*******************************************************************************
   StartApplet class
*******************************************************************************/
import java.applet.*;
import java.awt.*;

//******************************************************************************
// StartApplet class provides the front-end where user can invoke TestWizard.
//******************************************************************************

public class StartApplet
       extends Applet
{
   //---------------------------------------------------------------------------
   // Button for starting TestWizard.
   protected Button  startButton;

   //***************************************************************************
   // init() displays a button in the applet that will start TestWizard. 
   //---------------------------------------------------------------------------
   public void init() {
      startButton = new Button("Start "+TestWizard.APP_NAME);

      startButton.setFont(new Font("Courier", Font.ITALIC, 16));
      startButton.setBackground(Color.pink);

      add(startButton);
   }

   //---------------------------------------------------------------------------
   // action() catches events on the button, then starts TestWizard.
   //---------------------------------------------------------------------------
   public boolean action(Event eve, Object arg) {
      if (eve.target == startButton) {
         TestWizard  twiz = new TestWizard();

         twiz.setFont(new Font("Courier", Font.PLAIN, 16));
         twiz.resize(500,200);
         twiz.show();

         return true;
      }

      return super.action(eve, arg);
   }
}
