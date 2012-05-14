import javax.swing.JOptionPane;
 
public class WrapJOptionPane{
        public static JOptionPane JOptionPane2(Object[] message, int msgtype, int buttons){       
                JOptionPane pane = new JOptionPane( message, msgtype, buttons);
                return pane;
        }
}


