import java.lang.System;

import org.python.util.PythonInterpreter;

public class Main {
    public static void main(String[] args) {
        System.setProperty("python.home","/");
        System.setProperty("python.import.site","false");
        //Put the command line arguments ( args ) into sys.argv
        PythonInterpreter.initialize(System.getProperties(), System.getProperties(), args);
        // Now create an interpreter
        PythonInterpreter interp = new PythonInterpreter();
        interp.exec("try:\n import visualitzador\n visualitzador.main() \nexcept SystemExit: pass");
    }
}
