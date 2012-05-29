import java.lang.System;

import org.python.util.PythonInterpreter;

public class Main {
    public static void main(String[] args) {
        //Put the command line arguments ( args ) into sys.argv
        PythonInterpreter.initialize(System.getProperties(), System.getProperties(), args);
        // Now create an interpreter
        PythonInterpreter interp = new PythonInterpreter();
        interp.exec("try:\n import visualitzador\n visualitzador.main() \nexcept SystemExit: import sys\n sys.exit(0)");
    }
}
