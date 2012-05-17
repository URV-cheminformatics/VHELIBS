import java.lang.System;
import java.util.Properties;

import org.python.core.PySystemState;
import org.python.core.PyException;
import org.python.util.PythonInterpreter;

public class Main {
    public static void main(String[] args) throws PyException {
        //Put the command line arguments ( args ) into sys.argv
        PySystemState.initialize(PySystemState.getBaseProperties(), new Properties(), args);
        // Now create an interpreter
        PythonInterpreter interp = new PythonInterpreter();
        interp.exec("try:\n import visualitzador\n visualitzador.main() \nexcept SystemExit: pass");
        System.exit(0);
    }
}
