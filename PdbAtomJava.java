
public class PdbAtomJava{
    //Represents an atom from a PDB file
    public String residue;
    public String hetid;
    public float[] xyz;
    public PdbAtomJava(String record){
        //Needs an ATOM or HETATM record
        residue = record.substring(17,27);
        hetid = residue.substring(0,3);
        hetid = hetid.trim();
        xyz =  new float[3];
        xyz[0] = Float.parseFloat( record.substring(30,38));
        xyz[1] = Float.parseFloat( record.substring(38,46));
        xyz[2] = Float.parseFloat( record.substring(46,54));
    }
    public double __or__(PdbAtomJava other){
        //Return squared distance
        double distance = Math.pow(this.xyz[0] - other.xyz[0],2) + Math.pow(this.xyz[1] - other.xyz[1],2) + Math.pow(this.xyz[2] - other.xyz[2],2);
        return distance;
    }
}
