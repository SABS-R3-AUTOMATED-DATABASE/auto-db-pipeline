'''
    Created on 17 Apr 2013
    
    @author: dunbar

    @change: Modified the visualise function to take a list of entities as input and show them all in the same session. 
    The session will be associated with the first entity in the list 
    
'''

from subprocess import PIPE, Popen
import tempfile, os
from functools import partial

class Pymol:
    """
    A pymol instance class. 
    Adapted from:
     http://svn.gna.org/svn/relax/1.3/generic_fns/pymol_control.py
    Relax NMR program
    """

    class Cmd:
        """
        Class to allow the use of the pymol api in the subprocess.
        It mimicks pymol.cmd.
        """
        def __init__(self,pipe):
            self.pipe=pipe
        
        def __getattr__(self,f):
            # partial allows you to pass the function name to the _interpreter
            return partial(self._interpreter, function=f)

        def _interpreter(self,*args,**kwargs):

            cmd_str = "%s " % kwargs["function"]
            
            if args:
                cmd_str += " , ".join(args)
            
            for kw in kwargs:
                if kw == "function":continue
                cmd_str += " , "
                cmd_str += "%s=%s" % (str(kw), str(kwargs[kw]))
            self.pipe.write(cmd_str+"\n")
        
    def exec_cmd(self, command=None):
        """
        Execute a PyMOL command.
        @param command:         The PyMOL command to send into the program.
        @type command:          str
        """

        # Reopen the GUI if needed.
        if not self.running():
            self.open_gui()
        self.pymol.write(command + '\n')


    def open_gui(self):
        """Open the PyMOL GUI."""

        # Open PyMOL as a pipe.
        self.pymol = Popen(['pymol', '-qpK'], stdin=PIPE).stdin
        # allow for cmd to still be use
        self.cmd = self.Cmd(self.pymol)

    def load_pdb(self, file):
        """
        Load a pdb into the session.
        """

        # Test if PyMOL is running.
        if not self.running():
            return

        # Open the file in PyMOL.
        self.exec_cmd("load " + file)


    def running(self):
        """
        Test if PyMOL is running.

        @return:    Whether the Pymol pipe
        @rtype:     bool
        """

        # Test if command line PyMOL is already running.
        if not hasattr(self, 'pymol'):
            return False

        # Test if the pipe has been broken.
        try:
            self.pymol.write('\n')
        except IOError:
            return False

        # PyMOL is running.
        return True

    def quit(self):
        """
        Quit pymol.
        """
        self.exec_cmd("quit")
        del self.pymol
    
   
    
def visualise(ent, return_session=False):
    """
    Function to visualise an entity or entities.
    
    @param ent: A AbPDB entity (e.g. structure, ABchain, Fab, Fragment...) or a list of entities.
    """
    pymol_open =False
    if not isinstance(ent,list):
        entities = [ent]
    else:
        entities = ent
    for entity in entities:
        entity_name=""
        p=entity.get_parent()
        while 1:
            if p:
                entity_name=str(p.id)+"_"
                p=p.get_parent()
            else:
                break
        if entity.level == "AS":
            entity_name += entity.id
        elif entity.level == "C":
            entity_name += "chain_" + entity.id    
        elif entity.level == "R":
            entity_name += "residue_" + entity.parent.id + str(entity.id[1]) + entity.id[2]    
        elif entity.level == "F":
            try:
                entity_name += "fragment_" + entity.id + "_chain_" + entity.parent.parent.id
            except AttributeError:
                entity_name += "fragment_" + entity.id 
        elif entity.level == "H":
            entity_name += "chains_" + "_".join(list(entity.child_dict.keys()))
        elif entity.level == "M":
            entity_name += "model_"+str(entity.id)
    
        fd, file_name = tempfile.mkstemp('.pdb', entity_name ) # somehow have to clear up.
        
        entity.save( file_name )
            
        if not pymol_open:
            entity.pymol = Pymol()
            pymol_open = entity.pymol
            entity.pymol.open_gui()
        
        pymol_open.cmd.load(file_name, object=entity_name)
        #entity.pymol.cmd.load(file_name, object=entity_name)
    if return_session:
        return pymol_open
    
    
if __name__ == "__main__":
    p = Pymol()
    p.open_gui()
    p.cmd.load("/homes/dunbar/Documents/Current/2bdl.pdb")
    p.cmd.color("green","2bdl")
    p.cmd.hide("everything")
    p.cmd.show("cartoon")
    