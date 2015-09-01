"""
workflow_tensor.py

Defines a class **workflow_tensor** that acts as a data container and workflow 
guideline to calculate a non-adiabatic friction tensor and 
vibrational lifetimes.
"""
import coolvib
import coolvib.parser as parser
import coolvib.routines.friction_tensor as tensor
import coolvib.routines.spectral_function as spectral
import coolvib.routines.output as output
from ase.atoms import Atoms

class workflow_tensor():
    def __init__(self, atoms, code='aims', active_atoms=None):
        """
        workflow_tensor object.
        
        The workflow_tensor objet is a basic data container and workflow guideline 
        for all calculations necessary to generate a friction tensor and 
        calculate vibrational lifetimes due to non-adiabatic couplings. 

        Parameters:

        atoms: ase.atoms class 
            defines the atoms, molecules and or surfaces for 
            which the calculation is performed.

        code: str 
            A string that defines from which quantum chemistry package the input 
            is received.

        active_atoms: list
            A list of atom indices indicating for which atoms in the system the 
            friction tensor will be calculated and for which atoms first order 
            matrices exist.
           
        Examples:

            This is how this class is initialized::

                from ase.all import Atoms
                a = Atoms("CO", [(0.,0.,0.), (0.,0.,1.0)])
                model = workflow_tensor(a, code='siesta', active_atoms=[0, 1])

        """
        if isinstance(atoms, Atoms):
            self.atoms = atoms
        else:
            raise ValueError('model needs an instance of an ase Atoms class')
        self.code = code
        self.active_atoms = active_atoms

        self.input_read=False

    def set_code(self, code):
        """
        Sets the code keyword
        """

        if code in codes:
            self.code = code
        else:
            print 'This code is not supported, choose one of ', codes

    def set_active_atoms(self, active_atoms=None):
        """
        Sets the active_atoms list
        """

        self.active_atoms = active_atoms

    def read_input_data(self, **kwargs):
        """
        Reads all necessary input data including first order matrices
        wavefunctions, eigenvalues, fermi_level, etc.

        This is a wrapper function for specific parser routines, depending 
        on which code is used. For detailed input explanations please see 
        the documentation in the :py:coolvib.parser: package.

        """

        if self.code is None:
            print 'Please set code first using .set_code'

        method_to_call =getattr(parser,parser.codes[self.code]+'_tensor')
        self = method_to_call(self, **kwargs)

        self.input_read = True

    def calculate_friction_tensor(self,mode='default',**kwargs):
        """
        Calculates the friction tensor for a given window and delta function, 
        assumes that the spectral function has already been calculated, otherwise 
        calculates it.

        Parameters:

        mode: str

        This is a wrapper function for :py:coolvib.routines.friction_tensor.calculate_tensor:
        For detailed documentation please look there.
        """

        calculate_spec = True
        if hasattr(self, 'spectral_function'):
            calculate_spec = False
        
        #first calculate spectral function
        if calculate_spec:
            self.calculate_spectral_function(mode,**kwargs)

        ###calculate the tensor by applying the delta function
        x_axis = self.x_axis
        n_dim = len(self.active_atoms)*3
        spectral_function = self.spectral_function

        self.friction_tensor = tensor.calculate_tensor(
                n_dim,
                x_axis,
                spectral_function,
                **kwargs)

    def calculate_spectral_function(self,mode='default', **kwargs):
        """
        initiates calculation of the spectral functions and 
        takes the same arguments as coolvib.routines.calculate_spectral_function
     
        Parameters:

        mode : str
            'default' or 'momentum' for q=0 or q not 0 versions
        
        """

        masses = self.atoms.get_masses()[self.active_atoms]
        #CASE CODE
        if parser.code_type[self.code] is 'local':
            #CASE mode
            if mode is 'default':
                self.x_axis, self.spectral_function = spectral.calculate_spectral_function_tensor(
                        self.fermi_energy,
                        self.eigenvalues,
                        self.kpoints,
                        self.psi,
                        self.first_order_H,
                        self.first_order_S,
                        masses,
                        **kwargs)
            elif mode is 'momentum':
                self.x_axis, self.spectral_function = spectral.calculate_spectral_function_tensor_q(
                        self.fermi_energy,
                        self.eigenvalues,
                        self.kpoints,
                        self.psi,
                        self.first_order_H,
                        self.first_order_S,
                        masses,
                        self.basis_pos,
                        **kwargs)

        else:
            #plane wave stuff
            raise NotImplementedError('calculate_spectral_function: code is not supported!')

    def analyse_friction_tensor(self):
        """
        diagonalizes friction tensor and does spectral analysis
        """

        if hasattr(self, 'friction_tensor'):
            pass
        else:
            self.calculate_friction_tensor()

        self.friction_eigenvectors, self.friction_eigenvalues = \
                tensor.analyse_tensor(
                self.friction_tensor)


    def print_spectral_function(self, filename='nacs-spectrum.out'):
        """
        prints spectral function to a file

        Parameters:

        filename: str
            name of the output file, default is nacs-spectrum.out
        """

        output.print_spectral_function(
                self.x_axis, self.spectral_function,len(self.active_atoms)*3,
                filename)
    

    def print_jmol_friction_eigenvectors(self, filename='friction-eigenvectors.jmol'):
        """
        prints friction eigenvectors in jmol format
        
        Parameters:

        filename: str
            name of the output file, default is friction-eigenvectors.jmol 
        """
        
        if hasattr(self, 'friction_eigenvectors'):
            pass
        else:
            self.analyse_friction_tensor()

        output.print_jmol_friction_eigenvectors(
                self.atoms, self.active_atoms, 
                self.friction_eigenvectors,
                self.friction_eigenvalues,
                filename)

