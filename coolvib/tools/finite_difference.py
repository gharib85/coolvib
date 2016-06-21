import os
from ase.vibrations import Vibrations
from ase.utils import opencew
from os import remove
from os.path import isfile, getsize

class finite_difference(Vibrations):
    """
    This class is derived from the ASE Vibrations class and only 
    slightly modified to output Hamiltonian and overlap matrix output into 
    specific folders. With this data coolvib can calculate the finite difference 
    nonadiabatic coupling elements.
    """


    def run(self, clean=False):
        """Run the vibration calculations.

        This will calculate the forces for 6 displacements per atom +/-x,
        +/-y, +/-z. Only those calculations that are not already done will be
        started. Be aware that an interrupted calculation may produce an empty
        file (ending with .pckl), which must be deleted before restarting the
        job. Otherwise the forces will not be calculated for that
        displacement.

        Note that the calculations for the different displacements can be done
        simultaneously by several independent processes. This feature relies
        on the existence of files and the subsequent creation of the file in
        case it is not found.
        """

        filename = self.name + '.eq.pckl'
        fd = opencew(filename)
        if fd is not None:
            self.calculate(filename, fd)
            try:
                os.mkdir('./eq')
            except:
                pass
            os.system('mv *.out ./eq/')

        p = self.atoms.positions.copy()
        for a in self.indices:
            for i in range(3):
                for sign in [-1, 1]:
                    for ndis in range(1, self.nfree // 2 + 1):
                        filename = ('%s.%d%s%s.pckl' %
                                    (self.name, a, 'xyz'[i],
                                     ndis * ' +-'[sign]))
                        filename2 = ('a%dc%d%s' %
                                    (a, i,
                                     ndis * ' +-'[sign]))

                        if (isfile(filename) and getsize(filename) == 0):
                            remove(filename)
                        fd = opencew(filename)
                        if fd is not None:
                            disp = ndis * sign * self.delta
                            self.atoms.positions[a, i] = p[a, i] + disp
                            self.calculate(filename, fd)
                            try:
                                os.mkdir(filename2)
                            except:
                                pass
                            os.system('mv *.out '+filename2+'/')
                            self.atoms.positions[a, i] = p[a, i]

        if clean:
            self.clean()
