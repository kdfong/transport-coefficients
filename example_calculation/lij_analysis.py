import numpy as np
from scipy import integrate
import os


class TransportCoefficients:
    """
    Class for using Green-Kubo relations to compute transport coefficients L++, L+-, and L-- 
    from LAMMPS simulation data.
    """

    def __init__(self, v_cation_filename, v_anion_filename, V, times, T = 298):
        """
        :param v_anion_filename: string, filename for dump file containing anion velocities
        :param v_cation_filename: string, filename for dump file containing cation velocities
        :param V: float, system volume (cubic Angstroms)
        :param times: [float], times (in femtoseconds) at which velocity data was collected in dump files
        :param T: absolute temperature
        """
        self.V = V
        self.times = times
        self.T = T
        self.v_cation_filename = v_cation_filename
        self.v_anion_filename = v_anion_filename

        # Unit Conversions and constants
        global A2m  # Angstroms to m
        A2m = 1e-10
        global fs2s  # femtosecond to second
        fs2s = 1e-15
        global kb  # Boltzmann constant, J/k
        kb = 1.3806504e-23
        global F  # Faraday's constant, C/mol electron
        F = 96485

    def parse_dump_ions(self, filename):
        """
        Converts LAMMPS dump file to list of sum_i(atom_i velocity) over time.

        In LAMMPS input script, the dump file should be created as follows:
        1) Create a group for the anion or cation of interest
        2) Dump the velocities of each of those ions every time step.

        Example code:
            # Calculate velocities of anions
            group anion type 2  # creat group containing all anions
            dump anion_dump anion custom 1 v_anion.out id vx vy vz  # dump velocities every time step

            # Calculate velocities of cations
            group cation type 3
            dump cation_dump cation custom 1 v_cation.out id vx vy vz

        :param filename: string, the name of the dump file with the ion velocities (e.g. v_anion.out)
        :return: v = [vx, vy, vz], array[float, float, float]: velocity components summed over all ions over time
        """

        v = []  # vx, vy, vx summed over all atoms
        vx = 0
        vy = 0
        vz = 0
        flag = -1  # -1 in between, 1 data collection

        with open(filename, 'r') as f:
            for line in (f):
                if 'ITEM: ATOMS id vx vy vz' in line:  # start of new data collection
                    flag = 1
                    v.append([vx, vy, vz])
                    vx = 0
                    vy = 0
                    vz = 0
                elif 'ITEM: TIMESTEP' in line:  # start of new timestep
                    flag = -1
                elif flag == 1:
                    data = line.split('\n')[0].split(' ')[1:]
                    vx += float(data[0])
                    vy += float(data[1])
                    vz += float(data[2])
        v = np.array(v)

        return v

    def cross_corr(self, x, y):
        """
        Computes cross-correlation function of x and y.
        :param x: array[float], data set 1
        :param y: array[float], data set 2
        :return: array[float], cross-correlation of x and y.
        """
        N = len(x)
        F1 = np.fft.fft(x, n=2 ** (N * 2 - 1).bit_length())  # 2*N because of zero-padding, use next highest power of 2
        F2 = np.fft.fft(y, n=2 ** (N * 2 - 1).bit_length())
        PSD = F1 * F2.conjugate()
        res = np.fft.ifft(PSD)
        res = (res[:N]).real
        n = N * np.ones(N) - np.arange(0, N)
        return res / n

    def compute_acf(self, v_cation, v_anion):
        """
        Computes the correlation function between sum(velocity) of different types of ions.
        :param v_cation: array[float], velocity summed over all cations
        :param v_anion: array[float], velocity summed over all anions
        :return: acf_plusplus, acf_plusminus, acf_minusminus: array[float, float, float], correlation function of
        sum(velocity) for cation-cation, cation-anion, and anion-anion
        """
        n_times = min(len(v_cation), len(v_anion))  # number of times for which data is collected
        # Initialize outputs
        acf_plusplus = np.zeros([n_times, 3])
        acf_plusminus = np.zeros([n_times, 3])
        acf_minusminus = np.zeros([n_times, 3])
        # Compute correlation functions
        for i in range(3):  # three spatial directions, xyz
            acf_plusplus[:, i] = self.cross_corr(v_cation[:, i], v_cation[:, i])
            acf_plusminus[:, i] = self.cross_corr(v_cation[:, i], v_anion[:, i])
            acf_minusminus[:, i] = self.cross_corr(v_anion[:, i], v_anion[:, i])

        return acf_plusplus, acf_plusminus, acf_minusminus

    def compute_lij(self, acf_plusplus, acf_plusminus, acf_minusminus):
        """
        Computes transport coefficients Lij, in units of 1/(J s m) (averaged over all three spatial directions)
        by integrating the associated correlation functions
        :param acf_plusplus: array[float, float, float], correlation function of sum(velocity) for
        cation-cation
        :param acf_plusminus: array[float, float, float], correlation function of sum(velocity) for
        cation-anion
        :param acf_minusminus: array[float, float, float], correlation function of sum(velocity) for
        anion-anion
        :return: L_plusplus, L_plusminus, L_minusminus: array[float], transport coefficients in units of 1/(J s m)
        over time
        """
        convert_lij = 1 / (A2m * fs2s)  # final units 1/(J s m)

        n_times = len(acf_plusplus[:,0])
        L_plusplus = np.zeros([n_times - 1, 3])
        L_plusminus = np.zeros([n_times - 1, 3])
        L_minusminus = np.zeros([n_times - 1, 3])

        for i in range(3):  # three spatial directions, xyz
            L_plusplus[:, i] = 1/(kb*self.T*self.V)*integrate.cumtrapz(acf_plusplus[:,i], self.times)*convert_lij
            L_plusminus[:, i] = 1/(kb*self.T*self.V)*integrate.cumtrapz(acf_plusminus[:,i], self.times)*convert_lij
            L_minusminus[:, i] = 1/(kb*self.T*self.V)*integrate.cumtrapz(acf_minusminus[:,i], self.times)*convert_lij

        # average over all three spatial directions
        L_plusplus_avg = np.mean(L_plusplus, axis=1)
        L_plusminus_avg = np.mean(L_plusminus, axis=1)
        L_minusminus_avg = np.mean(L_minusminus, axis=1)

        return L_plusplus_avg, L_plusminus_avg, L_minusminus_avg

    def compute_all(self):
        """
        Compute all correlation functions and their corresponding transport coefficients (cation-cation, cation-anion,
        anion-anion). Saves each of these as a numpy array (.npy) file.
        :param v_anion_filename: string, filename for dump file containing anion velocities
        :param v_cation_filename: string, filename for dump file containing cation velocities
        :return: acf_plusplus, acf_plusminus, acf_minusminus, L_plusplus_avg, L_plusminus_avg, L_minusminus_avg:
        array[float, float, float],array[float, float, float],array[float, float, float],array[float],array[float],
        array[float]
        """
        v_cation = self.parse_dump_ions(self.v_cation_filename)
        np.save('v_cation.npy', v_cation)
        v_anion = self.parse_dump_ions(self.v_anion_filename)
        np.save('v_anion.npy', v_anion)

        acf_plusplus, acf_plusminus, acf_minusminus = self.compute_acf(v_cation, v_anion)
        np.save('acf_plusplus.npy', acf_plusplus)
        np.save('acf_plusminus.npy', acf_plusminus)
        np.save('acf_minusminus.npy', acf_minusminus)

        L_plusplus_avg, L_plusminus_avg, L_minusminus_avg = self.compute_lij(acf_plusplus, acf_plusminus, acf_minusminus)
        np.save('L_plusplus.npy', L_plusplus_avg)
        np.save('L_plusminus.npy', L_plusminus_avg)
        np.save('L_minusminus.npy', L_minusminus_avg)

        return acf_plusplus, acf_plusminus, acf_minusminus, L_plusplus_avg, L_plusminus_avg, L_minusminus_avg