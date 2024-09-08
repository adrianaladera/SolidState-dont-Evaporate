import numpy as np

def make_kpoints(path):
    with open("{}/POSCAR".format(path), 'r') as f:

        file = f.readlines()
        a1 = np.array([float(file[2].split()[0]), float(file[2].split()[1]), float(file[2].split()[2])])
        a2 = np.array([float(file[3].split()[0]), float(file[3].split()[1]), float(file[3].split()[2])])
        a3 = np.array([float(file[4].split()[0]), float(file[4].split()[1]), float(file[4].split()[2])])

        V = np.linalg.norm(np.dot(a1, np.cross(a2, a3)))
        len_b1 = np.linalg.norm(((2*np.pi) * np.cross(a2, a3)) / V)
        len_b2 = np.linalg.norm(((2*np.pi) * np.cross(a3, a1)) / V)
        len_b3 = np.linalg.norm(((2*np.pi) * np.cross(a1, a2)) / V)

        N1 = round(len_b1 * (1 / min([len_b1, len_b2, len_b3])))
        N2 = round(len_b2 * (1 / min([len_b1, len_b2, len_b3])))
        N3 = round(len_b3 * (1 / min([len_b1, len_b2, len_b3])))

        with open("{}/KPOINTS".format(path), 'w') as kp:
            kp.write("Regular k-point mesh\n")
            kp.write("0\n")
            kp.write("Gamma\n")
            kp.write("{} {} {}\n".format(N1, N2, N3))
        kp.close()
    f.close()

x = input("Enter path for POSCAR: ")
make_kpoints(x)
