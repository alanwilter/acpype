"""
Store extra information for Atom's elements
"""
# copied from cpptraj's code "Atom.cpp"
# /** Values taken from 'http://www.webelements.com/' */

__all__ = ['mass_element_dict', 'mass_atomic_number_dict']

mass_arr = [
    1.0, 1.00794, 10.811, 12.0107, 14.0067, 15.9994, 18.9984032, 30.973762,
    32.065, 35.453, 79.904, 55.845, 40.078, 126.90447, 24.3050, 63.546, 6.941,
    39.0983, 85.4678, 132.9054519, 65.38, 22.98976928, 26.9815386, 39.948,
    74.92160, 107.8682, 196.966569, 210, 9.012182, 137.327, 208.98040, 51.9961,
    58.933195, 112.411, 223, 69.723, 72.64, 4.002602, 178.49, 200.59, 114.818,
    192.217, 83.798, 54.938045, 95.96, 20.1797, 58.6934, 92.90638, 190.23,
    106.42, 195.084, 207.2, 209, 101.07, 102.90550, 186.207, 222, 226, 28.0855,
    44.955912, 78.96, 87.62, 118.710, 121.760, 47.867, 98, 127.60, 180.94788,
    204.3833, 50.9415, 183.84, 131.293, 91.224, 88.90585, 174.9668, 0.0
]

atomic_numbers = [
    0, 1, 5, 6, 7, 8, 9, 15, 16, 17, 35, 26, 20, 53, 12, 29, 3, 19, 37, 55, 30,
    11, 13, 18, 33, 47, 79, 85, 4, 56, 83, 24, 27, 48, 87, 31, 32, 2, 72, 80,
    49, 77, 36, 25, 42, 10, 28, 41, 76, 46, 78, 82, 84, 44, 45, 75, 86, 88, 14,
    21, 34, 38, 50, 51, 22, 43, 52, 73, 81, 23, 74, 54, 40, 39, 71, 0
]

# copied from cpptraj
# Atom names corresponding to AtomicElementType.
atom_elements = [
    "??", "H", "B", "C", "N", "O", "F", "P", "S", "CL", "BR", "FE", "CA", "I",
    "MG", "CU", "LI", "K", "RB", "CS", "ZN", "NA", "AL", "AR", "AS", "AG",
    "AU", "AT", "BE", "BA", "BI", "CR", "CO", "CD", "FR", "GA", "GE", "HE",
    "HF", "HG", "IN", "IR", "KR", "MN", "MO", "NE", "NI", "NB", "OS", "PD",
    "PT", "PB", "PO", "RU", "RH", "RE", "RN", "RA", "SI", "SC", "SE", "SR",
    "SN", "SB", "TI", "TC", "TE", "TA", "TL", "V", "W", "XE", "ZR", "Y", "LU",
    "XP"
]

mass_atomic_number_dict = dict(zip(atomic_numbers, mass_arr))

mass_element_dict = dict(zip(atom_elements, mass_arr))

atomic_number_element_dict = dict(zip(atomic_numbers, atom_elements))

# Element is a dict with key=atomic_number
Element = dict(zip(atomic_numbers, zip(atom_elements, mass_arr)))
