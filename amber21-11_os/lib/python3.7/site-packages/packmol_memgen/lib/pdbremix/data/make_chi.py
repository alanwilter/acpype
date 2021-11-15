import pprint

chi_topology_str = """
: CYS
 N CA CB SG 
: ASP
 N CA CB CG 
 CA CB CG OD1 
: GLU
 N CA CB CG 
 CA CB CG CD 
 CB CG CD OE1 
: PHE
 N CA CB CG 
 CA CB CG CD1 
: HIS
 N CA CB CG 
 CA CB CG ND1 
: ILE
 N CA CB CG1 
 CA CB CG1 CD1 
: LYS
 N CA CB CG 
 CA CB CG CD 
 CB CG CD CE 
 CG CD CE NZ 
: LYN
 N CA CB CG 
 CA CB CG CD 
 CB CG CD CE 
 CG CD CE NZ 
: LYP
 N CA CB CG 
 CA CB CG CD 
 CB CG CD CE 
 CG CD CE NZ 
: LEU
 N CA CB CG 
 CA CB CG CD1 
: MET
 N CA CB CG 
 CA CB CG SD 
 CB CG SD CE 
: ASN
 N CA CB CG 
 CA CB CG OD1 
: PRO
 N CA CB CG 
 CA CB CG CD 
 CB CG CD N 
 CG CD N CA 
: GLN
 N CA CB CG 
 CA CB CG CD 
 CB CG CD OE1 
: ARG
 N CA CB CG 
 CA CB CG CD 
 CB CG CD NE 
 CG CD NE CZ 
: SER
 N CA CB OG 
: THR
 N CA CB OG1 
: VAL
 N CA CB CG1 
: TRP
 N CA CB CG 
 CA CB CG CD1 
: TYR
 N CA CB CG 
 CA CB CG CD1 

: PHD
 N CA CB CG 
 CA CB CG OD1 
"""

chi_topology = {}

curr_res_type = ""
for line in chi_topology_str.splitlines():
  word_list = line.split()
  if word_list:
    if ":" in word_list[0]:
      curr_res_type = word_list[1]
      chi_topology[curr_res_type] = []
    elif curr_res_type:
      chi_topology[curr_res_type].append(word_list[0:4])

pprint.pprint(chi_topology, indent=2)

