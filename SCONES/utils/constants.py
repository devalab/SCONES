AA_ID_DICT = {
    "A": 0,
    "C": 1,
    "D": 2,
    "E": 3,
    "F": 4,
    "G": 5,
    "H": 6,
    "I": 7,
    "K": 8,
    "L": 9,
    "M": 10,
    "N": 11,
    "P": 12,
    "Q": 13,
    "R": 14,
    "S": 15,
    "T": 16,
    "V": 17,
    "W": 18,
    "Y": 19,
}

# https://github.com/minghuilab/PremPS/blame/main/v1.0.0/PremPS.py#L58-L61
AA_SASA = {'A':118.1,'R':256.0,'N':165.5,'D':158.7,'C':146.1,'Q':193.2,
           'E':186.2,'G':88.1,'H':202.5,'I':181.0,'L':193.1,'K':225.8,
           'M':203.4,'F':222.8,'P':146.8,'S':129.8,'T':152.5,'W':266.3,
           'Y':236.8,'V':164.5,'X':88.1}

AA_SC_PKA = {
    "A" : 0,
    "C" : 8.14,
    "D" : 3.71,
    "E" : 4.15,
    "F" : 0,
    "G" : 0,
    "H" : 6.04,
    "I" : 0,
    "K" : 10.67,
    "L" : 0,
    "M" : 0,
    "N" : 0,
    "P" : 0,
    "Q" : 0,
    "R" : 12.10, 
    "S" : 0,
    "T" : 0,
    "V" : 0,
    "W" : 0,
    "Y" : 10.10,
}

ID_AA_DICT = {v: k for k, v in AA_ID_DICT.items()}

SS_ID_DICT = {
    "H" : 0, # alpha helix
    "B" : 1, # isolated beta-bridge residue
    "E" : 2, # strand
    "G" : 3, # 3-10 helix
    "I" : 4, # pi helix
    "T" : 5, # turn
    "S" : 6, # bend
	"-" : 7, # coil
}

SS_NAME_DICT = {
    "H" : "alpha helix",
    "B" : "isolated beta-bridge residue",
    "E" : "strand",
    "G" : "3-10 helix",
    "I" : "pi helix",
    "T" : "turn",
    "S" : "bend",
	"-" : "coil"
}

AA_PROPERTY_DICT = {
    # [hydrophobic, polar, charge, acidic, basic, neutral]
    "A": [1,    0,  0],
    "C": [0.5,  0,  0],
    "D": [0,    0,  -1],
    "E": [0,    0,  -1],
    "F": [1,    0,  0],
    "G": [1,    0,  0],
    "H": [0,    0,  1],
    "I": [1,    0,  0],
    "K": [0,    0,  1],
    "L": [1,    0,  0],
    "M": [1,    0,  0],
    "N": [0,    1,  0],
    "P": [1,    0,  0],
    "Q": [0,    1,  0],
    "R": [0,    0,  1],
    "S": [0,    1,  0],
    "T": [0,    1,  0],
    "V": [1,    0,  0],
    "W": [1,    0,  0],
    "Y": [1,    0,  0]
}

AA_NAME_DICT = {
    "A": "Alanine",
    "C": "Cysteine",
    "D": "Aspartic Acid",
    "E": "Glutamic Acid",
    "F": "Phenylalanine",
    "G": "Glycine",
    "H": "Histidine",
    "I": "Isoleucine",
    "K": "Lysine",
    "L": "Leucine",
    "M": "Methionine",
    "N": "Asparagine",
    "P": "Proline",
    "Q": "Glutamine",
    "R": "Arginine",
    "S": "Serine",
    "T": "Threonine",
    "V": "Valine",
    "W": "Tryptophan",
    "Y": "Tyrosine",
}

AA_COLOR_DICT = {
    "A": "#BBBBBB",
    "C": "#D0FF00",
    "D": "#FFBF00",
    "E": "#FFBF00",
    "F": "#BBBBBB",
    "G": "#BBBBBB",
    "H": "#4A4AB0",
    "I": "#BBBBBB",
    "K": "#4A4AB0",
    "L": "#BBBBBB",
    "M": "#D0FF00",
    "N": "#EEEEEE",
    "P": "#BBBBBB",
    "Q": "#EEEEEE",
    "R": "#4A4AB0",
    "S": "#96E617",
    "T": "#96E617",
    "V": "#BBBBBB",
    "W": "#BBBBBB",
    "Y": "#BBBBBB",
}