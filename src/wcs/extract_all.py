import re

def extract_function_def(filename, line_no):
    """Extract function definition starting at line_no (1-indexed)"""
    with open(filename, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
    
    idx = line_no - 1
    sig_line = lines[idx].rstrip()
    
    # Look for parameter declarations in next lines
    params = []
    j = idx + 1
    while j < len(lines):
        l = lines[j].strip()
        if l == '' or l.startswith('/*') or l.startswith('*'):
            j += 1
            continue
        if l.startswith('{'):
            break
        # Parameter declaration
        params.append(lines[j].rstrip())
        j += 1
    
    return sig_line, params

files_funcs = {
    'wcstrig.c': [84, 107, 130, 151, 169, 187, 205],
    'proj.c': [275, 305, 338, 408, 456, 534, 585, 681, 701, 728, 761, 782, 804, 841, 1077, 1174, 1335, 1434, 1505, 1585, 1683, 1797, 1899, 2026, 2149, 2248, 2393, 2468, 2560, 2657, 2787, 3087, 3458, 3658],
    'cel.c': [311, 568, 597],
    'lin.c': [201, 238, 269, 301],
    'sph.c': [143, 200],
    'wcs.c': [456],
}

for fname in ['wcstrig.c', 'proj.c', 'cel.c', 'lin.c', 'sph.c', 'wcs.c']:
    if fname not in files_funcs:
        continue
    for line_no in sorted(files_funcs[fname]):
        sig, params = extract_function_def(fname, line_no)
        print(f"\n{fname}:{line_no}: {sig}")
        for p in params:
            print(f"  {p}")
