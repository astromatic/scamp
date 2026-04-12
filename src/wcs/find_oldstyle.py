import re
import sys

def find_old_style_functions(filename):
    """Find K&R style function definitions"""
    with open(filename, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
    
    results = []
    i = 0
    while i < len(lines):
        line = lines[i]
        # Match function-like pattern: return_type name(params) with { on later line
        # or empty parens ()
        
        # Pattern 1: Function with empty parens: "type name()"
        match = re.match(r'^\s*([a-zA-Z_][a-zA-Z0-9_*\s]+?)\s+([a-zA-Z_][a-zA-Z0-9_]*)\s*\(\s*\)\s*$', line)
        if match and i + 1 < len(lines) and '{' in lines[i+1]:
            results.append({
                'type': 'empty_parens',
                'line': i + 1,
                'definition': line.rstrip(),
                'filename': filename,
                'return_type': match.group(1).strip(),
                'func_name': match.group(2),
            })
            i += 1
            continue
        
        # Pattern 2: K&R style - function name with params but no types: func(a, b)
        # followed by type declarations
        match = re.match(r'^\s*([a-zA-Z_][a-zA-Z0-9_*\s]+?)\s+([a-zA-Z_][a-zA-Z0-9_]*)\s*\(\s*([a-zA-Z_][a-zA-Z0-9_,\s]*)\s*\)\s*$', line)
        if match:
            return_type = match.group(1).strip()
            func_name = match.group(2)
            params_str = match.group(3).strip()
            
            # Check if next line(s) have type declarations
            j = i + 1
            type_decls = []
            while j < len(lines) and not '{' in lines[j]:
                decl_line = lines[j].strip()
                if decl_line and not decl_line.startswith('/*') and not decl_line.startswith('*'):
                    # Check if it looks like a type declaration
                    if re.match(r'^(const\s+)?[a-zA-Z_][a-zA-Z0-9_\s*]+\s*[a-zA-Z_][a-zA-Z0-9_]*\s*(=.*)?;?\s*$', decl_line):
                        type_decls.append((j + 1, decl_line))
                        j += 1
                    elif decl_line == '':
                        j += 1
                    else:
                        break
                else:
                    j += 1
            
            if type_decls:
                results.append({
                    'type': 'k_and_r',
                    'line': i + 1,
                    'definition': line.rstrip(),
                    'filename': filename,
                    'return_type': return_type,
                    'func_name': func_name,
                    'params': params_str,
                    'type_decls': type_decls,
                })
                i = j
                continue
        
        i += 1
    
    return results

# Check all files
files = [
    'wcstrig.c', 'poly.c', 'wcs.c', 'proj.c', 'tnx.c', 'cel.c', 'lin.c', 'sph.c'
]

for f in files:
    try:
        results = find_old_style_functions(f)
        for r in results:
            print(f"\n{f}:{r['line']}: {r['definition']}")
            if r['type'] == 'k_and_r' and 'type_decls' in r:
                for line_no, decl in r['type_decls']:
                    print(f"{f}:{line_no}: {decl}")
    except FileNotFoundError:
        pass
