import re

def scan_file(filename):
    with open(filename, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
    
    results = []
    i = 0
    while i < len(lines):
        line = lines[i]
        
        # Skip comments and empty lines
        stripped = line.strip()
        if not stripped or stripped.startswith('/*') or stripped.startswith('*') or stripped.startswith('//'):
            i += 1
            continue
        
        # Pattern 1: empty parens like "void foo()"
        match = re.match(r'^\s*([a-zA-Z_][a-zA-Z0-9_*\s]+?)\s+([a-zA-Z_][a-zA-Z0-9_]*)\s*\(\s*\)\s*$', line)
        if match and i + 1 < len(lines):
            next_line = lines[i+1].strip()
            if next_line.startswith('{'):
                sig = line.rstrip()
                results.append(('empty_parens', i+1, sig))
                i += 1
                continue
        
        # Pattern 2: K&R with params on declaration line: func(a, b)
        match = re.match(r'^\s*([a-zA-Z_][a-zA-Z0-9_*\s]+?)\s+([a-zA-Z_][a-zA-Z0-9_]*)\s*\(\s*([a-zA-Z_][a-zA-Z0-9_,\s]*)\s*\)\s*$', line)
        if match:
            return_type = match.group(1).strip()
            func_name = match.group(2)
            params = match.group(3).strip()
            
            # Look ahead for parameter type declarations
            j = i + 1
            decls = []
            blank_lines = 0
            while j < len(lines):
                next_line = lines[j]
                next_stripped = next_line.strip()
                
                if next_stripped == '':
                    blank_lines += 1
                    j += 1
                    continue
                    
                if next_stripped.startswith('{'):
                    # Found opening brace - this is K&R style
                    if decls or (blank_lines > 0 and j - i > 1):
                        sig = line.rstrip()
                        results.append(('k_and_r', i+1, sig, decls))
                    break
                
                if next_stripped.startswith('/*') or next_stripped.startswith('*'):
                    j += 1
                    continue
                
                # Check if it looks like a parameter type declaration
                if re.match(r'^(const\s+)?(struct\s+)?[a-zA-Z_][a-zA-Z0-9_\s*]+[a-zA-Z_][a-zA-Z0-9_]*\s*(=.*)?;?\s*$', next_stripped):
                    decls.append((j+1, next_line.rstrip()))
                    j += 1
                    continue
                else:
                    break
            
            i = j + 1 if j < len(lines) else len(lines)
            continue
        
        i += 1
    
    return results

files = ['wcstrig.c', 'poly.c', 'wcs.c', 'proj.c', 'tnx.c', 'cel.c', 'lin.c', 'sph.c']

all_results = []
for fname in files:
    try:
        results = scan_file(fname)
        for r in results:
            all_results.append((fname,) + r)
    except:
        pass

# Sort by filename and line number
all_results.sort(key=lambda x: (x[0], x[2]))

for fname, kind, line_no, sig, *rest in all_results:
    print(f"{fname}:{line_no}: {sig}")
    if rest and rest[0]:
        for decl_line, decl_text in rest[0]:
            print(f"  {fname}:{decl_line}: {decl_text}")
