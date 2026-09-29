import zipfile, sys, csv
from xml.etree.ElementTree import iterparse

NS = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'

def col_index(ref):
    col = 0
    for ch in ref:
        if ch.isalpha():
            col = col * 26 + (ord(ch.upper()) - 64)
        else:
            break
    return col - 1

def parse_xlsx(path, out_csv):
    z = zipfile.ZipFile(path)
    shared = []
    if 'xl/sharedStrings.xml' in z.namelist():
        for ev, el in iterparse(z.open('xl/sharedStrings.xml'), events=('end',)):
            if el.tag == NS + 'si':
                shared.append(''.join(t.text or '' for t in el.iter(NS + 't')))
                el.clear()
    sheet = [n for n in z.namelist() if n.startswith('xl/worksheets/sheet')][0]
    rows_out = []
    for ev, row in iterparse(z.open(sheet), events=('end',)):
        if row.tag != NS + 'row':
            continue
        cells = {}
        for c in row.iter(NS + 'c'):
            ref = c.get('r'); idx = col_index(ref)
            v = c.find(NS + 'v'); t = c.get('t')
            if t == 's':
                val = shared[int(v.text)] if v is not None else ''
            elif t == 'inlineStr':
                is_el = c.find(NS + 'is')
                val = ''.join(t2.text or '' for t2 in is_el.iter(NS + 't')) if is_el is not None else ''
            else:
                val = v.text if v is not None else ''
            cells[idx] = val
        if cells:
            n = max(cells) + 1
            rows_out.append([cells.get(i, '') for i in range(n)])
        row.clear()
    with open(out_csv, 'w', newline='') as f:
        csv.writer(f).writerows(rows_out)
    return len(rows_out)

if __name__ == '__main__':
    print(parse_xlsx(sys.argv[1], sys.argv[2]))
