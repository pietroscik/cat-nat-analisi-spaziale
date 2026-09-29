#!/usr/bin/env python3
"""Recovery generalizzato dei .xls SA_*.xls INGV: DIFAT sfasato di 1 e continuazione azzerata.
Ricostruzione: blocchi FAT ai settori [difat0+1 .. difat0+n_fat-1], stream Workbook con
off-by-one (il primo settore e' l'header OLE). Estrae tutti i fogli in CSV."""
import struct, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parse_biff import parse_workbook, parse_sheet

ENDOFCHAIN = 0xFFFFFFFE
FREE = 0xFFFFFFFF

def recover(path):
    data = open(path, 'rb').read()
    nsec = len(data) // 512
    n_fat = struct.unpack_from('<I', data, 0x2C)[0]
    dir_start = struct.unpack_from('<I', data, 0x30)[0]
    difat = struct.unpack_from('<109I', data, 0x4C)
    d0 = next((x for x in difat if x < 0xFFFFFFFA), None)
    if d0 is None:
        raise ValueError('difat vuoto')
    base = d0 + 1                      # sfasamento di uno
    nblocks = n_fat - 1
    fat = []
    for j in range(nblocks):
        s = base + j
        if s >= nsec:
            break
        fat.extend(struct.unpack_from('<128I', data, s * 512))
    # directory: settori letti direttamente (il FAT ricostruito non copre il fondo)
    dir_secs = [s for s in (dir_start, dir_start + 1, dir_start + 2) if s < nsec]
    dirdata = b''.join(data[s * 512:(s + 1) * 512] for s in dir_secs)
    wb = None
    for off in range(0, len(dirdata), 128):
        e = dirdata[off:off + 128]
        if len(e) < 128:
            break
        etype = e[0x42]
        if etype != 2:
            continue
        nlen = struct.unpack_from('<H', e, 0x40)[0]
        name = e[:max(0, nlen - 2)].decode('utf-16-le', 'replace')
        s0 = struct.unpack_from('<I', e, 0x74)[0]
        sz = struct.unpack_from('<I', e, 0x78)[0]
        if name not in ('Workbook', 'Book'):
            continue
        secs, cur, seen = [], s0, set()
        while cur not in (ENDOFCHAIN, FREE) and cur not in seen and cur < len(fat):
            seen.add(cur)
            secs.append(cur)
            cur = fat[cur]
        stream = b''.join(data[s * 512:(s + 1) * 512] for s in secs)
        if stream[:8] == b'\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1':
            stream = stream[512:]   # scarta il settore-header OLE duplicato (off-by-one)
        wb = stream                # catena completa del Workbook
        break
    if wb is None or wb[:2] != b'\x09\x08':
        raise ValueError(f'BIF non valido in {path}: prime {wb[:8].hex() if wb else None}')
    return wb

def main():
    src, out_prefix = sys.argv[1], sys.argv[2]
    wb = recover(src)
    print(f'{os.path.basename(src)}: stream BIFF {len(wb)} byte')
    boundsheets, sst = parse_workbook(wb)
    print('  fogli:', [n for n, _ in boundsheets], '| sst:', len(sst))
    import csv
    for name, lbpl in boundsheets:
        cells = parse_sheet(wb, lbpl, sst)
        maxr = max(r for r, c in cells)
        maxc = max(c for r, c in cells)
        rows = [[cells.get((r, c), '') for c in range(maxc + 1)] for r in range(maxr + 1)]
        safe = name.replace('/', '_')
        out = f'{out_prefix}_{safe}.csv'
        with open(out, 'w', newline='', encoding='utf-8') as f:
            csv.writer(f).writerows(rows)
        print(f'  -> {out}: {len(rows)} righe x {maxc + 1} colonne')

if __name__ == '__main__':
    main()
