#!/usr/bin/env python3
"""Parser minimale OLE2 + BIFF8 per .xls: estrae le celle numeriche/stringa di tutti i fogli.
Stdlib only. Gestisce MULRK, RK, NUMBER, LABELSST (SST con CONTINUE), LABEL, MULBLANK."""
import struct, sys, io

ENDOFCHAIN = 0xFFFFFFFE

def read_ole_stream(path, stream_name=b'Workbook'):
    with open(path, 'rb') as f:
        data = f.read()
    if data[:8] != b'\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1':
        raise ValueError('non e un OLE2')
    sector_shift = struct.unpack_from('<H', data, 0x1E)[0]
    mini_shift = struct.unpack_from('<H', data, 0x20)[0]
    sec_size = 1 << sector_shift
    n_fat = struct.unpack_from('<I', data, 0x2C)[0]
    dir_start = struct.unpack_from('<I', data, 0x30)[0]
    mini_cutoff = struct.unpack_from('<I', data, 0x38)[0]
    # DIFAT: primi 109 settori FAT
    difat = list(struct.unpack_from('<109I', data, 0x4C))
    # DIFAT dalla testata; se la continuazione e' difettosa, ricostruisci per contiguita'
    fat_sects = []
    for s in difat:
        if s < 0xFFFFFFFA:
            fat_sects.append(s)
    if len(fat_sects) < n_fat and fat_sects:
        last = fat_sects[-1]
        need = n_fat - len(fat_sects)
        max_sec = len(data) // sec_size
        for s in range(last + 1, last + 1 + need):
            if s <= max_sec - 1:
                fat_sects.append(s)
    fat_sects = fat_sects[:n_fat]
    fat = []
    for s in fat_sects:
        fat.extend(struct.unpack_from(f'<{sec_size // 4}I', data, s * sec_size))
    def chain(start):
        out, cur, seen = [], start, set()
        while cur not in (ENDOFCHAIN, 0xFFFFFFFF) and cur not in seen and cur < len(fat):
            seen.add(cur)
            out.append(cur)
            cur = fat[cur]
        return out
    def read_sectors(sectors):
        return b''.join(data[s * sec_size:(s + 1) * sec_size] for s in sectors)
    # directory
    dir_data = read_sectors(chain(dir_start))
    entries = []
    for off in range(0, len(dir_data), 128):
        e = dir_data[off:off + 128]
        if len(e) < 128:
            break
        nlen = struct.unpack_from('<H', e, 0x40)[0]
        etype = e[0x42]
        if etype == 0:
            continue
        name = e[:max(0, nlen - 2)].decode('utf-16-le', 'replace')
        start = struct.unpack_from('<I', e, 0x74)[0]
        size = struct.unpack_from('<I', e, 0x78)[0]
        entries.append((name, etype, start, size))
    target = None
    for name, etype, start, size in entries:
        if etype == 2 and (name.encode() == stream_name or name.encode() == b'Book'):
            target = (start, size)
            break
    if target is None:
        raise ValueError('stream Workbook non trovato: ' + repr(entries))
    start, size = target
    if size >= mini_cutoff:
        return read_sectors(chain(start))[:size]
    # mini stream (radice)
    root = None
    for name, etype, s0, sz in entries:
        if etype == 5:
            root = (s0, sz)
    rstart, rsize = root
    mini_data = read_sectors(chain(rstart))
    # miniFAT
    mf_start = struct.unpack_from('<I', data, 0x40)[0]
    if mf_start in (ENDOFCHAIN, 0xFFFFFFFF):
        return b''
    mfat = []
    for s in chain(mf_start):
        mfat.extend(struct.unpack_from(f'<{sec_size // 4}I', data, s * sec_size))
    msize = 1 << mini_shift
    out, cur, seen = [], start, set()
    while cur not in (ENDOFCHAIN, 0xFFFFFFFF) and cur not in seen and cur < len(mfat):
        seen.add(cur)
        out.append(mini_data[cur * msize:(cur + 1) * msize])
        cur = mfat[cur]
    return b''.join(out)[:size]


def rk_to_float(rk):
    if rk & 0x02:  # intero a 30 bit con segno
        val = rk >> 2
        if val >= 1 << 29:
            val -= 1 << 30
        val = float(val)
    else:
        # i 30 bit alti (bit 31..2) sono i bit 63..34 del double IEEE:
        # ricostruzione = shift << 32, nessuna correzione aggiuntiva
        bits = (rk & 0xFFFFFFFC) << 32
        val = struct.unpack('<d', struct.pack('<Q', bits))[0]
    if rk & 0x01:
        val /= 100.0
    return val

def records(stream):
    pos, n = 0, len(stream)
    while pos + 4 <= n:
        rtype, rlen = struct.unpack_from('<HH', stream, pos)
        pos += 4
        yield rtype, stream[pos:pos + rlen]
        pos += rlen

def parse_workbook(stream):
    # 1) SST e BOUNDSHEET dal substream globale
    sst = []
    boundsheets = []
    pending_sst = None
    pos = 0
    for rtype, rdata in records(stream):
        if rtype == 0x0085:  # BOUNDSHEET
            lbpl = struct.unpack_from('<I', rdata, 0)[0]
            cch = rdata[6]
            grbit = rdata[7]
            off = 8
            if grbit & 0x01:
                name = rdata[off:off + cch * 2].decode('utf-16-le', 'replace')
            else:
                name = rdata[off:off + cch].decode('latin-1', 'replace')
            boundsheets.append((name, lbpl))
        elif rtype == 0x00FC:  # SST inizio
            pending_sst = bytearray(struct.pack('<HH', 0x00FC, len(rdata))) + rdata
        elif rtype == 0x003C and pending_sst is not None:  # CONTINUE dell'SST
            pending_sst += b'\x3C' + struct.pack('<H', len(rdata)) + rdata
        elif rtype == 0x000A and pending_sst is not None:  # EOF del globals
            break
    if pending_sst is not None:
        # ricostruisci la lista dei record (per gestire i CONTINUE come stream)
        raw = bytes(pending_sst)
        sst = parse_sst_records(raw)
    return boundsheets, sst

def parse_sst_records(raw):
    # raw = concatenazione dei record SST+CONTINUE (con header)
    # semplificazione: unisci i payload, ma un record CONTINUE puo' ripetere
    # il flag di compressione per una stringa spezzata.
    chunks = []
    pos = 0
    while pos + 4 <= len(raw):
        rtype, rlen = struct.unpack_from('<HH', raw, pos)
        pos += 4
        chunks.append(raw[pos:pos + rlen])
        pos += rlen
    total, unique = struct.unpack_from('<II', chunks[0], 0)
    unique = total  # il campo unique in questi file e' corrotto; il numero di stringhe e' total
    # ricostruzione sequenziale con gestione della compressione per chunk
    out = []
    ci = 0
    off = 8
    cur = chunks[0]
    state_16 = False
    remaining_in_chunk = len(cur)
    def take(nbytes_or_chars, is_str_tail, grbit_holder):
        nonlocal ci, off, cur
        # gestisce anche il cambio di codifica all'inizio di un CONTINUE
        pass
    # implementazione classica con cursore virtuale:
    def vread(n):
        nonlocal ci, off, cur
        buf = b''
        while n > 0:
            avail = len(cur) - off
            if avail == 0:
                ci += 1
                if ci >= len(chunks):
                    raise EOFError
                cur = chunks[ci]
                off = 0
                avail = len(cur)
                if avail == 0:
                    continue
            take = min(n, avail)
            buf += cur[off:off + take]
            off += take
            n -= take
        return buf
    # il flag di compressione per stringhe spezzate si ribadisce all'inizio del CONTINUE:
    # gestito da vread se il chiamante controlla il boundary -> approssimazione:
    # leggiamo flag per stringa e cambiamo codifica se il cursore era resettato.
    try:
        for _ in range(unique):
            cch = struct.unpack('<H', vread(2))[0]
            flags = vread(1)[0]
            is16 = bool(flags & 0x01)
            if flags & 0x04:
                cfmt = struct.unpack('<H', vread(2))[0]
            else:
                cfmt = 0
            if flags & 0x08:
                vread(4)
            chars = []
            for _i in range(cch):
                # il flag si ribadisce a inizio CONTINUE: se off==0 siamo appena passati a un nuovo chunk
                if off == 0 and len(chars) > 0:
                    b2 = vread(1)
                    is16 = bool(b2[0] & 0x01)
                    chars.append(chr(struct.unpack('<H' if is16 else '<B', vread(2 if is16 else 1))[0]) if is16 else chr(b2[0]))
                    continue
                if is16:
                    chars.append(chr(struct.unpack('<H', vread(2))[0]))
                else:
                    chars.append(chr(vread(1)[0]))
            out.append(''.join(chars))
            for _f in range(cfmt * 4):
                vread(2)
    except (EOFError, Exception) as ex:
        sys.stderr.write(f'[sst] interrotto: {type(ex).__name__}: {ex}\n')
    return out

def parse_sheet(stream, offset, sst):
    cells = {}  # (row, col) -> value
    pos = offset
    n = len(stream)
    # salta fino al BOF del foglio
    while pos + 4 <= n:
        rtype, rlen = struct.unpack_from('<HH', stream, pos)
        pos += 4
        data = stream[pos:pos + rlen]
        pos += rlen
        if rtype == 0x0809:
            break
    while pos + 4 <= n:
        rtype, rlen = struct.unpack_from('<HH', stream, pos)
        pos += 4
        data = stream[pos:pos + rlen]
        pos += rlen
        if rtype == 0x000A:  # EOF
            break
        elif rtype == 0x00BD:  # MULRK: row(2) colFirst(2) [XF(2) RK(4)]* colLast(2)
            row, colF = struct.unpack_from('<HH', data, 0)
            ncols = (rlen - 6) // 6
            for j in range(ncols):
                off = 6 + j * 6   # salta i 2 byte di XF di ogni cella
                if off + 4 > len(data):
                    break  # record troncato al termine dello stream
                rk, = struct.unpack_from('<I', data, off)
                cells[(row, colF + j)] = rk_to_float(rk)
        elif rtype == 0x027E:  # RK
            row, col = struct.unpack_from('<HH', data, 0)
            rk, = struct.unpack_from('<I', data, 6)
            cells[(row, col)] = rk_to_float(rk)
        elif rtype == 0x0203:  # NUMBER
            row, col = struct.unpack_from('<HH', data, 0)
            cells[(row, col)] = struct.unpack_from('<d', data, 6)[0]
        elif rtype == 0x00FD:  # LABELSST
            row, col = struct.unpack_from('<HH', data, 0)
            isst, = struct.unpack_from('<I', data, 6)
            cells[(row, col)] = sst[isst] if isst < len(sst) else f'#{isst}'
        elif rtype == 0x0204:  # LABEL (inline)
            row, col = struct.unpack_from('<HH', data, 0)
            cch = struct.unpack_from('<H', data, 6)[0]
            g = data[8]
            if g & 0x01:
                cells[(row, col)] = data[9:9 + cch * 2].decode('utf-16-le', 'replace')
            else:
                cells[(row, col)] = data[9:9 + cch].decode('latin-1', 'replace')
    return cells

def xls_to_csv(path, out_prefix):
    stream = read_ole_stream(path)
    boundsheets, sst = parse_workbook(stream)
    sys.stderr.write(f'boundsheets: {boundsheets} | sst: {len(sst)} stringhe\n')
    if sst:
        sys.stderr.write('  sst[0:12]: ' + repr(sst[:12]) + '\n')
    allrows = []
    for idx, (name, lbpl) in enumerate(boundsheets):
        cells = parse_sheet(stream, lbpl, sst)
        maxr = max(r for r, c in cells) if cells else -1
        maxc = max(c for r, c in cells) if cells else -1
        rows = []
        for r in range(maxr + 1):
            rows.append([cells.get((r, c), '') for c in range(maxc + 1)])
        out = f'{out_prefix}_{idx}.csv'
        with open(out, 'w', newline='', encoding='utf-8') as f:
            import csv as csvmod
            csvmod.writer(f).writerows(rows)
        sys.stderr.write(f'foglio "{name}": {len(rows)} righe x {maxc + 1} colonne -> {out}\n')
        allrows.append((name, out))
    return allrows

if __name__ == '__main__':
    xls_to_csv(sys.argv[1], sys.argv[2])
