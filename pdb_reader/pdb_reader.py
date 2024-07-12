import pathlib
import re
import numpy as np
import os
import urllib3
# g_starts=split("0 6 12 16 17 21 22 26 30 38 46 54 60 72 76 78")
# g_lenghts=split("6 5 4 1 3 1 4 1 8 8 8 6 6 4 2 2")

PATH_SCRIPT = pathlib.Path(__file__).parent


def _interpret_line(line):
    try:
        atomType = line[12:14].strip().capitalize()
        if atomType == "F":
            atomType = "Fe"
    except IndexError:
        return None
    if atomType[0].isdigit():
        atomType = atomType[1]
    if len(atomType) > 1 and atomType[1].isdigit():
        atomType = atomType[0]

    try:
        name = line[13:15].strip()
    except Exception:
        name = atomType
    try:
        (x, y, z) = (float(line[30:38]), float(line[38:46]), float(line[46:55]))
    except IndexError:
        return None
    except ValueError:
        return None
    try:
        res = line[17:20]
    except Exception:
        res = ""
    try:
        resN = int(line[22:27])
    except Exception:
        resN = 0
    try:
        Bfactor = float(line[60:66])
    except Exception:
        Bfactor = 0.0
    return (atomType, name, res, resN, x, y, z, Bfactor)


def download_pdb(pdb_id, out_fname=None):
    out_fname = PATH_SCRIPT / "pdb" / f"{pdb_id}.pdb"
    if out_fname.is_file():
        print(f"File {out_fname} already available, reading from local folder")
        return PDB(out_fname)
    address = (
        "http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=%s"
        % pdb_id
    )
    p = urllib3.urlopen(address)
    lines = p.readlines()
    if out_fname is None:
        out_fname = PATH_SCRIPT / "pdb" / f"{pdb_id}.pdb"
    # try to save file
    try:
        folder = out_fname.parent
        if not folder.is_dir():
            os.makedirs(folder)
        with open(out_fname, "w") as f:
            f.write("".join(lines))
    except OSError as e:
        print(f"Could not save {out_fname}, error was {e}")
    pdb = pdb_reader.PDB(lines)
    return pdb


class PDB(object):
    def __init__(self, fname_or_lines):
        if isinstance(fname_or_lines,(str,pathlib.Path)):
            self.fname = os.path.basename(fname_or_lines)
        else:
            self.fname = "not_known"
        self.lines = self._read_lines_from_file(fname_or_lines)
    
    def _read_lines_from_file(self, fname_or_lines):
        if isinstance(fname_or_lines,(str,pathlib.Path)):
            with open(fname_or_lines, "r") as f:
                lines = f.readlines()
        else:
            lines = fname_or_lines
        info = []
        data = []
        reg_data = re.compile("^ATOM|^HETATM")
        for line in lines:
            if reg_data.match(line) is not None:
                temp = _interpret_line(line)
                if temp is None:
                    info.append(line)
                else:
                    data.append(temp)
            else:
                info.append(line)

        self.info = info
        try:
            dt = np.format_parser(
                "a3,a3,a3,i,f,f,f,f",
                ["atomType", "name", "res", "resNum", "x", "y", "z", "Bfactor"],
                [],
            ).dtype
            self.data = np.asarray(data, dtype=dt)
        except NameError:
            self.data = data

    def __str__(self):
        aType = np.unique(self.get_atom_types())
        N = {}
        for a in aType:
            N[a] = np.count_nonzero(self.get_atoms() == a)
        return "%s %s" % (self.fname, str(N))

    def __repr__(self):
        return self.__str__()

    def get_atom_types(self):
        return self.data["atomType"].astype(str)

    def get_atoms(self):
        return self.data["name"].astype(str)

    def get_residues(self):
        return self.data["res"].astype(str)

    def get_coords(self):
        return np.transpose(
            np.asarray((self.data["x"], self.data["y"], self.data["z"]))
        )

    def get_resnum(self):
        return self.data["resNum"]

    def write_pdb(self, fname, coords=None, idxFilter=None):
        name = self.get_atoms()
        res = self.get_residues()
        resNum = self.get_resnum()
        pos = self.get_coords()
        if coords is not None:
            pos = coords
        write_pdb(fname, name, res, pos, resNum=resNum, idxFilter=idxFilter)


def read_pdb(fname):
    return PDB(fname)


def write_pdb(fname, name, res, pos, resNum="auto", idxFilter=None):
    name = name[idxFilter]
    res = res[idxFilter]
    pos = pos[idxFilter]
    if resNum != "auto":
        resNum[idxFilter]
    pos = np.asarray(pos)
    nAtoms = pos.shape[0]
    chain = "A"
    f = open(fname, "w")
    _resNum = 1
    for i in range(0, nAtoms):
        if len(name[i]) < 4:
            s = "ATOM   %4d  %-3s" % (i + 1, name[i])
        else:
            s = "ATOM   %4d %-4s" % (i + 1, name[i])
        if resNum != "auto":
            s += " %s %s %3d    " % (res[i], chain, resNum[i])
        else:
            s += " %s %s %3d    " % (res[i], chain, _resNum)
        s += " %+7.3f %+7.3f %+7.3f\n" % (pos[i, 0], pos[i, 1], pos[i, 2])
        if (i > 0) & (res[i] != res[i - 1]):
            _resNum += 1
        f.write(s)
    f.close()
