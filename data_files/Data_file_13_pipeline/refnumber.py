"""Renumber bracketed citations by order of first appearance (Vancouver / BMC style) for every rendered output.
The Markdown sources keep the stable numbers of references_verified.md; only rendered DOCX/PDF/preprint get sequential numbers.
    rn = Renumberer(refs)               # refs: {old_number: text}
    body = rn.convert(body)             # rewrites [65,66,71] -> [1,2,3], [61–63] -> [4–6]; extends the map as new numbers appear
    tables = rn.convert(tables)         # same map, extended for citations that occur only in tables
    rn.reference_list()                 # [(new, text), ...] in final order"""
import re
CIT = re.compile(r'\[(\d+(?:\s*[,–-]\s*\d+)*)\]')
class Renumberer:
    def __init__(self, refs): self.refs = refs; self.map = {}
    def _new(self, old):
        if old not in self.map: self.map[old] = len(self.map) + 1
        return self.map[old]
    def _expand(self, group):
        out = []
        for part in re.split(r'\s*,\s*', group):
            a, *b = re.split(r'\s*[–-]\s*', part); a = int(a); b = int(b[-1]) if b else a
            out += list(range(a, b + 1))
        return out
    @staticmethod
    def _compress(nums):
        nums = sorted(set(nums)); parts = []; i = 0
        while i < len(nums):
            j = i
            while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1: j += 1
            parts.append(f'{nums[i]}–{nums[j]}' if j - i >= 2 else (f'{nums[i]},{nums[j]}' if j == i + 1 else str(nums[i]))); i = j + 1
        return '[' + ','.join(parts) + ']'
    def convert(self, text):
        def sub(m):
            olds = [o for o in self._expand(m.group(1)) if o in self.refs]
            if not olds: return m.group(0)           # not a citation (e.g. a bracketed number in a table) — leave untouched
            return self._compress([self._new(o) for o in olds])
        return CIT.sub(sub, text)
    def reference_list(self):
        return [(new, self.refs[old]) for old, new in sorted(self.map.items(), key=lambda kv: kv[1])]
