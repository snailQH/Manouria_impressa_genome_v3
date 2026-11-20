#!/usr/bin/env python3
"""Correct the metadata of the Zenodo record minted from the GitHub release.

Zenodo's GitHub integration fills the record from the repository: the title gets the
"owner/repo:" prefix, the only creator is the GitHub account owner (so Jian Wang is
missing), there are no keywords and the only related identifier is the GitHub tree.
This script sets the intended metadata (the deposit procedure) on a
*published* record: open edit -> PUT metadata -> publish. Files are never touched.

Credentials: the token is read from the environment only (never a file, never argv):
    export ZENODO_TOKEN=...        # zenodo.org/account/settings/applications/tokens/new
                                   # scopes: deposit:write deposit:actions
    python3 zenodo_fix_metadata.py            # dry run: prints the diff
    python3 zenodo_fix_metadata.py --apply    # writes and re-publishes
    python3 zenodo_fix_metadata.py --record <record-id> --apply
"""
import json
import os
import sys
import urllib.error
import urllib.request

RECORD = '24387888'          # v3.0.1; concept 10.5281/zenodo.24056325
API = 'https://zenodo.org/api'

META = {
    'title': 'Manouria impressa genome assembly Mimp_v3 — pipeline and Data files 1–13',
    'upload_type': 'dataset',
    'version': 'v3.0',
    'language': 'eng',
    'creators': [
        {'name': 'Wang, Jian', 'affiliation': 'Ministry of Education Key Laboratory for Ecology of Tropical Islands, College of Life Sciences, Hainan Normal University, Haikou, China', 'orcid': '0000-0002-6164-4390'},
        {'name': 'Li, Qinghui', 'affiliation': 'Shanghai ImmunoBiotech Co., Ltd., Shanghai, China', 'orcid': '0000-0002-8494-0547'},
    ],
    'license': 'cc-by-4.0',      # data files; the code in the archive is MIT (stated in the description)
    'keywords': ['Manouria impressa', 'impressed tortoise', 'Testudines', 'genome assembly',
                 'haplotype-resolved', 'PacBio HiFi', 'Hi-C', 'chromosome-level', 'genome annotation',
                 'mitochondrial genome', 'Data Note'],
    'related_identifiers': [
        {'identifier': 'https://github.com/snailQH/Manouria_impressa_genome_v3/tree/v3.0', 'relation': 'isSupplementTo', 'scheme': 'url'},
        {'identifier': '10.1101/2025.03.10.642296', 'relation': 'isSupplementTo', 'scheme': 'doi'},
        {'identifier': 'PRJCA038476', 'relation': 'isDerivedFrom', 'scheme': 'other'},
        {'identifier': 'GWHHIVU00000000.2', 'relation': 'isDerivedFrom', 'scheme': 'other'},
        {'identifier': 'C_AA510199.1', 'relation': 'isDerivedFrom', 'scheme': 'other'},
    ],
    'notes': ('Data files 1–13 of the Data Note are inside the archive under data_files/. '
              'Code (pipeline, scripts, environments) is MIT-licensed (REPO_LICENSE); the data files are CC BY 4.0. '
              'Accessions: NGDC GWH GWHHIVU00000000.2 (Mimp_v3), GWHKKMH00000000.1 / GWHKKMI00000000.1 (haplotypes), '
              'GSA CRA024571 and CRA024586, GenBase C_AA510199.1 (mitogenome), BioProject PRJCA038476.'),
}


def call(method, url, token, payload=None):
    data = json.dumps(payload).encode if payload is not None else None
    req = urllib.request.Request(url, data=data, method=method,
                                 headers={'Content-Type': 'application/json',
                                          'Authorization': f'Bearer {token}'})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            body = r.read.decode
            return json.loads(body) if body.strip else {}
    except urllib.error.HTTPError as e:
        sys.exit(f'{method} {url} -> HTTP {e.code}: {e.read.decode[:400]}')


def main:
    apply_ = '--apply' in sys.argv
    rec = sys.argv[sys.argv.index('--record') + 1] if '--record' in sys.argv else RECORD
    token = os.environ.get('ZENODO_TOKEN')
    if not token:
        sys.exit('ZENODO_TOKEN is not set. export ZENODO_TOKEN=... (scopes deposit:write, deposit:actions), then re-run.')

    cur = call('GET', f'{API}/deposit/depositions/{rec}', token)
    old = cur.get('metadata', {})
    print(f'record {rec}  doi {cur.get("doi")}  concept {cur.get("conceptdoi")}  files {len(cur.get("files", []))}')
    for k, v in META.items:
        o = old.get(k)
        if k == 'creators':
            o = [c.get('name') for c in (o or [])]; v2 = [c['name'] for c in v]
        elif k == 'license':
            o = (o or {}).get('id') if isinstance(o, dict) else o; v2 = v
        else:
            v2 = v
        if o != v2:
            print(f'  {k}:\n    old: {str(o)[:150]}\n    new: {str(v2)[:150]}')
    if not apply_:
        print('\ndry run — nothing written. Re-run with --apply to update and re-publish.')
        return

    merged = dict(old); merged.update(META)
    merged.pop('doi', None); merged.pop('prereserve_doi', None)
    call('POST', f'{API}/deposit/depositions/{rec}/actions/edit', token)          # published -> editable
    call('PUT', f'{API}/deposit/depositions/{rec}', token, {'metadata': merged})
    out = call('POST', f'{API}/deposit/depositions/{rec}/actions/publish', token)
    print(f'published: {out.get("doi_url") or out.get("doi")}  (concept {out.get("conceptdoi")})')


if __name__ == '__main__':
    main
