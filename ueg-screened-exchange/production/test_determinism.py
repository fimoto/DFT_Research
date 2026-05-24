"""Regression checks for refusing changed results; no expensive integrals."""
import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import runtime
if __name__=='__main__':runtime.enter_cli()
import common
import render
import compute


class DeterminismGuards(unittest.TestCase):
    def test_registered_source_revision_preserves_payload(self):
        rel=Path('map/rs_1_z_0.json')
        row=json.loads((common.ROOT/'production/raw'/rel).read_text())
        row['metadata']['sources']=common.source_hashes()
        common.check_reference(rel,row)
        row['settings']['m']+=1
        with self.assertRaisesRegex(RuntimeError,'differs'):
            common.check_reference(rel,row)

    def test_unknown_source_revision_is_rejected(self):
        row=json.loads((common.ROOT/'production/raw/map/rs_1_z_0.json').read_text())
        changed=common.source_hashes();changed['src/ueg_exchange/lindhard.py']='0'*64
        with patch.object(common,'source_hashes',return_value=changed):
            with self.assertRaisesRegex(RuntimeError,'without an audited'):
                common.check_source_compatibility(row['metadata']['sources'])

    def test_verified_force_preserves_historical_canonical_record(self):
        path=common.ROOT/'production/raw/map/rs_1_z_0.json'
        original=path.read_bytes();row=json.loads(original)
        with patch.object(sys,'argv',['compute.py','--dataset','map','--rs','1',
                                     '--zeta','0','--force']), \
             patch.object(compute,'direct',return_value=row),patch.object(compute,'save') as save:
            compute.main()
        save.assert_not_called()
        self.assertEqual(path.read_bytes(),original)

    def test_altered_replicate_is_rejected(self):
        p=next((common.ROOT/'production/raw/tableII').glob('*.json'))
        row=json.loads(p.read_text());common.check_reference(p.relative_to(common.ROOT/'production/raw'),row)
        row['term_replicates_ry'][0][0]+=1e-9
        with self.assertRaisesRegex(RuntimeError,'differs'):
            common.check_reference(p.relative_to(common.ROOT/'production/raw'),row)

    def test_changed_cached_record_is_rejected_without_write(self):
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'tableII/rs_5_z_0.json';path.parent.mkdir(parents=True)
            row=json.loads((common.ROOT/'production/raw/tableII/rs_5_z_0.json').read_text())
            row['mean_mry']+=0.01;path.write_bytes(common.encoded(row));before=path.read_bytes()
            r=subprocess.run([sys.executable,str(common.ROOT/'production/compute.py'),
                '--dataset','tableII','--rs','5','--output',td],text=True,capture_output=True)
            self.assertNotEqual(r.returncode,0)
            self.assertIn('differs from the released reference',r.stderr)
            self.assertEqual(before,path.read_bytes())

    def test_runtime_mismatch_is_rejected(self):
        expected=copy.deepcopy(runtime.contract());expected['python']='0.0.0'
        with patch.object(runtime,'contract',return_value=expected):
            with self.assertRaisesRegex(RuntimeError,'Runtime differs'):runtime.validate()

    def test_release_bypass_cannot_target_released_raw(self):
        r=subprocess.run([sys.executable,str(common.ROOT/'production/compute.py'),
            '--dataset','tableII','--prepare-release'],text=True,capture_output=True)
        self.assertNotEqual(r.returncode,0)
        self.assertIn('requires a separate --output',r.stderr)

    def test_changed_plot_does_not_replace_any_output(self):
        with tempfile.TemporaryDirectory() as td:
            target=Path(td)/'spin/figs/fig_map.pdf';target.parent.mkdir(parents=True)
            target.write_bytes(b'original output')
            def changed_plot(summary,stage):
                f=stage/'spin/figs/fig_map.pdf';f.parent.mkdir(parents=True)
                f.write_bytes(b'changed output')
            with patch.object(sys,'argv',['render.py','--output',td]), \
                 patch.object(render,'summarize',return_value={}), \
                 patch.object(render,'write_tables'),patch.object(render,'plots',side_effect=changed_plot):
                with self.assertRaisesRegex(RuntimeError,'no output was replaced'):render.main()
            self.assertEqual(target.read_bytes(),b'original output')


if __name__=='__main__':unittest.main()
