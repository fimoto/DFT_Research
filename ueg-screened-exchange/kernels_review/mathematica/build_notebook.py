"""Packaging-only helper; reviewers need Mathematica, not Python."""
from pathlib import Path
import json
root=Path(__file__).resolve().parent
source=(root/'kernels_checks.wl').read_text()
quote=lambda s:json.dumps(s,ensure_ascii=True)
# A string BoxData is a linear input representation. Keep the exact .wl
# program in one input cell so evaluation order is unambiguous.
cells=[
 'Cell['+quote('Symbolic and numerical checks: screened second-order exchange')+', "Title"]',
 'Cell['+quote('Unexecuted Wolfram Language companion. Start a fresh kernel and evaluate the input cell. This notebook embeds the full script and requires no Python or external code file. Timeouts and unresolved symbolic results are not PASS. See README.md for equation correspondence and analytic obligations.')+', "Text"]',
 'Cell[BoxData['+quote(source)+'], "Input"]',
 'Cell['+quote('After execution inspect JPAReview`results and JPAReview`summary. Numeric quadrature is not interval arithmetic; analytic limiting arguments and diagram identification remain in the manuscript.')+', "Text"]'
]
(root/'kernels_mathematica_review.nb').write_text('Notebook[{\n'+',\n'.join(cells)+'\n}, WindowTitle->"JPA Mathematica checks", StyleDefinitions->"Default.nb"]\n')
# Verify exact code embedding, not Wolfram syntax or kernel evaluation.
nb=(root/'kernels_mathematica_review.nb').read_text()
encoded=nb.split('Cell[BoxData[',1)[1].split('], "Input"]',1)[0]
assert json.loads(encoded)==source
print('PASS: notebook embeds exact Wolfram source. No kernel execution claimed.')
