from pathlib import Path
import hashlib, json, html
from PIL import Image, ImageDraw, ImageFont

root = Path(__file__).resolve().parents[1]
out = root / 'docs'
out.mkdir(parents=True, exist_ok=True)
W,H,S = 1200,960,3
im = Image.new('RGB',(W*S,H*S),'white'); d = ImageDraw.Draw(im)
svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="12in" height="9.6in" viewBox="0 0 {W} {H}">', '<rect width="1200" height="960" fill="white"/>']
fontdir=Path('C:/Windows/Fonts')
def txt(x,y,t,size=18,bold=False,color='#172b3a'):
    f=ImageFont.truetype(str(fontdir/('arialbd.ttf' if bold else 'arial.ttf')),size*S)
    d.text((x*S,y*S),t,font=f,fill=color)
    svg.append(f'<text x="{x}" y="{y+size*.92}" font-family="Arial, sans-serif" font-size="{size}" font-weight="{700 if bold else 400}" fill="{color}">{html.escape(t)}</text>')
def rect(x,y,w,h,fill,stroke='#b6c6cf'):
    d.rounded_rectangle((x*S,y*S,(x+w)*S,(y+h)*S),radius=10*S,fill=fill,outline=stroke,width=2*S)
    svg.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="10" fill="{fill}" stroke="{stroke}" stroke-width="2"/>')
def box(x,y,w,h,title,lines,fill='#eff6f8'):
    rect(x,y,w,h,fill); txt(x+18,y+13,title,20,True)
    for i,line in enumerate(lines):txt(x+18,y+45+i*25,line,17)
def arrow(x1,y1,x2,y2):
    d.line((x1*S,y1*S,x2*S,y2*S),fill='#466577',width=3*S)
    if y2>y1: pts=[(x2,y2),(x2-6,y2-9),(x2+6,y2-9)]
    else: pts=[(x2,y2),(x2-9,y2-6),(x2-9,y2+6)]
    d.polygon([(x*S,y*S) for x,y in pts],fill='#466577')
    svg.append(f'<path d="M{x1},{y1} L{x2},{y2}" stroke="#466577" stroke-width="3"/><polygon points="'+ ' '.join(f'{x},{y}' for x,y in pts)+'" fill="#466577"/>')

txt(40,25,'Repaired evaluation workflow',32,True)
txt(40,70,'Implemented software core; biological evaluation has not been run',21,False,'#805019')
box(40,115,1120,107,'Required before a biological run',[
    'Review labels, sequence integrity, terminal completeness and relatedness; declare the target population.',
    'Supply aligned sequences, binary labels and reviewed group IDs. These reviews remain pending.'
], '#fff7e9')
arrow(330,222,330,251)
box(40,252,560,104,'Validate input and define outer folds',[
    'Reject invalid residues and sequences shorter than window W.',
    'Use stratified group folds; fail if a fold lacks either class.'
])
arrow(320,356,320,384)
box(40,385,560,101,'Within each outer training partition',[
    'Create group-preserving fitting / calibration partitions.',
    'Keep the outer test partition outside all fitting steps.'
])
arrow(320,486,320,514)
box(40,515,560,128,'Tune using fitting records only',[
    'Inner group folds: fit positive / negative sequence profiles',
    'inside the pipeline, then fit and tune an RBF SVM.',
    'Refit the selected pipeline on its fitting partition.'
])
arrow(320,643,320,671)
box(40,672,560,103,'Fit the calibration component',[
    'Fit a sigmoid on the held-out calibration records.',
    'Apply pipeline and sigmoid to the untouched outer test.'
])
arrow(600,723,645,723)
box(646,643,514,132,'Combine and save outer-test outputs',[
    'Average calibrated component scores for each record.',
    'Classify at the fixed score threshold of 0.5.',
    'Return classes, scores, fold IDs and fitting traces.'
])
box(646,252,514,138,'Partition boundary',[
    'Outer test labels are not used to fit profiles,',
    'select hyperparameters or fit the sigmoid.',
    'Group separation applies at every split level.'
], '#eaf4ee')
box(646,421,514,163,'Still required after a biological run',[
    'Check matched-window performance and uncertainty.',
    'Run justified permutation controls; reconcile counts.',
    'Assess calibration and limits of generalization.',
    'Generate result figures from saved record outputs.'
], '#fff7e9')
rect(40,818,1120,94,'#f3f5f7')
txt(58,834,'Scope: this schematic describes src/repair_model.py, not the deployed web predictor.',17,True)
txt(58,864,'Synthetic software tests do not establish biological accuracy or calibrated localization probabilities.',17)
svg.append('</svg>')
(out/'Workflow_Figure.svg').write_text('\n'.join(svg),encoding='utf-8')
im.save(out/'Workflow_Figure_300dpi.png',dpi=(300,300))
caption='''**Workflow Figure. Repaired evaluation design (biological run pending).** The schematic describes the supplied software core, not completed biological analyses or the deployed website. Input and provenance review precede use. Group boundaries are preserved in outer evaluation, calibration and inner tuning partitions. Supervised profiles are fitted within the classifier pipeline on training records at each tuning step. Each calibration component tunes on its fitting partition and fits a sigmoid using held-out calibration records; its outputs are averaged for the untouched outer test partition. The core returns predicted classes, continuous scores, fold identifiers and fitting traces. Biological performance, uncertainty, permutation controls and calibration assessment remain unexecuted. Source: `src/repair_model.py`; implementation inspected during this audit. This methods schematic does not replace the missing data-dependent Figures 1–3.
'''
(out/'Workflow_Figure_caption.md').write_text(caption,encoding='utf-8')
manifest={'artifact_type':'methods_schematic','status':'EXECUTED_FIGURE_RENDER_ONLY','biological_analysis':'NOT_RUN','source':'src/repair_model.py','source_sha256':hashlib.sha256((root/'src/repair_model.py').read_bytes()).hexdigest(),'png_pixels':[W*S,H*S],'png_dpi':300,'claims':['validate_sequences rejects invalid residues and short sequences','checked_splits preserves groups and checks classes','nested_cv uses separate fitting calibration and outer test indices','GridSearchCV fits BiProfileEncoder inside Pipeline','sigmoid fitting uses calibration records','outer test component scores are averaged and thresholded at 0.5'],'figure_not_a_performance_result':True}
(out/'Workflow_Figure_provenance.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print('Rendered methods schematic; no biological analysis performed.')
