"""Build the source-grounded patient and segmentation-record flow diagram."""
from pathlib import Path
import argparse, collections, csv, datetime, hashlib, json, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from PIL import Image, ImageChops

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def rows(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def write_csv(path, data):
    with Path(path).open('w', encoding='utf-8-sig', newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)

def main(argv=None):
    parser=argparse.ArgumentParser()
    parser.add_argument('--run-root',required=True,type=Path)
    parser.add_argument('--output-root',required=True,type=Path)
    args=parser.parse_args(argv)
    run=args.run_root.resolve();out=args.output_root.resolve();out.mkdir(parents=True,exist_ok=True)
    rev=run/'revision_R1_20261014'
    inputs={
        'candidate_disposition':run/'FINAL_CANDIDATE_DISPOSITION.csv',
        'image_roster':run/'FINAL_IMAGE_AND_ALIAS_ROSTER.csv',
        'representative_roster':run/'FINAL_PATIENT_AND_REPRESENTATIVE_ROSTER.csv',
        'code_mapping':rev/'01_data_audit/ALL_PAIR_REVIEW_CODE_MAPPING_KEEP_SEPARATE.csv',
        'source_eligibility':rev/'01_data_audit/phase_analysis_eligibility.csv',
        'all_source_pairs':rev/'01_data_audit/phase_pairing_audit.csv',
        'source_clinical':rev/'01_data_audit/source_clinical_233.csv',
        'actual_summary':rev/'03_model_development/final_analysis_summary.json',
        'actual_matrix':rev/'03_model_development/analysis_matrix.csv',
    }
    before={str(p):sha(p) for p in inputs.values()}
    candidates=rows(inputs['candidate_disposition']);images=rows(inputs['image_roster'])
    pairmap=rows(inputs['code_mapping']);eligibility=rows(inputs['source_eligibility'])
    source_pairs=rows(inputs['all_source_pairs']);clinical=rows(inputs['source_clinical'])
    summary=json.loads(inputs['actual_summary'].read_text(encoding='utf-8-sig'))
    code={r['patient_id']:r['review_code'].split('L')[0] for r in pairmap}
    included={r['patient_id'] for r in candidates if r['eligibility']=='include'}
    excluded={r['patient_id'] for r in candidates if r['eligibility']=='exclude'}
    active=[r for r in images if r['active_analysis_record'].lower()=='true']
    inactive=[r for r in images if r['active_analysis_record'].lower()!='true']
    unused_excluded=[r for r in inactive if r['patient_id'] in excluded]
    unused_retained=[r for r in inactive if r['patient_id'] in included]
    recovered=[r for r in eligibility if r['in_old_195'].lower()=='false']
    arterial_pairs=[r for r in source_pairs if int(float(r['ct_phase']))==1]
    expected_delay={'36'};expected_unknown={'83','94'}
    no_roi=excluded-expected_delay-expected_unknown
    assert len(candidates)==126 and len(included)==118 and len(excluded)==8
    assert no_roi=={'338','360','399','92','96'}
    assert len(images)==230 and len(active)==202 and len(inactive)==28
    assert len(unused_excluded)==16 and len(unused_retained)==12
    assert {r['patient_id'] for r in active}==included
    assert len(rows(inputs['representative_roster']))==len(rows(inputs['actual_matrix']))==summary['n']==118
    assert summary['events']==83 and summary['censored']==35 and summary['lesions']==202
    assert len(clinical)==len(eligibility)==233 and len(source_pairs)==377
    assert len(arterial_pairs)==230 and len({r['patient_id'] for r in arterial_pairs})==126
    other_pairs=[r for r in source_pairs if int(float(r['ct_phase']))!=1]
    assert len({r['patient_id'] for r in other_pairs})==107
    assert not ({r['patient_id'] for r in arterial_pairs} & {r['patient_id'] for r in other_pairs})
    assert len(recovered)==38 and sum(int(r['source_masks']) for r in recovered)==59
    assert sum(r['in_old_195'].lower()=='true' for r in eligibility)==195
    assert sum(int(r['source_masks']) for r in eligibility if r['in_old_195'].lower()=='true')==318
    assert len(source_pairs)-len(arterial_pairs)==147 and len(clinical)-len(candidates)==107
    excluded_rows=[]
    for pid in sorted(excluded,key=lambda p:code[p]):
        reason='Delayed phase' if pid in expected_delay else 'Indeterminate phase' if pid in expected_unknown else 'No valid ROI'
        excluded_rows.append({'review_code':code[pid],'patient_id':pid,'patient_exclusion_reason':reason,
                              'unused_source_records':sum(r['patient_id']==pid for r in unused_excluded)})
    write_csv(out/'Figure_1_excluded_patient_source.csv',excluded_rows)
    write_csv(out/'Figure_1_unused_records_source.csv',[
        {'sample_id':r['sample_id'],'patient_id':r['patient_id'],'patient_retained':r['patient_id'] in included,
         'record_disposition':'Unused source segmentation record'} for r in inactive])
    write_csv(out/'Figure_1_recovered_source.csv',[
        {'patient_id':r['patient_id'],'recovered_source_records':int(r['source_masks'])} for r in recovered])
    counts=[('source_patients',233),('source_records',377),('earlier_patients',195),('earlier_pairs',318),
            ('recovered_patients',38),('recovered_pairs',59),('arterial_candidates',126),('arterial_source_records',230),
            ('other_source_phase_patients',107),('other_source_phase_records',147),('included_patients',118),
            ('excluded_patients',8),('delayed_phase_patients',1),('indeterminate_phase_patients',2),('no_valid_ROI_patients',5),
            ('active_records',202),('unused_records',28),('unused_records_excluded_patients',16),
            ('unused_records_retained_patients',12),('actual_deaths',83),('actual_censored',35)]
    write_csv(out/'Figure_1_counts_source.csv',[{'quantity':k,'n':v} for k,v in counts])
    plt.rcParams.update({'font.family':'Arial','font.size':8.5,'svg.fonttype':'none','pdf.fonttype':42,
                         'axes.linewidth':1.1,'text.color':'#172B35'})
    fig,ax=plt.subplots(figsize=(7.09,4.8))
    fig.subplots_adjust(left=.015,right=.985,bottom=.015,top=.985)
    ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
    boxdata=[]
    def box(x,y,w,h,text,fill='#F0F5F7',fs=8.5,bold=False):
        ax.add_patch(Rectangle((x,y),w,h,facecolor=fill,edgecolor='#617D8B',linewidth=1.1))
        artist=ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=fs,
                       fontweight='bold' if bold else 'normal',linespacing=1.25)
        boxdata.append((artist,(x,y,w,h)))
    def arrow(start,end):
        ax.annotate('',xy=end,xytext=start,arrowprops={'arrowstyle':'-|>','lw':1.15,'color':'#617D8B','shrinkA':0,'shrinkB':2})
    box(.02,.85,.96,.13,'Reconstructed source: 233 patients / 377 CT–mask records\nEarlier available subset: 195 patients / 318 pairs\nArchive recovery: +38 patients / +59 pairs',fs=7.6)
    box(.02,.70,.60,.10,'Source-labelled arterial-phase candidates\n126 patients / 230 segmentation records',fs=8.0,bold=True)
    box(.69,.70,.29,.10,'Other source phases\n107 patients / 147 records',fill='#F3F3F3',fs=7.4)
    arrow((.32,.85),(.32,.80));arrow((.835,.85),(.835,.80))
    ax.text(.02,.625,'a  Patient eligibility',fontsize=9.0,fontweight='bold',ha='left',va='center')
    box(.02,.43,.24,.10,'126 candidate\npatients',fs=8.0)
    box(.34,.43,.24,.10,'118 included patients\n83 deaths / 35 censored',fill='#E7F0F4',fs=7.4,bold=True)
    box(.67,.37,.31,.16,'8 patients excluded\n1 delayed phase\n2 indeterminate phase\n5 without a valid ROI',fill='#F3F3F3',fs=7.4)
    arrow((.26,.48),(.34,.48))
    ax.plot([.14,.14,.825],[.53,.57,.57],color='#617D8B',lw=1.15)
    arrow((.825,.57),(.825,.53))
    ax.text(.02,.300,'b  Segmentation-record eligibility',fontsize=9.0,fontweight='bold',ha='left',va='center')
    box(.02,.11,.24,.10,'230 source\nsegmentation records',fs=7.5)
    box(.34,.11,.24,.10,'202 active records\n118 included patients',fill='#E7F0F4',fs=7.4,bold=True)
    box(.67,.065,.31,.145,'28 unused source records\n16 from the 8 excluded patients\n12 from patients still included',fill='#F3F3F3',fs=7.0)
    arrow((.26,.16),(.34,.16))
    ax.plot([.14,.14,.825],[.21,.245,.245],color='#617D8B',lw=1.15)
    arrow((.825,.245),(.825,.21))
    ax.text(.02,.010,'Segmentation records do not redefine the source clinical lesion counts.',fontsize=6.9,ha='left',va='bottom')
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer();bounds=[]
    for text,b in boxdata:
        bbox=text.get_window_extent(renderer).transformed(ax.transAxes.inverted())
        x,y,w,h=b;ok=bbox.x0>=x+.006 and bbox.x1<=x+w-.006 and bbox.y0>=y+.004 and bbox.y1<=y+h-.004
        assert ok, text.get_text()
        bounds.append({'text':text.get_text(),'within_box':bool(ok)})
    for ext in ['svg','pdf','png','tiff']:
        fig.savefig(out/('Figure_1.'+ext),dpi=600,facecolor='white',**({'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}))
    plt.close(fig)
    with Image.open(out/'Figure_1.png') as a,Image.open(out/'Figure_1.tiff') as b:
        assert a.size==b.size and ImageChops.difference(a.convert('RGB'),b.convert('RGB')).getbbox() is None
        pixels=list(a.size);dpi=list(a.info.get('dpi',()))
    assert '<text' in (out/'Figure_1.svg').read_text(encoding='utf-8')
    assert before=={str(p):sha(p) for p in inputs.values()}
    legend=('Source reconstruction and final patient and segmentation-record eligibility. '
            'The earlier available subset of 195 patients and 318 CT–mask pairs was expanded by recovery of 38 patients and 59 pairs, yielding 233 patients and 377 source records. '
            'Of these, 126 patients with 230 source-labelled arterial-phase segmentation records formed the candidate pool; the remaining 107 patients and 147 records had other source phase labels. '
            '(a) Eight candidate patients were excluded: one with a delayed phase, two with an indeterminate phase, and five without a valid ROI; 118 patients remained (83 deaths and 35 censored). '
            '(b) Of the 230 candidate source records, 202 were active and 28 were unused: 16 belonged to the eight excluded patients and 12 belonged to patients who remained included. '
            'Archive recovery is distinct from medical eligibility, and segmentation-record counts do not redefine source clinical lesion counts. ROI, region of interest.')
    (out/'Figure_1_legend.txt').write_text(legend+'\n',encoding='utf-8')
    audit={'status':'PASS_SOURCE_COUNTS_EXPORTS_AND_TEXT_BOUNDS_PENDING_VISUAL_REVIEW',
           'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
           'backend':'Python matplotlib','python':sys.version,'command_argv':sys.argv,
           'run_root':str(run),'counts':dict(counts),'excluded_patients':excluded_rows,
           'source_hashes_before':before,'source_hashes_after':{str(p):sha(p) for p in inputs.values()},
           'source_files_unchanged':True,'text_box_checks':bounds,
           'png_tiff_identical':True,'png_pixels':pixels,'png_dpi':dpi,'svg_editable_text':True,'figsize_inches':[7.09,4.8],
           'new_clinical_judgments':False,'new_model_runs':0,
           'exports':[{'path':str(out/('Figure_1.'+ext)),'sha256':sha(out/('Figure_1.'+ext)),'bytes':(out/('Figure_1.'+ext)).stat().st_size} for ext in ['svg','pdf','png','tiff']]}
    (out/'Figure_1_count_and_export_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({'status':audit['status'],'counts':dict(counts),'directory':str(out)},ensure_ascii=False))

if __name__=='__main__':main()
