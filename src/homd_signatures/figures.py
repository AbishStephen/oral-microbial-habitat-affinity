"""Rebuild the agreed five main figures from recovered plotting code and frozen inputs."""
from pathlib import Path
import argparse, json, os, shutil, subprocess, sys, time
from .reproduce import REPO

def validate_export(path):
    """Validate one immutable byte snapshot of a raster/vector export."""
    import io
    data=path.read_bytes()
    if path.suffix == '.pdf':
        if not (data.startswith(b'%PDF-') and data.rstrip().endswith(b'%%EOF')):
            raise ValueError(f'Incomplete PDF export: {path.name}')
    elif path.suffix == '.svg':
        import xml.etree.ElementTree as ET
        ET.fromstring(data)
    elif path.suffix == '.png':
        from PIL import Image
        with Image.open(io.BytesIO(data)) as im:im.verify()
        with Image.open(io.BytesIO(data)) as im:im.load()
    return data

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=REPO/'build/main_figures');a=ap.parse_args()
    a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=True);report=[];exports={}
    for num in range(1,6):
        target=a.output/f'Figure{num}';shutil.copytree(REPO/'figure_sources'/f'Figure{num}',target,dirs_exist_ok=True)
        script=target/f'build_figure{num}.py'
        start=time.time();env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',MPLBACKEND='Agg')
        with (target/'rebuild.log').open('w') as f:subprocess.run([sys.executable,str(script)],cwd=target,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
        for p in target.rglob('*'):
            if p.suffix in ['.png','.pdf','.svg']:exports[str(p.relative_to(a.output))]=validate_export(p)
        outputs=[str(p.relative_to(a.output)) for p in target.rglob('*') if p.suffix in ['.png','.pdf','.svg']]
        assert any(p.endswith('.pdf') for p in outputs) and any(p.endswith('.svg') for p in outputs)
        report.append(dict(figure=num,seconds=round(time.time()-start,2),outputs=outputs));print(json.dumps(report[-1]),flush=True)
    import zipfile
    with zipfile.ZipFile(a.output/'verified_figure_exports.zip','w',zipfile.ZIP_DEFLATED) as z:
        for name,data in exports.items():z.writestr(name,data)
    (a.output/'figure_rebuild_report.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
