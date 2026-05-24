"""Execute cells in a fresh IPython process without sockets.
Usage: python run_inprocess.py
The captured outputs are actual execution results; failure stops the run.
"""
import os,time,traceback,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"production"))
import runtime
if __name__ == "__main__": runtime.enter_cli()
from pathlib import Path
import nbformat
from IPython.core.interactiveshell import InteractiveShell
from IPython.utils.capture import capture_output

def main():
    started=time.perf_counter()
    root=Path(__file__).resolve().parent;os.chdir(root)
    path=root/'spin_reviewer_verification.ipynb';book=nbformat.read(path,as_version=4)
    shell=InteractiveShell.instance()
    count=0
    for c in book.cells:
        if c.cell_type=='code':c.outputs=[];c.execution_count=None
    book.metadata['execution_method']='Fresh-process IPython in-process execution; no socket kernel'
    for i,c in enumerate(book.cells):
        if c.cell_type!='code':continue
        count+=1;tic=time.perf_counter()
        with capture_output() as captured:
            result=shell.run_cell(c.source,store_history=True)
        c.execution_count=count;c.outputs=[]
        if captured.stdout:c.outputs.append(nbformat.v4.new_output('stream',name='stdout',text=captured.stdout))
        if captured.stderr:c.outputs.append(nbformat.v4.new_output('stream',name='stderr',text=captured.stderr))
        for o in captured.outputs:c.outputs.append(nbformat.v4.new_output('display_data',data=o.data,metadata=o.metadata))
        err=result.error_before_exec or result.error_in_exec
        if err:c.outputs.append(nbformat.v4.new_output('error',ename=type(err).__name__,evalue=str(err),traceback=traceback.format_exception(type(err),err,err.__traceback__)))
        nbformat.write(book,path)
        print('Cell %d, execution %d, %.1fs'%(i,count,time.perf_counter()-tic),flush=True)
        if err:raise RuntimeError('Cell %d failed: %s'%(i,err)) from err
    book.metadata['execution_seconds']=round(time.perf_counter()-started,3)
    nbformat.write(book,path)
    print('Notebook execution complete.',flush=True)
if __name__=='__main__':main()
