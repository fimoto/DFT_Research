import contextlib,io,json,sys,time,traceback,platform
import nbformat
src,dst=sys.argv[1:3]; n=nbformat.read(src,as_version=4)
ns={'__name__':'__main__'}; start=time.time(); count=0
for c in n.cells:
 if c.cell_type=='code': c.outputs=[];c.execution_count=None
try:
 for i,c in enumerate(n.cells):
  if c.cell_type!='code':continue
  count+=1;c.execution_count=count;buf=io.StringIO()
  print('CELL',i,round(time.time()-start,1),flush=True)
  try:
   with contextlib.redirect_stdout(buf),contextlib.redirect_stderr(buf):exec(compile(c.source,f'cell_{i}','exec'),ns)
  except Exception as e:
   c.outputs.append(nbformat.v4.new_output('error',ename=type(e).__name__,evalue=str(e),traceback=traceback.format_exc().splitlines()))
   raise
  finally:
   if buf.getvalue():c.outputs.insert(0,nbformat.v4.new_output('stream',name='stdout',text=buf.getvalue()))
   nbformat.write(n,dst)
 n.metadata['execution_method']='Fresh Python process; code cells executed in document order with captured stdout/stderr (no live Jupyter socket in this environment).'
 n.metadata['execution_seconds']=round(time.time()-start,3)
finally:nbformat.write(n,dst)
print('Failed labels:',[label for label,ok in ns.get('RESULTS',[]) if not ok],flush=True)
print('DONE',round(time.time()-start,1),flush=True)
