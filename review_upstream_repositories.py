"""Record both upstream histories and executable-code changes without modifying the clones."""
from pathlib import Path
import ast
import hashlib
import json
import subprocess
import pandas as pd
from build_submission import write_json, sha256
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'analysis/upstream_review'
REPOS=['Virtual-Cell-Challenge-2026','arc-virtual-cell-atlas']


def git(repo,*args):
    safe=str((ROOT/repo).resolve()).replace('\\','/')
    return subprocess.run(['git','-c','safe.directory='+safe,'-C',str(ROOT/repo),*args],
                          capture_output=True,text=True,encoding='utf-8',check=True).stdout


class WithoutDocs(ast.NodeTransformer):
    def generic_visit(self,node):
        node=super().generic_visit(node)
        if isinstance(node,(ast.Module,ast.FunctionDef,ast.AsyncFunctionDef,ast.ClassDef)) and node.body:
            first=node.body[0]
            if isinstance(first,ast.Expr) and isinstance(first.value,ast.Constant) and isinstance(first.value.value,str):
                node.body=node.body[1:]
        return node


def executable_digest(source):
    return hashlib.sha256(ast.dump(WithoutDocs().visit(ast.parse(source)),include_attributes=False).encode()).hexdigest()


def main():
    OUT.mkdir(parents=True,exist_ok=True);summaries=[]
    for repo in REPOS:
        assert not git(repo,'status','--porcelain').strip(), 'Review requires a clean upstream clone'
        records=[]
        log=git(repo,'log','--all','--reverse','--format=%H|%aI|%s')
        for line in log.splitlines():
            commit,date,message=line.split('|',2)
            paths=git(repo,'diff-tree','--root','--no-commit-id','--name-only','-r',commit).splitlines()
            parent=git(repo,'rev-list','--parents','-n','1',commit).split()[1:]
            changed=[]
            for path in paths:
                if not path.endswith('.py'):continue
                try:after=git(repo,'show',commit+':'+path)
                except subprocess.CalledProcessError:after=None
                try:before=git(repo,'show',parent[0]+':'+path) if parent else None
                except subprocess.CalledProcessError:before=None
                if after is None or before is None or executable_digest(after)!=executable_digest(before):changed.append(path)
            records.append({'commit':commit,'author_date':date,'subject':message,'changed_paths':json.dumps(paths),
                            'changed_python_executable_paths':json.dumps(changed),'python_executable_changed':bool(changed)})
        pd.DataFrame(records).to_csv(OUT/(repo+'_commits.csv'),index=False)
        files=[]
        for rel in git(repo,'ls-files').splitlines():
            path=ROOT/repo/rel
            record={'path':rel,'bytes':path.stat().st_size,'sha256':sha256(path)}
            if rel.endswith('.ipynb'):
                notebook=json.loads(path.read_text(encoding='utf-8'));cells=[c for c in notebook['cells'] if c['cell_type']=='code']
                record['code_cells']=len(cells);record['code_sha256']=hashlib.sha256(''.join(''.join(c['source']) for c in cells).encode()).hexdigest()
            files.append(record)
        summary={'repository':repo,'origin':git(repo,'remote','get-url','origin').strip(),'head':git(repo,'rev-parse','HEAD').strip(),
                 'origin_main':git(repo,'rev-parse','origin/main').strip(),'commits_all_refs':len(records),
                 'commits_main':int(git(repo,'rev-list','--count','origin/main').strip()),
                 'python_executable_change_commits':sum(r['python_executable_changed'] for r in records),
                 'tracked_files':files,'branches':git(repo,'branch','-a').splitlines(),'working_tree_clean':True}
        if repo=='arc-virtual-cell-atlas':
            summary['egress_guidance_fix_branch']=git(repo,'branch','-a','--contains','8fc7a9e').splitlines()
            summary['egress_guidance_fix_merged_into_main']=subprocess.run(['git','-c','safe.directory='+str((ROOT/repo).resolve()).replace('\\','/'),'-C',str(ROOT/repo),'merge-base','--is-ancestor','8fc7a9e','origin/main'],capture_output=True).returncode==0
        summaries.append(summary)
        print(repo,summary['head'],len(records),'commits across refs;',summary['python_executable_change_commits'],'executable Python change commits',flush=True)
    write_json(OUT/'repository_review.json',{'script_sha256':sha256(Path(__file__)),'repositories':summaries,
        'scope':'All reachable commit metadata and file changes; current tracked files and notebook code hashes; AST comparison ignores docstrings. Notebook tutorial diffs and selected historical methods reviewed separately.'})


if __name__=='__main__':main()
