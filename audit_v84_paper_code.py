"""Numerically test attention-mask semantics in the pinned author release only."""
import ast,types,json
from pathlib import Path
import torch
from torch import nn
from torch.nn import functional as F
from v57_pie import write,sha
root=Path.cwd();source=root/'experiments/v83_context_expansion/analysis/author_code/Finetune_perturbseq__finetune_perturbseq.py';tree=ast.parse(source.read_text());names={'FlashAttentionEncoderLayer','FlashAttentionDecoderLayer'};nodes=[n for n in tree.body if isinstance(n,ast.ClassDef) and n.name in names];assert len(nodes)==2;module=types.ModuleType('audited_author_attention');module.__dict__.update(nn=nn,F=F,torch=torch);exec(compile(ast.Module(body=nodes,type_ignores=[]),str(source),'exec'),module.__dict__)
torch.manual_seed(84);x=torch.randn(1,4,8);zero=torch.zeros((1,4),dtype=torch.bool);enc=module.FlashAttentionEncoderLayer(8,2,16,0).eval();dec=module.FlashAttentionDecoderLayer(8,2,16,0).eval()
with torch.inference_mode():
 y=enc(x,src_key_padding_mask=zero);changed=x.clone();changed[:,0]+=torch.randn(8);y2=enc(changed,src_key_padding_mask=zero);off=float((y[:,1:]-y2[:,1:]).abs().max());assert off==0
 m1=torch.randn_like(x);m2=torch.randn_like(x);z1=dec(x,m1,tgt_key_padding_mask=zero,memory_key_padding_mask=zero);z2=dec(x,m2,tgt_key_padding_mask=zero,memory_key_padding_mask=zero);memory_effect=float((z1-z2).abs().max());assert memory_effect==0
 z3=dec(changed,m1,tgt_key_padding_mask=zero,memory_key_padding_mask=zero);forward_effect=float((z1[:,1:]-z3[:,1:]).abs().max());assert forward_effect==0
 later=x.clone();later[:,-1]+=torch.randn(8);z4=dec(later,m1,tgt_key_padding_mask=zero,memory_key_padding_mask=zero);backward_effect=float((z1[:,:-1]-z4[:,:-1]).abs().max());assert backward_effect>0
record=dict(pinned_repository_commit='8996c7ca8c336d636de4ba06a4d5d97b2bbe2ce8',author_source_sha256=sha(source),torch_version=torch.__version__,tests_passed=True,encoder_off_token_effect_with_all_false_padding=off,decoder_memory_effect_with_all_false_padding=memory_effect,first_token_to_later_decoder_token_effect=forward_effect,last_token_to_earlier_decoder_token_effect=backward_effect,explanation='The released code applies masked_fill(-inf) to boolean masks, leaving boolean dtype. PyTorch SDPA treats True as allowed. Thus all-false encoder/cross masks block all attention; upper-triangular decoder masks allow later tokens. This verifies the inspected branch behavior, not the architecture used to train published checkpoints or the validity of the paper.',missing_architecture_configuration='Training parser never declares embedding_dim,nhead,num_encoder_layer,num_decoder_layer,dropout; checked DeepSpeed configs do not define them either.',checkpoint_executed=False,claim_about_paper_benchmark=False,code_sha256=sha(__file__))
write(root/'experiments/v84_intervention_operator/analysis/scdebart_code_audit.json',record);print(json.dumps(record,indent=2))
