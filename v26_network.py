"""A full-axis conditional logistic-normal multinomial VAE for V26."""
import torch
from torch import nn
class Decoder(nn.Module):
 def __init__(self,inputs,genes):
  super().__init__();self.q1=nn.Linear(inputs+24,128);self.q2=nn.Linear(128,32);self.d1=nn.Linear(40,64);self.d2=nn.Linear(64,128);self.out=nn.Linear(128,genes)
 def posterior(self,x,s,c):
  v=self.q2(torch.nn.functional.silu(self.q1(torch.cat([x,s,c],1))));mu,lv=v.chunk(2,1);return mu,lv.clamp(-6,3)
 def residual(self,z,s,c):
  v=torch.nn.functional.silu(self.d1(torch.cat([z,s,c],1)));v=torch.nn.functional.silu(self.d2(v));return 4*torch.tanh(self.out(v)/4)
 def forward(self,x,s,c,base,mask):
  mu,lv=self.posterior(x,s,c);z=mu+torch.exp(.5*lv)*torch.randn_like(mu);delta=self.residual(z,s,c)*mask;logits=torch.log(base.clamp_min(1e-12))+delta;logits=logits.masked_fill(base<=0,-1e9);lp=torch.log_softmax(logits,1);kl=.5*(mu.square()+lv.exp()-lv-1).sum(1);return lp,kl
