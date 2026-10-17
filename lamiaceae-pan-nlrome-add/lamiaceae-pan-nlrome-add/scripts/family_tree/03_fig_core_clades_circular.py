#!/usr/bin/env python3
"""Circular NB-ARC clade figure. Usage: 03_fig_core_clades_circular.py PARTITION_DIR CLADE_CONSERVATION.tsv CLADE_GENES.tsv OUT_PREFIX [midpoint|CladeNN]
PARTITION_DIR holds *_collapsed.treefile, *_clade_info.csv, *_anchors.tsv from 01_partition_refs.py (prefix "refs").
"""
import csv, numpy as np, matplotlib, pandas as pd
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib import cm
from Bio import Phylo
import sys; F=sys.argv[1].rstrip("/")+"/"; CONS=sys.argv[2]; GENES=sys.argv[3]; OUT=sys.argv[4]
ROOT=sys.argv[5] if len(sys.argv)>5 else "midpoint"  # "midpoint" or a clade id to use as outgroup, e.g. Clade18
t=Phylo.read(F+"refs_collapsed.treefile","newick"); (t.root_at_midpoint() if ROOT=="midpoint" else t.root_with_outgroup({"name":ROOT})); t.ladderize()
info={r["clade_id"]:r for r in csv.DictReader(open(F+"refs_clade_info.csv"))}
cons=pd.read_csv(CONS,sep="\t").set_index("clade")
exp=pd.read_csv(GENES,sep="\t").set_index("clade")
COL=dict(lam="#0F8F7E",ath="#C9184A",cnl="#0072B2",tnl="#D17A00",rnl="#7A5195")
INK="#222"; MUTED="#6b6b6b"; BG="#fcfcfb"
NOTE={"Clade18":"ADR1/NRG1","Clade59":"ZAR1 refs","Clade10":"RPP13","Clade06":"MLA","Clade11":"RPS2/RPS5","Clade14":"Ptr1","Clade27":"R3a/I2","Clade25":"Mi-1/Prf","Clade57":"NRC1-4","Clade71":"NRC0"}
leaves=t.get_terminals(); N=len(leaves)
GAP=0.14; step=(2*np.pi-GAP)/N
th={l:GAP/2+step*(i+0.5) for i,l in enumerate(leaves)}
for c in t.get_nonterminals(order="postorder"):
    th[c]=np.mean([th[k] for k in c.clades])
x={t.root:0}
for c in t.get_nonterminals(order="preorder"):
    for k in c.clades: x[k]=x[c]+(k.branch_length or 0)
xm=max(x[l] for l in leaves); RT=0.36
r=lambda n: RT*x[n]/xm
def core(cid): return cid in cons.index and cons.loc[cid,"n_species"]>=35 and cons.loc[cid,"copies_per_species"]<=3
def stable(cid): return core(cid) and cons.loc[cid,"median_copies"]<=1   # single-copy-like core (tier 1); other core clades = widespread, locally expanded (tier 2)
fig=plt.figure(figsize=(12.5,12.5),facecolor=BG)
ax=fig.add_axes([0.0,0.06,1.0,0.88],projection="polar",facecolor=BG)
ax.set_theta_zero_location("N"); ax.set_theta_direction(-1); ax.set_ylim(0,0.80); ax.axis("off")
R1a,R1b=0.40,0.455   # composition ring
R2a,R2b=0.465,0.505  # species ring
R3a=0.515; R3max=0.17 # genes bars
RL=R3a+R3max+0.015
# core wedges
for l in leaves:
    if core(l.name):
        ax.bar(th[l],RL-R1a+0.115,width=step,bottom=R1a-0.005,color="#CBD2D9" if stable(l.name) else "#ECEEF1",zorder=0,linewidth=0)
# tree
for c in t.get_nonterminals():
    a=[th[k] for k in c.clades]; rr=r(c)
    ax.plot(np.linspace(min(a),max(a),60),[rr]*60,c="#555",lw=.7,zorder=2)
    for k in c.clades: ax.plot([th[k]]*2,[rr,r(k)],c="#555",lw=.7,zorder=2)
spcm=plt.get_cmap("Greys")
for l in leaves:
    cid=l.name; ri=info[cid]; T=th[l]
    ax.plot([T]*2,[r(l),R1a-0.004],c="#cfcfcf",lw=.35,ls=":",zorder=1)
    rc=ri["ref_classes"]; rcd=dict(p.split(":") for p in rc.split(";")) if rc else {}
    parts=[("lam",int(ri["n_lamiaceae"])),("ath",int(ri["n_arabidopsis"])),("cnl",int(rcd.get("CNL",0))),("tnl",int(rcd.get("TNL",0))),("rnl",int(rcd.get("RNL",0)))]
    tot=sum(n for _,n in parts); bot=R1a
    for k,n in parts:
        if n:
            h=(R1b-R1a)*n/tot
            ax.bar(T,h,width=step*0.88,bottom=bot,color=COL[k],linewidth=0,zorder=3); bot+=h
    if cid in cons.index:
        ns=cons.loc[cid,"n_species"]
        ax.bar(T,R2b-R2a,width=step*0.88,bottom=R2a,color=spcm(0.12+0.88*ns/46),linewidth=0,zorder=3)
        g=exp.loc[cid,"lamiaceae_genes_expanded"]
        h=R3max*np.log10(g+1)/np.log10(900)
        ax.bar(T,h,width=step*0.6,bottom=R3a,color="#222" if core(cid) else "#8a8a8a",linewidth=0,zorder=3)
    else:
        ax.bar(T,R2b-R2a,width=step*0.88,bottom=R2a,color="none",edgecolor="#d5d5d5",linewidth=.4,zorder=3)
    deg=np.degrees(T); rot=90-deg; left=np.sin(T)<0
    if left: rot+=180
    lab=cid.replace("Clade","C")+(f"  {NOTE[cid]}" if cid in NOTE else "")
    if left: lab=(f"{NOTE[cid]}  " if cid in NOTE else "")+cid.replace("Clade","C")
    ax.text(T,RL,lab,rotation=rot,rotation_mode="anchor",ha="right" if left else "left",va="center",fontsize=6.3,
            fontweight="bold" if core(cid) else "normal",color=INK,zorder=4)
# guide rings for gene scale
for g in (10,100,900):
    rr=R3a+R3max*np.log10(g+1)/np.log10(900)
    ax.plot(np.linspace(GAP/2,2*np.pi-GAP/2,200),[rr]*200,c="#bbb",lw=.4,ls="--",zorder=1)
    ax.text(GAP/2*0.5,rr,str(g),fontsize=5.5,color=MUTED,ha="center",va="center",rotation=0)
fig.suptitle("NB-ARC family tree of Lamiaceae NLRs anchored with RefPlantNLR (85 clades)",fontsize=12,y=0.975)
H=[Patch(color=COL["lam"],label="Lamiaceae"),Patch(color=COL["ath"],label="Arabidopsis"),Patch(color=COL["cnl"],label="reference CNL"),Patch(color=COL["tnl"],label="reference TNL"),Patch(color=COL["rnl"],label="reference RNL"),Patch(color="#CBD2D9",label="conserved, single-copy-like"),Patch(color="#ECEEF1",label="widespread, locally expanded")]
fig.legend(handles=H,loc="lower center",bbox_to_anchor=(0.5,0.035),ncol=7,fontsize=7.6,frameon=False)
fig.text(0.5,0.012,"Rings, inside to outside: composition of tree leaves per clade (share of leaves); Lamiaceae species present (grey scale, 1–46); Lamiaceae genes per clade (log scale, bars). "
         "Shaded clades: ≥35 species and ≤3 genes/species on average; darker = median 1 copy per species (single-copy-like), lighter = widespread with local expansions.\ncd-hit (85%) representatives with clusters expanded to genes. Root: "+("midpoint (display only)" if ROOT=="midpoint" else "between the TNL/RNL clade ("+ROOT.replace("Clade","C")+", with the Klebsormidium sequence) and the CNL clades")+".",
         ha="center",fontsize=6.8,color=MUTED)
cax=fig.add_axes([0.80,0.945,0.15,0.012]); sm=cm.ScalarMappable(cmap=spcm,norm=plt.Normalize(1,46)); cb=fig.colorbar(sm,cax=cax,orientation="horizontal"); cb.ax.tick_params(labelsize=6,length=2); cb.set_label("Lamiaceae species (of 46)",fontsize=6.5,labelpad=1)
for e in ("png","pdf"): fig.savefig(f"{OUT}.{e}",dpi=200,facecolor=BG)
print("ok")
